# Riggs Engine — Architecture

> **Status**: Phase 3 (TUI Dashboard) — complete.
> This document describes the system as it currently exists, not as planned.
> Updated: 2026-08-06

## Overview

The Riggs Engine (codename for the Monod Engine) simulates **epigenetic methylation
as a volatile biological memory system** using the Gillespie Stochastic Simulation
Algorithm (SSA). DNA methylation state is treated as switchable RAM bits layered on
an immutable DNA sequence (ROM). The engine never mutates base pairs — only the
methylation state vector changes.

## Layer Architecture

```
Phase 4: Distributed Workers (gRPC + failure detector)   [not started]
Phase 3: TUI (Bubble Tea dashboard)                      [<-- COMPLETE]
Phase 2: Biological Instruction Set                      [COMPLETE]
Phase 1: Gillespie SSA Core                              [COMPLETE]
```

**Phase 1** is the pure math tier. It implements the Gillespie SSA Direct Method
as a general-purpose stochastic reaction simulator. It knows nothing about biology —
`Reaction`, `State`, and `Run()` operate on abstract species counts and rate constants.

**Phase 2** adds biological semantics. The `bio` package maps real biological entities
(genomes, CpG sites, dCas9 complexes, enzymes, environmental triggers) onto Phase 1
types. The core pattern is **reaction generation**: `System.Build()` produces
`[]Reaction` + `State` + `[]ScheduledEvent` values that the Phase 1 engine simulates.

**Phase 3** provides real-time visualization. The `tui` package bridges the simulation
goroutine to a Bubble Tea terminal dashboard via tick-driven coalescing. The simulation
runs at thousands of events/sec; the TUI reads a coalesced snapshot at 10 FPS.

## Package Architecture

```
cmd/riggs/main.go          Presentation Tier — CLI + TUI entrypoint
        │
        ├──→ (genome/toy)     Direct ensemble runner output
        │
        └──→ (tui)            Bubble Tea program
                │
                ▼
tui/                        Visualization Tier — terminal dashboard
  bridge.go                   SimRunner, SimSnapshot, RingBuffer, WALEntry
  model.go                    Bubble Tea Model (Init/Update/View)
  views.go                    renderMemoryMap, renderWALTail, renderStats
        │
        ▼
bio/                        Logic Tier — biology-specific types
  genome.go                   Genome, Site, SiteContext (the "disk")
  complex.go                  TargetingComplex (dCas9 pointer)
  trigger.go                  EnvironmentalTrigger (WAL mechanic)
  system.go                   System.Build() → reactions + state + schedule
        │
        ▼
gillespie/                  Math Tier — pure SSA engine
  types.go                    Reaction, State, TrajectoryRecord, ScheduledEvent
  ssa.go                      Step(), Run(), RunWithSchedule(), TimeWeightedFraction()
  ensemble.go                 RunEnsemble() — worker-pool + Welford aggregation
```

### Data flow

```
bio.System.Build()
    │
    ├──→ []gillespie.Reaction       (2N + K stochastic channels)
    ├──→ gillespie.State            (N + K species: methylation + bound states)
    └──→ []gillespie.ScheduledEvent (deterministic trigger events)
            │
            ▼
    gillespie.RunWithSchedule()     (hybrid SSA: stochastic + deterministic)
            │
            ▼
    gillespie.TrajectoryRecord      (sparse event stream)
            │
            ▼
    gillespie.TimeWeightedFraction() → per-trajectory scalar
            │
            ▼
    gillespie.RunEnsemble()         (Welford aggregation → EnsembleResult)
```

## Phase 1 Components

### Package: `gillespie`

The raw stochastic simulation engine. Contains no biology-specific logic.

| File | Responsibility |
|------|---------------|
| `types.go` | Core data types: `Reaction`, `State`, `TrajectoryRecord`, `EnsembleResult`, `ScheduledEvent` |
| `ssa.go` | SSA Direct Method: `Step()`, `Run()`, `RunWithSchedule()` (hybrid SSA), `TimeWeightedFraction()` |
| `ensemble.go` | Worker-pool concurrent ensemble runner with deterministic seeding, optional `Schedule` support |
| `ssa_test.go` | Table-driven tests including statistical validation against analytical CTMC results |
| `ssa_bench_test.go` | Allocation and throughput benchmarks for the inner loop |

## Phase 2 Components

### Package: `bio`

Biological instruction set — reaction generators over the Phase 1 core.

| File | Responsibility |
|------|---------------|
| `genome.go` | `Genome`, `Site`, `SiteContext` — the immutable DNA layer |
| `complex.go` | `TargetingComplex` — dCas9-guide RNA pointer with enhanced write/erase |
| `trigger.go` | `EnvironmentalTrigger` — deterministic stimulus that binds complexes |
| `system.go` | `System.Build()` — maps biology types → Gillespie reactions/state/schedule |
| `bio_test.go` | Statistical validation: targeted vs background, trigger timing, unbinding |
| `bio_bench_test.go` | Scaling benchmarks at N=10/100/1000 sites |

### Package: `cmd/riggs`

CLI entrypoint with two modes:
- `toy`: Phase 1 single-site, two-state CTMC validation
- `genome`: Phase 2 multi-site genome with targeted methylation and triggers

## Key Design Decisions

### Flat slices, not maps or interfaces

Species counts are `[]int` indexed by species ID. Reaction propensities are `func`
fields, not interface methods. Both choices eliminate allocation and dispatch overhead
in the inner SSA loop, which matters on a 4-core / 8GB target with no GPU.

### Pre-allocated scratch buffers

The `Step()` function accepts a pre-allocated `[]float64` scratch buffer for propensity
computation. The ensemble runner allocates one per worker at startup. Target:
**0 allocs/op** in the inner loop.

### Reaction generation, not hardcoded channels

The `bio` package generates `[]Reaction` programmatically — it doesn't register
reactions one-by-one. `System.Build()` loops over sites and complexes to produce
the full reaction set in one pass. The propensity closures capture only immutable
data (site indices, rate constants, pre-computed complex-target mappings).

### Conditional propensity gating (not separate channels)

Enhanced write/erase rates from a bound complex are folded into the *same* reaction
channel as the background rate. The propensity closure checks the complex's
bound-state species. This avoids channel-count explosion: 2N channels for N sites,
not 4N (background + enhanced × write + erase).

### Hybrid SSA for deterministic triggers

Environmental triggers fire at known times, not stochastically. `RunWithSchedule()`
interleaves deterministic events with SSA steps: at each iteration, compare next
stochastic time τ with next scheduled event time, fire whichever is sooner.
Preserves SSA exactness for stochastic channels.

### Pluggable, seedable RNG

Each trajectory (and each worker) receives its own `*rand.Rand` with a deterministic
seed derived from a base seed + trajectory index. No shared RNG, no mutex contention,
fully reproducible output.

### Worker-pool, not goroutine-per-trajectory

Ensemble simulation uses a fixed-size worker pool (default: `GOMAXPROCS`) consuming
trajectory jobs from a channel. Prevents unbounded goroutine spawn and controls
memory footprint to `O(workers)`, not `O(trajectories)`.

## Phase 3 Components

### Package: `tui`

Terminal dashboard bridging the simulation to Bubble Tea.

| File | Responsibility |
|------|---------------|
| `bridge.go` | `SimRunner` (simulation goroutine), `SimSnapshot` (shared state), `RingBuffer` (bounded WAL history), `WALEntry` |
| `model.go` | Bubble Tea `Model` with `Init`/`Update`/`View`, 10 FPS tick, keybinding dispatch |
| `views.go` | `renderMemoryMap` (genome bitmap), `renderWALTail` (scrolling log), `renderStats` (per-site bars + ensemble reference) |
| `bridge_test.go` | 12 tests: ring buffer correctness, snapshot deep copy, pause/resume/step-once lifecycle, statistical MethFracs validation, WAL population, scheduled event interleaving |

### Phase 3 Design Decisions

#### Tick-driven coalescing (not push-driven)

The simulation fires thousands of events/sec. The TUI renders at 10 FPS. The
bridge uses **tick-driven coalescing**: the simulation goroutine writes to a shared
`SimSnapshot` behind a `sync.Mutex` after every SSA step. The Bubble Tea event
loop reads the snapshot every 100ms (10 FPS tick). Between ticks, all intermediate
states are coalesced — the TUI only sees the latest state.

Why not push individual events as `tea.Msg`? Because:
- The Bubble Tea message queue would grow unboundedly at 10K+ events/sec.
- The terminal can't usefully render faster than ~30 FPS anyway.
- Coalescing is correct because the memory-map view shows *current* state, not history.

#### Backpressure policy

| Data | Policy | Rationale |
|------|--------|-----------|
| `State.Counts` | Last-writer-wins | TUI always sees latest state; intermediate states don't matter for the bitmap |
| Sim time, event count | Last-writer-wins | Monotonically increasing; latest is always correct |
| WAL tail | Ring buffer (256 entries, 8KB) | Fixed RAM bound; view only shows most recent events |
| Per-site methylation fractions | Running integrator (exact) | Simulation goroutine maintains Σ(state × dt) / elapsed; TUI reads the result directly |

The key invariant: **aggregate statistics (MethFracs) are always exact**, even though
the event log view is bounded. The running integrator accumulates (state × duration)
at every SSA step — no events are ever skipped for the integrator, only for the
ring buffer display.

#### Ring buffer for WAL tail

The WAL tail view uses a fixed-size `[256]WALEntry` array (not a slice). When full,
the oldest entry is overwritten. This is a hard RAM constraint — with event rates of
10K+/sec, an unbounded slice would consume megabytes per second. The 256-entry buffer
stores ~5 seconds of history at typical event rates, which is more than enough for
the scrolling `tail -f` view.

#### SimRunner owns the hybrid-SSA loop

The `tui` package inlines the hybrid-SSA interleaving logic (stochastic τ vs. next
scheduled event time) in its own `SimRunner.doOneStep()`. This avoids modifying the
`gillespie/` package — it uses the existing `Step()` propensity computation directly.
The `gillespie/` and `bio/` packages remain completely unchanged in Phase 3.

## Dependencies

The `gillespie/` and `bio/` packages use Go standard library only — **zero external
dependencies**. The `tui/` and `cmd/riggs/` packages add Charm libraries for the TUI.

| Package | Usage |
|---------|-------|
| `math`, `math/rand/v2` | SSA numerics, seedable PRNG (PCG) |
| `sync` | `WaitGroup`, `Mutex`, `Cond` for worker pool and TUI bridge |
| `runtime` | `GOMAXPROCS` for default pool size |
| `sort` | Sorting scheduled events by time |
| `fmt`, `flag`, `os`, `time` | CLI argument parsing, output, process control |
| `github.com/charmbracelet/bubbletea` | TUI framework (Phase 3 only) |
| `github.com/charmbracelet/lipgloss` | Terminal styling (Phase 3 only) |
