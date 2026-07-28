# Riggs Engine — Architecture

> **Status**: Phase 1 (Gillespie SSA Core) — in progress.
> This document describes the system as it currently exists, not as planned.
> Updated: 2026-07-28

## Overview

The Riggs Engine (codename for the Monod Engine) simulates **epigenetic methylation
as a volatile biological memory system** using the Gillespie Stochastic Simulation
Algorithm (SSA). DNA methylation state is treated as switchable RAM bits layered on
an immutable DNA sequence (ROM). The engine never mutates base pairs — only the
methylation state vector changes.

## Layer Architecture

```
Phase 4: Distributed Workers (gRPC + failure detector)   [not started]
Phase 3: TUI (Bubble Tea dashboard)                      [not started]
Phase 2: Biological Instruction Set                      [not started]
Phase 1: Gillespie SSA Core                              [<-- YOU ARE HERE]
```

**Phase 1** is the pure math tier. It implements the Gillespie SSA Direct Method
as a general-purpose stochastic reaction simulator. It knows nothing about biology —
`Reaction`, `State`, and `Run()` operate on abstract species counts and rate constants.
Phase 2 will add biological semantics (genome, dCas9, enzymes) as reaction *generators*
that produce `Reaction` values for the Phase 1 core to simulate.

## Phase 1 Components

### Package: `gillespie`

The raw stochastic simulation engine. Contains no biology-specific logic.

| File | Responsibility |
|------|---------------|
| `types.go` | Core data types: `Reaction`, `State`, `TrajectoryRecord`, `EnsembleResult` |
| `ssa.go` | SSA Direct Method: `Step()` (single event) and `Run()` (full trajectory) |
| `ensemble.go` | Worker-pool concurrent ensemble runner with deterministic seeding |
| `ssa_test.go` | Table-driven tests including statistical validation against analytical CTMC results |
| `ssa_bench_test.go` | Allocation and throughput benchmarks for the inner loop |

### Package: `cmd/riggs`

Minimal CLI entrypoint. Runs a toy single-CpG-site system and prints statistics.

## Key Design Decisions

### Flat slices, not maps or interfaces

Species counts are `[]int` indexed by species ID. Reaction propensities are `func`
fields, not interface methods. Both choices eliminate allocation and dispatch overhead
in the inner SSA loop, which matters on a 4-core / 8GB target with no GPU.

### Pre-allocated scratch buffers

The `Step()` function accepts a pre-allocated `[]float64` scratch buffer for propensity
computation. The ensemble runner allocates one per worker at startup. Target:
**0 allocs/op** in the inner loop.

### Pluggable, seedable RNG

Each trajectory (and each worker) receives its own `*rand.Rand` with a deterministic
seed derived from a base seed + trajectory index. No shared RNG, no mutex contention,
fully reproducible output.

### Worker-pool, not goroutine-per-trajectory

Ensemble simulation uses a fixed-size worker pool (default: `GOMAXPROCS`) consuming
trajectory jobs from a channel. Prevents unbounded goroutine spawn and controls
memory footprint to `O(workers)`, not `O(trajectories)`.

## Dependencies

Phase 1 uses **only Go standard library** — no external dependencies.

| stdlib package | Usage |
|---------------|-------|
| `math` | `Log`, `Inf`, `Sqrt` |
| `math/rand/v2` | Per-trajectory seedable PRNG |
| `sync` | `WaitGroup` for worker pool |
| `runtime` | `GOMAXPROCS` for default pool size |
| `fmt`, `flag` | CLI argument parsing and output |
