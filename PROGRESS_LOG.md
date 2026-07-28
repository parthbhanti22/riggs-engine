# Riggs Engine — Progress Log

> First-person, timestamped entries documenting what was built, why, and what to
> understand as a result. This is the learning journal for the project.

---

## Entry 1 — 2026-07-28T01:30+05:30 — Project Kickoff & Scaffolding

### What was built

Initialized the Go module (`github.com/pxrth9/riggs`) and created the `gillespie`
package scaffolding:

- **`gillespie/types.go`**: Core types for the SSA engine — `Reaction`, `State`,
  `TrajectoryRecord`, `EnsembleResult`.
- **`ARCHITECTURE.md`**: Documents the four-phase layer architecture and Phase 1's
  design decisions.
- **`BIO_MAPPING.md`**: Maps the CpG toy system's biology to Go constructs, including
  the analytical validation target.

### Why these type designs

**`Reaction.Propensity` is a `func` field, not an interface method.** The SSA inner
loop calls `Propensity()` for every reaction channel at every step — potentially
millions of times per trajectory. An interface method would add vtable dispatch
overhead (indirect call through the interface table). A `func` field is a direct
function pointer — the Go compiler can inline small closures and avoid the indirection.
On a CPU-only box where we can't hide latency behind GPU parallelism, this matters.

**`State.Counts` is `[]int`, not `map[string]int`.** Maps have allocation and hashing
overhead per access. For the hot state vector that gets read/written every step, a
flat slice indexed by integer species ID is the right choice. The tradeoff is that
species are identified by index rather than name, so we need the BIO_MAPPING.md
table to keep track of what index 0 means ("methylation state of CpG site").

**`TrajectoryRecord` stores sparse events, not fixed-timestep snapshots.** The SSA
produces events at irregular intervals (the exponential waiting times). Recording
every event is the natural representation — no interpolation needed, no wasted
storage for quiescent periods. For analysis (e.g. time-weighted averages), we
integrate over the inter-event intervals.

### What to understand

The Gillespie algorithm produces exact sample paths of a continuous-time Markov chain.
"Exact" means there is no time-discretization error — unlike Euler or Runge-Kutta
methods for ODEs, or tau-leaping approximations. Each step computes exactly which
reaction fires next and exactly when. The price is that each step advances by only
one reaction event, so wall-clock cost is proportional to the total number of events
in the trajectory. For our toy system (two reaction channels, counts ∈ {0,1}), this
is very fast. For Phase 2 with hundreds of CpG sites and many reaction channels, we
may need tau-leaping (batch multiple reactions per step) or the Gibson-Bruck Next
Reaction Method (maintain a priority queue of next-firing-times to avoid recomputing
all propensities). Both are documented future optimizations — not implemented now.

### Design note: memory footprint

On the 8GB target machine (3.7GB visible to WSL), the main memory concern is
per-trajectory allocation during ensemble runs. Each trajectory produces a
`TrajectoryRecord` with event timestamps and state snapshots. For the toy system
with `tMax=1000` and an average event rate of ~1.0/time-unit (k_write + k_erase),
a trajectory has ~1000 events × (8 bytes timestamp + 8 bytes count + 8 bytes
reaction index) ≈ 24KB. With 4 concurrent workers, that's ~96KB live. Even 10,000
total trajectories produce only ~240MB of results, well within budget.

The hot-path concern is the propensity scratch buffer in `Step()`. If allocated
per-call, that's a heap allocation every step. The plan is to pre-allocate one
per worker and pass it through — targeting 0 allocs/op in the inner loop.

---

## Open Questions (for review before Phase 2)

- **Environmental trigger schema**: The WAL analogy implies recording *what* caused
  each methylation event (temperature spike? signal molecule?). What schema should
  triggers follow? A simple enum? A structured event with timestamp + type + magnitude?

- **Target reaction-channel count for Phase 2**: How many CpG sites should the
  initial genome model support? 10? 100? 1000? This affects whether we need sparse
  delta vectors and optimized propensity updates (Next Reaction Method) vs. the
  simple linear scan we use now.

- **State persistence**: Should trajectories be streamable to disk during simulation
  (for very long runs), or is in-memory accumulation sufficient for Phase 2?

---

## Entry 2 — 2026-07-28T01:42+05:30 — SSA Core Loop & Numerical Stability

### What was built

`gillespie/ssa.go` — the complete SSA Direct Method implementation:

- **`Step()`**: a single Gillespie step. Computes propensities, draws exponential
  waiting time, selects which reaction fires, applies state change. Accepts a
  pre-allocated scratch buffer for propensity computation — **achieves 0 allocs/op**
  in benchmarks.
- **`Run()`**: drives `Step()` in a loop until `tMax` or absorbing state, recording
  every event into a `TrajectoryRecord`.
- **`TimeWeightedFraction()`**: computes the fraction of time a species was in a
  nonzero state by integrating the piecewise-constant step function defined by the
  event stream.

### The exponential waiting-time derivation (what to understand)

The core mathematical insight of the Gillespie algorithm is that if you have N
reaction channels with propensities a_1, a_2, ..., a_N, then the time until
*any one of them* fires is exponentially distributed with rate a_0 = Σ a_i.

Why? Each channel independently has an exponential clock with its own rate. The
minimum of N independent exponentials is itself exponential with the sum rate.
This is the "memoryless race" property of exponential distributions.

So if r1 ~ Uniform(0,1), then τ = -ln(r1) / a_0 samples from Exp(a_0) via
inverse CDF. And the probability that reaction j wins the race is a_j / a_0
(proportional to its firing rate). We implement this by picking a uniform point
r2 * a_0 in [0, a_0) and scanning the cumulative propensity sum.

The combined statement is: *given the current state x, draw (τ, j) jointly from
the distribution that says "the next event happens at time τ from now, and it's
reaction j."* Gillespie proved this gives the exact distribution of a continuous-
time Markov chain — no time-discretization error, no approximation.

### Numerical stability issues encountered

1. **ln(0) prevention**: `rand.Float64()` returns `[0, 1)`. If it returns exactly
   0, `ln(0) = -Inf` → infinite time step. We use `r1 = 1.0 - rng.Float64()` to
   map `[0, 1)` → `(0, 1]`, guaranteeing `r1 > 0`. This is standard SSA practice.

2. **Cumulative sum floating-point undershoot**: When selecting which reaction fires,
   rounding can cause the cumulative sum to never exceed `r2 * a0`. The fallback
   `fired = len(reactions) - 1` catches this. In practice with 2 channels it's nearly
   impossible, but it's essential for correctness with many channels in Phase 2.

3. **Propensity underflow / absorbing state**: If all propensities are zero, `a0 = 0`
   and no reaction can fire. We detect this and return `fired = -1, dt = +Inf` rather
   than dividing by zero. For the CpG toy system this can't happen (one channel is
   always active), but general systems can have absorbing states.

### Test results

7 tests pass, including 4 statistical validation cases against the analytical steady-
state probability P(methylated) = k_write / (k_write + k_erase):

| k_write | k_erase | Expected | Observed | Deviation |
|---------|---------|----------|----------|-----------|
| 0.3 | 0.7 | 0.3000 | 0.2998 | 0.0002 |
| 0.5 | 0.5 | 0.5000 | 0.4996 | 0.0004 |
| 0.9 | 0.1 | 0.9000 | 0.8992 | 0.0008 |
| 0.1 | 0.9 | 0.1000 | 0.1001 | 0.0001 |

All deviations are well within the ±0.02 tolerance (3000 trajectories, tMax=1000).

---

## Entry 3 — 2026-07-28T01:50+05:30 — Concurrent Ensemble Runner

### What was built

`gillespie/ensemble.go` — a worker-pool ensemble runner that launches N independent
SSA trajectories concurrently and aggregates results.

### Goroutine / worker-pool design

```
                 ┌─ worker 0 ─→ Run() → fraction → results ─┐
  jobs (indices) ├─ worker 1 ─→ Run() → fraction → results ─┤→ Welford aggregator
                 ├─ worker 2 ─→ Run() → fraction → results ─┤
                 └─ worker 3 ─→ Run() → fraction → results ─┘
```

**Why a fixed pool, not goroutine-per-trajectory?** If we launched 10,000 goroutines
(one per trajectory), the Go scheduler would handle it, but:
- Each goroutine has an 8KB initial stack that can grow. 10K goroutines = 80MB stacks
  before any work happens.
- All goroutines would race to allocate `TrajectoryRecord`s simultaneously, creating
  GC pressure spikes.
- On 4 cores, only 4 can run at once anyway — the rest just add scheduling overhead.

The worker pool keeps exactly `GOMAXPROCS` goroutines alive. Each worker:
1. Pulls a trajectory index from the buffered job channel.
2. Creates a deterministic `*rand.Rand` seeded from `baseSeed + index`.
3. Runs the full trajectory.
4. Computes the time-weighted fraction.
5. Sends just the scalar result (not the whole TrajectoryRecord) to the result channel.

**The key memory insight**: the `TrajectoryRecord` is discarded after computing the
fraction. So live memory is O(workers), not O(trajectories). With 4 workers and ~24KB
per trajectory, live footprint is ~96KB. Even the aggregated scalar results (8 bytes
per trajectory × 10K) is only 80KB. Total: under 200KB regardless of trajectory count.

### Aggregation: Welford's online algorithm

For computing mean and variance, we use Welford's online algorithm rather than the
naive "accumulate sum and sum-of-squares" approach. The naive method suffers from
catastrophic cancellation when the mean is large relative to the standard deviation —
`E[X²] - (E[X])²` subtracts two nearly-equal large numbers, losing precision. Welford's
method accumulates deviations from a running mean, which is numerically stable.

The 95% confidence interval uses `mean ± 1.96 * stddev / √N`, valid by CLT for large N.

### Deterministic seeding

Each trajectory gets `rand.NewPCG(baseSeed + trajectoryIndex, 0)`. This means:
- Same baseSeed → same results (reproducibility).
- No shared RNG → no mutex contention between workers.
- The second PCG parameter (stream) is 0 because each trajectory already has a unique
  seed; we don't need the stream dimension.

With Workers=1, the ordering is fully deterministic (FIFO through the job channel).
With Workers>1, the *set* of per-trajectory results is identical (each trajectory's
RNG is independent), but the *aggregation order* may vary — Welford's algorithm
produces the same result regardless of ordering, so the ensemble result is still
deterministic in exact arithmetic. In floating-point, there can be sub-ulp differences
due to addition ordering. The reproducibility test uses Workers=1 to avoid this.

---

## Entry 4 — 2026-07-28T01:55+05:30 — Phase 1 Complete: Final Report

### Summary

Phase 1 of the Riggs Engine is complete. The `gillespie` package implements the SSA
Direct Method as a general-purpose stochastic reaction simulator, validated against
the analytical two-state CTMC.

### What was built

| Component | File | Description |
|-----------|------|-------------|
| Core types | `gillespie/types.go` | Reaction, State, TrajectoryRecord, EnsembleResult |
| SSA engine | `gillespie/ssa.go` | Step(), Run(), TimeWeightedFraction() |
| Ensemble runner | `gillespie/ensemble.go` | Worker-pool concurrency, Welford aggregation |
| Tests | `gillespie/ssa_test.go` | 10 tests (7 unit + 3 ensemble), 4 statistical validations |
| Benchmarks | `gillespie/ssa_bench_test.go` | 3 benchmarks (Step, Run, Ensemble) |
| CLI | `cmd/riggs/main.go` | Toy CpG system runner |
| Docs | ARCHITECTURE.md, BIO_MAPPING.md, PROGRESS_LOG.md | Living documentation |

### Test coverage

10 tests, all passing:
- `TestStep_AbsorbingState` — deterministic edge case
- `TestStep_SingleActiveReaction` — deterministic correctness
- `TestStep_Deterministic` — seed reproducibility
- `TestRun_TwoState_SteadyState` — 4 statistical validations (k_w/k_e ratios)
- `TestRun_Reproducibility` — trajectory-level seed reproducibility
- `TestTimeWeightedFraction_EdgeCases` — 4 boundary conditions
- `TestClone` — deep copy independence
- `TestRunEnsemble_SteadyState` — concurrent statistical validation
- `TestRunEnsemble_Reproducibility` — ensemble seed reproducibility
- `TestRunEnsemble_SingleWorker` — degenerate pool case

### Benchmark results (Intel i3-1005G1 @ 1.20GHz, 4 cores)

| Benchmark | ns/op | B/op | allocs/op |
|-----------|------:|-----:|----------:|
| Step_TwoState | ~43 | 0 | **0** |
| Run_1000Events | ~49,000 | 68,921 | 425 |
| RunEnsemble_1000Traj | ~45M | 68.9M | 427K |

The critical result: **Step() achieves 0 allocs/op** with the pre-allocated scratch
buffer. The Run() allocations (425/trajectory) come from TrajectoryRecord growth and
per-event state snapshots — acceptable and expected for sparse event recording.

CLI throughput: **~16,000 trajectories/sec** with default settings (1000 trajectories,
tMax=1000, 4 workers).

### How Phase 2 plugs in (preview)

The Phase 1 core is deliberately biology-agnostic. Phase 2 will add biological
semantics as *reaction generators* — types that produce `[]Reaction` values for the
Phase 1 engine to simulate:

- **`Genome`**: a fixed DNA sequence with multiple CpG sites → the "disk/ROM". Each
  CpG site becomes a species in the `State.Counts` vector.
- **`dCas9` pointer**: navigates to a specific CpG site (index) without cutting →
  determines which reaction channels are active.
- **Methyltransferase / Demethylase**: each produces a `Reaction` for every targeted
  CpG site, with the appropriate Deltas and Propensity closures.
- **Environmental triggers**: events that modulate rate constants over time → will
  need an event-driven rate-update mechanism, or a time-dependent propensity function.

The key architectural decision Phase 2 needs to make: whether to keep the linear
propensity scan (O(R) per step, R = number of reaction channels) or switch to the
Gibson-Bruck Next Reaction Method (O(log R) per step via a priority queue of
next-firing-times). For 10-100 CpG sites (20-200 reaction channels), the linear
scan is fine. For 1000+ sites, the Next Reaction Method would help.

### Future optimizations (documented, not implemented)

- **Tau-leaping**: batch multiple reactions per step when individual event resolution
  isn't needed. Trades exactness for speed when propensities change slowly.
- **Gibson-Bruck Next Reaction Method**: maintain a priority queue of next-firing-times.
  Only recompute propensities for reactions whose reactants were affected by the last
  event. O(log R) per step instead of O(R).
- **Partial propensity recalculation**: instead of recomputing all propensities after
  each event, use the dependency graph to update only affected ones.

---

## Open Questions (for review before Phase 2)

- **Environmental trigger schema**: The WAL analogy implies recording *what* caused
  each methylation event (temperature spike? signal molecule?). What schema should
  triggers follow? A simple enum? A structured event with timestamp + type + magnitude?

- **Target reaction-channel count for Phase 2**: How many CpG sites should the
  initial genome model support? 10? 100? 1000? This affects whether we need sparse
  delta vectors and optimized propensity updates (Next Reaction Method) vs. the
  simple linear scan we use now.

- **State persistence**: Should trajectories be streamable to disk during simulation
  (for very long runs), or is in-memory accumulation sufficient for Phase 2?

- **Rate modulation**: Should environmental triggers change rate constants dynamically
  (time-inhomogeneous process), or should they be modeled as additional reaction
  channels that change state and indirectly affect propensities?
