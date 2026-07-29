package gillespie

import (
	"math"
	"math/rand/v2"
	"runtime"
	"sync"
)

// EnsembleConfig configures the concurrent ensemble runner.
type EnsembleConfig struct {
	// Initial is the starting state for every trajectory (cloned per trajectory).
	Initial State

	// Reactions is the set of reaction channels (shared read-only across workers).
	Reactions []Reaction

	// TMax is the maximum simulation time for each trajectory.
	TMax float64

	// NumTrajectories is the total number of independent trajectories to simulate.
	NumTrajectories int

	// Workers is the number of concurrent worker goroutines.
	// Defaults to runtime.GOMAXPROCS(0) if zero.
	Workers int

	// BaseSeed is the base seed for deterministic PRNG seeding.
	// Each trajectory i gets its own RNG seeded from (BaseSeed + i).
	// This ensures full reproducibility and zero contention (no shared RNG).
	BaseSeed uint64

	// ObserveSpecies is the index of the species to compute time-weighted
	// fraction for in the ensemble result.
	ObserveSpecies int

	// Schedule is an optional list of deterministic events (sorted by time)
	// to interleave with stochastic SSA steps. If non-nil and non-empty,
	// workers call RunWithSchedule() instead of Run().
	// Phase 1 callers leave this nil for the original behavior.
	Schedule []ScheduledEvent
}

// trajectoryResult holds the per-trajectory output sent back on the result channel.
type trajectoryResult struct {
	fraction float64 // time-weighted fraction of the observed species
}

// RunEnsemble executes NumTrajectories independent SSA trajectories using a
// fixed-size worker pool and aggregates the results.
//
// Architecture:
//
//	                ┌─ worker 0 ─→ Run() → fraction → results ─┐
//	 jobs (indices) ├─ worker 1 ─→ Run() → fraction → results ─┤→ aggregator
//	                ├─ worker 2 ─→ Run() → fraction → results ─┤
//	                └─ worker 3 ─→ Run() → fraction → results ─┘
//
// The job channel carries trajectory indices. Each worker:
//  1. Receives a trajectory index from the job channel.
//  2. Creates a deterministic RNG seeded from (BaseSeed + index).
//  3. Runs a full trajectory with Run().
//  4. Computes the time-weighted fraction and sends it on the result channel.
//
// Memory implications on 8GB (3.7GB WSL):
//   - Each worker holds one live State + one TrajectoryRecord at a time.
//   - TrajectoryRecord is discarded after computing the fraction — only the
//     scalar result survives. So live memory is O(workers), not O(trajectories).
//   - With 4 workers and ~1000 events/trajectory at ~24 bytes/event, live
//     footprint is roughly 4 × 24KB ≈ 96KB. Trivial.
func RunEnsemble(cfg EnsembleConfig) EnsembleResult {
	workers := cfg.Workers
	if workers <= 0 {
		workers = runtime.GOMAXPROCS(0)
	}
	// Don't spawn more workers than trajectories.
	if workers > cfg.NumTrajectories {
		workers = cfg.NumTrajectories
	}

	jobs := make(chan int, cfg.NumTrajectories)
	results := make(chan trajectoryResult, cfg.NumTrajectories)

	// Launch worker pool.
	var wg sync.WaitGroup
	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for trajIdx := range jobs {
				// Deterministic per-trajectory seed: baseSeed + trajectory index.
				// Using PCG with (seed, 0) — the second parameter is the stream,
				// which we leave at 0 since each trajectory already has a unique seed.
				rng := rand.New(rand.NewPCG(cfg.BaseSeed+uint64(trajIdx), 0))
				var rec TrajectoryRecord
				if len(cfg.Schedule) > 0 {
					rec = RunWithSchedule(cfg.Initial, cfg.Reactions, cfg.Schedule, cfg.TMax, rng)
				} else {
					rec = Run(cfg.Initial, cfg.Reactions, cfg.TMax, rng)
				}
				frac := TimeWeightedFraction(rec, cfg.ObserveSpecies, cfg.Initial.Counts[cfg.ObserveSpecies], cfg.TMax)
				results <- trajectoryResult{fraction: frac}
			}
		}()
	}

	// Feed jobs.
	for i := 0; i < cfg.NumTrajectories; i++ {
		jobs <- i
	}
	close(jobs)

	// Wait for all workers to finish, then close results.
	go func() {
		wg.Wait()
		close(results)
	}()

	// Aggregate results: compute mean, stddev, 95% CI.
	// We use Welford's online algorithm for numerical stability —
	// it avoids catastrophic cancellation that can occur with the
	// naive (sum-of-squares - square-of-sum) formula when the mean
	// is large relative to the variance.
	n := 0
	mean := 0.0
	m2 := 0.0 // sum of squared deviations from the running mean

	for r := range results {
		n++
		delta := r.fraction - mean
		mean += delta / float64(n)
		delta2 := r.fraction - mean
		m2 += delta * delta2
	}

	var stddev float64
	if n > 1 {
		stddev = math.Sqrt(m2 / float64(n-1)) // sample standard deviation
	}

	// 95% confidence interval: mean ± 1.96 * stddev / sqrt(N)
	// This assumes the per-trajectory fractions are approximately normally
	// distributed (justified by CLT for large N).
	se := stddev / math.Sqrt(float64(n)) // standard error of the mean
	return EnsembleResult{
		MeanFraction: mean,
		StdDev:       stddev,
		ConfLo:       mean - 1.96*se,
		ConfHi:       mean + 1.96*se,
		N:            n,
	}
}
