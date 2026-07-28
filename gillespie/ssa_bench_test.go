package gillespie

import (
	"math/rand/v2"
	"testing"
)

// BenchmarkStep_TwoState measures the cost of a single SSA step for the
// two-state CpG system. This is the hot inner loop — the main target for
// allocation optimization.
//
// We pre-allocate the scratch buffer outside the loop to measure the
// steady-state cost (which is what the ensemble runner achieves).
func BenchmarkStep_TwoState(b *testing.B) {
	_, reactions := newCpGSystem(0.3, 0.7)
	state := State{Counts: []int{0}, Time: 0}
	rng := rand.New(rand.NewPCG(42, 0))
	scratch := make([]float64, len(reactions))

	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		Step(&state, reactions, rng, scratch)
	}
}

// BenchmarkRun_1000Events measures a full trajectory run of approximately
// 1000 events (tMax=1000 with total rate ~1.0).
func BenchmarkRun_1000Events(b *testing.B) {
	initial, reactions := newCpGSystem(0.3, 0.7)
	rng := rand.New(rand.NewPCG(42, 0))

	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		Run(initial, reactions, 1000.0, rng)
	}
}

// BenchmarkRunEnsemble_1000Traj measures the full ensemble runner with
// 1000 trajectories across 4 workers.
func BenchmarkRunEnsemble_1000Traj(b *testing.B) {
	initial, reactions := newCpGSystem(0.3, 0.7)

	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		RunEnsemble(EnsembleConfig{
			Initial:         initial,
			Reactions:       reactions,
			TMax:            1000.0,
			NumTrajectories: 1000,
			Workers:         4,
			BaseSeed:        42,
			ObserveSpecies:  0,
		})
	}
}
