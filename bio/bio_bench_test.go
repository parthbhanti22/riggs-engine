package bio

import (
	"fmt"
	"math/rand/v2"
	"testing"

	"github.com/pxrth9/riggs/gillespie"
)

// --- Phase 2 Scaling Benchmarks ---
//
// These benchmarks verify that reaction-channel count and simulation time
// scale linearly with site count N, confirming the O(N+K) channel design.
//
// Benchmark matrix:
//   N=10    → 21 channels  (2*10 + 1 complex)
//   N=100   → 201 channels (2*100 + 1 complex)
//   N=1000  → 2001 channels (2*1000 + 1 complex)

func makeBenchSystem(nSites int) (*System, BuildResult) {
	g := NewGenome(nSites, 100, nil)
	sys := &System{
		Genome: g,
		Complexes: []TargetingComplex{
			{Index: 0, TargetSite: 0, KOff: 0.1, EnhWrite: 0.1, EnhErase: 0.5},
		},
		Triggers: []EnvironmentalTrigger{
			{Name: "bench_trigger", ComplexIdx: 0, FireTimes: []float64{10.0}},
		},
		KBgWrite: 0.001,
		KBgErase: 0.01,
	}
	result := sys.Build()
	return sys, result
}

// BenchmarkBuild measures the time to construct reactions from the bio system.
func BenchmarkBuild(b *testing.B) {
	for _, nSites := range []int{10, 100, 1000} {
		g := NewGenome(nSites, 100, nil)
		sys := &System{
			Genome: g,
			Complexes: []TargetingComplex{
				{Index: 0, TargetSite: 0, KOff: 0.1, EnhWrite: 0.1, EnhErase: 0.5},
			},
			KBgWrite: 0.001,
			KBgErase: 0.01,
		}
		b.Run(fmt.Sprintf("N=%d", nSites), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				result := sys.Build()
				_ = result.Reactions // prevent optimization
			}
		})
	}
}

// BenchmarkRunTrajectory measures per-trajectory time with RunWithSchedule.
// Uses a short tMax to keep event counts manageable at large N.
func BenchmarkRunTrajectory(b *testing.B) {
	for _, nSites := range []int{10, 100, 1000} {
		_, result := makeBenchSystem(nSites)

		// Scale tMax so each benchmark takes a comparable amount of time.
		// With k_bg ≈ 0.011 per site, total rate ≈ 0.022*N.
		// Want ~500 events per trajectory: tMax ≈ 500 / (0.022*N)
		tMax := 500.0 / (0.022 * float64(nSites))
		if tMax < 10 {
			tMax = 10
		}

		b.Run(fmt.Sprintf("N=%d", nSites), func(b *testing.B) {
			b.ReportMetric(float64(len(result.Reactions)), "channels")
			b.ReportMetric(float64(len(result.InitialState.Counts)), "species")
			for i := 0; i < b.N; i++ {
				rng := rand.New(rand.NewPCG(uint64(i), 0))
				rec := gillespie.RunWithSchedule(
					result.InitialState, result.Reactions,
					result.Schedule, tMax, rng,
				)
				_ = rec.Times
			}
		})
	}
}

// BenchmarkEnsemble measures ensemble throughput at different scales.
func BenchmarkEnsemble(b *testing.B) {
	for _, nSites := range []int{10, 100, 1000} {
		_, result := makeBenchSystem(nSites)
		tMax := 500.0 / (0.022 * float64(nSites))
		if tMax < 10 {
			tMax = 10
		}

		b.Run(fmt.Sprintf("N=%d", nSites), func(b *testing.B) {
			nTraj := 100
			b.ReportMetric(float64(len(result.Reactions)), "channels")
			for i := 0; i < b.N; i++ {
				ens := gillespie.RunEnsemble(gillespie.EnsembleConfig{
					Initial:         result.InitialState,
					Reactions:       result.Reactions,
					TMax:            tMax,
					NumTrajectories: nTraj,
					Workers:         4,
					BaseSeed:        uint64(i * nTraj),
					ObserveSpecies:  0,
					Schedule:        result.Schedule,
				})
				_ = ens.MeanFraction
			}
		})
	}
}
