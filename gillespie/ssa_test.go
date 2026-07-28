package gillespie

import (
	"math"
	"math/rand/v2"
	"testing"
)

// newCpGSystem creates the toy two-state CpG methylation system.
// Species 0: methylation state (0 = unmethylated, 1 = methylated).
// Reaction 0: methylation (write), rate kWrite. Fires only when unmethylated.
// Reaction 1: demethylation (erase), rate kErase. Fires only when methylated.
func newCpGSystem(kWrite, kErase float64) (State, []Reaction) {
	initial := State{
		Counts: []int{0}, // start unmethylated
		Time:   0,
	}
	reactions := []Reaction{
		{
			Deltas:    []int{1},
			RateConst: kWrite,
			Propensity: func(counts []int) float64 {
				// a_write = k_write * (1 - counts[0])
				// Fires only when unmethylated (counts[0] == 0).
				return kWrite * float64(1-counts[0])
			},
		},
		{
			Deltas:    []int{-1},
			RateConst: kErase,
			Propensity: func(counts []int) float64 {
				// a_erase = k_erase * counts[0]
				// Fires only when methylated (counts[0] == 1).
				return kErase * float64(counts[0])
			},
		},
	}
	return initial, reactions
}

// TestStep_AbsorbingState verifies that Step returns fired=-1 and dt=+Inf
// when all propensities are zero (no reaction can fire).
func TestStep_AbsorbingState(t *testing.T) {
	state := State{Counts: []int{5}, Time: 0}
	// A reaction that requires counts[0] == 0, but counts[0] is 5.
	reactions := []Reaction{
		{
			Deltas:    []int{1},
			RateConst: 1.0,
			Propensity: func(counts []int) float64 {
				if counts[0] == 0 {
					return 1.0
				}
				return 0.0
			},
		},
	}
	rng := rand.New(rand.NewPCG(42, 0))
	scratch := make([]float64, len(reactions))

	fired, dt := Step(&state, reactions, rng, scratch)
	if fired != -1 {
		t.Errorf("expected fired=-1 for absorbing state, got %d", fired)
	}
	if !math.IsInf(dt, 1) {
		t.Errorf("expected dt=+Inf for absorbing state, got %f", dt)
	}
	// State should be unchanged
	if state.Counts[0] != 5 {
		t.Errorf("state was modified in absorbing state: counts[0]=%d", state.Counts[0])
	}
	if state.Time != 0 {
		t.Errorf("time was modified in absorbing state: time=%f", state.Time)
	}
}

// TestStep_SingleActiveReaction verifies that when only one reaction channel
// has nonzero propensity, it always fires.
func TestStep_SingleActiveReaction(t *testing.T) {
	// CpG system starting unmethylated: only the methylation reaction can fire.
	state := State{Counts: []int{0}, Time: 0}
	_, reactions := newCpGSystem(0.5, 0.5)

	rng := rand.New(rand.NewPCG(123, 0))
	scratch := make([]float64, len(reactions))

	// Run 100 steps. The first step must always fire reaction 0 (methylation),
	// since counts[0]=0 means only the write reaction has nonzero propensity.
	fired, dt := Step(&state, reactions, rng, scratch)
	if fired != 0 {
		t.Errorf("expected reaction 0 (methylation) to fire, got %d", fired)
	}
	if dt <= 0 {
		t.Errorf("expected positive dt, got %f", dt)
	}
	if state.Counts[0] != 1 {
		t.Errorf("expected counts[0]=1 after methylation, got %d", state.Counts[0])
	}
}

// TestStep_Deterministic verifies that the same seed produces identical results.
func TestStep_Deterministic(t *testing.T) {
	state1 := State{Counts: []int{0}, Time: 0}
	state2 := State{Counts: []int{0}, Time: 0}
	_, reactions := newCpGSystem(0.3, 0.7)

	rng1 := rand.New(rand.NewPCG(999, 0))
	rng2 := rand.New(rand.NewPCG(999, 0))
	scratch1 := make([]float64, len(reactions))
	scratch2 := make([]float64, len(reactions))

	for i := 0; i < 50; i++ {
		f1, dt1 := Step(&state1, reactions, rng1, scratch1)
		f2, dt2 := Step(&state2, reactions, rng2, scratch2)
		if f1 != f2 || dt1 != dt2 {
			t.Fatalf("step %d: seed reproducibility failed: (%d,%f) vs (%d,%f)", i, f1, dt1, f2, dt2)
		}
	}
}

// TestRun_TwoState_SteadyState is the core statistical validation test.
//
// For the two-state CTMC:
//
//	unmethylated --[k_write]--> methylated
//	methylated   --[k_erase]--> unmethylated
//
// The steady-state probability of being methylated is:
//
//	P(methylated) = k_write / (k_write + k_erase)
//
// We run many independent trajectories and check that the ensemble average
// of the time-weighted methylation fraction converges to this value.
func TestRun_TwoState_SteadyState(t *testing.T) {
	tests := []struct {
		name      string
		kWrite    float64
		kErase    float64
		tMax      float64
		nTraj     int
		tolerance float64 // allowed deviation from analytical result
	}{
		{
			name:      "k_w=0.3, k_e=0.7, expect P=0.3",
			kWrite:    0.3,
			kErase:    0.7,
			tMax:      1000.0,
			nTraj:     3000,
			tolerance: 0.02,
		},
		{
			name:      "k_w=0.5, k_e=0.5, expect P=0.5",
			kWrite:    0.5,
			kErase:    0.5,
			tMax:      1000.0,
			nTraj:     3000,
			tolerance: 0.02,
		},
		{
			name:      "k_w=0.9, k_e=0.1, expect P=0.9",
			kWrite:    0.9,
			kErase:    0.1,
			tMax:      1000.0,
			nTraj:     3000,
			tolerance: 0.02,
		},
		{
			name:      "asymmetric: k_w=0.1, k_e=0.9, expect P=0.1",
			kWrite:    0.1,
			kErase:    0.9,
			tMax:      1000.0,
			nTraj:     3000,
			tolerance: 0.02,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			initial, reactions := newCpGSystem(tt.kWrite, tt.kErase)
			expected := tt.kWrite / (tt.kWrite + tt.kErase)

			sumFrac := 0.0
			for i := 0; i < tt.nTraj; i++ {
				rng := rand.New(rand.NewPCG(uint64(i), 0))
				rec := Run(initial, reactions, tt.tMax, rng)
				frac := TimeWeightedFraction(rec, 0, initial.Counts[0], tt.tMax)
				sumFrac += frac
			}
			meanFrac := sumFrac / float64(tt.nTraj)

			deviation := math.Abs(meanFrac - expected)
			if deviation > tt.tolerance {
				t.Errorf("mean methylation fraction = %.4f, expected %.4f (deviation %.4f > tolerance %.4f)",
					meanFrac, expected, deviation, tt.tolerance)
			} else {
				t.Logf("PASS: mean=%.4f, expected=%.4f, deviation=%.4f (within ±%.4f)",
					meanFrac, expected, deviation, tt.tolerance)
			}
		})
	}
}

// TestRun_Reproducibility verifies that Run produces identical trajectories
// with the same seed.
func TestRun_Reproducibility(t *testing.T) {
	initial, reactions := newCpGSystem(0.3, 0.7)
	tMax := 100.0

	rng1 := rand.New(rand.NewPCG(42, 0))
	rng2 := rand.New(rand.NewPCG(42, 0))

	rec1 := Run(initial, reactions, tMax, rng1)
	rec2 := Run(initial, reactions, tMax, rng2)

	if len(rec1.Times) != len(rec2.Times) {
		t.Fatalf("trajectory lengths differ: %d vs %d", len(rec1.Times), len(rec2.Times))
	}
	for i := range rec1.Times {
		if rec1.Times[i] != rec2.Times[i] {
			t.Fatalf("event %d: times differ: %f vs %f", i, rec1.Times[i], rec2.Times[i])
		}
		if rec1.FiredReaction[i] != rec2.FiredReaction[i] {
			t.Fatalf("event %d: reactions differ: %d vs %d", i, rec1.FiredReaction[i], rec2.FiredReaction[i])
		}
	}
}

// TestTimeWeightedFraction_EdgeCases tests the observable computation for
// boundary conditions.
func TestTimeWeightedFraction_EdgeCases(t *testing.T) {
	tests := []struct {
		name         string
		rec          TrajectoryRecord
		initialCount int
		tMax         float64
		expected     float64
	}{
		{
			name:         "empty trajectory, start unmethylated",
			rec:          TrajectoryRecord{},
			initialCount: 0,
			tMax:         100.0,
			expected:     0.0, // never methylated
		},
		{
			name:         "empty trajectory, start methylated",
			rec:          TrajectoryRecord{},
			initialCount: 1,
			tMax:         100.0,
			expected:     1.0, // always methylated
		},
		{
			name: "single event at midpoint",
			rec: TrajectoryRecord{
				Times:         []float64{50.0},
				States:        [][]int{{1}},
				FiredReaction: []int{0},
			},
			initialCount: 0,
			tMax:         100.0,
			expected:     0.5, // unmethylated for 50, methylated for 50
		},
		{
			name:         "tMax = 0",
			rec:          TrajectoryRecord{},
			initialCount: 1,
			tMax:         0,
			expected:     0.0, // no time window
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := TimeWeightedFraction(tt.rec, 0, tt.initialCount, tt.tMax)
			if math.Abs(got-tt.expected) > 1e-10 {
				t.Errorf("got %.6f, expected %.6f", got, tt.expected)
			}
		})
	}
}

// TestClone verifies that State.Clone produces an independent deep copy.
func TestClone(t *testing.T) {
	original := State{Counts: []int{1, 2, 3}, Time: 5.5}
	clone := original.Clone()

	// Modify clone
	clone.Counts[0] = 99
	clone.Time = 99.9

	// Original should be unchanged
	if original.Counts[0] != 1 {
		t.Errorf("Clone aliased Counts: original.Counts[0]=%d", original.Counts[0])
	}
	if original.Time != 5.5 {
		t.Errorf("Clone aliased Time: original.Time=%f", original.Time)
	}
}

// --- Ensemble runner tests ---

// TestRunEnsemble_SteadyState validates the concurrent ensemble runner against
// the same analytical result as the sequential test.
func TestRunEnsemble_SteadyState(t *testing.T) {
	initial, reactions := newCpGSystem(0.3, 0.7)
	expected := 0.3

	result := RunEnsemble(EnsembleConfig{
		Initial:         initial,
		Reactions:       reactions,
		TMax:            1000.0,
		NumTrajectories: 5000,
		Workers:         4,
		BaseSeed:        42,
		ObserveSpecies:  0,
	})

	deviation := math.Abs(result.MeanFraction - expected)
	if deviation > 0.02 {
		t.Errorf("ensemble mean=%.4f, expected=%.4f, deviation=%.4f > 0.02",
			result.MeanFraction, expected, deviation)
	}

	// The analytical result should fall within the 95% CI
	if expected < result.ConfLo || expected > result.ConfHi {
		t.Errorf("analytical value %.4f outside 95%% CI [%.4f, %.4f]",
			expected, result.ConfLo, result.ConfHi)
	}

	t.Logf("Ensemble: mean=%.4f, stddev=%.4f, CI=[%.4f, %.4f], N=%d",
		result.MeanFraction, result.StdDev, result.ConfLo, result.ConfHi, result.N)
}

// TestRunEnsemble_Reproducibility verifies that the same base seed produces
// identical ensemble results.
func TestRunEnsemble_Reproducibility(t *testing.T) {
	initial, reactions := newCpGSystem(0.3, 0.7)

	cfg := EnsembleConfig{
		Initial:         initial,
		Reactions:       reactions,
		TMax:            100.0,
		NumTrajectories: 500,
		Workers:         1, // single worker for deterministic ordering
		BaseSeed:        99,
		ObserveSpecies:  0,
	}

	r1 := RunEnsemble(cfg)
	r2 := RunEnsemble(cfg)

	if r1.MeanFraction != r2.MeanFraction {
		t.Errorf("non-reproducible: mean1=%.6f, mean2=%.6f", r1.MeanFraction, r2.MeanFraction)
	}
	if r1.StdDev != r2.StdDev {
		t.Errorf("non-reproducible: stddev1=%.6f, stddev2=%.6f", r1.StdDev, r2.StdDev)
	}
}

// TestRunEnsemble_SingleWorker verifies the ensemble runner works with a
// single worker (degenerate case, no concurrency).
func TestRunEnsemble_SingleWorker(t *testing.T) {
	initial, reactions := newCpGSystem(0.5, 0.5)

	result := RunEnsemble(EnsembleConfig{
		Initial:         initial,
		Reactions:       reactions,
		TMax:            100.0,
		NumTrajectories: 100,
		Workers:         1,
		BaseSeed:        7,
		ObserveSpecies:  0,
	})

	if result.N != 100 {
		t.Errorf("expected N=100, got %d", result.N)
	}
	// Just a sanity check that the result is in a reasonable range.
	if result.MeanFraction < 0.3 || result.MeanFraction > 0.7 {
		t.Errorf("mean=%.4f outside reasonable range for k_w=k_e=0.5", result.MeanFraction)
	}
}
