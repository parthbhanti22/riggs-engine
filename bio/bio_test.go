package bio

import (
	"math"
	"math/rand/v2"
	"testing"

	"github.com/pxrth9/riggs/gillespie"
)

// --- Checkpoint 1: Genome/Site scaffolding + System.Build() basics ---

func TestNewGenome(t *testing.T) {
	t.Run("basic creation", func(t *testing.T) {
		g := NewGenome(10, 100, []SiteContext{ContextPromoter, ContextGeneBody})
		if g.NumSites() != 10 {
			t.Fatalf("expected 10 sites, got %d", g.NumSites())
		}
		for i, s := range g.Sites {
			if s.Index != i {
				t.Errorf("site %d: Index=%d, expected %d", i, s.Index, i)
			}
			if s.Coordinate != i*100 {
				t.Errorf("site %d: Coordinate=%d, expected %d", i, s.Coordinate, i*100)
			}
			expectedCtx := ContextPromoter
			if i%2 == 1 {
				expectedCtx = ContextGeneBody
			}
			if s.Context != expectedCtx {
				t.Errorf("site %d: Context=%v, expected %v", i, s.Context, expectedCtx)
			}
		}
	})

	t.Run("nil context pattern defaults to intergenic", func(t *testing.T) {
		g := NewGenome(3, 50, nil)
		for _, s := range g.Sites {
			if s.Context != ContextIntergenic {
				t.Errorf("site %d: expected ContextIntergenic, got %v", s.Index, s.Context)
			}
		}
	})
}

func TestSiteContextString(t *testing.T) {
	tests := []struct {
		ctx  SiteContext
		want string
	}{
		{ContextPromoter, "promoter"},
		{ContextGeneBody, "gene-body"},
		{ContextIntergenic, "intergenic"},
		{SiteContext(99), "unknown"},
	}
	for _, tt := range tests {
		if got := tt.ctx.String(); got != tt.want {
			t.Errorf("SiteContext(%d).String() = %q, want %q", tt.ctx, got, tt.want)
		}
	}
}

func TestBuild_ChannelCount(t *testing.T) {
	tests := []struct {
		name             string
		nSites           int
		nComplexes       int
		expectedChannels int
	}{
		{"10 sites, 0 complexes", 10, 0, 20},     // 2*10 + 0
		{"10 sites, 1 complex", 10, 1, 21},        // 2*10 + 1
		{"10 sites, 5 complexes", 10, 5, 25},      // 2*10 + 5
		{"100 sites, 0 complexes", 100, 0, 200},   // 2*100 + 0
		{"100 sites, 10 complexes", 100, 10, 210}, // 2*100 + 10
		{"1000 sites, 0 complexes", 1000, 0, 2000},
		{"1000 sites, 50 complexes", 1000, 50, 2050},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			g := NewGenome(tt.nSites, 100, nil)
			complexes := make([]TargetingComplex, tt.nComplexes)
			for i := range complexes {
				complexes[i] = TargetingComplex{
					Index:      i,
					TargetSite: i % tt.nSites,
					KOff:       1.0,
					EnhWrite:   0.1,
					EnhErase:   0.5,
				}
			}
			sys := &System{
				Genome:    g,
				Complexes: complexes,
				KBgWrite:  0.001,
				KBgErase:  0.01,
			}
			result := sys.Build()
			if len(result.Reactions) != tt.expectedChannels {
				t.Errorf("expected %d channels, got %d", tt.expectedChannels, len(result.Reactions))
			}
		})
	}
}

func TestBuild_SpeciesCount(t *testing.T) {
	tests := []struct {
		name    string
		nSites  int
		nComplx int
		want    int
	}{
		{"10 sites, 0 complexes", 10, 0, 10},
		{"10 sites, 3 complexes", 10, 3, 13},
		{"50 sites, 5 complexes", 50, 5, 55},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			g := NewGenome(tt.nSites, 100, nil)
			complexes := make([]TargetingComplex, tt.nComplx)
			for i := range complexes {
				complexes[i] = TargetingComplex{
					Index: i, TargetSite: 0, KOff: 1.0,
					EnhWrite: 0.1, EnhErase: 0.5,
				}
			}
			sys := &System{Genome: g, Complexes: complexes, KBgWrite: 0.001, KBgErase: 0.01}
			result := sys.Build()
			if len(result.InitialState.Counts) != tt.want {
				t.Errorf("expected %d species, got %d", tt.want, len(result.InitialState.Counts))
			}
		})
	}
}

func TestBuild_InitialStateAllZero(t *testing.T) {
	g := NewGenome(10, 100, nil)
	sys := &System{
		Genome:   g,
		KBgWrite: 0.001,
		KBgErase: 0.01,
	}
	result := sys.Build()
	for i, c := range result.InitialState.Counts {
		if c != 0 {
			t.Errorf("species %d: expected 0 (all unmethylated/unbound), got %d", i, c)
		}
	}
}

// TestBackgroundOnly_SteadyState validates that a genome with no targeting
// complexes reaches the expected background steady-state methylation at each site.
//
// For background-only sites, each site is an independent two-state CTMC:
//
//	P(methylated) = k_bg_write / (k_bg_write + k_bg_erase)
//
// With k_bg_write=0.01, k_bg_erase=0.02: expected = 0.01/0.03 ≈ 0.333
func TestBackgroundOnly_SteadyState(t *testing.T) {
	nSites := 10
	kBgW := 0.01
	kBgE := 0.02
	tMax := 5000.0
	nTraj := 2000
	expected := kBgW / (kBgW + kBgE) // 0.333...

	g := NewGenome(nSites, 100, nil)
	sys := &System{
		Genome:   g,
		KBgWrite: kBgW,
		KBgErase: kBgE,
	}
	result := sys.Build()

	// Run ensemble for each site and check convergence.
	// Since all sites are independent and identical (no complexes), we can
	// aggregate across both trajectories and sites for a stronger test.
	totalFrac := 0.0
	totalSamples := 0

	for traj := 0; traj < nTraj; traj++ {
		rng := rand.New(rand.NewPCG(uint64(traj), 0))
		rec := gillespie.Run(result.InitialState, result.Reactions, tMax, rng)

		for site := 0; site < nSites; site++ {
			frac := gillespie.TimeWeightedFraction(rec, site, 0, tMax)
			totalFrac += frac
			totalSamples++
		}
	}

	meanFrac := totalFrac / float64(totalSamples)
	deviation := math.Abs(meanFrac - expected)
	tolerance := 0.02

	if deviation > tolerance {
		t.Errorf("background steady-state: mean=%.4f, expected=%.4f, deviation=%.4f > %.4f",
			meanFrac, expected, deviation, tolerance)
	} else {
		t.Logf("PASS: mean=%.4f, expected=%.4f, deviation=%.4f (within ±%.4f, N=%d samples)",
			meanFrac, expected, deviation, tolerance, totalSamples)
	}
}

// --- Checkpoint 2: TargetingComplex binding + conditional propensity ---

// TestTargetedVsBackground is the key Phase 2 statistical validation test.
//
// Setup: 10 CpG sites. One targeting complex permanently bound at site 0
// (we set the bound species to 1 in the initial state, simulating a pre-bound
// complex — full trigger-driven binding is tested in Checkpoint 3).
//
// Expected behavior:
//   - Site 0 (targeted): effective write rate = k_bg_write + k_enh_write = 0.001 + 0.1 = 0.101
//     effective erase rate = k_bg_erase + k_enh_erase = 0.01 + 0.5 = 0.51
//     steady-state = 0.101 / (0.101 + 0.51) ≈ 0.1653
//   - Sites 1-9 (background only): write = 0.001, erase = 0.01
//     steady-state = 0.001 / (0.001 + 0.01) ≈ 0.0909
//
// The complex unbinding channel (KOff) would normally cause the complex to
// dissociate, but by starting with it bound and setting KOff = 0, we get a
// permanently bound complex for clean validation of the propensity logic.
func TestTargetedVsBackground(t *testing.T) {
	nSites := 10
	kBgW := 0.001
	kBgE := 0.01
	kEnhW := 0.1
	kEnhE := 0.5
	tMax := 5000.0
	nTraj := 2000

	g := NewGenome(nSites, 100, nil)
	sys := &System{
		Genome: g,
		Complexes: []TargetingComplex{
			{
				Index:      0,
				TargetSite: 0,
				KOff:       0, // permanently bound for this test
				EnhWrite:   kEnhW,
				EnhErase:   kEnhE,
			},
		},
		KBgWrite: kBgW,
		KBgErase: kBgE,
	}
	result := sys.Build()

	// Set complex 0 as initially bound
	initialWithBound := result.InitialState.Clone()
	initialWithBound.Counts[nSites] = 1 // species N+0 = bound state of complex 0

	expectedTargeted := (kBgW + kEnhW) / ((kBgW + kEnhW) + (kBgE + kEnhE))
	expectedBackground := kBgW / (kBgW + kBgE)

	var sumTargeted, sumBackground float64
	nBgSamples := 0

	for traj := 0; traj < nTraj; traj++ {
		rng := rand.New(rand.NewPCG(uint64(traj), 0))
		rec := gillespie.Run(initialWithBound, result.Reactions, tMax, rng)

		// Site 0: targeted
		sumTargeted += gillespie.TimeWeightedFraction(rec, 0, 0, tMax)

		// Sites 1-9: background only
		for site := 1; site < nSites; site++ {
			sumBackground += gillespie.TimeWeightedFraction(rec, site, 0, tMax)
			nBgSamples++
		}
	}

	meanTargeted := sumTargeted / float64(nTraj)
	meanBackground := sumBackground / float64(nBgSamples)

	t.Logf("Targeted site 0:   mean=%.4f, expected=%.4f", meanTargeted, expectedTargeted)
	t.Logf("Background sites:  mean=%.4f, expected=%.4f", meanBackground, expectedBackground)

	// Validate targeted site converges to enhanced rate steady-state
	if math.Abs(meanTargeted-expectedTargeted) > 0.02 {
		t.Errorf("targeted site: deviation %.4f > 0.02",
			math.Abs(meanTargeted-expectedTargeted))
	}

	// Validate background sites converge to background rate steady-state
	if math.Abs(meanBackground-expectedBackground) > 0.02 {
		t.Errorf("background sites: deviation %.4f > 0.02",
			math.Abs(meanBackground-expectedBackground))
	}

	// THE KEY ASSERTION: targeted site must be measurably different from background
	if meanTargeted <= meanBackground {
		t.Errorf("targeted site (%.4f) should have higher methylation than background (%.4f)",
			meanTargeted, meanBackground)
	}
	t.Logf("Targeted/background ratio: %.2fx", meanTargeted/meanBackground)
}

// TestComplexUnbinding validates that a bound complex unbinds stochastically
// with the expected mean residence time ≈ 1/KOff.
func TestComplexUnbinding(t *testing.T) {
	nSites := 1
	kOff := 0.5 // mean residence time = 2.0 time units

	g := NewGenome(nSites, 100, nil)
	sys := &System{
		Genome: g,
		Complexes: []TargetingComplex{
			{Index: 0, TargetSite: 0, KOff: kOff, EnhWrite: 0.1, EnhErase: 0.5},
		},
		KBgWrite: 0.001,
		KBgErase: 0.01,
	}
	result := sys.Build()

	// Start with complex bound
	initialBound := result.InitialState.Clone()
	initialBound.Counts[nSites] = 1 // complex 0 bound

	tMax := 100.0
	nTraj := 3000

	// Measure time-weighted fraction the complex is bound.
	// Expected: kOff unbinding competes with no binding (no trigger here),
	// so the complex should unbind once and never rebind.
	// Expected fraction bound = (1/kOff) / tMax for a single unbinding event
	// averaged over many trajectories... actually, it's just the exponential
	// wait. Mean time bound ≈ 1/kOff = 2.0, then unbound for the rest.
	// So fraction bound ≈ 2.0 / 100.0 = 0.02
	var sumBoundFrac float64
	for traj := 0; traj < nTraj; traj++ {
		rng := rand.New(rand.NewPCG(uint64(traj)+10000, 0))
		rec := gillespie.Run(initialBound, result.Reactions, tMax, rng)
		frac := gillespie.TimeWeightedFraction(rec, nSites, 1, tMax)
		sumBoundFrac += frac
	}

	meanBoundFrac := sumBoundFrac / float64(nTraj)
	expectedFrac := (1.0 / kOff) / tMax // ≈ 0.02

	t.Logf("Mean bound fraction: %.4f, expected ≈ %.4f", meanBoundFrac, expectedFrac)

	// Allow generous tolerance since the exponential distribution has high variance
	if math.Abs(meanBoundFrac-expectedFrac) > 0.01 {
		t.Errorf("bound fraction deviation %.4f > 0.01",
			math.Abs(meanBoundFrac-expectedFrac))
	}
}

// TestPropensityFastPath verifies that sites with no targeting complexes use
// the fast-path propensity (pure background rate).
func TestPropensityFastPath(t *testing.T) {
	g := NewGenome(5, 100, nil)
	sys := &System{
		Genome: g,
		Complexes: []TargetingComplex{
			{Index: 0, TargetSite: 2, KOff: 1.0, EnhWrite: 0.1, EnhErase: 0.5},
		},
		KBgWrite: 0.01,
		KBgErase: 0.02,
	}
	result := sys.Build()

	// Site 0 (untargeted): propensity should be pure background
	// When unmethylated (counts[0]=0): write propensity = 0.01
	counts := make([]int, len(result.InitialState.Counts))
	writeP := result.Reactions[0].Propensity(counts)   // site 0 write
	if math.Abs(writeP-0.01) > 1e-10 {
		t.Errorf("site 0 write propensity: got %.6f, expected 0.01", writeP)
	}

	// Site 2 (targeted, complex unbound): should still be background only
	writeP2 := result.Reactions[4].Propensity(counts) // site 2 write (reactions are: 0w,0e,1w,1e,2w,2e,...)
	if math.Abs(writeP2-0.01) > 1e-10 {
		t.Errorf("site 2 write propensity (complex unbound): got %.6f, expected 0.01", writeP2)
	}

	// Site 2 (targeted, complex bound): should be enhanced
	counts[5] = 1 // species 5 = complex 0 bound state (nSites=5, complex index=0)
	writeP2Enh := result.Reactions[4].Propensity(counts)
	expectedEnh := 0.01 + 0.1 // bg + enhanced
	if math.Abs(writeP2Enh-expectedEnh) > 1e-10 {
		t.Errorf("site 2 write propensity (complex bound): got %.6f, expected %.6f", writeP2Enh, expectedEnh)
	}
}

// --- Checkpoint 3: EnvironmentalTrigger + RunWithSchedule ---

// TestTriggerTiming verifies that a scheduled trigger fires at the correct
// time and binds the complex (not before, not after).
func TestTriggerTiming(t *testing.T) {
	nSites := 5
	triggerTime := 50.0

	g := NewGenome(nSites, 100, nil)
	sys := &System{
		Genome: g,
		Complexes: []TargetingComplex{
			{Index: 0, TargetSite: 0, KOff: 0.01, EnhWrite: 0.1, EnhErase: 0.5},
		},
		Triggers: []EnvironmentalTrigger{
			{Name: "test_trigger", ComplexIdx: 0, FireTimes: []float64{triggerTime}},
		},
		KBgWrite: 0.001,
		KBgErase: 0.01,
	}
	result := sys.Build()

	rng := rand.New(rand.NewPCG(42, 0))
	rec := gillespie.RunWithSchedule(result.InitialState, result.Reactions, result.Schedule, 200.0, rng)

	// Find the scheduled event in the trajectory record.
	// Scheduled events have FiredReaction = -(tag + 1) = -(0 + 1) = -1
	foundTrigger := false
	for i, fired := range rec.FiredReaction {
		if fired == -1 {
			foundTrigger = true
			if math.Abs(rec.Times[i]-triggerTime) > 1e-10 {
				t.Errorf("trigger fired at t=%.4f, expected t=%.4f", rec.Times[i], triggerTime)
			}
			// Verify the complex is now bound
			complexSpecies := nSites + 0 // species index for complex 0
			if rec.States[i][complexSpecies] != 1 {
				t.Errorf("complex not bound after trigger: counts[%d]=%d",
					complexSpecies, rec.States[i][complexSpecies])
			}
			break
		}
	}
	if !foundTrigger {
		t.Error("trigger event not found in trajectory record")
	}

	// Verify complex was NOT bound before the trigger time
	complexSpecies := nSites + 0
	for i, tTime := range rec.Times {
		if tTime < triggerTime && rec.States[i][complexSpecies] != 0 {
			t.Errorf("complex was bound at t=%.4f (before trigger at t=%.4f)",
				tTime, triggerTime)
			break
		}
	}
}

// TestTriggerIdempotency verifies that a trigger firing while the complex
// is already bound does not corrupt the state (clamping to [0,1]).
func TestTriggerIdempotency(t *testing.T) {
	nSites := 1

	g := NewGenome(nSites, 100, nil)
	sys := &System{
		Genome: g,
		Complexes: []TargetingComplex{
			{Index: 0, TargetSite: 0, KOff: 0.001, EnhWrite: 0.1, EnhErase: 0.5},
		},
		Triggers: []EnvironmentalTrigger{
			{
				Name:       "rapid_fire",
				ComplexIdx: 0,
				// Fire twice in quick succession — second should be a no-op
				FireTimes: []float64{10.0, 10.1},
			},
		},
		KBgWrite: 0.001,
		KBgErase: 0.01,
	}
	result := sys.Build()

	rng := rand.New(rand.NewPCG(42, 0))
	rec := gillespie.RunWithSchedule(result.InitialState, result.Reactions, result.Schedule, 100.0, rng)

	// After both triggers, bound state should still be exactly 1 (not 2)
	complexSpecies := nSites + 0
	for i, tTime := range rec.Times {
		if rec.States[i][complexSpecies] > 1 {
			t.Errorf("complex bound state=%d at t=%.4f — clamping failed",
				rec.States[i][complexSpecies], tTime)
		}
	}
}

// TestRunWithSchedule_Deterministic verifies that RunWithSchedule produces
// identical results with the same seed.
func TestRunWithSchedule_Deterministic(t *testing.T) {
	nSites := 5
	g := NewGenome(nSites, 100, nil)
	sys := &System{
		Genome: g,
		Complexes: []TargetingComplex{
			{Index: 0, TargetSite: 2, KOff: 0.5, EnhWrite: 0.1, EnhErase: 0.5},
		},
		Triggers: []EnvironmentalTrigger{
			{Name: "det_test", ComplexIdx: 0, FireTimes: []float64{50.0, 150.0}},
		},
		KBgWrite: 0.01,
		KBgErase: 0.02,
	}
	result := sys.Build()

	tMax := 200.0

	rng1 := rand.New(rand.NewPCG(123, 0))
	rec1 := gillespie.RunWithSchedule(result.InitialState, result.Reactions, result.Schedule, tMax, rng1)

	rng2 := rand.New(rand.NewPCG(123, 0))
	rec2 := gillespie.RunWithSchedule(result.InitialState, result.Reactions, result.Schedule, tMax, rng2)

	if len(rec1.Times) != len(rec2.Times) {
		t.Fatalf("different event counts: %d vs %d", len(rec1.Times), len(rec2.Times))
	}
	for i := range rec1.Times {
		if rec1.Times[i] != rec2.Times[i] {
			t.Errorf("event %d: time %.6f vs %.6f", i, rec1.Times[i], rec2.Times[i])
			break
		}
		if rec1.FiredReaction[i] != rec2.FiredReaction[i] {
			t.Errorf("event %d: fired %d vs %d", i, rec1.FiredReaction[i], rec2.FiredReaction[i])
			break
		}
	}
}

// TestFullTriggerDrivenMethylation is the end-to-end test for the full
// Phase 2 pipeline: genome → complexes → triggers → RunWithSchedule → ensemble.
//
// Setup: 10 sites. Complex targets site 0, trigger fires at t=0 (immediate).
// Complex has high enhanced write, low enhanced erase. KOff means it will
// unbind eventually, but the trigger keeps re-binding it.
//
// We compare site 0 methylation (targeted + repeatedly triggered) against
// untargeted sites (background only). Site 0 should be significantly higher.
func TestFullTriggerDrivenMethylation(t *testing.T) {
	nSites := 10
	kBgW := 0.001
	kBgE := 0.01
	tMax := 500.0
	nTraj := 1000

	g := NewGenome(nSites, 100, nil)
	sys := &System{
		Genome: g,
		Complexes: []TargetingComplex{
			{
				Index:      0,
				TargetSite: 0,
				KOff:       0.1, // mean residence = 10 time units
				EnhWrite:   0.5,
				EnhErase:   0.05,
			},
		},
		Triggers: []EnvironmentalTrigger{
			{
				Name:       "periodic_stimulus",
				ComplexIdx: 0,
				// Trigger fires every 20 time units — complex keeps getting re-bound
				FireTimes: func() []float64 {
					times := make([]float64, 0)
					for t := 0.0; t < tMax; t += 20.0 {
						times = append(times, t)
					}
					return times
				}(),
			},
		},
		KBgWrite: kBgW,
		KBgErase: kBgE,
	}
	result := sys.Build()

	// Use ensemble runner with schedule
	ensResult := gillespie.RunEnsemble(gillespie.EnsembleConfig{
		Initial:         result.InitialState,
		Reactions:       result.Reactions,
		TMax:            tMax,
		NumTrajectories: nTraj,
		Workers:         4,
		BaseSeed:        42,
		ObserveSpecies:  0, // observe site 0 (targeted)
		Schedule:        result.Schedule,
	})

	// Also compute background site methylation
	bgResult := gillespie.RunEnsemble(gillespie.EnsembleConfig{
		Initial:         result.InitialState,
		Reactions:       result.Reactions,
		TMax:            tMax,
		NumTrajectories: nTraj,
		Workers:         4,
		BaseSeed:        42,
		ObserveSpecies:  1, // observe site 1 (background only)
		Schedule:        result.Schedule,
	})

	t.Logf("Site 0 (targeted, triggered): mean=%.4f, stddev=%.4f", ensResult.MeanFraction, ensResult.StdDev)
	t.Logf("Site 1 (background):          mean=%.4f, stddev=%.4f", bgResult.MeanFraction, bgResult.StdDev)
	t.Logf("Targeted/background ratio:    %.2fx", ensResult.MeanFraction/bgResult.MeanFraction)

	// The targeted site should be measurably higher than background
	if ensResult.MeanFraction <= bgResult.MeanFraction {
		t.Errorf("targeted site (%.4f) should have higher methylation than background (%.4f)",
			ensResult.MeanFraction, bgResult.MeanFraction)
	}

	// The targeted site should be notably higher (not just marginally)
	ratio := ensResult.MeanFraction / bgResult.MeanFraction
	if ratio < 2.0 {
		t.Errorf("expected targeted/background ratio > 2.0, got %.2f", ratio)
	}
}
