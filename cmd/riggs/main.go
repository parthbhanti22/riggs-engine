// Command riggs runs stochastic simulations of CpG methylation systems using
// the Gillespie SSA engine.
//
// Phase 1 mode (toy): single CpG site, two-state CTMC.
// Phase 2 mode (genome): multi-site genome with targeted methylation via
// dCas9 complexes and environmental triggers.
// Phase 3 mode (tui): interactive Bubble Tea dashboard with live simulation.
//
// Usage:
//
//	riggs [-mode toy|genome|tui] [-trajectories N] [-tmax T] [-seed S] [-workers W]
//	      [-sites N] [-kw rate] [-ke rate]
package main

import (
	"flag"
	"fmt"
	"os"
	"runtime"
	"time"

	tea "github.com/charmbracelet/bubbletea"

	"github.com/pxrth9/riggs/bio"
	"github.com/pxrth9/riggs/gillespie"
	riggstui "github.com/pxrth9/riggs/tui"
)

func main() {
	// --- CLI flags ---
	mode := flag.String("mode", "genome", "simulation mode: 'toy' (Phase 1), 'genome' (Phase 2), or 'tui' (Phase 3 dashboard)")
	trajectories := flag.Int("trajectories", 1000, "number of independent trajectories")
	tMax := flag.Float64("tmax", 1000.0, "maximum simulation time per trajectory")
	kWrite := flag.Float64("kw", 0.3, "methyltransferase rate constant (toy mode only)")
	kErase := flag.Float64("ke", 0.7, "demethylase rate constant (toy mode only)")
	seed := flag.Uint64("seed", 42, "base seed for deterministic PRNG")
	workers := flag.Int("workers", runtime.GOMAXPROCS(0), "number of worker goroutines")
	nSites := flag.Int("sites", 50, "number of CpG sites (genome mode)")
	flag.Parse()

	switch *mode {
	case "toy":
		runToy(*trajectories, *tMax, *kWrite, *kErase, *seed, *workers)
	case "genome":
		runGenome(*trajectories, *tMax, *seed, *workers, *nSites)
	case "tui":
		runTUI(*tMax, *seed, *workers, *nSites)
	default:
		fmt.Fprintf(os.Stderr, "error: unknown mode %q (use 'toy', 'genome', or 'tui')\n", *mode)
		os.Exit(1)
	}
}

// runToy is the Phase 1 single-site simulation.
func runToy(trajectories int, tMax, kw, ke float64, seed uint64, workers int) {
	if kw <= 0 || ke <= 0 {
		fmt.Fprintf(os.Stderr, "error: rate constants must be positive (kw=%.4f, ke=%.4f)\n", kw, ke)
		os.Exit(1)
	}

	initial := gillespie.State{Counts: []int{0}, Time: 0}
	reactions := []gillespie.Reaction{
		{
			Deltas: []int{1}, RateConst: kw,
			Propensity: func(counts []int) float64 { return kw * float64(1-counts[0]) },
		},
		{
			Deltas: []int{-1}, RateConst: ke,
			Propensity: func(counts []int) float64 { return ke * float64(counts[0]) },
		},
	}

	fmt.Println("╔══════════════════════════════════════════════════════════╗")
	fmt.Println("║           RIGGS ENGINE — Phase 1: Gillespie SSA        ║")
	fmt.Println("║         Toy CpG Methylation System (2-state CTMC)      ║")
	fmt.Println("╚══════════════════════════════════════════════════════════╝")
	fmt.Println()
	fmt.Printf("  System:       1 CpG site, 2 reaction channels\n")
	fmt.Printf("  k_write:      %.4f (methyltransferase rate)\n", kw)
	fmt.Printf("  k_erase:      %.4f (demethylase rate)\n", ke)
	fmt.Printf("  Analytical:   P(methylated) = %.4f\n", kw/(kw+ke))
	fmt.Println()

	start := time.Now()
	result := gillespie.RunEnsemble(gillespie.EnsembleConfig{
		Initial: initial, Reactions: reactions, TMax: tMax,
		NumTrajectories: trajectories, Workers: workers,
		BaseSeed: seed, ObserveSpecies: 0,
	})
	elapsed := time.Since(start)

	expected := kw / (kw + ke)
	fmt.Println("  ── Results ──────────────────────────────────────────────")
	fmt.Printf("  Mean methylation fraction:  %.6f\n", result.MeanFraction)
	fmt.Printf("  Analytical expectation:     %.6f\n", expected)
	fmt.Printf("  Deviation:                  %.6f\n", result.MeanFraction-expected)
	fmt.Printf("  95%% CI:                     [%.6f, %.6f]\n", result.ConfLo, result.ConfHi)
	fmt.Printf("  Throughput:                 %.0f trajectories/sec\n", float64(result.N)/elapsed.Seconds())
	fmt.Println()

	if expected >= result.ConfLo && expected <= result.ConfHi {
		fmt.Println("  ✓ Analytical value falls within 95% CI — PASS")
	} else {
		fmt.Println("  ✗ Analytical value outside 95% CI — FAIL")
		os.Exit(1)
	}
}

// runGenome is the Phase 2 multi-site simulation with targeting and triggers.
func runGenome(trajectories int, tMax float64, seed uint64, workers, nSites int) {
	if nSites <= 0 {
		fmt.Fprintf(os.Stderr, "error: site count must be positive (%d)\n", nSites)
		os.Exit(1)
	}

	// --- Build genome system ---
	// CpG island: alternating promoter/gene-body context
	g := bio.NewGenome(nSites, 100, []bio.SiteContext{bio.ContextPromoter, bio.ContextGeneBody})

	// One targeting complex at site 0 (a methyltransferase-fused dCas9)
	complexes := []bio.TargetingComplex{
		{
			Index:      0,
			TargetSite: 0,
			KOff:       0.1,  // mean residence = 10 time units
			EnhWrite:   0.5,  // strong writer
			EnhErase:   0.05, // weak eraser (biased toward writing)
		},
	}

	// Periodic trigger every 50 time units
	triggerTimes := make([]float64, 0)
	for t := 0.0; t < tMax; t += 50.0 {
		triggerTimes = append(triggerTimes, t)
	}
	triggers := []bio.EnvironmentalTrigger{
		{
			Name:       "periodic_signal",
			ComplexIdx: 0,
			FireTimes:  triggerTimes,
		},
	}

	sys := &bio.System{
		Genome:    g,
		Complexes: complexes,
		Triggers:  triggers,
		KBgWrite:  0.001,
		KBgErase:  0.01,
	}
	buildResult := sys.Build()

	fmt.Println("╔══════════════════════════════════════════════════════════╗")
	fmt.Println("║           RIGGS ENGINE — Phase 2: Multi-Site Genome    ║")
	fmt.Println("║         CpG Island with Targeted Methylation           ║")
	fmt.Println("╚══════════════════════════════════════════════════════════╝")
	fmt.Println()
	fmt.Printf("  Genome:       %d CpG sites, %d reaction channels\n", buildResult.NumSites, len(buildResult.Reactions))
	fmt.Printf("  Species:      %d (%d sites + %d complexes)\n",
		len(buildResult.InitialState.Counts), buildResult.NumSites, buildResult.NumComplexes)
	fmt.Printf("  Complexes:    %d (targeting site 0)\n", buildResult.NumComplexes)
	fmt.Printf("  Triggers:     %d events (every 50 time units)\n", len(buildResult.Schedule))
	fmt.Printf("  k_bg_write:   %.4f\n", sys.KBgWrite)
	fmt.Printf("  k_bg_erase:   %.4f\n", sys.KBgErase)
	fmt.Printf("  k_enh_write:  %.4f (when complex bound)\n", complexes[0].EnhWrite)
	fmt.Printf("  k_enh_erase:  %.4f (when complex bound)\n", complexes[0].EnhErase)
	fmt.Printf("  k_off:        %.4f (complex unbinding)\n", complexes[0].KOff)
	fmt.Println()
	fmt.Printf("  Trajectories: %d\n", trajectories)
	fmt.Printf("  tMax:         %.1f\n", tMax)
	fmt.Printf("  Workers:      %d\n", workers)
	fmt.Printf("  Base seed:    %d\n", seed)
	fmt.Println()

	// --- Run for targeted site (0) ---
	start := time.Now()
	targetResult := gillespie.RunEnsemble(gillespie.EnsembleConfig{
		Initial:         buildResult.InitialState,
		Reactions:       buildResult.Reactions,
		TMax:            tMax,
		NumTrajectories: trajectories,
		Workers:         workers,
		BaseSeed:        seed,
		ObserveSpecies:  0,
		Schedule:        buildResult.Schedule,
	})

	// --- Run for a background site (1) ---
	bgResult := gillespie.RunEnsemble(gillespie.EnsembleConfig{
		Initial:         buildResult.InitialState,
		Reactions:       buildResult.Reactions,
		TMax:            tMax,
		NumTrajectories: trajectories,
		Workers:         workers,
		BaseSeed:        seed,
		ObserveSpecies:  1,
		Schedule:        buildResult.Schedule,
	})
	elapsed := time.Since(start)

	bgExpected := sys.KBgWrite / (sys.KBgWrite + sys.KBgErase)

	fmt.Println("  ── Results ──────────────────────────────────────────────")
	fmt.Println()
	fmt.Println("  Site 0 (targeted by dCas9-methyltransferase):")
	fmt.Printf("    Mean methylation:    %.6f\n", targetResult.MeanFraction)
	fmt.Printf("    Std dev:             %.6f\n", targetResult.StdDev)
	fmt.Printf("    95%% CI:              [%.6f, %.6f]\n", targetResult.ConfLo, targetResult.ConfHi)
	fmt.Println()
	fmt.Println("  Site 1 (background only, no targeting):")
	fmt.Printf("    Mean methylation:    %.6f\n", bgResult.MeanFraction)
	fmt.Printf("    Std dev:             %.6f\n", bgResult.StdDev)
	fmt.Printf("    95%% CI:              [%.6f, %.6f]\n", bgResult.ConfLo, bgResult.ConfHi)
	fmt.Printf("    Background expected: %.6f\n", bgExpected)
	fmt.Println()

	ratio := targetResult.MeanFraction / bgResult.MeanFraction
	fmt.Printf("  Targeted/background ratio: %.2fx\n", ratio)
	fmt.Printf("  Wall-clock time:           %s\n", elapsed.Round(time.Millisecond))
	fmt.Printf("  Throughput:                %.0f trajectories/sec (×2 sites)\n",
		float64(trajectories*2)/elapsed.Seconds())
	fmt.Println()

	// --- Validation ---
	if ratio > 2.0 {
		fmt.Printf("  ✓ Targeted site shows %.1fx enhancement — PASS\n", ratio)
	} else {
		fmt.Printf("  ✗ Enhancement ratio %.2f too low (expected >2x) — FAIL\n", ratio)
		os.Exit(1)
	}
}

// runTUI launches the Phase 3 interactive Bubble Tea dashboard.
func runTUI(tMax float64, seed uint64, workers, nSites int) {
	if nSites <= 0 {
		fmt.Fprintf(os.Stderr, "error: site count must be positive (%d)\n", nSites)
		os.Exit(1)
	}

	// --- Build genome system (same as genome mode) ---
	g := bio.NewGenome(nSites, 100, []bio.SiteContext{bio.ContextPromoter, bio.ContextGeneBody})

	complexes := []bio.TargetingComplex{
		{
			Index:      0,
			TargetSite: 0,
			KOff:       0.1,
			EnhWrite:   0.5,
			EnhErase:   0.05,
		},
	}

	triggerTimes := make([]float64, 0)
	for t := 0.0; t < tMax; t += 50.0 {
		triggerTimes = append(triggerTimes, t)
	}
	triggers := []bio.EnvironmentalTrigger{
		{
			Name:       "periodic_signal",
			ComplexIdx: 0,
			FireTimes:  triggerTimes,
		},
	}

	sys := &bio.System{
		Genome:    g,
		Complexes: complexes,
		Triggers:  triggers,
		KBgWrite:  0.001,
		KBgErase:  0.01,
	}
	buildResult := sys.Build()

	// --- Pre-compute ensemble reference stats ---
	// Run a quick ensemble to provide reference statistics in the dashboard.
	ensembleN := 200
	ensembleMean := make([]float64, nSites)
	ensembleStd := make([]float64, nSites)
	for site := 0; site < nSites && site < 2; site++ {
		result := gillespie.RunEnsemble(gillespie.EnsembleConfig{
			Initial:         buildResult.InitialState,
			Reactions:       buildResult.Reactions,
			TMax:            tMax,
			NumTrajectories: ensembleN,
			Workers:         workers,
			BaseSeed:        seed + 1000,
			ObserveSpecies:  site,
			Schedule:        buildResult.Schedule,
		})
		ensembleMean[site] = result.MeanFraction
		ensembleStd[site] = result.StdDev
	}

	// --- Launch TUI ---
	runner := riggstui.NewSimRunner(buildResult, g, complexes, tMax, seed)
	model := riggstui.NewModel(runner, g)
	model.SetEnsembleStats(ensembleMean, ensembleStd, ensembleN)

	p := tea.NewProgram(model, tea.WithAltScreen())
	if _, err := p.Run(); err != nil {
		fmt.Fprintf(os.Stderr, "error running TUI: %v\n", err)
		os.Exit(1)
	}
}
