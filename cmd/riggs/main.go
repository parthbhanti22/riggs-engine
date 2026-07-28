// Command riggs runs a toy single-CpG-site stochastic simulation to demonstrate
// the Gillespie SSA engine.
//
// This is the Phase 1 entrypoint — it models one CpG dinucleotide that toggles
// between methylated and unmethylated states, driven by two reaction channels:
//   - Methylation (write): rate k_write, fires when unmethylated
//   - Demethylation (erase): rate k_erase, fires when methylated
//
// The ensemble average should converge to k_write / (k_write + k_erase).
//
// Usage:
//
//	riggs [-trajectories N] [-tmax T] [-kw rate] [-ke rate] [-seed S] [-workers W]
package main

import (
	"flag"
	"fmt"
	"os"
	"runtime"
	"time"

	"github.com/pxrth9/riggs/gillespie"
)

func main() {
	// --- CLI flags ---
	trajectories := flag.Int("trajectories", 1000, "number of independent trajectories")
	tMax := flag.Float64("tmax", 1000.0, "maximum simulation time per trajectory")
	kWrite := flag.Float64("kw", 0.3, "methyltransferase rate constant (write)")
	kErase := flag.Float64("ke", 0.7, "demethylase rate constant (erase)")
	seed := flag.Uint64("seed", 42, "base seed for deterministic PRNG")
	workers := flag.Int("workers", runtime.GOMAXPROCS(0), "number of worker goroutines")
	flag.Parse()

	// --- Validate inputs ---
	if *kWrite <= 0 || *kErase <= 0 {
		fmt.Fprintf(os.Stderr, "error: rate constants must be positive (kw=%.4f, ke=%.4f)\n", *kWrite, *kErase)
		os.Exit(1)
	}
	if *trajectories <= 0 {
		fmt.Fprintf(os.Stderr, "error: trajectory count must be positive (%d)\n", *trajectories)
		os.Exit(1)
	}
	if *tMax <= 0 {
		fmt.Fprintf(os.Stderr, "error: tmax must be positive (%.4f)\n", *tMax)
		os.Exit(1)
	}

	// --- Build the toy CpG system ---
	kw, ke := *kWrite, *kErase
	initial := gillespie.State{
		Counts: []int{0}, // start unmethylated
		Time:   0,
	}
	reactions := []gillespie.Reaction{
		{
			Deltas:    []int{1},
			RateConst: kw,
			Propensity: func(counts []int) float64 {
				return kw * float64(1-counts[0])
			},
		},
		{
			Deltas:    []int{-1},
			RateConst: ke,
			Propensity: func(counts []int) float64 {
				return ke * float64(counts[0])
			},
		},
	}

	// --- Run ensemble ---
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
	fmt.Printf("  Trajectories: %d\n", *trajectories)
	fmt.Printf("  tMax:         %.1f\n", *tMax)
	fmt.Printf("  Workers:      %d\n", *workers)
	fmt.Printf("  Base seed:    %d\n", *seed)
	fmt.Println()

	start := time.Now()
	result := gillespie.RunEnsemble(gillespie.EnsembleConfig{
		Initial:         initial,
		Reactions:       reactions,
		TMax:            *tMax,
		NumTrajectories: *trajectories,
		Workers:         *workers,
		BaseSeed:        *seed,
		ObserveSpecies:  0,
	})
	elapsed := time.Since(start)

	// --- Report ---
	expected := kw / (kw + ke)
	fmt.Println("  ── Results ──────────────────────────────────────────────")
	fmt.Printf("  Mean methylation fraction:  %.6f\n", result.MeanFraction)
	fmt.Printf("  Analytical expectation:     %.6f\n", expected)
	fmt.Printf("  Deviation:                  %.6f\n", result.MeanFraction-expected)
	fmt.Printf("  Std dev (per-trajectory):   %.6f\n", result.StdDev)
	fmt.Printf("  95%% CI:                     [%.6f, %.6f]\n", result.ConfLo, result.ConfHi)
	fmt.Printf("  N trajectories:             %d\n", result.N)
	fmt.Printf("  Wall-clock time:            %s\n", elapsed.Round(time.Millisecond))
	fmt.Printf("  Throughput:                 %.0f trajectories/sec\n", float64(result.N)/elapsed.Seconds())
	fmt.Println()

	// --- Validation ---
	if expected >= result.ConfLo && expected <= result.ConfHi {
		fmt.Println("  ✓ Analytical value falls within 95% CI — PASS")
	} else {
		fmt.Println("  ✗ Analytical value outside 95% CI — FAIL")
		os.Exit(1)
	}
}
