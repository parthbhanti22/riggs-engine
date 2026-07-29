package bio

import (
	"sort"

	"github.com/pxrth9/riggs/gillespie"
)

// System assembles a Genome, TargetingComplexes, and EnvironmentalTriggers
// into a set of Gillespie reactions, an initial state, and a schedule of
// deterministic events.
//
// This is the core "reaction generator" pattern: the bio layer knows about
// biology; it produces gillespie.Reaction values that the math tier simulates.
// The System struct is the bridge between the two.
//
// Species layout in State.Counts:
//
//	[0, N)       methylation state of each CpG site (0=unmethylated, 1=methylated)
//	[N, N+K)     bound state of each targeting complex (0=unbound, 1=bound)
//
// Where N = genome.NumSites(), K = len(Complexes).
//
// Reaction channel layout:
//
//	[0, 2N)      per-site write + erase (background + conditional enhancement)
//	[2N, 2N+K)   per-complex unbinding
//	Total: 2N + K = O(N) channels (linear scaling, not combinatorial)
type System struct {
	// Genome is the immutable DNA with CpG sites — the "disk."
	Genome *Genome

	// Complexes is the set of dCas9-effector targeting complexes.
	// Each complex targets one site and contributes enhanced write/erase
	// rates when bound.
	Complexes []TargetingComplex

	// Triggers is the set of environmental triggers that bind complexes
	// at scheduled times.
	Triggers []EnvironmentalTrigger

	// KBgWrite is the background (spontaneous) methylation rate constant
	// per site. Represents passive/maintenance methyltransferase activity
	// that occurs independently of any targeting complex.
	KBgWrite float64

	// KBgErase is the background (spontaneous) demethylation rate constant
	// per site. Represents passive demethylation (e.g., replication-dependent
	// dilution of methylation marks).
	KBgErase float64
}

// BuildResult holds the output of System.Build().
type BuildResult struct {
	// Reactions is the set of Gillespie reaction channels.
	Reactions []gillespie.Reaction

	// InitialState is the starting state (all unmethylated, all unbound).
	InitialState gillespie.State

	// Schedule is the list of deterministic events (from triggers),
	// sorted by time. Nil if no triggers are configured.
	Schedule []gillespie.ScheduledEvent

	// NumSites is the number of CpG sites (for interpreting State.Counts).
	NumSites int

	// NumComplexes is the number of targeting complexes.
	NumComplexes int
}

// Build produces the Gillespie reactions, initial state, and scheduled events
// from the biological system configuration.
//
// The reaction generation follows this pattern:
//
//  1. For each CpG site: generate a write channel and an erase channel.
//     The propensity function checks both the site's methylation state AND
//     the bound state of any targeting complexes aimed at that site.
//     This is how enhanced rates work — not separate channels, but conditional
//     propensity within the same channel.
//
//  2. For each targeting complex: generate an unbinding channel.
//     Binding is handled by the deterministic trigger schedule, not a
//     stochastic channel.
//
//  3. For each trigger: generate ScheduledEvents at the specified fire times.
//
// Total channels: 2N + K where N = sites, K = complexes. This is O(N+K),
// satisfying the linear-scaling constraint.
func (sys *System) Build() BuildResult {
	N := sys.Genome.NumSites()
	K := len(sys.Complexes)
	numSpecies := N + K

	// --- Initial state: all unmethylated, all complexes unbound ---
	initialState := gillespie.State{
		Counts: make([]int, numSpecies),
		Time:   0,
	}

	// Pre-compute: for each site, which complexes target it?
	// This avoids repeated linear scans during propensity evaluation.
	siteComplexes := make([][]siteComplexInfo, N)
	for _, c := range sys.Complexes {
		if c.TargetSite >= 0 && c.TargetSite < N {
			siteComplexes[c.TargetSite] = append(siteComplexes[c.TargetSite], siteComplexInfo{
				speciesIdx: N + c.Index,
				enhWrite:   c.EnhWrite,
				enhErase:   c.EnhErase,
			})
		}
	}

	reactions := make([]gillespie.Reaction, 0, 2*N+K)

	// --- Per-site methylation channels (2N total) ---
	for _, site := range sys.Genome.Sites {
		siteIdx := site.Index
		kBgW := sys.KBgWrite
		kBgE := sys.KBgErase
		complexes := siteComplexes[siteIdx]

		// Write channel: propensity = (k_bg + Σ(k_enh_j * bound_j)) * (1 - meth_i)
		//
		// The closure captures siteIdx (which species to check), kBgW (background rate),
		// and the complexes slice (which bound-state species to check for enhancement).
		// All captures are immutable — safe for concurrent use.
		deltasW := make([]int, numSpecies)
		deltasW[siteIdx] = 1
		reactions = append(reactions, gillespie.Reaction{
			Deltas:    deltasW,
			RateConst: kBgW,
			Propensity: makeWritePropensity(siteIdx, kBgW, complexes),
		})

		// Erase channel: propensity = (k_bg + Σ(k_enh_j * bound_j)) * meth_i
		deltasE := make([]int, numSpecies)
		deltasE[siteIdx] = -1
		reactions = append(reactions, gillespie.Reaction{
			Deltas:    deltasE,
			RateConst: kBgE,
			Propensity: makeErasePropensity(siteIdx, kBgE, complexes),
		})
	}

	// --- Per-complex unbinding channels (K total) ---
	for _, c := range sys.Complexes {
		specIdx := N + c.Index
		kOff := c.KOff
		deltas := make([]int, numSpecies)
		deltas[specIdx] = -1
		reactions = append(reactions, gillespie.Reaction{
			Deltas:    deltas,
			RateConst: kOff,
			Propensity: func(counts []int) float64 {
				return kOff * float64(counts[specIdx])
			},
		})
	}

	// --- Scheduled events from triggers ---
	var schedule []gillespie.ScheduledEvent
	for _, trig := range sys.Triggers {
		specIdx := N + trig.ComplexIdx
		for _, t := range trig.FireTimes {
			deltas := make([]int, numSpecies)
			deltas[specIdx] = 1
			schedule = append(schedule, gillespie.ScheduledEvent{
				Time:   t,
				Deltas: deltas,
				Tag:    trig.ComplexIdx,
			})
		}
	}
	if len(schedule) > 0 {
		sort.Slice(schedule, func(i, j int) bool {
			return schedule[i].Time < schedule[j].Time
		})
	}

	return BuildResult{
		Reactions:    reactions,
		InitialState: initialState,
		Schedule:     schedule,
		NumSites:     N,
		NumComplexes: K,
	}
}

// siteComplexInfo holds pre-computed information about a targeting complex
// that targets a specific site. Used during reaction generation to build
// propensity closures without repeated lookups.
type siteComplexInfo struct {
	speciesIdx int     // index into State.Counts for this complex's bound state
	enhWrite   float64 // enhanced write rate when bound
	enhErase   float64 // enhanced erase rate when bound
}

// propensityTarget is a compact struct captured by propensity closures.
// Contains only the species index and rate relevant to one direction
// (write or erase), to minimize closure capture size.
type propensityTarget struct {
	specIdx int
	rate    float64
}

// makeWritePropensity creates the propensity closure for a site's write channel.
// Extracted as a named function to avoid closure-variable capture bugs in loops.
func makeWritePropensity(siteIdx int, kBg float64, complexes []siteComplexInfo) func([]int) float64 {
	if len(complexes) == 0 {
		// Fast path: no complexes target this site.
		// Pure background rate — no need to check bound states.
		return func(counts []int) float64 {
			return kBg * float64(1-counts[siteIdx])
		}
	}

	// Slow path: check bound states of targeting complexes.
	// Pre-copy the relevant data to avoid capturing the full complexes slice.
	infos := make([]propensityTarget, len(complexes))
	for i, c := range complexes {
		infos[i] = propensityTarget{specIdx: c.speciesIdx, rate: c.enhWrite}
	}

	return func(counts []int) float64 {
		if counts[siteIdx] == 1 {
			return 0 // already methylated
		}
		rate := kBg
		for _, info := range infos {
			rate += info.rate * float64(counts[info.specIdx])
		}
		return rate
	}
}

// makeErasePropensity creates the propensity closure for a site's erase channel.
func makeErasePropensity(siteIdx int, kBg float64, complexes []siteComplexInfo) func([]int) float64 {
	if len(complexes) == 0 {
		return func(counts []int) float64 {
			return kBg * float64(counts[siteIdx])
		}
	}

	infos := make([]propensityTarget, len(complexes))
	for i, c := range complexes {
		infos[i] = propensityTarget{specIdx: c.speciesIdx, rate: c.enhErase}
	}

	return func(counts []int) float64 {
		if counts[siteIdx] == 0 {
			return 0 // already unmethylated
		}
		rate := kBg
		for _, info := range infos {
			rate += info.rate * float64(counts[info.specIdx])
		}
		return rate
	}
}
