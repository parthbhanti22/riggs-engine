package bio

// TargetingComplex represents a catalytically-dead dCas9 protein fused to a
// guide RNA, targeting one specific CpG site in the genome.
//
// In the biological system:
//   - dCas9 (catalytically dead Cas9) navigates to a specific genomic address
//     without cutting the DNA — it's the "pointer" in the ROM/RAM analogy.
//   - A guide RNA determines which site the dCas9 binds to (the "address").
//   - An effector domain (methyltransferase or demethylase) is fused to the
//     dCas9, providing enhanced write/erase capability at the target site.
//
// In the Gillespie state space:
//   - The complex's bound state is a species: State.Counts[genome.NumSites + Index].
//   - Bound = 1 (complex is attached to its target site), Unbound = 0.
//   - Binding is driven by EnvironmentalTrigger (deterministic schedule).
//   - Unbinding is a stochastic reaction channel with rate KOff.
//   - While bound, the complex enhances the write/erase propensity at its target
//     site — the propensity function checks the bound-state species, not a
//     separate reaction type.
type TargetingComplex struct {
	// Index is this complex's identifier.
	// Its bound-state species index = genome.NumSites() + Index.
	Index int

	// TargetSite is the Site.Index that this complex targets.
	TargetSite int

	// KOff is the unbinding rate constant (units: 1/time).
	// The complex stays bound for an exponentially-distributed duration
	// with mean 1/KOff before spontaneously dissociating.
	KOff float64

	// EnhWrite is the additional write (methylation) rate constant contributed
	// by this complex when it is bound at the target site.
	// Added to the background rate in the propensity function.
	EnhWrite float64

	// EnhErase is the additional erase (demethylation) rate constant contributed
	// by this complex when it is bound at the target site.
	// Added to the background rate in the propensity function.
	EnhErase float64
}
