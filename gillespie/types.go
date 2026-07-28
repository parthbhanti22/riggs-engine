// Package gillespie implements the Stochastic Simulation Algorithm (SSA)
// Direct Method (Gillespie 1977) for exact Monte Carlo simulation of
// coupled chemical reaction systems modeled as continuous-time Markov chains.
//
// Design priorities for an 8GB RAM, 4-core, CPU-only target:
//   - Flat int slices for species counts (struct-of-arrays, not slice-of-pointers)
//   - Func-typed propensity fields to avoid interface dispatch in the inner loop
//   - Pre-allocatable scratch buffers to eliminate per-step heap allocations
package gillespie

// Reaction represents a single reaction channel in the SSA.
//
// In the Gillespie framework, a reaction channel is defined by:
//   - A stoichiometry vector (Deltas): how species counts change when this reaction fires.
//   - A stochastic rate constant (RateConst, often written c_i): encodes the intrinsic
//     probability per unit time that a particular combination of reactant molecules will react.
//   - A propensity function (Propensity, often written a_i(x)): the rate at which this
//     reaction fires given the current state x. For mass-action kinetics this is
//     c_i * (combinatorial product of reactant counts), but we keep it as a general
//     function to support non-mass-action channels in Phase 2.
type Reaction struct {
	// Deltas is the state-change vector: Deltas[s] is the change to species s
	// when this reaction fires. Length must equal the number of species.
	// Example: methylation of a single CpG site → Deltas = [+1].
	Deltas []int

	// RateConst is the stochastic rate constant c_i.
	// Units: 1/time for unimolecular, 1/(time·count) for bimolecular, etc.
	// For the CpG toy system, this is k_write or k_erase directly.
	RateConst float64

	// Propensity computes a_i(x) given the current species count vector.
	// This function MUST be pure (no side effects, no captured mutable state)
	// and MUST NOT modify the passed-in counts slice.
	//
	// We use a func field rather than an interface to avoid virtual dispatch
	// overhead in the inner SSA loop. The closure captures only the rate constant
	// (immutable) — the counts slice is passed per call.
	Propensity func(counts []int) float64
}

// State holds the mutable simulation state: the species population vector
// and the continuous simulation clock.
//
// The species vector is a flat []int indexed by species ID (0-based).
// For the CpG toy system: Counts[0] = 1 means methylated, 0 means unmethylated.
type State struct {
	// Counts is the species population vector x(t).
	// Flat int slice — no pointers, no maps, no boxing.
	Counts []int

	// Time is the current simulation clock in continuous time.
	// Advances by exponentially-distributed increments τ at each step.
	Time float64
}

// Clone returns a deep copy of the State. Used to give each trajectory
// its own independent starting state without aliasing the Counts slice.
func (s State) Clone() State {
	c := make([]int, len(s.Counts))
	copy(c, s.Counts)
	return State{Counts: c, Time: s.Time}
}

// TrajectoryRecord stores the event stream from a single SSA trajectory.
// Each entry corresponds to one reaction firing event.
//
// We record sparse events (only when state changes), not fixed-timestep
// snapshots — this is the natural output of the SSA and avoids interpolation.
type TrajectoryRecord struct {
	// Times contains the simulation clock value after each reaction event.
	Times []float64

	// States contains a snapshot of the species count vector after each event.
	// States[i] corresponds to the state at Times[i].
	// Each inner slice is an independent copy (no aliasing).
	States [][]int

	// FiredReaction contains the index of the reaction that fired at each event.
	FiredReaction []int
}

// EnsembleResult aggregates statistics over many independent trajectories
// of the same system.
type EnsembleResult struct {
	// MeanFraction is the mean of the per-trajectory time-weighted observable.
	// For the CpG toy system: mean fraction of time spent methylated.
	MeanFraction float64

	// StdDev is the sample standard deviation of per-trajectory fractions.
	StdDev float64

	// ConfLo and ConfHi are the 95% confidence interval bounds on MeanFraction.
	// Computed as mean ± 1.96 * stddev / sqrt(N).
	ConfLo, ConfHi float64

	// N is the number of trajectories that contributed to this result.
	N int
}
