package bio

// EnvironmentalTrigger represents an external stimulus that instantiates a
// TargetingComplex at its target site at predetermined times.
//
// This is the "write-ahead log" mechanic: when the trigger fires, it binds
// a targeting complex to a CpG site, which enhances methylation/demethylation
// activity at that site. The trigger's occurrence becomes chronologically
// recorded as a methylation mark — the biological equivalent of appending
// an entry to a write-ahead log.
//
// The trigger uses a fixed external schedule (deterministic times), not a
// stochastic Poisson process. This models a researcher or engineered system
// applying a stimulus at known times. In the Gillespie simulation, triggers
// are implemented as scheduled interrupts interleaved with stochastic SSA
// steps (hybrid SSA technique).
//
// Idempotency: if a trigger fires while its complex is already bound, the
// event is a no-op (the bound state stays at 1). This means the WAL records
// "attempted triggers chronologically" — not "guaranteed fresh bindings."
type EnvironmentalTrigger struct {
	// Name is a human-readable label for this trigger (e.g., "heat_shock").
	Name string

	// ComplexIdx is the index of the TargetingComplex this trigger activates.
	// When the trigger fires, it sets State.Counts[genome.NumSites + ComplexIdx] = 1.
	ComplexIdx int

	// FireTimes is the deterministic schedule: absolute simulation times at
	// which this trigger fires. Must be sorted in ascending order.
	FireTimes []float64
}
