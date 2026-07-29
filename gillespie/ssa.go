package gillespie

import (
	"math"
	"math/rand/v2"
)

// Step executes one step of the SSA Direct Method (Gillespie 1977).
//
// Given the current state and a set of reaction channels, it:
//  1. Computes the propensity a_i for each reaction channel.
//  2. Sums them to get total propensity a0.
//  3. Draws the time to next event from an exponential distribution: τ = -ln(r1)/a0.
//  4. Selects which reaction fires by partitioning [0, a0) by cumulative propensity.
//  5. Applies the selected reaction's state-change vector.
//
// Parameters:
//   - state: mutable simulation state (modified in place).
//   - reactions: the reaction channels (read-only).
//   - rng: a seedable PRNG — must not be shared across goroutines.
//   - scratch: pre-allocated buffer for propensity values. Must have len >= len(reactions).
//     Passing a pre-allocated scratch buffer eliminates the main source of per-step
//     heap allocation in the inner loop. The ensemble runner allocates one per worker.
//
// Returns:
//   - fired: index of the reaction that fired, or -1 if the system is in an absorbing
//     state (all propensities are zero — no further reactions can occur).
//   - dt: the time elapsed until this event. +Inf if absorbing.
func Step(state *State, reactions []Reaction, rng *rand.Rand, scratch []float64) (fired int, dt float64) {
	// --- 1. Compute propensities ---
	// For each reaction channel i, a_i(x) = rate_i * f(current_counts).
	// Total propensity a0 = Σ a_i is the overall rate of "something happens".
	a0 := 0.0
	for i, r := range reactions {
		scratch[i] = r.Propensity(state.Counts)
		a0 += scratch[i]
	}

	// --- Absorbing state check ---
	// If a0 = 0, no reaction can fire. The system is stuck.
	// For the CpG toy system this can't happen (one of the two channels is always
	// active), but it's essential for correctness in general systems.
	if a0 <= 0 {
		return -1, math.Inf(1)
	}

	// --- 2. Time to next event ---
	// The waiting time τ between events in a Poisson process with rate a0 is
	// exponentially distributed: P(τ > t) = exp(-a0 * t).
	//
	// To sample from Exp(a0), we use the inverse-CDF method:
	//   If U ~ Uniform(0,1), then τ = -ln(U) / a0 ~ Exp(a0).
	//
	// Subtlety: rand.Float64() returns values in [0, 1). If it returns exactly 0,
	// ln(0) = -Inf, which would produce an infinite time step — nonsensical.
	// We use (1.0 - rng.Float64()) to map [0,1) → (0,1], guaranteeing r1 > 0.
	// This is a standard numerical hygiene trick in SSA implementations.
	r1 := 1.0 - rng.Float64()
	dt = -math.Log(r1) / a0

	// --- 3. Select which reaction fires ---
	// We partition [0, a0) into intervals [0, a_0), [a_0, a_0+a_1), ... and
	// pick a uniform random point r2*a0 in [0, a0). The interval it falls in
	// determines which reaction fires. This gives each reaction a probability
	// proportional to its propensity — exactly as the theory requires.
	//
	// Floating-point edge case: due to rounding, the cumulative sum might not
	// quite reach a0, so r2*a0 could fall past all intervals. The fallback
	// assignment (fired = last reaction) handles this gracefully.
	r2 := rng.Float64()
	target := r2 * a0
	cum := 0.0
	fired = len(reactions) - 1 // fallback for floating-point edge case
	for i := range reactions {
		cum += scratch[i]
		if cum > target {
			fired = i
			break
		}
	}

	// --- 4. Update state ---
	// Apply the stoichiometry vector: for each species s, counts[s] += delta[s].
	for s, delta := range reactions[fired].Deltas {
		state.Counts[s] += delta
	}
	state.Time += dt

	return fired, dt
}

// Run executes a complete SSA trajectory from the given initial state until
// the simulation clock reaches tMax or the system enters an absorbing state.
//
// It records every reaction event in a TrajectoryRecord (sparse event stream).
// The initial state is NOT recorded — only post-event snapshots.
//
// Parameters:
//   - initial: the starting state (cloned internally — caller's copy is not modified).
//   - reactions: the reaction channels.
//   - tMax: maximum simulation time.
//   - rng: seedable PRNG for this trajectory.
//
// Returns a TrajectoryRecord containing the full event history.
func Run(initial State, reactions []Reaction, tMax float64, rng *rand.Rand) TrajectoryRecord {
	state := initial.Clone()
	scratch := make([]float64, len(reactions)) // one allocation for the whole trajectory

	// Pre-allocate with a reasonable capacity estimate.
	// For the two-state CpG system with total rate ~1.0, expect ~tMax events.
	// Over-estimating slightly is cheaper than repeated slice growth.
	estEvents := int(tMax * 1.5)
	if estEvents < 64 {
		estEvents = 64
	}
	rec := TrajectoryRecord{
		Times:         make([]float64, 0, estEvents),
		States:        make([][]int, 0, estEvents),
		FiredReaction: make([]int, 0, estEvents),
	}

	for state.Time < tMax {
		fired, _ := Step(&state, reactions, rng, scratch)
		if fired < 0 {
			break // absorbing state
		}
		// Don't record events past tMax — the last Step may have pushed
		// the clock beyond tMax. We truncate to keep the trajectory clean.
		if state.Time > tMax {
			break
		}

		rec.Times = append(rec.Times, state.Time)
		// Snapshot the counts — must copy to avoid aliasing.
		snap := make([]int, len(state.Counts))
		copy(snap, state.Counts)
		rec.States = append(rec.States, snap)
		rec.FiredReaction = append(rec.FiredReaction, fired)
	}

	return rec
}

// TimeWeightedFraction computes the fraction of total time that species at
// the given index was in a nonzero state, using time-weighted integration
// over the trajectory's event stream.
//
// For the CpG toy system (species 0, count ∈ {0,1}), this gives the fraction
// of time spent methylated — the observable we validate against the analytical
// result k_write / (k_write + k_erase).
//
// The calculation works by integrating the step function defined by the event
// stream: between events, the state is constant, so the integral is just
// (value × duration) summed over all intervals.
//
// Parameters:
//   - rec: the trajectory event record.
//   - speciesIdx: which species to measure.
//   - initialCount: the count of that species at t=0 (before any events).
//   - tMax: the end time of the simulation window.
func TimeWeightedFraction(rec TrajectoryRecord, speciesIdx int, initialCount int, tMax float64) float64 {
	if tMax <= 0 {
		return 0
	}

	totalWeighted := 0.0
	prevTime := 0.0
	currentCount := initialCount

	for i, t := range rec.Times {
		// Accumulate time spent in current state
		dt := t - prevTime
		if currentCount > 0 {
			totalWeighted += dt
		}
		// Transition to new state
		currentCount = rec.States[i][speciesIdx]
		prevTime = t
	}

	// Account for the final interval: from last event to tMax
	dt := tMax - prevTime
	if currentCount > 0 {
		totalWeighted += dt
	}

	return totalWeighted / tMax
}

// RunWithSchedule executes a complete SSA trajectory with deterministic
// scheduled events interleaved with stochastic steps (hybrid SSA).
//
// At each iteration, the algorithm compares:
//   - The next stochastic event time (τ, from the SSA step)
//   - The next scheduled event time
//
// Whichever is sooner fires first. If the scheduled event fires, its Deltas
// are applied with clamping (binary species are clamped to [0,1] to enforce
// idempotency — a trigger firing while its complex is already bound is a no-op).
//
// This is a standard technique for handling deterministic events in SSA
// simulations (Anderson 2007). It preserves the exactness of the SSA for
// stochastic channels while correctly interleaving deterministic interventions.
//
// The schedule must be sorted by Time in ascending order.
//
// Parameters:
//   - initial: starting state (cloned internally).
//   - reactions: stochastic reaction channels.
//   - schedule: deterministic events sorted by time.
//   - tMax: maximum simulation time.
//   - rng: seedable PRNG.
//
// FiredReaction values in the returned TrajectoryRecord use negative indices
// for scheduled events: -(tag + 1), where tag is the ScheduledEvent.Tag.
// This distinguishes them from stochastic reaction indices (which are >= 0).
func RunWithSchedule(initial State, reactions []Reaction, schedule []ScheduledEvent,
	tMax float64, rng *rand.Rand) TrajectoryRecord {

	state := initial.Clone()
	scratch := make([]float64, len(reactions))
	schedIdx := 0 // pointer into the sorted schedule

	// Pre-allocate record with a reasonable capacity estimate.
	estEvents := int(tMax*1.5) + len(schedule)
	if estEvents < 64 {
		estEvents = 64
	}
	rec := TrajectoryRecord{
		Times:         make([]float64, 0, estEvents),
		States:        make([][]int, 0, estEvents),
		FiredReaction: make([]int, 0, estEvents),
	}

	for state.Time < tMax {
		// --- Compute next stochastic event ---
		// We need to know the stochastic τ WITHOUT applying it yet,
		// so we can compare with the next scheduled event time.
		a0 := 0.0
		for i, r := range reactions {
			scratch[i] = r.Propensity(state.Counts)
			a0 += scratch[i]
		}

		var stochTime float64
		if a0 <= 0 {
			stochTime = math.Inf(1) // no stochastic events possible
		} else {
			r1 := 1.0 - rng.Float64()
			stochTime = state.Time + (-math.Log(r1) / a0)
		}

		// --- Check next scheduled event ---
		var nextSched *ScheduledEvent
		if schedIdx < len(schedule) {
			nextSched = &schedule[schedIdx]
		}

		// --- Determine which fires first ---
		if nextSched != nil && nextSched.Time <= stochTime && nextSched.Time <= tMax {
			// Scheduled event fires first.
			state.Time = nextSched.Time

			// Apply deltas with clamping for binary species.
			// Clamping prevents counts from going below 0 or above 1 for
			// binary-state species (methylation, bound state).
			for s, delta := range nextSched.Deltas {
				if delta != 0 {
					newVal := state.Counts[s] + delta
					if newVal < 0 {
						newVal = 0
					} else if newVal > 1 {
						newVal = 1
					}
					state.Counts[s] = newVal
				}
			}
			schedIdx++

			// Record this event with a negative tag to distinguish
			// from stochastic reactions.
			if state.Time <= tMax {
				rec.Times = append(rec.Times, state.Time)
				snap := make([]int, len(state.Counts))
				copy(snap, state.Counts)
				rec.States = append(rec.States, snap)
				rec.FiredReaction = append(rec.FiredReaction, -(nextSched.Tag + 1))
			}
		} else if stochTime <= tMax {
			// Stochastic event fires first — standard SSA selection.
			state.Time = stochTime

			// Select which reaction fires.
			r2 := rng.Float64()
			target := r2 * a0
			cum := 0.0
			fired := len(reactions) - 1
			for i := range reactions {
				cum += scratch[i]
				if cum > target {
					fired = i
					break
				}
			}

			// Apply deltas.
			for s, delta := range reactions[fired].Deltas {
				state.Counts[s] += delta
			}

			rec.Times = append(rec.Times, state.Time)
			snap := make([]int, len(state.Counts))
			copy(snap, state.Counts)
			rec.States = append(rec.States, snap)
			rec.FiredReaction = append(rec.FiredReaction, fired)
		} else {
			// Both next stochastic and next scheduled are past tMax — done.
			break
		}
	}

	return rec
}
