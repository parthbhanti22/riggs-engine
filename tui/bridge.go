// Package tui implements the Bubble Tea terminal dashboard for the Riggs Engine.
//
// Architecture: the simulation runs on its own goroutine, stepping the SSA in
// place and updating a shared SimSnapshot behind a mutex. The Bubble Tea event
// loop reads snapshots at 10 FPS via a tick command — no individual events are
// ever pushed into the TUI message queue.
//
// This design solves the producer/consumer mismatch: the simulation fires
// thousands of events/sec while the terminal renders at 10 FPS. The snapshot
// mutex is held for ~200ns per update (memcpy of a 50-int slice + ring buffer
// append), so contention is negligible.
package tui

import (
	"math"
	"math/rand/v2"
	"sync"

	"github.com/pxrth9/riggs/bio"
	"github.com/pxrth9/riggs/gillespie"
)

// WALEntry records a single reaction event for the WAL tail view.
// This is the "write-ahead log" made visible — each entry is one line
// in the scrolling tail -f style log.
type WALEntry struct {
	Time          float64 // simulation time when the event fired
	FiredReaction int     // >=0 stochastic channel, <0 scheduled (-(tag+1))
	SiteAffected  int     // which site's methylation changed (-1 for complex-only)
	NewValue      int     // new value of the affected species after the event
}

// ringCapacity is the fixed size of the WAL ring buffer.
// 256 entries ≈ 8KB. At 50 sites with ~1 event/time-unit, stores ~5 sec of history.
const ringCapacity = 256

// RingBuffer is a fixed-size circular buffer of WALEntry values.
// When full, the oldest entry is overwritten by the newest.
// This is a hard requirement for the RAM budget — no unbounded slices.
type RingBuffer struct {
	entries [ringCapacity]WALEntry
	head    int // next write position (wraps around)
	count   int // number of valid entries (0..ringCapacity)
}

// Push appends an entry to the ring buffer, overwriting the oldest if full.
func (rb *RingBuffer) Push(e WALEntry) {
	rb.entries[rb.head] = e
	rb.head = (rb.head + 1) % ringCapacity
	if rb.count < ringCapacity {
		rb.count++
	}
}

// Len returns the number of valid entries in the buffer.
func (rb *RingBuffer) Len() int {
	return rb.count
}

// Get returns the i-th entry in chronological order (0 = oldest valid entry).
// Panics if i is out of range.
func (rb *RingBuffer) Get(i int) WALEntry {
	if i < 0 || i >= rb.count {
		panic("RingBuffer.Get: index out of range")
	}
	// The oldest entry is at (head - count + capacity) % capacity
	idx := (rb.head - rb.count + i + ringCapacity) % ringCapacity
	return rb.entries[idx]
}

// Tail returns the last n entries in chronological order (oldest first).
// If n > count, returns all valid entries.
func (rb *RingBuffer) Tail(n int) []WALEntry {
	if n > rb.count {
		n = rb.count
	}
	result := make([]WALEntry, n)
	start := rb.count - n
	for i := 0; i < n; i++ {
		result[i] = rb.Get(start + i)
	}
	return result
}

// SimSnapshot holds the coalesced simulation state shared between the
// simulation goroutine and the Bubble Tea event loop.
//
// The simulation goroutine writes to this struct after every SSA step.
// The TUI reads from it every 100ms (10 FPS tick). Both sides acquire
// the mutex, which is held for ~200-300ns per access.
//
// Invariant: MethFracs contains exact time-weighted methylation fractions
// computed by the simulation goroutine's running integrator. The TUI never
// computes its own averages — it reads these directly.
type SimSnapshot struct {
	mu sync.Mutex

	// Counts is a copy of the current State.Counts (N+K species).
	Counts []int

	// SimTime is the current simulation clock.
	SimTime float64

	// EventCount is the total number of reaction events fired.
	EventCount int64

	// Paused indicates whether the simulation is currently paused.
	Paused bool

	// Done indicates the simulation has completed (tMax reached or absorbing state).
	Done bool

	// WAL is the bounded ring buffer of recent events.
	WAL RingBuffer

	// MethFracs holds the exact time-weighted methylation fraction for each
	// CpG site. Updated by the simulation goroutine's running integrator.
	// Length = NumSites. These are exact, not sampled.
	MethFracs []float64

	// NumSites is the number of CpG sites (for partitioning Counts).
	NumSites int

	// NumComplexes is the number of targeting complexes.
	NumComplexes int
}

// CopyTo performs a deep copy of the snapshot into dst.
// The caller must hold s.mu (or call this from within a locked section).
func (s *SimSnapshot) CopyTo(dst *SimSnapshot) {
	dst.SimTime = s.SimTime
	dst.EventCount = s.EventCount
	dst.Paused = s.Paused
	dst.Done = s.Done
	dst.WAL = s.WAL // array copy (value type, not pointer)
	dst.NumSites = s.NumSites
	dst.NumComplexes = s.NumComplexes

	// Deep copy slices
	if len(dst.Counts) != len(s.Counts) {
		dst.Counts = make([]int, len(s.Counts))
	}
	copy(dst.Counts, s.Counts)

	if len(dst.MethFracs) != len(s.MethFracs) {
		dst.MethFracs = make([]float64, len(s.MethFracs))
	}
	copy(dst.MethFracs, s.MethFracs)
}

// SimRunner owns the simulation goroutine and the shared snapshot.
// It runs the hybrid SSA loop (stochastic + scheduled events) in a
// background goroutine and updates the snapshot after each step.
type SimRunner struct {
	// Simulation configuration (immutable after Start())
	reactions []gillespie.Reaction
	schedule  []gillespie.ScheduledEvent
	tMax      float64
	numSites  int
	numCplx   int
	genome    *bio.Genome

	// Mutable simulation state (owned by the simulation goroutine only)
	state     gillespie.State
	rng       *rand.Rand
	scratch   []float64
	schedIdx  int       // pointer into sorted schedule
	methInteg []float64 // per-site integral of methylation state over time
	lastTime  float64   // last time the integrator was updated

	// Shared state (protected by snapshot.mu)
	snapshot SimSnapshot

	// Control signals
	pauseCh  chan struct{} // send to request pause toggle
	stepCh   chan struct{} // send to request single-step (when paused)
	stopCh   chan struct{} // close to stop the simulation goroutine
	pauseMu  sync.Mutex
	isPaused bool
	pauseCnd *sync.Cond
}

// NewSimRunner creates a SimRunner from a bio.System build result.
func NewSimRunner(buildResult bio.BuildResult, genome *bio.Genome, tMax float64, seed uint64) *SimRunner {
	r := &SimRunner{
		reactions: buildResult.Reactions,
		schedule:  buildResult.Schedule,
		tMax:      tMax,
		numSites:  buildResult.NumSites,
		numCplx:   buildResult.NumComplexes,
		genome:    genome,
		state:     buildResult.InitialState.Clone(),
		rng:       rand.New(rand.NewPCG(seed, 0)),
		scratch:   make([]float64, len(buildResult.Reactions)),
		methInteg: make([]float64, buildResult.NumSites),
		lastTime:  0,
		pauseCh:   make(chan struct{}, 1),
		stepCh:    make(chan struct{}, 1),
		stopCh:    make(chan struct{}),
	}
	r.pauseCnd = sync.NewCond(&r.pauseMu)

	// Initialize snapshot
	r.snapshot.Counts = make([]int, len(r.state.Counts))
	copy(r.snapshot.Counts, r.state.Counts)
	r.snapshot.MethFracs = make([]float64, buildResult.NumSites)
	r.snapshot.NumSites = buildResult.NumSites
	r.snapshot.NumComplexes = buildResult.NumComplexes

	return r
}

// Start launches the simulation goroutine.
func (r *SimRunner) Start() {
	go r.runLoop()
}

// TogglePause toggles the pause state.
func (r *SimRunner) TogglePause() {
	r.pauseMu.Lock()
	r.isPaused = !r.isPaused
	paused := r.isPaused
	if !r.isPaused {
		r.pauseCnd.Signal()
	}
	r.pauseMu.Unlock()

	r.snapshot.mu.Lock()
	r.snapshot.Paused = paused
	r.snapshot.mu.Unlock()
}

// StepOnce advances the simulation by exactly one event (only when paused).
// Signals the condition variable to wake the goroutine from its Wait().
func (r *SimRunner) StepOnce() {
	select {
	case r.stepCh <- struct{}{}:
	default:
	}
	// Wake the goroutine so it can check stepCh
	r.pauseCnd.Signal()
}

// Stop terminates the simulation goroutine.
func (r *SimRunner) Stop() {
	close(r.stopCh)
	// Wake up if paused
	r.pauseMu.Lock()
	r.isPaused = false
	r.pauseCnd.Signal()
	r.pauseMu.Unlock()
}

// ReadSnapshot copies the current snapshot into dst.
// Thread-safe — acquires the snapshot mutex.
func (r *SimRunner) ReadSnapshot(dst *SimSnapshot) {
	r.snapshot.mu.Lock()
	r.snapshot.CopyTo(dst)
	r.snapshot.mu.Unlock()
}

// runLoop is the simulation goroutine's main loop.
// It implements the hybrid SSA: at each iteration, compare the next
// stochastic event time with the next scheduled event time, fire
// whichever is sooner. This is identical logic to RunWithSchedule()
// but without building a TrajectoryRecord.
func (r *SimRunner) runLoop() {
	var eventCount int64

	for r.state.Time < r.tMax {
		// --- Check stop signal ---
		select {
		case <-r.stopCh:
			return
		default:
		}

		// --- Check pause state ---
		r.pauseMu.Lock()
		for r.isPaused {
			// Check for step-once requests before and after waiting
			select {
			case <-r.stepCh:
				r.pauseMu.Unlock()
				r.doOneStep(&eventCount)
				r.pauseMu.Lock()
				continue
			default:
			}
			// Block until signaled (by TogglePause, StepOnce, or Stop)
			r.pauseCnd.Wait()
			// After waking, check stop
			select {
			case <-r.stopCh:
				r.pauseMu.Unlock()
				return
			default:
			}
			// After waking, check step-once (may have been signaled by StepOnce())
			select {
			case <-r.stepCh:
				r.pauseMu.Unlock()
				r.doOneStep(&eventCount)
				r.pauseMu.Lock()
				continue
			default:
			}
		}
		r.pauseMu.Unlock()

		r.doOneStep(&eventCount)
	}

	// Simulation complete — mark as done
	r.snapshot.mu.Lock()
	r.snapshot.Done = true
	r.snapshot.mu.Unlock()
}

// doOneStep executes one hybrid SSA step (stochastic or scheduled event).
func (r *SimRunner) doOneStep(eventCount *int64) {
	// Compute next stochastic event time
	a0 := 0.0
	for i, rx := range r.reactions {
		r.scratch[i] = rx.Propensity(r.state.Counts)
		a0 += r.scratch[i]
	}

	var stochTime float64
	if a0 <= 0 {
		stochTime = math.Inf(1)
	} else {
		r1 := 1.0 - r.rng.Float64()
		stochTime = r.state.Time + (-math.Log(r1) / a0)
	}

	// Check next scheduled event
	var nextSched *gillespie.ScheduledEvent
	if r.schedIdx < len(r.schedule) {
		nextSched = &r.schedule[r.schedIdx]
	}

	// Determine which fires first
	if nextSched != nil && nextSched.Time <= stochTime && nextSched.Time <= r.tMax {
		// --- Scheduled event fires ---
		r.updateIntegrator(nextSched.Time)
		r.state.Time = nextSched.Time

		// Apply deltas with clamping
		siteAffected := -1
		newValue := 0
		for s, delta := range nextSched.Deltas {
			if delta != 0 {
				nv := r.state.Counts[s] + delta
				if nv < 0 {
					nv = 0
				} else if nv > 1 {
					nv = 1
				}
				r.state.Counts[s] = nv
				siteAffected = s
				newValue = nv
			}
		}

		firedTag := -(nextSched.Tag + 1)
		r.schedIdx++
		*eventCount++

		r.pushSnapshot(*eventCount, firedTag, siteAffected, newValue)

	} else if stochTime <= r.tMax {
		// --- Stochastic event fires ---
		r.updateIntegrator(stochTime)
		r.state.Time = stochTime

		// Select which reaction fires
		r2 := r.rng.Float64()
		target := r2 * a0
		cum := 0.0
		fired := len(r.reactions) - 1
		for i := range r.reactions {
			cum += r.scratch[i]
			if cum > target {
				fired = i
				break
			}
		}

		// Apply deltas and find affected site
		siteAffected := -1
		newValue := 0
		for s, delta := range r.reactions[fired].Deltas {
			if delta != 0 {
				r.state.Counts[s] += delta
				siteAffected = s
				newValue = r.state.Counts[s]
			}
		}

		*eventCount++
		r.pushSnapshot(*eventCount, fired, siteAffected, newValue)

	} else {
		// Both past tMax — advance time to tMax and finish
		r.updateIntegrator(r.tMax)
		r.state.Time = r.tMax
	}
}

// updateIntegrator updates the per-site running integral of methylation state.
// This is how we maintain exact time-weighted fractions without recording
// every event — the integrator accumulates (state × duration) continuously.
func (r *SimRunner) updateIntegrator(newTime float64) {
	dt := newTime - r.lastTime
	if dt <= 0 {
		return
	}
	for i := 0; i < r.numSites; i++ {
		r.methInteg[i] += float64(r.state.Counts[i]) * dt
	}
	r.lastTime = newTime
}

// pushSnapshot updates the shared snapshot with the current state.
// Called after every event from the simulation goroutine.
func (r *SimRunner) pushSnapshot(eventCount int64, firedReaction, siteAffected, newValue int) {
	r.snapshot.mu.Lock()
	copy(r.snapshot.Counts, r.state.Counts)
	r.snapshot.SimTime = r.state.Time
	r.snapshot.EventCount = eventCount

	// Update exact methylation fractions
	if r.state.Time > 0 {
		for i := 0; i < r.numSites; i++ {
			r.snapshot.MethFracs[i] = r.methInteg[i] / r.state.Time
		}
	}

	// Append to WAL ring buffer
	r.snapshot.WAL.Push(WALEntry{
		Time:          r.state.Time,
		FiredReaction: firedReaction,
		SiteAffected:  siteAffected,
		NewValue:      newValue,
	})

	r.snapshot.mu.Unlock()
}
