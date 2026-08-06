package tui

import (
	"math"
	"testing"
	"time"

	"github.com/pxrth9/riggs/bio"
)

// --- RingBuffer tests ---

func TestRingBuffer_PushAndLen(t *testing.T) {
	var rb RingBuffer
	if rb.Len() != 0 {
		t.Fatalf("empty buffer: Len() = %d, want 0", rb.Len())
	}

	for i := 0; i < 10; i++ {
		rb.Push(WALEntry{Time: float64(i)})
	}
	if rb.Len() != 10 {
		t.Fatalf("after 10 pushes: Len() = %d, want 10", rb.Len())
	}
}

func TestRingBuffer_Overflow(t *testing.T) {
	var rb RingBuffer

	// Fill beyond capacity — oldest should be overwritten
	for i := 0; i < ringCapacity+50; i++ {
		rb.Push(WALEntry{Time: float64(i)})
	}

	if rb.Len() != ringCapacity {
		t.Fatalf("after overflow: Len() = %d, want %d", rb.Len(), ringCapacity)
	}

	// Oldest should be entry 50 (entries 0-49 overwritten)
	oldest := rb.Get(0)
	if oldest.Time != 50.0 {
		t.Errorf("oldest entry: Time=%.0f, want 50", oldest.Time)
	}

	// Newest should be entry (capacity + 49)
	newest := rb.Get(rb.Len() - 1)
	expected := float64(ringCapacity + 49)
	if newest.Time != expected {
		t.Errorf("newest entry: Time=%.0f, want %.0f", newest.Time, expected)
	}
}

func TestRingBuffer_Tail(t *testing.T) {
	var rb RingBuffer
	for i := 0; i < 20; i++ {
		rb.Push(WALEntry{Time: float64(i), SiteAffected: i})
	}

	tail := rb.Tail(5)
	if len(tail) != 5 {
		t.Fatalf("Tail(5): len=%d, want 5", len(tail))
	}

	// Should be entries 15, 16, 17, 18, 19
	for i, e := range tail {
		expected := float64(15 + i)
		if e.Time != expected {
			t.Errorf("Tail[%d]: Time=%.0f, want %.0f", i, e.Time, expected)
		}
	}
}

func TestRingBuffer_TailMoreThanCount(t *testing.T) {
	var rb RingBuffer
	rb.Push(WALEntry{Time: 1.0})
	rb.Push(WALEntry{Time: 2.0})

	tail := rb.Tail(100) // request more than available
	if len(tail) != 2 {
		t.Fatalf("Tail(100) with 2 entries: len=%d, want 2", len(tail))
	}
}

func TestRingBuffer_ChronologicalOrder(t *testing.T) {
	var rb RingBuffer
	// Push more than capacity to test wrap-around ordering
	for i := 0; i < ringCapacity*3; i++ {
		rb.Push(WALEntry{Time: float64(i)})
	}

	// All entries should be in chronological order
	for i := 1; i < rb.Len(); i++ {
		prev := rb.Get(i - 1)
		curr := rb.Get(i)
		if curr.Time <= prev.Time {
			t.Errorf("entry %d (t=%.0f) <= entry %d (t=%.0f) — not chronological",
				i, curr.Time, i-1, prev.Time)
			break
		}
	}
}

// --- SimSnapshot tests ---

func TestSimSnapshot_CopyTo(t *testing.T) {
	src := &SimSnapshot{
		Counts:       []int{1, 0, 1, 0, 1},
		SimTime:      123.456,
		EventCount:   999,
		Paused:       true,
		Done:         false,
		NumSites:     4,
		NumComplexes: 1,
		MethFracs:    []float64{0.5, 0.1, 0.9, 0.3},
	}
	src.WAL.Push(WALEntry{Time: 100.0, SiteAffected: 2})

	var dst SimSnapshot
	src.CopyTo(&dst)

	// Verify deep copy
	if dst.SimTime != 123.456 {
		t.Errorf("SimTime: got %f, want 123.456", dst.SimTime)
	}
	if dst.EventCount != 999 {
		t.Errorf("EventCount: got %d, want 999", dst.EventCount)
	}
	if !dst.Paused {
		t.Error("Paused: got false, want true")
	}

	// Modify source — destination should be independent
	src.Counts[0] = 99
	src.MethFracs[0] = 99.9
	if dst.Counts[0] == 99 {
		t.Error("Counts not deep-copied — still aliased to source")
	}
	if dst.MethFracs[0] == 99.9 {
		t.Error("MethFracs not deep-copied — still aliased to source")
	}

	// WAL should be copied
	if dst.WAL.Len() != 1 {
		t.Errorf("WAL.Len(): got %d, want 1", dst.WAL.Len())
	}
}

// --- SimRunner tests ---

func makeTestSystem() (*bio.System, bio.BuildResult) {
	g := bio.NewGenome(5, 100, nil)
	sys := &bio.System{
		Genome: g,
		Complexes: []bio.TargetingComplex{
			{Index: 0, TargetSite: 0, KOff: 0.1, EnhWrite: 0.5, EnhErase: 0.05},
		},
		Triggers: []bio.EnvironmentalTrigger{
			{Name: "test", ComplexIdx: 0, FireTimes: []float64{10.0, 50.0}},
		},
		KBgWrite: 0.001,
		KBgErase: 0.01,
	}
	return sys, sys.Build()
}

func TestSimRunner_StartAndStop(t *testing.T) {
	sys, result := makeTestSystem()
	runner := NewSimRunner(result, sys.Genome, sys.Complexes, 100.0, 42)
	runner.Start()

	// Wait briefly for some events to fire
	time.Sleep(50 * time.Millisecond)

	var snap SimSnapshot
	runner.ReadSnapshot(&snap)

	if snap.EventCount == 0 {
		t.Error("no events fired after 50ms — simulation may not be running")
	}
	if snap.SimTime <= 0 {
		t.Error("simulation time not advancing")
	}

	runner.Stop()
}

func TestSimRunner_PauseResume(t *testing.T) {
	sys, result := makeTestSystem()
	// Use very large tMax so the simulation is definitely still running
	// during the pause/resume window (with these rates, tMax=1000 completes
	// in ~1ms — far too fast for the test's timing windows).
	runner := NewSimRunner(result, sys.Genome, sys.Complexes, 1e9, 42)
	runner.Start()

	// Let it run briefly
	time.Sleep(20 * time.Millisecond)

	// Pause
	runner.TogglePause()
	time.Sleep(10 * time.Millisecond)

	var snap1 SimSnapshot
	runner.ReadSnapshot(&snap1)

	// Wait while paused — event count should not change
	time.Sleep(50 * time.Millisecond)

	var snap2 SimSnapshot
	runner.ReadSnapshot(&snap2)

	if snap2.EventCount != snap1.EventCount {
		t.Errorf("events fired while paused: %d → %d", snap1.EventCount, snap2.EventCount)
	}
	if !snap2.Paused {
		t.Error("snapshot should show Paused=true")
	}

	// Resume
	runner.TogglePause()
	time.Sleep(100 * time.Millisecond)

	var snap3 SimSnapshot
	runner.ReadSnapshot(&snap3)

	if snap3.EventCount <= snap2.EventCount {
		t.Error("no events fired after resume")
	}

	runner.Stop()
}

func TestSimRunner_StepOnce(t *testing.T) {
	sys, result := makeTestSystem()
	runner := NewSimRunner(result, sys.Genome, sys.Complexes, 1000.0, 42)

	// Start paused
	runner.Start()
	runner.TogglePause()
	time.Sleep(10 * time.Millisecond)

	var before SimSnapshot
	runner.ReadSnapshot(&before)

	// Step once
	runner.StepOnce()
	time.Sleep(10 * time.Millisecond)

	var after SimSnapshot
	runner.ReadSnapshot(&after)

	if after.EventCount != before.EventCount+1 {
		t.Errorf("StepOnce: events went from %d to %d, expected +1",
			before.EventCount, after.EventCount)
	}

	runner.Stop()
}

func TestSimRunner_MethFracsExact(t *testing.T) {
	// Validate that the running integrator in SimRunner converges to the
	// analytical CTMC steady-state. We run multiple short trajectories
	// and check the mean per-site fraction.
	//
	// For background-only: P(methylated) = k_w / (k_w + k_e) = 0.01/0.03 ≈ 0.3333
	kW := 0.01
	kE := 0.02
	expected := kW / (kW + kE) // 0.3333

	g := bio.NewGenome(5, 100, nil)
	sys := &bio.System{
		Genome:   g,
		KBgWrite: kW,
		KBgErase: kE,
	}
	result := sys.Build()

	nRuns := 20
	tMax := 2000.0
	siteFracSum := make([]float64, 5)

	for run := 0; run < nRuns; run++ {
		runner := NewSimRunner(result, g, nil, tMax, uint64(run))
		runner.Start()

		// Wait for completion
		for {
			time.Sleep(20 * time.Millisecond)
			var snap SimSnapshot
			runner.ReadSnapshot(&snap)
			if snap.Done {
				for site := 0; site < 5; site++ {
					siteFracSum[site] += snap.MethFracs[site]
				}
				break
			}
		}
		runner.Stop()
	}

	// Check mean across runs
	for site := 0; site < 5; site++ {
		mean := siteFracSum[site] / float64(nRuns)
		if math.Abs(mean-expected) > 0.05 {
			t.Errorf("site %d: mean MethFrac=%.4f, expected=%.4f (deviation > 0.05)",
				site, mean, expected)
		}
	}
}

func TestSimRunner_WALPopulated(t *testing.T) {
	sys, result := makeTestSystem()
	runner := NewSimRunner(result, sys.Genome, sys.Complexes, 100.0, 42)
	runner.Start()

	time.Sleep(50 * time.Millisecond)

	var snap SimSnapshot
	runner.ReadSnapshot(&snap)

	if snap.WAL.Len() == 0 {
		t.Error("WAL ring buffer is empty after simulation ran")
	}

	// Entries should be in chronological order
	for i := 1; i < snap.WAL.Len(); i++ {
		prev := snap.WAL.Get(i - 1)
		curr := snap.WAL.Get(i)
		if curr.Time < prev.Time {
			t.Errorf("WAL entry %d (t=%.4f) < entry %d (t=%.4f) — not chronological",
				i, curr.Time, i-1, prev.Time)
			break
		}
	}

	runner.Stop()
}

func TestSimRunner_ScheduledEventsInterleaved(t *testing.T) {
	// Verify that scheduled events (triggers) appear in the WAL with
	// negative FiredReaction values.
	sys, result := makeTestSystem()
	// Trigger fires at t=10.0 and t=50.0
	runner := NewSimRunner(result, sys.Genome, sys.Complexes, 100.0, 42)
	runner.Start()

	// Wait for simulation to complete
	for {
		time.Sleep(50 * time.Millisecond)
		var snap SimSnapshot
		runner.ReadSnapshot(&snap)
		if snap.Done || snap.SimTime > 60.0 {
			// Check for scheduled events in WAL
			foundScheduled := false
			for i := 0; i < snap.WAL.Len(); i++ {
				entry := snap.WAL.Get(i)
				if entry.FiredReaction < 0 {
					foundScheduled = true
					break
				}
			}
			if !foundScheduled {
				t.Error("no scheduled events found in WAL — triggers may not be interleaved")
			}
			break
		}
	}

	runner.Stop()
}
