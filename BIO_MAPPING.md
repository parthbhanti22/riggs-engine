# Riggs Engine — Biology-to-Code Mapping

> This document maps biological concepts to their Go representations in the Riggs
> codebase. It exists because the whole point of this project is learning the
> biology-to-code translation — not just getting simulations to run.
>
> **Key invariant**: The DNA sequence is never mutated. Only the methylation state
> vector changes. This is the ROM/RAM distinction that defines the project.

## Phase 1: Single CpG Site (Toy System)

The simplest possible epigenetic memory: one CpG dinucleotide that is either
methylated (bit = 1) or unmethylated (bit = 0). This is a two-state continuous-time
Markov chain (CTMC).

### State Mapping

| Biology | Go Representation | Notes |
|---------|------------------|-------|
| CpG site methylation state | `State.Counts[0]` | `0` = unmethylated, `1` = methylated. Single species, count is always 0 or 1 because we model *one* site. |
| "Unmethylated" | `Counts[0] == 0` | The default / ground state. |
| "Methylated" | `Counts[0] == 1` | A methyl group (CH₃) is attached to the cytosine. |
| Simulation clock | `State.Time` | Continuous time in arbitrary units. Advances by exponentially-distributed intervals. |

### Reaction Channel Mapping

| Biology | Go `Reaction` | Propensity `a_i(x)` | Rate Constant |
|---------|--------------|---------------------|---------------|
| **Methylation** by methyltransferase (DNMT-like enzyme) | `Deltas: [+1]`, `RateConst: k_write` | `k_write * (1 - Counts[0])` — fires only when unmethylated | `k_write = 0.3` (default) |
| **Demethylation** by demethylase (TET-like enzyme) | `Deltas: [-1]`, `RateConst: k_erase` | `k_erase * Counts[0]` — fires only when methylated | `k_erase = 0.7` (default) |

### Why the propensity functions look like that

In the Gillespie framework, propensity = rate constant × number of ways the reaction
can occur given the current state. For this two-state system:

- **Methylation**: can only happen if the site is unmethylated. "Number of ways" =
  `(1 - Counts[0])` — this is 1 when unmethylated (count=0) and 0 when already
  methylated (count=1). So `a_write = k_write * (1 - Counts[0])`.

- **Demethylation**: can only happen if the site is methylated. "Number of ways" =
  `Counts[0]` — this is 1 when methylated and 0 when unmethylated. So
  `a_erase = k_erase * Counts[0]`.

This is a degenerate case of mass-action kinetics where counts are always 0 or 1.
In Phase 2 with multiple CpG sites, counts can be larger and the combinatorial
factors become real combinatorics.

### Analytical Validation Target

This two-state system is exactly solvable. The long-run (steady-state) fraction of
time spent in the methylated state is:

```
P(methylated) = k_write / (k_write + k_erase)
```

With `k_write = 0.3`, `k_erase = 0.7`:

```
P(methylated) = 0.3 / (0.3 + 0.7) = 0.3
```

This is the target our ensemble average must converge to. The unit tests validate
this with a statistical tolerance of ±0.02 over 5000 trajectories.

### Biological Analogy Table

| Computer Analogy | Biology | Phase 1 Code |
|-----------------|---------|-------------|
| ROM / Disk | DNA sequence (immutable) | Not modeled in Phase 1 — no genome struct yet |
| RAM bit | Methylation state at one CpG | `State.Counts[0]` ∈ {0, 1} |
| Write operation | Methyltransferase attaches CH₃ | Reaction channel 0 (Deltas: [+1]) |
| Erase operation | Demethylase removes CH₃ | Reaction channel 1 (Deltas: [-1]) |
| Pointer / address | dCas9 navigating to target site | Not modeled in Phase 1 |
| Write-ahead log | Chronological methylation events | `TrajectoryRecord.Times` + `FiredReaction` |

## Phase 2: Multi-Site Genome with Targeted Methylation

Phase 2 extends the single-site toy to a genome with N CpG sites, programmable
targeting via dCas9 complexes, and environmentally-triggered methylation events.

### Entity Mapping

| Biology | Go Type | Package | Notes |
|---------|---------|---------|-------|
| Genome (DNA sequence) | `Genome` struct | `bio` | Ordered, immutable collection of `Site` structs — the "disk" |
| CpG dinucleotide | `Site` struct | `bio` | Index, coordinate, context tag. Methylation state tracked externally in `State.Counts[site.Index]` |
| Genomic context | `SiteContext` enum | `bio` | `ContextPromoter`, `ContextGeneBody`, `ContextIntergenic`. Phase 2: metadata only (does not affect rates) |
| dCas9 pointer | `TargetingComplex` struct | `bio` | Catalytically dead Cas9 + guide RNA targeting one site. Bound state = species in `State.Counts` |
| Methyltransferase effector | `TargetingComplex.EnhWrite` field | `bio` | Enhanced write rate when bound. Folds into the propensity function, not a separate channel |
| Demethylase effector | `TargetingComplex.EnhErase` field | `bio` | Enhanced erase rate when bound |
| Environmental trigger | `EnvironmentalTrigger` struct | `bio` | Deterministic scheduled event that binds a complex at its target site |
| System builder | `System.Build()` method | `bio` | Maps biology types → `[]gillespie.Reaction` + `gillespie.State` + `[]gillespie.ScheduledEvent` |

### Species Layout in `State.Counts`

```
Index:    [0]  [1]  [2]  ...  [N-1]   [N]    [N+1]  ...  [N+K-1]
Species:  m₀   m₁   m₂  ...  mₙ₋₁    b₀     b₁    ...   bₖ₋₁
          ↑ methylation states ↑       ↑ complex bound states ↑
          (0=unmethylated, 1=methylated)  (0=unbound, 1=bound)
```

Where N = number of CpG sites, K = number of targeting complexes.

### Reaction Channel Mapping

| Channel | Count | Deltas | Propensity `a_i(x)` |
|---------|------:|--------|---------------------|
| Site write (per site i) | N | `Δ[i] = +1` | `(k_bg_w + Σ{k_enh_w_j · bound_j}) · (1 - meth_i)` |
| Site erase (per site i) | N | `Δ[i] = -1` | `(k_bg_e + Σ{k_enh_e_j · bound_j}) · meth_i` |
| Complex unbind (per complex k) | K | `Δ[N+k] = -1` | `k_off · bound_k` |
| **Total** | **2N + K** | | **O(N + K) — linear scaling** |

The sum `Σ{k_enh_w_j · bound_j}` is over all complexes targeting site i. For
most sites this is empty (fast path: pure background rate). For targeted sites,
the closure checks the bound-state species of each complex aimed at that site.

### Propensity Function Design

The propensity for writing at site i demonstrates the "conditional enhancement"
pattern — the same reaction channel handles both background and enhanced rates:

```go
// Fast path (no complex targets this site):
propensity = k_bg_write * (1 - counts[siteIdx])

// Slow path (one or more complexes target this site):
rate := k_bg_write
for _, info := range complexesTargetingThisSite {
    rate += info.enhWriteRate * float64(counts[info.boundSpeciesIdx])
}
propensity = rate * (1 - counts[siteIdx])
```

This is evaluated ~2000 times per SSA step at N=1000 (once per reaction channel).
The fast path avoids the inner loop entirely for untargeted sites.

### Environmental Trigger: The Write-Ahead Log Mechanic

The EnvironmentalTrigger implements the WAL (write-ahead log) analogy:

1. **Trigger fires** at a scheduled time t (deterministic, not Poisson).
2. **Complex binds** at its target site: `State.Counts[N + complexIdx] = 1`.
3. **Enhanced propensity** activates: the bound complex gates higher write/erase
   rates at the target site via the propensity function.
4. **Methylation mark** accumulates stochastically while the complex is bound.
5. **Complex unbinds** stochastically (rate = KOff), ending the enhanced period.

The trigger's occurrence is chronologically recorded in the trajectory event stream
(FiredReaction = -(tag + 1) for scheduled events). The methylation mark that
results from the trigger is the biological "WAL entry" — it records that the
environmental event happened, readable from the methylation state long after the
trigger and the complex are gone.

**Arrival process assumption**: Environmental triggers use a **fixed external
schedule** (deterministic times), not a stochastic Poisson process. This models
a researcher or engineered system applying a stimulus at known times. If the
trigger fires while the complex is already bound, the event is a no-op (clamped
to [0,1]).

**Hybrid SSA technique**: Deterministic triggers are interleaved with stochastic
SSA steps via `RunWithSchedule()`. At each iteration, the next stochastic event
time τ is compared with the next scheduled event time. Whichever is sooner fires
first. This preserves SSA exactness for stochastic channels.

### Rate Constants (Default Values)

| Parameter | Value | Biological interpretation |
|-----------|------:|--------------------------|
| `k_bg_write` | 0.001 | Spontaneous methylation (rare without enzyme) |
| `k_bg_erase` | 0.01 | Passive demethylation (replication dilution) |
| `k_enh_write` | 0.5 | Active methyltransferase (DNMT-like) when bound |
| `k_enh_erase` | 0.05 | Active demethylase (TET-like) when bound |
| `k_off` | 0.1 | Complex dissociation (mean residence = 10 time units) |

Background steady-state: `0.001 / (0.001 + 0.01) ≈ 0.091` (9.1% methylated).
With bound complex: effective write = 0.501, effective erase = 0.06, steady-state ≈ 0.893 (89.3% methylated).

### Biological Analogy Table (Extended)

| Computer Analogy | Biology | Phase 2 Code |
|-----------------|---------|-------------|
| ROM / Disk | DNA sequence (immutable) | `Genome` struct with `[]Site` |
| RAM bit array | Methylation state vector | `State.Counts[0:N]` ∈ {0, 1}ᴺ |
| Pointer / address | dCas9 guide RNA | `TargetingComplex.TargetSite` → site index |
| Pointer register | Complex bound state | `State.Counts[N:N+K]` ∈ {0, 1}ᴷ |
| Write operation | Methyltransferase writes CH₃ | Write reaction channel, enhanced by bound complex |
| Erase operation | Demethylase removes CH₃ | Erase reaction channel, enhanced by bound complex |
| WAL entry | Trigger → bind → methylate | `EnvironmentalTrigger` → `ScheduledEvent` → propensity gating |
| Write-ahead log | Chronological trigger/event stream | `TrajectoryRecord.FiredReaction` (negative = scheduled, positive = stochastic) |
| Memory bus / data path | Reaction generation pipeline | `System.Build()` → `[]Reaction` + `State` + `[]ScheduledEvent` |

## Phase 3 Preview (do not implement yet)

The following will be needed when Phase 3 (TUI) lands:

| Feature | Expected | Notes |
|---------|----------|-------|
| Per-site methylation heatmap | TUI grid: sites × time | Color-coded 0→1 per CpG site |
| Complex binding timeline | TUI strip | Shows when/where complexes are bound |
| Trigger markers | TUI annotations | Vertical lines at trigger fire times |
| Genome context coloring | Site-level tags | Promoter vs gene-body visual distinction |
| Aggregated statistics | Rolling mean | Ensemble mean + CI for selected sites |
