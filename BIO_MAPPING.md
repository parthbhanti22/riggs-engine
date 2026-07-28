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

## Phase 2 Preview (do not implement yet)

The following will extend this table when Phase 2 lands:

| Biology | Expected Go Type | Notes |
|---------|-----------------|-------|
| Genome (DNA sequence) | `Genome` struct | Fixed sequence, many CpG sites — the "disk" |
| dCas9 pointer | `Pointer` or `Guide` | Navigates to specific CpG addresses without cutting |
| Methyltransferase | Reaction generator | Produces `Reaction` values for each targetable CpG |
| Demethylase | Reaction generator | Same, for erasure |
| Environmental trigger | TBD | Schema for what the WAL records — open question |
