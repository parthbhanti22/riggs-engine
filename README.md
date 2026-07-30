# Riggs Engine

The Riggs Engine is a high-performance stochastic simulation engine written in Go. It is designed to model epigenetic methylation as a volatile biological memory system. By treating DNA sequences as read-only memory (ROM) and chemical methylation tags as switchable random-access memory (RAM), the engine simulates how environmental stimuli and targeted biological effectors (such as dCas9-fused enzymes) write and erase information in the genome over time.

## Overview

Biological processes are inherently noisy and stochastic. To accurately model how a single epigenetic "bit" flips on or off, the Riggs Engine utilizes the Gillespie Stochastic Simulation Algorithm (SSA). This allows the simulation of continuous-time Markov chains (CTMC) to produce exact trajectories of the system's state.

Because a single trajectory is highly random, the engine utilizes a concurrent Monte Carlo approach, running thousands of independent trajectories in parallel to calculate precise statistical probabilities of methylation at any given genomic coordinate.

## Core Features

*   **Zero-Allocation Hot Paths:** The inner simulation loop is heavily optimized, utilizing pre-allocated buffers and flat array data structures to achieve zero heap allocations per step.
*   **Deterministic Reproducibility:** The concurrent worker pool relies on independent, seed-derived pseudo-random number generators (PRNGs). This ensures zero mutex contention and guarantees that multi-threaded simulations are 100% reproducible.
*   **Hybrid SSA Implementation:** The engine seamlessly interleaves deterministic scheduled events (environmental triggers) with stochastic reaction channels, allowing researchers to model specific interventions occurring at precise times.
*   **O(N) Scaling:** The biological logic tier utilizes a conditional propensity gating architecture, ensuring that reaction channel complexity scales linearly with the size of the simulated genome.

## Architecture

The project is structured into distinct, decoupled tiers:

1.  **Math Tier (`gillespie/`):** The pure mathematical core. It implements the Gillespie SSA Direct Method, knowing nothing of biology. It operates entirely on abstract species counts, continuous reaction rates, and state vectors.
2.  **Logic Tier (`bio/`):** The biological instruction set. It maps biological entities (genomes, CpG sites, dCas9 complexes) into reaction generators that produce the mathematical channels consumed by the Math Tier.
3.  **Presentation Tier (`cmd/riggs/`):** A command-line interface for executing ensemble simulations and observing statistical outputs.

## Biological Model

*   **The Genome:** Modeled as an ordered, immutable collection of sites (bits). Base pairs are never mutated; only the overlying epigenetic state vector changes.
*   **Targeting Complexes:** Represents biological navigation tools, such as catalytically dead Cas9 (dCas9) fused with methyltransferase or demethylase effectors. These complexes can conditionally bind to specific sites, drastically altering the local reaction rates.
*   **Environmental Triggers:** Deterministic events that mimic a Write-Ahead Log (WAL). A trigger binds a complex to the genome, initiating targeted epigenetic writing, leaving a durable chemical mark long after the biological effector has unbound.

## Getting Started

### Prerequisites

*   Go 1.22 or higher.

### Installation

Clone the repository and build the binary:

```bash
git clone https://github.com/pxrth9/riggs.git
cd riggs
go build -o riggs ./cmd/riggs
```

### Usage

The engine currently supports two simulation modes via the command-line interface.

**1. Toy Mode (Phase 1 Validation)**
Simulates a single-site system to validate the mathematical core against analytical CTMC predictions.

```bash
./riggs -mode toy -trajectories 1000
```

**2. Genome Mode (Phase 2 Simulation)**
Simulates a multi-site genome (default: 50 sites) subject to targeted methylation and periodic environmental triggers.

```bash
./riggs -mode genome -sites 50 -trajectories 500 -tmax 500
```

### Command Line Flags

*   `-mode`: Simulation mode (`toy` or `genome`).
*   `-trajectories`: Number of independent trajectories to simulate in the ensemble (default: 1000).
*   `-tmax`: Maximum simulation time per trajectory (default: 1000.0).
*   `-seed`: Base seed for deterministic PRNGs (default: 42).
*   `-workers`: Number of concurrent worker threads (defaults to all available logical CPUs).
*   `-sites`: Number of CpG sites to simulate in genome mode (default: 50).
