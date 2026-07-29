// Package bio provides the biological instruction set for the Riggs Engine.
//
// It maps real biological entities — genomes, CpG sites, dCas9 targeting complexes,
// methyltransferases, demethylases, and environmental triggers — onto the
// gillespie.Reaction / gillespie.State types from Phase 1.
//
// The key architectural pattern is "reaction generation": bio types know how to
// produce []gillespie.Reaction slices that the SSA engine can simulate. The bio
// package contains all biology-specific logic; the gillespie package remains a
// pure math tier with no biological semantics.
//
// The core invariant: the DNA sequence is never mutated. Only the methylation
// state vector changes. This is the ROM/RAM distinction that defines the project.
package bio

// SiteContext tags the genomic context of a CpG site.
//
// In real biology, methylation has different downstream effects depending on where
// it occurs: promoter methylation tends to silence gene transcription, while
// gene-body methylation typically doesn't. In Phase 2, this is metadata only —
// it does not affect simulation rates. Future phases may use it to compute
// transcriptional consequences of the methylation state.
type SiteContext int

const (
	// ContextPromoter indicates a CpG site in a gene's promoter region.
	// Methylation here typically silences transcription of the downstream gene.
	ContextPromoter SiteContext = iota

	// ContextGeneBody indicates a CpG site within a gene's coding/intronic region.
	// Methylation here is generally tolerated and doesn't silence transcription.
	ContextGeneBody

	// ContextIntergenic indicates a CpG site between genes.
	// Methylation here is often part of repeat-element silencing.
	ContextIntergenic
)

// String returns a human-readable name for the site context.
func (c SiteContext) String() string {
	switch c {
	case ContextPromoter:
		return "promoter"
	case ContextGeneBody:
		return "gene-body"
	case ContextIntergenic:
		return "intergenic"
	default:
		return "unknown"
	}
}

// Site represents a single CpG dinucleotide in the genome.
//
// Each Site maps to one species in the Gillespie State.Counts vector:
//
//	State.Counts[site.Index] = 0  →  unmethylated
//	State.Counts[site.Index] = 1  →  methylated
//
// The Site itself is immutable once created — only the methylation state
// (tracked in the Gillespie state vector) changes during simulation.
type Site struct {
	// Index is this site's position in the Genome's site list.
	// Also serves as the species index in State.Counts.
	Index int

	// Coordinate is the genomic position in base pairs.
	// Used for display, sorting, and (in future phases) distance-dependent effects.
	Coordinate int

	// Context is the genomic context tag (promoter, gene-body, intergenic).
	// Metadata only in Phase 2 — does not affect rates.
	Context SiteContext
}

// Genome is an ordered, immutable collection of CpG sites.
// This is the "disk" — the DNA sequence is fixed; only methylation state changes.
//
// The Genome struct itself contains no mutable state. Methylation state is tracked
// externally in the Gillespie State.Counts vector, where each site's Index
// corresponds to a species.
type Genome struct {
	// Sites is the ordered list of CpG sites in this genome.
	// Site indices are assigned sequentially: Sites[i].Index == i.
	Sites []Site
}

// NumSites returns the number of CpG sites in the genome.
func (g *Genome) NumSites() int {
	return len(g.Sites)
}

// NewGenome creates a genome with n CpG sites spaced at regular intervals.
//
// The contextPattern parameter assigns genomic context to sites in a repeating
// pattern. For example, if contextPattern is [ContextPromoter, ContextGeneBody],
// then site 0 gets ContextPromoter, site 1 gets ContextGeneBody, site 2 gets
// ContextPromoter, and so on. If contextPattern is nil or empty, all sites
// default to ContextIntergenic.
//
// Parameters:
//   - n: number of CpG sites.
//   - spacing: distance in base pairs between consecutive sites.
//   - contextPattern: repeating pattern of site contexts.
func NewGenome(n int, spacing int, contextPattern []SiteContext) *Genome {
	sites := make([]Site, n)
	for i := range sites {
		ctx := ContextIntergenic
		if len(contextPattern) > 0 {
			ctx = contextPattern[i%len(contextPattern)]
		}
		sites[i] = Site{
			Index:      i,
			Coordinate: i * spacing,
			Context:    ctx,
		}
	}
	return &Genome{Sites: sites}
}
