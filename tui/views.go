package tui

import (
	"fmt"
	"strings"

	"github.com/charmbracelet/lipgloss"

	"github.com/pxrth9/riggs/bio"
)

// renderMemoryMap renders the genome memory-map view — a horizontal bitmap
// strip where each CpG site is one character cell, making the "DNA as
// addressable memory" metaphor visually literal.
//
// Layout:
//
//	GENOME MEMORY MAP ──────────────────────────────
//	Addr: 0    5    10   15   20   25   30   35   40
//	      ████░░░░████░░░█░░░░░░░░░░░░░░░░░░░░████░
//	Ctx:  P P P G G G G G P P P G G G G G P P P G G
//	Bind: ◆──────────────────────────────────────────
func renderMemoryMap(snap *SimSnapshot, genome *bio.Genome, selected int, width int) string {
	if snap.NumSites == 0 {
		return ""
	}

	var b strings.Builder

	// Section title
	b.WriteString(titleStyle.Render("GENOME MEMORY MAP"))
	b.WriteString("\n")

	// Address ruler
	b.WriteString(dimStyle.Render("  Addr: "))
	for i := 0; i < snap.NumSites; i++ {
		if i%5 == 0 {
			label := fmt.Sprintf("%-5d", i)
			if i+5 > snap.NumSites {
				label = fmt.Sprintf("%d", i)
			}
			b.WriteString(dimStyle.Render(label))
			i += len(label) - 1
		}
	}
	b.WriteString("\n")

	// Memory cells: █ = methylated, ░ = unmethylated
	b.WriteString(dimStyle.Render("        "))
	for i := 0; i < snap.NumSites && i < len(snap.Counts); i++ {
		ch := "░"
		style := unmethylatedStyle
		if snap.Counts[i] == 1 {
			ch = "█"
			style = methylatedStyle
		}
		if i == selected {
			b.WriteString(selectedStyle.Render(ch))
		} else {
			b.WriteString(style.Render(ch))
		}
	}
	b.WriteString("\n")

	// Context row: P = promoter, G = gene body, I = intergenic
	b.WriteString(dimStyle.Render("  Ctx:  "))
	if genome != nil {
		for i := 0; i < snap.NumSites && i < len(genome.Sites); i++ {
			var ch string
			var style lipgloss.Style
			switch genome.Sites[i].Context {
			case bio.ContextPromoter:
				ch = "P"
				style = promoterStyle
			case bio.ContextGeneBody:
				ch = "G"
				style = geneBodyStyle
			default:
				ch = "·"
				style = intergenicStyle
			}
			if i == selected {
				b.WriteString(selectedStyle.Render(ch))
			} else {
				b.WriteString(style.Render(ch))
			}
		}
	}
	b.WriteString("\n")

	// Binding row: shows where complexes are bound
	b.WriteString(dimStyle.Render("  Bind: "))
	bindMap := make(map[int]bool) // site index → bound?
	if snap.NumComplexes > 0 && snap.ComplexTargets != nil {
		for k := 0; k < snap.NumComplexes; k++ {
			specIdx := snap.NumSites + k
			if specIdx < len(snap.Counts) && snap.Counts[specIdx] == 1 {
				if k < len(snap.ComplexTargets) {
					bindMap[snap.ComplexTargets[k]] = true
				}
			}
		}
	}

	for i := 0; i < snap.NumSites; i++ {
		if bindMap[i] {
			b.WriteString(boundStyle.Render("◆"))
		} else {
			b.WriteString(dimStyle.Render("─"))
		}
	}
	b.WriteString("\n")

	return b.String()
}

// renderWALTail renders the WAL (write-ahead log) tail view — a scrolling
// log of recent reaction events styled like `tail -f`.
func renderWALTail(snap *SimSnapshot, height int, width int) string {
	var b strings.Builder

	b.WriteString(titleStyle.Render("WAL TAIL"))
	b.WriteString(dimStyle.Render("  (most recent events)"))
	b.WriteString("\n")

	entries := snap.WAL.Tail(height)
	if len(entries) == 0 {
		b.WriteString(dimStyle.Render("  (no events yet)"))
		b.WriteString("\n")
		return b.String()
	}

	for _, e := range entries {
		ts := dimStyle.Render(fmt.Sprintf("  t=%-10.3f ", e.Time))
		b.WriteString(ts)

		if e.FiredReaction < 0 {
			// Scheduled event (trigger)
			tag := -(e.FiredReaction + 1)
			b.WriteString(triggerStyle.Render(fmt.Sprintf("*TRIGGER bind[%d]", tag)))
		} else {
			// Stochastic event
			site := e.SiteAffected
			if site >= 0 && site < snap.NumSites {
				if e.NewValue == 1 {
					b.WriteString(methylatedStyle.Render(fmt.Sprintf(" WRITE  site[%d]  0→1", site)))
				} else {
					b.WriteString(dimStyle.Render(fmt.Sprintf(" ERASE  site[%d]  1→0", site)))
				}
			} else if site >= snap.NumSites {
				// Complex unbind event
				cplxIdx := site - snap.NumSites
				b.WriteString(dimStyle.Render(fmt.Sprintf(" UNBIND cplx[%d]", cplxIdx)))
			} else {
				b.WriteString(dimStyle.Render(fmt.Sprintf(" EVENT  rxn[%d]", e.FiredReaction)))
			}
		}
		b.WriteString("\n")
	}

	return b.String()
}

// renderStats renders the aggregate statistics panel with per-site
// methylation fractions, bar charts, and ensemble reference data.
func renderStats(snap *SimSnapshot, genome *bio.Genome, selected int,
	ensembleMean, ensembleStd []float64, ensembleN int, width int) string {

	var b strings.Builder

	b.WriteString(titleStyle.Render("STATISTICS"))
	b.WriteString("\n")

	// Simulation clock and event count
	b.WriteString(statsLabelStyle.Render("  Sim Time:"))
	b.WriteString(statsValueStyle.Render(fmt.Sprintf("%.2f", snap.SimTime)))
	b.WriteString("    ")
	b.WriteString(statsLabelStyle.Render("Events:"))
	b.WriteString(statsValueStyle.Render(fmt.Sprintf("%d", snap.EventCount)))
	b.WriteString("\n")

	// Selected site detail
	if selected >= 0 && selected < snap.NumSites {
		b.WriteString("\n")

		ctx := ""
		if genome != nil && selected < len(genome.Sites) {
			switch genome.Sites[selected].Context {
			case bio.ContextPromoter:
				ctx = promoterStyle.Render(" [P]")
			case bio.ContextGeneBody:
				ctx = geneBodyStyle.Render(" [G]")
			default:
				ctx = intergenicStyle.Render(" [·]")
			}
		}

		b.WriteString(statsLabelStyle.Render(fmt.Sprintf("  Site %d:", selected)))
		b.WriteString(ctx)
		b.WriteString("\n")

		// Current state
		state := "unmethylated"
		if selected < len(snap.Counts) && snap.Counts[selected] == 1 {
			state = methylatedStyle.Render("METHYLATED")
		} else {
			state = unmethylatedStyle.Render("unmethylated")
		}
		b.WriteString(statsLabelStyle.Render("    State:"))
		b.WriteString(state)
		b.WriteString("\n")

		// Time-weighted fraction with bar
		if selected < len(snap.MethFracs) {
			frac := snap.MethFracs[selected]
			bar := renderBar(frac, 20)
			b.WriteString(statsLabelStyle.Render("    Meth frac:"))
			b.WriteString(statsValueStyle.Render(fmt.Sprintf("%.4f ", frac)))
			b.WriteString(bar)
			b.WriteString("\n")
		}

		// Ensemble reference (if available)
		if ensembleMean != nil && selected < len(ensembleMean) {
			b.WriteString(statsLabelStyle.Render("    Ensemble:"))
			b.WriteString(dimStyle.Render(fmt.Sprintf("%.4f ± %.4f  (N=%d)",
				ensembleMean[selected], ensembleStd[selected], ensembleN)))
			b.WriteString("\n")
		}
	}

	// Complex status
	if snap.NumComplexes > 0 {
		b.WriteString("\n")
		for k := 0; k < snap.NumComplexes; k++ {
			specIdx := snap.NumSites + k
			status := dimStyle.Render("UNBOUND")
			if specIdx < len(snap.Counts) && snap.Counts[specIdx] == 1 {
				status = boundStyle.Render("BOUND")
			}
			b.WriteString(statsLabelStyle.Render(fmt.Sprintf("  Complex %d:", k)))
			b.WriteString(status)
			b.WriteString("\n")
		}
	}

	return b.String()
}

// renderBar renders a horizontal bar chart of a fraction (0.0 to 1.0).
func renderBar(frac float64, barWidth int) string {
	if frac < 0 {
		frac = 0
	}
	if frac > 1 {
		frac = 1
	}
	filled := int(frac * float64(barWidth))
	empty := barWidth - filled

	bar := barFilledStyle.Render(strings.Repeat("█", filled)) +
		barEmptyStyle.Render(strings.Repeat("░", empty))
	return bar
}
