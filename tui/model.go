package tui

import (
	"fmt"
	"time"

	"github.com/charmbracelet/bubbletea"
	"github.com/charmbracelet/lipgloss"

	"github.com/pxrth9/riggs/bio"
)

// tickInterval is the target refresh rate for the TUI (10 FPS).
const tickInterval = 100 * time.Millisecond

// tickMsg is sent by the tick command every tickInterval.
type tickMsg time.Time

// Model is the Bubble Tea model for the Riggs Engine dashboard.
type Model struct {
	runner   *SimRunner
	snapshot *SimSnapshot // local copy, updated every tick
	genome   *bio.Genome

	// UI state
	selectedSite int  // which site is highlighted by arrow keys
	width        int  // terminal width
	height       int  // terminal height
	quitting     bool // user pressed Q

	// Ensemble reference (pre-computed, static)
	ensembleMean []float64 // per-site mean methylation from ensemble
	ensembleStd  []float64 // per-site std dev
	ensembleN    int       // number of trajectories
}

// NewModel creates a new TUI model from a bio.System build result.
func NewModel(runner *SimRunner, genome *bio.Genome) Model {
	return Model{
		runner:       runner,
		snapshot:     &SimSnapshot{},
		genome:       genome,
		selectedSite: 0,
	}
}

// SetEnsembleStats sets pre-computed ensemble statistics for the reference display.
func (m *Model) SetEnsembleStats(mean, std []float64, n int) {
	m.ensembleMean = mean
	m.ensembleStd = std
	m.ensembleN = n
}

// Init starts the simulation goroutine and the tick timer.
func (m Model) Init() tea.Cmd {
	m.runner.Start()
	return tea.Tick(tickInterval, func(t time.Time) tea.Msg {
		return tickMsg(t)
	})
}

// Update handles incoming messages (ticks, key presses, window resize).
func (m Model) Update(msg tea.Msg) (tea.Model, tea.Cmd) {
	switch msg := msg.(type) {

	case tickMsg:
		// Read the latest snapshot from the simulation goroutine
		m.runner.ReadSnapshot(m.snapshot)

		// Schedule next tick
		return m, tea.Tick(tickInterval, func(t time.Time) tea.Msg {
			return tickMsg(t)
		})

	case tea.KeyMsg:
		switch msg.String() {
		case "q", "ctrl+c":
			m.quitting = true
			m.runner.Stop()
			return m, tea.Quit

		case " ":
			// Toggle pause/resume
			m.runner.TogglePause()
			// Immediately update snapshot to reflect pause state
			m.runner.ReadSnapshot(m.snapshot)
			return m, nil

		case "n":
			// Step one event (only effective when paused)
			m.runner.StepOnce()
			return m, nil

		case "left", "h":
			if m.selectedSite > 0 {
				m.selectedSite--
			}
			return m, nil

		case "right", "l":
			if m.selectedSite < m.snapshot.NumSites-1 {
				m.selectedSite++
			}
			return m, nil
		}

	case tea.WindowSizeMsg:
		m.width = msg.Width
		m.height = msg.Height
		return m, nil
	}

	return m, nil
}

// View renders the complete dashboard.
func (m Model) View() string {
	if m.quitting {
		return ""
	}

	if m.width == 0 {
		return "Initializing..."
	}

	var s string

	// Header
	s += renderHeader(m.width, m.snapshot)
	s += "\n"

	// Memory Map (genome bitmap)
	s += renderMemoryMap(m.snapshot, m.genome, m.selectedSite, m.width)
	s += "\n"

	// WAL Tail (scrolling event log)
	walHeight := m.height - 20 // reserve space for header, map, stats, help
	if walHeight < 5 {
		walHeight = 5
	}
	if walHeight > 20 {
		walHeight = 20
	}
	s += renderWALTail(m.snapshot, walHeight, m.width)
	s += "\n"

	// Stats Panel
	s += renderStats(m.snapshot, m.genome, m.selectedSite, m.ensembleMean, m.ensembleStd, m.ensembleN, m.width)
	s += "\n"

	// Help bar
	s += renderHelp(m.snapshot, m.width)

	return s
}

// --- Styles ---

var (
	headerStyle = lipgloss.NewStyle().
			Bold(true).
			Foreground(lipgloss.Color("#7C3AED")).
			BorderStyle(lipgloss.NormalBorder()).
			BorderBottom(true).
			BorderForeground(lipgloss.Color("#4C1D95"))

	titleStyle = lipgloss.NewStyle().
			Bold(true).
			Foreground(lipgloss.Color("#A78BFA"))

	dimStyle = lipgloss.NewStyle().
			Foreground(lipgloss.Color("#6B7280"))

	methylatedStyle = lipgloss.NewStyle().
			Foreground(lipgloss.Color("#10B981")).
			Bold(true)

	unmethylatedStyle = lipgloss.NewStyle().
				Foreground(lipgloss.Color("#374151"))

	selectedStyle = lipgloss.NewStyle().
			Background(lipgloss.Color("#4C1D95")).
			Foreground(lipgloss.Color("#F9FAFB")).
			Bold(true)

	promoterStyle = lipgloss.NewStyle().
			Foreground(lipgloss.Color("#F59E0B"))

	geneBodyStyle = lipgloss.NewStyle().
			Foreground(lipgloss.Color("#6366F1"))

	intergenicStyle = lipgloss.NewStyle().
			Foreground(lipgloss.Color("#6B7280"))

	boundStyle = lipgloss.NewStyle().
			Foreground(lipgloss.Color("#EF4444")).
			Bold(true)

	triggerStyle = lipgloss.NewStyle().
			Foreground(lipgloss.Color("#F59E0B")).
			Bold(true)

	statsLabelStyle = lipgloss.NewStyle().
			Foreground(lipgloss.Color("#9CA3AF")).
			Width(18)

	statsValueStyle = lipgloss.NewStyle().
			Foreground(lipgloss.Color("#E5E7EB")).
			Bold(true)

	barFilledStyle = lipgloss.NewStyle().
			Foreground(lipgloss.Color("#10B981"))

	barEmptyStyle = lipgloss.NewStyle().
			Foreground(lipgloss.Color("#374151"))

	helpStyle = lipgloss.NewStyle().
			Foreground(lipgloss.Color("#6B7280"))

	pausedStyle = lipgloss.NewStyle().
			Foreground(lipgloss.Color("#EF4444")).
			Bold(true)

	doneStyle = lipgloss.NewStyle().
			Foreground(lipgloss.Color("#10B981")).
			Bold(true)

	sectionStyle = lipgloss.NewStyle().
			BorderStyle(lipgloss.NormalBorder()).
			BorderBottom(true).
			BorderForeground(lipgloss.Color("#374151"))
)

// renderHeader renders the top status bar.
func renderHeader(width int, snap *SimSnapshot) string {
	title := titleStyle.Render("RIGGS ENGINE")

	status := ""
	if snap.Done {
		status = doneStyle.Render(" DONE")
	} else if snap.Paused {
		status = pausedStyle.Render(" PAUSED")
	} else {
		status = dimStyle.Render(" RUNNING")
	}

	clock := fmt.Sprintf("t=%.2f", snap.SimTime)
	events := fmt.Sprintf("events=%d", snap.EventCount)
	right := dimStyle.Render(clock + "  " + events)

	// Lay out: title + status on left, clock on right
	left := title + status
	gap := width - lipgloss.Width(left) - lipgloss.Width(right)
	if gap < 2 {
		gap = 2
	}
	padding := ""
	for i := 0; i < gap; i++ {
		padding += " "
	}

	return headerStyle.Width(width).Render(left + padding + right)
}

// renderHelp renders the bottom help bar.
func renderHelp(snap *SimSnapshot, width int) string {
	help := "Space: pause/resume  N: step  ←/→: select site  Q: quit"
	return helpStyle.Width(width).Render(help)
}
