# Context Transfer Document — RLC Circuit Simulator

**Generated:** 2026-07-09, end of session, immediately before a context-window
reset. Read this first in the new session — it's a dense snapshot, not a
tutorial. The living, forward-looking plan is [ROADMAP.md](ROADMAP.md);
this document is a point-in-time handoff.

---

## 1. Project goals and architecture

### Origin and current scope
Started as a single-purpose simulator for one physics problem (series RLC:
R=1000 Ω, L=3.5 H, C=2 µF, E(t)=120·sin(377t) V, Q(0)=0, I(0)=0). Over four
milestones it grew into a **three-app suite sharing a physics engine**:

1. **`rlc_simulator.py`** — fixed-preset interactive simulator. Two circuit
   families (Series: R/RC/RL/LC/RLC; Parallel: R∥C/R∥L/R∥L∥C/Tank), each
   drivable by AC (E₀·sin ωt) or DC (step). 9 topologies × 2 sources = 18
   exact closed-form solutions.
2. **`rlc_builder.py`** — free-form circuit builder. Place R/L/C/VSRC/ISRC
   on a snapped grid, wire them, probe any node/branch, simulate via the
   general MNA engine. Not limited to the 9 fixed presets.
3. Both sit on top of a **general netlist engine** (`rlc_netlist.py` +
   `rlc_mna.py`) that can simulate *any* R/L/C/source network via
   trapezoidal Modified Nodal Analysis (the SPICE method).

### File map

| File | Role | Depends on |
|---|---|---|
| `rlc_config.py` | Shared constants: defaults, input bounds, theme colors, topology configs (`TOPOS`, `PARALLEL_TOPOS`) | — |
| `rlc_solver.py` | Exact closed-form solutions for all 9 fixed presets × AC/DC + RK4 verification. Pure numpy. | `rlc_config` |
| `rlc_schematic.py` | Matplotlib circuit drawing: `Schematic` (series), `ParallelSchematic` (parallel) | `rlc_config` |
| `rlc_app.py` | `RLCApp` — the interactive UI for `rlc_simulator.py`. Dispatches on `self.family` ("Series"/"Parallel") | `rlc_config`, `rlc_solver`, `rlc_schematic` |
| `rlc_simulator.py` | Entry point for `RLCApp` + comprehensive `--test` suite | `rlc_app`, `rlc_solver`, `rlc_config`, `rlc_netlist`, `rlc_mna` |
| `rlc_netlist.py` | General netlist data model: `Component`, `Netlist` (JSON round-trip). Pure Python, no numpy. | — |
| `rlc_mna.py` | General trapezoidal-MNA transient solver + probes API (`probe_voltage/current/charge/energy`). Pure numpy, optional scipy for speed. | `rlc_netlist`-shaped input (duck-typed, no import) |
| `rlc_builder.py` | `BuilderApp` — free-form circuit editor UI + its own `--test` suite | `rlc_config`, `rlc_netlist`, `rlc_mna` |
| `README.md` | User-facing docs for all 3 entry points | — |
| `ROADMAP.md` | Living dev plan: 4 milestones (all done), design-decision write-ups, backlog, architecture notes for contributors | — |

### Environment
Windows 11, PowerShell primary shell (Bash tool also available, POSIX
syntax). Python 3.14, numpy 2.3.5, matplotlib 3.10.7, scipy 1.16.3
(optional — soft dependency, used for speed in `rlc_mna.py` via raw LAPACK
`getrs`, falls back to `numpy.linalg.solve` if absent).

### Git
Repo exists, remote `https://github.com/BudimanApri/RLC-circuit-simulator.git`,
user identity `Archleon <budimanapriu999@gmail.com>`. **The user commits
their own work — do not `git commit` unless explicitly asked.** As of this
handoff, `README.md`, `ROADMAP.md`, `rlc_simulator.py` are modified and
`rlc_builder.py`, `rlc_mna.py`, `rlc_netlist.py` are new/untracked (all from
Milestones 3–4, not yet committed by the user). There's also an untracked
`Folder circuit/Circuit1.json` — **the user's own save from testing
`rlc_builder.py` interactively**, not something Claude created. Don't touch
or "clean up" that file without asking.

---

## 2. Current implementation progress

All four originally-planned roadmap milestones are **done**:

- ✅ **Milestone 1** — AC/DC source toggle. Every topology has exact DC
  step-response formulas alongside the AC ones.
- ✅ **Milestone 2** — Parallel circuit family (R∥C, R∥L, R∥L∥C current-
  driven; Tank voltage-driven). Multi-trace charts, amplitude slider that
  relabels V↔A, impedance/resonance panel for the two resonant presets.
- ✅ **Milestone 3** — General MNA netlist engine (`rlc_netlist.py` +
  `rlc_mna.py`). Cross-validated against all 9 fixed presets × AC/DC to
  ~10⁻⁴ relative error.
- ✅ **Milestone 4** — Free-form circuit builder (`rlc_builder.py`). Grid-
  based click-click placement/wiring, live simulation, probing, value/
  source editing, JSON save/load.

**Verification state:** both test suites pass cleanly as of the last run
this session:
```
python rlc_simulator.py --test out.png     # ~40+ assertions, all green
python rlc_builder.py --test out.png       # 6 test groups, all green
```
MNA performance benchmark in the last run: 20.0 recomputes/sec for a
19-component circuit at 3000 samples (right at the informal 20fps target).

**What's NOT done** (see §4 and ROADMAP.md's Backlog section) — these are
deliberate scope cuts, not oversights:
- Parallel family in `rlc_app.py` has no steady-state/RK4/envelope overlay
  (series-only for now); parallel mode only shows the 5τ transient marker.
- `rlc_builder.py` has no undo/redo, no short-circuit detection (only
  floating-node validation), and no UI toggle for charge/energy probing
  (the engine supports `probe_charge`/`probe_energy`, just not wired to a
  click target).

---

## 3. Key technical decisions

### Physics
- **Parallel R∥C/R∥L/R∥L∥C are current-source driven**, not voltage —
  verified experimentally that voltage-forcing decouples parallel branches
  entirely (kills resonance). This is the standard circuit-theory dual of
  the series (voltage-driven) case. **Tank** (R + L∥C) stays voltage-driven
  since R is genuinely in series there.
- `rlc_solver.solve_second_general()` is a reusable damped-2nd-order ODE
  engine (arbitrary AC forcing phase, or arbitrary DC post-step initial
  derivative) — backs both R∥L∥C and Tank. Reach for it before writing a
  new bespoke 2nd-order solver.
- **MNA method**: trapezoidal integration. R → conductance stamp. L/C →
  companion model (conductance + history current source), `Geq_C = 2C/h`,
  `Geq_L = h/(2L)`. VSRC → extra branch-current unknown (the "modified" in
  MNA). ISRC → RHS injection only.
- **MNA t=0 is a SEPARATE solve**, not part of the regular trapezoidal
  loop — uses the `hstep → 0` limit (capacitors become 0V shorts, inductors
  become opens), matching this app's universal Q(0)=0/I(0)=0 convention.
  This was a real bug in the first draft (see §4).
- **MNA performance**: bottleneck was Python/dict overhead per timestep,
  not the linear solve itself. Fixed via (a) vectorizing per-*kind*-of-
  component (grouped numpy arrays) instead of per-Component-object Python
  loops, (b) calling LAPACK's `getrs` directly instead of scipy's
  `lu_solve` wrapper (~4-5x faster for matrices this small).

### UI
- Amplitude slider in `rlc_app.py` dynamically relabels/rescales V↔A
  depending on family/preset (`RLCApp._configure_amplitude_slider()`),
  confirmed matplotlib `Slider.valmin/valmax/label` can be changed in
  place post-construction.
- `rlc_app.py` dispatches on `self.family` at a handful of chokepoints
  (`recompute()`, `_draw_upto()`, `_place_cursor()`) into `_series`/
  `_parallel` sibling methods, with parallel-only artists as a second,
  initially-hidden set on the *same* 3 axes (not new axes) — keeps the
  transient-marker machinery shared between families.
- **`rlc_builder.py` uses click-click placement, not drag** — matplotlib
  has no native drag-and-drop primitive; two clicks (anchor, then an
  *adjacent* grid point) gets nearly all the usability without custom
  motion-tracking state. Orientation is automatic from click direction —
  no separate rotate step exists or is needed.
- **Wires merge grid points into one netlist node** via union-find, not a
  near-zero-resistance component. Ground designation renames that node's
  group to `"0"`.
- **Source polarity convention**: whichever grid point is clicked *first*
  is `node_a` = VSRC's "+" terminal / the node ISRC current flows *from*.
  Same convention `rlc_netlist.series_netlist()`/`parallel_netlist()` use.
  Getting this backwards doesn't crash — it silently negates every result
  (I and Q both flip sign). Bit me once while writing tests; worth
  remembering when hand-building test circuits.

### Testing philosophy (both `--test` suites)
- Cross-validate new physics against whatever already-exact reference
  exists — RK4 for closed-form solvers, closed-form for MNA, MNA-via-real-
  clicks for the builder. Never trust "it looks right."
- **Builder tests must dispatch through the real `_on_click` handler**
  (via a fake event object with `.inaxes`/`.xdata`/`.ydata`), not call
  internal methods like `_place`/`_select` directly — the dispatch logic
  itself (tool routing, grid snapping, adjacency checks, hit-testing) is
  where bugs actually live.
- Always render and **visually inspect** `--test` screenshots after any
  layout change — cramped card title/content overlap is the single most
  common regression and is invisible to assertions.

---

## 4. Known constraints / bugs

### From the user, not yet specified
The user said *"a few things need to fix"* immediately before requesting
this handoff document, but did not enumerate them before the context reset.
**This is the first thing to ask about / listen for in the new session.**
Given the untracked `Folder circuit/Circuit1.json`, the user has been
actively testing `rlc_builder.py` interactively (not just via `--test`) —
the fixes likely came from that hands-on session, possibly UI/UX issues
that don't show up in scripted tests (real mouse clicks, real file dialogs,
visual layout at actual window size, etc.).

### Already-known limitations (deliberate scope cuts, documented in ROADMAP.md backlog)
- No undo/redo in the builder.
- No short-circuit detection (`Netlist.validate()` only catches floating
  nodes with <2 connections — a wire directly shorting a source's two
  terminals would just solve with a large current, no warning).
- No UI toggle for charge/energy probes in the builder (engine has
  `probe_charge`/`probe_energy`, no click target exposes them).
- MNA performance degrades for high sample counts (~1 recompute/sec at
  48000 samples) though it's fine at the app's typical 3000-sample floor
  (~20/sec). Not re-profiled under sustained *interactive* builder use.
- Parallel family in `rlc_app.py`: no steady-state overlay, RK4 overlay, or
  transient envelope (series-only). The DISPLAY & CONTROLS checkboxes for
  these exist but are inert when `family == "Parallel"`.

### Real bugs already fixed this project (context, not open issues)
- MNA t=0 nonzero-Q bug (see §3) — fixed via separate operating-point solve.
- Builder auto-probe bug: `_place()` used to call `_select()` (which both
  selects for editing *and* toggles the probe) on every newly-placed
  component, so by the time a circuit was built every component ended up
  probed as an unwanted side effect. Fixed by splitting `_select_for_edit()`
  (editing only, used by `_place()`) from `_select()` (editing + probe
  toggle, used by SELECT-tool clicks). **If touching `rlc_builder.py`
  again, double-check which of these two methods any new call site should
  use** — picking the wrong one silently changes probe state.

---

## 5. Immediate next steps

1. **Get the specifics of "a few things need to fix" from the user** — this
   is the actual next task, deferred past the context reset. Likely UI/UX
   friction found while manually testing `rlc_builder.py` (see the
   untracked save file above).
2. Fix whatever's identified, following the established verification
   pattern: reproduce, fix, re-run both `--test` suites, visually inspect
   any changed screenshots, live TkAgg smoke test if it's an interaction
   bug.
3. Beyond the immediate fixes, the ROADMAP.md backlog (§4 above) is the
   next source of work if the user wants to keep going: undo/redo,
   short-circuit warnings, charge/energy probe UI, parallel-family
   overlay parity, interactive-latency profiling.
4. Remember the uncommitted state (§1) — the user commits manually; don't
   `git add`/`git commit` proactively.

---

## How to verify the current state on resume

```bash
cd "D:\Google antigravity\RLC circuit"
python rlc_simulator.py --test out1.png    # fixed-preset engine + UI
python rlc_builder.py --test out2.png      # free-form builder
```
Both should print all-green assertion messages and save screenshots with
no exceptions. If either fails, that's a regression from whatever the user
changed by hand between sessions (they may have been editing/testing
directly) — diff against this document's understanding of "last known good."
