# Development Roadmap

Where this project is heading, split into milestones that are each usable and
testable on their own. Update the checkboxes as work lands.

## Current state (v7 — July 2026)

- [x] Exact analytic simulation of the series **RLC / RL / RC / LC / R**
      topologies with zero initial state
- [x] **AC/DC source toggle** (Milestone 1, v4): sinusoidal E0·sin(ωt) or a
      DC step E0 at t=0, with exact closed-form step responses for every
      topology (classic RC/RL charging curves, LC step ringing, RLC step
      response); ω slider disables in DC mode, equation text and impedance
      panel adapt
- [x] **Parallel circuit family** (Milestone 2, v5): a Series/Parallel
      toggle unlocks 4 presets — **R∥C, R∥L, R∥L∥C** (current-driven, the
      dual of the series case) and **Tank** (R + an L∥C tank, voltage-
      driven). R∥L∥C shows genuine parallel (anti)resonance; Tank shows the
      notch/band-stop character of a resonant tank. The charts become
      multi-trace (voltage, per-branch currents, charge on C), the
      amplitude slider relabels itself V↔A automatically, and the
      impedance panel shows |Zp|/phase for the two resonant presets
- [x] Three synced charts: charge Q(t), current I(t), and **voltages**
      E(t), V_R, V_L, V_C (KVL-exact, color-matched to the schematic)
- [x] Transient → steady-state marker (5τ, slowest decay mode), decay
      envelope, steady-state overlay, RK4 verification overlay (series
      family; parallel family has the 5τ marker only — see Milestone 2 notes)
- [x] Impedance triangle, phase angle, resonance gauge, damping badge
      (AC-only; replaced by an explanatory note in DC mode)
- [x] Sliders + exact numeric input, animation with 1×/2×/4× speed,
      hover readout, adaptive schematic
- [x] English UI (translated from Indonesian in v3)
- [x] Headless test suite: RK4 cross-check for all series topologies and
      all parallel presets × both source types, KVL/KCL checks, V_L =
      L·dI/dt check, widget round-trips (AC/DC toggle, family toggle),
      per-topology and per-preset screenshots
- [x] **General MNA netlist engine** (Milestone 3, v6): `rlc_netlist.py` +
      `rlc_mna.py` — an arbitrary R/L/C/VSRC/ISRC netlist can now be
      simulated by trapezoidal Modified Nodal Analysis instead of only the
      9 fixed presets. Cross-validated against every one of those 9 presets
      (max relative error ~10⁻⁴, consistent with 2nd-order trapezoidal
      accuracy).
- [x] **Free-form circuit builder** (Milestone 4, v7): `rlc_builder.py` —
      a standalone app (`python rlc_builder.py`) that puts a UI on top of
      the Milestone 3 engine. Place R/L/C/VSRC/ISRC on a snapped grid,
      wire them with click-click placement (no dragging), pick a ground,
      probe any node's voltage or component's current, edit values and
      source type/frequency live, save/load as JSON. Resolves
      automatically after every edit. See its own section below.

## Design position: should we build a free-form circuit builder?

**Yes, and it's now done (Milestone 4, v7).** Building it in one jump would
have meant replacing the physics engine, the data model, and the UI at the
same time, which is where projects stall — so it happened in stages instead,
each one usable and tested on its own:

- The original engine solved **one fixed ODE family** per topology in
  closed form. Exact and fast, but it couldn't describe an arbitrary R/L/C
  network.
- A free-form builder needed a **netlist** data model ("component X connects
  node a to node b") and a **general solver** — Modified Nodal Analysis
  (MNA) with a trapezoidal time-stepper, how SPICE-class simulators work.
  `rlc_netlist.py` and `rlc_mna.py` (Milestone 3) are exactly that, and are
  validated against every fixed preset the app knows. "Parallel circuits"
  (or any other wiring) is no longer a special case requiring bespoke
  algebra — it's just a different netlist.
- Per-branch currents and per-node voltages fall out of MNA naturally,
  which is exactly the "probe every point of the circuit" feature that
  motivated this — the probes API existed before any picking UI did, and
  `rlc_builder.py` (Milestone 4) is that UI: a grid-based editor where you
  place, wire, ground, and probe an arbitrary circuit, with results
  updating live.

What's left (see Milestone 4's "deferred to backlog" note and the backlog
section below) is refinement — undo/redo, short-circuit detection, a UI
toggle for charge/energy probes — not anything architectural.

## Milestone 1 — AC/DC source toggle *(small step, high value)* — DONE (v4)

- [x] Source model abstraction: `E(t) = E0·sin(ωt)` **or** `E(t) = E0` (DC
      step at t = 0; ω slider disabled in DC mode)
- [x] Closed-form DC solutions per topology (classic charging curves:
      RC charge-up, RL current rise, LC ringing, RLC step response)
- [x] UI: AC/DC buttons under the title; equation text, analysis panel, and
      phasor/resonance panel adapt (phasors are AC-only concepts — replaced
      by a note in DC mode)
- [x] Tests: RK4 cross-check of every topology × source combination, incl.
      the pure-inductor (R=0) DC edge case

## Milestone 2 — Fixed parallel topologies — DONE (v5)

- [x] Preset list of common circuits, chosen from the same topology-button
      row (relabeled per family):
      - R∥C and R∥L — **current-driven** (see note below), 1st order, no
        resonance
      - R∥L∥C — current-driven, 2nd order, **true parallel antiresonance**
      - Tank: R in series with an L∥C tank — voltage-driven, 2nd order,
        notch/band-stop character
- [x] Solver: closed-form for all four (no numeric-integration fallback was
      needed — `solve_second_general` in `rlc_solver.py` generalizes the
      2nd-order engine to arbitrary sinusoidal phase (AC) or an arbitrary
      post-step initial derivative (DC), which both R∥L∥C and Tank reduce
      to); RK4 cross-check still backs every preset in `--test`
- [x] UI: charts become **multi-trace** — voltage (chart 1), branch
      currents I_R/I_L/I_C/I_total with a legend (chart 2), charge on C
      (chart 3, with a note when there's no capacitor); schematic gained
      `ParallelSchematic` with vertical-branch drawings and a proper
      current-source symbol
- [x] Analysis panel: parallel impedance |Zp| and phase (via complex
      admittance Yp = 1/R + j(ωC − 1/ωL), Zp = 1/Yp), resonance gauge reused
      from the series panel, shown for R∥L∥C and Tank only

**Design decision — current-source driving:** R∥C, R∥L, and R∥L∥C are
driven by a current source, not the voltage source used everywhere else in
this app. An ideal voltage source forced directly across parallel branches
decouples them completely — no branch affects another, so R∥L∥C would show
no resonance at all (verified this experimentally before committing to the
current-source design; see the physics notes in `rlc_solver.py`). A current
source is the standard dual of the series case and is what real "parallel
resonance" circuit theory uses. Tank keeps a voltage source since R is
genuinely in series with it there. This is also why the amplitude slider
relabels itself (V↔A) when switching family/preset — implemented via
`RLCApp._configure_amplitude_slider()`, which reconfigures a matplotlib
`Slider`'s valmin/valmax/label in place (confirmed to work cleanly).

**Deferred to a later pass:** steady-state overlay, RK4 verification
overlay, and the transient-envelope band are series-family-only for now —
the parallel `recompute()` path only wires up the 5τ marker (which is
generic on `alpha_settle` and was extracted into a shared
`_update_transient_marker()` helper). Re-enable these for parallel once
there's a concrete request; the DISPLAY & CONTROLS checkboxes already exist
but are inert while `family == "Parallel"`.

## Milestone 3 — General circuit engine (MNA netlist solver) — DONE (v6)

The enabler for everything after it.

- [x] Data model (`rlc_netlist.py`): `Component(kind, node_a, node_b, value,
      name, source_type, freq)` for R/L/C/VSRC/ISRC + `Netlist(components,
      ground)`; JSON round-trip via `to_json()`/`from_json()`; a cheap
      `validate()` catches floating nodes / self-shorts before they'd
      otherwise surface as a singular matrix; `series_netlist()` /
      `parallel_netlist()` build netlists matching every existing preset,
      used by the cross-validation tests
- [x] MNA assembly (`rlc_mna.py`): conductance stamps for R plus trapezoidal
      companion models for L and C (`Geq_C=2C/h, Geq_L=h/(2L)` with a
      history current source each), an extra branch-current unknown per
      voltage source (the "modified" part of MNA); a dedicated t=0 initial
      solve (`h → 0` limit: capacitors become 0V shorts, inductors become
      opens) matching this app's Q(0)=0/I(0)=0 convention everywhere else —
      this was the one real bug caught during development (see below)
- [x] Probes API: `probe_voltage(result, node)`, `probe_current(result,
      component_name)`, `probe_charge(result, capacitor_name)`,
      `probe_energy(result, component_name)` (stored energy for L/C,
      cumulative trapezoidal-integrated dissipation for R)
- [x] Cross-validated against all 9 existing closed-form presets (5 series +
      4 parallel) × AC/DC = 18 combinations, permanently in
      `rlc_simulator.py --test`; max relative error ~1.4×10⁻⁴ (the LC/R∥L∥C
      undamped cases, expected — no damping to smooth out discretization
      noise), typically ~10⁻⁶. A 4-panel overlay plot (MNA dots on the
      exact curves, indistinguishable by eye) served as the visual proof.
- [x] Performance: vectorized per-kind-of-component (not per-component)
      RHS assembly + raw LAPACK `getrs` (bypassing scipy's `lu_solve`
      wrapper, ~4-5x faster for matrices this small) got a 19-component
      circuit to **~20 recomputes/sec at 3000 samples** (this app's
      minimum resolution) — at or near the 20fps target for typical use.
      Degrades for high-resolution sweeps (~1/sec at 48000 samples, the
      app's max) since per-step cost turned out to be dominated by fixed
      Python/numpy dispatch overhead rather than component count — profiled
      and documented in `rlc_mna.py`'s module docstring. Good enough to
      proceed; revisit only once Milestone 4 actually wires this into live
      slider dragging and it's a felt problem, not a hypothetical one.

**The one real bug found:** the first implementation solved every sample
including t=0 through the same trapezoidal step, which for a source that's
already "on" at t=0 (DC step, or `sin(0)=0` for AC — no issue there)
produced a nonzero Q(0) instead of the exact 0 every closed-form solver
assumes. Standard SPICE behavior — and the fix — is to treat t=0 as a
separate operating-point solve, not one more trapezoidal step from a
fictitious t=-h.

## Milestone 4 — Free-form circuit builder UI — DONE (v7)

The engine (`rlc_netlist.Netlist` + `rlc_mna.simulate`) already existed and
was tested — this milestone was purely UI: `rlc_builder.py` produces a
`Netlist` from user interaction and hands it to `simulate()`.

- [x] Grid-based editor canvas (`BuilderApp` in `rlc_builder.py`): a palette
      of 9 tools (R, L, C, VSRC, ISRC, Wire, Ground, Select, Delete); place
      and delete components on a snapped grid. **No separate rotate step**
      — orientation is implicit from which adjacent grid point you click
      second (a deliberate simplification, see design note below). Node
      detection is automatic: components sharing a grid point automatically
      share a netlist node, and the Wire tool explicitly merges two grid
      points (via union-find) for routing around
- [x] Probe tool: with the Select tool, clicking a component toggles its
      current into the results (and lets you edit its value/source type/
      frequency); clicking a bare grid point toggles that node's voltage.
      Two charts (voltage, current) update live with a color-coded legend
      per probe. Charge probing (`probe_charge`) exists in the engine but
      has no UI toggle yet — see backlog
- [x] Component value editing: a `TextBox` for the value (not the slider
      pattern from the fixed-preset apps — free-form component values don't
      fit a bounded slider range), plus AC/DC buttons and a frequency box
      for sources, shown/hidden contextually
- [x] Validation: `Netlist.validate()` surfaces as a plain-English status
      message ("Add a ground reference…", "Not solvable yet: floating
      node(s)…") instead of a crash or a silent wrong answer. Save/load via
      `Netlist.to_json()`/`from_json()`, using a real OS file dialog
      (`tkinter.filedialog`) with a headless-safe fallback for `--test`
- [x] Tests: `rlc_builder.py --test` drives the *real* click-dispatch code
      path (`_on_click`, not internal methods directly) to build a series
      RLC circuit from scratch and checks it against the exact closed-form
      solver (~5×10⁻⁵ relative error — consistent with the coarser default
      2000-sample grid), plus probe toggling, source editing, delete,
      save/load round-trip, and validation-message coverage for every
      incomplete-circuit state

**Design decision — click-click wiring, not click-drag:** true continuous
dragging needs a fair amount of custom motion-tracking state in matplotlib
(which has no native drag-and-drop primitive) for a payoff that a two-click
model gets almost all of: click a grid point, click an *adjacent* one, done
— and a non-adjacent second click just relocates the anchor rather than
erroring, which turned out forgiving enough in testing that a drag would
add complexity without adding much usability. This also sidesteps rotation
entirely: whichever adjacent point you pick determines horizontal vs.
vertical automatically.

**Bug caught during testing:** the first version of `_place()` auto-selected
*and* auto-probed every newly-placed component (selecting for editing was
supposed to be a convenience; toggling the probe was meant only for a
deliberate click on an existing component). The result: by the time a
circuit was fully built, every single component was probed as a side effect
— not wrong data, just a confusing default. Fixed by splitting
`_select_for_edit()` (editing only, used right after placement) from
`_select()` (editing *and* probe toggle, used by Select-tool clicks).

**Deferred to backlog:** undo/redo, short-circuit warnings (`validate()`
only catches floating nodes today, not e.g. a wire directly shorting a
source), a UI toggle for charge/energy probes, and a from-scratch
performance re-check under sustained interactive use (each click-triggered
resolve is a single small circuit, not the ~20-component stress case from
Milestone 3, so no slowdown was observed in testing — but it wasn't
rigorously profiled here either).

## Backlog / nice-to-haves (any time)

- [ ] Steady-state overlay, RK4 verification overlay, and transient envelope
      for the parallel family (deferred from Milestone 2 — see its notes)
- [ ] Builder: undo/redo for placements/wires/deletes
- [ ] Builder: short-circuit warnings (`Netlist.validate()` only catches
      floating nodes today, not e.g. a wire directly shorting a source)
- [ ] Builder: a UI toggle to probe charge (`probe_charge`) and energy
      (`probe_energy`) — both already exist in `rlc_mna.py`, just not
      exposed as a click target yet, unlike voltage/current
- [ ] Builder: rigorous interactive-latency profiling under sustained real
      use (not just the single-resolve-per-click cases exercised by
      `--test`) — revisit `rlc_mna.py`'s per-step cost (Milestone 3 notes)
      if this ever becomes a felt problem
- [ ] Export: CSV of traces, one-click PNG of the figure
- [ ] Phasor diagram with rotating vectors (animated at reduced speed)
- [ ] Preset save/load for parameter sets
- [ ] Dark theme toggle

## Architecture notes for contributors

- `rlc_solver.py`, `rlc_netlist.py`, and `rlc_mna.py` are all pure numpy (no
  matplotlib) — keep it that way so engines can be tested headlessly and
  reused by a future web UI.
- Every new physics path must ship with a cross-check test (closed-form vs
  RK4 in `rlc_simulator.py --test`; closed-form vs MNA there too; and
  click-built-circuit vs closed-form in `rlc_builder.py --test`). The
  KVL/KCL identity and V_L = L·dI/dt style checks are cheap and catch sign
  errors early; for MNA specifically, remember the t=0 initial-condition
  solve is a *separate* system from the ongoing trapezoidal sweep (see
  Milestone 3 notes) — a new component kind needs its own `hstep -> 0`
  limiting behavior worked out, not just its steady-state companion model.
- UI layout is hand-tuned figure coordinates in both `rlc_app.py` and
  `rlc_builder.py`; after any layout change, re-render the `--test`
  screenshots and inspect them — cramped card spacing (title text
  overlapping content) is the most common regression, easy to miss without
  actually looking at the PNG.
- `rlc_app.py` dispatches on `self.family` ("Series"/"Parallel") at a small
  number of chokepoints — `recompute()`, `_draw_upto()`, `_place_cursor()` —
  each delegating to a `_series`/`_parallel` sibling method. Series-mode
  code is otherwise untouched by Milestone 2; new parallel-only artists
  were added as a second, initially-hidden set on the *same* 3 axes rather
  than new axes, so the two families share the transient-marker machinery.
  Keep following this pattern rather than branching deep inside a shared
  method — it made Milestone 2 safe to add without re-testing Milestone 1.
- When a physics case needs a second-order response to an arbitrary
  sinusoidal phase (not just sin) or an arbitrary post-step derivative (not
  just 0), reach for `solve_second_general()` before writing a bespoke
  solver — it already backs both R∥L∥C and Tank.
- `rlc_builder.py` tests should drive `BuilderApp._on_click` (via a fake
  event object with `.inaxes`/`.xdata`/`.ydata`) rather than calling
  internal methods like `_place`/`_select` directly — the dispatch logic
  itself (tool routing, grid snapping, adjacency checks, hit-testing) is
  exactly what's most likely to have bugs, and calling internals only
  bypasses it. When toggling something as a side effect of another action
  (e.g. `_place()` selecting the new component), think about whether that
  side effect should compose with *every* caller — it shouldn't always
  (see the auto-probe bug in Milestone 4's notes).
