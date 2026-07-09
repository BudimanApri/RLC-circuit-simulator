# Development Roadmap

Where this project is heading, split into milestones that are each usable and
testable on their own. Update the checkboxes as work lands.

## Current state (v5 — July 2026)

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

## Design position: should we build a free-form circuit builder?

**Yes — as the end goal, reached in stages.** Building it in one jump would
mean replacing the physics engine, the data model, and the UI at the same
time, which is where projects stall. The staged plan below keeps the app
working at every step:

- The current engine solves **one fixed ODE family** in closed form. That is
  exact and fast, but it cannot describe an arbitrary R/L/C network.
- A free-form builder needs a **netlist** data model ("component X connects
  node a to node b") and a **general solver** — the standard approach is
  Modified Nodal Analysis (MNA) with a small implicit time-stepper, which is
  how SPICE-class simulators work. Once that engine exists, "parallel
  circuits" is not a special case anymore: *any* wiring works.
- Per-branch currents and per-node voltages fall out of MNA naturally, which
  is exactly the "probe every point of the circuit" feature you described —
  in a parallel circuit each branch has its own I, V, and (for capacitors) Q,
  so the charts must become probe-driven rather than single-trace.

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

## Milestone 3 — General circuit engine (MNA netlist solver)

The enabler for everything after it.

- [ ] Data model: `Component(kind, value, node_a, node_b)` + `Netlist`
      (components, ground node); serializable to/from JSON so circuits can
      be saved and shared
- [ ] MNA assembly (conductance matrix + sources) with companion models for
      L and C; trapezoidal integration (stable, 2nd-order accurate)
- [ ] Probes API: voltage at any node, current through any branch, charge on
      any capacitor, energy per component
- [ ] Keep the existing closed-form engines (5 series + 4 parallel presets;
      exact formulas remain the "Solution" panel feature) and use them +
      RK4 to cross-validate the MNA engine in tests
- [ ] Performance target: interactive slider dragging at ≥ 20 fps for
      circuits up to ~20 components (vectorized numpy, LU factorization
      reuse while the timestep is constant)

## Milestone 4 — Free-form circuit builder UI

- [ ] Grid-based editor canvas: place / rotate / delete components from a
      palette (R, L, C, AC/DC source, wire), click-drag wiring with
      automatic node detection
- [ ] Probe tool: click any node or component to add it to the charts; each
      probe gets a color and a legend entry; charts show V/I/Q per probe
      instead of one global trace
- [ ] Component value editing via the existing slider + numeric box pattern
      (select a component → its value binds to the controls)
- [ ] Validation & niceties: dangling-wire detection, short-circuit warnings,
      undo/redo, save/load circuit JSON
- [ ] UI technology decision: matplotlib widgets can carry an MVP (the
      current card layout already proves it), but if the editor outgrows it,
      the solver/netlist modules stay unchanged and only the view layer is
      swapped (e.g. a small web front-end talking to the same Python engine)

## Backlog / nice-to-haves (any time)

- [ ] Steady-state overlay, RK4 verification overlay, and transient envelope
      for the parallel family (deferred from Milestone 2 — see its notes)
- [ ] Energy view: energy stored in L and C, energy dissipated in R over time
- [ ] Export: CSV of traces, one-click PNG of the figure
- [ ] Phasor diagram with rotating vectors (animated at reduced speed)
- [ ] Preset save/load for parameter sets
- [ ] Dark theme toggle

## Architecture notes for contributors

- `rlc_solver.py` is pure numpy (no matplotlib) — keep it that way so engines
  can be tested headlessly and reused by a future web UI.
- Every new physics path must ship with a cross-check test in
  `rlc_simulator.py --test` (analytic vs RK4 today; MNA vs RK4 later).
  The KVL/KCL identity and V_L = L·dI/dt style checks are cheap and catch
  sign errors early.
- UI layout is hand-tuned figure coordinates in `rlc_app.py`; after any
  layout change, re-render the `--test` screenshots and inspect them.
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
