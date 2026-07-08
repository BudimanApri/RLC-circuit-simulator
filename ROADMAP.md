# Development Roadmap

Where this project is heading, split into milestones that are each usable and
testable on their own. Update the checkboxes as work lands.

## Current state (v4 — July 2026)

- [x] Exact analytic simulation of the series **RLC / RL / RC / LC / R**
      topologies with zero initial state
- [x] **AC/DC source toggle** (Milestone 1, v4): sinusoidal E0·sin(ωt) or a
      DC step E0 at t=0, with exact closed-form step responses for every
      topology (classic RC/RL charging curves, LC step ringing, RLC step
      response); ω slider disables in DC mode, equation text and impedance
      panel adapt
- [x] Three synced charts: charge Q(t), current I(t), and **voltages**
      E(t), V_R, V_L, V_C (KVL-exact, color-matched to the schematic)
- [x] Transient → steady-state marker (5τ, slowest decay mode), decay
      envelope, steady-state overlay, RK4 verification overlay
- [x] Impedance triangle, phase angle, resonance gauge, damping badge
      (AC-only; replaced by an explanatory note in DC mode)
- [x] Sliders + exact numeric input, animation with 1×/2×/4× speed,
      hover readout, adaptive schematic
- [x] English UI (translated from Indonesian in v3)
- [x] Headless test suite: RK4 cross-check for all topologies × both source
      types, KVL check, V_L = L·dI/dt check, widget round-trips (incl. the
      AC/DC toggle), per-topology and per-source screenshots

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

## Milestone 2 — Fixed parallel topologies

- [ ] Preset list of common two-branch circuits, chosen from a second row of
      topology buttons, e.g.:
      - R ∥ C and R ∥ L (driven directly by the source)
      - R in series with an L ∥ C tank (the classic resonant tank)
      - R ∥ L ∥ C (parallel resonance / antiresonance)
- [ ] Solver: these are still small linear ODE systems — keep closed-form
      where practical, otherwise integrate numerically (this is a good
      dress rehearsal for Milestone 3)
- [ ] UI: charts become **multi-trace** — one current per branch with a
      legend (I_total, I_R, I_L, I_C), per-branch voltage equal across
      parallel branches; schematic gains parallel branch drawings
- [ ] Analysis panel: parallel impedance Z_p = (1/Z₁ + 1/Z₂)⁻¹,
      antiresonance indicator

## Milestone 3 — General circuit engine (MNA netlist solver)

The enabler for everything after it.

- [ ] Data model: `Component(kind, value, node_a, node_b)` + `Netlist`
      (components, ground node); serializable to/from JSON so circuits can
      be saved and shared
- [ ] MNA assembly (conductance matrix + sources) with companion models for
      L and C; trapezoidal integration (stable, 2nd-order accurate)
- [ ] Probes API: voltage at any node, current through any branch, charge on
      any capacitor, energy per component
- [ ] Keep the existing closed-form engine for the five series presets
      (exact formulas remain the "assignment answer" feature) and use it +
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
  The KVL identity and V_L = L·dI/dt style checks are cheap and catch sign
  errors early.
- UI layout is hand-tuned figure coordinates in `rlc_app.py`; after any
  layout change, re-render the `--test` screenshots and inspect them.
