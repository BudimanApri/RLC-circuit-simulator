# Series Circuit Simulator — R / RC / RL / LC / RLC

Besides the full RLC circuit (the default, matching the assignment), the
simulator can also run the **RL**, **RC**, **LC**, and **R-only** variants —
absent components are genuinely removed from the equation and the schematic,
not just set to zero. Each topology can also be driven by an **AC** source
(E₀·sin(ωt), the assignment default) or a **DC** step (E₀ applied at t = 0),
producing the classic charging/discharging curves from circuit theory.

## Running

```
python rlc_simulator.py
```

Requires `numpy` and `matplotlib` (`pip install numpy matplotlib`).

Headless test mode (numeric verification + screenshots of every topology):

```
python rlc_simulator.py --test output.png
```

## Controls

| Control | Function |
|---|---|
| **RLC / RL / RC / LC / R** buttons (CIRCUIT card) | Choose the topology; sliders of absent components are disabled automatically |
| **AC / DC** buttons (under the equation) | Choose the source: sinusoidal E₀·sin(ωt), or a DC step E₀ at t = 0. The ω slider disables in DC mode, and the impedance/resonance panel (an AC-only concept) is replaced by a note |
| **Play** button / space bar | Animate the Q(t), I(t), and V(t) curves progressively |
| **1× / 2× / 4×** buttons (or keys `1`/`2`/`4`) | Playback speed |
| **Reset** button / key `r` | Restore every parameter to the assignment values |
| Sliders R, L, C, E₀, ω | Explore the effect of each component |
| **Numeric box** next to each slider | Type an exact value and press Enter (e.g. `1234.5` or `3,75`) — not limited by slider resolution or range |
| "Time span" slider | Change the displayed time window |
| **Steady-state** checkbox | Show the steady-state solution (dashed) |
| **RK4 (numeric)** checkbox | Verify the analytic solution with numerical integration |
| **Transient envelope** checkbox | Orange envelope: steady-state amplitude + transient decay bound |
| Hover the charts | Read t, Q, I and all component voltages at the cursor (while paused) |

## Charts

1. **Charge Q(t)** — capacitor charge (or *charge delivered* for topologies
   without a capacitor, where Q keeps a permanent offset).
2. **Current I(t)** — with initial-condition markers and RK4/steady-state
   overlays.
3. **Voltages** — source E(t) plus the per-component voltages V_R = R·I,
   V_L = L·dI/dt, and V_C = Q/C, color-matched to the schematic. Kirchhoff's
   voltage law V_R + V_L + V_C = E(t) holds exactly at every instant. Near
   resonance the plot shows the classic voltage magnification: V_L and V_C
   each exceed the 120 V source amplitude (~158 V) while nearly cancelling
   each other.

## Circuit variants

| Topology | Equation | Character |
|---|---|---|
| **RLC** | L·Q″ + R·Q′ + Q/C = E(t) | 2nd order: underdamped / critical / overdamped |
| **RL** | L·Q″ + R·Q′ = E(t) | 1st-order in current, τ = L/R; Q = charge delivered (permanent offset) |
| **RC** | R·Q′ + Q/C = E(t) | 1st order, τ = RC, pure exponential decay |
| **LC** | L·Q″ + Q/C = E(t) | Undamped: the transient never decays (beats near resonance) |
| **R** | R·I = E(t) | Algebraic: current is immediately steady-state, no transient |

Every topology is solved with its own exact closed-form solution (RC and R
as true 1st-order systems, not numeric limits of the 2nd-order case), and all
of them are verified against RK4 integration in `--test` mode, including a
KVL check and an independent V_L = L·dI/dt check.

## Source types

| Source | E(t) | What it shows |
|---|---|---|
| **AC** (default) | E₀·sin(ωt) | Default with default values |
| **DC** | E₀ for t ≥ 0 (step) | Classic charging/discharging curves: RC capacitor charge-up, RL current rise, LC step ringing (never settles), RLC step response. The steady state is a constant (or a ramp for topologies without a capacitor) instead of a sinusoid |

The damping classification (underdamped/critical/overdamped/1st-order,
α, ω_d, the characteristic roots) is identical between AC and DC — it only
depends on R, L, C. Only the steady-state / particular solution and the
"Solution" formulas differ. Switching source type keeps the same topology,
component values, and time span.

## Display elements

- **Transient → steady-state marker**: the transient region is shaded orange
  and bounded by a dashed line at t = 5τ ("transient practically over"), with
  TRANSIENT / STEADY-STATE labels above the chart. τ uses the slowest decay
  mode (for overdamped: 1/|r₁|). If 5τ exceeds the time window the whole
  chart is shaded; with no transient at all (R-only) the whole chart is
  labeled STEADY-STATE.
- **Damping badge** in the ANALYSIS panel: underdamped / critically damped /
  overdamped / 1st-order / undamped / no transient, colored by case.
- **IMPEDANCE & RESONANCE panel**: the impedance triangle (R, X = X_L−X_C, Z)
  drawn to true scale, |Z| and φ values, the circuit character
  (inductive/capacitive/resistive), and a gauge of ω against ω₀ (only when
  both L and C are present; the marker turns green near resonance).
- **Adaptive resolution**: the number of curve points follows ω, ω₀, and the
  time constants so curves stay smooth even for extreme typed values.

## File structure

| File | Contents |
|---|---|
| `rlc_simulator.py` | Entry point + `--test` mode |
| `rlc_app.py` | Matplotlib UI (charts, cards, widgets, animation) |
| `rlc_solver.py` | Exact analytic solutions per topology + RK4 check (pure numpy) |
| `rlc_schematic.py` | Circuit schematic that follows the topology |
| `rlc_config.py` | Shared constants: assignment values, input bounds, theme, topologies |

See [ROADMAP.md](ROADMAP.md) for the development plan (parallel circuits,
free-form circuit builder, AC/DC sources, per-branch probes).
