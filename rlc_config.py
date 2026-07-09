# -*- coding: utf-8 -*-
"""Shared configuration: defaults, input bounds, theme, topologies."""

# Default values: C in µF, T (time span) in ms
DEF = dict(R=1000.0, L=3.5, C=2.0, E0=120.0, W=377.0, T=80.0)

# Physical bounds for typed input (wider than the slider ranges)
BOUNDS = dict(R=(0.0, 1e6), L=(1e-3, 1e3), C=(1e-4, 1e4),
              E0=(0.0, 1e5), W=(0.1, 1e5), T=(1.0, 5000.0))

# (key, label, lo, hi, step) — visual slider ranges
SLIDERS = [("R", "R  (Ω)", 0.0, 4000.0, 10.0),
           ("L", "L  (H)", 0.5, 10.0, 0.05),
           ("C", "C  (µF)", 0.2, 10.0, 0.05),
           ("E0", "E₀  (V)", 0.0, 300.0, 5.0),
           ("W", "ω  (rad/s)", 50.0, 800.0, 1.0),
           ("T", "Time span (ms)", 10.0, 300.0, 5.0)]

# ---- color theme -------------------------------------------------------------
TH = dict(bg="#eef2f7", card="#ffffff", edge="#d9e0ec",
          text="#0f172a", sub="#5b6779", accent="#2563eb",
          q="#2563eb", i="#dc2626", ss="#64748b", rk="#0f172a",
          env="#f59e0b", trans="#b45309", steady="#047857")

# Component / source colors (shared by schematic and voltage traces)
C_SRC, C_R, C_L, C_C = "#0369a1", "#b45309", "#047857", "#6d28d9"

# Damping badge colors: substring of damping label → (text color, bg color)
BADGE = {"underdamped": ("#1d4ed8", "#dbeafe"),
         "critically": ("#15803d", "#dcfce7"),
         "overdamped": ("#6d28d9", "#ede9fe"),
         "1st-order": ("#0f766e", "#ccfbf1"),
         "no transient": ("#475569", "#e2e8f0"),
         "undamped": ("#475569", "#e2e8f0")}

# ---- source types --------------------------------------------------------------
SOURCE_ORDER = ["AC", "DC"]

# ---- circuit topologies -------------------------------------------------------
TOPO_ORDER = ["RLC", "RL", "RC", "LC", "R"]
TOPOS = {
    "RLC": dict(R=True, L=True, C=True,
                eq="L·Q″ + R·Q′ + Q/C = E₀·sin(ωt)      "
                   "Q(0) = 0,   I(0) = 0",
                eq_dc="L·Q″ + R·Q′ + Q/C = E₀ (step at t=0)      "
                      "Q(0) = 0,   I(0) = 0"),
    "RL": dict(R=True, L=True, C=False,
               eq="L·Q″ + R·Q′ = E₀·sin(ωt)      Q(0) = 0,   I(0) = 0      "
                  "(no capacitor: Q = charge delivered)",
               eq_dc="L·Q″ + R·Q′ = E₀ (step at t=0)      Q(0) = 0,   "
                     "I(0) = 0      (no capacitor: Q = charge delivered)"),
    "RC": dict(R=True, L=False, C=True,
               eq="R·Q′ + Q/C = E₀·sin(ωt)      Q(0) = 0      "
                  "(1st-order, no inductor)",
               eq_dc="R·Q′ + Q/C = E₀ (step at t=0)      Q(0) = 0      "
                     "(1st-order, no inductor)"),
    "LC": dict(R=False, L=True, C=True,
               eq="L·Q″ + Q/C = E₀·sin(ωt)      Q(0) = 0,   I(0) = 0      "
                  "(no resistance: undamped oscillation)",
               eq_dc="L·Q″ + Q/C = E₀ (step at t=0)      Q(0) = 0,   "
                     "I(0) = 0      (no resistance: undamped oscillation)"),
    "R": dict(R=True, L=False, C=False,
              eq="R·I = E₀·sin(ωt)   →   I(t) = E(t)/R      "
                 "(algebraic, no transient)",
              eq_dc="R·I = E₀ (step at t=0)   →   I(t) = E₀/R      "
                    "(algebraic, no transient)"),
}

# ---- circuit family: series (one shared current) vs. parallel (branch currents) --
FAMILY_ORDER = ["Series", "Parallel"]

# R∥C, R∥L, R∥L∥C are driven by a CURRENT source I(t) — the dual of the
# series case — because an ideal *voltage* source forced directly across
# parallel branches would decouple them entirely (no interaction, no
# resonance). Tank (R in series with an L∥C tank) stays voltage-driven
# since R is genuinely in series with the source there.
PARALLEL_ORDER = ["RC_P", "RL_P", "RLC_P", "TANK"]
PARALLEL_TOPOS = {
    "RC_P": dict(L=False, C=True, src="I", label="R∥C",
                 eq="C·V′ + V/R = I₀·sin(ωt)      V(0) = 0      "
                    "(current-driven, 1st-order)",
                 eq_dc="C·V′ + V/R = I₀ (step at t=0)      V(0) = 0      "
                       "(current-driven, 1st-order)"),
    "RL_P": dict(L=True, C=False, src="I", label="R∥L",
                 eq="(L/R)·I_L′ + I_L = I₀·sin(ωt)      I_L(0) = 0      "
                    "(current-driven, 1st-order in I_L)",
                 eq_dc="(L/R)·I_L′ + I_L = I₀ (step at t=0)      "
                       "I_L(0) = 0      (current-driven, 1st-order)"),
    "RLC_P": dict(L=True, C=True, src="I", label="R∥L∥C",
                  eq="C·V″ + V′/R + V/L = I₀·ω·cos(ωt)      "
                     "V(0) = 0,   V′(0) = 0      (current-driven)",
                  eq_dc="C·V″ + V′/R + V/L = 0  for t>0      V(0) = 0,   "
                        "V′(0) = I₀/C      (current-driven step)"),
    "TANK": dict(L=True, C=True, src="E", label="Tank",
                 eq="C·V_t″ + V_t′/R + V_t/L = (E₀ω/R)·cos(ωt)      "
                    "V_t(0) = 0,   V_t′(0) = 0      (R + L∥C tank)",
                 eq_dc="C·V_t″ + V_t′/R + V_t/L = 0  for t>0      "
                       "V_t(0) = 0,   V_t′(0) = E₀/(RC)      "
                       "(R + L∥C tank, step)"),
}

# Default/slider spec for the amplitude control in Parallel family: current-
# driven presets need a current (A) in a much smaller, more realistic range
# than the series voltage slider; the Tank preset keeps the voltage slider
# unchanged since it is voltage-driven like the series topologies.
PARALLEL_I0_DEFAULT = 0.1
PARALLEL_I0_BOUNDS = (0.0, 1e4)
PARALLEL_I0_SLIDER = (0.0, 0.5, 0.005)
