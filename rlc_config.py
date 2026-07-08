# -*- coding: utf-8 -*-
"""Shared configuration: assignment defaults, input bounds, theme, topologies."""

# Assignment values: C in µF, T (time span) in ms
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
