# -*- coding: utf-8 -*-
"""
Interactive Circuit Simulator  —  Series (R/RC/RL/LC/RLC) and
Parallel (R∥C, R∥L, R∥L∥C, Tank) families
=======================================================================
Series family (Kirchhoff, voltage-driven):
    L·Q'' + R·Q' + Q/C = E(t) = E0·sin(ω·t) or a DC step, zero initial
    state; absent components are truly removed from the equation (RC and
    R are exact 1st-order systems, not numeric limits of the 2nd-order one).

Parallel family (Milestone 2): R∥C, R∥L, and R∥L∥C are driven by a
CURRENT source (the dual of the series case — an ideal voltage source
forced directly across parallel branches would decouple them entirely).
Tank (R in series with an L∥C tank) stays voltage-driven since R is
genuinely in series with the source there.

Assignment defaults (series RLC topology):
    R = 1000 Ω,  L = 3.5 H,  C = 2×10⁻⁶ F,  E(t) = 120·sin(377·t) Volt

How to run:
    python rlc_simulator.py

Controls:
    [Series][Parallel]   : choose the circuit family
    [RLC][RL][RC][LC][R] or [R∥C][R∥L][R∥L∥C][Tank] : choose the topology
    [AC][DC]         : choose the source type (sinusoidal vs. a step at t=0)
    [Play] / space   : animate the curves (screen-recording)
    [1×] [2×] [4×]   : animation speed (or keyboard keys 1 / 2 / 4)
    [Reset] / 'r'    : restore every parameter to the assignment values
    Sliders          : vary R, L, C, the amplitude, ω, and the time span
                       (the amplitude slider is volts for Series/Tank and
                       amps for the current-driven parallel presets)
    Numeric boxes    : type an exact value (e.g. 1234.5 or 3,75) then Enter —
                       not limited by the slider resolution/range
    Checkboxes       : steady-state, RK4 verification, transient envelope
                       (series family only, for now)
    Hover the charts : read the numeric values at the cursor

File structure:
    rlc_config.py    — shared constants (defaults, theme, topologies)
    rlc_solver.py    — exact analytic solutions (series + parallel) + RK4
    rlc_schematic.py — circuit schematics (series and parallel)
    rlc_app.py       — matplotlib UI (charts, cards, widgets)
    rlc_simulator.py — entry point + test mode (--test [output.png])
"""

import os
import sys

import matplotlib

if "--test" in sys.argv:
    matplotlib.use("Agg")

import numpy as np                              # noqa: E402
import matplotlib.pyplot as plt                 # noqa: E402

from rlc_config import (DEF, TOPO_ORDER, TOPOS, PARALLEL_ORDER,   # noqa: E402
                        PARALLEL_TOPOS, PARALLEL_I0_DEFAULT)
from rlc_solver import (solve, solve_rk4, solve_parallel,         # noqa: E402
                        solve_parallel_rk4)
from rlc_app import RLCApp                                         # noqa: E402


def _self_test(app):
    R, L, C, E0, w = DEF["R"], DEF["L"], DEF["C"], DEF["E0"], DEF["W"]
    t = np.linspace(0, 0.08, 3000)

    # 1) analytic vs RK4 for ALL topologies, both AC and DC sources
    for mode in ("AC", "DC"):
        for topo in TOPO_ORDER:
            sol = solve(topo, R, L, C, E0, w, t, mode=mode)
            Qr, Ir = solve_rk4(topo, R, L, C, E0, w, t, mode=mode)
            rq = float(np.max(np.abs(Qr - sol["Q"]))) / \
                max(float(np.max(np.abs(sol["Q"]))), 1e-15)
            ri = float(np.max(np.abs(Ir - sol["I"]))) / \
                max(float(np.max(np.abs(sol["I"]))), 1e-15)
            print(f"{mode} {topo:3s} | {sol['damping']:38s} | "
                  f"rel|dQ| = {rq:.1e}   rel|dI| = {ri:.1e}")
            assert rq < 1e-6 and ri < 1e-6, (mode, topo, rq, ri)

    # RL with R=0 (pure inductor) is an edge case worth its own DC check
    sol = solve("RL", 0.0, L, C, E0, w, t, mode="DC")
    Qr, Ir = solve_rk4("RL", 0.0, L, C, E0, w, t, mode="DC")
    rq = float(np.max(np.abs(Qr - sol["Q"]))) / max(float(np.max(np.abs(sol["Q"]))), 1e-15)
    ri = float(np.max(np.abs(Ir - sol["I"]))) / max(float(np.max(np.abs(sol["I"]))), 1e-15)
    assert rq < 1e-6 and ri < 1e-6, ("RL R=0 DC", rq, ri)
    print("DC edge : pure inductor (R=0) step response OK")

    # 2) voltages: KVL holds and V_L matches L·dI/dt independently
    for mode in ("AC", "DC"):
        for topo in TOPO_ORDER:
            sol = solve(topo, R, L, C, E0, w, t, mode=mode)
            kvl = float(np.max(np.abs(sol["VR"] + sol["VL"] + sol["VC"]
                                      - sol["E"])))
            assert kvl < 1e-9 * E0, (mode, topo, kvl)
            if TOPOS[topo]["L"]:
                dIdt = np.gradient(sol["I"], t)
                err = float(np.max(np.abs(sol["VL"] - L * dIdt)[2:-2])) / \
                    max(float(np.max(np.abs(sol["VL"]))), 1e-12)
                assert err < 2e-3, (mode, topo, err)
    print("voltage : KVL + V_L = L·dI/dt OK  (AC and DC)")

    # 2b) parallel circuits (Milestone 2): analytic vs RK4 + exact KCL
    t2 = np.linspace(0, 0.08, 4000)
    for mode in ("AC", "DC"):
        for preset in PARALLEL_ORDER:
            amp = E0 if preset == "TANK" else PARALLEL_I0_DEFAULT
            psol = solve_parallel(preset, R, L, C, amp, w, t2, mode=mode)
            rk = solve_parallel_rk4(preset, R, L, C, amp, w, t2, mode=mode)
            # RL_P's RK4 state variable is I_L directly; every other preset
            # integrates V (or V_tank).
            ref = psol["I_L"] if preset == "RL_P" else psol["V"]
            rv = float(np.max(np.abs(rk[0] - ref))) / \
                max(float(np.max(np.abs(ref))), 1e-15)
            print(f"{mode} {preset:5s} | {psol['damping']:22s} | "
                  f"rel|dV| = {rv:.1e}")
            assert rv < 1e-6, (mode, preset, rv)

            branches = np.zeros_like(t2)
            for key in ("I_R", "I_L", "I_C"):
                if psol[key] is not None:
                    branches = branches + psol[key]
            if preset == "TANK":
                kcl = float(np.max(np.abs(psol["I_L"] + psol["I_C"]
                                          - psol["I_total"])))
            else:
                kcl = float(np.max(np.abs(branches - psol["I_total"])))
            assert kcl < 1e-9 * max(amp, 1.0), (mode, preset, kcl)
    print("parallel: RK4 cross-check + exact KCL OK  (AC and DC, "
          "all 4 presets)")

    # R=0 edge case for R∥L (current-driven): R=0 is a dead short across L,
    # so I_L stays ~0 for any practical time window (all current takes the
    # R=0 path) — I_L itself is tiny here (~1e-12), so check against I0's
    # scale rather than I_L's own near-zero scale, which would blow up a
    # relative-error check on noise.
    psol = solve_parallel("RL_P", 0.0, L, C, PARALLEL_I0_DEFAULT, w, t2,
                          mode="DC")
    rk = solve_parallel_rk4("RL_P", 0.0, L, C, PARALLEL_I0_DEFAULT, w, t2,
                            mode="DC")
    rv = float(np.max(np.abs(rk[0] - psol["I_L"]))) / PARALLEL_I0_DEFAULT
    assert rv < 1e-9, rv
    print("parallel: R∥L with R=0 edge case OK")

    # 3) assignment answer numbers (default RLC)
    sol = solve("RLC", R, L, C, E0, w, t)
    print("damping :", sol["damping"])
    print("alpha   :", sol["alpha"], " wd:", sol.get("wd"))
    print("A, B    :", sol["A"], sol["B"])
    print("Z, Iamp :", sol["Z"], sol["Iamp"])
    print(sol["fq"])
    print(sol["fi"])

    # 4) typed-input path: exact value, decimal comma, invalid input,
    #    value beyond the slider range
    app.boxes["R"].set_val("1234.5")
    assert abs(app.vals["R"] - 1234.5) < 1e-9, app.vals["R"]
    app.boxes["L"].set_val("3,75")
    assert abs(app.vals["L"] - 3.75) < 1e-9, app.vals["L"]
    app.boxes["C"].set_val("not-a-number")
    assert abs(app.vals["C"] - 2.0) < 1e-12, app.vals["C"]
    app.boxes["W"].set_val("2500")
    assert abs(app.vals["W"] - 2500.0) < 1e-9, app.vals["W"]
    assert abs(app.sliders["W"].val - 800.0) < 1e-9  # thumb pinned at the end
    app.reset()
    assert all(abs(app.vals[kk] - DEF[kk]) < 1e-12 for kk in DEF)
    print("widget  : typed input + reset OK")

    # 5) animation speed: one 4× tick advances four times as far
    app._set_speed(4)
    app.idx = 0
    app._tick()
    assert app.idx == min(app.step * 4, app.npts), app.idx
    app._set_speed(1)
    app.idx = 0
    app._tick()
    assert app.idx == min(app.step, app.npts), app.idx
    print("widget  : speed 1x/2x/4x OK")

    # 6) topology selection: absent-component sliders disabled + screenshots
    out = sys.argv[sys.argv.index("--test") + 1] if \
        len(sys.argv) > sys.argv.index("--test") + 1 else "rlc_test.png"
    base, ext = os.path.splitext(out)
    app._set_topo("RC")
    assert app.sliders["L"].active is False
    assert app.sliders["R"].active is True
    assert app.ln_vl.get_visible() is False
    assert "1st-order" in app.sol["damping"]
    app.fig.savefig(base + "_rc" + ext, dpi=110)
    app._set_topo("LC")
    assert app.sliders["R"].active is False
    assert "undamped" in app.sol["damping"]
    app.fig.savefig(base + "_lc" + ext, dpi=110)
    app._set_topo("RL")
    assert app.sliders["C"].active is False
    assert "1st-order" in app.sol["damping"]
    app.fig.savefig(base + "_rl" + ext, dpi=110)
    app._set_topo("R")
    assert "no transient" in app.sol["damping"]
    app.fig.savefig(base + "_r" + ext, dpi=110)
    app._set_topo("RLC")
    assert all(app.sliders[kk].active for kk in ("R", "L", "C"))
    assert all(ln.get_visible() for ln in
               (app.ln_vr, app.ln_vl, app.ln_vc, app.ln_e))
    print("widget  : topology selection OK")

    # 7) AC/DC source toggle: ω disabled in DC, impedance panel hidden,
    #    equation text switches, screenshots of a few DC step responses
    assert app.source_mode == "AC"
    assert app.sliders["W"].active is True
    assert app.ax_ph.get_visible() is True
    app._set_source("DC")
    assert app.sliders["W"].active is False
    assert app.ax_ph.get_visible() is False
    assert app.txt_dc_note.get_visible() is True
    assert "step" in app.txt_eq.get_text()
    assert "underdamped" in app.sol["damping"]      # RLC DC still 2nd order
    app.fig.savefig(base + "_dc_rlc" + ext, dpi=110)
    app._set_topo("RC")
    assert "1st-order" in app.sol["damping"]
    app.fig.savefig(base + "_dc_rc" + ext, dpi=110)
    app._set_topo("LC")
    assert "undamped" in app.sol["damping"]
    app.fig.savefig(base + "_dc_lc" + ext, dpi=110)
    app._set_source("AC")
    assert app.sliders["W"].active is True
    assert app.ax_ph.get_visible() is True
    assert "sin" in app.txt_eq.get_text()
    app._set_topo("RLC")
    print("widget  : AC/DC source toggle OK")

    # 7b) parallel family (Milestone 2): family toggle, preset switching,
    #     amplitude-slider unit swap, impedance panel behaviour, screenshots
    assert app.family == "Series"
    assert app.ln_q.get_visible() is True
    assert app.ln_pv.get_visible() is False
    app._set_family("Parallel")
    assert app.ln_q.get_visible() is False
    assert app.ln_pv.get_visible() is True
    assert app.parallel_preset == "RC_P"
    assert app._amp_is_voltage is False
    assert abs(app.vals["E0"] - PARALLEL_I0_DEFAULT) < 1e-9
    assert app.ln_pil.get_visible() is False       # RC_P has no L branch
    assert app.txt_no_c.get_visible() is False      # RC_P has a capacitor
    app.fig.savefig(base + "_par_rc" + ext, dpi=110)

    app._set_parallel_preset("RL_P")
    assert app.ln_pic.get_visible() is False        # RL_P has no C branch
    assert app.txt_no_c.get_visible() is True        # nothing to plot on ax_v
    app.fig.savefig(base + "_par_rl" + ext, dpi=110)

    app._set_parallel_preset("RLC_P")
    assert app.ln_pir.get_visible() and app.ln_pil.get_visible() \
        and app.ln_pic.get_visible()
    assert app.psol["Zp"] is not None                # AC: has antiresonance
    app.fig.savefig(base + "_par_rlc" + ext, dpi=110)

    app._set_parallel_preset("TANK")
    assert app._amp_is_voltage is True               # Tank is voltage-driven
    assert abs(app.vals["E0"] - DEF["E0"]) < 1e-9
    assert app.ln_pe.get_visible() is True            # Tank shows E(t) too
    app.fig.savefig(base + "_par_tank" + ext, dpi=110)

    app._set_source("DC")
    assert app.ax_res.get_visible() is False          # DC hides the gauge
    app.fig.savefig(base + "_par_tank_dc" + ext, dpi=110)
    app._set_source("AC")

    app._set_family("Series")
    assert app.ln_q.get_visible() is True
    assert app.ln_pv.get_visible() is False
    assert app._amp_is_voltage is True
    assert abs(app.vals["E0"] - DEF["E0"]) < 1e-9
    print("widget  : parallel family (Milestone 2) OK")

    # 8) animation + overlay path, then the main screenshot
    app.chk.set_active(0)
    app.chk.set_active(1)
    app.idx = app.npts // 2
    app._draw_upto(app.idx)
    app.fig.savefig(out, dpi=110)
    print("saved   :", out)


if __name__ == "__main__":
    app = RLCApp()
    if "--test" in sys.argv:
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
        _self_test(app)
    else:
        print("Circuit Simulator (Series R/RC/RL/LC/RLC + Parallel "
              "R∥C/R∥L/R∥L∥C/Tank) - close the window to quit.")
        print("Controls: family buttons (Series/Parallel), topology "
              "buttons, AC/DC, Play (space),")
        print("          Reset ('r'), sliders + numeric boxes (Enter), "
              "speed 1x/2x/4x, checkboxes,")
        print("          hover the charts to read values.")
        plt.show()
