# -*- coding: utf-8 -*-
"""
Interactive Series Circuit Simulator  —  R / RC / RL / LC / RLC
===============================================================
General equation (Kirchhoff):  L·Q'' + R·Q' + Q/C = E(t) = E0·sin(ω·t)
with zero initial state; the topology is chosen with the buttons in the
CIRCUIT card, and absent components are truly removed from the equation
(RC and R are solved as exact 1st-order systems, not as numeric limits).

Assignment defaults (RLC topology):
    R = 1000 Ω,  L = 3.5 H,  C = 2×10⁻⁶ F,  E(t) = 120·sin(377·t) Volt

How to run:
    python rlc_simulator.py

Controls:
    [RLC][RL][RC][LC][R] : choose the circuit topology (CIRCUIT card)
    [AC][DC]         : choose the source type (sinusoidal vs. a step at t=0)
    [Play] / space   : animate the Q(t), I(t), V(t) curves (screen-recording)
    [1×] [2×] [4×]   : animation speed (or keyboard keys 1 / 2 / 4)
    [Reset] / 'r'    : restore every parameter to the assignment values
    Sliders          : vary R, L, C, E0, ω, and the time span
    Numeric boxes    : type an exact value (e.g. 1234.5 or 3,75) then Enter —
                       not limited by the slider resolution/range
    Checkboxes       : steady-state, RK4 verification, transient envelope
    Hover the charts : read t, Q, I and the component voltages

File structure:
    rlc_config.py    — shared constants (assignment values, theme, topologies)
    rlc_solver.py    — exact analytic solutions per topology + RK4 check
    rlc_schematic.py — circuit schematic that follows the topology
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

from rlc_config import DEF, TOPO_ORDER, TOPOS   # noqa: E402
from rlc_solver import solve, solve_rk4         # noqa: E402
from rlc_app import RLCApp                      # noqa: E402


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
        print("Series Circuit Simulator (R/RC/RL/LC/RLC) - close the window "
              "to quit.")
        print("Controls: topology buttons (RLC/RL/RC/LC/R), Play (space), "
              "Reset ('r'),")
        print("          sliders + numeric boxes (Enter), speed 1x/2x/4x, "
              "checkboxes,")
        print("          hover the charts to read values.")
        plt.show()
