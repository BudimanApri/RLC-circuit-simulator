# -*- coding: utf-8 -*-
"""
Exact analytic solutions + RK4 verification for the series circuit family
R / RC / RL / LC / RLC with source E(t) = E0·sin(ωt) and zero initial state.

Every solver returns a dict `sol` with keys:
    Q, I            : full solution (arrays, units C and A)
    Qss, Iss        : effective long-term behaviour — includes the permanent
                      charge offset on topologies without a capacitor
    env_q_lo/hi     : lower/upper envelope of Q
    env_i_lo/hi     : lower/upper envelope of I
    E, VR, VL, VC   : source and per-component voltages (arrays, Volt)
    VRm, VLm, VCm   : steady-state voltage amplitudes per component
    alpha           : displayed damping coefficient (R/2L for 2nd order, 1/RC…)
    alpha_settle    : decay rate of the slowest transient mode;
                      0 = transient never decays, inf = no transient at all
    damping, fq, fi : damping label + answer formulas for Q(t) and I(t)
    XL, XC, Z, phi, w0, Qamp, Iamp, A, B, and wd / r1 / r2 when relevant
"""

import math

import numpy as np

from rlc_config import TOPOS


def _g(x):
    return f"{x:.4g}"


def _pm(x):
    """Format ' + 1.23e-04' / ' - 1.23e-04' for chaining formula terms."""
    return (" + " if x >= 0 else " - ") + f"{abs(x):.4g}"


def _pack(sol, Qh, Qss, Ih, Iss, center, envh_q, envh_i, alpha_settle, fq, fi):
    """Assemble the uniform output for all solvers."""
    spread_q = sol["Qamp"] + envh_q
    spread_i = sol["Iamp"] + envh_i
    sol.update(Q=Qh + Qss, I=Ih + Iss,
               Qss=Qss + center, Iss=Iss,
               env_q_hi=center + spread_q, env_q_lo=center - spread_q,
               env_i_hi=spread_i, env_i_lo=-spread_i,
               alpha_settle=alpha_settle, fq=fq, fi=fi)
    return sol


def solve_second(R, L, k, E0, w, t):
    """2nd order:  L·Q'' + R·Q' + k·Q = E0·sin(wt),  Q(0)=I(0)=0.

    k = 1/C; k = 0 means no capacitor (RL / pure L) — the root s = 0
    yields a constant charge offset that never decays.
    """
    D = max((k - L * w * w) ** 2 + (R * w) ** 2, 1e-30)
    A = E0 * (k - L * w * w) / D
    B = -E0 * R * w / D
    Qss = A * np.sin(w * t) + B * np.cos(w * t)
    Iss = A * w * np.cos(w * t) - B * w * np.sin(w * t)

    # Homogeneous part with Qh(0) = -Qp(0), Qh'(0) = -Qp'(0) so Q(0)=I(0)=0
    q0, i0 = -B, -A * w
    alpha = R / (2.0 * L)
    disc = R * R - 4.0 * L * k
    scale = R * R + 4.0 * L * k
    center = 0.0

    sol = dict(A=A, B=B, alpha=alpha, XL=w * L, XC=k / w,
               w0=math.sqrt(k / L), Qamp=math.hypot(A, B))
    sol["Z"] = max(math.hypot(R, sol["XL"] - sol["XC"]), 1e-12)
    sol["phi"] = math.degrees(math.atan2(sol["XL"] - sol["XC"], R))
    sol["Iamp"] = E0 / sol["Z"]

    if disc < -1e-12 * scale:                   # underdamped
        wd = math.sqrt(-disc) / (2.0 * L)
        c1 = q0
        c2 = (i0 + alpha * c1) / wd
        e = np.exp(-alpha * t)
        Qh = e * (c1 * np.cos(wd * t) + c2 * np.sin(wd * t))
        Ih = e * ((c2 * wd - c1 * alpha) * np.cos(wd * t)
                  - (c1 * wd + c2 * alpha) * np.sin(wd * t))
        d1, d2 = c2 * wd - c1 * alpha, -(c1 * wd + c2 * alpha)
        envh_q = math.hypot(c1, c2) * e
        envh_i = math.hypot(d1, d2) * e
        alpha_settle = alpha
        sol.update(damping="underdamped", wd=wd)
        if alpha == 0:
            fq = (f"Q(t) = {_g(c1)}·cos({_g(wd)} t){_pm(c2)}·sin({_g(wd)} t)"
                  f"{_pm(A)}·sin({_g(w)} t){_pm(B)}·cos({_g(w)} t)   [Coulomb]")
            fi = (f"I(t) = {_g(d1)}·cos({_g(wd)} t){_pm(d2)}·sin({_g(wd)} t)"
                  f"{_pm(A * w)}·cos({_g(w)} t){_pm(-B * w)}·sin({_g(w)} t)"
                  f"   [Ampere]")
        else:
            fq = (f"Q(t) = e^(-{_g(alpha)} t)·[{_g(c1)}·cos({_g(wd)} t)"
                  f"{_pm(c2)}·sin({_g(wd)} t)]"
                  f"{_pm(A)}·sin({_g(w)} t){_pm(B)}·cos({_g(w)} t)   [Coulomb]")
            fi = (f"I(t) = e^(-{_g(alpha)} t)·[{_g(d1)}·cos({_g(wd)} t)"
                  f"{_pm(d2)}·sin({_g(wd)} t)]"
                  f"{_pm(A * w)}·cos({_g(w)} t){_pm(-B * w)}·sin({_g(w)} t)"
                  f"   [Ampere]")
    elif disc > 1e-12 * scale:                  # overdamped, incl. RL (k = 0)
        sq = math.sqrt(disc) / (2.0 * L)
        r1, r2 = -alpha + sq, -alpha - sq
        if k == 0.0:
            r1, r2 = 0.0, -R / L                # exact, avoids rounding error
        c1 = (i0 - r2 * q0) / (r1 - r2)
        c2 = q0 - c1
        Qh = c1 * np.exp(r1 * t) + c2 * np.exp(r2 * t)
        Ih = c1 * r1 * np.exp(r1 * t) + c2 * r2 * np.exp(r2 * t)
        sol.update(r1=r1, r2=r2)
        if k == 0.0:
            # e^(0·t) term = permanent charge offset; current decays ~ e^(-R/L t)
            center = c1
            envh_q = abs(c2) * np.exp(r2 * t)
            envh_i = abs(c2 * r2) * np.exp(r2 * t)
            alpha_settle = -r2
            sol["damping"] = "1st-order in current (τ = L/R)"
            fq = (f"Q(t) = {_g(c1)}{_pm(c2)}·e^({_g(r2)} t)"
                  f"{_pm(A)}·sin({_g(w)} t){_pm(B)}·cos({_g(w)} t)   [Coulomb]")
            fi = (f"I(t) = {_g(c2 * r2)}·e^({_g(r2)} t)"
                  f"{_pm(A * w)}·cos({_g(w)} t){_pm(-B * w)}·sin({_g(w)} t)"
                  f"   [Ampere]")
        else:
            envh_q = abs(c1) * np.exp(r1 * t) + abs(c2) * np.exp(r2 * t)
            envh_i = (abs(c1 * r1) * np.exp(r1 * t)
                      + abs(c2 * r2) * np.exp(r2 * t))
            alpha_settle = -r1                  # slowest root
            sol["damping"] = "overdamped"
            fq = (f"Q(t) = {_g(c1)}·e^({_g(r1)} t){_pm(c2)}·e^({_g(r2)} t)"
                  f"{_pm(A)}·sin({_g(w)} t){_pm(B)}·cos({_g(w)} t)   [Coulomb]")
            fi = (f"I(t) = {_g(c1 * r1)}·e^({_g(r1)} t)"
                  f"{_pm(c2 * r2)}·e^({_g(r2)} t)"
                  f"{_pm(A * w)}·cos({_g(w)} t){_pm(-B * w)}·sin({_g(w)} t)"
                  f"   [Ampere]")
    else:                                       # critically damped
        r = -alpha
        c1 = q0
        c2 = i0 - r * c1
        e = np.exp(r * t)
        Qh = (c1 + c2 * t) * e
        Ih = (c2 + r * (c1 + c2 * t)) * e
        envh_q = (abs(c1) + abs(c2) * t) * e
        envh_i = (abs(c2 + r * c1) + abs(r * c2) * t) * e
        alpha_settle = alpha
        sol.update(damping="critically damped", r1=r)
        fq = (f"Q(t) = ({_g(c1)}{_pm(c2)}·t)·e^({_g(r)} t)"
              f"{_pm(A)}·sin({_g(w)} t){_pm(B)}·cos({_g(w)} t)   [Coulomb]")
        fi = (f"I(t) = ({_g(c2 + r * c1)}{_pm(r * c2)}·t)·e^({_g(r)} t)"
              f"{_pm(A * w)}·cos({_g(w)} t){_pm(-B * w)}·sin({_g(w)} t)"
              f"   [Ampere]")

    if alpha == 0 and "1st-order" not in sol["damping"]:
        sol["damping"] = "undamped (no damping)"
    return _pack(sol, Qh, Qss, Ih, Iss, center, envh_q, envh_i,
                 alpha_settle, fq, fi)


def solve_first(R, k, E0, w, t):
    """1st order (no inductor):  R·Q' + k·Q = E0·sin(wt),  Q(0)=0.

    k = 1/C; k = 0 means a pure resistor — I(t) = E(t)/R with no transient.
    I(0) = 0 is satisfied automatically because E(0) = 0.
    """
    R = max(R, 1e-9)
    D = max(k * k + (R * w) ** 2, 1e-30)
    A = E0 * k / D
    B = -E0 * R * w / D
    Qss = A * np.sin(w * t) + B * np.cos(w * t)
    Iss = A * w * np.cos(w * t) - B * w * np.sin(w * t)

    K = -B                                      # so that Q(0) = 0
    a1 = k / R
    e = np.exp(-a1 * t)
    Qh = K * e
    Ih = -a1 * K * e

    sol = dict(A=A, B=B, alpha=a1, XL=0.0, XC=k / w, w0=0.0,
               Qamp=math.hypot(A, B))
    sol["Z"] = max(math.hypot(R, sol["XC"]), 1e-12)
    sol["phi"] = math.degrees(math.atan2(-sol["XC"], R))
    sol["Iamp"] = E0 / sol["Z"]

    if k > 0:                                   # RC
        center = 0.0
        envh_q = abs(K) * e
        envh_i = a1 * abs(K) * e
        alpha_settle = a1
        sol["damping"] = "1st-order: exponential decay (τ = RC)"
        fq = (f"Q(t) = {_g(K)}·e^(-{_g(a1)} t)"
              f"{_pm(A)}·sin({_g(w)} t){_pm(B)}·cos({_g(w)} t)   [Coulomb]")
        fi = (f"I(t) = {_g(-a1 * K)}·e^(-{_g(a1)} t)"
              f"{_pm(A * w)}·cos({_g(w)} t){_pm(-B * w)}·sin({_g(w)} t)"
              f"   [Ampere]")
    else:                                       # pure R: Qh = K is constant
        center = K
        envh_q = np.zeros_like(t)
        envh_i = np.zeros_like(t)
        alpha_settle = math.inf
        sol["damping"] = "no transient (purely resistive)"
        fq = f"Q(t) = {_g(K)}{_pm(B)}·cos({_g(w)} t)   [Coulomb]"
        fi = f"I(t) = {_g(E0 / R)}·sin({_g(w)} t)   [Ampere]"
    return _pack(sol, Qh, Qss, Ih, Iss, center, envh_q, envh_i,
                 alpha_settle, fq, fi)


def solve_second_dc(R, L, k, E0, t):
    """2nd order DC step:  L·Q'' + R·Q' + k·Q = E0 for t>=0,  Q(0)=I(0)=0.

    k = 1/C; k = 0 means no capacitor (RL / pure inductor). Unlike the AC
    case the steady state is a constant (or, without a capacitor, a ramp) —
    the damping classification (roots, α, ω_d) is identical to the AC
    solver since it only depends on R, L, C, not on the source shape.
    """
    alpha = R / (2.0 * L)
    sol = dict(alpha=alpha, w0=(math.sqrt(k / L) if k > 0 else 0.0))

    if k > 0:
        disc = R * R - 4.0 * L * k
        scale = R * R + 4.0 * L * k
        Qp = E0 / k                              # final capacitor charge
        q0, i0 = -Qp, 0.0                         # so that Q(0)=0, I(0)=0

        if disc < -1e-12 * scale:                # underdamped
            wd = math.sqrt(-disc) / (2.0 * L)
            c1 = q0
            c2 = (i0 + alpha * c1) / wd
            e = np.exp(-alpha * t)
            Qh = e * (c1 * np.cos(wd * t) + c2 * np.sin(wd * t))
            d1, d2 = c2 * wd - c1 * alpha, -(c1 * wd + c2 * alpha)
            Ih = e * (d1 * np.cos(wd * t) + d2 * np.sin(wd * t))
            envh_q = math.hypot(c1, c2) * e
            envh_i = math.hypot(d1, d2) * e
            alpha_settle = alpha
            sol.update(damping="underdamped", wd=wd)
            if alpha == 0:
                fq = (f"Q(t) = {_g(Qp)}{_pm(c1)}·cos({_g(wd)} t)"
                      f"{_pm(c2)}·sin({_g(wd)} t)   [Coulomb]")
                fi = (f"I(t) = {_g(d1)}·cos({_g(wd)} t)"
                      f"{_pm(d2)}·sin({_g(wd)} t)   [Ampere]")
            else:
                fq = (f"Q(t) = {_g(Qp)} + e^(-{_g(alpha)} t)·"
                      f"[{_g(c1)}·cos({_g(wd)} t){_pm(c2)}·sin({_g(wd)} t)]"
                      f"   [Coulomb]")
                fi = (f"I(t) = e^(-{_g(alpha)} t)·[{_g(d1)}·cos({_g(wd)} t)"
                      f"{_pm(d2)}·sin({_g(wd)} t)]   [Ampere]")
        elif disc > 1e-12 * scale:               # overdamped
            sq = math.sqrt(disc) / (2.0 * L)
            r1, r2 = -alpha + sq, -alpha - sq
            c1 = (i0 - r2 * q0) / (r1 - r2)
            c2 = q0 - c1
            Qh = c1 * np.exp(r1 * t) + c2 * np.exp(r2 * t)
            Ih = c1 * r1 * np.exp(r1 * t) + c2 * r2 * np.exp(r2 * t)
            envh_q = abs(c1) * np.exp(r1 * t) + abs(c2) * np.exp(r2 * t)
            envh_i = (abs(c1 * r1) * np.exp(r1 * t)
                      + abs(c2 * r2) * np.exp(r2 * t))
            alpha_settle = -r1
            sol.update(damping="overdamped", r1=r1, r2=r2)
            fq = (f"Q(t) = {_g(Qp)}{_pm(c1)}·e^({_g(r1)} t)"
                  f"{_pm(c2)}·e^({_g(r2)} t)   [Coulomb]")
            fi = (f"I(t) = {_g(c1 * r1)}·e^({_g(r1)} t)"
                  f"{_pm(c2 * r2)}·e^({_g(r2)} t)   [Ampere]")
        else:                                     # critically damped
            r = -alpha
            c1 = q0
            c2 = i0 - r * c1
            e = np.exp(r * t)
            Qh = (c1 + c2 * t) * e
            Ih = (c2 + r * (c1 + c2 * t)) * e
            envh_q = (abs(c1) + abs(c2) * t) * e
            envh_i = (abs(c2 + r * c1) + abs(r * c2) * t) * e
            alpha_settle = alpha
            sol.update(damping="critically damped", r1=r)
            fq = (f"Q(t) = {_g(Qp)} + ({_g(c1)}{_pm(c2)}·t)·e^({_g(r)} t)"
                  f"   [Coulomb]")
            fi = (f"I(t) = ({_g(c2 + r * c1)}{_pm(r * c2)}·t)·e^({_g(r)} t)"
                  f"   [Ampere]")

        Q = Qh + Qp
        I = Ih
        Qss = np.full_like(t, Qp)
        Iss = np.zeros_like(t)
    else:
        # k == 0: no capacitor — RL (R>0) or a pure inductor (R=0)
        if R > 1e-9:
            a = E0 / R                            # final current
            c2 = E0 * L / (R * R)
            c1 = -c2
            e = np.exp(-(R / L) * t)
            Q = a * t + c1 + c2 * e
            I = a * (1.0 - e)
            Qss = a * t
            Iss = np.full_like(t, a)
            envh_q = abs(c2) * e
            envh_i = a * e
            alpha_settle = R / L
            sol.update(damping="1st-order in current (τ = L/R)")
            fq = (f"Q(t) = {_g(a)}·t {_pm(c1)} + {_g(c2)}·"
                  f"e^(-{_g(R / L)} t)   [Coulomb]")
            fi = f"I(t) = {_g(a)}·(1 − e^(-{_g(R / L)} t))   [Ampere]"
        else:                                     # pure inductor: ramps forever
            I = (E0 / L) * t
            Q = (E0 / (2.0 * L)) * t * t
            Qss, Iss = Q, I
            envh_q = np.zeros_like(t)
            envh_i = np.zeros_like(t)
            alpha_settle = 0.0
            sol.update(damping="undamped (no damping)")
            fq = f"Q(t) = {_g(E0 / (2 * L))}·t²   [Coulomb]"
            fi = f"I(t) = {_g(E0 / L)}·t   [Ampere]"

    if alpha == 0 and "1st-order" not in sol["damping"]:
        sol["damping"] = "undamped (no damping)"
    sol.update(Q=Q, I=I, Qss=Qss, Iss=Iss,
               env_q_hi=Qss + envh_q, env_q_lo=Qss - envh_q,
               env_i_hi=Iss + envh_i, env_i_lo=Iss - envh_i,
               alpha_settle=alpha_settle, fq=fq, fi=fi)
    return sol


def solve_first_dc(R, k, E0, t):
    """1st order DC step (no inductor):  R·Q' + k·Q = E0 for t>=0,  Q(0)=0.

    k = 1/C; k = 0 means a pure resistor — I(t) = E0/R with no transient.
    """
    R = max(R, 1e-9)
    sol = dict(alpha=(k / R if k > 0 else 0.0), w0=0.0)

    if k > 0:                                    # RC charge-up
        Qp = E0 / k
        a1 = k / R
        e = np.exp(-a1 * t)
        Q = Qp * (1.0 - e)
        I = (E0 / R) * e
        Qss = np.full_like(t, Qp)
        Iss = np.zeros_like(t)
        envh_q = Qp * e
        envh_i = (E0 / R) * e
        alpha_settle = a1
        sol["damping"] = "1st-order: exponential decay (τ = RC)"
        fq = f"Q(t) = {_g(Qp)}·(1 − e^(-{_g(a1)} t))   [Coulomb]"
        fi = f"I(t) = {_g(E0 / R)}·e^(-{_g(a1)} t)   [Ampere]"
    else:                                         # pure resistor
        I = np.full_like(t, E0 / R)
        Q = (E0 / R) * t
        Qss, Iss = Q, I
        envh_q = np.zeros_like(t)
        envh_i = np.zeros_like(t)
        alpha_settle = math.inf
        sol["damping"] = "no transient (purely resistive)"
        fq = f"Q(t) = {_g(E0 / R)}·t   [Coulomb]"
        fi = f"I(t) = {_g(E0 / R)}   [Ampere]"

    sol.update(Q=Q, I=I, Qss=Qss, Iss=Iss,
               env_q_hi=Qss + envh_q, env_q_lo=Qss - envh_q,
               env_i_hi=Iss + envh_i, env_i_lo=Iss - envh_i,
               alpha_settle=alpha_settle, fq=fq, fi=fi)
    return sol


def solve(topo, R, L, Cuf, E0, w, t, mode="AC"):
    """Pick the solver for the topology; absent components are removed.

    `mode` selects the source: "AC" for E0·sin(ωt) (ω used), or "DC" for a
    constant step E0 applied at t=0 (ω ignored). Also derives the source
    and per-component voltages:
        VR = R·I,   VC = Q/C,   VL = E − VR − VC  (Kirchhoff, exact)
    """
    cfg = TOPOS[topo]
    Reff = R if cfg["R"] else 0.0
    k = 1.0 / (Cuf * 1e-6) if cfg["C"] else 0.0

    if mode == "DC":
        if cfg["L"]:
            sol = solve_second_dc(Reff, L, k, E0, t)
        else:
            sol = solve_first_dc(Reff, k, E0, t)
        E = np.full_like(t, E0)
    else:
        if cfg["L"]:
            sol = solve_second(Reff, L, k, E0, w, t)
        else:
            sol = solve_first(Reff, k, E0, w, t)
        E = E0 * np.sin(w * t)

    Rv = max(Reff, 1e-9) if cfg["R"] else 0.0
    sol["E"] = E
    sol["VR"] = Rv * sol["I"]
    sol["VC"] = k * sol["Q"]
    sol["VL"] = E - sol["VR"] - sol["VC"]
    if mode == "AC":
        sol["VRm"] = Rv * sol["Iamp"]
        sol["VLm"] = sol["XL"] * sol["Iamp"]
        sol["VCm"] = sol["XC"] * sol["Iamp"]
    return sol


def solve_rk4(topo, R, L, Cuf, E0, w, t, mode="AC"):
    """RK4 numerical integration as an independent check of the analytics."""
    cfg = TOPOS[topo]
    Reff = R if cfg["R"] else 0.0
    k = 1.0 / (Cuf * 1e-6) if cfg["C"] else 0.0
    dt = t[1] - t[0]
    n_sub = max(2, int(math.ceil(dt / 8e-6)))
    h = dt / n_sub
    Q = np.empty_like(t)
    I = np.empty_like(t)
    Q[0] = I[0] = 0.0

    def esrc(tt):
        return E0 * math.sin(w * tt) if mode == "AC" else E0

    if cfg["L"]:
        q = i = 0.0

        def f(tt, q, i):
            return i, (esrc(tt) - Reff * i - k * q) / L

        for n in range(1, len(t)):
            tt = t[n - 1]
            for _ in range(n_sub):
                k1q, k1i = f(tt, q, i)
                k2q, k2i = f(tt + h / 2, q + h / 2 * k1q, i + h / 2 * k1i)
                k3q, k3i = f(tt + h / 2, q + h / 2 * k2q, i + h / 2 * k2i)
                k4q, k4i = f(tt + h, q + h * k3q, i + h * k3i)
                q += h / 6 * (k1q + 2 * k2q + 2 * k3q + k4q)
                i += h / 6 * (k1i + 2 * k2i + 2 * k3i + k4i)
                tt += h
            Q[n], I[n] = q, i
    else:
        Rf = max(Reff, 1e-9)

        def g(tt, q):
            return (esrc(tt) - k * q) / Rf

        q = 0.0
        for n in range(1, len(t)):
            tt = t[n - 1]
            for _ in range(n_sub):
                k1 = g(tt, q)
                k2 = g(tt + h / 2, q + h / 2 * k1)
                k3 = g(tt + h / 2, q + h / 2 * k2)
                k4 = g(tt + h, q + h * k3)
                q += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
                tt += h
            Q[n] = q
        Earr = E0 * np.sin(w * t) if mode == "AC" else np.full_like(t, E0)
        I = (Earr - k * Q) / Rf
    return Q, I


# =============================================================================
# Parallel circuits (Milestone 2)
#
# R∥C, R∥L, R∥L∥C are driven by a CURRENT source I(t) — the dual of the
# series case. An ideal *voltage* source forced directly across parallel
# branches would decouple them completely (no branch affects another), which
# is both non-standard and pedagogically flat, so a current source is used
# instead — this is also what makes the true "antiresonance" of R∥L∥C show
# up. Tank (R in series with an L∥C tank) stays voltage-driven since R is
# genuinely in series with the source there.
#
# Every parallel solver returns a dict with keys:
#     V, Vp           : node/tank voltage and its time-derivative (arrays)
#     Src, src_kind    : the prescribed source waveform and "I" or "E"
#     I_R, I_L, I_C    : branch currents (None where the branch is absent)
#     I_total          : total current supplied by the source
#     Q_C              : capacitor charge (None if no capacitor)
#     alpha, alpha_settle, damping, wd / r1, r2 : same vocabulary as solve()
#     Zp, phi_p (AC only) : magnitude/phase of the impedance the source sees
#     w0 (2nd-order presets only) : natural frequency 1/√(LC)
#     fV, fI           : formula strings for the "Solution" panel
# =============================================================================

def solve_second_general(Leff, Reff, keff, t, mode="AC", w=0.0, F0=0.0,
                          phase=0.0, y0=0.0, yp0=0.0, yname="y"):
    """General damped 2nd-order response:  Leff·y'' + Reff·y' + keff·y = F(t).

    AC (mode="AC"): F(t) = F0·sin(wt + phase), with y(0)=0, y'(0)=0 forced
    (the usual "switch closes at t=0" convention).

    DC (mode="DC"): F(t) = 0 for t>0 (a pure homogeneous response), with the
    given initial state y(0)=y0, y'(0)=yp0 — used to represent the
    instantaneous KCL/KVL jump immediately after a step at t=0, which for
    these parallel/tank circuits is generally nonzero even though y itself
    starts at 0 (e.g. a capacitor voltage can't jump, but its slope can).

    The damping classification (roots, α, ω_d) only depends on Leff, Reff,
    keff — identical in structure to solve_second / solve_second_dc.
    """
    alpha = Reff / (2.0 * Leff)
    disc = Reff * Reff - 4.0 * Leff * keff
    scale = Reff * Reff + 4.0 * Leff * keff if keff > 0 else max(Reff * Reff, 1.0)

    if mode == "AC":
        D = max((keff - Leff * w * w) ** 2 + (Reff * w) ** 2, 1e-30)
        cp, sp = math.cos(phase), math.sin(phase)
        A = F0 * ((keff - Leff * w * w) * cp + Reff * w * sp) / D
        B = F0 * ((keff - Leff * w * w) * sp - Reff * w * cp) / D
        yss = A * np.sin(w * t) + B * np.cos(w * t)
        ypss = A * w * np.cos(w * t) - B * w * np.sin(w * t)
        q0, i0 = -B, -A * w
    else:
        A = B = 0.0
        yss = np.zeros_like(t)
        ypss = np.zeros_like(t)
        q0, i0 = y0, yp0

    sol = dict(alpha=alpha, w0=math.sqrt(keff / Leff) if keff > 0 else 0.0)

    if disc < -1e-12 * scale:                        # underdamped
        wd = math.sqrt(-disc) / (2.0 * Leff)
        c1 = q0
        c2 = (i0 + alpha * c1) / wd
        e = np.exp(-alpha * t)
        yh = e * (c1 * np.cos(wd * t) + c2 * np.sin(wd * t))
        d1, d2 = c2 * wd - c1 * alpha, -(c1 * wd + c2 * alpha)
        yph = e * (d1 * np.cos(wd * t) + d2 * np.sin(wd * t))
        envh_y = math.hypot(c1, c2) * e
        envh_yp = math.hypot(d1, d2) * e
        alpha_settle = alpha
        sol.update(damping="underdamped", wd=wd)
        fy = (f"{yname}(t) = e^(-{_g(alpha)} t)·[{_g(c1)}·cos({_g(wd)} t)"
              f"{_pm(c2)}·sin({_g(wd)} t)]"
              + (f"{_pm(A)}·sin({_g(w)} t){_pm(B)}·cos({_g(w)} t)"
                 if mode == "AC" else ""))
    elif disc > 1e-12 * scale:                        # overdamped
        sq = math.sqrt(disc) / (2.0 * Leff)
        r1, r2 = -alpha + sq, -alpha - sq
        c1 = (i0 - r2 * q0) / (r1 - r2)
        c2 = q0 - c1
        yh = c1 * np.exp(r1 * t) + c2 * np.exp(r2 * t)
        yph = c1 * r1 * np.exp(r1 * t) + c2 * r2 * np.exp(r2 * t)
        envh_y = abs(c1) * np.exp(r1 * t) + abs(c2) * np.exp(r2 * t)
        envh_yp = abs(c1 * r1) * np.exp(r1 * t) + abs(c2 * r2) * np.exp(r2 * t)
        alpha_settle = -r1
        sol.update(damping="overdamped", r1=r1, r2=r2)
        fy = (f"{yname}(t) = {_g(c1)}·e^({_g(r1)} t){_pm(c2)}·e^({_g(r2)} t)"
              + (f"{_pm(A)}·sin({_g(w)} t){_pm(B)}·cos({_g(w)} t)"
                 if mode == "AC" else ""))
    else:                                             # critically damped
        r = -alpha
        c1 = q0
        c2 = i0 - r * c1
        e = np.exp(r * t)
        yh = (c1 + c2 * t) * e
        yph = (c2 + r * (c1 + c2 * t)) * e
        envh_y = (abs(c1) + abs(c2) * t) * e
        envh_yp = (abs(c2 + r * c1) + abs(r * c2) * t) * e
        alpha_settle = alpha
        sol.update(damping="critically damped", r1=r)
        fy = (f"{yname}(t) = ({_g(c1)}{_pm(c2)}·t)·e^({_g(r)} t)"
              + (f"{_pm(A)}·sin({_g(w)} t){_pm(B)}·cos({_g(w)} t)"
                 if mode == "AC" else ""))

    if alpha == 0:
        sol["damping"] = "undamped (no damping)"

    y = yh + yss
    yp = yph + ypss
    spread_y = (math.hypot(A, B) if mode == "AC" else 0.0) + envh_y
    spread_yp = (math.hypot(A * w, B * w) if mode == "AC" else 0.0) + envh_yp
    sol.update(y=y, yp=yp, env_y_hi=spread_y, env_y_lo=-spread_y,
               env_yp_hi=spread_yp, env_yp_lo=-spread_yp,
               alpha_settle=alpha_settle, fy=fy + f"   [Volt]")
    return sol


def solve_parallel_rc(R, Cuf, I0, w, t, mode="AC"):
    """R∥C driven by a current source:  C·V' + V/R = I(t),  V(0)=0."""
    Cf = Cuf * 1e-6
    Rp, kp = Cf, 1.0 / max(R, 1e-9)
    a1 = kp / Rp
    if mode == "AC":
        D = max(kp * kp + (Rp * w) ** 2, 1e-30)
        A = I0 * kp / D
        B = -I0 * Rp * w / D
        Vss = A * np.sin(w * t) + B * np.cos(w * t)
        Vpss = A * w * np.cos(w * t) - B * w * np.sin(w * t)
        K = -B
        e = np.exp(-a1 * t)
        V = Vss + K * e
        Vp = Vpss - K * a1 * e
        envh = abs(K) * e
        spread = math.hypot(A, B) + envh
        fV = (f"V(t) = {_g(K)}·e^(-{_g(a1)} t){_pm(A)}·sin({_g(w)} t)"
              f"{_pm(B)}·cos({_g(w)} t)   [Volt]")
    else:
        Vp_final = I0 / kp
        e = np.exp(-a1 * t)
        V = Vp_final * (1.0 - e)
        Vp = Vp_final * a1 * e
        envh = Vp_final * e
        spread = envh
        fV = f"V(t) = {_g(Vp_final)}·(1 − e^(-{_g(a1)} t))   [Volt]"

    I_R = V / max(R, 1e-9)
    I_C = Cf * Vp
    Src = I0 * np.sin(w * t) if mode == "AC" else np.full_like(t, I0)
    fI = f"I_C(t) = C·V′(t)   [Ampere]   (I_R(t) = V(t)/R)"
    return dict(V=V, Vp=Vp, Src=Src, src_kind="I", I_R=I_R, I_L=None,
                I_C=I_C, I_total=Src, Q_C=Cf * V,
                alpha=a1, alpha_settle=a1, damping="1st-order (τ = RC)",
                w0=0.0, env_y_hi=spread, env_y_lo=-spread,
                fV=fV, fI=fI)


def solve_parallel_rl(R, L, I0, w, t, mode="AC"):
    """R∥L driven by a current source:  (L/R)·I_L' + I_L = I(t),  I_L(0)=0."""
    Rr = max(R, 1e-9)
    Rp, kp = L / Rr, 1.0
    a1 = kp / Rp
    if mode == "AC":
        D = max(kp * kp + (Rp * w) ** 2, 1e-30)
        A = I0 * kp / D
        B = -I0 * Rp * w / D
        ILss = A * np.sin(w * t) + B * np.cos(w * t)
        ILpss = A * w * np.cos(w * t) - B * w * np.sin(w * t)
        K = -B
        e = np.exp(-a1 * t)
        I_L = ILss + K * e
        ILp = ILpss - K * a1 * e
        envh = abs(K) * e
        spread = math.hypot(A, B) + envh
        fI = (f"I_L(t) = {_g(K)}·e^(-{_g(a1)} t){_pm(A)}·sin({_g(w)} t)"
              f"{_pm(B)}·cos({_g(w)} t)   [Ampere]")
    else:
        IL_final = I0 / kp
        e = np.exp(-a1 * t)
        I_L = IL_final * (1.0 - e)
        ILp = IL_final * a1 * e
        envh = IL_final * e
        spread = envh
        fI = f"I_L(t) = {_g(IL_final)}·(1 − e^(-{_g(a1)} t))   [Ampere]"

    V = L * ILp
    I_R = V / Rr
    Src = I0 * np.sin(w * t) if mode == "AC" else np.full_like(t, I0)
    fV = f"V(t) = L·I_L′(t)   [Volt]"
    return dict(V=V, Vp=None, Src=Src, src_kind="I", I_R=I_R, I_L=I_L,
                I_C=None, I_total=Src, Q_C=None,
                alpha=a1, alpha_settle=a1, damping="1st-order (τ = L/R)",
                w0=0.0, env_y_hi=spread, env_y_lo=-spread,
                fV=fV, fI=fI)


def solve_parallel_rlc(R, L, Cuf, I0, w, t, mode="AC"):
    """R∥L∥C driven by a current source (true parallel resonance):

        C·V'' + V'/R + V/L = I'(t)
    """
    Cf = Cuf * 1e-6
    Rr = max(R, 1e-9)
    if mode == "AC":
        gs = solve_second_general(Cf, 1.0 / Rr, 1.0 / L, t, mode="AC", w=w,
                                  F0=I0 * w, phase=math.pi / 2, yname="V")
    else:
        gs = solve_second_general(Cf, 1.0 / Rr, 1.0 / L, t, mode="DC",
                                  y0=0.0, yp0=I0 / Cf, yname="V")
    V, Vp = gs["y"], gs["yp"]
    Src = I0 * np.sin(w * t) if mode == "AC" else np.full_like(t, I0)
    I_C = Cf * Vp
    I_R = V / Rr
    I_L = Src - I_R - I_C

    Zp = phi_p = None
    if mode == "AC":
        Yp = complex(1.0 / Rr, w * Cf - 1.0 / (w * L))
        Zpc = 1.0 / Yp
        Zp, phi_p = abs(Zpc), math.degrees(math.atan2(Zpc.imag, Zpc.real))

    fI = f"I_L(t) = I(t) − I_R(t) − I_C(t)   [Ampere]   (KCL)"
    out = dict(V=V, Vp=Vp, Src=Src, src_kind="I", I_R=I_R, I_L=I_L, I_C=I_C,
               I_total=Src, Q_C=Cf * V, alpha=gs["alpha"],
               alpha_settle=gs["alpha_settle"], damping=gs["damping"],
               w0=gs["w0"], env_y_hi=gs["env_y_hi"], env_y_lo=gs["env_y_lo"],
               fV=gs["fy"], fI=fI, Zp=Zp, phi_p=phi_p)
    if "wd" in gs:
        out["wd"] = gs["wd"]
    if "r1" in gs:
        out["r1"] = gs["r1"]
    if "r2" in gs:
        out["r2"] = gs["r2"]
    return out


def solve_tank(R, L, Cuf, E0, w, t, mode="AC"):
    """Tank: R in series with an L∥C tank, driven by a voltage source:

        C·V_t'' + V_t'/R + V_t/L = E'(t)/R
    """
    Cf = Cuf * 1e-6
    Rr = max(R, 1e-9)
    if mode == "AC":
        gs = solve_second_general(Cf, 1.0 / Rr, 1.0 / L, t, mode="AC", w=w,
                                  F0=E0 * w / Rr, phase=math.pi / 2,
                                  yname="V_t")
    else:
        gs = solve_second_general(Cf, 1.0 / Rr, 1.0 / L, t, mode="DC",
                                  y0=0.0, yp0=E0 / (Rr * Cf), yname="V_t")
    Vt, Vtp = gs["y"], gs["yp"]
    Esrc = E0 * np.sin(w * t) if mode == "AC" else np.full_like(t, E0)
    I_C = Cf * Vtp
    I_total = (Esrc - Vt) / Rr                # current through R = into the tank
    I_L = I_total - I_C

    Ztot = phi_tot = None
    if mode == "AC":
        Ztank = 1.0 / complex(0.0, w * Cf - 1.0 / (w * L)) if \
            abs(w * Cf - 1.0 / (w * L)) > 1e-30 else complex(1e18, 0.0)
        Zt = complex(Rr, 0.0) + Ztank
        Ztot, phi_tot = abs(Zt), math.degrees(math.atan2(Zt.imag, Zt.real))

    fI = f"I(t) = (E(t) − V_t(t))/R   [Ampere]   (through R, into the tank)"
    return dict(V=Vt, Vp=Vtp, Src=Esrc, src_kind="E", I_R=I_total, I_L=I_L,
                I_C=I_C, I_total=I_total, Q_C=Cf * Vt, alpha=gs["alpha"],
                alpha_settle=gs["alpha_settle"], damping=gs["damping"],
                w0=gs["w0"], env_y_hi=gs["env_y_hi"], env_y_lo=gs["env_y_lo"],
                fV=gs["fy"], fI=fI, Zp=Ztot, phi_p=phi_tot,
                **({"wd": gs["wd"]} if "wd" in gs else {}),
                **({"r1": gs["r1"], "r2": gs["r2"]} if "r2" in gs else
                   ({"r1": gs["r1"]} if "r1" in gs else {})))


def solve_parallel(preset, R, L, Cuf, E0, w, t, mode="AC"):
    """Dispatch to the right parallel-preset solver."""
    if preset == "RC_P":
        return solve_parallel_rc(R, Cuf, E0, w, t, mode)
    if preset == "RL_P":
        return solve_parallel_rl(R, L, E0, w, t, mode)
    if preset == "RLC_P":
        return solve_parallel_rlc(R, L, Cuf, E0, w, t, mode)
    if preset == "TANK":
        return solve_tank(R, L, Cuf, E0, w, t, mode)
    raise ValueError(f"unknown parallel preset: {preset!r}")


def solve_parallel_rk4(preset, R, L, Cuf, E0, w, t, mode="AC"):
    """RK4 numerical integration as an independent check of the parallel
    analytics. Returns (V, I_total_or_secondary) matching the state
    variable each preset actually solves for, for cross-checking against
    the analytic `V`/`I_L` arrays."""
    Cf = Cuf * 1e-6
    Rr = max(R, 1e-9)
    dt = t[1] - t[0]
    n_sub = max(2, int(math.ceil(dt / 8e-6)))
    h = dt / n_sub

    def src(tt):
        return E0 * math.sin(w * tt) if mode == "AC" else E0

    if preset == "RC_P":
        V = np.empty_like(t)
        V[0] = 0.0
        v = 0.0

        def g(tt, v):
            return (src(tt) - v / Rr) / Cf

        for n in range(1, len(t)):
            tt = t[n - 1]
            for _ in range(n_sub):
                k1 = g(tt, v)
                k2 = g(tt + h / 2, v + h / 2 * k1)
                k3 = g(tt + h / 2, v + h / 2 * k2)
                k4 = g(tt + h, v + h * k3)
                v += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
                tt += h
            V[n] = v
        return V, None

    if preset == "RL_P":
        IL = np.empty_like(t)
        IL[0] = 0.0
        i_l = 0.0

        def g(tt, i_l):
            return (src(tt) - i_l) * Rr / L

        for n in range(1, len(t)):
            tt = t[n - 1]
            for _ in range(n_sub):
                k1 = g(tt, i_l)
                k2 = g(tt + h / 2, i_l + h / 2 * k1)
                k3 = g(tt + h / 2, i_l + h / 2 * k2)
                k4 = g(tt + h, i_l + h * k3)
                i_l += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
                tt += h
            IL[n] = i_l
        return IL, None

    if preset == "RLC_P":
        V = np.empty_like(t)
        Vp = np.empty_like(t)
        V[0] = 0.0
        Vp[0] = E0 / Cf if mode == "DC" else 0.0
        v, vp = V[0], Vp[0]

        def dsrc(tt):
            return E0 * w * math.cos(w * tt) if mode == "AC" else 0.0

        def f(tt, v, vp):
            return vp, (dsrc(tt) - vp / Rr - v / L) / Cf

        for n in range(1, len(t)):
            tt = t[n - 1]
            for _ in range(n_sub):
                k1v, k1p = f(tt, v, vp)
                k2v, k2p = f(tt + h / 2, v + h / 2 * k1v, vp + h / 2 * k1p)
                k3v, k3p = f(tt + h / 2, v + h / 2 * k2v, vp + h / 2 * k2p)
                k4v, k4p = f(tt + h, v + h * k3v, vp + h * k3p)
                v += h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
                vp += h / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
                tt += h
            V[n], Vp[n] = v, vp
        return V, Vp

    if preset == "TANK":
        Vt = np.empty_like(t)
        Vtp = np.empty_like(t)
        Vt[0] = 0.0
        Vtp[0] = E0 / (Rr * Cf) if mode == "DC" else 0.0
        v, vp = Vt[0], Vtp[0]

        def dsrc(tt):
            return (E0 * w * math.cos(w * tt) if mode == "AC" else 0.0) / Rr

        def f(tt, v, vp):
            return vp, (dsrc(tt) - vp / Rr - v / L) / Cf

        for n in range(1, len(t)):
            tt = t[n - 1]
            for _ in range(n_sub):
                k1v, k1p = f(tt, v, vp)
                k2v, k2p = f(tt + h / 2, v + h / 2 * k1v, vp + h / 2 * k1p)
                k3v, k3p = f(tt + h / 2, v + h / 2 * k2v, vp + h / 2 * k2p)
                k4v, k4p = f(tt + h, v + h * k3v, vp + h * k3p)
                v += h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
                vp += h / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
                tt += h
            Vt[n], Vtp[n] = v, vp
        return Vt, Vtp

    raise ValueError(f"unknown parallel preset: {preset!r}")
