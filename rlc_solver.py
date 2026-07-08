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
