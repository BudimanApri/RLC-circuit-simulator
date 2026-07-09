# -*- coding: utf-8 -*-
"""
General transient circuit solver via trapezoidal Modified Nodal Analysis
(Milestone 3) — the engine that makes an arbitrary `rlc_netlist.Netlist`
simulatable, instead of only the fixed presets `rlc_solver.py` knows in
closed form.

Method
------
Unknowns: one voltage per non-ground node, plus one branch-current unknown
per voltage source (the "modified" part of MNA). Resistors, capacitors, and
inductors are all reduced to a conductance stamp between two nodes — L and
C via their trapezoidal *companion models*, which turn "dV/dt" / "dI/dt"
into an equivalent (conductance, history current source) pair that only
depends on the previous timestep:

    Capacitor (I = C·dV/dt), branch current a→b:
        Geq = 2C/h,   Ihist = Geq·V_prev + I_prev
        I(t+h) = Geq·V(t+h) − Ihist

    Inductor (V = L·dI/dt), branch current a→b:
        Geq = h/(2L), Ihist = I_prev + Geq·V_prev
        I(t+h) = Geq·V(t+h) + Ihist

Both reduce dynamics to "solve one linear system per step", which is why
they can share the exact same conductance-stamping code path as a plain
resistor — only the history term differs, and it only touches the RHS
vector, never the matrix itself. Since R/L/C values and h are fixed for an
entire sweep, the system matrix G is CONSTANT across all steps: it is
LU-factorized once and reused (via scipy if available, else a fresh
`numpy.linalg.solve` per step — slower, still correct).

Branch current convention: for every 2-terminal component between node_a
and node_b, I is the current flowing node_a -> node_b (associated
reference, V = Va - Vb). This is the same convention as SPICE.

Validated against the closed-form solutions in `rlc_solver.py` for all 9
existing presets (5 series + 4 parallel) in `rlc_simulator.py --test`.
"""

import math
from dataclasses import dataclass, field

import numpy as np

try:
    from scipy.linalg import lu_factor
    from scipy.linalg.lapack import get_lapack_funcs
    _HAVE_SCIPY = True
except ImportError:
    _HAVE_SCIPY = False


@dataclass
class MnaResult:
    t: np.ndarray
    node_v: dict = field(default_factory=dict)   # node name -> array
    comp_i: dict = field(default_factory=dict)   # component name -> array
    comp_v: dict = field(default_factory=dict)   # component name -> array
    netlist: object = None

    def component(self, name):
        return next(c for c in self.netlist.components if c.name == name)


class _Group:
    """Same-kind components' node indices/values as plain numpy arrays, so
    the per-timestep work is a handful of vectorized numpy calls instead of
    a Python-level loop with per-component dict lookups (profiling showed
    that loop — not the linear solve itself — was >75% of the wall time)."""

    __slots__ = ("comp", "ia", "ib", "gidx", "val")

    def __init__(self, comps, node_idx, size, all_names):
        self.comp = comps
        dump = size                     # ground -> this row/col, discarded
        self.ia = np.array([node_idx.get(c.node_a, dump) for c in comps],
                           dtype=np.intp)
        self.ib = np.array([node_idx.get(c.node_b, dump) for c in comps],
                           dtype=np.intp)
        self.val = np.array([c.value for c in comps], dtype=float)
        self.gidx = np.array([all_names[c.name] for c in comps],
                             dtype=np.intp)


def _stamp_conductance(G, ia, ib, g, size):
    """Add conductance g[k] between ia[k]/ib[k] for every k, vectorized.
    Ground terminals point at row/col `size` (discarded via slicing after)."""
    n = size + 1
    np.add.at(G.reshape(-1), ia * n + ia, g)
    np.add.at(G.reshape(-1), ib * n + ib, g)
    np.add.at(G.reshape(-1), ia * n + ib, -g)
    np.add.at(G.reshape(-1), ib * n + ia, -g)


def simulate(netlist, t):
    """Run a trapezoidal-MNA transient simulation over the time array `t`
    (must be uniformly spaced, t[0] == 0 — matches this app's "switch
    closes at t=0, everything starts at zero" convention)."""
    netlist.validate()
    h = float(t[1] - t[0])
    t = np.asarray(t, dtype=float)
    nT = len(t)

    nodes = netlist.non_ground_nodes
    node_idx = {n: i for i, n in enumerate(nodes)}
    N = len(nodes)
    all_comps = netlist.components
    all_names = {c.name: i for i, c in enumerate(all_comps)}
    n_comp = len(all_comps)

    Rs = [c for c in all_comps if c.kind == "R"]
    Cs = [c for c in all_comps if c.kind == "C"]
    Ls = [c for c in all_comps if c.kind == "L"]
    Vs = [c for c in all_comps if c.kind == "VSRC"]
    Is = [c for c in all_comps if c.kind == "ISRC"]
    vs_idx = {c.name: N + i for i, c in enumerate(Vs)}
    size = N + len(Vs)                  # unknowns: node voltages + I(VSRC)

    gR = _Group(Rs, node_idx, size, all_names)
    gC = _Group(Cs, node_idx, size, all_names)
    gL = _Group(Ls, node_idx, size, all_names)
    gI = _Group(Is, node_idx, size, all_names)
    # VSRC needs its own extra-unknown row/col too, built like a Group plus k
    v_ia = np.array([node_idx.get(c.node_a, size) for c in Vs], dtype=np.intp)
    v_ib = np.array([node_idx.get(c.node_b, size) for c in Vs], dtype=np.intp)
    v_k = np.array([vs_idx[c.name] for c in Vs], dtype=np.intp)
    v_gidx = np.array([all_names[c.name] for c in Vs], dtype=np.intp)

    # every source's waveform is a pure function of t (no solver-state
    # dependency) -> precompute it once as a full array, not per-step
    def wave_matrix(comps):
        if not comps:
            return np.zeros((0, nT))
        return np.array([np.asarray(c.waveform(t), dtype=float)
                         for c in comps])

    I_wave = wave_matrix(Is)
    V_wave = wave_matrix(Vs)

    def build_matrix(hstep):
        """G (size+1 square, last row/col a discarded ground dump) and the
        per-L/C companion conductance for step size hstep. Called twice:
        the real h for the ongoing sweep, and hstep -> 0 for the t=0
        initial solve (capacitors -> 0V shorts, inductors -> opens, i.e.
        this app's usual Q(0)=0 / I(0)=0 convention)."""
        G = np.zeros((size + 1, size + 1))
        _stamp_conductance(G, gR.ia, gR.ib, 1.0 / gR.val, size)
        geqC = 2.0 * gC.val / hstep
        _stamp_conductance(G, gC.ia, gC.ib, geqC, size)
        geqL = hstep / (2.0 * gL.val)
        _stamp_conductance(G, gL.ia, gL.ib, geqL, size)
        if len(Vs):
            G[v_ia, v_k] += 1.0
            G[v_k, v_ia] += 1.0
            G[v_ib, v_k] -= 1.0
            G[v_k, v_ib] -= 1.0
        return G[:size, :size], geqC, geqL

    def make_solver(G):
        if _HAVE_SCIPY:
            # Call LAPACK's getrs directly instead of scipy.linalg.lu_solve:
            # for a system this small (a handful to a few dozen unknowns),
            # lu_solve's own argument-checking wrapper costs more than the
            # solve itself — raw getrs measured ~4-5x faster per call.
            lu, piv = lu_factor(G, check_finite=False)
            getrs, = get_lapack_funcs(("getrs",), (lu,))
            return lambda b: getrs(lu, piv, b, 0, 0)[0]
        return lambda b: np.linalg.solve(G, b)

    node_v_arr = np.empty((N, nT))
    comp_i_arr = np.empty((n_comp, nT))
    comp_v_arr = np.empty((n_comp, nT))

    b_buf = np.zeros(size + 1)           # reused every step (no realloc)
    xg_buf = np.zeros(size + 1)          # xg_buf[size] stays 0 = ground

    def step(solve_fn, geqC, geqL, hist_vC, hist_iC, hist_vL, hist_iL, n):
        b_buf[:] = 0.0
        if len(Is):
            vals = I_wave[:, n]
            np.add.at(b_buf, gI.ia, -vals)
            np.add.at(b_buf, gI.ib, vals)
        if len(Vs):
            b_buf[v_k] = V_wave[:, n]
        ihC = geqC * hist_vC + hist_iC
        if len(Cs):
            np.add.at(b_buf, gC.ia, ihC)
            np.add.at(b_buf, gC.ib, -ihC)
        ihL = hist_iL + geqL * hist_vL
        if len(Ls):
            np.add.at(b_buf, gL.ia, -ihL)
            np.add.at(b_buf, gL.ib, ihL)

        x = solve_fn(b_buf[:size])
        xg_buf[:size] = x                # xg_buf[size] (ground) stays 0.0
        xg = xg_buf
        node_v_arr[:, n] = x[:N]

        if len(Rs):
            v = xg[gR.ia] - xg[gR.ib]
            comp_v_arr[gR.gidx, n] = v
            comp_i_arr[gR.gidx, n] = v / gR.val
        v_c = xg[gC.ia] - xg[gC.ib] if len(Cs) else hist_vC
        i_c = geqC * v_c - ihC
        if len(Cs):
            comp_v_arr[gC.gidx, n] = v_c
            comp_i_arr[gC.gidx, n] = i_c
        v_l = xg[gL.ia] - xg[gL.ib] if len(Ls) else hist_vL
        i_l = geqL * v_l + ihL
        if len(Ls):
            comp_v_arr[gL.gidx, n] = v_l
            comp_i_arr[gL.gidx, n] = i_l
        if len(Vs):
            comp_v_arr[v_gidx, n] = xg[v_ia] - xg[v_ib]
            comp_i_arr[v_gidx, n] = x[v_k]
        if len(Is):
            comp_v_arr[gI.gidx, n] = xg[gI.ia] - xg[gI.ib]
            comp_i_arr[gI.gidx, n] = I_wave[:, n]

        return v_c, i_c, v_l, i_l

    zeros_c = np.zeros(len(Cs))
    zeros_l = np.zeros(len(Ls))

    # ---- t=0 initial solve (see build_matrix docstring) --------------------
    G0, geqC0, geqL0 = build_matrix(h * 1e-9)
    hv_c, hi_c, hv_l, hi_l = step(make_solver(G0), geqC0, geqL0,
                                  zeros_c, zeros_c, zeros_l, zeros_l, 0)

    # ---- the ongoing trapezoidal sweep, t[1] .. t[-1] ------------------------
    G, geqC, geqL = build_matrix(h)
    solve_fn = make_solver(G)
    for n in range(1, nT):
        hv_c, hi_c, hv_l, hi_l = step(solve_fn, geqC, geqL,
                                      hv_c, hi_c, hv_l, hi_l, n)

    node_v = {nodes[i]: node_v_arr[i] for i in range(N)}
    node_v[netlist.ground] = np.zeros(nT)
    comp_i = {all_comps[i].name: comp_i_arr[i] for i in range(n_comp)}
    comp_v = {all_comps[i].name: comp_v_arr[i] for i in range(n_comp)}
    return MnaResult(t=t, node_v=node_v, comp_i=comp_i, comp_v=comp_v,
                     netlist=netlist)


# ---- probes: thin, self-documenting accessors ------------------------------
def probe_voltage(result, node):
    """Voltage at any node, including the ground node (always 0)."""
    return result.node_v[node]


def probe_current(result, comp_name):
    """Branch current through any component, node_a -> node_b."""
    return result.comp_i[comp_name]


def probe_charge(result, comp_name):
    """Charge on a capacitor: Q = C*V."""
    c = result.component(comp_name)
    if c.kind != "C":
        raise ValueError(f"{comp_name} is a {c.kind}, not a capacitor")
    return c.value * result.comp_v[comp_name]


def probe_energy(result, comp_name):
    """Energy in Joules: instantaneous stored energy for L/C
    (½LI² / ½CV²), cumulative dissipated energy for R (∫I²R dt, via the
    trapezoidal rule so it uses the same time discretization as the sweep
    itself)."""
    c = result.component(comp_name)
    if c.kind == "C":
        return 0.5 * c.value * result.comp_v[comp_name] ** 2
    if c.kind == "L":
        return 0.5 * c.value * result.comp_i[comp_name] ** 2
    if c.kind == "R":
        p = result.comp_i[comp_name] ** 2 * c.value
        dt = np.diff(result.t)
        inc = 0.5 * (p[:-1] + p[1:]) * dt
        return np.concatenate(([0.0], np.cumsum(inc)))
    raise ValueError(f"no energy defined for source component {comp_name}")
