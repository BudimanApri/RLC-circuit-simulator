# -*- coding: utf-8 -*-
"""
Netlist data model (Milestone 3).

A `Netlist` is just a bag of `Component`s connecting named nodes, plus a
designated ground node (fixed at V=0). This is the same representation
SPICE-class simulators use, and it is what lets `rlc_mna.py` solve *any*
wiring instead of the five fixed series presets and four fixed parallel
presets `rlc_solver.py` already knows in closed form.

Component kinds:
    "R"    resistor,        value = ohms
    "L"    inductor,        value = henries
    "C"    capacitor,       value = farads
    "VSRC" voltage source,  value = amplitude (volts), node_a = "+" terminal
    "ISRC" current source,  value = amplitude (amps),  flows node_a -> node_b

Sources carry `source_type` ("AC" or "DC") and, for AC, `freq` (rad/s):
    AC:  x(t) = value * sin(freq * t)
    DC:  x(t) = value                 for all t >= 0 (switch-closure step)

This module is pure Python/dataclasses — no numpy, no matplotlib — so it
stays trivially testable and JSON-serializable on its own.
"""

import json
import math
from dataclasses import dataclass, field, asdict

GROUND = "0"
_KINDS = ("R", "L", "C", "VSRC", "ISRC")
_SOURCE_KINDS = ("VSRC", "ISRC")


@dataclass
class Component:
    kind: str            # "R" | "L" | "C" | "VSRC" | "ISRC"
    node_a: str
    node_b: str
    value: float
    name: str = ""
    source_type: str = "DC"     # "AC" or "DC"; only meaningful for sources
    freq: float = 0.0           # rad/s; only meaningful for AC sources

    def __post_init__(self):
        if self.kind not in _KINDS:
            raise ValueError(f"unknown component kind: {self.kind!r}")
        if self.node_a == self.node_b:
            raise ValueError(f"{self.name or self.kind}: node_a == node_b "
                             f"({self.node_a!r}) — a component can't short "
                             f"a node to itself")
        if self.is_source() and self.source_type not in ("AC", "DC"):
            raise ValueError(f"unknown source_type: {self.source_type!r}")

    def is_source(self):
        return self.kind in _SOURCE_KINDS

    def waveform(self, t):
        """Evaluate this source at time(s) t (scalar or numpy array).
        Non-sources raise — only call this on VSRC/ISRC components."""
        if not self.is_source():
            raise ValueError(f"{self.kind} is not a source")
        if self.source_type == "AC":
            try:
                import numpy as np
                return self.value * np.sin(self.freq * t)
            except ImportError:
                return self.value * math.sin(self.freq * t)
        try:
            import numpy as np
            return self.value * np.ones_like(t, dtype=float) \
                if hasattr(t, "__len__") else self.value
        except ImportError:
            return self.value


@dataclass
class Netlist:
    components: list = field(default_factory=list)
    ground: str = GROUND

    def add(self, kind, node_a, node_b, value, name="", source_type="DC",
            freq=0.0):
        """Convenience constructor: build, append, and return a Component."""
        if not name:
            n_existing = sum(1 for c in self.components if c.kind == kind)
            name = f"{kind}{n_existing + 1}"
        c = Component(kind, node_a, node_b, value, name, source_type, freq)
        self.components.append(c)
        return c

    @property
    def nodes(self):
        """Sorted set of every node name referenced, ground included."""
        ns = {self.ground}
        for c in self.components:
            ns.add(c.node_a)
            ns.add(c.node_b)
        return sorted(ns)

    @property
    def non_ground_nodes(self):
        return [n for n in self.nodes if n != self.ground]

    def sources(self):
        return [c for c in self.components if c.is_source()]

    def validate(self):
        """Cheap sanity checks — catches the mistakes that would otherwise
        surface as a singular MNA matrix. Not full DRC (dangling-wire /
        short-circuit detection is a Milestone 4 UI concern); this only
        guards against inputs the solver flat-out cannot handle."""
        if not self.components:
            raise ValueError("netlist has no components")
        if self.ground not in self.nodes:
            raise ValueError(f"ground node {self.ground!r} is never used "
                             f"by any component")
        touches = {n: 0 for n in self.non_ground_nodes}
        for c in self.components:
            for n in (c.node_a, c.node_b):
                if n in touches:
                    touches[n] += 1
        floating = [n for n, cnt in touches.items() if cnt < 2]
        if floating:
            raise ValueError(f"floating node(s) with only one connection: "
                             f"{floating} — every node needs >=2 component "
                             f"terminals (or a path to ground) to be solvable")

    # ---- JSON round-trip ---------------------------------------------------
    def to_dict(self):
        return dict(ground=self.ground,
                   components=[asdict(c) for c in self.components])

    def to_json(self, indent=2):
        return json.dumps(self.to_dict(), indent=indent)

    @classmethod
    def from_dict(cls, d):
        nl = cls(ground=d.get("ground", GROUND))
        nl.components = [Component(**c) for c in d["components"]]
        return nl

    @classmethod
    def from_json(cls, s):
        return cls.from_dict(json.loads(s))


# ---- convenience builders for the topologies this app already knows -------
# These mirror rlc_config.TOPOS / PARALLEL_TOPOS exactly, and exist mainly so
# rlc_mna.py's cross-validation tests (and any future UI code) don't have to
# hand-wire nodes every time. "n1"/"n2" are internal nodes; "0" is ground.

def series_netlist(topo, R, L, Cuf, E0, w, mode="AC"):
    """Series R[-L][-C] loop driven by a voltage source, matching the
    topology semantics of rlc_solver.solve()."""
    from rlc_config import TOPOS
    cfg = TOPOS[topo]
    Cf = Cuf * 1e-6
    nl = Netlist()
    nl.add("VSRC", "n0", "0", E0, name="E", source_type=mode, freq=w)

    chain = []
    if cfg["R"]:
        chain.append(("R", R))
    if cfg["L"]:
        chain.append(("L", L))
    if cfg["C"]:
        chain.append(("C", Cf))
    # every topology has at least one of R/L/C, so chain is never empty

    node = "n0"
    for i, (kind, val) in enumerate(chain):
        nxt = "0" if i == len(chain) - 1 else f"n{i + 1}"
        nl.add(kind, node, nxt, val)
        node = nxt
    return nl


def parallel_netlist(preset, R, L, Cuf, E0, w, mode="AC"):
    """Parallel preset driven the same way rlc_solver.solve_parallel() does:
    R∥C / R∥L / R∥L∥C by a current source, Tank by a voltage source."""
    from rlc_config import PARALLEL_TOPOS
    cfg = PARALLEL_TOPOS[preset]
    Cf = Cuf * 1e-6
    nl = Netlist()
    if preset == "TANK":
        nl.add("VSRC", "n0", "0", E0, name="E", source_type=mode, freq=w)
        nl.add("R", "n0", "n1", R)
        nl.add("L", "n1", "0", L)
        nl.add("C", "n1", "0", Cf)
    else:
        nl.add("ISRC", "0", "n0", E0, name="I", source_type=mode, freq=w)
        nl.add("R", "n0", "0", R)
        if cfg["L"]:
            nl.add("L", "n0", "0", L)
        if cfg["C"]:
            nl.add("C", "n0", "0", Cf)
    return nl
