# -*- coding: utf-8 -*-
"""Series circuit schematic: the R, L, C components can be hidden per
topology (replaced by a straight wire) so the drawing always matches the
equation being solved."""

import numpy as np

from rlc_config import TH, C_SRC, C_R, C_L, C_C


class Schematic:
    def __init__(self, ax):
        self.ax = ax
        ax.set_xlim(-0.4, 10.4)
        ax.set_ylim(-0.7, 6.4)
        ax.set_aspect("equal")
        ax.axis("off")
        wire = dict(color="#334155", lw=1.6, solid_capstyle="round")

        # -- always visible: AC source, frame wires, switch --------------------
        th = np.linspace(0, 2 * np.pi, 80)
        ax.plot(1 + 0.45 * np.cos(th), 3 + 0.45 * np.sin(th),
                color=C_SRC, lw=1.6)
        ts = np.linspace(-0.3, 0.3, 60)
        ax.plot(1 + ts, 3 + 0.16 * np.sin(ts / 0.3 * np.pi),
                color=C_SRC, lw=1.2)
        ax.plot([1, 1], [1, 2.55], **wire)
        ax.plot([1, 1], [3.45, 5], **wire)
        ax.plot([1, 2.2], [5, 5], **wire)
        ax.plot([4.2, 5.2], [5, 5], **wire)
        ax.plot([7.6, 9], [5, 5], **wire)
        ax.plot([9, 9], [5, 3.25], **wire)
        ax.plot([9, 9], [2.75, 1], **wire)
        ax.plot([9, 5.6], [1, 1], **wire)
        ax.plot([4.4, 1], [1, 1], **wire)
        ax.plot([4.4, 5.45], [1, 1.45], color="#334155", lw=1.5)
        ax.plot([4.4, 5.6], [1, 1], "o", color="#334155", ms=3.5)
        ax.text(5.0, 1.78, "switch: t = 0", fontsize=7.2, ha="center",
                color=TH["sub"])

        # -- resistor (zigzag) + bypass wire ------------------------------------
        xs = np.concatenate(([2.2], np.linspace(2.35, 4.05, 7), [4.2]))
        ys = [5, 5.28, 4.72, 5.28, 4.72, 5.28, 4.72, 5.28, 5]
        zig, = ax.plot(xs, ys, color=C_R, lw=1.7, solid_capstyle="round")
        self.grp_R = [zig]
        self.wire_R, = ax.plot([2.2, 4.2], [5, 5], **wire)

        # -- inductor (four arcs) + bypass wire ----------------------------------
        tha = np.linspace(0, np.pi, 40)
        self.grp_L = []
        for kk in range(4):
            cx = 5.5 + kk * 0.6
            arc, = ax.plot(cx + 0.3 * np.cos(tha[::-1]),
                           5 + 0.35 * np.sin(tha),
                           color=C_L, lw=1.7, solid_capstyle="round")
            self.grp_L.append(arc)
        self.wire_L, = ax.plot([5.2, 7.6], [5, 5], **wire)

        # -- capacitor (two plates) + gap-bridging wire ---------------------------
        p1, = ax.plot([8.55, 9.45], [3.25, 3.25], color=C_C, lw=2.0)
        p2, = ax.plot([8.55, 9.45], [2.75, 2.75], color=C_C, lw=2.0)
        self.grp_C = [p1, p2]
        self.wire_C, = ax.plot([9, 9], [3.25, 2.75], **wire)

        # -- component value labels ------------------------------------------------
        self.labels = dict(
            R=ax.text(3.2, 5.66, "", fontsize=8.4, ha="center",
                      color=C_R, weight="bold"),
            L=ax.text(6.4, 5.66, "", fontsize=8.4, ha="center",
                      color=C_L, weight="bold"),
            C=ax.text(8.35, 3.0, "", fontsize=8.4, ha="right", va="center",
                      color=C_C, weight="bold"),
            E=ax.text(5.0, -0.25, "", fontsize=8.4, ha="center",
                      color=C_SRC, weight="bold"))

    def set_topology(self, has_R, has_L, has_C):
        """Show the present components; replace absent ones with plain wire."""
        for art in self.grp_R:
            art.set_visible(has_R)
        self.wire_R.set_visible(not has_R)
        self.labels["R"].set_visible(has_R)
        for art in self.grp_L:
            art.set_visible(has_L)
        self.wire_L.set_visible(not has_L)
        self.labels["L"].set_visible(has_L)
        for art in self.grp_C:
            art.set_visible(has_C)
        self.wire_C.set_visible(not has_C)
        self.labels["C"].set_visible(has_C)


# ---- vertical-branch component primitives (parallel schematic) ----------------
_WIRE = dict(color="#334155", lw=1.6, solid_capstyle="round")


def _vert_resistor(ax, x, ytop, ybot):
    span = ytop - ybot
    y_hi, y_lo = ytop - 0.15 * span, ytop - 0.85 * span
    ax.plot([x, x], [ytop, y_hi], **_WIRE)
    yz = np.linspace(y_hi, y_lo, 7)
    xz = [x + (0.28 if i % 2 == 0 else -0.28) for i in range(7)]
    ax.plot([x] + xz + [x], [y_hi] + list(yz) + [y_lo],
            color=C_R, lw=1.7, solid_capstyle="round")
    ax.plot([x, x], [y_lo, ybot], **_WIRE)


def _vert_inductor(ax, x, ytop, ybot):
    span = ytop - ybot
    y_hi = ytop - 0.15 * span
    y_lo = y_hi - 2.4
    ax.plot([x, x], [ytop, y_hi], **_WIRE)
    tha = np.linspace(0, np.pi, 40)
    for kk in range(4):
        cy = y_hi - 0.3 - kk * 0.6
        ax.plot(x + 0.32 * np.sin(tha), cy + 0.3 * np.cos(tha[::-1]),
                color=C_L, lw=1.7, solid_capstyle="round")
    ax.plot([x, x], [y_lo, ybot], **_WIRE)


def _vert_capacitor(ax, x, ytop, ybot):
    ymid = (ytop + ybot) / 2.0
    yt, yb = ymid + 0.25, ymid - 0.25
    ax.plot([x, x], [ytop, yt], **_WIRE)
    ax.plot([x - 0.45, x + 0.45], [yt, yt], color=C_C, lw=2.0)
    ax.plot([x - 0.45, x + 0.45], [yb, yb], color=C_C, lw=2.0)
    ax.plot([x, x], [yb, ybot], **_WIRE)


def _current_source(ax, cx=1.0, cy=3.0, r=0.45):
    th = np.linspace(0, 2 * np.pi, 80)
    ax.plot(cx + r * np.cos(th), cy + r * np.sin(th), color=C_SRC, lw=1.6)
    ax.annotate("", xy=(cx, cy + 0.28), xytext=(cx, cy - 0.28),
                arrowprops=dict(arrowstyle="-|>", color=C_SRC, lw=1.6,
                                mutation_scale=13))


def _voltage_source(ax, cx=1.0, cy=3.0, r=0.45):
    th = np.linspace(0, 2 * np.pi, 80)
    ax.plot(cx + r * np.cos(th), cy + r * np.sin(th), color=C_SRC, lw=1.6)
    ts = np.linspace(-0.3, 0.3, 60)
    ax.plot(cx + ts, cy + 0.16 * np.sin(ts / 0.3 * np.pi), color=C_SRC, lw=1.2)


class ParallelSchematic:
    """Schematic for the parallel presets (Milestone 2). Since the branch
    layout differs a lot between a 2-branch, 3-branch, and the series-R-
    into-a-tank drawing, this redraws the small schematic axes from scratch
    on every preset switch rather than trying to hide/show one fixed set of
    artists — `set_preset` is only called on a topology-button click, so the
    extra draw cost is irrelevant."""

    def __init__(self, ax):
        self.ax = ax
        self.labels = {}

    def _reset_axes(self):
        ax = self.ax
        ax.clear()
        ax.set_xlim(-0.4, 10.4)
        ax.set_ylim(-0.7, 6.4)
        ax.set_aspect("equal")
        ax.axis("off")

    def _rails_and_switch(self, ax, x_last_branch=9.0):
        """Bottom rail (with switch) common to every parallel preset."""
        ax.plot([1, 1], [1, 2.55], **_WIRE)
        ax.plot([x_last_branch, 5.6], [1, 1], **_WIRE)
        ax.plot([4.4, 1], [1, 1], **_WIRE)
        ax.plot([4.4, 5.45], [1, 1.45], color="#334155", lw=1.5)
        ax.plot([4.4, 5.6], [1, 1], "o", color="#334155", ms=3.5)
        ax.text(5.0, 1.78, "switch: t = 0", fontsize=7.2, ha="center",
                color=TH["sub"])

    def set_preset(self, preset):
        self._reset_axes()
        ax = self.ax
        if preset == "TANK":
            self._draw_tank(ax)
        else:
            self._draw_direct(ax, preset)

    def _draw_direct(self, ax, preset):
        """R∥C, R∥L, R∥L∥C: current source with vertical branches across it."""
        _current_source(ax)
        ax.plot([1, 1], [3.45, 5], **_WIRE)
        branches = {"RC_P": ["R", "C"], "RL_P": ["R", "L"],
                    "RLC_P": ["R", "L", "C"]}[preset]
        n = len(branches)
        xs = (np.linspace(3.4, 8.6, n) if n > 1 else [6.0])
        ax.plot([1, xs[0]], [5, 5], **_WIRE)
        for a, b in zip(xs[:-1], xs[1:]):
            ax.plot([a, b], [5, 5], **_WIRE)
        self._rails_and_switch(ax, xs[-1])

        draw_fn = {"R": _vert_resistor, "L": _vert_inductor,
                  "C": _vert_capacitor}
        colors = {"R": C_R, "L": C_L, "C": C_C}
        self.labels = {}
        for kind, x in zip(branches, xs):
            draw_fn[kind](ax, x, 5.0, 1.0)
            self.labels[kind] = ax.text(x, 5.66, "", fontsize=8.2,
                                        ha="center", va="bottom",
                                        color=colors[kind], weight="bold")
        self.labels["E"] = ax.text(1.0, -0.25, "", fontsize=8.4, ha="center",
                                   color=C_SRC, weight="bold")

    def _draw_tank(self, ax):
        """R in series, then the top rail splits into an L∥C tank."""
        _voltage_source(ax, cx=1.0, cy=3.0)
        ax.plot([1, 1], [3.45, 5], **_WIRE)
        ax.plot([1, 2.2], [5, 5], **_WIRE)
        xs_r = np.concatenate(([2.2], np.linspace(2.35, 4.05, 7), [4.2]))
        ys_r = [5, 5.28, 4.72, 5.28, 4.72, 5.28, 4.72, 5.28, 5]
        ax.plot(xs_r, ys_r, color=C_R, lw=1.7, solid_capstyle="round")
        ax.plot([4.2, 6.5], [5, 5], **_WIRE)
        ax.plot([6.5, 8.5], [5, 5], **_WIRE)          # bridges the tank tops
        self._rails_and_switch(ax, 8.5)
        _vert_inductor(ax, 6.5, 5.0, 1.0)
        _vert_capacitor(ax, 8.5, 5.0, 1.0)
        self.labels = dict(
            R=ax.text(3.2, 5.66, "", fontsize=8.2, ha="center", color=C_R,
                      weight="bold"),
            L=ax.text(6.5, 5.66, "", fontsize=8.2, ha="center", color=C_L,
                      weight="bold"),
            C=ax.text(8.5, 5.66, "", fontsize=8.2, ha="center", color=C_C,
                      weight="bold"),
            E=ax.text(1.0, -0.25, "", fontsize=8.4, ha="center",
                      color=C_SRC, weight="bold"))
