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
