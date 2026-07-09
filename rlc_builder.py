# -*- coding: utf-8 -*-
"""
Free-form circuit builder (Milestone 4) — a matplotlib UI on top of the
Milestone 3 engine (rlc_netlist.Netlist + rlc_mna.simulate). Place R, L, C,
and voltage/current sources on a snapped grid, wire them together, pick a
ground reference, probe any node's voltage or any component's current, and
watch the results update live.

How to run:
    python rlc_builder.py

Interaction model (click-click, not drag — far more robust in matplotlib
than continuous dragging, and orientation falls out for free: whichever
adjacent grid point you click second determines the direction):

    1. Pick a tool from the palette row: R, L, C, VSRC, ISRC, Wire, Ground,
       Select, Delete.
    2. Components (R/L/C/VSRC/ISRC): click a grid point, then click an
       *adjacent* grid point — horizontal, vertical, or 45-degree diagonal —
       to place it between them. Clicking a non-adjacent point instead just
       moves your anchor there.
       Wires: click a start point, then click *any other* point in a
       straight line from it (not necessarily adjacent) — every grid point
       in between is wired automatically, so a long run takes two clicks
       instead of one per hop.
    3. Ground: click any grid point to designate it "0" (moves if you click
       elsewhere).
    4. Select: click a component to edit its value (and, for sources,
       source type/frequency, and polarity via "flip") in the PROPERTIES
       card, and to toggle it into the results as a probed current. Click a
       bare grid point to toggle its voltage into the results.
    5. Delete: click a component or wire to remove it.

Source polarity: whichever grid point you click *first* when placing a
VSRC/ISRC is node_a — the "+" terminal for VSRC, the tail of the current
arrow for ISRC. This is marked on the schematic (+/− labels, arrow) and can
be reversed after the fact with the "flip polarity" button in PROPERTIES.
DC sources draw as a battery symbol (long/short bars); AC sources draw as a
circle with a sine wiggle. When every source in the circuit is DC, small
grey arrows show the actual (steady-state) current direction through each
component — omitted for AC since direction reverses every half-cycle.

The circuit resolves automatically (via rlc_mna.simulate) after every edit,
as soon as it has a ground and at least one source and validates.
"""

import json
import math
import sys

import numpy as np
import matplotlib

if "--test" in sys.argv:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.widgets import Button, TextBox

from rlc_config import TH, C_R, C_L, C_C, C_SRC
from rlc_netlist import Netlist
from rlc_mna import (simulate, probe_voltage, probe_current, probe_charge,
                     probe_energy)

GX, GY = 13, 9                       # grid columns/rows (0-indexed)
TOOLS = ["R", "L", "C", "VSRC", "ISRC", "WIRE", "GROUND", "SELECT", "DELETE"]
KIND_COLOR = {"R": C_R, "L": C_L, "C": C_C, "VSRC": C_SRC, "ISRC": C_SRC}
PROBE_COLORS = ["#2563eb", "#dc2626", "#059669", "#d97706", "#7c3aed",
               "#0891b2", "#be185d", "#65a30d"]


class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def _draw_component(ax, kind, p1, p2, color, source_type="DC", lw_lead=1.6):
    """Draw one component symbol on the single grid edge p1->p2 (adjacent
    grid points, one grid step apart — horizontal, vertical, or diagonal).
    Works for any direction via a local (along, perpendicular) basis, so
    there is one code path instead of separate per-orientation routines.

    p1 is always node_a: the "+" terminal for VSRC, the tail of the current
    arrow for ISRC (the netlist's own polarity convention). `source_type`
    ("AC"/"DC") only affects VSRC/ISRC symbols: AC draws the classic
    circle-with-sine-wiggle, DC draws a battery bar pair (long/thin near
    node_a "+", short/thick near node_b "−") so the two aren't visually
    confusable — a DC source drawn with an AC sine symbol was a reported
    point of confusion."""
    a, b = np.array(p1, dtype=float), np.array(p2, dtype=float)
    d = b - a
    length = np.hypot(*d)
    u = d / length
    perp = np.array([-u[1], u[0]])
    wire = dict(color="#334155", lw=lw_lead, solid_capstyle="round", zorder=2)

    def pt(t, w=0.0):
        return a + u * t * length + perp * w

    if kind == "WIRE":
        ax.plot([a[0], b[0]], [a[1], b[1]], color="#334155", lw=1.8,
                solid_capstyle="round", zorder=2)
        return

    lead = 0.22
    p_lo, p_hi = pt(lead), pt(1 - lead)
    ax.plot([a[0], p_lo[0]], [a[1], p_lo[1]], **wire)
    ax.plot([p_hi[0], b[0]], [p_hi[1], b[1]], **wire)

    if kind == "R":
        n = 6
        pts = [p_lo]
        for i in range(1, n):
            tt = lead + (1 - 2 * lead) * i / n
            w = 0.15 if i % 2 else -0.15
            pts.append(pt(tt, w))
        pts.append(p_hi)
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=color,
                lw=1.8, solid_capstyle="round", zorder=3)
    elif kind == "L":
        n_arcs = 3
        seg = (1 - 2 * lead) / n_arcs
        th = np.linspace(0, np.pi, 20)
        for k in range(n_arcs):
            t0 = lead + k * seg
            ts = t0 + seg * (th / np.pi)
            ws = 0.16 * np.sin(th)
            xs = [pt(tt, ww)[0] for tt, ww in zip(ts, ws)]
            ys = [pt(tt, ww)[1] for tt, ww in zip(ts, ws)]
            ax.plot(xs, ys, color=color, lw=1.8, solid_capstyle="round",
                    zorder=3)
    elif kind == "C":
        g = 0.07
        ca, cb = pt(0.5 - g), pt(0.5 + g)
        pw = 0.22
        for c in (ca, cb):
            p1_ = c + perp * pw
            p2_ = c - perp * pw
            ax.plot([p1_[0], p2_[0]], [p1_[1], p2_[1]], color=color, lw=2.4,
                    zorder=3)
        ax.plot([p_lo[0], ca[0]], [p_lo[1], ca[1]], **wire)
        ax.plot([p_hi[0], cb[0]], [p_hi[1], cb[1]], **wire)
    elif kind in ("VSRC", "ISRC"):
        c = pt(0.5)
        r = 0.28
        th = np.linspace(0, 2 * np.pi, 40)
        ax.plot(c[0] + r * np.cos(th), c[1] + r * np.sin(th), color=color,
                lw=1.6, zorder=3)
        is_ac = source_type == "AC"
        if kind == "VSRC":
            if is_ac:
                ts = np.linspace(-0.18, 0.18, 30)
                wig = (c + np.outer(ts, u)
                      + np.outer(0.10 * np.sin(ts / 0.18 * np.pi), perp))
                ax.plot(wig[:, 0], wig[:, 1], color=color, lw=1.2, zorder=3)
            else:
                # battery bars: long/thin toward node_a ("+"), short/thick
                # toward node_b ("-") -- distinct at a glance from the AC
                # sine symbol above
                p_plus = c - u * 0.09
                p_minus = c + u * 0.09
                for cc, w, lw_ in ((p_plus, 0.17, 1.4), (p_minus, 0.09, 3.4)):
                    e1, e2 = cc + perp * w, cc - perp * w
                    ax.plot([e1[0], e2[0]], [e1[1], e2[1]], color=color,
                            lw=lw_, solid_capstyle="butt", zorder=3)
            lbl_plus = pt(0.08, 0.42)
            lbl_minus = pt(0.92, 0.42)
            ax.text(*lbl_plus, "+", color=color, fontsize=7.5,
                   weight="bold", ha="center", va="center", zorder=5)
            ax.text(*lbl_minus, "−", color=color, fontsize=8.5,
                   weight="bold", ha="center", va="center", zorder=5)
        else:
            tail, head = c - u * 0.16, c + u * 0.16
            ax.annotate("", xy=tuple(head), xytext=tuple(tail),
                       arrowprops=dict(arrowstyle="-|>", color=color,
                                       lw=1.6, mutation_scale=11), zorder=4)
            if is_ac:
                # small sine tick tucked inside the circle, above the arrow
                # -- kept inside the circle radius (not above it) since the
                # component's name label sits right at the circle's top edge
                ts = np.linspace(-0.11, 0.11, 20)
                base = c + perp * 0.15
                wig = (base + np.outer(ts, u)
                      + np.outer(0.045 * np.sin(ts / 0.11 * np.pi), perp))
                ax.plot(wig[:, 0], wig[:, 1], color=color, lw=1.0, zorder=4)
        ax.plot([p_lo[0], (c - u * r)[0]], [p_lo[1], (c - u * r)[1]], **wire)
        ax.plot([(c + u * r)[0], p_hi[0]], [(c + u * r)[1], p_hi[1]], **wire)


def _draw_current_arrow(ax, p1, p2, sign, color="#334155"):
    """Small directional arrow offset off to the side of the p1->p2 edge,
    showing which way current actually flows (sign is the last-sample sign
    of the component's node_a->node_b current). Only meaningful for DC —
    callers should not call this for circuits with any AC source, since the
    direction reverses every half-cycle there."""
    a, b = np.array(p1, dtype=float), np.array(p2, dtype=float)
    d = b - a
    length = np.hypot(*d)
    u = d / length
    perp = np.array([-u[1], u[0]])
    mid = a + u * length * 0.5 + perp * 0.34
    s = 1.0 if sign >= 0 else -1.0
    tail, head = mid - u * 0.14 * s, mid + u * 0.14 * s
    ax.annotate("", xy=tuple(head), xytext=tuple(tail),
               arrowprops=dict(arrowstyle="-|>", color=color, lw=1.5,
                               mutation_scale=10), zorder=5)


class BuilderApp:
    def __init__(self):
        self.placements = []          # dicts: kind,p1,p2,value,name,source_type,freq
        self.wires = []                # (p1, p2) pairs
        self.ground = None             # grid point tuple or None
        self.tool = "R"
        self.anchor = None             # pending first click for 2-click placement
        self.selected = None           # index into self.placements
        self.probes_v = set()          # grid points
        self.probes_i = set()          # indices into self.placements
        self.probes_q = set()          # indices into self.placements (C only)
        self.probes_e = set()          # indices into self.placements (R/L/C)
        self.T_ms = 80.0
        self.result = None
        self.status = ""
        self._next_id = {k: 1 for k in ("R", "L", "C", "VSRC", "ISRC")}

        self.fig = plt.figure(figsize=(14.5, 9.0))
        self.fig.patch.set_facecolor(TH["bg"])
        # Work around a matplotlib TextBox bug: clicking into a TextBox
        # (e.g. the value box) grabs the canvas's mouse; clicking a
        # *different* widget's axes afterwards makes TextBox stop accepting
        # keystrokes but never actually releases that grab, so the next
        # widget's own grab_mouse() raises "Another Axes already grabs
        # mouse input" and the click is lost. Connected here, before any
        # widget is constructed, so it runs first for every click and
        # releases any stale grab before the clicked widget tries to
        # acquire its own.
        self.fig.canvas.mpl_connect("button_press_event",
                                    self._release_stale_grab)
        try:
            self.fig.canvas.manager.set_window_title(
                "Circuit Builder — free-form (Milestone 4)")
        except Exception:
            pass

        self.fig.text(0.02, 0.975, "Circuit Builder", fontsize=15,
                      weight="bold", color=TH["text"], va="center")
        self.fig.text(0.02, 0.952,
                      "components: click 2 points, adjacent incl. diagonal "
                      "· wires: click 2 points on a line, any distance",
                      fontsize=8.5, color=TH["sub"], va="center")

        # ---- palette --------------------------------------------------------
        self.tool_btns = {}
        for j, name in enumerate(TOOLS):
            axb = self.fig.add_axes([0.02 + j * 0.075, 0.905, 0.070, 0.032])
            b = Button(axb, name, color="#e2e8f0", hovercolor="#cbd5e1")
            b.label.set_fontsize(8.5)
            for sp in axb.spines.values():
                sp.set_visible(False)
            b.on_clicked(lambda ev, nm=name: self._set_tool(nm))
            self.tool_btns[name] = b

        # ---- canvas card ------------------------------------------------------
        self._card([0.02, 0.335, 0.575, 0.555], "CANVAS")
        self.ax_canvas = self.fig.add_axes([0.03, 0.345, 0.555, 0.53])
        self.ax_canvas.set_xlim(-0.7, GX - 1 + 0.7)
        self.ax_canvas.set_ylim(-0.7, GY - 1 + 0.7)
        self.ax_canvas.set_aspect("equal")
        self.ax_canvas.axis("off")
        gx, gy = np.meshgrid(range(GX), range(GY))
        self.grid_dots, = self.ax_canvas.plot(gx.ravel(), gy.ravel(), "o",
                                              ms=3, color="#c7d0de",
                                              zorder=1, picker=False)
        self.anchor_dot, = self.ax_canvas.plot([], [], "o", ms=11,
                                                mfc="none",
                                                mec=TH["accent"], mew=2,
                                                zorder=6)

        # ---- properties card ----------------------------------------------------
        self._card([0.605, 0.730, 0.375, 0.160], "PROPERTIES")
        self.txt_sel = self.fig.text(0.615, 0.845, "(nothing selected)",
                                     fontsize=9, weight="bold",
                                     color=TH["text"])
        self.fig.text(0.615, 0.816, "value", fontsize=7.5, color=TH["sub"])
        axv = self.fig.add_axes([0.615, 0.782, 0.10, 0.028])
        self.box_value = TextBox(axv, "", initial="", color="#f8fafc")
        self.box_value.on_submit(self._on_value_submit)

        self.txt_source_lbl = self.fig.text(0.735, 0.816, "source",
                                            fontsize=7.5, color=TH["sub"])
        self.src_btns = {}
        for j, name in enumerate(("AC", "DC")):
            axb = self.fig.add_axes([0.735 + j * 0.058, 0.782, 0.054, 0.028])
            b = Button(axb, name, color="#e2e8f0", hovercolor="#cbd5e1")
            b.label.set_fontsize(8)
            for sp in axb.spines.values():
                sp.set_visible(False)
            b.on_clicked(lambda ev, nm=name: self._on_source_type(nm))
            self.src_btns[name] = b

        self.txt_freq_lbl = self.fig.text(0.860, 0.816, "ω (rad/s)",
                                          fontsize=7.5, color=TH["sub"])
        axf = self.fig.add_axes([0.860, 0.782, 0.10, 0.028])
        self.box_freq = TextBox(axf, "", initial="", color="#f8fafc")
        self.box_freq.on_submit(self._on_freq_submit)

        # probe toggles for non-sources: Q (charge, capacitors only) and E
        # (energy, R/L/C) share the same footprint as the source/ω controls
        # above -- the two groups are mutually exclusive (a component is
        # either a source or a passive element) so only one is ever visible
        self.txt_probe_lbl = self.fig.text(0.735, 0.816, "probe",
                                           fontsize=7.5, color=TH["sub"])
        self.probe_btns = {}
        for j, (key, label) in enumerate((("Q", "Q  (charge)"),
                                          ("E", "E  (energy)"))):
            axb = self.fig.add_axes([0.735 + j * 0.125, 0.782, 0.115, 0.028])
            b = Button(axb, label, color="#e2e8f0", hovercolor="#cbd5e1")
            b.label.set_fontsize(7.5)
            for sp in axb.spines.values():
                sp.set_visible(False)
            b.on_clicked(lambda ev, k=key: self._on_probe_toggle(k))
            self.probe_btns[key] = b

        self.txt_hint = self.fig.text(
            0.615, 0.750, "select a component (SELECT tool) to edit it here",
            fontsize=7.2, color=TH["sub"])
        axflip = self.fig.add_axes([0.615, 0.745, 0.235, 0.026])
        self.btn_flip = Button(axflip, "⇄ flip polarity (+ / −)",
                               color="#e2e8f0", hovercolor="#cbd5e1")
        self.btn_flip.label.set_fontsize(7.3)
        for sp in axflip.spines.values():
            sp.set_visible(False)
        self.btn_flip.on_clicked(lambda ev: self._flip_selected())
        self.btn_flip.ax.set_visible(False)
        self.btn_flip.set_active(False)

        # ---- circuit card (time span, save/load, status) -------------------------
        self._card([0.605, 0.335, 0.375, 0.375], "CIRCUIT")
        axt = self.fig.add_axes([0.680, 0.645, 0.27, 0.020])
        from matplotlib.widgets import Slider
        self.sl_T = Slider(axt, "Time span (ms)", 5.0, 500.0,
                           valinit=self.T_ms, valstep=5.0,
                           color=TH["accent"])
        self.sl_T.label.set_fontsize(8)
        self.sl_T.valtext.set_fontsize(8)
        self.sl_T.on_changed(self._on_T_changed)

        axsave = self.fig.add_axes([0.615, 0.590, 0.115, 0.032])
        self.btn_save = Button(axsave, "Save JSON", color="#e2e8f0",
                               hovercolor="#cbd5e1")
        self.btn_save.label.set_fontsize(8)
        self.btn_save.on_clicked(lambda ev: self.save_json())

        axload = self.fig.add_axes([0.740, 0.590, 0.115, 0.032])
        self.btn_load = Button(axload, "Load JSON", color="#e2e8f0",
                               hovercolor="#cbd5e1")
        self.btn_load.label.set_fontsize(8)
        self.btn_load.on_clicked(lambda ev: self.load_json())

        axclear = self.fig.add_axes([0.865, 0.590, 0.10, 0.032])
        self.btn_clear = Button(axclear, "Clear all", color="#fee2e2",
                                hovercolor="#fecaca")
        self.btn_clear.label.set_fontsize(8)
        self.btn_clear.on_clicked(lambda ev: self.clear_all())

        for axx in (axsave, axload, axclear):
            for sp in axx.spines.values():
                sp.set_visible(False)

        self.txt_status = self.fig.text(0.615, 0.545, "", fontsize=8,
                                        color=TH["sub"], va="top",
                                        linespacing=1.5, wrap=True)
        self.txt_netlist = self.fig.text(0.615, 0.345, "", fontsize=7.3,
                                         family="monospace", color=TH["sub"],
                                         va="bottom", linespacing=1.4)

        # ---- results charts --------------------------------------------------------
        self._card([0.02, 0.02, 0.955, 0.29], "RESULTS  —  probed nodes "
                   "(voltage), components (current, charge, energy)")
        self.ax_v = self.fig.add_axes([0.045, 0.045, 0.27, 0.235])
        self.ax_i = self.fig.add_axes([0.363, 0.045, 0.27, 0.235])
        self.ax_q = self.fig.add_axes([0.681, 0.045, 0.225, 0.235])
        self.ax_e = self.ax_q.twinx()          # shares x, own y-axis (Joules)
        for ax, lab in ((self.ax_v, "Voltage (V)"), (self.ax_i, "Current (A)"),
                        (self.ax_q, "Charge (C)")):
            ax.set_facecolor("white")
            ax.set_xlabel("t (ms)", fontsize=8)
            ax.set_ylabel(lab, fontsize=8)
            ax.tick_params(labelsize=7.5, colors=TH["sub"])
            ax.grid(True, ls="--", lw=0.5, alpha=0.4)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
        self.ax_e.set_ylabel("Energy (J)", fontsize=8)
        self.ax_e.tick_params(labelsize=7.5, colors=TH["sub"])
        self.ax_e.spines["top"].set_visible(False)

        self.fig.canvas.mpl_connect("button_press_event", self._on_click)
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)
        self._set_tool("R")
        self._show_properties(None)
        self._resolve()
        self._set_status("Pick a tool, place components on the grid, set "
                         "a ground, then add a source.")

    # ------------------------------------------------------------------ card
    def _card(self, rect, title=None):
        x, y, w, h = rect
        box = FancyBboxPatch((x, y), w, h, transform=self.fig.transFigure,
                             boxstyle="round,pad=0,rounding_size=0.008",
                             mutation_aspect=1.62, fc=TH["card"],
                             ec=TH["edge"], lw=1.1, zorder=-1)
        self.fig.add_artist(box)
        if title:
            self.fig.text(x + 0.011, y + h - 0.016, title, fontsize=7.8,
                          weight="bold", color=TH["sub"], va="top")

    # -------------------------------------------------------------- palette
    def _set_tool(self, name):
        self.tool = name
        self.anchor = None
        self.anchor_dot.set_data([], [])
        for nm, b in self.tool_btns.items():
            active = (nm == name)
            b.color = TH["accent"] if active else "#e2e8f0"
            b.hovercolor = "#1e40af" if active else "#cbd5e1"
            b.ax.set_facecolor(b.color)
            b.label.set_color("white" if active else TH["text"])
            b.label.set_fontweight("bold" if active else "normal")
        self.fig.canvas.draw_idle()

    def _set_status(self, msg):
        self.status = msg
        self.txt_status.set_text(msg)
        self.fig.canvas.draw_idle()

    # ------------------------------------------------------------- geometry
    def _snap_point(self, x, y, tol=0.42):
        col, row = round(x), round(y)
        if 0 <= col < GX and 0 <= row < GY and \
                math.hypot(x - col, y - row) <= tol:
            return (col, row)
        return None

    def _edge_key(self, p1, p2):
        return tuple(sorted((p1, p2)))

    def _edge_occupied(self, p1, p2):
        key = self._edge_key(p1, p2)
        for pl in self.placements:
            if self._edge_key(pl["p1"], pl["p2"]) == key:
                return True
        for w in self.wires:
            if self._edge_key(*w) == key:
                return True
        return False

    def _hit_component(self, x, y, tol=0.32):
        best, best_d = None, tol
        for i, pl in enumerate(self.placements):
            mx = (pl["p1"][0] + pl["p2"][0]) / 2.0
            my = (pl["p1"][1] + pl["p2"][1]) / 2.0
            d = math.hypot(x - mx, y - my)
            if d < best_d:
                best, best_d = i, d
        return best

    def _hit_wire(self, x, y, tol=0.28):
        best, best_d = None, tol
        for i, (p1, p2) in enumerate(self.wires):
            mx, my = (p1[0] + p2[0]) / 2.0, (p1[1] + p2[1]) / 2.0
            d = math.hypot(x - mx, y - my)
            if d < best_d:
                best, best_d = i, d
        return best

    def _release_stale_grab(self, ev):
        # NB: ev.inaxes is not trustworthy here -- while a grab is active,
        # matplotlib's own LocationEvent.__init__ sets .inaxes to the
        # grabber itself regardless of where the click actually landed
        # (that's the whole point of a mouse grab). Check the actual pixel
        # position against the grabber's bounding box instead, the same
        # way TextBox._click itself decides whether a click was "outside".
        grabber = self.fig.canvas.mouse_grabber
        if grabber is not None and not grabber.contains(ev)[0]:
            self.fig.canvas.release_mouse(grabber)

    # -------------------------------------------------------------- clicks
    def _on_click(self, ev):
        if ev.inaxes is not self.ax_canvas or ev.xdata is None:
            return
        x, y = ev.xdata, ev.ydata

        if self.tool == "GROUND":
            p = self._snap_point(x, y)
            if p is not None:
                self.ground = p
                self._set_status(f"Ground set at {p}.")
                self._resolve()
            return

        if self.tool == "DELETE":
            ci = self._hit_component(x, y)
            if ci is not None:
                self._delete_component(ci)
                return
            wi = self._hit_wire(x, y)
            if wi is not None:
                del self.wires[wi]
                self._set_status("Wire removed.")
                self._resolve()
            return

        if self.tool == "SELECT":
            ci = self._hit_component(x, y)
            if ci is not None:
                self._select(ci)
                return
            p = self._snap_point(x, y)
            if p is not None:
                if p in self.probes_v:
                    self.probes_v.discard(p)
                    self._set_status(f"Voltage probe removed at {p}.")
                else:
                    self.probes_v.add(p)
                    self._set_status(f"Probing voltage at {p}.")
                self._redraw_canvas()
                self._redraw_charts()
            return

        # placement tools: R, L, C, VSRC, ISRC, WIRE
        p = self._snap_point(x, y)
        if p is None:
            return
        if self.anchor is None:
            self.anchor = p
            self.anchor_dot.set_data([p[0]], [p[1]])
            self.fig.canvas.draw_idle()
            return
        if p == self.anchor:
            self.anchor = None
            self.anchor_dot.set_data([], [])
            self.fig.canvas.draw_idle()
            return

        if self.tool == "WIRE":
            path = self._wire_path(self.anchor, p)
            if path is None:
                self.anchor = p            # not a straight line -> move anchor
                self.anchor_dot.set_data([p[0]], [p[1]])
                self.fig.canvas.draw_idle()
                return
            self._place_wire_path(path)
            self.anchor = None
            self.anchor_dot.set_data([], [])
            return

        dx = abs(p[0] - self.anchor[0])
        dy = abs(p[1] - self.anchor[1])
        if max(dx, dy) != 1:
            self.anchor = p                    # not adjacent -> move anchor
            self.anchor_dot.set_data([p[0]], [p[1]])
            self.fig.canvas.draw_idle()
            return
        self._place(self.anchor, p)
        self.anchor = None
        self.anchor_dot.set_data([], [])

    def _on_key(self, ev):
        if ev.key == "escape":
            self.anchor = None
            self.anchor_dot.set_data([], [])
            self.fig.canvas.draw_idle()

    # ---------------------------------------------------------- placement
    def _wire_path(self, a, b):
        """Grid points from a to b inclusive, straight-line only (0/45/90
        degrees). Returns None if a->b isn't a straight 8-direction line —
        the caller then treats it as a normal "move the anchor" click. This
        is what lets a long wire run be placed with two clicks instead of
        one click per grid hop."""
        dx, dy = b[0] - a[0], b[1] - a[1]
        if dx == 0 and dy == 0:
            return None
        if not (dx == 0 or dy == 0 or abs(dx) == abs(dy)):
            return None
        steps = max(abs(dx), abs(dy))
        sx = 0 if dx == 0 else dx // abs(dx)
        sy = 0 if dy == 0 else dy // abs(dy)
        return [(a[0] + sx * k, a[1] + sy * k) for k in range(steps + 1)]

    def _place_wire_path(self, path):
        placed = skipped = 0
        for p1, p2 in zip(path, path[1:]):
            if self._edge_occupied(p1, p2):
                skipped += 1
                continue
            self.wires.append((p1, p2))
            placed += 1
        if skipped:
            self._set_status(f"Wire {path[0]}-{path[-1]}: {placed} "
                             f"segment(s) placed, {skipped} skipped "
                             f"(already occupied).")
        else:
            self._set_status(f"Wire placed {path[0]}-{path[-1]} "
                             f"({placed} segment(s)).")
        self._resolve()

    def _place(self, p1, p2):
        if self._edge_occupied(p1, p2):
            self._set_status("That edge already has something on it — "
                             "delete it first.")
            return
        kind = self.tool
        n = self._next_id[kind]
        self._next_id[kind] += 1
        name = f"{kind}{n}"
        default_val = {"R": 1000.0, "L": 1.0, "C": 1e-6,
                       "VSRC": 120.0, "ISRC": 0.1}[kind]
        pl = dict(kind=kind, p1=p1, p2=p2, value=default_val, name=name,
                 source_type="AC", freq=377.0)
        self.placements.append(pl)
        self._set_status(f"Placed {name} between {p1} and {p2}.")
        self._select_for_edit(len(self.placements) - 1)
        self._resolve()

    def _delete_component(self, idx):
        pl = self.placements.pop(idx)
        self.probes_i = {j if j < idx else j - 1 for j in self.probes_i
                         if j != idx}
        self.probes_q = {j if j < idx else j - 1 for j in self.probes_q
                         if j != idx}
        self.probes_e = {j if j < idx else j - 1 for j in self.probes_e
                         if j != idx}
        if self.selected == idx:
            self.selected = None
            self._show_properties(None)
        elif self.selected is not None and self.selected > idx:
            self.selected -= 1
        self._set_status(f"Deleted {pl['name']}.")
        self._resolve()

    def clear_all(self):
        self.placements.clear()
        self.wires.clear()
        self.ground = None
        self.probes_v.clear()
        self.probes_i.clear()
        self.probes_q.clear()
        self.probes_e.clear()
        self.selected = None
        self.anchor = None
        self._next_id = {k: 1 for k in ("R", "L", "C", "VSRC", "ISRC")}
        self._show_properties(None)
        self._set_status("Cleared.")
        self._resolve()

    # ----------------------------------------------------------- selection
    def _select_for_edit(self, idx):
        """Select a component for value/source editing only — does NOT
        touch its probe status. Used right after placing a new component
        (so its value is immediately editable) without silently probing
        every component you place."""
        self.selected = idx
        self._show_properties(self.placements[idx])
        self._redraw_canvas()

    def _select(self, idx):
        """SELECT-tool click on a component: select it for editing *and*
        toggle it as a current probe — the two are combined here because
        that's what a deliberate click on an already-placed component
        means; placement itself uses `_select_for_edit` instead."""
        self._select_for_edit(idx)
        if idx in self.probes_i:
            self.probes_i.discard(idx)
        else:
            self.probes_i.add(idx)
        self._redraw_canvas()
        self._redraw_charts()

    def _show_widget(self, widget, show):
        """Toggle a Button/TextBox's visibility AND activate/deactivate it
        together. Several widgets in PROPERTIES occupy the *same* screen
        rectangle (source/ω controls vs. Q/E probe controls are mutually
        exclusive, so their axes overlap to save space) and are switched by
        visibility alone. `ax.set_visible(False)` only stops it from being
        *drawn* -- matplotlib's own hit-testing (`Widget.ignore()`) only
        checks `.active`, not visibility, so a hidden-but-still-active
        Button sitting under a visible one still receives every click on
        that shared rectangle and tries to grab the mouse for itself too,
        colliding with the visible widget's own grab and raising
        'Another Axes already grabs mouse input'. Always use this instead
        of touching `.ax.set_visible()` directly for any toggled widget."""
        widget.ax.set_visible(show)
        widget.set_active(show)

    def _show_properties(self, pl):
        if pl is None:
            self.txt_sel.set_text("(nothing selected)")
            self.box_value.set_val("")
            self.box_freq.set_val("")
            self._show_widget(self.box_freq, False)
            self._show_widget(self.btn_flip, False)
            self.txt_hint.set_visible(True)
            self.txt_source_lbl.set_visible(False)
            self.txt_freq_lbl.set_visible(False)
            self.txt_probe_lbl.set_visible(False)
            for b in self.src_btns.values():
                b.color = "#e2e8f0"
                b.ax.set_facecolor(b.color)
                self._show_widget(b, False)
            for b in self.probe_btns.values():
                self._show_widget(b, False)
            self.fig.canvas.draw_idle()
            return
        is_src = pl["kind"] in ("VSRC", "ISRC")
        label = "+" if pl["kind"] == "VSRC" else "current from"
        self.txt_sel.set_text(f"{pl['name']}  ({pl['kind']})" +
                              (f"   [{label} @ {pl['p1']}]" if is_src
                               else ""))
        self.box_value.set_val(f"{pl['value']:.6g}")
        show_freq = is_src and pl["source_type"] == "AC"
        self._show_widget(self.box_freq, show_freq)
        self.txt_freq_lbl.set_visible(show_freq)
        self.box_freq.set_val(f"{pl['freq']:.6g}")
        self._show_widget(self.btn_flip, is_src)
        self.txt_hint.set_visible(not is_src)
        self.txt_source_lbl.set_visible(is_src)
        for nm, b in self.src_btns.items():
            self._show_widget(b, is_src)
            active = is_src and nm == pl["source_type"]
            b.color = TH["accent"] if active else "#e2e8f0"
            b.ax.set_facecolor(b.color)
            b.label.set_color("white" if active else TH["text"])

        can_q = pl["kind"] == "C"                       # probe_charge: C only
        can_e = pl["kind"] in ("R", "L", "C")            # probe_energy: not sources
        self.txt_probe_lbl.set_visible(not is_src)
        idx = self.selected
        for key, can in (("Q", can_q), ("E", can_e)):
            b = self.probe_btns[key]
            self._show_widget(b, not is_src and can)
            active = can and idx is not None and \
                idx in (self.probes_q if key == "Q" else self.probes_e)
            b.color = TH["accent"] if active else "#e2e8f0"
            b.ax.set_facecolor(b.color)
            b.label.set_color("white" if active else TH["text"])
        self.fig.canvas.draw_idle()

    def _on_probe_toggle(self, key):
        """Toggle the selected component into/out of the charge (Q, C only)
        or energy (E, R/L/C) probe set -- kept as separate buttons rather
        than folded into the SELECT-tool click (which already toggles the
        current probe) since a component can be probed for any combination
        of I/Q/E at once and a single click can't express that."""
        if self.selected is None:
            return
        pl = self.placements[self.selected]
        target = self.probes_q if key == "Q" else self.probes_e
        allowed = pl["kind"] == "C" if key == "Q" else \
            pl["kind"] in ("R", "L", "C")
        if not allowed:
            return
        kind_name = "charge" if key == "Q" else "energy"
        if self.selected in target:
            target.discard(self.selected)
            self._set_status(f"{kind_name.capitalize()} probe removed "
                             f"from {pl['name']}.")
        else:
            target.add(self.selected)
            self._set_status(f"Probing {kind_name} of {pl['name']}.")
        self._show_properties(pl)
        self._redraw_canvas()
        self._redraw_charts()

    def _flip_selected(self):
        """Swap node_a/node_b on the selected source, reversing its
        polarity (VSRC "+"/"-") or current direction (ISRC) without having
        to delete and re-place it in the opposite click order."""
        if self.selected is None:
            return
        pl = self.placements[self.selected]
        if pl["kind"] not in ("VSRC", "ISRC"):
            return
        pl["p1"], pl["p2"] = pl["p2"], pl["p1"]
        self._set_status(f"Flipped {pl['name']} — now referenced from "
                         f"{pl['p1']}.")
        self._show_properties(pl)
        self._resolve()

    def _on_value_submit(self, text):
        if self.selected is None:
            return
        try:
            v = float(text.strip().replace(",", "."))
        except ValueError:
            self._show_properties(self.placements[self.selected])
            return
        if v <= 0:
            self._show_properties(self.placements[self.selected])
            return
        self.placements[self.selected]["value"] = v
        self._resolve()

    def _on_freq_submit(self, text):
        if self.selected is None:
            return
        try:
            v = float(text.strip().replace(",", "."))
        except ValueError:
            self._show_properties(self.placements[self.selected])
            return
        self.placements[self.selected]["freq"] = max(v, 0.0)
        self._resolve()

    def _on_source_type(self, name):
        if self.selected is None or \
                self.placements[self.selected]["kind"] not in \
                ("VSRC", "ISRC"):
            return
        self.placements[self.selected]["source_type"] = name
        self._show_properties(self.placements[self.selected])
        self._resolve()

    def _on_T_changed(self, val):
        self.T_ms = float(val)
        self._resolve()

    # --------------------------------------------------------------- netlist
    def _canonical_names(self):
        uf = UnionFind()
        for p1, p2 in self.wires:
            uf.union(p1, p2)
        ground_root = uf.find(self.ground) if self.ground is not None else None

        def name_of(p):
            root = uf.find(p)
            if ground_root is not None and root == ground_root:
                return "0"
            return f"n{root[0]}_{root[1]}"
        return name_of

    def build_netlist(self):
        """Raises ValueError (with a plain-English message) if a wire
        connects both terminals of the same component to one node -- a
        short circuit. This is checked here, before constructing any
        Component, because Component.__post_init__ also rejects
        node_a == node_b but with a message aimed at library callers, not
        UI users, and only for the first offender found."""
        name_of = self._canonical_names()
        shorted = [pl["name"] for pl in self.placements
                  if name_of(pl["p1"]) == name_of(pl["p2"])]
        if shorted:
            raise ValueError(
                f"{', '.join(shorted)} — a wire connects both of that "
                f"component's terminals to the same node. Remove the "
                f"shorting wire, or move the component.")
        nl = Netlist(ground="0")
        for pl in self.placements:
            nl.add(pl["kind"], name_of(pl["p1"]), name_of(pl["p2"]),
                  pl["value"], name=pl["name"],
                  source_type=pl.get("source_type", "DC"),
                  freq=pl.get("freq", 0.0))
        return nl, name_of

    def _resolve(self):
        self.result = None
        try:
            nl, name_of = self.build_netlist()
        except ValueError as e:
            self._name_of = self._canonical_names()
            self.txt_netlist.set_text("(short circuit -- see status)")
            self._set_status(f"Short circuit: {e}")
            self._redraw_canvas()
            self._redraw_charts()
            return
        self._name_of = name_of
        self.txt_netlist.set_text(self._netlist_summary(nl))

        if self.ground is None:
            self._set_status("Add a ground reference (GROUND tool) to "
                             "simulate.")
            self._redraw_canvas()
            self._redraw_charts()
            return
        if not any(c.is_source() for c in nl.components):
            self._set_status("Add a voltage or current source to "
                             "simulate.")
            self._redraw_canvas()
            self._redraw_charts()
            return
        try:
            nl.validate()
        except ValueError as e:
            self._set_status(f"Not solvable yet: {e}")
            self._redraw_canvas()
            self._redraw_charts()
            return

        t = np.linspace(0.0, self.T_ms * 1e-3, 2000)
        try:
            self.result = simulate(nl, t)
            self._set_status(f"Solved: {len(nl.components)} components, "
                             f"{len(nl.non_ground_nodes)} nodes.")
        except Exception as e:                        # noqa: BLE001
            self._set_status(f"Solve failed: {e}")
        self._redraw_canvas()
        self._redraw_charts()

    def _netlist_summary(self, nl):
        if not nl.components:
            return "(empty circuit)"
        lines = [f"{c.name:6s} {c.node_a:>6s} - {c.node_b:<6s} "
                f"{c.value:.4g}" + (f"  {c.source_type}" if c.is_source()
                                    else "")
                for c in nl.components]
        return "\n".join(lines[-10:])

    # --------------------------------------------------------------- drawing
    def _redraw_canvas(self):
        ax = self.ax_canvas
        for art in list(ax.patches) + [a for a in ax.lines
                                       if a not in (self.grid_dots,
                                                    self.anchor_dot)] + \
                ax.texts + ax.collections:
            art.remove()

        sources = [pl for pl in self.placements
                  if pl["kind"] in ("VSRC", "ISRC")]
        show_flow = (self.result is not None and sources
                    and all(pl["source_type"] == "DC" for pl in sources))

        for p1, p2 in self.wires:
            _draw_component(ax, "WIRE", p1, p2, "#334155")
        for i, pl in enumerate(self.placements):
            color = KIND_COLOR[pl["kind"]]
            _draw_component(ax, pl["kind"], pl["p1"], pl["p2"], color,
                            source_type=pl.get("source_type", "DC"))
            mx = (pl["p1"][0] + pl["p2"][0]) / 2.0
            my = (pl["p1"][1] + pl["p2"][1]) / 2.0
            ax.text(mx, my + 0.30, pl["name"], fontsize=6.5, ha="center",
                   color=color, zorder=5)
            if show_flow:
                try:
                    i_arr = probe_current(self.result, pl["name"])
                    if len(i_arr):
                        _draw_current_arrow(ax, pl["p1"], pl["p2"],
                                           float(i_arr[-1]))
                except KeyError:
                    pass
            if i == self.selected:
                ax.plot(mx, my, "s", ms=14, mfc="none", mec=TH["accent"],
                       mew=1.6, zorder=4)
            if i in self.probes_i:
                pc = PROBE_COLORS[list(self.probes_i).index(i)
                                  % len(PROBE_COLORS)] \
                    if False else self._probe_color_i(i)
                ax.plot(mx, my, "o", ms=16, mfc="none", mec=pc, mew=1.8,
                       zorder=4)
            if i in self.probes_q:
                ax.plot(mx, my, "D", ms=9, mfc="none",
                       mec=self._probe_color_q(i), mew=1.3, zorder=4)
            if i in self.probes_e:
                ax.plot(mx, my, "^", ms=11, mfc="none",
                       mec=self._probe_color_e(i), mew=1.3, zorder=4)

        if self.ground is not None:
            gx, gy = self.ground
            ax.plot([gx, gx], [gy, gy - 0.22], color="#334155", lw=1.6,
                   zorder=3)
            for k, wgt in enumerate((0.20, 0.13, 0.06)):
                yy = gy - 0.22 - k * 0.07
                ax.plot([gx - wgt, gx + wgt], [yy, yy], color="#334155",
                       lw=1.4, zorder=3)

        for p in self.probes_v:
            pc = self._probe_color_v(p)
            ax.plot(p[0], p[1], "o", ms=12, mfc="none", mec=pc, mew=1.8,
                   zorder=4)

        self.fig.canvas.draw_idle()

    def _probe_color_v(self, p):
        ordered = sorted(self.probes_v)
        return PROBE_COLORS[ordered.index(p) % len(PROBE_COLORS)]

    def _probe_color_i(self, idx):
        ordered = sorted(self.probes_i)
        return PROBE_COLORS[ordered.index(idx) % len(PROBE_COLORS)]

    def _probe_color_q(self, idx):
        ordered = sorted(self.probes_q)
        return PROBE_COLORS[ordered.index(idx) % len(PROBE_COLORS)]

    def _probe_color_e(self, idx):
        ordered = sorted(self.probes_e)
        return PROBE_COLORS[ordered.index(idx) % len(PROBE_COLORS)]

    def _redraw_charts(self):
        for ax, lab in ((self.ax_v, "Voltage (V)"), (self.ax_i, "Current (A)"),
                        (self.ax_q, "Charge (C)")):
            ax.cla()
            ax.set_facecolor("white")
            ax.set_xlabel("t (ms)", fontsize=8)
            ax.set_ylabel(lab, fontsize=8)
            ax.tick_params(labelsize=7.5, colors=TH["sub"])
            ax.grid(True, ls="--", lw=0.5, alpha=0.4)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
        self.ax_e.cla()
        # cla() resets a twin axes back to default (left-side) tick/label
        # placement, undoing what twinx() set up at construction -- have to
        # re-assert "this axis lives on the right" on every redraw, not
        # just once at __init__ time.
        self.ax_e.yaxis.tick_right()
        self.ax_e.yaxis.set_label_position("right")
        self.ax_e.set_facecolor("none")
        self.ax_e.set_ylabel("Energy (J)", fontsize=8)
        self.ax_e.tick_params(labelsize=7.5, colors=TH["sub"])
        self.ax_e.spines["top"].set_visible(False)

        if self.result is None:
            for ax in (self.ax_v, self.ax_i, self.ax_q):
                ax.text(0.5, 0.5, "no solved result yet",
                       transform=ax.transAxes, ha="center",
                       color=TH["sub"], fontsize=9)
            self.fig.canvas.draw_idle()
            return

        tms = self.result.t * 1e3
        for p in sorted(self.probes_v):
            name = self._name_of(p)
            try:
                v = probe_voltage(self.result, name)
            except KeyError:
                continue
            self.ax_v.plot(tms, v, color=self._probe_color_v(p), lw=1.5,
                           label=f"V{p}")
        for i in sorted(self.probes_i):
            if i >= len(self.placements):
                continue
            pl = self.placements[i]
            try:
                iv = probe_current(self.result, pl["name"])
            except KeyError:
                continue
            self.ax_i.plot(tms, iv, color=self._probe_color_i(i), lw=1.5,
                           label=f"I({pl['name']})")
        q_lines, e_lines = [], []
        for i in sorted(self.probes_q):
            if i >= len(self.placements) or self.placements[i]["kind"] != "C":
                continue
            pl = self.placements[i]
            try:
                qv = probe_charge(self.result, pl["name"])
            except KeyError:
                continue
            ln, = self.ax_q.plot(tms, qv, color=self._probe_color_q(i),
                                 lw=1.5, ls="-", label=f"Q({pl['name']})")
            q_lines.append(ln)
        for i in sorted(self.probes_e):
            if i >= len(self.placements) or \
                    self.placements[i]["kind"] not in ("R", "L", "C"):
                continue
            pl = self.placements[i]
            try:
                ev = probe_energy(self.result, pl["name"])
            except (KeyError, ValueError):
                continue
            ln, = self.ax_e.plot(tms, ev, color=self._probe_color_e(i),
                                 lw=1.5, ls="--", label=f"E({pl['name']})")
            e_lines.append(ln)
        if self.probes_v:
            self.ax_v.legend(fontsize=7, frameon=False)
        if self.probes_i:
            self.ax_i.legend(fontsize=7, frameon=False)
        if q_lines or e_lines:
            lines = q_lines + e_lines
            self.ax_q.legend(lines, [ln.get_label() for ln in lines],
                            fontsize=6, frameon=True, framealpha=0.8,
                            facecolor="white", edgecolor="none",
                            loc="upper left", labelspacing=0.25,
                            handlelength=1.4, borderpad=0.3)
        self.fig.canvas.draw_idle()

    # ------------------------------------------------------------- save/load
    def save_json(self):
        # placements/wires/ground below are the actual round-trip source of
        # truth (see load_json) -- `netlist` is just a convenience snapshot
        # for anyone inspecting the file by hand, so saving must still work
        # even mid-short-circuit (build_netlist() raises ValueError then)
        try:
            nl, _ = self.build_netlist()
            netlist_dict = nl.to_dict()
        except ValueError:
            netlist_dict = None
        data = dict(netlist=netlist_dict,
                   placements=[{**pl, "p1": list(pl["p1"]),
                               "p2": list(pl["p2"])}
                              for pl in self.placements],
                   wires=[[list(p1), list(p2)] for p1, p2 in self.wires],
                   ground=list(self.ground) if self.ground else None,
                   T_ms=self.T_ms)
        path = self._file_dialog(save=True)
        if path is None:
            return
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
        self._set_status(f"Saved to {path}")

    def load_json(self):
        path = self._file_dialog(save=False)
        if path is None:
            return
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        self.clear_all()
        for pl in data["placements"]:
            pl = dict(pl)
            pl["p1"], pl["p2"] = tuple(pl["p1"]), tuple(pl["p2"])
            self.placements.append(pl)
        self.wires = [(tuple(p1), tuple(p2)) for p1, p2 in data["wires"]]
        self.ground = tuple(data["ground"]) if data.get("ground") else None
        self.T_ms = data.get("T_ms", 80.0)
        self.sl_T.eventson = False
        self.sl_T.set_val(self.T_ms)
        self.sl_T.eventson = True
        for k in self._next_id:
            used = [int(pl["name"][len(k):]) for pl in self.placements
                   if pl["kind"] == k and pl["name"][len(k):].isdigit()]
            self._next_id[k] = (max(used) + 1) if used else 1
        self._set_status(f"Loaded from {path}")
        self._resolve()

    def _file_dialog(self, save):
        try:
            import tkinter as tk
            from tkinter import filedialog
            root = tk.Tk()
            root.withdraw()
            if save:
                path = filedialog.asksaveasfilename(
                    defaultextension=".json",
                    filetypes=[("Circuit JSON", "*.json")])
            else:
                path = filedialog.askopenfilename(
                    filetypes=[("Circuit JSON", "*.json")])
            root.destroy()
            return path or None
        except Exception:
            path = "circuit.json"
            self._set_status(f"(no file dialog available — using "
                             f"./{path})")
            return path


class _FakeEvent:
    """A minimal stand-in for a matplotlib MouseEvent, just enough for
    BuilderApp._on_click to read .inaxes/.xdata/.ydata — used so --test can
    exercise the *real* click-dispatch code path (tool routing, grid
    snapping, adjacency checks, hit-testing) instead of calling internal
    methods directly and missing bugs those methods' callers would catch."""

    def __init__(self, ax, x, y):
        self.inaxes = ax
        self.xdata = x
        self.ydata = y


def _click(app, x, y):
    app._on_click(_FakeEvent(app.ax_canvas, x, y))


def _self_test(app, out):
    import os
    from rlc_solver import solve

    base, ext = os.path.splitext(out)

    # 1) build a series RLC entirely through the real click dispatcher and
    #    verify it matches the exact closed-form solution
    app._set_tool("VSRC")
    _click(app, 0, 1)
    _click(app, 0, 0)
    assert app.placements[-1]["kind"] == "VSRC", app.placements
    app.placements[-1].update(value=120.0, freq=377.0, source_type="AC")

    app._set_tool("R")
    _click(app, 0, 1)
    _click(app, 1, 1)
    app.placements[-1]["value"] = 1000.0

    app._set_tool("L")
    _click(app, 1, 1)
    _click(app, 2, 1)
    app.placements[-1]["value"] = 3.5

    app._set_tool("C")
    _click(app, 2, 1)
    _click(app, 2, 0)
    app.placements[-1]["value"] = 2e-6

    app._set_tool("WIRE")
    _click(app, 2, 0)
    _click(app, 1, 0)
    _click(app, 1, 0)
    _click(app, 0, 0)

    app._set_tool("GROUND")
    _click(app, 0, 0)

    assert app.result is not None, app.status
    nl, _ = app.build_netlist()
    assert len(nl.components) == 4
    t = app.result.t
    sol = solve("RLC", 1000.0, 3.5, 2.0, 120.0, 377.0, t, mode="AC")
    I_r1 = probe_current(app.result, "R1")
    Q_c1 = probe_charge(app.result, "C1")
    rI = float(np.max(np.abs(I_r1 - sol["I"]))) / \
        max(float(np.max(np.abs(sol["I"]))), 1e-15)
    rQ = float(np.max(np.abs(Q_c1 - sol["Q"]))) / \
        max(float(np.max(np.abs(sol["Q"]))), 1e-15)
    print(f"builder : click-built series RLC vs exact solver — "
          f"I relerr={rI:.1e}  Q relerr={rQ:.1e}")
    assert rI < 5e-3 and rQ < 5e-3, (rI, rQ)

    # 2) probing: select R1 for a current probe, a node for a voltage probe
    app._set_tool("SELECT")
    r_idx = next(i for i, pl in enumerate(app.placements)
                if pl["kind"] == "R")
    app._select(r_idx)
    assert r_idx in app.probes_i
    app.probes_v.add((1, 1))
    app._redraw_charts()
    print("builder : probe toggling OK (placement no longer auto-probes)")

    # 2b) charge/energy probing: Q is C-only, E is R/L/C (not sources); both
    #     toggle independently of the I probe and of each other. Toggling
    #     these doesn't touch app.result, so it's safe to leave them set for
    #     the rest of the suite -- later tests probe rlc_mna directly.
    c_idx = next(i for i, pl in enumerate(app.placements)
                if pl["kind"] == "C")
    v_idx0 = next(i for i, pl in enumerate(app.placements)
                 if pl["kind"] == "VSRC")

    app._select_for_edit(c_idx)
    app._on_probe_toggle("Q")
    app._on_probe_toggle("E")
    assert c_idx in app.probes_q and c_idx in app.probes_e

    app._select_for_edit(r_idx)
    app._on_probe_toggle("Q")                 # R is not a capacitor -> no-op
    assert r_idx not in app.probes_q
    app._on_probe_toggle("E")
    assert r_idx in app.probes_e

    app._select_for_edit(v_idx0)
    app._on_probe_toggle("E")                 # sources have no energy probe
    assert v_idx0 not in app.probes_e

    app._redraw_canvas()                      # Q/E canvas badges must not raise
    app._redraw_charts()
    qv = probe_charge(app.result, app.placements[c_idx]["name"])
    ev = probe_energy(app.result, app.placements[r_idx]["name"])
    assert len(qv) == len(app.result.t) and len(ev) == len(app.result.t)
    print("builder : charge/energy probe toggling (Q=C-only, E=not-sources) OK")

    # 3) source editing: AC/DC toggle + frequency
    v_idx = next(i for i, pl in enumerate(app.placements)
                if pl["kind"] == "VSRC")
    app._select_for_edit(v_idx)
    assert app.box_freq.ax.get_visible()
    app._on_source_type("DC")
    assert app.placements[v_idx]["source_type"] == "DC"
    app._on_freq_submit("500")
    assert app.placements[v_idx]["freq"] == 500.0
    app._on_source_type("AC")           # restore for the screenshot
    app._on_freq_submit("377")
    print("builder : source type/frequency editing OK")

    app.fig.savefig(base + "_rlc" + ext, dpi=110)

    # 4) flip polarity: reversing a source's node_a/node_b must negate every
    #    downstream result (I and Q both flip sign) -- this is the concrete
    #    check that polarity is both physically wired up and user-
    #    controllable, not just silently determined by click order
    pre_I = probe_current(app.result, "R1").copy()
    pre_Q = probe_charge(app.result, "C1").copy()
    app._set_tool("SELECT")
    app._select_for_edit(v_idx)
    app._flip_selected()
    post_I = probe_current(app.result, "R1")
    post_Q = probe_charge(app.result, "C1")
    assert np.allclose(post_I, -pre_I, atol=1e-9), \
        float(np.max(np.abs(post_I + pre_I)))
    assert np.allclose(post_Q, -pre_Q, atol=1e-9), \
        float(np.max(np.abs(post_Q + pre_Q)))
    app._flip_selected()                # flip back for the tests below
    print("builder : flip polarity negates I/Q as expected")

    # 5) diagonal placement: a component can span a 45-degree grid edge, not
    #    just horizontal/vertical (needed for e.g. delta-star wiring)
    app4 = BuilderApp()
    app4._set_tool("R")
    _click(app4, 4, 4)
    _click(app4, 5, 5)
    assert app4.placements and app4.placements[-1]["p1"] == (4, 4) and \
        app4.placements[-1]["p2"] == (5, 5), app4.placements
    print("builder : diagonal component placement OK")

    # 6) multi-hop wire: two clicks on a straight (horizontal/vertical/
    #    diagonal) line auto-fill every intermediate grid point instead of
    #    needing one click per hop
    app4._set_tool("WIRE")
    _click(app4, 0, 8)
    _click(app4, 4, 8)                  # 4 grid steps away, same row
    assert len(app4.wires) == 4, app4.wires
    assert app4.wires[0] == ((0, 8), (1, 8))
    assert app4.wires[-1] == ((3, 8), (4, 8))
    app4._set_tool("WIRE")
    _click(app4, 6, 0)
    _click(app4, 9, 3)                  # diagonal, 3 grid steps away
    assert len(app4.wires) == 4 + 3, app4.wires
    print("builder : multi-hop wire auto-fill (straight + diagonal) OK")

    # 7) DC current-flow arrows render without raising; the redraw path
    #    that computes them (show_flow) only engages when every source in
    #    the circuit is DC
    app5 = BuilderApp()
    app5._set_tool("VSRC")
    _click(app5, 0, 0)
    _click(app5, 0, 1)
    app5.placements[-1].update(value=10.0, source_type="DC")
    app5._set_tool("R")
    _click(app5, 0, 1)
    _click(app5, 1, 1)
    app5._set_tool("WIRE")
    _click(app5, 1, 1)
    _click(app5, 1, 0)
    _click(app5, 1, 0)
    _click(app5, 0, 0)
    app5._set_tool("GROUND")
    _click(app5, 0, 0)
    assert app5.result is not None, app5.status
    app5._redraw_canvas()               # must not raise for an all-DC circuit
    print("builder : DC current-flow arrows render without error")

    # 8) delete via the real click dispatcher
    n_before = len(app.placements)
    app._set_tool("DELETE")
    _click(app, 1.5, 1.0)               # L1's midpoint
    assert len(app.placements) == n_before - 1
    print("builder : delete via click OK")

    # 9) save/load JSON round-trip (bypass the tkinter file dialog)
    path = base + "_saved.json"
    data_before = app.build_netlist()[0].to_dict()
    app._file_dialog = lambda save: path
    app.save_json()
    app2 = BuilderApp()
    app2._file_dialog = lambda save: path
    app2.load_json()
    data_after = app2.build_netlist()[0].to_dict()
    assert data_before == data_after
    print("builder : save/load JSON round-trip OK")

    # 10) validation feedback for incomplete circuits never crashes
    app3 = BuilderApp()
    assert "Pick a tool" in app3.status
    app3._set_tool("R")
    app3._place((5, 5), (5, 6))
    assert "ground" in app3.status.lower()
    app3.ground = (5, 5)
    app3._resolve()
    assert "source" in app3.status.lower()
    app3._set_tool("VSRC")
    app3._place((5, 6), (6, 6))          # (6,6) only touched once -> floating
    assert "not solvable" in app3.status.lower()
    print("builder : validation feedback (no ground/source, floating "
          "node) OK")

    # 10b) short circuit: a wire connecting both of a component's own
    #      terminals must never crash the app (it used to, uncaught, before
    #      this fix) and must give a plain-English status instead
    app6 = BuilderApp()
    app6._set_tool("VSRC")
    app6._place((0, 0), (1, 0))
    app6._set_tool("GROUND")
    app6.ground = (0, 0)
    app6._set_tool("WIRE")
    app6._place_wire_path([(0, 0), (0, 3)])
    app6._place_wire_path([(0, 3), (1, 3)])
    app6._place_wire_path([(1, 3), (1, 0)])   # closes a loop shorting VSRC1
    assert app6.result is None, app6.status
    assert "short circuit" in app6.status.lower() and "VSRC1" in app6.status
    app6._file_dialog = lambda save: base + "_short.json"
    app6.save_json()                          # must not raise mid-short
    assert "Saved" in app6.status
    print("builder : short-circuit detection (no crash, clear status) OK")

    # 11) clear_all leaves a consistent empty state
    app3.clear_all()
    assert not app3.placements and not app3.wires and app3.ground is None
    app3.fig.savefig(base + "_empty" + ext, dpi=100)

    # 12) a widget (e.g. the value TextBox) left mid-edit must not break the
    #     next click on a different widget -- matplotlib's TextBox grabs the
    #     canvas's mouse on click-in and only releases it via a full click
    #     elsewhere; get this wrong and every subsequent widget click raises
    #     "Another Axes already grabs mouse input" and is silently lost.
    #     Uses real MouseEvents (not the _FakeEvent/_click helper above,
    #     which only exercises our own canvas dispatch) since the bug lives
    #     inside matplotlib's own widget internals.
    from matplotlib.backend_bases import MouseEvent
    app7 = BuilderApp()
    app7._set_tool("C")
    app7._place((0, 0), (1, 0))
    app7._set_tool("SELECT")
    app7._select_for_edit(0)
    canvas7 = app7.fig.canvas

    def _press(ax, frac=0.5):
        x, y = ax.transAxes.transform((frac, 0.5))
        canvas7.callbacks.process(
            "button_press_event",
            MouseEvent("button_press_event", canvas7, x, y, button=1))

    def _release(ax, frac=0.5):
        x, y = ax.transAxes.transform((frac, 0.5))
        canvas7.callbacks.process(
            "button_release_event",
            MouseEvent("button_release_event", canvas7, x, y, button=1))

    _press(app7.box_value.ax)                 # click into value box, don't
    assert canvas7.mouse_grabber is app7.box_value.ax    # press Enter
    _press(app7.probe_btns["Q"].ax)           # must not raise RuntimeError
    _release(app7.probe_btns["Q"].ax)
    assert 0 in app7.probes_q, "Q click was lost to the stale grab"
    print("builder : stale TextBox mouse-grab doesn't break the next "
          "widget click OK")

    # 12b) the actual bug the user hit: select a component (SELECT tool,
    #     real canvas click, no TextBox involved at all), then click Q/E --
    #     still raised the same RuntimeError, for an unrelated reason: the
    #     Q/E buttons' axes geometrically overlap the AC/DC buttons' and
    #     the freq box's axes (both groups reuse the same screen rectangle
    #     since they're mutually exclusive), and a *hidden* Button is still
    #     hit-testable in matplotlib (Widget.ignore() only checks .active,
    #     not visibility) -- so one real click could reach two overlapping
    #     Buttons, both try to grab_mouse, and the second one collides.
    app8 = BuilderApp()
    app8._set_tool("C")
    app8._place((0, 0), (1, 0))
    app8._set_tool("SELECT")
    _click(app8, 0.5, 0)                      # real canvas-click selection,
    assert app8.selected == 0                 # like the user's own repro
    assert app8.src_btns["AC"].active is False and \
        app8.box_freq.active is False, \
        "overlapping source/freq widgets must be deactivated, not just hidden"

    canvas8 = app8.fig.canvas

    def _click_widget(btn, frac=0.15):
        # frac=0.15 (near the button's left edge, not dead center) matters:
        # Q's axes rect is wide enough that its exact center falls in the
        # gap between the (overlapping, hidden) AC/DC buttons underneath
        # it, which would make this test pass even without the fix -- an
        # off-center click is what actually lands inside the overlap.
        x, y = btn.ax.transAxes.transform((frac, 0.5))
        canvas8.callbacks.process(
            "button_press_event",
            MouseEvent("button_press_event", canvas8, x, y, button=1))
        canvas8.callbacks.process(
            "button_release_event",
            MouseEvent("button_release_event", canvas8, x, y, button=1))

    _click_widget(app8.probe_btns["Q"])       # must not raise
    assert 0 in app8.probes_q
    _click_widget(app8.probe_btns["E"])       # must not raise
    assert 0 in app8.probes_e
    print("builder : overlapping hidden widgets don't intercept clicks on "
          "the visible one OK")

    # 13) the Energy axis (a twinx() of the charge chart) must stay on the
    #     right after a redraw -- cla() resets a twin axes back to the
    #     default left-side tick/label position, undoing what twinx() set
    #     up at construction time, unless every redraw re-asserts it
    app._redraw_charts()
    app._redraw_charts()                      # twice: catches a fix that
    assert app.ax_e.yaxis.label_position == "right"   # only survives once
    assert app.ax_e.yaxis.get_ticks_position() == "right"
    print("builder : Energy axis stays right-aligned across redraws OK")

    app.fig.savefig(out, dpi=110)
    print("saved   :", out)


if __name__ == "__main__":
    app = BuilderApp()
    if "--test" in sys.argv:
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
        out_arg = sys.argv[sys.argv.index("--test") + 1] if \
            len(sys.argv) > sys.argv.index("--test") + 1 else \
            "rlc_builder_test.png"
        _self_test(app, out_arg)
    else:
        print("Circuit Builder — close the window to quit.")
        print("Pick a tool from the palette, click grid points to place "
              "and wire components,")
        print("set a ground, add a source, and probe nodes/components to "
              "see results.")
        plt.show()
