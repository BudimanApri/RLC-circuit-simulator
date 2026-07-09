# -*- coding: utf-8 -*-
"""Interactive matplotlib UI for the series R / RC / RL / LC / RLC circuit
simulator. The physics lives in rlc_solver, the schematic in rlc_schematic,
and shared constants in rlc_config."""

import math
import time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.transforms import blended_transform_factory
from matplotlib.widgets import Slider, Button, CheckButtons, TextBox

from rlc_config import (DEF, BOUNDS, SLIDERS, TH, BADGE, TOPOS, TOPO_ORDER,
                        SOURCE_ORDER, FAMILY_ORDER, PARALLEL_ORDER,
                        PARALLEL_TOPOS, PARALLEL_I0_DEFAULT,
                        PARALLEL_I0_BOUNDS, PARALLEL_I0_SLIDER,
                        C_SRC, C_R, C_L, C_C)
from rlc_schematic import Schematic, ParallelSchematic
from rlc_solver import solve, solve_rk4, solve_parallel


class RLCApp:
    def __init__(self):
        self.playing = False
        self.vals = dict(DEF)
        self.source_mode = "AC"
        self.family = "Series"
        self.parallel_preset = "RC_P"
        self._amp_is_voltage = True
        self._amplitude_bounds = BOUNDS["E0"]
        self.rk4_data = None
        self._spans = []
        self._last_hover = 0.0
        self.hover_on = False
        self.v_leg = None
        self.p_leg = None
        self.pv_leg = None

        self.fig = plt.figure(figsize=(13.8, 8.5))
        self.fig.patch.set_facecolor(TH["bg"])
        try:
            self.fig.canvas.manager.set_window_title(
                "Series Circuit Simulator — R/RC/RL/LC/RLC")
        except Exception:
            pass

        # ---- title + equation (follows the selected topology) -----------------
        self.fig.text(0.035, 0.972, "Series Circuit Simulator",
                      fontsize=15, weight="bold", color=TH["text"],
                      va="center")
        self.txt_eq = self.fig.text(0.035, 0.941, "", fontsize=9,
                                    color=TH["sub"], va="center")

        # ---- AC / DC source toggle ---------------------------------------------
        self.fig.text(0.035, 0.9175, "Source:", fontsize=7.8,
                      weight="bold", color=TH["sub"], va="center")
        self.source_btns = {}
        for j, name in enumerate(SOURCE_ORDER):
            axb = self.fig.add_axes([0.100 + j * 0.048, 0.906, 0.043, 0.023])
            b = Button(axb, name, color="#e2e8f0", hovercolor="#cbd5e1")
            b.label.set_fontsize(7.5)
            for sp in axb.spines.values():
                sp.set_visible(False)
            b.on_clicked(lambda ev, nm=name: self._set_source(nm))
            self.source_btns[name] = b

        # ---- circuit family toggle: series vs. parallel -------------------------
        self.fig.text(0.225, 0.9175, "Family:", fontsize=7.8, weight="bold",
                      color=TH["sub"], va="center")
        self.family_btns = {}
        for j, name in enumerate(FAMILY_ORDER):
            axb = self.fig.add_axes([0.283 + j * 0.066, 0.906, 0.061, 0.023])
            b = Button(axb, name, color="#e2e8f0", hovercolor="#cbd5e1")
            b.label.set_fontsize(7.5)
            for sp in axb.spines.values():
                sp.set_visible(False)
            b.on_clicked(lambda ev, nm=name: self._set_family(nm))
            self.family_btns[name] = b

        # ---- chart axes: charge, current, voltages -----------------------------
        self.ax_q = self.fig.add_axes([0.055, 0.725, 0.585, 0.175])
        self.ax_i = self.fig.add_axes([0.055, 0.535, 0.585, 0.165],
                                      sharex=self.ax_q)
        self.ax_v = self.fig.add_axes([0.055, 0.345, 0.585, 0.165],
                                      sharex=self.ax_q)
        self.ax_q.set_ylabel("Charge  Q  (mC)", fontsize=9, color=TH["text"])
        self.ax_i.set_ylabel("Current  I  (A)", fontsize=9, color=TH["text"])
        self.ax_v.set_ylabel("Voltage  (V)", fontsize=9, color=TH["text"])
        self.ax_v.set_xlabel("t  (ms)", fontsize=9, color=TH["text"])
        for ax in (self.ax_q, self.ax_i, self.ax_v):
            ax.set_facecolor("white")
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            for sp in ("left", "bottom"):
                ax.spines[sp].set_color("#b8c2d4")
            ax.tick_params(colors=TH["sub"], labelsize=8.5)
            ax.grid(True, ls="--", lw=0.6, alpha=0.45, color="#c3cddd")
        plt.setp(self.ax_q.get_xticklabels(), visible=False)
        plt.setp(self.ax_i.get_xticklabels(), visible=False)

        env_style = dict(color=TH["env"], lw=1.0, ls=(0, (5, 2, 1, 2)),
                         alpha=0.85, zorder=1.8)
        self.ln_qenv_p, = self.ax_q.plot([], [], **env_style)
        self.ln_qenv_m, = self.ax_q.plot([], [], **env_style)
        self.ln_ienv_p, = self.ax_i.plot([], [], **env_style)
        self.ln_ienv_m, = self.ax_i.plot([], [], **env_style)
        self.ln_qss, = self.ax_q.plot([], [], "--", color=TH["ss"], lw=1.1,
                                      zorder=2.2)
        self.ln_iss, = self.ax_i.plot([], [], "--", color=TH["ss"], lw=1.1,
                                      zorder=2.2)
        self.ln_q, = self.ax_q.plot([], [], color=TH["q"], lw=1.7, zorder=3)
        self.ln_i, = self.ax_i.plot([], [], color=TH["i"], lw=1.7, zorder=3)
        self.ln_qrk, = self.ax_q.plot([], [], ":", color=TH["rk"], lw=1.4,
                                      zorder=3.5)
        self.ln_irk, = self.ax_i.plot([], [], ":", color=TH["rk"], lw=1.4,
                                      zorder=3.5)

        # voltage traces, colored to match the schematic components
        self.ln_e, = self.ax_v.plot([], [], color=C_SRC, lw=1.0, alpha=0.75,
                                    zorder=2.5)
        self.ln_vr, = self.ax_v.plot([], [], color=C_R, lw=1.3, zorder=3)
        self.ln_vl, = self.ax_v.plot([], [], color=C_L, lw=1.3, zorder=3)
        self.ln_vc, = self.ax_v.plot([], [], color=C_C, lw=1.3, zorder=3)

        self.fig_legend = self.fig.legend(
            [self.ln_q, self.ln_i, self.ln_qss, self.ln_qrk, self.ln_qenv_p],
            ["Q(t) analytic", "I(t) analytic", "steady-state",
             "RK4 (numeric)", "transient envelope"],
            loc="upper right", bbox_to_anchor=(0.992, 0.998), ncol=5,
            fontsize=8, frameon=False, handlelength=1.7,
            columnspacing=1.1, handletextpad=0.5, labelcolor=TH["text"])

        # initial-condition markers Q(0)=0, I(0)=0 (series-mode only)
        self.dot_q0, = self.ax_q.plot(0, 0, "o", color=TH["q"], ms=6,
                                      zorder=5)
        self.dot_i0, = self.ax_i.plot(0, 0, "o", color=TH["i"], ms=6,
                                      zorder=5)
        self.ann_q0 = self.ax_q.annotate("Q(0) = 0", (0, 0),
                                         textcoords="offset points",
                                         xytext=(8, -14), fontsize=8,
                                         color=TH["q"])
        self.ann_i0 = self.ax_i.annotate("I(0) = 0", (0, 0),
                                         textcoords="offset points",
                                         xytext=(8, -14), fontsize=8,
                                         color=TH["i"])

        # ---- parallel-mode chart artists (Milestone 2) -------------------------
        # Share the same 3 axes as series mode but repurposed: ax_q -> voltage,
        # ax_i -> multi-trace branch currents, ax_v -> charge on C. All hidden
        # by default; toggled on when family == "Parallel".
        self.ln_pv, = self.ax_q.plot([], [], color=TH["q"], lw=1.7, zorder=3,
                                     visible=False)
        self.ln_pe, = self.ax_q.plot([], [], "--", color=C_SRC, lw=1.2,
                                     zorder=2.5, visible=False)
        self.ln_pir, = self.ax_i.plot([], [], color=C_R, lw=1.5, zorder=3,
                                      visible=False)
        self.ln_pil, = self.ax_i.plot([], [], color=C_L, lw=1.5, zorder=3,
                                      visible=False)
        self.ln_pic, = self.ax_i.plot([], [], color=C_C, lw=1.5, zorder=3,
                                      visible=False)
        self.ln_pitot, = self.ax_i.plot([], [], "--", color=TH["text"],
                                        lw=1.1, alpha=0.8, zorder=2.5,
                                        visible=False)
        self.ln_pqc, = self.ax_v.plot([], [], color=C_C, lw=1.7, zorder=3,
                                      visible=False)
        self.txt_no_c = self.ax_v.text(
            0.5, 0.5, "No capacitor in this circuit — nothing to plot here.",
            transform=self.ax_v.transAxes, ha="center", va="center",
            fontsize=9, color=TH["sub"], visible=False)
        self.txt_p = self.ax_q.text(0.987, 0.94, "", fontsize=8,
                                    transform=self.ax_q.transAxes,
                                    linespacing=1.3, ha="right", va="top",
                                    color="white", zorder=7, visible=False,
                                    bbox=dict(boxstyle="round,pad=0.45",
                                              fc="#0f172a", ec="none",
                                              alpha=0.82))

        # transient / steady-state markers (5τ)
        tr_q = blended_transform_factory(self.ax_q.transData,
                                         self.ax_q.transAxes)
        tr_v = blended_transform_factory(self.ax_v.transData,
                                         self.ax_v.transAxes)
        vln = dict(color="#d97706", ls=(0, (4, 3)), lw=1.2, alpha=0.9,
                   zorder=1.6, visible=False)
        self.vln_q = self.ax_q.axvline(0, **vln)
        self.vln_i = self.ax_i.axvline(0, **vln)
        self.vln_v = self.ax_v.axvline(0, **vln)
        self.lab_trans = self.ax_q.text(0, 0.95, "TRANSIENT", transform=tr_q,
                                        ha="center", va="top", fontsize=7.3,
                                        weight="bold", color=TH["trans"],
                                        alpha=0.9, zorder=4, visible=False)
        self.lab_steady = self.ax_q.text(0, 0.95, "STEADY-STATE",
                                         transform=tr_q, ha="center", va="top",
                                         fontsize=7.3, weight="bold",
                                         color=TH["steady"], alpha=0.9,
                                         zorder=4, visible=False)
        self.lab_5t = self.ax_v.text(0, 0.08, "", transform=tr_v, ha="left",
                                     va="bottom", fontsize=7.3,
                                     color="#9a3412", zorder=5, visible=False,
                                     bbox=dict(boxstyle="round,pad=0.32",
                                               fc="#fff7ed", ec="#fdba74",
                                               lw=0.8))

        # animation / hover cursor + running dots + readout boxes
        cur = dict(color="0.35", lw=0.9, visible=False, zorder=4.5)
        self.cur_q = self.ax_q.axvline(0, **cur)
        self.cur_i = self.ax_i.axvline(0, **cur)
        self.cur_v = self.ax_v.axvline(0, **cur)
        self.dot_q, = self.ax_q.plot([], [], "o", color=TH["q"], ms=7.5,
                                     mec="white", mew=1.3, visible=False,
                                     zorder=6)
        self.dot_i, = self.ax_i.plot([], [], "o", color=TH["i"], ms=7.5,
                                     mec="white", mew=1.3, visible=False,
                                     zorder=6)
        tip = dict(ha="right", va="top", color="white", zorder=7,
                   visible=False,
                   bbox=dict(boxstyle="round,pad=0.45", fc="#0f172a",
                             ec="none", alpha=0.82))
        self.txt_t = self.ax_q.text(0.987, 0.94, "", fontsize=8.5,
                                    transform=self.ax_q.transAxes,
                                    linespacing=1.35, **tip)
        self.txt_vt = self.ax_v.text(0.987, 0.93, "", fontsize=7.8,
                                     transform=self.ax_v.transAxes, **tip)

        # ---- card: circuit schematic + topology buttons -------------------------
        # The 5 button slots are index-addressable and relabeled/rebound on
        # every family switch: Series uses all 5 (TOPO_ORDER), Parallel uses
        # the first 4 (PARALLEL_ORDER) and hides the 5th slot.
        self._card([0.665, 0.700, 0.325, 0.235], "CIRCUIT")
        self.topo = "RLC"
        self.topo_btns = []
        for j in range(5):
            axb = self.fig.add_axes([0.775 + j * 0.0425, 0.906, 0.040, 0.023])
            b = Button(axb, "", color="#e2e8f0", hovercolor="#cbd5e1")
            b.label.set_fontsize(7.5)
            for sp in axb.spines.values():
                sp.set_visible(False)
            b.on_clicked(lambda ev, idx=j: self._topo_button_clicked(idx))
            self.topo_btns.append(b)
        self.ax_sch = self.fig.add_axes([0.673, 0.700, 0.309, 0.196])
        self.sch = Schematic(self.ax_sch)
        self.ax_psch = self.fig.add_axes([0.673, 0.700, 0.309, 0.196])
        self.ax_psch.set_visible(False)
        self.p_sch = ParallelSchematic(self.ax_psch)
        self.p_sch.set_preset(self.parallel_preset)

        # ---- card: analysis -------------------------------------------------------
        self._card([0.665, 0.487, 0.325, 0.198], "ANALYSIS")
        self.txt_badge = self.fig.text(0.978, 0.671, "", ha="right", va="top",
                                       fontsize=7.6, weight="bold",
                                       bbox=dict(boxstyle="round,pad=0.35",
                                                 fc="#dbeafe", ec="none"))
        self.txt_info = self.fig.text(0.676, 0.648, "", va="top", ha="left",
                                      fontsize=8.2, linespacing=1.45,
                                      color=TH["text"])

        # ---- card: impedance & resonance --------------------------------------------
        self._card([0.665, 0.325, 0.325, 0.145], "IMPEDANCE  &  RESONANCE")
        self.ax_ph = self.fig.add_axes([0.672, 0.330, 0.105, 0.112])
        self.ax_ph.axis("off")
        self.ph_R, = self.ax_ph.plot([], [], color="#64748b", lw=1.8,
                                     solid_capstyle="round")
        self.ph_X, = self.ax_ph.plot([], [], color="#7c3aed", lw=1.8,
                                     solid_capstyle="round")
        self.ph_Z, = self.ax_ph.plot([], [], color=TH["accent"], lw=2.4,
                                     solid_capstyle="round")
        self.ax_ph.plot(0, 0, "o", ms=3.5, color="#334155", zorder=5)
        an = dict(textcoords="offset points", fontsize=7, weight="bold")
        self.an_R = self.ax_ph.annotate("R", (0, 0), xytext=(0, -8),
                                        ha="center", color="#64748b", **an)
        self.an_X = self.ax_ph.annotate("X", (0, 0), xytext=(5, 0),
                                        ha="left", va="center",
                                        color="#7c3aed", **an)
        self.an_Z = self.ax_ph.annotate("Z", (0, 0), xytext=(-4, 5),
                                        ha="right", color=TH["accent"], **an)

        self.txt_z = self.fig.text(0.795, 0.428, "", fontsize=9.5,
                                   weight="bold", color=TH["text"])
        self.txt_phi = self.fig.text(0.908, 0.428, "", fontsize=9.5,
                                     weight="bold", color=TH["text"])
        self.txt_lag = self.fig.text(0.795, 0.407, "", fontsize=6.8,
                                     color=TH["sub"])
        self.txt_ratio = self.fig.text(0.795, 0.387, "", fontsize=7.8,
                                       weight="bold", color=TH["sub"])
        self.ax_res = self.fig.add_axes([0.795, 0.330, 0.185, 0.050])
        self.ax_res.set_xlim(-2.3, 2.3)
        self.ax_res.set_ylim(0, 1)
        self.ax_res.axis("off")
        self.ax_res.add_patch(Rectangle((-2.05, 0.40), 4.1, 0.30,
                                        fc="#e2e8f0", ec="none"))
        self.ax_res.plot([0, 0], [0.30, 0.80], color="#94a3b8", lw=1.1)
        for xv, lab in [(-2, "ω₀/4"), (-1, "ω₀/2"), (0, "ω₀"),
                        (1, "2ω₀"), (2, "4ω₀")]:
            self.ax_res.text(xv, 0.24, lab, ha="center", va="top",
                             fontsize=6.2, color=TH["sub"])
        self.res_dot, = self.ax_res.plot([0], [0.86], marker="v", ms=6.5,
                                         mfc="#475569", mec="white", mew=0.8,
                                         ls="none", zorder=5, clip_on=False)
        self.txt_dc_note = self.fig.text(
            0.8275, 0.400, "Impedance and resonance are AC-only\n"
            "concepts. Switch the source to AC\nto see them here.",
            fontsize=8, ha="center", va="center", color=TH["sub"],
            linespacing=1.6, visible=False)

        # ---- card: solution formulas (assignment answer) -----------------------------
        self._card([0.030, 0.195, 0.955, 0.105], "SOLUTION  (ASSIGNMENT ANSWER)")
        self.fig.text(0.041, 0.2635, "●", fontsize=6, color=TH["q"])
        self.fig.text(0.041, 0.2375, "●", fontsize=6, color=TH["i"])
        self.txt_fq = self.fig.text(0.052, 0.260, "", fontsize=8.3,
                                    family="monospace", color=TH["text"])
        self.txt_fi = self.fig.text(0.052, 0.234, "", fontsize=8.3,
                                    family="monospace", color=TH["text"])
        self.txt_note = self.fig.text(0.052, 0.207, "", fontsize=7.5,
                                      color=TH["sub"])

        # ---- card: parameters (sliders + numeric boxes) --------------------------------
        self._card([0.030, 0.022, 0.535, 0.150],
                   "PARAMETERS   —   drag a slider, or type a value and "
                   "press Enter")
        self.sliders, self.boxes = {}, {}
        for j, (key, lab, lo, hi, st) in enumerate(SLIDERS):
            y = 0.128 - j * 0.0185
            axs = self.fig.add_axes([0.125, y, 0.29, 0.011])
            try:
                s = Slider(axs, lab, lo, hi, valinit=DEF[key], valstep=st,
                           color=TH["accent"], track_color="#dbe3ef",
                           initcolor="#b6c2d6",
                           handle_style=dict(facecolor="white",
                                             edgecolor=TH["accent"], size=9))
            except TypeError:
                s = Slider(axs, lab, lo, hi, valinit=DEF[key], valstep=st,
                           color=TH["accent"])
            s.label.set_fontsize(8.6)
            s.label.set_color(TH["text"])
            s.valtext.set_visible(False)
            s.on_changed(self._slider_changed(key))
            self.sliders[key] = s

            axb = self.fig.add_axes([0.432, y - 0.0035, 0.075, 0.0175])
            try:
                tb = TextBox(axb, "", initial=f"{DEF[key]:.6g}",
                             textalignment="center", color="#f8fafc",
                             hovercolor="#eef2f7")
            except TypeError:
                tb = TextBox(axb, "", initial=f"{DEF[key]:.6g}",
                             color="#f8fafc", hovercolor="#eef2f7")
            tb.text_disp.set_fontsize(8.4)
            tb.text_disp.set_color(TH["text"])
            for sp in axb.spines.values():
                sp.set_color("#c3cddd")
                sp.set_linewidth(1.0)
            tb.on_submit(self._box_submitted(key))
            self.boxes[key] = tb

        # ---- card: display & controls -----------------------------------------------
        self._card([0.585, 0.022, 0.405, 0.150], "DISPLAY  &  CONTROLS")
        ax_chk = self.fig.add_axes([0.595, 0.030, 0.155, 0.112])
        ax_chk.axis("off")
        labels = ["Steady-state", "RK4 (numeric)", "Transient envelope"]
        actives = [False, False, True]
        try:
            self.chk = CheckButtons(
                ax_chk, labels, actives,
                frame_props=dict(sizes=[52] * 3, facecolor="white",
                                 edgecolor="#94a3b8", linewidth=1.1),
                check_props=dict(sizes=[52] * 3, facecolor=TH["accent"]),
                label_props=dict(fontsize=[8.6] * 3, color=[TH["text"]] * 3))
        except TypeError:
            self.chk = CheckButtons(ax_chk, labels, actives)
        self.chk.on_clicked(self._on_change)

        axp = self.fig.add_axes([0.775, 0.098, 0.095, 0.048])
        self.btn_play = Button(axp, "▶  Play", color=TH["accent"],
                               hovercolor="#1e40af")
        self.btn_play.label.set_color("white")
        self.btn_play.label.set_fontsize(9.5)
        self.btn_play.label.set_fontweight("bold")
        self.btn_play.on_clicked(lambda ev: self.toggle_play())

        axr = self.fig.add_axes([0.882, 0.098, 0.095, 0.048])
        self.btn_reset = Button(axr, "Reset", color="#e2e8f0",
                                hovercolor="#cbd5e1")
        self.btn_reset.label.set_fontsize(9.5)
        self.btn_reset.label.set_color(TH["text"])
        self.btn_reset.on_clicked(lambda ev: self.reset())
        for axx in (axp, axr):
            for sp in axx.spines.values():
                sp.set_visible(False)

        self.fig.text(0.775, 0.070, "speed:", fontsize=7.5, color=TH["sub"])
        self.speed = 1
        self.speed_btns = {}
        for mult, xx in ((1, 0.836), (2, 0.874), (4, 0.912)):
            axb = self.fig.add_axes([xx, 0.060, 0.034, 0.030])
            b = Button(axb, f"{mult}×", color="#e2e8f0", hovercolor="#cbd5e1")
            b.label.set_fontsize(8)
            for sp in axb.spines.values():
                sp.set_visible(False)
            b.on_clicked(lambda ev, m=mult: self._set_speed(m))
            self.speed_btns[mult] = b
        self._set_speed(1)

        self.fig.text(0.775, 0.033,
                      "space = Play/Pause    r = Reset    "
                      "hover the charts to read values",
                      fontsize=6.7, color=TH["sub"])

        self.fig.canvas.mpl_connect("key_press_event", self._on_key)
        self.fig.canvas.mpl_connect("motion_notify_event", self._on_move)
        self.fig.canvas.mpl_connect("figure_leave_event", self._on_leave)
        self.timer = self.fig.canvas.new_timer(interval=28)
        self.timer.add_callback(self._tick)

        self._set_source("AC")
        self._set_family("Series")
        self._set_topo("RLC")

    # ------------------------------------------------------------------ card
    def _card(self, rect, title=None):
        x, y, w, h = rect
        box = FancyBboxPatch((x, y), w, h, transform=self.fig.transFigure,
                             boxstyle="round,pad=0,rounding_size=0.008",
                             mutation_aspect=1.62, fc=TH["card"],
                             ec=TH["edge"], lw=1.1, zorder=-1)
        self.fig.add_artist(box)
        if title:
            self.fig.text(x + 0.011, y + h - 0.0135, title, fontsize=7.8,
                          weight="bold", color=TH["sub"], va="top")

    # ------------------------------------------------------------------ util
    def params(self):
        v = self.vals
        return v["R"], v["L"], v["C"], v["E0"], v["W"]

    def _sync_box(self, key):
        tb = self.boxes[key]
        tb.eventson = False
        tb.set_val(f"{self.vals[key]:.6g}")
        tb.eventson = True

    def _slider_changed(self, key):
        def cb(val):
            self.vals[key] = float(val)
            self._sync_box(key)
            self._stop_anim()
            self.recompute()
        return cb

    def _bounds_for(self, key):
        return self._amplitude_bounds if key == "E0" else BOUNDS[key]

    def _box_submitted(self, key):
        def cb(text):
            try:
                v = float(text.strip().replace(",", "."))
            except ValueError:
                self._sync_box(key)          # invalid input → restore display
                return
            lo, hi = self._bounds_for(key)
            v = min(max(v, lo), hi)
            if v == self.vals[key]:
                self._sync_box(key)
                return
            self.vals[key] = v
            s = self.sliders[key]
            s.eventson = False
            s.set_val(min(max(v, s.valmin), s.valmax))
            s.eventson = True
            self._sync_box(key)
            self._stop_anim()
            self.recompute()
        return cb

    # ------------------------------------------------------------- topology
    def _enable_param(self, key, on):
        s = self.sliders[key]
        tb = self.boxes[key]
        s.set_active(on)
        tb.set_active(on)
        a_main = 1.0 if on else 0.25
        a_txt = 1.0 if on else 0.35
        s.poly.set_alpha(a_main)
        s.label.set_alpha(a_txt)
        for attr in ("track", "_handle"):
            art = getattr(s, attr, None)
            if art is not None:
                art.set_alpha(a_main)
        tb.text_disp.set_alpha(a_txt)
        tb.ax.set_facecolor("#f8fafc" if on else "#edf1f7")

    def _style_button(self, b, active):
        b.color = TH["accent"] if active else "#e2e8f0"
        b.hovercolor = "#1e40af" if active else "#cbd5e1"
        b.ax.set_facecolor(b.color)
        b.label.set_color("white" if active else TH["text"])
        b.label.set_fontweight("bold" if active else "normal")

    def _refresh_topo_buttons(self):
        """Relabel/restyle the 5 shared button slots for the active family
        (Series uses all 5; Parallel uses 4 and hides the 5th slot)."""
        if self.family == "Series":
            order, active_name = TOPO_ORDER, self.topo
        else:
            order, active_name = PARALLEL_ORDER, self.parallel_preset
        for j, b in enumerate(self.topo_btns):
            if j < len(order):
                name = order[j]
                label = name if self.family == "Series" \
                    else PARALLEL_TOPOS[name]["label"]
                b.label.set_text(label)
                b.ax.set_visible(True)
                self._style_button(b, name == active_name)
            else:
                b.ax.set_visible(False)

    def _topo_button_clicked(self, idx):
        if self.family == "Series":
            if idx < len(TOPO_ORDER):
                self._set_topo(TOPO_ORDER[idx])
        elif idx < len(PARALLEL_ORDER):
            self._set_parallel_preset(PARALLEL_ORDER[idx])

    def _set_topo(self, name):
        self.topo = name
        cfg = TOPOS[name]
        for key, on in (("R", cfg["R"]), ("L", cfg["L"]), ("C", cfg["C"])):
            self._enable_param(key, on)
        self._refresh_topo_buttons()
        self.sch.set_topology(cfg["R"], cfg["L"], cfg["C"])
        self._update_equation_text()
        self.ax_q.set_ylabel(
            "Charge  Q  (mC)" if cfg["C"] else "Charge delivered  Q  (mC)",
            fontsize=9, color=TH["text"])
        self._stop_anim()
        self.recompute()

    def _set_parallel_preset(self, name):
        self.parallel_preset = name
        self._refresh_topo_buttons()
        self.p_sch.set_preset(name)
        self._update_equation_text()
        self._configure_amplitude_slider()
        self._stop_anim()
        self.recompute()

    def _set_family(self, name):
        self.family = name
        for m, b in self.family_btns.items():
            self._style_button(b, m == name)

        is_series = (name == "Series")
        self.ax_sch.set_visible(is_series)
        self.ax_psch.set_visible(not is_series)
        self._set_series_artists_visible(is_series)
        self._set_parallel_artists_visible(not is_series)
        self._refresh_topo_buttons()
        self._update_equation_text()
        self._configure_amplitude_slider()

        if is_series:
            cfg = TOPOS[self.topo]
            self.ax_q.set_ylabel(
                "Charge  Q  (mC)" if cfg["C"] else "Charge delivered  Q  (mC)",
                fontsize=9, color=TH["text"])
            self.ax_i.set_ylabel("Current  I  (A)", fontsize=9,
                                 color=TH["text"])
            self.ax_v.set_ylabel("Voltage  (V)", fontsize=9,
                                 color=TH["text"])
        else:
            self.ax_q.set_ylabel("Voltage  V  (V)", fontsize=9,
                                 color=TH["text"])
            self.ax_i.set_ylabel("Branch currents  (A)", fontsize=9,
                                 color=TH["text"])
            self.ax_v.set_ylabel("Charge on C   Q_C  (mC)", fontsize=9,
                                 color=TH["text"])

        self._stop_anim()
        self.recompute()

    def _set_series_artists_visible(self, flag):
        for art in (self.ln_q, self.ln_i, self.ln_e, self.ln_vr, self.ln_vl,
                    self.ln_vc, self.ln_qss, self.ln_iss, self.ln_qrk,
                    self.ln_irk, self.ln_qenv_p, self.ln_qenv_m,
                    self.ln_ienv_p, self.ln_ienv_m, self.dot_q0,
                    self.dot_i0, self.ann_q0, self.ann_i0):
            art.set_visible(flag)
        self.fig_legend.set_visible(flag)
        if not flag and self.v_leg is not None:
            self.v_leg.set_visible(False)

    def _set_parallel_artists_visible(self, flag):
        for art in (self.ln_pv, self.ln_pir, self.ln_pil, self.ln_pic,
                    self.ln_pitot):
            art.set_visible(flag)
        if not flag:
            for art in (self.ln_pe, self.ln_pqc, self.txt_no_c, self.txt_p):
                art.set_visible(False)
            if self.p_leg is not None:
                self.p_leg.set_visible(False)
            if self.pv_leg is not None:
                self.pv_leg.set_visible(False)

    def _configure_amplitude_slider(self):
        """The amplitude slider is volts for Series and for the (voltage-
        driven) Tank preset, and amps for the current-driven parallel
        presets. Reconfigure its range/label, and reset its value only when
        the underlying *unit* actually changes (not on every preset click
        within the same unit)."""
        want_voltage = (self.family == "Series"
                        or self.parallel_preset == "TANK")
        s, tb = self.sliders["E0"], self.boxes["E0"]
        if want_voltage == self._amp_is_voltage:
            return
        self._amp_is_voltage = want_voltage
        if want_voltage:
            lo, hi, step = 0.0, 300.0, 5.0
            label, default = "E₀  (V)", DEF["E0"]
            self._amplitude_bounds = BOUNDS["E0"]
        else:
            lo, hi, step = PARALLEL_I0_SLIDER
            label, default = "I₀  (A)", PARALLEL_I0_DEFAULT
            self._amplitude_bounds = PARALLEL_I0_BOUNDS
        s.valmin, s.valmax, s.valstep = lo, hi, step
        s.ax.set_xlim(lo, hi)
        s.label.set_text(label)
        self.vals["E0"] = default
        s.eventson = False
        s.set_val(default)
        s.eventson = True
        self._sync_box("E0")

    def _set_source(self, mode):
        self.source_mode = mode
        for m, b in self.source_btns.items():
            self._style_button(b, m == mode)
        self._enable_param("W", mode == "AC")
        self._update_equation_text()
        self._stop_anim()
        self.recompute()

    def _update_equation_text(self):
        cfg = TOPOS[self.topo] if self.family == "Series" \
            else PARALLEL_TOPOS[self.parallel_preset]
        self.txt_eq.set_text(cfg["eq"] if self.source_mode == "AC"
                             else cfg["eq_dc"])

    def _set_speed(self, mult):
        self.speed = mult
        for m, b in self.speed_btns.items():
            active = (m == mult)
            b.color = TH["accent"] if active else "#e2e8f0"
            b.hovercolor = "#1e40af" if active else "#cbd5e1"
            b.ax.set_facecolor(b.color)
            b.label.set_color("white" if active else TH["text"])
            b.label.set_fontweight("bold" if active else "normal")
        self.fig.canvas.draw_idle()

    def _typing(self):
        return any(getattr(tb, "capturekeystrokes", False)
                   for tb in self.boxes.values())

    def _on_change(self, _=None):
        self._stop_anim()
        self.recompute()

    def _on_key(self, ev):
        if self._typing():
            return
        if ev.key == " ":
            self.toggle_play()
        elif ev.key == "r":
            self.reset()
        elif ev.key in ("1", "2", "4"):
            self._set_speed(int(ev.key))

    def reset(self):
        self._stop_anim()
        self.vals = dict(DEF)
        if not self._amp_is_voltage:
            self.vals["E0"] = PARALLEL_I0_DEFAULT
        for key, s in self.sliders.items():
            s.eventson = False
            s.set_val(self.vals[key])
            s.eventson = True
            self._sync_box(key)
        self.recompute()

    # ------------------------------------------------------------- computation
    def recompute(self):
        if self.family == "Series":
            self._recompute_series()
        else:
            self._recompute_parallel()

    def _recompute_series(self):
        R, L, Cuf, E0, w = self.params()
        Tms = self.vals["T"]
        mode = self.source_mode
        cfg = TOPOS[self.topo]
        Cf = Cuf * 1e-6
        k = 1.0 / Cf if cfg["C"] else 0.0
        Reff = R if cfg["R"] else 0.0

        # adaptive samples: ≥ ~36 points per fastest cycle / time constant
        # (a DC source has no oscillation of its own, so w only matters in AC)
        w_eff = w if mode == "AC" else 0.0
        if cfg["L"]:
            wfast = max(w_eff, math.sqrt(k / L), Reff / L)
        else:
            wfast = max(w_eff, k / max(Reff, 1e-9))
        self.npts = int(np.clip(Tms * 1e-3 * max(wfast, 1.0) / (2 * math.pi)
                                * 36.0, 3000, 48000))
        self.step = max(1, self.npts // 420)
        self.t = np.linspace(0.0, Tms * 1e-3, self.npts)
        self.tms = self.t * 1e3

        sol = solve(self.topo, R, L, Cuf, E0, w, self.t, mode=mode)
        self.sol = sol
        self.Q, self.I = sol["Q"], sol["I"]
        self.Qss, self.Iss = sol["Qss"], sol["Iss"]
        self.E, self.VR = sol["E"], sol["VR"]
        self.VL, self.VC = sol["VL"], sol["VC"]

        show_ss, show_rk, show_env = self.chk.get_status()
        self.ln_qss.set_visible(show_ss)
        self.ln_iss.set_visible(show_ss)

        for ln, dat in ((self.ln_qenv_p, sol["env_q_hi"] * 1e3),
                        (self.ln_qenv_m, sol["env_q_lo"] * 1e3),
                        (self.ln_ienv_p, sol["env_i_hi"]),
                        (self.ln_ienv_m, sol["env_i_lo"])):
            ln.set_data(self.tms, dat)
            ln.set_visible(show_env)

        # voltage traces: source always, components per topology
        self.ln_vr.set_visible(cfg["R"])
        self.ln_vl.set_visible(cfg["L"])
        self.ln_vc.set_visible(cfg["C"])
        if self.v_leg is not None:
            self.v_leg.remove()
        handles, labels = [self.ln_e], ["E(t)"]
        for flag, ln, lab in ((cfg["R"], self.ln_vr, "V_R"),
                              (cfg["L"], self.ln_vl, "V_L"),
                              (cfg["C"], self.ln_vc, "V_C")):
            if flag:
                handles.append(ln)
                labels.append(lab)
        self.v_leg = self.ax_v.legend(
            handles, labels, loc="lower right", bbox_to_anchor=(1.0, 1.0),
            ncol=4, fontsize=6.4, frameon=False, handlelength=1.2,
            columnspacing=0.9, handletextpad=0.4, borderaxespad=0.0,
            labelcolor=TH["sub"])

        note = ("Solid lines = exact analytic solution.   Orange envelope = "
                "steady-state amplitude + transient decay.   Enable "
                "'RK4 (numeric)' for independent verification.")
        self.ln_qrk.set_visible(show_rk)
        self.ln_irk.set_visible(show_rk)
        self.rk4_data = None
        if show_rk:
            Qr, Ir = solve_rk4(self.topo, R, L, Cuf, E0, w, self.t, mode=mode)
            self.rk4_data = (Qr, Ir)
            dq = float(np.max(np.abs(Qr - self.Q)))
            di = float(np.max(np.abs(Ir - self.I)))
            note = (f"RK4 verification:  max|ΔQ| = {dq:.2e} C,  "
                    f"max|ΔI| = {di:.2e} A   →   numeric ≡ analytic ✓")
        self.txt_note.set_text(note)

        tmax = self.tms[-1]
        self._update_transient_marker(sol["alpha_settle"], tmax)

        # ---- axis limits (fixed during animation; asymmetric when offset) -------
        def limits(data, lo_arr, hi_arr):
            lo = float(np.min(data))
            hi = float(np.max(data))
            if show_env:
                lo = min(lo, float(np.min(lo_arr)))
                hi = max(hi, float(np.max(hi_arr)))
            span = max(hi - lo, 1e-12)
            return lo - 0.18 * span, hi + 0.18 * span

        qlo, qhi = limits(self.Q, sol["env_q_lo"], sol["env_q_hi"])
        ilo, ihi = limits(self.I, sol["env_i_lo"], sol["env_i_hi"])
        self.ax_q.set_xlim(0, tmax)
        self.ax_q.set_ylim(qlo * 1e3, qhi * 1e3)
        self.ax_i.set_ylim(ilo, ihi)
        vm = float(np.max(np.abs(self.E)))
        for flag, arr in ((cfg["R"], self.VR), (cfg["L"], self.VL),
                          (cfg["C"], self.VC)):
            if flag:
                vm = max(vm, float(np.max(np.abs(arr))))
        vm = max(vm, 1e-9)
        self.ax_v.set_ylim(-1.15 * vm, 1.15 * vm)

        # ---- schematic labels follow the parameters ------------------------------
        self.sch.labels["R"].set_text(f"R = {R:.4g} Ω")
        self.sch.labels["L"].set_text(f"L = {L:.4g} H")
        self.sch.labels["C"].set_text(f"C = {Cuf:.4g} µF")
        self.sch.labels["E"].set_text(
            f"E(t) = {E0:.4g}·sin({w:.4g}·t)  V" if mode == "AC"
            else f"E(t) = {E0:.4g} V  (DC step)")

        # ---- analysis panel ---------------------------------------------------------
        for key, (fg, bg) in BADGE.items():
            if key in sol["damping"]:
                self.txt_badge.set_color(fg)
                self.txt_badge.get_bbox_patch().set_facecolor(bg)
                break
        self.txt_badge.set_text(sol["damping"])
        if mode == "AC":
            self.txt_info.set_text(
                self._analysis_text_ac(sol, R, L, Cuf, E0, w))
        else:
            self.txt_info.set_text(
                self._analysis_text_dc(sol, R, L, Cuf, E0))

        # ---- impedance & resonance panel (AC-only concept) ----------------------------
        show_ac_panel = (mode == "AC")
        self.ax_ph.set_visible(show_ac_panel)
        for t_ in (self.txt_z, self.txt_phi, self.txt_lag, self.txt_ratio):
            t_.set_visible(show_ac_panel)
        self.txt_dc_note.set_visible(not show_ac_panel)
        self.txt_dc_note.set_text(
            "Impedance and resonance are AC-only\nconcepts. Switch the "
            "source to AC\nto see them here.")
        show_gauge = show_ac_panel and cfg["L"] and cfg["C"]
        self.ax_res.set_visible(show_gauge)

        if show_ac_panel:
            X = sol["XL"] - sol["XC"]
            self.ph_R.set_data([0, Reff], [0, 0])
            self.ph_X.set_data([Reff, Reff], [0, X])
            self.ph_Z.set_data([0, Reff], [0, X])
            self.an_R.xy = (Reff / 2.0, 0)
            self.an_X.xy = (Reff, X / 2.0)
            self.an_Z.xy = (Reff / 2.0, X / 2.0)
            self.an_R.set_visible(cfg["R"])
            self.an_X.set_visible(cfg["L"] or cfg["C"])
            # manual limits so x-scale = y-scale (an honest triangle)
            pos = self.ax_ph.get_position()
            wa = (pos.width * self.fig.get_figwidth()) / \
                 (pos.height * self.fig.get_figheight())
            ylo, yhi = min(0.0, X), max(0.0, X)
            dx = 1.30 * max(max(Reff, 1e-9), wa * max(yhi - ylo, 1e-9))
            dy = dx / wa
            self.ax_ph.set_xlim(Reff / 2.0 - dx / 2.0, Reff / 2.0 + dx / 2.0)
            self.ax_ph.set_ylim((ylo + yhi) / 2.0 - dy / 2.0,
                                (ylo + yhi) / 2.0 + dy / 2.0)

            self.txt_z.set_text(f"|Z| = {sol['Z']:.6g} Ω")
            self.txt_phi.set_text(f"φ = {sol['phi']:+.3g}°")
            if sol["phi"] > 0.5:
                tag = "inductive  —  current lags the voltage"
            elif sol["phi"] < -0.5:
                tag = "capacitive  —  current leads the voltage"
            elif show_gauge:
                tag = "≈ purely resistive  (very near resonance)"
            else:
                tag = "resistive  —  current in phase with the voltage"
            self.txt_lag.set_text(tag)

            if show_gauge:
                ratio = w / sol["w0"]
                near = abs(ratio - 1.0) < 0.05
                self.res_dot.set_data(
                    [float(np.clip(math.log2(ratio), -2.05, 2.05))], [0.86])
                self.res_dot.set_markerfacecolor(
                    "#16a34a" if near else "#475569")
                self.txt_ratio.set_text(
                    f"ω/ω₀ = {ratio:.3f}"
                    + ("    —  near resonance!" if near else ""))
                self.txt_ratio.set_color("#15803d" if near else TH["sub"])
            else:
                self.txt_ratio.set_text("resonance requires both L and C")
                self.txt_ratio.set_color(TH["sub"])

        self.txt_fq.set_text(sol["fq"])
        self.txt_fi.set_text(sol["fi"])

        self.idx = self.npts
        self._draw_upto(self.npts, final=True)

    def _update_transient_marker(self, alpha_settle, tmax):
        """Shade the transient region and place the 5τ marker on ax_q/ax_i/
        ax_v. Generic on alpha_settle alone, so both the series and the
        parallel recompute() paths share this."""
        for p in self._spans:
            p.remove()
        self._spans = []
        a_s = alpha_settle
        if a_s == math.inf:
            t5 = 0.0                      # no transient at all
        elif a_s > 0:
            t5 = 5.0 / a_s * 1e3
        else:
            t5 = math.inf                 # transient never decays

        self.lab_steady.set_text("STEADY-STATE")
        if t5 == 0.0:
            for ln in (self.vln_q, self.vln_i, self.vln_v):
                ln.set_visible(False)
            self.lab_trans.set_visible(False)
            self.lab_steady.set_visible(True)
            self.lab_steady.set_position((tmax / 2.0, 0.95))
            self.lab_steady.set_text("STEADY-STATE   (no transient)")
            self.lab_5t.set_visible(False)
        elif math.isfinite(t5):
            end = min(t5, tmax)
            for ax in (self.ax_q, self.ax_i, self.ax_v):
                self._spans.append(ax.axvspan(0.0, end, fc=TH["env"],
                                              alpha=0.07, lw=0, zorder=0.6))
            inside = t5 < 0.985 * tmax
            for ln in (self.vln_q, self.vln_i, self.vln_v):
                ln.set_visible(inside)
                if inside:
                    ln.set_xdata([t5, t5])
            if inside:
                self.lab_trans.set_visible(t5 > 0.10 * tmax)
                self.lab_trans.set_position((t5 / 2.0, 0.95))
                self.lab_trans.set_text("TRANSIENT")
                self.lab_steady.set_visible((tmax - t5) > 0.14 * tmax)
                self.lab_steady.set_position(((t5 + tmax) / 2.0, 0.95))
                ha = "left" if t5 <= 0.70 * tmax else "right"
                dx = 0.012 * tmax if ha == "left" else -0.012 * tmax
                self.lab_5t.set_visible(True)
                self.lab_5t.set_ha(ha)
                self.lab_5t.set_position((t5 + dx, 0.08))
                self.lab_5t.set_text(
                    f"transient practically over:  5τ ≈ {t5:.3g} ms")
            else:
                self.lab_trans.set_visible(True)
                self.lab_trans.set_position((tmax / 2.0, 0.95))
                self.lab_trans.set_text(
                    f"TRANSIENT   (5τ ≈ {t5:.3g} ms  >  time window)")
                self.lab_steady.set_visible(False)
                self.lab_5t.set_visible(False)
        else:
            for ln in (self.vln_q, self.vln_i, self.vln_v):
                ln.set_visible(False)
            self.lab_trans.set_visible(False)
            self.lab_steady.set_visible(False)
            self.lab_5t.set_visible(False)

    # --------------------------------------------------------- parallel compute
    def _recompute_parallel(self):
        R, L, Cuf, E0, w = self.params()
        Tms = self.vals["T"]
        mode = self.source_mode
        preset = self.parallel_preset
        pcfg = PARALLEL_TOPOS[preset]
        Cf = Cuf * 1e-6
        Rr = max(R, 1e-9)

        w_eff = w if mode == "AC" else 0.0
        wfast = w_eff
        if pcfg["C"] and pcfg["L"]:
            wfast = max(wfast, math.sqrt(1.0 / (L * Cf)))
        elif pcfg["C"]:
            wfast = max(wfast, 1.0 / (Rr * Cf))
        elif pcfg["L"]:
            wfast = max(wfast, Rr / L)
        self.npts = int(np.clip(Tms * 1e-3 * max(wfast, 1.0) / (2 * math.pi)
                                * 36.0, 3000, 48000))
        self.step = max(1, self.npts // 420)
        self.t = np.linspace(0.0, Tms * 1e-3, self.npts)
        self.tms = self.t * 1e3

        sol = solve_parallel(preset, R, L, Cuf, E0, w, self.t, mode=mode)
        self.psol = sol

        # -- chart 1 (ax_q): voltage ---------------------------------------------
        self.ln_pv.set_data(self.tms, sol["V"])
        show_e = (sol["src_kind"] == "E")
        self.ln_pe.set_visible(show_e)
        if show_e:
            self.ln_pe.set_data(self.tms, sol["Src"])
        if self.pv_leg is not None:
            self.pv_leg.remove()
            self.pv_leg = None
        if show_e:
            self.pv_leg = self.ax_q.legend(
                [self.ln_pv, self.ln_pe], ["V_tank(t)", "E(t) source"],
                loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2,
                fontsize=6.8, frameon=False, handlelength=1.4,
                columnspacing=1.0, handletextpad=0.4, borderaxespad=0.0,
                labelcolor=TH["sub"])

        # -- chart 2 (ax_i): branch currents -------------------------------------
        for ln, arr in ((self.ln_pir, sol["I_R"]), (self.ln_pil, sol["I_L"]),
                        (self.ln_pic, sol["I_C"])):
            present = arr is not None
            ln.set_visible(present)
            if present:
                ln.set_data(self.tms, arr)
        self.ln_pitot.set_data(self.tms, sol["I_total"])
        if self.p_leg is not None:
            self.p_leg.remove()
        handles, labels = [], []
        if sol["I_R"] is not None:
            handles.append(self.ln_pir)
            labels.append("I_R")
        if sol["I_L"] is not None:
            handles.append(self.ln_pil)
            labels.append("I_L")
        if sol["I_C"] is not None:
            handles.append(self.ln_pic)
            labels.append("I_C")
        handles.append(self.ln_pitot)
        labels.append("I_total" if pcfg["src"] == "I" else "I (through R)")
        self.p_leg = self.ax_i.legend(
            handles, labels, loc="lower right", bbox_to_anchor=(1.0, 1.0),
            ncol=4, fontsize=6.4, frameon=False, handlelength=1.2,
            columnspacing=0.9, handletextpad=0.4, borderaxespad=0.0,
            labelcolor=TH["sub"])

        # -- chart 3 (ax_v): charge on C ------------------------------------------
        has_c = sol["Q_C"] is not None
        self.ln_pqc.set_visible(has_c)
        self.txt_no_c.set_visible(not has_c)
        if has_c:
            self.ln_pqc.set_data(self.tms, sol["Q_C"] * 1e3)

        self.txt_note.set_text(
            "Parallel mode: steady-state overlay, RK4 verification, and "
            "the transient envelope are not available yet for parallel "
            "circuits (planned for a future update).")

        tmax = self.tms[-1]
        self._update_transient_marker(sol["alpha_settle"], tmax)

        # -- axis limits ------------------------------------------------------------
        def limits(data):
            lo, hi = float(np.min(data)), float(np.max(data))
            span = max(hi - lo, 1e-12)
            return lo - 0.18 * span, hi + 0.18 * span

        vlo, vhi = limits(np.concatenate([sol["V"], sol["Src"]])
                          if show_e else sol["V"])
        self.ax_q.set_xlim(0, tmax)
        self.ax_q.set_ylim(vlo, vhi)

        i_arrays = [sol["I_total"]]
        for arr in (sol["I_R"], sol["I_L"], sol["I_C"]):
            if arr is not None:
                i_arrays.append(arr)
        ilo, ihi = limits(np.concatenate(i_arrays))
        self.ax_i.set_ylim(ilo, ihi)

        if has_c:
            qlo, qhi = limits(sol["Q_C"] * 1e3)
            self.ax_v.set_ylim(qlo, qhi)

        # -- schematic labels ---------------------------------------------------------
        self.p_sch.labels["R"].set_text(f"R = {R:.4g} Ω")
        if pcfg["L"]:
            self.p_sch.labels["L"].set_text(f"L = {L:.4g} H")
        if pcfg["C"]:
            self.p_sch.labels["C"].set_text(f"C = {Cuf:.4g} µF")
        unit, sym = ("V", "E") if pcfg["src"] == "E" else ("A", "I")
        if mode == "AC":
            self.p_sch.labels["E"].set_text(
                f"{sym}(t) = {E0:.4g}·sin({w:.4g}·t)  {unit}")
        else:
            self.p_sch.labels["E"].set_text(
                f"{sym}(t) = {E0:.4g} {unit}  (DC step)")

        # -- analysis panel -----------------------------------------------------------
        for key, (fg, bg) in BADGE.items():
            if key in sol["damping"]:
                self.txt_badge.set_color(fg)
                self.txt_badge.get_bbox_patch().set_facecolor(bg)
                break
        self.txt_badge.set_text(sol["damping"])
        self.txt_info.set_text(
            self._analysis_text_parallel(sol, R, L, Cuf, E0, mode))

        # -- impedance & resonance panel ------------------------------------------------
        # The R-X-Z triangle is a series-circuit visual; skip it here and show
        # Zp/phi_p numerically (only meaningful for R∥L∥C and Tank in AC mode).
        show_ac_panel = (mode == "AC" and sol.get("Zp") is not None)
        self.ax_ph.set_visible(False)
        for t_ in (self.txt_z, self.txt_phi, self.txt_lag, self.txt_ratio):
            t_.set_visible(show_ac_panel)
        self.ax_res.set_visible(show_ac_panel)
        self.txt_dc_note.set_visible(not show_ac_panel)
        if show_ac_panel:
            self.txt_z.set_text(f"|Zp| = {sol['Zp']:.6g} Ω")
            self.txt_phi.set_text(f"φ = {sol['phi_p']:+.3g}°")
            if sol["phi_p"] > 0.5:
                tag = "inductive-leaning  —  V lags the source current"
            elif sol["phi_p"] < -0.5:
                tag = "capacitive-leaning  —  V leads the source current"
            else:
                tag = "≈ resistive  (very near resonance)"
            self.txt_lag.set_text(tag)
            ratio = w / sol["w0"] if sol["w0"] else 1.0
            near = abs(ratio - 1.0) < 0.05
            self.res_dot.set_data(
                [float(np.clip(math.log2(ratio), -2.05, 2.05))], [0.86])
            self.res_dot.set_markerfacecolor("#16a34a" if near else "#475569")
            self.txt_ratio.set_text(
                f"ω/ω₀ = {ratio:.3f}"
                + ("    —  near resonance!" if near else ""))
            self.txt_ratio.set_color("#15803d" if near else TH["sub"])
        elif mode == "AC":
            self.txt_dc_note.set_text(
                f"{pcfg['label']} has no resonance — only R∥L∥C\nand Tank "
                "do. Pick one of those presets to\nsee impedance & "
                "antiresonance here.")
        else:
            self.txt_dc_note.set_text(
                "Impedance and resonance are AC-only\nconcepts. Switch "
                "the source to AC\nto see them here.")

        self.txt_fq.set_text(sol["fV"])
        self.txt_fi.set_text(sol["fI"])

        self.idx = self.npts
        self._draw_upto(self.npts, final=True)

    def _analysis_text_ac(self, sol, R, L, Cuf, E0, w):
        cfg = TOPOS[self.topo]
        Cf = Cuf * 1e-6
        a_s = sol["alpha_settle"]
        if a_s == math.inf:
            settle = "no transient: current is immediately steady-state"
        elif a_s > 0:
            settle = f"transient < 1% after 5τ ≈ {5.0 / a_s * 1e3:.3g} ms"
        else:
            settle = "transient never decays (no damping)"

        vparts = []
        if cfg["R"]:
            vparts.append(f"V_R ≈ {sol['VRm']:.4g} V")
        if cfg["L"]:
            vparts.append(f"V_L ≈ {sol['VLm']:.4g} V")
        if cfg["C"]:
            vparts.append(f"V_C ≈ {sol['VCm']:.4g} V")
        vline = "V_max (steady):  " + ",   ".join(vparts)

        if cfg["L"] and cfg["C"]:                # RLC / LC (full 2nd order)
            Rterm = R if cfg["R"] else 0.0
            tau_s = "∞" if a_s == 0 else f"{1.0 / a_s * 1e3:.3g} ms"
            if "wd" in sol:
                mid = f"ω_d (transient oscillation) = {sol['wd']:.4g} rad/s"
                tau_lab = f"τ = 1/α = {tau_s}"
            elif "r2" in sol:
                mid = (f"real roots:  r₁ = {sol['r1']:.4g},   "
                       f"r₂ = {sol['r2']:.4g}  s⁻¹")
                tau_lab = f"τ_slow = 1/|r₁| = {tau_s}"
            else:
                mid = f"double root:  r = {sol['r1']:.4g}  s⁻¹"
                tau_lab = f"τ = 1/α = {tau_s}"
            return (f"R²  vs  4L/C :   {Rterm * Rterm:.3g}   vs   "
                    f"{4 * L / Cf:.3g}\n"
                    f"α = R/(2L) = {sol['alpha']:.4g} s⁻¹       {tau_lab}\n"
                    f"{mid}\n"
                    f"ω₀ = 1/√(LC) = {sol['w0']:.4g} rad/s\n"
                    f"X_L = ωL = {sol['XL']:.4g} Ω      "
                    f"X_C = 1/(ωC) = {sol['XC']:.4g} Ω\n"
                    f"Steady-state:  I_max = E₀/|Z| = {sol['Iamp']:.4g} A,"
                    f"   Q_max = {sol['Qamp'] * 1e3:.4g} mC\n"
                    f"{vline}\n"
                    f"{settle}")
        if cfg["L"]:                             # RL
            if R <= 0:
                return ("R = 0  →  pure inductor:   L·I′ = E(t)\n"
                        f"I(t) = (E₀/(ωL))·(1 − cos ωt):  "
                        f"current offset never decays\n"
                        f"X_L = ωL = {sol['XL']:.4g} Ω\n"
                        f"I_max = 2·E₀/(ωL) = {2 * E0 / (w * L):.4g} A\n"
                        f"{vline}\n"
                        f"{settle}")
            return (f"Effective equation (1st-order in I):   "
                    f"L·I′ + R·I = E(t)\n"
                    f"τ = L/R = {L / R * 1e3:.3g} ms       "
                    f"α = R/L = {R / L:.4g} s⁻¹\n"
                    f"roots:  r₁ = 0 (charge offset),   "
                    f"r₂ = −R/L = {sol['r2']:.4g} s⁻¹\n"
                    f"X_L = ωL = {sol['XL']:.4g} Ω      "
                    f"X_C = 0  (no capacitor)\n"
                    f"Steady-state:  I_max = E₀/|Z| = {sol['Iamp']:.4g} A,"
                    f"   ΔQ = {sol['Qamp'] * 1e3:.4g} mC\n"
                    f"{vline}\n"
                    f"{settle}")
        if cfg["C"]:                             # RC
            return (f"Equation (1st-order):   R·Q′ + Q/C = E(t)\n"
                    f"τ = RC = {max(R, 1e-9) * Cf * 1e3:.3g} ms       "
                    f"α = 1/(RC) = {sol['alpha']:.4g} s⁻¹\n"
                    f"X_C = 1/(ωC) = {sol['XC']:.4g} Ω      "
                    f"X_L = 0  (no inductor)\n"
                    f"Steady-state:  I_max = E₀/|Z| = {sol['Iamp']:.4g} A,"
                    f"   Q_max = {sol['Qamp'] * 1e3:.4g} mC\n"
                    f"{vline}\n"
                    f"{settle}")
        return (f"Algebraic:   I(t) = E(t)/R    (no ODE)\n"
                f"Z = R = {sol['Z']:.6g} Ω       φ = 0°  "
                f"(current in phase with the voltage)\n"
                f"I_max = E₀/R = {sol['Iamp']:.4g} A\n"
                f"Q(t) = charge delivered:  oscillates "
                f"0 … {2 * sol['Qamp'] * 1e3:.4g} mC\n"
                f"{vline}\n"
                f"{settle}")

    def _analysis_text_dc(self, sol, R, L, Cuf, E0):
        cfg = TOPOS[self.topo]
        Cf = Cuf * 1e-6
        a_s = sol["alpha_settle"]
        if a_s == math.inf:
            settle = "no transient: current is instantly at its final value"
        elif a_s == 0:
            settle = "never settles (undamped oscillation or unbounded ramp)"
        else:
            settle = f"transient < 1% after 5τ ≈ {5.0 / a_s * 1e3:.3g} ms"

        if cfg["L"] and cfg["C"]:                # RLC / LC (full 2nd order)
            Rterm = R if cfg["R"] else 0.0
            Qp = E0 * Cf
            if "wd" in sol:
                mid = f"ω_d (oscillation) = {sol['wd']:.4g} rad/s"
            elif "r2" in sol:
                mid = (f"real roots:  r₁ = {sol['r1']:.4g},   "
                       f"r₂ = {sol['r2']:.4g}  s⁻¹")
            else:
                mid = f"double root:  r = {sol['r1']:.4g}  s⁻¹"
            if a_s > 0:
                final = (f"Final (t→∞):  I → 0 A,   V_C → E₀ = {E0:.4g} V,"
                         f"   Q → {Qp * 1e3:.4g} mC")
            else:
                final = (f"No true steady state (R = 0): Q oscillates "
                         f"0 … {2 * Qp * 1e3:.4g} mC forever")
            return (f"R²  vs  4L/C :   {Rterm * Rterm:.3g}   vs   "
                    f"{4 * L / Cf:.3g}\n"
                    f"α = R/(2L) = {sol['alpha']:.4g} s⁻¹\n"
                    f"{mid}\n"
                    f"ω₀ = 1/√(LC) = {sol['w0']:.4g} rad/s\n"
                    f"Q(t) = {Qp * 1e3:.4g} mC·(1 − decaying/oscillating "
                    f"term)\n"
                    f"{final}\n"
                    f"{settle}")
        if cfg["L"]:                             # RL
            if R <= 0:
                return ("R = 0  →  pure inductor:   L·I′ = E₀\n"
                        f"I(t) = (E₀/L)·t,   Q(t) = (E₀/2L)·t²   "
                        f"(both ramp forever, no steady state)\n"
                        f"{settle}")
            tau = L / R
            return (f"Effective equation (1st-order in I):   "
                    f"L·I′ + R·I = E₀\n"
                    f"τ = L/R = {tau * 1e3:.3g} ms       "
                    f"α = R/L = {R / L:.4g} s⁻¹\n"
                    f"I(t) = (E₀/R)·(1 − e^(−t/τ))\n"
                    f"Final (t→∞):  I → E₀/R = {E0 / R:.4g} A,   "
                    f"V_R → {E0:.4g} V,   V_L → 0 V\n"
                    f"{settle}")
        if cfg["C"]:                             # RC
            Rr = max(R, 1e-9)
            tau = Rr * Cf
            return (f"Equation (1st-order):   R·Q′ + Q/C = E₀\n"
                    f"τ = RC = {tau * 1e3:.3g} ms       "
                    f"α = 1/(RC) = {sol['alpha']:.4g} s⁻¹\n"
                    f"Q(t) = E₀C·(1 − e^(−t/τ)),   "
                    f"I(t) = (E₀/R)·e^(−t/τ)\n"
                    f"Final (t→∞):  Q → {E0 * Cf * 1e3:.4g} mC,   "
                    f"V_C → {E0:.4g} V,   I → 0 A\n"
                    f"{settle}")
        Rr = max(R, 1e-9)
        return (f"Algebraic:   I(t) = E₀/R = {E0 / Rr:.4g} A    "
                f"(no ODE, instant response)\n"
                f"Q(t) = (E₀/R)·t: charge delivered grows without bound\n"
                f"{settle}")

    def _analysis_text_parallel(self, sol, R, L, Cuf, E0, mode):
        preset = self.parallel_preset
        pcfg = PARALLEL_TOPOS[preset]
        Cf = Cuf * 1e-6
        Rr = max(R, 1e-9)
        a_s = sol["alpha_settle"]
        if a_s == math.inf:
            settle = "no transient: instantly at its final value"
        elif a_s == 0:
            settle = "never settles (undamped oscillation or unbounded ramp)"
        else:
            settle = f"transient < 1% after 5τ ≈ {5.0 / a_s * 1e3:.3g} ms"
        src_line = ("Source: current-driven" if pcfg["src"] == "I"
                    else "Source: voltage-driven")

        if preset == "RC_P":
            tau = Rr * Cf
            return (f"{src_line}      τ = RC = {tau * 1e3:.3g} ms      "
                    f"α = 1/(RC) = {sol['alpha']:.4g} s⁻¹\n"
                    f"KCL:  I(t) = I_R(t) + I_C(t)      "
                    f"(I_R = V/R,   I_C = C·V′)\n"
                    f"{settle}")
        if preset == "RL_P":
            tau = L / Rr
            return (f"{src_line}      τ = L/R = {tau * 1e3:.3g} ms      "
                    f"α = R/L = {sol['alpha']:.4g} s⁻¹\n"
                    f"KCL:  I(t) = I_R(t) + I_L(t)      "
                    f"(V = L·I_L′,   I_R = V/R)\n"
                    f"{settle}")

        if "wd" in sol:
            mid = f"ω_d (oscillation) = {sol['wd']:.4g} rad/s"
        elif "r2" in sol:
            mid = (f"real roots:  r₁ = {sol['r1']:.4g},   "
                   f"r₂ = {sol['r2']:.4g}  s⁻¹")
        else:
            mid = f"double root:  r = {sol['r1']:.4g}  s⁻¹"

        if preset == "RLC_P":
            return (f"{src_line}   (true parallel resonance)\n"
                    f"α = 1/(2RC) = {sol['alpha']:.4g} s⁻¹\n"
                    f"{mid}\n"
                    f"ω₀ = 1/√(LC) = {sol['w0']:.4g} rad/s\n"
                    f"KCL:  I(t) = I_R(t) + I_L(t) + I_C(t)\n"
                    f"{settle}")
        return (f"{src_line}   (R limits current into the L∥C tank)\n"
                f"α = 1/(2RC) = {sol['alpha']:.4g} s⁻¹\n"
                f"{mid}\n"
                f"ω₀ = 1/√(LC) = {sol['w0']:.4g} rad/s      — the tank "
                f"blocks current hardest near ω₀\n"
                f"KVL:  E(t) = I(t)·R + V_t(t)\n"
                f"{settle}")

    # --------------------------------------------------------------- animation
    def _place_cursor(self, i):
        if self.family == "Series":
            self._place_cursor_series(i)
        else:
            self._place_cursor_parallel(i)

    def _place_cursor_series(self, i):
        cfg = TOPOS[self.topo]
        tm = self.tms[i - 1]
        for ln in (self.cur_q, self.cur_i, self.cur_v):
            ln.set_xdata([tm, tm])
        self.dot_q.set_data([tm], [self.Q[i - 1] * 1e3])
        self.dot_i.set_data([tm], [self.I[i - 1]])
        self.txt_t.set_text(f"t = {tm:.2f} ms\n"
                            f"Q = {self.Q[i - 1] * 1e3:+.4f} mC\n"
                            f"I = {self.I[i - 1]:+.4f} A")
        parts = [f"E {self.E[i - 1]:+.1f}"]
        if cfg["R"]:
            parts.append(f"V_R {self.VR[i - 1]:+.1f}")
        if cfg["L"]:
            parts.append(f"V_L {self.VL[i - 1]:+.1f}")
        if cfg["C"]:
            parts.append(f"V_C {self.VC[i - 1]:+.1f}")
        self.txt_vt.set_text("   ".join(parts) + "   V")
        for art in (self.cur_q, self.cur_i, self.cur_v, self.dot_q,
                    self.dot_i, self.txt_t, self.txt_vt):
            art.set_visible(True)

    def _place_cursor_parallel(self, i):
        sol = self.psol
        tm = self.tms[i - 1]
        for ln in (self.cur_q, self.cur_i, self.cur_v):
            ln.set_xdata([tm, tm])
        self.dot_q.set_data([tm], [sol["V"][i - 1]])
        self.dot_i.set_data([tm], [sol["I_total"][i - 1]])
        lines = [f"t = {tm:.2f} ms", f"V = {sol['V'][i - 1]:+.4f} V"]
        if sol["I_R"] is not None:
            lines.append(f"I_R = {sol['I_R'][i - 1]:+.4f} A")
        if sol["I_L"] is not None:
            lines.append(f"I_L = {sol['I_L'][i - 1]:+.4f} A")
        if sol["I_C"] is not None:
            lines.append(f"I_C = {sol['I_C'][i - 1]:+.4f} A")
        lines.append(f"I_total = {sol['I_total'][i - 1]:+.4f} A")
        self.txt_p.set_text("\n".join(lines))
        for art in (self.cur_q, self.cur_i, self.cur_v, self.dot_q,
                    self.dot_i, self.txt_p):
            art.set_visible(True)

    def _hide_cursor(self):
        for art in (self.cur_q, self.cur_i, self.cur_v, self.dot_q,
                    self.dot_i, self.txt_t, self.txt_vt, self.txt_p):
            art.set_visible(False)

    def _draw_upto(self, idx, final=False):
        if self.family == "Series":
            self._draw_upto_series(idx, final)
        else:
            self._draw_upto_parallel(idx, final)

    def _draw_upto_series(self, idx, final=False):
        i = max(1, min(idx, self.npts))
        self.ln_q.set_data(self.tms[:i], self.Q[:i] * 1e3)
        self.ln_i.set_data(self.tms[:i], self.I[:i])
        self.ln_e.set_data(self.tms[:i], self.E[:i])
        for ln, arr in ((self.ln_vr, self.VR), (self.ln_vl, self.VL),
                        (self.ln_vc, self.VC)):
            if ln.get_visible():
                ln.set_data(self.tms[:i], arr[:i])
        if self.ln_qss.get_visible():
            self.ln_qss.set_data(self.tms[:i], self.Qss[:i] * 1e3)
            self.ln_iss.set_data(self.tms[:i], self.Iss[:i])
        if self.ln_qrk.get_visible() and self.rk4_data is not None:
            Qr, Ir = self.rk4_data
            self.ln_qrk.set_data(self.tms[:i], Qr[:i] * 1e3)
            self.ln_irk.set_data(self.tms[:i], Ir[:i])
        if final:
            self._hide_cursor()
        else:
            self._place_cursor(i)
        self.fig.canvas.draw_idle()

    def _draw_upto_parallel(self, idx, final=False):
        i = max(1, min(idx, self.npts))
        sol = self.psol
        self.ln_pv.set_data(self.tms[:i], sol["V"][:i])
        if self.ln_pe.get_visible():
            self.ln_pe.set_data(self.tms[:i], sol["Src"][:i])
        for ln, arr in ((self.ln_pir, sol["I_R"]), (self.ln_pil, sol["I_L"]),
                        (self.ln_pic, sol["I_C"])):
            if ln.get_visible() and arr is not None:
                ln.set_data(self.tms[:i], arr[:i])
        self.ln_pitot.set_data(self.tms[:i], sol["I_total"][:i])
        if self.ln_pqc.get_visible():
            self.ln_pqc.set_data(self.tms[:i], sol["Q_C"][:i] * 1e3)
        if final:
            self._hide_cursor()
        else:
            self._place_cursor(i)
        self.fig.canvas.draw_idle()

    def _tick(self):
        self.idx = min(self.idx + self.step * self.speed, self.npts)
        if self.idx >= self.npts:
            self._stop_anim()
            self._draw_upto(self.npts, final=True)
        else:
            self._draw_upto(self.idx)

    def _stop_anim(self):
        if self.playing:
            self.timer.stop()
            self.playing = False
            self.btn_play.label.set_text("▶  Play")

    def toggle_play(self):
        if self.playing:
            self.timer.stop()
            self.playing = False
            self.btn_play.label.set_text("▶  Play")
        else:
            if self.idx >= self.npts:
                self.idx = 0
            self.playing = True
            self.btn_play.label.set_text("❚❚  Pause")
            self.timer.start()
        self.fig.canvas.draw_idle()

    # ----------------------------------------------------------- hover readout
    def _on_move(self, ev):
        if self.playing:
            return
        if ev.inaxes not in (self.ax_q, self.ax_i, self.ax_v) \
                or ev.xdata is None:
            self._on_leave(ev)
            return
        now = time.monotonic()
        if now - self._last_hover < 0.04:
            return
        self._last_hover = now
        i = int(np.clip(np.searchsorted(self.tms, ev.xdata), 1,
                        self.npts - 1))
        self.hover_on = True
        self._place_cursor(i + 1)
        self.fig.canvas.draw_idle()

    def _on_leave(self, ev):
        if self.hover_on and not self.playing:
            self.hover_on = False
            self._hide_cursor()
            self.fig.canvas.draw_idle()
