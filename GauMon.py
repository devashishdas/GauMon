#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gaussian Monitor (PyQt5)

"""

import os, sys, re, time, json, signal
from pathlib import Path
from typing import Optional, List, Tuple, Dict

import numpy as np
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QTabWidget, QFileDialog, QHBoxLayout, QVBoxLayout,
    QPushButton, QLabel, QLineEdit, QTextEdit, QTableWidget, QTableWidgetItem,
    QGroupBox, QFormLayout, QMessageBox, QSplitter, QSizePolicy, QDoubleSpinBox
)

# ---------- Defaults (user can change in GUI) ----------
SIZE_THRESHOLD_DEFAULT = 8.0      # now "Max displacement" threshold
RMSD_THRESHOLD_DEFAULT = 2.0

# ------------------------ pyqtgraph ------------------------
try:
    import pyqtgraph as pg
    from pyqtgraph import DateAxisItem, TextItem
    pg.setConfigOptions(antialias=True)
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')
    PEN_LINE = pg.mkPen(width=2)
    PEN_DASH_RED = pg.mkPen((200, 0, 0), width=2, style=Qt.DashLine)
    PEN_DASH_BLUE = pg.mkPen((0, 70, 200), width=2, style=Qt.DashLine)
except Exception:
    pg = None
    DateAxisItem = None
    TextItem = None
    PEN_LINE = PEN_DASH_RED = PEN_DASH_BLUE = None

# ------------------------ Link map (optional json override) ------------------------
def load_link_map(json_path: Path) -> Dict[int, str]:
    try:
        data = json.loads(json_path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    out: Dict[int, str] = {}
    for k, v in data.items():
        try: out[int(str(k).upper().lstrip("L"))] = str(v)
        except Exception: continue
    return out

_json_file = Path(__file__).with_name("gaussian_links.json")
LINK_MAP: Dict[int, str] = load_link_map(_json_file) or {
    101: "Input parser & molecule setup",
    103: "Berny geometry optimizer",
    202: "Orientation / symmetry / stoichiometry",
    301: "Basis set initialization",
    302: "Integral/overlap preparation",
    303: "DFT grid / XC setup",
    401: "Initial guess & spin handling",
    502: "SCF/Kohn–Sham at current geometry",
    601: "Energy gradient (forces)",
    716: "Vibrational analysis",
    721: "Thermochemistry",
}

# ------------------------ Regex ------------------------
ENTER_LINK_RE = re.compile(r"\(Enter .*?[\\/](?:g16|G16)[\\/]+l(\d+)\.exe\)")
STEP_RE = re.compile(r"\bStep number\s+(\d+)\s+out of a maximum of\s+(\d+)")
RMS_FORCE_RE = re.compile(r"RMS Force\s*=\s*([0-9.D+\-]+)")
SCF_DONE_RE = re.compile(
    r"SCF Done:\s+E\([RU]?[A-Z0-9]+\)\s*=\s*([\-\d\.D\+]+)\s+A\.U\.\s+after\s+(\d+)\s+cycles?",
    re.IGNORECASE,
)
FREQ_RE = re.compile(r"Frequencies --\s+")
NORMAL_TERM_RE = re.compile(r"Normal termination of Gaussian 16", re.IGNORECASE)
ERROR_TERM_RE = re.compile(r"Error termination", re.IGNORECASE)
CHARGE_MULT_RE = re.compile(r"Charge\s*=\s*([-+]?\d+)\s+Multiplicity\s*=\s*(\d+)")

# match the exact Gaussian coordinate table header line
COLUMN_HEADER_RE = re.compile(r"^\s*Number\s+Number\s+Type\s+X\s+Y\s+Z", re.MULTILINE)
HEADERS = ("Input orientation:", "Standard orientation:", "Z-Matrix orientation:", "Symbolic Z-matrix:")

# ------------------------ Small utilities ------------------------
def pretty_secs(s: float) -> str:
    s = int(max(0, s)); h = s // 3600; m = (s % 3600) // 60; x = s % 60
    return f"{h:02d}:{m:02d}:{x:02d}"

def one_line_readonly() -> QLineEdit:
    w = QLineEdit(); w.setReadOnly(True); w.setCursorPosition(0)
    w.setMinimumHeight(24); w.setMaximumHeight(24)
    w.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
    return w

def make_table(cols: List[str]) -> QTableWidget:
    tbl = QTableWidget(0, len(cols))
    tbl.setHorizontalHeaderLabels(cols)
    tbl.horizontalHeader().setStretchLastSection(True)
    return tbl

# ------------------------ Data container ------------------------
class LogStatus:
    def __init__(self):
        self.last_enter_link: Optional[int] = None
        self.step: Optional[int] = None
        self.step_max: Optional[int] = None
        self.rms_force: Optional[str] = None
        self.last_scf_e: Optional[float] = None
        self.last_scf_cycles: Optional[int] = None
        self.freq_blocks: int = 0
        self.termination: Optional[str] = None
        self.charge: Optional[int] = None
        self.mult: Optional[int] = None
        self.tail_text: str = ""
        self.full_text: str = ""

# ------------------------ Geometry parsing & metrics ------------------------
def _parse_orientation_table(lines: List[str], i: int) -> Tuple[List[Tuple[float,float,float]], int]:
    i += 1
    while i < len(lines) and not lines[i].strip().startswith("Center"):
        i += 1
    i += 2  # header + dashed line
    coords: List[Tuple[float,float,float]] = []
    while i < len(lines):
        s = lines[i].strip()
        if s == '' or s.startswith('-----'): break
        parts = lines[i].split()
        if len(parts) >= 6 and parts[0].isdigit():
            try:
                x, y, z = map(float, parts[-3:])
                coords.append((x, y, z))
            except Exception:
                pass
        i += 1
    return coords, i

def _parse_symbolic_zmatrix(lines: List[str], i: int) -> Tuple[List[Tuple[float,float,float]], int]:
    coords: List[Tuple[float,float,float]] = []
    i += 1
    while i < len(lines):
        s = lines[i].strip()
        if s == '': break
        got = []
        for tok in reversed(s.replace(',', ' ').split()):
            try:
                got.append(float(tok))
                if len(got) == 3: break
            except Exception:
                continue
        if len(got) == 3:
            z, y, x = got  # reversed order
            coords.append((x, y, z))
        i += 1
    return coords, i

def extract_geometries(text: str) -> List[List[Tuple[float,float,float]]]:
    lines = text.splitlines()
    blocks: List[List[Tuple[float,float,float]]] = []
    i = 0
    while i < len(lines):
        head = lines[i].strip()
        if head.startswith(HEADERS[0]) or head.startswith(HEADERS[1]) or head.startswith(HEADERS[2]):
            coords, i = _parse_orientation_table(lines, i)
            if coords: blocks.append(coords)
        elif head.startswith(HEADERS[3]):
            coords, i = _parse_symbolic_zmatrix(lines, i)
            if coords: blocks.append(coords)
        i += 1
    return blocks

def extract_geometries_by_header(text: str) -> List[List[Tuple[float,float,float]]]:
    lines = text.splitlines()
    blocks: List[List[Tuple[float,float,float]]] = []
    i = 0
    while i < len(lines):
        if COLUMN_HEADER_RE.search(lines[i]):
            # print(f"[DEBUG -vvv] Found coordinate header at line {i}")
            i += 1
            while i < len(lines) and lines[i].strip().startswith('-'):
                i += 1
            coords: List[Tuple[float,float,float]] = []
            while i < len(lines):
                s = lines[i].strip()
                if s == '' or s.startswith('-----'): break
                parts = lines[i].split()
                if len(parts) >= 6 and parts[0].isdigit():
                    try:
                        x, y, z = map(float, parts[-3:])
                        coords.append((x, y, z))
                    except Exception:
                        pass
                i += 1
            if coords:
                # print(f"[DEBUG -vvv] Parsed {len(coords)} atoms")
                blocks.append(coords)
        i += 1
    return blocks

def kabsch_rmsd(P: List[Tuple[float,float,float]], Q: List[Tuple[float,float,float]]) -> float:
    Pm = np.asarray(P, float); Qm = np.asarray(Q, float)
    if Pm.shape != Qm.shape or Pm.size == 0: return float('nan')
    Pc = Pm - Pm.mean(0); Qc = Qm - Qm.mean(0)
    V, _, Wt = np.linalg.svd(Pc.T @ Qc)
    if (np.linalg.det(V) * np.linalg.det(Wt)) < 0.0: V[:, -1] *= -1
    U = V @ Wt
    diff = (Pc @ U) - Qc
    return float(np.sqrt((diff * diff).sum() / len(P)))

def compute_max_displacement(ref: List[Tuple[float, float, float]],
                             coords: List[Tuple[float, float, float]]) -> float:
    """Maximum atom displacement relative to the reference geometry."""
    ref_arr = np.asarray(ref, float)
    arr = np.asarray(coords, float)
    if ref_arr.shape != arr.shape or len(arr) < 1:
        return float('nan')
    disp = np.linalg.norm(arr - ref_arr, axis=1)
    return float(disp.max())

# ------------------------ Parsing one log pass ------------------------
def parse_log_file(path: Path, tail_lines: int = 4000) -> LogStatus:
    st = LogStatus()
    if not path.exists(): return st
    try:
        text = path.read_text(errors="ignore")
    except Exception:
        return st
    st.full_text = text
    lines = text.splitlines()
    tail = lines[-tail_lines:]
    st.tail_text = "\n".join(lines[-800:])

    for ln in reversed(tail):
        if NORMAL_TERM_RE.search(ln): st.termination = 'normal'; break
        if ERROR_TERM_RE.search(ln):  st.termination = 'error';  break

    for ln in tail:
        m = ENTER_LINK_RE.search(ln)
        if m: st.last_enter_link = int(m.group(1))
        m = STEP_RE.search(ln)
        if m: st.step, st.step_max = int(m.group(1)), int(m.group(2))
        m = RMS_FORCE_RE.search(ln)
        if m: st.rms_force = m.group(1)
        m = SCF_DONE_RE.search(ln)
        if m:
            try:
                st.last_scf_e = float(m.group(1).replace('D','E'))
                st.last_scf_cycles = int(m.group(2))
            except Exception:
                pass
        if FREQ_RE.search(ln): st.freq_blocks += 1
        m = CHARGE_MULT_RE.search(ln)
        if m:
            st.charge = int(m.group(1)); st.mult = int(m.group(2))
    return st

# ======================== SHARED MONITOR (logic + plots) ========================
class BaseMonitor(QWidget):
    """Holds state, shared polling, plotting, and table refresh for both tabs."""
    def __init__(self, parent=None):
        super().__init__(parent)
        # file / growth
        self.current_log: Optional[Path] = None
        self._last_size = 0
        self._last_growth_ts = time.time()

        # histories
        self.scf_history: List[Tuple[float,int,float]] = []  # (E, cycles, ts)
        self.step_rms: Dict[int, float] = {}
        self.step_times: Dict[int, float] = {}
        self.link_history: List[Tuple[float,int,str]] = []
        self._last_link_seen: Optional[int] = None

        # geometry tracking
        self.geom_disp: List[float] = []   # max displacement vs initial
        self.rmsd_series: List[float] = []
        self.prev_rmsd_series: List[float] = []
        self._geom_count_seen: int = 0
        self._first_geom: Optional[List[Tuple[float,float,float]]] = None
        self._last_geom_for_prev: Optional[List[Tuple[float,float,float]]] = None
        self.rmsd_current: Optional[float] = None

        # thresholds / alarms
        self.size_alarm_threshold: float = SIZE_THRESHOLD_DEFAULT
        self.rmsd_alarm_threshold: float = RMSD_THRESHOLD_DEFAULT
        self.alarm_active: bool = False
        self._banner_on: bool = False

        # opts
        self.max_steps_default = 150
        self._last_detected_max: Optional[int] = None

        # timers
        self._epoch_for_timer: Optional[float] = None
        self._last_scf_count_seen: int = 0

        # --- UI placeholders (wired by subclasses) ---
        self.alarm_banner: QLabel = QLabel()
        self.timer_lbl: QLabel = QLabel()
        self.route_line = self.cm_line = self.link_line = one_line_readonly()
        self.opt_line = self.scf_line = self.freq_line = one_line_readonly()
        self.disp_line = self.rmsd_line = self.eta_line = one_line_readonly()
        self.energy_table: QTableWidget = make_table(['#', 'E (a.u.)', 'Cycles'])
        self.energy_header: QLabel = QLabel('All SCFs (0)')
        self.history_table: QTableWidget = make_table(['Time', 'Link', 'Description'])
        self.tail_view: QTextEdit = QTextEdit()
        self.plot_scf = self.curve_scf = None
        self.plot_scf_time = self.curve_scf_time = None
        self.plot_rms = self.curve_rms = None
        self.plot_disp = self.curve_disp = None
        self.disp_threshold_line = self.disp_overlay = None
        self.plot_rmsd = self.curve_rmsd_plot = self.curve_rmsd_prev = None
        self.rmsd_threshold_line = self.rmsd_overlay = None

        # controls for thresholds (wired in tabs)
        self.spin_disp_thr: Optional[QDoubleSpinBox] = None
        self.spin_rmsd_thr: Optional[QDoubleSpinBox] = None

        # timers (constructed here so both tabs share logic)
        self.poll_timer = QTimer(self); self.poll_timer.setInterval(800); self.poll_timer.timeout.connect(self.poll)
        self.clock_timer = QTimer(self); self.clock_timer.setInterval(1000); self.clock_timer.timeout.connect(self.update_timer)
        self.alarm_timer = QTimer(self); self.alarm_timer.setInterval(450); self.alarm_timer.timeout.connect(self._tick_alarm_blink)

    # ---------- Shared parsing helpers ----------
    def _update_scf_history(self, log_path: Path) -> bool:
        try:
            lines = log_path.read_text(errors='ignore').splitlines()
        except Exception:
            return False
        scfs: List[Tuple[float,int]] = []
        for ln in lines:
            m = SCF_DONE_RE.search(ln)
            if m:
                try: scfs.append((float(m.group(1).replace('D','E')), int(m.group(2))))
                except Exception: pass
        appended = False
        if len(scfs) > len(self.scf_history):
            now = time.time()
            for e, cyc in scfs[len(self.scf_history):]:
                self.scf_history.append((e, cyc, now))
                appended = True
        return appended

    def _update_step_stats(self, st: LogStatus):
        if st.step is not None:
            if st.step not in self.step_times: self.step_times[st.step] = time.time()
            if st.rms_force:
                try: self.step_rms[st.step] = float(st.rms_force.replace('D','E'))
                except Exception: pass
        if st.step_max: self._last_detected_max = st.step_max

    def _append_history_row(self, ts: float, link_id: int, desc: str):
        r = self.history_table.rowCount()
        self.history_table.insertRow(r)
        self.history_table.setItem(r, 0, QTableWidgetItem(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(ts))))
        self.history_table.setItem(r, 1, QTableWidgetItem(f"l{link_id:03d}.exe"))
        self.history_table.setItem(r, 2, QTableWidgetItem(desc))

    def _update_link_history(self, st: LogStatus):
        if st.last_enter_link and st.last_enter_link != self._last_link_seen:
            self._last_link_seen = st.last_enter_link
            desc = LINK_MAP.get(st.last_enter_link, 'Unknown')
            self.link_history.append((time.time(), st.last_enter_link, desc))
            self._append_history_row(time.time(), st.last_enter_link, desc)

    def _update_geometry_metrics(self, st: LogStatus, force: bool = False):
        # Always scan the full log for geometries
        # print(f"[DEBUG -vvv] Geometry scan triggered (force={force}, link={st.last_enter_link})")

        # Decide if we should scan at all
        should_scan = (
            force or (st.last_enter_link == 202)
            or COLUMN_HEADER_RE.search(st.tail_text)
            or any(h in st.tail_text for h in HEADERS)
        )
        if not should_scan:
            return

        # Parse geometries from full text
        blocks = extract_geometries_by_header(st.full_text) or extract_geometries(st.full_text)

        if not blocks:
            if force:
                print("[ERROR] New energy detected but NO coordinate blocks found in log!")
            return

        # Initialize reference geometry
        if self._first_geom is None:
            self._first_geom = blocks[0]

        # Process only new blocks
        for b in blocks[self._geom_count_seen:]:
            try:
                disp = compute_max_displacement(self._first_geom, b) if self._first_geom else float('nan')
            except Exception:
                disp = float('nan')
            self.geom_disp.append(disp)

            # RMSD vs first
            if self._first_geom is not None and len(b) == len(self._first_geom):
                rmsd = kabsch_rmsd(self._first_geom, b)
            else:
                rmsd = float('nan')
            self.rmsd_series.append(rmsd)

            # RMSD vs previous
            if self._last_geom_for_prev is not None and len(b) == len(self._last_geom_for_prev):
                prev_r = kabsch_rmsd(self._last_geom_for_prev, b)
            else:
                prev_r = 0.0
            self.prev_rmsd_series.append(prev_r)

            self._last_geom_for_prev = b

        # Update current RMSD
        last = None
        if self._first_geom:
            for b in reversed(blocks):
                if len(b) == len(self._first_geom):
                    last = b
                    break
        self.rmsd_current = (
            kabsch_rmsd(self._first_geom, last) if (last and self._first_geom) else None
        )

        self._geom_count_seen = len(blocks)

    def _compute_eta_text(self) -> str:
        if len(self.step_times) < 2:
            return f"Max time (≤{self._last_detected_max or self.max_steps_default}): —"
        steps = sorted(self.step_times.keys())
        t0, tN = self.step_times[steps[0]], self.step_times[steps[-1]]
        avg = max(1e-6, (tN - t0) / max(1, steps[-1] - steps[0]))
        remaining = max(0, (self._last_detected_max or self.max_steps_default) - steps[-1])
        return f"Max time (≤{self._last_detected_max or self.max_steps_default}): ~{pretty_secs(remaining*avg)}  (avg {avg:.1f}s/step)"

    # ---------- Shared timers / UI updates ----------
    def start_monitors(self):
        self._epoch_for_timer = time.time()
        self._last_size = 0; self._last_growth_ts = time.time()
        self.clock_timer.start(); self.poll_timer.start()

    def stop_monitors(self):
        self.clock_timer.stop(); self.alarm_timer.stop(); self.poll_timer.stop()

    def update_timer(self):
        t0 = self._epoch_for_timer or time.time()
        self.timer_lbl.setText(f"⏱ {pretty_secs(time.time() - t0)}")

    def _tick_alarm_blink(self):
        self._banner_on = not self._banner_on
        self.alarm_banner.setStyleSheet(
            "padding:6px; border-radius:6px; font-weight:800;"
            + ("color:white; background:#c21807;" if self._banner_on else "color:black; background:#ffd54f;")
        )

    # ---------- The single shared poll ----------
    def poll(self):
        if not self.current_log: return
        st = parse_log_file(self.current_log)

        try:
            size = self.current_log.stat().st_size
            if size > self._last_size:
                self._last_size = size; self._last_growth_ts = time.time()
        except Exception:
            pass

        new_energy = self._update_scf_history(self.current_log)
        self._update_step_stats(st)
        self._update_link_history(st)
        self._update_geometry_metrics(st, force=new_energy)

        if st.last_enter_link:
            self.link_line.setText(f"l{st.last_enter_link:03d}.exe — {LINK_MAP.get(st.last_enter_link,'Unknown')}")
        self.opt_line.setText(
            f"step {st.step}/{st.step_max or '—'} | RMS Force: {st.rms_force or '—'}" if st.step is not None else '—'
        )
        self.scf_line.setText(
            f"E={st.last_scf_e:.9f} a.u. | cycles={st.last_scf_cycles}" if st.last_scf_e is not None else '—'
        )
        self.freq_line.setText(f"blocks={st.freq_blocks}")
        if st.charge is not None and st.mult is not None:
            self.cm_line.setText(f"{st.charge}/{st.mult}")
        self.eta_line.setText(self._compute_eta_text())

        if self.rmsd_current is not None and np.isfinite(self.rmsd_current):
            self.rmsd_line.setText(f"{self.rmsd_current:.3f} Å")
        else:
            self.rmsd_line.setText('—')

        latest_disp = self.geom_disp[-1] if self.geom_disp else None
        if latest_disp is not None and np.isfinite(latest_disp):
            self.disp_line.setText(f"{latest_disp:.2f} Å (thr={self.size_alarm_threshold:.1f} Å)")
            if latest_disp > self.size_alarm_threshold:
                self.disp_line.setStyleSheet("color: #b00020; font-weight: 700;")
                if not self.alarm_active:
                    self.alarm_active = True; self.alarm_banner.setVisible(True); self.alarm_timer.start()
            else:
                self.disp_line.setStyleSheet("")
                if self.alarm_active:
                    self.alarm_active = False; self.alarm_banner.setVisible(False); self.alarm_timer.stop()
        else:
            self.disp_line.setText('—'); self.disp_line.setStyleSheet("")
            self.alarm_active = False; self.alarm_banner.setVisible(False); self.alarm_timer.stop()

        self._refresh_energy_table()
        self._refresh_plots()
        self.tail_view.setPlainText(st.tail_text)

    # ---------- Table + plots ----------
    def _refresh_energy_table(self):
        total = len(self.scf_history)
        self.energy_header.setText(f'All SCFs ({total})')
        self.energy_table.setRowCount(total)
        for i, (e, c, _) in enumerate(self.scf_history, start=1):
            self.energy_table.setItem(i-1, 0, QTableWidgetItem(str(i)))
            self.energy_table.setItem(i-1, 1, QTableWidgetItem(f"{e:.9f}"))
            self.energy_table.setItem(i-1, 2, QTableWidgetItem(str(c)))
        if total: self.energy_table.scrollToBottom()

    def _snap_small(self, arr: List[float], tol: float = 1e-2) -> List[float]:
        out = []
        for v in arr:
            if not np.isfinite(v): out.append(v)
            else: out.append(0.0 if abs(v) < tol else v)
        return out

    def _finite_xy(self, x, y):
        x = np.asarray(x, float); y = np.asarray(y, float)
        m = np.isfinite(x) & np.isfinite(y)
        return (x[m], y[m]) if np.any(m) else (np.array([]), np.array([]))

    def _plot_set(self, plot_widget, curve, x, y, x_default=(0,1), y_default=(0,1)):
        if plot_widget is None or curve is None: return False
        x2, y2 = self._finite_xy(x, y)
        curve.setData(x2, y2)
        try:
            if x2.size == 0: plot_widget.setRange(xRange=x_default, yRange=y_default, padding=0.05)
            else: plot_widget.enableAutoRange(x=True, y=True)
        except Exception:
            pass
        return x2.size > 0

    def _refresh_plots(self):
        if pg is None: return

        # SCF vs index / time
        xs = list(range(len(self.scf_history)))
        ys = [e for (e,_,__) in self.scf_history]
        self._plot_set(self.plot_scf, self.curve_scf, xs, ys, x_default=(0,1), y_default=(-1,1))

        ts = [t for (_,_,t) in self.scf_history]
        self._plot_set(self.plot_scf_time, self.curve_scf_time, ts, ys,
                       x_default=(time.time()-1.0, time.time()+1.0), y_default=(-1,1))

        # RMS Force vs step
        if self.step_rms:
            steps = sorted(self.step_rms.keys())
            rvals = [self.step_rms[s] for s in steps]
            self._plot_set(self.plot_rms, self.curve_rms, steps, rvals)
        else:
            self._plot_set(self.plot_rms, self.curve_rms, [], [], x_default=(0,1), y_default=(0,1))

        # Max displacement vs geometry index
        gi = list(range(len(self.geom_disp)))
        disp_vals = self._snap_small(self.geom_disp, tol=1e-2)
        shown_disp = self._plot_set(self.plot_disp, self.curve_disp, gi, disp_vals)
        if shown_disp and disp_vals:
            last_x = gi[-1]; last_y = disp_vals[-1]
            if np.isfinite(last_y) and last_y > self.size_alarm_threshold:
                self.disp_overlay.setText('💥 DRIFT!'); self.disp_overlay.setPos(last_x, float(last_y))
            else:
                self.disp_overlay.setText('')
        else:
            self.disp_overlay.setText('')

        # RMSD vs first + vs previous
        gi_r = list(range(len(self.rmsd_series)))
        rmsd_vals = self._snap_small(self.rmsd_series, tol=1e-2)
        shown_r = self._plot_set(self.plot_rmsd, self.curve_rmsd_plot, gi_r, rmsd_vals)

        gi_prev = list(range(len(self.prev_rmsd_series)))
        rmsd_prev_vals = self._snap_small(self.prev_rmsd_series, tol=1e-2)
        self._plot_set(self.plot_rmsd, self.curve_rmsd_prev, gi_prev, rmsd_prev_vals)

        if shown_r and rmsd_vals:
            last_r = rmsd_vals[-1]
            if np.isfinite(last_r) and last_r > self.rmsd_alarm_threshold:
                self.rmsd_overlay.setText('⚠️ RMSD HIGH'); self.rmsd_overlay.setPos(gi_r[-1], float(last_r))
            else:
                self.rmsd_overlay.setText('')
        else:
            self.rmsd_overlay.setText('')

    # ---------- threshold wiring ----------
    def _init_threshold_controls(self, form: QFormLayout):
        # Max displacement threshold
        self.spin_disp_thr = QDoubleSpinBox()
        self.spin_disp_thr.setRange(0.0, 1e6)
        self.spin_disp_thr.setDecimals(3)
        self.spin_disp_thr.setSingleStep(0.1)
        self.spin_disp_thr.setValue(self.size_alarm_threshold)
        self.spin_disp_thr.valueChanged.connect(self._on_disp_thr_changed)

        # RMSD threshold
        self.spin_rmsd_thr = QDoubleSpinBox()
        self.spin_rmsd_thr.setRange(0.0, 1e6)
        self.spin_rmsd_thr.setDecimals(3)
        self.spin_rmsd_thr.setSingleStep(0.1)
        self.spin_rmsd_thr.setValue(self.rmsd_alarm_threshold)
        self.spin_rmsd_thr.valueChanged.connect(self._on_rmsd_thr_changed)

        form.addRow('Max displacement thr (Å):', self.spin_disp_thr)
        form.addRow('RMSD threshold (Å):', self.spin_rmsd_thr)

    def _on_disp_thr_changed(self, v: float):
        self.size_alarm_threshold = float(v)
        if self.disp_threshold_line is not None:
            try: self.disp_threshold_line.setPos(self.size_alarm_threshold)
            except Exception: pass

    def _on_rmsd_thr_changed(self, v: float):
        self.rmsd_alarm_threshold = float(v)
        if self.rmsd_threshold_line is not None:
            try: self.rmsd_threshold_line.setPos(self.rmsd_alarm_threshold)
            except Exception: pass

# ======================== TAB: RUN GAUSSIAN ========================
class RunGaussianTab(BaseMonitor):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._proc_pid: Optional[int] = None
        self._proc_pgid: Optional[int] = None

        # --- controls ---
        self.exe_edit = QLineEdit('g16')
        self.gjf_edit = QLineEdit()
        self.gjf_btn = QPushButton('Browse .gjf…')
        self.work_edit = QLineEdit('')
        self.log_edit = QLineEdit('')
        self.start_btn = QPushButton('Start')
        self.stop_btn = QPushButton('Stop'); self.stop_btn.setEnabled(False)
        self.kill_btn = QPushButton('Kill lXXX.exe')

        # --- left panel (Run Info / Energy / History) ---
        self.alarm_banner = QLabel('🚨 DRIFT DETECTED'); self.alarm_banner.setVisible(False)
        self.alarm_banner.setAlignment(Qt.AlignCenter)
        self.alarm_banner.setStyleSheet("padding:6px; border-radius:6px; font-weight:800;")
        self.timer_lbl = QLabel('⏱ 00:00:00'); self.timer_lbl.setStyleSheet("font-weight:600;")
        self.route_line = one_line_readonly()
        self.cm_line = one_line_readonly()
        self.link_line = one_line_readonly()
        self.opt_line = one_line_readonly()
        self.scf_line = one_line_readonly()
        self.freq_line = one_line_readonly()
        self.disp_line = one_line_readonly()
        self.rmsd_line = one_line_readonly()
        self.eta_line = one_line_readonly()

        runinfo_form = QFormLayout()
        runinfo_form.addRow('Timer:', self.timer_lbl)
        runinfo_form.addRow('Route:', self.route_line)
        runinfo_form.addRow('Q/M:', self.cm_line)
        runinfo_form.addRow('Link:', self.link_line)
        runinfo_form.addRow('Opt:', self.opt_line)
        runinfo_form.addRow('SCF:', self.scf_line)
        runinfo_form.addRow('Freq:', self.freq_line)
        runinfo_form.addRow('Max displacement:', self.disp_line)
        runinfo_form.addRow('RMSD:', self.rmsd_line)
        runinfo_form.addRow('ETA:', self.eta_line)
        # thresholds (user adjustable)
        self._init_threshold_controls(runinfo_form)

        runinfo_wrap = QWidget(); v = QVBoxLayout(); v.addWidget(self.alarm_banner)
        ri = QWidget(); ri.setLayout(runinfo_form); v.addWidget(ri); runinfo_wrap.setLayout(v)

        energy = QWidget(); self.energy_table = make_table(['#', 'E (a.u.)', 'Cycles'])
        self.energy_header = QLabel('All SCFs (0)')
        ev = QVBoxLayout(); ev.addWidget(self.energy_header); ev.addWidget(self.energy_table); energy.setLayout(ev)

        history = QWidget(); self.history_table = make_table(['Time', 'Link', 'Description'])
        hv = QVBoxLayout(); hv.addWidget(self.history_table); history.setLayout(hv)

        left_tabs = QTabWidget()
        left_tabs.addTab(runinfo_wrap, 'Run Info')
        left_tabs.addTab(energy, 'Energy')
        left_tabs.addTab(history, 'History')

        # --- right plots ---
        self._build_plots()

        # --- top (run controls) ---
        top = QGroupBox('Run Gaussian')
        form = QFormLayout()
        form.addRow('Gaussian exec:', self.exe_edit)
        gjf_row = QHBoxLayout(); gjf_row.addWidget(self.gjf_edit); gjf_row.addWidget(self.gjf_btn)
        gjf_holder = QWidget(); gjf_holder.setLayout(gjf_row)
        form.addRow('Input .gjf:', gjf_holder)
        form.addRow('Working dir (optional):', self.work_edit)
        form.addRow('Log file name (auto if blank):', self.log_edit)
        btn_row = QHBoxLayout()
        btn_row.addWidget(self.start_btn); btn_row.addWidget(self.stop_btn); btn_row.addWidget(self.kill_btn)
        holder2 = QWidget(); holder2.setLayout(btn_row); form.addRow('', holder2)
        top.setLayout(form)

        # --- bottom split ---
        bottom = QSplitter(Qt.Horizontal)
        left_box = QGroupBox('Info / Energy / History'); lv = QVBoxLayout(); lv.addWidget(left_tabs); left_box.setLayout(lv)
        right_box = QGroupBox('Plots'); rv = QVBoxLayout(); rv.addWidget(self.graph_tabs); right_box.setLayout(rv)
        bottom.addWidget(left_box); bottom.addWidget(right_box); bottom.setSizes([520, 700])

        # --- log tab ---
        self.tail_view = QTextEdit(); self.tail_view.setReadOnly(True)
        log_tab = QWidget(); log_layout = QVBoxLayout(); log_layout.addWidget(self.tail_view); log_tab.setLayout(log_layout)

        # --- root tabs ---
        dash = QWidget(); dash_layout = QVBoxLayout(); dash_layout.addWidget(top); dash_layout.addWidget(bottom); dash.setLayout(dash_layout)
        self.inner_tabs = QTabWidget(); self.inner_tabs.addTab(dash, 'Dashboard'); self.inner_tabs.addTab(log_tab, 'Log')
        root = QVBoxLayout(); root.addWidget(self.inner_tabs); self.setLayout(root)

        # signals
        self.gjf_btn.clicked.connect(self.choose_gjf)
        self.start_btn.clicked.connect(self.start_gaussian)
        self.stop_btn.clicked.connect(self.stop_gaussian)
        self.kill_btn.clicked.connect(self.kill_links)

    def _build_plots(self):
        self.graph_tabs = QTabWidget()
        if pg is None:
            self.graph_tabs.addTab(QLabel('pyqtgraph not installed — run: pip install pyqtgraph'), 'Graphs')
            return
        self.plot_scf = pg.PlotWidget(title='SCF Energy (a.u.) vs index'); self.plot_scf.showGrid(x=True, y=True)
        self.curve_scf = self.plot_scf.plot([], [], pen=PEN_LINE, symbol='o')
        self.graph_tabs.addTab(self.plot_scf, 'SCF vs Index')

        axis = DateAxisItem(orientation='bottom') if DateAxisItem else None
        self.plot_scf_time = pg.PlotWidget(title='SCF Energy (a.u.) vs wall-time', axisItems={'bottom': axis} if axis else None)
        self.plot_scf_time.showGrid(x=True, y=True)
        self.curve_scf_time = self.plot_scf_time.plot([], [], pen=PEN_LINE, symbol='o')
        self.graph_tabs.addTab(self.plot_scf_time, 'SCF vs Time')

        self.plot_rms = pg.PlotWidget(title='RMS Force vs Optimization Step'); self.plot_rms.showGrid(x=True, y=True)
        self.curve_rms = self.plot_rms.plot([], [], pen=PEN_LINE, symbol='o')
        self.graph_tabs.addTab(self.plot_rms, 'RMS Force')

        self.plot_disp = pg.PlotWidget(title='Max atom displacement vs initial (Å)'); self.plot_disp.showGrid(x=True, y=True)
        self.curve_disp = self.plot_disp.plot([], [], pen=PEN_LINE, symbol='o')
        self.disp_threshold_line = pg.InfiniteLine(pos=self.size_alarm_threshold, angle=0, pen=PEN_DASH_RED); self.plot_disp.addItem(self.disp_threshold_line)
        self.disp_overlay = TextItem('', color=(200, 0, 0), anchor=(0.5, 1.2)); self.plot_disp.addItem(self.disp_overlay)
        self.graph_tabs.addTab(self.plot_disp, 'Displacement')

        self.plot_rmsd = pg.PlotWidget(title='RMSD vs geometry (Å)'); self.plot_rmsd.showGrid(x=True, y=True)
        self.curve_rmsd_plot = self.plot_rmsd.plot([], [], pen=PEN_LINE, symbol='o')
        self.curve_rmsd_prev = self.plot_rmsd.plot([], [], pen=PEN_DASH_BLUE, symbol='o')
        self.rmsd_threshold_line = pg.InfiniteLine(pos=self.rmsd_alarm_threshold, angle=0, pen=PEN_DASH_BLUE); self.plot_rmsd.addItem(self.rmsd_threshold_line)
        self.rmsd_overlay = TextItem('', color=(0, 70, 200), anchor=(0.5, 1.2)); self.plot_rmsd.addItem(self.rmsd_overlay)
        self.graph_tabs.addTab(self.plot_rmsd, 'RMSD')

    # ---------- launch / stop ----------
    def _spawn_gaussian(self, exe: str, gjf: Path, work: Path) -> bool:
        try:
            os.makedirs(work, exist_ok=True)
            pid = os.fork()
            if pid == 0:
                os.setsid()
                os.chdir(str(work))
                devnull = os.open(os.devnull, os.O_RDWR)
                os.dup2(devnull, 0); os.dup2(devnull, 1); os.dup2(devnull, 2)
                os.execvp(exe, [exe, str(gjf)])
                os._exit(127)
            else:
                self._proc_pid = pid
                try: self._proc_pgid = os.getpgid(pid)
                except Exception: self._proc_pgid = None
                return True
        except Exception as e:
            QMessageBox.critical(self, "Launch error", f"Failed to start Gaussian: {e}")
            return False

    def choose_gjf(self):
        start_dir = self.work_edit.text().strip() or (str(Path(self.gjf_edit.text()).parent) if self.gjf_edit.text().strip() else os.getcwd())
        path, _ = QFileDialog.getOpenFileName(self, 'Select Gaussian Input (.gjf)', start_dir,
                                              'Gaussian Input (*.gjf *.com);;All Files (*)')
        if not path: return
        p = Path(path)
        self.gjf_edit.setText(str(p))
        self.log_edit.setText(p.with_suffix('.log').name)
        if not self.work_edit.text().strip(): self.work_edit.setText(str(p.parent))
        route, ch, mu = read_route_charge_from_gjf(p)
        self.route_line.setText(route or '—')
        self.cm_line.setText(f"{ch}/{mu}" if ch is not None else '—')

    def start_gaussian(self):
        exe = self.exe_edit.text().strip() or 'g16'
        gjf = Path(self.gjf_edit.text().strip())
        if not gjf.exists():
            QMessageBox.warning(self, 'Missing input', 'Please select a valid .gjf file.'); return
        work = Path(self.work_edit.text().strip()) if self.work_edit.text().strip() else gjf.parent
        log_name = self.log_edit.text().strip() or gjf.with_suffix('.log').name
        self.current_log = work / log_name

        # reset state
        self.scf_history.clear(); self.step_rms.clear(); self.step_times.clear()
        self.link_history.clear(); self._last_link_seen = None
        self.history_table.setRowCount(0); self.energy_table.setRowCount(0)
        self._last_detected_max = None
        self.geom_disp.clear(); self.rmsd_series.clear(); self.prev_rmsd_series.clear()
        self._last_geom_for_prev = None; self._geom_count_seen = 0
        self.disp_line.setText('—'); self.disp_line.setStyleSheet("")
        self.rmsd_line.setText('—'); self._first_geom = None
        self.alarm_active = False; self.alarm_banner.setVisible(False); self.alarm_timer.stop()
        self._last_scf_count_seen = 0

        if not self._spawn_gaussian(exe, gjf, work): return
        self.start_btn.setEnabled(False); self.stop_btn.setEnabled(True)
        self.start_monitors()

    def stop_gaussian(self):
        ok = False
        try:
            if self._proc_pgid is not None:
                os.killpg(self._proc_pgid, signal.SIGTERM); time.sleep(1.0)
                os.killpg(self._proc_pgid, signal.SIGKILL); ok = True
        except ProcessLookupError:
            ok = True
        except Exception:
            pass
        try:
            if self._proc_pid is not None and not ok:
                os.system(f"pkill -TERM -P {self._proc_pid}"); time.sleep(1.0)
                os.system(f"pkill -KILL -P {self._proc_pid}"); ok = True
        except Exception:
            pass

        self.stop_btn.setEnabled(False); self.start_btn.setEnabled(True)
        self.stop_monitors()

    def kill_links(self):
        os.system(r"pkill -9 -f 'l[0-9][0-9][0-9]\.exe'")

# ======================== TAB: PARSE LOG ONLY ========================
class ParseOnlyTab(BaseMonitor):
    def __init__(self, parent=None):
        super().__init__(parent)

        # controls
        self.log_edit = QLineEdit()
        self.log_btn = QPushButton('Browse .log…')
        self.attach_btn = QPushButton('Attach')

        # left panel
        self.alarm_banner = QLabel('🚨 DRIFT DETECTED'); self.alarm_banner.setVisible(False)
        self.alarm_banner.setAlignment(Qt.AlignCenter)
        self.alarm_banner.setStyleSheet("padding:6px; border-radius:6px; font-weight:800;")
        self.timer_lbl = QLabel('⏱ 00:00:00'); self.timer_lbl.setStyleSheet("font-weight:600;")
        self.route_line = one_line_readonly()
        self.cm_line = one_line_readonly()
        self.link_line = one_line_readonly()
        self.opt_line = one_line_readonly()
        self.scf_line = one_line_readonly()
        self.freq_line = one_line_readonly()
        self.disp_line = one_line_readonly()
        self.rmsd_line = one_line_readonly()
        self.eta_line = one_line_readonly()

        runinfo_form = QFormLayout()
        runinfo_form.addRow('Timer:', self.timer_lbl)
        runinfo_form.addRow('Route:', self.route_line)
        runinfo_form.addRow('Q/M:', self.cm_line)
        runinfo_form.addRow('Link:', self.link_line)
        runinfo_form.addRow('Opt:', self.opt_line)
        runinfo_form.addRow('SCF:', self.scf_line)
        runinfo_form.addRow('Freq:', self.freq_line)
        runinfo_form.addRow('Max displacement:', self.disp_line)
        runinfo_form.addRow('RMSD:', self.rmsd_line)
        runinfo_form.addRow('ETA:', self.eta_line)
        # thresholds (user adjustable)
        self._init_threshold_controls(runinfo_form)

        runinfo_wrap = QWidget(); v = QVBoxLayout(); v.addWidget(self.alarm_banner)
        ri = QWidget(); ri.setLayout(runinfo_form); v.addWidget(ri); runinfo_wrap.setLayout(v)

        energy = QWidget(); self.energy_table = make_table(['#', 'E (a.u.)', 'Cycles'])
        self.energy_header = QLabel('All SCFs (0)')
        ev = QVBoxLayout(); ev.addWidget(self.energy_header); ev.addWidget(self.energy_table); energy.setLayout(ev)

        history = QWidget(); self.history_table = make_table(['Time', 'Link', 'Description'])
        hv = QVBoxLayout(); hv.addWidget(self.history_table); history.setLayout(hv)

        left_tabs = QTabWidget()
        left_tabs.addTab(runinfo_wrap, 'Run Info')
        left_tabs.addTab(energy, 'Energy')
        left_tabs.addTab(history, 'History')

        # plots
        self._build_plots()

        # top row (attach)
        form = QFormLayout()
        lr = QHBoxLayout(); lr.addWidget(self.log_edit); lr.addWidget(self.log_btn); lr.addWidget(self.attach_btn)
        holder = QWidget(); holder.setLayout(lr)
        form.addRow('Log file:', holder)

        # bottom
        bottom = QSplitter(Qt.Horizontal)
        left_box = QGroupBox('Info / Energy / History'); lv = QVBoxLayout(); lv.addWidget(left_tabs); left_box.setLayout(lv)
        right_box = QGroupBox('Plots'); rv = QVBoxLayout(); rv.addWidget(self.graph_tabs); right_box.setLayout(rv)
        bottom.addWidget(left_box); bottom.addWidget(right_box); bottom.setSizes([520, 700])

        # log tab
        self.tail_view = QTextEdit(); self.tail_view.setReadOnly(True)
        log_tab = QWidget(); log_layout = QVBoxLayout(); log_layout.addWidget(self.tail_view); log_tab.setLayout(log_layout)

        # root tabs
        dash = QWidget(); dash_layout = QVBoxLayout(); dash_layout.addLayout(form); dash_layout.addWidget(bottom); dash.setLayout(dash_layout)
        self.inner_tabs = QTabWidget(); self.inner_tabs.addTab(dash, 'Dashboard'); self.inner_tabs.addTab(log_tab, 'Log')
        root = QVBoxLayout(); root.addWidget(self.inner_tabs); self.setLayout(root)

        # signals
        self.log_btn.clicked.connect(self.choose_log)
        self.attach_btn.clicked.connect(self.attach)

    def _build_plots(self):
        self.graph_tabs = QTabWidget()
        if pg is None:
            self.graph_tabs.addTab(QLabel('pyqtgraph not installed — run: pip install pyqtgraph'), 'Graphs')
            return
        self.plot_scf = pg.PlotWidget(title='SCF Energy (a.u.) vs index'); self.plot_scf.showGrid(x=True, y=True)
        self.curve_scf = self.plot_scf.plot([], [], pen=PEN_LINE, symbol='o')
        self.graph_tabs.addTab(self.plot_scf, 'SCF vs Index')

        axis = DateAxisItem(orientation='bottom') if DateAxisItem else None
        self.plot_scf_time = pg.PlotWidget(title='SCF Energy (a.u.) vs wall-time', axisItems={'bottom': axis} if axis else None)
        self.plot_scf_time.showGrid(x=True, y=True)
        self.curve_scf_time = self.plot_scf_time.plot([], [], pen=PEN_LINE, symbol='o')
        self.graph_tabs.addTab(self.plot_scf_time, 'SCF vs Time')

        self.plot_rms = pg.PlotWidget(title='RMS Force vs Optimization Step'); self.plot_rms.showGrid(x=True, y=True)
        self.curve_rms = self.plot_rms.plot([], [], pen=PEN_LINE, symbol='o')
        self.graph_tabs.addTab(self.plot_rms, 'RMS Force')

        self.plot_disp = pg.PlotWidget(title='Max atom displacement vs initial (Å)'); self.plot_disp.showGrid(x=True, y=True)
        self.curve_disp = self.plot_disp.plot([], [], pen=PEN_LINE, symbol='o')
        self.disp_threshold_line = pg.InfiniteLine(pos=self.size_alarm_threshold, angle=0, pen=PEN_DASH_RED); self.plot_disp.addItem(self.disp_threshold_line)
        self.disp_overlay = TextItem('', color=(200, 0, 0), anchor=(0.5, 1.2)); self.plot_disp.addItem(self.disp_overlay)
        self.graph_tabs.addTab(self.plot_disp, 'Displacement')

        self.plot_rmsd = pg.PlotWidget(title='RMSD vs geometry (Å)'); self.plot_rmsd.showGrid(x=True, y=True)
        self.curve_rmsd_plot = self.plot_rmsd.plot([], [], pen=PEN_LINE, symbol='o')
        self.curve_rmsd_prev = self.plot_rmsd.plot([], [], pen=PEN_DASH_BLUE, symbol='o')
        self.rmsd_threshold_line = pg.InfiniteLine(pos=self.rmsd_alarm_threshold, angle=0, pen=PEN_DASH_BLUE); self.plot_rmsd.addItem(self.rmsd_threshold_line)
        self.rmsd_overlay = TextItem('', color=(0, 70, 200), anchor=(0.5, 1.2)); self.plot_rmsd.addItem(self.rmsd_overlay)
        self.graph_tabs.addTab(self.plot_rmsd, 'RMSD')

    # attach / choose
    def choose_log(self):
        start_dir = str(Path(self.log_edit.text().strip()).parent) if self.log_edit.text().strip() else os.getcwd()
        path, _ = QFileDialog.getOpenFileName(self, 'Select Gaussian Log (.log)', start_dir,
                                              'Gaussian Log (*.log *.out *.log.txt);;All Files (*)')
        if path: self.log_edit.setText(path)

    def attach(self):
        p = Path(self.log_edit.text().strip())
        if not p.exists():
            QMessageBox.warning(self, 'Missing log', 'Please choose an existing .log file.'); return
        self.current_log = p
        # reset
        self.scf_history.clear(); self.step_rms.clear(); self.step_times.clear()
        self.link_history.clear(); self._last_link_seen = None
        self.history_table.setRowCount(0); self.energy_table.setRowCount(0)
        self._last_detected_max = None
        self.geom_disp.clear(); self.rmsd_series.clear(); self.prev_rmsd_series.clear()
        self._last_geom_for_prev = None; self._geom_count_seen = 0
        self.disp_line.setText('—'); self.disp_line.setStyleSheet("")
        self.rmsd_line.setText('—'); self._first_geom = None
        self.alarm_active = False; self.alarm_banner.setVisible(False); self.alarm_timer.stop()
        self._last_scf_count_seen = 0
        self.start_monitors()

# ------------------------ Utility: read route/charge from GJF ------------------------
def read_route_charge_from_gjf(path: Path) -> Tuple[str, Optional[int], Optional[int]]:
    route_lines: List[str] = []; charge = mult = None
    try: lines = path.read_text(errors="ignore").splitlines()
    except Exception: return "", None, None
    in_route = False; state = "start"
    for ln in lines:
        s = ln.strip()
        if s.startswith('#'):
            in_route = True; route_lines.append(s); state = "route"; continue
        if in_route:
            if s == "" or s.startswith("--"):
                in_route = False; state = "after_route_blank"; continue
            route_lines.append(s)
        if state == "after_route_blank" and s != "": state = "title"; continue
        if state == "title" and s == "": state = "charge_mult"; continue
        if state == "charge_mult" and s:
            parts = s.split()
            if len(parts) >= 2 and all(p.lstrip('+-').isdigit() for p in parts[:2]):
                try: charge, mult = int(parts[0]), int(parts[1])
                except Exception: pass
            break
    return " ".join(route_lines), charge, mult

# ------------------------ Main window ------------------------
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Gaussian Monitor — Robust Stop, Geometry Refresh, Displacement thresholds')
        self.resize(1280, 940)
        tabs = QTabWidget()
        tabs.addTab(RunGaussianTab(self), 'Run Gaussian')
        tabs.addTab(ParseOnlyTab(self), 'Parse Log Only')
        self.setCentralWidget(tabs)

def main():
    app = QApplication(sys.argv)
    w = MainWindow(); w.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()

