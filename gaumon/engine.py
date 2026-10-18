"""CPU-only Gaussian log parsing and owned-process lifecycle."""
from __future__ import annotations
import math
import os
from pathlib import Path
import re
import json
from datetime import datetime
import signal
import subprocess
import threading
import time
import uuid
from collections import deque
import numpy as np
from .deep_scf import DeepSCFParser
from .deep_qcnr import DeepQCNRParser
from .deep_freq import DeepFreqParser

SCF = re.compile(r'SCF Done:\s+E\([^)]*\)\s*=\s*([-+\d.EeDd]+)\s+A\.U\.\s+after\s+(\d+)\s+cycles?', re.I)
STEP = re.compile(r'Step number\s+(\d+)\s+out of a maximum of\s+(\d+)', re.I)
FORCE = re.compile(r'RMS\s+Force\s*=\s*([-+\d.EeDd]+)', re.I)
CHARGE = re.compile(r'Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)')
LINK = re.compile(r'(?:Enter\s+[^\n]*?[/\\])l(\d{2,4})\.exe', re.I)
NORMAL = re.compile(r'Normal termination of Gaussian\s+\d+', re.I)
ERROR = re.compile(r'Error termination', re.I)
FREQ = re.compile(r'Frequencies\s+--\s*(.*)', re.I)
OPT_COMPLETE = re.compile(r'Optimization completed\.', re.I)
LEAVE_TIME = re.compile(r'Leave Link\s+\d+\s+at\s+([^,]+),', re.I)
ERROR_EVIDENCE = re.compile(r'FormBX had a problem|Bend failed for angle|Tors failed for dihedral|Error termination|End of file in ZSymb|Convergence failure|Linear angle in (?:Bend|Tors)', re.I)
HEADER = re.compile(r'^\s*Number\s+Number\s+Type\s+X\s+Y\s+Z', re.I)
MAX_SERIES = 12000

ROOT = Path(__file__).resolve().parents[1]
LINK_MAP_PATH = ROOT / 'assets' / 'gaussian_links.json'
try:
    LINK_MAP = json.loads(LINK_MAP_PATH.read_text(encoding='utf-8')) if LINK_MAP_PATH.is_file() else {}
except Exception:
    LINK_MAP = {}


def format_seconds(seconds: float | None) -> str:
    if seconds is None or not math.isfinite(seconds) or seconds < 0:
        return '—'
    seconds = int(round(seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f'{h:02d}:{m:02d}:{s:02d}'


def kabsch(reference, current):
    a, b = np.asarray(reference, float), np.asarray(current, float)
    if a.shape != b.shape or a.ndim != 2 or a.shape[1] != 3 or len(a) == 0:
        return None, None
    x, y = a-a.mean(0), b-b.mean(0)
    u, _, vt = np.linalg.svd(y.T @ x)
    d = np.linalg.det(u @ vt)
    diag = np.diag([1, 1, -1 if d < 0 else 1])
    aligned = y @ (u @ diag @ vt)
    distances = np.linalg.norm(aligned-x, axis=1)
    return float(np.sqrt(np.mean(distances**2))), float(distances.max())


class Parser:
    def __init__(self):
        self.reset()

    def reset(self):
        self.sizes = 0
        self.buffer = ''
        self.rows = deque(maxlen=MAX_SERIES)
        self.forces = deque(maxlen=MAX_SERIES)
        self.rmsds = deque(maxlen=MAX_SERIES)
        self.displacements = deque(maxlen=MAX_SERIES)
        self.prev_rmsds = deque(maxlen=MAX_SERIES)
        self.geometries = 0
        self.initial = None
        self.previous = None
        self.latest_atoms = []
        self.step = None
        self.step_max = None
        self.charge = None
        self.mult = None
        self.link = None
        self.freq_blocks = 0
        self.frequency_values = []
        self.opt_completed = False
        self.termination = None
        self.route = ''
        self._route_pending = False
        self._geom_stage = 0
        self._geom_rows = []
        self._last_force = None
        self._last_energy = None
        self._last_cycles = None
        self.tail = deque(maxlen=130)
        self.step_events = deque(maxlen=MAX_SERIES)
        self.link_events = deque(maxlen=MAX_SERIES)
        self.step_seq = 0
        self.link_seq = 0
        self.scf_seq = 0
        self.revision = 0
        self.last_log_timestamp = None
        self.diagnostic_evidence = deque(maxlen=70)
        self.deep = DeepSCFParser()
        self.qcnr = DeepQCNRParser()
        self.freq_deep = DeepFreqParser()

    def feed(self, data: str):
        self.buffer += data
        lines = self.buffer.split('\n')
        self.buffer = lines.pop()
        for line in lines:
            self._line(line.rstrip('\r'))
        if len(self.buffer) > 2_000_000:
            self.buffer = self.buffer[-4096:]

    def _geometry_line(self, line):
        if 'Standard orientation:' in line or 'Input orientation:' in line or 'Z-Matrix orientation:' in line:
            self._geom_stage = 1
            self._geom_rows = []
            return
        if self._geom_stage == 1:
            if HEADER.search(line):
                self._geom_stage = 2
            return
        if self._geom_stage == 2:
            if line.strip().startswith('---'):
                self._geom_stage = 3
            return
        if self._geom_stage == 3:
            bits = line.split()
            if len(bits) >= 6 and bits[0].isdigit() and bits[1].isdigit():
                try:
                    self._geom_rows.append((int(bits[1]), [float(v.replace('D', 'E')) for v in bits[-3:]]))
                except ValueError:
                    pass
                return
            if self._geom_rows:
                self._accept_geometry(self._geom_rows)
            self._geom_stage = 0
            self._geom_rows = []

    def _accept_geometry(self, atoms):
        numbers = [a[0] for a in atoms]
        coords = [a[1] for a in atoms]
        if self.initial is None:
            self.initial = (numbers, coords)
        if self.initial[0] != numbers:
            return  # atom identities/order must match
        rms, disp = kabsch(self.initial[1], coords)
        prev, _ = kabsch(self.previous[1], coords) if self.previous and self.previous[0] == numbers else (0.0, 0.0)
        if rms is not None:
            self.geometries += 1
            x = self.geometries
            self.rmsds.append([x, rms]); self.displacements.append([x, disp])
            self.prev_rmsds.append([x, prev])
            self.latest_atoms = [{'element': n, 'x': round(c[0], 5), 'y': round(c[1], 5), 'z': round(c[2], 5)} for n,c in atoms]
            self.previous = (numbers, coords)

    def _line(self, line):
        self.tail.append(line)
        self.revision += 1
        self.deep.feed_line(line)
        self.qcnr.feed_line(line)
        self.freq_deep.feed_line(line)
        if ERROR_EVIDENCE.search(line):
            self.diagnostic_evidence.append(line)
        stamp = LEAVE_TIME.search(line)
        if stamp:
            try:
                parsed = datetime.strptime(' '.join(stamp.group(1).split()), '%a %b %d %H:%M:%S %Y')
                self.last_log_timestamp = parsed.timestamp()
            except (ValueError, OverflowError, OSError):
                pass
        self._geometry_line(line)
        if line.lstrip().startswith('#'):
            self.route = (self.route + ' ' + line.strip()).strip()[:600]
        m = SCF.search(line)
        if m:
            self._last_energy = float(m[1].replace('D', 'E'))
            self._last_cycles = int(m[2])
            self.scf_seq += 1
            self.rows.append([self.scf_seq, self._last_energy, self._last_cycles])
        m = STEP.search(line)
        if m:
            new_step, new_max = int(m[1]), int(m[2])
            # Record each distinct observed optimizer step, including when a complete
            # existing log is parsed in one read. Do not invent cycle timestamps.
            if new_step != self.step or not self.step_events:
                self.step_seq += 1
                self.step_events.append({'seq': self.step_seq, 'step': new_step,
                                         'step_max': new_max, 'log_timestamp': self.last_log_timestamp,
                                         'scf_cycles': self._last_cycles, 'energy': self._last_energy,
                                         'force': None})
            self.step, self.step_max = new_step, new_max
        m = FORCE.search(line)
        if m:
            self._last_force = float(m[1].replace('D', 'E'))
            self.forces.append([self.step or len(self.forces)+1, self._last_force])
            if self.step_events and self.step_events[-1]['step'] == self.step:
                self.step_events[-1]['force'] = self._last_force
        m = CHARGE.search(line)
        if m: self.charge, self.mult = int(m[1]), int(m[2])
        m = LINK.search(line)
        if m:
            new_link = 'L' + m[1]
            self.link = new_link
            self.link_seq += 1
            self.link_events.append({'seq': self.link_seq, 'link': new_link,
                                     'log_timestamp': self.last_log_timestamp})
        match_freq = FREQ.search(line)
        if match_freq:
            self.freq_blocks += 1
            for tok in match_freq.group(1).split():
                try: self.frequency_values.append(float(tok.replace('D', 'E')))
                except ValueError: pass
        if OPT_COMPLETE.search(line): self.opt_completed = True
        if NORMAL.search(line): self.termination = 'normal'
        if ERROR.search(line): self.termination = 'error'
        if 'Entering Gaussian System' in line and self.termination:
            self.termination = None  # new Link1 job

    def snapshot(self):
        def last(series): return series[-1][1] if series else None
        return dict(energy=self._last_energy, cycles=self._last_cycles,
                    step=self.step, step_max=self.step_max, force=self._last_force,
                    rmsd=last(self.rmsds), displacement=last(self.displacements),
                    prev_rmsd=last(self.prev_rmsds), geometries=self.geometries,
                    atoms=self.latest_atoms[:400], charge=self.charge, mult=self.mult,
                    link=self.link, route=self.route, frequencies=self.freq_blocks,
                    frequency_values=self.frequency_values[-100:],
                    imaginary_modes=sum(f < -1e-7 for f in self.frequency_values),
                    opt_completed=self.opt_completed,
                    termination=self.termination, scf=list(self.rows), forces=list(self.forces),
                    rmsds=list(self.rmsds), prev_rmsds=list(self.prev_rmsds),
                    displacements=list(self.displacements), tail='\n'.join(self.tail),
                    diagnostic_text='\n'.join(self.diagnostic_evidence), revision=self.revision,
                    deep_scf=self.deep.overview(), deep_qcnr=self.qcnr.overview(), deep_freq=self.freq_deep.overview())


class Monitor:
    TERMINAL = {'COMPLETED', 'FAILED', 'STOPPED', 'INCOMPLETE'}
    def __init__(self):
        self.lock = threading.RLock()
        self.proc = None
        self.pgid = None
        self.log_path = None
        self.mode = None
        self.parser = Parser()
        self.state = 'READY'
        self.detail = 'Choose a Gaussian input or attach an existing log'
        self.start_time = None
        self.end_time = None
        self.last_growth = None
        self.offset = 0
        self.file_identity = None
        self.stop_requested = False
        self.external_pid = None
        self.external_birth = None
        self.session_id = None
        self._last_poll = 0.0
        self.wall_start = None
        self.step_history = deque(maxlen=400)
        self.cycle_durations = deque(maxlen=150)
        self.last_seen_step = None
        self.last_step_elapsed = None
        self.link_history = deque(maxlen=120)
        self.last_seen_link = None
        self._step_event_cursor = 0
        self._link_event_cursor = 0
        self._prev_cycle_event = None
        self._last_live_step_time = None
        self._initial_backfill = False
        self._timer_origin = 'none'
        self._has_post_attach_growth = False

    def _begin(self, path, mode):
        self.parser.reset()
        self.log_path = Path(path)
        self.mode = mode
        self.state = 'STARTING'
        self.detail = 'Waiting for output'
        self.start_time = time.monotonic()
        self.end_time = None
        self.last_growth = self.start_time
        self.offset = 0
        self.file_identity = None
        self.stop_requested = False
        self.session_id = uuid.uuid4().hex
        self.wall_start = time.time()
        self.step_history.clear(); self.cycle_durations.clear(); self.link_history.clear()
        self.last_seen_step = None; self.last_step_elapsed = None; self.last_seen_link = None
        self._step_event_cursor = self._link_event_cursor = 0
        self._prev_cycle_event = None
        self._last_live_step_time = None
        self._initial_backfill = mode == 'external'
        self._timer_origin = 'none'
        self._has_post_attach_growth = False
        self._last_poll = 0.0

    def launch(self, executable, input_path, workdir, log_name):
        with self.lock:
            if self.proc and self.proc.poll() is None:
                raise ValueError('A GauMon-managed calculation is already running')
            input_path = Path(input_path).expanduser().resolve()
            if not input_path.is_file() or input_path.suffix.lower() not in ('.gjf', '.com'):
                raise ValueError('Choose an existing .gjf or .com input file')
            work = Path(workdir).expanduser().resolve() if workdir else input_path.parent
            if not work.is_dir(): raise ValueError('Working directory does not exist')
            name = log_name.strip() or input_path.with_suffix('.log').name
            if Path(name).name != name or name in ('.', '..'):
                raise ValueError('Log filename must be a simple filename, not a path')
            path = work / name
            if path.resolve() == input_path or (path.exists() and path.stat().st_size > 0):
                raise ValueError('Output file exists and is not empty; choose a new filename')
            if not executable.strip(): raise ValueError('Gaussian executable is required')
            self._begin(path, 'managed')
            self.external_pid = self.external_birth = None
            try:
                with input_path.open('rb') as src, path.open('wb') as out:
                    self.proc = subprocess.Popen([executable], stdin=src, stdout=out,
                                                 stderr=subprocess.STDOUT, cwd=work,
                                                 start_new_session=(os.name == 'posix'))
                self.pgid = self.proc.pid if os.name == 'posix' else None
                self.detail = f'Process {self.proc.pid} started'
            except Exception:
                self.proc = None
                self.state = 'FAILED'; self.end_time = time.monotonic()
                raise
            return self.snapshot()

    def attach(self, path, pid=None):
        with self.lock:
            if self.proc and self.proc.poll() is None:
                raise ValueError('A managed calculation is still running; cannot replace it')
            path = Path(path).expanduser().resolve()
            if not path.is_file(): raise ValueError('Log file does not exist')
            self.external_pid = self.external_birth = None
            if pid is not None:
                import psutil
                proc = psutil.Process(int(pid))
                if proc.uids().real != os.getuid():
                    raise ValueError('Can only attach to your own process')
                self.external_pid, self.external_birth = proc.pid, proc.create_time()
            self.proc = None
            self._begin(path, 'external')
            self.detail = 'Attached to existing log (read-only)'
            return self.snapshot()

    def _external_alive(self):
        if self.external_pid is None: return None
        try:
            import psutil
            proc = psutil.Process(self.external_pid)
            return abs(proc.create_time()-self.external_birth) < 0.01 and proc.is_running() and proc.status() != psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess: return False
        except (psutil.AccessDenied, OSError): return None

    def stop(self, force=False):
        with self.lock:
            if self.mode != 'managed' or self.proc is None or self.proc.poll() is not None:
                raise ValueError('No GauMon-owned running process to stop')
            self.stop_requested = True
            if self.pgid is not None:
                os.killpg(self.pgid, signal.SIGKILL if force else signal.SIGTERM)
            elif force: self.proc.kill()
            else: self.proc.terminate()
            self.detail = 'Force stop requested' if force else 'Graceful stop requested'
            if not force:
                proc, pgid, session = self.proc, self.pgid, self.session_id
                def escalate():
                    time.sleep(4)
                    with self.lock:
                        if self.proc is proc and self.session_id == session and proc.poll() is None:
                            try:
                                if pgid is not None: os.killpg(pgid, signal.SIGKILL)
                                else: proc.kill()
                            except ProcessLookupError: pass
                threading.Thread(target=escalate, daemon=True).start()
            return self.snapshot()

    def _ingest_events(self, now: float, *, initial_attachment: bool = False):
        """Replay all parsed events, including old records, without faking timestamps.

        Cycle durations use timestamps printed by Gaussian when present and
        monotonic. Otherwise they are inferred *only* from newly observed live
        steps, one step per poll, after a baseline observation has been made.
        """
        steps = [e for e in self.parser.step_events if e['seq'] > self._step_event_cursor]
        links = [e for e in self.parser.link_events if e['seq'] > self._link_event_cursor]
        if steps:
            self._step_event_cursor = steps[-1]['seq']
        if links:
            self._link_event_cursor = links[-1]['seq']
        # An attached backlog is NOT N cycles completed in zero seconds.
        for idx, event in enumerate(steps):
            timestamp = event.get('log_timestamp')
            prev = self._prev_cycle_event
            cycle_time = None
            source = 'unavailable'
            if prev and event['step'] > prev['step']:
                prev_t = prev.get('log_timestamp')
                if prev_t is not None and timestamp is not None and timestamp > prev_t:
                    gap = event['step'] - prev['step']
                    # Gaps indicate incomplete logs, resets, or missing optimizer events:
                    # do not claim one cycle time from a multi-step jump.
                    if gap == 1:
                        cycle_time = timestamp - prev_t
                        source = 'Gaussian link-time proxy'
                elif (not initial_attachment and len(steps) == 1 and
                      self._last_live_step_time is not None and event['step'] == prev['step']+1):
                    cycle_time = now - self._last_live_step_time
                    source = 'monitor observation'
            if initial_attachment:
                # Only explicit log timestamps can validate historical cycles.
                self._last_live_step_time = None
            else:
                # First post-attach update establishes a baseline, not a full cycle.
                if len(steps) == 1:
                    self._last_live_step_time = now
                else:
                    self._last_live_step_time = None
            if cycle_time is not None and 0 < cycle_time < 14 * 86400:
                self.cycle_durations.append(cycle_time)
                self._timer_origin = source
            else:
                cycle_time = None
            self.step_history.append({'step': event['step'],
                                      'observed_s': None if initial_attachment else round(max(0, now-self.start_time), 1),
                                      'cycle_time_s': round(cycle_time, 2) if cycle_time is not None else None,
                                      'timing_source': source if cycle_time is not None else 'unavailable',
                                      'log_timestamp': timestamp,
                                      'scf_cycles': event['scf_cycles'],
                                      'energy': event['energy'], 'force': event['force']})
            self._prev_cycle_event = event
            self.last_seen_step = event['step']
        for event in links:
            self.link_history.append({'link': event['link'],
                                      'observed_s': None if initial_attachment else round(max(0, now-self.start_time),1),
                                      'log_timestamp': event.get('log_timestamp'),
                                      'description': LINK_MAP.get(event['link'], '')})
            self.last_seen_link = event['link']

    def poll(self):
        with self.lock:
            now = time.monotonic()
            if self.state in self.TERMINAL or self.state == 'READY': return self.snapshot()
            if now-self._last_poll < 0.35: return self.snapshot()
            self._last_poll = now
            try:
                stat = self.log_path.stat()
                identity = (stat.st_dev, stat.st_ino)
                if identity != self.file_identity or stat.st_size < self.offset:
                    if self.file_identity is not None:
                        self.step_history.clear(); self.cycle_durations.clear(); self.link_history.clear()
                        self._step_event_cursor = self._link_event_cursor = 0
                        self._prev_cycle_event = None; self._last_live_step_time = None
                    self.parser.reset(); self.offset = 0; self.file_identity = identity
                if stat.st_size > self.offset:
                    if self.mode == 'external' and self.offset > 0:
                        self._has_post_attach_growth = True
                    with self.log_path.open('rb') as f:
                        f.seek(self.offset)
                        while True:
                            data = f.read(1024*1024)
                            if not data: break
                            self.parser.feed(data.decode('utf-8', errors='replace'))
                        self.offset = f.tell()
                    self.last_growth = now
                    self._ingest_events(now, initial_attachment=self._initial_backfill)
                    self._initial_backfill = False
            except OSError:
                pass
            alive = self.proc.poll() is None if self.mode == 'managed' and self.proc else self._external_alive()
            code = self.proc.poll() if self.mode == 'managed' and self.proc else None
            marker = self.parser.termination
            if alive is False:
                self.state = ('STOPPED' if self.stop_requested else
                              'FAILED' if marker == 'error' or (code is not None and code != 0) else
                              'COMPLETED' if marker == 'normal' else 'INCOMPLETE')
            elif alive is None and marker in ('normal', 'error'):
                self.state = 'COMPLETED' if marker == 'normal' else 'FAILED'
            elif alive is True and now-self.last_growth > 180:
                self.state = 'STALLED'
            elif self.offset:
                self.state = ('ATTACHED' if self.mode == 'external' and alive is None
                              and not self._has_post_attach_growth else 'RUNNING')
            else:
                self.state = 'STARTING'
            if self.state in self.TERMINAL:
                self.end_time = now
                self.detail = {'COMPLETED':'Gaussian terminated normally', 'FAILED':'Gaussian failed or exited with an error',
                               'STOPPED':'Owned calculation stopped', 'INCOMPLETE':'Process exited without final Gaussian termination'}[self.state]
            elif self.state == 'STALLED': self.detail = 'No new output for 180 seconds; process may still be working'
            elif self.state == 'RUNNING': self.detail = ('Verified process running' if alive is True else 'Log changed after attachment; process not verified')
            elif self.state == 'ATTACHED': self.detail = 'Existing log loaded; process state unverified (provide PID)'
            return self.snapshot()

    def snapshot(self):
        t = self.end_time if self.end_time is not None else time.monotonic()
        seconds = max(0, t-self.start_time) if self.start_time is not None else 0
        data = self.parser.snapshot()
        cycle_avg = (sum(self.cycle_durations) / len(self.cycle_durations)) if self.cycle_durations else None
        recent_avg = (sum(list(self.cycle_durations)[-5:]) / min(5, len(self.cycle_durations))) if self.cycle_durations else None
        base_epoch = time.time()
        eta_pace = recent_avg if recent_avg is not None else cycle_avg
        if self.state in self.TERMINAL:
            eta_pace = None
        eta_records = []
        for extra in (10, 50, 100, 150):
            if eta_pace is None:
                eta_records.append({'extra_cycles': extra, 'seconds': None, 'duration': '—', 'datetime': '—'})
            else:
                future_seconds = extra * eta_pace
                eta_records.append({
                    'extra_cycles': extra,
                    'seconds': round(future_seconds, 1),
                    'duration': format_seconds(future_seconds),
                    'datetime': datetime.fromtimestamp(base_epoch + future_seconds).strftime('%Y-%m-%d %H:%M:%S'),
                })
        remaining = None
        eta_to_max = None
        if data.get('step') is not None and data.get('step_max') is not None:
            remaining = max(0, data['step_max'] - data['step'])
            if eta_pace is not None:
                eta_to_max = {
                    'remaining_cycles': remaining,
                    'seconds': round(remaining * eta_pace, 1),
                    'duration': format_seconds(remaining * eta_pace),
                    'datetime': datetime.fromtimestamp(base_epoch + remaining * eta_pace).strftime('%Y-%m-%d %H:%M:%S'),
                }
        data.update(state=self.state, detail=self.detail, elapsed=round(seconds, 1),
                    pid=self.proc.pid if self.proc else self.external_pid,
                    mode=self.mode, log=str(self.log_path) if self.log_path else None,
                    owned=bool(self.proc and self.proc.poll() is None),
                    verified=bool((self.proc and self.proc.poll() is None) or self.external_pid),
                    session=self.session_id,
                    step_history=list(self.step_history),
                    cycle_durations=list(self.cycle_durations),
                    cycle_average_s=round(cycle_avg, 2) if cycle_avg is not None else None,
                    cycle_timing_source=self._timer_origin,
                    cycle_sample_count=len(self.cycle_durations),
                    recent_cycle_average_s=round(recent_avg, 2) if recent_avg is not None else None,
                    link_history=list(self.link_history),
                    link_description=LINK_MAP.get(data.get('link') or '', ''),
                    eta_records=eta_records,
                    eta_to_max=eta_to_max,
                    observed_started_at=datetime.fromtimestamp(self.wall_start).strftime('%Y-%m-%d %H:%M:%S') if self.wall_start else None)
        return data

    def get_deep_block(self, block_id=None):
        with self.lock:
            return self.parser.deep.get_block(block_id)

    def get_qcnr_solve(self, solve_id=None):
        with self.lock:
            return self.parser.qcnr.get_solve(solve_id)

    def get_qcnr_outer_history(self):
        with self.lock:
            return self.parser.qcnr.outer_history()

    def get_freq_group(self, group_id=None):
        with self.lock:
            return self.parser.freq_deep.get_group(group_id)

    def jobs(self):
        try:
            import psutil
            out = []
            for proc in psutil.process_iter(['name', 'cmdline', 'pid']):
                try:
                    if proc.uids().real != os.getuid(): continue
                    name = proc.info['name'] or ''
                    if name in ('g16','g09','g03') or re.fullmatch(r'l\d{3,4}\.exe', name, re.I):
                        out.append({'pid':proc.pid, 'name':name, 'command':' '.join(proc.info['cmdline'] or [])[:200]})
                except (psutil.AccessDenied, psutil.NoSuchProcess): pass
            return out
        except ImportError: return []
