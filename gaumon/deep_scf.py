"""Incremental deep Gaussian SCF-cycle diagnostics.

One block represents one SCF procedure at one geometry; it is NOT an
optimization cycle. Records are kept in bounded, opt-in-viewable history.
Original Gaussian logs are never stored in the distributed package.
"""
from __future__ import annotations
from collections import deque
from copy import deepcopy
import re
from typing import Any

# Supports Gaussian's Fortran D exponents, including adjacent fixed-width values.
_NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DEde][+-]?\d+)?'
NUMBER_RE = re.compile(_NUMBER)
CYCLE_RE = re.compile(r'^\s*Cycle\s+(\d+)\s+Pass\s+(\d+)\s+IDiag\s+(\d+)\s*:', re.I)
SCF_DONE_RE = re.compile(r'\bSCF Done\s*:', re.I)
STEP_RE = re.compile(r'\bStep number\s+(\d+)\s+out of a maximum of\s+(\d+)', re.I)
FIELD_RE = re.compile(r'(?<![A-Za-z])(?P<key>RMSU|RMSDP|MaxDP|OVMax|Delta-E|DE|EnMin|ErrMin|ErrMax|BMatC|BMatP|WtCom|WtEn|GapD|Gap|RMSU|NSaved|IDIUse|EMaxC|Rises|Damp|E)\s*=\s*(?P<value>'+_NUMBER+r'|[TF])\b', re.I)
DIIS_RE = re.compile(r'\bDIIS:\s*error\s*=\s*('+_NUMBER+r')', re.I)
COEFF_RE = re.compile(r'^\s*(Coeff-Com|Coeff-En|Coeff|CP):\s*(.*)', re.I)
CP_INLINE_RE = re.compile(r'\bCP:\s*(.*)', re.I)
RMSU_RE = re.compile(r'\bRMSU\s*=\s*('+_NUMBER+r')', re.I)
MAX_BLOCKS = 120
MAX_ROWS_PER_BLOCK = 400


def real_number(value: str):
    """Fortran-style double precision converted without changing units."""
    try:
        return float(value.replace('d', 'E').replace('D', 'E'))
    except (ValueError, TypeError):
        return None


def _fresh_block(index: int):
    return {'id': index, 'label': f'SCF #{index}', 'optimization_step': None,
            'status': 'running', 'rows': [], 'final_energy': None,
            'reported_cycles': None, 'truncated': False}


class DeepSCFParser:
    """Receives complete lines from the parent incremental parser.

    Starts a new block on each new Cycle 1 (or a cycle index reset). On SCF
    completion the last block freezes, and subsequent Cycle 1 begins a new
    blank live block. A later Berny Step number labels the just-finished block.
    """
    def __init__(self):
        self.reset()

    def reset(self):
        self.blocks = deque(maxlen=MAX_BLOCKS)
        self.current = None
        self.last_completed = None
        self.next_id = 1
        self.revision = 0
        self.row = None
        self.coeff_kind = None
        self.last_step = None

    def _start(self, cycle: int, pass_n: int, idiag: int):
        block = _fresh_block(self.next_id)
        self.next_id += 1
        self.blocks.append(block)
        self.current = block
        self.row = None
        self.coeff_kind = None
        self._new_row(cycle, pass_n, idiag)

    def _new_row(self, cycle: int, pass_n: int, idiag: int):
        row = {'cycle': cycle, 'pass': pass_n, 'idiag': idiag,
               'energy': None, 'delta_e': None, 'diis_error': None,
               'rmsu': None, 'rmsdp': None, 'maxdp': None, 'ovmax': None,
               'gap': None, 'gaps': [], 'cp': [], 'coeff': [],
               'coeff_com': [], 'coeff_en': [], 'nsaved': None,
               'wtcom': None, 'wten': None, 'rises': None, 'damp': None}
        if len(self.current['rows']) >= MAX_ROWS_PER_BLOCK:
            self.current['rows'].pop(0)
            self.current['truncated'] = True
        self.current['rows'].append(row)
        self.row = row
        self.coeff_kind = None
        self.revision += 1

    def feed_line(self, line: str):
        m = CYCLE_RE.search(line)
        if m:
            cycle, pass_n, idiag = map(int, m.groups())
            if (self.current is None or self.current['status'] != 'running' or
                    (self.row is not None and cycle < self.row['cycle'])):
                self._start(cycle, pass_n, idiag)
            elif (self.row is not None and (cycle, pass_n) == (self.row['cycle'], self.row['pass'])):
                # Repeated header, no duplicate data point.
                pass
            else:
                self._new_row(cycle, pass_n, idiag)
            return
        if SCF_DONE_RE.search(line):
            if self.current is not None and self.current['status'] == 'running':
                self.current['status'] = 'complete'
                self.last_completed = self.current
                self.current['reported_cycles'] = (int(re.search(r'after\s+(\d+)\s+cycles?', line, re.I).group(1))
                    if re.search(r'after\s+(\d+)\s+cycles?', line, re.I) else None)
                energy_match = re.search(r'\bE\([^)]*\)\s*=\s*('+_NUMBER+r')', line, re.I)
                self.current['final_energy'] = real_number(energy_match.group(1)) if energy_match else None
                self.revision += 1
            return
        m = STEP_RE.search(line)
        if m:
            step = int(m[1])
            self.last_step = step
            # Gaussian prints optimizer step after the SCF Done at that geometry.
            if self.last_completed is not None and self.last_completed['optimization_step'] is None:
                self.last_completed['optimization_step'] = step
                self.last_completed['label'] = f'Optimization {step} · SCF #{self.last_completed["id"]}'
                self.revision += 1
            return
        if self.current is None or self.current['status'] != 'running' or self.row is None:
            return
        row = self.row
        # Gaussian uses continuation lines for coefficients and CP. A label
        # line changes the active continuation, any other line ends it.
        coeff_match = COEFF_RE.match(line)
        if coeff_match:
            kind, payload = coeff_match.groups()
            key = kind.casefold().replace('-', '_')
            target = {'coeff_com':'coeff_com', 'coeff_en':'coeff_en',
                      'coeff':'coeff', 'cp':'cp'}[key]
            row[target].extend(v for t in NUMBER_RE.findall(payload) if (v := real_number(t)) is not None)
            self.coeff_kind = target
            self.revision += 1
            return
        if self.coeff_kind is not None and not line.lstrip().startswith(('E=', 'DIIS:', 'NSaved=', 'ErrMax=', 'Gap=', 'GapD=', 'RMSDP=', 'RMSU=', 'IDIUse=')):
            if line.startswith(' ' * 15) and not re.search(r'[A-Za-z]', line):
                numbers = NUMBER_RE.findall(line)
                if numbers:
                    row[self.coeff_kind].extend(v for t in numbers if (v := real_number(t)) is not None)
                    self.revision += 1
                    return
        self.coeff_kind = None
        # CP often shares the RMSU line, e.g. RMSU= ... CP: ...
        cp_inline = CP_INLINE_RE.search(line)
        if cp_inline:
            row['cp'].extend(v for t in NUMBER_RE.findall(cp_inline.group(1))
                             if (v := real_number(t)) is not None)
            self.coeff_kind = 'cp'
            self.revision += 1
        values = list(FIELD_RE.finditer(line))
        for match in values:
            key, val = match.group('key').lower(), match.group('value')
            target = {'e':'energy', 'delta-e':'delta_e', 'de':'de',
                      'rmsu':'rmsu', 'rmsdp':'rmsdp', 'maxdp':'maxdp',
                      'ovmax':'ovmax', 'gap':'gap', 'rises':'rises',
                      'damp':'damp', 'nsaved':'nsaved', 'wtcom':'wtcom',
                      'wten':'wten', 'enmin':'enmin', 'errmin':'errmin',
                      'errmax':'errmax', 'bmatc':'bmatc', 'bmatp':'bmatp',
                      'idiuse':'idiuse', 'emaxc':'emaxc', 'gapd':'gapd'}.get(key)
            if target is None:
                continue
            num = val if target in ('rises','damp') else real_number(val)
            if target == 'gap':
                if num is not None: row['gaps'].append(num); row['gap'] = num
            else:
                row[target] = num
            self.revision += 1
        diis = DIIS_RE.search(line)
        if diis:
            row['diis_error'] = real_number(diis.group(1))
            self.revision += 1
        # Some Gaussian printouts omit Delta-E, but include DE on the summary.
        if row['delta_e'] is None and row.get('de') is not None:
            row['delta_e'] = row['de']

    def overview(self):
        blocks = list(self.blocks)
        current = self.current
        return {
            'revision': self.revision,
            'active_block_id': current['id'] if current else None,
            'active_block_status': current['status'] if current else None,
            'active_cycle': self.row['cycle'] if current and current['status'] == 'running' and self.row else None,
            'blocks': [{'id': b['id'], 'label': b['label'], 'status': b['status'],
                        'optimization_step': b['optimization_step'],
                        'cycles': len(b['rows']), 'reported_cycles': b['reported_cycles'],
                        'truncated': b['truncated']} for b in blocks],
        }

    def get_block(self, block_id=None):
        block = self.current if block_id is None else next((b for b in self.blocks if b['id'] == block_id), None)
        return deepcopy(block) if block is not None else None
