"""Incremental Gaussian quadratic-convergence SCF / LinEq1 diagnostics.

LinEq1 belongs to the coupled/quadratic SCF linear solver (typically L508),
NOT vibrational frequencies. All history is bounded; no external files copied.
"""
from collections import deque
from copy import deepcopy
import re

NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DdEe][+-]?\d+)?'
QITER = re.compile(r'^\s*Iteration\s+(\d+)\s+EE\s*=\s*('+NUMBER+r')', re.I)
LIN = re.compile(r'^\s*LinEq1:\s*Iter\s*=\s*(\d+)\s+NonCon\s*=\s*(\d+)\s+RMS\s*=\s*('+NUMBER+r')\s+Max\s*=\s*('+NUMBER+r')\s+NDo\s*=\s*(\d+)', re.I)
LIMITS = re.compile(r'QCNR:\s*CnvC1\s*=\s*('+NUMBER+r')\s+CnvC2\s*=\s*('+NUMBER+r')',re.I)
CONVERGED = re.compile(r'Linear equations converged to\s*('+NUMBER+r')\s+('+NUMBER+r')\s+after\s+(\d+)\s+iterations',re.I)
FIELDS = re.compile(r'(?<![\w-])(?P<key>Delta-E|Grad|Max\s+rot)\s*=\s*(?P<value>'+NUMBER+r')',re.I)
LINK = re.compile(r'^\s*\(Enter\s+.*[/\\]l(\d+)\.exe\)',re.I)
FINISH = re.compile(r'\b(?:SCF Done|Error termination|Normal termination)\b',re.I)
MAX_SOLVES = 120
MAX_ITERS = 3000
MAX_OUTER = 400

def _number(value):
    return float(value.replace('D','E').replace('d','E'))

class DeepQCNRParser:
    """Consumes complete lines, retaining successive nested linear solves.

    New SCF/QCNR orbital iteration may create one LinEq1 solve. A new
    Iter=0 (or decreasing iteration index) starts a fresh inner solve.
    A solver is 'converged' only on NonCon=0 or explicit solver message.
    """
    def __init__(self):
        self.reset()

    def reset(self):
        self.solves = deque(maxlen=MAX_SOLVES)
        self.outer = deque(maxlen=MAX_OUTER)
        self.current = None
        self.outer_current = None
        self.current_link = None
        self.in_qcnr = False
        self.next_id = 1
        self.revision = 0
        self.epoch = 0

    def _close_active(self, status='interrupted'):
        if self.current and self.current['status'] == 'running':
            self.current['status'] = status
            self.revision += 1

    def _begin_solve(self, iteration):
        self._close_active()
        data = dict(id=self.next_id, label=f'QCNR {self.epoch} · outer {self.outer_current["iteration"] if self.outer_current else "—"} · solve {self.next_id}',
                    outer_iteration=self.outer_current['iteration'] if self.outer_current else None,
                    epoch=self.epoch, rows=[], status='running',
                    target_rms=self.outer_current.get('target_rms') if self.outer_current else None,
                    target_max=self.outer_current.get('target_max') if self.outer_current else None,
                    actual_convergence=None, truncated=False)
        self.next_id += 1
        self.current = data
        self.solves.append(data)
        self.revision += 1

    def feed_line(self, line):
        link = LINK.match(line)
        if link:
            new_link = f'L{link[1]}'
            if self.in_qcnr and new_link != self.current_link:
                self._close_active()
                self.in_qcnr=False
            self.current_link = new_link
        if 'Quadratic Convergence SCF Method' in line:
            self._close_active()
            self.epoch += 1
            self.in_qcnr = True
            self.outer_current = None
            self.revision += 1
            return
        if FINISH.search(line):
            self._close_active()
            self.in_qcnr=False
            return
        outer = QITER.match(line)
        if outer and (self.in_qcnr or self.current_link == 'L508'):
            self._close_active()
            r=dict(iteration=int(outer[1]),energy=_number(outer[2]),delta_e=None,
                   grad=None,max_rot=None,target_rms=None,target_max=None,epoch=self.epoch)
            for match in FIELDS.finditer(line):
                key=match['key'].lower().replace(' ','_').replace('-','_')
                r[key]=_number(match['value'])
            self.outer.append(r)
            self.outer_current = r
            self.revision += 1
            return
        limits = LIMITS.search(line)
        if limits and self.outer_current is not None:
            self.outer_current['target_rms']=_number(limits[1])
            self.outer_current['target_max']=_number(limits[2])
            self.revision += 1
            return
        lin = LIN.match(line)
        if lin:
            iteration,noncon,rms,mx,ndo=int(lin[1]),int(lin[2]),_number(lin[3]),_number(lin[4]),int(lin[5])
            if (self.current is None or self.current['status'] != 'running' or not self.current['rows'] or
                iteration <= self.current['rows'][-1]['iteration']):
                # The first line is a new solve; duplicate row retransmissions replace instead.
                if self.current and self.current['status']=='running' and self.current['rows'] and iteration==self.current['rows'][-1]['iteration']:
                    self.current['rows'][-1].update(noncon=noncon,rms=rms,max=mx,ndo=ndo)
                    self.revision += 1
                    return
                self._begin_solve(iteration)
            r=dict(iteration=iteration,noncon=noncon,rms=rms,max=mx,ndo=ndo)
            if len(self.current['rows']) >= MAX_ITERS:
                self.current['rows'].pop(0)
                self.current['truncated'] = True
            self.current['rows'].append(r)
            if noncon==0:
                self.current['status']='converged'
            self.revision += 1
            return
        done = CONVERGED.search(line)
        if done and self.current:
            self.current['status']='converged'
            self.current['actual_convergence'] = {'rms':_number(done[1]), 'max':_number(done[2]),'iterations_reported':int(done[3])}
            self.revision += 1

    def overview(self):
        return dict(revision=self.revision,
                    active_solve_id=self.current['id'] if self.current else None,
                    active_iteration=self.current['rows'][-1]['iteration'] if self.current and self.current['rows'] else None,
                    current_link=self.current_link,
                    in_qcnr=self.in_qcnr,
                    outer_iterations=len(self.outer),
                    solves=[{'id':b['id'],'label':b['label'],'status':b['status'],
                             'outer_iteration':b['outer_iteration'],'iterations':len(b['rows']),
                             'epoch':b['epoch'],'truncated':b['truncated']} for b in self.solves])

    def get_solve(self, solve_id=None):
        b=self.current if solve_id is None else next((s for s in self.solves if s['id']==solve_id),None)
        return deepcopy(b) if b is not None else None

    def outer_history(self):
        return deepcopy(list(self.outer))
