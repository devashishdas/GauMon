"""Gaussian log event and timing analytics. No Gaussian executable required.

Important: link completion durations are Gaussian's *reported* elapsed times;
unfinished links and gaps are not silently assigned a duration.
"""
from __future__ import annotations
from collections import defaultdict, deque
from itertools import pairwise
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
import os
import re
import statistics
import time

LINK_DESC={1:'Setup / job instructions',101:'Molecule and title',103:'Geometry optimizer',202:'Molecule orientation',301:'Basis setup',302:'One-electron integrals',401:'Initial orbitals',502:'Standard SCF',508:'Quadratic SCF',601:'Electron density and properties',701:'Derivative integrals',702:'Integral derivatives',703:'Two-electron derivatives',716:'Forces and frequencies',9999:'Job finishing'}
ENTER=re.compile(r'\(Enter\s+.*?[\\/]l(\d+)\.exe\)', re.I)
LEAVE=re.compile(r'Leave Link\s+(\d+)\s+at\s+(.+?),\s+MaxMem=.*?cpu:\s*([\d.]+)\s+elap:\s*([\d.]+)',re.I)
STEP=re.compile(r'\bStep number\s+(\d+)\s+out of a maximum of\s+(\d+)',re.I)
SCF=re.compile(r'SCF Done:\s*E\([^)]*\)\s*=\s*(-?[\d.]+(?:[DEde][+-]?\d+)?)\s+A\.U\.\s+after\s+(\d+)\s+cycles',re.I)
TERM_OK=re.compile(r'Normal termination of Gaussian',re.I)
TERM_BAD=re.compile(r'Error termination(?:\s+via\s+Lnk1e)?',re.I)
FREQ=re.compile(r'Frequencies\s*--',re.I)
RUN_BANNER=re.compile(r'Entering Gaussian System',re.I)
NUM_D=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DEde][+-]?\d+)?'
ERROR_RULES=(
 ('error','Internal-coordinate problem (FormBX)',r'FormBX had a problem|Bend failed for angle|Tors failed for dihedral','Geometry optimizer could not construct/use an angle or torsion. Inspect reported atoms, geometry and constraints; consider Cartesian optimization for an unconstrained job.'),
 ('error','Electronic SCF did not converge',r'Convergence failure|SCF has not converged|SCF failed to converge','Electronic iterations did not converge. Inspect spin, method, initial guess and DIIS/XQC behavior.'),
 ('error','Gaussian terminated with an error',r'Error termination(?:\s+via\s+Lnk1e)?','Look above this line for the first specific Gaussian error. This line alone does not explain the cause.'),
 ('warning','Optimization step limit',r'Number of steps exceeded|Optimization stopped','Optimization did not meet its convergence criteria within the allowed steps.'),
 ('warning','Internal angle near linear',r'Linear angle in Bend|Linear angle in Tors','Check whether nearly collinear atoms make the internal-coordinate description unstable.'),
)

def adjacent_step_durations(pairs):
    """Accept deques or any iterable; do not slice a deque in the GUI refresh.

    Only strictly consecutive, positive-duration observations count. The
    observations must be live timestamps, not historical backfilled steps.
    """
    result=[]
    for (a,t),(b,u) in pairwise(pairs):
        if b==a+1 and u>t:
            result.append((b,u-t))
    return result


@dataclass
class Problem:
    level:str
    title:str
    line:int
    evidence:str
    explanation:str

@dataclass
class LinkEvent:
    link:int
    elapsed:float
    cpu:float
    line:int
    date:str
    description:str

@dataclass
class Snapshot:
    path:str=''
    line_count:int=0
    events:deque[LinkEvent]=field(default_factory=lambda:deque(maxlen=2000))
    link_aggregate:dict[int,list]=field(default_factory=dict)
    event_count:int=0
    steps:deque[tuple[int,int,int]]=field(default_factory=lambda:deque(maxlen=6000)) # step/max/line
    problems:list[Problem]=field(default_factory=list)
    scf:deque[tuple[float,int,int]]=field(default_factory=lambda:deque(maxlen=15000))
    frequencies:int=0
    termination:str=''
    current_link:int|None=None
    last_growth:float|None=None
    observed_step_times:deque[tuple[int,float]]=field(default_factory=lambda:deque(maxlen=500))
    # Clock times printed alongside Leave Link are stage completion stamps,
    # not exact times of optimization step boundaries.

    def totals(self):
        return {k:list(v) for k,v in self.link_aggregate.items()}

    def eta(self, horizons=(10,50,100,150)):
        """Measured active-tail interval estimate; never treats historical backfill as timing."""
        gaps=[duration for step,duration in adjacent_step_durations(self.observed_step_times)]
        if self.termination:
            return {'available':False,'reason':'Job has reached a terminal marker; no future ETA.','horizons':[]}
        if len(gaps)<2:
            return {'available':False,'reason':'Waiting for ≥2 complete, newly observed optimization-step intervals. Historical steps are not timed by attachment.','horizons':[]}
        recent=gaps[-7:]
        med=statistics.median(recent)
        # Spread here is descriptive, not a calibrated prediction interval.
        lo=min(recent);hi=max(recent)
        return {'available':True,'pace_s':med,'sample_count':len(gaps),'min_observed_s':lo,'max_observed_s':hi,
                'uncertainty':'Uncalibrated observed pace range, NOT a probabilistic confidence interval.',
                'horizons':[{'steps':h,'seconds':med*h,'earliest_rough_s':lo*h,'latest_rough_s':hi*h} for h in horizons]}

class IncrementalLog:
    """Read append-only log chunks without re-reading previously parsed lines."""
    def __init__(self,path:Path):
        self.path=Path(path).expanduser().resolve()
        self.snapshot=Snapshot(path=str(self.path))
        self.offset=0
        self.partial=b''
        self.identity=None
        self.last_step=None
        self.backfill=True
        self.last_size=0
        self.last_growth_mono=None
        self.last_poll_error=''

    def poll(self,max_bytes=512*1024):
        if not self.path.is_file():
            self.last_poll_error='Log not found or not readable.'
            return self.snapshot
        try:
            st=self.path.stat()
            identity=(st.st_dev,st.st_ino)
            if self.identity is not None and (self.identity!=identity or st.st_size<self.offset):
                self.snapshot=Snapshot(path=str(self.path));self.offset=0;self.partial=b'';self.backfill=True;self.last_step=None
            self.identity=identity
            initial=self.backfill
            with self.path.open('rb') as fh:
                fh.seek(self.offset)
                data=fh.read(max_bytes)
            self.offset+=len(data)
            if data:
                if not initial:self.last_growth_mono=time.monotonic()
                buf=self.partial+data
                *lines,self.partial=buf.split(b'\n')
                for line in lines:
                    self._feed(line.decode('utf-8','replace').rstrip('\r'),historical=initial)
            if self.offset>=st.st_size:
                self.backfill=False
            self.last_poll_error=''
        except OSError as exc:
            self.last_poll_error=str(exc)
        return self.snapshot

    def _feed(self,line,historical):
        s=self.snapshot;s.line_count+=1;n=s.line_count
        m=ENTER.search(line)
        if m:s.current_link=int(m[1])
        m=LEAVE.search(line)
        if m:
            lid=int(m[1]);elapsed=float(m[4]);cpu=float(m[3]);date=m[2].strip()
            s.events.append(LinkEvent(lid,elapsed,cpu,n,date,LINK_DESC.get(lid,'Other Gaussian module')))
            s.event_count+=1
            counts=s.link_aggregate.setdefault(lid,[0.0,0.0,0])
            counts[0]+=elapsed;counts[1]+=cpu;counts[2]+=1
        m=STEP.search(line)
        if m:
            step,maximum=int(m[1]),int(m[2]);s.steps.append((step,maximum,n))
            if not historical and self.last_step!=step:
                s.observed_step_times.append((step,time.monotonic()))
            self.last_step=step
        m=SCF.search(line)
        if m:
            s.scf.append((float(m[1].replace('D','E').replace('d','E')),int(m[2]),n))
        if FREQ.search(line):s.frequencies+=1
        if TERM_OK.search(line):s.termination='COMPLETED'
        if TERM_BAD.search(line):s.termination='FAILED'
        for level,title,pattern,description in ERROR_RULES:
            if re.search(pattern,line,re.I):
                # Bounded findings, but retain all distinct matched error evidence
                if len(s.problems)<200:
                    s.problems.append(Problem(level,title,n,line.strip()[:500],description))
                break


def inspect_text(text):
    """For deterministic unit tests and headless analysis."""
    obj=IncrementalLog(Path('not-used.log'))
    for line in text.splitlines():obj._feed(line,historical=True)
    return obj.snapshot
