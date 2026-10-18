"""Incremental parsing of Gaussian printed harmonic frequencies.

Frequency mode output is separate from SCF and LinEq1 linear equation solves.
No eigensolutions are inferred: we display only quantities in the log.
"""
from collections import deque
from copy import deepcopy
import re

NUMBER=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DdEe][+-]?\d+)?'
NUM=re.compile(NUMBER)
FREQ=re.compile(r'^\s*Frequencies\s*--\s*(.*)',re.I)
PROPERTY=re.compile(r'^\s*(Red\. masses|Frc consts|IR Inten|Raman Activ|Depolar \(P\)|Depolar \(U\))\s*--\s*(.*)',re.I)
LINK=re.compile(r'^\s*\(Enter\s+.*[/\\]l(\d+)\.exe\)',re.I)
MAX_GROUPS=40
MAX_MODES=2000
FIELDS={'red. masses':'red_mass', 'frc consts':'force_constant',
        'ir inten':'ir_intensity','raman activ':'raman_activity',
        'depolar (p)':'depolar_p','depolar (u)':'depolar_u'}

class DeepFreqParser:
    def __init__(self):self.reset()

    def reset(self):
        self.groups=deque(maxlen=MAX_GROUPS)
        self.current=None
        self.next_id=1
        self.revision=0
        self.last_batch=[]
        self.next_mode=1
        self.current_link=None
        self.phase=0

    def feed_line(self, line):
        link=LINK.match(line)
        if link:
            self.current_link=f'L{link[1]}'
        if 'Harmonic frequencies (cm**-1)' in line:
            self.phase+=1
            if self.current and self.current['status']=='running':self.current['status']='printed'
            self.current=None
            self.last_batch=[]
            self.next_mode=1
            self.revision+=1
            return
        if 'Normal termination of Gaussian' in line or 'Error termination' in line:
            if self.current and self.current['status']=='running':
                self.current['status']='printed' if 'Normal termination' in line else 'interrupted'
                self.revision+=1
            return
        m=FREQ.match(line)
        if m:
            if self.current is None or self.current['status'] != 'running':
                self.current={'id':self.next_id,'label':f'Frequency set {self.next_id}',
                    'modes':[],'status':'running','phase':self.phase,'truncated':False}
                self.next_id+=1
                self.groups.append(self.current)
                self.next_mode=1
            vals=[]
            for x in NUM.findall(m[1]):
                try:vals.append(float(x.replace('D','E').replace('d','E')))
                except ValueError:pass
            self.last_batch=[]
            for v in vals:
                row={'mode':self.next_mode,'frequency':v,'red_mass':None,
                     'force_constant':None,'ir_intensity':None,'raman_activity':None,
                     'depolar_p':None,'depolar_u':None}
                if len(self.current['modes'])>=MAX_MODES:
                    self.current['modes'].pop(0)
                    self.current['truncated']=True
                self.current['modes'].append(row)
                self.last_batch.append(row)
                self.next_mode+=1
            self.revision+=1
            return
        p=PROPERTY.match(line)
        if p and self.current is not None and self.last_batch:
            key=FIELDS.get(p[1].strip().lower())
            if not key:return
            vals=[]
            for x in NUM.findall(p[2]):
                try:vals.append(float(x.replace('D','E').replace('d','E')))
                except ValueError:pass
            for row,val in zip(self.last_batch,vals):row[key]=val
            self.revision+=1

    def overview(self):
        return {'revision':self.revision,
                'active_group_id':self.current['id'] if self.current else None,
                'groups':[{'id':g['id'],'label':g['label'],'status':g['status'],
                           'modes':len(g['modes']),'imaginary':sum(m['frequency']<0 for m in g['modes'])}
                          for g in self.groups]}

    def get_group(self,group_id=None):
        g=self.current if group_id is None else next((x for x in self.groups if x['id']==group_id),None)
        return deepcopy(g) if g else None
