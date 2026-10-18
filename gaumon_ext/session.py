"""Explicitly enrolled, local-only, resumable job records. Never auto-run jobs."""
from __future__ import annotations
import json,os,uuid,time
from datetime import datetime,timezone
from pathlib import Path

class SessionStore:
    def __init__(self,path=None):
        self.path=Path(path).expanduser() if path else Path.home()/'.local/state/gaumon/session_jobs.json'
        self.items=[]
        self.load()
    def load(self):
        try:
            data=json.loads(self.path.read_text(encoding='utf-8'))
            if isinstance(data,dict) and data.get('schema')==1:
                self.items=[x for x in data.get('jobs',[]) if isinstance(x,dict) and x.get('path')]
        except (OSError,ValueError,TypeError):self.items=[]
    def add(self,path,source='opened',input_path=None):
        p=str(Path(path).expanduser().resolve())
        for x in self.items:
            if x['path']==p:
                if input_path and str(Path(input_path).expanduser().resolve())!=x.get('input_path'):
                    x['input_path']=str(Path(input_path).expanduser().resolve());self.save()
                return x
        entry=dict(id=uuid.uuid4().hex,path=p,source=source,added=datetime.now(timezone.utc).isoformat(),last_known='UNVERIFIED',pid=None,input_path=str(Path(input_path).expanduser().resolve()) if input_path else None)
        self.items.append(entry);self.save();return entry
    def remove(self,job_id):
        self.items=[x for x in self.items if x['id']!=job_id];self.save()
    def save(self):
        self.path.parent.mkdir(parents=True,exist_ok=True)
        target=self.path.with_name(self.path.name+'.tmp')
        data=json.dumps({'schema':1,'jobs':self.items},indent=2)
        fd=os.open(str(target),os.O_WRONLY|os.O_CREAT|os.O_TRUNC,0o600)
        try:
            with os.fdopen(fd,'w',encoding='utf-8') as out:out.write(data)
            os.replace(target,self.path)
        finally:
            if target.exists():target.unlink()


def describe_observed_job(engine, now=None):
    """Human-readable session-list status; never assumes a saved job is attached.

    Loading a SessionStore restores paths only. A missing engine means the log
    has NOT been opened in this GUI session. A missing parser snapshot means
    there is no safe termination state to report yet. Neither means the
    external Gaussian job has finished, failed, or is running.
    """
    if engine is None:
        return 'SAVED · RESUME TO ATTACH'
    if getattr(engine, 'last_poll_error', ''):
        return 'UNREADABLE'
    snap = getattr(engine, 'snapshot', None)
    if snap is None:
        return 'ATTACHING · WAITING FOR LOG'
    term = getattr(snap, 'termination', '')
    if term == 'FAILED':
        return 'FAILED'
    if term == 'COMPLETED':
        return 'COMPLETED'
    growth = getattr(engine, 'last_growth_mono', None)
    if growth is not None:
        clock = time.monotonic() if now is None else now
        if clock - growth < 30 and clock >= growth:
            return 'LOG GROWING · PROCESS UNVERIFIED'
    return 'LOG ATTACHED · PROCESS UNVERIFIED'
