"""Conservative syntax/structure checks: warning != Gaussian execution guarantee."""
from dataclasses import dataclass
import re

@dataclass
class Issue:
    severity:str
    line:int
    column:int
    message:str
    fix:str

CHARGE=re.compile(r'^\s*([+-]?\d+)\s+(\d+)\s*$')
ROUTE=re.compile(r'^\s*#')
LINK1=re.compile(r'^\s*--Link1--\s*$',re.I)


def check_input(text):
    lines=text.splitlines()
    issues=[]
    # Iterate sections explicitly; blank line separators are significant.
    sections=[];beg=0
    for i,line in enumerate(lines):
        if LINK1.match(line):sections.append((beg,i,lines[beg:i]));beg=i+1
    sections.append((beg,len(lines),lines[beg:]))
    for start,end,segment in sections:
        if not any(line.strip() for line in segment):continue
        routes=[i for i,line in enumerate(segment) if ROUTE.match(line)]
        if not routes:
            issues.append(Issue('error',start+1,1,'Missing Gaussian route section beginning with #','Add a route such as #p hf/3-21g sp.'))
            continue
        route=routes[0]
        if route>0 and not all(x.lstrip().startswith('%') or not x.strip() for x in segment[:route]):
            issues.append(Issue('warning',start+1,1,'Unexpected content before the route section','Check Link0 keywords and the route placement.'))
        k=route+1
        while k<len(segment) and segment[k].strip() and not LINK1.match(segment[k]):k+=1
        if k>=len(segment):
            issues.append(Issue('error',start+len(segment),1,'Missing blank line after route section','Insert an empty line after the route card.'))
            continue
        k+=1
        title_start=k
        while k<len(segment) and segment[k].strip():k+=1
        if k==title_start:
            issues.append(Issue('warning',start+title_start+1,1,'Empty Gaussian title card','Add a descriptive title.'))
        if k>=len(segment):
            issues.append(Issue('error',start+len(segment),1,'Missing blank line after title','Add an empty line between title and charge/multiplicity.'))
            continue
        k+=1
        while k<len(segment) and not segment[k].strip():k+=1
        if k>=len(segment):
            issues.append(Issue('error',start+end-start,1,'Missing charge and multiplicity line','Add e.g. 0 1 and geometry or appropriate geom=check options.'))
            continue
        match=CHARGE.match(segment[k])
        if not match:
            issues.append(Issue('error',start+k+1,1,'Expected charge and multiplicity (e.g. 0 1)','Insert a valid charge and a positive integer multiplicity.'))
        elif int(match[2])<1:
            issues.append(Issue('error',start+k+1,1,'Spin multiplicity must be ≥1','Use a physically justified positive multiplicity.'))
        route_text=' '.join(x.strip().lower() for x in segment[route:route+5] if x.strip() and not x.startswith('%'))
        if 'geom=connectivity' in route_text and not any(re.match(r'^\s*\d+\s+\d+\s+[\d.]+',x) for x in segment[k+1:]):
            issues.append(Issue('warning',start+route+1,1,'geom=connectivity set, but no connectivity bonds were recognized','Inspect the geometry and connectivity block; some valid formats may not match this heuristic.'))
        if end==len(lines) and text and not text.endswith('\n\n'):
            issues.append(Issue('warning',len(lines),max(1,len(lines[-1]) if lines else 1),'No final blank line after input','Ensure required Gaussian blocks are separated and terminated by blank lines.'))
    for i,line in enumerate(lines):
        if '\t' in line:
            issues.append(Issue('info',i+1,line.index('\t')+1,'Tab character detected','Use visible whitespace mode to inspect formatting; tabs are not automatically invalid.'))
        if line.rstrip(' ')!=line:
            issues.append(Issue('info',i+1,len(line.rstrip(' '))+1,'Trailing spaces','Consider removing trailing spaces where unintended.'))
    return issues
