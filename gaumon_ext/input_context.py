"""Gaussian input discovery and route-keyword help, without executing Gaussian.

Auto-discovery is deliberately conservative: do not pretend a printed route
section is equivalent to the original .gjf file.
"""
from __future__ import annotations
from pathlib import Path
import re

# (meaning, plain-language caution). Text is advisory, never a validator.
KEYWORDS = {
    '#p': ('Verbose Gaussian output', 'Requests more detailed diagnostic printing, useful for GauMon and SCF diagnostics but increases log size.'),
    '#n': ('Normal Gaussian print level', 'Typical amount of printed output.'),
    '#t': ('Terse Gaussian output', 'Reduces log text and may omit details needed for advanced live diagnostics.'),
    'opt': ('Optimize geometry', 'Move atoms to find a stationary structure; the allowed maximum steps is not an ETA.'),
    'freq': ('Calculate vibrational frequencies', 'Analyze curvature/vibration at the present geometry; negative frequencies require interpretation.'),
    'sp': ('Single-point energy', 'Calculate electronic energy without moving atomic positions.'),
    'force': ('Compute atomic forces', 'Evaluate energy gradients; not a per-atom CPU-time measurement.'),
    'irc': ('Follow a reaction path', 'Trace an intrinsic reaction coordinate; may require many separate calculations.'),
    'scan': ('Scan a geometry coordinate', 'Compute a sequence of structures; total time depends on all scanned points.'),
    'td': ('Excited-state calculation', 'Time-dependent method; number of requested states strongly affects cost.'),
    'stable': ('Check wavefunction stability', 'A converged SCF solution need not be the lowest-energy stable solution.'),
    'pop': ('Population analysis', 'Mulliken/NBO-related options depend on which population method is requested.'),
    'density': ('Specify electronic density', 'Selects how a density matrix is used or obtained.'),
    'guess': ('Initial orbital guess', 'Controls the starting wavefunction; guess=mix helps break spin symmetry.'),
    'scf': ('Electronic convergence settings', 'scf=xqc can switch from conventional iterations to quadratic SCF; may take longer.'),
    'xqc': ('Extra quadratic SCF convergence', 'A fallback iterative method for challenging electronic structures.'),
    'qc': ('Quadratically convergent SCF', 'Alternative SCF solver; not a vibrational-frequency calculation.'),
    'diis': ('DIIS convergence acceleration', 'Mixes previous iteration information to help electronic convergence.'),
    'maxcycle': ('Maximum iteration count', 'A safety limit, NOT the number of iterations the job is expected to need.'),
    'nosymm': ('Disable use of molecular symmetry', 'Can simplify interpretation; may cost more computational time.'),
    'geom': ('Geometry input settings', 'geom=connectivity asks Gaussian to read the connectivity block.'),
    'connectivity': ('Bond connectivity data', 'Connectivity describes atom bonds; it does not assign computational time to bonds.'),
    'checkpoint': ('Read saved checkpoint information', 'Checkpoint-based options require compatible, existing checkpoint data.'),
    'scrf': ('Implicit solvent model', 'Continuum solvent treatment; may add iterations and cost.'),
    'pcm': ('Polarizable continuum model', 'An implicit solvation treatment.'),
    'solvent': ('Select the solvent', 'Name and availability depend on the chosen solvation model.'),
    'oniom': ('Multi-layer QM/MM or QM:QM model', 'Different molecular regions are treated at different levels of theory.'),
    'empiricaldispersion': ('Dispersion correction', 'Adds a dispersion treatment to the selected method.'),
    'int': ('Numerical integration settings', 'DFT integration grid choices affect precision and computational cost.'),
    'ultrafine': ('Fine numerical DFT grid', 'Usually more computationally demanding than coarser DFT grids.'),
    'grid': ('Integration-grid settings', 'Affects accuracy and CPU cost of numerical integration.'),
    'ub3lyp': ('Unrestricted B3LYP density functional', 'Allows different spin-up and spin-down orbitals, important for open-shell systems.'),
    'b3lyp': ('B3LYP density functional', 'A common hybrid density-functional approximation.'),
    'hf': ('Hartree–Fock electronic structure', 'An approximate electronic structure method.'),
    'uhf': ('Unrestricted Hartree–Fock', 'Allows separate alpha and beta molecular orbitals.'),
    'def2svp': ('def2-SVP basis set', 'A moderate-sized double-zeta basis; larger bases generally cost more.'),
    'def2tzvp': ('def2-TZVP basis set', 'A larger basis set with more functions and generally greater cost.'),
    'sdd': ('Stuttgart/Dresden basis and core potential', 'Often used for heavier elements; suitability depends on the atom and method.'),
    '6-31g': ('Pople split-valence basis family', 'The detailed polarization/diffuse-function specification matters.'),
    'nprocshared': ('Shared-memory processor count', 'Maximum cores requested, not a guarantee of proportional speedup.'),
    'mem': ('Memory limit', 'Requests a memory budget; insufficient memory can increase disk use.'),
    'chk': ('Checkpoint file', 'Stores state for later continuation or analysis; avoid overwriting unrelated jobs.'),
    'oldchk': ('Read a previous checkpoint', 'Requires an existing compatible checkpoint file.'),
    'rwf': ('Read/write scratch file', 'Scratch I/O may strongly affect performance.'),
}

LINK0_RE=re.compile(r'^\s*%\s*([A-Za-z][A-Za-z0-9]*)\s*=',re.I)
ROUTE_START=re.compile(r'^\s*#(?:[NpT])?(?:\s|$)',re.I)
TOKENS=re.compile(r'[A-Za-z][A-Za-z0-9]*(?:-[0-9A-Za-z]+)*|6-31g',re.I)


def keyword_info(token):
    special=token.strip().lower()
    if special in ('#p','#n','#t'):return KEYWORDS[special]
    key=token.lower().strip('%#')
    if key in KEYWORDS:return KEYWORDS[key]
    if key.startswith('6-31g'):return KEYWORDS['6-31g']
    return None


def explain_input(text):
    """Return meaningful occurrences in Link0/route sections only.

    Results are (line_number, keyword, meaning, caution) in source order.
    Does not claim all valid Gaussian grammar is recognized.
    """
    results=[]
    route=False
    for lineno, line in enumerate(text.splitlines(),1):
        if line.strip().lower()=='--link1--':route=False;continue
        m=LINK0_RE.match(line)
        if m:
            info=keyword_info(m[1]);
            if info:results.append((lineno,'%'+m[1],*info))
            continue
        if line.lstrip().startswith('#'):
            route=True
            mprint=re.match(r'^\s*#([NnPpTt])(?=\s|$)',line)
            if mprint:
                token='#'+mprint[1].lower()
                results.append((lineno,token,*KEYWORDS[token]))
        elif route and not line.strip():
            route=False;continue
        if not route:continue
        # Remove #p / #n/#t print-level marker, which is not part of the method.
        content=re.sub(r'^\s*#(?:[NnPpTt])?\s*', '', line, count=1)
        for m in TOKENS.finditer(content):
            info=keyword_info(m[0])
            if info:
                results.append((lineno,m[0],*info))
    return results


def matching_input_for_log(log_path, cwd=None):
    """Return only an existing same-stem .gjf/.com adjacent to the log.

    No arbitrary recursive search and no use of basename-only matches elsewhere:
    those could silently associate the wrong scientific input with a job.
    """
    log=Path(log_path).expanduser().resolve()
    stem=log.stem
    if stem.lower().endswith('.log'): stem=Path(stem).stem
    matches=[log.with_name(stem+ext) for ext in ('.gjf','.com','.GJF','.COM')]
    for p in matches:
        if p.is_file():return p
    return None
