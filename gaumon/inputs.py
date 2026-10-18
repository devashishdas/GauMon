"""Gaussian input validation and small, reproducible test templates.

Input generation uses public Gaussian input syntax and contains no Gaussian
proprietary binaries or code. Test outcomes are not guaranteed by templates.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import re

@dataclass(frozen=True)
class Case:
    key: str
    label: str
    formula: str
    description: str
    expected: str
    filename: str

CASES = (
    Case('sp', 'Single-point energy', 'H₂', 'Minimal H₂ electronic energy', 'SCF energy; normal termination', '01_energy_H2.gjf'),
    Case('optfreq', 'Optimization + frequency', 'H₂O', 'Optimize water, then calculate frequencies', 'Optimization complete; frequency lines; no imaginary modes expected', '02_opt_freq_H2O.gjf'),
    Case('freq', 'Frequency only', 'H₂O', 'Harmonic frequencies at a supplied near-equilibrium geometry', 'Frequency lines; normal termination', '03_freq_H2O.gjf'),
    Case('ts', 'Transition-state search', 'NH₃', 'Planar ammonia umbrella inversion saddle point', 'TS search, then one imaginary mode if convergence succeeds', '04_ts_NH3_inversion.gjf'),
)
ADVANCED_CASES = (
    Case('trimer', 'Advanced: water trimer Opt+Freq', '(H₂O)₃',
         'Hydrogen-bonded 9-atom optimization and harmonic frequency',
         'Several optimization cycles; possible structural rearrangement',
         '05_advanced_water_trimer_optfreq.gjf'),
    Case('stretched', 'Advanced: stretched H₂ UHF', 'H₂',
         'Long-bond unrestricted SCF with mixed initial guess and XQC fallback',
         'Check iterative SCF stability, spin contamination, and convergence',
         '06_advanced_stretched_H2_UHF.gjf'),
    Case('fe2s2', 'Advanced: Fe₂S₂(SH)₄ toy SP', 'Fe₂S₆H₄',
         'Illustrative open-shell metal cluster (not a validated broken-symmetry model)',
         'SCF stability and convergence; investigate spin-state sensitivity',
         '07_advanced_Fe2S2_toy_SP.gjf'),
)
ALL_CASES = CASES + ADVANCED_CASES
BY_KEY = {case.key: case for case in ALL_CASES}
ROUTE = re.compile(r'^\s*#', re.I)
CM = re.compile(r'^\s*(-?\d+)\s+(\d+)\s*$')
XYZ = re.compile(r'^\s*([A-Z][a-z]?|\d+)\s+([-+\d.EeDd]+)\s+([-+\d.EeDd]+)\s+([-+\d.EeDd]+)(?:\s|$)')


def validate_input(text: str) -> list[str]:
    """Validate the common Gaussian route/title/charge/Cartesian layout.

    A conservative check for the 4 shipped templates; it does not replace
    Gaussian's parser. Existing advanced inputs may use ModRedundant/ChkAll etc.
    """
    lines = text.replace('\r', '').split('\n')
    errors = []
    route_start = next((i for i, s in enumerate(lines) if ROUTE.match(s)), None)
    if route_start is None:
        return ['Missing Gaussian route section (line starting with #).']
    i = route_start + 1
    while i < len(lines) and lines[i].strip():
        i += 1
    if i == len(lines):
        return ['Route section not followed by an empty separator line.']
    i += 1
    title = []
    while i < len(lines) and lines[i].strip():
        title.append(lines[i].strip()); i += 1
    if not title:
        errors.append('Title Card Required: missing non-empty title between route and charge/multiplicity.')
    if i == len(lines):
        errors.append('Missing blank separator after title.'); return errors
    i += 1
    if i >= len(lines) or not CM.match(lines[i]):
        errors.append('Missing charge/multiplicity line such as 0 1.'); return errors
    i += 1
    if not any(XYZ.match(s) for s in lines[i:]):
        errors.append('No recognizable Cartesian atomic coordinates after charge/multiplicity.')
    if not text.replace('\r\n','\n').endswith('\n\n'):
        errors.append('End of file in ZSymb risk: add an empty line after the final input section (including coordinates).')
    route = ' '.join(lines[route_start:i]).lower()
    if 'geom=connectivity' in route and not any(re.match(r'^\s*\d+\s+\d+(?:\s+\d+(?:\.\d+)?)?', line)
                                                 for line in lines[i:] if line.strip()):
        # This is only a warning: connectivity could be provided indirectly.
        errors.append('geom=connectivity requested but no connectivity block was recognized; remove the keyword or supply bonds.')
    return errors


def load_case(case_key: str, root: Path | None = None) -> tuple[Case, str]:
    case = BY_KEY[case_key]
    root = root or Path(__file__).resolve().parents[1] / 'examples' / 'inputs'
    text = (root / case.filename).read_text(encoding='utf-8')
    return case, text


def write_case(case_key: str, output_dir: Path, *, overwrite: bool = False) -> Path:
    case, content = load_case(case_key)
    issues = validate_input(content)
    if issues:
        raise ValueError(f'Template invalid: {"; ".join(issues)}')
    dest = Path(output_dir).expanduser().resolve()
    dest.mkdir(parents=True, exist_ok=True)
    outfile = dest / case.filename
    if outfile.exists() and not overwrite:
        raise FileExistsError(f'Input already exists: {outfile}')
    outfile.write_text(content, encoding='utf-8')
    return outfile
