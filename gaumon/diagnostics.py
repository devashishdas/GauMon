"""Conservative and actionable diagnoses of Gaussian output.

References: Zhe Wang, Gaussian Common Errors and Solutions,
https://wongzit.github.io/gaussian-common-errors-and-solutions/
The recommendations are independently phrased; no source content is copied.
"""
from __future__ import annotations
import re
from dataclasses import dataclass

GUIDE_URL = 'https://wongzit.github.io/gaussian-common-errors-and-solutions/'
GUIDE_CREDIT = 'Zhe Wang — Gaussian Common Errors and Solutions'

@dataclass(frozen=True)
class Diagnosis:
    severity: str
    title: str
    description: str
    action: str
    source_name: str = ''
    source_url: str = ''


def diagnose_log(text: str) -> list[Diagnosis]:
    """Prioritize specific failures before generic Gaussian termination lines."""
    found: list[Diagnosis] = []
    angle_matches = re.findall(r'Bend failed for angle\s+(\d+)\s*-\s*(\d+)\s*-\s*(\d+)', text, re.I)
    torsion_matches = re.findall(r'Tors failed for dihedral\s+(\d+)\s*-\s*(\d+)\s*-\s*(\d+)\s*-\s*(\d+)', text, re.I)
    formbx = re.search(r'FormBX had a problem|Linear angle in (?:Bend|Tors)|Error in internal coordinate system', text, re.I)
    if formbx:
        angles = ['–'.join(t) for t in dict.fromkeys(angle_matches)]
        torsions = ['–'.join(t) for t in dict.fromkeys(torsion_matches)]
        evidence = []
        if angles:
            evidence.append('Failed angle(s): ' + ', '.join(angles[:12]))
        if torsions:
            evidence.append('Failed torsion(s): ' + ', '.join(torsions[:12]))
        evidence.append('Berny could not construct/use one or more redundant internal coordinates (typically l103.exe).'
                        ' Near-collinear atoms or unsuitable coordinate definitions are plausible; the exact geometry must be checked.')
        found.append(Diagnosis('ERROR', 'Internal-coordinate failure — FormBX (L103)',
            ' '.join(evidence),
            'Inspect the latest valid structure, especially the numbered atoms and any constraints/scan definitions. '
            'For an unconstrained optimization, restart from a sensible geometry using Opt=Cartesian; '
            'this may take more steps. Do not blindly combine Opt=Cartesian with Opt=ModRedundant. '
            'Alternatively, regenerate the internal coordinates from the last valid structure.',
            GUIDE_CREDIT, GUIDE_URL))
    if 'End of file in ZSymb' in text:
        found.append(Diagnosis('ERROR', 'Incomplete geometry / input (ZSymb, L101)',
            'Gaussian reached end of input while reading its molecular specification.',
            'Check the terminating blank line after coordinates; when geom=connectivity is used, '
            'provide connectivity or remove the option.', GUIDE_CREDIT, GUIDE_URL))
    if 'Title Card Required' in text:
        found.append(Diagnosis('INFO', 'Title Card Required may be literal text',
            'Gaussian may print this as an input title, not necessarily an error.',
            'Look for explicit error lines before diagnosing missing input sections.'))
    if re.search(r'Convergence failure|SCF has not converged|SCF failed to converge', text, re.I):
        found.append(Diagnosis('ERROR', 'Electronic convergence problem',
            'A self-consistent-field convergence failure was reported.',
            'Inspect the initial guess, spin state, method and convergence settings.', GUIDE_CREDIT, GUIDE_URL))
    if re.search(r'Optimization stopped|Number of steps exceeded', text, re.I):
        found.append(Diagnosis('WARNING', 'Optimization step limit reached',
            'The optimizer did not satisfy its convergence criteria.',
            'Review geometry and gradients before increasing the cycle limit.', GUIDE_CREDIT, GUIDE_URL))
    if 'segmentation violation' in text.lower() or 'segmentation fault' in text.lower():
        found.append(Diagnosis('ERROR', 'Segmentation fault reported',
            'A segmentation-violation traceback is present; it could follow an earlier failure.',
            'Fix the earliest explicit Gaussian error first. If a valid minimal input still crashes, examine the installation.'))
    if re.search(r'Error termination via Lnk1e', text, re.I):
        found.append(Diagnosis('ERROR', 'Gaussian error termination',
            'A Gaussian link reported unsuccessful termination.',
            'Inspect the first specific error (above the termination message).'))
    return found
