#!/usr/bin/env python3
"""Synthetic Gaussian-style stream for GauMon testing only. NOT quantum chemistry."""
from __future__ import annotations
import os
import re
import sys
import time


def emit(s: str = ''):
    print(s, flush=True)


def geometry(step: int, is_nh3: bool):
    z = 0.08 * step
    atoms = ([(7, 0.000, 0.000, 0.000), (1, 1.000, 0., z),
              (1, -0.500, .866, 0.), (1, -0.500, -.866, 0.)]
             if is_nh3 else
             [(8, 0., 0., 0.), (1, 0., .72, .55+z), (1, 0., -.72, .55-z)])
    emit(' Standard orientation:')
    emit(' ---------------------------------------------------------------------')
    emit(' Center     Atomic      Atomic             Coordinates (Angstroms)')
    emit(' Number     Number       Type             X           Y           Z')
    emit(' ---------------------------------------------------------------------')
    for i, (n, x, y, zz) in enumerate(atoms, start=1):
        emit(f' {i:>5d} {n:>10d}           0    {x:>12.6f} {y:>11.6f} {zz:>11.6f}')
    emit(' ---------------------------------------------------------------------')


def mock_scf_cycles(step: int, count: int, delay: float, *, tough=False):
    """Synthetic SCF cycles to exercise Deep SCF plots; not numerical chemistry."""
    import math
    for k in range(1, count+1):
        emit(f' Cycle {k:>3d}  Pass {0 if k < 3 else 1}  IDiag  1:')
        energy = -76.0 - 0.01*step - 0.04*(1-math.exp(-0.34*k))
        if tough:
            energy += 0.03*math.sin(k*1.7) * math.exp(-k/15)
        delta = -(0.04*0.34)*math.exp(-0.34*k)
        diis = 0.04*math.exp(-0.24*k)
        rmsu = 0.009*math.exp(-0.26*k)
        rmsdp = 0.008*math.exp(-0.28*k)
        maxdp = 0.15*math.exp(-0.20*k)
        emit(f' RMSU= {rmsu:.3E}    CP: 9.57D-01 1.12D-02 7.42D-01')
        emit(f' E= {energy: .12f} Delta-E= {delta: .9f} Rises=F Damp=F')
        emit(f' DIIS: error= {diis:.3E} at cycle {k:3d} NSaved= {min(k,9):3d}.')
        emit(f' NSaved= {min(k,9)} IEnMin= {k} EnMin= {energy: .12f} IErMin= {k} ErrMin= {diis:.3E}')
        emit(f' ErrMax= {diis:.3E} 0.00D+00 EMaxC= 1.00D-01 BMatC= 2.22D-02 BMatP= 1.17D+00')
        emit(' IDIUse=3 WtCom= 9.80D-01 WtEn= 2.00D-02')
        emit(' Coeff-Com: -0.165D+00 0.225D+00 0.274D-01-0.804D-01')
        emit(' Coeff-En:   0.000D+00 0.000D+00 0.100D+01')
        emit(' Coeff: -0.148D+00 0.202D+00 0.246D-01-0.721D-01')
        emit(' Gap= 0.284 Goal= None Shift= 0.000')
        emit(f' RMSDP={rmsdp:.3E} MaxDP={maxdp:.3E} DE={delta:.3E} OVMax= 5.82D-02')
        if delay > 0:
            time.sleep(delay / count)
    return energy


def main():
    text = sys.stdin.read()
    delay = float(os.environ.get('GAUMON_MOCK_STEP_SECONDS', '0.85'))
    lo = text.lower()
    case = ('trimer' if 'advanced test - hydrogen bonded water trimer' in lo else
            'stretched' if 'advanced test - stretched h2' in lo else
            'fe2s2' if 'advanced test - illustrative fe2s2' in lo else
            'ts' if 'opt=(ts' in lo else
            'optfreq' if 'opt freq' in lo else
            'freq' if re.search(r'#.*\bfreq\b', lo) else 'sp')
    bad = 'GAUMON_BAD_TITLE' in text
    emit(' Gaussian 16: MOCK TEST OUTPUT -- NOT A REAL GAUSSIAN CALCULATION')
    emit(' (Enter /mock/g16/l101.exe)')
    if bad:
        emit(' Title Card Required')
        emit(' Error termination via Lnk1e in /mock/g16/l101.exe')
        return 2
    emit(' Charge = 0 Multiplicity = 1')
    emit(' #p ' + ('b3lyp/3-21g sp' if case in ('sp','stretched','fe2s2') else 'b3lyp/3-21g opt freq'))
    steps = 7 if case == 'trimer' else 5 if case in ('optfreq','ts') else 1
    for j in range(1, steps+1):
        emit(' (Enter /mock/g16/l502.exe)')
        cycles = (24 if case in ('stretched','fe2s2') else 16 if case=='trimer' else 9)
        energy = mock_scf_cycles(j, cycles, delay, tough=case in ('trimer','stretched','fe2s2'))
        emit(f' SCF Done:  E(RB3LYP) =  {energy:.9f}     A.U. after   {cycles} cycles')
        if case in ('optfreq','ts','trimer'):
            emit(' (Enter /mock/g16/l103.exe)')
            emit(f' Step number {j} out of a maximum of 35')
        emit(f' RMS Force = {0.01/(j*j):.5E}')
        geometry(j, case == 'ts')
    if case in ('optfreq', 'ts','trimer'):
        emit(' Optimization completed.')
    if case in ('freq','optfreq','ts','trimer'):
        emit(' (Enter /mock/g16/l716.exe)')
        if case == 'ts':
            emit(' Frequencies -- -850.0000   450.0000   750.0000')
            emit(' Frequencies --  900.0000  1200.0000  1600.0000')
        else:
            emit(' Frequencies --  1600.0000   3700.0000  3800.0000')
        time.sleep(delay)
    emit(' Normal termination of Gaussian 16 at mock time')
    return 0

if __name__ == '__main__':
    sys.exit(main())
