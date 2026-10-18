#!/usr/bin/env python3
"""Generate synthetic high-difficulty Gaussian-style SCF logs (NO chemistry)."""
from __future__ import annotations
import argparse
import logging
import math
from pathlib import Path

LOG=logging.getLogger('gaumon.stress')


def cycle(i: int, offset: float, *, swing: bool = False) -> str:
    de=(-0.09*math.exp(-i/5.0)) + ((-1)**i*0.04*math.exp(-i/18.0) if swing else 0)
    energy=offset + (-0.21*(1-math.exp(-i/4.5)))+(0.015*math.sin(i*1.8) if swing else 0)
    err=(0.035*math.exp(-i/13.0))*(3 if swing and i%7==0 else 1)
    lines=[f' Cycle {i:>3}  Pass {0 if i<4 else 1}  IDiag  1:',
           f' RMSU=  {0.004*math.exp(-i/11):.3E} CP: 9.56D-01 1.03D-02 7.27D-01',
           '                    CP:  8.35D-01',
           f' E= {energy:.14f} Delta-E= {de:.11f} Rises={"T" if de>0 else "F"} Damp=F',
           f' DIIS: error= {err:.3E} at cycle {i:3} NSaved= {min(9,i):3}.',
           f' NSaved= {min(9,i)} IEnMin= {i} EnMin= {energy:.14f} IErMin= {i} ErrMin= {err:.3E}',
           f' ErrMax= {err:.3E} 0.00D+00 EMaxC= 1.00D-01 BMatC= 3.03D-01 BMatP= 3.37D+00',
           ' IDIUse=3 WtCom= 8.97D-01 WtEn= 1.03D-01',
           ' Coeff-Com: -0.165D+00 0.225D+00 0.274D-01-0.804D-01',
           ' Coeff-Com:  0.767D+00',
           ' Coeff-En:  0.000D+00 0.000D+00 0.100D+01',
           ' Coeff: -0.148D+00 0.202D+00 0.246D-01-0.721D-01',
           ' Gap= 0.330 Goal= None Shift= 0.000',
           ' Gap= 0.286 Goal= None Shift= 0.000',
           f' RMSDP={0.006*math.exp(-i/10):.3E} MaxDP={0.15*math.exp(-i/8):.3E} DE={de:.3E} OVMax= 5.82D-02']
    return '\n'.join(lines)+'\n'


def build(case: str) -> str:
    if case not in {'many_steps','oscillation','incomplete','chunked'}:
        raise ValueError(f'Unknown stress case {case}')
    out=['GauMon synthetic SCF stress case; NOT Gaussian output\n']
    if case=='many_steps':
        for step in range(1,9):
            out.append(' (Enter /mock/g16/l502.exe)\n')
            for i in range(1,10+step):
                out.append(cycle(i,-100.0-step))
            out.append(f' SCF Done: E(RHF) = -{100.0+step:.10f} A.U. after {9+step} cycles\n')
            out.append(' (Enter /mock/g16/l103.exe)\n')
            out.append(f' Step number {step} out of a maximum of 40\n')
        out.append(' Normal termination of Gaussian 16 at mock time\n')
    elif case=='oscillation':
        out.append(' (Enter /mock/g16/l502.exe)\n')
        for i in range(1,96):
            out.append(cycle(i,-345.0,swing=True))
        out.append(' Convergence failure -- run terminated.\nError termination via Lnk1e in /mock/g16/l502.exe\n')
    elif case=='incomplete':
        out.append(' (Enter /mock/g16/l502.exe)\n')
        for i in range(1,12):out.append(cycle(i,-4200.0))
        out.append(' Cycle  12 Pass 1 IDiag 1:\n RMSU=1.00D-04 CP: 9.5D-01\n')
        # Deliberate abrupt end without SCF Done.
    else:
        out.append(' (Enter /mock/g16/l502.exe)\n')
        for i in range(1,5):out.append(cycle(i,-50.0))
        out.append(' SCF Done: E(RHF) = -50.0000 A.U. after 4 cycles\n')
        out.append(' (Enter /mock/g16/l103.exe)\n Step number 1 out of a maximum of 25\n')
        out.append(' Cycle  1 Pass 1 IDiag 1:\n E= -50.1000 Delta-E= -0.10000')  # partial line
    return ''.join(out)


def main(argv=None):
    parser=argparse.ArgumentParser(description='Create stress logs without running Gaussian')
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'examples'/'stress_logs')
    parser.add_argument('--case',choices=('many_steps','oscillation','incomplete','chunked','all'),default='all')
    args=parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO,format='%(levelname)s: %(message)s')
    choices=('many_steps','oscillation','incomplete','chunked') if args.case=='all' else (args.case,)
    try:
        args.output.mkdir(parents=True,exist_ok=True)
        for key in choices:
            file=args.output/f'{key}.log'
            file.write_text(build(key),encoding='utf-8')
            LOG.info('Created synthetic %s, bytes=%s',file,file.stat().st_size)
    except OSError as exc:
        LOG.error('Failed to write fixture: %s',exc)
        return 1
    return 0

if __name__=='__main__':
    raise SystemExit(main())
