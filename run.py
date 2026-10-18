#!/usr/bin/env python3
"""One entry point for native Qt GUI, validation, synthetic tests, real runs.

All real Gaussian executions require explicit --run-real CASE. Unit tests and
simulations never launch Gaussian or modify existing Gaussian outputs.
"""
from __future__ import annotations
import argparse
import importlib.util
import logging
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time
import unittest

from gaumon import __version__
from gaumon.engine import Monitor
from gaumon.inputs import BY_KEY, CASES, ADVANCED_CASES, ALL_CASES, load_case, write_case, validate_input
from gaumon.diagnostics import diagnose_log
from gaumon_ext.qtcheck import qt_environment_report
from gaumon_ext.checker import check_input as detailed_input_check
from gaumon_ext.analytics import IncrementalLog, LINK_DESC

LOG = logging.getLogger('gaumon')
ROOT = Path(__file__).resolve().parent


def env_check() -> int:
    ok, report=qt_environment_report()
    print('GauMon Desktop '+__version__)
    print(report)
    return 0 if ok else 2


def check_input(path: Path) -> int:
    if not path.is_file():
        LOG.error('Input not found: %s',path);return 2
    issues=validate_input(path.read_text(encoding='utf-8',errors='replace'))
    if not issues:
        print(f'PASS: common input sections present: {path}')
        return 0
    for issue in issues:print('INPUT PROBLEM:',issue)
    return 1


def list_cases() -> int:
    for case in ALL_CASES:
        print(f'{case.key:8} {case.filename:32} {case.formula:5} {case.description}')
    return 0


def export_cases(output_dir: Path) -> int:
    errors=0
    for case in ALL_CASES:
        try:print(f'SAVED: {write_case(case.key,output_dir)}')
        except Exception as exc:
            LOG.error('%s: %s',case.filename,exc);errors+=1
    return 1 if errors else 0


def wait_for_case(monitor: Monitor, timeout: float, expected: str | None) -> int:
    deadline=time.monotonic()+timeout
    prev=None
    while time.monotonic()<deadline:
        snap=monitor.poll()
        stage=snap['state']
        identity=(stage,snap['step'],len(snap['scf']),snap['freq_blocks'] if 'freq_blocks' in snap else snap['frequencies'])
        if identity!=prev:
            print(f'[{snap["elapsed"]:5.1f}s] state={stage:<11} step={str(snap["step"]):<4} '
                  f'SCF={len(snap["scf"]):<3} FREQ={snap["frequencies"]} '
                  f'IMAG={snap["imaginary_modes"]} {snap["detail"]}',flush=True)
            prev=identity
        if stage in Monitor.TERMINAL:
            print(f'FINAL: {stage}; elapsed observation {snap["elapsed"]:.1f}s; '
                  f'energy={snap["energy"]}; opt_completed={snap["opt_completed"]}; '
                  f'frequencies={snap["frequencies"]}; imaginary_modes={snap["imaginary_modes"]}; '
                  f'log={snap["log"]}')
            return 0 if stage=='COMPLETED' else 1
        time.sleep(.12)
    if snap['owned']:
        try:monitor.stop(force=False)
        except Exception:pass
    LOG.error('Timed out after %.1f seconds; no successful completion reported',timeout)
    return 124


def execute_case(key: str, output_dir: Path | None, *, real: bool, executable: str, timeout: float) -> int:
    if key not in BY_KEY:
        raise ValueError(f'Unknown case: {key}')
    if real:
        if not executable.strip():raise ValueError('Gaussian executable must not be empty')
        executable=shutil.which(executable) or executable
        if '/' in executable and not Path(executable).is_file():
            raise FileNotFoundError(executable)
        output_dir=Path(output_dir or Path.cwd()/'gaumon_real_results').resolve()
        output_dir.mkdir(parents=True,exist_ok=True)
        input_path=write_case(key,output_dir)
        LOG.warning('Launching a real licensed Gaussian calculation. It may take time and consume resources.')
    else:
        output_dir=Path(output_dir).resolve() if output_dir else Path(tempfile.mkdtemp(prefix='gaumon_mock_'))
        output_dir.mkdir(parents=True,exist_ok=True)
        input_path=write_case(key,output_dir)
        executable=str(ROOT/'examples'/'mock_g16')
        LOG.info('Running synthetic simulation; NO scientific calculation is performed.')
    monitor=Monitor()
    snap=monitor.launch(executable,str(input_path),str(output_dir),
                        f'{key}_{"real" if real else "synthetic"}.log')
    print(f'INPUT: {input_path}\nLOG: {snap["log"]}\nPID: {snap["pid"]}')
    return wait_for_case(monitor,timeout,key)


def main(argv=None) -> int:
    cli=argparse.ArgumentParser(description='GauMon standalone CPU monitor, tests, and Gaussian case validator')
    modes=cli.add_mutually_exclusive_group()
    modes.add_argument('--gui',action='store_true',help='launch Qt standalone desktop (default)')
    modes.add_argument('--demo-ui',action='store_true',help='Qt desktop and mock opt+freq run (no Gaussian)')
    modes.add_argument('--test',action='store_true',help='run headless automated tests')
    modes.add_argument('--env-check',action='store_true',help='verify Python/Conda environment')
    modes.add_argument('--list-cases',action='store_true',help='show seven basic/advanced Gaussian test templates')
    modes.add_argument('--export-cases',type=Path,metavar='DIR',help='save all seven .gjf templates')
    modes.add_argument('--validate',type=Path,metavar='INPUT',help='validate Gaussian input title/charge/geometry')
    modes.add_argument('--check-input',type=Path,metavar='INPUT',help='show line/column diagnostic checks')
    modes.add_argument('--timing-report',type=Path,metavar='LOG',help='measured Gaussian Link timing table')
    modes.add_argument('--warnings-report',type=Path,metavar='LOG',help='warning/error evidence with line numbers')
    modes.add_argument('--diagnose',type=Path,metavar='LOG',help='diagnose common Gaussian errors')
    modes.add_argument('--attach',type=Path,metavar='LOG',help='open GUI attached to an existing log from anywhere')
    modes.add_argument('--inspect',type=Path,metavar='LOG',help='headless summary and ETA for an existing log')
    modes.add_argument('--list-jobs',action='store_true',help='show visible Gaussian processes for this user')
    modes.add_argument('--deep-report',type=Path,metavar='LOG',help='inspect SCF-iteration blocks without GUI')
    modes.add_argument('--qcnr-report',type=Path,metavar='LOG',help='inspect QCNR / LinEq1 linear-solver iteration history')
    modes.add_argument('--freq-report',type=Path,metavar='LOG',help='inspect Gaussian vibrational frequency modes')
    modes.add_argument('--stress-case',choices=('many_steps','oscillation','incomplete','chunked'),help='parse a bundled synthetic stress log')
    cli.add_argument('--pid',type=int,help='external Gaussian PID to validate with --attach or --inspect')
    modes.add_argument('--mock-case',choices=BY_KEY.keys(),help='simulate one case with mock executable')
    modes.add_argument('--run-real',choices=BY_KEY.keys(),help='explicitly launch licensed Gaussian test')
    cli.add_argument('--output',type=Path,help='output directory for mock or real case')
    cli.add_argument('--gaussian-exe',default='g16',help='real Gaussian executable or absolute path')
    cli.add_argument('--timeout',type=float,default=120.0,help='test timeout in seconds (real jobs may require more)')
    cli.add_argument('--debug',action='store_true',help='verbose logs and traceback')
    cli.add_argument('--version',action='version',version=__version__)
    a=cli.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if a.debug else logging.INFO,
        format='%(asctime)s | %(levelname)-7s | %(message)s')
    try:
        if a.timeout<=0:cli.error('--timeout must be >0')
        if a.env_check:return env_check()
        if a.list_cases:return list_cases()
        if a.list_jobs:
            import json
            print(json.dumps(Monitor().jobs(),indent=2))
            return 0
        if a.export_cases:return export_cases(a.export_cases)
        if a.validate:return check_input(a.validate)
        if a.check_input:
            if not a.check_input.is_file():raise FileNotFoundError(a.check_input)
            issues=detailed_input_check(a.check_input.read_text(errors='replace'))
            for item in issues:
                print(f'{item.severity.upper()} line {item.line}:{item.column} — {item.message}\n    {item.fix}')
            return 1 if any(x.severity=='error' for x in issues) else 0
        if a.timing_report or a.warnings_report:
            filename=a.timing_report or a.warnings_report
            if not filename.is_file():raise FileNotFoundError(filename)
            reader=IncrementalLog(filename)
            while reader.offset<filename.stat().st_size:
                prior=reader.offset
                reader.poll(max_bytes=1048576)
                if reader.offset==prior:break
            snap=reader.snapshot
            if a.timing_report:
                print(f'{filename.name}: {snap.event_count} completed Link calls')
                print('Link | Calls | Elapsed s | CPU s | Purpose')
                for lid,(elapsed,cpu,n) in sorted(snap.totals().items(),key=lambda pair:pair[1][0],reverse=True):
                    print(f'L{lid:<4} {n:>7} {elapsed:>12.2f} {cpu:>10.2f}  {LINK_DESC.get(lid,"Other Gaussian stage")}')
                print('Unfinished stages and per-atom runtimes cannot be measured from these records.')
            else:
                print(f'{filename.name}: {len(snap.problems)} recognized warnings/errors')
                for issue in snap.problems:
                    print(f'[{issue.level}] line {issue.line}: {issue.title}\n    {issue.evidence}\n    {issue.explanation}')
            return 0
        if a.diagnose:
            if not a.diagnose.is_file():raise FileNotFoundError(a.diagnose)
            issues=diagnose_log(a.diagnose.read_text(errors='replace')[-200000:])
            if issues:
                for item in issues:
                    print(f'[{item.severity}] {item.title}\n  Evidence: {item.description}\n  Action: {item.action}')
                    if item.source_url:
                        print(f'  Credit: {item.source_name} — {item.source_url}')
            else:print('No recognized diagnostics; this does NOT verify scientific correctness.')
            return 0
        if a.pid is not None and (a.pid<=0 or (not a.attach and not a.inspect)):
            cli.error('--pid requires a positive integer and --attach or --inspect')
        if a.deep_report or a.stress_case:
            from gaumon.engine import Parser
            filename=(ROOT/'examples'/'stress_logs'/f'{a.stress_case}.log') if a.stress_case else a.deep_report
            if not filename.is_file():
                raise FileNotFoundError(str(filename))
            parser=Parser()
            with filename.open('r',encoding='utf-8',errors='replace') as f:
                while chunk:=f.read(1048576):
                    parser.feed(chunk)
            # Preserve partial line, if any, for summary only: a live log may be
            # mid-write, so don't fabricate fields from unfinished text.
            metadata=parser.deep.overview()
            print(f'Deep SCF report: {filename.name}  |  {len(metadata["blocks"])} procedures')
            print('Procedure | optimizer step | status | SCF iterations | final energy (Ha) | last DIIS error')
            for b in metadata['blocks']:
                data=parser.deep.get_block(b['id'])
                last=data['rows'][-1] if data['rows'] else {}
                print(f'{b["id"]:>9} | {str(b["optimization_step"] or "—"):>14} | {b["status"]:<8} | '
                      f'{b["cycles"]:>14} | {str(data["final_energy"] or "—"):>17} | {last.get("diis_error")}')
            return 0
        if a.qcnr_report or a.freq_report:
            from gaumon.engine import Parser
            filename=a.qcnr_report if a.qcnr_report else a.freq_report
            if not filename.is_file():raise FileNotFoundError(filename)
            parser=Parser()
            with filename.open('r',encoding='utf-8',errors='replace') as stream:
                while chunk:=stream.read(1048576):parser.feed(chunk)
            if a.qcnr_report:
                overview=parser.qcnr.overview()
                print(f'QCNR / LinEq1 report: {filename.name}; outer QCNR iterations: {overview["outer_iterations"]}; solves: {len(overview["solves"])}')
                print('Solve | QCNR outer | status | LinEq1 rows | final residual RMS | final residual max')
                for entry in overview['solves']:
                    block=parser.qcnr.get_solve(entry['id'])
                    end=block['rows'][-1] if block['rows'] else {}
                    print(f'{entry["id"]:>5} | {str(entry["outer_iteration"]):>10} | {entry["status"]:<11} | '
                          f'{len(block["rows"]):>11} | {str(end.get("rms")):>18} | {end.get("max")}')
                print('Note: LinEq1 is an SCF linear equation solve, not a vibrational frequency.')
            else:
                overview=parser.freq_deep.overview()
                print(f'Deep Frequencies report: {filename.name}; printed groups: {len(overview["groups"])}')
                for entry in overview['groups']:
                    group=parser.freq_deep.get_group(entry['id'])
                    masses=sum(row['red_mass'] is not None for row in group['modes'])
                    print(f'Group {entry["id"]}: {entry["modes"]} modes, {entry["imaginary"]} imaginary, reduced masses for {masses} modes')
                if not overview['groups']:
                    print('No Frequencies -- records. LinEq1 / QCNR output is not frequency data.')
            return 0
        if a.inspect:
            import json
            monitor=Monitor()
            monitor.attach(a.inspect.resolve(),a.pid)
            # Initial file processing; do not report a fake historical runtime.
            state=monitor.poll()
            info={key:state.get(key) for key in (
                'state','detail','step','step_max','energy','cycles','link',
                'link_description','scf','step_history','link_history','eta_records',
                'eta_to_max','cycle_sample_count','cycle_timing_source')}
            info['diagnoses']=[vars(d) for d in diagnose_log(state.get('diagnostic_text','')+'\n'+state['tail'])]
            print(json.dumps(info,indent=2))
            return 0
        if a.mock_case:return execute_case(a.mock_case,a.output,real=False,executable='',timeout=a.timeout)
        if a.run_real:return execute_case(a.run_real,a.output,real=True,executable=a.gaussian_exe,timeout=a.timeout)
        if a.test:
            suite=unittest.defaultTestLoader.discover(str(ROOT/'tests'),pattern='test_*.py')
            result=unittest.TextTestRunner(verbosity=2).run(suite)
            return 0 if result.wasSuccessful() else 1
        try:
            from gaumon.ui import run_gui
        except ImportError as exc:
            LOG.error('Qt GUI dependencies unavailable: %s\nRun: python run.py --env-check',exc)
            return 2
        return run_gui(mock_on_start=a.demo_ui, attach_log=a.attach.resolve() if a.attach else None, attach_pid=a.pid)
    except (OSError,ValueError,RuntimeError) as exc:
        LOG.exception('GauMon operation failed') if a.debug else LOG.error('%s',exc)
        return 2

if __name__=='__main__':sys.exit(main())
