#!/usr/bin/env python3
"""Check and screenshot GauMon layouts on a computer with PyQt5 installed.

No Gaussian is executed, no GPU or internet is needed. Outputs screenshots
for QA and reports obvious clipped/overlapping cards and plots.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main() -> int:
    parser = argparse.ArgumentParser(description='Save native GauMon UI layout verification screenshots')
    parser.add_argument('--output', type=Path, default=Path('gaumon_gui_qa'))
    parser.add_argument('--offscreen', action='store_true', help='use Qt offscreen platform for CI / remote nodes')
    args = parser.parse_args()
    if args.offscreen:
        os.environ['QT_QPA_PLATFORM'] = 'offscreen'
    from PyQt5.QtCore import QPoint, QRect
    from PyQt5.QtWidgets import QApplication
    import pyqtgraph as pg
    from gaumon.ui import MainWindow
    pg.setConfigOptions(useOpenGL=False, antialias=True)
    app = QApplication.instance() or QApplication([])
    app.setStyle('Fusion')
    args.output.mkdir(parents=True, exist_ok=True)
    report = []
    for width, height in ((1366,768),(1920,1080)):
        window = MainWindow()
        window.show()
        window.resize(width,height)
        app.processEvents()
        issues = []
        bounds={}
        for name, card in window.plots.items():
            pos=card.mapTo(window, QPoint(0,0))
            rect=QRect(pos,card.size())
            bounds[name]=rect
            if rect.width() < 210 or rect.height() < 115:
                issues.append(f'Plot {name}: undersized ({rect.width()}x{rect.height()})')
        for a in bounds:
            for b in bounds:
                if a>=b:continue
                if bounds[a].intersects(bounds[b]):
                    issues.append(f'Overlap: {a} with {b}')
        for name,metric in window.metrics.items():
            if metric.width() < 115:
                issues.append(f'Metric {name}: undersized ({metric.width()}px)')
            # Compare current content against available label width only; fonts can vary.
            text=metric.value.text()
            if metric.value.fontMetrics().horizontalAdvance(text) > metric.value.width()+2:
                issues.append(f'Metric {name}: value text clipped')
        path=args.output/f'overview_{width}x{height}.png'
        if not window.grab().save(str(path)):
            issues.append('Qt screenshot could not be saved')
        window.show_page(2)
        app.processEvents()
        path2=args.output/f'timeline_{width}x{height}.png'
        window.grab().save(str(path2))
        window.open_deep()
        app.processEvents()
        deep_bounds={}
        for name,card in window.deep_charts.items():
            rect=QRect(card.mapTo(window,QPoint(0,0)),card.size())
            deep_bounds[name]=rect
            if rect.width()<190 or rect.height()<125:
                issues.append(f'Deep plot {name}: undersized ({rect.width()}x{rect.height()})')
        for a in deep_bounds:
            for b in deep_bounds:
                if a>=b:continue
                if deep_bounds[a].intersects(deep_bounds[b]):
                    issues.append(f'Deep overlap: {a} with {b}')
        path3=args.output/f'deep_scf_{width}x{height}.png'
        if not window.grab().save(str(path3)):
            issues.append('Deep Qt screenshot could not be saved')
        window.close_deep()
        window.open_qcnr()
        app.processEvents()
        for name,card in window.qcnr_plots.items():
            rect=QRect(card.mapTo(window,QPoint(0,0)),card.size())
            if rect.width()<185 or rect.height()<110:
                issues.append(f'QCNR plot {name}: undersized ({rect.width()}x{rect.height()})')
        path4=args.output/f'deep_qcnr_{width}x{height}.png'
        if not window.grab().save(str(path4)):
            issues.append('QCNR screenshot could not be saved')
        window.close_qcnr()
        window.open_freq()
        app.processEvents()
        path5=args.output/f'deep_frequencies_{width}x{height}.png'
        if not window.grab().save(str(path5)):
            issues.append('Frequency screenshot could not be saved')
        window.close_freq()
        extra_screens={}
        for index,title in [(10,'time_graphs'),(11,'input_checker'),(12,'advanced_eta'),
                            (13,'warnings_errors'),(14,'session_jobs')]:
            window.show_page(index)
            app.processEvents()
            path_extra=args.output/f'{title}_{width}x{height}.png'
            if not window.grab().save(str(path_extra)):
                issues.append(f'{title} screenshot could not be saved')
            extra_screens[title+'_screenshot']=str(path_extra)
        report.append({'resolution':[width,height], 'issues':issues,
                       **extra_screens, 'deep_qcnr_screenshot':str(path4),'deep_frequencies_screenshot':str(path5),
                       'overview_screenshot':str(path),'timeline_screenshot':str(path2),
                       'deep_scf_screenshot':str(path3)})
        window.close()
    dest=args.output/'layout_report.json'
    dest.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))
    print('Layout screenshots and report:',args.output.resolve())
    return 1 if any(r['issues'] for r in report) else 0


if __name__=='__main__':
    sys.exit(main())
