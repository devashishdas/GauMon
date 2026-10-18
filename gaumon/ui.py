"""Standalone, CPU-only Qt dashboard for GauMon.

No Chromium/WebEngine, no HTTP server, no CUDA, and no OpenGL plot backend.
The dashboard itself has no scroll area; data tables and logs scroll internally.
"""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor
import math
import csv
import logging
import os
from pathlib import Path
import tempfile
import time
from typing import Optional

from PyQt5.QtCore import Qt, QTimer, QPointF, QUrl
from PyQt5.QtGui import QColor, QFont, QPainter, QPen, QBrush, QDesktopServices, QKeySequence
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QFrame, QLabel,
    QVBoxLayout, QHBoxLayout, QGridLayout, QStackedWidget, QLineEdit,
    QPlainTextEdit, QTableWidget, QTableWidgetItem, QAbstractItemView, QHeaderView,
    QFormLayout, QFileDialog, QMessageBox, QPushButton, QComboBox, QCheckBox, QDoubleSpinBox,
    QButtonGroup, QSplitter, QSizePolicy, QInputDialog, QShortcut)
import pyqtgraph as pg

from .engine import Monitor
from .inputs import CASES, ADVANCED_CASES, ALL_CASES, BY_KEY, validate_input, write_case, load_case
from .diagnostics import diagnose_log, GUIDE_CREDIT, GUIDE_URL
from gaumon_ext.gui import AnalysisHub, TimePage, InputPage, ETAPage, WarningsPage, JobsPage
from gaumon_ext.input_context import matching_input_for_log

LOG = logging.getLogger('gaumon.ui')

C = {
    'bg': '#F6F7F3', 'panel': '#FFFFFF', 'sidebar': '#FCFCF9', 'text': '#24352D',
    'sub': '#6D7F75', 'muted': '#94A49A', 'border': '#E4EBE4',
    'green': '#28966C', 'lightgreen': '#E5F4EA', 'orange': '#DC9450',
    'coral': '#D76567', 'violet': '#9872AF', 'ink': '#23372F',
}

STYLE = """
QMainWindow, QWidget#root { background: #F6F7F3; }
QWidget { color: #24352D; font-family: 'Noto Sans', 'DejaVu Sans', sans-serif; font-size: 13px; }
QFrame#panel, QFrame#metric, QFrame#header, QFrame#side, QFrame#plotCard, QFrame#subPanel {
 background: #FFFFFF; border: 1px solid #E4EBE4; border-radius: 12px;
}
QFrame#side { background: #FCFCF9; }
QLabel#brand { font-size: 23px; font-weight: 800; color: #23382E; }
QLabel#small { font-size: 11px; color: #6D7F75; }
QLabel#eyebrow { font-size: 11px; font-weight: 700; color: #28966C; }
QLabel#section { font-size: 15px; font-weight: 700; color: #23382E; }
QLabel#cardTitle { font-size: 11px; font-weight: 700; color: #6D7F75; }
QLabel#cardValue { font-size: 20px; font-weight: 800; color: #24352D; }
QLabel#detailKey { font-size: 12px; color: #718278; }
QLabel#detailValue { font-size: 12px; font-weight: 650; color: #324438; }
QLabel#pill { border-radius: 9px; background: #E5F4EA; color: #21885E;
 padding: 7px 11px; font-size: 11px; font-weight: 750; }
QPushButton { background: #F7FAF7; border: 1px solid #DFE8DE; border-radius: 8px;
 padding: 7px 12px; min-height: 18px; color: #354A3B; font-weight: 650; }
QPushButton:hover { border-color: #94CDB1; background: #F0F7F0; }
QPushButton:disabled { color: #ACB9B1; background: #F7F8F7; border-color: #EBEFEB; }
QPushButton#primary { background: #28966C; color: white; border-color: #28966C; }
QPushButton#primary:hover { background: #207C5A; }
QPushButton#danger { background: #FFF5F4; color: #C65B5D; border-color: #F2D8D8; }
QPushButton#nav { text-align: left; padding: 5px 10px; min-height: 16px; background: transparent;
 border: none; font-weight: 650; color: #67786D; }
QPushButton#nav:hover { background: #F0F6F0; }
QPushButton#nav:checked { background: #E6F4EC; color: #1C865D; }
QLineEdit, QComboBox, QDoubleSpinBox { background: #FFFFFF; color: #24352D; border: 1px solid #DDE5DC;
 border-radius: 8px; padding: 8px 10px; selection-background-color: #BDE4CE; }
QLineEdit:focus, QComboBox:focus { border-color: #28966C; }
QPlainTextEdit { background: #FBFCFA; color: #2E4735; border: 1px solid #E4EBE4;
 border-radius: 9px; padding: 8px; font-family: 'DejaVu Sans Mono',monospace; font-size: 11px; }
QTableWidget { background: white; border: 1px solid #E4EBE4; border-radius: 9px;
 gridline-color: #EFF2ED; alternate-background-color: #FAFCFA; selection-background-color: #E5F4EA;
 selection-color: #214332; }
QHeaderView::section { background: #F6F9F6; border: none; border-bottom: 1px solid #E4EBE4;
 padding: 8px 9px; font-size: 11px; font-weight: 700; color: #64776A; }
QCheckBox { spacing: 6px; }
QFrame#insight { background: #EFF7F0; border: 1px solid #D5EBDC; border-radius: 11px; }
QFrame#warningBanner { background: #FFF6ED; border: 1px solid #F4DBC7; border-radius: 11px; }
QLabel#calloutTitle { font-size: 12px; font-weight: 750; color: #276644; }
QLabel#calloutValue { font-size: 22px; font-weight: 800; color: #24352D; }
QLabel#sectionHeading { font-size: 14px; font-weight: 750; color: #24352D; }
QPushButton#subtle { background: #FFFFFF; border: 1px solid #E1EAE1; }

"""


def row_layout(gap: int = 8) -> QHBoxLayout:
    layout = QHBoxLayout(); layout.setSpacing(gap); layout.setContentsMargins(0,0,0,0)
    return layout


def col_layout(gap: int = 8) -> QVBoxLayout:
    layout = QVBoxLayout(); layout.setSpacing(gap); layout.setContentsMargins(0,0,0,0)
    return layout


def label(text: str, style: str = 'small', *, wrap: bool = False) -> QLabel:
    result = QLabel(text); result.setObjectName(style); result.setWordWrap(wrap)
    return result


class Metric(QFrame):
    def __init__(self, caption: str, initial: str = '—'):
        super().__init__(); self.setObjectName('metric')
        self.setMinimumWidth(0)
        v = QVBoxLayout(self); v.setContentsMargins(13, 10, 10, 9); v.setSpacing(3)
        self.caption = label(caption.upper(), 'cardTitle')
        self.value = label(initial, 'cardValue')
        self.value.setTextInteractionFlags(Qt.TextSelectableByMouse)
        v.addWidget(self.caption); v.addWidget(self.value)


class PlotCard(QFrame):
    def __init__(self, title: str, color: str, *, is_log: bool = False, paired: bool = False):
        super().__init__(); self.setObjectName('plotCard')
        wrap = col_layout(4); wrap.setContentsMargins(12, 9, 12, 10); self.setLayout(wrap)
        name = label(title, 'detailValue'); wrap.addWidget(name)
        self.plot = pg.PlotWidget(background='w')
        self.plot.setMinimumSize(180, 115)
        self.plot.setMouseEnabled(x=True,y=True)
        self.plot.showGrid(x=True,y=True,alpha=0.13)
        self.plot.getPlotItem().setMenuEnabled(False)
        for axis in ('left','bottom'):
            ax = self.plot.getAxis(axis)
            ax.setPen(pg.mkPen(C['border']))
            ax.setTextPen(pg.mkPen(C['sub']))
            ax.setStyle(tickFont=QFont('DejaVu Sans', 8), tickLength=-3)
        if is_log:
            self.plot.setLogMode(y=True)
        self.curve = self.plot.plot(pen=pg.mkPen(color,width=2.4), symbol='o',
                                    symbolSize=4, symbolBrush=color, symbolPen=None)
        self.curve2 = self.plot.plot(pen=pg.mkPen(C['violet'], width=2.1),
                                     symbol=None) if paired else None
        wrap.addWidget(self.plot,stretch=1)

    def refresh(self, data, second=None):
        if data:
            self.curve.setData([v[0] for v in data], [v[1] for v in data])
        else:
            self.curve.setData([],[])
        if self.curve2 is not None:
            if second: self.curve2.setData([v[0] for v in second], [v[1] for v in second])
            else: self.curve2.setData([],[])


class MoleculeCanvas(QWidget):
    """2D XY projection only; not a 3D molecular renderer."""
    def __init__(self):
        super().__init__(); self.atoms = []; self.setMinimumSize(240,180)

    def set_atoms(self, atoms):
        self.atoms = atoms or []; self.update()

    def paintEvent(self, event):
        p=QPainter(self); p.setRenderHint(QPainter.Antialiasing,True)
        p.fillRect(self.rect(), QColor('#FAFCF9'))
        if not self.atoms:
            p.setPen(QColor(C['sub'])); p.drawText(self.rect(), Qt.AlignCenter,
                'Geometry will appear here\nafter the first coordinate snapshot')
            return
        atoms=self.atoms[:200]
        xs=[a['x'] for a in atoms]; ys=[a['y'] for a in atoms]
        cx=(max(xs)+min(xs))/2;cy=(max(ys)+min(ys))/2
        extent=max(max(xs)-min(xs),max(ys)-min(ys),1.8)
        scale=min(self.width()-68,self.height()-55)/extent
        def point(a):
            return QPointF(self.width()/2+(a['x']-cx)*scale,
                           self.height()/2-(a['y']-cy)*scale)
        from math import dist
        palette={1:('#FFFFFF','H'),6:('#525D58','C'),7:('#8C6EB0','N'),8:('#DB7068','O'),
                 16:('#D9BB58','S'),26:('#A57B53','Fe')}
        for i,a in enumerate(atoms):
            for b in atoms[i+1:]:
                if dist((a['x'],a['y'],a['z']),(b['x'],b['y'],b['z'])) < (1.27 if 1 in (a['element'],b['element']) else 1.92):
                    p.setPen(QPen(QColor('#C1CDC3'),3,Qt.SolidLine,Qt.RoundCap))
                    p.drawLine(point(a),point(b))
        p.setFont(QFont('DejaVu Sans',8,QFont.Bold))
        for a in atoms:
            pt=point(a); symbol=palette.get(a['element'],('#B79E81',str(a['element'])))
            radius=12 if a['element']!=1 else 8
            p.setPen(QPen(QColor('#C8D2C9'),1));p.setBrush(QBrush(QColor(symbol[0])))
            p.drawEllipse(pt,radius,radius)
            p.setPen(QColor('#38483E') if a['element']==1 else QColor('#FFFFFF'))
            p.drawText(int(pt.x()-radius),int(pt.y()-radius),radius*2,radius*2,Qt.AlignCenter,symbol[1])


class MainWindow(QMainWindow):
    def __init__(self, *, mock_on_start=False):
        super().__init__()
        self.monitor=Monitor(); self.mock_dir: Optional[tempfile.TemporaryDirectory]=None
        self.last_session=None
        self.last_plot_sig=None
        self.last_energy_sig=None
        self.last_cycles_sig=None
        self.last_links_sig=None
        self.last_deep_sig=None
        self._deep_page_open=False
        self._deep_history_selection='follow'
        self._deep_block_for_export=None
        self.last_qcnr_sig=None
        self._qcnr_page_open=False
        self._qcnr_selection='follow'
        self._qcnr_export=None
        self.last_freq_sig=None
        self._freq_page_open=False
        self._freq_selection='follow'
        self._freq_export=None
        # This is a record of ONLY explicitly enrolled jobs. The store is loaded
        # at startup; saved logs are NOT opened or restarted without user action.
        self.hub=AnalysisHub(self)
        self._pending_enroll_log=None
        self._pending_enroll_input=None
        self._poll_executor=ThreadPoolExecutor(max_workers=1,thread_name_prefix='gaumon-log')
        self._poll_future=None
        self.disp_limit=8.0
        self.rmsd_limit=2.0
        self.setWindowTitle('GauMon Desktop')
        self.setMinimumSize(1160, 680)
        available=QApplication.primaryScreen().availableGeometry()
        self.resize(min(1400,available.width()),min(820,available.height()))
        self.setStyleSheet(STYLE)
        self._build()
        self.work_entry.setText(os.getcwd())
        self.shortcut_full=QShortcut(QKeySequence('F11'),self)
        self.shortcut_full.activated.connect(self.toggle_fullscreen)
        self.shortcut_exit=QShortcut(QKeySequence('Escape'),self)
        self.shortcut_exit.activated.connect(self.exit_fullscreen)
        self.timer=QTimer(self);self.timer.setInterval(750);self.timer.timeout.connect(self.tick);self.timer.start()
        if mock_on_start: QTimer.singleShot(200,lambda:self.simulate('optfreq'))

    def toggle_fullscreen(self):
        if self.isFullScreen():self.showMaximized()
        else:self.showFullScreen()
        self.full_btn.setText('Exit full screen  Esc' if self.isFullScreen() else 'Full screen  F11')

    def exit_fullscreen(self):
        if self.isFullScreen():
            self.showMaximized()
            self.full_btn.setText('Full screen  F11')

    def _panel(self) -> tuple[QFrame,QVBoxLayout]:
        panel=QFrame();panel.setObjectName('panel')
        v=QVBoxLayout(panel);v.setSpacing(10);v.setContentsMargins(14,12,14,12)
        return panel,v

    def _build(self):
        root=QWidget();root.setObjectName('root');self.setCentralWidget(root)
        outer=QVBoxLayout(root);outer.setContentsMargins(14,12,14,8);outer.setSpacing(10)
        header=QFrame();header.setObjectName('header'); header.setFixedHeight(59)
        hl=QHBoxLayout(header);hl.setContentsMargins(16,8,16,8);hl.setSpacing(10)
        hl.addWidget(label('GauMon','brand'))
        hl.addWidget(label('Computational chemistry monitor','small'))
        hl.addStretch(1)
        self.header_state=label('● READY','pill'); hl.addWidget(self.header_state)
        hl.addWidget(label('Live calculation insights','small'))
        self.full_btn=QPushButton('Full screen  F11');self.full_btn.setObjectName('subtle')
        self.full_btn.setToolTip('F11 switches to immersive full-screen view. Press Esc to restore the maximized window.')
        self.full_btn.clicked.connect(self.toggle_fullscreen);hl.addWidget(self.full_btn)
        outer.addWidget(header)
        self.failure_banner=QFrame();self.failure_banner.setObjectName('warningBanner')
        alert=row_layout(9);alert.setContentsMargins(14,7,14,7);self.failure_banner.setLayout(alert)
        self.failure_title=label('Calculation diagnostic','detailValue')
        alert.addWidget(self.failure_title,1)
        see_failure=QPushButton('Inspect failure  →');see_failure.setObjectName('subtle')
        see_failure.clicked.connect(lambda:self.show_page(5))
        alert.addWidget(see_failure)
        self.failure_banner.setVisible(False)
        outer.addWidget(self.failure_banner)
        metrics=row_layout(8)
        self.metrics={}
        for key,title in [('state','Calculation'),('elapsed','Elapsed'),('energy','Last SCF energy  ·  Ha'),
                          ('step','Optimization step'),('force','RMS force'),('rmsd','Aligned RMSD  ·  Å')]:
            card=Metric(title);metrics.addWidget(card,1); self.metrics[key]=card
        self.metric_bar=QWidget()
        self.metric_bar.setLayout(metrics)
        outer.addWidget(self.metric_bar)
        content=row_layout(11)
        side=QFrame();side.setObjectName('side');side.setFixedWidth(210)
        side_v=col_layout(4);side_v.setContentsMargins(9,13,9,10);side.setLayout(side_v)
        side_v.addWidget(label('WORKSPACE','eyebrow'))
        self.navgroup=QButtonGroup(self);self.navgroup.setExclusive(True)
        self.pages=QStackedWidget()
        names=['Overview','Run calculation','Timeline & ETA','Analysis','Geometry','Live log','Test cases','Deep SCF','Deep QCNR','Deep Frequencies','Time graphs','Input checker / editor','Advanced time estimator','Warnings & errors','Session jobs']
        self.nav_buttons=[]
        for index,name in enumerate(names):
            button=QPushButton(name);button.setObjectName('nav');button.setCheckable(True)
            button.setCursor(Qt.PointingHandCursor)
            button.clicked.connect(lambda checked,ind=index:self.show_page(ind))
            self.navgroup.addButton(button);side_v.addWidget(button);self.nav_buttons.append(button)
            if name in ('Deep SCF','Deep QCNR','Deep Frequencies'):
                button.setVisible(False)
        self.nav_buttons[0].setChecked(True)
        side_v.addStretch(1)
        side_v.addWidget(label('Local desktop · CPU renderer', 'small',wrap=True))
        content.addWidget(side)
        content.addWidget(self.pages,1)
        outer.addLayout(content,1)
        foot=row_layout()
        self.footer=label('Ready · Choose an input, attach a log, or run a mock test','small');foot.addWidget(self.footer)
        foot.addStretch(1)
        foot.addWidget(label('Renderer: CPU / Qt','small'))
        outer.addLayout(foot)
        self._overview()
        self._run_page()
        self._timeline_page()
        self._analysis_page()
        self._geometry_page()
        self._log_page()
        self._tests_page()
        self._deep_page()
        self._qcnr_page()
        self._freq_page()
        self._utility_pages()

    def _utility_pages(self):
        """Port the five v6.1 utilities into the existing v5.2 native window."""
        self.time_page=TimePage(self.hub)
        self.input_checker=InputPage(self.hub)
        self.advanced_eta=ETAPage(self.hub)
        self.warning_page=WarningsPage(self.hub)
        self.jobs_page=JobsPage(self.hub,self._job_selected)
        for page in (self.time_page,self.input_checker,self.advanced_eta,
                     self.warning_page,self.jobs_page):
            self.pages.addWidget(page)

    def _job_selected(self):
        """Show selected enrolled log everywhere. Reattach remains read-only."""
        job=self.hub.current()
        if job is None:
            return
        path=Path(job.path)
        if path.resolve() != (self.monitor.log_path.resolve() if self.monitor.log_path else None):
            if self.monitor.snapshot()['owned']:
                self.warn('This GauMon-owned Gaussian calculation is still running. To retain safe stop controls, switching the main monitor is disabled until it finishes.')
                owner=next((x['id'] for x in self.hub.store.items if x['path']==str(self.monitor.log_path.resolve())),None)
                if owner in self.hub.jobs:
                    self.hub.choose(owner)
                return
            try:
                self.monitor.attach(path)
            except Exception as exc:
                self.warn(str(exc)); return
            self._clear_visual_caches()
        record=next((r for r in self.hub.store.items if r['path']==str(path.resolve())),None)
        self._autoload_input(path,record.get('input_path') if record else None)
        self.show_page(10)

    def _clear_visual_caches(self):
        self.last_session=None
        self.last_plot_sig=None
        self.last_energy_sig=None
        self.last_cycles_sig=None
        self.last_links_sig=None
        self.last_deep_sig=None
        self.last_qcnr_sig=None
        self.last_freq_sig=None

    def _enroll_log(self,path,source='opened',input_path=None):
        """Never enroll arbitrary processes or files unless user selected them."""
        try:
            matched=Path(input_path) if input_path else matching_input_for_log(path)
            self.hub.open_path(path,source,input_path=matched)
        except (OSError,RuntimeError) as exc:
            LOG.warning('Could not add log to session list: %s',exc)

    def _autoload_input(self,log_path,explicit_input=None):
        """Load a verified original input, never a guessed/reconstructed file."""
        result=self.input_checker.autoload_for_log(log_path,explicit_input)
        if result and self.input_checker.path is not None:
            self.input_entry.setText(str(self.input_checker.path))
        return result

    def show_page(self, n):
        if n == 7 and not self._deep_page_open:
            self.open_deep()
            return
        if n == 8 and not self._qcnr_page_open:
            self.open_qcnr()
            return
        if n == 9 and not self._freq_page_open:
            self.open_freq()
            return
        self.pages.setCurrentIndex(n); self.nav_buttons[n].setChecked(True)
        self.metric_bar.setVisible(n not in (7,8,9,10,11,12,13,14))
        self.tick()

    def open_deep(self):
        # A deliberate user action reveals the advanced page; ordinary monitoring
        # never draws the expensive deep diagnostic charts.
        self._deep_page_open=True
        self.nav_buttons[7].setVisible(True)
        self._deep_history_selection='follow'
        self.deep_choose.setCurrentIndex(0)
        self.last_deep_sig=None
        self.pages.setCurrentIndex(7)
        self.metric_bar.setVisible(False)
        self.nav_buttons[7].setChecked(True)
        self.tick()

    def open_freq(self):
        self._freq_page_open=True
        self.nav_buttons[9].setVisible(True)
        self._freq_selection='follow'
        self.last_freq_sig=None
        self.pages.setCurrentIndex(9)
        self.metric_bar.setVisible(False)
        self.nav_buttons[9].setChecked(True)
        self.tick()

    def close_freq(self):
        self._freq_page_open=False
        self.nav_buttons[9].setVisible(False)
        self.show_page(0)

    def open_qcnr(self):
        self._qcnr_page_open=True
        self.nav_buttons[8].setVisible(True)
        self._qcnr_selection='follow'
        self.last_qcnr_sig=None
        self.pages.setCurrentIndex(8)
        self.metric_bar.setVisible(False)
        self.nav_buttons[8].setChecked(True)
        self.tick()

    def close_qcnr(self):
        self._qcnr_page_open=False
        self.nav_buttons[8].setVisible(False)
        self.show_page(0)

    def close_deep(self):
        self._deep_page_open=False
        self.nav_buttons[7].setVisible(False)
        self.metric_bar.setVisible(True)
        self.show_page(0)


    def _overview(self):
        page=QWidget();lay=row_layout(11);page.setLayout(lay)
        plots=QWidget(); grid=QGridLayout(plots);grid.setSpacing(10);grid.setContentsMargins(0,0,0,0)
        self.plots={}
        specs=[('energy','SCF energy  ·  Hartree',C['green'],False,False),
               ('force','RMS force  ·  Hartree / Bohr',C['orange'],True,False),
               ('disp','Max atomic displacement  ·  Å',C['coral'],False,False),
               ('rmsd','RMSD  ·  initial / previous  ·  Å',C['green'],False,True)]
        for i,(key,title,color,log,paired) in enumerate(specs):
            widget=PlotCard(title,color,is_log=log,paired=paired)
            grid.addWidget(widget,i//2,i%2); self.plots[key]=widget
        self.disp_line_item=pg.InfiniteLine(pos=self.disp_limit,angle=0,
            pen=pg.mkPen(C['coral'],width=1.5,style=Qt.DashLine))
        self.rmsd_line_item=pg.InfiniteLine(pos=self.rmsd_limit,angle=0,
            pen=pg.mkPen(C['violet'],width=1.5,style=Qt.DashLine))
        self.plots['disp'].plot.addItem(self.disp_line_item)
        self.plots['rmsd'].plot.addItem(self.rmsd_line_item)
        grid.setRowStretch(0,1);grid.setRowStretch(1,1)
        grid.setColumnStretch(0,1);grid.setColumnStretch(1,1)
        lay.addWidget(plots,1)
        right,v=self._panel(); right.setFixedWidth(278)
        v.addWidget(label('Run details','section'))
        self.info={}
        for key,title in [('status','Status'),('pid','Process ID'),('mode','Connection'),
                          ('link','Current link'),('link_desc','Link purpose'),('scf_cycles','SCF iterations'),('charge','Charge / spin'),('route','Method / route'),
                          ('freq','Frequency blocks'),('geom','Geometry snapshots'),
                          ('disp','Max displacement'),('cycle_avg','Cycle pace'),('eta_max','ETA to max')]:
            row=row_layout(5);row.addWidget(label(title,'detailKey'))
            row.addStretch(1)
            value=label('—','detailValue');value.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
            value.setMaximumWidth(142);value.setToolTip(title)
            row.addWidget(value);v.addLayout(row);self.info[key]=value
        v.addStretch(1)
        divider=QFrame();divider.setFrameShape(QFrame.HLine); divider.setStyleSheet('background: #E4EBE4; max-height: 1px;')
        v.addWidget(divider)
        self.quick_state=label('Select a calculation to start monitoring.', 'small',wrap=True)
        v.addWidget(self.quick_state)
        eta_tile=QFrame();eta_tile.setObjectName('insight')
        eta_v=col_layout(3);eta_v.setContentsMargins(12,9,12,9);eta_tile.setLayout(eta_v)
        eta_v.addWidget(label('ESTIMATED NEXT 10 CYCLES','eyebrow'))
        self.eta_peek=label('Awaiting cycle data','calloutValue')
        eta_v.addWidget(self.eta_peek)
        self.eta_peek_time=label('Requires at least two timed steps','small',wrap=True)
        eta_v.addWidget(self.eta_peek_time)
        v.addWidget(eta_tile)
        link_title=label('Current program stage','section'); v.addWidget(link_title)
        self.link_peek=label('No link records yet','small',wrap=True)
        self.link_peek.setMaximumHeight(47);v.addWidget(self.link_peek)
        timeline_btn=QPushButton('Open timeline & ETA  →');timeline_btn.setObjectName('subtle')
        timeline_btn.clicked.connect(lambda:self.show_page(2));v.addWidget(timeline_btn)
        deep_btn=QPushButton('Deep SCF diagnostics  →');deep_btn.setObjectName('subtle')
        deep_btn.clicked.connect(self.open_deep);v.addWidget(deep_btn)
        qcnr_btn=QPushButton('Deep QCNR / LinEq1  →');qcnr_btn.setObjectName('subtle')
        qcnr_btn.clicked.connect(self.open_qcnr);v.addWidget(qcnr_btn)
        freq_btn=QPushButton('Deep Frequencies  →');freq_btn.setObjectName('subtle')
        freq_btn.clicked.connect(self.open_freq);v.addWidget(freq_btn)
        nav=row_layout(5)
        run=QPushButton('New job');run.setObjectName('primary');run.clicked.connect(lambda:self.show_page(1))
        log=QPushButton('Log');log.clicked.connect(lambda:self.show_page(5))
        nav.addWidget(run,1);nav.addWidget(log,1)
        v.addLayout(nav)
        lay.addWidget(right)
        self.pages.addWidget(page)

    def _field(self,label_text,default='',browse=False):
        row=row_layout(7)
        entry=QLineEdit(default);entry.setMinimumHeight(29)
        row.addWidget(entry,1)
        if browse:
            b=QPushButton('Browse');b.clicked.connect(lambda:self.browse(entry,browse))
            row.addWidget(b)
        return label(label_text,'detailKey'),row,entry

    def browse(self, edit, mode):
        current=edit.text().strip() or os.getcwd()
        if mode=='input':
            path,_=QFileDialog.getOpenFileName(self,'Gaussian input',current,'Gaussian Input (*.gjf *.com);;All files (*)')
        elif mode=='log':
            path,_=QFileDialog.getOpenFileName(self,'Gaussian output',current,'Gaussian Logs (*.log *.out *.txt);;All files (*)')
        else:
            path=QFileDialog.getExistingDirectory(self,'Working directory',current)
        if path:
            edit.setText(path)
            if mode=='input':
                p=Path(path);self.work_entry.setText(str(p.parent))
                self.log_entry.setText(p.with_suffix('.log').name)
                if not self.input_checker.load(p,automatic=True):
                    # Do not point the Start button at a different file from
                    # the source the user elected to keep in the editor.
                    self.input_entry.setText(str(self.input_checker.path) if self.input_checker.path else '')
                    if self.input_checker.path:
                        self.work_entry.setText(str(self.input_checker.path.parent))
                        self.log_entry.setText(self.input_checker.path.with_suffix('.log').name)

    def _run_page(self):
        page=QWidget();cols=row_layout(12);page.setLayout(cols)
        left,v=self._panel();v.addWidget(label('Run a Gaussian calculation','section'))
        v.addWidget(label('A licensed Gaussian executable is required only for actual calculations.','small',wrap=True))
        form=QFormLayout();form.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)
        form.setVerticalSpacing(13);form.setHorizontalSpacing(16)
        specs=[('Gaussian executable','g16',None,'exec_entry'),
               ('Input .gjf / .com','', 'input','input_entry'),
               ('Working directory',os.getcwd(), 'directory','work_entry'),
               ('Output log filename','',None,'log_entry')]
        for name,default,browse,attr in specs:
            text=row_layout();edit=QLineEdit(default);edit.setMinimumHeight(31)
            text.addWidget(edit,1)
            if browse:
                btn=QPushButton('Browse…');btn.clicked.connect(lambda _=False,e=edit,m=browse:self.browse(e,m))
                text.addWidget(btn)
            holder=QWidget();holder.setLayout(text)
            form.addRow(name,holder);setattr(self,attr,edit)
        v.addLayout(form)
        actions=row_layout(8)
        self.start_btn=QPushButton('Start Gaussian');self.start_btn.setObjectName('primary')
        self.start_btn.clicked.connect(self.launch)
        self.stop_btn=QPushButton('Stop gracefully');self.stop_btn.clicked.connect(lambda:self.stop(False))
        self.kill_btn=QPushButton('Force stop');self.kill_btn.setObjectName('danger')
        self.kill_btn.clicked.connect(lambda:self.stop(True))
        actions.addWidget(self.start_btn,2);actions.addWidget(self.stop_btn,1);actions.addWidget(self.kill_btn,1)
        v.addLayout(actions)
        v.addStretch(1)
        cols.addWidget(left,3)
        right,r=self._panel();r.addWidget(label('Attach an existing log','section'))
        r.addWidget(label('An attached external calculation is read-only and cannot be stopped by GauMon.','small',wrap=True))
        field=row_layout(4);self.attach_entry=QLineEdit();self.attach_entry.setPlaceholderText('/path/to/job.log')
        field.addWidget(self.attach_entry,1)
        btn=QPushButton('Browse');btn.clicked.connect(lambda:self.browse(self.attach_entry,'log'));field.addWidget(btn)
        r.addLayout(field)
        r.addWidget(label('Optional process PID','detailKey'))
        self.pid_entry=QLineEdit();self.pid_entry.setPlaceholderText('Leave empty for log-only observation')
        r.addWidget(self.pid_entry)
        list_jobs=QPushButton('Find running Gaussian processes')
        list_jobs.setObjectName('subtle');list_jobs.clicked.connect(self.pick_running_job)
        r.addWidget(list_jobs)
        attach=QPushButton('Attach and inspect');attach.setObjectName('primary');attach.clicked.connect(self.attach)
        r.addWidget(attach)
        self.input_hint=label('Input validation will identify common missing title-card errors before launching.','small',wrap=True)
        r.addStretch(1);r.addWidget(self.input_hint)
        cols.addWidget(right,2)
        self.pages.addWidget(page)

    def _timeline_page(self):
        page=QWidget();outer=col_layout(12);page.setLayout(outer)
        heading=row_layout(8)
        heading.addWidget(label('Optimization timeline & forecasts','section'))
        heading.addStretch(1)
        self.timeline_quality=label('Timing confidence: waiting for data','small')
        heading.addWidget(self.timeline_quality)
        outer.addLayout(heading)
        projections=row_layout(10)
        self.eta_tiles={}
        for n in (10,50,100,150):
            tile=QFrame();tile.setObjectName('insight')
            v=col_layout(3);v.setContentsMargins(15,12,15,12);tile.setLayout(v)
            v.addWidget(label(f'NEXT +{n} CYCLES','eyebrow'))
            duration=label('—','calloutValue');v.addWidget(duration)
            due=label('No reliable estimate yet','small',wrap=True);v.addWidget(due)
            projections.addWidget(tile,1)
            self.eta_tiles[n]=(duration,due)
        outer.addLayout(projections)
        body=row_layout(12)
        left,l=self._panel()
        h=row_layout();h.addWidget(label('Cycle history','section'));h.addStretch(1)
        self.cycle_count_label=label('0 steps','small');h.addWidget(self.cycle_count_label)
        l.addLayout(h)
        self.cycle_table=QTableWidget(0,5)
        self.cycle_table.setHorizontalHeaderLabels(['Step','SCF cycles','Cycle time','Source','Energy (Ha)'])
        self._style_table(self.cycle_table)
        self.cycle_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        l.addWidget(self.cycle_table,1)
        body.addWidget(left,5)
        right,r=self._panel()
        r.addWidget(label('Gaussian link history','section'))
        self.link_table=QTableWidget(0,3)
        self.link_table.setHorizontalHeaderLabels(['Link','Log time','Function'])
        self._style_table(self.link_table)
        self.link_table.horizontalHeader().setSectionResizeMode(0,QHeaderView.ResizeToContents)
        self.link_table.horizontalHeader().setSectionResizeMode(1,QHeaderView.ResizeToContents)
        self.link_table.horizontalHeader().setSectionResizeMode(2,QHeaderView.Stretch)
        r.addWidget(self.link_table,1)
        link_current=row_layout(8)
        self.current_link_label=label('Current: —','detailValue')
        link_current.addWidget(self.current_link_label);link_current.addStretch(1)
        r.addLayout(link_current)
        self.link_explanation=label('Link descriptions from gaussian_links.json','small',wrap=True)
        r.addWidget(self.link_explanation)
        body.addWidget(right,4)
        outer.addLayout(body,1)
        note=label('Forecasts are extrapolations, not guaranteed Gaussian completion times. '
                   'Only measured cycle intervals are used. Earlier events in attached logs '
                   'may lack timestamps; finished jobs have no remaining ETA.', 'small',wrap=True)
        outer.addWidget(note)
        self.pages.addWidget(page)

    def _analysis_page(self):
        page=QWidget();lay=row_layout(10);page.setLayout(lay)
        left,v=self._panel()
        topbar=row_layout(8);topbar.addWidget(label('SCF history','section'));topbar.addStretch(1)
        deep=QPushButton('Deep SCF  →');deep.setObjectName('subtle')
        deep.clicked.connect(self.open_deep);topbar.addWidget(deep)
        v.addLayout(topbar)
        self.energy_table=QTableWidget(0,3); self.energy_table.setHorizontalHeaderLabels(['Index','Energy  (Ha)','SCF cycles'])
        self._style_table(self.energy_table)
        v.addWidget(self.energy_table,1); lay.addWidget(left,3)
        right,r=self._panel();r.addWidget(label('Calculation summary','section'))
        self.summary=QPlainTextEdit();self.summary.setReadOnly(True);r.addWidget(self.summary,1)
        r.addWidget(label('Geometry alerts  ·  observation only', 'section'))
        thresholds=QFormLayout(); thresholds.setVerticalSpacing(6)
        self.disp_spin=QDoubleSpinBox();self.disp_spin.setRange(.01,99999.0);self.disp_spin.setDecimals(2)
        self.disp_spin.setSuffix(' Å');self.disp_spin.setValue(self.disp_limit)
        self.rmsd_spin=QDoubleSpinBox();self.rmsd_spin.setRange(.01,99999.0);self.rmsd_spin.setDecimals(2)
        self.rmsd_spin.setSuffix(' Å');self.rmsd_spin.setValue(self.rmsd_limit)
        self.disp_spin.valueChanged.connect(self._change_disp_limit)
        self.rmsd_spin.valueChanged.connect(self._change_rmsd_limit)
        thresholds.addRow('Max displacement',self.disp_spin)
        thresholds.addRow('RMSD vs initial',self.rmsd_spin)
        r.addLayout(thresholds)
        lay.addWidget(right,2);self.pages.addWidget(page)

    def _change_disp_limit(self,value):
        self.disp_limit=float(value)
        self.disp_line_item.setPos(self.disp_limit)

    def _change_rmsd_limit(self,value):
        self.rmsd_limit=float(value)
        self.rmsd_line_item.setPos(self.rmsd_limit)

    def _geometry_page(self):
        page=QWidget();v=col_layout(8);page.setLayout(v)
        panel,inner=self._panel()
        head=row_layout();head.addWidget(label('Geometry inspection','section'))
        head.addStretch(1);head.addWidget(label('XY projection · not a 3D molecular model','small'))
        inner.addLayout(head)
        self.molecule=MoleculeCanvas();inner.addWidget(self.molecule,1)
        self.geom_note=label('No geometry available. Coordinate data is parsed from Gaussian orientation tables.','small',wrap=True)
        inner.addWidget(self.geom_note)
        v.addWidget(panel,1);self.pages.addWidget(page)

    def _log_page(self):
        page=QWidget();v=col_layout(8);page.setLayout(v)
        panel,inner=self._panel()
        h=row_layout();h.addWidget(label('Live Gaussian output','section'))
        h.addStretch(1);self.follow=QCheckBox('Follow latest');self.follow.setChecked(True);h.addWidget(self.follow)
        inner.addLayout(h)
        self.log_box=QPlainTextEdit();self.log_box.setReadOnly(True);self.log_box.setLineWrapMode(QPlainTextEdit.NoWrap)
        inner.addWidget(self.log_box,1)
        diag_head=row_layout()
        diag_head.addWidget(label('Diagnostics','section'))
        diag_head.addStretch(1)
        self.source_btn=QPushButton('Error guide · Zhe Wang ↗')
        self.source_btn.setObjectName('subtle')
        self.source_btn.clicked.connect(lambda: QDesktopServices.openUrl(QUrl(GUIDE_URL)))
        diag_head.addWidget(self.source_btn)
        inner.addLayout(diag_head)
        self.diagnostic_label=QPlainTextEdit()
        self.diagnostic_label.setReadOnly(True)
        self.diagnostic_label.setLineWrapMode(QPlainTextEdit.WidgetWidth)
        self.diagnostic_label.setFixedHeight(138)
        self.diagnostic_label.setPlainText('No detected errors · monitoring output')
        inner.addWidget(self.diagnostic_label)
        v.addWidget(panel,1);self.pages.addWidget(page)

    def _tests_page(self):
        page=QWidget();v=col_layout(8);page.setLayout(v)
        panel,inner=self._panel()
        inner.addWidget(label('Gaussian test laboratory','section'))
        inner.addWidget(label('Four basic inputs and three advanced stress inputs. Synthetic runs do not perform quantum chemistry.','small',wrap=True))
        self.test_table=QTableWidget(len(ALL_CASES),4)
        self.test_table.setHorizontalHeaderLabels(['Type','System','What it checks','Expected evidence'])
        self._style_table(self.test_table)
        for i,case in enumerate(ALL_CASES):
            for j,text in enumerate((case.label,case.formula,case.description,case.expected)):
                item=QTableWidgetItem(text);item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                item.setToolTip(text);self.test_table.setItem(i,j,item)
        self.test_table.setCurrentCell(0,0)
        self.test_table.horizontalHeader().setSectionResizeMode(0,QHeaderView.ResizeToContents)
        self.test_table.horizontalHeader().setSectionResizeMode(1,QHeaderView.ResizeToContents)
        self.test_table.horizontalHeader().setSectionResizeMode(2,QHeaderView.Stretch)
        self.test_table.horizontalHeader().setSectionResizeMode(3,QHeaderView.Stretch)
        self.test_table.verticalHeader().setDefaultSectionSize(52)
        inner.addWidget(self.test_table,1)
        actions=row_layout(9)
        mock=QPushButton('Run synthetic test');mock.setObjectName('primary');mock.clicked.connect(self.simulate_selected)
        export=QPushButton('Export selected .gjf');export.clicked.connect(self.export_selected)
        export_all=QPushButton('Export all seven');export_all.clicked.connect(self.export_all)
        actions.addWidget(mock);actions.addWidget(export);actions.addWidget(export_all);actions.addStretch(1)
        inner.addLayout(actions)
        stressbar=row_layout(8)
        stressbar.addWidget(label('SCF stress fixture:', 'detailKey'))
        self.stress_choose=QComboBox()
        for case in ('many_steps','oscillation','incomplete','chunked'):
            self.stress_choose.addItem(case.replace('_',' ').title(),case)
        stressbar.addWidget(self.stress_choose)
        open_stress=QPushButton('Inspect in Deep SCF →')
        open_stress.clicked.connect(self.open_stress_fixture)
        stressbar.addWidget(open_stress)
        stressbar.addStretch(1)
        inner.addLayout(stressbar)
        self.test_note=label('Select any row; use the synthetic run to verify GUI and parser without Gaussian.','small',wrap=True)
        inner.addWidget(self.test_note)
        v.addWidget(panel,1);self.pages.addWidget(page)

    def _deep_page(self):
        """On-demand scientific SCF plots with history and reset-per-SCF semantics."""
        page=QWidget();outer=col_layout(9);page.setLayout(outer)
        header=row_layout(10)
        header.addWidget(label('Deep SCF convergence','section'))
        header.addStretch(1)
        header.addWidget(label('SCF procedure:', 'detailKey'))
        self.deep_choose=QComboBox()
        self.deep_choose.setMinimumWidth(255)
        self.deep_choose.addItem('Follow newest · reset at next SCF', 'follow')
        self.deep_choose.currentIndexChanged.connect(self._deep_changed)
        header.addWidget(self.deep_choose)
        close=QPushButton('Close deep view');close.setObjectName('subtle')
        close.clicked.connect(self.close_deep)
        header.addWidget(close)
        outer.addLayout(header)
        summary=row_layout(10)
        self.deep_kpis={}
        for key, title in [('step','OPTIMIZATION'),('iter','SCF CYCLE'),
                           ('diis','DIIS ERROR'),('delta','ΔE · HARTREE')]:
            metric=Metric(title)
            metric.value.setText('—')
            summary.addWidget(metric,1);self.deep_kpis[key]=metric
        outer.addLayout(summary)
        chart_area=QWidget();grid=QGridLayout(chart_area);grid.setSpacing(9);grid.setContentsMargins(0,0,0,0)
        self.deep_charts={}
        entries=[('energy','Electronic energy · Hartree',C['green'],False),
                 ('delta','Absolute energy change · Hartree',C['orange'],True),
                 ('diis','DIIS error / RMSU · diagnostic scale',C['violet'],True),
                 ('density','RMSDP / MaxDP · density change',C['coral'],True)]
        for idx,(key, title, color, log_mode) in enumerate(entries):
            card=PlotCard(title,color,is_log=log_mode,paired=key in ('diis','density'))
            card.plot.setMinimumSize(160,105)
            grid.addWidget(card,idx//2,idx%2)
            self.deep_charts[key]=card
        grid.setRowStretch(0,1);grid.setRowStretch(1,1)
        grid.setColumnStretch(0,1);grid.setColumnStretch(1,1)
        outer.addWidget(chart_area,4)
        details=row_layout(10)
        table_panel,v=self._panel()
        top=row_layout(6);top.addWidget(label('SCF iteration records','section'));top.addStretch(1)
        save=QPushButton('Export CSV');save.setObjectName('subtle');save.clicked.connect(self.export_deep_csv)
        top.addWidget(save);v.addLayout(top)
        self.deep_table=QTableWidget(0,8)
        self.deep_table.setHorizontalHeaderLabels(['Cycle','Pass','Energy (Ha)','ΔE (Ha)','DIIS error','RMSU','RMSDP','MaxDP'])
        self._style_table(self.deep_table)
        self.deep_table.setMinimumHeight(110)
        self.deep_table.verticalHeader().setDefaultSectionSize(27)
        self.deep_table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.deep_table.horizontalHeader().setMinimumSectionSize(68)
        self.deep_table.itemSelectionChanged.connect(self._deep_row_selected)
        v.addWidget(self.deep_table,1);details.addWidget(table_panel,3)
        readout,rv=self._panel()
        rv.addWidget(label('Selected iteration','section'))
        self.deep_detail=QPlainTextEdit();self.deep_detail.setReadOnly(True)
        self.deep_detail.setMinimumHeight(105)
        rv.addWidget(self.deep_detail,1)
        details.addWidget(readout,2)
        outer.addLayout(details,2)
        self.deep_note=label('The view follows one SCF procedure at a time. When Gaussian begins a new SCF, the live plots reset. Older procedures remain selectable. Log-scale plots omit zero or missing points.','small',wrap=True)
        outer.addWidget(self.deep_note)
        self.pages.addWidget(page)

    def _deep_changed(self):
        if not hasattr(self,'deep_choose'):
            return
        self._deep_history_selection=self.deep_choose.currentData()
        self.last_deep_sig=None
        if self._deep_page_open:
            self._update_deep(self.monitor.snapshot())

    @staticmethod
    def _deep_fmt(value, digits=4):
        if value is None:
            return '—'
        if isinstance(value, (int,float)):
            return f'{value:.{digits}g}'
        return str(value)

    def _deep_row_selected(self):
        if not hasattr(self,'deep_table'):
            return
        selection=self.deep_table.currentRow()
        block=self._deep_block_for_export
        if not block or selection < 0 or selection >= len(block['rows']):
            return
        rec=block['rows'][selection]
        lines=[f'Cycle {rec["cycle"]} · Pass {rec["pass"]} · IDiag {rec["idiag"]}']
        for key in ('energy','delta_e','diis_error','rmsu','rmsdp','maxdp','ovmax',
                    'nsaved','wtcom','wten','bmatc','bmatp','errmax','errmin',
                    'enmin','idiuse','rises','damp'):
            if rec.get(key) is not None:
                lines.append(f'{key:12} {self._deep_fmt(rec[key],10)}')
        for key in ('gaps','cp','coeff','coeff_com','coeff_en'):
            if rec.get(key):
                values=', '.join(self._deep_fmt(v,6) for v in rec[key][:20])
                if len(rec[key])>20:
                    values+=' …'
                lines.append(f'{key}: {values}')
        self.deep_detail.setPlainText('\n'.join(lines))

    def _update_deep(self, state):
        if not self._deep_page_open or self.pages.currentIndex() != 7:
            return
        meta=state.get('deep_scf') or {}
        blocks=meta.get('blocks') or []
        selected=self._deep_history_selection
        stamp=(state.get('session'), meta.get('revision'), selected)
        if stamp == self.last_deep_sig:
            return
        self.last_deep_sig=stamp
        # Selectors for archive; users can inspect prior SCF procedures without
        # changing which calculation is monitored.
        old=self.deep_choose.blockSignals(True)
        self.deep_choose.clear()
        self.deep_choose.addItem('Follow newest · reset when next SCF begins','follow')
        for item in reversed(blocks):
            name=item['label']+' · '+item['status']+' · '+str(item['cycles'])+' cycles'
            self.deep_choose.addItem(name,item['id'])
        index=self.deep_choose.findData(selected)
        if index<0:
            index=0
            self._deep_history_selection='follow'
        self.deep_choose.setCurrentIndex(index)
        self.deep_choose.blockSignals(old)
        target=None if self._deep_history_selection=='follow' else self._deep_history_selection
        block=self.monitor.get_deep_block(target)
        self._deep_block_for_export=block
        for card in self.deep_charts.values():
            card.refresh([])
        self.deep_table.setRowCount(0)
        if not block:
            self.deep_note.setText('No Gaussian SCF iteration records found yet. Attach a log with Cycle/Pass entries or run a calculation.')
            for m in self.deep_kpis.values(): m.value.setText('—')
            self.deep_detail.clear()
            return
        rows=block['rows']
        self.deep_kpis['step'].value.setText(str(block.get('optimization_step') or '—'))
        self.deep_kpis['iter'].value.setText(str(rows[-1]['cycle']) if rows else '—')
        last=rows[-1] if rows else {}
        self.deep_kpis['diis'].value.setText(self._deep_fmt(last.get('diis_error')))
        self.deep_kpis['delta'].value.setText(self._deep_fmt(last.get('delta_e'),7))
        def points(key, absolute=False):
            return [(r['cycle'],abs(r[key]) if absolute else r[key]) for r in rows
                    if r.get(key) is not None and math.isfinite(r[key]) and
                    (not absolute or r[key] != 0)]
        # The three magnitude/error plots use log axes: don't draw zeros.
        self.deep_charts['energy'].refresh(points('energy'))
        self.deep_charts['delta'].refresh(points('delta_e',True))
        self.deep_charts['diis'].refresh(points('diis_error',True),points('rmsu',True))
        self.deep_charts['density'].refresh(points('rmsdp',True),points('maxdp',True))
        self.deep_table.setRowCount(len(rows))
        fields=('cycle','pass','energy','delta_e','diis_error','rmsu','rmsdp','maxdp')
        for i,rec in enumerate(rows):
            for j,key in enumerate(fields):
                self.deep_table.setItem(i,j,QTableWidgetItem(self._deep_fmt(rec.get(key),10)))
        if rows:
            self.deep_table.selectRow(len(rows)-1)
        self.deep_note.setText(f'{block["label"]} · {block["status"].upper()} · {len(rows)} recorded cycles'
                               + (' · HISTORY TRUNCATED' if block.get('truncated') else '')
                               + ' · SCF iteration ≠ geometry-optimization step')

    def export_deep_csv(self):
        block=self._deep_block_for_export
        if not block or not block['rows']:
            self.warn('No deep SCF data to export.');return
        file,_=QFileDialog.getSaveFileName(self,'Export SCF iterations',
                    str(Path.cwd()/f'GauMon_SCF_{block["id"]}.csv'),'CSV files (*.csv)')
        if not file:return
        try:
            fields=('cycle','pass','idiag','energy','delta_e','diis_error','rmsu',
                    'rmsdp','maxdp','ovmax','gap','nsaved','wtcom','wten','rises','damp',
                    'cp','coeff','coeff_com','coeff_en','gaps')
            with open(file,'w',newline='',encoding='utf-8') as out:
                writer=csv.DictWriter(out,fieldnames=fields,extrasaction='ignore')
                writer.writeheader()
                for rec in block['rows']:
                    writer.writerow(rec)
        except OSError as exc:
            self.warn(f'CSV export failed: {exc}')

    def _qcnr_page(self):
        """Opt-in QCNR orbital / nested LinEq1 solver convergence workspace."""
        page=QWidget();v=col_layout(9);page.setLayout(v)
        head=row_layout(8)
        head.addWidget(label('Deep QCNR  /  LinEq1','section'))
        head.addStretch(1)
        self.qcnr_choose=QComboBox();self.qcnr_choose.setMinimumWidth(265)
        self.qcnr_choose.addItem('Follow newest linear solve','follow')
        self.qcnr_choose.currentIndexChanged.connect(self._qcnr_changed)
        head.addWidget(self.qcnr_choose)
        close=QPushButton('Close deep view');close.setObjectName('subtle')
        close.clicked.connect(self.close_qcnr);head.addWidget(close)
        v.addLayout(head)
        kpis=row_layout(9);self.qcnr_kpis={}
        for key,title in [('outer','QCNR ITERATION'),('inner','LINEAR ITERATION'),('rms','LINEAR RMS'),('status','SOLVER STATUS')]:
            card=Metric(title);card.value.setText('—');kpis.addWidget(card,1);self.qcnr_kpis[key]=card
        v.addLayout(kpis)
        canvas=QWidget();grid=QGridLayout(canvas);grid.setSpacing(9);grid.setContentsMargins(0,0,0,0)
        self.qcnr_plots={}
        cards=[('rms','LinEq1 · RMS residual',C['green'],True),
               ('max','LinEq1 · Max residual',C['coral'],True),
               ('energy','QCNR outer · EE / Hartree',C['violet'],False),
               ('grad','QCNR outer · Gradient',C['orange'],False)]
        for i,(key,caption,color,log_mode) in enumerate(cards):
            chart=PlotCard(caption,color,is_log=log_mode,paired=False)
            chart.plot.setMinimumSize(160,103)
            grid.addWidget(chart,i//2,i%2);self.qcnr_plots[key]=chart
        grid.setColumnStretch(0,1);grid.setColumnStretch(1,1)
        grid.setRowStretch(0,1);grid.setRowStretch(1,1)
        v.addWidget(canvas,4)
        bottom=row_layout(9)
        panel, pv=self._panel();h=row_layout(8)
        h.addWidget(label('LinEq1 iteration history','section'))
        h.addStretch(1)
        btn=QPushButton('Export CSV');btn.setObjectName('subtle')
        btn.clicked.connect(self.export_qcnr_csv);h.addWidget(btn)
        pv.addLayout(h)
        self.qcnr_table=QTableWidget(0,5)
        self.qcnr_table.setHorizontalHeaderLabels(['Iter','NonCon','RMS','Max','NDo'])
        self._style_table(self.qcnr_table)
        self.qcnr_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.qcnr_table.verticalHeader().setDefaultSectionSize(27)
        pv.addWidget(self.qcnr_table,1);bottom.addWidget(panel,3)
        detail, dv=self._panel();dv.addWidget(label('Solver interpretation','section'))
        self.qcnr_details=QPlainTextEdit();self.qcnr_details.setReadOnly(True)
        dv.addWidget(self.qcnr_details,1)
        bottom.addWidget(detail,2)
        v.addLayout(bottom,2)
        self.qcnr_note=label('LinEq1 is the nested linear-equation solver in quadratic-convergence SCF (L508); it is NOT a frequency calculation. The live residual plots reset for each new linear solve.','small',wrap=True)
        v.addWidget(self.qcnr_note)
        self.pages.addWidget(page)

    def _qcnr_changed(self):
        if not hasattr(self,'qcnr_choose'):return
        self._qcnr_selection=self.qcnr_choose.currentData()
        self.last_qcnr_sig=None
        if self._qcnr_page_open:self._update_qcnr(self.monitor.snapshot())

    def _update_qcnr(self,state):
        if not self._qcnr_page_open or self.pages.currentIndex()!=8:return
        meta=state.get('deep_qcnr') or {}
        selected=self._qcnr_selection
        signature=(state.get('session'),meta.get('revision'),selected)
        if signature==self.last_qcnr_sig:return
        self.last_qcnr_sig=signature
        old=self.qcnr_choose.blockSignals(True)
        self.qcnr_choose.clear()
        self.qcnr_choose.addItem('Follow newest · reset on next LinEq1','follow')
        for entry in reversed(meta.get('solves') or []):
            self.qcnr_choose.addItem(f"{entry['label']} · {entry['status']} · {entry['iterations']} points",entry['id'])
        target_index=self.qcnr_choose.findData(selected)
        if target_index<0:target_index=0;self._qcnr_selection='follow'
        self.qcnr_choose.setCurrentIndex(target_index)
        self.qcnr_choose.blockSignals(old)
        which=None if self._qcnr_selection=='follow' else self._qcnr_selection
        solve=self.monitor.get_qcnr_solve(which)
        self._qcnr_export=solve
        outer=self.monitor.get_qcnr_outer_history()
        # Filter outer iterations to the same QCNR epoch as the chosen solve.
        outer=[p for p in outer if solve is not None and p['epoch']==solve['epoch']]
        for chart in self.qcnr_plots.values():chart.refresh([])
        self.qcnr_table.setRowCount(0)
        if not solve:
            for card in self.qcnr_kpis.values():card.value.setText('—')
            self.qcnr_details.setPlainText('No LinEq1 records in this log. A frequency or ordinary DIIS calculation may not execute quadratic SCF.')
            return
        rows=solve['rows'];last=rows[-1] if rows else {}
        self.qcnr_kpis['outer'].value.setText(str(solve.get('outer_iteration') or '—'))
        self.qcnr_kpis['inner'].value.setText(str(last.get('iteration','—')))
        self.qcnr_kpis['rms'].value.setText(self._deep_fmt(last.get('rms'),5))
        self.qcnr_kpis['status'].value.setText(solve['status'].upper())
        def coords(vals,key,positive=False):
            return [(x,y) for rec in vals if (x:=rec.get('iteration')) is not None and
                    (y:=rec.get(key)) is not None and math.isfinite(y) and (not positive or y>0)]
        self.qcnr_plots['rms'].refresh(coords(rows,'rms',True))
        self.qcnr_plots['max'].refresh(coords(rows,'max',True))
        self.qcnr_plots['energy'].refresh(coords(outer,'energy'))
        # Plot gradient only: max rotation is measured in degrees, not gradient units.
        self.qcnr_plots['grad'].refresh(coords(outer,'grad'))
        self.qcnr_table.setRowCount(len(rows))
        for i,rec in enumerate(rows):
            for j,key in enumerate(('iteration','noncon','rms','max','ndo')):
                self.qcnr_table.setItem(i,j,QTableWidgetItem(self._deep_fmt(rec.get(key),7)))
        if rows:self.qcnr_table.selectRow(len(rows)-1)
        convergence=solve.get('actual_convergence') or {}
        self.qcnr_details.setPlainText('\n'.join([
            f"Link: {meta.get('current_link') or '—'} · Quadratic SCF / Newton-Raphson",
            f"QCNR outer iteration: {solve.get('outer_iteration') or '—'}",
            f"Inner solver status: {solve['status']}",
            f"LinEq1 records: {len(rows)}",
            f"Latest RMS residual: {self._deep_fmt(last.get('rms'),7)}",
            f"Latest maximum residual: {self._deep_fmt(last.get('max'),7)}",
            f"NonCon (unconverged RHS): {last.get('noncon','—')}",
            f"NDo (systems processed): {last.get('ndo','—')}",
            f"Requested RMS / max thresholds: {self._deep_fmt(solve.get('target_rms'))} / {self._deep_fmt(solve.get('target_max'))}",
            f"Solver-reported iterations: {convergence.get('iterations_reported','—')}",
            f"QCNR orbital max rotation: {self._deep_fmt(next((r.get('max_rot') for r in reversed(outer) if r.get('iteration')==solve.get('outer_iteration')), None))} degrees",
            'A completed inner solve does not prove overall SCF convergence.',
        ]))
        self.qcnr_note.setText(f"{solve['label']}  ·  {solve['status'].upper()}  ·  {len(rows)} linear iterations"
                               + (' · HISTORY TRUNCATED' if solve.get('truncated') else '')
                               + ' · Selecting another solve does not change the monitored job.')

    def export_qcnr_csv(self):
        solve=self._qcnr_export
        if not solve or not solve['rows']:
            self.warn('No LinEq1 iterations to export.');return
        filename,_=QFileDialog.getSaveFileName(self,'Export LinEq1 data',
            str(Path.cwd()/f'GauMon_LinEq1_{solve["id"]}.csv'),'CSV files (*.csv)')
        if not filename:return
        try:
            with open(filename,'w',newline='',encoding='utf-8') as output:
                writer=csv.DictWriter(output,fieldnames=('iteration','noncon','rms','max','ndo'),extrasaction='ignore')
                writer.writeheader();writer.writerows(solve['rows'])
        except OSError as error:self.warn(f'CSV export failed: {error}')

    def _freq_page(self):
        """Opt-in harmonic mode properties (not LinEq1)."""
        page=QWidget();v=col_layout(9);page.setLayout(v)
        head=row_layout(8)
        head.addWidget(label('Deep Frequencies','section'))
        head.addStretch(1)
        self.freq_choose=QComboBox();self.freq_choose.setMinimumWidth(250)
        self.freq_choose.addItem('Latest printed spectrum','follow')
        self.freq_choose.currentIndexChanged.connect(self._freq_changed)
        head.addWidget(self.freq_choose)
        close=QPushButton('Close deep view');close.setObjectName('subtle')
        close.clicked.connect(self.close_freq);head.addWidget(close)
        v.addLayout(head)
        kpis=row_layout(9);self.freq_kpis={}
        for key,title in [('total','MODES PRINTED'),('imaginary','IMAGINARY MODES'),('lowest','LOWEST · CM⁻¹'),('highest','HIGHEST · CM⁻¹')]:
            card=Metric(title);card.value.setText('—');kpis.addWidget(card,1);self.freq_kpis[key]=card
        v.addLayout(kpis)
        area=QWidget();charts=QGridLayout(area);charts.setContentsMargins(0,0,0,0);charts.setSpacing(9)
        self.freq_charts={}
        for i,(key,title,color) in enumerate((('modes','Harmonic frequency by mode · cm⁻¹',C['green']),
                                              ('ir','IR intensity vs frequency · printed units',C['coral']))):
            plot=PlotCard(title,color)
            charts.addWidget(plot,0,i);self.freq_charts[key]=plot
            charts.setColumnStretch(i,1)
        v.addWidget(area,3)
        b=row_layout(9)
        panel,pv=self._panel()
        h=row_layout(9);h.addWidget(label('Printed normal modes','section'))
        h.addStretch(1);btn=QPushButton('Export CSV');btn.setObjectName('subtle')
        btn.clicked.connect(self.export_freq_csv);h.addWidget(btn);pv.addLayout(h)
        self.freq_table=QTableWidget(0,5)
        self.freq_table.setHorizontalHeaderLabels(['Mode','Freq cm⁻¹','Red. mass','Force const','IR inten'])
        self._style_table(self.freq_table)
        self.freq_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.freq_table.verticalHeader().setDefaultSectionSize(27)
        pv.addWidget(self.freq_table,1);b.addWidget(panel,3)
        right,rv=self._panel();rv.addWidget(label('Scientific context','section'))
        self.freq_details=QPlainTextEdit();self.freq_details.setReadOnly(True)
        rv.addWidget(self.freq_details,1);b.addWidget(right,2)
        v.addLayout(b,2)
        self.freq_note=label('Only actual Gaussian Frequencies -- records are shown. Negative frequencies indicate imaginary modes; inspect normal-mode displacement before claiming a transition state.','small',wrap=True)
        v.addWidget(self.freq_note)
        self.pages.addWidget(page)

    def _freq_changed(self):
        if not hasattr(self,'freq_choose'):return
        self._freq_selection=self.freq_choose.currentData()
        self.last_freq_sig=None
        if self._freq_page_open:self._update_freq(self.monitor.snapshot())

    def _update_freq(self,state):
        if not self._freq_page_open or self.pages.currentIndex()!=9:return
        meta=state.get('deep_freq') or {}
        selected=self._freq_selection
        sig=(state.get('session'),meta.get('revision'),selected)
        if sig==self.last_freq_sig:return
        self.last_freq_sig=sig
        old=self.freq_choose.blockSignals(True)
        self.freq_choose.clear()
        self.freq_choose.addItem('Latest printed spectrum','follow')
        for ent in reversed(meta.get('groups') or []):
            self.freq_choose.addItem(f"{ent['label']} · {ent['modes']} modes",ent['id'])
        ix=self.freq_choose.findData(selected)
        if ix<0:ix=0;self._freq_selection='follow'
        self.freq_choose.setCurrentIndex(ix)
        self.freq_choose.blockSignals(old)
        which=None if self._freq_selection=='follow' else self._freq_selection
        group=self.monitor.get_freq_group(which)
        self._freq_export=group
        self.freq_table.setRowCount(0)
        for chart in self.freq_charts.values():chart.refresh([])
        if not group:
            for m in self.freq_kpis.values():m.value.setText('—')
            self.freq_details.setPlainText('No vibrational frequencies printed in this log. A single-point SCF / LinEq1 record does not contain harmonic modes.')
            return
        rows=group['modes'];freqs=[r['frequency'] for r in rows]
        self.freq_kpis['total'].value.setText(str(len(rows)))
        self.freq_kpis['imaginary'].value.setText(str(sum(f<-1e-7 for f in freqs)))
        self.freq_kpis['lowest'].value.setText(self._deep_fmt(min(freqs),7) if freqs else '—')
        self.freq_kpis['highest'].value.setText(self._deep_fmt(max(freqs),7) if freqs else '—')
        self.freq_charts['modes'].refresh([(r['mode'],r['frequency']) for r in rows])
        self.freq_charts['ir'].refresh(sorted([(r['frequency'],r['ir_intensity']) for r in rows
                  if r.get('ir_intensity') is not None],key=lambda v:v[0]))
        self.freq_table.setRowCount(len(rows))
        for i,r in enumerate(rows):
            for j,key in enumerate(('mode','frequency','red_mass','force_constant','ir_intensity')):
                self.freq_table.setItem(i,j,QTableWidgetItem(self._deep_fmt(r.get(key),8)))
        self.freq_details.setPlainText('\n'.join([
            'Frequencies: Gaussian harmonic normal modes (cm⁻¹)',
            f"Number of printed modes: {len(rows)}",
            f"Imaginary (negative) modes: {sum(f<0 for f in freqs)}",
            f"IR intensities available: {sum(r.get('ir_intensity') is not None for r in rows)}",
            'IR peak heights are unbroadened intensity data, not a simulated spectrum.',
            'Frequency-only calculations at nonstationary structures require caution.',
            'These values are separate from L508 / QCNR solver convergence.',
        ]))
        self.freq_note.setText(f"{group['label']} · {len(rows)} harmonic modes"
                     + (' · HISTORY TRUNCATED' if group.get('truncated') else '')
                     + ' · Negative cm⁻¹ values are imaginary frequencies.')

    def export_freq_csv(self):
        group=self._freq_export
        if not group or not group['modes']:
            self.warn('No vibrational modes to export.');return
        filename,_=QFileDialog.getSaveFileName(self,'Export printed Gaussian frequencies',
            str(Path.cwd()/f'GauMon_Frequencies_{group["id"]}.csv'),'CSV files (*.csv)')
        if not filename:return
        try:
            fields=('mode','frequency','red_mass','force_constant','ir_intensity',
                    'raman_activity','depolar_p','depolar_u')
            with open(filename,'w',newline='',encoding='utf-8') as output:
                writer=csv.DictWriter(output,fieldnames=fields,extrasaction='ignore')
                writer.writeheader();writer.writerows(group['modes'])
        except OSError as error:self.warn(f'CSV export failed: {error}')

    @staticmethod
    def _style_table(table):
        table.setAlternatingRowColors(True)
        table.setSelectionBehavior(QAbstractItemView.SelectRows)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        table.horizontalHeader().setStretchLastSection(True)
        table.verticalHeader().hide()

    def selected_case(self):
        return ALL_CASES[max(0,self.test_table.currentRow())]

    def simulate_selected(self):
        self.simulate(self.selected_case().key)

    def simulate(self,case_key: str):
        if self.monitor.snapshot()['owned']:
            self.warn('A managed process is already running. Finish or stop it first.');return
        try:
            if self.mock_dir: self.mock_dir.cleanup()
            self.mock_dir=tempfile.TemporaryDirectory(prefix='gaumon_cpu_test_')
            dirpath=Path(self.mock_dir.name)
            gjf=write_case(case_key,dirpath)
            executable=Path(__file__).resolve().parents[1]/'examples'/'mock_g16'
            snap=self.monitor.launch(str(executable),str(gjf),str(dirpath),f'{case_key}_mock.log')
            self.test_note.setText(f'{case_key.upper()} synthetic calculation running (not Gaussian): {snap["log"]}')
            self.last_session=None
            self.last_plot_sig=None
            self.last_energy_sig=None
            self.last_cycles_sig=None
            self.last_links_sig=None
            self.show_page(0)
            self.tick()
        except Exception as exc:
            LOG.exception('Cannot run synthetic case')
            self.warn(str(exc))

    def open_stress_fixture(self):
        key=self.stress_choose.currentData()
        file=Path(__file__).resolve().parents[1]/'examples'/'stress_logs'/f'{key}.log'
        if not file.is_file():
            self.warn(f'Synthetic stress log unavailable: {file}')
            return
        try:
            self.attach_path(file)
            self.open_deep()
        except Exception as exc:
            LOG.exception('Cannot load synthetic stress fixture')
            self.warn(str(exc))

    def export_selected(self):
        dest=QFileDialog.getExistingDirectory(self,'Export test input',str(Path.cwd()))
        if not dest: return
        try:
            filename=write_case(self.selected_case().key,Path(dest))
            self.input_entry.setText(str(filename));self.work_entry.setText(str(filename.parent))
            self.log_entry.setText(filename.with_suffix('.log').name)
            self.show_page(1)
        except Exception as exc: self.warn(str(exc))

    def export_all(self):
        dest=QFileDialog.getExistingDirectory(self,'Export 7 test inputs',str(Path.cwd()))
        if not dest: return
        errors=[]
        for case in ALL_CASES:
            try:write_case(case.key,Path(dest))
            except Exception as exc:errors.append(f'{case.filename}: {exc}')
        if errors:self.warn('\n'.join(errors))
        else:QMessageBox.information(self,'Export complete',f'Seven Gaussian inputs saved to:\n{dest}')

    def launch(self):
        inp=Path(self.input_entry.text().strip()).expanduser()
        if not inp.is_file():self.warn('Choose a valid Gaussian .gjf/.com input first.');return
        issues=validate_input(inp.read_text(encoding='utf-8',errors='replace'))
        if issues:
            self.warn('Input validation failed:\n\n'+'\n'.join(issues));return
        try:
            self.monitor.launch(self.exec_entry.text().strip(),str(inp),self.work_entry.text().strip(),
                                self.log_entry.text().strip())
            self._clear_visual_caches()
            self._pending_enroll_log=self.monitor.snapshot()['log']
            self._pending_enroll_input=str(inp.resolve())
            self._autoload_input(self._pending_enroll_log,inp.resolve())
            self.show_page(0);self.tick()
        except Exception as exc:
            LOG.exception('Gaussian launch failed');self.warn(str(exc))

    def pick_running_job(self):
        jobs=self.monitor.jobs()
        if not jobs:
            QMessageBox.information(self,'Gaussian processes',
                'No Gaussian processes owned by this user were detected. '                'You can still attach directly to a log file without a PID.')
            return
        choices=[f"PID {j['pid']}   ·   {j['name']}   ·   {j['command'][:55]}" for j in jobs]
        item,ok=QInputDialog.getItem(self,'Running Gaussian processes',
            'Select a process (attachment is read-only):',choices,0,False)
        if ok and item:
            index=choices.index(item)
            self.pid_entry.setText(str(jobs[index]['pid']))
            self.show_page(1)
            self.input_hint.setText('PID selected. Choose its actual log file before attaching.')

    def attach(self):
        value=self.pid_entry.text().strip()
        try:
            pid=int(value) if value else None
            if pid is not None and pid<=0:raise ValueError('PID must be positive')
            self.monitor.attach(self.attach_entry.text().strip(),pid)
            self._clear_visual_caches()
            log_path=self.monitor.snapshot()['log']
            self._enroll_log(log_path)
            self._autoload_input(log_path)
            self.show_page(0);self.tick()
        except Exception as exc:
            LOG.exception('Attach failed');self.warn(str(exc))

    def attach_path(self, path: Path, pid: Optional[int] = None):
        self.attach_entry.setText(str(path))
        self.pid_entry.setText(str(pid) if pid is not None else '')
        self.attach()

    def stop(self,force: bool):
        if not self.monitor.snapshot()['owned']:
            self.warn('No GauMon-owned active Gaussian process. Attached external processes cannot be terminated.');return
        message='Force-kill only this GauMon-owned calculation?' if force else 'Request a graceful stop of this GauMon-owned calculation?'
        if QMessageBox.question(self,'Confirm stop',message,QMessageBox.Yes|QMessageBox.No,QMessageBox.No)!=QMessageBox.Yes:
            return
        try:self.monitor.stop(force=force)
        except Exception as exc:self.warn(str(exc))

    def warn(self,message):QMessageBox.warning(self,'GauMon',message)

    def tick(self):
        if self._pending_enroll_log:
            path=Path(self._pending_enroll_log)
            if path.is_file():
                self._enroll_log(path,'started-in-GauMon',self._pending_enroll_input)
                self._pending_enroll_log=None
                self._pending_enroll_input=None
        # A log attached after hours of calculation can be very large.
        # Keep parsing off the Qt event loop; only paint completed snapshots.
        if self._poll_future is None:
            self._poll_future=self._poll_executor.submit(self.monitor.poll)
            return
        if not self._poll_future.done():
            return
        try:
            state=self._poll_future.result()
        except Exception as exc:
            LOG.exception('Monitor poll failure')
            self.footer.setText(f'Monitor error: {exc}')
            self._poll_future=None
            return
        self._poll_future=self._poll_executor.submit(self.monitor.poll)
        if state.get('session') != self.monitor.session_id:
            return  # old completed future after the user selected another log
        active=state['state']
        self.header_state.setText('● '+active)
        self.header_state.setStyleSheet('color: #C55E5B; background: #FFF1EF;' if active in ('FAILED','INCOMPLETE') else
                                        'color: #B67937; background: #FFF4E7;' if active=='STALLED' else
                                        'color: #566B5E; background: #F0F2EE;' if active in ('READY','STARTING','ATTACHED') else
                                        'color: #258A62; background: #E5F4EA;')
        seconds=round(state['elapsed']);duration=f'{seconds//3600:02d}:{(seconds%3600)//60:02d}:{seconds%60:02d}'
        values={'state':active,'elapsed':duration,
                'energy':f'{state["energy"]:.6f}' if state['energy'] is not None else '—',
                'step':f'{state["step"]} / {state["step_max"] or "—"}' if state['step'] is not None else '—',
                'force':f'{state["force"]:.2e}' if state['force'] is not None else '—',
                'rmsd':f'{state["rmsd"]:.3f}' if state['rmsd'] is not None else '—'}
        for k,t in values.items():self.metrics[k].value.setText(t)
        self.info['status'].setText(active)
        self.info['pid'].setText(str(state['pid'] or '—'))
        self.info['mode'].setText(state['mode'] or '—')
        self.info['link'].setText(state['link'] or '—')
        self.info['link_desc'].setText((state.get('link_description') or '—')[:20] + ('…' if len(state.get('link_description') or '') > 20 else ''))
        self.info['link_desc'].setToolTip(state.get('link_description') or '—')
        self.link_peek.setText((state.get('link') or '—') + ' · ' + (state.get('link_description') or 'No description available'))
        self.current_link_label.setText('Current: ' + (state.get('link') or '—'))
        self.link_explanation.setText(state.get('link_description') or 'No link entered yet')
        self.info['scf_cycles'].setText(str(state['cycles']) if state['cycles'] is not None else '—')
        self.info['charge'].setText(f'{state["charge"]} / {state["mult"]}' if state['charge'] is not None else '—')
        route=state['route'] or '—';self.info['route'].setText(route[:23]+('…' if len(route)>23 else ''))
        self.info['route'].setToolTip(route)
        self.info['freq'].setText(str(state['frequencies']))
        self.info['geom'].setText(str(state['geometries']))
        self.info['disp'].setText(f'{state["displacement"]:.3f} Å' if state['displacement'] is not None else '—')
        self.info['cycle_avg'].setText((f'{state.get("recent_cycle_average_s", state.get("cycle_average_s")):.1f} s' if state.get('recent_cycle_average_s') is not None else '—'))
        eta_max = state.get('eta_to_max')
        self.info['eta_max'].setText((eta_max['duration'] if eta_max else '—'))
        self.info['eta_max'].setToolTip((eta_max['datetime'] if eta_max else 'Not enough completed optimization cycles yet'))
        records=state.get('eta_records',[])
        if records:
            r10=records[0]
            self.eta_peek.setText(r10['duration'])
            self.eta_peek_time.setText('Likely by '+r10['datetime'] if r10['datetime']!='—' else 'Awaiting reliable cycle timing')
            for rec in records:
                if rec['extra_cycles'] in self.eta_tiles:
                    duration_widget,due_widget=self.eta_tiles[rec['extra_cycles']]
                    duration_widget.setText(rec['duration'])
                    due_widget.setText(rec['datetime'] if rec['datetime']!='—' else 'Awaiting measured intervals')
        sample_n=state.get('cycle_sample_count',0)
        source=state.get('cycle_timing_source') or 'unknown'
        self.timeline_quality.setText(f'{sample_n} timed intervals · {source}' if sample_n else 'No verified cycle timings')
        self.cycle_count_label.setText(str(len(state.get('step_history',[])))+' observed steps')
        above_disp=(state['displacement'] is not None and state['displacement']>self.disp_limit)
        above_rmsd=(state['rmsd'] is not None and state['rmsd']>self.rmsd_limit)
        warning=(' • Displacement above threshold' if above_disp else '') + \
                (' • RMSD above threshold' if above_rmsd else '')
        self.quick_state.setText(state['detail'] + warning)
        all_diagnoses=diagnose_log(state.get('diagnostic_text','')+'\n'+state.get('tail',''))
        specific=next((d for d in all_diagnoses if d.severity=='ERROR' and
                       'Gaussian error termination' not in d.title),None)
        self.failure_banner.setVisible(active=='FAILED' and specific is not None)
        if specific and active=='FAILED':
            self.failure_title.setText('Failure identified · '+specific.title)
            self.failure_title.setToolTip(specific.description)
        self.metrics['rmsd'].value.setStyleSheet('color: #D76567;' if above_rmsd else '')
        self.footer.setText(f'{state["detail"]}  ·  {state["log"] or "No output file selected"}')
        self.stop_btn.setEnabled(state['owned']);self.kill_btn.setEnabled(state['owned'])
        self.start_btn.setEnabled(not state['owned'])
        plot_stamp=(state['session'],len(state['scf']),state['scf'][-1][0] if state['scf'] else 0,
                    len(state['forces']),len(state['displacements']),len(state['rmsds']),
                    len(state['prev_rmsds']),state['displacements'][-1][0] if state['displacements'] else 0)
        if plot_stamp != self.last_plot_sig:
            self.last_plot_sig=plot_stamp
            self.plots['energy'].refresh([[r[0],r[1]] for r in state['scf']])
            self.plots['force'].refresh(state['forces'])
            self.plots['disp'].refresh(state['displacements'])
            self.plots['rmsd'].refresh(state['rmsds'],state['prev_rmsds'])
        self.molecule.set_atoms(state['atoms'])
        self.geom_note.setText(f'{state["geometries"]} coordinate snapshots · RMSD vs initial {values["rmsd"]} Å · 2D projection only')
        # Avoid rebuilding analysis/log widgets when the data did not change.
        stamp=(state['session'],state.get('revision'),state['frequencies'],active)
        if stamp!=self.last_session:
            self.last_session=stamp
            energy_sig=(state['session'],len(state['scf']),state['scf'][-1][0] if state['scf'] else None)
            if energy_sig != self.last_energy_sig:
                self.last_energy_sig=energy_sig
                self.energy_table.setRowCount(len(state['scf']))
                for i,row in enumerate(state['scf']):
                    for j,val in enumerate((str(row[0]),f'{row[1]:.9f}',str(row[2]))):
                        self.energy_table.setItem(i,j,QTableWidgetItem(val))
            cycle_sig=(state['session'],len(state.get('step_history',[])),
                       state['step_history'][-1]['step'] if state.get('step_history') else None)
            if cycle_sig != self.last_cycles_sig:
                self.last_cycles_sig=cycle_sig
                steps=state.get('step_history',[])
                self.cycle_table.setRowCount(len(steps))
                for i,row in enumerate(steps):
                    fields=(str(row['step']), str(row['scf_cycles'] or '—'),
                            '—' if row['cycle_time_s'] is None else f"{row['cycle_time_s']:.1f} s",
                            row.get('timing_source','—'),
                            '—' if row['energy'] is None else f"{row['energy']:.7f}")
                    for j,val in enumerate(fields):
                        self.cycle_table.setItem(i,j,QTableWidgetItem(val))
            links=state.get('link_history',[])
            links_sig=(state['session'],len(links),links[-1]['link'] if links else None)
            if links_sig != self.last_links_sig:
                self.last_links_sig=links_sig
                self.link_table.setRowCount(len(links))
                for i,entry in enumerate(links):
                    logtime=entry.get('log_timestamp')
                    # Timestamp printed by Gaussian can be missing. Do not fabricate it.
                    time_label=time.strftime('%H:%M:%S',time.localtime(logtime)) if logtime else '—'
                    for j,val in enumerate((entry['link'],time_label,entry['description'] or '—')):
                        self.link_table.setItem(i,j,QTableWidgetItem(val))
            self.summary.setPlainText('\n'.join([
                f'Status:            {active}',
                f'Connection:        {state["mode"] or "not selected"}',
                f'Process PID:       {state["pid"] or "unavailable"}',
                f'Charge / mult:     {state["charge"]} / {state["mult"]}',
                f'Route:             {state["route"]}',
                f'Optimization:      {state["step"]} / {state["step_max"]}',
                f'Opt completed:     {state["opt_completed"]}',
                f'SCF cycles:        {state["cycles"]}',
                f'Frequency blocks:  {state["frequencies"]}',
                f'Imaginary modes:   {state["imaginary_modes"]}',
                f'Frequencies cm⁻¹:   {state["frequency_values"]}',
                f'Coordinates seen:  {state["geometries"]}',
                f'Gaussian finish:   {state["termination"]}',
                f'Observed elapsed:  {duration}',
            ]))
            if self.log_box.toPlainText()!=state['tail']:
                self.log_box.setPlainText(state['tail'])
                if self.follow.isChecked():
                    scroll=self.log_box.verticalScrollBar();scroll.setValue(scroll.maximum())
            diagnoses=diagnose_log(state.get('diagnostic_text','')+'\n'+state['tail'])
            if diagnoses:
                top=diagnoses[0]
                self.diagnostic_label.setPlainText(top.title+'\n'+top.description+'\nSuggested action: '+top.action)
                self.diagnostic_label.setToolTip(top.source_name or 'Automated rule, verify against actual Gaussian output')
            else:
                self.diagnostic_label.setPlainText('No recognized Gaussian errors in monitored output.')

        self._update_deep(state)
        self._update_qcnr(state)
        self._update_freq(state)

    def closeEvent(self,event):
        if self.mock_dir and self.monitor.snapshot()['owned']:
            QMessageBox.warning(self,'Mock calculation active',
                'Allow the synthetic test to finish or stop it before closing this window.')
            event.ignore();return
        if self.monitor.snapshot()['owned']:
            response=QMessageBox.question(self,'Calculation still active',
                'A GauMon-owned calculation is still running.\n'
                'Closing the monitor will NOT stop the Gaussian process. Close anyway?',
                QMessageBox.Yes|QMessageBox.No,QMessageBox.No)
            if response!=QMessageBox.Yes:
                event.ignore();return
        if self.mock_dir and not self.monitor.snapshot()['owned']:
            self.mock_dir.cleanup()
        self.timer.stop()
        self.hub.timer.stop()
        self._poll_executor.shutdown(wait=False, cancel_futures=True)
        event.accept()


def run_gui(*,mock_on_start=False,attach_log=None,attach_pid=None) -> int:
    # Force pyqtgraph's 2D QPainter path rather than any OpenGL backend.
    pg.setConfigOptions(useOpenGL=False, antialias=True, background='w', foreground=C['sub'])
    app=QApplication.instance() or QApplication([])
    app.setApplicationName('GauMon Desktop CPU')
    app.setStyle('Fusion')
    window=MainWindow(mock_on_start=mock_on_start)
    if os.environ.get('GAUMON_START_MAXIMIZED','').strip() == '1':
        window.showMaximized()  # XRDP fallback if full-screen integration is troublesome.
    else:
        window.showFullScreen()  # Full screen by default; F11/Esc restores maximized mode.
        window.full_btn.setText('Exit full screen  Esc')
    if attach_log is not None:
        QTimer.singleShot(100,lambda:window.attach_path(Path(attach_log),attach_pid))
    return app.exec_()
