"""Five opt-in scientific utility pages for the source-bundled GauMon desktop.
CPU-only Qt / pyqtgraph. Views observe logs; they never start Gaussian jobs.
"""
from __future__ import annotations
from datetime import datetime,timedelta
from pathlib import Path
import csv
import math
import os
import re
import statistics
import time
from PyQt5.QtCore import Qt,QTimer,QRect,QSize,pyqtSignal
from PyQt5.QtGui import QColor,QFont,QPainter,QTextCharFormat,QSyntaxHighlighter,QTextFormat,QTextOption
from PyQt5.QtWidgets import (QWidget,QVBoxLayout,QHBoxLayout,QLabel,QPushButton,
 QLineEdit,QTableWidget,QTableWidgetItem,QHeaderView,QFileDialog,QMessageBox,
 QPlainTextEdit,QSplitter,QCheckBox,QTextEdit,QComboBox,QSpinBox,QFrame,QAbstractItemView,QTabWidget)
try:
 import pyqtgraph as pg
except ImportError:
 pg=None
from .analytics import IncrementalLog,LINK_DESC,adjacent_step_durations
from .session import SessionStore,describe_observed_job
from .checker import check_input
from .input_context import explain_input,keyword_info,matching_input_for_log,TOKENS

GREEN='#276447'; MUTED='#66756A'; WHITE='#FFFFFF'; PALE='#F0F6F1'; AMBER='#A06A27'; RED='#B44740'


def subtitle(parent,title,description):
    parent.addWidget(label(title,20,True))
    note=label(description,12);note.setWordWrap(True);note.setStyleSheet('color:'+MUTED);parent.addWidget(note)

def label(content,size=12,bold=False):
    x=QLabel(str(content));f=x.font();f.setPointSize(size);f.setBold(bold);x.setFont(f);return x

def tbl(names):
    t=QTableWidget(0,len(names));t.setHorizontalHeaderLabels(names)
    t.setAlternatingRowColors(True);t.setSelectionBehavior(QAbstractItemView.SelectRows)
    t.setEditTriggers(QAbstractItemView.NoEditTriggers)
    t.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    t.verticalHeader().setVisible(False)
    return t

def fill(table,rows):
    table.setRowCount(len(rows))
    for i,row in enumerate(rows):
        for j,value in enumerate(row):
            table.setItem(i,j,QTableWidgetItem(str(value)))


def fmt(t):
    if t is None:return '—'
    if t<60:return f'{t:.1f} s'
    if t<3600:return f'{t/60:.1f} min'
    return f'{t/3600:.2f} h'

class PlotView(QWidget):
    def __init__(self,caption):
        super().__init__()
        v=QVBoxLayout(self);v.setContentsMargins(0,0,0,0)
        v.addWidget(label(caption,13,True))
        if pg:
            self.plot=pg.PlotWidget(background=WHITE)
            self.plot.showGrid(x=True,y=True,alpha=.12)
            self.plot.setMinimumHeight(155)
            self.plot.getAxis('left').setTextPen(MUTED)
            self.plot.getAxis('bottom').setTextPen(MUTED)
            v.addWidget(self.plot)
        else:
            self.plot=None;v.addWidget(label('Install pyqtgraph for interactive plotting.'))
    def bars(self,values,names):
        if not self.plot:return
        self.plot.clear()
        if not values:return
        bars=pg.BarGraphItem(x=list(range(len(values))),height=values,width=.65,brush=GREEN)
        self.plot.addItem(bars)
        self.plot.getAxis('bottom').setTicks([list(enumerate(names))]);self.plot.enableAutoRange()
    def line(self,x,y):
        if not self.plot:return
        self.plot.clear()
        if x:self.plot.plot(x,y,pen=pg.mkPen(GREEN,width=2.4),symbol='o',symbolSize=5,symbolBrush=GREEN)

class AnalysisHub(QWidget):
    """Manages only explicitly opened / launched jobs. Resumption reads logs only."""
    changed=pyqtSignal()
    def __init__(self,parent=None,history_path=None):
        super().__init__(parent)
        self.store=SessionStore(history_path)
        self.jobs={}
        self.selected=None
        self.timer=QTimer(self);self.timer.setInterval(1200);self.timer.timeout.connect(self.refresh)
        self.timer.start()
    def open_path(self,path,source='opened',input_path=None):
        p=Path(path).expanduser().resolve()
        if not p.is_file():raise FileNotFoundError(str(p))
        existing=next((r for r in self.store.items if r['path']==str(p)),None)
        if (not existing or existing['id'] not in self.jobs) and len(self.jobs)>=20:
            raise RuntimeError('At most 20 simultaneously attached logs are supported in this desktop session.')
        record=self.store.add(p,source,input_path=input_path)
        job_id=record['id']
        if job_id not in self.jobs:
            self.jobs[job_id]=IncrementalLog(p)
        self.selected=job_id
        self.refresh()
        self.changed.emit()
        return record
    def resume(self,job_id):
        record=next((r for r in self.store.items if r['id']==job_id),None)
        if record is None:raise ValueError('Unknown session record')
        return self.open_path(record['path'],record.get('source','opened'),record.get('input_path'))
    def choose(self,job_id):
        if job_id not in self.jobs:raise ValueError('Job is not active; select Resume first')
        self.selected=job_id;self.changed.emit()
    def current(self):
        job=self.jobs.get(self.selected)
        return job.snapshot if job else None
    def refresh(self):
        changed=False
        for ident,job in list(self.jobs.items()):
            before=(job.offset,job.snapshot.termination,job.last_poll_error)
            job.poll(max_bytes=256*1024)
            after=(job.offset,job.snapshot.termination,job.last_poll_error)
            if before!=after:changed=True
        if changed:self.changed.emit()

class TimePage(QWidget):
    def __init__(self,hub):
        super().__init__();self.hub=hub
        v=QVBoxLayout(self);v.setSpacing(9)
        subtitle(v,'Time graphs','Where did the computer actually spend time? Gaussian reports elapsed (“wall-clock”) and CPU seconds when each Link finishes. These are measurements for calculation modules—not for individual atoms or bonds. An unfinished Link has no completed-time record yet.')
        self.header=label('Open a log to see completed Gaussian Link timings.',13,True);v.addWidget(self.header)
        charts=QHBoxLayout();self.total=PlotView('Elapsed time by module (seconds)');self.sequence=PlotView('Completed Link calls in order (seconds)')
        charts.addWidget(self.total,1);charts.addWidget(self.sequence,1);v.addLayout(charts,2)
        self.table=tbl(['Link','What it does','Calls','Elapsed s','CPU s','Avg elapsed s','Share of measured elapsed'])
        v.addWidget(self.table,3)
        v.addWidget(label('Interpretation: Elapsed ≠ CPU. Elapsed sums cover completed Link calls only; parallel CPU seconds can exceed wall-clock time. “Slowest Link” identifies a computation stage, not a chemically problematic atom.',11))
        hub.changed.connect(self.redraw)
    def redraw(self):
        s=self.hub.current()
        if not s:self.header.setText('No session log selected');return
        totals=s.totals();measured=sum(x[0] for x in totals.values())
        keys=sorted(totals,key=lambda k:totals[k][0],reverse=True)
        self.header.setText(f'{Path(s.path).name} · {s.event_count} completed link calls · measured {fmt(measured)} · current L{s.current_link or "—"}')
        rows=[]
        for k in keys:
            elapsed,cpu,n=totals[k]
            rows.append((f'L{k}',LINK_DESC.get(k,'Other Gaussian module'),n,f'{elapsed:.2f}',f'{cpu:.2f}',f'{elapsed/n:.2f}',f'{elapsed/max(measured,1e-20)*100:.1f}%'))
        fill(self.table,rows)
        self.total.bars([totals[k][0] for k in keys[:10]],[f'L{k}' for k in keys[:10]])
        recent=list(s.events)[-35:]
        self.sequence.bars([e.elapsed for e in recent],[str(i+1) for i,e in enumerate(recent)])

class InputHighlighter(QSyntaxHighlighter):
    def __init__(self,doc):
        super().__init__(doc)
        self.styles=[]
        for pattern,color,bold in [(r'^\s*%[^\n]*',GREEN,True),(r'^\s*#[^\n]*', '#86562F',True),(r'^\s*--Link1--\s*$',RED,True), (r'^\s*[+-]?\d+\s+\d+\s*$', '#8B4674',True)]:
            fmt=QTextCharFormat();fmt.setForeground(QColor(color));fmt.setFontWeight(QFont.Bold if bold else QFont.Normal)
            self.styles.append((re.compile(pattern,re.I),fmt))
    def highlightBlock(self,text):
        for pattern,fmt in self.styles:
            m=pattern.search(text)
            if m:self.setFormat(m.start(),m.end()-m.start(),fmt)
        is_route=text.lstrip().startswith('#') or self.previousBlockState()==1
        if not text.strip() or text.strip().lower()=='--link1--':
            self.setCurrentBlockState(0)
        elif is_route:
            self.setCurrentBlockState(1)
            fmt=QTextCharFormat();fmt.setForeground(QColor(GREEN));fmt.setFontWeight(QFont.Bold)
            for match in TOKENS.finditer(text):
                if keyword_info(match[0]):
                    self.setFormat(match.start(),len(match[0]),fmt)
        else:
            self.setCurrentBlockState(0)

class CodeEditor(QPlainTextEdit):
    def __init__(self):
        super().__init__();self.gutter=QWidget(self);self.gutter.paintEvent=self.paint_gutter
        self.blockCountChanged.connect(self.update_gutter_width)
        self.updateRequest.connect(self.update_gutter)
        self.update_gutter_width();self.cursorPositionChanged.connect(self.highlight_cursor)
        self.setTabStopDistance(32)
        self.setFont(QFont('DejaVu Sans Mono',11))
        self.highlighter=InputHighlighter(self.document());self.error_lines=set()
        self.highlight_cursor()
    def gutter_width(self):return 14+self.fontMetrics().horizontalAdvance('9')*max(3,len(str(self.blockCount())))
    def update_gutter_width(self,*_):self.setViewportMargins(self.gutter_width(),0,0,0)
    def resizeEvent(self,e):
        super().resizeEvent(e);r=self.contentsRect();self.gutter.setGeometry(QRect(r.left(),r.top(),self.gutter_width(),r.height()))
    def update_gutter(self,rect,dy):
        if dy:self.gutter.scroll(0,dy)
        else:self.gutter.update(0,rect.y(),self.gutter.width(),rect.height())
    def paint_gutter(self,e):
        painter=QPainter(self.gutter);painter.fillRect(e.rect(),QColor('#EEF3EE'))
        block=self.firstVisibleBlock();num=block.blockNumber();top=int(self.blockBoundingGeometry(block).translated(self.contentOffset()).top())
        while block.isValid() and top<=e.rect().bottom():
            height=int(self.blockBoundingRect(block).height())
            if block.isVisible() and top+height>=e.rect().top():
                painter.setPen(QColor(RED if num+1 in self.error_lines else MUTED))
                painter.drawText(0,top,self.gutter.width()-5,height,Qt.AlignRight|Qt.AlignVCenter,str(num+1))
            block=block.next();top+=height;num+=1
    def highlight_cursor(self):
        extras=[]
        sel=QTextEdit.ExtraSelection();sel.format.setBackground(QColor('#EFF6EF'))
        sel.format.setProperty(QTextFormat.FullWidthSelection,True);sel.cursor=self.textCursor();sel.cursor.clearSelection();extras.append(sel)
        for l in list(self.error_lines)[:200]:
            if l<1 or l>self.blockCount():continue
            bl=self.document().findBlockByNumber(l-1)
            err=QTextEdit.ExtraSelection();err.format.setBackground(QColor('#FFF0ED'))
            err.format.setProperty(QTextFormat.FullWidthSelection,True)
            err.cursor=self.textCursor();err.cursor.setPosition(bl.position());err.cursor.clearSelection();extras.append(err)
        self.setExtraSelections(extras)
    def show_whitespace(self,show):
        opts=self.document().defaultTextOption()
        flags=opts.flags()
        if show:flags|=QTextOption.ShowTabsAndSpaces|QTextOption.ShowLineAndParagraphSeparators
        else:flags &= ~(QTextOption.ShowTabsAndSpaces|QTextOption.ShowLineAndParagraphSeparators)
        opts.setFlags(flags);self.document().setDefaultTextOption(opts)

class InputPage(QWidget):
    def __init__(self,hub):
        super().__init__();self.hub=hub;self.path=None
        v=QVBoxLayout(self);subtitle(v,'Input checker & editor','Check the Gaussian .gjf/.com structure before starting a job. The checker highlights suspicious lines, missing blank separators, and whitespace; it does NOT replace Gaussian’s own input validator. Editing never launches a calculation.')
        row=QHBoxLayout();self.path_lbl=label('Unsaved input',11);row.addWidget(self.path_lbl,1)
        for caption,action in [('Open .gjf',self.open_file),('Check input',self.check),('Save as…',self.save_as)]:
            b=QPushButton(caption);b.clicked.connect(action);row.addWidget(b)
        v.addLayout(row)
        self.whitespace=QCheckBox('Show spaces, tabs and line endings');self.whitespace.toggled.connect(self.toggle_ws);v.addWidget(self.whitespace)
        self.source=label('Source: no input selected. Auto-detection uses only same-name .gjf/.com next to an attached log.',11)
        self.source.setWordWrap(True);v.addWidget(self.source)
        split=QSplitter(Qt.Vertical);self.edit=CodeEditor();split.addWidget(self.edit)
        self.edit.cursorPositionChanged.connect(self._cursor_keyword)

        panel=QWidget();p=QVBoxLayout(panel);p.setContentsMargins(0,0,0,0)
        self.guides=QTabWidget()
        issue_tab=QWidget();issue_v=QVBoxLayout(issue_tab);issue_v.setContentsMargins(0,0,0,0)
        issue_v.addWidget(label('Select an issue to jump to that line',11))
        self.issues=tbl(['Severity','Line','Column','Explanation','Suggested correction'])
        self.issues.verticalHeader().setDefaultSectionSize(42)
        self.issues.setWordWrap(True)
        self.issues.cellDoubleClicked.connect(self.jump);issue_v.addWidget(self.issues)
        self.guides.addTab(issue_tab,'Input issues')
        keyword_tab=QWidget();key_v=QVBoxLayout(keyword_tab);key_v.setContentsMargins(0,0,0,0)
        key_v.addWidget(label('Recognized Gaussian keywords · double-click to jump to source',11))
        self.keyword_table=tbl(['Line','Keyword','Plain-English meaning','Important details'])
        self.keyword_table.setWordWrap(True)
        self.keyword_table.verticalHeader().setDefaultSectionSize(54)
        self.keyword_table.cellDoubleClicked.connect(self.jump_keyword)
        key_v.addWidget(self.keyword_table)
        self.guides.addTab(keyword_tab,'Keyword explanations')
        p.addWidget(self.guides)
        split.addWidget(panel);split.setSizes([430,210]);v.addWidget(split,1)
        self.current_keyword=label('Place the cursor on a recognized route keyword to see its explanation here.',11)
        self.current_keyword.setWordWrap(True);v.addWidget(self.current_keyword)
        self.status=label('Waiting for an input file.',11);v.addWidget(self.status)
        self._check_timer=QTimer(self);self._check_timer.setInterval(450)
        self._check_timer.setSingleShot(True);self._check_timer.timeout.connect(self.check)
        self.edit.textChanged.connect(self._check_timer.start)
    def toggle_ws(self,b):self.edit.show_whitespace(b)
    def open_file(self):
        file,_=QFileDialog.getOpenFileName(self,'Open Gaussian input',str(Path.cwd()),'Gaussian input (*.gjf *.com);;All files (*)')
        if file:self.load(file)
    def load(self,file,*,automatic=False):
        f=Path(file).expanduser().resolve()
        # Never discard unsaved scientific input when switching jobs.
        if self.edit.document().isModified():
            if automatic and self.path and f==self.path:
                return True
            choice=QMessageBox.question(self,'Unsaved input edits',
                'The editor contains unsaved changes. Replace them with the selected input?',
                QMessageBox.Yes|QMessageBox.No,QMessageBox.No)
            if choice!=QMessageBox.Yes:return False
        try:
            content=f.read_text(encoding='utf-8',errors='replace')
        except OSError as exc:
            if not automatic:QMessageBox.warning(self,'Input read failed',str(exc))
            self.source.setText(f'Input could not be loaded: {exc}')
            return False
        self.edit.setPlainText(content)
        self.edit.document().setModified(False)
        self.path=f;self.path_lbl.setText(str(f))
        self.source.setText(('Automatically loaded original input' if automatic else 'Opened input')+f': {f}. Editing does not modify the original unless you explicitly save to that path.')
        self.check()
        return True
    def autoload_for_log(self,log_path,explicit_input=None):
        """Choose explicit provenance or a same-stem adjacent file, never guess globally."""
        source=Path(explicit_input).expanduser().resolve() if explicit_input else matching_input_for_log(log_path)
        if source and source.is_file():
            return self.load(source,automatic=True)
        self.source.setText(f'No original .gjf/.com found beside {Path(log_path).name}. Gaussian logs may contain a route and geometry but not the original input file path. Use Open .gjf to select it; the currently shown editor text is not verified as belonging to this log.')
        return False
    def check(self):
        issues=check_input(self.edit.toPlainText())
        fill(self.issues,[(p.severity.upper(),p.line,p.column,p.message,p.fix) for p in issues])
        self.edit.error_lines={x.line for x in issues if x.severity=='error'}
        self.edit.highlight_cursor();self.edit.gutter.update()
        self._keyword_rows=explain_input(self.edit.toPlainText())
        fill(self.keyword_table,self._keyword_rows)
        for row,info in enumerate(self._keyword_rows):
            for col in (2,3):
                item=self.keyword_table.item(row,col)
                if item:item.setToolTip(info[col])
        self._cursor_keyword()
        e=sum(i.severity=='error' for i in issues);w=sum(i.severity=='warning' for i in issues)
        self.status.setText(f'{e} structural errors · {w} warnings · {len(issues)} total findings (heuristic checks only)')
        return issues
    def jump(self,row,col):
        cell=self.issues.item(row,1)
        if cell:
            block=self.edit.document().findBlockByNumber(max(0,int(cell.text())-1))
            cur=self.edit.textCursor();cur.setPosition(block.position());self.edit.setTextCursor(cur);self.edit.setFocus()
    def jump_keyword(self,row,col):
        if hasattr(self,'_keyword_rows') and 0<=row<len(self._keyword_rows):
            line=self._keyword_rows[row][0]
            block=self.edit.document().findBlockByNumber(line-1)
            cur=self.edit.textCursor();cur.setPosition(block.position())
            self.edit.setTextCursor(cur);self.edit.setFocus()
    def _cursor_keyword(self):
        if not hasattr(self,'current_keyword'):return
        cursor=self.edit.textCursor();line=cursor.block().text()
        if not any(row[0]==cursor.blockNumber()+1 for row in getattr(self,'_keyword_rows',())):
            self.current_keyword.setText('Select a green route keyword or open the Keyword explanations tab to see what it means.')
            return
        pos=max(0,cursor.positionInBlock())
        for m in TOKENS.finditer(line):
            if m.start()<=pos<=m.end():
                info=keyword_info(m[0])
                if info:
                    self.current_keyword.setText(f'{m[0]} — {info[0]}. {info[1]}')
                    return
        self.current_keyword.setText('Place the cursor on a green highlighted Gaussian keyword to see its meaning here. Full keyword list is in the Keyword explanations tab.')
    def save_as(self):
        default=(str(self.path.with_name(self.path.stem+'.checked'+self.path.suffix)) if self.path else str(Path.cwd()/'input.checked.gjf'))
        fname,_=QFileDialog.getSaveFileName(self,'Save checked input (new file)',default,'Gaussian input (*.gjf *.com)')
        if not fname:return
        f=Path(fname)
        if self.path and f.resolve()==self.path.resolve():
            if QMessageBox.question(self,'Overwrite original?','This overwrites the original Gaussian input. Continue?',QMessageBox.Yes|QMessageBox.No)!=QMessageBox.Yes:return
        try:
            f.write_text(self.edit.toPlainText(),encoding='utf-8');self.path=f;self.path_lbl.setText(str(f));self.edit.document().setModified(False);self.status.setText(f'Saved {f}')
        except OSError as exc:QMessageBox.warning(self,'Save failed',str(exc))

class ETAPage(QWidget):
    def __init__(self,hub):
        super().__init__();self.hub=hub
        v=QVBoxLayout(self);subtitle(v,'Advanced time estimator','Measured time and predicted time are different. GauMon measures completed Gaussian Link calls; next-step forecasts require new optimization steps observed while this session is open. The maximum allowed optimization steps is NOT the expected number of steps to convergence.')
        self.now_lbl=label('No log selected.',14,True);v.addWidget(self.now_lbl)
        self.note=label('Attach a log and wait for newly observed optimization steps.',12);self.note.setWordWrap(True);v.addWidget(self.note)
        self.stage_lbl=label('Current Gaussian module duration estimate: —',12);self.stage_lbl.setWordWrap(True);v.addWidget(self.stage_lbl)
        self.chart=PlotView('Recent observed time per optimization step');v.addWidget(self.chart,2)
        self.table=tbl(['Additional steps','Provisional duration','Local projected time','Observed pace range (not calibrated)']);v.addWidget(self.table,2)
        v.addWidget(label('This is not a guaranteed job-finish prediction. Frequency analysis, Link1 stages, XQC transitions and early convergence can change total runtime substantially.',11))
        hub.changed.connect(self.redraw)
    def redraw(self):
        s=self.hub.current()
        if not s:return
        eta=s.eta();self.now_lbl.setText(f'{Path(s.path).name} · optimizer step {s.steps[-1][0] if s.steps else "—"}/{s.steps[-1][1] if s.steps else "—"}')
        self.note.setText(eta.get('reason') or f'Median observed step duration: {fmt(eta["pace_s"])} from {eta["sample_count"]} intervals. No calibrated completion probability is available.')
        current=s.current_link
        elapsed=[e.elapsed for e in s.events if e.link==current]
        if current and elapsed:
            self.stage_lbl.setText(f'Current L{current}: recent completed calls took a median {fmt(statistics.median(elapsed))} (sample n={len(elapsed)}). This is observed WHOLE-MODULE duration, not time remaining.')
        elif current:
            self.stage_lbl.setText(f'Current L{current}: no completed calls of this module to estimate its duration.')
        else:
            self.stage_lbl.setText('Current module: unknown from the attached log.')
        durations=adjacent_step_durations(s.observed_step_times)
        self.chart.line([step for step,duration in durations],
                        [duration for step,duration in durations])
        out=[]
        if eta['available']:
            for item in eta['horizons']:
                t=datetime.now()+timedelta(seconds=item['seconds'])
                out.append((f'+{item["steps"]}',fmt(item['seconds']),t.strftime('%Y-%m-%d %H:%M'),f'{fmt(item["earliest_rough_s"])} to {fmt(item["latest_rough_s"])}'))
        fill(self.table,out)

class WarningsPage(QWidget):
    def __init__(self,hub):
        super().__init__();self.hub=hub
        v=QVBoxLayout(self)
        subtitle(v,'Warnings & errors','Understand what Gaussian actually printed. RED = explicit error evidence; AMBER = warning or suspected issue. Recommendations are suggestions to verify—not automatic changes to a scientific calculation.')
        self.summary=label('No job selected.');v.addWidget(self.summary)
        self.table=tbl(['Level','Log line','Gaussian evidence','Meaning / next check']);self.table.setWordWrap(True)
        self.table.verticalHeader().setDefaultSectionSize(55)
        v.addWidget(self.table,1)
        guide=label('Example: “FormBX had a problem” generally indicates an internal-coordinate optimization problem; inspect nearby failed angles and torsions before changing the optimization method. A generic “Error termination” only confirms that the job stopped abnormally.',11)
        guide.setWordWrap(True);v.addWidget(guide)
        hub.changed.connect(self.redraw)
    def redraw(self):
        s=self.hub.current()
        if not s:return
        self.summary.setText(f'{Path(s.path).name} · {len(s.problems)} recognized findings · status {s.termination or "not confirmed"}')
        fill(self.table,[(p.level.upper(),p.line,p.evidence,p.explanation) for p in s.problems])

class JobsPage(QWidget):
    def __init__(self,hub,selection_changed):
        super().__init__();self.hub=hub;self.selection_changed=selection_changed
        v=QVBoxLayout(self)
        subtitle(v,'Session jobs','Only jobs explicitly opened or launched through GauMon appear here. Saved paths survive a desktop restart; select Resume to reattach read-only. No external calculation is restarted and no unrelated process is scanned or stopped.')
        row=QHBoxLayout()
        for name,fn in [('Open existing log…',self.open),('Resume selected',self.resume),('Select for graphs',self.select),('Remove record',self.remove)]:
            b=QPushButton(name);b.clicked.connect(fn);row.addWidget(b)
        row.addStretch(1);v.addLayout(row)
        self.table=tbl(['File','Session status','Source','Completed links','Optimization steps','Saved date'])
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.cellDoubleClicked.connect(lambda *_: self.select());v.addWidget(self.table,1)
        self.footer=label('Saved records are not proof that a process is still running.',11);v.addWidget(self.footer)
        hub.changed.connect(self.redraw);self.redraw()
    def active_id(self):
        row=self.table.currentRow()
        return self.ids[row] if row>=0 and row<len(self.ids) else None
    def open(self):
        p,_=QFileDialog.getOpenFileName(self,'Open Gaussian log',str(Path.cwd()),'Gaussian logs (*.log *.out *.txt);;All files (*)')
        if p:
            try:self.hub.open_path(p);self.selection_changed()
            except (OSError,RuntimeError) as err:QMessageBox.warning(self,'Cannot read log',str(err))
    def resume(self):
        ident=self.active_id()
        if ident:
            try:self.hub.resume(ident);self.selection_changed()
            except (OSError,RuntimeError) as err:QMessageBox.warning(self,'Cannot resume',str(err))
    def select(self):
        ident=self.active_id()
        if ident in self.hub.jobs:self.hub.choose(ident);self.selection_changed()
        elif ident:QMessageBox.information(self,'Resume first','This saved record is inactive. Click Resume selected to reattach read-only.')
    def remove(self):
        ident=self.active_id()
        if ident:
            self.hub.jobs.pop(ident,None)
            if self.hub.selected==ident:self.hub.selected=None
            self.hub.store.remove(ident);self.hub.changed.emit()
    def redraw(self):
        ident=self.active_id() if hasattr(self,'ids') else None
        self.ids=[x['id'] for x in self.hub.store.items]
        rows=[]
        for x in self.hub.store.items:
            engine=self.hub.jobs.get(x['id'])
            st=getattr(engine,'snapshot',None)
            status=describe_observed_job(engine)
            rows.append((Path(x['path']).name,status,x.get('source','opened'),
                         st.event_count if st else '—',len(st.steps) if st else '—',x['added'][:16]))
        selected=self.table.currentRow();fill(self.table,rows)
        for idx,record in enumerate(self.hub.store.items):
            cell=self.table.item(idx,0)
            if cell:cell.setToolTip(record['path'])
        try:
            ix=self.ids.index(ident);self.table.selectRow(ix)
        except (ValueError,TypeError):
            if selected>=0 and selected<len(rows):self.table.selectRow(selected)
        self.footer.setText(f'{len(rows)} explicitly enrolled jobs · {len(self.hub.jobs)} attached log readers · session file: {self.hub.store.path}')
