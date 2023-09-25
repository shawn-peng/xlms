import copy
import sys
import time
import threading
from typing import *
import multiprocessing as mp

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5 import QtCore

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from mixture_model import MixtureModel
from xlms import XLMS_Dataset
# from run import get_cons_str
from myutils import *


# from plotter import ProcessPlotter


class MplCanvas(FigureCanvasQTAgg):
    updatePlot = pyqtSignal(np.ndarray, list, np.ndarray, AttrObj)

    def __init__(self, parent=None, width=16, height=9, dpi=100):
        fig = Figure(figsize=(width, height), layout='tight')
        self.fig = fig
        # self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)
        # self.resize(width*dpi, height*dpi)
        self.updatePlot.connect(self._plot)
        self.parent = parent
        # self.inPlotting = False
        self.stopPlotting = False

    def plot(self, *args):
        if not self.stopPlotting:
            print('emitting plot signal')
            self.updatePlot.emit(*args)

    def _plot(self, *args):
        with QSignalBlocker(self):
            print('updating plot...')
            self.stopPlotting = True
            self.parent.model._plot(*args, fig=self.fig)
            self.draw()
            self.stopPlotting = False
            print('plot updated')


class Slider(QWidget):
    valueChanged = pyqtSignal()

    def __init__(self, lower=0, upper=100, val=50, step=1):
        super().__init__()
        self.stepsize = step
        lower = int(lower / step)
        upper = int(upper / step)

        b = QHBoxLayout()
        b.setAlignment(Qt.AlignRight)
        self.setLayout(b)
        self.l = QLineEdit()
        self.l.setFixedWidth(40)
        self.l.textChanged.connect(self.valueinput)
        # self.l.setFixedHeight(14)
        s = QSlider(Qt.Horizontal)
        self.s = s
        s.setValue(val)
        # s.setSingleStep(step)
        # s.setPageStep(step * 10)
        self.valuechange()
        s.setMinimum(lower)
        s.setMaximum(upper)
        s.setFixedWidth(100)
        s.valueChanged.connect(self.valuechange)
        s.setTickPosition(QSlider.TicksAbove)
        b.addWidget(self.l)
        b.addWidget(self.s)

    def valuechange(self):
        v = self.s.value() * self.stepsize
        self.l.setText(('%8.2f' % v))
        self.valueChanged.emit()

    def valueinput(self):
        v = float(self.l.text()) / self.stepsize
        self.s.setValue(int(v))
        self.valueChanged.emit()

    @property
    def val(self):
        # return self.s.value() * self.stepsize
        return float(self.l.text())

    def set_value(self, v, direction='input'):
        v = v / self.stepsize
        if direction == 'input':
            self.l.setText(('%8.2f' % v))
            self.s.setValue(int(v))
        elif direction == 'output':
            with QSignalBlocker(self):
                self.l.setText(('%8.2f' % v))
                self.s.setValue(int(v))


class Worker(QObject):
    target: Callable
    finished = pyqtSignal()
    progress = pyqtSignal(int)

    def run(self):
        """Long-running task."""
        self.target()
        self.finished.emit()


class WorkerThread(QThread):
    finished = pyqtSignal()

    def __init__(self):
        self.model = None
        self.X = None
        super().__init__()

    def run(self):
        print('worker thread started')
        self.model.fit(self.X)
        print('model fit done')
        self.finished.emit()


datasets = [
    'alban',
    'Alinden',
    'ALott',
    'CPSF',
    'D1810',
    'ecoli_xl',
    # 'KKT4',
    'MS2000225',
    'peplib',
    'QE',
    'RPA',
    # 'Gordon',
]


class ParamWindow(QWidget):
    model_update_signal = pyqtSignal(AttrObj)
    model_stop_signal = pyqtSignal()

    def __init__(self, param_settings):
        global app
        super().__init__()
        w = self
        w.setWindowTitle('Parameters')

        self.model_update_signal.connect(self._model_updates)
        self.model_stop_signal.connect(self._model_stopped)
        self.model_event_handlers = {
            # 'plot':    self.update_plot,
            'update':  self.model_updates,
            'stopped': self.model_stopped,
        }

        mainbox = QVBoxLayout()
        w.setLayout(mainbox)

        settings_widget = QWidget()
        mainbox.addWidget(settings_widget)
        settings_box = QHBoxLayout()
        settings_widget.setLayout(settings_box)

        btnlist = QWidget()
        btnbox = QVBoxLayout()
        btnlist.setLayout(btnbox)
        settings_box.addWidget(btnlist)

        iconholder = QLabel()
        iconholder.setAlignment(Qt.AlignCenter)
        iconholder.setPixmap(app.style().standardPixmap(QStyle.SP_MessageBoxWarning))
        iconholder.hide()
        self.iconholder = iconholder
        btnbox.addWidget(iconholder)

        comboBox = QComboBox(self)
        for dataset in datasets:
            comboBox.addItem(dataset)
        btnbox.addWidget(comboBox)
        comboBox.setCurrentText('Alinden')
        comboBox.activated[str].connect(self.dataset_choice)

        vis_btn = QPushButton('Visualize')
        vis_btn.clicked.connect(self.visualize_model)
        btnbox.addWidget(vis_btn)

        run_btn = QPushButton('Run')
        run_btn.clicked.connect(self.run_model)
        self.run_btn = run_btn
        btnbox.addWidget(run_btn)

        stop_btn = QPushButton('Stop')
        stop_btn.clicked.connect(self.stop_model)
        btnbox.addWidget(stop_btn)

        btnbox.addStretch(1)

        condlist = QWidget()
        condbox = QVBoxLayout()
        condlist.setLayout(condbox)
        settings_box.addWidget(condlist)

        self.options = [
            'gaussian_model',
            'ic2_comp',
            'show_plotting',
            'weights',
            'mode',
            'pdf',
            'cdf',
            'weighted_pdf',
        ]
        self.options = {k: QCheckBox(k) for k in self.options}
        for btn in self.options.values():
            condbox.addWidget(btn)

        self.options['show_plotting'].setChecked(True)
        self.options['ic2_comp'].setChecked(True)
        self.options['mode'].setChecked(True)
        self.options['pdf'].setChecked(True)

        canvas = MplCanvas(self)
        self.canvas = canvas
        mainbox.addWidget(canvas)

        alpha_base = 1
        sls = {'C': alpha_base, 'IC': alpha_base, 'IC2': alpha_base, 'I1': -alpha_base, 'I2': -alpha_base}
        self.model = MixtureModel(sls, ic2_comp=True, plot_func=self.canvas.plot,
                                  event_notify_func=self.model_events)

        grid = QGridLayout()

        i = 3
        for _, params in param_settings.items():
            for pa in params.keys():
                grid.addWidget(QLabel(pa), i + 1, 1)
                i += 1
            break

        j = 1
        self.sliders = {}
        for i in range(self.model.n_samples):
            grid.addWidget(QLabel(f'S{i + 1}'), i + 2, 1)
        self.canvas.stopPlotting = True
        for c, params in param_settings.items():
            grid.addWidget(QLabel(c), 1, j + 1)
            for i in range(self.model.n_samples):
                if c in self.model.weights[i]:
                    s = Slider(0, 1, step=0.01)
                    pa = f'w{i + 1}'
                    self.sliders[(c, pa)] = s
                    s.valueChanged.connect((lambda c, pa: lambda: self.param_change(c, pa))(c, pa))
                    grid.addWidget(s, i + 2, j + 1)
            i += 2
            for pa, (v, lo, hi) in params.items():
                s = Slider(lo, hi, v, step=0.05)
                self.sliders[(c, pa)] = s
                s.valueChanged.connect((lambda c, pa: lambda: self.param_change(c, pa))(c, pa))
                grid.addWidget(s, i + 1, j + 1)
                i += 1
            j += 1
        self.canvas.stopPlotting = False

        param_grid = QWidget()
        param_grid.setLayout(grid)
        settings_box.addWidget(param_grid)

        self.thread = WorkerThread()

        self.thread.finished.connect(self.model_finished)


    def dataset_choice(self, dataset_name):
        dataset = XLMS_Dataset(dataset_name)
        dataset.mat = dataset.mat[:, dataset.mat[1, :] != 0]
        self.dataset = dataset

        n = dataset.mat.shape[1]
        SAMPLE_SIZE = 20000
        if SAMPLE_SIZE > 0:
            nsample = min(SAMPLE_SIZE, n)
            np.random.seed(42)
            rind = np.random.choice(np.arange(n), nsample, replace=False)
            dataset.mat = dataset.mat[:, rind]

        config = self.get_config()
        sls = {}
        alpha_base = 2
        sls = {'C': alpha_base, 'IC': alpha_base, 'IC2': alpha_base, 'I1': -alpha_base, 'I2': -alpha_base}
        self.model = MixtureModel(sls, ic2_comp=True, plot_func=self.canvas.plot,
                                  event_notify_func=self.model_events)
        X = self.dataset.mat.T
        self.X = X
        self.model.init_model(X)

    def showEvent(self, event):
        self.canvas.resize(1600, 700)
        dataset_name = 'Alinden'
        dataset = XLMS_Dataset(dataset_name)
        dataset.mat = dataset.mat[:, dataset.mat[1, :] != 0]
        self.dataset = dataset

        n = dataset.mat.shape[1]
        SAMPLE_SIZE = 20000
        if SAMPLE_SIZE > 0:
            nsample = min(SAMPLE_SIZE, n)
            np.random.seed(42)
            rind = np.random.choice(np.arange(n), nsample, replace=False)
            dataset.mat = dataset.mat[:, rind]

        config = self.get_config()
        sls = {}
        alpha_base = 2
        sls = {'C': alpha_base, 'IC': alpha_base, 'IC2': alpha_base, 'I1': -alpha_base, 'I2': -alpha_base}
        self.model = MixtureModel(sls, ic2_comp=True, plot_func=self.canvas.plot, show_plotting=True,
                                  event_notify_func=self.model_events)
        X = self.dataset.mat.T
        self.X = X
        self.model.init_model(X)

        event.accept()
        pass


    def model_stopped(self):
        self.model_stop_signal.emit()

    def _model_stopped(self):
        self.run_btn.setEnabled(True)
        msg_box = QMessageBox(QMessageBox.Information, 'Info', 'Model stopped')
        msg_box.exec()
        self.reset_worker_thread()

    def model_finished(self):
        self.run_btn.setEnabled(True)
        msg_box = QMessageBox(QMessageBox.Information, 'Info', 'Model converged')
        msg_box.exec()
        self.reset_worker_thread()

    def reset_worker_thread(self):
        self.thread.wait()
        print('worker thread ended')
        self.thread = WorkerThread()
        self.thread.finished.connect(self.model_finished)

    def weight_change(self, c, si):
        v = self.get_weight()

    def param_change(self, c, pa):
        v = self.get_param(c, pa)
        if pa[0] == 'w':
            si = int(pa[1]) - 1
            self.model.weights[si][c] = v
        else:
            dist = self.model.all_comps[c]
            dist.__setattr__(pa, v)
            dist.calc_alt_params()
        self.visualize_model()
        if not self.model.check_constraints():
            self.iconholder.show()
            self.run_btn.setEnabled(False)
        else:
            self.iconholder.hide()
            self.run_btn.setEnabled(True)

    def get_option(self, opt):
        return self.options[opt].isChecked()

    def get_param(self, c, pa):
        return float(self.sliders[(c, pa)].val)

    def run(self):
        self.show()


    def run_model(self):
        print('setting model...')
        alpha_base = 2
        sls = {'C': alpha_base, 'IC': alpha_base, 'IC2': alpha_base, 'I1': -alpha_base, 'I2': -alpha_base}
        self.model = MixtureModel(sls, **self.get_config())
        self.model.init_range(self.X)
        for (c, pa), s in self.sliders.items():
            if pa[0] == 'w':
                i = int(pa[1]) - 1
                self.model.weights[i][c] = s.val
            else:
                dist = self.model.all_comps[c]
                dist.__setattr__(pa, s.val)
                dist.calc_alt_params()
        self.model.starting_pos = self.model.frozen()
        self.model.create_constraints()
        self.model.initialized = True
        print('frozen_model', self.model.frozen())

        self.thread.model = self.model
        self.thread.X = self.X
        print('thread', self.thread)
        print('thread.model', self.thread.model)
        print('thread.X', self.thread.X)

        print('starting model')
        self.thread.start()
        self.run_btn.setEnabled(False)

    def visualize_model(self):
        config = self.get_config()
        alpha_base = 1
        sls = {'C': alpha_base, 'IC': alpha_base, 'IC2': alpha_base, 'I1': -alpha_base, 'I2': -alpha_base}
        model = MixtureModel(sls, config)
        for c, dist in model.all_comps.items():
            for pa in ['mu', 'sigma', 'alpha']:
                dist.__setattr__(pa, self.get_param(c, pa))
            dist.calc_alt_params()
        X = self.X
        self.model.plot(X, [self.model.log_likelihood(X)], self.model.sep_log_likelihood(X))

    def get_config(self):
        config = {
            'gaussian_model':    self.get_option('gaussian_model'),
            'ic2_comp':          self.get_option('ic2_comp'),
            'show_plotting':     self.get_option('show_plotting'),
            'plot_interval':     10,
            'tolerance':         1e-8,
            'constraints':       [
                cons for cons in ['weights', 'mode', 'pdf', 'cdf', 'weighted_pdf'] if self.get_option(cons)
            ],
            'plot_func':         self.canvas.plot,
            'event_notify_func': self.model_events,
        }
        return config

    def model_events(self, *args):
        # print('model event', args)
        event = args[0]
        self.model_event_handlers[event](*args[1:])

    def model_updates(self, frozen_model):
        # print('sending update signal')
        self.model_update_signal.emit(frozen_model)

    def _model_updates(self, frozen_model):
        for c, dist in frozen_model.all_comps.items():
            for i in range(frozen_model.n_samples):
                if c in frozen_model.weights[i]:
                    self.sliders[(c, f'w{i + 1}')].set_value(frozen_model.weights[i][c], direction='output')
            self.sliders[(c, 'mu')].set_value(dist.mu, direction='output')
            self.sliders[(c, 'sigma')].set_value(dist.sigma, direction='output')
            self.sliders[(c, 'alpha')].set_value(dist.alpha, direction='output')

    def _model_run(self):
        alpha_base = 2
        sls = {'C': alpha_base, 'IC': alpha_base, 'IC2': alpha_base, 'I1': -alpha_base, 'I2': -alpha_base}
        self.model = MixtureModel(sls, **self.get_config())
        for (c, pa), s in self.sliders.items():
            dist = self.model.all_comps[c]
            dist.__setattr__(pa, s.val)
            dist.calc_alt_params()
        self.model.initialized = True
        print('fitting model to data')
        ll, lls = self.model.fit(self.dataset.mat.T)
        self.run_btn.setEnabled(True)

    def stop_model(self):
        self.model.stop()
        self.run_btn.setEnabled(True)


if __name__ == '__main__':
    QApplication.setStyle('fusion')
    app = QApplication(sys.argv)

    param_ranges = {'mu': (150, -100, 300), 'sigma': (50, 1, 200), 'alpha': (2, -10, 10)}
    w = ParamWindow({c: param_ranges for c in ['C', 'IC', 'I1', 'IC2', 'I2']})
    w.resize(1200, 1000)
    w.run()

    app.exec()
    print('app ended somehow')
