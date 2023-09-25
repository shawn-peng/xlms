import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import skew_normal
import normal_gpu
from constraints import *
from param_binary_search import *
from kanren import *
import sympy
import time
from myutils import *
from plotter import *

from collections import defaultdict
from typing import *
import multiprocessing as mp

SN = skew_normal.SkewNormal
N = normal_gpu.Normal


class PlotWrapper:
    def __init__(self, func):
        # self.static = {}
        self.func = func

    def __call__(self, *args):
        cmd = args[0]
        if cmd == 'set_static':
            name, val = args[1:]
            self.static[name] = val
        elif cmd == 'plot':
            X = self.static['X']
            self.func(X, *args[1:])


class MixtureModelBase:
    def __init__(self, plot_func=None, title=None, **kwargs):
        self.plot_func = plot_func
        # if not plot_pipe:
        #     plot_pipe, plotter_pipe = mp.Pipe()
        #     self.plotter = ProcessPlotter((self._plot))
        #     self.plot_process = mp.Process(
        #         target=self.plotter, args=(plotter_pipe,), daemon=True)
        #     self.plot_process.start()
        # self.plot_pipe = plot_pipe
        self.title = title
        self.plotting_X = None
        self.fig = None

    def frozen(self):
        return AttrObj({
            **{
                k: deepcopy(self.__getattribute__(k))
                for k in ['binwidth', 'plotstep', 'n_samples', 'weights', 'comps', 'all_comps', 'starting_pos', 'title']
            }})

    def __del__(self):
        # self.plot_pipe.send(None)
        pass

    def init_range(self, X):
        xmax = np.max(X)
        xmin = np.min(X)
        self.n_xsamples = 200
        xstep = (xmax - xmin) / self.n_xsamples
        # x = np.arange(xmin, xmax + self.plotstep, self.plotstep)
        self.xrange = np.arange(xmin, xmax + xstep, xstep)

    # def plot(self, X, lls, slls, finished=False):
    #     send = self.plot_pipe.send
    #
    #     m = TimeMeter()
    #     if X is not self.plotting_X:
    #         self.plotting_X = X
    #         data = ('set_static_pos_arg', 0, copy(X))
    #         send(data)
    #     data = ('plot', None, lls, slls, self.frozen())  # None is a placeholder for static arg
    #     print('Sending plot data')
    #     send(data)
    #     print('done', m.read())
    def plot(self, X, lls, slls):
        if self.plot_func:
            return self.plot_func(X, lls, slls, self.frozen())
        if not self.fig:
            self.fig = plt.figure(figsize=(16, 9))
            plt.ion()
            plt.show()
        fig = self.fig
        return self._plot(X, lls, slls, self.frozen(), fig)

    @staticmethod
    def _plot(X, lls, slls, frozen_model, fig=None):
        # print('frozen_model', frozen_model)
        # print('frozen_model starting_pos', frozen_model.starting_pos)
        # frozen_model = AttrObj(frozen_model)
        xmax = np.max(X)
        xmin = -50
        # xmin = x.min()
        # xmax = x.max()
        x = np.arange(xmin, xmax + frozen_model.plotstep, frozen_model.plotstep)
        bins = np.arange(xmin, xmax + frozen_model.binwidth, frozen_model.binwidth)
        # distc = {}
        axs = []
        idmap = defaultdict(lambda: len(idmap))
        nsub = frozen_model.n_samples + 1
        # fig.clf()
        if not fig.axes:
            axs = fig.subplots(nsub, 1)
        else:
            axs = fig.axes
        for i in range(len(frozen_model.comps)):
            # print(f'plotting score {i + 1}')
            # xi = X[:, i]
            # ax = fig.subplot(nsub, 1, i + 1, label=f's{i + 1}')
            ax = axs[i]
            ax.cla()
            yi = np.zeros_like(x, dtype=float)
            legends = []
            # cnts, _ = np.histogram(X[:, i], bins=bins, density=True)
            # ymax = cnts.max()
            for j, (cname, cdist) in enumerate(frozen_model.comps[i].items()):
                # ax.plot(cdist.pdf(x), linestyle='--')
                yj = frozen_model.weights[i][cname] * cdist.pdf(x)
                yi += yj
                ax.plot(x, yj, c='C%d' % idmap[cname])
                legends.append(cname)
            ymax = yi.max()
            for j, (cname, cdist) in enumerate(frozen_model.comps[i].items()):
                scdist = frozen_model.starting_pos.comps[i][cname]
                # print(f'plotting component {cname}')
                str_params = f'{frozen_model.weights[i][cname]:.2f}' \
                             f' {cname}: m={cdist.mu:.2f}({scdist.mu:.2f}),' \
                             f' s={cdist.sigma:.2f}({scdist.sigma:.2f}),' \
                             f' a={cdist.alpha:.2f}({scdist.alpha:.2f})'
                ax.text(xmin, ymax * (1 - 0.1 * (j + 1)), str_params)
            # print('draw mixture')
            ax.plot(x, yi, c='C%d' % idmap['mixture'])
            legends.append(f'mixture{i + 1}')
            # print('draw hist')
            ax.hist(X[:, i], bins=bins, density=True, facecolor='w', ec='k', lw=1)
            # fig.bars()
            # h = info.hist[i]
            # fig.bar(h[1][:-1], h[0], facecolor='none', edgecolor='k')
            # print('draw texts')
            ax.text(xmax * 0.76, ymax * 0.9, f'll = {slls[i]:.5f}')
            # print('draw legends')
            ax.legend(legends)
        # ax = fig.subplot(nsub, 1, nsub, label=f'll_curve')
        ax = axs[-1]
        ax.cla()
        ax.plot(lls[1:])
        if lls:
            ax.text(0, lls[-1] - 0.1, f'll = {lls[-1]:.5f}')
        ax = fig.axes[0]
        # plt.axes(ax)
        if frozen_model.title:
            ax.set_title(frozen_model.title)
        plt.pause(0.01)

    def log_likelihood(self, X):
        ll = np.log(self.likelihood(X))
        ll[np.isinf(ll)] = ll[~np.isinf(ll)].min()
        return ll.mean()

    def sep_log_likelihood(self, X):
        ll = np.log(self.likelihood(X))
        ll[np.isinf(ll)] = ll[~np.isinf(ll)].min()
        return ll.mean(1)

    def likelihood(self, X):
        X = np.array(X)
        n, d = X.shape
        assert d == 2
        pj = list(map(lambda c: np.zeros((n, len(c))), self.comps))
        p = []
        for i in range(len(self.comps)):
            ws = self.weights[i]
            for j, (cname, cdist) in enumerate(self.comps[i].items()):
                pj[i][:, j] = ws[cname] * cdist.pdf(X[:, i])
            p0 = np.sum(pj[i], 1)
            p.append(p0)
        return p
