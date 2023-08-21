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

from collections import defaultdict
from typing import *

SN = skew_normal.SkewNormal
N = normal_gpu.Normal


class MixtureModelBase:

    def plot(self, X, lls, slls):
        xmax = np.max(X)
        xmin = -100
        x = np.arange(xmin, xmax + self.plotstep, self.plotstep)
        bins = np.arange(xmin, xmax + self.binwidth, self.binwidth)
        # distc = {}
        axs = []
        idmap = defaultdict(lambda: len(idmap))
        nsub = self.n_samples + 1
        for i in range(len(self.comps)):
            # xi = X[:, i]
            ax = plt.subplot(nsub, 1, i + 1, label=f's{i + 1}')
            ax.cla()
            axs.append(ax)
            yi = np.zeros_like(x, dtype=float)
            legends = []
            cnts, _ = np.histogram(X[:, i], bins=bins, density=True)
            ymax = cnts.max()
            for j, (cname, cdist) in enumerate(self.comps[i].items()):
                scdist = self.starting_pos.comps[i][cname]
                yj = self.weights[i][cname] * cdist.pdf(x)
                plt.plot(x, yj, c='C%d' % idmap[cname])
                legends.append(cname)
                str_params = f'{cname}: mu={cdist.mu:.2f}({scdist.mu:.2f}),' \
                             f' sigma={cdist.sigma:.2f}({scdist.sigma:.2f}),' \
                             f' alpha={cdist.alpha:.2f}({scdist.alpha:.2f})'
                plt.text(xmin, ymax * (1 - 0.1 * (j + 1)), str_params)
                yi += yj
            plt.plot(x, yi, c='C%d' % idmap['mixture'])
            legends.append(f'mixture{i + 1}')
            plt.hist(X[:, i], bins=bins, density=True, facecolor='w', ec='k', lw=1)
            plt.text(xmax * 0.76, ymax * 0.9, f'll = {slls[i]:.5f}')
            plt.legend(legends)
        ax = plt.subplot(nsub, 1, nsub, label=f'll_curve')
        ax.cla()
        plt.plot(lls[1:])
        if lls:
            plt.text(0, lls[-1] - 0.1, f'll = {lls[-1]:.5f}')
        ax = plt.gcf().axes[0]
        plt.axes(ax)
        if self.title:
            plt.title(self.title)
        plt.pause(0.5)

    def log_likelihood(self, X):
        return np.log(self.likelihood(X)).mean()

    def sep_log_likelihood(self, X):
        return np.log(self.likelihood(X)).mean(1)

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
