import matplotlib.pyplot as plt
import scipy
import scipy.integrate
from scipy import stats
import numpy as np
import normal
from constraints import *
from param_binary_search import *
from plotting import *
from kanren import *
import sympy

from collections import defaultdict
from typing import *

Norm = normal.Normal
inf = np.inf

factors = []


class OrderStat2c:
    def __init__(self, dist1, dist2):
        self.dist1 = dist1
        self.dist2 = dist2
        self.dist1_higher = self._dist1_higher()
        self.dist2_higher = 1 - self.dist1_higher
        # self.mt11 = self._mt11()
        # self.mt21 = self._mt21()
        # self.mt12 = self._mt12()
        # self.mt22 = self._mt22()

    def __repr__(self):
        return str(self.__dict__)

    @classmethod
    def from_complementary(cls, other: 'OrderStat2c'):
        new = copy(other)
        new.dist1, new.dist2 = other.dist2, other.dist1
        new.dist1_higher, new.dist2_higher = other.dist2_higher, other.dist1_higher
        # new.mt11, new.mt12 = other.mt21, other.mt22
        # new.mt21, new.mt22 = other.mt11, other.mt12
        return new

    def _dist1_higher(self):
        unnormalized = lambda x: self.dist1.pdf(x) * self.dist2.cdf(x)
        factor_1 = scipy.integrate.quad(unnormalized, -inf, self.dist1.mu)[0]
        factor_2 = scipy.integrate.quad(unnormalized, self.dist1.mu, inf)[0]
        return factor_1 + factor_2

    def _dist2_higher(self):
        return 1 - self.dist1_higher

    def _mt11(self):
        f = lambda x: x * self.dist1.pdf(x) * self.dist2.cdf(x) / self.dist1_higher
        return scipy.integrate.quad(f, -inf, inf)[0] - self.dist1.mu

    def _mt21(self):
        f = lambda x: x * self.dist2.pdf(x) * (1 - self.dist1.cdf(x)) / self.dist1_higher
        return scipy.integrate.quad(f, -inf, inf)[0] - self.dist2.mu

    def _mt12(self):
        f = lambda x: (x - self.dist1.mu)**2 * self.dist1.pdf(x) * self.dist2.cdf(x) / self.dist1_higher
        return scipy.integrate.quad(f, -inf, inf)[0] - self.dist1.sigma**2

    def _mt22(self):
        f = lambda x: (x - self.dist2.mu)**2 * self.dist2.pdf(x) * (1 - self.dist1.cdf(x)) / self.dist1_higher
        return scipy.integrate.quad(f, -inf, inf)[0] - self.dist2.sigma**2

    def pdf(self, X):
        p1 = self.dist1.pdf(X[:, 0])
        p2 = self.dist2.pdf(X[:, 1])
        return p1 * p2 / self.dist1_higher

    def top_pdf1(self, x):
        return self.dist1.pdf(x) * self.dist2.cdf(x) / self.dist1_higher

    def bot_pdf2(self, x):
        return self.dist2.pdf(x) * (1 - self.dist1.cdf(x)) / self.dist1_higher

    # @staticmethod
    # def top_pdf1(dist1, dist2, X):
    #     unnormalized = lambda x: dist1.pdf(x) * dist2.cdf(x)
    #     factor = scipy.integrate.quad(unnormalized, -inf, inf)[0]
    #     factors.append(factor)
    #     return unnormalized(X) / factor
    #
    # @staticmethod
    # def bot_pdf2(dist1, dist2, X):
    #     unnormalized = lambda x: dist2.pdf(x) * (1 - dist1.cdf(x))
    #     factor = scipy.integrate.quad(unnormalized, -inf, inf)[0]
    #     return unnormalized(X) / factor


class OrderStatMixtureModel:
    def __init__(self,
                 tolerance=1e-4,
                 binwidth=2.5,
                 plotstep=0.1,
                 show_plotting=False,
                 **kwargs):
        # comps = {cname: Norm() for cname in ['C', 'IC', 'I']},
        self.comps = {cname: Norm() for cname in ['C', 'IC', 'I']}
        # self.join_comps = [
        #     ['C', 'C', 'IC', 'IC', 'I'],
        #     ['IC', 'I', 'IC', 'I', 'I'],
        # ]
        self.join_comps = [
            ['C', 'IC', 'C', 'I', 'IC', 'IC', 'IC', 'I', 'I', 'I'],
            ['IC', 'C', 'I', 'C', 'IC', 'IC', 'I', 'IC', 'I', 'I'],
        ]
        # self.join_comps = [
        #     ['C', 'C', 'IC', 'IC', 'IC', 'IC', 'I', 'I', 'I', 'I'],
        #     ['IC', 'I', 'C', 'IC', 'IC', 'I', 'C', 'IC', 'I', 'I'],
        # ]
        self.tolerance = tolerance
        self.binwidth = binwidth
        self.plotstep = plotstep
        self.show_plotting = show_plotting
        self.m = len(self.join_comps[0])
        self.weights = np.ones(self.m)
        self.doss = {}
        self.pcomp = {'C': 0.1, 'IC': 0.2}

    def likelihood(self, X):
        n = X.shape[0]
        # m = len(self.join_comps[0])
        m = self.m
        p = np.zeros((n, m))
        for j in range(m):
            cname1 = self.join_comps[0][j]
            cname2 = self.join_comps[1][j]
            # p1 = self.weights[j] * OrderStat2c.top_pdf1(self.comps[cname1], self.comps[cname2], X[:, 0])
            # p2 = self.weights[j] * OrderStat2c.bot_pdf2(self.comps[cname1], self.comps[cname2], X[:, 1])
            dos = self.doss[(cname1, cname2)]
            # dist1 = self.comps[cname1]
            # dist2 = self.comps[cname2]
            p[:, j] = self.weights[j] * dos.pdf(X)  # OrderStat2c(dist1, dist2).pdf(X)
        return p

    def log_likelihood(self, X):
        return np.log(self.likelihood(X).sum(1)).mean()

    def sep_log_likelihood(self, X):
        return np.log(self.likelihood(X).sum(1)).mean()

    def pred(self, X):
        X = np.array(X)
        n, d = X.shape
        assert d == 2
        m = self.m
        # p = np.zeros((n, m))
        pj = self.likelihood(X)
        p0 = np.sum(pj, 1).reshape((-1, 1))
        p = pj / p0
        return p

    def fit(self, X):
        X = X[X[:, 1] != 0]

        n, d = X.shape
        m = self.m

        x1 = X[:, 0]
        x2 = X[:, 1]

        sigma = np.sqrt(X[:, 0].var())
        for j, (cname, _) in enumerate(self.comps.items()):
            # mu = X[:, 0].mean() + 0.5 * sigma - j * 0.5 * sigma
            # self.comps[cname].mu = mu
            self.comps[cname].sigma = sigma / 2
        self.comps['C'].mu = 200
        self.comps['IC'].mu = 180
        self.comps['I'].mu = 150

        self.doss = {}
        for j in range(m):
            cname1 = self.join_comps[0][j]
            cname2 = self.join_comps[1][j]
            d1 = self.comps[cname1]
            d2 = self.comps[cname2]
            if (cname2, cname1) in self.doss:
                self.doss[(cname1, cname2)] = OrderStat2c.from_complementary(self.doss[(cname2, cname1)])
            else:
                self.doss[(cname1, cname2)] = OrderStat2c(d1, d2)

        for j in range(m):
            # self.weights[j] = 1/m
            cname1 = self.join_comps[0][j]
            cname2 = self.join_comps[1][j]
            d1 = self.comps[cname1]
            d2 = self.comps[cname2]
            dos = self.doss[(cname1, cname2)]

            p = 1
            if cname1 == 'C' or cname2 == 'C':
                p *= self.pcomp['C']
                if cname1 == 'IC' or cname2 == 'IC':
                    p *= self.pcomp['IC']
                else:
                    p *= 1 - self.pcomp['IC']
            else:
                p *= 1 - self.pcomp['C']
                if cname1 == 'IC':
                    p *= self.pcomp['IC']
                else:
                    p *= 1 - self.pcomp['IC']
                if cname2 == 'IC':
                    p *= self.pcomp['IC']
                else:
                    p *= 1 - self.pcomp['IC']
                if cname1 != cname2:
                    p *= 2
            p *= dos.dist1_higher
            self.weights[j] = p

        prev_ll = -np.inf
        ll = self.log_likelihood(X)
        lls = [ll]
        slls = self.sep_log_likelihood(X)
        plt.figure(figsize=[18, 9])
        if self.show_plotting:
            plt.ion()
            plt.show()
            self.plot(X, lls, slls)
        while abs(ll - prev_ll) > self.tolerance:
            prev_ll = ll
            rs = self.pred(X)
            # sumRs = np.sum(rs, 0)

            # for j in range(m):
            #     self.weights[j] = sumRs[j] / n

            d = defaultdict(lambda: np.zeros_like(X))
            d_weights = defaultdict(float)
            # d_mu = defaultdict(lambda: [0] * 2)
            # d_sigma = defaultdict(lambda: [0] * 2)
            for j in range(m):
                cname1 = self.join_comps[0][j]
                cname2 = self.join_comps[1][j]
                d1 = self.comps[cname1]
                d2 = self.comps[cname2]
                dos = self.doss[(cname1, cname2)]
                # pj = dos.pdf(X)
                r = rs[:, j]

                # plt.figure()
                # weighted_hist2d(X, r)
                # plt.title(f"{cname1, cname2}")

                if (cname2, cname1) in d_weights:
                    d_weights[(cname2, cname1)] += r.mean()
                else:
                    d_weights[(cname1, cname2)] += r.mean()

                d[cname1][:, 0] += r
                d[cname2][:, 1] += r

                # plt.figure()
                # weighted_hist(np.vstack([x1, x2]), np.vstack([d[cname1][:, 0], d[cname1][:, 1]]))
                # plt.title(cname1)
                #
                # plt.figure()
                # weighted_hist(np.vstack([x1, x2]), np.vstack([d[cname2][:, 0], d[cname2][:, 1]]))
                # plt.title(cname2)

                # prsum = np.sum(pr)
                # d_mu[cname1][0] += np.sum((x1 - dos.mt11) * pr)
                # d_mu[cname1][1] += prsum
                #
                # d_sigma[cname1][0] += np.sum(((x1 - d1.mu) ** 2 - dos.mt12) * pr)
                # d_sigma[cname1][1] += prsum
                # print(cname1, np.sqrt(d_sigma[cname1][0] / d_sigma[cname1][1]))
                #
                # d_mu[cname2][0] += np.sum((x2 - dos.mt21) * pr)
                # d_mu[cname2][1] += prsum
                #
                # d_sigma[cname2][0] += np.sum(((x2 - d2.mu) ** 2 - dos.mt22) * pr)
                # d_sigma[cname2][1] += prsum
                # print(cname2, np.sqrt(d_sigma[cname2][0] / d_sigma[cname2][1]))

            # for cname, cdist in self.comps.items():
            #     cdist.mu = d_mu[cname][0] / d_mu[cname][1]
            #     cdist.sigma = np.sqrt(d_sigma[cname][0] / d_sigma[cname][1])
            for cname, r in d.items():
                # plt.figure()
                # weighted_hist(np.vstack([x1, x2]), np.vstack([r[:, 0], r[:, 1]]))
                # plt.title(cname)
                self.comps[cname].iterfit(X, r)

            self.doss = {}
            for j in range(m):
                cname1 = self.join_comps[0][j]
                cname2 = self.join_comps[1][j]
                d1 = self.comps[cname1]
                d2 = self.comps[cname2]

                if (cname2, cname1) in self.doss:
                    self.doss[(cname1, cname2)] = OrderStat2c.from_complementary(self.doss[(cname2, cname1)])
                else:
                    self.doss[(cname1, cname2)] = OrderStat2c(d1, d2)
                dos = self.doss[(cname1, cname2)]

                if (cname1, cname2) in d_weights:
                    self.weights[j] = d_weights[(cname1, cname2)] * dos.dist1_higher
                else:
                    self.weights[j] = d_weights[(cname2, cname1)] * dos.dist1_higher

            ll = self.log_likelihood(X)
            lls.append(ll)
            slls = self.sep_log_likelihood(X)
            if self.show_plotting:
                self.plot(X, lls, slls)

        slls = self.sep_log_likelihood(X)
        self.plot(X, lls, slls)
        return ll, lls

    def plot(self, X, lls, slls):
        n, d = X.shape
        m = self.m
        xmax = np.max(X)
        x = np.arange(0, xmax + self.plotstep, self.plotstep)
        bins = np.arange(0, xmax + self.binwidth, self.binwidth)
        axs = []
        idmap = defaultdict(lambda: len(idmap))
        for i in range(d):
            ax = plt.subplot(3, 1, i + 1, label=f'{i}')
            ax.cla()
            axs.append(ax)
            yi = np.zeros_like(x, dtype=float)
            legends = []
            ys = defaultdict(lambda: np.zeros_like(x))
            for j in range(m):
                cname1 = self.join_comps[0][j]
                cname2 = self.join_comps[1][j]
                dos = self.doss[(cname1, cname2)]
                if i == 0:
                    ys[cname1] += self.weights[j] * dos.top_pdf1(x)
                elif i == 1:
                    ys[cname2] += self.weights[j] * dos.bot_pdf2(x)
            for cname, yj in ys.items():
                plt.plot(x, yj, c='C%d' % idmap[cname])
                yi += yj
                legends.append(cname)
            plt.plot(x, yi, c='C%d' % idmap['mixture'])
            legends.append(f'mixture{i + 1}')
            plt.hist(X[:, i], bins=bins, density=True, facecolor='w', ec='k', lw=1)
            plt.text(225, 0.015, f'll = {slls:.5f}')
            plt.legend(legends)
        ax = plt.subplot(3, 1, 3)
        ax.cla()
        plt.plot(lls[1:])
        plt.pause(0.001)
