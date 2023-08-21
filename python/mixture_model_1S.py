import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import skew_normal
from constraints import *
from param_binary_search import *
from mixture_model_base import MixtureModelBase
from kanren import *
import sympy

from collections import defaultdict
from typing import *

SN = skew_normal.SkewNormal


def eval_eqs(eqs: List[str]):
    def eval_eq(eq: str):
        return sympy.parse_expr(eq)

    return list(map(eval_eq, eqs))


comp_id = {'C': 1, 'IC': 2, 'IC2': 3, 'I1': 4, 'I2': 5}


class MixtureModel1S(MixtureModelBase):
    def __init__(self, skew_dirs,
                 tolerance=1e-4,
                 binwidth=2.5,
                 plotstep=0.1,
                 show_plotting=False):
        """

        :param skew_dirs: skew directions
        :param ic2_comp: model IC2 component separately
        :param constraint_types: options [ weights, pdf, cdf ]
        :param tolerance: convergence threshold
        :param binwidth: binwidth for histogram plotting
        :param plotstep: step for the grid when plotting curves (e.g. pdf, cdf)
        :param show_plotting: plot each step during running of the program
        """
        self.tolerance = tolerance
        self.binwidth = binwidth
        self.plotstep = plotstep
        self.show_plotting = show_plotting
        self.join_comps = [
            ['C', 'IC', 'I1'],
        ]
        self.n_samples = len(self.join_comps)
        # self.comp_rels = {
        #     'C': 'IC',
        #     'IC': 'I1',
        # }

        self.comps = [{
            cname: SN(skew_dirs[cname]) for cname in sorted(set(comps), key=lambda x: comp_id[x])
        } for comps in self.join_comps]
        self.weights = [{
            cname: 1 for cname in sorted(set(comps), key=lambda x: comp_id[x])
        } for comps in self.join_comps]

        self.all_comps = {}
        for i in range(self.n_samples):
            for cname, cdist in self.comps[i].items():
                self.all_comps[cname] = cdist

    def pred(self, X):
        X = np.array(X)
        n, d = X.shape
        assert d == 2
        pj = list(map(lambda c: np.zeros((n, len(c))), self.comps))
        p = [np.zeros((n, len(self.comps[i]))) for i in range(self.n_samples)]
        for i in range(len(self.comps)):
            ws = self.weights[i]
            for j, (cname, cdist) in enumerate(self.comps[i].items()):
                pj[i][:, j] = ws[cname] * cdist.pdf(X[:, i])
            p0 = np.sum(pj[i], 1).reshape((-1, 1))
            p[i] = pj[i] / p0
        return p

    # def get_fitting_params
    def get_constraint(self, k):
        # if k == 'C_IC':
        pass

    def solve_weights(self, sumRs, n):
        ws = {}
        res = []
        for i in range(self.n_samples):
            res.append({})
            for j, (cname, cdist) in enumerate(self.comps[i].items()):
                # self.weights[i][cname] = sumRs[i][j] / n
                res[i][cname] = sumRs[i][j] / n
                ws[f'w{i + 1}{cname}'] = sumRs[i][j] / n

        return res

    def fit(self, X):
        prev_ll = -np.inf
        X = X[X[:, 1] != 0]

        n, d = X.shape

        xmax = np.max(X)
        xmin = np.min(X)
        self.n_xsamples = 200
        xstep = (xmax - xmin) / self.n_xsamples
        # x = np.arange(xmin, xmax + self.plotstep, self.plotstep)
        x = np.arange(xmin, xmax + xstep, xstep)
        sigma = np.sqrt(X[:, 0].var())

        comp_constraints = {}

        for i in range(len(self.comps)):
            for j, (cname, _) in enumerate(self.comps[i].items()):
                # mu = xmax - (xmax - xmin) / len(self.comps[i]) * (1 / 2 + j)
                mu = X[:, 0].mean() + 1 * sigma - j * 1 * sigma
                # self.comps[i][cname] = SN()
                self.weights[i][cname] = 1 / len(self.comps[i])
                self.comps[i][cname].mu = mu
                self.comps[i][cname].sigma = sigma / 3

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

            # sum_rs = []
            new_weights = self.solve_weights([np.sum(r, 0) for r in rs], n)
            self.weights = new_weights

            d = defaultdict(lambda: [None, [], []])
            for i, r in enumerate(rs):
                for j, (cname, cdist) in enumerate(self.comps[i].items()):
                    # self.weights[i][cname] = sum_rs[i][j] / n
                    d[cname][0] = cdist
                    d[cname][1].append(X[:, i])
                    d[cname][2].append(r[:, j])
                    # self.comps[i][cname].iterfit(X[:, i], r[:, j])

            for cname, (cdist, xs, rs) in d.items():

                # new_dist = cdist.iterfit(np.array(xs), np.array(rs), cons)
                cdist.iterfit(np.array(xs), np.array(rs))
                # self.all_comps[cname] = new_dist
                # for i in range(len(self.comps)):
                #     if cname in self.comps[i]:
                #         self.comps[i][cname] = new_dist

            ll = self.log_likelihood(X)
            lls.append(ll)
            slls = self.sep_log_likelihood(X)
            if self.show_plotting:
                self.plot(X, lls, slls)

        slls = self.sep_log_likelihood(X)
        self.plot(X, lls, slls)
        return ll, lls
