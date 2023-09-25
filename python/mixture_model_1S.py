import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import skew_normal
from constraints import *
from param_binary_search import *
from mixture_model_base import MixtureModelBase
from kanren import *
import sympy
from myutils import *

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
                 constraints=(),
                 tolerance=1e-4,
                 binwidth=2.5,
                 plotstep=0.1,
                 show_plotting=False,
                 plot_interval=1,
                 **kwargs):
        """

        :param skew_dirs: skew directions
        :param ic2_comp: model IC2 component separately
        :param constraint_types: options [ weights, pdf, cdf ]
        :param tolerance: convergence threshold
        :param binwidth: binwidth for histogram plotting
        :param plotstep: step for the grid when plotting curves (e.g. pdf, cdf)
        :param show_plotting: plot each step during running of the program
        """
        super().__init__(**kwargs)
        self.tolerance = tolerance
        self.binwidth = binwidth
        self.plotstep = plotstep
        self.show_plotting = show_plotting
        self.plot_interval = plot_interval
        self.join_comps = [
            ['C', 'IC', 'I1'],
        ]
        self.n_samples = len(self.join_comps)
        self.constraints = constraints
        self.starting_pos = None
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

        self.weight_constrained = 'weights' in constraints
        self.mode_constrained = 'mode' in constraints
        self.pdf_constrained = 'pdf' in constraints
        self.cdf_constrained = 'cdf' in constraints
        self.weighted_pdf_constrained = 'weighted_pdf' in constraints

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
            truncate_zero(p0)
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

    def create_constraints(self):
        x = self.xrange
        relative_constraints = {}
        self.relative_constraints = relative_constraints

        relative_constraints['C_IC'] = RelativeConstraint(self.all_comps['C'], self.all_comps['IC'],
                                                          x_range=x, mode=self.mode_constrained,
                                                          pdf=self.pdf_constrained,
                                                          cdf=self.cdf_constrained)
        relative_constraints['IC_I1'] = RelativeConstraint(self.all_comps['IC'], self.all_comps['I1'],
                                                           x_range=x, mode=self.mode_constrained,
                                                           pdf=self.pdf_constrained,
                                                           cdf=self.cdf_constrained)
        comp_constraints = {}
        self.comp_constraints = comp_constraints
        # if cname == 'C':
        comp_constraints['C'] = ComposedChecker(
            relative_constraints['C_IC'].getDistChecker('right'),
        )
        # elif cname == 'IC':
        comp_constraints['IC'] = ComposedChecker(
            relative_constraints['C_IC'].getDistChecker('left'),
            relative_constraints['IC_I1'].getDistChecker('right'),
        )
        # elif cname == 'I1':
        comp_constraints['I1'] = ComposedChecker(
            relative_constraints['IC_I1'].getDistChecker('left'),
        )


    def fit(self, X):
        prev_ll = -np.inf
        X = X[X[:, 1] != 0]

        n, d = X.shape

        self.init_range(X)

        xmax = np.max(X)
        xmin = np.min(X)
        self.n_xsamples = 200
        xstep = (xmax - xmin) / self.n_xsamples
        # x = np.arange(xmin, xmax + self.plotstep, self.plotstep)
        x = np.arange(xmin, xmax + xstep, xstep)
        sigma = np.sqrt(X[:, 0].var())

        for i in range(len(self.comps)):
            for j, (cname, _) in enumerate(self.comps[i].items()):
                # mu = xmax - (xmax - xmin) / len(self.comps[i]) * (1 / 2 + j)
                mu = X[:, 0].mean() + 1 * sigma - j * 1 * sigma
                # self.comps[i][cname] = SN()
                self.weights[i][cname] = 1 / len(self.comps[i])
                self.comps[i][cname].mu = mu
                self.comps[i][cname].sigma = sigma / 3
                self.comps[i][cname].calc_alt_params()

        self.create_constraints()
        self.starting_pos = self.frozen()

        self.ll = self.log_likelihood(X)
        self.lls = [self.ll]
        self.slls = self.sep_log_likelihood(X)
        # plt.figure(figsize=[18, 9])
        if self.show_plotting:
            # plt.ion()
            # plt.show()
            self.plot(X, self.lls, self.slls)
        prev_t = time.time()
        while abs(self.ll - prev_ll) > self.tolerance:
            prev_ll = self.ll
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
                # cdist.iterfit(np.array(xs), np.array(rs))
                cdist.iterfit(np.array(xs), np.array(rs), self.comp_constraints[cname])
                # self.all_comps[cname] = new_dist
                # for i in range(len(self.comps)):
                #     if cname in self.comps[i]:
                #         self.comps[i][cname] = new_dist

            self.ll = self.log_likelihood(X)
            self.lls.append(self.ll)
            self.slls = self.sep_log_likelihood(X)
            meter = TimeMeter()
            cur_t = time.time()
            if self.show_plotting and cur_t - prev_t >= self.plot_interval:
                # thread = threading.Thread(target=lambda : self.plot(X, self.lls, self.slls))
                # thread = threading.Thread(target=update_fig)
                # thread.start()
                print('plotting...')
                self.plot(X, self.lls, self.slls)
                print('|', meter.read())
                print('plot finished')

                prev_t = time.time()

        self.slls = self.sep_log_likelihood(X)
        self.plot(X, self.lls, self.slls)
        return self.ll, self.lls
