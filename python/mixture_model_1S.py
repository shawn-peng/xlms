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
                 init_strategy=None,
                 mu_strategy=None,
                 seedoff=0,
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
        self.init_strategy = init_strategy
        self.mu_strategy = mu_strategy
        self.seedoff = seedoff
        # self.comp_rels = {
        #     'C': 'IC',
        #     'IC': 'I1',
        # }

        self.comps = [{
            cname: SN(skew_dirs[cname]) for cname in sorted(set(comps), key=lambda x: comp_id[x])
        } for comps in self.join_comps]
        self.weights = [{
            cname: DynamicParam(1) for cname in sorted(set(comps), key=lambda x: comp_id[x])
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
        self.initialized = False
        self.starting_pos = None
        self.lls = []
        self.ll = None
        self.stopped = False

    def log(self, *args):
        print(f'{self.title} id {self.seedoff}:', *args)

    def rand_sigmas(self, sigma, slow=0.5, shigh=1.0):
        sigmas = {}
        for cname, cdist in self.all_comps.items():
            sigma_scale = np.random.uniform(slow, shigh)
            # sigma_scale = 1.0
            sigmas[cname] = np.float64(sigma * sigma_scale)
        return sigmas

    def rand_alphas(self, alpha_scale):
        alphas = {}
        for cname, cdist in self.all_comps.items():
            # alpha_scale = np.random.uniform(slow, shigh)
            alpha_scale = np.random.uniform(1 / alpha_scale, alpha_scale)
            # alpha_scale = 1.0
            alphas[cname] = np.float64(1.0 * np.sign(cdist.alpha) * alpha_scale)
        return alphas

    def rand_mus_distance(self, xmax, sigma, scale=2.5):
        mus = {}
        mu = xmax
        j = 0
        for cname, cdist in self.all_comps.items():
            mu_offset = np.random.uniform(0, 1)
            mu -= scale * mu_offset * sigma
            mus[cname] = np.float64(mu)
            j += 1
        return mus

    def rand_mus_disturb(self, mu, sigma, scale=1):
        mus = {}
        mu = mu + 2 * sigma
        j = 0
        for cname, cdist in self.all_comps.items():
            mu_offset = np.random.uniform(0, 1)
            mus[cname] = np.float64(mu + scale * mu_offset * sigma)
            mu -= sigma
            j += 1
        return mus

    def rand_mus_split(self, xmax, xmin, mu, sigma, scale=1.5):
        mus = {}
        mus['IC'] = np.random.uniform(mu - scale * sigma, mu + scale * sigma)
        mus['C'] = np.random.uniform(mus['IC'], xmax)
        mus['I1'] = np.random.uniform(xmin, mus['IC'])
        return mus

    def mus_from_sample(self, sample):
        sample = np.sort(sample)[::-1]
        mus = {}
        j = 0
        for cname, cdist in self.all_comps.items():
            mus[cname] = np.float64(sample[j])
            j += 1
        return mus

    def rand_mus_uniform(self, xmax, xmin):
        sample = np.random.uniform(xmin, xmax, len(self.all_comps))
        mus = self.mus_from_sample(sample)
        return mus

    def rand_mus_gaussian(self, mu, sigma):
        sample = np.random.normal(mu, sigma, len(self.all_comps))
        mus = self.mus_from_sample(sample)
        return mus

    def init_model(self, X):
        X = X.astype(np.float64)
        self.log('start init ...')
        self.init_range(X)

        self.create_constraints()

        def plot():
            self.plot(X, [], self.sep_log_likelihood(X))

        sigma = np.sqrt(X[:, 0].var())
        mu = X[:, 0].mean()
        xmin = X.min()
        xmax = X.max()
        xmin = -50.0
        if self.init_strategy == 'random':
            # seed = int(time.time()) + self.seedoff
            seed = self.seedoff
            # seed = 31
            # seed = self.seedoff + 8
            # seed = 4
            print(f'seed {seed}')
            np.random.seed(seed)

            for i in range(len(self.comps)):
                for j, (cname, _) in enumerate(self.comps[i].items()):
                    self.weights[i][cname].set(np.float64(1 / len(self.comps[i])))

            frozen_model = self.frozen()
            # self.starting_pos = self.frozen()
            # plt.ion()
            while True:
                if self.mu_strategy == 'distance':
                    mus = self.rand_mus_distance(xmax, sigma)
                elif self.mu_strategy == 'disturb':
                    mus = self.rand_mus_disturb(mu, sigma)
                elif self.mu_strategy == 'split':
                    mus = self.rand_mus_split(xmax, xmin, mu, sigma)
                elif self.mu_strategy == 'uniform':
                    mus = self.rand_mus_uniform(xmax, 50)
                elif self.mu_strategy == 'gaussian':
                    mus = self.rand_mus_gaussian(mu, sigma)
                else:
                    mus = self.rand_mus_uniform(xmax, 50)
                sigmas = self.rand_sigmas(sigma, 0.25, 1.0)
                alphas = self.rand_alphas(frozen_model.all_comps['C'].alpha)
                for cname, cdist in self.all_comps.items():
                    cdist.mu = mus[cname]
                    cdist.sigma = sigmas[cname]
                    cdist.alpha = alphas[cname]
                # self.log(self.comps)
                # self.starting_pos = self.frozen()
                # self.plot(X, [], self.sep_log_likelihood(X))
                if self.check_constraints():
                    break
                self.log('resample params')
            self.log(self.comps)
        # elif self.init_strategy == 'one_sample':
        #     model1s = MixtureModel1S(self.skew_dirs, self.constraints, self.tolerance, self.binwidth, self.plotstep,
        #                              False, self.plot_interval, **self.kwargs)
        #     model1s.fit(X)
        else:
            for i in range(len(self.comps)):
                for j, (cname, _) in enumerate(self.comps[i].items()):
                    mu = X[:, 0].mean() + 0.5 * sigma - j * 0.5 * sigma
                    self.weights[i][cname].set(np.float64(1 / len(self.comps[i])))
                    self.comps[i][cname].mu = np.float64(mu)
                    self.comps[i][cname].sigma = np.float64(sigma)
                    self.comps[i][cname].calc_alt_params()

        self.ll = self.log_likelihood(X)
        self.lls = [self.ll]
        self.slls = self.sep_log_likelihood(X)

        self.starting_pos = self.frozen()
        self.initialized = True

        if self.show_plotting:
            self.plot(X, self.lls, self.slls)
        self.log('init finished ...')

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

    def check_constraints(self):
        for c, dist in self.all_comps.items():
            if not self.comp_constraints[c](dist, 2e-7):
                return False
        return True

    @property
    def cons_satisfied(self):
        return self.check_constraints()

    @property
    def fdr_thres(self):
        return self.fdr_thres_score(self.xrange, 0.01)

    def fdr_thres_score(self, x, fdr_thres):
        fdr = self.fdr(x)
        cond = fdr <= fdr_thres
        inds = np.argwhere(cond)
        if not len(inds):
            return np.inf
        return x[inds[0]][0]

    @property
    def fdr_curve(self):
        return self.fdr(self.plotting_xrange)

    def fdr(self, x):
        tp = self.weights[0]['C'] * (1 - self.all_comps['C'].cdf(x))
        fpic = self.weights[0]['IC'] * (1 - self.all_comps['IC'].cdf(x))
        fpi1 = self.weights[0]['I1'] * (1 - self.all_comps['I1'].cdf(x))
        fp = fpic + fpi1
        fdr = fp / (tp + fp)
        return fdr

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

    def update_weights(self, new_weights):
        for i in range(self.n_samples):
            for cname, w in self.weights[i].items():
                w.set(new_weights[i][cname])

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

        print(f'model initialized {self.initialized}')
        if not self.initialized:
            self.init_model(X)

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
            # self.weights = new_weights
            self.update_weights(new_weights)

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
