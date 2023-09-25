import threading

import matplotlib.pyplot as plt
# from scipy import stats
import numpy as np
import skew_normal as skew_normal
# import skew_normal_cupy as skew_normal
import normal_gpu
from constraints import *
from param_binary_search import *
from mixture_model_1S import MixtureModel1S
from mixture_model_base import MixtureModelBase
from myutils import *

from kanren import *
import sympy
import time
from itertools import starmap

from collections import defaultdict
from typing import *

SN = skew_normal.SkewNormal
N = normal_gpu.Normal


def eval_eqs(eqs: List[str]):
    def eval_eq(eq: str):
        return sympy.parse_expr(eq)

    return list(map(eval_eq, eqs))


comp_id = {'C': 1, 'IC': 2, 'IC2': 3, 'I1': 4, 'I2': 5}


class MixtureModel(MixtureModelBase):
    def __init__(self, skew_dirs,
                 ic2_comp=False,
                 constraints=('weights',),
                 tolerance=1e-8,
                 binwidth=2.5,
                 plotstep=0.1,
                 show_plotting=False,
                 plot_interval=1,
                 initialized=False,
                 # title=None,
                 event_notify_func=None,
                 init_strategy=None,
                 **kwargs):
        """

        :param skew_dirs: skew directions
        :param ic2_comp: model IC2 component separately
        :param constraints: options [ weights, mode, pdf, cdf, weighted_pdf ]
        :param tolerance: convergence threshold
        :param binwidth: binwidth for histogram plotting
        :param plotstep: step for the grid when plotting curves (e.g. pdf, cdf)
        :param show_plotting: plot each step during running of the program
        :param title: title for the figure
        """
        super().__init__(**kwargs)
        self.skew_dirs = skew_dirs
        self.tolerance = tolerance
        self.binwidth = binwidth
        self.plotstep = plotstep
        self.show_plotting = show_plotting
        self.plot_interval = plot_interval
        self.initialized = initialized
        self.event_notify_func = event_notify_func
        self.ic2_comp = ic2_comp
        if ic2_comp:
            self.join_comps = [
                ['C', 'C', 'IC', 'IC', 'IC', 'I1', 'I1', 'I1'],
                ['IC', 'I1', 'C', 'IC2', 'I1', 'C', 'IC', 'I2'],
            ]
        else:
            self.join_comps = [
                ['C', 'C', 'IC', 'IC', 'IC', 'I1', 'I1', 'I1'],
                ['IC', 'I1', 'C', 'IC', 'I1', 'C', 'IC', 'I1'],
            ]
        self.n_samples = len(self.join_comps)
        self.constraints = constraints
        self.init_strategy = init_strategy
        self.kwargs = kwargs
        # constraint_types = constraints
        # self.comp_rels = {
        #     'C': 'IC',
        #     'IC': 'I1',
        # }

        self.all_comps = {}
        for i in range(self.n_samples):
            for cname in self.join_comps[i]:
                if cname not in self.all_comps:
                    self.all_comps[cname] = SN(skew_dirs[cname])

        # self.comps = [{
        #     cname: SN(skew_dirs[cname]) for cname in sorted(set(comps), key=lambda x: comp_id[x])
        # } for comps in self.join_comps]
        self.comps = [{
            cname: self.all_comps[cname] for cname in sorted(set(comps), key=lambda x: comp_id[x])
        } for comps in self.join_comps]
        self.weights = [{
            cname: 1 for cname in sorted(set(comps), key=lambda x: comp_id[x])
        } for comps in self.join_comps]
        # self.weights = [
        #     NamedArray(sorted(set(comps), key=lambda x: comp_id[x]), 1)
        #     for comps in self.join_comps
        # ]

        self.weight_constrained = 'weights' in constraints
        self.mode_constrained = 'mode' in constraints
        self.pdf_constrained = 'pdf' in constraints
        self.cdf_constrained = 'cdf' in constraints
        self.weighted_pdf_constrained = 'weighted_pdf' in constraints
        self.starting_pos = None
        self.lls = []
        self.ll = None
        self.stopped = False

        # self.comp_constraints = {
        #     'IC':  (self.all_comps['C'], self.all_comps['IC']),
        #     'IC2': (self.all_comps['IC'], self.all_comps['IC2']),
        #     'I1':  (self.all_comps['IC'], self.all_comps['I1']),
        #     'I2':  (self.all_comps['I1'], self.all_comps['I2']),
        # }
        # self.comp_constraints = defaultdict(lambda: {'pdf': False, 'cdf': False})
        # for cons in constraints:
        #     if cons == 'weights':
        #         continue
        #     ct, c1, rel, c2 = cons.split()
        #     if rel == '<':
        #         c1, c2 = c2, c1
        #         rel = '>'
        #         self.comp_constraints[c2]['right_comp'] = c1
        #         self.comp_constraints[c2][ct] = True

        # self.comps
        self.sym_lambdas = [sympy.Symbol(f'lambda_{i + 1}', real=True) for i in range(self.n_samples)]
        self.sym_etas = []
        self.sym_Rs = []
        self.sym_ws = []
        self.syms = ['lambda1', 'lambda2']
        for i in range(self.n_samples):
            # sym_etas = {}
            # sym_Rs = {}
            # sym_ws = {}
            for j, cname in enumerate(self.comps[i].keys()):
                # sym_etas[cname] = sympy.Symbol(f'eta{i + 1}{cname}', real=True)
                # sym_Rs[cname] = sympy.Symbol(f'R{i + 1}{cname}', real=True)
                # sym_ws[cname] = sympy.Symbol(f'w{i + 1}{cname}', real=True)
                spart = f'{i + 1}{cname}'
                self.syms += list(map(lambda v: v + spart, ['eta', 'R', 'w']))
            # self.sym_etas.append(sym_etas)
            # self.sym_Rs.append(sym_Rs)
            # self.sym_ws.append(sym_ws)
        self.syms = {name: sympy.Symbol(name, real=True) for name in self.syms}

        for i in range(self.n_samples):
            sym_ws = {}
            for j, cname in enumerate(self.comps[i].keys()):
                spart = f'{i + 1}{cname}'
                sym_ws[cname] = self.syms[f'w{spart}']
            # print(self.sym_ws)
            self.sym_ws.append(sym_ws)

        self.Lag_sys_eqs = []
        ind = {}
        cnames = []
        for i in range(self.n_samples):
            for cname in self.comps[i].keys():
                ind[(i, cname)] = len(ind)
                cnames.append(cname)
        g = np.zeros([len(ind), len(ind)])
        eqs = [0] * len(ind)
        for i, cname in ind.keys():
            # eqs[ind[(i, cname)]] = sympy.parse_expr(f'R{i + 1}{cname} / w{i + 1}{cname} + lambda{i + 1}')
            eqs[ind[(i, cname)]] = (self.syms[f'R{i + 1}{cname}'] / self.syms[f'w{i + 1}{cname}']
                                    + self.syms[f'lambda{i + 1}'])
        for cs in zip(*self.join_comps):
            t = tuple(map(lambda c: ind[c], enumerate(cs)))
            g[t] = 1
            cname = cnames[t[1]]
            eqs[t[0]] += self.syms[f'eta{2}{cname}']
        for cname in self.comps[1]:
            eqs[ind[(1, cname)]] -= self.syms[f'eta{2}{cname}']

        # ind1 = {}
        # for cname in self.comps[0].keys():
        #     ind1[cname] = len(ind1)
        # weight_cons1 = [self.syms[f'w1{cname}'] for cname in self.comps[0].keys()]
        # for cs in zip(*self.join_comps):
        #     i = ind1[cs[0]]
        #     weight_cons1[i] -= self.syms[f'w2{cs[1]}']
        #
        # ind2 = {}
        # for cname in self.comps[1].keys():
        #     ind2[cname] = len(ind2)
        # eqs2 = [self.syms[f'w2{cname}'] for cname in self.comps[1].keys()]
        # weight_cons2 = [self.syms[f'w2{cname}'] for cname in self.comps[1].keys()]
        # for cs in zip(*self.join_comps):
        #     i = ind2[cs[1]]
        #     eqs2[i] -= self.syms[f'w1{cs[0]}']
        #     weight_cons2[i] -= self.syms[f'w1{cs[0]}']
        # self.weight_cons = weight_cons1 + weight_cons2
        self.weight_cons = {}
        for i in range(self.n_samples):
            for cs in zip(*self.join_comps):
                self.weight_cons[f'{i + 1}{cs[i]}'] = self.syms[f'w{i + 1}{cs[i]}']
        for i in range(self.n_samples):
            for cs in zip(*self.join_comps):
                self.weight_cons[f'{i + 1}{cs[i]}'] -= self.syms[f'w{2 - i}{cs[1 - i]}']

        # for cname in self.comps[1].keys():
        #     i = ind2[cname]
        #     # eqs2[i] = sympy.Or(sympy.Eq(eqs2[i], 0), sympy.Eq(self.syms[f'eta2{cname}'], 0))
        #     eqs2[i] *= self.syms[f'eta2{cname}']
        #     # self.weight_cons[i] = self.weight_cons[i] <= 0
        eqs2 = []
        for cname in self.comps[1].keys():
            eqs2.append(self.syms[f'eta2{cname}'] * self.weight_cons[f'2{cname}'])

        # for cname in self.comps[1].keys():
        #     weight_cons.append(self.syms[f'eta{2}{cname}'] <= 0)

        eqs3 = []
        for i in range(self.n_samples):
            eqs3.append(sum(self.sym_ws[i].values()) - 1)
            # weight_cons.append(self.syms[f'lambda{i+1}'] <= 0)

        ineqs = []
        for i in range(self.n_samples):
            for cname in self.comps[i].keys():
                ineqs.append(self.syms[f'w{i + 1}{cname}'] >= 0)
        self.pos_ws = ineqs
        self.Lag_sys_eqs = eqs + eqs3  # + ineqs + eqs2 + weight_cons

        pass

    def stop(self):
        self.stopped = True

    def pred(self, X):
        X = np.array(X)
        n, d = X.shape
        assert d == 2
        pj = list(map(lambda c: np.zeros((n, len(c))), self.comps))
        p = list(map(lambda c: np.zeros((n, len(c))), self.comps))
        # p = [np.zeros((n, 3)), np.zeros((n, 4))]
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
                ws[self.syms[f'w{i + 1}{cname}']] = sumRs[i][j] / n

        if not self.weight_constrained:
            return res

        sys_eqs = []
        vals = {}
        for i in range(self.n_samples):
            for j, cname in enumerate(self.comps[i].keys()):
                # sys_eqs.append(sympy.Eq(self.syms[f'R{i + 1}{cname}'], sumRs[i][j]))
                vals[self.syms[f'R{i + 1}{cname}']] = sumRs[i][j]

        flag = False
        active_cons = []
        for j, cname in enumerate(self.comps[1].keys()):
            # if self.weight_cons[j].subs(ws) <= 0:
            if self.weight_cons[f'2{cname}'].subs(ws) <= 0:
                # sys_eqs.append(sympy.Eq(self.syms[f'eta{2}{cname}'], 0))
                # sys_eqs.append(self.syms[f'eta{2}{cname}'])
                vals[self.syms[f'eta{2}{cname}']] = 0
            else:
                flag = True
                active_cons.append(cname)
                # active_cons += 1
                # sys_eqs.append(self.weight_cons[j])
                # sys_eqs.append(self.weight_cons[f'2{cname}'])
        if not flag:
            return res
        if len(active_cons) > 1:
            print(f'got {active_cons} active constraints, might unable to find solution')
            return res

        def all_sys_eqs():
            for i in range(len(active_cons)):
                newsys, newvals = deepcopy(sys_eqs), deepcopy(vals)
                for j in range(len(active_cons)):
                    cname = active_cons[j]
                    if i == j:
                        newsys.append(self.weight_cons[f'2{cname}'])
                    else:
                        newvals[self.syms[f'eta{2}{cname}']] = 0
                yield newsys, newvals

        def solve(sys_eqs, vals):
            sys_eqs += self.Lag_sys_eqs
            sys_eqs = list(map(lambda x: x.subs(vals), sys_eqs))
            solutions = sympy.solve(sys_eqs)
            solutions = [s for s in solutions if all(map(lambda x: x.subs(s), self.pos_ws))]
            if len(solutions) < 1:
                return
            solution = solutions[0]
            return solution

        flag = False
        for solution in starmap(solve, all_sys_eqs()):
            if not solution:
                continue
            flag = True
            break

        # assert flag
        if not flag:
            return res


        # if len(solutions) != 1:
        #     print(f'failed to solve weights, got {len(solutions)} solutions, return unconstrained weights')
        #     return res

        # print(solution)
        for i in range(self.n_samples):
            for j, cname in enumerate(self.comps[i].keys()):
                # self.weights[i][cname] = float(solution[self.syms[f'w{i + 1}{cname}']])
                res[i][cname] = float(solution[self.syms[f'w{i + 1}{cname}']].evalf(chop=True))
                ws[self.syms[f'w{i + 1}{cname}']] = float(solution[self.syms[f'w{i + 1}{cname}']].evalf(chop=True))

        # print(self.weights)
        for j, cname in enumerate(self.comps[1].keys()):
            r = float(self.weight_cons[f'2{cname}'].subs(ws))
            assert r < 0 or np.isclose(r, 0)

        if len(active_cons) > 1:
            print('Issue solved')

        return res

    def putdata(self, X):
        X = X[X[:, 1] != 0].astype(np.float32)
        self.X = X
        xmax = np.max(X)
        xmin = -100
        self.xplot = np.arange(xmin, xmax + self.plotstep, self.plotstep)
        bins = np.arange(xmin, xmax + self.binwidth, self.binwidth)
        self.hist = [
            np.histogram(X[:, 0], bins, density=True),
            np.histogram(X[:, 1], bins, density=True),
        ]

    def init_model(self, X):
        self.init_range(X)

        sigma = np.sqrt(X[:, 0].var())
        xmax = X.max()
        mu = xmax
        if self.init_strategy == 'random':
            np.random.seed(int(time.time()))
            for i in range(len(self.comps)):
                for j, (cname, _) in enumerate(self.comps[i].items()):
                    mu_offset = np.random.uniform(0, 1)
                    print(f'{cname} mu_offset {mu_offset}')
                    sigma_scale = np.random.random()
                    mu -= j * mu_offset * sigma
                    self.weights[i][cname] = np.float32(1 / len(self.comps[i]))
                    self.comps[i][cname].mu = np.float32(mu)
                    self.comps[i][cname].sigma = np.float32(sigma)
                    self.comps[i][cname].calc_alt_params()
            self.weights[1]['C'] *= np.float32(0.05)
        # elif self.init_strategy == 'one_sample':
        #     model1s = MixtureModel1S(self.skew_dirs, self.constraints, self.tolerance, self.binwidth, self.plotstep,
        #                              False, self.plot_interval, **self.kwargs)
        #     model1s.fit(X)
        else:
            for i in range(len(self.comps)):
                for j, (cname, _) in enumerate(self.comps[i].items()):
                    mu = X[:, 0].mean() + 0.5 * sigma - j * 0.5 * sigma
                    self.weights[i][cname] = np.float32(1 / len(self.comps[i]))
                    self.comps[i][cname].mu = np.float32(mu)
                    self.comps[i][cname].sigma = np.float32(sigma)
                    self.comps[i][cname].calc_alt_params()
            self.weights[1]['C'] *= np.float32(0.05)

        self.ll = self.log_likelihood(X)
        self.lls = [self.ll]
        self.slls = self.sep_log_likelihood(X)

        self.create_constraints()

        self.starting_pos = self.frozen()
        self.initialized = True

        if self.show_plotting:
            self.plot(X, self.lls, self.slls)
        if self.event_notify_func:
            self.event_notify_func('update', self.starting_pos)

    def create_constraints(self):
        x = self.xrange
        relative_constraints = {}
        self.relative_constraints = relative_constraints

        relative_constraints['C_IC'] = RelativeConstraint(self.all_comps['C'], self.all_comps['IC'],
                                                          x_range=x, mode=self.mode_constrained,
                                                          pdf=self.pdf_constrained,
                                                          cdf=self.cdf_constrained)
        relative_constraints['IC_C_w2'] = RelativeConstraint(self.all_comps['IC'], self.all_comps['C'],
                                                             weights=(self.weights[1]['IC'],
                                                                      self.weights[1]['C']),
                                                             x_range=x, pdf=self.weighted_pdf_constrained,
                                                             cdf=False, mode=False)
        relative_constraints['IC_I1'] = RelativeConstraint(self.all_comps['IC'], self.all_comps['I1'],
                                                           x_range=x, mode=self.mode_constrained,
                                                           pdf=self.pdf_constrained,
                                                           cdf=self.cdf_constrained)
        if self.ic2_comp:
            relative_constraints['IC_IC2'] = RelativeConstraint(self.all_comps['IC'], self.all_comps['IC2'],
                                                                x_range=x, mode=self.mode_constrained,
                                                                pdf=self.pdf_constrained,
                                                                cdf=self.cdf_constrained)
            relative_constraints['IC_IC2_w2'] = RelativeConstraint(self.all_comps['IC'], self.all_comps['IC2'],
                                                                   weights=(self.weights[1]['IC'],
                                                                            self.weights[1]['IC2']),
                                                                   x_range=x, pdf=self.weighted_pdf_constrained,
                                                                   cdf=False, mode=False)
            relative_constraints['IC2_I1'] = RelativeConstraint(self.all_comps['IC2'], self.all_comps['I1'],
                                                                x_range=x, mode=self.mode_constrained,
                                                                pdf=self.pdf_constrained,
                                                                cdf=self.cdf_constrained)
            relative_constraints['I1_I2'] = RelativeConstraint(self.all_comps['I1'], self.all_comps['I2'],
                                                               x_range=x, mode=self.mode_constrained,
                                                               pdf=self.pdf_constrained,
                                                               cdf=self.cdf_constrained)
        # cons = None
        comp_constraints = {}
        self.comp_constraints = comp_constraints
        # if cname == 'C':
        comp_constraints['C'] = ComposedChecker(
            relative_constraints['C_IC'].getDistChecker('right'),
            relative_constraints['IC_C_w2'].getDistChecker('left'),
        )
        # elif cname == 'IC':
        if self.ic2_comp:
            comp_constraints['IC'] = ComposedChecker(
                relative_constraints['C_IC'].getDistChecker('left'),
                relative_constraints['IC_C_w2'].getDistChecker('right'),
                relative_constraints['IC_IC2'].getDistChecker('right'),
                relative_constraints['IC_IC2_w2'].getDistChecker('right'),
                relative_constraints['IC_I1'].getDistChecker('right'),
            )
        else:
            comp_constraints['IC'] = ComposedChecker(
                relative_constraints['C_IC'].getDistChecker('left'),
                relative_constraints['IC_C_w2'].getDistChecker('right'),
                relative_constraints['IC_I1'].getDistChecker('right'),
            )
        # elif cname == 'IC2':
        if self.ic2_comp:
            comp_constraints['IC2'] = ComposedChecker(
                relative_constraints['IC_IC2'].getDistChecker('left'),
                relative_constraints['IC_IC2_w2'].getDistChecker('left'),
                relative_constraints['IC2_I1'].getDistChecker('right'),
            )
        # elif cname == 'I1':
        if self.ic2_comp:
            comp_constraints['I1'] = ComposedChecker(
                relative_constraints['IC_I1'].getDistChecker('left'),
                relative_constraints['IC2_I1'].getDistChecker('left'),
                relative_constraints['I1_I2'].getDistChecker('right'),
            )
        else:
            comp_constraints['I1'] = ComposedChecker(
                relative_constraints['IC_I1'].getDistChecker('left'),
            )
        # elif cname == 'I2':
        if self.ic2_comp:
            comp_constraints['I2'] = relative_constraints['I1_I2'].getDistChecker('left')

    def check_constraints(self):
        for c, dist in self.all_comps.items():
            comp = f'{c}'
            cons = self.comp_constraints[comp]
            if not cons(dist):
                return False
        return True

    def fit(self, X):
        prev_ll = -np.inf
        # X = X[X[:, 1] != 0].astype(np.float32)
        # xmax = np.max(X)
        # xmin = -100
        # xplot = np.arange(xmin, xmax + self.plotstep, self.plotstep)
        # bins = np.arange(xmin, xmax + self.binwidth, self.binwidth)
        # self.hist = [
        #     np.histogram(X[:, 0], bins, density=True),
        #     np.histogram(X[:, 1], bins, density=True),
        # ]

        n, d = X.shape

        print(f'model initialized {self.initialized}')
        if not self.initialized:
            self.init_model(X)

        # self.starting_pos = deepcopy(self)
        self.starting_pos = self.frozen()
        # self.putdata(X)

        self.ll = self.log_likelihood(X)
        self.lls = [self.ll]
        self.slls = self.sep_log_likelihood(X)
        if self.show_plotting:
            self.plot(X, self.lls, self.slls)
            pass
        prev_t = time.time()
        while abs(self.ll - prev_ll) > self.tolerance and not self.stopped:
            # print(f'iteration {len(self.lls)}')
            prev_ll = self.ll
            rs = self.pred(X)

            # print('ll', self.ll)
            # sum_rs = []
            new_weights = self.solve_weights([np.sum(r, 0) for r in rs], n)

            # print(f'solved weights')

            def proj_weights():
                cons = self.relative_constraints['IC_C_w2'].getWeightChecker('left')
                old_w2c = self.weights[1]['C']

                w2c = new_weights[1]['C']
                w2ic = new_weights[1]['IC']
                old_sum = w2c + w2ic
                new_w2c = param_binary_search(old_w2c, w2c, cons)
                rw = new_w2c / (new_w2c + w2ic)
                restw = w2c - new_w2c
                w2c = new_w2c + rw * restw
                w2ic += (1 - rw) * restw
                new_sum = w2c + w2ic
                w2ic += old_sum - new_sum

                new_weights[1]['C'] = w2c
                new_weights[1]['IC'] = w2ic

            proj_weights()
            self.weights = new_weights
            # print('projected weights', self.weights)

            d = defaultdict(lambda: [None, [], []])
            for i, r in enumerate(rs):
                for j, (cname, cdist) in enumerate(self.comps[i].items()):
                    # self.weights[i][cname] = sum_rs[i][j] / n
                    d[cname][0] = cdist
                    d[cname][1].append(X[:, i])
                    d[cname][2].append(r[:, j])
                    # self.comps[i][cname].iterfit(X[:, i], r[:, j])

            for cname, (cdist, xs, rs) in d.items():
                cdist.iterfit(np.array(xs), np.array(rs), self.comp_constraints[cname])

            # ll = self.log_likelihood(X)
            self.ll = self.log_likelihood(X)
            # assert ll >= self.lls[-1]
            if self.ll < self.lls[-1]:
                print('ll decreased', f'{self.lls[-1]} -> {self.ll}')
            self.lls.append(self.ll)
            self.slls = self.sep_log_likelihood(X)

            isplotting = False

            # def update_fig():
            #     nonlocal isplotting
            #     if isplotting:
            #         return
            #     isplotting = True
            #     self.plot(X, self.lls, self.slls)
            #     isplotting = False

            meter = TimeMeter()
            cur_t = time.time()
            # print(self.show_plotting, self.plot_interval, cur_t - prev_t)
            if self.event_notify_func:
                self.event_notify_func('update', self.frozen())
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
        if self.stopped and self.event_notify_func:
            self.event_notify_func('stopped')
        return self.ll, self.lls
