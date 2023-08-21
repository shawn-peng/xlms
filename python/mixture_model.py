import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import skew_normal
import normal_gpu
from constraints import *
from param_binary_search import *
from mixture_model_base import MixtureModelBase
from kanren import *
import sympy
import time

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
    def __init__(self, skew_dirs, ic2_comp=False,
                 constraints=('weights',),
                 tolerance=1e-4,
                 binwidth=2.5,
                 plotstep=0.1,
                 show_plotting=False,
                 plot_interval=1,
                 initialized=False,
                 title=None):
        """

        :param skew_dirs: skew directions
        :param ic2_comp: model IC2 component separately
        :param constraint_types: options [ weights, pdf, cdf ]
        :param tolerance: convergence threshold
        :param binwidth: binwidth for histogram plotting
        :param plotstep: step for the grid when plotting curves (e.g. pdf, cdf)
        :param show_plotting: plot each step during running of the program
        :param title: title for the figure
        """
        self.tolerance = tolerance
        self.binwidth = binwidth
        self.plotstep = plotstep
        self.show_plotting = show_plotting
        self.plot_interval = plot_interval
        self.initialized = initialized
        self.title = title
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
        constraint_types = constraints
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

        self.weight_constrained = 'weights' in constraint_types
        self.pdf_constrained = 'pdf' in constraint_types
        self.cdf_constrained = 'cdf' in constraint_types
        self.weighted_pdf_constrained = 'weighted_pdf' in constraint_types
        self.starting_pos = None
        self.lls = []
        self.ll = None

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
                self.weight_cons[f'{i+1}{cs[i]}'] = self.syms[f'w{i+1}{cs[i]}']
        for i in range(self.n_samples):
            for cs in zip(*self.join_comps):
                self.weight_cons[f'{i+1}{cs[i]}'] -= self.syms[f'w{2-i}{cs[1-i]}']

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
                ineqs.append(self.syms[f'w{i+1}{cname}'] >= 0)
        self.pos_ws = ineqs
        self.Lag_sys_eqs = eqs + eqs3  # + ineqs + eqs2 + weight_cons

        pass

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
        for j, cname in enumerate(self.comps[1].keys()):
            # if self.weight_cons[j].subs(ws) <= 0:
            if self.weight_cons[f'2{cname}'].subs(ws) <= 0:
                # sys_eqs.append(sympy.Eq(self.syms[f'eta{2}{cname}'], 0))
                # sys_eqs.append(self.syms[f'eta{2}{cname}'])
                vals[self.syms[f'eta{2}{cname}']] = 0
            else:
                flag = True
                # sys_eqs.append(self.weight_cons[j])
                sys_eqs.append(self.weight_cons[f'2{cname}'])
        if not flag:
            return res
        sys_eqs += self.Lag_sys_eqs
        sys_eqs = list(map(lambda x: x.subs(vals), sys_eqs))
        solutions = sympy.solve(sys_eqs)
        solutions = [s for s in solutions if all(map(lambda x: x.subs(s), self.pos_ws))]
        if len(solutions) != 1:
            print(f'failed to solve weights, got {len(solutions)} solutions')
            assert 0
        solution = solutions[0]
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

        if not self.initialized:
            for i in range(len(self.comps)):
                for j, (cname, _) in enumerate(self.comps[i].items()):
                    # mu = xmax - (xmax - xmin) / len(self.comps[i]) * (1 / 2 + j)
                    mu = X[:, 0].mean() + 0.5 * sigma - j * 0.5 * sigma
                    # self.comps[i][cname] = SN()
                    self.weights[i][cname] = 1 / len(self.comps[i])
                    self.comps[i][cname].mu = mu
                    self.comps[i][cname].sigma = sigma

            self.weights[1]['C'] *= 0.05

        self.starting_pos = deepcopy(self)

        comp_constraints['C_IC'] = RelativeConstraint(self.all_comps['C'], self.all_comps['IC'],
                                                      x_range=x, pdf=self.pdf_constrained,
                                                      cdf=self.cdf_constrained)
        comp_constraints['IC_C_w2'] = RelativeConstraint(self.all_comps['IC'], self.all_comps['C'],
                                                         weights=(self.weights[1]['IC'], self.weights[1]['C']),
                                                         x_range=x, pdf=self.weighted_pdf_constrained,
                                                         cdf=False)
        if self.ic2_comp:
            comp_constraints['IC_IC2'] = RelativeConstraint(self.all_comps['IC'], self.all_comps['IC2'],
                                                            x_range=x, pdf=self.pdf_constrained,
                                                            cdf=self.cdf_constrained)
            comp_constraints['I1_I2'] = RelativeConstraint(self.all_comps['I1'], self.all_comps['I2'],
                                                           x_range=x, pdf=self.pdf_constrained,
                                                           cdf=self.cdf_constrained)
        comp_constraints['IC_I1'] = RelativeConstraint(self.all_comps['IC'], self.all_comps['I1'],
                                                       x_range=x, pdf=self.pdf_constrained,
                                                       cdf=self.cdf_constrained)

        self.ll = self.log_likelihood(X)
        self.lls = [self.ll]
        self.slls = self.sep_log_likelihood(X)
        plt.figure(figsize=[18, 9])
        if self.show_plotting:
            plt.ion()
            plt.show()
            self.plot(X, self.lls, self.slls)
        prev_t = time.time()
        while abs(self.ll - prev_ll) > self.tolerance:
            prev_ll = self.ll
            rs = self.pred(X)

            # sum_rs = []
            new_weights = self.solve_weights([np.sum(r, 0) for r in rs], n)

            def proj_weights():
                # cons = RelativeConstraint(self.all_comps['IC'], self.all_comps['C'],
                #                           weights=(self.weights[1]['IC'], self.weights[1]['C']),
                #                           x_range=x, pdf=self.weighted_pdf_constrained,
                #                           cdf=False).getWeightChecker('left')
                cons = comp_constraints['IC_C_w2'].getWeightChecker('left')
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

            d = defaultdict(lambda: [None, [], []])
            for i, r in enumerate(rs):
                for j, (cname, cdist) in enumerate(self.comps[i].items()):
                    # self.weights[i][cname] = sum_rs[i][j] / n
                    d[cname][0] = cdist
                    d[cname][1].append(X[:, i])
                    d[cname][2].append(r[:, j])
                    # self.comps[i][cname].iterfit(X[:, i], r[:, j])

            for cname, (cdist, xs, rs) in d.items():
                cons = None
                if cname == 'C':
                    cons = ComposedChecker(
                        comp_constraints['C_IC'].getDistChecker('right'),
                        comp_constraints['IC_C_w2'].getDistChecker('left'),
                    )
                elif cname == 'IC':
                    if self.ic2_comp:
                        cons = ComposedChecker(
                            comp_constraints['C_IC'].getDistChecker('left'),
                            comp_constraints['IC_C_w2'].getDistChecker('right'),
                            comp_constraints['IC_IC2'].getDistChecker('right'),
                            comp_constraints['IC_I1'].getDistChecker('right'),
                        )
                    else:
                        cons = ComposedChecker(
                            comp_constraints['C_IC'].getDistChecker('left'),
                            comp_constraints['IC_C_w2'].getDistChecker('right'),
                            comp_constraints['IC_I1'].getDistChecker('right'),
                        )
                elif cname == 'IC2':
                    cons = comp_constraints['IC_IC2'].getDistChecker('left')
                elif cname == 'I1':
                    if self.ic2_comp:
                        cons = ComposedChecker(
                            comp_constraints['IC_I1'].getDistChecker('left'),
                            comp_constraints['I1_I2'].getDistChecker('right'),
                        )
                    else:
                        cons = ComposedChecker(
                            comp_constraints['IC_I1'].getDistChecker('left'),
                        )
                elif cname == 'I2':
                    cons = comp_constraints['I1_I2'].getDistChecker('left')

                # new_dist = cdist.iterfit(np.array(xs), np.array(rs), cons)
                cdist.iterfit(np.array(xs), np.array(rs), cons)
                # self.all_comps[cname] = new_dist
                # for i in range(len(self.comps)):
                #     if cname in self.comps[i]:
                #         self.comps[i][cname] = new_dist

            # ll = self.log_likelihood(X)
            self.ll = self.log_likelihood(X)
            # assert ll >= self.lls[-1]
            if self.ll < self.lls[-1]:
                print('ll decreased', f'{self.lls[-1]} -> {self.ll}')
            self.lls.append(self.ll)
            self.slls = self.sep_log_likelihood(X)

            cur_t = time.time()
            if self.show_plotting and cur_t - prev_t >= self.plot_interval:
                self.plot(X, self.lls, self.slls)
                prev_t = time.time()

        self.slls = self.sep_log_likelihood(X)
        self.plot(X, self.lls, self.slls)
        return self.ll, self.lls
