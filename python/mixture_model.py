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
import traceback as tb

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
                 max_iteration=10000,
                 binwidth=2.5,
                 plotstep=0.1,
                 show_plotting=False,
                 plot_interval=1,
                 initialized=False,
                 # title=None,
                 event_notify_func=None,
                 init_strategy=None,
                 seedoff=0,
                 **kwargs):
        """

        :param skew_dirs: skew directions
        :param ic2_comp: model IC2 component separately
        :param constraints: options [ weights, mode, pdf, cdf, weighted_pdf ]
        :param tolerance: convergence threshold
        :param max_iteration: maximum iteration will go through
        :param binwidth: binwidth for histogram plotting
        :param plotstep: step for the grid when plotting curves (e.g. pdf, cdf)
        :param show_plotting: plot each step during running of the program
        :param title: title for the figure
        """
        super().__init__(**kwargs)
        self.skew_dirs = skew_dirs
        self.tolerance = tolerance
        self.max_iteration = max_iteration
        self.binwidth = binwidth
        self.plotstep = plotstep
        self.show_plotting = show_plotting
        self.plot_interval = plot_interval
        self.initialized = initialized
        self.event_notify_func = event_notify_func
        self.ic2_comp = ic2_comp
        self.seedoff = seedoff
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
                    self.all_comps[cname] = SN(skew_dirs[cname], name=cname)

        self.all_comps = dict(sorted(self.all_comps.items(), key=lambda x: comp_id[x[0]]))

        # self.comps = [{
        #     cname: SN(skew_dirs[cname]) for cname in sorted(set(comps), key=lambda x: comp_id[x])
        # } for comps in self.join_comps]
        self.comps = [{
            cname: self.all_comps[cname] for cname in sorted(set(comps), key=lambda x: comp_id[x])
        } for comps in self.join_comps]
        self.weights = [{
            cname: DynamicParam(1) for cname in sorted(set(comps), key=lambda x: comp_id[x])
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
        # for i in range(self.n_samples):
        #     for cname in self.comps[i].keys():
        #         eqs[ind[(i, cname)]] *= self.syms[f'w{i + 1}{cname}']

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
                # truncate_zero(pj[i][:, j])
            p0 = np.sum(pj[i], 1).reshape((-1, 1))
            truncate_zero(p0)
            # assert not np.any(p0)
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
        violated_cons = []
        for j, cname in enumerate(self.comps[1].keys()):
            # if self.weight_cons[j].subs(ws) <= 0:
            if self.weight_cons[f'2{cname}'].subs(ws) <= 0:
                # constraint is satified, set it to inactive
                pass
                # sys_eqs.append(sympy.Eq(self.syms[f'eta{2}{cname}'], 0))
                # sys_eqs.append(self.syms[f'eta{2}{cname}'])
                # vals[self.syms[f'eta{2}{cname}']] = 0
            else:
                flag = True
                violated_cons.append(cname)
                # violated_cons += 1
                # sys_eqs.append(self.weight_cons[j])
                # sys_eqs.append(self.weight_cons[f'2{cname}'])

        if not flag:
            return res

        if len(violated_cons) > 1:
            pass
            # print(f'got {violated_cons} active constraints, might unable to find solution')
            # issue solved
            # return res

        # usually we get two at most, if we have two violated cons
        # try two cons separately, check wether all cons are satisfied

        def choose_n(l, n):
            if n == 0:
                yield ()
                return
            if n > len(l):
                return
            x = l[0]
            for r in choose_n(l[1:], n - 1):
                yield (x,) + r
            yield from choose_n(l[1:], n)

        def gen_sys_eq(active_cnames):
            # print(active_cnames)
            newsys, newvals = deepcopy(sys_eqs), deepcopy(vals)
            # for cname in violated_cons:
            for cname in self.comps[1].keys():
                if cname in active_cnames:
                    newsys.append(self.weight_cons[f'2{cname}'])
                else:
                    # newvals[self.syms[f'eta{2}{cname}']] = 0
                    newsys.append(self.syms[f'eta{2}{cname}'])
            return newsys, newvals, active_cnames

        def all_sys_eqs():
            for i in range(len(violated_cons)):
                for selected_cnames in choose_n(violated_cons, i + 1):
                    yield gen_sys_eq(selected_cnames)

        def get_var_syms_n_guess():
            syms = []
            guess = []
            for i in range(self.n_samples):
                for cname in self.comps[i].keys():
                    syms.append(self.syms[f'w{i + 1}{cname}'])
                    guess.append(res[i][cname])

            for i in range(self.n_samples):
                syms.append(self.syms[f'lambda{i + 1}'])
                guess.append(n)

            for cname in self.comps[1].keys():
                syms.append(self.syms[f'eta2{cname}'])
                guess.append(0)
            return syms, guess

        def find_neg_w(solution):
            for i in range(self.n_samples):
                for cname in self.comps[i].keys():
                    sym = self.syms[f'w{i + 1}{cname}']
                    if solution[sym] < 0:
                        return sym

        def solve_2IC2_2I1(sys_eqs, vals):
            w2 = self.syms['w2I1']
            eta2 = self.syms['eta2I1']
            lambda2 = self.syms['lambda2']
            R1C = self.syms['R1C']
            R1IC = self.syms['R1IC']
            R1I1 = self.syms['R1I1']
            R2C = self.syms['R2C']
            R2IC = self.syms['R2IC']
            R2I1 = self.syms['R2I1']
            R2IC2 = self.syms['R2IC2']
            R2I2 = self.syms['R2I2']
            a = R1C + R1IC + R2I1
            partial_sys_eqs = [
                w2 ** 2 - w2 + R2I1 / n,
                eta2 - (w2 * (n - a) - R2I1) / (w2 * (1 - w2)),
                lambda2 + n - w2 * eta2 - 2 * w2 * n + a,
            ]  # + sys_eqs
            sym_vars, _ = get_var_syms_n_guess()
            solutions = sympy.solve(partial_sys_eqs, [w2, eta2, lambda2])
            new_sys = list(unify_vec(sys_eqs, [w2, eta2, lambda2], solutions))
            pass

        def unify(sys_eqs, solutions):
            for solution in solutions:
                sym_vars = list(solution.keys())
                vals = list(solution.values())
                yield sympy.substitution(sys_eqs, sym_vars, vals)

        def unify_vec(sys_eqs, syms, solutions):
            for solution in solutions:
                yield sympy.substitution(sys_eqs, syms, solution)

        def solve(sys_eqs, vals, active_cons):
            # if active_cons == ('IC2', 'I1'):
            #     return solve_2IC2_2I1(sys_eqs, vals)
            nactive = len(active_cons)
            use_nsolve = False
            # use_nsolve = True
            if nactive > 1:
                use_nsolve = True
            sys_eqs += self.Lag_sys_eqs
            sys_eqs = list(map(lambda x: x.subs(vals), sys_eqs))
            # sys_eqs = [x for x in sys_eqs if x != 0]
            try:
                # if active_cons == ('IC2', 'I1'):
                #     solutions = solve_2IC2_2I1(sys_eqs, vals)
                syms, guess = get_var_syms_n_guess()
                if nactive > 1:
                    # solutions = sympy.nsolve(sys_eqs, syms, guess)
                    # solutions = [dict(zip(syms, solutions[:, j]))
                    #              for j in range(solutions.shape[1])]
                    solutions = sympy.nonlinsolve(sys_eqs, *syms)
                    solutions = [dict(zip(syms, s)) for s in solutions]
                    # print(solutions)
                else:
                    # print(sys_eqs)
                    solutions = sympy.nsolve(sys_eqs, syms, guess)
                    solutions = [dict(zip(syms, solutions[:, j]))
                                 for j in range(solutions.shape[1])]
                    # solutions = sympy.solve(sys_eqs)
            except Exception as e:
                # if e is not ZeroDivisionError:
                #     tb.print_exc()
                return

            solutions = [s for s in solutions if all(map(lambda x: x.subs(s), self.pos_ws))]
            if len(solutions) < 1:
                return
            solution = solutions[0]
            return solution

        def extract_res(solution):
            for i in range(self.n_samples):
                for j, cname in enumerate(self.comps[i].keys()):
                    # self.weights[i][cname] = float(solution[self.syms[f'w{i + 1}{cname}']])
                    w = float(solution[self.syms[f'w{i + 1}{cname}']].evalf(chop=True))
                    res[i][cname] = w
                    ws[self.syms[f'w{i + 1}{cname}']] = w

        def check_cons(ws):
            for j, cname in enumerate(self.comps[1].keys()):
                r = float(self.weight_cons[f'2{cname}'].subs(ws))
                # assert r < 0 or np.isclose(r, 0)
                # print(r)
                if not (r < 0 or np.isclose(r, 0)):
                    return False
            return True

        flag = False
        # solutions = []

        # Start from one cons, check if any solution for one cons satifies all cons
        # if not, try two cons, three cons, etc.
        for solution in starmap(solve, all_sys_eqs()):
            if not solution:
                continue

            extract_res(solution)
            if not check_cons(ws):
                continue

            flag = True
            break
            # solutions.append(solution)

        # assert flag
        if not flag:
            return res
            assert False

        # if len(solutions) != 1:
        #     print(f'failed to solve weights, got {len(solutions)} solutions, return unconstrained weights')
        #     return res

        return res

    def putdata(self, X):
        X = X[X[:, 1] != 0].astype(np.float128)
        self.X = X
        xmax = np.max(X)
        xmin = -100
        self.xplot = np.arange(xmin, xmax + self.plotstep, self.plotstep)
        bins = np.arange(xmin, xmax + self.binwidth, self.binwidth)
        self.hist = [
            np.histogram(X[:, 0], bins, density=True),
            np.histogram(X[:, 1], bins, density=True),
        ]

    def check_constraints(self):
        for c, dist in self.all_comps.items():
            if not self.comp_constraints[c](dist, 2e-7):
                return False
        return True

    def log(self, *args):
        print(f'{self.title} id {self.seedoff}:', *args)

    def init_model(self, X):
        X = X.astype(np.float128)
        self.log('start init ...')
        self.init_range(X)

        self.create_constraints()

        sigma = np.sqrt(X[:, 0].var())
        mu = X[:, 0].mean()
        xmin = X.min()
        xmax = X.max()
        xmin = -50.0
        if self.init_strategy == 'random':
            # seed = int(time.time()) + self.seedoff
            seed = self.seedoff + 8
            # seed = 4
            print(f'seed {seed}')
            np.random.seed(seed)

            for i in range(len(self.comps)):
                for j, (cname, _) in enumerate(self.comps[i].items()):
                    self.weights[i][cname].set(np.float128(1 / len(self.comps[i])))
            self.weights[1]['C'] *= np.float128(0.001)

            frozen_model = self.frozen()
            while True:
                # mu = xmax
                j = 0
                # mus = np.random.uniform(xmin, xmax, len(self.all_comps))
                # mus[::-1].sort()
                IC = self.all_comps['IC']
                IC.mu = np.random.uniform(mu - 1.5 * sigma, mu + 1.5 * sigma)
                self.all_comps['C'].mu = np.random.uniform(IC.mu, xmax)
                self.all_comps['I1'].mu = np.random.uniform(xmin, IC.mu)
                self.all_comps['IC2'].mu = np.random.uniform(self.all_comps['I1'].mu, IC.mu)
                self.all_comps['I2'].mu = np.random.uniform(xmin, self.all_comps['I1'].mu)
                for cname, cdist in self.all_comps.items():
                    # mu_offset = np.random.uniform(0, 1)
                    # print(f'{cname} mu_offset {mu_offset}')
                    sigma_scale = np.random.uniform(0.5, 1.0)
                    # sigma_scale = 1.0
                    alpha_scale = np.random.uniform(0.0, 2.0)
                    # alpha_scale = 1.0
                    # mu -= 2 * mu_offset * sigma
                    # cdist.mu = np.float128(mu)
                    # cdist.mu = np.float128(mu - (j + mu_offset) * sigma)
                    cdist.sigma = np.float128(sigma)
                    cdist.sigma *= np.float128(sigma_scale)
                    cdist.alpha = frozen_model.all_comps[cname].alpha * np.float128(alpha_scale)
                    cdist.calc_alt_params()
                    j += 1
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
                    self.weights[i][cname].set(np.float128(1 / len(self.comps[i])))
                    self.comps[i][cname].mu = np.float128(mu)
                    self.comps[i][cname].sigma = np.float128(sigma)
                    self.comps[i][cname].calc_alt_params()
            self.weights[1]['C'] *= np.float128(0.001)

        self.ll = self.log_likelihood(X)
        self.lls = [self.ll]
        self.slls = self.sep_log_likelihood(X)

        self.starting_pos = self.frozen()
        self.initialized = True

        if self.show_plotting:
            self.plot(X, self.lls, self.slls)
        if self.event_notify_func:
            self.event_notify_func('update', self.starting_pos)
        self.log('init finished ...')

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
                                                             # weights=(lambda: self.weights[1]['IC'],
                                                             #          lambda: self.weights[1]['C']),
                                                             x_range=x, pdf=self.weighted_pdf_constrained,
                                                             cdf=False, mode=False, left_tail=False)
        relative_constraints['IC_I1'] = RelativeConstraint(self.all_comps['IC'], self.all_comps['I1'],
                                                           x_range=x, mode=self.mode_constrained,
                                                           pdf=self.pdf_constrained,
                                                           cdf=self.cdf_constrained)
        if self.ic2_comp:
            relative_constraints['IC_IC2'] = RelativeConstraint(self.all_comps['IC'], self.all_comps['IC2'],
                                                                x_range=x, mode=self.mode_constrained,
                                                                pdf=self.pdf_constrained,
                                                                cdf=self.cdf_constrained)
            # relative_constraints['IC_IC2_w2'] = RelativeConstraint(self.all_comps['IC'], self.all_comps['IC2'],
            #                                                        weights=(self.weights[1]['IC'],
            #                                                                 self.weights[1]['IC2']),
            #                                                        x_range=x, pdf=self.weighted_pdf_constrained,
            #                                                        cdf=False, mode=False, left_tail=False)
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
                # relative_constraints['IC_IC2_w2'].getDistChecker('right'),
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
                # relative_constraints['IC_IC2_w2'].getDistChecker('left'),
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

    def update_weights(self, new_weights):
        for i in range(self.n_samples):
            for cname, w in self.weights[i].items():
                w.set(new_weights[i][cname])

    def fit(self, X):
        prev_ll = -np.inf
        # X = X[X[:, 1] != 0].astype(np.float128)
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
            plt.ion()
            self.plot(X, self.lls, self.slls)
            pass
        prev_t = time.time()
        while abs(self.ll - prev_ll) > self.tolerance and not self.stopped and len(self.lls) <= self.max_iteration:
            # print(f'iteration {len(self.lls)}')
            prev_ll = self.ll
            rs = self.pred(X)

            # print('ll', self.ll)
            # sum_rs = []
            new_weights = self.solve_weights([np.sum(r, 0) for r in rs], n)

            # print(f'solved weights')

            def proj_weights():
                cons = self.relative_constraints['IC_C_w2'].getWeightChecker('left')
                old_w2c = self.weights[1]['C'].get()

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
            self.update_weights(new_weights)
            # self.weights = new_weights
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
                self.log('ll decreased', f'{self.lls[-1]} -> {self.ll}')
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
            # print(cur_t)
            if cur_t - prev_t >= self.plot_interval:
                self.log(f'{len(self.lls)} iterations')
                if self.show_plotting:
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
