from copy import copy, deepcopy
import matplotlib.pyplot as plt


class RelativeConstraint:
    def __init__(self, right_dist, left_dist, x_range, weights=(1, 1), pdf=True, cdf=False):
        self.right_dist = right_dist
        self.left_dist = left_dist
        self.x_range = x_range
        self.weights = list(weights)
        self.pdf = pdf
        self.cdf = cdf

    def __repr__(self):
        return f'{self.right_dist, self.left_dist, self.weights}'

    def debug(self):
        plt.figure()
        mode_right = self.right_dist.mode
        x_pdf = self.x_range[self.x_range > mode_right]
        pr = self.right_dist.pdf(x_pdf)
        pl = self.left_dist.pdf(x_pdf)
        plt.plot(x_pdf, pr)
        plt.plot(x_pdf, pl)
        plt.plot(x_pdf, pr - pl)
        ax2 = plt.twinx()
        cr = self.right_dist.cdf(self.x_range)
        cl = self.left_dist.cdf(self.x_range)
        ax2.plot(self.x_range, cr)
        ax2.plot(self.x_range, cl)
        ax2.plot(self.x_range, cl - cr)

    def check(self, fuzzy=0):
        if self.pdf:
            mode_right = self.right_dist.mode
            x_pdf = self.x_range[self.x_range > mode_right]
            pr = self.right_dist.pdf(x_pdf)
            pl = self.left_dist.pdf(x_pdf)
            # for wr, wl in self.weights:
            wr, wl = self.weights
            # cond = wr * self.right_dist.pdf(x_pdf) >= wl * self.left_dist.pdf(x_pdf)
            cond = (wr * pr - wl * pl) >= -fuzzy
            if not cond.all():
                return False

        if self.cdf:
            x_cdf = self.x_range
            cr = self.right_dist.cdf(x_cdf)
            cl = self.left_dist.cdf(x_cdf)
            cond = (cl - cr) >= -fuzzy * 5e4
            if not cond.all():
                return False

        return True

    def getDistChecker(self, comp):
        return DistChecker(self, comp)

    def getWeightChecker(self, comp):
        return WeightChecker(self, comp)

    # def dist_check_func(self, dist):
    #     def param_check_func(param_name):
    #         def _check_func(param_val):
    #             setattr(dist, param_name, param_val)
    #             return self.check()
    #         return _check_func
    #     return param_check_func


# def compose_constraints(constraints):
# class ComposedConstraint:
#     def __init__(self, constraints):
#         self.constraints = constraints
#
#     def check(self):
#         for cons in self.constraints:
#             if not cons.check():
#                 return False
#         return True


class WeightChecker:
    def __init__(self, constraint, comp):
        self.constraint = constraint
        self.comp = comp
        assert self.comp in ['left', 'right']

    def __repr__(self):
        return str(self.__dict__)

    def __call__(self, w, fuzzy=0):
        # cons = copy(self.constraint)
        cons = self.constraint
        if self.comp == 'right':
            cons.weights[0] = w
        elif self.comp == 'left':
            cons.weights[1] = w
        return cons.check(fuzzy)


class DistChecker:
    def __init__(self, constraint, comp):
        self.constraint = constraint
        self.comp = comp
        assert self.comp in ['left', 'right']

    def __repr__(self):
        return str(self.__dict__)

    def __call__(self, dist, fuzzy=0):
        # cons = copy(self.constraint)
        cons = self.constraint
        if self.comp == 'left':
            return dist is cons.left_dist and cons.check(fuzzy)
            # cons.left_dist = dist
        elif self.comp == 'right':
            return dist is cons.right_dist and cons.check(fuzzy)
            # cons.right_dist = dist

        # return cons.check(fuzzy)

    def debug(self):
        self.constraint.debug()

    def get_param_checker(self, param_name):
        # cons = deepcopy(self.constraint)
        cons = self.constraint
        if self.comp == 'left':
            dist = cons.left_dist
        elif self.comp == 'right':
            dist = cons.right_dist

        def _check_func(param_val, fuzzy=0):
            setattr(dist, param_name, param_val)
            return cons.check(fuzzy)

        return _check_func


class ComposedChecker:
    def __init__(self, *checkers):
        self.checkers = checkers

    def __repr__(self):
        return str(self.__dict__)

    def __call__(self, dist, fuzzy=0):
        flags = [checker(dist, fuzzy) for checker in self.checkers]
        return all(flags)

    def debug(self):
        for checker in self.checkers:
            checker.constraint.debug()

    def get_param_checker(self, param_name):
        def _check_func(param_val, fuzzy=0):
            return all(checker.get_param_checker(param_name)(param_val, fuzzy) for checker in self.checkers)

        return _check_func