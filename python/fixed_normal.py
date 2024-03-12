from normal import Normal


class FixedNormal(Normal):
    def iterfit(self, X, weights, cons=None):
        """same as Normal, except do nothing on fitting"""
