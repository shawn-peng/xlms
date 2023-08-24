import numpy as np
from scipy import stats
from param_binary_search import *
from normal_gpu import Normal


def trunc_norm_moments(mu, sigma):
    """Array trunc norm moments"""
    cdf = stats.norm.cdf(mu / sigma)
    flags = cdf == 0
    pdf = stats.norm.pdf(mu / sigma)
    p = pdf / cdf
    p[flags] = abs(mu[flags] / sigma)

    m1 = mu + sigma * p
    m2 = mu ** 2 + sigma ** 2 + sigma * mu * p
    return m1, m2


class SkewNormal:
    def __init__(self, alpha=1, mu=0, sigma=1):
        self.mu = mu
        self.sigma = sigma
        self.alpha = alpha
        self.delta = None
        self.Delta = None
        self.Gamma = None
        self.calc_alt_params()

    @property
    def mode(self):
        self.calc_alt_params()
        u_z = np.sqrt(2/np.pi) * self.delta
        sigma_z = np.sqrt(1 - u_z ** 2)
        gamma_1 = (4 - np.pi) / 2 * (self.delta * np.sqrt(2 / np.pi)) ** 3 / (1 - 2 * self.delta ** 2 / np.pi) ** (3/2)
        m0 = u_z - gamma_1 * sigma_z / 2 + np.sign(self.alpha) / 2 * np.exp(-2 * np.pi / np.abs(self.alpha))
        mode = self.mu + self.sigma * m0
        return mode

    def calc_alt_params(self):
        self.delta = self.alpha / np.sqrt(1 + self.alpha ** 2)
        self.Delta = self.sigma * self.delta
        self.Gamma = self.sigma ** 2 * (1 - self.delta ** 2)

    def from_alt_params(self):
        self.alpha = np.sign(self.Delta) * np.sqrt(self.Delta ** 2 / self.Gamma)
        self.sigma = np.sqrt(self.Gamma + self.Delta ** 2)

    def iterfit(self, X, weights, dist_cons=None):
        if dist_cons:
            assert dist_cons(self, 2e-7)
        # new_sn = self
        tn_mu = self.delta / self.sigma * (X - self.mu)
        tn_sigma = np.sqrt(1 - self.delta ** 2)
        # v = stats.truncnorm.moment(1, 0, np.Inf, tn_mu, tn_sigma)
        v, w = trunc_norm_moments(tn_mu, tn_sigma)
        # new_sn.mu = np.sum(weights * (X - v * self.Delta)) / np.sum(weights)
        mu = np.sum(weights * (X - v * self.Delta)) / np.sum(weights)
        if dist_cons:
            mu_test_func = dist_cons.get_param_checker('mu')
            mu = param_binary_search(self.mu, mu, mu_test_func)
        self.mu = mu
        # new_sn.mu = self.mu
        new_sn = SkewNormal(self.alpha, self.mu, self.sigma)
        new_sn.Delta = np.sum(weights * v * (X - self.mu)) / np.sum(weights * w)
        new_sn.Gamma = np.sum(weights * (
                (X - new_sn.mu)**2
                - 2 * v * (X - new_sn.mu) * new_sn.Delta
                + w * new_sn.Delta**2)
            ) / np.sum(weights)
        new_sn.from_alt_params()
        if dist_cons:
            sigma_test_func = dist_cons.get_param_checker('sigma')
            new_sn.sigma = param_binary_search(self.sigma, new_sn.sigma, sigma_test_func)
            self.sigma = new_sn.sigma
            alpha_test_func = dist_cons.get_param_checker('alpha')
            new_sn.alpha = param_binary_search(self.alpha, new_sn.alpha, alpha_test_func)
            self.alpha = new_sn.alpha
        # self.mu = new_sn.mu
        # self.sigma = new_sn.sigma
        # self.alpha = new_sn.alpha
        self.calc_alt_params()
        if dist_cons:
            assert dist_cons(self, 2e-7)
        return new_sn

    def cdf(self, x):
        return stats.skewnorm.cdf(x, self.alpha, self.mu, self.sigma)

    def pdf(self, x):
        # return stats.skewnorm.pdf(x, self.alpha, self.mu, self.sigma)
        r = 2 * Normal(self.mu, self.sigma).pdf(x) * Normal().cdf(self.alpha * (x - self.mu) / self.sigma)
        return r

    def __repr__(self):
        return f'<SkewNormal alpha={self.alpha}, mu={self.mu}, sigma={self.sigma}>'

