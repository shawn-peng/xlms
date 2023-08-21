import torch.distributions
from scipy import stats
import numpy as np
from param_binary_search import *
import matplotlib.pyplot as plt
# from torchquad import Simpson, MonteCarlo, set_up_backend
import torchquad
from numbers import Real

import gpu

# integrator = Simpson()
# integrator = MonteCarlo()
# integrator = torchquad.Gaussian()
# integrator = torchquad.Boole()
integrator = torchquad.VEGAS()

inf = np.inf


class Normal:
    def __init__(self, mu=0.0, sigma=1.0):
        self.mu = mu
        self.sigma = sigma

    def iterfit(self, X, weights):
        self.mu = np.sum(weights * X) / np.sum(weights)
        self.sigma = np.sqrt(np.sum(weights * (X - self.mu) ** 2) / np.sum(weights))
        pass

    def pdf(self, X):
        return stats.norm.pdf(X, self.mu, self.sigma)

    def cdf(self, X):
        return stats.norm.pdf(X, self.mu, self.sigma)
        # # f = lambda x: torch.exp(torch.distributions.Normal(self.mu, self.sigma).log_prob(x))
        # var = (self.sigma ** 2)
        # # log_scale = torch.math.log(self.sigma) if isinstance(self.sigma, Real) else self.sigma.log()
        # # -((X - self.mu) ** 2) / (2 * var) - log_scale - torch.math.log(torch.math.sqrt(2 * torch.math.pi))
        # new_mus = - (X - self.mu) / self.sigma
        # # new_mus = [0]
        # new_mus = torch.Tensor(new_mus)
        # lower_bound = new_mus.abs().max() * -100
        # f = lambda x: torch.exp(torch.Tensor(-((x - new_mus) ** 2 / 2))) / torch.math.sqrt(2 * torch.math.pi)
        # f2 = lambda x: torch.exp(torch.Tensor(-((x) ** 2 / 2))) / torch.math.sqrt(2 * torch.math.pi)
        # # f = lambda x: torch.exp(torch.distributions.Normal(self.mu, self.sigma).log_prob(x))
        # res = integrator.integrate(f, dim=1, N=9999, integration_domain=[[-1e3, 0]], backend="torch", )
        # res2 = integrator.integrate(f, dim=1, N=9999, integration_domain=[[lower_bound, 0]], backend="torch", )
        # print(res)
        # return np.array(res)

    def plot(self, X, weights):
        d = X.shape[1]
        plt.figure()

        xmin = X.min()
        xmax = X.max()
        step = 0.1
        x = np.arange(xmin, xmax + step, step)
        p = self.pdf(x)
        for i in range(d):
            ax = plt.subplot(d, 2, i * 2 + 1)
            ax.hist(X[:, i], 100, weights=weights[:, i], density=True)
            ax.plot(x, p)
            plt.xlim([0, 300])
            ax = plt.subplot(d, 2, i * 2 + 2)
            ax.scatter(X[:, 0], X[:, 1], c=weights[:, i], s=1)

    def __repr__(self):
        return f'<Normal mu={self.mu}, sigma={self.sigma}>'
