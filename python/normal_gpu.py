import torch
import torch.distributions
import matplotlib.pyplot as plt
import math

import gpu


class Normal:
    def __init__(self, mu=0.0, sigma=1.0):
        self.mu = mu
        self.sigma = sigma

    def iterfit(self, X, weights):
        self.mu = torch.sum(weights * X) / torch.sum(weights)
        self.sigma = torch.sqrt(torch.sum(weights * (X - self.mu) ** 2) / torch.sum(weights))
        pass

    def pdf(self, X):
        X = torch.tensor(X)
        p = torch.exp(-((X - self.mu) / self.sigma) ** 2 / 2) / (self.sigma * math.sqrt(2 * math.pi))
        return p.cpu().numpy()

    def cdf(self, X):
        res = ((1 + torch.erf((torch.Tensor(X) - self.mu) / (self.sigma * math.sqrt(2)))) / 2)
        return res.cpu().numpy()

    def plot(self, X, weights):
        d = X.shape[1]
        plt.figure()

        xmin = X.min()
        xmax = X.max()
        step = 0.1
        x = torch.arange(xmin, xmax + step, step)
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
