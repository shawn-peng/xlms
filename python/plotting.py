import matplotlib.pyplot as plt
from scipy import stats
import numpy as np


def weighted_hist(x, w=None, bins=100):
    h = np.histogram(x, bins, density=True, weights=w)
    if w is not None:
        plt.bar(h[1][:-1], h[0] * w.mean())
    else:
        plt.bar(h[1][:-1], h[0])


def weighted_hist2d(X, w=None, bins=100):
    h = np.histogram2d(X[:, 0], X[:, 1], bins=bins, density=True, weights=w)
    H, xedges, yedges = h
    if w is not None:
        H *= w.mean()
    plt.imshow(H, interpolation='nearest', origin='lower',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
    plt.colorbar()
