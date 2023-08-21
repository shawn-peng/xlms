from scipy import stats

import matplotlib.pyplot as plt
import numpy as np

from skew_normal import SkewNormal
from xlms import XLMS_Dataset


def skew_normal_test():
    d = SkewNormal(5, 1, 2)
    # d = stats.skewnorm(5, 1, 2)
    x = np.arange(-1, 8, 0.01)
    y = d.pdf(x)
    plt.plot(x, y)
    y = d.cdf(x)
    ax2 = plt.twinx()
    ax2.plot(x, y)
    plt.show()


# skew_normal_test()


def xlms_dataset_test():
    ds = XLMS_Dataset('alban')


xlms_dataset_test()
