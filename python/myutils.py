import time
import json
import pprint

import numpy as np


class TimeMeter:
    def __init__(self):
        self.s = time.time_ns()

    def read(self):
        t = time.time_ns()
        return (t - self.s) / 1e6


def truncate_zero(m):
    flags = m == 0
    m[flags] = m[~flags].min()


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


class NamedArray(object):
    def __init__(self, fields, val=0):
        self.idx = {f: i for i, f in enumerate(fields)}
        self.arr = np.zeros(len(self.idx))
        self.arr[:] = val

    def __getitem__(self, item):
        idx = self.idx[item]
        return self.arr[idx:idx + 1].squeeze()

    def __setitem__(self, item, val):
        idx = self.idx[item]
        self.arr[idx] = val


class AttrObj(object):
    def __init__(self, attrs):
        self.attrs = attrs

    def __getstate__(self):
        return self.attrs

    def __setstate__(self, state):
        self.attrs = state

    def __getattr__(self, item):
        if item in self.attrs:
            return self.attrs[item]
        # print(f'searching {item} in super')
        return super().__getattr__(item)

    def __str__(self):
        return pprint.pformat(self.attrs, indent=2)

    def __repr__(self):
        return str(self)


class DynamicParam:
    def __init__(self, val, getter=None, setter=None):
        self.val = val
        # self.getter = getter
        # self.setter = setter

    def __mul__(self, other):
        return self.val * other

    def __str__(self):
        return str(self.val)

    def __repr__(self):
        return str(self.val)

    def __format__(self, format_spec):
        return self.val.__format__(format_spec)

    def __float__(self):
        return self.val

    def get(self):
        return self.val

    def set(self, value):
        self.val = value

    def __imul__(self, other):
        self.val *= other
        return self
