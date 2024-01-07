import os
import pickle

import numpy as np

from xlms import XLMS_Dataset

datasets = [
    'alban',
    'Alinden',
    'ALott',
    'CPSF',
    'D1810',
    'ecoli_xl',
    # 'KKT4',
    'MS2000225',
    'peplib',
    'QE',
    'RPA',
    # 'Gordon',
]

MAX_SAMPLE_SIZE = 50000

NUM_BOOTSTRAP = 50


def gen_bootstrap_index(dataset_name):
    # dataset = XLMS_Dataset(dataset_name)
    dataset = XLMS_Dataset(dataset_name, nodup=True)
    bs_indices = []
    for i in range(NUM_BOOTSTRAP):
        rind = gen_bootstrap_index_i(dataset, rseed=i)
        bs_indices.append(rind)
    bs_indices = np.vstack(bs_indices)
    return bs_indices


def gen_bootstrap_index_i(dataset, rseed=0):
    n = dataset.mat.shape[1]
    nsample = min(MAX_SAMPLE_SIZE, n)
    np.random.seed(rseed)
    rind = np.random.choice(np.arange(n), nsample, replace=True)

    return rind


if __name__ == '__main__':
    bootstrap_dir = '../bootstrap/'
    if not os.path.exists(bootstrap_dir):
        os.makedirs(bootstrap_dir)

    index = {}
    for dataset_name in datasets:
        dataset_indices = gen_bootstrap_index(dataset_name)
        index[dataset_name] = dataset_indices

    pickle.dump(index, open(f'{bootstrap_dir}/index.pickle', 'wb'))



