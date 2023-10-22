import json
import multiprocessing
# multiprocessing.set_start_method('spawn')
import traceback

import os

os.environ['OPENBLAS_NUM_THREADS'] = '1'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
import sys
import glob
import subprocess
from functools import *
from itertools import *
import argparse

from mixture_model_base import MixtureModelBase
from mixture_model import MixtureModel
from mixture_model_1S import MixtureModel1S
from order_stat_mixture_model import OrderStatMixtureModel
from order_stat_mixture_model_v2 import OrderStatMixtureModelV2
from xlms import XLMS_Dataset
from myutils import *

import normal
import skew_normal

Norm = normal.Normal
SN = skew_normal.SkewNormal

map_cons_str = {
    'weights':      'w',
    'mode':         'mode',
    'pdf':          'pdf',
    'cdf':          'cdf',
    'weighted_pdf': 'wpdf',
}


def get_cons_str(constraints):
    return ';'.join(map(lambda x: map_cons_str[x], constraints))


config = sys.argv[1]

model_samples = 2
gaussian_model = False
ic2_comp = True
dir_suffix = ''
alpha_bases = [0.] if gaussian_model else [2.]

model_class = f'{model_samples}S{"g" if gaussian_model else ""}{"2" if ic2_comp else ""}'

base_figure_dir = f'../figures_python_diffsign_{model_class}_{config}{dir_suffix}' \
                  f'_initskew_{"_".join([f"{alpha_base:.0f}" for alpha_base in alpha_bases])}'

basic_settings = {}

settings = {
    'no_constraint':                            {**basic_settings, 'constraints': []},
    'weight_constraints':                       {**basic_settings, 'constraints': ['weights']},
    'unweighted_pdf_no_weight_constraints':     {**basic_settings, 'constraints': ['pdf']},
    'unweighted_cdf_no_weight_constraints':     {**basic_settings, 'constraints': ['cdf']},
    'unweighted_pdf_cdf_no_weight_constraints': {**basic_settings, 'constraints': ['pdf', 'cdf']},
    'unweighted_pdf':                           {**basic_settings, 'constraints': ['weights', 'pdf']},
    'unweighted_cdf':                           {**basic_settings, 'constraints': ['weights', 'cdf']},
    'unweighted_pdf_cdf':                       {**basic_settings, 'constraints': ['weights', 'pdf', 'cdf']},
    'all_constraints':                          {**basic_settings,
                                                 'constraints': ['weights', 'pdf', 'cdf', 'weighted_pdf']},
    'unweighted_pdf_cdf_mode':                  {**basic_settings, 'constraints': ['weights', 'mode', 'pdf', 'cdf']},
    'unweighted_pdf_mode':                      {**basic_settings, 'constraints': ['weights', 'mode', 'pdf']},
    'unweighted_cdf_mode':                      {**basic_settings, 'constraints': ['weights', 'mode', 'cdf']},
    'unweighted_pdf_cdf_weighted_pdf_mode':     {**basic_settings,
                                                 'constraints': ['weights', 'mode', 'pdf', 'cdf', 'weighted_pdf']},
    'unweighted_pdf_weighted_pdf_mode':         {**basic_settings,
                                                 'constraints': ['weights', 'mode', 'pdf', 'weighted_pdf']},
    'unweighted_cdf_weighted_pdf_mode':         {**basic_settings,
                                                 'constraints': ['weights', 'mode', 'cdf', 'weighted_pdf']},
}


def plot_dataset(dataset_name):
    tda_info = json.load(open(f'../results/info/{dataset_name}.json'))
    dataset = XLMS_Dataset(dataset_name)
    dataset.mat = dataset.mat[:, dataset.mat[1, :] != 0]

    res_dir = f'{base_figure_dir}/{dataset_name}'
    models = pickle.load(open(f'{res_dir}/models.pickle', 'rb'))

    best = max(models, key=lambda x: x['ll'])
    # best = max(models, key=lambda x: x['slls'][0])
    if 'slls' not in best:
        best['slls'] = best['model'].sep_log_likelihood(dataset.mat.T)
    if best['ll'] == -np.inf:
        print('no solution for', dataset_name)
        return
    tda_fdr1 = tda_info['fdr_thres']
    # best['model'].plot(dataset.mat.T, best['lls'], best['model'].sep_log_likelihood(dataset.mat.T))
    fig = plt.figure(figsize=(16, 9))
    MixtureModelBase._plot(dataset.mat.T, best['lls'], best['slls'], best['model'], fig)
    # plt.subplot(3, 1, 1)
    ax = plt.gcf().axes[0]
    plt.axes(ax)
    plt.axvline(tda_fdr1, linestyle='--')
    plt.text(tda_fdr1, 0.003, '$\leftarrow$ TDA 1% FDR threshold')
    if 'cons_sat' in best:
        plt.title(
            f"{dataset_name} {best['sls']} ll={best['ll']:.05f}"
            f" constraints={get_cons_str(settings[config]['constraints'])}"
            f" {'Y' if best['cons_sat'] else 'N'}")
    else:
        plt.title(
            f"{dataset_name} {best['sls']} ll={best['ll']:.05f}"
            f" constraints={get_cons_str(settings[config]['constraints'])}")
    plt.savefig(f'{res_dir}/best.png')


print(base_figure_dir)
for d in glob.glob(base_figure_dir + '/*'):
    if not os.path.isdir(d):
        continue
    print(d)
    dataset_name = os.path.basename(d)
    print(dataset_name)
    plot_dataset(dataset_name)
