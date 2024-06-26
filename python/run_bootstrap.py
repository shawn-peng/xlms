import json
import multiprocessing
import shutil
# multiprocessing.set_start_method('spawn')
import traceback

import os

os.environ['OPENBLAS_NUM_THREADS'] = '1'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
import sys
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

SAMPLE_SIZE = 50000
tolerance = 1e-8
max_iteration = 5000
show_plotting = True
# show_plotting = False
plot_interval = 5
# gaussian_model = True
gaussian_model = False
model_samples = 2
ic2_comp = True
# ic2_comp = False
# init_strategy = None
init_strategy = 'random'
mu_strategy = 'gaussian'
random_size = 2
parallel = True
# parallel = False
# inner_parallel = True
inner_parallel = False
num_workers = 20
if init_strategy == 'random' and parallel and inner_parallel:
    show_plotting = False

if os.path.exists('on_server'):
    show_plotting = False

# alpha_base = 0. if gaussian_model else 2.
alpha_bases = [0.] if gaussian_model else [1., 2., 5.]
# alpha_bases = [0.] if gaussian_model else [2.]

run_all = True
# run_all = False
# dataset_to_run = 'ecoli_xl'
# dataset_to_run = 'alban'
dataset_to_run = 'Alinden'
# dataset_to_run = 'D1810'
# dataset_to_run = 'MS2000225'
if run_all:
    show_plotting = False

common_settings_2s = {}

# config = 'common'
# config = 'no_constraint'
# config = 'weight_constraints'
# config = 'unweighted_pdf'
# config = 'unweighted_cdf'
# config = 'unweighted_pdf_cdf'
# config = 'unweighted_pdf_no_weight_constraints'
# config = 'unweighted_cdf_no_weight_constraints'
# config = 'unweighted_pdf_cdf_no_weight_constraints'
config = 'unweighted_pdf_mode'
# config = 'unweighted_cdf_mode'
# config = 'unweighted_pdf_cdf_mode'
# config = 'unweighted_pdf_weighted_pdf_mode'
# config = 'unweighted_cdf_weighted_pdf_mode'
# config = 'unweighted_pdf_cdf_weighted_pdf_mode'
# if len(sys.argv) > 1:
#     config = sys.argv[1]

dir_suffix = '_4'

parser = argparse.ArgumentParser(prog='XLMS')
parser.add_argument('-d', '--dataset', default=dataset_to_run)
parser.add_argument('-c', '--config', default=config)
parser.add_argument('-s', '--model_samples', type=int, default=model_samples)
parser.add_argument('--suffix', default=dir_suffix)
parser.add_argument('-j', '--jobs', type=int, default=num_workers)
parser.add_argument('-r', '--random_size', type=int, default=random_size)
parser.add_argument('-q', '--part', type=int, default=-1)
parser.add_argument('-o', '--random_i', type=int, default=-1)
parser.add_argument('-t', '--tolerance', type=float, default=tolerance)
parser.add_argument('-a', '--all', action='store_true', default=False)
parser.add_argument('-p', '--parallel', action='store_true', default=False)
parser.add_argument('-i', '--inner_parallel', action='store_true', default=False)
parser.add_argument('--mu_strategy', default=mu_strategy)
parser.add_argument('--bootstrap_i', type=int, default=0)

args = parser.parse_args()

print(args)

dataset_to_run = args.dataset
config = args.config
model_samples = args.model_samples
dir_suffix = args.suffix
tolerance = args.tolerance
run_all = args.all
parallel = args.parallel
inner_parallel = args.inner_parallel
num_workers = args.jobs
random_size = args.random_size
part = args.part
random_i = args.random_i
mu_strategy = args.mu_strategy
bootstrap_i = args.bootstrap_i

map_cons_str = {
    'weights':      'w',
    'mode':         'mode',
    'pdf':          'pdf',
    'cdf':          'cdf',
    'weighted_pdf': 'wpdf',
}

basic_settings = {'ic2_comp':      ic2_comp,
                  'tolerance':     tolerance,
                  'show_plotting': show_plotting,
                  'plot_interval': plot_interval,
                  'init_strategy': init_strategy,
                  'mu_strategy':   mu_strategy,
                  'max_iteration': max_iteration,
                  }

settings = {
    'common':                                   {'tolerance': tolerance, 'show_plotting': show_plotting},
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


def get_cons_str(constraints):
    return ';'.join(map(lambda x: map_cons_str[x], constraints))


bootstrap_dir = '../bootstrap/'

model_class = f'{model_samples}S{"g" if gaussian_model else ""}{"2" if ic2_comp else ""}'
base_figure_dir = None


def capture_args(locals):
    # return {k: locals[k] for k in ['dataset_name', 'dataset', 'res_dir']}
    return [locals[k] for k in ['dataset_name', 'dataset', 'res_dir']]


def run_model(sls, dataset_name, dataset, tda_info, res_dir, modelid=0):
    title = f"({dataset.mat.shape[1] / 1000:.1f}k) {dataset_name} constraints={get_cons_str(settings[config]['constraints'])}"
    print('model:', modelid)
    if model_samples == 1:
        model = MixtureModel1S(sls, **settings[config], title=title, seedoff=modelid)
    elif model_samples == 2:
        model = MixtureModel(sls, **settings[config], title=title, seedoff=modelid)
    else:
        assert False, 'Invalid model samples'
    try:
        ll, lls = model.fit(dataset.mat.T)
        pass
    except Exception as e:
        print(f'Exception {traceback.format_exc()} happened for {dataset_name}, stopped at middle')
        dump_pickle = f'{res_dir}/model_{modelid}.pickle'
        pickle.dump(model.frozen(), open(dump_pickle, 'wb'))
        ll = model.ll
        lls = model.lls
        # model.plot(dataset.mat.T, model.lls, model.slls)
    # ax = plt.gcf().axes[0]
    # plt.axes(ax)
    # tda_fdr1 = tda_info['fdr_thres']
    # plt.axvline(tda_fdr1, linestyle='--')
    # plt.text(tda_fdr1, 0.003, r'$\leftarrow$ TDA 1% FDR threshold')
    # plt.title(
    #     f"({dataset.mat.shape[1] / 1000:.1f}k) {dataset_name} {sls} ll={ll:.05f}"
    #     f" constraints={get_cons_str(settings[config]['constraints'])}"
    #     f" {'Y' if model.cons_satisfied else 'N'}")
    # plt.savefig(res_dir + '_'.join(map(str, sls.values())) + '.png')
    return {
        'll':       ll,
        'lls':      lls,
        'sls':      sls,
        'slls':     model.slls,
        'model':    model.frozen(),
        'cons_sat': model.cons_satisfied,
    }


def run_rand_models(n, sls, dataset_name, dataset, tda_info, res_dir):
    rand_dir = f"{res_dir}/random_{'_'.join(map(str, sls.values()))}/"
    if not os.path.exists(rand_dir):
        os.makedirs(rand_dir)
    models_pickle = f'{rand_dir}/models.pickle'
    if os.path.exists(models_pickle):
        models = pickle.load(open(models_pickle, 'rb'))
        print(models_pickle, 'loaded')
    else:
        if parallel and inner_parallel:
            with multiprocessing.get_context('spawn').Pool(num_workers) as pool:
                models = pool.starmap(run_model,
                                      [(sls, dataset_name, dataset, tda_info, rand_dir, i) for i in range(n)])
        else:
            models = list(starmap(run_model, [(sls, dataset_name, dataset, tda_info, rand_dir, i) for i in range(n)]))
        print('sorting models')
        models = list(sorted(models, key=lambda x: x['ll'], reverse=True))
        if not os.path.exists(rand_dir):
            os.makedirs(rand_dir)
        print('dumping models')
        pickle.dump(models, open(models_pickle, 'wb'))

    print('plotting models')
    fig = plt.figure(figsize=(16, 9))
    for i, model in enumerate(models):
        fname = f'rank_{i + 1}.png'
        MixtureModelBase._plot(dataset.mat.T, model['lls'], model['slls'], model['model'], fig)
        ax = plt.gcf().axes[0]
        plt.axes(ax)
        tda_fdr1 = tda_info['fdr_thres']
        plt.axvline(tda_fdr1, linestyle='--')
        plt.text(tda_fdr1, 0.003, r'$\leftarrow$ TDA 1% FDR threshold')
        plt.title(
            f"({dataset.mat.shape[1] / 1000:.1f}k) {dataset_name} {model['sls']} ll={model['ll']:.05f}"
            f" constraints={get_cons_str(settings[config]['constraints'])}"
            f" {'Y' if model['cons_sat'] else 'N'}")
        plt.savefig(f'{rand_dir}/{fname}')
        break
    shutil.copyfile(f'{rand_dir}/rank_1.png', res_dir + '_'.join(map(str, sls.values())) + '.png')

    """ Plot ll hist """
    print('plotting ll hist')
    lls = [model['ll'] for model in models]
    plt.figure()
    plt.hist(lls)
    plt.title(mu_strategy)
    plt.savefig(f'{rand_dir}/ll_hist.png')
    return models[0]


def run_bootstrap_i(dataset_name, bootstrap_i):
    """ Bootstrap for a dataset with the same configuration

    """
    # dataset = XLMS_Dataset(dataset_to_run)
    dataset = XLMS_Dataset(dataset_name, nodup=True)

    bs_indices = pickle.load(open(f'{bootstrap_dir}/index.pickle', 'rb'))

    bs_index = bs_indices[dataset_name][bootstrap_i]

    dataset.mat = dataset.mat[:, bs_index]

    global base_figure_dir
    base_dir = f'{bootstrap_dir}/diffsign_{model_class}_{config}{dir_suffix}' \
               f'_initskew_{"_".join([f"{alpha_base:.0f}" for alpha_base in alpha_bases])}'

    res_dir = f'{base_dir}/{dataset_name}/{bootstrap_i}/'
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)
    tda_info = json.load(open(f'../results/info/{dataset_name}.json'))

    dataset.mat = dataset.mat[:, dataset.mat[1, :] != 0]

    choices = sum([[
        # {'C': alpha_base, 'IC': alpha_base, 'IC2': alpha_base, 'I1': -alpha_base, 'I2': -alpha_base},
        # {'C': alpha_base, 'IC': alpha_base, 'IC2': alpha_base, 'I1': alpha_base, 'I2': alpha_base},
        {'C': alpha_base, 'IC': alpha_base, 'IC2': -alpha_base, 'I1': -alpha_base, 'I2': -alpha_base},
        {'C': alpha_base, 'IC': -alpha_base, 'IC2': -alpha_base, 'I1': -alpha_base, 'I2': -alpha_base},
    ] for alpha_base in alpha_bases], [])

    if part >= 0:
        choices = choices[part:part + 1]

    args = capture_args(locals())

    if init_strategy == 'random':
        models = list(
            map(lambda sls: run_rand_models(random_size, sls, dataset_name, dataset, tda_info, res_dir), choices))
    else:
        models = list(starmap(run_model, choices, *args))
    lock = acquireLock(f'{res_dir}/models_pickle.lock')
    if os.path.exists(f'{res_dir}/models.pickle'):
        models = pickle.load(open(f'{res_dir}/models.pickle', 'rb')) + models
    pickle.dump(models, open(f'{res_dir}models.pickle', 'wb'))
    releaseLock(lock)

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
    plt.text(tda_fdr1, 0.003, r'$\leftarrow$ TDA 1% FDR threshold')
    plt.title(
        f"({dataset.mat.shape[1] / 1000:.1f}k) {dataset_name} {best['sls']} ll={best['ll']:.05f}"
        f" constraints={get_cons_str(settings[config]['constraints'])}"
        f" {'Y' if best['cons_sat'] else 'N'}")
    plt.savefig(f'{res_dir}/best.png')


if __name__ == '__main__':
    __spec__ = None
    run_bootstrap_i(dataset_to_run, bootstrap_i)


