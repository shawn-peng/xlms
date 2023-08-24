import multiprocessing
import traceback

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
import subprocess

from mixture_model import MixtureModel
from mixture_model_1S import MixtureModel1S
from order_stat_mixture_model import OrderStatMixtureModel
from order_stat_mixture_model_v2 import OrderStatMixtureModelV2
from xlms import XLMS_Dataset

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

SAMPLE_SIZE = 20000
# ic2_comp = True
ic2_comp = False
tolerance = 1e-8
show_plotting = True
# show_plotting = False
plot_interval = 5
# gaussian_model = True
gaussian_model = False
model_samples = 2

alpha_base = 0. if gaussian_model else 2.

# run_all = True
run_all = False
if run_all:
    show_plotting = False

common_settings_2s = {}

settings = {
    'common':                                   {'tolerance': tolerance, 'show_plotting': show_plotting},
    'no_constraint':                            {'ic2_comp':      ic2_comp, 'tolerance': tolerance,
                                                 'show_plotting': show_plotting,
                                                 'plot_interval': plot_interval, 'constraints': []},
    'weight_constraints':                       {'ic2_comp':      ic2_comp, 'tolerance': tolerance,
                                                 'show_plotting': show_plotting,
                                                 'plot_interval': plot_interval, 'constraints': ['weights']},
    'unweighted_pdf_no_weight_constraints':     {'ic2_comp':      ic2_comp, 'tolerance': tolerance,
                                                 'show_plotting': show_plotting,
                                                 'plot_interval': plot_interval, 'constraints': ['pdf']},
    'unweighted_cdf_no_weight_constraints':     {'ic2_comp':      ic2_comp, 'tolerance': tolerance,
                                                 'show_plotting': show_plotting,
                                                 'plot_interval': plot_interval, 'constraints': ['cdf']},
    'unweighted_pdf_cdf_no_weight_constraints': {'ic2_comp':      ic2_comp, 'tolerance': tolerance,
                                                 'show_plotting': show_plotting,
                                                 'plot_interval': plot_interval, 'constraints': ['pdf', 'cdf']},
    'unweighted_pdf':                           {'ic2_comp':      ic2_comp, 'tolerance': tolerance,
                                                 'show_plotting': show_plotting,
                                                 'plot_interval': plot_interval, 'constraints': ['weights', 'pdf']},
    'unweighted_cdf':                           {'ic2_comp':      ic2_comp, 'tolerance': tolerance,
                                                 'show_plotting': show_plotting,
                                                 'plot_interval': plot_interval, 'constraints': ['weights', 'cdf']},
    'unweighted_pdf_cdf':                       {'ic2_comp':      ic2_comp, 'tolerance': tolerance,
                                                 'show_plotting': show_plotting,
                                                 'plot_interval': plot_interval,
                                                 'constraints':   ['weights', 'pdf', 'cdf']},
    'all_constraints':                          {'ic2_comp':      ic2_comp, 'tolerance': tolerance,
                                                 'show_plotting': show_plotting,
                                                 'plot_interval': plot_interval,
                                                 'constraints':   ['weights', 'pdf', 'cdf', 'weighted_pdf']},
}

# config = 'common'
# config = 'no_constraint'
# config = 'weight_constraints'
# config = 'unweighted_pdf'
# config = 'unweighted_cdf'
config = 'unweighted_pdf_cdf'
# config = 'unweighted_pdf_no_weight_constraints'
# config = 'unweighted_cdf_no_weight_constraints'
# config = 'unweighted_pdf_cdf_no_weight_constraints'

map_cons_str = {
    'weights':      'w',
    'pdf':          'pdf',
    'cdf':          'cdf',
    'weighted_pdf': 'wpdf',
}


def get_cons_str(constraints):
    return ';'.join(map(lambda x: map_cons_str[x], constraints))


model_class = f'{model_samples}S{"g" if gaussian_model else ""}{"2" if ic2_comp else ""}'
base_figure_dir = f'figures_python_{model_class}_{config}_initskew_{alpha_base:.0f}'


# base_figure_dir = f'figures_python_2S_{config}'


def run_dataset(dataset_name):
    def enum_signs(comps):
        def get_signs(n, prefix):
            if not n:
                yield {cname: s for s, cname in zip(prefix, comps)}
                return
            yield from get_signs(n - 1, prefix + (+1,))
            yield from get_signs(n - 1, prefix + (-1,))

        yield from get_signs(len(comps), ())

    # res_dir = f'../figures_python_1S_{config}/{dataset_name}/'
    res_dir = f'../{base_figure_dir}/{dataset_name}/'
    # res_dir = f'../figures_python_order_stats_skewnorm_IC_I/{dataset_name}/'
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)
    dataset = XLMS_Dataset(dataset_name)
    dataset.mat = dataset.mat[:, dataset.mat[1, :] != 0]

    n = dataset.mat.shape[1]
    if SAMPLE_SIZE > 0:
        nsample = min(SAMPLE_SIZE, n)
        np.random.seed(42)
        rind = np.random.choice(np.arange(n), nsample, replace=False)
        dataset.mat = dataset.mat[:, rind]

    # choices = [sls for sls in enum_signs(['C', 'IC', 'IC2', 'I1', 'I2']) if
    #            sls['C'] > 0 > sls['I1'] == sls['I2']]
    # models = [{} for _ in range(len(choices))]

    # for i, sls in enumerate(choices):
    def run_model(sls):
        title = f"{dataset_name} constraints={get_cons_str(settings[config]['constraints'])}"
        if model_samples == 1:
            model = MixtureModel1S(sls, **settings[config], title=title)
        elif model_samples == 2:
            model = MixtureModel(sls, **settings[config], title=title)
        # comps = {'C': Norm(), 'IC': SN(1), 'I': SN(-1)}
        comps = {'C': Norm(), 'IC': SN(1), 'I': SN(-1)}
        # comps = {'C': Norm(), 'IC': Norm(), 'I': Norm()}
        # model = OrderStatMixtureModel(**settings['common'])
        # ll, lls = model.fit(dataset.mat.T)
        # model = pickle.load(open('temp_model.pickle', 'rb'))
        # model.initialized = True
        try:
            ll, lls = model.fit(dataset.mat.T)
        except Exception as e:
            print(f'Exception {traceback.format_exc()} happened for {dataset_name}, stopped at middle')
            ll = model.ll
            lls = model.lls
            model.plot(dataset.mat.T, model.lls, model.slls)
        # models[i]['ll'] = ll
        # models[i]['lls'] = lls
        # models[i]['sls'] = sls
        # models[i]['model'] = model
        # plt.subplot(3, 1, 1)
        ax = plt.gcf().axes[0]
        plt.axes(ax)
        plt.title(f"{dataset_name} {sls} ll={ll:.05f} constraints={get_cons_str(settings[config]['constraints'])}")
        plt.savefig(res_dir + '_'.join(map(str, sls.values())) + '.png')
        return {
            'll':    ll,
            'lls':   lls,
            'sls':   sls,
            'model': model,
        }

    choices = [{'C': alpha_base, 'IC': alpha_base, 'IC2': alpha_base, 'I1': -alpha_base, 'I2': -alpha_base}]
    # choices = [{'C': 0, 'IC': 0, 'IC2': 0, 'I1': 0, 'I2': 0}]
    # pool = multiprocessing.Pool(32)
    models = list(map(run_model, choices))
    pickle.dump(models, open(f'{res_dir}models.pickle', 'wb'))

    best = max(models, key=lambda x: x['ll'])
    if best['ll'] == -np.inf:
        print('no solution for', dataset_name)
        return
    best['model'].plot(dataset.mat.T, best['lls'], best['model'].sep_log_likelihood(dataset.mat.T))
    # plt.subplot(3, 1, 1)
    ax = plt.gcf().axes[0]
    plt.axes(ax)
    plt.title(
        f"{dataset_name} {best['sls']} ll={best['ll']:.05f} constraints={get_cons_str(settings[config]['constraints'])}")
    plt.savefig(f'{res_dir}/best.png')


# run_dataset('alban')
# for dataset_name in datasets:
#     run_dataset(dataset_name)
# pool = multiprocessing.Pool(10)
if __name__ == '__main__':
    if run_all:
        multiprocessing.set_start_method('spawn')
        with multiprocessing.Pool(10) as pool:
            res = list(pool.map(run_dataset, datasets))
    else:
        # run_dataset('KKT4')
        run_dataset('peplib')
        # run_dataset('Alinden')
        # run_dataset('ALott')
        # run_dataset('MS2000225')
    # subprocess.call(['C:\\cygwin64\\bin\\bash.exe', '-l', '../copy_best_figures.sh', base_figure_dir],
    #                 cwd='../')

    print(base_figure_dir)
    os.system(f'cd .. && ./copy_best_figures.sh {base_figure_dir}')

