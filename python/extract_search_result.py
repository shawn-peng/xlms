import os
import glob

import h5py
import hdf5storage

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import json
from collections import defaultdict

from load_search_result import load_idxmls, load_mzids
from plot_search_result import plot_result

info_dir = './results/info/'


def extract_res(name, suffix=''):
    data_tab = load_idxmls(f'{name}/*.idXML', 'results/openpepxllf/knime4.6/')
    if not os.path.exists(info_dir):
        os.makedirs(info_dir)
    data_tab.to_csv(f'{info_dir}/{name}_res.csv')


def extract_decoy_res(name):
    data_tab = load_idxmls(f'{name}_decoy/*.idXML', 'results/openpepxllf/knime4.6/')
    if not os.path.exists(info_dir):
        os.makedirs(info_dir)
    data_tab.to_csv(f'{info_dir}/{name}_decoy_res.csv')


def load_res(name):
    res_csv = f'{info_dir}/{name}_res.csv'
    res_csv = csv.DictReader(open(res_csv))
    # data_tab = pd.read_csv(res_csv, index_col=0)
    data_tab = pd.DataFrame(res_csv)
    data_tab = data_tab.set_index('')
    data_tab['xl_rank'] = data_tab['xl_rank'].astype(int)
    data_tab['OpenPepXL:score'] = data_tab['OpenPepXL:score'].astype(float)
    return data_tab


def load_decoy_res(name):
    res_csv = f'{info_dir}/{name}_decoy_res.csv'
    res_csv = csv.DictReader(open(res_csv))
    # data_tab = pd.read_csv(res_csv, index_col=0)
    data_tab = pd.DataFrame(res_csv)
    data_tab = data_tab.set_index('')
    data_tab['xl_rank'] = data_tab['xl_rank'].astype(int)
    data_tab['OpenPepXL:score'] = data_tab['OpenPepXL:score'].astype(float)
    return data_tab


def extract_mzid_smat(tab: pd.DataFrame, score_column='OpenPepXL:score', noisotope=False):
    spec_matches = defaultdict(list)
    for i, row in tab.iterrows():
        specid = row['spectrumID']
        rank = row['rank']
        s = row[score_column]
        print(specid, rank, row['xl_type'], s)
        #         for k, v in row.items():
        #             print(k, v)
        if not s:
            continue
        if noisotope and spec_matches[specid] and spec_matches[specid][-1] == s:
            continue
        if len(spec_matches[specid]) >= 2:
            continue
        spec_matches[specid].append(row[score_column])
    for l in spec_matches.values():
        if len(l) == 1:
            l.append(0)
    return pd.DataFrame(spec_matches).transpose().rename(columns={0: 's1', 1: 's2'})


def extract_idxml_smat(tab: pd.DataFrame, score_column='OpenPepXL:score', noisotope=False, keep_diff_xl_pos=False):
    spec_matches = defaultdict(list)
    spec_peptides = defaultdict(set)
    for i, row in tab.iterrows():
        specid = row['spectrum_reference']
        rank = row['xl_rank']
        s = row[score_column]
        # print(specid, rank, row['xl_type'], s)
        # for k, v in row.items():
        #     print(k, v)
        if not s:
            continue
        if noisotope and spec_matches[specid] and spec_matches[specid][-1] == s:
            continue
        if len(spec_matches[specid]) >= 2:
            continue
        pep = (row['sequence'], row['sequence_beta'])
        if pep in spec_peptides[specid]:
            # print('skipping same peptides with different xl_pos')
            # print(pep)
            continue
        if not keep_diff_xl_pos:
            spec_peptides[specid].add(pep)
        spec_matches[specid].append(row[score_column])
    for l in spec_matches.values():
        if len(l) == 1:
            l.append(-10000)
    return pd.DataFrame(spec_matches).transpose().rename(columns={0: 's1', 1: 's2'})


def extract_to_matfile(tab, filename, keep_diff_xl_pos=False):
    matfiledata = {}  # make a dictionary to store the MAT data in
    mat = extract_idxml_smat(tab, keep_diff_xl_pos=keep_diff_xl_pos)
    matfiledata[u'mat'] = mat.to_numpy().T
    # *** u prefix for variable name = unicode format, no issues thru Python 3.5;
    # advise keeping u prefix indicator format based on feedback despite docs ***

    #     matfiledata[u'variable2'] = np.ones(300)
    matdata_dir = './results/matdata/scoremats/'
    if not os.path.exists(matdata_dir):
        os.makedirs(matdata_dir)
    hdf5storage.writes(matfiledata, os.path.join(matdata_dir, filename), matlab_compatible=True)


def extract_to_matfile_for_dataset(name):
    info_dir = './results/info/'
    res_csv = f'{info_dir}/{name}_res.csv'
    res_csv = csv.DictReader(open(res_csv))
    # data_tab = pd.read_csv(res_csv, index_col=0)
    data_tab = pd.DataFrame(res_csv)
    data_tab = data_tab.set_index('')
    data_tab['xl_rank'] = data_tab['xl_rank'].astype(int)
    data_tab['OpenPepXL:score'] = data_tab['OpenPepXL:score'].astype(float) * 300

    plot_result(data_tab)
    # extract_to_matfile(data_tab, f'{name}.mat', keep_diff_xl_pos=True)
    extract_to_matfile(data_tab, f'{name}_nodup.mat', keep_diff_xl_pos=False)


def extract_TDA_info(data_tab, output_file):
    data_tab = data_tab[(data_tab['xl_rank'] == 1) | (data_tab['xl_rank'] == '1')]
    data_tab = data_tab.sort_values('OpenPepXL:score', ascending=False)
    data_tab['OpenPepXL:score'] *= 300

    print(len(data_tab))

    a_decoy = data_tab['xl_target_decoy_alpha'] == 'decoy'
    b_decoy = data_tab['xl_target_decoy_beta'] == 'decoy'
    ab_decoy = a_decoy & b_decoy
    decoy = a_decoy | b_decoy
    print(decoy)
    num_dd = (a_decoy & b_decoy).sum()
    num_d = (a_decoy | b_decoy).sum()
    print(num_dd, num_d)

    ndecoy = 0
    ndd = 0
    ntarget = 0
    fdr_thres = 0
    qvals = np.zeros(len(data_tab))
    scores = []
    curve_fdr = []
    curve_decoy = []
    curve_dd = []
    curve_matches = []

    specs = set()

    qstart = 0
    i = 0
    for ind, row in data_tab.iterrows():
        spec = row['spectrum_reference']
        if spec in specs:
            print('dup spec', spec)
            continue
        specs.add(spec)
        scores.append(data_tab['OpenPepXL:score'].iloc[i])
        if decoy.iloc[i]:
            if ntarget == 0:
                qval = 1
            else:
                qval = (ndecoy - 2 * ndd) / ntarget
            qvals[qstart:i] = qval
            qstart = i
            if qval >= 0.01 and fdr_thres == 0:
                fdr_thres = data_tab['OpenPepXL:score'].iloc[i]
            ndecoy += 1
        else:
            ntarget += 1
        if ab_decoy.iloc[i]:
            ndd += 1

        curve_decoy.append(ndecoy)
        curve_dd.append(ndd)
        curve_matches.append(ntarget + ndecoy)
        i += 1

    print('num spec', len(specs))

    if ntarget == 0:
        qval = 1
    else:
        qval = (ndecoy - 2 * ndd) / ntarget
    qvals[qstart:i + 1] = qval
    print(i)
    curve_fdr = list(qvals)

    plt.figure()
    plt.plot(curve_fdr)

    plt.figure()
    plt.plot(curve_dd)
    plt.plot(curve_decoy)
    plt.plot(curve_matches)

    print(len(data_tab['OpenPepXL:score']))

    #     plt.figure()
    #     plt.plot(data_tab['OpenPepXL:score'], curve_fdr)

    #     plt.figure()
    #     plt.plot(data_tab['OpenPepXL:score'], curve_matches)
    #     plt.plot(data_tab['OpenPepXL:score'], curve_dd)
    #     plt.plot(data_tab['OpenPepXL:score'], curve_decoy)
    #     plt.plot(data_tab['OpenPepXL:score'], np.array(curve_decoy) - 2*np.array(curve_dd))

    info = {
        'ndecoy':        ndecoy,
        'ndd':           ndd,
        'ntarget':       ntarget,
        'fdr_thres':     fdr_thres,
        'scores':        scores,
        'curve_fdr':     curve_fdr,
        'curve_decoy':   curve_decoy,
        'curve_dd':      curve_dd,
        'curve_matches': curve_matches,
    }
    json.dump(info, open(output_file, 'w'))
    return info


def extract_TDA_info_for_dataset(name):
    info_dir = './results/info/'
    data_tab = load_res(name)
    # print(data_tab)

    if not os.path.exists(info_dir):
        os.makedirs(info_dir)
    info = extract_TDA_info(data_tab, f'{info_dir}/{name}.json')
    return info
