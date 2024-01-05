import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_result(tab):
    sp_ncol = 5
    sp_nrow = (len(tab.columns)+4) // sp_ncol


    plt.figure(figsize=[4*sp_ncol, 4*sp_nrow])

    i = 0
    for col in (tab.columns):
        dtype = tab[col].dtype
        if dtype not in [np.float64, np.int32, np.int64]:
            continue
        i += 1
        plt.subplot(sp_nrow, sp_ncol, i)
        plt.hist(tab[col], bins=50)
        plt.title(col)


def plot_idxml_scores(tab, score_col='OpenPepXL:score', showdecoy=False):
    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(2, hspace=0)
    axs = gs.subplots(sharex=True, sharey=True)
    rank1 = tab[tab['xl_rank']==1]
    rank2 = tab[tab['xl_rank']==2]
    rank1_target = rank1[~rank1['is decoy']]
    rank1_decoy  = rank1[ rank1['is decoy']]
    if showdecoy:
        axs[0].hist(rank1_target[score_col], 50)
        axs[0].hist(rank1_decoy[score_col], 50)
    else:
        axs[0].hist(rank1[score_col], 50)

    axs[1].hist(rank2[score_col], 50)
    plt.title('Score')

    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(2, hspace=0)
    axs = gs.subplots(sharex=True, sharey=True)
    if showdecoy:
        axs[0].hist(np.log(rank1_target[score_col]), 50)
        axs[0].hist(np.log(rank1_decoy[score_col]), 50)
    else:
        axs[0].hist(np.log(rank1[score_col]), 50)

    axs[1].hist(np.log(tab[score_col][(tab['xl_rank']==2) & (tab[score_col]!=0)]), 50)
    plt.title('log(Score)')

    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(2, hspace=0)
    axs = gs.subplots(sharex=True, sharey=True)
#     n = tab['pos2'] - tab['pos1']
    n = tab['end'].apply(lambda x: x[0]) - tab['start'].apply(lambda x: x[0])
    axs[0].hist(rank1[score_col] / n, 50)
    axs[1].hist(rank2[score_col] / n, 50)
    plt.title('Score/PepL')

    fig = plt.figure(figsize=(8,6))
#     s1 = tab[score_col][tab['xl_rank']==1].set_index['spectrum_reference']
#     s2 = tab[score_col][tab['xl_rank']==2].set_index['spectrum_reference']

    fig = plt.figure(figsize=(8,6))
    smat = extract_idxml_smat(tab)
    smat = smat[smat['s2']!=0]
    plt.scatter(smat['s1'], smat['s2'])
    plt.xlabel('s1')
    plt.ylabel('s2')

    fig = plt.figure(figsize=(8,6))
    smat = extract_idxml_smat(tab, noisotope=True)
    smat = smat[smat['s2']!=0]
    plt.scatter(smat['s1'], smat['s2'])
    plt.xlabel('s1')
    plt.ylabel('s2')
#     plt.plot(s1, s2)

def plot_mzid_scores(tab, score_col='OpenPepXL:score'):
    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(2, hspace=0)
    axs = gs.subplots(sharex=True, sharey=True)
    rank1 = tab[tab['xl_rank']==1]
    rank2 = tab[tab['xl_rank']==2]
    rank1_target = rank1[~rank1['is decoy']]
    rank1_decoy  = rank1[ rank1['is decoy']]

    axs[0].hist(tab[score_col][tab['rank']==1], 50)
    axs[1].hist(tab[score_col][tab['rank']==2], 50)
    plt.title('Score')

    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(2, hspace=0)
    axs = gs.subplots(sharex=True, sharey=True)
    axs[0].hist(np.log(tab[score_col][(tab['rank']==1) & (tab[score_col]!=0)]), 50)
    axs[1].hist(np.log(tab[score_col][(tab['rank']==2) & (tab[score_col]!=0)]), 50)
    plt.title('log(Score)')

    fig = plt.figure(figsize=(8,6))
    gs = fig.add_gridspec(2, hspace=0)
    axs = gs.subplots(sharex=True, sharey=True)
    n = tab['end'].apply(lambda x: x[0]) - tab['start'].apply(lambda x: x[0])
    axs[0].hist(tab[score_col][tab['rank']==1] / n, 50)
    axs[1].hist(tab[score_col][tab['rank']==2] / n, 50)
    plt.title('Score/PepL')


    fig = plt.figure(figsize=(8,6))
    smat = extract_mzid_smat(tab)
    smat = smat[smat['s2']!=0]
    plt.scatter(smat['s1'], smat['s2'])
    plt.xlabel('s1')
    plt.ylabel('s2')

    fig = plt.figure(figsize=(8,6))
    smat = extract_mzid_smat(tab, noisotope=True)
    smat = smat[smat['s2']!=0]
    plt.scatter(smat['s1'], smat['s2'])
    plt.xlabel('s1')
    plt.ylabel('s2')

