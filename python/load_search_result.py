import os
import glob

from pyteomics import mzid
# from pyteomics.openms import idxml
import idxml
import pandas as pd


def parse_attrs(s):
    return {k: v for k, v in [attrstr.split('=') for attrstr in s.split()]}


def load_mzids(pattern, directory='./results/openpepxllf/knime4.6/'):
    tab = None
    for f in glob.glob('%s/%s' % (directory, pattern)):
        #         for item in mzid.read(f):
        #             print(item)
        t = mzid.DataFrame(f)

        def get_spec_ref(row):
            ref = parse_attrs(row['spectrumID'])
            return 'file=%s scan=%s' % (os.path.basename(f), ref['scan'])

        t['spectrumID'] = t.apply(get_spec_ref, axis=1)
        if tab is None:
            tab = t.sort_values(['spectrumID', 'rank'])
        else:
            tab = pd.concat((tab, t.sort_values(['spectrumID', 'rank'])), ignore_index=True)
    tab = pd.DataFrame(tab)
    return tab


def load_idxmls(pattern, res_dir='XLSearch/openpepxllf_results'):
    tab = None

    for f in glob.glob(os.path.join(res_dir, pattern)):
        print(f)
        t = idxml.DataFrame(f)

        def get_spec_ref(row):
            ref = parse_attrs(row['spectrum_reference'])
            return 'file=%s scan=%s' % (os.path.basename(f), ref['scan'])

        t['spectrum_reference'] = t.apply(get_spec_ref, axis=1)
        if tab is None:
            tab = t.sort_values(['spectrum_reference', 'xl_rank'])
        else:
            tab = pd.concat((tab, t.sort_values(['spectrum_reference', 'xl_rank'])), ignore_index=True)
    tab = pd.DataFrame(tab)

    return tab
