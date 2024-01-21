import xml.etree.ElementTree as ET
from collections import defaultdict
import pymzml
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import traceback as tb


def get_user_param_value(param):
    if param.attrib['type'] == 'string':
        return param.attrib['value']
    elif param.attrib['type'] == 'int':
        return int(param.attrib['value'])
    elif param.attrib['type'] == 'float':
        return float(param.attrib['value'])
    else:
        assert False, param.attrib['type']


def from_user_param(node):
    r = {}
    for x in node.findall('UserParam'):
        r[x.attrib['name']] = get_user_param_value(x)
    return r


class XMLNode:
    def __init__(self, node):
        self.node = node
        self.d = {}
        self.d.update(node.attrib)
        self.d.update(from_user_param(node))

    def __getattr__(self, attr):
        if attr in self.d:
            return self.d[attr]
        print(self.d)
        return super().__getattr__(attr)

    def __getitem__(self, item):
        return self.d[item]


class XLMS_res:
    def __init__(self, filename=None):
        self.filename = filename
        if filename:
            self.xml = ET.parse(filename)

    @staticmethod
    def from_old(old):
        new = XLMS_res()
        new.filename = old.filename
        new.xml = old.xml
        return new

    def pep_matches(self):
        return self.xml.findall('IdentificationRun/PeptideIdentification')

    def all_spectrum(self):
        r = defaultdict(list)
        for pepid in self.pep_matches():
            # pephits = pepid.findall('PeptideHit')
            # # assert 1 <= len(pephits) <= 2, pephits
            # pephits = list(map(PepHit, pephits))
            r[pepid.attrib['spectrum_reference']].append(pepid)
        return r

    def DataFrame(self):
        res = []

        # chain_cols = ['sequence', 'charge', 'aa_before', 'aa_after', 'start', 'end', 'protein_refs']
        # chain_cols = ['sequence']
        for psm_node in self.pep_matches():
            row = {}
            alpha = XMLNode(psm_node.find('PeptideHit[1]'))
            # row['alpha_sequence'] = alpha.sequence
            psm = XMLNode(psm_node)
            # print(psm.spectrum_reference)
            if 'xl_type' not in psm.d:
                print(psm.spectrum_reference, f'wrong format in file {self.filename}, skipping')
                continue
                # assert False
            row.update(psm.d)
            for k, v in alpha.d.items():
                if k in row:
                    assert v == row[k]
            row.update(alpha.d)
            try:
                if psm.xl_type == 'cross-link':
                    beta = XMLNode(psm_node.find('PeptideHit[2]'))
                elif psm.xl_type == 'mono-link':
                    pass
                elif psm.xl_type == 'loop-link':
                    pass
                else:
                    print(psm.xl_type)
                    assert False
                    break
            except e:
                tb.print_exc()
            res.append(row)
        return pd.DataFrame(res)


def DataFrame(filename):
    xlms_res = XLMS_res(filename)
    return xlms_res.DataFrame()


def apply_pepid_filter(d, f):
    r = {k: [x for x in l if f(x)] for k, l in d.items()}
    return r


def filter_top2(x):
    return XMLNode(x).xl_rank <= 2


def apply_spectrum_filter(d, f):
    r = {k: l for k, l in d.items() if f(l)}
    return r


def filter_high_s2(l):
    if len(l) < 2:
        return False
    p2 = l[1]
    return XMLNode(p2.find('PeptideHit'))['OpenPepXL:score'] * 300 > 175


def filter_all_xl(l):
    for x in l:
        if XMLNode(x).xl_type != 'cross-link':
            return False
    return True
