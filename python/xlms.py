import json
import multiprocessing
import hdf5storage

results_dir = '../results'

def register_results_dir(d):
    global results_dir
    results_dir = d


class XLMS_Dataset:
    def __init__(self, dataset_name, nodup=False):
        suffix = '_nodup' if nodup else ''
        self.data_file = f'{results_dir}/matdata/scoremats/{dataset_name}{suffix}.mat'
        self.mat_vars = hdf5storage.loadmat(self.data_file)
        self.mat = self.mat_vars['mat']
        self.info_file = f'{results_dir}/info/{dataset_name}.json'
        self.info = json.load(open(self.info_file))
        print(self.mat_vars)
        print(self.info.keys())
