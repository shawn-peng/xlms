# Decoy-free FDR Estimation for Cross-linked Peptide Mass Spectrometry

## Instructions

### Clone this repository
```git clone https://github.com/shawn-peng/xlms.git```

### Dependencies
Install pytorch with cuda support.

Install other dependencies,
```sh
pip install -r python/requirements.txt
```

### Run Decoy & No-decoy MS searches using OpenPepXLLF (from OpenMS software).
Install OpenPepXLLF software. Please refer to this install instruction of OpenMS,
https://openms.readthedocs.io/en/latest/openms-applications-and-tools/installation.html

Run the search following OpenPepXLLF document and output search results in idXML format.

For decoy search, run the search algorithm with a target-decoy protein sequence database that has decoy sequences appended.

For no-decoy search, just run with target database.

If you want to use scripts provided to parse the results, please put the search results under directory
```results/openpepxllf/knime4.6/<dataset_name>/```,
and for decoy search results, put them under ```results/openpepxllf/knime4.6/<dataset_name>_decoy/```

### Extract Matching Information
Use the following code to extract search results
```Python 3
from extract_search_result import *

dataset_name = '<dataset_name>'

extract_res(dataset_name)
extract_to_matfile_for_dataset(dataset_name)
extract_decoy_res(dataset_name)
extract_TDA_info_for_dataset(dataset_name)
```

The scores matrix for top two matches will be put in ```results/matdata/scoremats```


### Run the mixture model to fit to data

Go into ```python/``` directory

Run the algorithm with default settings using
```python run.py --dataset <dataset_name>```

### Advanced command line options
```
usage: run.py [-h] [-d DATASET] [-c CONFIG] [-s MODEL_SAMPLES] [--suffix SUFFIX] [-j JOBS] [-r RANDOM_SIZE] [-q PART] [-o RANDOM_I] [-t TOLERANCE] [-a] [-p] [-i] [--show_plotting]
              [--clear_results] [--mu_strategy MU_STRATEGY]

options:
  -h, --help            show this help message and exit
  -d DATASET, --dataset DATASET
                        the dataset name to estimate FDR
  -c CONFIG, --config CONFIG
                        constrains setting, [no_constraint, unweighted_pdf_mode]
  -s MODEL_SAMPLES, --model_samples MODEL_SAMPLES
                        number of samples of the mixture model, 1 or 2
  --suffix SUFFIX       add a suffix to the output directory to avoid overwrite previous results
  -j JOBS, --jobs JOBS  number of threads
  -r RANDOM_SIZE, --random_size RANDOM_SIZE
                        number of random starts for each skewness setting
  -q PART, --part PART  the index of skewness setting to run, [0 - 5]
  -o RANDOM_I, --random_i RANDOM_I
                        the index of random start to run
  -t TOLERANCE, --tolerance TOLERANCE
                        the tolerance to check the convergence of the algorithm
  -a, --all             run all the datasets mentioned in the paper
  -p, --parallel        run algorithm in parallel for different starts
  -i, --inner_parallel  run algorithm in parallel for random starts
  --show_plotting       show plotting while fitting the model (this will slow down the process largely)
  --mu_strategy MU_STRATEGY
                        the strategy for generating the initial mu parameters for components, [gaussian, gamma, uniform]

```

