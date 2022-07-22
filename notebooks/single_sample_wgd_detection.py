#%% Imports
import os
import sys

import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import medicc

#%% Load data
data_folder = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           "../examples/gundem_et_al_2015")               
input_df = medicc.io.read_and_parse_input_data(os.path.join(data_folder, "PTX011_input_df.tsv"))
samples = input_df.index.get_level_values('sample_id').unique()

#%% calculate WGD status for each sample
results = pd.Series(index=samples, dtype=float)
for sample in samples:
    has_wgd = medicc.core.detect_wgd(input_df, sample)
    results[sample] = has_wgd

results

#%% Bootstrap WGD detection for each sample
N_bootstrap = 100
results_bootstrap = pd.Series(0, index=samples, dtype=float)

for _ in range(N_bootstrap):
    cur_bootstrap_df = medicc.bootstrap.chr_wise_bootstrap_df(input_df)
    for sample in samples:
        has_wgd = medicc.core.detect_wgd(cur_bootstrap_df, sample)
        results_bootstrap[sample] += int(has_wgd)

results_bootstrap
#%% WGD assignment based on a 5% threshold
bootstrap_threshold = 5
results_bootstrap > bootstrap_threshold
