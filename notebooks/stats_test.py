#%% imports and configurations
import imp
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns
import fstlib
import medicc
import medicc.stats

# %%
#file = 'example4_wgd_pdm_total.tsv'
file = 'example7_pdm_total.tsv'
#file = 'example_2_desc_pdm_total.tsv'
#file = 'OV03-04_descr_pdm_total.tsv'
df = pd.read_csv(os.path.join('../examples/test', file), sep='\t').set_index('sample_id')
D = df.values#[0:19,0:19]

# %% test star topology
imp.reload(medicc.stats)
p_star = medicc.stats.star_topology_test(D)
p_clock = medicc.stats.molecular_clock_test(D)
print("P star: %.4f, p clock: %.4f" % (p_star, p_clock))
# %%
