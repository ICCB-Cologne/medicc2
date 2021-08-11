#%% Import modules
import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from plotting_params import plotting_params, set_plotting_params

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import fstlib
import medicc

SEED = 42
set_plotting_params()

CALC_NEW = False
DATA_FILE = 'data/Fig_2A.tsv'

#%% read fsts
T_no_wgd_asymm = medicc.io.read_fst('../objects/no_wgd_asymm.fst')
T_no_wgd_symm = ~T_no_wgd_asymm * T_no_wgd_asymm
T_wgd_asymm = medicc.io.read_fst('../objects/wgd_asymm.fst')


# %% time composition function
def time_distance(fst, s1, s2, n_samples=10, max_length_factor=10, type='legacy'):

    elapsed = np.zeros((n_samples, max_length_factor))
    distances = np.zeros((n_samples, max_length_factor))
    fsas1 = [fstlib.factory.from_string(s1 * k, isymbols=fst.input_symbols(), osymbols=fst.output_symbols(
    ), arc_type=fstlib.Semiring.TROPICAL) for k in range(1, max_length_factor+1)]
    fsas2 = [fstlib.factory.from_string(s2 * k, isymbols=fst.input_symbols(), osymbols=fst.output_symbols(
    ), arc_type=fstlib.Semiring.TROPICAL) for k in range(1, max_length_factor+1)]
    for j in range(n_samples):
        for i in range(max_length_factor):
            t1 = fsas1[i]
            t2 = fsas2[i]
            start = time.time()
            if type == 'legacy':
                c1 = fstlib.compose(t1.arcsort('olabel'), fst)
                c2 = fstlib.compose(c1, t2.arcsort('ilabel'))
                d = float(fstlib.shortestdistance(c2, reverse=True)[c2.start()])
            elif type == 'kernel':
                d = float(fstlib.kernel_score(fst, t1, t2))
            else:
                d = np.nan
            end = time.time()
            elapsed[j, i] = end-start
            distances[j, i] = d

    return elapsed, distances


#%% 
if CALC_NEW or not os.path.isfile(DATA_FILE):
    print('Calculating data')
    # Simulate
    sequence_len = 20
    diploid = '1' * sequence_len
    sequence, nloss, ngains, nwgd = medicc.sim.evolve(diploid, mu=20, pgain=0.8, pwgd=0, seed=SEED)

    # run timing experiment
    n_samples = 50
    mlf = 10
    elapsed_S_no_wgd_legacy = time_distance(T_no_wgd_symm, diploid, sequence, n_samples, mlf, type='legacy')
    elapsed_T_wgd_kernel = time_distance(T_wgd_asymm, diploid, sequence, n_samples, mlf, type='kernel')
    elapsed_T_no_wgd_kernel = time_distance(T_no_wgd_asymm, diploid, sequence, n_samples, mlf, type='kernel')

    df_S_no_wgd_legacy = pd.DataFrame(elapsed_S_no_wgd_legacy[0])
    df_T_wgd_kernel = pd.DataFrame(elapsed_T_wgd_kernel[0])
    df_T_no_wgd_kernel = pd.DataFrame(elapsed_T_no_wgd_kernel[0])

    # Create DF
    timing = pd.concat([df_S_no_wgd_legacy, df_T_wgd_kernel, df_T_no_wgd_kernel],
                    keys=['no_wgd_legacy', 'wgd_kernel', 'no_wgd_kernel'],
                    names=['tag', 'rep'])
    timing.columns = (timing.columns+1) * sequence_len
    timing.columns.name = 'length'
    timing = timing.stack().to_frame('Time').reset_index('length')
    timing['Algorithm'] = ''
    timing['WGD'] = ''
    timing.loc['no_wgd_legacy', 'Algorithm'] = 'Legacy composition\n(Schwarz et al. 2014)'
    timing.loc['wgd_kernel', 'Algorithm'] = 'Lazy composition'
    timing.loc['no_wgd_kernel', 'Algorithm'] = 'Lazy composition'
    timing.loc['no_wgd_legacy', 'WGD'] = 'No WGD'
    timing.loc['wgd_kernel', 'WGD'] = 'WGD'
    timing.loc['no_wgd_kernel', 'WGD'] = 'No WGD'

    timing.to_csv(DATA_FILE, sep='\t')
else:
    print('Loading data')
    timing = pd.read_csv(DATA_FILE, sep='\t')

#%% Plot Figure
fig, ax = plt.subplots(figsize=(plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']))
sns.lineplot(data=timing, x='length', y='Time', hue='Algorithm', style='WGD',
             markers=True, legend=True, ci='sd', ax=ax)
ax.set_xlabel('CNP length (#segments)')
ax.set_ylabel('Time (s)')

fig.savefig('final_figures/Fig_2A.pdf', bbox_inches='tight')
fig.savefig('final_figures/Fig_2A.png', bbox_inches='tight', dpi=600)
