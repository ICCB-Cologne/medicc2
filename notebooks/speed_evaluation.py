#%% imports and configurations
import os
import sys
import time
import numpy as np
import pandas as pd
from plotnine.themes.elements import element_blank
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import plotnine as gg
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
%matplotlib inline
import fstlib
import medicc

#%% read fsts 
T_no_wgd_asymm = fstlib.read('../objects/no_wgd_asymm.fst')
T_no_wgd_symm = ~T_no_wgd_asymm * T_no_wgd_asymm
T_wgd_asymm = fstlib.read('../objects/wgd_asymm.fst')
T_wgd_symm = ~T_wgd_asymm * T_wgd_asymm
symbol_table = T_wgd_asymm.input_symbols()

# %% time composition function and create random sequneces
sequence_len = 20
seq1 = '1' * sequence_len
seq2, nloss, ngains, nwgd = medicc.sim.evolve(seq1, mu=20, pgain=0.8, pwgd=0)

def time_distance(fst, s1, s2, replicates=10, max_length_factor=10, type='legacy'):
    elapsed = np.zeros((replicates, max_length_factor))
    distances = np.zeros((replicates, max_length_factor))
    fsas1 = [fstlib.factory.from_string(s1 * k, isymbols=fst.input_symbols(), osymbols=fst.output_symbols(), arc_type=fstlib.Semiring.TROPICAL) for k in range(1, max_length_factor+1)]
    fsas2 = [fstlib.factory.from_string(s2 * k, isymbols=fst.input_symbols(), osymbols=fst.output_symbols(), arc_type=fstlib.Semiring.TROPICAL) for k in range(1, max_length_factor+1)]
    for j in range(replicates):
        for i in range(max_length_factor):
            t1 = fsas1[i]
            t2 = fsas2[i]
            start = time.time()
            if type=='legacy':
                c1 = fstlib.compose(t1.arcsort('olabel'), fst)
                c2 = fstlib.compose(c1, t2.arcsort('ilabel'))
                d = float(fstlib.shortestdistance(c2, reverse=True)[c2.start()])
            elif type == 'score':
                d = float(fstlib.score(fst, t1, t2))
            elif type == 'kernel':
                d = float(fstlib.kernel_score(fst, t1, t2))
            else:
                d = np.nan
            end = time.time()
            elapsed[j, i] = end-start
            distances[j, i] = d

    return elapsed, distances

#%% run timing experiment
nsamp = 10
mlf = 10
elapsed_S_no_wgd_legacy = time_distance(T_no_wgd_symm, seq1, seq2, nsamp, mlf, type='legacy')
elapsed_S_no_wgd_score = time_distance(T_no_wgd_symm, seq1, seq2, nsamp, mlf, type='score')
elapsed_S_wgd_score = time_distance(T_wgd_symm, seq1, seq2, nsamp, mlf, type='score')
elapsed_T_wgd_kernel = time_distance(T_wgd_asymm, seq1, seq2, nsamp, mlf, type='kernel')
elapsed_T_no_wgd_kernel = time_distance(T_no_wgd_asymm, seq1, seq2, nsamp, mlf, type='kernel')

df_S_no_wgd_legacy = pd.DataFrame(elapsed_S_no_wgd_legacy[0])
df_S_no_wgd_score = pd.DataFrame(elapsed_S_no_wgd_score[0])
df_S_wgd_score = pd.DataFrame(elapsed_S_wgd_score[0])
df_T_wgd_kernel = pd.DataFrame(elapsed_T_wgd_kernel[0])
df_T_no_wgd_kernel = pd.DataFrame(elapsed_T_no_wgd_kernel[0])

#%% plot
timing = pd.concat([df_S_no_wgd_legacy, df_S_no_wgd_score, df_S_wgd_score, df_T_wgd_kernel, df_T_no_wgd_kernel], 
    keys=['no_wgd_legacy', 'no_wgd_lazy', 'wgd_lazy', 'wgd_kernel', 'no_wgd_kernel'],
    names=['tag', 'rep'])
timing.columns = (timing.columns+1) * sequence_len
timing.columns.name = 'length'
timing = timing.stack().to_frame('Time').reset_index('length')
timing['Algorithm'] = ''
timing['WGD'] = ''
timing.loc['no_wgd_legacy', 'Algorithm'] = 'Legacy composition\n(Schwarz et al. 2014)'
timing.loc['no_wgd_lazy', 'Algorithm'] = 'Lazy composition'
timing.loc['wgd_lazy', 'Algorithm'] = 'Lazy composition'
timing.loc['wgd_kernel', 'Algorithm'] = 'Lazy kernel composition'
timing.loc['no_wgd_kernel', 'Algorithm'] = 'Lazy kernel composition'
timing.loc['no_wgd_legacy', 'WGD'] = 'No WGD'
timing.loc['no_wgd_lazy', 'WGD'] = 'No WGD'
timing.loc['wgd_lazy', 'WGD'] = 'WGD'
timing.loc['wgd_kernel', 'WGD'] = 'WGD'
timing.loc['no_wgd_kernel', 'WGD'] = 'No WGD'
scale = 0.8

#%% Joint panel
fig, ax = plt.subplots(figsize=(8*scale, 6*scale))
sns.lineplot(data=timing, x='length', y='Time', hue='Algorithm', style='WGD', 
    palette=sns.color_palette('Set1', 3), markers=True, legend=True, ci='sd', ax=ax)
ax.set_xlabel('CNP length (#segments)')
ax.set_ylabel('Time (s)')
fig.show()
fig.savefig('figures/speed_joint.pdf', bbox_inches='tight')

#%% WGD panel
fig, ax = plt.subplots(figsize=(8*scale, 6*scale))
sns.lineplot(data=timing.query("WGD=='WGD'"), x='length', y='Time', hue='Algorithm', 
    palette=sns.color_palette(None, 3)[1:], style=1, markers=True, legend=False, ax=ax)
ax.set_xlabel('CNP length (#segments)')
ax.set_ylabel('Time (s)')
ax.set_title('WGD')
ax.legend(handles=ax.lines, labels=["Lazy composition","Lazy kernel composition"])
ax.set_ylim([-0.05, timing.Time.max()+0.05])
fig.show()
fig.savefig('figures/speed_wgd.pdf', bbox_inches='tight')

#%% No WGD panel
fig, ax = plt.subplots(figsize=(8*scale, 6*scale))
sns.lineplot(data=timing.query("WGD=='No WGD'"), x='length', y='Time', hue='Algorithm', style=1, markers=True, legend=False, ax=ax)
ax.set_xlabel('CNP length (#segments)')
ax.set_ylabel('Time (s)')
ax.set_title('No WGD')
ax.legend(handles=ax.lines, labels=["Legacy composition", "Lazy composition","Lazy kernel composition"])
ax.set_ylim([-0.05, timing.Time.max()+0.05])
fig.show()
fig.savefig('figures/speed_no_wgd.pdf', bbox_inches='tight')

# %%
t1 = fstlib.factory.from_string(seq1 * 10, isymbols=symbol_table, osymbols=symbol_table, arc_type=fstlib.Semiring.TROPICAL)
t2 = fstlib.factory.from_string(seq2 * 10, isymbols=symbol_table, osymbols=symbol_table, arc_type=fstlib.Semiring.TROPICAL)

#%% SYMM classic
%%timeit
c1 = fstlib.compose(t1.arcsort('olabel'), T_wgd_symm)
c2 = fstlib.compose(c1, t2.arcsort('ilabel'))
d1 = fstlib.shortestdistance(c2,reverse=True)[c2.start()]

#%% SYMM lazy
%%timeit
d2 = fstlib.score(T_wgd_symm, t1, t2)

# %% ASYMM kernel classic
%%timeit
c1 = fstlib.compose(T_wgd_asymm, t1.arcsort('ilabel'))
c2 = fstlib.compose(T_wgd_asymm, t2.arcsort('ilabel')).arcsort('ilabel')
inter = fstlib.compose(~c1, c2)
d3 = fstlib.shortestdistance(inter, reverse=True)[inter.start()]

#%% ASYMM kernel lazy
%%timeit
d4 = fstlib.kernel_score(T_wgd_asymm, t1, t2)

# %%
