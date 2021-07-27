#%% imports and configurations
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
import seaborn as sns

from plotting_params import plotting_params, set_plotting_params

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import fstlib
import medicc

set_plotting_params()
SEED = 42
# %matplotlib inline
#%%
def fst_dist(row, fst, seq_in='in', seq_out='out', kernel=False):
    fsa_in = fstlib.factory.from_string(row[seq_in], isymbols=fst.input_symbols(
    ), osymbols=fst.output_symbols(), arc_type=fst.arc_type())
    fsa_out = fstlib.factory.from_string(row[seq_out], isymbols=fst.input_symbols(
    ), osymbols=fst.output_symbols(), arc_type=fst.arc_type())
    if kernel:
        dist = float(fstlib.kernel_score(fst, fsa_in, fsa_out))
    else:
        dist = float(fstlib.score(fst, fsa_in, fsa_out))
    return dist


def plot_variables(df, variables):

    descr = {'dist_fst_wgd_asymm': 'MED (with WGD)',
             'dist_fst_no_wgd_asymm': 'MED (without WGD)',
             'nevents': 'True distance'}

    fig, ax = plt.subplots(figsize=(plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']))

    xname, yname = variables
    plotdat = df.groupby([xname, yname, 'is_wgd']).sum()['Count'].reset_index()
    plotdat.rename({'is_wgd': 'WGD'}, axis=1, inplace=True)
    sns.scatterplot(x=xname, y=yname, size="Count", sizes=(50, 250),
                    alpha=0.5, hue='WGD', data=plotdat, axes=ax)
    ax.set_xlabel(descr[xname])
    ax.set_ylabel(descr[yname])
    ax.axline((0, 0), slope=1, alpha=0.5, ls='--', color='grey')
    r2 = sp.stats.pearsonr(df[xname], df[yname])[0]**2
    ax.text(0.98, 0.02, "$r^2 = {:.2f}$".format(r2), transform=ax.transAxes, 
            ha='right', va='bottom', fontsize=plotting_params['FONTSIZE_MEDIUM'])
    return fig

#%% Load FSTs
T_wgd_asymm = fstlib.Fst.read('../objects/wgd_asymm.fst')
T_no_wgd_asymm = fstlib.Fst.read('../objects/no_wgd_asymm.fst')

# %% generate samples
chr_length = 10
nchr = 5
nsamples = 1000
diploid = "X".join(['1' * chr_length, ] * nchr)
results = pd.DataFrame([(diploid,) +
                     medicc.sim.evolve(diploid, mu=10, verbose=False, mincn=0, seed=i) +
                     medicc.sim.evolve(diploid, mu=5, verbose=False, mincn=0, seed=i) for i in range(nsamples)],
                    columns=['in', 'out', 'nloss', 'ngain', 'nwgd', 'out2', 'nloss2', 'ngain2', 'nwgd2'])
results['nevents'] = results['ngain'] + results['nloss'] + results['nwgd']
results['is_wgd'] = results['nwgd']>0
#%% annotate
results['dist_fst_wgd_asymm'] = results.apply(fst_dist, fst=T_wgd_asymm, axis=1)
results['dist_fst_no_wgd_asymm'] = results.apply(fst_dist, fst=T_no_wgd_asymm, axis=1)
results['Count'] = 1


#%% Create Figures
fig_main = plot_variables(results, ('nevents', 'dist_fst_wgd_asymm'))
fig_main.savefig('final_figures/Fig_2B.pdf', bbox_inches='tight')
fig_main.savefig('final_figures/Fig_2B.png', bbox_inches='tight', dpi=600)

fig_main = plot_variables(results, ('nevents', 'dist_fst_no_wgd_asymm'))
fig_main.savefig('final_figures/Supp_1A.pdf', bbox_inches='tight')
fig_main.savefig('final_figures/Supp_1A.png', bbox_inches='tight', dpi=600)
