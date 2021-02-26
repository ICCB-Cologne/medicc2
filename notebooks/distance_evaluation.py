#%% imports and configurations
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import imp
import numpy as np
import pandas as pd
import scipy as sp
import sklearn as skl
import sklearn.metrics
import matplotlib as mpl
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns
import fstlib
import medicc

def fst_dist(row, fst, seq_in='in', seq_out='out', kernel=False):
    fsa_in = fstlib.factory.from_string(row[seq_in], isymbols=fst.input_symbols(), osymbols=fst.output_symbols(), arc_type=fst.arc_type())
    fsa_out = fstlib.factory.from_string(row[seq_out], isymbols=fst.input_symbols(), osymbols=fst.output_symbols(), arc_type=fst.arc_type())
    if kernel:
        dist = float(fstlib.kernel_score(fst, fsa_in, fsa_out))    
    else:
        dist = float(fstlib.score(fst, fsa_in, fsa_out))
    return dist

def euclidean_dist(row, seq_in='in', seq_out='out'):
    chrs_in = row[seq_in].split('X')
    chrs_out = row[seq_out].split('X')
    dist = 0
    for pair in zip(chrs_in, chrs_out):
        vec1 = np.array(list(pair[0])).astype('int')
        vec2 = np.array(list(pair[1])).astype('int')
        dist += np.sum((vec1-vec2)**2)
    dist = np.sqrt(dist)
    return dist

def to_major_minor(row, seq_in1='out', seq_in2='out2', major=True, zero_wins=False):
    str1 = row[seq_in1]
    str2 = row[seq_in2]
    str_out = []
    for pair in zip(str1, str2):
        a, b = pair
        if a != 'X' and b != 'X':        
            cna=medicc.tools.hex2int(a)
            cnb=medicc.tools.hex2int(b)
            if zero_wins and cna==0:
                majorstr = medicc.tools.int2hex(cna)
                minorstr = medicc.tools.int2hex(cnb)
            elif zero_wins and cnb == 0:
                majorstr = medicc.tools.int2hex(cnb)
                minorstr = medicc.tools.int2hex(cna)
            else:
                majorstr = medicc.tools.int2hex(max(cna, cnb))
                minorstr = medicc.tools.int2hex(min(cna, cnb))

            if major:
                str_out.append(majorstr)
            else:
                str_out.append(minorstr)
        else:
            str_out.append('X')
    return "".join(str_out)

#%%
T_wgd_asymm = fstlib.Fst.read('../objects/wgd_asymm.fst')
T_no_wgd_asymm = fstlib.Fst.read('../objects/no_wgd_asymm.fst')

# %% generate samples
chr_length = 10
nchr =5
nsamples = 1000
diploid = "X".join(['1' * chr_length,] * nchr)
rand = pd.DataFrame([(diploid,) + 
    medicc.sim.evolve(diploid, mu=10, verbose=False, mincn=0) + 
    medicc.sim.evolve(diploid, mu=5, verbose=False, mincn=0) for i in range(nsamples)], 
    columns=['in', 'out', 'nloss', 'ngain', 'nwgd', 'out2', 'nloss2', 'ngain2', 'nwgd2'])
rand['is_wgd'] = rand['nwgd']>0
rand['is_wgd2'] = rand['nwgd2']>0
rand['nevents'] = rand['ngain'] + rand['nloss'] + rand['nwgd']
rand['nevents2'] = rand['ngain2'] + rand['nloss2'] + rand['nwgd2']
rand['nevents_tree'] = rand['nevents'] + rand['nevents2']
rand['out_major'] = rand.apply(to_major_minor, major=True, axis=1)
rand['out_minor'] = rand.apply(to_major_minor, major=False, axis=1)

#%% annotate
imp.reload(medicc.stats)
## branch distances left
rand['dist_fst_wgd_asymm'] = rand.apply(fst_dist, fst=T_wgd_asymm, axis=1)
rand['dist_fst_no_wgd_asymm'] = rand.apply(fst_dist, fst=T_no_wgd_asymm, axis=1)
rand['dist_fst_wgd_symm'] = rand.apply(fst_dist, fst=T_wgd_asymm, kernel=True, axis=1)
rand['dist_fst_no_wgd_symm'] = rand.apply(fst_dist, fst=T_no_wgd_asymm, kernel=True, axis=1)
rand['dist_euclidean'] = rand.apply(euclidean_dist, axis=1)

## branch distances right
rand['dist_fst_wgd_asymm2'] = rand.apply(fst_dist, fst=T_wgd_asymm, seq_out='out2', axis=1)
rand['dist_fst_no_wgd_asymm2'] = rand.apply(fst_dist, fst=T_no_wgd_asymm, seq_out='out2', axis=1)
rand['dist_fst_wgd_symm2'] = rand.apply(fst_dist, fst=T_wgd_asymm, seq_out='out2', kernel=True, axis=1)
rand['dist_fst_no_wgd_symm2'] = rand.apply(fst_dist, fst=T_no_wgd_asymm, seq_out='out2', kernel=True, axis=1)
rand['dist_euclidean2'] = rand.apply(euclidean_dist, seq_out='out2', axis=1)

## symmetric distances
rand['tree_dist_fst_wgd_symm'] = rand.apply(fst_dist, fst=T_wgd_asymm, seq_in='out', seq_out='out2', kernel=True, axis=1)
rand['tree_dist_euclidean'] = rand.apply(euclidean_dist, seq_in='out', seq_out='out2', axis=1)

## major/minor branch distances
rand['major_dist_fst_wgd_asymm'] = rand.apply(fst_dist, fst=T_wgd_asymm, seq_in='in', seq_out='out_major', axis=1)
rand['minor_dist_fst_wgd_asymm'] = rand.apply(fst_dist, fst=T_wgd_asymm, seq_in='in', seq_out='out_minor', axis=1)
rand['mm_dist'] = rand['major_dist_fst_wgd_asymm'] + rand['minor_dist_fst_wgd_asymm']

rand['Count'] = 1
#rand['p_wgd'] = rand.apply(lambda x: medicc.stats.wgd_test(x['dist_fst_wgd_asymm'], x['dist_fst_no_wgd_asymm'], 5, 10), axis=1)
#rand.query('~is_wgd').eval('dist_fst_no_wgd_asymm-dist_fst_wgd_asymm').hist(bins=20)
#rand.query('~is_wgd').eval('dist_fst_no_wgd_asymm-dist_fst_wgd_asymm').mean()

#%%
descr = {'dist_fst_wgd_asymm':'MED (WGD)',
    'dist_fst_no_wgd_asymm':'MED (no WGD)',
    'dist_fst_wgd_symm':'Symmetric MED (WGD)',
    'dist_fst_no_wgd_symm':'Symmetric MED (no WGD)',
    'nevents_tree':'True tree distance',
    'tree_dist_fst_wgd_symm':'Symmetric MED (WGD)',
    'dist_euclidean':'Euclidean distance',
    'tree_dist_euclidean':'Euclidean tree distance',
    'nevents': 'True distance'}


# %% final plot
scale = 0.8
distnames = [('nevents', c) for c in rand.columns if c.startswith('dist') and c.endswith('asymm')]
distnames.append(('nevents', 'dist_euclidean'))
distnames.append(('dist_fst_wgd_asymm', 'dist_fst_wgd_symm'))
distnames.append(('nevents_tree', 'tree_dist_fst_wgd_symm'))
distnames.append(('nevents_tree', 'tree_dist_euclidean'))
distnames = np.array(distnames)

def plot_variables(df, variables):
    nplots = len(variables)
    if nplots>=2:
        ncol = 2
    else:
        ncol = 1
    nrow = np.ceil(nplots / ncol)
    fig = plt.figure(figsize=(6*ncol*scale, 6*nrow*scale))
    ax = [None,] * nplots
    for i, names in enumerate(variables):
        xname, yname = names
        ax[i] = fig.add_subplot(nrow, ncol, i+1)
        plotdat = df.groupby([xname, yname, 'is_wgd']).sum()['Count'].reset_index()
        plotdat.rename({'is_wgd':'WGD'}, axis=1, inplace=True)
        sns.scatterplot(x=xname, y=yname, size="Count", sizes=(20,500), alpha=0.5, hue='WGD', data=plotdat, axes=ax[i])
        ax[i].set_xlabel(descr[xname])
        ax[i].set_ylabel(descr[yname])
        ax[i].axline((0,0), slope=1, alpha=0.5, ls='--')
        r2 = sp.stats.pearsonr(df[xname], df[yname])[0]**2
        ax[i].text(0.98, 0.02, r"$r^2 = %.2f$" % r2, transform=ax[i].transAxes, ha='right', va='bottom')
    return fig

with plt.rc_context({"legend.labelspacing": 0.5}):
    fig_main = plot_variables(rand, [distnames[0],])
    fig_main.savefig('figures/distances_main_wgd.pdf', bbox_inches='tight')

    fig_main = plot_variables(rand, [distnames[1],])
    fig_main.savefig('figures/distances_main_no_wgd.pdf', bbox_inches='tight')

    fig_main = plot_variables(rand, [distnames[-1],])
    fig_main.savefig('figures/distances_main_euclid.pdf', bbox_inches='tight')





# %% difference
diff = ((rand['dist_fst_no_wgd_asymm'] + rand['dist_fst_no_wgd_asymm2']) - (rand['dist_fst_wgd_asymm'] + rand['dist_fst_wgd_asymm2']))
rand['diff_total'] = diff
rand['nwgd_total'] = rand['nwgd'] + rand['nwgd2']
rand['is_wgd_total'] = rand['nwgd_total']>0

# %%
sns.displot(x='diff_total', col='is_wgd_total', data=rand)

# %%
fig, ax = plt.subplots(figsize=(8,6))
plotdat = rand.groupby(['dist_fst_wgd_asymm', 'dist_fst_no_wgd_asymm', 'is_wgd']).sum()['Count'].reset_index()
plotdat.rename({'is_wgd':'WGD'}, axis=1, inplace=True)
sns.scatterplot(x='dist_fst_wgd_asymm', y='dist_fst_no_wgd_asymm', size="Count", sizes=(20,500), alpha=0.5, hue='WGD', data=plotdat, ax=ax)
ax.set_xlabel("Minimum event distance (with WGD)")
ax.set_ylabel("Minimum event distance (without WGD)")
fig.show()
# %%
