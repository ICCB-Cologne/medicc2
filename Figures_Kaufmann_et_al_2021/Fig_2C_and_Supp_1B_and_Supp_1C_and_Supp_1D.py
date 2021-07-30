# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from plotting_params import plotting_params, set_plotting_params

set_plotting_params()

#%% load results
# Due to the number of simulations the results were computed on a computing cluster, see Script TODO
results_grf = pd.read_csv('data/Fig_2C_grf.tsv', sep='\t', index_col=0)
results_rf = pd.read_csv('data/Fig_2C_rf.tsv', sep='\t', index_col=0)
results_quartet = pd.read_csv('data/Fig_2C_quartet.tsv', sep='\t', index_col=0)

results_grf = results_grf.loc[results_grf['Number of Leaves'] > 3]
results_rf = results_rf.loc[results_rf['Number of Leaves'] > 3]
results_quartet = results_quartet.loc[results_quartet['Number of Leaves'] > 3]


#%% Figure 2C
print('Fig 2C')
fig, ax = plt.subplots(figsize=(plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']))

cur_results = results_rf.loc[results_grf['WGD'] == 'Low WGD']
cur_results = cur_results.loc[cur_results['Method'] != 'MEDALT']
cur_results = cur_results.loc[cur_results['Rate'].isin([0.01, 0.025, 0.05])]

nr = cur_results.groupby('Number of Leaves').count().iloc[0, 0]

for i, method in enumerate(['Euclidean Min. Ev.', 'Euclidean NJ', 'Manhattan Min. Ev.', 'Manhattan NJ', 'MEDICC2']):
    cur_method_results = cur_results.loc[cur_results['Method'] == method]
    n_leaves = np.sort(np.unique(cur_method_results['Number of Leaves']))
    ax.plot(n_leaves, cur_method_results.groupby('Number of Leaves')['Distance'].mean().values)

    ax.errorbar(x=n_leaves, 
                y=cur_method_results.groupby('Number of Leaves')['Distance'].mean().values.astype(float),
    #             yerr=cur_method_results.groupby('Number of Leaves')['Distance'].std() / np.sqrt(nr),
                capsize=5, capthick=2, ms=5, marker='o', label=method, color=sns.color_palette('husl', 6)[i])
ax.set_xlabel('Number of leaves')
ax.set_ylabel('Generalized RF distance')

plt.legend(loc='upper right')

fig.savefig('final_figures/Fig_2C.pdf', bbox_inches='tight')
fig.savefig('final_figures/Fig_2C.png', bbox_inches='tight', dpi=600)
#%% Figure Supp 1B Full ranges
print('Supp 1B')
fig, axs = plt.subplots(ncols=3, nrows=2, sharey=True,
                        figsize=(plotting_params['WIDTH_FULL'], plotting_params['WIDTH_HALF']))

cur_results = results_grf.loc[results_grf['Method'] != 'MEDALT']
cur_results = cur_results.loc[cur_results['Number of Leaves'].isin([5, 10, 15, 20])]

for wgd, ax in zip(['No WGD', 'Low WGD', 'High WGD'], axs.T):
    cur_wgd_results = cur_results.loc[cur_results['WGD'] == wgd]

    n_leaves = np.sort(np.unique(cur_wgd_results['Number of Leaves']))
    rates = np.sort(np.unique(cur_wgd_results['Rate']))

    nr = cur_wgd_results.groupby('Number of Leaves').count().iloc[0, 0]

    for i, method in enumerate(['Euclidean Min. Ev.', 'Euclidean NJ', 'Manhattan Min. Ev.', 'Manhattan NJ', 'MEDICC2']):
        cur_method_results = cur_wgd_results.loc[cur_wgd_results['Method'] == method]

        ax[0].errorbar(x=n_leaves, y=cur_method_results.groupby('Number of Leaves')['Distance'].mean(),
                       yerr=cur_method_results.groupby('Number of Leaves')['Distance'].std() / np.sqrt(nr),
                       capsize=5, capthick=2, ms=5, marker='o', label=method, color=sns.color_palette('husl', 6)[i])
        ax[1].errorbar(x=rates, y=cur_method_results.groupby('Rate')['Distance'].mean(),
                       yerr=cur_method_results.groupby('Rate')['Distance'].std() / np.sqrt(nr),
                       capsize=5, capthick=2, ms=5, marker='o', label=method, color=sns.color_palette('husl', 6)[i])

axs[0, 0].set_title('No WGD', fontsize=plotting_params['FONTSIZE_MEDIUM'])
axs[0, 1].set_title('Low WGD', fontsize=plotting_params['FONTSIZE_MEDIUM'])
axs[0, 2].set_title('High WGD', fontsize=plotting_params['FONTSIZE_MEDIUM'])
axs[0, 1].set_xlabel('Number of leaves')
axs[1, 1].set_xlabel('Mutation Rate')
axs[0, 0].set_ylabel('Generalized RF')
axs[1, 0].set_ylabel('Generalized RF')

axs[0, 2].legend(bbox_to_anchor=(1, 1.))
# axs[0, 0].legend(loc='upper left')
plt.tight_layout()

fig.savefig('final_figures/Supp_1B.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_1B.png', bbox_inches='tight', dpi=600)


#%% Figure Supp 1C (Include MEDALT)
print('Supp 1C')
fig, ax = plt.subplots(figsize=(plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']))

cur_results = results_rf.loc[results_grf['WGD'] == 'Low WGD']
cur_results = cur_results.loc[cur_results['Rate'].isin([0.01, 0.025, 0.05])]

nr = cur_results.groupby('Number of Leaves').count().iloc[0, 0]

for i, method in enumerate(['Euclidean Min. Ev.', 'Euclidean NJ', 'Manhattan Min. Ev.', 'Manhattan NJ', 'MEDICC2', 'MEDALT']):
    cur_method_results = cur_results.loc[cur_results['Method'] == method]
    n_leaves = np.sort(np.unique(cur_method_results['Number of Leaves']))
    ax.plot(n_leaves, cur_method_results.groupby('Number of Leaves')['Distance'].mean().values)

    ax.errorbar(x=n_leaves,
                y=cur_method_results.groupby('Number of Leaves')['Distance'].mean().values.astype(float),
                yerr=cur_method_results.groupby('Number of Leaves')['Distance'].std() / np.sqrt(nr),
                capsize=5, capthick=2, ms=5, marker='o', label=method, color=sns.color_palette('husl', 6)[i])
ax.set_xlabel('Number of leaves')
ax.set_ylabel('Generalized RF distance')

plt.legend(loc='upper right')
fig.savefig('final_figures/Supp_1C.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_1C.png', bbox_inches='tight', dpi=600)

#%% Figure Supp 1D (Quartet and RF)
print('Supp 1D')
fig, axs = plt.subplots(ncols=2, figsize=(2*plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']))


for ax, results, name in zip(axs, [results_quartet, results_rf], ['Quartet', 'Robinson-Foulds']):
    cur_results = results.loc[results['WGD'] == 'Low WGD']
    cur_results = cur_results.loc[cur_results['Method'] != 'MEDALT']
    cur_results = cur_results.loc[cur_results['Rate'].isin([0.01, 0.025, 0.05])]
    nr = cur_results.groupby('Number of Leaves').count().iloc[0, 0]

    for i, method in enumerate(['Euclidean Min. Ev.', 'Euclidean NJ', 'Manhattan Min. Ev.', 'Manhattan NJ', 'MEDICC2']):
        cur_method_results = cur_results.loc[cur_results['Method'] == method]
        
        ax.errorbar(x=n_leaves, y=cur_method_results.groupby('Number of Leaves')['Distance'].mean(),
                    yerr=cur_method_results.groupby('Number of Leaves')['Distance'].std() / np.sqrt(nr),
                    capsize=5, capthick=2, ms=5, marker='o', label=method, color=sns.color_palette('husl', 6)[i])
    ax.set_xlabel('Number of leaves')
    ax.set_ylabel('{} distance'.format(name))

axs[1].legend(loc='upper right')

fig.savefig('final_figures/Supp_1D.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_1D.png', bbox_inches='tight', dpi=600)
