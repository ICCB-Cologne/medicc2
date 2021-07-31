# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from plotting_params import plotting_params, set_plotting_params

set_plotting_params()
SEED = 42

#%% load results
# Due to the number of simulations the results were computed on a computing cluster, see Script TODO
results = pd.read_csv('data/Fig_2C_and_Supp_1B_and_Supp_1C.tsv', sep='\t', index_col=False)

#%% Figure 2C
fig, ax = plt.subplots(figsize=(plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']))

cur_results = results.loc[results['WGD'] == 'Low WGD']
cur_results = cur_results.loc[cur_results['Method'] != 'MEDALT']
cur_results = cur_results.loc[cur_results['Rate'].isin([0.01, 0.025, 0.05])]

n_leaves = np.sort(np.unique(cur_results['Number of Leaves']))
nr = cur_results.groupby('Number of Leaves').count().iloc[0, 0]

for i, method in enumerate(['Euclidean Min. Ev.', 'Euclidean NJ', 'Hamming Min. Ev.', 'Hamming NJ', 'MEDICC2']):
    cur_method_results = results.loc[results['Method'] == method]
    ax.errorbar(x=n_leaves, y=cur_method_results.groupby('Number of Leaves')['Distance'].mean(),
                yerr=cur_method_results.groupby('Number of Leaves')['Distance'].std() / np.sqrt(nr),
                capsize=5, capthick=2, ms=5, marker='o', label=method, color=sns.color_palette('Set1', 6)[i])
ax.set_xlabel('Number of leaves')
ax.set_ylabel('Generalized RF distance')

ax.set_yticks([0.0, 0.05, 0.1, 0.15, 0.2, 0.25])
ax.set_ylim(0, 0.25)

plt.legend(loc='lower right')

fig.savefig('final_figures/Fig_2C.pdf', bbox_inches='tight')
fig.savefig('final_figures/Fig_2C.png', bbox_inches='tight', dpi=600)

#%% Figure Supp 1B
fig, axs = plt.subplots(ncols=3, nrows=2, sharey=True,
                        figsize=(plotting_params['WIDTH_FULL'], plotting_params['WIDTH_HALF']))

cur_results = results.loc[results['Method'] != 'MEDALT']
cur_results = cur_results.loc[cur_results['Number of Leaves'].isin([5, 10, 15, 20])]


for wgd, ax in zip(['No WGD', 'Low WGD', 'High WGD'], axs.T):
    cur_wgd_results = cur_results.loc[cur_results['WGD'] == wgd]

    n_leaves = np.sort(np.unique(cur_wgd_results['Number of Leaves']))
    rates = np.sort(np.unique(cur_wgd_results['Rate']))

    nr = cur_wgd_results.groupby('Number of Leaves').count().iloc[0, 0]

    for i, method in enumerate(['Euclidean Min. Ev.', 'Euclidean NJ', 'Hamming Min. Ev.', 'Hamming NJ', 'MEDICC2']):
        cur_method_results = cur_wgd_results.loc[results['Method'] == method]

        ax[0].errorbar(x=n_leaves, y=cur_method_results.groupby('Number of Leaves')['Distance'].mean(),
                       yerr=cur_method_results.groupby('Number of Leaves')[
            'Distance'].std() / np.sqrt(nr),
            capsize=5, capthick=2, ms=5, marker='o', label=method, color=sns.color_palette('Set1', 6)[i])
        ax[1].errorbar(x=rates, y=cur_method_results.groupby('Rate')['Distance'].mean(),
                       yerr=cur_method_results.groupby('Rate')['Distance'].std() / np.sqrt(nr),
                       capsize=5, capthick=2, ms=5, marker='o', label=method, color=sns.color_palette('Set1', 6)[i])

axs[0, 0].set_title('No WGD', fontsize=plotting_params['FONTSIZE_MEDIUM'])
axs[0, 1].set_title('Low WGD', fontsize=plotting_params['FONTSIZE_MEDIUM'])
axs[0, 2].set_title('High WGD', fontsize=plotting_params['FONTSIZE_MEDIUM'])
axs[0, 1].set_xlabel('Number of leaves')
axs[1, 1].set_xlabel('Mutation Rate')
axs[0, 0].set_ylabel('Generalized RF')
axs[1, 0].set_ylabel('Generalized RF')

for ax in axs[:, 0]:
    ax.set_yticks([0.0, 0.05, 0.1, 0.15, 0.2])
    ax.set_ylim(0, 0.225)

axs[0, 2].legend(bbox_to_anchor=(1, 1.))
# axs[0, 0].legend(loc='upper left')
plt.tight_layout()

fig.savefig('final_figures/Supp_1B.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_1B.png', bbox_inches='tight', dpi=600)


#%% Figure Supp 1C
# TODO
# fig, ax = plt.subplots(figsize=(plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']))

# cur_results = results.loc[results['WGD'] == 'Low WGD']
# cur_results = cur_results.loc[cur_results['Rate'].isin([0.01, 0.025, 0.05])]

# n_leaves = np.sort(np.unique(cur_results['Number of Leaves']))
# nr = cur_results.groupby('Number of Leaves').count().iloc[0, 0]

# for i, method in enumerate(['Euclidean Min. Ev.', 'Euclidean NJ', 'Hamming Min. Ev.', 'Hamming NJ', 'MEDICC2', 'MEDALT']):
#     cur_method_results = results.loc[results['Method'] == method]
#     ax.errorbar(x=n_leaves, y=cur_method_results.groupby('Number of Leaves')['Distance'].mean(),
#                 yerr=cur_method_results.groupby('Number of Leaves')['Distance'].std() / np.sqrt(nr),
#                 capsize=5, capthick=2, ms=5, marker='o', label=method, color=sns.color_palette('Set1', 6)[i])
# ax.set_xlabel('Number of leaves')
# ax.set_ylabel('Generalized RF distance')

# ax.set_yticks([0.0, 0.05, 0.1, 0.15, 0.2, 0.25])
# ax.set_ylim(0, 0.25)

# plt.legend(loc='lower right')

# fig.savefig('final_figures/Supp_1C.pdf', bbox_inches='tight')
# fig.savefig('final_figures/Supp_1C.png', bbox_inches='tight', dpi=600)
