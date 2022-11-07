# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from plotting_params import plotting_params, set_plotting_params

set_plotting_params()

#%% load results
# Due to the number of simulations the results were computed on a computing cluster, see Script TODO
results = pd.read_csv('data/Fig_2B_and_Supp_3.tsv', sep='\t', index_col=0)
random_results = pd.read_csv('data/Fig_2B_and_Supp_3_random_baseline.tsv', sep='\t', index_col=0)
method_rename = {'euclidean_nj': 'Euclidean NJ',
          'manhattan_nj': 'Hamming NJ',
          'euclidean_fastme': 'Euclidean Min. Ev.',
          'manhattan_fastme': 'Hamming Min. Ev.',
          'MEDALT': 'MEDALT',
          'medicc': 'MEDICC2',
          'sitka': 'Sitka'}
results['Method'] = results['Method'].apply(lambda x: method_rename[x])
wgd_rename = {'nowgd': 'No WGD',
          'lowwgd': 'Low WGD',
          'highwgd': 'High WGD'}
results['WGD'] = results['WGD'].apply(lambda x: wgd_rename[x])

#%% Figure 2B
print('Fig 2B')
fig, ax = plt.subplots(figsize=(plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']))

cur_results = results.copy()
cur_results['Distance'] = cur_results['grf']
cur_results = cur_results.loc[cur_results['Rate'].isin([0.01, 0.025, 0.05])]

for i, method in enumerate(['MEDICC2', 'Euclidean Min. Ev.', 'Euclidean NJ', 'Hamming Min. Ev.', 'Hamming NJ', 'MEDALT']):
    cur_method_results = cur_results.loc[cur_results['Method'] == method]
    n_leaves = np.sort(np.unique(cur_method_results['Number of Leaves']))
    nr = cur_method_results.groupby('Number of Leaves').count().iloc[0, 0]

    ax.plot(n_leaves, cur_method_results.groupby('Number of Leaves')['Distance'].mean().values)

    ax.errorbar(x=n_leaves, 
                y=cur_method_results.groupby('Number of Leaves')['Distance'].mean().values.astype(float),
                yerr=cur_method_results.groupby('Number of Leaves')['Distance'].std() / np.sqrt(nr),
                capsize=5, capthick=2, ms=5, marker='o', label=method, color='C' + str(i))
ax.set_xlabel('Number of leaves')
ax.set_ylabel('Generalized RF distance')

plt.legend(facecolor='white', frameon=True, framealpha=0.75)

plt.tight_layout()

fig.savefig('final_figures/Fig_2B.pdf', bbox_inches='tight')
fig.savefig('final_figures/Fig_2B.png', bbox_inches='tight', dpi=300)
plt.close()

#%% Figure Supp 3A Full ranges
print('Supp 3A')
fig, axs = plt.subplots(ncols=3, nrows=2, sharey=True,
                        figsize=(plotting_params['WIDTH_FULL'], plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']))

cur_results = results.copy()
cur_results = cur_results.loc[cur_results['Method'] != 'Sitka']
cur_results['Distance'] = cur_results['grf']

cur_results = cur_results.loc[cur_results['Number of Leaves'].isin([5, 10, 15, 20])]

for wgd, ax in zip(['No WGD', 'Low WGD', 'High WGD'], axs.T):
    cur_wgd_results = cur_results.loc[cur_results['WGD'] == wgd]

    n_leaves = np.sort(np.unique(cur_wgd_results['Number of Leaves']))
    rates = np.sort(np.unique(cur_wgd_results['Rate']))

    for i, method in enumerate(['MEDICC2', 'Euclidean Min. Ev.', 'Euclidean NJ', 'Hamming Min. Ev.', 'Hamming NJ', 'MEDALT']):
        cur_method_results = cur_wgd_results.loc[cur_wgd_results['Method'] == method]
        nr = cur_method_results.groupby('Number of Leaves').count().iloc[0, 0]

        ax[0].errorbar(x=n_leaves, y=cur_method_results.groupby('Number of Leaves')['Distance'].mean(),
                       yerr=cur_method_results.groupby('Number of Leaves')['Distance'].std() / np.sqrt(nr),
                       capsize=5, capthick=2, ms=5, marker='o', label=method, color='C' + str(i))
        ax[1].errorbar(x=rates, y=cur_method_results.groupby('Rate')['Distance'].mean(),
                       yerr=cur_method_results.groupby('Rate')['Distance'].std() / np.sqrt(nr),
                       capsize=5, capthick=2, ms=5, marker='o', label=method, color='C' + str(i))

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
fig.savefig('final_figures/Supp_3A.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_3A.png', bbox_inches='tight', dpi=300)
plt.close()


#%% Figure Supp 3B (Include Sitka)
print('Supp 3B')
fig, ax = plt.subplots(figsize=(0.8*plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']))

cur_results = results.copy()
cur_results['Distance'] = cur_results['grf']
cur_results = cur_results.loc[cur_results['WGD'] == 'Low WGD']
cur_results = cur_results.loc[cur_results['Rate'].isin([0.01, 0.025, 0.05])]

nr = cur_results.groupby('Number of Leaves').count().iloc[0, 0]

# for i, method in enumerate(['Euclidean Min. Ev.', 'Euclidean NJ', 'Manhattan NJ', 'MEDICC2']):
for i, method in enumerate(['MEDICC2', 'Euclidean Min. Ev.', 'Euclidean NJ', 'Hamming Min. Ev.', 'Hamming NJ', 'MEDALT', 'Sitka']):
    cur_method_results = cur_results.loc[cur_results['Method'] == method]
    n_leaves = np.sort(np.unique(cur_method_results['Number of Leaves']))
    nr = cur_method_results.groupby('Number of Leaves').count().iloc[0, 0]
    ax.plot(n_leaves, cur_method_results.groupby('Number of Leaves')['Distance'].mean().values)

    ax.errorbar(x=n_leaves,
                y=cur_method_results.groupby('Number of Leaves')['Distance'].mean().values.astype(float),
                yerr=cur_method_results.groupby('Number of Leaves')['Distance'].std() / np.sqrt(nr),
                capsize=5, capthick=2, ms=5, marker='o', label=method, color='C' + str(i))
ax.set_xlabel('Number of leaves')
ax.set_ylabel('Generalized RF distance')

random = random_results.groupby('leaves').mean()['grf']
ax.plot(random.index, random.values, '--', color='grey', lw=3, label='random tree')

plt.legend(facecolor='white', frameon=True, framealpha=0.75)

plt.tight_layout()
fig.savefig('final_figures/Supp_3B.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_3B.png', bbox_inches='tight', dpi=300)
plt.close()

#%% Figure Supp 3C (Quartet and RF)
print('Supp 3C')
fig, axs = plt.subplots(ncols=2, figsize=(0.8*2*plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']))

for ax, distance, name in zip(axs, ['quartet', 'rf'], ['Quartet', 'Robinson-Foulds']):
    cur_results = results.copy()
    cur_results['Distance'] = cur_results[distance]
    cur_results = cur_results.loc[cur_results['WGD'] == 'Low WGD']
    cur_results = cur_results.loc[cur_results['Rate'].isin([0.01, 0.025, 0.05])]
    
    if distance == 'quartet':
        cur_results = cur_results.loc[cur_results['Number of Leaves'] < 500]

    for i, method in enumerate(['MEDICC2', 'Euclidean Min. Ev.', 'Euclidean NJ', 'Hamming Min. Ev.', 'Hamming NJ', 'MEDALT', 'Sitka']):
        cur_method_results = cur_results.loc[cur_results['Method'] == method]
        nr = cur_method_results.groupby('Number of Leaves').count().iloc[0, 0]
        n_leaves = np.sort(np.unique(cur_method_results['Number of Leaves']))

        ax.errorbar(x=n_leaves, y=cur_method_results.groupby('Number of Leaves')['Distance'].mean(),
                    yerr=cur_method_results.groupby('Number of Leaves')['Distance'].std() / np.sqrt(nr),
                    capsize=5, capthick=2, ms=5, marker='o', label=method, color='C' + str(i))
    ax.set_xlabel('Number of leaves')
    ax.set_ylabel('{} distance'.format(name))

axs[1].legend(loc='upper left')

plt.tight_layout()

fig.savefig('final_figures/Supp_3C.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_3C.png', bbox_inches='tight', dpi=300)
plt.close()
