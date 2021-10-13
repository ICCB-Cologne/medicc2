# imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from plotting_params import plotting_params, set_plotting_params

set_plotting_params()

plotting_data = pd.read_csv('data/Fig_2A.tsv', index_col=0, sep='\t')

# rename MEDICC to MEDICC2 entries
plotting_data['method'] = plotting_data['method'].apply(
    lambda x: {'MEDICC': 'MEDICC2', 'next seg': 'next segment', 'random chrom': 'random segment'}.get(x, x))


# Fig 2A: SV validation filtered at 10^7 and only for duplications and deletions
cur_plotting_data = plotting_data.loc[plotting_data['filter']=='filter_7_dup_del']
max_y = cur_plotting_data['count'].max()

fig, ax = plt.subplots(figsize=(
    plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']))

sns.barplot(data=cur_plotting_data.loc[cur_plotting_data['start_both'] == 'start matches'],
            x='method', y='count', ax=ax, color='C0', label='start matches')
sns.barplot(data=cur_plotting_data.loc[cur_plotting_data['start_both'] == 'start/end matches'],
            x='method', y='count', ax=ax, color='C1', label='start and end matches')

for n, method in enumerate(['MEDICC2', 'next segment', 'random segment']):
    cur_count = cur_plotting_data.loc[(cur_plotting_data['method'] == method) & (
        cur_plotting_data['start_both'] == 'start/end matches'), 'count'].values[0]
    cur_fraction = cur_count / cur_plotting_data.loc[(cur_plotting_data['method'] == method) & (cur_plotting_data['start_both'] == 'start matches'), 'count'].values[0]
    ax.text(n, cur_count-(max_y/50), f"{int(np.round(100*cur_fraction, 0))}%",
            ha='center', va='top', fontsize=plotting_params['FONTSIZE_MEDIUM'])

ax.set_xticklabels(ax.get_xticklabels(), fontsize=plotting_params['FONTSIZE_MEDIUM'])
ax.set_ylabel('count')
ax.set_xlabel('')
ax.legend(loc='best')
fig.suptitle('Overlap of MEDICC2 events with structural variants\nAggregated over all PCAWG samples\nFiltered at 10Mbp, only dup/del', fontsize=plotting_params['FONTSIZE_LARGE'])
plt.tight_layout()

fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Fig_2A.pdf', bbox_inches='tight')
fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Fig_2A.png',
            bbox_inches='tight', dpi=600)
plt.close()

# Supp TODO: SV validation for all filters
fig, axs = plt.subplots(ncols=6, figsize=(1.5*plotting_params['WIDTH_FULL'], 
                                          plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']))
for ax, filter in zip(axs,
                      np.sort(np.unique(plotting_data['filter']))):
    cur_plotting_data = plotting_data.loc[plotting_data['filter'] == filter]
    max_y = cur_plotting_data['count'].max()

    sns.barplot(data=cur_plotting_data.loc[cur_plotting_data['start_both'] == 'start matches'],
                x='method', y='count', ax=ax, color='C0', label='start matches')
    sns.barplot(data=cur_plotting_data.loc[cur_plotting_data['start_both'] == 'start/end matches'],
                x='method', y='count', ax=ax, color='C1', label='start and end matches')
    ax.set_title(f'Filtered at $10^{filter.split("_")[1]}$bp{", only dup/del" if filter.split("_")[-1]=="del" else ""}')
    for n, method in enumerate(['MEDICC2', 'next segment', 'random segment']):
        cur_count = cur_plotting_data.loc[(cur_plotting_data['method'] == method) & (
            cur_plotting_data['start_both'] == 'start/end matches'), 'count'].values[0]
        cur_fraction = cur_count / cur_plotting_data.loc[(cur_plotting_data['method'] == method) & (
            cur_plotting_data['start_both'] == 'start matches'), 'count'].values[0]
        ax.text(n, cur_count-(max_y/50), f"{int(np.round(100*cur_fraction, 0))}%",
                ha='center', va='top', fontsize=plotting_params['FONTSIZE_MEDIUM'])

    ax.set_xticklabels(['MEDICC2', 'next\nsegment', 'random\nsegment'], fontsize=plotting_params['FONTSIZE_MEDIUM'])
    ax.set_ylabel('')
    ax.set_xlabel('')

axs[-1].set_ylabel('count')
axs[-1].legend(bbox_to_anchor=(1, 1))
fig.suptitle('Overlap of MEDICC2 events with structural variants\nAggregated over all PCAWG samples',
             fontsize=plotting_params['FONTSIZE_LARGE'])
plt.tight_layout()

fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Supp_TODO_SV_validation.pdf', bbox_inches='tight')
fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Supp_TODO_SV_validation.png',
            bbox_inches='tight', dpi=600)
plt.close()
