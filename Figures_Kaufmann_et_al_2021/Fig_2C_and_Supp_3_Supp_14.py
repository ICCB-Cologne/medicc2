# imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import norm

from plotting_params import plotting_params, set_plotting_params

set_plotting_params()


def z_score(cur_plotting_data, method1, method2):

    x1 = cur_plotting_data.loc[(cur_plotting_data['method'] == method1) & (
        cur_plotting_data['start_both'] == 'start/end matches'), 'count'].values[0]
    x2 = cur_plotting_data.loc[(cur_plotting_data['method'] == method2) & (
        cur_plotting_data['start_both'] == 'start/end matches'), 'count'].values[0]
    n1 = cur_plotting_data.loc[(cur_plotting_data['method'] == method1) & (
        cur_plotting_data['start_both'] == 'one matches'), 'count'].values[0]
    n2 = cur_plotting_data.loc[(cur_plotting_data['method'] == method2) & (
        cur_plotting_data['start_both'] == 'one matches'), 'count'].values[0]

    p1 = x1/n1
    p2 = x2/n2

    p_dash = (n1*p1 + n2*p2) / (n1+n2)
    z = (p1 - p2) / np.sqrt(p_dash * (1 - p_dash) * (1/n1 + 1/n2))
    p = norm.cdf(-z)

    return p


plotting_data = pd.read_csv('data/Fig_2C.tsv', index_col=0, sep='\t')

# rename MEDICC to MEDICC2 entries
plotting_data['method'] = plotting_data['method'].apply(
    lambda x: {'MEDICC': 'MEDICC2', 'next seg': 'next segment', 'random chrom': 'random segment'}.get(x, x))


# Fig 2C: SV validation filtered at 10^7 and only for duplications and deletions
cur_plotting_data = plotting_data.loc[plotting_data['filter']=='filter_7_dup_del']
cur_plotting_data.loc[cur_plotting_data['start_both'] == 'one matches', 'count'] = cur_plotting_data.loc[cur_plotting_data['start_both'] ==
                                                                                                         'one matches', 'count'].values + cur_plotting_data.loc[cur_plotting_data['start_both'] == 'start/end matches', 'count'].values

max_y = cur_plotting_data['count'].max()

fig, ax = plt.subplots(figsize=(
    plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']))

sns.barplot(data=cur_plotting_data.loc[cur_plotting_data['start_both'] == 'one matches'],
            x='method', y='count', ax=ax, color='C0', label='_nolegend_') # label was: "One segment boundary matches"
sns.barplot(data=cur_plotting_data.loc[cur_plotting_data['start_both'] == 'start/end matches'],
            x='method', y='count', ax=ax, color='C1', label='_nolegend_') # label was: "Both segment boundaries match"

for n, method in enumerate(['MEDICC2', 'neighbor', 'random']):
    cur_count = cur_plotting_data.loc[(cur_plotting_data['method'] == method) & (
        cur_plotting_data['start_both'] == 'start/end matches'), 'count'].values[0]
    cur_fraction = cur_count / cur_plotting_data.loc[(cur_plotting_data['method'] == method) & (cur_plotting_data['start_both'] == 'one matches'), 'count'].values[0]
    ax.text(n, cur_count+(max_y/50), f"{int(np.round(100*cur_fraction, 0))}%",
            ha='center', va='bottom', color='white', fontsize=plotting_params['FONTSIZE_MEDIUM'])


# statistical annotation, note that actual calculation is done in external notebook using a z-test
yticks = ax.get_yticks()
x1, x2 = 0, 1
h = (max_y/50)
hl = (max_y/25)
y = cur_plotting_data['count'].max() + h
p_value = z_score(cur_plotting_data, 'MEDICC2', 'neighbor')
cur_text = "$P$<.0001" if p_value < 0.0001 else f"P={p_value:.2g}"
ax.plot([x1, x1, x2, x2], [y, y+hl, y+hl, y], lw=1.5, c='black')
ax.text((x1+x2)*.5, y+hl*1.15, cur_text, ha='center', va='bottom', color='black', fontsize=plotting_params['FONTSIZE_MEDIUM'])

x1, x2 = 0, 2
h = 4*(max_y/50)
y = cur_plotting_data['count'].max() + 1.5*h
p_value = z_score(cur_plotting_data, 'MEDICC2', 'random')
cur_text = "$P$<.0001" if p_value < 0.0001 else f"P={p_value:.2g}"
ax.plot([x1, x1, x2, x2], [y, y+hl, y+hl, y], lw=1.5, c='black')
ax.text((x1+x2)*.5, y+hl*1.15, cur_text, ha='center', va='bottom', color='black', fontsize=plotting_params['FONTSIZE_MEDIUM'])

ax.set_ylim(0, ax.get_ylim()[1] + hl)
ax.set_yticks(yticks[:-1])
ax.set_xticklabels(['MEDICC2 event\nboundary', 'Neighboring segment\nboundary', 'Random segment\nboundary'],
                   fontsize=plotting_params['FONTSIZE_MEDIUM'])
ax.set_ylabel('Frequency')
ax.set_xlabel('')
ax.legend(loc='center right')
#fig.suptitle('Overlap of MEDICC2 events with\nstructural variants in all PCAWG samples', fontsize=plotting_params['FONTSIZE_LARGE'])
plt.tight_layout()

fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Fig_2C.pdf', bbox_inches='tight')
fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Fig_2C.png',
            bbox_inches='tight', dpi=600)
plt.close()

# Supp 3: SV validation for all filters for all PCAWG data
fig, axs = plt.subplots(ncols=3, nrows=2, sharex=True, figsize=(plotting_params['WIDTH_FULL'], 
                                          plotting_params['WIDTH_FULL']/plotting_params['ASPECT_RATIO']))
axs = axs.ravel()
for ax, filter in zip(axs, np.sort(np.unique(plotting_data['filter']))):
    cur_plotting_data = plotting_data.loc[plotting_data['filter'] == filter]
    cur_plotting_data.loc[cur_plotting_data['start_both'] == 'one matches', 'count'] = cur_plotting_data.loc[cur_plotting_data['start_both'] ==
                                                                                                            'one matches', 'count'].values + cur_plotting_data.loc[cur_plotting_data['start_both'] == 'start/end matches', 'count'].values
    max_y = cur_plotting_data['count'].max()

    sns.barplot(data=cur_plotting_data.loc[cur_plotting_data['start_both'] == 'one matches'],
                x='method', y='count', ax=ax, color='C0', label='One segment boundary matches')
    sns.barplot(data=cur_plotting_data.loc[cur_plotting_data['start_both'] == 'start/end matches'],
                x='method', y='count', ax=ax, color='C1', label='Both segment boundaries match')
    ax.set_title(f'Filtered at $10^{filter.split("_")[1]}$bp\n{"only duplications and deletions" if filter.split("_")[-1]=="del" else ""}')
    for n, method in enumerate(['MEDICC2', 'neighbor', 'random']):
        cur_count = cur_plotting_data.loc[(cur_plotting_data['method'] == method) & (
            cur_plotting_data['start_both'] == 'start/end matches'), 'count'].values[0]
        cur_fraction = cur_count / cur_plotting_data.loc[(cur_plotting_data['method'] == method) & (
            cur_plotting_data['start_both'] == 'one matches'), 'count'].values[0]
        ax.text(n, cur_count+(max_y/50), f"{int(np.round(100*cur_fraction, 0))}%",
                ha='center', va='bottom', fontsize=plotting_params['FONTSIZE_MEDIUM'])

    ax.set_xticklabels(['MEDICC2', 'Neighboring\nsegment', 'Random\nsegment'],
                       fontsize=plotting_params['FONTSIZE_MEDIUM'], rotation=45)
    ax.set_ylabel('')
    ax.set_xlabel('')

    # statistical annotation
    yticks = ax.get_yticks()
    x1, x2 = 0, 1
    h = (max_y/50)
    hl = (max_y/50)
    y = cur_plotting_data['count'].max() + h
    p_value = z_score(cur_plotting_data, 'MEDICC2', 'neighbor')
    cur_text = "$P$<.0001" if p_value < 0.0001 else f"P={p_value:.2g}"
    ax.plot([x1, x1, x2, x2], [y, y+hl, y+hl, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+hl*1.15, cur_text, ha='center', va='bottom',
            color='black', fontsize=plotting_params['FONTSIZE_MEDIUM'])

    x1, x2 = 0, 2
    h = 4*(max_y/50)
    y = cur_plotting_data['count'].max() + 1.5*h
    p_value = z_score(cur_plotting_data, 'MEDICC2', 'random')
    cur_text = "$P$<.0001" if p_value < 0.0001 else f"P={p_value:.2g}"
    ax.plot([x1, x1, x2, x2], [y, y+hl, y+hl, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+hl*1.15, cur_text, ha='center', va='bottom',
            color='black', fontsize=plotting_params['FONTSIZE_MEDIUM'])

    ax.set_ylim(0, ax.get_ylim()[1] + hl)
    ax.set_yticks(yticks[:-1])

axs[0].set_ylabel('Frequency')
axs[3].set_ylabel('Frequency')
#axs[-1].legend(bbox_to_anchor=(1, 1))
#fig.suptitle('Overlap of MEDICC2 events with structural variants in all PCAWG samples',
#             fontsize=plotting_params['FONTSIZE_LARGE'])
plt.tight_layout()

fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Supp_3.pdf', bbox_inches='tight')
fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Supp_3.png',
            bbox_inches='tight', dpi=600)
plt.close()


# Supp 14: SV validation for all filters for all Gundem data
plotting_data = pd.read_csv('data/Fig_2C_Gundem.tsv', index_col=0, sep='\t')
plotting_data['method'] = plotting_data['method'].apply(
    lambda x: {'MEDICC': 'MEDICC2', 'next seg': 'next segment', 'random chrom': 'random segment'}.get(x, x))

fig, axs = plt.subplots(ncols=3, nrows=2, sharex=True, figsize=(plotting_params['WIDTH_FULL'], 
                                          plotting_params['WIDTH_FULL']/plotting_params['ASPECT_RATIO']))
axs = axs.ravel()

for ax, filter in zip(axs,
                      np.sort(np.unique(plotting_data['filter']))):
    cur_plotting_data = plotting_data.loc[plotting_data['filter'] == filter]
    cur_plotting_data.loc[cur_plotting_data['start_both'] == 'one matches', 'count'] = cur_plotting_data.loc[cur_plotting_data['start_both'] ==
                                                                                                             'one matches', 'count'].values + cur_plotting_data.loc[cur_plotting_data['start_both'] == 'start/end matches', 'count'].values
    max_y = cur_plotting_data['count'].max()

    sns.barplot(data=cur_plotting_data.loc[cur_plotting_data['start_both'] == 'one matches'],
                x='method', y='count', ax=ax, color='C0', label='One segment boundary matches')
    sns.barplot(data=cur_plotting_data.loc[cur_plotting_data['start_both'] == 'start/end matches'],
                x='method', y='count', ax=ax, color='C1', label='Both segment boundaries match')
    ax.set_title(
        f'Filtered at $10^{filter.split("_")[1]}$bp\n{"only duplications and deletions" if filter.split("_")[-1]=="del" else ""}')
    for n, method in enumerate(['MEDICC2', 'neighbor', 'random']):
        cur_count = cur_plotting_data.loc[(cur_plotting_data['method'] == method) & (
            cur_plotting_data['start_both'] == 'start/end matches'), 'count'].values[0]
        cur_fraction = cur_count / cur_plotting_data.loc[(cur_plotting_data['method'] == method) & (
            cur_plotting_data['start_both'] == 'one matches'), 'count'].values[0]
        ax.text(n, cur_count+(max_y/50), f"{int(np.round(100*cur_fraction, 0))}%",
                ha='center', va='bottom', fontsize=plotting_params['FONTSIZE_MEDIUM'])

    ax.set_xticklabels(['MEDICC2', 'Neighboring\nsegment', 'Random\nsegment'],
                       fontsize=plotting_params['FONTSIZE_MEDIUM'], rotation=45)
    ax.set_ylabel('')
    ax.set_xlabel('')

    # statistical annotation
    yticks = ax.get_yticks()
    x1, x2 = 0, 1
    h = (max_y/50)
    hl = (max_y/50)
    y = cur_plotting_data['count'].max() + h
    p_value = z_score(cur_plotting_data, 'MEDICC2', 'neighbor')
    cur_text = "$P$<.0001" if p_value < 0.0001 else f"P={p_value:.2g}"
    ax.plot([x1, x1, x2, x2], [y, y+hl, y+hl, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+hl*1.15, cur_text, ha='center', va='bottom',
            color='black', fontsize=plotting_params['FONTSIZE_MEDIUM'])

    x1, x2 = 0, 2
    h = 4*(max_y/50)
    y = cur_plotting_data['count'].max() + 1.5*h
    p_value = z_score(cur_plotting_data, 'MEDICC2', 'random')
    cur_text = "$P$<.0001" if p_value < 0.0001 else f"P={p_value:.2g}"
    ax.plot([x1, x1, x2, x2], [y, y+hl, y+hl, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+hl*1.15, cur_text, ha='center', va='bottom',
            color='black', fontsize=plotting_params['FONTSIZE_MEDIUM'])

    ax.set_ylim(0, ax.get_ylim()[1] + hl)
    ax.set_yticks(yticks[:-1])

axs[0].set_ylabel('Frequency')
axs[3].set_ylabel('Frequency')

#axs[2].legend(bbox_to_anchor=(1, 1))
#fig.suptitle('Overlap of MEDICC2 events with structural variants in all Gundem et al. samples',
#             fontsize=plotting_params['FONTSIZE_LARGE'])
plt.tight_layout()

fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Supp_14.pdf',
            bbox_inches='tight')
fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Supp_14.png',
            bbox_inches='tight', dpi=600)
plt.close()
