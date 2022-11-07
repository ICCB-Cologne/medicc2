# imports
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pickle
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

#%% Plotting data
with open('data/Fig_2C_and_Supp_7_Supp_8_Supp_9_Supp_21.pickle', 'rb') as f:
    plot_data = pickle.load(f)

plotting_data_pcawg = plot_data['plotting_data_pcawg']
plotting_data_gundem = plot_data['plotting_data_gundem']
sv_accuracy = plot_data['sv_accuracy']
all_distances = plot_data['all_distances']


#%% Figure 2C
print('Figure 2C')
plotting_data = plotting_data_pcawg

fig, ax = plt.subplots(figsize=(
    plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']))

cur_filter = 'filter_7_dup_del'

cur_plotting_data = plotting_data.loc[plotting_data['filter'] == cur_filter]
cur_plotting_data.loc[cur_plotting_data['start_both'] == 'one matches', 'count'] = cur_plotting_data.loc[cur_plotting_data['start_both'] ==
                                                                                                        'one matches', 'count'].values + cur_plotting_data.loc[cur_plotting_data['start_both'] == 'start/end matches', 'count'].values
max_y = cur_plotting_data['count'].max()

sns.barplot(data=cur_plotting_data.loc[cur_plotting_data['start_both'] == 'one matches'],
            x='method', y='count', ax=ax, color='C0', label='One segment\nboundary matches')
sns.barplot(data=cur_plotting_data.loc[cur_plotting_data['start_both'] == 'start/end matches'],
            x='method', y='count', ax=ax, color='C1', label='Both segment\nboundaries match')
# ax.set_title(f'Filtered at $10^{cur_filter.split("_")[1]}$ bp\n{"only duplications and deletions" if cur_filter.split("_")[-1]=="del" else ""}')
for n, method in enumerate(['MEDICC', 'neighbor', 'random']):
    cur_count = cur_plotting_data.loc[(cur_plotting_data['method'] == method) & (
        cur_plotting_data['start_both'] == 'start/end matches'), 'count'].values[0]
    cur_fraction = cur_count / cur_plotting_data.loc[(cur_plotting_data['method'] == method) & (
        cur_plotting_data['start_both'] == 'one matches'), 'count'].values[0]
    ax.text(n, cur_count+(max_y/50), f"{int(np.round(100*cur_fraction, 0))}%",
            ha='center', va='bottom', fontsize=10, color='white')

ax.set_xticklabels(['MEDICC2', 'Neighboring\nsegment', 'Random\nsegment'],
                   fontsize=10, rotation=45)
ax.set_ylabel('')
ax.set_xlabel('')

# statistical annotation
yticks = ax.get_yticks()
x1, x2 = 0, 1
h = (max_y/50)
hl = (max_y/50)
y = cur_plotting_data['count'].max() + h
p_value = z_score(cur_plotting_data, 'MEDICC', 'neighbor')
cur_text = "$P$<.0001" if p_value < 0.0001 else f"P={p_value:.2g}"
ax.plot([x1, x1, x2, x2], [y, y+hl, y+hl, y], lw=1.5, c='black')
ax.text((x1+x2)*.5, y+hl*1.15, cur_text, ha='center', va='bottom',
        color='black', fontsize=10)

x1, x2 = 0, 2
h = 4*(max_y/50)
y = cur_plotting_data['count'].max() + 1.5*h
p_value = z_score(cur_plotting_data, 'MEDICC', 'random')
cur_text = "$P$<.0001" if p_value < 0.0001 else f"P={p_value:.2g}"
ax.plot([x1, x1, x2, x2], [y, y+hl, y+hl, y], lw=1.5, c='black')
ax.text((x1+x2)*.5, y+hl*1.15, cur_text, ha='center', va='bottom',
        color='black', fontsize=10)

ax.set_ylim(0, ax.get_ylim()[1] + hl)
ax.set_yticks(yticks[:-1])

ax.spines.left.set_bounds((0, ax.get_yticks().max()))

ax.set_ylim(0, ax.get_ylim()[1] + hl)
ax.set_yticks(yticks[:-1])
ax.set_xticklabels(['MEDICC2 event\nboundary', 'Next segment\nboundary', 'Random segment\nboundary'],
                   fontsize=plotting_params['FONTSIZE_MEDIUM'], rotation=0)
ax.set_ylabel('Frequency')
ax.set_xlabel('')
plt.legend(facecolor='white', frameon=True, framealpha=0.75)

plt.tight_layout()
fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Fig_2C.pdf', bbox_inches='tight')
fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Fig_2C.png',
            bbox_inches='tight', dpi=600)
plt.close()


#%% Supp 7
print('Supp 7')

plotting_data = plotting_data_pcawg

fig, axs = plt.subplots(ncols=3, nrows=2, sharex=True, figsize=(plotting_params['WIDTH_FULL'], 
                                          plotting_params['WIDTH_FULL']/plotting_params['ASPECT_RATIO']))
axs = axs.T.ravel()
for ax, cur_filter in zip(axs, np.sort(np.unique(plotting_data['filter']))):
    cur_plotting_data = plotting_data.loc[plotting_data['filter'] == cur_filter]
    cur_plotting_data.loc[cur_plotting_data['start_both'] == 'one matches', 'count'] = cur_plotting_data.loc[cur_plotting_data['start_both'] ==
                                                                                                            'one matches', 'count'].values + cur_plotting_data.loc[cur_plotting_data['start_both'] == 'start/end matches', 'count'].values
    max_y = cur_plotting_data['count'].max()

    sns.barplot(data=cur_plotting_data.loc[cur_plotting_data['start_both'] == 'one matches'],
                x='method', y='count', ax=ax, color='C0', label='One segment boundary matches')
    sns.barplot(data=cur_plotting_data.loc[cur_plotting_data['start_both'] == 'start/end matches'],
                x='method', y='count', ax=ax, color='C1', label='Both segment boundaries match')
    ax.set_title(f'Filtered at $10^{cur_filter.split("_")[1]}$ bp\n{"only duplications and deletions" if cur_filter.split("_")[-1]=="del" else ""}')
    for n, method in enumerate(['MEDICC', 'neighbor', 'random']):
        cur_count = cur_plotting_data.loc[(cur_plotting_data['method'] == method) & (
            cur_plotting_data['start_both'] == 'start/end matches'), 'count'].values[0]
        cur_fraction = cur_count / cur_plotting_data.loc[(cur_plotting_data['method'] == method) & (
            cur_plotting_data['start_both'] == 'one matches'), 'count'].values[0]
        ax.text(n, cur_count+(max_y/50), f"{int(np.round(100*cur_fraction, 0))}%",
                ha='center', va='bottom', fontsize=10)

    ax.set_xticklabels(['MEDICC2', 'Neighboring\nsegment', 'Random\nsegment'],
                       fontsize=10, rotation=45)
    ax.set_ylabel('')
    ax.set_xlabel('')

    # statistical annotation
    yticks = ax.get_yticks()
    x1, x2 = 0, 1
    h = (max_y/50)
    hl = (max_y/50)
    y = cur_plotting_data['count'].max() + h
    p_value = z_score(cur_plotting_data, 'MEDICC', 'neighbor')
    cur_text = "$P$<.0001" if p_value < 0.0001 else f"P={p_value:.2g}"
    ax.plot([x1, x1, x2, x2], [y, y+hl, y+hl, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+hl*1.15, cur_text, ha='center', va='bottom',
            color='black', fontsize=10)

    x1, x2 = 0, 2
    h = 4*(max_y/50)
    y = cur_plotting_data['count'].max() + 1.5*h
    p_value = z_score(cur_plotting_data, 'MEDICC', 'random')
    cur_text = "$P$<.0001" if p_value < 0.0001 else f"P={p_value:.2g}"
    ax.plot([x1, x1, x2, x2], [y, y+hl, y+hl, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+hl*1.15, cur_text, ha='center', va='bottom',
            color='black', fontsize=10)

    ax.set_ylim(0, ax.get_ylim()[1] + hl)
    ax.set_yticks(yticks[:-1])

axs[0].set_ylabel('Frequency')
axs[1].set_ylabel('Frequency')
#axs[-1].legend(bbox_to_anchor=(1, 1))
#fig.suptitle('Overlap of MEDICC2 events with structural variants in all PCAWG samples',
#             fontsize=plotting_params['FONTSIZE_LARGE'])
for ax in axs:
    ax.spines.left.set_bounds((0, ax.get_yticks().max()))
    ax.set_xticklabels(['MEDICC2 event\nboundary', 'Next segment\nboundary', 'Random segment\nboundary'],
                   fontsize=plotting_params['FONTSIZE_MEDIUM'], rotation=45)
    ax.set_xlabel('')
    
axs[4].legend(bbox_to_anchor=(1, 1))
# ax[4].legend(facecolor='white', frameon=True, framealpha=0.75)


plt.tight_layout()
fig.savefig('final_figures/Supp_7.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_7.png', bbox_inches='tight', dpi=300)
plt.close()


#%% Supp 8
print('Supp 8')

fig, axs = plt.subplots(ncols=3, nrows=2, sharex=True, figsize=(plotting_params['WIDTH_FULL'], 
                                          plotting_params['WIDTH_FULL']/plotting_params['ASPECT_RATIO']))

for i, size_filter in enumerate([1e5, 1e6, 1e7]):
    for j, filter_dup_del in enumerate([False, True]):

        ax = axs[j, i]
        distance_to_next_sv = all_distances[(size_filter, filter_dup_del)]
        
        bins = np.arange(0, np.max(distance_to_next_sv), 1e7)
        count, bins = np.histogram(distance_to_next_sv, bins=bins)
        sns.histplot(distance_to_next_sv, bins=bins, kde_kws={'clip': (0.0, bins[-1])}, ax=ax, kde=True, color='C0')
        ax.set_title(
            f'Filtered at $10^{np.round(np.log10(size_filter), 0).astype(int)}$ bp\n{"only duplications and deletions" if filter_dup_del else ""}')
        # ax.set_title(f'filter {int(np.log10(size_filter))}{" (dup/del)" if filter_dup_del else ""}\n'
        #              f'{count[0]} ({count[0]/np.sum(count)*100:.1f}%) are within 10 Mbp');

for ax in axs[1, :]:
    ax.set_xlabel('Distance')

plt.suptitle('Distance between mismatched MEDICC2 and SV breakpoints')

plt.tight_layout()
fig.savefig('final_figures/Supp_8.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_8.png', bbox_inches='tight', dpi=300)
plt.close()


#%% Supp 9
print('Supp 9')
fig, ax = plt.subplots(figsize=(
    plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']))

cur_plotting_data = sv_accuracy.copy()
# cur_plotting_data.loc[cur_plotting_data['start_both'] == 'one matches', 'count'] = cur_plotting_data.loc[cur_plotting_data['start_both'] ==
#                                                                                                          'one matches', 'count'].values + cur_plotting_data.loc[cur_plotting_data['start_both'] == 'start/end matches', 'count'].values
max_y = cur_plotting_data['count'].max()

sns.barplot(data=cur_plotting_data.loc[cur_plotting_data['start_both'] == 'one matches'],
            x='filter', y='count', ax=ax, color='C0', label='One segment boundary matches')
sns.barplot(data=cur_plotting_data.loc[cur_plotting_data['start_both'] == 'start/end matches'],
            x='filter', y='count', ax=ax, color='C1', label='Both segment boundaries match')
ax.set_title('Overlap of structural variants with MEDICC2 events')
for n, filter in enumerate(cur_plotting_data['filter'].unique()):
    cur_count = cur_plotting_data.loc[(cur_plotting_data['filter'] == filter) & (
        cur_plotting_data['start_both'] == 'start/end matches'), 'count'].values[0]
    cur_fraction = cur_count / cur_plotting_data.loc[(cur_plotting_data['filter'] == filter) & (
        cur_plotting_data['start_both'] == 'one matches'), 'count'].values[0]
    ax.text(n, cur_count+(max_y/50), f"{int(np.round(100*cur_fraction, 0))}%",
            ha='center', va='bottom', fontsize=plotting_params['FONTSIZE_MEDIUM'])

# ax.set_xticklabels(['MEDICC2', 'Neighboring\nsegment', 'Random\nsegment'],
#                    fontsize=plotting_params['FONTSIZE_MEDIUM'], rotation=45)

ax.set_xticklabels([f'$10^{cur_filter.get_text().split("_")[1]}$ bp\n{"(dup/del)" if cur_filter.get_text().split("_")[-1]=="del" else ""}' for cur_filter in ax.get_xticklabels()])
                    
ax.set_ylim(0, ax.get_ylim()[1] + hl)
# ax.set_yticks(yticks[:-1])

ax.set_ylabel('Frequency')
ax.set_xlabel('Filter')

ax.legend(bbox_to_anchor=(1, 1))
#fig.suptitle('Overlap of MEDICC2 events with structural variants in all Gundem et al. samples',
#             fontsize=plotting_params['FONTSIZE_LARGE'])
plt.tight_layout()
fig.savefig('final_figures/Supp_9.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_9.png', bbox_inches='tight', dpi=300)
plt.close()


#%% Supp 21 (Gundem SV)
print('Supp 21')
plotting_data = plotting_data_gundem

plotting_data['method'] = plotting_data['method'].apply(
    lambda x: {'MEDICC': 'MEDICC2', 'next seg': 'next segment', 'random chrom': 'random segment'}.get(x, x))

fig, axs = plt.subplots(ncols=3, nrows=2, sharex=True, figsize=(plotting_params['WIDTH_FULL'], 
                                          plotting_params['WIDTH_FULL']/plotting_params['ASPECT_RATIO']))
axs = axs.T.ravel()

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
        f'Filtered at $10^{filter.split("_")[1]}$ bp\n{"only duplications and deletions" if filter.split("_")[-1]=="del" else ""}')
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
axs[1].set_ylabel('Frequency')
#axs[-1].legend(bbox_to_anchor=(1, 1))
#fig.suptitle('Overlap of MEDICC2 events with structural variants in all PCAWG samples',
#             fontsize=plotting_params['FONTSIZE_LARGE'])
for ax in axs:
    ax.spines.left.set_bounds((0, ax.get_yticks().max()))
    ax.set_xticklabels(['MEDICC2 event\nboundary', 'Next segment\nboundary', 'Random segment\nboundary'],
                   fontsize=plotting_params['FONTSIZE_MEDIUM'], rotation=45)
    ax.set_xlabel('')
    
axs[4].legend(bbox_to_anchor=(1, 1))
# ax[4].legend(facecolor='white', frameon=True, framealpha=0.75)

plt.tight_layout()
fig.savefig('final_figures/Supp_21.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_21.png', bbox_inches='tight', dpi=300)
plt.close()