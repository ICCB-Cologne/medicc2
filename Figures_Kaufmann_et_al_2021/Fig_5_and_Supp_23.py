# %%
import sys

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

sys.path.append('..')
import medicc

from plotting_params import plotting_params, set_plotting_params

set_plotting_params()
plt.rc('xtick', labelsize=plotting_params['FONTSIZE_TINY'])


# %% Define plotting function here so the plot doesn't change if plot_cn_heatmap is updated in the repo
def plot_cn_heatmap(input_df, final_tree=None, y_posns=None, cmax=8, 
                    alleles='total', tree_width_ratio=1, cbar_width_ratio=0.02, 
                    figsize=(20, 10), tree_line_width=0.5, tree_marker_size=0.5,
                    tree_label_colors=None, tree_label_func=None, cmap='coolwarm'):
    
    cur_sample_labels = np.unique(input_df.index.get_level_values('sample_id'))

    if not isinstance(alleles, list) and not isinstance(alleles, tuple):
        alleles = [alleles]
    nr_alleles = len(alleles)

    if final_tree is None:
        fig, axs = plt.subplots(figsize=figsize, ncols=1+nr_alleles, sharey=False,
                                gridspec_kw={'width_ratios': nr_alleles*[1] + [cbar_width_ratio]})

        if y_posns is None:
            y_posns = {s: i for i, s in enumerate(np.sort(cur_sample_labels))}

        cn_axes = axs[:-1]
    else:
        fig, axs = plt.subplots(figsize=figsize, ncols=2+nr_alleles, sharey=False,
                                gridspec_kw={'width_ratios': [tree_width_ratio] + nr_alleles*[1] + [cbar_width_ratio]})
        tree_ax = axs[0]
        cn_axes = axs[1:-1]

        y_posns = {k.name:v for k, v in medicc.plot._get_y_positions(final_tree, adjust=False).items()}
        
        _ = medicc.plot.plot_tree(final_tree, ax=tree_ax,
                      label_func=tree_label_func if tree_label_func is not None else lambda x: '',
                      hide_internal_nodes=True, show_branch_lengths=False, show_events=False,
                      line_width=tree_line_width, marker_size=tree_marker_size,
                      title='', label_colors=tree_label_colors)
        tree_ax.set_axis_off()
        tree_ax.set_axis_off()
        fig.set_constrained_layout_pads(w_pad=0, h_pad=0, hspace=0.0, wspace=100)

    cax = axs[-1]

    ind = [y_posns.get(x, -1) for x in cur_sample_labels]
    cur_sample_labels = cur_sample_labels[np.argsort(ind)]
    color_norm = mcolors.TwoSlopeNorm(vmin=0, vcenter=1, vmax=min(
        cmax, np.max(input_df.values.astype(int))))

    chr_ends = input_df.loc[cur_sample_labels[0]].copy()
    chr_ends['end_pos'] = np.cumsum([1]*len(chr_ends))
    chr_ends = chr_ends.reset_index().groupby('chrom').max()['end_pos']
    chr_ends.dropna(inplace=True)

    x_pos = np.append([0], np.cumsum(input_df.loc[cur_sample_labels].astype(int).unstack(
        'sample_id').loc[:, (alleles[0])].loc[:, cur_sample_labels].eval('end-start').values))
    y_pos = np.arange(len(cur_sample_labels)+1)+0.5

    for ax, allele in zip(cn_axes, alleles):
        im = ax.pcolormesh(x_pos, y_pos,
                        input_df.loc[cur_sample_labels].astype(int).unstack(
                            'sample_id').loc[:, (allele)].loc[:, cur_sample_labels].values.T,
                        cmap=cmap,
                        norm=color_norm)

        for _, line in chr_ends.iteritems():
            ax.axvline(x_pos[line], color='black', linewidth=0.75)
        
        xtick_pos = np.append([0], x_pos[chr_ends.values][:-1])
        xtick_pos = (xtick_pos + np.roll(xtick_pos, -1))/2
        xtick_pos[-1] += x_pos[-1]/2
        ax.set_xticks(xtick_pos)
        ax.set_xticklabels([x[3:] for x in chr_ends.index], ha='center', rotation=90, va='bottom')
        ax.tick_params(width=0)
        ax.xaxis.set_tick_params(labelbottom=False, labeltop=True, bottom=False)
        ax.set_yticks([])

    cax.pcolormesh([0, 1],
                   np.arange(0, cmax+2),
                   np.arange(0, cmax+1)[:, np.newaxis],
                   cmap=cmap,
                   norm=color_norm)
    cax.set_xticks([])
    cax.set_yticks(np.arange(0, cmax+1)+0.5)
    if np.max(input_df.values.astype(int)) > cmax:
        cax.set_yticklabels([str(x) + '+' if x == cmax else str(x)
                             for x in np.arange(0, cmax+1)], ha='left')
    else:
        cax.set_yticklabels(np.arange(0, cmax+1), ha='left')
    cax.yaxis.set_tick_params(left=False, labelleft=False, labelright=True)

    for ax in axs[:-1]:
        ax.set_ylim(len(cur_sample_labels)+1, 0)

    return fig


# %% Plot Fig 5
input_df = medicc.io.read_and_parse_input_data('data/Fig_5_df.tsv')
final_tree = medicc.io.import_tree("data/Fig_5_tree.new", 'diploid')

fig = plot_cn_heatmap(input_df, final_tree=final_tree, y_posns=None, cmax=6,
                             alleles=['cn_a', 'cn_b'], tree_width_ratio=1, cbar_width_ratio=0.05,
                             figsize=(plotting_params['WIDTH_FULL'], plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']), 
                             tree_line_width=0.5, tree_marker_size=0.5,
                             tree_label_colors=lambda x: 'grey', tree_label_func=None)

fig.savefig('final_figures/Fig_5.pdf', bbox_inches='tight')
fig.savefig('final_figures/Fig_5.png', bbox_inches='tight', dpi=600)


#%% Plot Supplement TN1 (Supp 23A)
input_df = medicc.io.read_and_parse_input_data('data/Supp_23A_df.tsv')
final_tree = medicc.io.import_tree("data/Supp_23A_tree.new", 'diploid')

fig = plot_cn_heatmap(input_df, final_tree=final_tree, y_posns=None, cmax=6,
                             alleles=['cn_a', 'cn_b'], tree_width_ratio=1, cbar_width_ratio=0.05,
                             figsize=(plotting_params['WIDTH_FULL'], plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']), 
                             tree_line_width=0.5, tree_marker_size=0.5,
                             tree_label_colors=lambda x: 'grey', tree_label_func=None)

fig.savefig('final_figures/Supp_23A.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_23A.png', bbox_inches='tight', dpi=600)

#%% Plot Supplement Fastme trees (Supp 23B)
input_df = medicc.io.read_and_parse_input_data('data/Supp_23A_df.tsv')
final_tree = medicc.io.import_tree("data/Supp_23B_tree.new", 'diploid')

fig = plot_cn_heatmap(input_df, final_tree=final_tree, y_posns=None, cmax=6,
                             alleles=['cn_a', 'cn_b'], tree_width_ratio=1, cbar_width_ratio=0.05,
                             figsize=(plotting_params['WIDTH_FULL'], plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']),
                             tree_line_width=0.5, tree_marker_size=0.5,
                             tree_label_colors=lambda x: 'grey', tree_label_func=None)
#fig.get_axes()[0].set_title('Minimum Evolution\nPatient TN1')

fig.savefig('final_figures/Supp_23B.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_23B.png', bbox_inches='tight', dpi=600)

#%% Plot Supplement Fastme trees (Supp 23C)
input_df = medicc.io.read_and_parse_input_data('data/Fig_5_df.tsv')
final_tree = medicc.io.import_tree("data/Supp_23C_tree.new", 'diploid')

fig = plot_cn_heatmap(input_df, final_tree=final_tree, y_posns=None, cmax=6,
                             alleles=['cn_a', 'cn_b'], tree_width_ratio=1, cbar_width_ratio=0.05,
                             figsize=(plotting_params['WIDTH_FULL'], plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']),
                             tree_line_width=0.5, tree_marker_size=0.5,
                             tree_label_colors=lambda x: 'grey', tree_label_func=None)
#fig.get_axes()[0].set_title('Minimum Evolution\nPatient TN2')

fig.savefig('final_figures/Supp_23C.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_23C.png', bbox_inches='tight', dpi=600)
