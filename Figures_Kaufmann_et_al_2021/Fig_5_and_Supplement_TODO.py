# %%
import sys

import matplotlib.pyplot as plt

sys.path.append('..')
import medicc
from medicc.plot import plot_cn_heatmap

from plotting_params import plotting_params, set_plotting_params

set_plotting_params()
plt.rc('xtick', labelsize=plotting_params['FONTSIZE_TINY'])

# %% Plot Fig 5
input_df = medicc.io.read_and_parse_input_data('data/Fig_5_df.tsv')
final_tree = medicc.io.import_tree("data/Fig_5_tree.new", 'diploid')

fig = plot_cn_heatmap(input_df, final_tree=final_tree, y_posns=None, cmax=8,
                      alleles=['cn_a', 'cn_b'], tree_width_ratio=1, cbar_width_ratio=0.05,
                      figsize=(plotting_params['WIDTH_FULL'], plotting_params['WIDTH_HALF']), 
                      tree_line_width=0.5, tree_marker_size=0.5,
                      tree_label_colors=lambda x: 'grey', tree_label_func=None)

fig.savefig('final_figures/Fig_5.pdf', pad_inches=0)
fig.savefig('final_figures/Fig_5.png', pad_inches=0, dpi=600)


#%% Plot Supplement TN1 TODO
input_df = medicc.io.read_and_parse_input_data('data/Supp_4_df.tsv')
final_tree = medicc.io.import_tree("data/Supp_4_tree.new", 'diploid')

fig = plot_cn_heatmap(input_df, final_tree=final_tree, y_posns=None, cmax=8,
                      alleles=['cn_a', 'cn_b'], tree_width_ratio=1, cbar_width_ratio=0.05,
                      figsize=(plotting_params['WIDTH_FULL'], plotting_params['WIDTH_HALF']), 
                      tree_line_width=0.5, tree_marker_size=0.5,
                      tree_label_colors=lambda x: 'grey', tree_label_func=None)

fig.savefig('final_figures/Supp_4.pdf', pad_inches=0)
fig.savefig('final_figures/Supp_4.png', pad_inches=0, dpi=600)

#%% Plot Supplement Fastme trees
