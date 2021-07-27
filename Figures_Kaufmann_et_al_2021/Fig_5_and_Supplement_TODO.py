# %%
import os
import sys

import matplotlib.pyplot as plt

sys.path.append('..')
import medicc
from medicc.plot import plot_cn_heatmap

from plotting_params import plotting_params, set_plotting_params

set_plotting_params()
plt.rc('xtick', labelsize=plotting_params['FONTSIZE_TINY'])

# %% Plot Fig 5
patient = 'TN20'

input_df = medicc.io.read_and_parse_input_data('data/Fig_5_df.tsv')
final_tree = medicc.io.import_tree("data/Fig_5_tree.new", 'diploid')

fig = plot_cn_heatmap(input_df, final_tree=final_tree, y_posns=None, cmax=8,
                      alleles=['cn_a', 'cn_b'], tree_width_ratio=1, cbar_width_ratio=0.05,
                      figsize=(plotting_params['WIDTH_FULL'], plotting_params['WIDTH_HALF']), 
                      title=None, tree_line_width=0.5, tree_marker_size=0.75)

fig.savefig('final_figures/Fig_5.pdf', pad_inches=0)
fig.savefig('final_figures/Fig_5.png', pad_inches=0, dpi=600)


#%% Plot Supplement TODO
