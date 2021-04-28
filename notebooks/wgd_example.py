#%% import and load data
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import fstlib
# %matplotlib inline
import medicc

fst = fstlib.Fst.read('../objects/wgd_asymm.fst')
#%% Load data and run MEDICC
input_df = medicc.io.read_and_parse_input_data(
    '../examples/WGD_example/example_input.tsv')

sample_labels, pdms, nj_tree, final_tree, output_df = medicc.main(
    input_df,
    fst)

#%% Save files
medicc.io.write_tree_files(final_tree, '../examples/output/WGD_example_final_tree')
output_df.to_csv('../examples/output/WGD_example_output_df.tsv', sep = '\t')

#%% Plot tree and copy number track
fig = medicc.plot.plot_cn_profiles(
	output_df,
	final_tree,
	title="WGD Example",
	normal_name='diploid',
	tree_width_scale=1,
	track_width_scale=0.75,
	height_scale=1,
	hide_normal_chromosomes=False,
	label_func=lambda x: x.replace('_', ' ').replace('taxon', 'sample '))
fig.savefig('../examples/output/WGD Example.pdf', bbox_inches='tight')
