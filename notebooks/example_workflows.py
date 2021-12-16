#%% import and load data
import os
import sys

import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import fstlib
import medicc

fst = fstlib.Fst.read('../objects/wgd_asymm.fst')

# 1) Simple example
#%% Load data from tsv table and run MEDICC
input_df = medicc.io.read_and_parse_input_data('../examples/simple_example/simple_example.tsv')

sample_labels, pdms, nj_tree, final_tree, output_df, events_df = medicc.main(
    input_df,
    fst)

#%% Save files and plot tree and copy number track
medicc.io.write_tree_files(final_tree, '../examples/output/simple_example_final_tree')
output_df.to_csv('../examples/output/simple_example_output_df.tsv', sep='\t')

fig = medicc.plot.plot_cn_profiles(
	output_df,
	final_tree,
	title="Simple Example",
	normal_name='diploid',
	tree_width_scale=1,
	track_width_scale=0.75,
	height_scale=1,
	hide_normal_chromosomes=False,
	label_func=lambda x: x.replace('_', ' ').replace('taxon', 'sample '))
fig.savefig('../examples/output/simple_example_cn_track.pdf', bbox_inches='tight')


# 2) OV03-04
#%% Alternatively one can input fasta files containing the copy numbers (old MEDICC format)
input_df = medicc.io.read_and_parse_input_data('../examples/OV03-04/OV03-04_descr.txt',
                                               input_type='fasta')


fasta_desc_file = pd.read_csv('../examples/OV03-04/OV03-04_descr.txt', sep=' ', header=None)
fasta_desc_file.head()

#%% Run MEDICC
sample_labels, pdms, nj_tree, final_tree, output_df, events_df = medicc.main(input_df, fst)

#%% Save files and plot tree and copy number track
medicc.io.write_tree_files(final_tree, '../examples/output/OV03-04_final_tree')
output_df.to_csv('../examples/output/OV03-04_output_df.tsv', sep='\t')

fig = medicc.plot.plot_cn_profiles(
	output_df,
	final_tree,
	title="OV03-04",
	normal_name='diploid',
	tree_width_scale=1,
	track_width_scale=0.75,
	height_scale=1,
	hide_normal_chromosomes=False,
	label_func=lambda x: x.replace('_', ' ').replace('taxon', 'sample '))
fig.savefig('../examples/output/OV03-04_cn_track.pdf', bbox_inches='tight')


#%% test star topology
p_star = medicc.stats.star_topology_test(pdms['total'].values)
p_clock = medicc.stats.molecular_clock_test(pdms['total'].values)
print("P star: %.4f, p clock: %.4f" % (p_star, p_clock))

# 3) WGD example
#%% Load data and run MEDICC
input_df = medicc.io.read_and_parse_input_data('../examples/WGD_example/WGD_example.tsv')
sample_labels, pdms, nj_tree, final_tree, output_df, events_df = medicc.main(input_df, fst)

#%% Save files and plot tree and copy number track
medicc.io.write_tree_files(final_tree, '../examples/output/WGD_example_final_tree')
output_df.to_csv('../examples/output/WGD_example_output_df.tsv', sep='\t')

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
fig.savefig('../examples/output/WGD_Example.pdf', bbox_inches='tight')

# %%
