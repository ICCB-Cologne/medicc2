#%% import and load data
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# %matplotlib inline
import medicc

#%%
evo001_df = medicc.io.read_and_parse_input_data('../examples/output/EVO001_EP_final_cn_profiles.tsv')
evo001_tree = medicc.io.import_tree('../examples/output/EVO001_EP_final_tree.new', 'diploid')
ptx011_df = medicc.io.read_and_parse_input_data('../examples/output/PTX011_final_cn_profiles.tsv')
ptx011_tree = medicc.io.import_tree('../examples/output/PTX011_final_tree.new', 'diploid')
ptx011_nowgd_df = medicc.io.read_and_parse_input_data('../examples/output/PTX011-noWGD_final_cn_profiles.tsv')
ptx011_nowgd_tree = medicc.io.import_tree('../examples/output/PTX011-noWGD_final_tree.new', 'diploid')

# %%
def plot_ptx011(df, tree, title='PTX011'):
	plotdat = df.reset_index()
	plotdat['chrom'] = plotdat.chrom.str.replace('chr','')
	plotdat.set_index(['sample_id', 'chrom', 'start', 'end'], inplace=True)
	labels = {'diploid':'Diploid', 
		'Primary':'Primary', 
		'LAdrenalMet':'Left\nAdrenal\nMetast.',
		'RRib7Met':'Right\nRib7\nMetast.',
		'RIngLNMet':'Right\nIngLN\nMetast.',
		'RSubduralMet':'Right\nSubdural\nMetast.',
		'internal_1':'MRCA\nAll',
		'internal_2':'MRCA\nMetast.',
		'internal_3':'',
		'internal_4':'MRCA\nRight\nside',
		}
	fig = medicc.plot.plot_cn_profiles(
		plotdat,
		tree, 
		title=title, 
		normal_name='diploid', 
		tree_width_scale=1,
		track_width_scale=0.6, 
		height_scale=1, 
		hide_normal_chromosomes=True,
		plot_clonal_summary=True,
		ignore_segment_lengths=False,
		label_func = lambda x:labels[x])
	return fig
fig = plot_ptx011(ptx011_df, ptx011_tree, title='PTX011')
fig.savefig('../examples/output/PTX011.pdf', bbox_inches='tight')
fig = plot_ptx011(ptx011_nowgd_df, ptx011_nowgd_tree, title='PTX011 w/o WGD')
fig.savefig('../examples/output/PTX011-noWGD.pdf', bbox_inches='tight')

# %%
fig = medicc.plot.plot_cn_profiles(
	evo001_df, 
	evo001_tree, 
	title="EVO001", 
	normal_name='diploid', 
	tree_width_scale=1.5,
	track_width_scale=1, 
	height_scale=1, 
	hide_normal_chromosomes=True,
    show_small_segments=True,
    plot_clonal_summary=True,
	label_func = lambda x:x.replace('_',' ').replace('G RLX001', 'G-RLX001\n'))
fig.savefig('../examples/output/EVO001.pdf', bbox_inches='tight')

# %%
