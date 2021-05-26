# %%
import os
import sys

import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import medicc

# %%
data_folder = "../examples/output_gundem_et_al_2015"
patients = [f.split('_')[2] for f in os.listdir(data_folder) if 'final_cn_profiles.tsv' in f]

for patient in patients:

    cur_df = medicc.io.read_and_parse_input_data(
        os.path.join(data_folder, "20210303_G_{}_gundem_phased_data_intersection_1mb_homdel_correct_df_final_cn_profiles.tsv".format(patient)))
    cur_tree = medicc.io.import_tree(
        os.path.join(data_folder, "20210303_G_{}_gundem_phased_data_intersection_1mb_homdel_correct_df_final_tree.new".format(patient)), 'diploid')

    labels = {'diploid': 'Diploid'}
    for label in cur_df.reset_index()['sample_id']:
        if 'diploid' not in label and 'internal' not in label:
            labels[label] = '_'.join([label.split('_')[1].split('-')[0], label.split('_')[-1]])

    fig = medicc.plot.plot_cn_profiles(
        cur_df,
        cur_tree,
        title=patient,
        normal_name='diploid',
        hide_normal_chromosomes=False,
        ignore_segment_lengths=False,
        label_func=lambda label: labels.get(label, label))

    for ax in fig.get_axes():
        ax.set_ylabel(ax.get_ylabel(), rotation=0, horizontalalignment='right')

    fig.savefig(os.path.join(data_folder, '{}.pdf'.format(patient)), bbox_inches='tight')

#%%
# Plot the branch-support of the trees. This might take some time
N_bootstrap = 100
for patient in patients:

    cur_tree = medicc.io.import_tree(
        os.path.join(data_folder, "20210303_G_{}_gundem_phased_data_intersection_1mb_homdel_correct_df_final_tree.new".format(patient)), 'diploid')
    cur_df = medicc.io.read_and_parse_input_data(
        os.path.join(data_folder, "20210303_G_{}_gundem_phased_data_intersection_1mb_homdel_correct_df_final_cn_profiles.tsv".format(patient)))
    # remove internal nodes from df
    cur_df = cur_df.loc[cur_df.index.get_level_values('sample_id').map(lambda x: 'internal' not in x)]

    labels = {'diploid': 'Diploid'}
    for label in cur_df.reset_index()['sample_id']:
        if 'diploid' not in label and 'internal' not in label:
            labels[label] = '_'.join([label.split('_')[1].split('-')[0], label.split('_')[-1]])

    trees_df, support_tree = medicc.bootstrap.run_bootstrap(cur_df, cur_tree, 
                                                            N_bootstrap=N_bootstrap, method='chr-wise')

    fig, ax = plt.subplots(figsize=(12, 12))
    fig = medicc.plot.plot_tree(support_tree,
                                title='support tree for {}'.format(patient),
                                label_func=lambda label: labels.get(label, label),
                                show_branch_support=True,
                                show_branch_lengths=True,
                                ax=ax)
    fig.savefig(os.path.join(data_folder, 'support_tree_{}.pdf'.format(patient)), bbox_inches='tight')
