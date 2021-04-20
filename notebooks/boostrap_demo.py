#%% import and load data
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import medicc
import medicc.bootstrap

# %matplotlib inline

#%%
# By setting the numpy seed, the created bootstrap datasets will always be the same
# Comment this line out if you want to create different datasets each time the notebook is run
np.random.seed(42)

#%%
data_folder = "../examples/output_gundem_et_al_2015"

main_tree = medicc.io.import_tree(os.path.join(data_folder, "20210303_G_PTX011_gundem_phased_data_intersection_1mb_homdel_correct_df_final_tree.new"), 'diploid')
main_df = medicc.io.read_tsv_as_dataframe(
    os.path.join(data_folder, "20210303_G_PTX011_gundem_phased_data_intersection_1mb_homdel_correct_df_final_cn_profiles.tsv"))
# remove internal nodes from df
main_df = main_df.loc[main_df.index.get_level_values('sample_id').map(lambda x: 'internal' not in x)]

labels = {'diploid': 'Diploid'}
for label in main_df.reset_index()['sample_id']:
    if 'diploid' not in label and 'internal' not in label:
        labels[label] = '_'.join([label.split('_')[1].split('-')[0], label.split('_')[-1]])

#%% Run chromosome-wise bootstrap
N_bootstrap = 20

bootstrap_trees_df, support_tree = medicc.bootstrap.run_bootstrap(main_df, main_tree,
                                                                  N_bootstrap=N_bootstrap, method='chr-wise')

print('{} chromosome-wise bootstrap datasets created'.format(N_bootstrap))
print('{} individual trees'.format(len(bootstrap_trees_df)))
print('The original tree was recreated {} times ({:.1f}%)'.format(bootstrap_trees_df.loc[bootstrap_trees_df['note'] == 'original', 'count'].iloc[0],
                                                                  100*bootstrap_trees_df.loc[bootstrap_trees_df['note'] == 'original', 'freq'].iloc[0]))

#%% Show the most prominent trees
fig, axs = plt.subplots(figsize=(21, 7), ncols=3)
for ax, (i, row) in zip(axs, bootstrap_trees_df.iterrows()):
    medicc.plot.plot_tree(row['tree'], 
                          label_func=lambda label: labels.get(label, label),
                          ax=ax,
                          title='{} ({:.0f}%)'.format(row['note'], 100*row['freq']))

fig.suptitle('Most common bootstrap trees', y=1.0, fontsize=25, weight='bold')
plt.show()

#%% Show the support for the individual branches of the main tree
fig, ax = plt.subplots(figsize=(12, 12))
fig = medicc.plot.plot_tree(support_tree,
                            title='support tree',
                            label_func=lambda label: labels.get(label, label),
                            show_branch_lengths=True,
                            show_tree_support=True,
                            ax=ax)
plt.show()


#%% Run segment-wise jackknife
N_bootstrap = 20

bootstrap_trees_df, support_tree = medicc.bootstrap.run_bootstrap(main_df, main_tree,
                                                                  N_bootstrap=N_bootstrap, method='segment-wise')

print('{} segment-wise jackknife datasets created'.format(N_bootstrap))
print('{} individual trees'.format(len(bootstrap_trees_df)))
print('The original tree was recreated {} times ({:.1f}%)'.format(bootstrap_trees_df.loc[bootstrap_trees_df['note'] == 'original', 'count'].iloc[0],
                                                                  100*bootstrap_trees_df.loc[bootstrap_trees_df['note'] == 'original', 'freq'].iloc[0]))
