# %%
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr

sys.path.append('..')
import medicc

from plotting_params import plotting_params, set_plotting_params

set_plotting_params()


#%%
print('plotting Fig 4A')
data_folder = "../examples/output_gundem_et_al_2015"
patients = np.sort([f.split('_')[2] for f in os.listdir(data_folder) if 'final_cn_profiles.tsv' in f])

all_changes = {}
WGD_results = pd.DataFrame(index=patients, columns=['nr_wgds', 'type'])
WGD_results['type'] = 'no WGD'

for patient in patients:
    cur_df = medicc.io.read_and_parse_input_data(
        os.path.join(data_folder, "20210303_G_{}_gundem_phased_data_intersection_1mb_homdel_correct_df_final_cn_profiles.tsv".format(patient)))
    cur_tree = medicc.io.import_tree(os.path.join(data_folder, "20210303_G_{}_gundem_phased_data_intersection_1mb_homdel_correct_df_final_tree.new".format(patient)),
                                     'diploid', file_format='newick')

    output_df = medicc.core.summarize_changes(cur_df[['cn_a', 'cn_b']], cur_tree, normal_name='diploid',
                                      allele_columns=['cn_a', 'cn_b'])

    WGD_nodes = np.unique(output_df.loc[output_df['is_wgd']].index.get_level_values('sample_id'))
    WGD_results.loc[patient, 'nr_wgds'] = len(WGD_nodes)
    if len(WGD_nodes):
        if list(cur_tree.find_clades(WGD_nodes[0]))[0] in cur_tree.root.clades:
            WGD_results.loc[patient, 'type'] = 'clonal'
        elif list(cur_tree.find_clades(WGD_nodes[0]))[0] in cur_tree.get_terminals():
            WGD_results.loc[patient, 'type'] = 'terminal'
        else:
            WGD_results.loc[patient, 'type'] = 'sub-clonal'
          
    all_changes[patient] = output_df
    
WGD_results['type'] = pd.Categorical(WGD_results['type'], ['no WGD','clonal','sub-clonal','terminal'])

plt.figure(figsize=(plotting_params['WIDTH_HALF']/2, plotting_params['WIDTH_HALF']))

sns.histplot(y=WGD_results['type'], color='C0', lw=1)
plt.ylabel('WGD type', rotation=90, labelpad=0)
# plt.ylabel('WGD type', rotation=0, labelpad=40)
plt.title('WGDs\n(10 patients)', fontsize=plotting_params['FONTSIZE_LARGE'])
plt.yticks(rotation=45, va='center')
plt.xticks(np.arange(0, 6))
plt.text(5+0.25, 0, '5', fontsize=plotting_params['FONTSIZE_LARGE'], va='center')
plt.text(3+0.25, 1, '3', fontsize=plotting_params['FONTSIZE_LARGE'], va='center')
plt.text(1+0.25, 2, '1', fontsize=plotting_params['FONTSIZE_LARGE'], va='center')
plt.text(1+0.25, 3, '1', fontsize=plotting_params['FONTSIZE_LARGE'], va='center')
plt.xlim(0, 5.75)

plt.tight_layout()
plt.savefig('final_figures/Fig_4A.pdf', pad_inches=0)
plt.savefig('final_figures/Fig_4A.png', pad_inches=0, dpi=600)

#%%
print('plotting Fig 4B')
data_folder = "../examples/output_gundem_et_al_2015"
patients = np.sort([f.split('_')[2] for f in os.listdir(data_folder) if 'final_cn_profiles.tsv' in f])

all_events = {}

for patient in patients:
    cur_df = medicc.io.read_and_parse_input_data(
        os.path.join(data_folder, "20210303_G_{}_gundem_phased_data_intersection_1mb_homdel_correct_df_final_cn_profiles.tsv".format(patient)))
    cur_tree = medicc.io.import_tree(os.path.join(data_folder, "20210303_G_{}_gundem_phased_data_intersection_1mb_homdel_correct_df_final_tree.new".format(patient)),
                                     'diploid', file_format='newick')

    cur_events = medicc.core.overlap_events(df=cur_df[['cn_a', 'cn_b']], 
                                            tree=cur_tree, 
                                            chromosome_bed='../objects/hg19_chromosome_arms.bed',
                                            regions_bed='data/Fig_4C.bed',
                                            replace_loss_with_loh=True)

    all_events[patient] = cur_events

combined_events = pd.concat(all_events.values()).reset_index().set_index('branch')
genes, counts = np.unique(combined_events['name'], return_counts=True)
#%%
gain_loss_diff = [((combined_events.loc[combined_events['name'] == arm, 'event'] == 'gain').sum() -
                   (combined_events.loc[combined_events['name'] == arm, 'event'] == 'loss').sum() -
                   (combined_events.loc[combined_events['name'] == arm, 'event'] == 'loh').sum()) for arm in [x for x in genes if 'chr' in x]]

arm_results = pd.DataFrame(columns=['chrom', 'gain_loss_diff'],
                           index=[x for x in genes if 'chr' in x])
arm_results['gain_loss_diff'] = gain_loss_diff
arm_results['chrom'] = arm_results.index.map(lambda x: x.replace('q', '').replace('p', ''))

arm_results['gain_loss_diff_wgd'] = arm_results['gain_loss_diff']
arm_results['gain_loss_diff_wgd'] += (combined_events['event'] == 'WGD').sum() * 2

for patient, cur_events in all_events.items():
    wgd_nodes = list(cur_events.loc[cur_events['name'] == 'WGD'].index)

    cur_tree = medicc.io.import_tree(os.path.join(data_folder, "20210303_G_{}_gundem_phased_data_intersection_1mb_homdel_correct_df_final_tree.new".format(patient)),
                                     'diploid', file_format='newick')

    for wgd_node in wgd_nodes:
        LOHs = np.unique(cur_events.loc[[x.name for x in cur_tree.get_path(wgd_node)]].loc[cur_events.loc[[
                         x.name for x in cur_tree.get_path(wgd_node)], 'event'] == 'loh', 'name'])
        LOHs = [x for x in LOHs if 'chr' in x]
        arm_results.loc[LOHs, 'gain_loss_diff_wgd'] -= 1

# Data from Davoli et al. 2013
armwise_scores = pd.DataFrame(columns=['arm', 'score'], 
                             data=[['chr1p',  -2.194424482], ['chr1q',  1.226691224], ['chr2p',  -0.058765864], ['chr2q',  0.722951857], ['chr3p',  -0.23241358], ['chr3q',  2.77507441], ['chr4p',  0.343354369], ['chr4q',  1.236932406], ['chr5p',  1.29642446], ['chr5q',  -0.728377682], ['chr6p',  -0.841818493], ['chr6q',  -0.783217918], ['chr7p',  5.195591398], ['chr7q',  4.588576125], ['chr8p',  1.391378151], ['chr8q',  1.26397449], ['chr9p',  -1.495356436], ['chr9q',  1.979101476], [
                                 'chr10p', -3.665105263], ['chr10q', -0.068322404], ['chr11p', 0.590530233], ['chr11q', -0.200811736], ['chr12p', 2.960242754], ['chr12q', 1.661808926], ['chr13q', -1.906547855], ['chr14q', -0.784903448], ['chr15q', -0.893062731], ['chr16p', 0.450231325], ['chr16q', -3.157515406], ['chr17p', -2.803987461], ['chr17q', 0.168961686], ['chr18p', 0], ['chr18q', -2.868015464], ['chr19p', 1.632886207], ['chr19q', -0.908804749], ['chr20p', 0.602858757], ['chr20q', 1.609560694], ['chr21q', -0.990394366], ['chr22q', -1.209528986]])
armwise_scores = armwise_scores.set_index('arm')

arm_results = arm_results.join(armwise_scores, how='inner')

#%%
plt.figure(figsize=(plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']))
sns.set_palette("tab20", plt.cm.tab20c.N)

for arm, row in arm_results.iterrows():
    plt.plot(row['score'], row['gain_loss_diff'], 'o', ms=plotting_params['MARKERSIZE_LARGE'])
    plt.text(row['score'], row['gain_loss_diff']-0.2, arm,
             fontsize=plotting_params['FONTSIZE_TINY'])
plt.title('Aggregated for all 10 patients\nPearson R={:.2f} (p={:.0e})'.format(*pearsonr(arm_results['score'],
                                                                                                              arm_results['gain_loss_diff'])),
          fontsize=plotting_params['FONTSIZE_LARGE'])

x = np.linspace(arm_results['score'].min(), arm_results['score'].max(), 100)
linear_model = np.polyfit(arm_results['score'], arm_results['gain_loss_diff'], 1)
linear_model_fn = np.poly1d(linear_model)
plt.plot(x, linear_model_fn(x), color="grey", lw=3, label='linear fit')

plt.legend()
plt.xlabel('OG-TSG score\n(Davoli 2013)')
plt.ylabel('#gains - #losses')
# plt.grid(False)
# plt.gca().set_facecolor('white')
plt.tight_layout()
plt.savefig('final_figures/Fig_4B.pdf', pad_inches=0)
plt.savefig('final_figures/Fig_4B.png', pad_inches=0, dpi=600)


#%%
print('plotting Fig 4C')
TSG_genes = pd.read_csv('data/Fig_4C_TSG.csv', header=2, index_col=0, sep=',')
OG_genes = pd.read_csv('data/Fig_4C_OG.csv', header=2, index_col=0, sep=',')

results_Davoli = pd.DataFrame(index=np.intersect1d(np.unique(combined_events['name']), list(TSG_genes.index) + list(OG_genes.index)),
                              columns=['type', 'gains', 'losses',
                                       'gains-losses', 'p-value', 'p-value_directed'],
                              dtype=float)


results_Davoli.loc[np.intersect1d(genes, list(TSG_genes.index)), 'type'] = 'TSG'
results_Davoli.loc[np.intersect1d(genes, list(OG_genes.index)), 'type'] = 'OG'
results_Davoli.loc[np.intersect1d(genes, list(
    TSG_genes.index)), 'p-value'] = TSG_genes.loc[np.intersect1d(genes, list(TSG_genes.index)), 'TUSON_p_value_TSG']
results_Davoli.loc[np.intersect1d(genes, list(
    OG_genes.index)), 'p-value'] = OG_genes.loc[np.intersect1d(genes, list(OG_genes.index)), 'TUSON_p_value_OG']
results_Davoli.loc[np.intersect1d(genes, list(
    TSG_genes.index)), 'p-value_directed'] = TSG_genes.loc[np.intersect1d(genes, list(TSG_genes.index)), 'TUSON_p_value_TSG']
results_Davoli.loc[np.intersect1d(genes, list(OG_genes.index)), 'p-value_directed'] = - \
    1*OG_genes.loc[np.intersect1d(genes, list(OG_genes.index)), 'TUSON_p_value_OG']
for gene in results_Davoli.index:
    results_Davoli.loc[gene, 'gains'] = len(combined_events.loc[np.logical_and(
        combined_events['name'] == gene,  combined_events['event'] == 'gain')])
    results_Davoli.loc[gene, 'losses'] = len(combined_events.loc[np.logical_and(
        combined_events['name'] == gene,  combined_events['event'].isin(['loss', 'loh']))])

results_Davoli['gains-losses'] = results_Davoli['gains'] - results_Davoli['losses']

#%%
p_min = 1
sns.set_palette("tab10",plt.cm.tab10.N)
fig, ax = plt.subplots(figsize=(plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']))

x1 = -1*np.log10(results_Davoli.loc[np.logical_and(results_Davoli['p-value'] < p_min, results_Davoli['type']=='OG'), 'p-value'].values.astype(float))
y1 = results_Davoli.loc[np.logical_and(results_Davoli['p-value'] < p_min, results_Davoli['type']=='OG'), 'gains-losses']
ax.plot(x1, y1, 'o', color='C2', ms=plotting_params['MARKERSIZE_SMALL'], label='Oncogenes')
x2 = np.log10(results_Davoli.loc[np.logical_and(
    results_Davoli['p-value'] < p_min, results_Davoli['type'] == 'TSG'), 'p-value'].values.astype(float))
y2 = results_Davoli.loc[np.logical_and(
    results_Davoli['p-value'] < p_min, results_Davoli['type'] == 'TSG'), 'gains-losses']
ax.plot(x2, y2, 'o', color='C3', ms=plotting_params['MARKERSIZE_SMALL'], label='Tumor-Suppressors')

for g, row in results_Davoli.loc[results_Davoli['p-value'] < 1e-20].iterrows():
    if row['type'] == 'OG':
        ax.text(-1*np.log10(row['p-value']), row['gains-losses']-0.7, g, 
                fontsize=plotting_params['FONTSIZE_TINY'])
    else:
        ax.text(np.log10(row['p-value']), row['gains-losses']-0.7, g, 
                fontsize=plotting_params['FONTSIZE_TINY'])

ax.axvline(0, color='grey')
ax.set_ylabel('#gains - #losses')

ax.legend(loc='upper left')
ax.set_title('{} TSG / OG genes\nPearson R = {:.2f} (p={:.0e})'.format(len(x1)+len(x2), 
                                                                       *pearsonr(np.append(x1, x2),
                                                                                 np.append(y1, y2))),
                fontsize=plotting_params['FONTSIZE_LARGE'])

x = np.linspace(np.append(x1, x2).min(), np.append(x1, x2).max(), 100)
linear_model = np.polyfit(np.append(x1, x2), np.append(y1, y2), 1)
linear_model_fn = np.poly1d(linear_model)
ax.plot(x, linear_model_fn(x), color="grey", lw=3, label='linear fit')
    
ax.set_xlabel('log10 of p-value for TSG/OG score \n(Davoli 2013)')

plt.savefig('final_figures/Fig_4C.pdf', pad_inches=0)
plt.savefig('final_figures/Fig_4C.png', pad_inches=0, dpi=600)

#%% Check correlation for filtered genes
for p_min in [1, 1e-5, 1e-10, 1e-15, 1e-20]:
    x1 = -1*np.log10(results_Davoli.loc[np.logical_and(results_Davoli['p-value']
                                                    < p_min, results_Davoli['type'] == 'OG'), 'p-value'].values.astype(float))
    y1 = results_Davoli.loc[np.logical_and(
        results_Davoli['p-value'] < p_min, results_Davoli['type'] == 'OG'), 'gains-losses']
    x2 = np.log10(results_Davoli.loc[np.logical_and(
        results_Davoli['p-value'] < p_min, results_Davoli['type'] == 'TSG'), 'p-value'].values.astype(float))
    y2 = results_Davoli.loc[np.logical_and(
        results_Davoli['p-value'] < p_min, results_Davoli['type'] == 'TSG'), 'gains-losses']

    print(p_min, len(x1)+len(x2), pearsonr(np.append(x1, x2), np.append(y1, y2)))

for top_N in [50, 100, 500]:
    cur_results = results_Davoli.sort_values('p-value').iloc[:top_N]
    x1 = -1*np.log10(cur_results.loc[cur_results['type'] == 'OG', 'p-value'].values.astype(float))
    y1 = cur_results.loc[cur_results['type'] == 'OG', 'gains-losses']
    x2 = np.log10(cur_results.loc[cur_results['type'] == 'TSG', 'p-value'].values.astype(float))
    y2 = cur_results.loc[cur_results['type'] == 'TSG', 'gains-losses']

    print(len(x1)+len(x2), pearsonr(np.append(x1, x2), np.append(y1, y2)))
