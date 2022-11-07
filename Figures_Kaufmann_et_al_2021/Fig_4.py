# %%
import os
import sys

import adjustText
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
data_folder = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           "../examples/output_gundem_et_al_2015")
patients = np.sort([f.split('_')[0] for f in os.listdir(data_folder) if 'final_cn_profiles.tsv' in f])

all_changes = {}
WGD_results = pd.DataFrame(index=patients, columns=['nr_wgds', 'type'])
WGD_results['type'] = 'No WGD'

for patient in patients:
    output_df = pd.read_csv(os.path.join(data_folder, "{}_final_cn_profiles.tsv".format(patient)),
                            sep='\t')
    cur_tree = medicc.io.import_tree(os.path.join(data_folder, "{}_final_tree.new".format(patient)),
                                     'diploid', file_format='newick')

    WGD_nodes = np.unique(output_df.loc[output_df['is_wgd'], 'sample_id'])
    WGD_results.loc[patient, 'nr_wgds'] = len(WGD_nodes)
    if len(WGD_nodes):
        if list(cur_tree.find_clades(WGD_nodes[0]))[0] in cur_tree.root.clades:
            WGD_results.loc[patient, 'type'] = 'Clonal'
        elif list(cur_tree.find_clades(WGD_nodes[0]))[0] in cur_tree.get_terminals():
            WGD_results.loc[patient, 'type'] = 'Terminal'
        else:
            WGD_results.loc[patient, 'type'] = 'Subclonal'
          
    all_changes[patient] = output_df
    
WGD_results['type'] = pd.Categorical(WGD_results['type'], ['No WGD','Clonal','Subclonal','Terminal'])

figwidth_panel1 = plotting_params['WIDTH_HALF']/2
plt.figure(figsize=(figwidth_panel1, plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']))
ax = sns.countplot(y='type', data=WGD_results)
for p in ax.patches:
    ax.annotate('{:d}'.format(int(p.get_width())), (p.get_width()-0.5, (p.get_y()+p.get_height()/2)), 
        va='center', 
        c='white', 
        fontsize=plotting_params['FONTSIZE_LARGE'])
ax.set_ylabel('WGD type', rotation=90, labelpad=0)
ax.set_xlabel('Frequency')
# plt.ylabel('WGD type', rotation=0, labelpad=40)
#ax.set_title('WGDs\n(10 patients)', fontsize=plotting_params['FONTSIZE_LARGE'])
plt.tight_layout()
plt.savefig('final_figures/Fig_4A.pdf', pad_inches=0)
plt.savefig('final_figures/Fig_4A.png', pad_inches=0, dpi=600)

#%%
print('plotting Fig 4B')
patients = np.sort([f.split('_')[0] for f in os.listdir(data_folder) if 'final_cn_profiles.tsv' in f])

all_region_overlaps = {}

for patient in patients:
    cur_df = medicc.io.read_and_parse_input_data(
        os.path.join(data_folder, "{}_final_cn_profiles.tsv".format(patient)))
    cur_tree = medicc.io.import_tree(os.path.join(data_folder, "{}_final_tree.new".format(patient)),
                                     'diploid', file_format='newick')

    cur_overlap = medicc.core.overlap_events(output_df=cur_df[['cn_a', 'cn_b']],
                                            tree=cur_tree,
                                            chromosome_bed='../medicc/objects/hg19_chromosome_arms.bed',
                                            regions_bed='../medicc/objects/Davoli_2013_TSG_OG_genes.bed',
                                            replace_loh_with_loss=True,
                                            alleles=['cn_a', 'cn_b'],
                                            replace_both_arms_with_chrom=False)

    cur_overlap['patient'] = patient
    all_region_overlaps[patient] = cur_overlap

combined_events = pd.concat(all_region_overlaps.values()).reset_index().set_index('branch')
regions = np.unique(combined_events['name'])
arms = [x for x in regions if 'chr' in x and ('p' in x or 'q' in x)]
genes = [x for x in regions if 'chr' not in x]
#%%
gain_loss_diff = [((combined_events.loc[combined_events['name'] == arm, 'event'] == 'gain').sum() -
                   (combined_events.loc[combined_events['name'] == arm, 'event'] == 'loss').sum() -
                   (combined_events.loc[combined_events['name'] == arm, 'event'] == 'loh').sum()) for arm in arms]

arm_results = pd.DataFrame(columns=['chrom', 'gain_loss_diff'],
                           index=arms)
arm_results['gain_loss_diff'] = gain_loss_diff
arm_results['chrom'] = arm_results.index.map(lambda x: x.replace('q', '').replace('p', ''))

# Data from Davoli et al. 2013
armwise_scores = pd.DataFrame(columns=['arm', 'score'], 
                             data=[['chr1p',  -2.194424482], ['chr1q',  1.226691224], ['chr2p',  -0.058765864], ['chr2q',  0.722951857], ['chr3p',  -0.23241358], ['chr3q',  2.77507441], ['chr4p',  0.343354369], ['chr4q',  1.236932406], ['chr5p',  1.29642446], ['chr5q',  -0.728377682], ['chr6p',  -0.841818493], ['chr6q',  -0.783217918], ['chr7p',  5.195591398], ['chr7q',  4.588576125], ['chr8p',  1.391378151], ['chr8q',  1.26397449], ['chr9p',  -1.495356436], ['chr9q',  1.979101476], [
                                 'chr10p', -3.665105263], ['chr10q', -0.068322404], ['chr11p', 0.590530233], ['chr11q', -0.200811736], ['chr12p', 2.960242754], ['chr12q', 1.661808926], ['chr13q', -1.906547855], ['chr14q', -0.784903448], ['chr15q', -0.893062731], ['chr16p', 0.450231325], ['chr16q', -3.157515406], ['chr17p', -2.803987461], ['chr17q', 0.168961686], ['chr18p', 0], ['chr18q', -2.868015464], ['chr19p', 1.632886207], ['chr19q', -0.908804749], ['chr20p', 0.602858757], ['chr20q', 1.609560694], ['chr21q', -0.990394366], ['chr22q', -1.209528986]])
armwise_scores = armwise_scores.set_index('arm')

arm_results = arm_results.join(armwise_scores, how='inner')

#%%
plt.figure(figsize=(plotting_params['WIDTH_HALF']-figwidth_panel1/4, plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']))

texts = []
for i, (arm, row) in enumerate(arm_results.iterrows()):
    plt.plot(row['score'], row['gain_loss_diff'], 
             'o', ms=plotting_params['MARKERSIZE_LARGE'], color=plt.cm.tab20c.colors[i%20])
    texts.append(plt.text(row['score'], row['gain_loss_diff']-0.2, arm.replace('chr', ''),
             fontsize=plotting_params['FONTSIZE_SMALL']))
adjustText.adjust_text(texts)
#plt.title('Aggregated for all 10 patients\nPearson R={:.2f} (p={:.0e})'.format(*pearsonr(arm_results['score'], arm_results['gain_loss_diff'])),
#          fontsize=plotting_params['FONTSIZE_LARGE'])

x = np.linspace(arm_results['score'].min(), arm_results['score'].max(), 100)
linear_model = np.polyfit(arm_results['score'], arm_results['gain_loss_diff'], 1)
linear_model_fn = np.poly1d(linear_model)
plt.plot(x, linear_model_fn(x), color="grey", lw=3, label='_nolegend_') # label was: "linear fit"

corr = pearsonr(arm_results['score'], arm_results['gain_loss_diff'])
r = corr[0]
p = corr[1]
print(f'correlation: r = {r}, p = {p}')
plt.text(0.98, 0.02, "$r = {:.2f}$".format(r), transform=plt.gca().transAxes, 
        ha='right', va='bottom', fontsize=plotting_params['FONTSIZE_MEDIUM'])

plt.xlabel('OG-TSG score')
plt.ylabel('#gains - #losses')

plt.tight_layout()
plt.savefig('final_figures/Fig_4B.pdf', pad_inches=0)
plt.savefig('final_figures/Fig_4B.png', pad_inches=0, dpi=600)


#%%
print('plotting Fig 4C')
TSG_genes = pd.read_csv('data/Fig_4C_TSG.csv', header=2, index_col=0, sep=',')
OG_genes = pd.read_csv('data/Fig_4C_OG.csv', header=2, index_col=0, sep=',')

results_Davoli = pd.DataFrame(index=np.intersect1d(genes, list(TSG_genes.index) + list(OG_genes.index)),
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
fig, ax = plt.subplots(figsize=(plotting_params['WIDTH_HALF']-figwidth_panel1/4, plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']))

x1 = -1*np.log10(results_Davoli.loc[np.logical_and(results_Davoli['p-value'] < p_min, results_Davoli['type']=='OG'), 'p-value'].values.astype(float))
y1 = results_Davoli.loc[np.logical_and(results_Davoli['p-value'] < p_min, results_Davoli['type']=='OG'), 'gains-losses']
ax.plot(x1, y1, 'o', color='C2', ms=plotting_params['MARKERSIZE_SMALL'], label='Oncogenes')

x2 = np.log10(results_Davoli.loc[np.logical_and(
    results_Davoli['p-value'] < p_min, results_Davoli['type'] == 'TSG'), 'p-value'].values.astype(float))
y2 = results_Davoli.loc[np.logical_and(
    results_Davoli['p-value'] < p_min, results_Davoli['type'] == 'TSG'), 'gains-losses']
ax.plot(x2, y2, 'o', color='C3', ms=plotting_params['MARKERSIZE_SMALL'], label='Tumour-suppressors')

texts.clear()
for g, row in results_Davoli.loc[results_Davoli['p-value'] < 1e-20].iterrows():
    if row['type'] == 'OG':
        texts.append(ax.text(-1*np.log10(row['p-value']), row['gains-losses']-0.7, g, 
                fontsize=plotting_params['FONTSIZE_SMALL']))
    else:
        texts.append(ax.text(np.log10(row['p-value']), row['gains-losses']-0.7, g, 
                fontsize=plotting_params['FONTSIZE_SMALL']))
adjustText.adjust_text(texts)

ax.axvline(0, color='grey')
ax.set_ylabel('#gains - #losses')

ax.legend(loc='upper left', frameon=False)
#ax.set_title('{} OG-TSG genes\nPearson R = {:.2f} (p={:.0e})'.format(len(x1)+len(x2), 
#                                                                       *pearsonr(np.append(x1, x2),
#                                                                                 np.append(y1, y2))),
#                fontsize=plotting_params['FONTSIZE_LARGE'])

x = np.linspace(np.append(x1, x2).min(), np.append(x1, x2).max(), 100)
linear_model = np.polyfit(np.append(x1, x2), np.append(y1, y2), 1)
linear_model_fn = np.poly1d(linear_model)
ax.plot(x, linear_model_fn(x), color="grey", lw=3, label='linear fit')
    
ax.set_xlabel('signed log10 p-value for OG-TSG score')

corr = pearsonr(np.append(x1, x2), np.append(y1, y2))
r = corr[0]
p = corr[1]
print(f'correlation: r = {r}, p = {p}')
ax.text(0.98, 0.02, "$r = {:.2f}$".format(r), transform=ax.transAxes, 
        ha='right', va='bottom', fontsize=plotting_params['FONTSIZE_MEDIUM'])

plt.tight_layout()
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
