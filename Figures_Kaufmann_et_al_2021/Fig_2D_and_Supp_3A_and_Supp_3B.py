#%% imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from plotting_params import plotting_params, set_plotting_params

set_plotting_params()

# Run the notebook 'pcawg.py' in the folder 'notebooks' to recreate the data from scratch
plotting_data = pd.read_csv('data/Fig_2D_and_Supp_3A_and_Supp_3B.tsv', index_col=None, sep='\t')
pcawg_color_palette = pd.read_csv('data/Supp_3A_and_Supp_3B_color_palette.tsv', index_col=0, sep='\t')
linex = np.linspace(0, plotting_data['hom'].max())
linea = -2
linec = 2.9
liney = linea * linex + linec

#%% Fig 2D: MEDICC predictions
fig, ax = plt.subplots(figsize=(plotting_params['WIDTH_HALF'], plotting_params['WIDTH_HALF']))
sns.scatterplot(x='hom', y='ploidy_pcawg', data=plotting_data, hue='pcawg_wgd', ax=ax, hue_order=['No WGD', 'WGD'])
sns.scatterplot(x='hom', y='ploidy_pcawg', data=plotting_data.loc[plotting_data['bootstrap_2_wgds']],
                color='C2', label='2 WGDs')
sns.scatterplot(x='hom', y='ploidy_pcawg', data=plotting_data.loc[plotting_data.eval('pcawg_wgd != wgd_status_medicc_bootstrap')],
                color='black', label='False Predictions')
sns.scatterplot(x='hom', y='ploidy_pcawg', data=plotting_data.loc[plotting_data.eval('pcawg_wgd != wgd_status_medicc_bootstrap') &
                                                           plotting_data['wgd_uncertain']],
                color='grey', label='False Predictions (uncertain)')
ax.set_title('False Predictions')
ax.set_xlabel('Fraction of genome with LOH')
ax.set_ylabel('Ploidy')
ax.plot(linex, liney, '--', color='grey')
fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Fig_2D.pdf', bbox_inches='tight')
fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Fig_2D.png', bbox_inches='tight', dpi=600)


# #%% Supp 3A: MEDICC MED distributions
palette = dict(zip(pcawg_color_palette.index, pcawg_color_palette['Hexadecimal']))
fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True,
                        figsize=(3*plotting_params['WIDTH_HALF'], plotting_params['WIDTH_FULL']))

cur_plotting_data = plotting_data.query("histology_useable").sort_values('histology_abbreviation')
sns.boxplot(x='dist_wgd', y='histology_abbreviation', palette=palette, data=cur_plotting_data, ax=axs[0])
axs[0].set_xlabel('MED (WGD) from diploid normal')
axs[0].set_ylabel('PCAWG histology')

#%% Supp 3B: MEDICC2 WGD score distributions
cur_plotting_data = plotting_data.query("histology_useable").sort_values('histology_abbreviation')
sns.boxplot(x='dist_diff', y='histology_abbreviation', palette=palette, data=cur_plotting_data, ax=axs[1])
axs[1].set_xlabel('MEDICC2 WGD score')
axs[1].set_ylabel('')
# axs[1].set_yticklabels([])

#%% Supp 3C: MEDICC2 WGD score distributions
plotting_data['has_wgd'] = plotting_data['wgd_status_medicc_bootstrap'] == 'WGD'
cur_plotting_data = plotting_data.query("histology_useable").groupby(
    'histology_abbreviation').mean()[['has_wgd']]

sns.barplot(x='has_wgd', y=cur_plotting_data.index, orient='h', data=cur_plotting_data, 
            palette=palette, ax=axs[2], edgecolor='black', linewidth=1)
axs[2].set_xlabel('Fraction of samples exhibiting WGD')
axs[2].set_ylabel('')
# axs[2].set_yticklabels([])

fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Supp_3A_3B_3C.pdf', bbox_inches='tight')
fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Supp_3A_3B_3C.png',
            bbox_inches='tight', dpi=600)
