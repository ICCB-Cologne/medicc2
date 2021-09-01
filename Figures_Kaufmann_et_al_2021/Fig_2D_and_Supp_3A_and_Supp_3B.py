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
# %% Figure with false MEDICC predictions
# plotting_data[['hom', 'ploidy_pcawg', 'pcawg_wgd']].to_csv('../Figures_Kaufmann_et_al_2021/data/Fig_2D.tsv',
#                                                     sep='\t', index=False)

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


# %% plot PCAWG MED distributions
palette = dict(zip(pcawg_color_palette.index, pcawg_color_palette['Hexadecimal']))
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(
    plotting_params['WIDTH_HALF'], plotting_params['WIDTH_FULL']))
cur_plotting_data = plotting_data.query("histology_useable").sort_values('histology_abbreviation')
sns.boxplot(x='dist_wgd', y='histology_abbreviation', palette=palette, data=cur_plotting_data, ax=ax)
ax.set_xlabel('MED (WGD) from diploid normal')
ax.set_ylabel('PCAWG histology')
fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Supp_3A.pdf', bbox_inches='tight')
fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Supp_3A.png',
            bbox_inches='tight', dpi=600)

#%%
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(
    plotting_params['WIDTH_HALF'], plotting_params['WIDTH_FULL']))
cur_plotting_data = plotting_data.query("histology_useable").sort_values('histology_abbreviation')
sns.boxplot(x='dist_diff', y='histology_abbreviation', palette=palette, data=cur_plotting_data, ax=ax)
ax.set_xlabel('MEDICC2 WGD score')
ax.set_ylabel('PCAWG histology')
fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Supp_3B.pdf', bbox_inches='tight')
fig.savefig('../Figures_Kaufmann_et_al_2021/final_figures/Supp_3B.png',
            bbox_inches='tight', dpi=600)
