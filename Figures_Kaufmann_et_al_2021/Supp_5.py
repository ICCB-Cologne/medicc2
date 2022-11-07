# %%
import os
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from plotting_params import plotting_params, set_plotting_params, label_axes

set_plotting_params()


with open('data/Supp_5.pickle', 'rb') as f:
    plot_data = pickle.load(f)

results = plot_data['medicc']
results_1_core = results.loc[1]
results_32_core = results.loc[32]
results_1_core['method'] = 'medicc_1'
results_32_core['method'] = 'medicc_32'
results_sitka = plot_data['sitka']
results_sitka['method'] = 'sitka'
results_medalt = plot_data['medalt']
results_medalt['method'] = 'medalt'

results = pd.concat([results_1_core.reset_index(), results_32_core.reset_index(), results_medalt.reset_index(), results_sitka.reset_index()]).set_index(['method', 'leaves', 'total_segs', 'n'])


fig, axs = plt.subplots(ncols=2, nrows=2, sharey=False, sharex=True,
                        figsize=(plotting_params['WIDTH_FULL'], plotting_params['WIDTH_FULL']/plotting_params['ASPECT_RATIO']))

for ax, method, title in zip(axs.flat,
                             ['medicc_1', 'medicc_32', 'medalt', 'sitka'],
                             ['MEDICC2 (single core)', 'MEDICC2 (32 cores)', 'MEDALT', 'Sitka']):
    for cur_segments in [100, 200, 400]:
    # for cur_segments in results.index.get_level_values('total_segs').unique():
        if cur_segments not in results.loc[method].index.get_level_values('total_segs').unique():
            ax.errorbar(x=[], y=[], yerr=[], marker='o', capsize=5, label=cur_segments, lw=3)
        else:
            cur = results.loc[(method, slice(None), cur_segments), 'time_hrs'].groupby('leaves')
            ax.errorbar(x=cur.count().index, y=cur.mean(), yerr=cur.std()/np.sqrt(cur.count()), marker='o', capsize=5, label=cur_segments, lw=3)
    ax.set_title(title, fontsize=plotting_params['FONTSIZE_LARGE'])
    ax.legend(title='Numer of\nsegments')

# axs[0, 1].legend(title='segments per chromosome', bbox_to_anchor=(1, 1))
axs[1, 0].set_xlabel('Number of samples')
axs[1, 1].set_xlabel('Number of samples')

for ax in axs.flat:
    # ax.legend(title='segments per chromosome')
    ax.set_ylabel(f'Time in hours')
    
label_axes(axs)

plt.tight_layout()
fig.savefig('final_figures/Supp_5.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_5.png', bbox_inches='tight', dpi=300)
plt.close()
