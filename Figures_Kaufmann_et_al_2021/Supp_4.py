import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import curve_fit

from plotting_params import plotting_params, set_plotting_params, label_axes, color_palette

set_plotting_params()


with open('data/Supp_4.pickle', 'rb') as f:
    plot_data = pickle.load(f)

results_combined = plot_data
results_combined = results_combined.loc[20]

method_rename = {'euclidean_nj_tree': 'Euclidean NJ',
              'hamming_nj_tree': 'Hamming NJ',
              'euclidean_fastme_tree': 'Euclidean Min. Ev.',
              'hamming_fastme_tree': 'Hamming Min. Ev.',
              'medalt_tree': 'MEDALT',
              'medicc_tree': 'MEDICC2',
              'sitka_tree': 'Sitka'}
results_combined = results_combined.rename(method_rename, axis=1)

results_long = pd.melt(results_combined.reset_index(), id_vars=['rate', 'n'])

leaves = 20
def f(x, p1):
    return p1 * x + leaves * 20

x = results_combined.groupby('rate').mean()['medicc_tree_size'].index
y = results_combined.groupby('rate').mean()['medicc_tree_size'].values
linear_fit = curve_fit(f, x, y)

residuals = y - f(x, linear_fit[0])
ss_res = np.sum(residuals**2)
ss_tot = np.sum((y-np.mean(y))**2)
R_squared = 1 - (ss_res / ss_tot)

fig, axs = plt.subplots(ncols=3, nrows=1,
                        figsize=(plotting_params['WIDTH_FULL'], plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']))


# panel A
ax = axs[0]
sns.violinplot(data=results_long.loc[results_long['variable']=='MEDICC2'], x='rate', y='value', ax=ax)
# sns.stripplot(data=results_long.loc[results_long['variable']=='medicc_tree'], x='rate', y='value', color='black', ax=ax)
ax.set_ylabel('Generalized RF distance')
# ax.set_title("effect of translocations/inversions\non MEDICC2's performance")

# panel B

ax = axs[1]
for col in ['MEDICC2', 'Euclidean Min. Ev.', 'Euclidean NJ', 'Hamming Min. Ev.', 'Hamming NJ', 'MEDALT', 'Sitka']:
    cur = results_combined.groupby('rate')        
    ax.errorbar(x=cur.count().index, y=cur.mean()[col], yerr=(cur.std()/np.sqrt(cur.count()))[col], label=col, capsize=5, marker='o')

ax.legend(facecolor='white', frameon=True, framealpha=0.75)
ax.set_title(f'')
ax.set_ylabel('Generalized RF distance')


# panel C
ax = axs[2]
cur_results = results_combined.copy()
# ax.axhline(leaves * 20, linestyle='--', c='k', label='number of gains/losses')
sns.scatterplot(data=cur_results, x='rate', y='medicc_tree_size', hue='rate', palette=color_palette,
                ax=ax, legend=False)

ax.plot(x, linear_fit[0][0]*x + leaves * 20, '--', color='k',
        label=f'{linear_fit[0][0]:.0f} * x + {leaves * 20}\n$R^2$ = {R_squared:.2f}')
# ax.set_title('')
ax.set_ylabel('Generalized RF distance')
ax.legend()

ax.set_ylabel('number of inferred\ngains/losses by MEDICC2')

# ax.set_title("effect of translocations/inversions\non MEDICC2's performance")
    
# all panels
for ax in axs:
    ax.set_xlabel('Ratio between occurences of\nCN neutral events and gains/losses')

label_axes(axs)

plt.tight_layout()
fig.savefig('final_figures/Supp_4.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_4.png', bbox_inches='tight', dpi=300)
plt.close()
