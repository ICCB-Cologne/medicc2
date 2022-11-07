import pickle
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import auc

from plotting_params import plotting_params, set_plotting_params, label_axes

set_plotting_params()

with open('data/Supp_11.pickle', 'rb') as f:
    plot_data = pickle.load(f)

results = plot_data

ths = np.arange(101)
accuracy = np.zeros(101)
precision = np.zeros(101)
recall = np.zeros(101)
f1 = np.zeros(101)

baseline = results['pcawg_wgd']=='WGD'
for i in range(101):
    prediction = results['nr_bootstraps_with_wgd'] >= i
    accuracy[i] = (prediction == baseline).sum() / len(baseline)
    precision[i] = np.logical_and(prediction, baseline).sum() / prediction.sum()
    recall[i] = np.logical_and(prediction, baseline).sum() / baseline.sum()
    f1[i] = 2*precision[i] * recall[i] / (precision[i] + recall[i])


fig, axs = plt.subplots(ncols=2,
                        figsize=(plotting_params['WIDTH_FULL'], plotting_params['WIDTH_HALF']/plotting_params['ASPECT_RATIO']))

axs[0].axhline(1, c='black', linestyle='--', zorder=0)
axs[0].plot(ths, accuracy, '-', label='accuracy')
axs[0].plot(ths, precision, '-', label='precision')
axs[0].plot(ths, recall, '-', label='recall')
axs[0].plot(ths, f1, '-', label='f1')
axs[0].legend(title='metric')
axs[0].set_xlim(-0.5, 100)
# axs[0].set_xlim(0.5, 100)

selection = np.arange(1, 11)
inlet_ax = axs[0].inset_axes([0.2, 0.07, 0.47, 0.47])
inlet_ax.axhline(1, c='black', linestyle='--')
inlet_ax.plot(ths[selection], accuracy[selection], 'o-', label='accuracy')
inlet_ax.plot(ths[selection], precision[selection], 'o-', label='precision')
inlet_ax.plot(ths[selection], recall[selection], 'o-', label='recall')
inlet_ax.plot(ths[selection], f1[selection], 'o-', label='f1-score')
inlet_ax.spines['top'].set_visible(True)
inlet_ax.spines['right'].set_visible(True)
# inlet_ax.legend()
inlet_ax.set_xticks(selection[::2])
axs[0].indicate_inset_zoom(inlet_ax, edgecolor="black", alpha=1)
axs[0].set_xlabel('Threshold for bootstrap WGD detection (out of 100 runs)')

axs[1].step(np.append(recall, 0), np.append(precision, 1), '-', where='mid', color='black',
            label=f'Precision-recall curve\nAUC = {auc(np.append(recall, 0), np.append(precision, 1)):.2f}')
axs[1].step(recall, precision, 'o', where='mid', color='black')

axs[1].set_ylabel('Precision')
axs[1].set_xlabel('Recall')
axs[1].legend()

label_axes(axs)

plt.tight_layout()
fig.savefig('final_figures/Supp_11.pdf', bbox_inches='tight')
fig.savefig('final_figures/Supp_11.png', bbox_inches='tight', dpi=300, facecolor='white')
plt.close()
