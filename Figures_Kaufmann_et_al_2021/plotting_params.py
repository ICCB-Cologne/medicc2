import matplotlib.pyplot as plt
import seaborn as sns

plotting_params = {
    'WIDTH_FULL': 12,
    'WIDTH_HALF': 6,
    'ASPECT_RATIO': 4/3,
    'FONTSIZE_LARGE': 12,
    'FONTSIZE_MEDIUM': 10,
    'FONTSIZE_SMALL': 8,
    'FONTSIZE_TINY': 5,
    'LINEWIDTH': 3,
    'MARKERSIZE_SMALL': 3,
    'MARKERSIZE_MEDIUM': 5,
    'MARKERSIZE_LARGE': 10,
    'LINEWIDTH_SMALL': 1}


def set_plotting_params():
    plt.rc('font', family='sans-serif')
    plt.rc('font', size=plotting_params['FONTSIZE_MEDIUM'])
    plt.rc('axes', titlesize=plotting_params['FONTSIZE_LARGE'])
    plt.rc('axes', labelsize=plotting_params['FONTSIZE_MEDIUM'])
    plt.rc('xtick', labelsize=plotting_params['FONTSIZE_SMALL'])
    plt.rc('ytick', labelsize=plotting_params['FONTSIZE_SMALL'])
    plt.rc('legend', fontsize=plotting_params['FONTSIZE_SMALL'])
    plt.rc('legend', frameon=False)
    plt.rc('figure', titlesize=plotting_params['FONTSIZE_LARGE'])

    # Seaborn standard palette
    sns.set_palette(sns.color_palette())
    # sns.set_palette(sns.color_palette("husl"))

    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
