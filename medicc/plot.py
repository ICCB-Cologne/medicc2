import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from medicc import core

COL_ALLELE_A = mpl.colors.to_rgba('orange')
COL_ALLELE_B = mpl.colors.to_rgba('teal')
COL_CLONAL = mpl.colors.to_rgba('lightgrey')
COL_NORMAL = mpl.colors.to_rgba('dimgray')
COL_GAIN = mpl.colors.to_rgba('red')
COL_LOSS = mpl.colors.to_rgba('blue')
COL_CHR_LABEL = mpl.colors.to_rgba('grey')
COL_VLINES = '#1f77b4'
COL_MARKER_INTERNAL = COL_VLINES
COL_MARKER_TERMINAL = 'black'
COL_MARKER_NORMAL = 'green'
COL_SUMMARY_LABEL = 'grey'
COL_BACKGROUND = 'white'
COL_BACKGROUND_HATCH = 'lightgray'
COL_PATCH_BACKGROUND = 'white'
ALPHA_PATCHES = 0.15
ALPHA_CLONAL = 0.3
BACKGROUND_HATCH_MARKER = '/////'
TREE_MARKER_SIZE = 40
YLABEL_FONT_SIZE = 8
YLABEL_TICK_SIZE = 6
XLABEL_FONT_SIZE = 10
XLABEL_TICK_SIZE = 8
CHR_LABEL_SIZE = 8


def plot_cn_profiles(
        input_df,
        input_tree=None,
        title=None,
        normal_name=None,
        mincn='auto',
        maxcn='auto',
        plot_summary=True,
        plot_subclonal_summary=True,
        plot_mrca_summary=True,
        hide_normal_chromosomes=False,
        ignore_segment_lengths=False, 
        tree_width_scale=1,
        track_width_scale=1, 
        height_scale=1,
        close_gaps=False,
        label_func = None):

    alleles = input_df.columns

    if input_tree is None or normal_name is None: 
        plot_summary = False
        plot_subclonal_summary = False

    if np.setdiff1d(['is_clonal', 'is_normal', 'is_gain', 'is_loss'], input_df.columns).size > 0:
        df = core.summarize_changes(input_df,
                                       input_tree,
                                       normal_name,
                                       ignore_segment_lengths=ignore_segment_lengths)

    if hide_normal_chromosomes:
        df = df.join(df.groupby('chrom')['is_normal'].all().to_frame('hide'))
        df = df.query("~hide").drop('hide', axis=1)

    if plot_summary or plot_subclonal_summary:
        agg_events = core.compute_change_events(df[alleles], input_tree)
        agg_events = agg_events.groupby(["chrom", "start", "end"], observed=True).sum()

    if mincn=='auto':
        mincn = df.min().min()
    else:
        mincn = int(mincn)        
    if maxcn=='auto':
        maxcn = df.max().max()
    else:
        maxcn = int(maxcn)

    samples = df.index.get_level_values('sample_id').unique()
    nsamp = len(samples)
    nsegs = df.loc[samples[0],:].groupby('chrom').size()

    df.reset_index(['start','end'], inplace=True)

    if close_gaps:
        cur_df = df.loc[samples[0], :].reset_index()[['chrom', 'start', 'end']]
        segment_lengths = cur_df['end'] - cur_df['start']
        cur_df['start_pos'] = np.cumsum(np.append([0], segment_lengths))[:-1]
        cur_df['end_pos'] = cur_df['start_pos'] + segment_lengths

        cur_df.set_index(['chrom', 'start'], inplace=True)
        cur_df.drop('end', axis=1, inplace=True)

        df = df.join(cur_df, on=['chrom', 'start'])

    else:
        chrom_offset = df.loc[samples[0],:].reset_index().groupby('chrom', sort=False).max()['end']
        chrom_offset.dropna(inplace=True)
        chrom_offset[:] = np.append(0, chrom_offset.cumsum().values[:-1])
        chrom_offset.name='offset'
        df = df.join(chrom_offset, on='chrom')
        df['start_pos'] = df['start'] + df['offset']
        df['end_pos'] = df['end'] + df['offset'] + 1

    ## determine clade colors
    clade_colors = {}
    for sample in samples:
        ## determine if sample is terminal
        is_terminal = True
        if input_tree is not None:
            matches = list(input_tree.find_clades(sample))
            if len(matches)>0:
                clade = matches[0]
                is_terminal = clade.is_terminal()
        ## determine if sample is normal
        is_normal = sample==normal_name
        clade_colors[sample] = COL_MARKER_TERMINAL
        if not is_terminal:
            clade_colors[sample] = COL_MARKER_INTERNAL
        if is_normal:
            clade_colors[sample] = COL_MARKER_NORMAL

    # %% define plot width and height and create figure
    track_width = nsegs.sum() * 0.2 * track_width_scale
    if input_tree is not None:
        tree_width = 2.5 * tree_width_scale ## in figsize
    else:
        tree_width = 0
    #plotheight =  abs(maxcn-mincn) * 0.2 * nsamp * height_scale
    plotheight =  4 * 0.2 * nsamp * height_scale
    plotwidth = tree_width + track_width
    tree_width_ratio = tree_width / plotwidth
    fig = plt.figure(figsize=(plotwidth, plotheight), constrained_layout=True)
    if input_tree is None:
        nrows = nsamp
        gs = fig.add_gridspec(nrows, 1)
        cn_axes = [fig.add_subplot(gs[i]) for i in range(0, nrows)]
        y_order = list(samples) ## as they appear
    else:
        nrows = nsamp + int(plot_summary) + int(plot_subclonal_summary) + int(plot_mrca_summary)
        gs = fig.add_gridspec(nrows, 2, width_ratios=[tree_width_ratio, 1-tree_width_ratio])
        tree_ax = fig.add_subplot(gs[0:nsamp, 0])
        cn_axes = [fig.add_subplot(gs[i]) for i in range(1,(2*(nrows))+1,2)]
        y_posns = _get_y_positions(input_tree, adjust=True)
        y_order = [x.name for x in y_posns if x.name is not None and x.name!='root'] ## as in tree
        plot_tree(input_tree, 
                  ax=tree_ax,
                  title=title,
                  label_func=lambda x: '',
                  label_colors=clade_colors,
                  branch_labels=lambda x: x.branch_length if x.name != 'root' and x.name is not None else None)
    
    # wspace=-0.03 corresponds ok with "internal_X" labels 
    fig.set_constrained_layout_pads(w_pad=0, h_pad=0, hspace=0.0, wspace=-0.03)
    ## iterate over samples and plot the track
    for sample, group in df.groupby('sample_id'):
        index_to_plot = y_order.index(sample)
        plot_axis_labels = (index_to_plot == (nsamp-1))
        _plot_cn_profile_for_sample(cn_axes[index_to_plot], 
            label_func(sample) if label_func is not None else sample, 
            group, 
            mincn-1, 
            maxcn+1,
            alleles,
            plot_xaxis_labels=plot_axis_labels if not (
                plot_summary + plot_subclonal_summary + plot_mrca_summary) else False,
            plot_yaxis_labels=True,
            yaxis_label_color=clade_colors[sample])
    
    if plot_mrca_summary:
        mrca = [x for x in input_tree.root.clades if x.name != normal_name][0].name    
        mrca_df = df.loc[df.index.get_level_values('sample_id') == mrca] - 1
        _plot_aggregated_events(mrca_df,
                                alleles, cn_axes[nsamp], 
                                close_gaps=close_gaps)
        cn_axes[nsamp].get_xaxis().set_visible(not (plot_summary or plot_subclonal_summary))
        cn_axes[nsamp].set_ylabel('MRCA')
    if plot_summary:
        _plot_aggregated_events(agg_events, alleles, cn_axes[-1], 
                                close_gaps=close_gaps)

    if plot_subclonal_summary:
        agg_events.loc[df.loc[df.index.get_level_values('sample_id') == 'diploid', 'is_clonal'].values] = 0
        _plot_aggregated_events(agg_events, alleles, cn_axes[-1 - int(plot_summary)], 
                                close_gaps=close_gaps)
        cn_axes[-1 - int(plot_summary)].get_xaxis().set_visible(not plot_summary)
        cn_axes[-1 - int(plot_summary)].set_ylabel('subclonal\nsummary')

    cn_axes[-1].set_xlabel('genome position', fontsize=XLABEL_FONT_SIZE)
    cn_axes[-1].xaxis.set_tick_params(labelsize=XLABEL_TICK_SIZE)
    cn_axes[-1].xaxis.offsetText.set_fontsize(XLABEL_TICK_SIZE)

    return fig


def _plot_cn_profile_for_sample(ax, sample_label, group, mincn, maxcn, alleles,
                                plot_xaxis_labels=True, plot_yaxis_labels=True, 
                                yaxis_label_color='black'):
    ## collect line segments and background patches
    event_patches = []
    bkg_patches = []
    lines_a = []
    lines_b = []
    alpha = []
    for idx, r in group.iterrows():
        lines_a.append([(r['start_pos'], r[alleles[0]]),(r['end_pos'], r[alleles[0]])])
        lines_b.append([(r['start_pos'], r[alleles[1]]),(r['end_pos'], r[alleles[1]])])
        rect = mpl.patches.Rectangle((r['start_pos'], mincn), r['end_pos']-r['start_pos'], maxcn-mincn, edgecolor=None, facecolor=COL_PATCH_BACKGROUND, alpha=1)
        alpha.append(ALPHA_CLONAL if r['is_clonal'] else 1.0)
        bkg_patches.append(rect)
        # Not used because clonal tracks are made transparent below
        # if r['is_clonal']:
        #     rect = mpl.patches.Rectangle((r['start_pos'], mincn), r['end_pos']-r['start_pos'], maxcn-mincn, edgecolor=None, facecolor=COL_CLONAL, alpha=0.1)
        #     event_patches.append(rect)
        if r['is_normal']:
            rect = mpl.patches.Rectangle((r['start_pos'], mincn), r['end_pos']-r['start_pos'], maxcn-mincn, edgecolor=None, facecolor=COL_NORMAL, alpha=ALPHA_PATCHES)
            event_patches.append(rect)
        if r['is_gain']:
            rect = mpl.patches.Rectangle((r['start_pos'], mincn), r['end_pos']-r['start_pos'], maxcn-mincn, edgecolor=None, facecolor=COL_GAIN, alpha=ALPHA_PATCHES)
            event_patches.append(rect)
        if r['is_loss']:
            rect = mpl.patches.Rectangle((r['start_pos'], mincn), r['end_pos']-r['start_pos'], maxcn-mincn, edgecolor=None, facecolor=COL_LOSS, alpha=ALPHA_PATCHES)
            event_patches.append(rect)

    #rc = mpl.collections.PatchCollection(rectangles, facecolors=rect_colors, alpha=0.1)
    events = mpl.collections.PatchCollection(event_patches, match_original=True, zorder=2)
    ax.add_collection(events)
    backgrounds = mpl.collections.PatchCollection(bkg_patches, match_original=True, zorder=1)
    ax.add_collection(backgrounds)
    plot_bkg = mpl.patches.Rectangle((0,0), 1, 1, transform=ax.transAxes, facecolor=COL_BACKGROUND, edgecolor=COL_BACKGROUND_HATCH, zorder=0, hatch=BACKGROUND_HATCH_MARKER)
    ax.add_patch(plot_bkg)
    colors_a = np.array([COL_ALLELE_A] * len(lines_a))
    colors_b = np.array([COL_ALLELE_B] * len(lines_b))
    # clonal mutations
    colors_a[:, 3] = np.array(alpha)
    colors_b[:, 3] = np.array(alpha)
    # a and b are overlapping
    colors_a[group[alleles[0]]==group[alleles[1]], 3] = 0.5
    colors_b[group[alleles[0]]==group[alleles[1]], 3] = 0.5
    colors = np.row_stack([colors_a, colors_b])
    lc = mpl.collections.LineCollection(lines_a + lines_b, colors=colors, linewidth=2)
    ax.add_collection(lc)
    ax.autoscale()

    ## draw segment boundaries
    seg_bound_first = group['start_pos'].values[0]
    seg_bounds = group['end_pos'].values
    ax.vlines(np.append(seg_bound_first, seg_bounds), ymin=mincn, ymax=maxcn, ls='--', alpha=0.25, color=COL_VLINES, linewidth=0.5)
    ## draw chromosome boundaries
    chr_ends = group.reset_index().groupby('chrom', sort=False).max()['end_pos']
    linex = chr_ends.values[:-1] ## don't plot last
    ax.vlines(linex, ymin=mincn, ymax=maxcn, color=COL_VLINES, linewidth=1)
    ## draw chromosome labels
    chr_label_pos = chr_ends
    chr_label_pos.loc[:] = np.roll(chr_label_pos.values, 1)
    chr_label_pos.iloc[0] = 0
    for chrom, pos in chr_label_pos.iteritems():
        ax.text(pos + 5e5, maxcn-0.35, chrom, va='top', color=COL_CHR_LABEL,
                fontweight='medium', fontsize=CHR_LABEL_SIZE)
    ## draw sample labels
    if plot_yaxis_labels:
        ax.set_ylabel(sample_label, fontsize=YLABEL_FONT_SIZE, rotation=0, ha='right', va='center')
        ax.yaxis.label.set_color(yaxis_label_color)

    ## axis modifications
    ax.get_xaxis().set_visible(plot_xaxis_labels)

    #ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True, prune='both', nbins=(maxcn-mincn)/2 + 1))
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True, prune='both', nbins=5))
    ax.yaxis.set_tick_params(labelsize=YLABEL_TICK_SIZE)
    ax.set_ylim(mincn, maxcn)
    ax.set_xlim(1, group['end_pos'].max())


def _plot_aggregated_events(agg_events_input, alleles, ax, close_gaps=False):

    agg_events = agg_events_input.copy()
    
    maxcn = agg_events[[alleles[0], alleles[1]]].max().max()+1
    mincn = agg_events[[alleles[0], alleles[1]]].min().min()-1
    
    if 'start_pos' not in agg_events.columns or 'end_pos' not in agg_events.columns:

        if close_gaps:
            cur_df = agg_events.reset_index()[['chrom', 'start', 'end']]
            segment_lengths = cur_df['end'] - cur_df['start']
            cur_df['start_pos'] = np.cumsum(np.append([0], segment_lengths))[:-1]
            cur_df['end_pos'] = cur_df['start_pos'] + segment_lengths

            cur_df.set_index(['chrom', 'start'], inplace=True)
            cur_df.drop('end', axis=1, inplace=True)

            agg_events = agg_events.join(cur_df, on=['chrom', 'start'])
        
        else:
            agg_events.reset_index(['start','end'], inplace=True)
            offset = agg_events.end.groupby('chrom').max()
            offset[:] = np.append(0, offset.cumsum().values[:-1])
            offset.name = 'offset'
            agg_events = agg_events.join(offset)
            agg_events['start_pos'] = agg_events['start'] + agg_events['offset']
            agg_events['end_pos'] = agg_events['end'] + agg_events['offset'] + 1

    # draw ractangles
    event_patches = []
    bkg_patches = []
    rect_colors = []
    lines_a = []
    lines_b = []
    
    for idx, r in agg_events.iterrows():
        rect = mpl.patches.Rectangle((r['start_pos'], mincn), r['end_pos']-r['start_pos'], maxcn-mincn, edgecolor=None, facecolor=COL_PATCH_BACKGROUND, alpha=1)
        bkg_patches.append(rect)
        lines_a.append([(r['start_pos'], r[alleles[0]]),(r['end_pos'], r[alleles[0]])])
        lines_b.append([(r['start_pos'], r[alleles[1]]),(r['end_pos'], r[alleles[1]])])
        
        rect = mpl.patches.Rectangle((r['start_pos'], mincn), 
                                     width = r['end_pos']-r['start_pos'], 
                                     height = maxcn - mincn, 
                                     facecolor ='grey', 
                                     edgecolor = None)
        event_patches.append(rect)
        rect_colors.append(COL_CLONAL)

    
    events = mpl.collections.PatchCollection(event_patches, facecolors = rect_colors, alpha = 0.1, zorder=2)
    ax.add_collection(events)
    backgrounds = mpl.collections.PatchCollection(bkg_patches, match_original=True, zorder=1)
    ax.add_collection(backgrounds)
    plot_bkg = mpl.patches.Rectangle((0,0), 1, 1, transform=ax.transAxes, facecolor=COL_BACKGROUND, edgecolor=COL_BACKGROUND_HATCH, hatch=BACKGROUND_HATCH_MARKER, zorder=0)
    ax.add_patch(plot_bkg)
    colors_a = np.array([COL_ALLELE_A] * len(lines_a))
    colors_b = np.array([COL_ALLELE_B] * len(lines_b))
    colors_a[agg_events[alleles[0]] == agg_events[alleles[1]], 3] = 0.5
    colors_b[agg_events[alleles[0]] == agg_events[alleles[1]], 3] = 0.5
    colors = np.row_stack([colors_a, colors_b])
    lc = mpl.collections.LineCollection(lines_a + lines_b, colors=colors)
    ax.add_collection(lc)
    ax.autoscale()
    
    ## draw segment boundaries
    seg_bound_first = agg_events['start_pos'].values[0]
    seg_bounds = agg_events['end_pos'].values
    ax.vlines(np.append(seg_bound_first, seg_bounds), 
              ymin = mincn, 
              ymax = maxcn, 
              ls = '--', 
              alpha = 0.25, 
              color = COL_VLINES, 
              linewidth = 0.5)
    
    ## draw chromosome boundaries
    chr_ends = agg_events.groupby('chrom').max()['end_pos']
    linex = chr_ends.values[:-1] ## don't plot last
    ax.vlines(linex, 
              ymin = mincn, 
              ymax = maxcn,
              color = COL_VLINES,
              linewidth = 1)
    
    ## draw chromosome labels
    chr_label_pos = chr_ends
    chr_label_pos.loc[:] = np.roll(chr_label_pos.values, 1)
    chr_label_pos.iloc[0] = 0
    for chrom, pos in chr_label_pos.iteritems():
        ax.text(pos + 5e5, maxcn-0.35, chrom, va='top', color=COL_CHR_LABEL,
                fontweight='medium', fontsize=CHR_LABEL_SIZE)

    ## axis and axis labels
    ax.set_ylabel("total\nsummary", fontsize=YLABEL_FONT_SIZE, rotation=0, ha='right', va='center')
    ax.yaxis.label.set_color(COL_SUMMARY_LABEL)
    
    #nbins = (agg_events[alleles].max().max()-agg_events[alleles].min().min()) / 2 + 1
    #ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True, prune='both', nbins=nbins))
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True, prune='both', nbins=5))
    ax.yaxis.set_tick_params(labelsize=YLABEL_TICK_SIZE)

    ax.set_ylim(mincn, maxcn)
    ax.set_xlim(1, agg_events['end_pos'].max())

def _get_x_positions(tree):
    """Create a mapping of each clade to its horizontal position.
    Dict of {clade: x-coord}
    """
    depths = tree.depths()
    # If there are no branch lengths, assume unit branch lengths
    if not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)
    return depths

def _get_y_positions(tree, adjust=False):
    """Create a mapping of each clade to its vertical position.
    Dict of {clade: y-coord}.
    Coordinates are negative, and integers for tips.
    """
    maxheight = tree.count_terminals()
    # Rows are defined by the tips
    heights = {
        tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))
    }

    # Internal nodes: place at midpoint of children
    def calc_row(clade):
        for subclade in clade:
            if subclade not in heights:
                calc_row(subclade)
        # Closure over heights
        heights[clade] = (
            heights[clade.clades[0]] + heights[clade.clades[-1]]
        ) / 2.0

    if tree.root.clades:
        calc_row(tree.root)
        
    if adjust:
        pos = pd.DataFrame([(clade, val) for clade, val in heights.items()], columns=['clade','pos']).sort_values('pos')
        pos['newpos'] = 0
        count = 0
        for i in pos.index:
            if pos.loc[i,'clade'].name is not None and pos.loc[i,'clade'].name != 'root':
                count = count+1
            pos.loc[i, 'newpos'] = count

        pos.set_index('clade', inplace=True)
        heights = pos.to_dict()['newpos']

    return heights


def plot_tree(input_tree,
              label_func=None,
              title='',
              ax=None,
              output_name=None,
              normal_name='diploid',
              show_branch_lengths=False,
              branch_labels=None,
              show_support=False,
              label_colors=None,
              **kwargs):
    """Plot the given tree using matplotlib (or pylab).
    The graphic is a rooted tree, drawn with roughly the same algorithm as
    draw_ascii.
    Additional keyword arguments passed into this function are used as pyplot
    options. The input format should be in the form of:
    pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict), or
    pyplot_option_name=(dict).
    Example using the pyplot options 'axhspan' and 'axvline'::
        from Bio import Phylo, AlignIO
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        constructor = DistanceTreeConstructor()
        aln = AlignIO.read(open('TreeConstruction/msa.phy'), 'phylip')
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        tree = constructor.upgma(dm)
        Phylo.draw(tree, axhspan=((0.25, 7.75), {'facecolor':'0.5'}),
        ... axvline={'x':0, 'ymin':0, 'ymax':1})
    Visual aspects of the plot can also be modified using pyplot's own functions
    and objects (via pylab or matplotlib). In particular, the pyplot.rcParams
    object can be used to scale the font size (rcParams["font.size"]) and line
    width (rcParams["lines.linewidth"]).
    :Parameters:
        label_func : callable
            A function to extract a label from a node. By default this is str(),
            but you can use a different function to select another string
            associated with each node. If this function returns None for a node,
            no label will be shown for that node.
        do_show : bool
            Whether to show() the plot automatically.
        show_support : bool
            Whether to display confidence values, if present on the tree.
        ax : matplotlib/pylab axes
            If a valid matplotlib.axes.Axes instance, the phylogram is plotted
            in that Axes. By default (None), a new figure is created.
        branch_labels : dict or callable
            A mapping of each clade to the label that will be shown along the
            branch leading to it. By default this is the confidence value(s) of
            the clade, taken from the ``confidence`` attribute, and can be
            easily toggled off with this function's ``show_support`` option.
            But if you would like to alter the formatting of confidence values,
            or label the branches with something other than confidence, then use
            this option.
        label_colors : dict or callable
            A function or a dictionary specifying the color of the tip label.
            If the tip label can't be found in the dict or label_colors is
            None, the label will be shown in black.
    """

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        try:
            import pylab as plt
        except ImportError:
            raise MEDICCPlotError(
                "Install matplotlib or pylab if you want to use draw."
            ) from None

    import matplotlib.collections as mpcollections

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 7))

    label_func=label_func if label_func is not None else lambda x: x

    # options for displaying label colors.
    if label_colors is not None:
        if callable(label_colors):
            def get_label_color(label):
                return label_colors(label)
        else:
            # label_colors is presumed to be a dict
            def get_label_color(label):
                return label_colors.get(label, "black")
    else:
        clade_colors = {}
        for sample in [x.name for x in list(input_tree.find_clades(''))]:
            ## determine if sample is terminal
            is_terminal = True
            matches = list(input_tree.find_clades(sample))
            if len(matches) > 0:
                clade = matches[0]
                is_terminal = clade.is_terminal()
            ## determine if sample is normal
            is_normal = sample == normal_name
            clade_colors[sample] = COL_MARKER_TERMINAL
            if not is_terminal:
                clade_colors[sample] = COL_MARKER_INTERNAL
            if is_normal:
                clade_colors[sample] = COL_MARKER_NORMAL
        
        def get_label_color(label):
            return clade_colors.get(label, "black")

    marker_func=lambda x: (TREE_MARKER_SIZE, get_label_color(x.name)) if x.name is not None else None

    ax.axes.get_yaxis().set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True, prune=None))
    ax.xaxis.set_tick_params(labelsize=XLABEL_TICK_SIZE)
    ax.xaxis.label.set_size(XLABEL_FONT_SIZE)
    ax.set_title(title, x=0.01, y=1.0, ha='left', va='bottom',
                fontweight='bold', fontsize=16, zorder=10)
    x_posns = _get_x_positions(input_tree)
    y_posns = _get_y_positions(input_tree, adjust=True)

    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []

    # Options for displaying branch labels / confidence
    def conf2str(conf):
        if conf is None or conf == 0:
            return None
        elif int(conf) == conf:
            return str(int(conf))
        else:
            return str(conf)

    if not branch_labels:
        if show_branch_lengths:
            def format_branch_label(x): 
                return np.round(x.branch_length, 1) if x.name != 'root' and x.name is not None else None
        elif show_support:
            def format_branch_label(clade):
                try:
                    confidences = clade.confidences
                    # phyloXML supports multiple confidences
                except AttributeError:
                    pass
                else:
                    return "/".join(conf2str(cnf.value) for cnf in confidences)
                if clade.confidence is not None:
                    return conf2str(clade.confidence)
                return None

        else:
            def format_branch_label(clade):
                return None

    elif isinstance(branch_labels, dict):
        def format_branch_label(clade):
            return branch_labels.get(clade)
    else:
        if not callable(branch_labels):
            raise TypeError(
                "branch_labels must be either a dict or a callable (function)"
            )
        def format_branch_label(clade):
            return conf2str(branch_labels(clade))

    # The function draw_clade closes over the axes object
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
    elif not isinstance(ax, plt.matplotlib.axes.Axes):
        raise ValueError("Invalid argument for ax: %s" % ax)

    def draw_clade_lines(
        use_linecollection=False,
        orientation="horizontal",
        y_here=0,
        x_start=0,
        x_here=0,
        y_bot=0,
        y_top=0,
        color="black",
        lw=".1",
    ):
        """Create a line with or without a line collection object.
        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if not use_linecollection and orientation == "horizontal":
            ax.hlines(y_here, x_start, x_here, color=color, lw=lw)
        elif use_linecollection and orientation == "horizontal":
            horizontal_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw
                )
            )
        elif not use_linecollection and orientation == "vertical":
            ax.vlines(x_here, y_bot, y_top, color=color)
        elif use_linecollection and orientation == "vertical":
            vertical_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw
                )
            )

    def draw_clade(clade, x_start, color, lw):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        # phyloXML-only graphics annotations
        if hasattr(clade, "color") and clade.color is not None:
            color = clade.color.to_hex()
        if hasattr(clade, "width") and clade.width is not None:
            lw = clade.width * plt.rcParams["lines.linewidth"]
        # Draw a horizontal line from start to here
        draw_clade_lines(
            use_linecollection=True,
            orientation="horizontal",
            y_here=y_here,
            x_start=x_start,
            x_here=x_here,
            color=color,
            lw=lw,
        )
        # Add node marker
        if marker_func is not None:
            marker = marker_func(clade)
            if marker is not None:
                marker_size, marker_col = marker_func(clade)
                ax.scatter(x_here, y_here, s=marker_size, c=marker_col, zorder=3)
        # Add node/taxon labels
        label = label_func(str(clade))
        if label not in (None, clade.__class__.__name__):
            ax.text(
                x_here + 1,
                y_here,
                " %s" % label,
                verticalalignment="center",
                color=get_label_color(label),
            )
        # Add label above the branch (optional)
        conf_label = format_branch_label(clade)
        if conf_label:
            ax.text(
                0.5 * (x_start + x_here),
                y_here - 0.15,
                conf_label,
                fontsize="small",
                horizontalalignment="center",
            )
        if clade.clades:
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.clades[0]]
            y_bot = y_posns[clade.clades[-1]]
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(
                use_linecollection=True,
                orientation="vertical",
                x_here=x_here,
                y_bot=y_bot,
                y_top=y_top,
                color=color,
                lw=lw,
            )
            # Draw descendents
            for child in clade:
                draw_clade(child, x_here, color, lw)

    draw_clade(input_tree.root, 0, "k", plt.rcParams["lines.linewidth"])

    # If line collections were used to create clade lines, here they are added
    # to the pyplot plot.
    for i in horizontal_linecollections:
        ax.add_collection(i)
    for i in vertical_linecollections:
        ax.add_collection(i)

    ax.set_xlabel("branch length")
    ax.set_ylabel("taxa")

    # Add margins around the `tree` to prevent overlapping the ax
    xmax = max(x_posns.values())
    #ax.set_xlim(-0.05 * xmax, 1.25 * xmax)
    ax.set_xlim(-0.05 * xmax, 1.05 * xmax)
    # Also invert the y-axis (origin at the top)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    ax.set_ylim(max(y_posns.values()) + 0.5, 0.5)

    # Parse and process key word arguments as pyplot options
    for key, value in kwargs.items():
        try:
            # Check that the pyplot option input is iterable, as required
            list(value)
        except TypeError:
            raise ValueError(
                'Keyword argument "%s=%s" is not in the format '
                "pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),"
                " or pyplot_option_name=(dict) " % (key, value)
            ) from None
        if isinstance(value, dict):
            getattr(plt, str(key))(**dict(value))
        elif not (isinstance(value[0], tuple)):
            getattr(plt, str(key))(*value)
        elif isinstance(value[0], tuple):
            getattr(plt, str(key))(*value[0], **dict(value[1]))

    if output_name is not None:
        plt.savefig(output_name + ".png")


class MEDICCPlotError(Exception):
    pass
