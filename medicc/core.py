import copy
import logging
import os
from itertools import combinations

import Bio
import fstlib
import numpy as np
import pandas as pd
import pyranges as pr
from joblib import Parallel, delayed

import medicc
from medicc import io, nj, stats, tools

# prepare logger 
logger = logging.getLogger(__name__)


def main(input_df,
         asymm_fst,
         normal_name='diploid',
         input_tree=None,
         ancestral_reconstruction=True,
         chr_separator='X',
         prune_weight=0,
         allele_columns=['cn_a', 'cn_b'],
         n_cores=None):
    """ MEDICC Main Method """

    symbol_table = asymm_fst.input_symbols()

    ## Validate input
    logger.info("Validating input.")
    io.validate_input(input_df, symbol_table)

    ## Compile input data into FSAs stored in dictionaries
    logger.info("Compiling input sequences into FSAs.")
    FSA_dict = create_standard_fsa_dict_from_data(input_df, symbol_table, chr_separator)
    sample_labels = input_df.index.get_level_values('sample_id').unique()

    ## Calculate pairwise distances
    logger.info("Calculating pairwise distance matrices")
    if n_cores is not None and n_cores > 1:
        pairwise_distances = parallelization_calc_pairwise_distance_matrix(sample_labels, 
                                                                       asymm_fst,
                                                                       FSA_dict,
                                                                       n_cores)
    else:
        pairwise_distances = calc_pairwise_distance_matrix(asymm_fst, FSA_dict)

    ## Reconstruct a tree
    if input_tree is None:
        logger.info("Inferring tree topology.")
        nj_tree = infer_tree_topology(
            pairwise_distances.values, pairwise_distances.index, normal_name=normal_name)
    else:
        logger.info("Tree provided, using it.")

        assert len([x for x in list(input_tree.find_clades()) if x.name is not None and 'internal' not in x.name]) == \
            len(np.unique(input_df.index.get_level_values('sample_id'))), \
            "Number of samples differs in input tree and input dataframe"
        assert np.all(
            np.sort([x.name for x in list(input_tree.find_clades()) if x.name is not None and 'internal' not in x.name]) ==
            np.sort(np.unique(input_df.index.get_level_values('sample_id')))), \
            "Input tree does not match input dataframe"
        
        # necessary for the way that reconstruct_ancestors is performed
        if ancestral_reconstruction:
            input_tree.root_with_outgroup([x for x in input_tree.root.clades if x.name != normal_name][0].name)

        nj_tree = input_tree


    tools.set_sequences_on_tree_from_df(nj_tree, input_df)

    final_tree = copy.deepcopy(nj_tree)

    if ancestral_reconstruction:
        logger.info("Reconstructing ancestors.")
        ancestors = medicc.reconstruct_ancestors(tree=final_tree,
                                                 samples_dict=FSA_dict,
                                                 fst=asymm_fst,
                                                 normal_name=normal_name,
                                                 prune_weight=prune_weight)

        ## Create and write output data frame with ancestors
        logger.info("Creating output copynumbers.")
        output_df = create_df_from_fsa(input_df, ancestors)

        ## Update branch lengths with ancestors
        logger.info("Updating branch lengths of final tree using ancestors.")
        tools.set_sequences_on_tree_from_df(final_tree, output_df)
        update_branch_lengths(final_tree, asymm_fst, ancestors, normal_name)

    nj_tree.root_with_outgroup(normal_name)
    final_tree.root_with_outgroup(normal_name)

    if ancestral_reconstruction:
        output_df, events_df = calculate_all_cn_events(
            final_tree, output_df, allele_columns, normal_name)
        if len(events_df) != final_tree.total_branch_length():
            logger.warn("Event recreation was faulty. Events in '_cn_events_df.tsv' might be incorrect")
    else:
        events_df = None
        output_df = input_df


    return sample_labels, pairwise_distances, nj_tree, final_tree, output_df, events_df


def create_standard_fsa_dict_from_data(input_data,
                                       symbol_table: fstlib.SymbolTable,
                                       separator: str = "X") -> dict:
    """ Creates a dictionary of FSAs from input DataFrame or Series.
    The keys of the dictionary are the sample/taxon names. 
    If the input is a DataFrame, the FSA will be the concatenated copy number profiles of all allele columns"""

    fsa_dict = {}
    if isinstance(input_data, pd.DataFrame):
        logger.info('Creating FSA for pd.DataFrame with the following data columns:\n{}'.format(
            input_data.columns))
        def aggregate_copy_number_profile(copy_number_profile):
            return separator.join([separator.join(["".join(x.astype('str'))
                                                   for _, x in cnp[allele].groupby('chrom')]) for allele in copy_number_profile.columns])

    elif isinstance(input_data, pd.Series):
        logger.info('Creating FSA for pd.Series with the name {}'.format(input_data.name))
        def aggregate_copy_number_profile(copy_number_profile):
            return separator.join(["".join(x.astype('str')) for _, x in copy_number_profile.groupby('chrom')])

    else:
        raise MEDICCError("Input to function create_standard_fsa_dict_from_data has to be either"
                          "pd.DataFrame or pd.Series. \n input provided was {}".format(type(input_data)))
    
    for taxon, cnp in input_data.groupby('sample_id'):
        cn_str = aggregate_copy_number_profile(cnp)
        fsa_dict[taxon] = fstlib.factory.from_string(cn_str,
                                                     arc_type="standard",
                                                     isymbols=symbol_table,
                                                     osymbols=symbol_table)

    return fsa_dict


def create_phasing_fsa_dict_from_df(input_df: pd.DataFrame, symbol_table: fstlib.SymbolTable, separator: str = "X") -> dict:
    """ Creates a dictionary of FSAs from two allele columns (Pandas DataFrame).
    The keys of the dictionary are the sample/taxon names. """
    allele_columns = input_df.columns
    if len(allele_columns) != 2:
        raise MEDICCError("Need exactly two alleles for phasing.")

    fsa_dict = {}
    for taxon, cnp in input_df.groupby('sample_id'):
        allele_a = cnp[allele_columns[0]]
        allele_b = cnp[allele_columns[1]]
        cn_str_a = separator.join(["".join(x) for _,x in allele_a.groupby(level='chrom', sort=False)])
        cn_str_b = separator.join(["".join(x) for _,x in allele_b.groupby(level='chrom', sort=False)])
        encoded = np.array([list(zip(cn_str_a, cn_str_b)), list(zip(cn_str_b, cn_str_a))])
        fsa_dict[taxon] = fstlib.factory.from_array(encoded, symbols=symbol_table, arc_type='standard')
        fsa_dict[taxon] = fstlib.determinize(fsa_dict[taxon]).minimize()

    return fsa_dict

def phase(input_df: pd.DataFrame, model_fst: fstlib.Fst, reference_sample='diploid', separator: str = 'X') -> pd.DataFrame:
    """ Phases every FST against the reference sample. 
    Returns two standard FSA dicts, one for each allele. """
    
    diploid_fsa = medicc.tools.create_diploid_fsa(model_fst)
    phasing_dict = medicc.create_phasing_fsa_dict_from_df(input_df, model_fst.input_symbols(), separator)
    fsa_dict_a, fsa_dict_b, _ = phase_dict(phasing_dict, model_fst, diploid_fsa)
    output_df = medicc.create_df_from_phasing_fsa(input_df, [fsa_dict_a, fsa_dict_b], separator)

    return output_df

def phase_dict(phasing_dict, model_fst, reference_fst):
    """ Phases every FST against the reference sample. 
    Returns two standard FSA dicts, one for each allele. """
    fsa_dict_a = {}    
    fsa_dict_b = {}
    scores = {}
    left = (reference_fst * model_fst).project('output')
    right = (~model_fst * reference_fst).project('input')
    for sample_id, sample_fst in phasing_dict.items():
        phased_fst = fstlib.align(sample_fst, left, right).topsort()
        score = fstlib.shortestdistance(phased_fst, reverse=True)[phased_fst.start()]
        scores[sample_id] = float(score)
        fsa_dict_a[sample_id] = fstlib.arcmap(phased_fst.copy().project('input'), map_type='rmweight')
        fsa_dict_b[sample_id] = fstlib.arcmap(phased_fst.project('output'), map_type='rmweight')
    
    return fsa_dict_a, fsa_dict_b, scores


def create_df_from_fsa(input_df: pd.DataFrame, fsa, separator: str = 'X'):
    """ 
    Takes a single FSA dict or a list of FSA dicts and extracts the copy number profiles.
    The allele names are taken from the input_df columns and the returned data frame has the same 
    number of rows and row index as the input_df. """

    alleles = input_df.columns
    if not isinstance(fsa, dict):
        raise MEDICCError("fsa input to create_df_from_fsa has to be a dict"
                          "Input type is {}".format(type(fsa)))

    nr_alleles = len(alleles)
    samples = input_df.index.get_level_values('sample_id').unique()
    output_df = input_df.unstack('sample_id')

    # Create dict and concat later to prevent pandas PerformanceWarning
    internal_cns = dict()
    for node in fsa:
        if node in samples:
            continue
        cns = tools.fsa_to_string(fsa[node]).split(separator)
        if len(cns) % nr_alleles != 0:
            raise MEDICCError('For sample {} we have {} haplotype-specific chromosomes for {} alleles'
                              '\nnumber of chromosomes has to be divisible by nr of alleles'.format(node,
                                                                                                    len(cns),
                                                                                                    nr_alleles))
        nr_chroms = int(len(cns) // nr_alleles)
        for i, allele in enumerate(alleles):
            cn = list(''.join(cns[(i*nr_chroms):((i+1)*nr_chroms)]))
            internal_cns[(allele, node)] = cn

    output_df = (pd.concat([output_df, pd.DataFrame(internal_cns, index=output_df.index)], axis=1)
                 .stack('sample_id')
                 .reorder_levels(['sample_id', 'chrom', 'start', 'end'])
                 .sort_index())

    return output_df


def create_df_from_phasing_fsa(input_df: pd.DataFrame, fsas, separator: str = 'X'):
    """ 
    Takes a two FSAs dicts from phasing and extracts the copy number profiles.
    The allele names are taken from the input_df columns and the returned data frame has the same 
    number of rows and row index as the input_df. """

    alleles = input_df.columns
    if len(fsas) != 2:
        raise MEDICCError("fsas has to be of length 2")
    if not all([isinstance(fsa, dict) for fsa in fsas]):
        raise MEDICCError("all fsas entries have to be dicts")
    if fsas[0].keys() != fsas[1].keys():
        raise MEDICCError("fsas keys have to be the same")


    output_df = input_df.copy()[[]]
    output_df[alleles] = ''

    for sample in fsas[0].keys():
        cns_a = tools.fsa_to_string(fsas[0][sample]).split(separator)
        cns_b = tools.fsa_to_string(fsas[1][sample]).split(separator)
        if len(cns_a) != len(cns_b):
            raise MEDICCError(f"length of alleles is not the same for sample {sample}")

        output_df.loc[sample, alleles[0]] = list(''.join(cns_a))
        output_df.loc[sample, alleles[1]] = list(''.join(cns_b))

    # output_df = output_df.stack('sample_id')
    # output_df = output_df.reorder_levels(['sample_id', 'chrom', 'start', 'end']).sort_index()
    
    return output_df


def parallelization_calc_pairwise_distance_matrix(sample_labels, asymm_fst, FSA_dict, n_cores):
    parallelization_groups = medicc.tools.create_parallelization_groups(len(sample_labels))
    parallelization_groups = [sample_labels[group] for group in parallelization_groups]
    logger.info("Running {} parallel runs on {} cores".format(len(parallelization_groups), n_cores))

    parallel_pairwise_distances = Parallel(n_jobs=n_cores)(delayed(calc_pairwise_distance_matrix)(
        asymm_fst, {key: val for key, val in FSA_dict.items() if key in cur_group}, True)
            for cur_group in parallelization_groups)

    pdm = medicc.tools.total_pdm_from_parallel_pdms(sample_labels, parallel_pairwise_distances)

    return pdm


def calc_pairwise_distance_matrix(model_fst, fsa_dict, parallel_run=True):
    '''Given a symmetric model FST and input FSAs in a form of a dictionary, output pairwise distance matrix'''

    samples = list(fsa_dict.keys())
    pdm = pd.DataFrame(0, index=samples, columns=samples, dtype=float)
    combs = list(combinations(samples, 2))
    ncombs = len(combs)

    for i, (sample_a, sample_b) in enumerate(combs):
        cur_dist = float(fstlib.kernel_score(model_fst, fsa_dict[sample_a], fsa_dict[sample_b]))
        pdm[sample_a][sample_b] = cur_dist
        pdm[sample_b][sample_a] = cur_dist

        if not parallel_run and (100*(i+1)/ncombs) % 10 == 0:  # log every 10%
            logger.info('%.2f%%', (i+1)/ncombs * 100)

    return pdm


def infer_tree_topology(pairwise_distances, labels, normal_name):
    if len(labels) > 2:
        tree = nj.NeighbourJoining(pairwise_distances, labels).tree

        tmpsearch = [c for c in tree.find_clades(name = normal_name)]
        normal_node = tmpsearch[0]
        root_path = tree.get_path(normal_node)[::-1]

        if len(root_path)>1:
            new_root = root_path[1]
            tree.root_with_outgroup(new_root)
    else:
        clade_ancestor = Bio.Phylo.PhyloXML.Clade(branch_length=0, name='internal_1')
        clade_ancestor.clades = [Bio.Phylo.PhyloXML.Clade(
            name=label, branch_length=0 if label == normal_name else 1) for label in labels]

        tree = Bio.Phylo.PhyloXML.Phylogeny(root=clade_ancestor)
        tree.root_with_outgroup(normal_name)

    return tree


def update_branch_lengths(tree, fst, ancestor_fsa, normal_name='diploid'):
    """ Updates the branch lengths in the tree using the internal nodes supplied in the FSA dict 
    """
    if len(ancestor_fsa) == 2:
        child_clade = [x for x in tree.find_clades() if x.name is not None and x.name != normal_name][0]
        child_clade.branch_length = float(fstlib.score(
            fst, ancestor_fsa[normal_name], ancestor_fsa[child_clade.name]))

    if not isinstance(ancestor_fsa, dict):
        raise MEDICCError("input ancestor_fsa to function update_branch_lengths has to be either a dict"
                          "provided type is {}".format(type(ancestor_fsa)))

    def _distance_to_child(fst, ancestor_fsa, sample_1, sample_2):
        return float(fstlib.score(fst, ancestor_fsa[sample_1], ancestor_fsa[sample_2]))

    for clade in tree.find_clades():
        if clade.name is None:
            continue
        children = clade.clades
        if len(children) != 0:
            for child in children:
                if child.name == normal_name:  # exception: evolution goes from diploid to internal node
                    brs = _distance_to_child(fst, ancestor_fsa, child.name, clade.name)
                else:
                    brs = _distance_to_child(fst, ancestor_fsa, clade.name, child.name)
                child.branch_length = brs


def calculate_all_cn_events(tree, cur_df, alleles=['cn_a', 'cn_b'], normal_name='diploid'):
    """Create a DataFrame containing all copy-number events in the current data

    Args:
        tree (Bio.Phylo.Tree): Phylogenetic tree created by MEDICC2's tree reconstruction
        cur_df (pandas.DataFrame): DataFrame containing the copy-numbers of the samples and internal nodes
        alleles (list, optional): List of alleles. Defaults to ['cn_a', 'cn_b'].
        normal_name (str, optional): Name of the normal sample. Defaults to 'diploid'.

    Returns:
        pandas.DataFrame: Updated copy-number DataFrame
        pandas.DataFrame: DataFrame of copy-number events
    """
    
    cur_df[['is_gain', 'is_loss', 'is_wgd']] = False
    cur_df[alleles] = cur_df[alleles].astype(int)
    if tree == None:
        cur_df[['is_normal', 'is_clonal']] = False
        events = None
    else:

        events = pd.DataFrame(columns=['sample_id', 'chrom', 'start',
                                    'end', 'allele', 'type', 'cn_child'])

        clades = [x for x in tree.find_clades()]

        for clade in clades:
            if not len(clade.clades):
                continue
            if clade.name is None:
                clade = copy.deepcopy(clade)
                clade.name = 'diploid'
            for child in clade.clades:
                if child.branch_length == 0:
                    continue

                cur_df, cur_events = calculate_cn_events_per_branch(
                    cur_df, clade.name, child.name, alleles=alleles)

                events = pd.concat([events, cur_events])

        events = events.reset_index(drop=True)

        is_normal = ~cur_df.unstack('sample_id')[['is_loss', 'is_gain', 'is_wgd']].any(axis=1)
        is_normal.name = 'is_normal'
        mrca = [x for x in tree.root.clades if x.name != normal_name][0].name
        is_clonal = ~cur_df.loc[cur_df.index.get_level_values('sample_id')!=mrca].unstack('sample_id')[['is_loss', 'is_gain', 'is_wgd']].any(axis=1)
        is_clonal.name = 'is_clonal'

        cur_df = cur_df.drop(['is_normal', 'is_clonal'], axis=1, errors='ignore')
        cur_df = (cur_df
                .join(is_normal, how='inner')
                .reorder_levels(['sample_id', 'chrom', 'start', 'end'])
                .sort_index()
                .join(is_clonal, how='inner')
                .reset_index())
        cur_df['chrom'] = tools.format_chromosomes(cur_df['chrom'])
        cur_df = (cur_df
                .set_index(['sample_id', 'chrom', 'start', 'end'])
                .sort_index())

        events = events.set_index(['sample_id', 'chrom', 'start', 'end'])

    return cur_df, events


def calculate_cn_events_per_branch(cur_df, parent_name, child_name, alleles=('cn_a', 'cn_b')):
    """Calculate copy-number events for a single branch. Used in calculate_all_cn_events

    Args:
        cur_df (pandas.DataFrame): DataFrame containing the copy-numbers of the samples and internal nodes
        parent_name (str): Name of the parent sample
        child_name (str): Name of the child sample
        alleles (list, optional): List of alleles. Defaults to ['cn_a', 'cn_b'].

    Returns:
        pandas.DataFrame: Updated copy-number DataFrame
        pandas.DataFrame: DataFrame of copy-number events
    """

    asymm_fst, asymm_fst_nowgd, asymm_fst_1_wgd, symbol_table = io.load_main_fsts(
        return_symbol_table=True)

    events_df = pd.DataFrame(columns=['sample_id', 'chrom', 'start', 'end', 'allele', 'type', 'cn_child'])

    cur_parent_cn = cur_df.loc[parent_name, alleles].astype(int)
    cur_child_cn = cur_df.loc[child_name, alleles].astype(int)
    cur_chroms = cur_df.loc['diploid'].index.get_level_values(
        'chrom').map(lambda x: int(x.split('chr')[-1])).values.astype(int)

    # 1. find total losss (loh)
    parent_loh = cur_parent_cn == 0
    for allele in alleles:

        cur_loh = cur_child_cn.loc[~parent_loh[allele], allele] == 0
        if cur_loh.sum() == 0:
            continue

        cur_df.loc[child_name, 'is_loss'] = (cur_df.loc[child_name, 'is_loss'].values 
                                             + np.logical_and(cur_child_cn[allele] == 0, 
                                                              cur_parent_cn[allele] != 0).values)

        max_previous_cn = np.max(
            np.unique(cur_parent_cn.loc[~parent_loh[allele], allele].loc[cur_loh]))

        for _ in np.arange(max_previous_cn):
            cur_loh_and_parental_val = np.logical_and(cur_loh.values, 
                                                        cur_parent_cn.loc[~parent_loh[allele], allele] > 0).values

            event_labels_ = ((np.cumsum(np.concatenate([[0], np.diff(
                (cur_loh_and_parental_val + cur_chroms[~parent_loh[allele]]))])
                * cur_loh_and_parental_val) + 1)
                * cur_loh_and_parental_val)

            event_labels = np.zeros_like(event_labels_)
            for i, j in enumerate(np.unique(event_labels_)):
                event_labels[event_labels_ == j] = i

            cur_parent_cn.loc[parent_loh.loc[~parent_loh[allele],
                                               allele].index[cur_loh_and_parental_val], allele] -= 1

            cur_events = (cur_parent_cn
                        .loc[~parent_loh[allele]]
                        .reset_index()
                        .loc[np.array([np.argmax(event_labels == ind) for ind in np.setdiff1d(np.unique(event_labels), [0])])]
                        [['chrom', 'start', 'end']].values)
            # adjust ends
            cur_events[:, 2] = (cur_parent_cn
                                .loc[~parent_loh[allele]]
                                .reset_index()
                                .loc[np.array([len(event_labels) - np.argmax(event_labels[::-1] == ind) - 1 for ind in np.setdiff1d(np.unique(event_labels), [0])])]
                                ['end'].values)

            cur_ind = np.arange(len(events_df), len(events_df)+len(cur_events))
            events_df = pd.concat([events_df, pd.DataFrame(index=cur_ind)])
            events_df.loc[cur_ind, 'sample_id'] = child_name
            events_df.loc[cur_ind, 'allele'] = allele
            events_df.loc[cur_ind, 'type'] = 'loh'
            events_df.loc[cur_ind, 'cn_child'] = 0
            events_df.loc[cur_ind, ['chrom', 'start', 'end']] = cur_events[:, :3]

            # recalculate parental_loss and cur_loh for next iteration
            parent_loh = cur_parent_cn <= 0
            cur_loh = cur_child_cn.loc[~parent_loh[allele], allele] == 0

        cur_parent_cn.loc[cur_parent_cn[allele] < 0, allele] = 0

    # 2. WGDs
    # only check if >30% of is gained
    wgd_candidate_threshold = 0.3

    widths = cur_df.loc[[child_name]].eval('end-start')
    fraction_gain = ((cur_df.loc[child_name, alleles] > 1).astype(int).sum(axis=1) * widths.loc[child_name]
                        ).sum() / (2 * widths.loc[child_name].sum())
    if fraction_gain > wgd_candidate_threshold:

        parent_fsa = fstlib.factory.from_string('X'.join(['X'.join(["".join(x.astype('str')) for _, x in cur_df.loc[parent_name, alleles][allele].groupby('chrom')]) for allele in alleles]),
                                                arc_type="standard",
                                                isymbols=symbol_table,
                                                osymbols=symbol_table)
        child_fsa = fstlib.factory.from_string('X'.join(['X'.join(["".join(x.astype('str')) for _, x in cur_df.loc[child_name, alleles][allele].groupby('chrom')]) for allele in alleles]),
                                                arc_type="standard",
                                                isymbols=symbol_table,
                                                osymbols=symbol_table)

        score_wgd = float(fstlib.score(asymm_fst, parent_fsa, child_fsa))
        fraction_double_gain = (((cur_df.loc[child_name, alleles] > 2)
                                 .astype(int)
                                 .sum(axis=1)
                                 * widths.loc[child_name]
                                 ).sum() / widths.loc[child_name].sum())

        # double wgd
        if fraction_double_gain and (float(fstlib.score(asymm_fst_1_wgd, parent_fsa, child_fsa)) != score_wgd):
            cur_parent_cn = 4 * cur_parent_cn
            events_df.loc[len(events_df.index)] = [child_name, 'chr0', cur_df.index.get_level_values('start').min(),
                                             cur_df.index.get_level_values('end').max(), 'both', 'wgd', 0]
            events_df.loc[len(events_df.index)] = [child_name, 'chr0', cur_df.index.get_level_values('start').min(),
                                                cur_df.index.get_level_values('end').max(), 'both', 'wgd', 0]
            cur_df.loc[child_name, 'is_wgd'] = True
        # single wgd
        elif float(fstlib.score(asymm_fst_nowgd, parent_fsa, child_fsa)) != score_wgd:
            cur_parent_cn = 2 * cur_parent_cn
            events_df.loc[len(events_df.index)] = [child_name, 'chr0', cur_df.index.get_level_values('start').min(),
                                             cur_df.index.get_level_values('end').max(), 'both', 'wgd', 0]
            cur_df.loc[child_name, 'is_wgd'] = True

    # 3. losses and gains
    loh_pos = (cur_parent_cn == 0)
    for allele in alleles:

        cn_changes = (cur_child_cn[allele] - cur_parent_cn[allele]).values
        all_cn_change_vals = np.unique(cn_changes)

        cur_df.loc[child_name, 'is_loss'] = np.logical_or(cur_df.loc[child_name, 'is_loss'].values,
                                                          cn_changes < 0)
        cur_df.loc[child_name, 'is_gain'] = np.logical_or(cur_df.loc[child_name, 'is_gain'].values,
                                                          cn_changes > 0)

        all_cn_change_vals = np.setdiff1d(np.arange(np.min(all_cn_change_vals), np.max(all_cn_change_vals)+1), [0])
        for cur_cn_change in all_cn_change_vals[np.argsort(np.abs(all_cn_change_vals))[::-1]]:
            cur_event = 'gain' if cur_cn_change > 0 else 'loss'

            cur_change_location = ((cur_child_cn.loc[~loh_pos[allele], allele] -
                        cur_parent_cn.loc[~loh_pos[allele], allele]) == cur_cn_change)

            event_labels_ = ((np.cumsum(np.concatenate([[0], np.diff(
                (cur_change_location.values + cur_chroms[~loh_pos[allele]]))])
                * cur_change_location.values) + 1)
                * cur_change_location.values)

            event_labels = np.zeros_like(event_labels_)
            for i, j in enumerate(np.unique(event_labels_)):
                event_labels[event_labels_ == j] = i

            cur_events = (cur_child_cn
                          .loc[~loh_pos[allele]]
                          .reset_index()
                          .loc[np.array([np.argmax(event_labels == val) for val in np.setdiff1d(np.unique(event_labels), [0])])]
                          [['chrom', 'start', 'end', allele]].values)
            # adjust ends
            cur_events[:, 2] = (cur_child_cn
                                .loc[~loh_pos[allele]]
                                .reset_index()
                                .loc[np.array([len(event_labels) - np.argmax(event_labels[::-1] == val) - 1 for val in np.setdiff1d(np.unique(event_labels), [0])])]
                                ['end'].values)

            cur_ind = np.arange(len(events_df), len(events_df)+len(cur_events))
            events_df = pd.concat([events_df, pd.DataFrame(index=cur_ind)])
            events_df.loc[cur_ind, 'sample_id'] = child_name
            events_df.loc[cur_ind, 'allele'] = allele
            events_df.loc[cur_ind, 'type'] = cur_event
            events_df.loc[cur_ind, 'cn_child'] = cur_events[:, 3]
            events_df.loc[cur_ind, ['chrom', 'start', 'end']] = cur_events[:, :3]

            cur_child_cn.loc[np.intersect1d(loh_pos.loc[~loh_pos[allele]].index,
                                            cur_change_location.loc[cur_change_location].index), allele] += (1 if (cur_cn_change < 0) else -1)

    events_df['chrom'] = tools.format_chromosomes(events_df['chrom'])
    events_df = (events_df[['sample_id', 'allele', 'chrom', 'start', 'end', 'type', 'cn_child']]
                 .reset_index(drop=True)
                 .sort_values(['sample_id', 'allele', 'chrom', 'start', 'end', 'type', 'cn_child']))

    return cur_df, events_df


def compute_cn_change(df, tree, normal_name='diploid'):
    """Compute the copy-number changes per segment in all branches

    Args:
        df (pandas.DataFrame): DataFrame containing the copy-numbers of samples and internal nodes
        tree (Bio.Phylo.Tree): Phylogenetic tree
        normal_name (str, optional): Name of normal sample. Defaults to 'diploid'.

    Returns:
        pandas.DataFrame: DataFrame containing the copy-number changes
    """    
    cn_change = df.copy()
    alleles = cn_change.columns
    for allele in alleles:
        cn_change[allele] = cn_change[allele].astype('int')

    clades = [clade for clade in tree.find_clades(order = "postorder") if clade.name is not None and clade.name != normal_name]
    for clade in clades:
        for child in clade.clades:
            cn_change.loc[child.name, alleles] = cn_change.loc[child.name, alleles].values - cn_change.loc[clade.name, alleles].values
    cn_change.loc[clades[-1].name, alleles] = cn_change.loc[clades[-1].name, alleles].values - cn_change.loc[normal_name, alleles].values
    cn_change.loc[normal_name, alleles] = 0

    return cn_change


def summarize_patient(tree, pdm, sample_labels, normal_name='diploid', events_df=None):
    """Calculate several summary values for the provided samples

    Args:
        tree (Bio.Phylo.Tree): Phylogenetic tree
        pdm (pandas.DataFrame): Pairwise distance matrix between the samples
        sample_labels (list): List of all samples
        normal_name (str, optional): Name of normal sample. Defaults to 'diploid'.
        events_df (pandas.DataFrame, optional): DataFrame containg all copy-number events. Defaults to None.

    Returns:
        pandas.DataFrame: Summary DataFrame
    """    
    branch_lengths = []
    for parent in tree.find_clades(terminal=False, order="level"):
        for child in parent.clades:
            if child.branch_length:
                branch_lengths.append(child.branch_length)
    
    nsamples = len(sample_labels)
    tree_length = np.sum(branch_lengths)
    avg_branch_length = np.mean(branch_lengths)
    min_branch_length = np.min(branch_lengths)
    max_branch_length = np.max(branch_lengths)
    median_branch_length = np.median(branch_lengths)
    # p_star = stats.star_topology_test(pdm)
    # p_clock = stats.molecular_clock_test(pdm,
    #                                      np.flatnonzero(np.array(sample_labels) == normal_name)[0])
    if events_df is None:
        wgd_status = "unknown"
    else:
        if "wgd" in events_df['type'].values:
            wgd_status = "WGD on branch " + \
                "and ".join(events_df.loc[events_df['type'] ==
                                          'wgd'].index.get_level_values('sample_id'))
        else:
            wgd_status = "no WGD"

    result = pd.Series({
        'nsamples': nsamples,
        'normal_name': normal_name,
        'tree_length': tree_length,
        'mean_branch_length': avg_branch_length,
        'median_branch_length': median_branch_length,
        'min_branch_length': min_branch_length,
        'max_branch_length': max_branch_length,
        # 'p_star': p_star,
        # 'p_clock': p_clock,
        'wgd_status': wgd_status,
    })
    
    return result


def overlap_events(events_df=None, output_df=None, tree=None, overlap_threshold=0.9,
                   chromosome_bed='default', regions_bed='default',
                   replace_loh_with_loss=True, alleles=['cn_a', 'cn_b'],
                   replace_both_arms_with_chrom=True, normal_name='diploid'):
    """Overlap copy-number events with regions of interest

    Args:
        events_df (pandas.DataFrame, optional): All copy-number events. Defaults to None.
        output_df (pandas.DataFrame, optional): DataFrame containing all copy-numbers. Defaults to None.
        tree (Bio.Phylo.Tree, optional): Phylogenetic tree. Defaults to None.
        overlap_threshold (float, optional): Threshold above which an overlap is considered. Defaults to 0.9.
        chromosome_bed (str, optional): Name of BED file containing chromosome arm information. Defaults to 'default'.
        regions_bed (str, optional): Name of BED file containing regions of interest. Defaults to 'default'.
        replace_loh_with_loss (bool, optional): If True, loh is considered like a normal loss. Defaults to True.
        alleles (list, optional): List of alleles. Defaults to ['cn_a', 'cn_b'].
        replace_both_arms_with_chrom (bool, optional): If True, an event in the p- and q-arm of a chromosome will be displayed as a single event. Defaults to True.
        normal_name (str, optional): Name of normal sample. Defaults to 'diploid'.

    Returns:
        pandas.DataFrame: DataFrame with events concerning the regions of interest
    """                   

    if chromosome_bed == 'default':
        chromosome_bed = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                      "objects", "hg19_chromosome_arms.bed")
    if regions_bed == 'default':
        regions_bed = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   "objects", "Davoli_2013_TSG_OG_genes.bed")

    all_events = pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'name', 'NumberOverlaps',
                                       'FractionOverlaps', 'event', 'branch']).set_index(['Chromosome', 'Start', 'End'])

    if events_df is None:
        if output_df is None or tree is None:
            raise MEDICCError("Either events_df or df and tree has to be specified")
        _, events_df = calculate_all_cn_events(tree, output_df, alleles=alleles, normal_name=normal_name)
    if replace_loh_with_loss:
        events_df.loc[events_df['type'] == 'loh', 'type'] = 'loss'

    # Read chromosome regions and other regions
    if chromosome_bed is None and regions_bed is None:
        raise MEDICCError("Either chromosome_bed or regions_bed has to be specified")

    chr_arm_regions = None
    if chromosome_bed is not None:
        chr_arm_regions = medicc.io.read_bed_file(chromosome_bed)
        whole_chromosome = chr_arm_regions.groupby('Chromosome').min()
        whole_chromosome['End'] = chr_arm_regions.groupby('Chromosome')['End'].max()
        whole_chromosome['name'] = whole_chromosome.index
        chr_arm_regions = pd.concat([chr_arm_regions, whole_chromosome.reset_index()]).sort_values('Chromosome')
        chr_arm_regions = pr.PyRanges(chr_arm_regions)

    regions = None
    if regions_bed is not None:
        regions = []
        if isinstance(regions_bed, list) or isinstance(regions_bed, tuple):
            for f in regions_bed:
                regions.append(pr.PyRanges(medicc.io.read_bed_file(f)))
        else:
            regions.append(pr.PyRanges(medicc.io.read_bed_file(regions_bed)))

    # Add WGD
    for ind, _ in events_df.loc[events_df['type'] == 'wgd'].iterrows():
        all_events.loc[('all', '0', '0')] = ['WGD', 1, 1., 'WGD', ind[0]]

    for cur_branch in events_df.index.get_level_values('sample_id').unique():
        for allele in alleles:
            cur_events_df = events_df.loc[cur_branch]
            cur_events_df = cur_events_df.loc[cur_events_df['allele']==allele]
            for event_type in ['gain', 'loss'] if replace_loh_with_loss else ['gain', 'loh', 'loss']:
                cur_events_ranges = pr.PyRanges(cur_events_df.loc[cur_events_df['type'] == event_type].reset_index(
                ).rename({'chrom': 'Chromosome', 'start': 'Start', 'end': 'End'}, axis=1))

                # Calculate chromosomal events
                if chr_arm_regions is not None:
                    chr_events = overlap_regions(
                        chr_arm_regions, cur_events_ranges, event_type, cur_branch, overlap_threshold)
                    # remove arms if the whole chromosome is in there
                    if replace_both_arms_with_chrom and len(chr_events) > 0:
                        chr_events = chr_events[~chr_events['name'].isin(np.concatenate(
                            [[name + 'p', name + 'q'] if ('q' not in name and 'p' not in name) else [] for name in chr_events['name']]))]
                    all_events = pd.concat([all_events, chr_events])

                # Calculate other events
                if regions is not None:
                    for region in regions:
                        chr_events = overlap_regions(
                            region, cur_events_ranges, event_type, cur_branch, overlap_threshold)
                        all_events = pd.concat([all_events, chr_events])

    all_events['final_name'] = all_events['name'].apply(lambda x: x.split(
        'chr')[-1]) + all_events['event'].apply(lambda x: ' +' if x == 'gain' else (' -' if x == 'loss' else (' 0' if x == 'loh' else '')))

    all_events = all_events.set_index(['branch', 'name']).drop('NumberOverlaps', axis=1)
    all_events = all_events.reset_index().set_index('branch')
    
    return all_events


def overlap_regions(region, cur_events_ranges, event, branch, overlap_threshold):

    cur_events_overlaps = region.coverage(cur_events_ranges).as_df()
    cur_events_overlaps = cur_events_overlaps.loc[cur_events_overlaps['FractionOverlaps']
                                                > overlap_threshold]
    cur_events_overlaps = cur_events_overlaps.set_index(['Chromosome', 'Start', 'End'])
    cur_events_overlaps['event'] = event
    cur_events_overlaps['branch'] = branch

    return cur_events_overlaps


class MEDICCError(Exception):
    pass
