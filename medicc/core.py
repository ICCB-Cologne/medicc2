import copy
import logging
import os
from itertools import combinations
from pathlib import Path

import Bio
import fstlib
import numpy as np
import pandas as pd
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
        pdms = {'total': parallelization_calc_pairwise_distance_matrix(sample_labels, 
                                                                       asymm_fst,
                                                                       FSA_dict,
                                                                       n_cores)}
    else:
        pdms = {'total': calc_pairwise_distance_matrix(asymm_fst, FSA_dict)}

    ## Reconstruct a tree
    if input_tree is None:
        logger.info("Inferring tree topology.")
        nj_tree = infer_tree_topology(
            pdms['total'].values, pdms['total'].index, normal_name=normal_name)
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
        logger.info("Creating output table.")
        output_df = create_df_from_fsa(input_df, ancestors)

        ## Update branch lengths with ancestors
        logger.info("Updating branch lengths of final tree using ancestors.")
        tools.set_sequences_on_tree_from_df(final_tree, output_df)
        update_branch_lengths(final_tree, asymm_fst, ancestors, normal_name)

    nj_tree.root_with_outgroup(normal_name)
    final_tree.root_with_outgroup(normal_name)

    if ancestral_reconstruction:
        output_df, events_df = calculate_all_events(
            final_tree, output_df, allele_columns, normal_name)
        if len(events_df) != final_tree.total_branch_length():
            logger.warn("Event recreation was faulty. Events in '_cn_events_df.tsv' might be incorrect")
    else:
        events_df = None
        output_df = input_df


    return sample_labels, pdms, nj_tree, final_tree, output_df, events_df


def main_legacy(input_df,
                asymm_fst,
                normal_name='diploid',
                input_tree=None,
                ancestral_reconstruction=True,
                chr_separator='X',
                prune_weight=0,
                allele_columns=['cn_a', 'cn_b'],
                n_cores=None):
    """ MEDICC Main Method 
    LEGACY VERSION: The alleles are treated separately in the WGD step"""

    symbol_table = asymm_fst.input_symbols()

    ## Validate input
    logger.info("Validating input.")
    io.validate_input(input_df, symbol_table)

    ## Compile input data into FSAs stored in dictionaries
    logger.info("Compiling input sequences into FSAs.")
    FSA_dicts = [create_standard_fsa_dict_from_data(input_df[c], symbol_table, chr_separator) 
                    for c in input_df]

    ## Calculate pairwise distances
    logger.info("Calculating pairwise distance matrices for both alleles")
    sample_labels = input_df.index.get_level_values('sample_id').unique()
    if n_cores is not None and n_cores > 1:
        pdms = {allele: parallelization_calc_pairwise_distance_matrix(sample_labels, asymm_fst, fsa_dict, n_cores)
                    for allele, fsa_dict in zip(input_df.columns, FSA_dicts)}
    else:
        pdms = {allele: calc_pairwise_distance_matrix(asymm_fst, fsa_dict) 
                    for allele, fsa_dict in zip(input_df.columns, FSA_dicts)}
    pdms['total'] = sum(pdms.values())

    ## Reconstruct a tree
    if input_tree is None:
        logger.info("Inferring tree topology.")
        nj_tree = infer_tree_topology(pdms['total'].values, sample_labels, normal_name=normal_name)
    else:
        logger.info("Tree provided, using it.")
        nj_tree = input_tree
    
    tools.set_sequences_on_tree(nj_tree, FSA_dicts, input_df.columns)

    final_tree = copy.deepcopy(nj_tree)

    if ancestral_reconstruction:
        logger.info("Reconstructing ancestors.")
        ancestors = [medicc.reconstruct_ancestors(tree=final_tree,
                                                  samples_dict=fsa_dict,
                                                  fst=asymm_fst,
                                                  normal_name=normal_name,
                                                  prune_weight=prune_weight)
                     for fsa_dict in FSA_dicts]

        ## Create and write output data frame with ancestors
        logger.info("Creating output table.")
        output_df = create_df_from_fsa(input_df, ancestors)

        ## Update branch lengths with ancestors
        logger.info("Updating branch lengths of final tree using ancestors.")
        tools.set_sequences_on_tree(final_tree, ancestors, input_df.columns)
        update_branch_lengths(final_tree, asymm_fst, ancestors, normal_name)
        
    else:
        output_df = input_df

    nj_tree.root_with_outgroup(normal_name)
    final_tree.root_with_outgroup(normal_name)

    if ancestral_reconstruction:
        output_df = summarize_changes_LEGACY(
            output_df, final_tree, normal_name=normal_name, allele_columns=allele_columns)

    return sample_labels, pdms, nj_tree, final_tree, output_df


# TODO replace with new events_df (also that doesn't use summarize_cn_changes)
def summarize_changes_LEGACY(input_df,
                      input_tree,
                      normal_name='diploid',
                      allele_specific=False,
                      calc_wgd=True,
                      chr_separator='X',
                      asymm_fst_nowgd=None,
                      allele_columns=['cn_a', 'cn_b']):
    df = input_df.copy()

    ## we're force converting to categoricals to always maintain the order of the chromosomes as given
    if not pd.api.types.is_categorical_dtype(df.index.get_level_values('chrom')):
        df.reset_index('chrom', inplace=True)
        df['chrom'] = pd.Categorical(df['chrom'], categories=df['chrom'].unique())
        df.set_index('chrom', inplace=True, append=True)
        df = df.reorder_levels(['sample_id', 'chrom', 'start', 'end'])

    ## test if region is fully conserved
    df.columns.name = 'allele'
    df = df.unstack('sample_id').stack('allele')
    is_normal = df.apply(lambda x: (x.loc[normal_name] == x).all(), axis=1).unstack(
        'allele').apply(lambda x: np.all(x), axis=1)
    is_clonal = df.drop(normal_name, axis=1).apply(lambda x: (x.iloc[0] == x).all(), axis=1).unstack(
        'allele').apply(lambda x: np.all(x), axis=1)

    for a in df:
        df[a] = df[a].astype(int)
    df = df.unstack('allele').stack('sample_id')
    df = df.reorder_levels(['sample_id', 'chrom', 'start', 'end'])
    df.sort_index(inplace=True)
    cats = df.index.get_level_values('chrom').categories
    df = df.join(is_clonal.to_frame('is_clonal').join(is_normal.to_frame('is_normal')))

    ## now work around pandas bug of dropping categoricals
    df.reset_index(inplace=True)
    df.loc[:, 'chrom'] = pd.Categorical(df['chrom'], categories=cats)
    df.set_index(input_df.index.names, inplace=True)
    df.sort_index(inplace=True)

    df['is_loh'] = False
    df['is_wgd'] = False

    if input_tree is not None:
        cn_changes = compute_cn_change(df[input_df.columns], input_tree, normal_name)
        df.loc[:, 'is_gain'] = np.any(cn_changes.values > 0, axis=1)
        df.loc[:, 'is_loss'] = np.any(cn_changes.values < 0, axis=1)
        df.loc[np.logical_and(cn_changes[allele_columns] < 0,
                              df[allele_columns] == 0).any(axis=1), 'is_loh'] = True
        if allele_specific:
            df.loc[:, 'is_gain_a'] = cn_changes['cn_a'].values > 0
            df.loc[:, 'is_loss_a'] = cn_changes['cn_a'].values < 0
            df.loc[:, 'is_gain_b'] = cn_changes['cn_b'].values > 0
            df.loc[:, 'is_loss_b'] = cn_changes['cn_b'].values < 0
    else:
        df['is_gain'] = False
        df['is_loss'] = False
        if allele_specific:
            df.loc[:, 'is_gain_a'] = False
            df.loc[:, 'is_loss_a'] = False
            df.loc[:, 'is_gain_b'] = False
            df.loc[:, 'is_loss_b'] = False

    if input_tree is not None and calc_wgd:

        if asymm_fst_nowgd is None:
            asymm_fst_nowgd = medicc.io.read_fst(no_wgd=True)

        # To save time only potential candidates are investigated further
        wgd_candidate_threshold = 0.3
        df['width'] = df.eval('end-start')
        fraction_gain = (df['is_gain'].astype(int) * df['width']
                         ).groupby('sample_id').sum() / df.loc[df.index.get_level_values('sample_id')[0], 'width'].sum()
        wgd_candidates = list(fraction_gain.index[(fraction_gain > wgd_candidate_threshold)])

        for candidate in wgd_candidates:
            if len(input_tree.get_path(candidate)) == 1:
                parent = normal_name
            else:
                parent = input_tree.get_path(candidate)[-2].name
            cur_FSA_dict = create_standard_fsa_dict_from_data(
                df.loc[[candidate, parent], allele_columns], asymm_fst_nowgd.input_symbols(), chr_separator)
            if float(fstlib.score(asymm_fst_nowgd, cur_FSA_dict[parent], cur_FSA_dict[candidate])) != list(input_tree.find_clades(candidate))[0].branch_length:
                df.loc[candidate, 'is_wgd'] = True
                df.loc[candidate, 'is_gain'] = False
                df.loc[candidate, 'is_loss'] = False
                if allele_specific:
                    df.loc[candidate, 'is_gain_a'] = False
                    df.loc[candidate, 'is_loss_a'] = False
                    df.loc[candidate, 'is_gain_b'] = False
                    df.loc[candidate, 'is_loss_b'] = False

                # calculate subsequent losses and gains
                cur_change = df.loc[candidate, allele_columns] - \
                    (df.loc[parent, allele_columns] + 1)
                df.loc[candidate, 'is_gain'] = np.any(cur_change.values > 0, axis=1)
                df.loc[candidate, 'is_loss'] = np.any(cur_change.values < 0, axis=1)
                if allele_specific:
                    df.loc[candidate, 'is_gain_a'] = cur_change['cn_a'].values > 0
                    df.loc[candidate, 'is_loss_a'] = cur_change['cn_a'].values < 0
                    df.loc[candidate, 'is_gain_b'] = cur_change['cn_b'].values > 0
                    df.loc[candidate, 'is_loss_b'] = cur_change['cn_b'].values < 0
        df.drop('width', axis=1, inplace=True)
    return df


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
    
    phasing_dict = medicc.create_phasing_fsa_dict_from_df(input_df, model_fst.input_symbols(), separator)
    fsa_dict_a, fsa_dict_b, _ = phase_dict(phasing_dict, model_fst, phasing_dict[reference_sample])
    output_df = medicc.create_df_from_fsa(input_df, [fsa_dict_a, fsa_dict_b], separator)

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


def create_df_from_fsa(input_df: pd.DataFrame,
                       fsa,
                       separator: str = 'X'):
    """ 
    Takes a single FSA dict or a list of FSA dicts and extracts the copy number profiles.
    The allele names are taken from the input_df columns and the returned data frame has the same 
    number of rows and row index as the input_df. """

    alleles = input_df.columns
    if isinstance(fsa, dict):
        nr_alleles = len(alleles)
        output_df = input_df.unstack('sample_id')

        for sample in fsa:
            cns = tools.fsa_to_string(fsa[sample]).split(separator)
            if len(cns) % nr_alleles != 0:
                raise MEDICCError('For sample {} we have {} haplotype-specific chromosomes for {} alleles'
                                  '\nnumber of chromosomes has to be divisible by nr of alleles'.format(sample,
                                                                                                        len(cns),
                                                                                                        nr_alleles))
            nr_chroms = int(len(cns) // nr_alleles)
            for i, allele in enumerate(alleles):
                cn = list(''.join(cns[(i*nr_chroms):((i+1)*nr_chroms)]))
                output_df.loc[:, (allele, sample)] = cn

    elif isinstance(fsa, list) and np.all([isinstance(fsa_item, dict) for fsa_item in fsa]):
        output_index = input_df.reset_index('sample_id').index.unique()

        result = {}
        for allele, fsa_dict in zip(alleles, fsa):
            for sample in fsa_dict:
                cn = list(tools.fsa_to_string(fsa_dict[sample]).replace(separator, ''))
                result[(allele, sample)] = cn

        output_df = pd.DataFrame(result, index=output_index)
        output_df.columns.names = ['allele', 'sample_id']

    else:
        raise MEDICCError("fsa input to create_df_from_fsa has to be either dict or list of dicts"
                          "Input type is {}".format(type(fsa)))

    output_df = output_df.stack('sample_id')
    output_df = output_df.reorder_levels(['sample_id', 'chrom', 'start', 'end']).sort_index()
    
    return output_df


def parallelization_calc_pairwise_distance_matrix(sample_labels, asymm_fst, FSA_dict, n_cores):
    parallelization_groups = medicc.tools.create_parallelization_groups(len(sample_labels))
    parallelization_groups = [sample_labels[group] for group in parallelization_groups]
    logger.info("Running {} parallel runs on {} cores".format(len(parallelization_groups), n_cores))

    parallel_pdms = Parallel(n_jobs=n_cores)(delayed(calc_pairwise_distance_matrix)(
        asymm_fst, {key: val for key, val in FSA_dict.items() if key in cur_group}, True)
            for cur_group in parallelization_groups)

    pdm = medicc.tools.total_pdm_from_parallel_pdms(sample_labels, parallel_pdms)

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


def infer_tree_topology(pdm, labels, normal_name):
    if len(labels) > 2:
        tree = nj.NeighbourJoining(pdm, labels).tree

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
    
    The two versions that are selected based on isinstance refer to the current and the legacy version 
    of MEDICC2. In the legacy version the two alleles are treated as separate FSTs and therefore 
    presented as a list.
    """
    if len(ancestor_fsa) == 2:
        child_clade = [x for x in tree.find_clades() if x.name is not None and x.name != normal_name][0]
        child_clade.branch_length = float(fstlib.score(
            fst, ancestor_fsa[normal_name], ancestor_fsa[child_clade.name]))

    if isinstance(ancestor_fsa, dict):
        def _distance_to_child(fst, ancestor_fsa, sample_1, sample_2):
            return float(fstlib.score(fst, ancestor_fsa[sample_1], ancestor_fsa[sample_2]))

    # legacy version
    elif isinstance(ancestor_fsa, list) and np.all([isinstance(ancestor_item, dict) for ancestor_item in ancestor_fsa]):
        def _distance_to_child(fst, ancestor_fsa_list, sample_1, sample_2):
            return np.sum([float(fstlib.score(fst, cur_ancestor_fsa[sample_1], cur_ancestor_fsa[sample_2]))
                           for cur_ancestor_fsa in ancestor_fsa_list])

    else:
        raise MEDICCError("input ancestor_fsa to function update_branch_lengths has to be either a dict"
                          "or a list of dicts\nprovided type is {}".format(type(ancestor_fsa)))

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


def calculate_all_events(tree, cur_df, alleles=['cn_a', 'cn_b'], normal_name='diploid'):
    
    cur_df[['is_gain', 'is_loss', 'is_wgd']] = False
    cur_df[alleles] = cur_df[alleles].astype(int)

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

    return cur_df, events


def calculate_cn_events_per_branch(cur_df, parent_name, child_name, alleles=('cn_a', 'cn_b')):

    asymm_fst, asymm_fst_nowgd, asymm_fst_1_wgd, symbol_table = io.load_main_fsts(
        return_symbol_table=True)

    events_df = pd.DataFrame(columns=['sample_id', 'chrom', 'start', 'end', 'allele', 'type', 'cn_child'])

    cur_parent_cn = cur_df.loc[parent_name, alleles].astype(int)
    cur_child_cn = cur_df.loc[child_name, alleles].astype(int)
    cur_chroms = cur_df.loc['diploid'].index.get_level_values(
        'chrom').map(lambda x: int(x.split('chr')[-1])).values.astype(int)

    # 1. find total losss (tloss)
    parent_tloss = cur_parent_cn == 0
    for allele in alleles:

        cur_tloss = cur_child_cn.loc[~parent_tloss[allele], allele] == 0
        if cur_tloss.sum() == 0:
            continue

        cur_df.loc[child_name, 'is_loss'] = (cur_df.loc[child_name, 'is_loss'].values 
                                             + np.logical_and(cur_child_cn[allele] == 0, 
                                                              cur_parent_cn[allele] != 0).values)

        max_previous_cn = np.max(
            np.unique(cur_parent_cn.loc[~parent_tloss[allele], allele].loc[cur_tloss]))

        for _ in np.arange(max_previous_cn):
            cur_tloss_and_parental_val = np.logical_and(cur_tloss.values, 
                                                        cur_parent_cn.loc[~parent_tloss[allele], allele] > 0).values

            event_labels_ = ((np.cumsum(np.concatenate([[0], np.diff(
                (cur_tloss_and_parental_val + cur_chroms[~parent_tloss[allele]]))])
                * cur_tloss_and_parental_val) + 1)
                * cur_tloss_and_parental_val)

            event_labels = np.zeros_like(event_labels_)
            for i, j in enumerate(np.unique(event_labels_)):
                event_labels[event_labels_ == j] = i

            cur_parent_cn.loc[parent_tloss.loc[~parent_tloss[allele],
                                               allele].index[cur_tloss_and_parental_val], allele] -= 1

            cur_events = (cur_parent_cn
                        .loc[~parent_tloss[allele]]
                        .reset_index()
                        .loc[np.array([np.argmax(event_labels == ind) for ind in np.setdiff1d(np.unique(event_labels), [0])])]
                        [['chrom', 'start', 'end']].values)
            # adjust ends
            cur_events[:, 2] = (cur_parent_cn
                                .loc[~parent_tloss[allele]]
                                .reset_index()
                                .loc[np.array([len(event_labels) - np.argmax(event_labels[::-1] == ind) - 1 for ind in np.setdiff1d(np.unique(event_labels), [0])])]
                                ['end'].values)

            cur_ind = np.arange(len(events_df), len(events_df)+len(cur_events))
            events_df = events_df.append(pd.DataFrame(index=cur_ind))
            events_df.loc[cur_ind, 'sample_id'] = child_name
            events_df.loc[cur_ind, 'allele'] = allele
            events_df.loc[cur_ind, 'type'] = 'tloss'
            events_df.loc[cur_ind, 'cn_child'] = 0
            events_df.loc[cur_ind, ['chrom', 'start', 'end']] = cur_events[:, :3]

            # recalculate parental_loss and cur_tloss for next iteration
            parent_tloss = cur_parent_cn <= 0
            cur_tloss = cur_child_cn.loc[~parent_tloss[allele], allele] == 0

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
    tloss_pos = (cur_parent_cn == 0)
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

            cur_change_location = ((cur_child_cn.loc[~tloss_pos[allele], allele] -
                        cur_parent_cn.loc[~tloss_pos[allele], allele]) == cur_cn_change)

            event_labels_ = ((np.cumsum(np.concatenate([[0], np.diff(
                (cur_change_location.values + cur_chroms[~tloss_pos[allele]]))])
                * cur_change_location.values) + 1)
                * cur_change_location.values)

            event_labels = np.zeros_like(event_labels_)
            for i, j in enumerate(np.unique(event_labels_)):
                event_labels[event_labels_ == j] = i

            cur_events = (cur_child_cn
                          .loc[~tloss_pos[allele]]
                          .reset_index()
                          .loc[np.array([np.argmax(event_labels == val) for val in np.setdiff1d(np.unique(event_labels), [0])])]
                          [['chrom', 'start', 'end', allele]].values)
            # adjust ends
            cur_events[:, 2] = (cur_child_cn
                                .loc[~tloss_pos[allele]]
                                .reset_index()
                                .loc[np.array([len(event_labels) - np.argmax(event_labels[::-1] == val) - 1 for val in np.setdiff1d(np.unique(event_labels), [0])])]
                                ['end'].values)

            cur_ind = np.arange(len(events_df), len(events_df)+len(cur_events))
            events_df = events_df.append(pd.DataFrame(index=cur_ind))
            events_df.loc[cur_ind, 'sample_id'] = child_name
            events_df.loc[cur_ind, 'allele'] = allele
            events_df.loc[cur_ind, 'type'] = cur_event
            events_df.loc[cur_ind, 'cn_child'] = cur_events[:, 3]
            events_df.loc[cur_ind, ['chrom', 'start', 'end']] = cur_events[:, :3]

            cur_child_cn.loc[np.intersect1d(tloss_pos.loc[~tloss_pos[allele]].index,
                                            cur_change_location.loc[cur_change_location].index), allele] += (1 if (cur_cn_change < 0) else -1)

    events_df['chrom'] = tools.format_chromosomes(events_df['chrom'])
    events_df = (events_df[['sample_id', 'allele', 'chrom', 'start', 'end', 'type', 'cn_child']]
                 .reset_index(drop=True)
                 .sort_values(['sample_id', 'allele', 'chrom', 'start', 'end', 'type', 'cn_child']))

    return cur_df, events_df


def compute_cn_change(df, tree, normal_name='diploid'):
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
    p_star = stats.star_topology_test(pdm)
    normal_index = np.flatnonzero(np.array(sample_labels) == normal_name)[0]
    p_clock = stats.molecular_clock_test(pdm, normal_index)
    if events_df is None:
        wgd_status = "unknown"
    else:
        if "wgd" in events_df['type'].values:
            wgd_status = "WGD on branch" + "and ".join(events_df.loc[events_df['type'] == 'wgd', 'sample_id'])
        else:
            wgd_status = "no WGD"

    result = pd.Series({
        'nsamples':nsamples,
        'normal_name':normal_name,
        'tree_length':tree_length,
        'mean_branch_length':avg_branch_length,
        'median_branch_length':median_branch_length,
        'min_branch_length':min_branch_length,
        'max_branch_length':max_branch_length,
        'p_star':p_star,
        'p_clock':p_clock,
        'wgd_status': wgd_status,
    })
    
    return result


# TODO: make compatible with new events_df
def overlap_events(OLD_events_df=None, df=None, tree=None, overlap_threshold=0.9,
                   chromosome_bed='../medicc/objects/hg19_chromosome_arms.bed', regions_bed=None,
                   replace_loss_with_loh=True, allele_specific=False,
                   replace_both_arms_with_chrom=True):

    # TODO move pyranges to main imports once included in main conda env (and yml file)
    import pyranges as pr

    all_events = pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'name', 'NumberOverlaps',
                                       'FractionOverlaps', 'event', 'branch']).set_index(['Chromosome', 'Start', 'End'])

    if OLD_events_df is None:
        if df is None or tree is None:
            raise MEDICCError("Either events_df or df and tree has to be specified")
        OLD_events_df = summarize_changes_LEGACY(df, tree, allele_specific=allele_specific)

    # Read chromosome regions and other regions
    if chromosome_bed is None and regions_bed is None:
        raise MEDICCError("Either chromosome_bed or regions_bed has to be specified")

    chr_arm_regions = None
    if chromosome_bed is not None:
        chr_arm_regions = medicc.io.read_bed_file(chromosome_bed)
        whole_chromosome = chr_arm_regions.groupby('Chromosome').min()
        whole_chromosome['End'] = chr_arm_regions.groupby('Chromosome')['End'].max()
        whole_chromosome['name'] = whole_chromosome.index
        chr_arm_regions = chr_arm_regions.append(whole_chromosome.reset_index()).sort_values('Chromosome')
        chr_arm_regions = pr.PyRanges(chr_arm_regions)

    regions = None
    if regions_bed is not None:
        regions = []
        if isinstance(regions_bed, list) or isinstance(regions_bed, tuple):
            for f in regions_bed:
                regions.append(pr.PyRanges(medicc.io.read_bed_file(f)))
        else:
            regions.append(pr.PyRanges(medicc.io.read_bed_file(regions_bed)))

    for branch in OLD_events_df.index.get_level_values('sample_id').unique():
        # add wgd
        if OLD_events_df.loc[branch, 'is_wgd'].any():
            all_events = all_events.append(pd.DataFrame([['all', '0', '0', 'WGD', len(OLD_events_df.loc[OLD_events_df.index.get_level_values('sample_id').unique()[0]]), 1.0, 'WGD', branch]],
                                                        columns=['Chromosome', 'Start', 'End', 'name', 'NumberOverlaps', 'FractionOverlaps', 'event', 'branch']).set_index(['Chromosome', 'Start', 'End']))
        
        for event in ['loh', 'gain', 'loss'] if replace_loss_with_loh else ['gain', 'loss']:

            cur_events_ranges = OLD_events_df.loc[branch].reset_index().rename(
                {'chrom': 'Chromosome', 'start': 'Start', 'end': 'End'}, axis=1)
            cur_events_ranges = cur_events_ranges.loc[cur_events_ranges['is_{}'.format(event)]]
            cur_events_ranges = pr.PyRanges(cur_events_ranges)

            # Calculate chromosomal events
            if chr_arm_regions is not None:
                chr_events = overlap_regions(
                    chr_arm_regions, cur_events_ranges, event, branch, overlap_threshold)
                # remove arms if the whole chromosome is in there
                if replace_both_arms_with_chrom and len(chr_events) > 0:
                    chr_events = chr_events[~chr_events['name'].isin(np.concatenate(
                        [[name + 'p', name + 'q'] if ('q' not in name and 'p' not in name) else [] for name in chr_events['name']]))]
                all_events = all_events.append(chr_events)

            # Calculate other events
            if regions is not None:
                for region in regions:
                    chr_events = overlap_regions(
                        region, cur_events_ranges, event, branch, overlap_threshold)
                    all_events = all_events.append(chr_events)

    all_events['final_name'] = all_events['name'].apply(lambda x: x.split(
        'chr')[-1]) + all_events['event'].apply(lambda x: ' +' if x == 'gain' else (' -' if x == 'loss' else (' loh' if x == 'loh' else '')))

    all_events.set_index(['branch', 'name'], inplace=True)

    if replace_loss_with_loh:
        all_events = all_events.loc[np.logical_or(all_events['event'] != 'loss',
                                                  ~all_events.index.isin(np.intersect1d(all_events.loc[all_events['event'] == 'loh'].index,
                                                                                        all_events.loc[all_events['event'] == 'loss'].index)))]

    all_events = all_events.reset_index().set_index('branch')
    
    return all_events


def overlap_regions(region, cur_events_ranges, event, branch, overlap_threshold):

    cur_events_overlaps = region.coverage(cur_events_ranges).as_df()
    cur_events_overlaps = cur_events_overlaps.loc[cur_events_overlaps['FractionOverlaps']
                                                > overlap_threshold]
    cur_events_overlaps.set_index(['Chromosome', 'Start', 'End'], inplace=True)
    cur_events_overlaps['event'] = event
    cur_events_overlaps['branch'] = branch

    return cur_events_overlaps

class MEDICCError(Exception):
    pass
