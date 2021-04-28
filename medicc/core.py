import copy
import logging
from itertools import combinations
from typing import Dict, List

import Bio
import fstlib
import numpy as np
import pandas as pd

import medicc
from medicc import io, nj, stats, tools

# prepare logger 
logger = logging.getLogger(__name__)


def main(input_df,
         asymm_fst,
         normal_name='diploid',
         input_tree=None,
         ancestral_reconstruction=True,
         chr_separator='X'):
    """ MEDICC Main Method """

    symbol_table = asymm_fst.input_symbols()

    ## Validate input
    logger.info("Validating input.")
    io.validate_input(input_df, symbol_table)

    ## Compile input data into FSAs stored in dictionaries
    logger.info("Compiling input sequences into FSAs.")
    FSA_dicts = [create_standard_fsa_dict_from_allele_column(input_df[c], symbol_table, chr_separator) for c in input_df]

    ## Calculate pairwise distances
    logger.info("Calculating pairwise distance matrices for both alleles")
    sample_labels = input_df.index.get_level_values('sample_id').unique()
    pdms = {allele:calc_pairwise_distance_matrix(asymm_fst, fsa_dict) for allele, fsa_dict in zip(input_df.columns, FSA_dicts)}
    pdms['total'] = sum(pdms.values())

    ## Reconstruct a tree
    if input_tree is None:
        logger.info("Inferring tree topology.")
        nj_tree = infer_tree_topology(pdms['total'], sample_labels, diploid = normal_name)
    else:
        logger.info("Tree provided, using it.")
        nj_tree = input_tree
    
    tools.set_sequences_on_tree(nj_tree, FSA_dicts, input_df.columns)

    final_tree = copy.deepcopy(nj_tree)

    if ancestral_reconstruction:
        logger.info("Reconstructing ancestors.")
        ancestors = [medicc.reconstruct_ancestors(tree=final_tree,
                                                    samples_dict=fsa_dict,
                                                    model=asymm_fst,
                                                    normal_name=normal_name)
                    for fsa_dict in FSA_dicts] 

        ## Update branch lengths with ancestors
        logger.info("Updating branch lengths of final tree using ancestors.")
        tools.set_sequences_on_tree(final_tree, ancestors, input_df.columns)
        update_branch_lengths(final_tree, asymm_fst, ancestors, normal_name)
        
        ## Create and write output data frame with ancestors
        logger.info("Creating output table.")
        output_df = create_df_from_fsa_dicts(input_df, ancestors)
    else:
        output_df = input_df

    nj_tree.root_with_outgroup(normal_name)
    final_tree.root_with_outgroup(normal_name)

    if ancestral_reconstruction:
        output_df = summarize_changes(output_df, final_tree, normal_name=normal_name)

    return sample_labels, pdms, nj_tree, final_tree, output_df


def summarize_changes(input_df, input_tree, normal_name=None,
                         ignore_segment_lengths=False):
    df = input_df.copy()

    ## we're force converting to categoricals to always maintain the order of the chromosomes as given
    if not pd.api.types.is_categorical_dtype(df.index.get_level_values('chrom')):
        df.reset_index('chrom', inplace=True)
        df['chrom'] = pd.Categorical(df['chrom'], categories=df['chrom'].unique())
        df.set_index('chrom', inplace=True, append=True)
        df = df.reorder_levels(['sample_id', 'chrom', 'start', 'end'])

    if ignore_segment_lengths:
        df.reset_index(['start', 'end'], inplace=True)

        def reset_start_end(x):
            x['start'] = np.arange(x.shape[0])+1
            x['end'] = np.arange(x.shape[0])+1
            return x

        df = df.groupby(['sample_id', 'chrom']).apply(reset_start_end)
        df.set_index(['start', 'end'], append=True, inplace=True)

    ## test if region is fully conserved
    df.columns.name = 'allele'
    df = df.unstack('sample_id').stack('allele')
    if normal_name is not None:
        is_normal = df.apply(lambda x: (x.loc[normal_name] == x).all(), axis=1).unstack(
            'allele').apply(lambda x: x[0] and x[1], axis=1)
    else:
        is_normal = df.apply(lambda x: (1 == x).all(), axis=1).unstack(
            'allele').apply(lambda x: x[0] and x[1], axis=1)
    is_clonal = df.drop(normal_name, axis=1).apply(lambda x: (x.iloc[0] == x).all(), axis=1).unstack(
        'allele').apply(lambda x: x[0] and x[1], axis=1)

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

    if input_tree is not None:
        dfderiv = compute_change_events(df[input_df.columns], input_tree)
        df.loc[:, 'is_gain'] = dfderiv.apply(lambda x: (x > 0).any(), axis=1)
        df.loc[:, 'is_loss'] = dfderiv.apply(lambda x: (x < 0).any(), axis=1)
    else:
        df['is_gain'] = False
        df['is_loss'] = False

    return df


def create_standard_fsa_dict_from_allele_column(input_column: pd.Series, symbol_table: fstlib.SymbolTable, separator: str = "X") -> dict:
    """ Creates a dictionary of FSAs from a single column/allele (Pandas Series) of the input data frame.
    The keys of the dictionary are the sample/taxon names. """
    fsa_dict = {}
    for taxon, cnp in input_column.groupby('sample_id'):
        cn_str = separator.join(["".join(x.astype('str')) for _,x in cnp.groupby('chrom')])
        fsa_dict[taxon] = fstlib.factory.from_string(cn_str, 
                                arc_type = "standard", 
                                isymbols = symbol_table, 
                                osymbols = symbol_table)

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
    output_df = medicc.create_df_from_fsa_dicts(input_df, [fsa_dict_a, fsa_dict_b], separator)

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

def create_df_from_fsa_dicts(input_df: pd.DataFrame, fsa_dicts: List[dict], separator: str = 'X'):
    """ Takes a list of FSA dicts where each entry corresponds to one allele and extracts the CNPs. 
    The allele names are taken from the input_df columns and the retured data frame has the same 
    number of rows and row index as the input_df. """
    alleles = input_df.columns
    output_index = input_df.reset_index('sample_id').index.unique()

    result = {}
    for allele, fsa_dict in zip(alleles, fsa_dicts):
        for sample in fsa_dict:
            cn = list(tools.fsa_to_string(fsa_dict[sample]).replace(separator, ''))
            result[(allele, sample)] = cn

    output_df = pd.DataFrame(result, index=output_index)
    output_df.columns.names = ['allele', 'sample_id']
    output_df = output_df.stack('sample_id').reorder_levels(['sample_id', 'chrom', 'start', 'end'])

    return output_df

# Given a symmetric model FST and input FSAs in a form of a dictionary, output pairwise distance matrix
def calc_pairwise_distance_matrix(model_fst, fsa_dict):
    samples = list(fsa_dict.keys())
    nsamples = len(samples)

    pwd = np.zeros((nsamples, nsamples))
    combs = list(combinations(list(range(nsamples)),2))

    ncombs = len(combs)
    for i, c in enumerate(combs):
        key1 = samples[c[0]]
        key2 = samples[c[1]]
        fsa1 = fsa_dict[key1]
        fsa2 = fsa_dict[key2]

        d = float(fstlib.kernel_score(model_fst, fsa1, fsa2))

        # Put into a pairwise distance matrix
        key1_idx = samples.index(key1)
        key2_idx = samples.index(key2)
        pwd[key1_idx][key2_idx] = d
        pwd[key2_idx][key1_idx] = d
        progress = (i+1)/ncombs * 100
        if i % 10 == 0 or i==(len(combs)-1): ## every 10 calculations or at the end
            logger.info('%.2f%%', progress)        

    return pwd

def infer_tree_topology(pdm, labels, diploid):
    tree = nj.NeighbourJoining(pdm, labels).tree
    
    input_tree = Bio.Phylo.BaseTree.copy.deepcopy(tree)
    tmpsearch = [c for c in input_tree.find_clades(name = diploid)]
    diploid = tmpsearch[0]
    root_path = input_tree.get_path(diploid)[::-1]

    if len(root_path)>1:
        new_root = root_path[1]
        input_tree.root_with_outgroup(new_root)

    ## from mythic: nj.tree.root_with_outgroup([{'name':s} for s in normal_samples], outgroup_branch_length=0)
    return input_tree

def update_branch_lengths(tree, fst, fsa_dicts, diploid):
    """ Updates the branch lengths in the tree using the internal nodes supplied in the FSA dicts """
    for clade in tree.find_clades():
        children = clade.clades
        if len(children) != 0:
            for child in children:
                if child.name == diploid: ## exception: evolution goes from diploid to internal node
                    brs = [float(fstlib.score(fst, fsa_dict[child.name], fsa_dict[clade.name]))
                        for fsa_dict in fsa_dicts]
                else:
                    brs = [float(fstlib.score(fst, fsa_dict[clade.name], fsa_dict[child.name]))
                        for fsa_dict in fsa_dicts]
                child.branch_length = np.sum(brs)

def compute_change_events(df, tree, normal_name='diploid'):
    dfderiv = df.copy()
    alleles = dfderiv.columns
    for c in alleles:
        dfderiv[c] = dfderiv[c].astype('int')

    clades = [clade for clade in tree.find_clades(order = "postorder") if clade.name is not None and clade.name != normal_name]
    for clade in clades:
        for child in clade.clades:
            dfderiv.loc[child.name, alleles] = dfderiv.loc[child.name, alleles].values - dfderiv.loc[clade.name, alleles].values
    dfderiv.loc[clades[-1].name, alleles] = dfderiv.loc[clades[-1].name, alleles].values - dfderiv.loc[normal_name, alleles].values
    dfderiv.loc[normal_name, alleles] = 0

    return dfderiv

def summarise_patient(tree, pdm, sample_labels, normal_name):
    branch_lengths = []
    for parent in tree.find_clades(terminal=False, order="level"):
        for child in parent.clades:
            if child.branch_length:
                branch_lengths.append(child.branch_length)
    
    nsamples=len(sample_labels)
    tree_length = np.sum(branch_lengths)
    avg_branch_length = np.mean(branch_lengths)
    min_branch_length = np.min(branch_lengths)
    max_branch_length = np.max(branch_lengths)
    median_branch_length = np.median(branch_lengths)
    p_star = stats.star_topology_test(pdm)
    normal_index = np.flatnonzero(np.array(sample_labels) == normal_name)[0]
    p_clock = stats.molecular_clock_test(pdm, normal_index)
    result = pd.Series({
        'nsamples':nsamples,
        'normal_name':normal_name,
        'tree_length':tree_length,
        'mean_branch_length':avg_branch_length,
        'median_branch_length':median_branch_length,
        'min_branch_length':min_branch_length,
        'max_branch_length':max_branch_length,
        'p_star':p_star,
        'p_clock':p_clock
    })
    
    return result

class MEDICCError(Exception):
    pass
