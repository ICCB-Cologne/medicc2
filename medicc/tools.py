import logging
import re

import Bio
import fstlib
import numpy as np
import pandas as pd


def set_sequences_on_tree_from_df(tree: Bio.Phylo.BaseTree, df: pd.DataFrame, clear_before=True):
    """Set sequences on tree from dataframe

    Args:
        tree (Bio.Phylo.BaseTree): Tree to set sequences on
        df (pd.DataFrame): DataFrame with copy number information
        clear_before (bool, optional): Clear old sequences. Defaults to True.
    """


    if not hasattr(tree.root, 'sequences'):
        tree = tree.as_phyloxml()
    for clade in tree.find_clades():
        if clear_before:
            clade.sequences.clear()
        for label, data in df.iteritems():
            try:
                clade.sequences.append(
                    Bio.Phylo.PhyloXML.Sequence(
                        name='X'.join(data.loc[clade.name].groupby('chrom').apply(lambda x: ''.join(x))),
                        symbol=label.upper())) 
            except KeyError:
                pass

def set_sequences_on_tree(tree, fsa_dicts, allele_labels, clear_before=True): 
    """LEGACY - treats alleles separately"""
    if not hasattr(tree.root, 'sequences'):
        tree = tree.as_phyloxml()
    for clade in tree.find_clades():
        if clear_before:
            clade.sequences.clear()
        for label, fsa_dict in zip(allele_labels, fsa_dicts):
            try:
                fsa = fsa_dict[clade.name]
                clade.sequences.append(
                    Bio.Phylo.PhyloXML.Sequence(
                        name=fsa_to_string(fsa),
                        symbol=label.upper())
                ) 
            except KeyError:
                pass


def fsa_to_string(fsa):
    fsa_string = fstlib.tools.strings(fsa).string[0]
    return fsa_string

def hex2int(x):
    return int(x, 16)

def int2hex(x):
    return hex(x)[2:]

def format_chromosomes(ds):
    """ Expects pandas Series with chromosome names. 
    The goal is to take recognisalbe chromosome names, i.e. chr4 or chrom3 and turn them into chr3 format.
    If the chromosomes names are not recognized, return them unchanged."""
    ds = ds.astype('str')
    pattern = re.compile(r"(chr|chrom)?((\d+)|X|Y)", flags=re.IGNORECASE)
    matches = ds.apply(pattern.match)
    matchable = ~matches.isnull().any()
    if matchable:
        newchr = matches.apply(lambda x:"chr%s" % x[2].upper())
        numchr = matches.apply(lambda x:int(x[3]) if x[3] is not None else -1)
        chrlevels = np.sort(numchr.unique())
        chrlevels = np.setdiff1d(chrlevels, [-1])
        chrcats = ["chr%d" % i for i in chrlevels]
        if 'chrX' in list(newchr):
            chrcats += ['chrX',]
        if 'chrY' in list(newchr):
            chrcats += ['chrY',]
        newchr = pd.Categorical(newchr, categories=chrcats)
    else:
        logger = logging.getLogger('medicc.io')
        logger.warn("Could not match the chromosome labels. Rename the chromosomes according chr1, "
                    "chr2, ... to avoid potential errors.")
        newchr = pd.Categorical(ds, categories=ds.unique())
    return newchr

def lcp(m):
    """ Returns the longest common prefix of a set of strings (in an iterable). """
    if m is None: 
        return None
    
    s1 = min(m)
    s2 = max(m)
    for i, c in enumerate(s1):
        if c != s2[i]:
            return s1[:i]
    return s1

def format_sample_ids(ds, underscore_replacement=None):
    """ Replaces _ with spaces and inserts \n in the sample id after the longest 
    common prefix. """
    if underscore_replacement is not None:
        output = ds.str.replace('_', underscore_replacement)
    else:
        output = ds.copy()
    prefix = lcp(output)
    if prefix!='':
        if underscore_replacement and prefix.endswith(underscore_replacement):
            output = output.str.replace(prefix, prefix[:-1] + '\n')
        else:
            output = output.str.replace(prefix, prefix + '\n')
    return output


def next_prime(N):
    '''calculate a prime number p with p >= N'''
    def is_prime(x):
        return all(x % i for i in range(3, x))

    cur_p = int(np.ceil(np.sqrt(N)))
    if cur_p % 2 == 0:
        cur_p += 1
    while not is_prime(cur_p):
        cur_p += 2
    return cur_p


def create_parallelization_groups(number_samples):
    '''Create subgroups of the samples to perform parallel runs of the pairwise distances matrix 
    calculations. 

    Args:
        number_samples (int): Total number of samples

    Returns:
        list of lists: Lists, each one containing up to p indices


    Method implemented as proposed in:
    Emmanuel Sapin, Matthew C Keller, Novel approach for parallelizing pairwise comparison problems as applied to detecting segments identical by decent in whole-genome data, 
    Bioinformatics, 2021;, btab084, https://doi.org/10.1093/bioinformatics/btab084
    '''
    p = next_prime(number_samples)
    p_matrix = np.arange(p**2).reshape((p, p)).T
    p_matrix[p_matrix >= number_samples] = -1

    p_matrix_permutations = [p_matrix]
    for _ in range(p-1):
        p_matrix_permutations.append(np.array([np.roll(p_matrix_permutations[-1][:, i], -i) for i in range(p)]).T)

    groups = [*p_matrix.T] + [x for cur_p in p_matrix_permutations for x in [*cur_p]]
    # remove -1 entries (numbers larger than number_samples)
    groups = [group[group >= 0] for group in groups]
    groups = [group for group in groups if len(group) > 1]

    return groups


def total_pdm_from_parallel_pdms(sample_labels, parallel_pdms):
    total_pdm = pd.DataFrame(index=sample_labels, columns=sample_labels, dtype=float)

    for cur_pdm in parallel_pdms:
        total_pdm.loc[cur_pdm.index, cur_pdm.index] = cur_pdm

    if total_pdm.isna().sum().sum() != 0:
        raise ValueError('Something went wrong with these indices:\n{}'.format(
            np.where(total_pdm.isna())))

    return total_pdm
