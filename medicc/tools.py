import re

import Bio
import fstlib
import numpy as np
import pandas as pd


def set_sequences_on_tree_from_df(tree, df, clear_before=True):
    if not hasattr(tree.root, 'sequences'):
        tree = tree.as_phyloxml()
    for clade in tree.find_clades():
        if clear_before:
            clade.sequences.clear()
        for label, data in df.iteritems():
            try:
                clade.sequences.append(
                    Bio.Phylo.PhyloXML.Sequence(
                        name='X'.join(data.groupby('chrom').apply(lambda x: ''.join(x))),
                        symbol=label.upper())) 
            except KeyError:
                pass

def set_sequences_on_tree(tree, fsa_dicts, allele_labels, clear_before=True): 
    """LEGCAY - treats alleles separately"""
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
        chrcats = ["chr%d" % i for i in chrlevels]
        if 'chrX' in newchr:
            chrcats += ['chrX',]
        if 'chrY' in newchr:
            chrcats += ['chrY',]
        newchr = pd.Categorical(newchr, categories=chrcats)
    else:
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
