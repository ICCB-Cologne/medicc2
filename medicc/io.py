import logging
import os
from typing import List

import Bio
import fstlib
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from medicc import core, io, tools

matplotlib.use("Agg")
logger = logging.getLogger(__name__)

def read_and_parse_input_data(filename, normal_name='diploid', input_type='tsv', separator='X', allele_columns=['cn_a','cn_b'], maxcn=8):
    ## Read in input data
    if input_type.lower() == "fasta" or input_type.lower() == 'f':
        logger.info("Reading FASTA input.")
        input_df = io.read_fasta_as_dataframe(filename, separator=separator, allele_columns=allele_columns, maxcn=maxcn)
    elif input_type.lower() == "tsv" or input_type.lower() == 't':
        logger.info("Reading Refphase TSV input.")
        input_df = io.read_tsv_as_dataframe(filename, allele_columns=allele_columns, maxcn=maxcn)
    else:
        raise MEDICCIOError("Unknown input type, possible options are FASTA or TSV.")

    ## Add normal sample if needed
    input_df = io.add_normal_sample(input_df, normal_name, allele_columns=allele_columns)
    nsamples = input_df.index.get_level_values('sample_id').unique().shape[0]
    nchr = input_df.index.get_level_values('chrom').unique().shape[0]
    nsegs = input_df.loc[normal_name,:].shape[0]
    logger.info("Read %d samples, %d chromosomes, %d segments per sample", nsamples, nchr, nsegs)
    return input_df

def validate_input(input_df, symbol_table):
    ## Check the number of alleles
    if len(input_df.columns)>2:
        raise MEDICCIOError("More than 2 alleles are currently not supported.")

    if len(input_df.columns)==0:
        raise MEDICCIOError("No alleles found.")

    ## Check if CN states are within bounds

    ## Check if symbols are in symbol table
    alphabet = {x[1] for x in symbol_table}
    data_chars = set(input_df.values.flatten())
    if not data_chars.issubset(alphabet):
        not_in_set = data_chars.difference(alphabet)
        raise MEDICCIOError("Not all input symbols are contained in symbol table. Offending symbols: %s" % str(not_in_set))

    ## Check data type start and end columns
    if (input_df.index.get_level_values('start').dtype != np.int or 
        input_df.index.get_level_values('end').dtype != np.int):
        raise MEDICCIOError("Start and end columns must be of type: integer.")

    ## Check data type payload columns - these should all be of type str (object)
    if not np.all([pd.api.types.is_string_dtype(x) for x in input_df.dtypes]): ## this shouldn't happen
        raise MEDICCIOError("Payload columns must be of type: string.")

    logger.info('Ok!')

def read_tsv_as_dataframe(path, allele_columns=['cn_a','cn_b'], maxcn=8):
    logger.info("Reading TSV file %s", path)
    input_file = pd.read_csv(path, sep = "\t")
    nexpected = 4 + len(allele_columns)
    input_columns = list(input_file.columns[0:4]) + allele_columns
    if input_file.shape[1] < nexpected:
        raise MEDICCIOError("TSV file needs at least %d columns (sample_id, chrom, start, end and the allele columns)!" % nexpected)

    ## check if allele_columns are present
    if not np.all(np.in1d(allele_columns, input_file.columns)):
        raise MEDICCIOError("Allele columns {%s} not found!" % ', '.join(allele_columns))

    logger.info("Successfully read input file. Using columns {%s}" % ', '.join(input_columns))
    new_columns = ['sample_id', 'chrom', 'start', 'end'] + allele_columns ## rename user input to what we know
    input_file = input_file[input_columns]
    input_file.columns = new_columns
    for c in allele_columns:
        if input_file[c].dtype in [np.dtype('float64'), np.dtype('float32')]:
            logger.warning("Floating point payload! I will round, but this might not be intended.")
            input_file[c] = input_file[c].round().astype('int')
        if input_file[c].dtype in [np.dtype('int64'), np.dtype('int32')]:
            if np.any(input_file[c]>maxcn):
                logger.warning("Integer CN > maxcn %d, capping.", maxcn)
                input_file[c] = np.fmin(input_file[c], maxcn)
    input_file['chrom'] = tools.format_chromosomes(input_file['chrom'])
    input_file.set_index(['sample_id', 'chrom', 'start', 'end'], inplace=True)
    input_file.sort_index(inplace=True)
    input_file[allele_columns] = input_file[allele_columns].astype(str)

    return input_file

def read_fasta_as_dataframe(infile: str, separator: str = 'X', allele_columns = ['cn_a','cn_b'], maxcn: int = 8):
    """ Reads FASTA decriptor file (old MEDICC input format) and reads the corresponding FASTA files to generate
    a data frame with the same format as the refphase input (TSV) format. """
    logger.info("Reading FASTA dataset from description file %s.", infile)
    description_file = pd.read_csv(infile,
                                    delim_whitespace = True,
                                    header = None,
                                    names = ["chrom",] + allele_columns,
                                    usecols = ["chrom",] + allele_columns)
    description_file.set_index('chrom', inplace=True)
    description_file.columns.name='allele'
    description_file = description_file.stack()
    inpath = os.path.dirname(infile)

    payload = []
    for filename in description_file:
        with open(os.path.join(inpath, filename), 'r') as fd:
            text = fd.read()
            df = pd.DataFrame([s.strip().split('\n') for s in text.split('>') if s.strip() !=''], columns=['sample_id', 'cnp'])
            df.set_index('sample_id', inplace=True)
            payload.append(df)

    payload = pd.concat(payload, keys=description_file.index, names=description_file.index.names)

    ## deal with the case where multiple chromosomes are given in one string
    payload=payload.cnp.str.split(separator, expand=True)
    if payload.shape[1]>1: ## multiple columns
        payload.columns.name='id'
        payload=payload.unstack('chrom').reorder_levels(['chrom','id'],axis=1)
        payload.columns = ["%s%d" % (c,i+1) for c,i in payload.columns]
        payload.columns.name = 'chrom'
        payload = payload.stack()
        payload = payload.reorder_levels([1,2,0]).sort_index()
    else:
        payload = payload.iloc[:,0]
        payload = payload.reorder_levels([2,0,1]).sort_index()

    payexpand = [pd.DataFrame(list(s), columns=['cn']) for s in payload]
    paylong = pd.concat(payexpand, keys=payload.index, names=payload.index.names + ['segment'])
    paylong['cn'] = paylong['cn'].apply(tools.hex2int)
    if np.any(paylong['cn'] > maxcn):
        logger.warning("Integer CN > maxcn %d, capping.", maxcn)
        paylong['cn'] = np.fmin(paylong['cn'], maxcn)
    paylong['cn'] = paylong['cn'].astype(str)
    result = paylong.unstack(['allele'])
    result = result.droplevel(0, axis=1).reset_index()
    result.columns.name = None
    result['start'] = result['segment']+1
    result['end'] = result['start']
    result = result[['sample_id','chrom','start','end'] + allele_columns]
    result['chrom'] = tools.format_chromosomes(result['chrom'])
    result.sort_values(['sample_id', 'chrom', 'start', 'end'], inplace=True)
    result.set_index(['sample_id', 'chrom', 'start', 'end'], inplace=True)
    result.sort_index(inplace=True)

    return result

def add_normal_sample(df, normal_name, allele_columns=['cn_a','cn_b']):
    """ Adds an artificial normal samples with the supplied name to the data frame.
    The normal sample has CN=1 on all supplied alleles. """
    samples = df.index.get_level_values('sample_id').unique()

    if normal_name is not None and normal_name not in samples:
        logger.info("Normal sample '%s' not found, adding artifical normal by the name: '%s'.", normal_name, normal_name)
        tmp=df.unstack('sample_id')
        for col in allele_columns:
            tmp.loc[:,(col, normal_name)] = '1'
        tmp = tmp.stack('sample_id')
        tmp = tmp.reorder_levels(['sample_id', 'chrom', 'start', 'end']).sort_index()
    else:
        tmp = df

    return tmp


def write_tree_files(tree, out_name: str, plot_tree=True, draw_ascii=True):
    """ Writes a Newick, PhyloXML, Ascii graphic and PNG grahic file of the tree. """
    Bio.Phylo.write(tree, out_name + ".new", "newick")
    Bio.Phylo.write(tree, out_name + ".xml", "phyloxml")

    if draw_ascii:
        with open(out_name + ".txt", "w") as f:
            Bio.Phylo.draw_ascii(tree, file = f)

    if plot_tree:
        fig = plt.figure(figsize = (7, 7))
        axes = fig.add_subplot(1, 1, 1)
        Bio.Phylo.draw(tree, axes = axes, do_show = False)
        plt.savefig(out_name + ".png")

def write_pdms(sample_labels, pdms, filename_prefix):
    """ Writes all PDMs in the pdms dictionary. """
    for allele, pdm in pdms.items():
        _write_pdm(sample_labels, pdm, "%s_%s.tsv" % (filename_prefix, allele))

def _write_pdm(labels, pdm, filename):
    """ Writes a single PDM to the given filename using the provided labels as row and column names. """
    pdm_df = pd.DataFrame(pdm, columns = labels, index = labels)
    pdm_df.to_csv(filename, sep='\t')
    return pdm_df

def import_tree(tree_file, normal_name, file_format='newick'):
    """ Loads a phylogenetic tree in the given format and roots it at the normal sample. """
    tree = Bio.Phylo.read(tree_file, file_format)
    input_tree = Bio.Phylo.BaseTree.copy.deepcopy(tree)
    tmpsearch = [c for c in input_tree.find_clades(name = normal_name)]
    normal_name = tmpsearch[0]
    root_path = input_tree.get_path(normal_name)[::-1]

    if len(root_path) > 1:
        new_root = root_path[1]
        print("Rooting with " + new_root.name)
        input_tree.root_with_outgroup(new_root)
    else:
        pass

    return input_tree

class MEDICCIOError(Exception):
    pass
