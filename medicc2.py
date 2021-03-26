import argparse
import logging
import os

import numpy as np

import fstlib
import medicc

logger = logging.getLogger('medicc-main')
objects_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "objects")

parser = argparse.ArgumentParser()
parser.add_argument("input_file", 
                    help = "a path to the input file")
parser.add_argument("output_dir", 
                    help ="a path to the output folder")
parser.add_argument("--input-type", "-i", type = str, dest = "input_type", default = "t", choices = ["f", "t"], required = False, 
                    help = "Choose the type of input: f for FASTA, t for TSV (default: TSV)")
parser.add_argument("--input-allele-columns", "-a",
                    type=str,
                    dest='input_allele_columns',
                    default='cn_a, cn_b',
                    required=False,
                    help = """Name of the CN columns (comma separated) if using TSV input format (default: 'cn_a, cn_b').
                    This also adjusts the number of alleles considered (min. 1, max. 2).""")
parser.add_argument("--input-chr-separator",
                    type=str,
                    dest='input_chr_separator',
                    default = 'X',
                    required=False,
                    help = 'Character used to separate chromosomes in the input data (condensed FASTA only, default: \"X\").')
parser.add_argument("--tree",
                    action = "store",
                    dest = "user_tree",
                    help = "Do not reconstruct tree, use provided tree instead (in newick format) and only perform ancestral reconstruction (default: None).",
                    required = False)
parser.add_argument("--topology-only", "-s",
                    action = "store_true", 
                    dest = "topology_only", 
                    help = "Output only tree topology, without reconstructing ancestors (default: false).",
                    required = False)
parser.add_argument("--normal-name", "-n",
                    default = "diploid", 
                    type = str,
                    dest = "normal_name", 
                    help = """ID of the sample to be treated as the normal sample. 
                    Trees are rooted at this sample for ancestral reconstruction (default: \"diploid\").
                    If the sample ID is not found, an artificial normal sample of the same name is created 
                    with CN states = 1 for each allele.""",
                    required = False)
parser.add_argument("--exclude-samples", "-x",
                    default = None,
                    type = str,
                    help = "Comma separated list of sample IDs to exclude.",
                    required=False)
parser.add_argument("--prefix", '-p', type=str, dest='prefix', default=None, 
                    help='Output prefix to be used (default: input filename).', required=False)
parser.add_argument("--no-wgd", action='store_true', default=False, 
                    help='Enable whole-genome doubling events (default: false).', required=False)
parser.add_argument("--no-plot", action='store_true', default=False, 
                    help='Disable plotting (default: false).', required=False)
parser.add_argument("--summary",
		            dest = "plot_summary",
		            action = "store_true", 
                    required = False, 
                    help = "Plot a track with event summary")
parser.add_argument("--fst", type=str, dest='fst', default=None,
                    help='Expert option: path to an alternative FST.')
parser.add_argument("--fst-chr-separator", type=str, dest='fst_chr_separator', default='X',
                    help = 'Expert option: character used to separate chromosomes in the FST (default: \"X\").')
parser.add_argument("--maxcn", type=int, dest='maxcn', default=8,
                    help='Expert option: maximum CN supported by the supplied FST.')
parser.add_argument("-v", "--verbose", action='store_true', default=False,
                    help='Enable versbose output (default: false).', required=False)

args = parser.parse_args()

if args.verbose:
    logging.getLogger('medicc').setLevel(logging.INFO)

output_dir = args.output_dir
normal_name = args.normal_name 
allele_columns = [x.strip() for x in args.input_allele_columns.split(',')]

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

## Determine prefix for output files 
if args.prefix is None:
    output_prefix = os.path.basename(os.path.splitext(args.input_file)[0])
else:
    output_prefix = args.prefix

## Read in symbol table and transducers
logger.info("Reading FST.")
if args.fst is not None:
    fst_path = args.fst
else:
    if args.no_wgd:
        fst_path = os.path.join(objects_dir, 'no_wgd_asymm.fst')
    else:
        fst_path = os.path.join(objects_dir, 'wgd_asymm.fst')
fst = fstlib.read(fst_path)

if args.user_tree is not None:
    logger.info("Importing user tree.")
    input_tree = medicc.io.import_tree(tree_file = args.user_tree, diploid = normal_name)
else:
    input_tree = None

## Load data
logger.info("Reading and parsing input data.")
input_df = medicc.io.read_and_parse_input_data(
    args.input_file,
    normal_name,
    args.input_type.strip(),
    args.input_chr_separator.strip(),
    allele_columns)

if args.exclude_samples is not None:
    exclude_samples = np.array([x.strip() for x in args.exclude_samples.split(',')])
    logger.info("Excluding samples {%s}." % ', '.join(exclude_samples))
    input_df = input_df.loc[~np.in1d(input_df.index.get_level_values('sample_id'), exclude_samples), :]


## Run main method
logger.info("Running main reconstruction routine.")
sample_labels, pdms, nj_tree, final_tree, output_df = medicc.main(
    input_df, 
    fst, 
    normal_name, 
    input_tree=input_tree, 
    ancestral_reconstruction=not args.topology_only,
    chr_separator=args.fst_chr_separator.strip())

## Output pairwise distance matrices
logger.info("Writing pairwise distance matrices.")
medicc.io.write_pdms(sample_labels, pdms, os.path.join(output_dir, output_prefix + "_pdm"))

## Write trees
logger.info("Writing trees.")
medicc.io.write_tree_files(tree = nj_tree, out_name = os.path.join(output_dir, output_prefix + "_nj_tree"))
medicc.io.write_tree_files(tree = final_tree, out_name = os.path.join(output_dir, output_prefix + "_final_tree"))

## Write ouput table
output_df.to_csv(os.path.join(output_dir, output_prefix + "_final_cn_profiles.tsv"), sep='\t')

## Summarise
logger.info("Writing patient summary.")
summary = medicc.summarise_patient(final_tree, pdms['total'], sample_labels, normal_name)
logger.info("Final tree length %d", summary.tree_length)
summary.to_csv(os.path.join(output_dir, output_prefix + "_summary.tsv"), index=True, header=False, sep='\t')

if not args.no_plot:
    logger.info("Plotting CN profiles.")
    plot_summary = args.plot_summary
    ##p = medicc.plot.plot_cn_profiles(output_df, sample_order=[c.name for c in final_tree.find_clades(order='postorder')])
    #gg.ggsave(p, filename=os.path.join(output_dir, output_prefix + '_summary.pdf'), limitsize=False)
    p = medicc.plot.plot_cn_profiles(
        output_df, 
        final_tree, 
        title=output_prefix, 
        normal_name=normal_name, 
        plot_summary = plot_summary,
        label_func = None)
    p.savefig(os.path.join(output_dir, output_prefix + '_cn_profiles.pdf'))
    
