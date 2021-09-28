# MEDICC2 - Whole-genome doubling-aware copy number phylogenies for cancer evolution

For more information see the accompanying  paper [Whole-genome doubling-aware copy number phylogenies for cancer evolution with MEDICC2](https://www.biorxiv.org/content/10.1101/2021.02.28.433227v2).

# Installation
Install MEDICC2 via conda (recommended), pip or from source. MEDICC2 was developed and tested on unix-built systems (Linux and MacOS). For Windows users we recommended WSL2.

Note that the notebooks and examples are not included when installing from conda or pip. When installing from pip or source, you need to make sure to have a working version of `gcc` and `gxx` installed.

## Installation via conda (recommended)
MEDICC2 can be installed via `conda install -c bioconda -c conda-forge medicc2`.

## Installation via pip
As MEDICC2 relies on OpenFST version 1.8.1 which is not packaged on PyPi you have to first install it using conda with `conda install -c conda-forge openfst`. Next you can install MEDICC2 via `pip install medicc2`.

## Installation from source
Clone the MEDICC2 repository and its submodules using `git clone --recursive https://bitbucket.org/schwarzlab/medicc2.git`. It is important to use the `--recursive` flag to also download the modified OpenFST submodule.

All dependencies including OpenFST (v1.8.1) should be directly installable via conda. A yaml file with a suggested MEDICC2 conda environment is provided in 'doc/medicc2.yml'. You can create a new conda environment with all requirements using `conda env create -f doc/medicc2.yml -n medicc_env`.

Then, inside the `medicc2` folder, run `pip install .` to install MEDICC2 to your environment. 

# Usage
After installing MEDICC2, you can use MEDICC2 functions in python scripts (through `import medicc`) and from the command line. General usage from the command line is `medicc2 path/to/input/file path/to/output/folder`. Run `medicc2 --help` for information on optional arguments.

Logging settings can be changed using the `medicc/logging_conf.yaml` file with the standard python logging syntax.

## Command line Flags

* `input_file`: path to the input file
* `output_dir`: path to the output folder
* `--input-type`, `-i`: Choose the type of input: f for FASTA, t for TSV. Default: 'TSV'
* `--input-allele-columns`, `-a`: Name of the CN columns (comma separated) if using TSV input format. This also adjusts the number of alleles considered (min. 1, max. 2). Default: 'cn_a, cn_b'
* `--input-chr-separator`: Character used to separate chromosomes in the input data (condensed FASTA only). Default: 'X'
* `--tree`: Do not reconstruct tree, use provided tree instead (in newick format) and only perform ancestral reconstruction. Default: None
* `--topology-only`, `-s`: Output only tree topology, without reconstructing ancestors. Default: False
* `--normal-name`, `-n`: ID of the sample to be treated as the normal sample. Trees are rooted at this sample for ancestral reconstruction. If the sample ID is not found, an artificial normal sample of the same name is created with CN states = 1 for each allele. Default: 'diploid'
* `--exclude-samples`, `-x`: Comma separated list of sample IDs to exclude. Default: None
* `--filter-segment-length`: Removes segments that are smaller than specified length. Default: None
* `--bootstrap-method`: Bootstrap method. Has to be either 'chr-wise' or 'segment-wise'. Default: 'chr-wise'
* `--bootstrap-nr`: Number of bootstrap runs to perform. Default: None
* `--prefix`, '-p': Output prefix to be used. None uses input filename. Default: None
* `--no-wgd`: Disable whole-genome doubling events. Default: False
* `--no-plot`: Disable plotting. Default: False
* `--legacy-version`: Use legacy version in which alleles are treated separately. Default: False
* `--total-copy-numbers`: Run for total copy number data instead of allele-specific data. Default: False
* `-j`, `--n-cores`: Number of cores to run on. Default: None
* `-v`, `--verbose`: Enable verbose output. Default: False
* `--maxcn`: Expert option: maximum CN at which the input is capped. Does not change FST. Default: 8
* `--prune-weight`: Expert option: Prune weight in ancestor reconstruction. Values >0 might result in more accurate ancestors but will require more time and memory. Default: 0
* `--fst`: Expert option: path to an alternative FST. Default: None
* `--fst-chr-separator`: Expert option: character used to separate chromosomes in the FST. Default: 'X'


## Input files
Input files can be either in fasta or tsv format:
* **fasta:** A description file should be provided to MEDICC. This file should include one line per file with the name of the chromosome and the corresponding file names. If fasta files are provided you have to use the flag `--input-type fasta`.
* **tsv:** Files should have the following columns: `sample_id`, `chrom`, `start`, `end` as well as columns for the copy numbers. MEDICC expects the copy number columns to be called `cn_a` and `cn_b`. Using the flag `--input-allele-columns` you can set your own copy number columns. If you want to use total copy numbers, make sure to use the flag `--total-copy-numbers`.

MEDICC2 follows the BED convention for segment coordinates, i.e. segment start is at 0 and the segment end is non-inclusive.

The folder `examples/simple_example` contains a simple example input both in fasta and tsv format.
The folder `examples/OV03-04` contains a larger example consisting of multiple fasta files. If you want to run MEDICC on this data run `medicc2 examples/OV03-04/OV03-04_descr.txt path/to/output/folder --input-type fasta`.

## Usage examples
For first time users we recommend to have a look at `examples/simple_example` to get an idea of how input data should look like. Then run `medicc2 examples/simple_example/simple_example.tsv path/to/output/folder` as an example of a standard MEDICC run. Finally, the notebook `notebooks/example_workflows.py` shows how the individual functions in the workflow are used.

The notebook `notebooks/bootstrap_demo.py` demonstrates how to use the bootstrapping routine and `notebooks/plot_demo.py` shows how to use the main plotting functions.

# Contact
Email questions, feature requests and bug reports to **Tom Kaufmann, tom.kaufmann@mdc-berlin.de**.

# License
MEDICC2 is available under [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html). It contains modified code of the *pywrapfst* Python module from [OpenFST](http://www.openfst.org/) as permitted by the [Apache 2](http://www.apache.org/licenses/LICENSE-2.0) license.

# Please cite
Kaufmann TL, Petkovic M, Watkins TBK, Colliver EC, Laskina S, Thapa N, Minussi DC, Navin N, Swanton C, Van Loo P, Haase K, Tarabichi M, Schwarz RF.
**MEDICC2: whole-genome doubling aware copy-number phylogenies for cancer evolution**  
bioRxiv 2021 Sep 6; doi: 10.1101/2021.02.28.433227 

Schwarz RF, Trinh A, Sipos B, Brenton JD, Goldman N, Markowetz F.  
**Phylogenetic quantification of intra-tumour heterogeneity.**  
PLoS Comput Biol. 2014 Apr 17;10(4):e1003535. doi: 10.1371/journal.pcbi.1003535.
