# MEDICC2 - Minimum Event Distance for Intra-tumour Copy number Comparisons

**Version 0.4.2 beta (May 2021)**

For more information see the accompanying  paper [Whole-genome doubling-aware copy number phylogenies for cancer evolution with MEDICC2](https://www.biorxiv.org/content/10.1101/2021.02.28.433227v1.abstract).

# Setup
## Download

Clone the MEDICC2 repository and its submodules using **git clone --recursive https://bitbucket.org/schwarzlab/medicc2.git**

**It is important to use the `--recursive` flag to also download the OpenFST submodule**
## Installation

Run **python setup.py build_ext --inplace** to compile the fstlib C extensions.

## Dependencies
All dependencies including OpenFST should be directly installable via conda.

**A YML file with a suggested MEDICC2 conda environment is provided in 'doc/medicc2.yml'.**

MEDICC2 contains modified code of the *pywrapfst* Python module from [OpenFST](http://www.openfst.org/) as permitted by the [Apache 2](http://www.apache.org/licenses/LICENSE-2.0) license.

## License
MEDICC2 is available under [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html). 

# Usage

General usage is `python medicc2.py path/to/input/file path/to/output/folder`. Run `python medicc2.py --help` for information on possible flags

Logging settings can be changed using the `logging_conf.yaml` file with the standard python logging syntax.

## Input files
Input files can be either in fasta or tsv format:
* **fasta:** A description file should be provided to MEDICC. This file should include one line per file with the name of the chromosome and the corresponding file names. If fasta files are provided you have to use the flag `--input-type fasta`.
* **tsv:** Files should have the following columns: `sample_id`, `chrom`, `start`, `end` as well as columns for the copy numbers. MEDICC expects the copy number columns to be called `cn_a` and `cn_b`. Using the flag `--input-allele-columns` you can set your own copy number columns. If you want to use total copy numbers, make sure to use the flag `--total-copy-numbers`.

The folder `examples/simple_example` contains a simple example input both in fasta and tsv format.
The folder `examples/OV03-04` contains a larger example consisting of multiple fasta files. If you want to run MEDICC on this data run `python medicc2.py examples/OV03-04/OV03-04_descr.txt path/to/output/folder --input-type fasta`.

## Usage examples
For first time users we recommend to have a look at `examples/simple_example` to get an idea of how input data should look like. Then run `python medicc2.py examples/simple_example/simple_example.tsv path/to/output/folder` as an example of a standard MEDICC run. Finally, the notebook `notebooks/example_workflows.py` shows how the individual functions in the workflow are used.

The notebook `notebooks/bootstrap_demo.py` demonstrates how to use the bootstrapping routine and `notebooks/plot_demo.py` shows how to use the main plotting functions.

# Contact
Email questions, feature requests and bug reports to **Tom Kaufmann, tom.kaufmann@mdc-berlin.de**.

# Please cite

Petkovic M, Watkins TBK, Colliver E, Laskina S, Swanton C, Haase K, Schwarz RF.  
**Whole-genome doubling-aware copy number phylogenies for cancer evolution with MEDICC2.**  
bioRxiv 2021 Feb 28; doi: 10.1101/2021.02.28.433227  

Schwarz RF, Trinh A, Sipos B, Brenton JD, Goldman N, Markowetz F.  
**Phylogenetic quantification of intra-tumour heterogeneity.**  
PLoS Comput Biol. 2014 Apr 17;10(4):e1003535. doi: 10.1371/journal.pcbi.1003535.
