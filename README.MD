# MEDICC2 - Minimum Event Distance for Intra-tumour Copy number Comparisons#

### Version 0.1 beta (Feb 2021). 

## Download: ##

* Clone the MEDICC2 repository using "git clone https://bitbucket.org/rfs/medicc2.git"

## Installation: ##

* Run **python setup.py build_ext --inplace** to compile the fstlib C extensions.

## Dependencies: ##
All dependencies including OpenFST should be directly installable via conda. There is typically no need for manual compilation. Tested on Unix and Windows WSL2.

* Python 3.9
* Mandatory: OpenFST 1.8.0, Biopython 1.78, Numpy 1.20.1, SciPy 1.6.0, Pandas 1.2.2, Matplotlib 3.3.4, Seaborn 0.11.1, Cython 0.29.21
* Optional: Scikit-learn 0.24.1

## Running & Testing: ##

Use **python medicc2.py --help** for general usage and info.

## Contact: ##
Email questions, feature requests and bug reports to **Roland Schwarz, roland.schwarz@mdc-berlin.de**.

## Please cite: ##

PLoS Comput Biol. 2014 Apr 17;10(4):e1003535. doi: 10.1371/journal.pcbi.1003535.
**Phylogenetic quantification of intra-tumour heterogeneity.**
Schwarz RF, Trinh A, Sipos B, Brenton JD, Goldman N, Markowetz F.

February 2021, Roland Schwarz
