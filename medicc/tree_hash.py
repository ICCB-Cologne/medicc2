import copy
from Bio import Phylo
import hashlib
import io


def canonicalize_clade(clade):
    """Recursively sort the clades to get a canonical form."""
    if clade.is_terminal():
        return
    clade.clades.sort(key=lambda x: str(x.name) if x.name else '')
    for subclade in clade.clades:
        canonicalize_clade(subclade)


def get_canonical_newick(tree):
    """Return a canonical Newick string for a Biopython tree."""
    # Deepcopy to avoid modifying the original tree
    import copy
    tree_copy = copy.deepcopy(tree)
    canonicalize_clade(tree_copy.root)
    handle = io.StringIO()
    Phylo.write(tree_copy, handle, 'newick')
    newick_str = handle.getvalue().strip()
    return newick_str


def tree_hash(tree):
    canon_newick = get_canonical_newick(tree)
    return hashlib.md5(canon_newick.encode()).hexdigest()





def strip_branch_lengths(tree):
    tree_no_branch_length = copy.deepcopy(tree)
    for clade in tree_no_branch_length.find_clades():
        clade.branch_length = None
    return tree_no_branch_length
