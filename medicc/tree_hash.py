import copy
from Bio import Phylo
import hashlib
import io


def canonicalize_clade(clade):
    """Recursively sort child clades by subtree structure."""
    if clade.is_terminal():
        return ("leaf", str(clade.name) if clade.name else "")

    # Build canonical signatures bottom-up, then sort children by signature.
    child_signatures = [(canonicalize_clade(subclade), subclade) for subclade in clade.clades]
    child_signatures.sort(key=lambda item: item[0])
    clade.clades = [subclade for _, subclade in child_signatures]
    return ("internal", tuple(signature for signature, _ in child_signatures))


def strip_internal_node_names(tree):
    """Return a copy of a tree where internal clade names are removed."""
    tree_no_internal_names = copy.deepcopy(tree)
    for clade in tree_no_internal_names.find_clades():
        if not clade.is_terminal():
            clade.name = None
    return tree_no_internal_names


def get_canonical_newick(tree):
    """Return a canonical Newick string for topology-based comparisons."""
    tree_copy = strip_internal_node_names(tree)
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
