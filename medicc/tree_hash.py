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



# ---------- Hashing with branch lengths (used for after SPR finished branch length calculations) ----------
# ---------- Prevents functioonally equivalent trees (trees are the same after collapsing zero-length branches) from being treated as different trees ----------

def get_collapsed_tree(tree):
    """
    Takes a Bio.Phylo tree, creates a deep copy, collapses internal branches 
    with length 0.0, and returns the newly modified tree object.
    """
    # 1. Create a deep copy so the original tree remains completely untouched
    new_tree = copy.deepcopy(tree)
    
    # 2. Safely collapse 0-length internal branches
    collapsed_something = True
    while collapsed_something:
        collapsed_something = False
        # Iterate over internal nodes only
        for clade in new_tree.get_nonterminals():
            # Skip the root node
            if clade == new_tree.root:
                continue
            
            # If the branch length is exactly 0.0, collapse it
            if clade.branch_length is not None and clade.branch_length == 0.0:
                new_tree.collapse(clade)
                collapsed_something = True
                break # Break out of the for-loop and restart to ensure safe traversal
                
    return new_tree


def get_topology_hash(tree):
    """
    Generates a unique, rotation-invariant MD5 hash for a Bio.Phylo tree's topology.
    Branch lengths and internal node names are ignored.
    """
    
    def get_canonical_string(clade):
        # 1. Base case: If it's a leaf, just return its name
        if clade.is_terminal():
            return str(clade.name)
        
        # 2. Recursive case: Get the canonical string for all children
        child_strings = [get_canonical_string(child) for child in clade.clades]
        
        # 3. Sort the children alphabetically to make it rotation-invariant
        child_strings.sort()
        
        # 4. Wrap them in parentheses to maintain the tree structure
        return "(" + ",".join(child_strings) + ")"

    # Generate the sorted, standardized string starting from the root
    canonical_representation = get_canonical_string(tree.root)
    
    # Hash the canonical string using MD5 (fast and consistent across Python runs)
    return hashlib.md5(canonical_representation.encode('utf-8')).hexdigest()



def hash_tree_with_collapse(tree):
    """
    Combines the collapsing and hashing steps into one function for convenience.
    """
    collapsed_tree = get_collapsed_tree(tree)
    return get_topology_hash(collapsed_tree)