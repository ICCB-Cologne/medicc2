"""Unit tests for medicc.nni — NNI move enumeration."""
import copy
import pytest
from Bio.Phylo.BaseTree import Clade, Tree

import medicc.nni as nni
import medicc.spr as spr
import medicc.tree_hash as tree_hash
from medicc.ancestors import _get_dirty_nodes_up


def _make_search_shape_tree(n_leaves):
    """Build a left-comb search-shape tree with n_leaves non-diploid leaves.

    Result shape:
      root("diploid") -> internal_0 -> internal_1 -> ... -> (s1, s2, ..., s_n)

    Specifically a caterpillar where each internal node has one leaf child and
    one internal child. The deepest internal node has two leaf children.
    """
    assert n_leaves >= 2
    leaves = [Clade(name=f"s{i+1}") for i in range(n_leaves)]
    # Build from the deepest internal upward
    deepest_internal_idx = n_leaves - 2  # internal_{n_leaves-2} has two leaves
    current = Clade(name=f"internal_{deepest_internal_idx}",
                    clades=[leaves[-2], leaves[-1]])
    # Walk outward: each step adds a new internal whose children are (next leaf, current)
    for i in range(deepest_internal_idx - 1, -1, -1):
        leaf = leaves[i]
        current = Clade(name=f"internal_{i}", clades=[leaf, current])
    root = Clade(name="diploid", clades=[current])
    return Tree(root=root, rooted=True)


def test_nni_neighbor_count_n3():
    tree = _make_search_shape_tree(3)
    # n=3 -> 2(n-3) = 0 neighbors
    assert list(nni.nni_neighbors(tree)) == []


def test_nni_neighbor_count_n4():
    tree = _make_search_shape_tree(4)
    # n=4 -> 2(n-3) = 2 neighbors
    assert len(list(nni.nni_neighbors(tree))) == 2


def test_nni_neighbor_count_n5():
    tree = _make_search_shape_tree(5)
    # n=5 -> 2(n-3) = 4 neighbors
    assert len(list(nni.nni_neighbors(tree))) == 4


def test_nni_neighbor_count_n8():
    tree = _make_search_shape_tree(8)
    # n=8 -> 2(n-3) = 10 neighbors
    assert len(list(nni.nni_neighbors(tree))) == 10


def test_nni_neighbors_are_distinct():
    tree = _make_search_shape_tree(6)
    hashes = set()
    input_hash = tree_hash.tree_hash(tree_hash.strip_branch_lengths(tree))
    for neighbor_tree, _ in nni.nni_neighbors(tree):
        h = tree_hash.tree_hash(tree_hash.strip_branch_lengths(neighbor_tree))
        assert h != input_hash, "NNI neighbor identical to input"
        assert h not in hashes, "NNI neighbors not pairwise distinct"
        hashes.add(h)
    # n=6 -> 2(n-3) = 6 neighbors
    assert len(hashes) == 6


def _child_sets_by_name(tree):
    return {c.name: frozenset(child.name for child in c.clades)
            for c in tree.find_clades() if len(c.clades) > 0}


def test_nni_neighbor_changes_exactly_two_internal_child_sets():
    tree = _make_search_shape_tree(5)
    input_csets = _child_sets_by_name(tree)
    for neighbor_tree, _ in nni.nni_neighbors(tree):
        nbr_csets = _child_sets_by_name(neighbor_tree)
        differing = [name for name in input_csets if input_csets[name] != nbr_csets[name]]
        assert len(differing) == 2, \
            f"NNI move should change exactly 2 internal child sets, got {len(differing)}"


def test_nni_preserves_internal_node_names():
    tree = _make_search_shape_tree(5)
    input_names = {c.name for c in tree.find_clades() if c.name is not None}
    for neighbor_tree, _ in nni.nni_neighbors(tree):
        nbr_names = {c.name for c in neighbor_tree.find_clades() if c.name is not None}
        assert input_names == nbr_names


def test_nni_skips_root_incident_edges():
    tree = _make_search_shape_tree(5)
    root = tree.root
    root_children = set(root.clades)
    for _, spr_result in nni.nni_neighbors(tree):
        # u (prune_grandparent) must not be root or a root-child
        assert spr_result.prune_grandparent != root.name
        assert spr_result.prune_grandparent not in {c.name for c in root_children}
