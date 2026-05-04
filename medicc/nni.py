"""NNI (Nearest Neighbor Interchange) move enumeration for MEDICC2.

Mirrors the shape and conventions of medicc/spr.py. Yields all NNI neighbors of
a rooted MEDICC tree along with a synthetic SPRResult for each, so they can be
fed into reconstruct_ancestors_incremental without any new infrastructure.

The NNI move at internal edge (u, v) — where u is the parent of v and v is itself
internal — picks the "uncle" A (the other child of u) and swaps it with one of v's
children (B or C). Two distinct moves per qualifying edge.

Edges incident to the root are excluded, matching the SPR convention in spr.py.
The MEDICC search-shape tree has the diploid as root with a single structural
child (the MRCA), so root and root-children are skipped.
"""
import copy

from medicc.spr import SPRResult


def _is_root_or_root_child(clade, root, root_children):
    return clade is root or clade in root_children


def _apply_nni_swap(tree, u_name, v_name, A_name, swap_target_name):
    """Deepcopy `tree` and apply an NNI swap at edge (u, v), exchanging A and swap_target.

    Pre: u is parent of v; A is u's other child; swap_target is one of v's children.
    Post: u's children are {swap_target, v}; v's children are {A, the_other_v_child}.

    Returns the new tree. No new internal nodes are created; all names preserved.
    """
    new_tree = copy.deepcopy(tree)
    # Locate u, v, A, swap_target by name in the copy
    name_to_clade = {c.name: c for c in new_tree.find_clades() if c.name is not None}
    u = name_to_clade[u_name]
    v = name_to_clade[v_name]
    A = name_to_clade[A_name]
    swap_target = name_to_clade[swap_target_name]

    # Rewire: remove A from u, append swap_target to u
    u.clades = [c for c in u.clades if c is not A]
    u.clades.append(swap_target)
    # Rewire: remove swap_target from v, append A to v
    v.clades = [c for c in v.clades if c is not swap_target]
    v.clades.append(A)

    return new_tree


def nni_neighbors(tree):
    """Yield (new_tree, synthetic_SPRResult) for every NNI neighbor of `tree`.

    Iterates over internal edges (u, v) skipping edges incident to the root,
    and for each qualifying edge yields the 2 uncle-swap variants.

    The synthetic SPRResult uses prune_grandparent=u.name and regraft_parent=v.name
    so reconstruct_ancestors_incremental._get_dirty_nodes_up walks
    {v, u, parent(u), ..., root} — exactly the dirty set for an NNI move.
    """
    root = tree.root
    root_children = set(root.clades)

    for u in tree.find_clades():
        if _is_root_or_root_child(u, root, root_children):
            continue
        # u must be an internal node with at least 2 children
        if len(u.clades) < 2:
            continue
        for v in u.clades:
            # v must itself be internal (have children to swap with A)
            if len(v.clades) < 2:
                continue
            siblings_of_v = [c for c in u.clades if c is not v]
            if not siblings_of_v:
                continue
            A = siblings_of_v[0]  # binary tree: exactly one sibling
            for swap_target in list(v.clades):
                new_tree = _apply_nni_swap(
                    tree=tree,
                    u_name=u.name,
                    v_name=v.name,
                    A_name=A.name,
                    swap_target_name=swap_target.name,
                )
                synthetic = SPRResult(
                    tree=new_tree,
                    pruned_subtree_root=swap_target.name,  # cosmetic only
                    prune_grandparent=u.name,
                    regraft_parent=v.name,
                )
                yield new_tree, synthetic
