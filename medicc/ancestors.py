import Bio
import Bio.Phylo
import fstlib
import numpy as np
import pandas as pd


def reconstruct_ancestors(tree, samples_dict, fst, normal_name):
    fsa_dict = samples_dict.copy()
    tree = Bio.Phylo.BaseTree.copy.deepcopy(tree)

    clade_list = [clade for clade in tree.find_clades(order="preorder") if clade.name != normal_name]
    # up the tree (leaf to root)
    for node in reversed(clade_list): 
        if len(node.clades) != 0:
            children = [item for item in node.clades if item.name != normal_name]
            left_name = children[0].name
            right_name = children[1].name

            ## project
            intersection = intersect_clades_detmin(fsa_dict[left_name], fsa_dict[right_name], fst, 
                                                   prune_weight=0, detmin_before_intersect=False, detmin_after_intersect=True)
            fsa_dict[node.name] = intersection

    # root node is calculated separately w.r.t. normal node
    root_name = clade_list[0].name 
    sp = fstlib.align(fst, fsa_dict[normal_name], fsa_dict[root_name])
    fsa_dict[root_name] = fstlib.arcmap(sp.copy().project('output'), map_type='rmweight')

    # down the tree (root to leaf)
    for node in clade_list:
        if len(node.clades) != 0:
            children = [q for q in node.clades if len(q.clades) != 0]
            for child in children:
                sp = fstlib.align(fst, fsa_dict[node.name], fsa_dict[child.name])
                fsa_dict[child.name] = fstlib.arcmap(sp.copy().project('output'), map_type='rmweight')

    return fsa_dict

def intersect_clades_detmin(left, right, fst, prune_weight=None, detmin_before_intersect=True, detmin_after_intersect=True):
    L = fstlib.compose(fst, left.arcsort('ilabel')).project('input')
    R = fstlib.compose(fst, right.arcsort('ilabel')).project('input')
    if detmin_before_intersect:
        L = fstlib.determinize(L).minimize()
        R = fstlib.determinize(R).minimize()
    intersection = fstlib.intersect(L.arcsort('olabel'), R)
    # For prune_weight=0, deletes all paths but the shortest one
    if prune_weight is not None:
        pruned = fstlib.prune(intersection, weight=prune_weight)
    else:
        pruned = intersection
    if detmin_after_intersect:
        pruned = fstlib.determinize(pruned).minimize()
    return pruned
