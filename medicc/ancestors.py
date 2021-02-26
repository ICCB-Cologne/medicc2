import numpy as np
import pandas as pd
import Bio
import Bio.Phylo

import fstlib

def reconstruct_ancestors(tree, samples_dict, model, normal_name):
    fsa_dict = samples_dict.copy()
    tree_allele = Bio.Phylo.BaseTree.copy.deepcopy(tree)

    depth_first = []
    for clade in tree_allele.find_clades(order = "preorder"):
        depth_first.append(clade)
    idx = [i for i, item in enumerate(depth_first) if item.name == normal_name]
    depth_first.pop(idx[0])
        
    # up the tree
    for node in reversed(depth_first): 
        if len(node.clades) != 0:
            children = [item for item in node.clades if item.name != normal_name]
            left_name = children[0].name
            right_name = children[1].name

            ## finalise
            #sp = fstlib.align(symm_fst, fsa_dict[left_name], fsa_dict[right_name])
            #fsa_dict[left_name] = fstlib.arcmap(sp.copy().project('input'), map_type='rmweight')
            #fsa_dict[right_name] = fstlib.arcmap(sp.copy().project('output'), map_type='rmweight')

            ## project
            intersection = intersect_clades_detmin(fsa_dict[left_name], fsa_dict[right_name], model, prune=0, detmin_before_intersect=False, detmin_after_intersect=True)
            fsa_dict[node.name] = intersection

    # root node separately
    root_name = depth_first[0].name 
    sp = fstlib.align(model, fsa_dict[normal_name], fsa_dict[root_name])
    fsa_dict[root_name] = fstlib.arcmap(sp.copy().project('output'), map_type='rmweight')

    # down the tree
    for node in depth_first:
        if len(node.clades) != 0:
            children = [q for q in node.clades if len(q.clades) != 0]
            for child in children:
                sp = fstlib.align(model, fsa_dict[node.name], fsa_dict[child.name])
                fsa_dict[child.name] = fstlib.arcmap(sp.copy().project('output'), map_type='rmweight')

    return fsa_dict

def intersect_clades_detmin(left, right, model, prune=None, detmin_before_intersect=True, detmin_after_intersect=True):
    L = fstlib.compose(model, left.arcsort('ilabel')).project('input')
    R = fstlib.compose(model, right.arcsort('ilabel')).project('input')
    if detmin_before_intersect:
        L = fstlib.determinize(L).minimize()
        R = fstlib.determinize(R).minimize()
    intersection = fstlib.intersect(L.arcsort('olabel'),R)
    if prune is not None:
        pruned = fstlib.prune(intersection, weight=prune)
    else:
        pruned = intersection
    if detmin_after_intersect:
        pruned = fstlib.determinize(pruned).minimize()
    return pruned
