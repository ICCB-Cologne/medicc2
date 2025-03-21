import copy

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree 
import random 

def get_random_candidate_clades(tree, prune = True, parent_name = None):
    # Get all clades
    all_clades = list(tree.find_clades())
    
    # Get the root and its direct children
    root = tree.root
    root_children = set(root.clades)
    
    # Filter out the root and its direct children
    if prune:
        # Find eligible clades for pruning
        eligible_clades = [clade for clade in all_clades 
                      if clade != root and clade not in root_children]
    else:
        # Find eligible clades for regrafting
        eligible_clades = [clade for clade in all_clades
                      if clade != root and clade.name != parent_name]
    
    # Check if we have any eligible clades
    if not eligible_clades:
        return None
    
    # Return a random clade from the eligible ones
    return random.choice(eligible_clades)

def find_sibling(tree, prune_clade):
    # Get the parent of the prune clade
    parent = tree.get_path(prune_clade)[-2]
    
    # Get the siblings of the prune clade
    siblings = [clade for clade in parent.clades if clade != prune_clade]
    
    return siblings

def get_parent_and_grandparent(tree, clade):
    path = tree.get_path(clade)
    path = [tree.root] + path
    parent = path[-2]
    grandparent = path[-3]
    return parent, grandparent

def connect_grandparent_to_sibling(tree, prune_clade):
    parent, grandparent = get_parent_and_grandparent(tree, prune_clade)
    siblings = find_sibling(tree, prune_clade)

    # Remove the parent from the grandparent's children
    grandparent.clades.remove(parent)

    # Connect the grandparent to the sibling
    grandparent.clades.extend(siblings)

    return tree 

def prune_tree(tree):
    prune_clade = get_random_candidate_clades(tree, prune = True)
    parent, grandparent = get_parent_and_grandparent(tree, prune_clade)
    tree = connect_grandparent_to_sibling(tree, prune_clade)

    return tree, prune_clade, parent.name 

def regraft_tree(pruned_tree, subtree, parent_name):
    regraft_clade = get_random_candidate_clades(pruned_tree, prune = False, parent_name=parent_name)

    # Find the parent of the regraft clade
    regraft_clade_path = [pruned_tree.root] + pruned_tree.get_path(regraft_clade)
    parent_reg = regraft_clade_path[-2]

    # detach the regraft clade to its parent
    parent_reg.clades.remove(regraft_clade)

    # Create a new clade with parent name and regraft clade and subtree as its children
    new_clade = Clade(name=parent_name, clades=[regraft_clade, subtree])

    # Attach this new clade to the parent of the regraft clade
    parent_reg.clades.append(new_clade)

    return pruned_tree

def spr_move(tree):
    # NOTE: This function modifies the tree in place!!!
    tree_copy = copy.deepcopy(tree)
    # Prune a subtree
    pruned_tree, subtree, parent_name = prune_tree(tree_copy)

    # Regraft the subtree
    tree_SPR = regraft_tree(pruned_tree, subtree, parent_name)

    return tree_SPR 