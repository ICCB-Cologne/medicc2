'''
Created on 09 Dec 2019

@author: mp 
'''

import logging

import Bio
import Bio.Nexus
import Bio.Phylo
import numpy as np

SS_TAXON_SEPERATOR = ","
CHR_SEPARATOR = "X"

logger = logging.getLogger(__name__)

class NeighbourJoining(object):
    """ computes a neighbour joining tree from a distance matrix """
    
    def __init__(self, distance_matrix, taxa):
        """ constructor:
        requires distance matrix as numpy array and list of taxa 
        names in the same order as in the matrix """
        
        self.distmat = distance_matrix
        self.taxa = taxa
        self.tree = None
        self.internal_nodes = []
        self.log = logger
        self._run_nj()
    
    def _run_nj(self):
        """ main algorithm routine """
        
        self.log.info("Neighbor joining:")
        name_clade_map = dict()
        internal_node_list = list()
        sequence_files = list(self.taxa) ## create mutable copy
        distances = self.distmat.copy()
        
        while len(sequence_files) > 2: ## do while still more than two files present
            q_matrix = self._nj_create_q_matrix(distances)
            (left, right) = self._nj_nearest_neighbours(q_matrix)
            left_name = sequence_files[left]
            right_name = sequence_files[right]
            self.log.debug("Choosing %s with distance %f" % ((left_name, right_name), distances[left,right]))
            new_distances = self._nj_get_distances_to_new_node((left, right), distances)
            distances = self._nj_forge_new_distance_matrix((left, right), distances, new_distances)
            ## ancestor_filename = "(%s-%s)" % (left_name, right_name)
            internal_node_count = len(internal_node_list) + 1
            ancestor_filename = "internal_%d" % internal_node_count
        
            clade_ancestor = Bio.Phylo.PhyloXML.Clade(branch_length=0)
            clade_ancestor.name = ancestor_filename
            
            try:
                clade_left = name_clade_map[left_name]
            except KeyError:
                clade_left = Bio.Phylo.PhyloXML.Clade()
                clade_left.name = left_name
                name_clade_map[left_name] = clade_left
            finally:
                clade_left.branch_length = new_distances[left]
            
            try:
                clade_right = name_clade_map[right_name]
            except KeyError:
                clade_right = Bio.Phylo.PhyloXML.Clade()
                clade_right.name = right_name
                name_clade_map[right_name] = clade_right        
            finally:
                clade_right.branch_length = new_distances[right]
            
            clade_ancestor.clades = [clade_left, clade_right]
            name_clade_map[ancestor_filename] = clade_ancestor
            internal_node_list.append(clade_ancestor)
            
            sequence_files.remove(left_name)
            sequence_files.remove(right_name)
            sequence_files.append(ancestor_filename)
            
            self.log.debug("Remaining number of sequences: %d" % len(sequence_files))
        
        ## here we have exactly 2 sequences left
        last_index = (sequence_files.index(ancestor_filename) + 1) % 2
        last_name = sequence_files[last_index]
        
        try:
            clade_last = name_clade_map[last_name]
        except KeyError:
            clade_last = Bio.Phylo.PhyloXML.Clade()
            clade_last.name = last_name
            name_clade_map[last_name] = clade_last        
        finally:
            clade_last.branch_length = new_distances[last_index]
        
        clade_ancestor.clades.append(clade_last) ## this one has 3 descendants
        
        self.tree = Bio.Phylo.PhyloXML.Phylogeny(root = clade_ancestor)
        self.internal_nodes = internal_node_list
    
    def _nj_create_q_matrix(self, D):
        """ computes the values of a neighbour-joining Q matrix out of the distance matrix D """
        Q = (D.shape[0]-2) * D - np.sum(D, axis=0, keepdims=True) - np.sum(D, axis=1, keepdims=True)
        
        return Q
    
    def _nj_nearest_neighbours(self, D):
        """ returns the nearest neighbours from the distance matrix D """
        np.fill_diagonal(D, np.inf)
        (min_i, min_j) = np.unravel_index(np.argmin(D), D.shape)
    
        return (min_i, min_j)
    
    def _nj_get_distances_to_new_node(self, left_right_tuple, D):
        (left, right) = left_right_tuple
        r = D.shape[0]
        distances = np.zeros((r))
        
        ## first the dist from the joined nodes to the new node
        d_left = 0.5 * D[left, right] + 1.0 / (2 * (r - 2)) * (np.sum(D[left, ]) - np.sum(D[right, ]))
        d_right = 0.5 * D[left, right] + 1.0 / (2 * (r - 2)) * (np.sum(D[right, ]) - np.sum(D[left, ]))

        ## now the dist from all other nodes to the new node
        distances = 0.5 * (D[left, :] - d_left) + 0.5 * (D[right, :] - d_right)

        distances[left] = d_left
        distances[right] = d_right

        return distances 
    
    def _nj_forge_new_distance_matrix(self, left_right_tuple, D, new_distances):
        (left, right) = left_right_tuple 
        r = D.shape[0]
        idx = np.repeat(True, r)
        idx[left] = False
        idx[right] = False
        
        E = D[idx][...,idx].copy()
    
        F = np.append(E, [new_distances[idx]], 0)
        G = np.append(F.transpose(), [np.append(new_distances[idx], [0])], 0).transpose()
        
        return G



