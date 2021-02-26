'''
Created on 17 May 2010

@author: rfs

FST algorithms
'''

from fstlib.core import *
import sys
import numpy as np
import scipy as sp
import fstlib.core
import array
import logging
import fstlib.tools
import functools

log = logging.getLogger(__name__)

class FstAlgoError(Exception):
	pass

class ConfusionNetworkFromLattice(object):
	""" Expects lattice in the form of a FST object. 
	Requires lattice to be topologically sorted
	(arcs are always directed from lower to higher state ids). It might require (not sure)
   the initial state to have index 0. """
	def __init__(self, lattice):
		self.lattice = lattice
		self.state2set = None ## state set ids for each state in lattice
		self.statesets = None ## list of lists of state ids contained in the state sets
		self.edgesets = None ## new edges between state sets
		
	def run(self, collapse_arcs=True):
		log.info('Constructing confusion network...')
		lattice = self.lattice
		outfst = fstlib.core.Fst(arc_type = self.lattice.arc_type())
		outfst.set_input_symbols(self.lattice.input_symbols())
		outfst.set_output_symbols(self.lattice.output_symbols())
		timeslices = fstlib.shortestdistance(fstlib.arcmap(fstlib.arcmap(fstlib.arcmap(self.lattice, map_type='rmweight'), map_type='times', weight=1), map_type='to_std'))
		statemap = np.array([int(float(t)) for t in timeslices])
		newstates = np.array(outfst.add_states(np.unique(statemap).shape[0]))
		outfst.set_start(newstates[0])
		for s in newstates:
			weights = [self.lattice.final(t) for t in np.flatnonzero(statemap==s) if self.lattice.is_final(t)]
			if len(weights)>0: ## newstate is final if any of the merged states were final (have non-inf weights)
				total_weight = functools.reduce(fstlib.times, weights)
				outfst.set_final(s, total_weight)

		for state in self.lattice.states():
			newstate = statemap[state]
			for arc in self.lattice.arcs(state):
				newnextstate = statemap[arc.nextstate]
				outfst.add_arc(newstate, (arc.ilabel, arc.olabel, arc.weight, newnextstate))

		if collapse_arcs:
			outfst=fstlib.determinize(outfst)
		return outfst

class PosteriorDecoding(object):
	""" expects pFSA in binary (i.e. openfst) format """
	def __init__(self, ifst):
		self.ifst = ifst

	def run(self):
		alpha = fstlib.shortestdistance(self.ifst)
		beta = fstlib.shortestdistance(self.ifst, reverse=True)
		
		outfst = fstlib.core.Fst(arc_type=self.ifst.arc_type())
		outfst.set_input_symbols(self.ifst.input_symbols())
		outfst.set_output_symbols(self.ifst.output_symbols())
		
		for state in self.ifst.states():
			newstate = outfst.add_state()
			if self.ifst.is_final(state):
				outfst.set_final(newstate, self.ifst.final(state))
			if state != newstate:
				raise FSTError('State mapping error') 
				## this must not happen; I'm not sure it can currently due to the way the states are handled
		initial_state = self.ifst.start()
		outfst.set_start(initial_state)
		for state in self.ifst.states():
			for arc in self.ifst.arcs(state):
				posterior = float(alpha[state]) + float(arc.weight) + float(beta[arc.nextstate]) - float(beta[initial_state])
				outfst.add_arc(state, (arc.ilabel, arc.olabel, posterior, arc.nextstate))
		return outfst

class AbstractCountDP(object):
	def __init__(self, fst, trans_cache=None, state_depth_cache=None):
		""" constructor """
		self.fst = fst
		
		## if transition cache was given use it, otherwise create it
		if trans_cache == None:
			self._build_trans_cache()
		else:
			self._trans_cache = trans_cache
		
		## if state depth cache was given use it, otherwise create it
		if state_depth_cache == None:
			self._build_state_depth_cache()
		else:
			self._state_depth_cache = state_depth_cache
			
		assert(fst != None)

	def _build_trans_cache(self):
		#self._trans_cache = numpy.empty( (len(self.fst.states), len(self.fst.states)), dtype="object")
		self._trans_cache = [0] * len(self.fst.states)
		self._trans_cache = [[i] * len(self.fst.states) for i in self._trans_cache]
		
		for i in range(0, len(self.fst.states)):
			from_transs = [t for t in self.fst.arcs if t.state_from.id == i]
			for j in range(0, len(self.fst.states)):
				transs = len([t for t in from_transs if t.state_to.id == j])
				self._trans_cache[i][j] = transs
				
	def _build_state_depth_cache(self):
		""" tries to speed up the algo """
		self._state_depth_cache=list()
		depth = 0
		self._build_state_depth_cache_recurse(self.fst.states[0], depth)
		
	def _build_state_depth_cache_recurse(self, state, depth):
		""" recursion function of the pre inspection """
		
		if len(self._state_depth_cache) <= depth:
			self._state_depth_cache.append(set())
		self._state_depth_cache[depth].add(state.id)
		
		targets = list(set([t.state_to for t in state.get_outgoing_transitions()]))
		for s in targets:
			self._build_state_depth_cache_recurse(s, depth + 1)

	def max_seq_len(self):
		return len(self._state_depth_cache)-1

	def get_transition_cache(self):
		return self._trans_cache
	
	def get_state_depth_cache(self):
		return self._state_depth_cache
				
		
class ForwardCountDP(AbstractCountDP):
	""" Dynamic programming forward algorithm.
	This is not in log space and very generic. """
	
	def __init__(self, fst, trans_cache=None, state_depth_cache=None):
		""" constructor """
		super(ForwardCountDP, self).__init__(fst, trans_cache, state_depth_cache)
		self.alpha = None
		self._run()
		
	def _run(self):
		""" main algorithm """
		
		nstates = len(self.fst.states)
		## initialize the forward matrix
		# the sequence index starts with 1, 0 means "no part of the sequence seen so far"
		##self.alpha = numpy.zeros((nstates, 1))
		self.alpha = list()
		for i in range(0, self.max_seq_len() + 1):
			self.alpha.append(array.array("L", [0] * nstates))
		
		# we assume the starting state is always state 0
		self.alpha[0][0] = 1 

		for seqpos in range(1, self.max_seq_len() + 1):
			for stateid_to in self._state_depth_cache[seqpos]:
				value = 0
				for stateid_from in self._state_depth_cache[seqpos-1]:
					value += self._trans_cache[stateid_from][stateid_to] * self.alpha[seqpos-1][stateid_from]
				self.alpha[seqpos][stateid_to] = value
	
class BackwardCountDP(AbstractCountDP):
	def __init__(self, fst, trans_cache=None, state_depth_cache=None):
		""" constructor """
		super(BackwardCountDP, self).__init__(fst, trans_cache, state_depth_cache)
		self.beta = None
		self._run()

	
	def _run(self):
		""" main algorithm """
		
		## initialize the forward matrix
		# the sequence index starts with 1, 0 means "no part of the sequence seen so far"
		nstates = len(self.fst.states)
		self.beta = list()
		for stateid_from in range(0, self.max_seq_len() + 1):
			self.beta.append(array.array("L", [0] * nstates))
		
		# all final states are assigned 1
		for id in self._state_depth_cache[self.max_seq_len()]:
			self.beta[self.max_seq_len()][id] = 1 

		for seqpos in range(self.max_seq_len()-1, -1, -1):
			for stateid_from in self._state_depth_cache[seqpos]:
				value = 0
				for stateid_to in self._state_depth_cache[seqpos + 1]:
					value += self.beta[seqpos+1][stateid_to] * self._trans_cache[stateid_from][stateid_to]  
				self.beta[seqpos][stateid_from] = value

class PosteriorDecodingCountDP(AbstractCountDP):
	""" Dynamic programming algorithm for posterior decoding of counts """
	
	def __init__(self, fst):
		
		""" constructor """
		self.forward_algo = ForwardCountDP(fst)
		self.backward_algo = BackwardCountDP(fst,  
								   self.forward_algo.get_transition_cache(),
								   self.forward_algo.get_state_depth_cache())
		super(PosteriorDecodingCountDP, self).__init__(fst, 
												  self.forward_algo.get_transition_cache(),
												  self.forward_algo.get_state_depth_cache())        
		self.posterior = [dict() for i in range(0, self.max_seq_len())]
		self._run()
		
	def _run(self):
		""" main algorithm """
		nstates = len(self.fst.states)
		
		for pos in range(0, self.max_seq_len()):
			for state_from in self._state_depth_cache[pos]:
				all_transitions = self.fst.states[state_from].get_outgoing_transitions()
				for state_to in self._state_depth_cache[pos+1]:
					transitions = [t for t in all_transitions if t.state_to.id == state_to]
					fw = self.forward_algo.alpha[pos][state_from]
					bw = self.backward_algo.beta[pos+1][state_to]
					for trans in transitions:
						try:
							value = self.posterior[pos][trans.isymbol]
						except KeyError:
							value = 0
						value += fw * numpy.exp(-trans.weight) * bw
						self.posterior[pos][trans.isymbol] = value



class AbstractDP(object):
	def __init__(self, fst, trans_cache=None):
		""" constructor """
		self.fst = fst
		if trans_cache == None:
			self._build_cache()
		else:
			self._trans_cache = trans_cache
		assert(fst != None)

	def _get_transitions(self, state_from, state_to):
		""" returns all transitions between states with caching.
		Expects numerical state ids. """
		
		return self._trans_cache[state_from][state_to]
	
	def _build_cache(self):
		#self._trans_cache = numpy.empty( (len(self.fst.states), len(self.fst.states)), dtype="object")
		self._trans_cache = [None] * len(self.fst.states)
		self._trans_cache = [[i] * len(self.fst.states) for i in self._trans_cache]
		
		for i in range(0, len(self.fst.states)):
			from_transs = [t for t in self.fst.arcs if t.state_from.id == i]
			for j in range(0, len(self.fst.states)):
				transs = array.array('f', [t.weight for t in from_transs if t.state_to.id == j])
				self._trans_cache[i][j] = transs
				
				
		
class PosteriorDecodingDP(AbstractDP):
	""" Dynamic programming algorithm for posterior decoding """
	
	def __init__(self, fst):
		""" constructor """
		super(PosteriorDecodingDP, self).__init__(fst)        
		self.forward_algo = ForwardDP(self.fst)
		self.backward_algo = BackwardDP(self.fst, 
								   self.forward_algo.get_sequence_length(), 
								   self.forward_algo.get_transition_cache())
		self._trans_cache = self.forward_algo.get_transition_cache()
		self.posterior = [dict() for i in range(0, self.forward_algo.get_sequence_length())]
		self._run()
		
	def _run(self):
		""" main algorithm """
		
		
		for pos in range(0, self.forward_algo.get_sequence_length()):
			for state_from in range(0, len(self.fst.states)):
				for state_to in range(0, len(self.fst.states)):
					transitions = self._get_transitions(state_from, state_to)
					fw = self.forward_algo.alpha[state_from, pos]
					bw = self.backward_algo.beta[state_to, pos+1]
					for trans in transitions:
						try:
							value = self.posterior[pos][trans.isymbol]
						except KeyError:
							value = 0
						value += fw * numpy.exp(-trans.weight) * bw
						self.posterior[pos][trans.isymbol] = value
						
		
class ForwardDP(AbstractDP):
	""" Dynamic programming forward algorithm.
	This is not in log space and very generic. """
	
	def __init__(self, fst, trans_cache=None):
		""" constructor """
		super(ForwardDP, self).__init__(fst, trans_cache)
		self.alpha = None        
		self._run()
		
	def _run(self):
		""" main algorithm """
		
		nstates = len(self.fst.states)
		## initialize the forward matrix
		# the sequence index starts with 1, 0 means "no part of the sequence seen so far"
		##self.alpha = numpy.zeros((nstates, 1))
		self.alpha = list()
		self.alpha.append(array.array("f", [0] * nstates))
		
		# we assume the starting state is always state 0
		self.alpha[0][0] = 1 

		j=0
		while j<len(self.alpha):
			j+=1
			column = array.array("f", [0] * nstates)
			any_weight = False        
			for i in range(0, nstates):
				value = 0
				for k in range(0, nstates):
					transitions = self._trans_cache[k][i]
					for t in transitions:
						value += self.alpha[j-1][k] * math.exp(-t)  
				column[i] = value
				if value > 0.0: any_weight = True
			if any_weight:
				self.alpha.append(column)
	
	def get_sequence_length(self):
		return len(self.alpha)
	
	def get_transition_cache(self):
		return self._trans_cache
				
class BackwardDP(AbstractDP):
	""" Dynamic programming backward algorithm.
	This is not in log space and very generic. """
	
	def __init__(self, fst, seqlen, trans_cache=None):
		""" constructor """
		super(BackwardDP, self).__init__(fst, trans_cache)
		self.seqlen = seqlen
		self.beta = None
		self._run()
		
	def _run(self):
		""" main algorithm """
		
		## initialize the forward matrix
		# the sequence index starts with 1, 0 means "no part of the sequence seen so far"
		self.beta = numpy.zeros((len(self.fst.states), self.seqlen + 1))
		
		# all final states are assigned 1
		for state in self.fst.states:
			if state.is_final:
				self.beta[state.id, self.seqlen] = 1 

		for j in range(self.seqlen-1, -1, -1):
			for i in range(0, len(self.fst.states)):
				value = 0
				for k in range(0, len(self.fst.states)):
					transitions = self._get_transitions(i, k)
					for t in transitions:
						value += self.beta[k, j+1] * numpy.exp(-t)  
				self.beta[i,j] = value
		
		

class PathDepthFirstSearch(object):
	""" recursively walks through the FST and prints the paths.
	Make sure the FST does not contain loops, otherwise the function will not terminate. """
	
	def __init__(self, fst):
		'''
		Constructor
		'''
		self.fst = fst
	
	def get_paths(self):
		""" generator that returns paths """

		if self.fst.num_states()==0:
			return

		state = self.fst.start()
		next_arc = {}
		stack = []
		path = fstlib.core.Path()

		while state != None:
			arcs = [a for a in self.fst.arcs(state)]
			i = 0
			try:
				i = next_arc[state]
			except KeyError:
				pass
			
			if self.fst.num_arcs(state) > i: ## move forward
				arc = arcs[i]  ## get forward transition
				path.append(arc) ## save transition
				stack.append(state) ## put state on stack
				next_arc[state] = i + 1 ## increase trans counter by 1
				state = arc.nextstate ## move to next state
			else: ## move back
				if self.fst.is_final(state): ## end node
					path.finalWeight = self.fst.final(state)
					yield path.copy()

				next_arc[state] = 0    
				try:
					state = stack.pop()
					path.pop()
				except IndexError:
					state = None
				
				
							
class FeatureDepthFirstSearch(object):
	""" recursively walks through the FST and prints the paths.
	Make sure the FST does not contain loops, otherwise the function will not terminate. """
	
	def __init__(self, fst):
		'''
		Constructor
		'''
		self.fst = fst
		self.paths = []
		self.path_limit = None
		assert(fst != None)
	
	def run(self):
		""" starts the recursion """
		self.paths = []

		self._iteration()
	
	def _iteration(self):
		""" interative walk through the DAG """
		
		state = self.fst.states[0]
		suffix_cache = {}
		next_trans = {}
		stack = []
		counter = 0
		max_depth = 0
		while state != None:
			transitions_out = state.get_outgoing_transitions()
			i = 0
			try:
				i = next_trans[state]
			except KeyError:
				pass
			
			if len(transitions_out) > i: ## move forward
				max_depth += 1
				trans = transitions_out[i]  ## get forward transition
				stack.append(state) ## put state on stack
				next_trans[state] = i + 1 ## increase trans counter by 1
				state = trans.state_to ## move to next state                
				
			else: ## generate suffixes and jump back
				## do the jump back
				#print counter
				#counter += 1
				
				previous_state = None
				try:
					previous_state = stack.pop()
				except IndexError:
					pass
				
				if previous_state != None:
					prev_trans_id = next_trans[previous_state] - 1
					prev_trans = previous_state.get_outgoing_transitions()[prev_trans_id]
					suffixes_current = [[max_depth - 0, -1, -1]]
					try:
						suffixes_current = suffix_cache[state]
						#print suffixes_current
					except KeyError:
						pass
					#print "word: %s" % str(word)
					#print "concatenating..."
					#new_suffixes = [[word] + s for s in suffixes]
					#print "new suffixes:"
					#print [[str(t) for t in s] for s in new_suffixes]
						
					suffixes_prev = None
					try:
						suffixes_prev = suffix_cache[previous_state]
					except KeyError:
						suffixes_prev = []
						suffix_cache[previous_state] = suffixes_prev

					for suf in suffixes_current:
						new_suffix = [suf[0] - 1, suf[1], suf[2]]
						if prev_trans.osymbol == "X":
							assert( not(new_suffix[1] > -1 and new_suffix[2] > -1) )
							if new_suffix[2] == -1:
								new_suffix[2] = suf[0] - 2
							elif new_suffix[1] == -1:
								new_suffix[1] = suf[0]
								
						suffixes_prev.append(new_suffix)
					
					#print suffixes_prev
						
						#sc = sc + new_suffixes
						#print "sc:"
						#print [[str(t) for t in s] for s in sc]
					
						
				
				state = previous_state
		
		self.paths = suffix_cache[self.fst.states[0]]

def posterior_decoding(ifst, collapse_arcs=True):
	pd = PosteriorDecoding(ifst).run()
	conf = ConfusionNetworkFromLattice(pd).run(collapse_arcs)
	return conf

def map_log_to_real(x):
	return np.exp(-float(x))

def map_real_to_log(x):
	return -np.log(float(x))