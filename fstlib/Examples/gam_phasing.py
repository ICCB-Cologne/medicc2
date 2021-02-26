import numpy as np
import imp
import os
import sys
import random
import pandas as pd
import logging
import matplotlib as mpl
import matplotlib.pyplot as pl
import time
import functools
import seaborn as sns
sys.path.append('../..')
import fstlib as fst


logging.basicConfig(level=logging.INFO)

def main():
	semiring = 'log'
	symbols = ['a','b','.','_','=','x']
	symbol_table = fst.SymbolTable()
	eps = '<eps>'
	symbol_table.add_symbol(eps)
	for s in symbols:
		symbol_table.add_symbol(s)

	symbol_file = 'gam.symbols'
	symbol_table.write_text(symbol_file)

	## some general FSTs used in this example
	## fst for transforming {a,b,.} sequence into {x,=}
	stayflip = fst.Fst(arc_type=semiring)
	stayflip.set_input_symbols(symbol_table)
	stayflip.set_output_symbols(symbol_table)
	stayflip.add_states(4)
	stayflip.set_start(0)
	stayflip.set_final(1,0)
	stayflip.set_final(2,0)
	stayflip.set_final(3,0)
	stayflip.add_arc(0, ('a', eps, 0, 1))
	stayflip.add_arc(0, ('b', eps, 0, 2))
	stayflip.add_arc(0, ('.', eps, 0, 3))
	stayflip.add_arc(1, ('a', '=', 0, 1))
	stayflip.add_arc(1, ('b', 'x', 0, 2))
	stayflip.add_arc(1, ('.', '_', 0, 3))
	stayflip.add_arc(2, ('b', '=', 0, 2))
	stayflip.add_arc(2, ('a', 'x', 0, 1))
	stayflip.add_arc(2, ('.', '_', 0, 3))
	stayflip.add_arc(3, ('.', '_', 0, 3))
	stayflip.add_arc(3, ('a', '_', 0, 1))
	stayflip.add_arc(3, ('b', '_', 0, 2))
	stayflip.verify()
	print(stayflip.to_real())

	## 'flip distance' between haplotypes 
	## (minimum number of a<->b flips of arbitrary length needed to transform one sequence into another)
	flipdist = fst.Fst(arc_type='standard')
	flipdist.set_input_symbols(symbol_table)
	flipdist.set_output_symbols(symbol_table)
	flipdist.add_states(2)
	flipdist.set_start(0)
	flipdist.set_final(0,0)
	flipdist.set_final(1,0)
	flipdist.add_arc(0, ('a', 'a', 0, 0))
	flipdist.add_arc(0, ('b', 'b', 0, 0))
	flipdist.add_arc(0, ('a', 'b', 1, 1))
	flipdist.add_arc(0, ('b', 'a', 1, 1))
	flipdist.add_arc(1, ('a', 'b', 0, 1))
	flipdist.add_arc(1, ('b', 'a', 0, 1))
	flipdist.add_arc(1, ('a', 'a', 0, 0))
	flipdist.add_arc(1, ('b', 'b', 0, 0))
	assert(flipdist.verify())

	## filter for inference (only allows paths starting with 'a')
	filterfst = fst.Fst(arc_type=semiring)
	filterfst.set_input_symbols(symbol_table)
	filterfst.set_output_symbols(symbol_table)
	filterfst.add_states(2)
	filterfst.set_start(0)
	filterfst.set_final(1,0)
	filterfst.add_arc(0, fst.Arc(symbol_table.find('a'), symbol_table.find('a'), fst.Weight.One(stayflip.weight_type()), 1))
	for (key, value) in symbol_table:
		if key!=0:
			filterfst.add_arc(1, fst.Arc(key, key, fst.Weight.One(stayflip.weight_type()), 1))
	filterfst.verify()

	#%% functions
	def create_gam_fst(p_switch, p_extend, p_drop, p_drop_extend, semiring='log'):
		""" creates the full gam and simple gam fsts from a set of drop and switch parameters """
		haplofst = fst.Fst(arc_type=semiring)
		haplofst.set_input_symbols(symbol_table)
		haplofst.set_output_symbols(symbol_table)
		haplofst.add_states(3)
		haplofst.set_start(0)
		haplofst.set_final(1,0)
		haplofst.set_final(2,0)
		haplofst.add_arc(0, (eps, eps, -np.log(0.5), 1))
		haplofst.add_arc(0, (eps, eps, -np.log(0.5), 2))
		haplofst.add_arc(1, ('a', 'a', -np.log(p_extend), 1))
		haplofst.add_arc(1, ('b', 'b', -np.log(p_extend), 1))
		haplofst.add_arc(1, ('a', 'b', -np.log(p_switch), 2))
		haplofst.add_arc(1, ('b', 'a', -np.log(p_switch), 2))
		haplofst.add_arc(2, ('a', 'b', -np.log(p_extend), 2))
		haplofst.add_arc(2, ('b', 'a', -np.log(p_extend), 2))
		haplofst.add_arc(2, ('a', 'a', -np.log(p_switch), 1))
		haplofst.add_arc(2, ('b', 'b', -np.log(p_switch), 1))
		fst.tools.normalize_alphabet(haplofst, inplace=True)
		haplofst.arcsort('olabel')
		assert(haplofst.verify())
		#print(haplofst.to_real())

		dropoutfst = fst.Fst(arc_type=semiring)
		dropoutfst.set_input_symbols(symbol_table)
		dropoutfst.set_output_symbols(symbol_table)
		dropoutfst.add_states(3)
		dropoutfst.set_start(0)
		dropoutfst.set_final(1,0)
		dropoutfst.set_final(2,0)
		dropoutfst.add_arc(0, (eps, eps, -np.log(0.5), 1))
		dropoutfst.add_arc(0, (eps, eps, -np.log(0.5), 2))
		dropoutfst.add_arc(1, ('a', 'a', -np.log(1-p_drop), 1))
		dropoutfst.add_arc(1, ('b', 'b', -np.log(1-p_drop), 1))
		dropoutfst.add_arc(1, ('a', '.', -np.log(p_drop), 2))
		dropoutfst.add_arc(1, ('b', '.', -np.log(p_drop), 2))
		dropoutfst.add_arc(2, ('a', '.', -np.log(p_drop_extend), 2))
		dropoutfst.add_arc(2, ('b', '.', -np.log(p_drop_extend), 2))
		dropoutfst.add_arc(2, ('a', 'a', -np.log(1-p_drop_extend), 1))
		dropoutfst.add_arc(2, ('b', 'b', -np.log(1-p_drop_extend), 1))
		fst.tools.normalize_alphabet(dropoutfst, inplace=True)
		assert(dropoutfst.verify())

		## create full gam model from haplotype and dropout transducers
		gam = (fst.compose(haplofst,dropoutfst)).rmepsilon().arcsort('olabel')

		## check that the gam fst is conditionally normalized
		gamdf = gam.to_dataframe(to_real=True)
		gamdf.groupby(['state_from','ilabel']).sum()

		## create simple gam model with an i.i.d dropout model
		## (seems to perform just as well as the full gam)
		simpledrop = fst.Fst(arc_type=semiring)
		simpledrop.set_input_symbols(symbol_table)
		simpledrop.set_output_symbols(symbol_table)
		start = simpledrop.add_state()
		simpledrop.set_start(start)
		simpledrop.add_arc(start, ('a', 'a', -np.log(0.5), start))
		simpledrop.add_arc(start, ('b', 'b', -np.log(0.5), start))
		simpledrop.add_arc(start, ('a', '.', -np.log(0.5), start))
		simpledrop.add_arc(start, ('b', '.', -np.log(0.5), start))
		simpledrop.set_final(start, 0)
		simplegam = fst.compose(haplofst, simpledrop).rmepsilon().arcsort('olabel')

		return gam, simplegam, haplofst, dropoutfst

	## create random samples
	def create_random_samples(truthfsa, modelfst, n):
		""" creates a set of n random samples from a given truth T using the model M provided (i.e. randgen((T o M)->) """
		composite = fst.compose(truthfsa, modelfst)
		samplepaths = fst.randgen(composite, select='log_prob', npath=n, weight=False).project('output').rmepsilon()
		#samplepaths = [fst.randgen(composite, select='log_prob', npath=1, weight=False).project("output").rmepsilon() for i in range(n)]
		#samplepaths = fst.determinize(samplepaths).minimize()
		samplestrings = fst.tools.strings(samplepaths)
		#samplestrings = pd.concat([fst.tools.strings(s) for s in samplepaths])
		samplefsas = [fst.factory.from_string(s, arc_type=semiring, isymbols=symbol_table, osymbols=symbol_table) for s in samplestrings.string]
		#samplefsas = [fst.arcmap(s, map_type='rmweight') for s in samplepaths]
		return samplefsas, samplestrings

	## reconstruction
	def evaluate(viterbi, truthfsa):
		""" evaluates the given result (viterbi) against the truth FSA using the flip distance and some text visualisation """
		truth = fst.tools.strings(truthfsa).string.iloc[0]
		viterbiseq = fst.tools.strings(viterbi).string.iloc[0]
		truthfsa_stayflip = fst.compose(truthfsa, stayflip).project('output').rmepsilon()
		truth_stayflip = fst.tools.strings(truthfsa_stayflip).string['path0']
		viterbi_stayflip = fst.arcmap(fst.compose(viterbi, fst.arcmap(stayflip, map_type='to_std')).project('output').rmepsilon(), map_type='rmweight')
		viterbiseq_stayflip = fst.tools.strings(viterbi_stayflip).string.iloc[0]
		truth2viterbi_sf_bool = np.array([truth_stayflip[i] == viterbiseq_stayflip[i] for i in range(len(truth_stayflip))])
		truth2viterbi_stayflip = ''.join(['.' if truth2viterbi_sf_bool[i] else 'X' for i in range(len(truth_stayflip))])
		#dist_sf = np.sum(~truth2viterbi_sf_bool)
		dist_ht = float(fst.shortestdistance(fst.compose(fst.compose(fst.arcmap(fst.arcmap(truthfsa, map_type='to_std'), map_type='rmweight'), flipdist).arcsort('olabel'), fst.arcmap(viterbi, map_type='rmweight')), reverse=True)[0])
		#dist_ht = np.sum(~truth2viterbi_bool)
		return '\n'.join([truth, ' '+truth_stayflip, ' '+truth2viterbi_stayflip, ' '+viterbiseq_stayflip, viterbiseq]), dist_ht

	def probabilistic_inference(ifst, samplefsas, delta=fst.DEF_DELTA, prior=None, parallel=False):
		""" general inference method for an arbitrary model fst (e.g. gam) based on projection and intersection of all samples """
		if prior is not None:
			ifst = fst.compose(prior, ifst).arcsort('olabel')

		if prior is not None:
			sampledists = [fst.determinize(fst.compose(ifst, r).project()).minimize(delta=delta).push(remove_total_weight=True).arcsort('olabel') for r in samplefsas]
		else:
			sampledists = [fst.determinize(fst.intersect(filterfst, fst.compose(ifst, r).project())).minimize(delta=delta).push(remove_total_weight=True).rmepsilon().arcsort('olabel') for r in samplefsas]
			#sampledists = [fst.determinize(fst.intersect(filterfst, fst.compose(ifst, r).project())).minimize(delta=delta).rmepsilon().arcsort('olabel') for r in samplefsas]
			#sampledists = [fst.determinize(fst.compose(ifst, r).project()).minimize(delta=delta).rmepsilon().arcsort('olabel') for r in samplefsas]

		#viterbi solution
		#joint_detmin = fst.tools.mass_intersect_quick(sampledists, detmin=True, delta_min=delta)
		if parallel:
			joint_detmin = fst.tools.multicommand_parallel(fst.intersect, sampledists, determinize=False, minimize=True, rmepsilon=False)
		else:
			joint_detmin = fst.tools.mass_intersect_quick(sampledists, determinize=False, minimize=False, rmepsilon=False)
			#joint_detmin = functools.reduce(fst.intersect, sampledists)

		viterbi = fst.determinize(fst.shortestpath(fst.arcmap(joint_detmin, map_type='to_std')).rmepsilon())
	
		return viterbi, sampledists

	def count_bigrams(samplefsas, delta=fst.DEF_DELTA):
		""" simple inference model which counts the number of bi-grams (adjacent stay/flip transitions) """
		switchback = fst.Fst(arc_type='standard')
		switchback.set_input_symbols(symbol_table)
		switchback.set_output_symbols(symbol_table)
		switchback.add_states(3)
		switchback.set_start(0)
		switchback.set_final(1,0)
		switchback.set_final(2,0)
		switchback.add_arc(0, ('a', eps, 0, 1))
		#switch.add_arc(0, ('b', eps, 0, 2))
		switchback.add_arc(1, ('a', '=', 0, 1))
		switchback.add_arc(1, ('b', 'x', 0, 2))
		switchback.add_arc(2, ('b', '=', 0, 2))
		switchback.add_arc(2, ('a', 'x', 0, 1))
		switchback.arcsort('ilabel')
		switchback.verify()

		count = fst.Fst(arc_type='log')
		count.set_input_symbols(symbol_table)
		count.set_output_symbols(symbol_table)
		i = count.add_state()
		count.add_arc(i, ('x','_',0,i))
		count.add_arc(i, ('=','_',0,i))
		count.add_arc(i, ('x','x',-1,i))
		count.add_arc(i, ('=','=',-1,i))
		count.add_arc(i, ('x','=',0,i))
		count.add_arc(i, ('=','x',0,i))
		count.set_start(i)
		count.set_final(i,0)
		count.arcsort('olabel')
		count.verify()
		
		resolved = [fst.compose(s, stayflip).project('output').rmepsilon() for s in samplefsas]
		resolved2 = [fst.compose(count, s).project().arcsort() for s in resolved]
		final = functools.reduce(fst.intersect, resolved2)

		for state in final.states():
			weightsum = np.sum([abs(float(arc.weight)) for arc in final.arcs(state)])+2
			mai = final.mutable_arcs(state)
			for arc in mai:
				arc.weight = fst.Weight(final.weight_type(), -np.log((abs(float(arc.weight))+1)/weightsum))
				mai.set_value(arc)

		finalmap = fst.determinize(fst.compose(fst.arcmap(switchback, map_type='to_log'), final).project()).minimize()
		viterbi = fst.determinize(fst.shortestpath(fst.arcmap(final, map_type='to_std')).rmepsilon())
		#viterbimap = fst.determinize(fst.compose(switchback, viterbi).project()).minimize().rmepsilon()
		viterbimap = fst.compose(switchback, viterbi).project()
		return viterbimap, resolved2, finalmap

	def mcmc_inference(ifst, samplefsas, nsamples=1000, startfsa=None, delta=fst.DEF_DELTA):
		""" attempt at implementing an MCMC sampler to infer the joint likelihood instead of many intersections """
		jump = fst.Fst(arc_type='log')
		jump.set_input_symbols(symbol_table)
		jump.set_output_symbols(symbol_table)
		i,j = jump.add_states(2)
		jump.set_start(i)
		jump.set_final(j, 0)
		jump.add_arc(i, ('a', 'a', 0, j))
		jump.add_arc(j, ('a', 'a', -np.log(0.9), j))
		jump.add_arc(j, ('b', 'b', -np.log(0.9), j))
		jump.add_arc(j, ('a', 'b', -np.log(0.1), j))
		jump.add_arc(j, ('b', 'a', -np.log(0.1), j))
		#jump.add_arc(j, (0, 0, 0, i))
		jump.arcsort()

		sampledists = [fst.determinize(fst.intersect(filterfst, fst.compose(ifst, r).project())).minimize(delta=delta).push(remove_total_weight=True).rmepsilon().arcsort('olabel') for r in samplefsas]
		confusions = [fst.algos.posterior_decoding(ifst) for ifst in sampledists]
		confusions_intersect = fst.tools.mass_intersect_quick(confusions)

		#longfsa = sampledists[0].copy()
		#for i in range(1, len(sampledists)):
		#	longfsa.concat(sampledists[i])

		## start point
		if startfsa is None:
			current = fst.arcmap(fst.arcmap(fst.determinize(fst.shortestpath(fst.arcmap(confusions_intersect, map_type='to_std'))), map_type='to_log'), map_type='rmweight')
		else:
			current = fst.arcmap(startfsa, map_type='rmweight')

		p_current = -np.sum([float(fst.shortestdistance(fst.intersect(current, s), reverse=True)[0]) for s in sampledists])
		#p_current = -float(fst.shortestdistance(fst.intersect(current.copy().closure().rmepsilon(), longfsa), reverse=True)[0])
		
		posterior = [current]
		posteriorll = [p_current]

		for i in range(nsamples):
			## draw random sample
			randsamp = np.random.randint(0, len(samplefsas))
			proposal = fst.arcmap(fst.randgen(fst.compose(current, jump), select='log_prob').project('output'), map_type='rmweight')
			#proposal = fst.arcmap(fst.randgen(fst.compose(fst.compose(current, jump).arcsort('olabel'),sampledists[randsamp]), select='log_prob').project("output"), map_type='rmweight')
			p_proposal = -np.sum([float(fst.shortestdistance(fst.intersect(proposal, s), reverse=True)[0]) for s in sampledists])
			#p_proposal = -float(fst.shortestdistance(fst.intersect(proposal.copy().closure().rmepsilon(), longfsa), reverse=True)[0])

			p_accept = np.exp(p_proposal - p_current)
			accept = np.random.rand() < p_accept
			
			#logging.warn(p_proposal)
			if accept:
				current = proposal
				p_current = p_proposal
				#logging.warn('accept')
			posterior.append(current)
			posteriorll.append(p_current)

		return posterior, posteriorll

	## make gam fsts
	p_drop_extend = 0.8
	p_drop = 0.5
	p_switch = 0.1
	p_extend = 0.3
	gam, simplegam, haplofst, dropoutfst = create_gam_fst(p_switch, p_extend, p_drop, p_drop_extend)

	truth = 'a'+''.join(np.random.choice(['a','b'], 50))
	truthfsa = fst.factory.from_string(truth, arc_type=semiring, isymbols=symbol_table, osymbols=symbol_table)

	nsamples = 50
	samplefsas, samplestrings = create_random_samples(truthfsa, gam, nsamples)
	samplestrings

	#samplestrings = ['a.b.b', '.a.b.', '.b..b']
	#samplefsas = [fst.factory.from_string(s, isymbols=symbol_table, osymbols=symbol_table) for s in samplestrings]

	print('Bi-grams:')
	start = time.time()
	viterbi, sampledists, prior = count_bigrams(samplefsas)
	end = time.time()
	evalstr, dist = evaluate(viterbi, truthfsa)
	print(evalstr)
	print('Distance: ' + str(dist))
	print('Time: ' + str(end-start))

	print('Simple GAM:')
	start = time.time()
	viterbi, sampledists = probabilistic_inference(simplegam, samplefsas)
	end = time.time()
	evalstr, dist = evaluate(viterbi, truthfsa)
	print(evalstr)
	print('Distance: ' + str(dist))
	print('Time: ' + str(end-start))

	print('Random:')
	start = time.time()
	viterbi = fst.factory.from_string('a'+''.join(np.random.choice(['a','b'], len(truth)-1)), isymbols=symbol_table, osymbols=symbol_table, arc_type='standard')
	end = time.time()
	evalstr, dist = evaluate(viterbi, truthfsa)
	print(evalstr)
	print('Distance: ' + str(dist))
	print('Time: ' + str(end-start))

	## accuracy / speed test
	nsteps = 400
	seqlen = 50

	length = np.zeros(nsteps)
	elapsed = np.zeros((3, nsteps))
	dist = np.zeros((3, nsteps))
	for i in range(nsteps):
		print(i)
		#truth = 'abbabaaaaabbababaababb' * 2
		truth = 'a'+''.join(np.random.choice(['a','b'], seqlen-1))
		#truth = truth*((i+1)*2)
		#truth = truth*(i+1)
		length[i] = len(truth)
		truthfsa = fst.factory.from_string(truth, arc_type=semiring, isymbols=symbol_table, osymbols=symbol_table)
		samplefsas, _ = create_random_samples(truthfsa, gam, 10)
		
		## bigrams
		start = time.time()
		viterbi, _, prior = count_bigrams(samplefsas)
		end = time.time()
		elapsed[2, i] = end-start
		_, d = evaluate(viterbi, truthfsa)
		dist[2, i] = d

		## simple gam
		start = time.time()
		viterbi, _ = probabilistic_inference(simplegam, samplefsas)
		end = time.time()
		elapsed[1, i] = end-start
		_, d = evaluate(viterbi, truthfsa)
		dist[1, i] = d

		## random control
		start = time.time()
		viterbi = fst.factory.from_string('a'+''.join(np.random.choice(['a','b'], seqlen-1)), isymbols=symbol_table, osymbols=symbol_table, arc_type='standard')
		end = time.time()
		elapsed[0, i] = end-start
		_, d = evaluate(viterbi, truthfsa)
		dist[0, i] = d

	models = ['control', 'gam', 'bi-gram']
	indexnames = ['model','sample']
	df = pd.DataFrame(elapsed, index=models).stack().to_frame('elapsed')
	df.index.names=indexnames
	df['dist'] = pd.DataFrame(dist, index=models).stack()
	lengths=pd.DataFrame(length, columns=['length'])
	lengths.index.name = 'sample'
	df = df.join(lengths)

	sns.boxplot(df.index.get_level_values('model'), df.dist, notch=True)
	g=sns.boxplot(df.index.get_level_values('model'), df.elapsed, notch=True)
	g.set_yscale('log')
	pl.hexbin(df.loc['bi-gram','dist'], df.loc['gam','dist'], gridsize=10)

if __name__ == "__main__":
	main()
