#%% imports
import sys
import numpy as np
sys.path.append('../..')
import fstlib

#%% create symbol table
symbols = ['H','T', 'f', 'u'] # head, tail, fair, unfair
symbol_table = fstlib.SymbolTable()
symbol_table.add_symbol('<eps>')
for s in symbols:
	symbol_table.add_symbol(s)


#%% create casino HMM
hmm = fstlib.Fst(arc_type='standard')
hmm.set_input_symbols(symbol_table)
hmm.set_output_symbols(symbol_table)
hmm.add_states(2)
hmm.set_start(0)
hmm.set_final(0,0)
hmm.set_final(1,0)
hmm.add_arc(0,('H','f',-np.log(0.35),0))
hmm.add_arc(0,('T','f',-np.log(0.35),0))
hmm.add_arc(0,('H','u',-np.log(0.15),1))
hmm.add_arc(0,('T','u',-np.log(0.15),1))
hmm.add_arc(1,('H','u',-np.log(0.48),1))
hmm.add_arc(1,('T','u',-np.log(0.12),1))
hmm.add_arc(1,('H','f',-np.log(0.2),0))
hmm.add_arc(1,('T','f',-np.log(0.2),0))
hmm.verify()
display(hmm)
print(hmm.to_real())
#hmm.view()

#%% create data
seq = 'HHTHTHTHHHHHHTHHHHHHHTHTHHTHT'
seqfsa = fstlib.factory.from_string(seq, arc_type='standard', isymbols=symbol_table, osymbols=symbol_table)
seqfsa.verify()
print(seqfsa)
display(seqfsa)

#%% Build HMM
unrolled = fstlib.compose(seqfsa, hmm)
print(unrolled)
display(unrolled)
#unrolled.view()

#%% Viterbi solution
sp = fstlib.determinize(fstlib.shortestpath(unrolled))
print(sp.to_real())
fstlib.tools.strings(sp, tape='input')
fstlib.tools.strings(sp, tape='output')

#%% MAP solution
posterior = fstlib.algos.posterior_decoding(unrolled.project('output'))

map = fstlib.determinize(fstlib.shortestpath(fstlib.arcmap(posterior, map_type='to_standard')))
fstlib.tools.strings(map, tape='input')
fstlib.tools.strings(map, tape='output')



# %%
