#%% imports
import sys
import numpy as np
from sklearn.manifold import MDS
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append('../..')
sys.path.append('../../../GenomicsToolkit')
import fstlib
import genometools as gtools

#%% symbol table
symbols = ['1','2', '3', '|'] # unmeth, partially meth, fully meth
symbol_table = fstlib.SymbolTable()
symbol_table.add_symbol('-')
for s in symbols:
	symbol_table.add_symbol(s)

#%% create FST
T1 = fstlib.Fst(arc_type='standard')
T1.set_input_symbols(symbol_table)
T1.set_output_symbols(symbol_table)
T1.add_states(1)
T1.set_start(0)
T1.set_final(0, 0)
T1.add_arc(0, ('1','1', 0, 0))
T1.add_arc(0, ('2','2', 0, 0))
T1.add_arc(0, ('3','3', 0, 0))
T1.add_arc(0, ('|','|', 0, 0))

T1.add_arc(0, ('1','2', 1, 0))
T1.add_arc(0, ('2','1', 1, 0))
T1.add_arc(0, ('1','3', 2, 0))
T1.add_arc(0, ('3','1', 2, 0))
T1.add_arc(0, ('2','3', 1, 0))
T1.add_arc(0, ('3','2', 1, 0))

T1.add_arc(0, ('-','1', 10, 0))
T1.add_arc(0, ('-','2', 10, 0))
T1.add_arc(0, ('-','3', 10, 0))
T1.add_arc(0, ('-','|', 11, 0))
T1.add_arc(0, ('1','-', 10, 0))
T1.add_arc(0, ('2','-', 10, 0))
T1.add_arc(0, ('3','-', 10, 0))
T1.add_arc(0, ('|','-', 11, 0))

s1 = fstlib.factory.from_string('11|33333|11|2222', isymbols=symbol_table, osymbols=symbol_table, arc_type=fstlib.Semiring.TROPICAL)
s2 = fstlib.factory.from_string('223|333|111111|1', isymbols=symbol_table, osymbols=symbol_table, arc_type=fstlib.Semiring.TROPICAL)

c1 = fstlib.compose(s1, T1)
final = fstlib.compose(c1, s2)
sp = fstlib.shortestpath(final)
sp
spin = fstlib.tools.strings(sp, tape='input')
spout = fstlib.tools.strings(sp, tape='output')
print('%s\n%s' % (spin.iloc[0,0], spout.iloc[0,0]))
print(spin.iloc[0,1])
# %%
