#%% imports
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import fstlib
import medicc

#%% example FST (Supp Fig 1a)
symbols = ['a','b'] 
st = fstlib.SymbolTable()
st.add_symbol('-')
for s in symbols:
        st.add_symbol(s)
fst = fstlib.Fst(arc_type='standard')
fst.set_input_symbols(st)
fst.set_output_symbols(st)
fst.add_states(2)
fst.set_start(0)
fst.set_final(0,0)
fst.set_final(1,0)
fst.add_arc(0, ('a', 'a', 0, 0))
fst.add_arc(0, ('a', 'b', 1, 1))
fst.add_arc(1, ('a', 'b', 1, 1))
fst

#%% manually create the fst for the example to keep it small
sep='X'
symbol_table = medicc.create_symbol_table(3, separator=sep)
n = len(medicc.factory._get_int_cns_from_symbol_table(symbol_table))
X1step = medicc.create_1step_del_fst(symbol_table, sep, exclude_zero=True)
X = medicc.create_nstep_fst(2, X1step)
XX = fstlib.encode_determinize_minimize(X*~X)
LOH = fstlib.encode_determinize_minimize(medicc.create_loh_fst(symbol_table, sep))
W1step = medicc.create_1step_WGD_fst(symbol_table, sep)
W = medicc.create_nstep_fst(3, W1step)
T = LOH * W * XX
T.info()

#%% Create genomes
g1 = fstlib.factory.from_string("1111X1111", arc_type="standard", isymbols=symbol_table, osymbols=symbol_table)
g2 = fstlib.factory.from_string("2033X2112", arc_type="standard", isymbols=symbol_table, osymbols=symbol_table)
# U is shown in Supp Fig 1h
U = g1 * T * g2
U

#%% final plotting
def fst_to_pdf(fst, acceptor=False, ranksep=0.4, nodesep=0.25, fontsize=14):
    g = fst.to_graphviz(ranksep=ranksep, nodesep=nodesep, acceptor=acceptor, fontsize=fontsize)
    g.format='pdf'
    return g.pipe()

with open('../Figures_Kaufmann_et_al_2021/final_figures/Supp_1_method_example_fst_gainloss.pdf', 'wb') as f:
    f.write(fst_to_pdf(XX, ranksep=0.25, nodesep=0.1, fontsize=12))
with open('../Figures_Kaufmann_et_al_2021/final_figures/Supp_1_method_example_fst_loh.pdf', 'wb') as f:
    f.write(fst_to_pdf(LOH, fontsize=12))
with open('../Figures_Kaufmann_et_al_2021/final_figures/Supp_1_method_example_fst_wgd.pdf', 'wb') as f:
    f.write(fst_to_pdf(W1step, fontsize=12))
with open('../Figures_Kaufmann_et_al_2021/final_figures/Supp_1_method_example_fsa1.pdf', 'wb') as f:
    f.write(fst_to_pdf(g1, acceptor=True, fontsize=10))
with open('../Figures_Kaufmann_et_al_2021/final_figures/Supp_1_method_example_fsa2.pdf', 'wb') as f:
    f.write(fst_to_pdf(g2, acceptor=True, fontsize=10))
with open('../Figures_Kaufmann_et_al_2021/final_figures/Supp_1_method_example_composition.pdf', 'wb') as f:
    f.write(fst_to_pdf(U, fontsize=12))
with open('../Figures_Kaufmann_et_al_2021/final_figures/Supp_1_method_example_simple_fst.pdf', 'wb') as f:
    f.write(fst_to_pdf(fst, fontsize=12))

