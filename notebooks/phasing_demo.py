#%% imports and configurations
import os
import sys

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# %matplotlib inline
import fstlib

#%%
T_wgd_asymm = fstlib.Fst.read('../objects/wgd_asymm.fst')
T_no_wgd_asymm = fstlib.Fst.read('../objects/no_wgd_asymm.fst')

# %% phasing solution
allele_a = '12102'
allele_b = '21013'
encoded = np.array([list(zip(allele_a, allele_b)), list(zip(allele_b, allele_a))])

td = fstlib.factory.from_array(np.array([['1',]*len(allele_a)]), symbols=T_wgd_asymm.input_symbols(), arc_type='standard')
tg = fstlib.factory.from_array(encoded, symbols=T_wgd_asymm.input_symbols(), arc_type='standard')

left = (td * T_wgd_asymm).project('output')
right = (~T_wgd_asymm * td).project('input')
phase = fstlib.shortestpath(left * tg * right)
phase_str = [''.join(x) for x in zip(*fstlib.tools.paths(phase, 'both')[0][0])]
print('final phasing: {} - {}'.format(*phase_str))
