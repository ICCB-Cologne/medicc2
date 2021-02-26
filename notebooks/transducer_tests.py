#%% imports and configurations
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import medicc
import medicc.test

#%% define alphabet and build WGD transducer
maxcn = 8 ## alphabet = 012345678; maxcn losses, maxcn-1 gains
maxwgd = 8
sep = "X"
symbol_table = medicc.create_symbol_table(maxcn, sep)

I = medicc.create_copynumber_fst(symbol_table, sep, enable_wgd=False)
T = medicc.create_copynumber_fst(symbol_table, sep, enable_wgd=True)
print(T.info('T').join(I.info('I')))

#%%
medicc.test.run_pair_tests(I, is_wgd=False)
medicc.test.run_pair_tests(T, is_wgd=True)

# %%
