# %%
import medicc
import fstlib


separator = "X"
max_cn = 8
max_pre_wgd_losses = 7

# FST weights
w_stay_gain = 0
w_open_gain = 1
w_extend_gain = 0.1

w_stay_loss = 0
w_open_loss = 1
w_extend_loss = 0.1

wgd_cost = 1

# %%
SYMBOL_TABLE = medicc.create_symbol_table(max_cn=max_cn, separator=separator)

n = len(medicc.factory._get_int_cns_from_symbol_table(SYMBOL_TABLE, separator))
G1step = ~medicc.factory.create_1step_del_fst(SYMBOL_TABLE,
                                              separator=separator,
                                              exclude_zero=True,
                                              w_stay=w_stay_gain,
                                              w_open=w_open_gain,
                                              w_extend=w_extend_gain)

L_LOH_1step = medicc.factory.create_1step_del_fst(SYMBOL_TABLE,
                                                  separator=separator,
                                                  exclude_zero=False,
                                                  w_stay=w_stay_loss,
                                                  w_open=w_open_loss,
                                                  w_extend=w_extend_loss)

G = medicc.create_nstep_fst(n - 1, G1step)
L_LOH = medicc.create_nstep_fst(n - 1, L_LOH_1step)

# %%
# XX = fstlib.encode_determinize_minimize(L_LOH * G) # NOT gonna work
XX = L_LOH * G


# %%
W_1step = medicc.factory.create_1step_WGD_fst(SYMBOL_TABLE, separator,
                                              wgd_cost=wgd_cost,
                                              minimize=False,
                                              wgd_x2=True,
                                              total_cn=False)

MEDICC2_FST_with_extend_modeling = W_1step * XX

MEDICC2_FST_with_extend_modeling.write("gain_loss_extend_0.1.fst")