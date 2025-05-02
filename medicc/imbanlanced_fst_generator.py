import medicc 
import medicc.factory
import numpy as np 
import fstlib
from fractions import Fraction 


# Parameters
output_fst_path_no_LOH = "objects/imbalanced_fst.fst"
output_fst_path_LOH = "objects/imbalanced_fst_LOH.fst"
p_gain = 0.8 # probability of gain
LOH = True


def convert_gain_loss_probability_to_integer(p_gain: float) -> tuple[int, int]:
    """
    Converts a probability of gain to an integer ratio of P(gain) : P(loss)
    that closely approximates the ratio p_gain : (1 - p_gain).
    
    Args:
        p_gain (float): The probability of gain (0 <= p_gain <= 1).
        
    Returns:
        (int, int): The integer ratio representing P(gain) : P(loss).
    """
    p_loss = 1 - p_gain 

    s_gain = 1 / p_gain 
    s_loss = 1 / p_loss 

    score_fraction = Fraction(s_gain / s_loss).limit_denominator()

    return score_fraction.numerator, score_fraction.denominator


separator = "X"
max_cn = 8
max_pre_wgd_losses = max_cn - 1

SYMBOL_TABLE = medicc.create_symbol_table(max_cn, separator)

n = len(medicc.factory._get_int_cns_from_symbol_table(SYMBOL_TABLE, separator))

s_gain, s_loss = convert_gain_loss_probability_to_integer(p_gain)
G1step_event_open_altered = ~medicc.factory.create_1step_del_fst(SYMBOL_TABLE, separator, exclude_zero=True, w_stay = 0, w_open = s_gain, w_extend = 0)
L_LOH_1step_event_open_altered = medicc.factory.create_1step_del_fst(SYMBOL_TABLE, separator, exclude_zero=False, w_stay = 0, w_open = s_loss, w_extend = 0)

if LOH:
    L_LOH_event_open_altered = medicc.create_nstep_fst(n - 1, L_LOH_1step_event_open_altered, minimize=True)
    G_event_open_altered = medicc.create_nstep_fst(n - 1, G1step_event_open_altered, minimize=True)
    L_LOH_G_open_altered = L_LOH_event_open_altered * G_event_open_altered # Does the order matter?
    L_LOH_G_open_altered_edm = fstlib.ext.encode_determinize_minimize(L_LOH_G_open_altered)
    L_LOH_G_open_altered_edm.write(output_fst_path_LOH)
else:
    l1_step_event_open_altered = medicc.factory.create_1step_del_fst(SYMBOL_TABLE, separator, exclude_zero=True, w_stay = 0, w_open = s_loss, w_extend = 0)
    G_event_open_altered = medicc.create_nstep_fst(n - 1, G1step_event_open_altered, minimize=True)
    L_event_open_altered = medicc.create_nstep_fst(n - 1, l1_step_event_open_altered, minimize=True)
    LG_event_open_altered = L_event_open_altered * G_event_open_altered # Does the order matter?
    LG_event_open_altered_edm = fstlib.ext.encode_determinize_minimize(LG_event_open_altered)
    LG_event_open_altered_edm.write(output_fst_path_no_LOH)

