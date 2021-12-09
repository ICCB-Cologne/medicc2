import logging

import fstlib
import numpy as np

from medicc import tools

SYMBOL_GAP = 'gap'

# prepare logger
logger = logging.getLogger(__name__)


def create_symbol_table(max_cn=8, separator='X'):
    alphabet = [tools.int2hex(x) for x in range(max_cn+1)]
    symbol_table = fstlib.SymbolTable()
    symbol_table.add_symbol(SYMBOL_GAP)
    for s in alphabet:
        symbol_table.add_symbol(s)

    if separator is not None and separator != '':
        symbol_table.add_symbol(separator)
    return symbol_table


def _get_int_cns_from_symbol_table(symbol_table, separator='X'):
    cns = {x[1]: tools.hex2int(x[1]) for x in symbol_table if x[0] != 0 and x[1] != separator}
    return cns


def create_1step_amp_fst(symbol_table, separator='X'):

    cns = _get_int_cns_from_symbol_table(symbol_table, separator)

    myfst = fstlib.Fst(arc_type='standard')
    myfst.set_input_symbols(symbol_table)
    myfst.set_output_symbols(symbol_table)

    myfst.add_states(2)
    myfst.set_start(0)
    myfst.set_final(0, 0)
    myfst.set_final(1, 0)

    myfst.add_arcs(0, [(s, s, 0, 0) for s in cns.keys()])
    if separator is not None and separator != '':
        myfst.add_arc(0, (separator, separator, 0, 0))  # add separator
    myfst.add_arcs(0, [(s, t, 1, 1) for s in cns.keys() for t in cns.keys()
                       if ((cns[t]-cns[s]) == 1) & (s != '0')])  # transitions from 0->1
    myfst.add_arcs(1, [(s, t, 0, 1) for s in cns.keys() for t in cns.keys() if (
        (cns[t]-cns[s]) == 1) & (s != '0')])  # extension of an open window
    if '0' in cns.keys():
        myfst.add_arc(1, ('0', '0', 0, 1))
    myfst.add_arcs(1, [(s, s, 0, 0) for s in cns.keys() if s != '0'])  # return
    if separator is not None and separator != '':
        myfst.add_arc(1, (separator, separator, 0, 0))

    return myfst


def create_1step_del_fst(symbol_table, separator='X', exclude_zero=False):

    cns = _get_int_cns_from_symbol_table(symbol_table, separator)

    myfst = fstlib.Fst(arc_type='standard')
    myfst.set_input_symbols(symbol_table)
    myfst.set_output_symbols(symbol_table)

    myfst.add_states(2)
    myfst.set_start(0)
    myfst.set_final(0, 0)
    myfst.set_final(1, 0)

    myfst.add_arcs(0, [(s, s, 0, 0) for s in cns.keys()])
    if separator is not None and separator != '':
        myfst.add_arc(0, (separator, separator, 0, 0))  # add separator
    if exclude_zero:
        myfst.add_arcs(0, [(s, t, 1, 1) for s in cns.keys()
                           for t in cns.keys() if (cns[s]-cns[t]) == 1 and t != '0'])  # 0->1
        myfst.add_arcs(1, [(s, t, 0, 1) for s in cns.keys() for t in cns.keys()
                           if (cns[s]-cns[t]) == 1 and t != '0'])  # extend an open window
    else:
        myfst.add_arcs(0, [(s, t, 1, 1) for s in cns.keys()
                           for t in cns.keys() if (cns[s]-cns[t]) == 1])  # 0->1
        myfst.add_arcs(1, [(s, t, 0, 1) for s in cns.keys()
                           for t in cns.keys() if (cns[s]-cns[t]) == 1])  # extend an open window
    if '0' in cns.keys():
        myfst.add_arc(1, ('0', '0', 0, 1))
    myfst.add_arcs(1, [(s, s, 0, 0) for s in cns.keys() if s != '0'])  # return
    if separator is not None and separator != '':
        myfst.add_arc(1, (separator, separator, 0, 0))

    return myfst


def create_loh_fst(symbol_table, separator='X'):

    cns = _get_int_cns_from_symbol_table(symbol_table, separator)

    myfst = fstlib.Fst(arc_type='standard')
    myfst.set_input_symbols(symbol_table)
    myfst.set_output_symbols(symbol_table)

    myfst.add_states(9)
    myfst.set_start(0)
    for i in range(9):
        myfst.set_final(i, 0)

    ## match 0 -> 0
    myfst.add_arcs(0, [(s, s, 0, 0) for s in cns.keys()])
    if separator is not None and separator != '':
        myfst.add_arc(0, (separator, separator, 0, 0))  # add separator

    ## others
    for state in range(1, 9):
        cost = state
        myfst.add_arcs(0, [(s, t, cost, state) for s in cns.keys()
                           for t in cns.keys() if (cns[s]-cns[t]) == cost and t == '0'])
        myfst.add_arcs(state, [(s, t, 0, state) for s in cns.keys()
                               for t in cns.keys() if (cns[s]-cns[t]) <= cost and t == '0'])
        myfst.add_arcs(state, [(s, t, cns[s]-cns[t]-cost, cns[s]-cns[t]) for s in cns.keys()
                               for t in cns.keys() if (cns[s]-cns[t]) > cost and t == '0'])
        myfst.add_arcs(state, [(s, s, 0, 0) for s in cns.keys()])
        if separator is not None and separator != '':
            myfst.add_arc(state, (separator, separator, 0, 0))

    return myfst


def create_1step_WGD_fst(symbol_table, separator='X', wgd_cost=1, minimize=True):

    cns = _get_int_cns_from_symbol_table(symbol_table, separator)

    W = fstlib.Fst()
    W.set_input_symbols(symbol_table)
    W.set_output_symbols(symbol_table)
    W.add_states(3)
    W.set_start(0)
    W.set_final(0, 0)
    W.set_final(1, 0)
    W.set_final(2, 0)
    W.add_arcs(0, [(s, t, wgd_cost, 1) for s in cns.keys()
                   for t in cns.keys() if (s != '0') and ((cns[t]/cns[s]) == 2.)])
    W.add_arc(1, ('0', '0', 0, 1))
    W.add_arcs(1, [(s, t, 0, 1) for s in cns.keys()
                   for t in cns.keys() if (s != '0') and ((cns[t]/cns[s]) == 2.)])
    W.add_arc(0, ('0', '0', 0, 0))
    if separator is not None and separator != '':
        W.add_arc(0, (separator, separator, 0, 0))
        W.add_arc(1, (separator, separator, 0, 1))
    W.add_arc(1, ('0', '0', 0, 1))

    W.add_arcs(0, [(s, s, 0, 2) for s in cns.keys() if s != '0'])
    W.add_arcs(2, [(s, s, 0, 2) for s in cns.keys()])
    W.add_arc(2, ('0', '0', 0, 2))
    if separator is not None and separator != '':
        W.add_arc(2, (separator, separator, 0, 2))
    W.arcsort('olabel')
    if minimize:
        W = fstlib.encode_determinize_minimize(W)

    return W


def create_filter_fst(allowed_symbols, symbol_table, separator='X'):
    F = fstlib.Fst()
    F.set_input_symbols(symbol_table)
    F.set_output_symbols(symbol_table)
    F.add_states(1)
    F.set_start(0)
    F.set_final(0, 0)
    F.add_arcs(0, [(s, s, 0, 0) for s in allowed_symbols])
    F.add_arc(0, (separator, separator, 0, 0))

    return F


def create_nstep_fst(n, one_step_fst, minimize=True):
    # Extend 1step FST
    nstep_fst = one_step_fst
    nstep_fst.arcsort(sort_type='olabel')

    for _ in range(n):
        nstep_fst = fstlib.compose(nstep_fst, one_step_fst)
        if minimize:
            nstep_fst = fstlib.encode_determinize_minimize(nstep_fst)
        nstep_fst.arcsort(sort_type='olabel')

    return nstep_fst


def create_copynumber_fst(symbol_table, sep='X', enable_wgd=False, wgd_cost=1, max_num_wgds=3):
    """ Creates the tree FST T which computes the asymmetric MED. """
    n = len(_get_int_cns_from_symbol_table(symbol_table, sep))
    X1step = create_1step_del_fst(symbol_table, sep, exclude_zero=True)
    X = create_nstep_fst(n-1, X1step)
    XX = fstlib.encode_determinize_minimize(X*~X)
    LOH = create_loh_fst(symbol_table, sep)

    if enable_wgd:
        W1step = create_1step_WGD_fst(symbol_table, sep, wgd_cost=wgd_cost,
                                      minimize=False)
        n = min(max_num_wgds, np.floor(np.log2(n)))
        if n > 1:
            W = create_nstep_fst(n-1, W1step)
        else:
            W = W1step
        T = LOH * W * XX
    else:
        T = LOH * XX

    return T


def create_phasing_fsa_from_strings(allele_a, allele_b, symbol_table, sep='X'):
    myfst = fstlib.Fst(arc_type='standard')
    myfst.set_input_symbols(symbol_table)
    myfst.set_output_symbols(symbol_table)
    prev_state = myfst.add_state()
    myfst.set_start(0)
    for a, b in zip(allele_a, allele_b):
        curr_state = myfst.add_state()
        myfst.add_arc(prev_state, (a, a, 0, curr_state))
        myfst.add_arc(prev_state, (b, b, 0, curr_state))
        prev_state = curr_state

    myfst.set_final(prev_state, 0)
    return myfst
