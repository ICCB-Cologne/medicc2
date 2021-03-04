import logging

import pandas as pd

import fstlib
import fstlib.cext.ops
from fstlib.cext import pywrapfst

logger = logging.getLogger(__name__)

try:
    from fstlib.cext.ops import *
    logger.info('Using accelerated C operations.')
except ImportError:
    logger.warn('Accelerated C operations not available.')

def align(model, ifst1, ifst2):
    ofst = fstlib.cext.ops.align(model.fst, ifst1.fst, ifst2.fst)
    return fstlib.Fst(ofst)

def encode_determinize_minimize(ifst, delta=1e-6):
    em = fstlib.EncodeMapper(arc_type=ifst.arc_type(), encode_labels=True, encode_weights=False)
    ofst = ifst.copy().encode(em)
    ofst = fstlib.determinize(ofst, delta=delta)
    ofst.minimize()
    ofst.decode(em)
    return ofst

def info(ifst, name=None):
    """ Outputs some info stats """
    info = {}
    info['fst_type'] = ifst.fst_type()
    info['arc_type'] = ifst.arc_type()
    info['nstates'] = ifst.num_states()
    info['narcs'] = 0
    info['ninputeps'] = 0
    info['noutputeps'] = 0
    for state in ifst.states():
        info['narcs'] += ifst.num_arcs(state)
        info['ninputeps'] += ifst.num_input_epsilons(state)
        info['noutputeps'] += ifst.num_output_epsilons(state)
    
    df = pd.DataFrame.from_dict(info, orient='index', columns=[name if name is not None else id(ifst)])
    return df

def kernel_score(model, ifst1, ifst2):
    if model.arc_type()=='standard':
        distance = fstlib.cext.ops.kernel_score_std(model.fst, ifst1.fst, ifst2.fst)
    elif model.arc_type()=='log':
        ##distance = fstlib.cext.ops.score_log(model.fst, ifst1.fst, ifst2.fst)
        distance = None
    else:
        distance = None
    return distance

def normalize_alphabet(ifst, inplace=True):
    """Normalizes fst so that outgoing transition weights of the same input symbol sum to 1"""
    
    ## convert to real if necessary
    if inplace:
        ofst = ifst
    else:
        ofst = ifst.copy()
        
    ofst.weight_map(fstlib.algos.map_log_to_real)

    for state in ofst.states():
        arcweights = [(arc.ilabel, float(arc.weight)) for arc in ofst.arcs(state)]
        labelweights = pd.DataFrame.from_records(arcweights, columns=['ilabel','weight']).groupby('ilabel').sum()
        mai = ofst.mutable_arcs(state)
        for arc in mai:
            labelsum = labelweights.weight.loc[arc.ilabel]
            arc.weight = fstlib.Weight(ofst.weight_type(), float(arc.weight) / labelsum)
            mai.set_value(arc)

    ## convert back
    ofst.weight_map(fstlib.algos.map_real_to_log)
    return ofst

def read(source):
    newfst = pywrapfst.Fst.read(source)
    return fstlib.Fst(newfst)

def score(model, ifst1, ifst2):
    if model.arc_type()=='standard':
        distance = fstlib.cext.ops.score_std(model.fst, ifst1.fst, ifst2.fst)
    elif model.arc_type()=='log':
        distance = fstlib.cext.ops.score_log(model.fst, ifst1.fst, ifst2.fst)
    else:
        pass
    return distance

def weight_map(ifst, func):
    newfst = ifst.copy()
    newfst.weight_map(func)
    return newfst
