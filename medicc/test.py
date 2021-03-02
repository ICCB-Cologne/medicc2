import logging
from dataclasses import dataclass

import fstlib
import numpy as np

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

@dataclass
class PairTest:
    seq_in: str
    seq_out: str
    expected_score_wgd: float
    expected_score_no_wgd: float

PAIR_TESTS = [
    PairTest('11111X1111', '22022X1111', 2, 2),
    PairTest('11111X1111', '22022X2222', 2, 3),
    PairTest('11101X1111', '10111X1111', 2, 2),
    PairTest('33233X1111', '00000X1111', 3, 3),
    PairTest('1111111111X1111111111', '2212222222X2222222222', 2, 3),
    PairTest('2222222222X2222222222', '3323333333X3333323333', 3, 4),
    PairTest('1111111X11X11X1111X1111', '3322112X22X23X2222X2200', 5, 9),
    PairTest('1111111111X111X111', '3332222221X333X333', 4, 6),
]

def run_pair_tests(fst: fstlib.Fst, is_wgd: bool):
    """ Runs all pair tests and reports the output. """
    test_results = np.array([_run_pair_test(test, fst, is_wgd) for test in PAIR_TESTS])
    return test_results

def _run_pair_test(pair_test: PairTest, fst: fstlib.Fst, is_wgd: bool) -> bool:
    """ Runs individual pair test and returns if passed. """
    td = fstlib.factory.from_string(pair_test.seq_in, isymbols=fst.input_symbols(), osymbols=fst.output_symbols(), arc_type=fst.arc_type())
    tg = fstlib.factory.from_string(pair_test.seq_out, isymbols=fst.input_symbols(), osymbols=fst.output_symbols(), arc_type=fst.arc_type())
    test_score = float(fstlib.score(fst, td, tg))
    if is_wgd:
        expected_score = pair_test.expected_score_wgd
    else:
        expected_score = pair_test.expected_score_no_wgd
    passed = False
    if test_score == expected_score:
        passed = True
    if passed:
        logfun = logger.info
    else:
        logfun = logger.warn
    logfun("Testing MED between %s and %s. Expected distance: %f (%s), returned distance: %f. Test %s!", 
        pair_test.seq_in, 
        pair_test.seq_out, 
        expected_score, 
        'WGD' if is_wgd else 'no WGD', 
        test_score, 
        'passed' if passed else 'failed')
    return passed

