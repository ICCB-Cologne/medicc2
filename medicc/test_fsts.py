from dataclasses import dataclass

import fstlib
import numpy as np
import pytest

import medicc


@dataclass
class PairTest:
    seq_in: str
    seq_out: str
    expected_score_wgd: float
    expected_score_no_wgd: float


PAIR_TESTS = [
    PairTest('11111X1111', '22022X1111', 2, 2),
    PairTest('11111X1111', '22022X2222', 2, 3),
    PairTest('11101X1111', '10111X1111', np.inf, np.inf), # expected distance 2 for symmetric metric
    PairTest('33233X1111', '00000X1111', 3, 3),
    PairTest('1111111111X1111111111', '2212222222X2222222222', 2, 3),
    PairTest('2222222222X2222222222', '3323333333X3333323333', 3, 4),
    PairTest('1111111X11X11X1111X1111', '3322112X22X23X2222X2200', 5, 9),
    PairTest('1111111111X111X111', '3332222221X333X333', 4, 6),
]

# TODO: symmetric, i.e. fstlib.kernel_score
def _run_pair_test(pair_test: PairTest, is_wgd: bool) -> bool:
    """Runs individual pair test"""
    maxcn = 8  # alphabet = 012345678; maxcn losses, maxcn-1 gains
    sep = "X"
    symbol_table = medicc.create_symbol_table(maxcn, sep)
    fst = medicc.create_copynumber_fst(symbol_table, sep, enable_wgd=is_wgd)

    td = fstlib.factory.from_string(pair_test.seq_in, isymbols=fst.input_symbols(), osymbols=fst.output_symbols(), arc_type=fst.arc_type())
    tg = fstlib.factory.from_string(pair_test.seq_out, isymbols=fst.input_symbols(), osymbols=fst.output_symbols(), arc_type=fst.arc_type())
    test_score = float(fstlib.score(fst, td, tg))
    if is_wgd:
        expected_score = pair_test.expected_score_wgd
    else:
        expected_score = pair_test.expected_score_no_wgd
    return test_score, expected_score


@pytest.mark.parametrize("pair", PAIR_TESTS)
def test_fstlib_with_wgd(pair: PairTest):

    test_score, expected_score = _run_pair_test(pair, is_wgd=True)
    assert test_score == expected_score

@pytest.mark.parametrize("pair", PAIR_TESTS)
def test_fstlib_without_wgd(pair: PairTest):

    test_score, expected_score = _run_pair_test(pair, is_wgd=False)
    assert test_score == expected_score
