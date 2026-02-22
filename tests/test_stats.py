from __future__ import annotations

import numpy as np

from src.utils.stats import coloc_posteriors, wakefield_abf


def test_wakefield_abf_zero_beta_matches_closed_form() -> None:
    beta = np.array([0.0, 0.0])
    se = np.array([0.1, 0.2])
    w = 0.04
    abf = wakefield_abf(beta, se, w)
    expected = np.sqrt(np.square(se) / (np.square(se) + w))
    assert np.allclose(abf, expected, atol=1e-12)


def test_coloc_posteriors_sum_to_one() -> None:
    abf1 = np.array([1.0, 2.0, 3.0])
    abf2 = np.array([1.5, 0.5, 2.0])
    pp = coloc_posteriors(abf1, abf2, 1e-4, 1e-4, 1e-5)
    total = pp["PP0"] + pp["PP1"] + pp["PP2"] + pp["PP3"] + pp["PP4"]
    assert abs(total - 1.0) < 1e-12
