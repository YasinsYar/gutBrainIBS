from __future__ import annotations

import math
from typing import Iterable

import numpy as np
from scipy.special import erfc


def normal_two_sided_pvalue(z: np.ndarray | float) -> np.ndarray | float:
    return erfc(np.abs(z) / np.sqrt(2.0))


def wakefield_abf(beta: np.ndarray, se: np.ndarray, w: float) -> np.ndarray:
    v = np.square(se)
    z2 = np.square(beta / se)
    r = w / (v + w)
    log_abf = 0.5 * (np.log1p(-r) + r * z2)
    return np.exp(log_abf)


def coloc_posteriors(abf1: np.ndarray, abf2: np.ndarray, p1: float, p2: float, p12: float) -> dict[str, float]:
    abf1 = np.asarray(abf1, dtype=float)
    abf2 = np.asarray(abf2, dtype=float)
    s1 = float(np.sum(abf1))
    s2 = float(np.sum(abf2))
    s12 = float(np.sum(abf1 * abf2))

    h0 = 1.0
    h1 = p1 * s1
    h2 = p2 * s2
    h3 = p1 * p2 * max(s1 * s2 - s12, 0.0)
    h4 = p12 * s12

    total = h0 + h1 + h2 + h3 + h4
    if total <= 0:
        return {"PP0": 1.0, "PP1": 0.0, "PP2": 0.0, "PP3": 0.0, "PP4": 0.0}
    return {
        "PP0": h0 / total,
        "PP1": h1 / total,
        "PP2": h2 / total,
        "PP3": h3 / total,
        "PP4": h4 / total,
    }


def fdr_bh(pvals: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(pvals), dtype=float)
    n = len(p)
    if n == 0:
        return np.array([])
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0.0, 1.0)
    return out


def ivw_mr(beta_exp: np.ndarray, se_exp: np.ndarray, beta_out: np.ndarray, se_out: np.ndarray) -> tuple[float, float, float]:
    ratio = beta_out / beta_exp
    var_ratio = np.square(se_out / beta_exp) + np.square(beta_out * se_exp / np.square(beta_exp))
    w = 1.0 / np.clip(var_ratio, 1e-12, None)
    beta = float(np.sum(w * ratio) / np.sum(w))
    se = math.sqrt(float(1.0 / np.sum(w)))
    z = beta / se
    p = float(math.erfc(abs(z) / math.sqrt(2.0)))
    return beta, se, p


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    idx = np.argsort(values)
    v = values[idx]
    w = weights[idx]
    cdf = np.cumsum(w) / np.sum(w)
    return float(v[np.searchsorted(cdf, 0.5)])


def egger_mr(beta_exp: np.ndarray, se_exp: np.ndarray, beta_out: np.ndarray, se_out: np.ndarray) -> tuple[float, float, float, float]:
    x = beta_exp
    y = beta_out
    w = 1.0 / np.clip(np.square(se_out), 1e-12, None)
    x_bar = float(np.sum(w * x) / np.sum(w))
    y_bar = float(np.sum(w * y) / np.sum(w))
    sxx = float(np.sum(w * np.square(x - x_bar)))
    sxy = float(np.sum(w * (x - x_bar) * (y - y_bar)))
    if sxx <= 0:
        return 0.0, 1.0, 0.0, 1.0
    slope = sxy / sxx
    intercept = y_bar - slope * x_bar

    resid = y - (intercept + slope * x)
    dof = max(len(x) - 2, 1)
    sigma2 = float(np.sum(w * np.square(resid)) / dof)
    se_slope = math.sqrt(sigma2 / sxx)
    se_intercept = math.sqrt(sigma2 * (1.0 / np.sum(w) + x_bar * x_bar / sxx))

    p_slope = float(math.erfc(abs(slope / max(se_slope, 1e-12)) / math.sqrt(2.0)))
    p_intercept = float(math.erfc(abs(intercept / max(se_intercept, 1e-12)) / math.sqrt(2.0)))
    return slope, p_slope, intercept, p_intercept
