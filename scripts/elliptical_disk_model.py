"""
Ryden (2004) / Padilla & Strauss (2008) elliptical-disk generative model.

Thickness gamma = C/A ~ truncated Gaussian (mu_g, sig_g).
Face-on ellipticity eps = 1 - B/A ~ lognormal in ln(eps) (mu, sig).
Projected axis ratio q from Binney (1985) / Ryden eqs. (12)-(15).

Padilla dust (optional): E(theta) with edge-on E0; sample cos(theta) with
weight psi ~ 10^{-0.4 E} (photometric adaptation; no 1/Vmax).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Optional

import numpy as np


@dataclass
class ShapeParams:
    mu_g: float  # mean thickness gamma = C/A
    sig_g: float
    mu: float  # mean ln(eps)
    sig: float
    E0: float = 0.0  # Padilla edge-on extinction (mag); 0 = Ryden

    def to_dict(self) -> dict:
        return asdict(self)


# Literature seeds
RYDEN_SEED = ShapeParams(mu_g=0.222, sig_g=0.057, mu=-1.85, sig=0.89, E0=0.0)
PADILLA_SEED = ShapeParams(mu_g=0.21, sig_g=0.05, mu=-2.33, sig=0.79, E0=0.45)


@dataclass
class CappedShapeParams:
    """Elliptical-disk shapes whose intrinsic FACE-ON b/a is capped near ``cap``.

    Motivation: Legacy Survey REX (round-object) typing removes disks with observed
    b/a above ~0.8 from the EXP sample (see REX_AND_ELLIPTICAL_DISK.md). Because a
    disk's projected b/a can never exceed its face-on b/a (= 1 - eps), we can encode
    that ceiling *physically* -- as a property of the intrinsic shape distribution --
    instead of via the ad-hoc post-hoc ``/0.8`` renormalisation used by ``scaled``.

    The intrinsic face-on axis ratio ``f = 1 - eps`` is drawn as:
        * uniform on ``[0, cap]``            (flat "otherwise")
        * exp(-(f - cap)/lam) on ``(cap, 1]``  (controllable fall-off past the cap)
    so ``lam -> 0`` is a hard cap at ``cap`` and larger ``lam`` softens it. Thickness
    ``gamma = C/A ~ N(mu_g, sig_g)`` (truncated) and orientation isotropic, exactly as
    Ryden (2004). This is a genuine generative Ryden distribution (real cos i spread),
    not the degenerate delta that ``scaled`` corresponds to.
    """

    mu_g: float          # mean thickness gamma = C/A
    sig_g: float
    cap: float = 0.8     # face-on b/a ceiling (REX edge)
    lam: float = 0.03    # exponential fall-off scale above the cap (0 => hard cap)
    E0: float = 0.0      # Padilla edge-on extinction; 0 = no dust in the shape model

    def to_dict(self) -> dict:
        return asdict(self)


def sample_faceon_ba_capped(
    n: int, cap: float, lam: float, rng: np.random.Generator
) -> np.ndarray:
    """Face-on b/a: uniform on [0, cap], exp tail exp(-(f-cap)/lam) on (cap, 1]."""
    cap = float(np.clip(cap, 1e-3, 1.0 - 1e-6))
    lam = float(max(lam, 1e-9))
    area_flat = cap
    denom = 1.0 - np.exp(-(1.0 - cap) / lam)  # tail normalisation on (cap, 1]
    area_tail = lam * denom
    p_flat = area_flat / (area_flat + area_tail)

    u = rng.uniform(0.0, 1.0, size=n)
    f = np.empty(n, dtype=float)
    flat = u < p_flat
    nflat = int(flat.sum())
    if nflat:
        f[flat] = rng.uniform(0.0, cap, size=nflat)
    ntail = n - nflat
    if ntail:
        v = rng.uniform(0.0, 1.0, size=ntail)
        # inverse-CDF of the truncated exponential on (cap, 1]
        x = cap - lam * np.log(1.0 - v * denom)
        f[~flat] = np.clip(x, cap, 1.0)
    return np.clip(f, 1e-4, 1.0 - 1e-9)


def generate_q_capped(
    params: CappedShapeParams,
    n: int,
    rng: Optional[np.random.Generator] = None,
    *,
    return_angles: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Same projection as ``generate_q`` but with the capped face-on b/a law."""
    rng = rng or np.random.default_rng()
    gamma = sample_gamma(n, params.mu_g, params.sig_g, rng)
    eps = 1.0 - sample_faceon_ba_capped(n, params.cap, params.lam, rng)
    cos_t = sample_cos_theta(n, gamma, params.E0, rng)
    theta = np.arccos(np.clip(cos_t, 0.0, 1.0))
    phi = rng.uniform(0.0, 2.0 * np.pi, size=n)
    q = projected_q(gamma, eps, theta, phi)
    if return_angles:
        return q, cos_t, theta
    return q


def sample_gamma(n: int, mu_g: float, sig_g: float, rng: np.random.Generator) -> np.ndarray:
    g = rng.normal(mu_g, max(sig_g, 1e-6), size=n)
    # reject/resample out of [0, 1]
    for _ in range(8):
        bad = (g < 0.0) | (g > 1.0)
        if not np.any(bad):
            break
        g[bad] = rng.normal(mu_g, max(sig_g, 1e-6), size=int(bad.sum()))
    return np.clip(g, 1e-4, 1.0 - 1e-6)


def sample_eps(n: int, mu: float, sig: float, rng: np.random.Generator) -> np.ndarray:
    # ln eps ~ N(mu, sig^2), eps in (0, 1)
    ln_e = rng.normal(mu, max(sig, 1e-6), size=n)
    eps = np.exp(ln_e)
    for _ in range(8):
        bad = (eps <= 0.0) | (eps >= 1.0)
        if not np.any(bad):
            break
        ln_e[bad] = rng.normal(mu, max(sig, 1e-6), size=int(bad.sum()))
        eps[bad] = np.exp(ln_e[bad])
    return np.clip(eps, 1e-6, 1.0 - 1e-6)


def projected_q(
    gamma: np.ndarray,
    eps: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
) -> np.ndarray:
    """Ryden (2004) eqs. 12-15. theta = inclination from face-on (0=face-on)."""
    ct = np.cos(theta)
    st = np.sin(theta)
    cp = np.cos(phi)
    sp = np.sin(phi)
    e2 = eps * (2.0 - eps)
    A = (1.0 - e2 * sp * sp) * ct * ct + gamma * gamma * st * st
    B = 4.0 * (eps * (2.0 - eps)) ** 2 * ct * ct * sp * sp * cp * cp
    C = 1.0 - e2 * cp * cp
    disc = (A - C) ** 2 + B
    disc = np.maximum(disc, 0.0)
    root = np.sqrt(disc)
    num = A + C - root
    den = A + C + root
    q2 = np.where(den > 0, num / den, 0.0)
    q2 = np.clip(q2, 0.0, 1.0)
    return np.sqrt(q2)


def padilla_E(cos_theta: np.ndarray, gamma: np.ndarray, E0: float) -> np.ndarray:
    """Padilla E(theta); use Ryden gamma≈C/A as path-length proxy y."""
    if E0 <= 0.0:
        return np.zeros_like(cos_theta)
    y = np.asarray(gamma, dtype=float)
    ct = np.asarray(cos_theta, dtype=float)
    E = np.where(ct > y, E0 * (1.0 + y - ct), E0)
    return np.maximum(E, 0.0)


def sample_cos_theta(
    n: int,
    gamma: np.ndarray,
    E0: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Uniform in cos(theta) if E0=0; else reject with psi ~ 10^{-0.4 E(gamma,theta)}."""
    if E0 <= 0.0:
        return rng.uniform(0.0, 1.0, size=n)

    out = np.empty(n, dtype=float)
    todo = np.ones(n, dtype=bool)
    # Cap iterations; leftover get uniform (should be rare)
    for _ in range(80):
        if not np.any(todo):
            break
        idx = np.flatnonzero(todo)
        k = len(idx)
        ct = rng.uniform(0.0, 1.0, size=k)
        E = padilla_E(ct, gamma[idx], E0)
        psi = 10.0 ** (-0.4 * E)
        accept = rng.uniform(0.0, 1.0, size=k) < psi
        out[idx[accept]] = ct[accept]
        todo[idx[accept]] = False
    if np.any(todo):
        out[todo] = rng.uniform(0.0, 1.0, size=int(todo.sum()))
    return out


def generate_q(
    params: ShapeParams,
    n: int,
    rng: Optional[np.random.Generator] = None,
    *,
    return_angles: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = rng or np.random.default_rng()
    gamma = sample_gamma(n, params.mu_g, params.sig_g, rng)
    eps = sample_eps(n, params.mu, params.sig, rng)
    cos_t = sample_cos_theta(n, gamma, params.E0, rng)
    theta = np.arccos(np.clip(cos_t, 0.0, 1.0))
    phi = rng.uniform(0.0, 2.0 * np.pi, size=n)
    q = projected_q(gamma, eps, theta, phi)
    if return_angles:
        return q, cos_t, theta
    return q


def histogram_q(
    q: np.ndarray,
    *,
    bin_width: float = 0.05,
    q_min: float = 0.0,
    q_max: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    edges = np.arange(q_min, q_max + 0.5 * bin_width, bin_width)
    if edges[-1] < q_max - 1e-12:
        edges = np.append(edges, q_max)
    counts, edges = np.histogram(q, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return counts.astype(float), edges, centers


def model_histogram(
    params: ShapeParams,
    n_model: int,
    data_counts: np.ndarray,
    edges: np.ndarray,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    q = generate_q(params, n_model, rng)
    model_counts, _ = np.histogram(q, bins=edges)
    model_counts = model_counts.astype(float)
    # scale to same total N as data
    n_data = float(np.sum(data_counts))
    n_mod = float(np.sum(model_counts))
    if n_mod > 0:
        model_counts *= n_data / n_mod
    return model_counts


def chi2_hist(
    data_counts: np.ndarray,
    model_counts: np.ndarray,
    *,
    min_count: float = 10.0,
) -> float:
    ok = data_counts >= min_count
    if not np.any(ok):
        return np.inf
    d = data_counts[ok]
    m = model_counts[ok]
    # Poisson-like variance from data
    var = np.maximum(d, 1.0)
    return float(np.sum((m - d) ** 2 / var))


class ConditionalCosiSampler:
    """Draw cos(i) from the model posterior P(cos i | q).

    A given projected b/a maps to a *distribution* of inclinations, not a single
    value (that is the whole point of an elliptical-disk model). We forward-model
    a large pool of (q, cos theta), bin it finely in q, and for each observed
    b/a draw a random model cos(theta) from the matching q-bin. This restores the
    full conditional spread; using the per-bin median instead collapses every
    galaxy onto a narrow band and yields a degenerate near-vertical "step" CDF.
    """

    def __init__(
        self,
        params: Optional[ShapeParams] = None,
        rng: Optional[np.random.Generator] = None,
        *,
        n_model: int = 2_000_000,
        n_bins: int = 100,
        pool: Optional[tuple[np.ndarray, np.ndarray]] = None,
    ) -> None:
        if pool is not None:
            qm, cm = pool  # (projected b/a, cos theta) from any generator
        else:
            rng = rng or np.random.default_rng(0)
            qm, cm, _ = generate_q(params, n_model, rng, return_angles=True)
        self._build(np.asarray(qm, float), np.asarray(cm, float), n_bins)

    @classmethod
    def from_pool(
        cls, qm: np.ndarray, cm: np.ndarray, *, n_bins: int = 100
    ) -> "ConditionalCosiSampler":
        """Build directly from a precomputed (b/a, cos theta) pool (e.g. capped model)."""
        return cls(pool=(qm, cm), n_bins=n_bins)

    def _build(self, qm: np.ndarray, cm: np.ndarray, n_bins: int) -> None:
        self.edges = np.linspace(0.0, 1.0, n_bins + 1)
        self.n_bins = int(n_bins)
        bm = np.clip(np.digitize(qm, self.edges) - 1, 0, self.n_bins - 1)
        order = np.argsort(bm, kind="stable")
        self._bm = bm[order]
        self._cm = cm[order]
        self._starts = np.searchsorted(self._bm, np.arange(self.n_bins), side="left")
        self._ends = np.searchsorted(self._bm, np.arange(self.n_bins), side="right")
        self._pop = np.flatnonzero(self._ends > self._starts)

    def sample(self, q_obs: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        q_obs = np.asarray(q_obs, dtype=float)
        out = np.full(q_obs.shape, np.nan, dtype=float)
        if q_obs.size == 0 or self._pop.size == 0:
            return out
        bo = np.clip(np.digitize(q_obs, self.edges) - 1, 0, self.n_bins - 1)
        for b in np.unique(bo):
            sel = np.flatnonzero(bo == b)
            lo, hi = int(self._starts[b]), int(self._ends[b])
            if hi <= lo:  # empty model bin: fall back to nearest populated bin
                b_use = int(self._pop[np.argmin(np.abs(self._pop - b))])
                lo, hi = int(self._starts[b_use]), int(self._ends[b_use])
            pick = rng.integers(lo, hi, size=sel.size)
            out[sel] = self._cm[pick]
        return out


def inclination_vs_q(
    params: ShapeParams,
    n: int = 400_000,
    bin_width: float = 0.05,
    rng: Optional[np.random.Generator] = None,
) -> dict[str, np.ndarray]:
    """Median / 16-84% of cos(theta) in bins of projected q."""
    rng = rng or np.random.default_rng(0)
    q, cos_t, _ = generate_q(params, n, rng, return_angles=True)
    edges = np.arange(0.0, 1.0 + 0.5 * bin_width, bin_width)
    centers = 0.5 * (edges[:-1] + edges[1:])
    med = np.full(len(centers), np.nan)
    lo = np.full(len(centers), np.nan)
    hi = np.full(len(centers), np.nan)
    for i in range(len(centers)):
        mask = (q >= edges[i]) & (q < edges[i + 1])
        if i == len(centers) - 1:
            mask = (q >= edges[i]) & (q <= edges[i + 1])
        if np.sum(mask) < 30:
            continue
        c = cos_t[mask]
        med[i] = float(np.median(c))
        lo[i] = float(np.percentile(c, 16))
        hi[i] = float(np.percentile(c, 84))
    return {"q": centers, "cosi_med": med, "cosi_lo": lo, "cosi_hi": hi}
