"""Shared helpers for the PSF b/a vs magnitude simulation."""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Any, Iterator

import numpy as np
import yaml

TOOL_DIR = Path(__file__).resolve().parent
REPO_ROOT = TOOL_DIR.parents[1]


def load_config(config_path: Path | None = None) -> dict[str, Any]:
    path = config_path or (TOOL_DIR / "config.yaml")
    with open(path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        raise ValueError(f"Invalid config: {path}")
    return cfg


def outputs_dir(cfg: dict[str, Any]) -> Path:
    rel = cfg.get("paths", {}).get("outputs_dir", "outputs")
    return (TOOL_DIR / rel).resolve()


def mag_grid(cfg: dict[str, Any]) -> list[float]:
    g = cfg["grid"]
    start = float(g["mag_start"])
    stop = float(g["mag_stop"])
    step = float(g["mag_step"])
    n = int(round((stop - start) / step)) + 1
    return [round(start + i * step, 6) for i in range(n)]


def realizations_for_mag(mag: float, cfg: dict[str, Any]) -> int:
    """Return number of noise realizations for a given magnitude."""
    r_cfg = cfg.get("realizations", {})
    if r_cfg.get("n_realizations") is not None:
        return max(1, int(r_cfg["n_realizations"]))
    by_mag = r_cfg.get("n_realizations_by_mag") or {"default": 1}
    default = int(by_mag.get("default", 1))
    best_threshold: float | None = None
    best_n = default
    for key, n in by_mag.items():
        if key == "default":
            continue
        if str(key).startswith(">="):
            thresh = float(str(key)[2:])
            if mag >= thresh and (best_threshold is None or thresh >= best_threshold):
                best_threshold = thresh
                best_n = int(n)
    return max(1, best_n)


def _fmt_float(val: float, prefix: str = "") -> str:
    s = f"{val:.1f}".replace(".", "p")
    return f"{prefix}{s}"


def make_galaxy_id(ba: float, re_arcsec: float, mag: float, realization: int) -> str:
    ba_tag = f"ba{int(round(ba * 100)):03d}"
    re_tag = _fmt_float(re_arcsec, "re")
    mag_tag = _fmt_float(mag, "m")
    return f"{ba_tag}_{re_tag}_{mag_tag}_r{realization:03d}"


def realization_seed(base_seed: int, galaxy_id: str, realization: int) -> int:
    digest = hashlib.sha256(f"{galaxy_id}:{realization}".encode()).hexdigest()
    return (base_seed + int(digest[:8], 16)) % (2**31 - 1)


def flux_e_from_mag(mag: float, zeropoint_e: float) -> float:
    return float(10.0 ** (-0.4 * (mag - zeropoint_e)))


def sky_e_per_pix(sky_mag_arcsec2: float, zeropoint_e: float, pixel_scale: float) -> float:
    sky_flux_per_arcsec2 = 10.0 ** (-0.4 * (sky_mag_arcsec2 - zeropoint_e))
    return float(sky_flux_per_arcsec2 * pixel_scale**2)


def galfit_center_1based(stamp_px: int) -> tuple[float, float]:
    """1-based GALFIT coordinates for GalSim-centered even-sized stamp."""
    c = stamp_px / 2.0 + 0.5
    return c, c


def re_pix_from_arcsec(re_arcsec: float, pixel_scale: float) -> float:
    return float(re_arcsec / pixel_scale)


def iter_grid_points(cfg: dict[str, Any]) -> Iterator[dict[str, Any]]:
    g = cfg["grid"]
    for ba in g["intrinsic_ba"]:
        for re_arcsec in g["intrinsic_re_arcsec"]:
            for mag in mag_grid(cfg):
                n_real = realizations_for_mag(mag, cfg)
                for r in range(n_real):
                    gid = make_galaxy_id(float(ba), float(re_arcsec), mag, r)
                    yield {
                        "galaxy_id": gid,
                        "ba_true": float(ba),
                        "re_arcsec_true": float(re_arcsec),
                        "mag_true": mag,
                        "realization": r,
                        "pa_true": float(g["pa_deg"]),
                        "sersic_n": float(g["sersic_n"]),
                    }


def ensure_output_layout(cfg: dict[str, Any]) -> dict[str, Path]:
    root = outputs_dir(cfg)
    layout = {
        "root": root,
        "mocks": root / "mocks",
        "fits": root / "fits",
        "catalogs": root / "catalogs",
        "plots": root / "plots",
    }
    for p in layout.values():
        p.mkdir(parents=True, exist_ok=True)
    return layout


def resolve_zeropoint_e(cfg: dict[str, Any]) -> float:
    zp = cfg.get("physics", {}).get("zeropoint_e")
    if zp is None:
        raise ValueError("physics.zeropoint_e must be set (fixed ZP; use coadd_exptime_sec for SNR tuning)")
    return float(zp)


def resolve_coadd_exptime(cfg: dict[str, Any], calibrate_fn) -> float:
    t = cfg.get("physics", {}).get("coadd_exptime_sec")
    if t is not None and math.isfinite(float(t)):
        return float(t)
    return float(calibrate_fn(cfg))
