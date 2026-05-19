"""
Parse GALFIT fit.log files.

Used by build_master_frb_galfit_from_logs.py and can replace brittle split()
logic in compile_galfit_logs.py (GALFIT writes full-width dash rules, not the
short '--------------------------------------------' substring).

Selection policy:
  - Split the log on lines that are only dashes (≥20).
  - Keep blocks that contain Chi^2/nu and at least one sersic component.
  - Prefer single-sersic blocks with sane Chi^2/nu (< 1e6) to skip numerical blow-ups.
  - If the last sane block fixes n (bracketed n on the sersic line), use the
    previous sane block for free-n summaries.
  - Otherwise use the last sane single-sersic block.
  - Multi-sersic fits: last qualifying block, then pick the Nth `sersic` line
    (default N=0).  For this repo's pipeline, `host_components.csv` is ordered
    with the FRB host first and `run_galfit_fitting.py` writes that row as
    GALFIT component 1, so N=0 is the host; neighbors are N=1, 2, ...
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import numpy as np

Q0_DEFAULT = 0.2

_DASH_SPLIT = re.compile(r"^\s*-{20,}\s*$", re.MULTILINE)
_SERSIC_LINE = re.compile(r"^\s*sersic\s*:", re.MULTILINE)


def split_fitlog_blocks(content: str) -> list[str]:
    parts = _DASH_SPLIT.split(content)
    return [p.strip() for p in parts if p.strip()]


def _count_sersic_lines(block: str) -> int:
    return len(_SERSIC_LINE.findall(block))


def _qualifying_blocks(content: str) -> list[str]:
    out = []
    for b in split_fitlog_blocks(content):
        if "Chi^2/nu" in b and "sersic" in b.lower():
            out.append(b)
    return out


def _chi2nu_from_block(block: str) -> float:
    for line in block.split("\n"):
        if "Chi^2/nu =" in line:
            try:
                return float(line.split("=")[1].split(",")[0].strip())
            except (ValueError, IndexError):
                return float("nan")
    return float("nan")


def _is_fixed_n_sersic_line(line: str) -> bool:
    """True if GALFIT printed n as constrained, e.g. [1.00] or [6.00]."""
    if not _SERSIC_LINE.match(line):
        return False
    return "[" in line and "]" in line


def _block_has_fixed_n(block: str) -> bool:
    for line in block.split("\n"):
        if _is_fixed_n_sersic_line(line):
            return True
    return False


def _select_block(blocks: list[str]) -> tuple[str, str]:
    """
    Returns (block, strategy) where strategy documents how the block was chosen.

    Policy:
      1. Among single-sersic iteration blocks, keep those with finite Chi^2/nu < 1e6
         (drops repeated numerical blow-ups).
      2. If the last sane block looks like a fixed-n tweak (bracketed n on the sersic line),
         use the previous sane block — matches free-n summaries when the final iteration
         fixes n (cf. progress.md).
      3. Otherwise use the last sane single-sersic block.
      4. If no sane blocks, fall back to legacy second-to-last / last among singles.
      5. Multi-sersic fits: last qualifying block, first sersic component.
    """
    if not blocks:
        return "", "empty"

    singles = [b for b in blocks if _count_sersic_lines(b) == 1]

    sane_singles: list[str] = []
    for b in singles:
        c2 = _chi2nu_from_block(b)
        if math.isfinite(c2) and abs(c2) < 1e6:
            sane_singles.append(b)

    if sane_singles:
        last_s = sane_singles[-1]
        if len(sane_singles) >= 2 and _block_has_fixed_n(last_s):
            return sane_singles[-2], "free_n_before_fixed_n_refine"
        return last_s, "last_sane_single_sersic"

    # Degenerate / failed logs: avoid picking a blow-up when middle iteration was OK
    if len(singles) >= 2:
        return singles[-2], "second_to_last_single_sersic_fallback"
    if len(singles) == 1:
        return singles[-1], "only_single_sersic_block"

    return blocks[-1], "last_block_first_sersic"


def _parse_float(tok: str) -> float:
    t = str(tok).replace("*", "").replace("[", "").replace("]", "").strip()
    return float(t)


def _parse_sersic_line_pair(line: str, err_line: str | None) -> dict[str, Any]:
    """Parse one `sersic : (...)` line and the following error line."""
    clean = (
        line.strip()
        .replace("(", " ")
        .replace(")", " ")
        .replace(",", " ")
        .replace("*", " ")
        .replace("[", " ")
        .replace("]", " ")
    )
    parts = clean.split()
    data: dict[str, Any] = {}
    if len(parts) < 9 or parts[0].strip() != "sersic":
        return data
    try:
        data["x"] = _parse_float(parts[2])
        data["y"] = _parse_float(parts[3])
        data["mag"] = _parse_float(parts[4])
        data["re"] = _parse_float(parts[5])
        data["n"] = _parse_float(parts[6])
        data["b_a"] = _parse_float(parts[7])
        data["pa"] = _parse_float(parts[8])
    except (ValueError, IndexError):
        return {}

    if err_line is None:
        return data
    if "sky" in err_line.lower():
        return data

    eclean = (
        err_line.replace("(", " ")
        .replace(")", " ")
        .replace(",", " ")
        .replace("*", " ")
        .replace("[", " ")
        .replace("]", " ")
    )
    ep = eclean.split()
    if len(ep) >= 7:
        try:
            data["x_err"] = _parse_float(ep[0])
            data["y_err"] = _parse_float(ep[1])
            data["mag_err"] = _parse_float(ep[2])
            data["re_err"] = _parse_float(ep[3])
            data["n_err"] = _parse_float(ep[4])
            data["b_a_err"] = _parse_float(ep[5])
            data["pa_err"] = _parse_float(ep[6])
        except (ValueError, IndexError):
            pass
    return data


def parse_fitlog_block(
    content: str, sersic_component_index: int = 0
) -> tuple[dict[str, Any], str]:
    """
    Parse selected iteration block from full fit.log text.
    Returns (dict of floats, strategy string).

    ``sersic_component_index`` selects which ``sersic :`` line to read when
    multiple components are present (0 = first = FRB host in our pipeline).
    """
    blocks = _qualifying_blocks(content)
    block, strategy = _select_block(blocks)
    if not block:
        return {}, strategy

    lines = block.split("\n")
    chi2nu: float | None = None
    for line in lines:
        if "Chi^2/nu =" in line:
            try:
                chi2nu = float(line.split("=")[1].split(",")[0].strip())
            except (ValueError, IndexError):
                pass
            break

    # Indices of sersic lines
    sidx = [i for i, ln in enumerate(lines) if _SERSIC_LINE.match(ln)]
    if not sidx:
        return {}, strategy

    k = int(sersic_component_index)
    if k < 0 or k >= len(sidx):
        k = 0
    i0 = sidx[k]
    err_line = None
    if i0 + 1 < len(lines):
        nxt = lines[i0 + 1]
        if "sky" not in nxt.lower():
            err_line = nxt

    data = _parse_sersic_line_pair(lines[i0], err_line)
    if chi2nu is not None:
        data["chi2nu"] = chi2nu
    return data, strategy


def parse_fitlog_file(
    path: str, sersic_component_index: int = 0
) -> tuple[dict[str, Any], str]:
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            content = f.read()
    except OSError:
        return {}, "missing_or_unreadable"

    if not content.strip():
        return {}, "empty_file"

    return parse_fitlog_block(content, sersic_component_index=sersic_component_index)


def inclination_from_b_a(b_a: float | None, q0: float = Q0_DEFAULT) -> float | None:
    if b_a is None or not math.isfinite(b_a):
        return None
    q = float(b_a)
    if q <= q0:
        return 90.0
    val = (q * q - q0 * q0) / (1.0 - q0 * q0)
    val = min(1.0, max(0.0, val))
    return round(float(np.degrees(np.arccos(np.sqrt(val)))), 4)


def inclination_err_from_b_a_err(
    b_a: float | None, b_a_err: float | None, q0: float = Q0_DEFAULT
) -> float | None:
    if b_a is None or not math.isfinite(b_a):
        return None
    if b_a_err is None or not math.isfinite(b_a_err) or b_a_err <= 0:
        return 0.0
    ba, e = float(b_a), float(b_a_err)
    q_up = min(1.0, ba + e)
    q_down = max(1e-9, ba - e)
    i_up = inclination_from_b_a(q_up, q0)
    i_down = inclination_from_b_a(q_down, q0)
    if i_up is None or i_down is None:
        return None
    err = abs(i_up - i_down) / 2.0
    if err == 0.0 and ba <= q0:
        i_boundary = inclination_from_b_a(min(1.0, q0 + e), q0)
        if i_boundary is not None:
            err = abs(90.0 - i_boundary)
    return round(float(err), 4)


# Backwards-compatible name for older scripts -------------------------------------------------

_LEGACY_LOG_KEYS = (
    "x",
    "y",
    "mag",
    "re",
    "n",
    "b_a",
    "pa",
    "chi2nu",
    "x_err",
    "y_err",
    "mag_err",
    "re_err",
    "n_err",
    "b_a_err",
    "pa_err",
)


_SKY_LEVEL_RE = re.compile(
    r"sky\s*:.*?\[.*?\]\s*([-\d.eE+]+)", re.IGNORECASE
)


def count_fitted_sersic_components(
    output_dir: str | Path,
    *,
    log_path: str | Path | None = None,
) -> int | None:
    """
    Number of Sérsic components in the pipeline fit (sky excluded).

    Prefer ``host_components.csv`` row count (Phase 3a fit list). Fall back to
    the last qualifying ``fit.log`` block's ``sersic`` line count.
    """
    odir = Path(output_dir)
    comp = odir / "host_components.csv"
    if comp.is_file():
        try:
            import pandas as pd

            n = len(pd.read_csv(comp))
            return int(n) if n > 0 else None
        except Exception:
            pass
    path = Path(log_path) if log_path is not None else odir / "fit.log"
    if not path.is_file():
        return None
    blocks = _qualifying_blocks(path.read_text(encoding="utf-8", errors="replace"))
    if not blocks:
        return None
    n = _count_sersic_lines(blocks[-1])
    return int(n) if n > 0 else None


def _sky_from_block(block: str) -> float | None:
    match = _SKY_LEVEL_RE.search(block)
    if not match:
        return None
    try:
        val = float(match.group(1))
    except ValueError:
        return None
    if not math.isfinite(val):
        return None
    return val


def _last_sky_line_value(content: str) -> float | None:
    """Last ``sky : [...] value`` line anywhere in the log (iteration or summary)."""
    last: float | None = None
    for line in content.splitlines():
        match = _SKY_LEVEL_RE.search(line)
        if not match:
            continue
        try:
            val = float(match.group(1))
        except ValueError:
            continue
        if math.isfinite(val):
            last = val
    return last


def sky_level_is_plausible(
    sky_adu: float | None,
    sky_ref: float | None = None,
    *,
    max_chi2nu: float = 1e6,
) -> bool:
    """Reject numerical blow-ups and skies far from the SExtractor seed."""
    if sky_adu is None or not math.isfinite(sky_adu):
        return False
    if sky_ref is not None and math.isfinite(sky_ref):
        ref = abs(float(sky_ref))
        # Typical pipeline backgrounds are ~1e-4 ADU; allow generous headroom.
        limit = max(1.0, ref * 1.0e4 + 0.05)
        if abs(sky_adu) > limit:
            return False
    elif abs(sky_adu) > 1.0:
        return False
    return True


def parse_fitlog_sky_level(
    log_path: str | Path,
    *,
    sky_ref: float | None = None,
    max_chi2nu: float = 1e6,
) -> float | None:
    """
    Fitted global sky level (ADU) from fit.log.

    Prefer the last qualifying block with sane Chi^2/nu (< ``max_chi2nu``).
    If GALFIT crashed before writing a summary block, fall back to the last
    ``sky :`` line in the file (iteration output). Values that fail
    ``sky_level_is_plausible`` are discarded.
    """
    path = Path(log_path) if not isinstance(log_path, Path) else log_path
    if not path.is_file():
        return None
    content = path.read_text(encoding="utf-8", errors="replace")
    blocks = _qualifying_blocks(content)
    for block in reversed(blocks):
        chi2 = _chi2nu_from_block(block)
        if not math.isfinite(chi2) or abs(chi2) >= max_chi2nu:
            continue
        val = _sky_from_block(block)
        if sky_level_is_plausible(val, sky_ref, max_chi2nu=max_chi2nu):
            return val

    val = _last_sky_line_value(content)
    if sky_level_is_plausible(val, sky_ref, max_chi2nu=max_chi2nu):
        return val
    return None


def parse_fitlog_full(
    log_path: str, sersic_component_index: int = 0
) -> dict[str, Any]:
    """
    Drop-in style replacement for compile_galfit_logs.parse_fitlog_full.
    Returns all legacy keys; missing values are ''.
    """
    data, _strategy = parse_fitlog_file(
        log_path, sersic_component_index=sersic_component_index
    )
    out: dict[str, Any] = {k: "" for k in _LEGACY_LOG_KEYS}
    for k in _LEGACY_LOG_KEYS:
        if k in data and data[k] is not None:
            out[k] = data[k]
    return out
