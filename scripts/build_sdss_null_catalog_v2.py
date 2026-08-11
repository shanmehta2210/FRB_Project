"""
Build SDSS DR16 null catalog (v2) — full imaging footprint, HTM random sampling.

Unlike v1 (``build_sdss_null_catalog.py``), v2:
  - Covers the full SDSS PhotoObj footprint (no joint Legacy Dec clip).
  - Uses SDSS HTM-hash strata for unbiased random draws (not TOP-N per RA).
  - Stores ``objID`` for deduplication.
  - Iterates until the production strict null pool at ``modelMag_r < 20`` reaches
    ``--min-strict-mag20-pool`` (default 50k).

Run from repo root::

    python scripts/build_sdss_null_catalog_v2.py
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from astroquery.sdss import SDSS

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (
    SDSS_CATALOG_V2_DEFAULT,
    SDSS_LNL_COVERAGE_MIN,
    SDSS_MIN_STRICT_MAG20_POOL,
    assign_sdss_profile_winner_columns,
    count_strict_mag20_pool,
    footprint_summary,
    sdss_htm_stratum_clause,
    sdss_htm_stratum_edges,
)

DEFAULT_OUT = SDSS_CATALOG_V2_DEFAULT
DEFAULT_CACHE_DIR = "SDSS_v2_fetch_cache"
INTERNAL_POOL_MARGIN = 1.1  # stop at 110% of target for buffer
DEFAULT_MAX_ROWS = 4_000_000
DEFAULT_BATCH_SIZE = 100_000
DEFAULT_N_STRATA = 32


def _sql_htm_batch(stratum_lo: int, stratum_hi: int, top_n: int) -> str:
    htm_clause = sdss_htm_stratum_clause(stratum_lo, stratum_hi)
    return f"""
    SELECT TOP {int(top_n)}
        p.objID AS objID,
        p.ra AS ra,
        p.dec AS dec,
        p.cmodelMag_r AS cmodelMag_r,
        p.petroMag_r AS petroMag_r,
        p.modelMag_r AS modelMag_r,
        p.modelMag_u AS modelMag_u,
        p.modelMag_g AS modelMag_g,
        p.deVMag_r AS deVMag_r,
        p.expMag_r AS expMag_r,
        p.lnLDeV_r AS lnLDeV_r,
        p.lnLExp_r AS lnLExp_r,
        p.deVAB_r AS deVAB_r,
        p.expAB_r AS expAB_r,
        p.fracDeV_r AS fracDeV_r,
        p.expRad_r AS expRad_r,
        p.deVRad_r AS deVRad_r,
        p.type AS type
    FROM PhotoObj AS p
    WHERE p.type = 3
      AND p.clean = 1
      AND p.mode = 1
      AND p.cmodelMag_r > 0 AND p.cmodelMag_r < 90
      AND p.deVAB_r > 0 AND p.deVAB_r <= 1
      AND p.expAB_r > 0 AND p.expAB_r <= 1
      AND {htm_clause}
    """


def _normalize_objid_column(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure lowercase ``objid`` column for deduplication."""
    out = df.copy()
    for col in list(out.columns):
        if col.lower() == "objid":
            if col != "objid":
                out = out.rename(columns={col: "objid"})
            break
    else:
        raise KeyError("Query result missing objID column")
    return out


def query_sdss_htm_batch(
    stratum_lo: int,
    stratum_hi: int,
    batch_size: int,
    retries: int,
) -> pd.DataFrame:
    sql = _sql_htm_batch(stratum_lo, stratum_hi, batch_size)
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            tbl = SDSS.query_sql(sql, timeout=900)
            df = tbl.to_pandas()
            if len(df) == 0:
                return pd.DataFrame()
            return _normalize_objid_column(df)
        except Exception as exc:
            last_err = exc
            if attempt < retries:
                time.sleep(2.0 * attempt)
    raise RuntimeError(
        f"HTM stratum [{stratum_lo}, {stratum_hi}] failed after {retries} attempts: {last_err}"
    )


def dedupe_raw(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    out = _normalize_objid_column(df)
    return out.drop_duplicates(subset=["objid"], keep="first").reset_index(drop=True)


def build_catalog(df: pd.DataFrame) -> pd.DataFrame:
    objid = pd.to_numeric(df["objid"], errors="coerce")
    rmag = pd.to_numeric(df["cmodelMag_r"], errors="coerce")
    expab = pd.to_numeric(df["expAB_r"], errors="coerce")
    devab = pd.to_numeric(df["deVAB_r"], errors="coerce")
    model_mag = pd.to_numeric(df["modelMag_r"], errors="coerce")
    model_mag_u = pd.to_numeric(df.get("modelMag_u"), errors="coerce")
    model_mag_g = pd.to_numeric(df.get("modelMag_g"), errors="coerce")
    dev_mag = pd.to_numeric(df["deVMag_r"], errors="coerce")
    exp_mag = pd.to_numeric(df["expMag_r"], errors="coerce")
    frac_dev = pd.to_numeric(df["fracDeV_r"], errors="coerce")

    d_dev = np.abs(dev_mag - model_mag)
    d_exp = np.abs(exp_mag - model_mag)
    use_dev = d_dev <= d_exp
    best_model_ba = np.where(use_dev, devab, expab)
    exp_re = pd.to_numeric(df.get("expRad_r"), errors="coerce")
    dev_re = pd.to_numeric(df.get("deVRad_r"), errors="coerce")
    best_model_re = np.where(use_dev, dev_re, exp_re)
    n_eff = 1.0 + 3.0 * frac_dev

    ln_dev = pd.to_numeric(df.get("lnLDeV_r"), errors="coerce")
    ln_exp = pd.to_numeric(df.get("lnLExp_r"), errors="coerce")

    out = pd.DataFrame(
        {
            "objID": objid.astype(np.int64),
            "RA_ICRS": pd.to_numeric(df["ra"], errors="coerce"),
            "DE_ICRS": pd.to_numeric(df["dec"], errors="coerce"),
            "rmag": rmag,
            "petroMag_r": pd.to_numeric(df.get("petroMag_r"), errors="coerce"),
            "modelMag_r": model_mag,
            "modelMag_u": model_mag_u,
            "modelMag_g": model_mag_g,
            "deVMag_r": dev_mag,
            "expMag_r": exp_mag,
            "lnLDeV_r": ln_dev,
            "lnLExp_r": ln_exp,
            "u_r": model_mag_u - model_mag,
            "g_r": model_mag_g - model_mag,
            "expAB_r": expab,
            "deVAB_r": devab,
            "best_model_ba_r": best_model_ba,
            "fracDeV_r": frac_dev,
            "expRad_r": exp_re,
            "deVRad_r": dev_re,
            "best_model_re_r": best_model_re,
            "n_eff_r": n_eff,
            "b_a": expab,
            "sdss_type": pd.to_numeric(df.get("type"), errors="coerce"),
        }
    )

    good = (
        np.isfinite(out["objID"])
        & np.isfinite(out["RA_ICRS"])
        & np.isfinite(out["DE_ICRS"])
        & np.isfinite(out["rmag"])
        & np.isfinite(out["expAB_r"])
        & np.isfinite(out["best_model_ba_r"])
        & (out["expAB_r"] >= 0.0)
        & (out["expAB_r"] <= 1.0)
        & (out["best_model_ba_r"] >= 0.0)
        & (out["best_model_ba_r"] <= 1.0)
    )
    out = out.loc[good].copy().reset_index(drop=True)
    out["objID"] = out["objID"].astype(np.int64)
    if out["lnLExp_r"].notna().any():
        out = assign_sdss_profile_winner_columns(out)
    return out


def _merge_catalog_batches(
    catalog: pd.DataFrame | None,
    batch_raw: pd.DataFrame,
) -> pd.DataFrame:
    """Incrementally append a raw batch to the built catalog (objID dedupe)."""
    batch_cat = build_catalog(batch_raw)
    if catalog is None or catalog.empty:
        return batch_cat
    out = pd.concat([catalog, batch_cat], ignore_index=True)
    return out.drop_duplicates(subset=["objID"], keep="first").reset_index(drop=True)


def query_sdss_htm_batched(
    *,
    min_strict_mag20: int,
    batch_size: int,
    n_strata: int,
    max_rows: int,
    retries: int,
    cache_dir: Path | None,
    resume: bool,
    start_stratum: int = 0,
    initial_catalog: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, list[dict]]:
    """
    Fetch HTM-stratified batches until strict mag20 pool target is met.

    Returns (built_catalog_df, batch_log_rows).
    """
    strata = sdss_htm_stratum_edges(n_strata)
    stop_target = int(min_strict_mag20 * INTERNAL_POOL_MARGIN)
    batch_log: list[dict] = []
    catalog: pd.DataFrame | None = initial_catalog
    frames: list[pd.DataFrame] = []

    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)

    for idx, (lo, hi) in enumerate(strata):
        if idx < start_stratum:
            continue
        cache_path = cache_dir / f"stratum_{idx:03d}_{lo}_{hi}.csv" if cache_dir else None
        batch_df: pd.DataFrame | None = None

        if resume and cache_path is not None and cache_path.is_file():
            print(f"  [cache] stratum {idx}: loading {cache_path.name}", flush=True)
            batch_df = pd.read_csv(cache_path)

        if batch_df is None:
            print(f"  [query] stratum {idx}/{n_strata - 1} HTM [{lo}, {hi}]...", flush=True)
            try:
                batch_df = query_sdss_htm_batch(lo, hi, batch_size, retries)
            except RuntimeError as exc:
                print(f"  [!] {exc}", flush=True)
                batch_log.append(
                    {
                        "stratum": idx,
                        "htm_lo": lo,
                        "htm_hi": hi,
                        "n_fetched": 0,
                        "status": "failed",
                    }
                )
                continue
            if cache_path is not None and len(batch_df):
                batch_df.to_csv(cache_path, index=False)

        n_batch = len(batch_df) if batch_df is not None else 0
        if n_batch:
            frames.append(batch_df)
            catalog = _merge_catalog_batches(catalog, batch_df)

        n_unique = len(catalog) if catalog is not None else 0
        pool_n = count_strict_mag20_pool(catalog) if n_unique else 0

        batch_log.append(
            {
                "stratum": idx,
                "htm_lo": lo,
                "htm_hi": hi,
                "n_fetched": n_batch,
                "n_unique_objid": n_unique,
                "strict_mag20_pool": pool_n,
                "status": "ok",
            }
        )
        print(
            f"  stratum {idx}: +{n_batch} rows -> "
            f"unique={n_unique}, strict_mag20={pool_n} (target {stop_target})",
            flush=True,
        )

        if pool_n >= stop_target:
            print(f"[*] Strict mag20 pool {pool_n} >= {stop_target}; stopping fetch.", flush=True)
            break
        if n_unique >= max_rows:
            print(f"[!] Reached --max-rows {max_rows} before pool target.", flush=True)
            break

    if catalog is None or catalog.empty:
        raise RuntimeError("No rows returned from any HTM stratum.")

    return catalog, batch_log


def assemble_catalog_from_cache(cache_dir: Path) -> pd.DataFrame:
    """Build catalog from all cached stratum CSVs (sorted by name)."""
    paths = sorted(cache_dir.glob("stratum_*.csv"))
    if not paths:
        raise FileNotFoundError(f"No stratum cache files in {cache_dir}")
    frames = [pd.read_csv(p) for p in paths]
    raw = dedupe_raw(pd.concat(frames, ignore_index=True))
    return build_catalog(raw)


def _max_cached_stratum_index(cache_dir: Path) -> int:
    paths = sorted(cache_dir.glob("stratum_*.csv"))
    if not paths:
        return -1
    # filename: stratum_013_26624_28671.csv
    return max(int(p.name.split("_")[1]) for p in paths)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=DEFAULT_OUT, help="Output CSV path.")
    parser.add_argument(
        "--min-strict-mag20-pool",
        type=int,
        default=SDSS_MIN_STRICT_MAG20_POOL,
        help="Minimum production strict null pool at modelMag_r < 20.",
    )
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument("--n-strata", type=int, default=DEFAULT_N_STRATA)
    parser.add_argument("--max-rows", type=int, default=DEFAULT_MAX_ROWS)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--cache-dir",
        default=DEFAULT_CACHE_DIR,
        help="Cache per-stratum CSV batches (enables resume).",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Do not read/write stratum cache files.",
    )
    parser.add_argument(
        "--no-resume",
        action="store_true",
        help="Ignore existing cache files.",
    )
    parser.add_argument(
        "--assemble-only",
        action="store_true",
        help="Bulk-load cache, fetch only missing strata, then write (faster resume).",
    )
    parser.add_argument(
        "--start-stratum",
        type=int,
        default=None,
        help="First HTM stratum index to query (default: 0, or after last cache with --assemble-only).",
    )
    args = parser.parse_args()

    cache_dir = None if args.no_cache else Path(args.cache_dir)

    start_stratum = args.start_stratum if args.start_stratum is not None else 0
    initial_catalog = None

    if args.assemble_only and cache_dir is not None and cache_dir.is_dir():
        print(f"[*] Bulk-assembling from {cache_dir}...", flush=True)
        initial_catalog = assemble_catalog_from_cache(cache_dir)
        pool_pre = count_strict_mag20_pool(initial_catalog)
        print(
            f"[*] Cache assembly: N={len(initial_catalog)}, strict_mag20={pool_pre}",
            flush=True,
        )
        stop_target = int(args.min_strict_mag20_pool * INTERNAL_POOL_MARGIN)
        if pool_pre >= stop_target:
            catalog = initial_catalog.sample(frac=1.0, random_state=args.seed).reset_index(drop=True)
            footprint_summary(catalog, "SDSS v2 full")
            if catalog["objID"].duplicated().any():
                raise SystemExit("[!] Duplicate objID in final catalog.")
            Path(args.out).parent.mkdir(parents=True, exist_ok=True)
            catalog.to_csv(args.out, index=False)
            print(f"Wrote: {args.out} ({len(catalog)} rows, strict_mag20_pool={pool_pre})")
            return
        start_stratum = _max_cached_stratum_index(cache_dir) + 1
        print(f"[*] Continuing fetch from stratum {start_stratum}", flush=True)

    print(
        "Querying SDSS PhotoObj v2 (full footprint, HTM random strata, "
        f"type=GALAXY, clean=1, mode=1)...",
        flush=True,
    )

    catalog, batch_log = query_sdss_htm_batched(
        min_strict_mag20=args.min_strict_mag20_pool,
        batch_size=args.batch_size,
        n_strata=args.n_strata,
        max_rows=args.max_rows,
        retries=args.retries,
        cache_dir=cache_dir,
        resume=not args.no_resume,
        start_stratum=start_stratum,
        initial_catalog=initial_catalog,
    )
    print(f"Built catalog: {len(catalog)} unique objID rows")

    log_path = Path(args.out).with_suffix(".batch_log.csv")
    if log_path.is_file() and batch_log:
        old = pd.read_csv(log_path)
        batch_log_df = pd.concat([old, pd.DataFrame(batch_log)], ignore_index=True)
        batch_log_df = batch_log_df.drop_duplicates(subset=["stratum"], keep="last")
    else:
        batch_log_df = pd.DataFrame(batch_log)
    batch_log_df.to_csv(log_path, index=False)
    print(f"Wrote batch log: {log_path}")

    catalog = catalog.drop_duplicates(subset=["objID"], keep="first").reset_index(drop=True)

    catalog = catalog.sample(frac=1.0, random_state=args.seed).reset_index(drop=True)
    footprint_summary(catalog, "SDSS v2 full")

    pool_n = count_strict_mag20_pool(catalog)
    print(f"[*] Strict mag20 production pool: N={pool_n}")

    lnl_frac = float(pd.to_numeric(catalog["lnLExp_r"], errors="coerce").notna().mean())
    print(f"[*] lnL coverage: {lnl_frac:.1%}")
    if lnl_frac < SDSS_LNL_COVERAGE_MIN:
        print(
            f"[!] lnL coverage {lnl_frac:.1%} < {SDSS_LNL_COVERAGE_MIN:.0%}; "
            "run patch_sdss_profile_winner.py --footprint full --in-csv ... --out-csv ..."
        )

    if pool_n < args.min_strict_mag20_pool:
        raise SystemExit(
            f"[!] Strict mag20 pool {pool_n} < {args.min_strict_mag20_pool}. "
            "Increase --max-rows, --n-strata, or --batch-size and re-run."
        )

    if catalog["objID"].duplicated().any():
        raise SystemExit("[!] Duplicate objID in final catalog.")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    catalog.to_csv(args.out, index=False)
    print(f"Wrote: {args.out} ({len(catalog)} rows, strict_mag20_pool={pool_n})")


if __name__ == "__main__":
    main()
