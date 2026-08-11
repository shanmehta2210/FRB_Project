"""
Build a random ~500k-row sample from the DES Y1 morphology catalog
(Tarsitano+2018), read directly from the public NCSA FITS files.

Why not TAP: the NOIRLAB Data Lab sync gateway 504s on morph queries and async
tile pulls are slow/flaky. The authoritative catalog is 6 FITS binary tables
(~18 GB total, 54,581,334 rows) on the DES public server. A full download at the
observed ~0.6 MB/s would take ~8 h, so instead we exploit the fixed-width FITS
binary-table layout and issue HTTP *range* requests for only the rows we sample.

Method: read each file's header (cheap, via fsspec) to get the row count, row
byte width (NAXIS1=335), the record dtype, and the data start offset (datLoc).
FITS binary data is big-endian, so we parse each range with a big-endian dtype.
We read many random contiguous row-blocks spread across all 6 files (rows are
ordered by tile, so random block starts sample random tiles across the
footprint), transferring only ~200 MB to net ~500k galaxies with fits.

ba_r = 1 - ELLIPTICITY_SERSIC_R   (epsilon = 1 - b/a; ALREADY calibrated)
mag_r = MAG_SERSIC_R              (ALREADY includes MAG_CAL; do NOT add again)
*_CAL_* columns store eta = truth_med - measured_med already applied (Appendix B)

Run from repo root::

    python scripts/build_des_y1_morph_sample.py                     # ~500k
    python scripts/build_des_y1_morph_sample.py --target-rows 3000  # quick test
"""

from __future__ import annotations

import argparse
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

BASE = "https://desdr-server.ncsa.illinois.edu/despublic/y1a1_files/morph_catalogs/"
FILES = [
    "Y1A1_morphological_000001.fits",
    "Y1A1_morphological_000002.fits",
    "Y1A1_morphological_000003.fits",
    "Y1A1_morphological_000004.fits",
    "Y1A1_morphological_000005.fits",
    "Y1A1_morphological_000006.fits",
]

KEEP = [
    "COADD_OBJECTS_ID", "TILENAME", "RA", "DEC", "SG", "FLAGS_BADREGION",
    "FIT_AVAILABLE_R", "SN_R", "MAG_SERSIC_R", "MAG_CAL_R", "RE_R",
    "N_SERSIC_R", "ELLIPTICITY_SERSIC_R", "ELLIPTICITY_SERSIC_CAL_R",
    "MAG_SERSIC_I", "ELLIPTICITY_SERSIC_I", "FIT_AVAILABLE_I", "SN_I",
    "MAG_SERSIC_G", "ELLIPTICITY_SERSIC_G",
]

# Defaults resolved relative to repo root (run from repo root).
DEFAULT_OUT = "catalog/DES_y1_morph_sample_500k.csv"
DEFAULT_CACHE = "DES_y1_morph_fetch_cache"
DEFAULT_TARGET = 500_000
DEFAULT_BLOCK = 4000
DEFAULT_OVERSAMPLE = 1.15
DEFAULT_RETRIES = 6
DEFAULT_SEED = 42
HTTP_TIMEOUT = 120


def get_meta(url: str, retries: int) -> tuple[int, int, np.dtype, int]:
    last = None
    for attempt in range(1, retries + 1):
        try:
            with fits.open(url, use_fsspec=True) as hdul:
                h = hdul[1].header
                nrows = int(h["NAXIS2"])
                rowbytes = int(h["NAXIS1"])
                dt_be = hdul[1].columns.dtype.newbyteorder(">")
                datloc = int(hdul.fileinfo(1)["datLoc"])
            return nrows, rowbytes, dt_be, datloc
        except Exception as exc:  # noqa: BLE001
            last = exc
            print(f"    [meta retry {attempt}/{retries}] {type(exc).__name__} {str(exc)[:90]}",
                  flush=True)
            time.sleep(min(30.0, 5.0 * attempt))
    raise RuntimeError(f"meta failed {url}: {last}")


def http_range(url: str, start: int, stop: int, retries: int) -> bytes:
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{stop-1}"})
            with urllib.request.urlopen(req, timeout=HTTP_TIMEOUT) as r:
                return r.read()
        except Exception as exc:  # noqa: BLE001
            last = exc
            print(f"      [range retry {attempt}/{retries}] {type(exc).__name__} {str(exc)[:80]}",
                  flush=True)
            time.sleep(min(30.0, 4.0 * attempt))
    raise RuntimeError(f"range [{start},{stop}) failed: {last}")


def block_to_df(buf: bytes, dt_be: np.dtype, nrows_block: int) -> pd.DataFrame:
    arr = np.frombuffer(buf, dtype=dt_be, count=nrows_block)
    data = {}
    for c in KEEP:
        col = arr[c]
        if col.dtype.kind == "S":
            col = np.char.decode(col, "ascii", "ignore")
            col = np.char.strip(col)
        else:
            # cast big-endian -> native for pandas
            col = col.astype(col.dtype.newbyteorder("="))
        data[c] = col
    return pd.DataFrame(data)


def build_catalog(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    fit = pd.to_numeric(out["FIT_AVAILABLE_R"], errors="coerce")
    e = pd.to_numeric(out["ELLIPTICITY_SERSIC_R"], errors="coerce")
    ba = 1.0 - e
    out["ba_r"] = ba.where((ba >= 0.0) & (ba <= 1.0))
    # Tarsitano+2018 Appendix B: MAG/RE/N/ELLIPTICITY_SERSIC already include
    # the eta corrections stored in *_CAL_*. Adding MAG_CAL again double-counts.
    mag = pd.to_numeric(out["MAG_SERSIC_R"], errors="coerce")
    out["mag_r"] = mag
    out["n_r"] = pd.to_numeric(out["N_SERSIC_R"], errors="coerce")
    out["re_r"] = pd.to_numeric(out["RE_R"], errors="coerce")
    out["RA_ICRS"] = pd.to_numeric(out["RA"], errors="coerce")
    out["DE_ICRS"] = pd.to_numeric(out["DEC"], errors="coerce")
    good = (
        (fit == 1)
        & out["ba_r"].notna()
        & out["mag_r"].notna()
        & out["RA_ICRS"].notna()
        & out["DE_ICRS"].notna()
        & (out["mag_r"] > 5) & (out["mag_r"] < 40)
    )
    return out.loc[good].reset_index(drop=True)


def plan_blocks(nrows: int, need_raw: int, block: int, rng) -> list[tuple[int, int]]:
    block = min(block, nrows)
    if nrows <= block:
        return [(0, nrows)]
    n_blocks = max(1, int(np.ceil(need_raw / block)))
    n_blocks = min(n_blocks, nrows // block)
    starts = rng.choice(nrows - block, size=n_blocks, replace=False)
    starts.sort()
    return [(int(s), int(s) + block) for s in starts]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--target-rows", type=int, default=DEFAULT_TARGET)
    p.add_argument("--block", type=int, default=DEFAULT_BLOCK)
    p.add_argument("--oversample", type=float, default=DEFAULT_OVERSAMPLE)
    p.add_argument("--retries", type=int, default=DEFAULT_RETRIES)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p.add_argument("--out", default=DEFAULT_OUT)
    p.add_argument("--cache-dir", default=DEFAULT_CACHE)
    p.add_argument("--no-cache", action="store_true")
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    cache_dir = None if args.no_cache else Path(args.cache_dir)
    if cache_dir:
        cache_dir.mkdir(parents=True, exist_ok=True)

    print("[*] Reading FITS headers ...", flush=True)
    meta = {}
    for fn in FILES:
        m = get_meta(BASE + fn, retries=args.retries)
        meta[fn] = m
        print(f"    {fn}: nrows={m[0]:,} rowbytes={m[1]} datLoc={m[3]}", flush=True)
    total = sum(m[0] for m in meta.values())
    print(f"[*] total morph rows={total:,}", flush=True)

    need_total_raw = int(args.target_rows * args.oversample)
    frames: list[pd.DataFrame] = []
    n_kept = 0

    for fn in FILES:
        if n_kept >= args.target_rows:
            break
        nrows, rowbytes, dt_be, datloc = meta[fn]
        url = BASE + fn
        cache_path = cache_dir / f"{fn}.csv" if cache_dir else None
        if cache_path is not None and cache_path.is_file():
            df = pd.read_csv(cache_path)
            frames.append(df)
            n_kept += len(df)
            print(f"[cache] {fn}: kept={len(df):,} -> total={n_kept:,}/{args.target_rows:,}",
                  flush=True)
            continue

        share = nrows / total
        need_raw = int(np.ceil(need_total_raw * share))
        blocks = plan_blocks(nrows, need_raw, args.block, rng)
        mb = len(blocks) * args.block * rowbytes / 1e6
        print(f"[*] {fn}: {len(blocks)} blocks x {args.block} rows (~{mb:.0f} MB)", flush=True)

        file_frames: list[pd.DataFrame] = []
        file_kept = 0
        t0 = time.time()
        for bi, (a, b) in enumerate(blocks):
            start = datloc + a * rowbytes
            stop = datloc + b * rowbytes
            buf = http_range(url, start, stop, retries=args.retries)
            n_block = len(buf) // rowbytes
            built = build_catalog(block_to_df(buf, dt_be, n_block))
            if len(built):
                file_frames.append(built)
                file_kept += len(built)
            if bi % 5 == 0 or bi == len(blocks) - 1:
                el = time.time() - t0
                print(f"    block {bi+1}/{len(blocks)} kept_file={file_kept:,} "
                      f"({(bi+1)/max(1e-9,el):.2f} blk/s, {el:.0f}s)", flush=True)

        if file_frames:
            df = pd.concat(file_frames, ignore_index=True)
            if cache_path is not None:
                df.to_csv(cache_path, index=False)
            frames.append(df)
            n_kept += len(df)
        print(f"[*] {fn}: kept={file_kept:,} in {time.time()-t0:.0f}s "
              f"-> total={n_kept:,}/{args.target_rows:,}", flush=True)

    if not frames:
        raise RuntimeError("No DES morph rows collected.")

    catalog = pd.concat(frames, ignore_index=True)
    catalog = catalog.drop_duplicates(subset=["COADD_OBJECTS_ID"], keep="first")
    if len(catalog) > args.target_rows:
        catalog = catalog.sample(n=args.target_rows, random_state=args.seed).reset_index(drop=True)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    catalog.to_csv(out_path, index=False)
    print(
        f"Wrote {out_path} ({len(catalog):,} rows)  "
        f"RA=[{catalog['RA_ICRS'].min():.2f},{catalog['RA_ICRS'].max():.2f}]  "
        f"Dec=[{catalog['DE_ICRS'].min():.2f},{catalog['DE_ICRS'].max():.2f}]  "
        f"median ba_r={catalog['ba_r'].median():.4f}  "
        f"median mag_r={catalog['mag_r'].median():.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
