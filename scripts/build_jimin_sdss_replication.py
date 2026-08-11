"""
Reproduce the advisor/Jimin SDSS DR16 SkyServer selection + Hubble cos(i) CDFs.

Two morphology variants (shared base WHERE):
  V1 (historical): fracDeV_r = 0 AND lnLDeV_r < lnLExp_r
  V2 (weaker):     lnLDeV_r < lnLExp_r only

Mag cut (``--mag-cut``):
  model (default): p.r BETWEEN 12 AND 21  [= modelMag_r, historical SQL]
  petro:           p.petroMag_r BETWEEN 12 AND 21  -> outputs under Jimin/petroMag/

Then: expAB_r > q0 (=0.2), Hubble cos(i). Outputs under
plots/plots_null/v2/sdss_audit/Jimin/.

Run from repo root::

    python scripts/build_jimin_sdss_replication.py
    python scripts/build_jimin_sdss_replication.py --mag-cut petro
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astroquery.sdss import SDSS

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import Q0, cosi_array_from_df  # noqa: E402
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

OUT_ROOT = REPO_ROOT / "plots" / "plots_null" / "v2" / "sdss_audit" / "Jimin"
DATA_RELEASE = 16
# Historical SkyServer SQL used TOP 10000; that truncates V2 (true N~17.7k).
# Default: no TOP so CasJobs-style full counts. Use --top 10000 to match SkyServer.
DEFAULT_TOP_N: int | None = None

SELECT_COLS = """
p.objid, p.ra, p.dec,
p.r AS modelMag_r,
p.cModelMag_u, p.cModelMag_g, p.cModelMag_r,
p.petroR90_r, p.petroR50_r,
p.petroMag_r, p.petroRad_r,
p.deVRad_r, p.deVAB_r, p.lnLDeV_r,
p.expRad_r, p.expAB_r, p.lnLExp_r,
p.lnLStar_r,
p.fracDeV_r, pz.nnAvgZ AS photz
"""

# Mag WHERE options. In SDSS PhotoObj, shorthand p.r ≡ modelMag_r (not petroMag_r).
MAG_CUTS = {
    "model": {
        "sql": "p.r BETWEEN 12 AND 21",
        "label": "modelMag_r (p.r)",
        "subdir": None,  # historical default -> Jimin/
    },
    "petro": {
        "sql": "p.petroMag_r BETWEEN 12 AND 21",
        "label": "petroMag_r",
        "subdir": "petroMag",
    },
}


def base_where(mag_sql: str) -> str:
    return f"""
p.ra BETWEEN 148.0 AND 152.0
AND p.dec BETWEEN 0.0 AND 4.0
AND {mag_sql}
AND p.mode = 1
AND p.clean = 1
AND p.type_r = 3
AND p.lnLStar_r < -10
AND pz.nnAvgZ > 0
AND p.score > 0.8
"""

VARIANTS = {
    "v1_fracDev0_and_lnL": {
        "label": "V1: fracDeV=0 AND lnLExp wins (historical)",
        "morph_sql": "AND p.fracDeV_r = 0 AND p.lnLDeV_r < p.lnLExp_r",
        "short": "fracDev0+lnL",
        "color": "#e41a1c",
    },
    "v2_lnL_exp": {
        "label": "V2: lnLExp wins only",
        "morph_sql": "AND p.lnLDeV_r < p.lnLExp_r",
        "short": "lnL only",
        "color": "#377eb8",
    },
}


def build_sql(morph_sql: str, *, top_n: int | None, mag_sql: str) -> str:
    top_clause = f"TOP {int(top_n)} " if top_n and top_n > 0 else ""
    return f"""
SELECT {top_clause}{SELECT_COLS}
FROM PhotoObj AS p
JOIN Photoz AS pz ON pz.objid = p.objid
WHERE
{base_where(mag_sql).strip()}
{morph_sql}
""".strip()


def empirical_cdf(vals: np.ndarray, x: np.ndarray) -> np.ndarray:
    s = np.sort(vals[np.isfinite(vals)])
    return np.searchsorted(s, x, side="right") / max(1, len(s))


def query_variant(
    key: str,
    morph_sql: str,
    *,
    top_n: int | None,
    mag_sql: str,
    timeout: int,
    retries: int,
) -> pd.DataFrame:
    sql = build_sql(morph_sql, top_n=top_n, mag_sql=mag_sql)
    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            print(f"[*] Querying DR{DATA_RELEASE} {key} (attempt {attempt}) ...", flush=True)
            tbl = SDSS.query_sql(sql, data_release=DATA_RELEASE, timeout=timeout)
            if tbl is None:
                return pd.DataFrame()
            return tbl.to_pandas()
        except Exception as exc:
            last_err = exc
            print(f"    failed: {exc}", flush=True)
            if attempt < retries:
                time.sleep(2.0 * attempt)
    raise RuntimeError(f"Query {key} failed after {retries} attempts: {last_err}")


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [c.lower() for c in out.columns]
    # unify names used later
    rename = {}
    if "expab_r" in out.columns:
        rename["expab_r"] = "expAB_r"
    if "fracdev_r" in out.columns:
        rename["fracdev_r"] = "fracDeV_r"
    if "lnldev_r" in out.columns:
        rename["lnldev_r"] = "lnLDeV_r"
    if "lnlexp_r" in out.columns:
        rename["lnlexp_r"] = "lnLExp_r"
    if "objid" in out.columns:
        rename["objid"] = "objID"
    return out.rename(columns=rename)


def strict_cosi(df: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame]:
    ba = pd.to_numeric(df["expAB_r"], errors="coerce")
    ok = ba.notna() & (ba > Q0) & (ba <= 1.0)
    sub = df.loc[ok].copy()
    cosi = cosi_array_from_df(sub, q_col="expAB_r", q0=Q0)
    return cosi, sub


def plot_cdf(cosi: np.ndarray, *, title: str, color: str, out: Path) -> float:
    med = float(np.median(cosi)) if len(cosi) else float("nan")
    x = np.linspace(0, 1, 401)
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot((0, 1), (0, 1), "k--", lw=1.2, label="Uniform")
    if len(cosi):
        ax.plot(x, empirical_cdf(cosi, x), color=color, lw=2.0, label="SDSS Jimin null")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(rf"$\cos(i)$ (Hubble, $q_0={Q0:g}$)")
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(f"{title}\nN = {len(cosi):,}  |  median cos(i) = {med:.3f}")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return med


def write_readme(
    out_dir: Path,
    *,
    rows: list[dict],
    truncated: dict[str, bool],
    top_n: int | None,
    mag_key: str,
    mag_label: str,
    mag_sql: str,
) -> None:
    r1 = next(r for r in rows if r["variant"].startswith("v1"))
    r2 = next(r for r in rows if r["variant"].startswith("v2"))
    top_desc = f"`TOP {top_n}`" if top_n and top_n > 0 else "no TOP (full COUNT)"
    trunc_note = ""
    if any(truncated.values()):
        trunc_note = (
            f"\n**Warning:** at least one query returned exactly {top_n} rows — "
            "truncation may have dropped objects.\n"
        )

    mag_note = ""
    if mag_key == "model":
        mag_note = (
            "\n**Mag cut confirmation:** SQL uses `p.r BETWEEN 12 AND 21`. In SDSS "
            "`PhotoObj`, shorthand `r` ≡ `modelMag_r` (exact match in DR16 probe: "
            "`max|r−modelMag_r|=0`). It is **not** `petroMag_r`. See sibling "
            "`petroMag/` for the Petrosian variant.\n"
        )
    else:
        mag_note = (
            "\n**Mag cut:** `p.petroMag_r BETWEEN 12 AND 21` (Petrosian). "
            "Historical Jimin SQL used `p.r` (= modelMag). Compare to parent "
            "`Jimin/` (model).\n"
        )

    text = f"""# Jimin / advisor SDSS DR16 query replication

Reproduction of the historical SkyServer SQL selection used by the advisor + Jimin,
plus Hubble cos(i) CDFs. Release: **DR16**. Mag cut: **{mag_label}** (`{mag_sql}`).

{trunc_note}{mag_note}
## Morphology cuts (two versions)

| Version | Morph SQL | Meaning |
|---------|-----------|---------|
| **V1** (historical) | `fracDeV_r = 0` **AND** `lnLDeV_r < lnLExp_r` | Pure cModel exponential **and** pure-exp likelihood wins |
| **V2** (weaker) | `lnLDeV_r < lnLExp_r` only | Pure-exp likelihood wins (SDSS `modelMag` rule) |

Shared base: RA 148–152, Dec 0–4, `{mag_sql}`, `mode=1`, `clean=1`, `type_r=3`,
`lnLStar_r < -10`, Photoz `nnAvgZ > 0`, `score > 0.8`, {top_desc}.

After fetch: keep **`expAB_r` > {Q0:g}**, Hubble cos(i) with `q0={Q0:g}`
(`cosi_array_from_df(..., q_col='expAB_r')`). Not `deVAB_r`.

## Why V1 is ~9k (not a fetch bug)

`COUNT(*)` **without TOP** on DR16 for the exact V1 WHERE with **model** `p.r` returns
**9,120** — identical to the fetched catalog. Historical `TOP 10000` never limited V1.

Ablated funnel (`scripts/_jimin_count_audit.py`, DR16, no TOP, **model** mag):

| Stage | N |
|-------|--:|
| Box + mode/clean + `type_r=3` | 144,070 |
| + `r` in [12, 21] | 47,351 |
| + `lnLStar_r < -10` | 33,045 |
| + Photoz `nnAvgZ > 0` | 32,978 |
| + `score > 0.8` (base, no morph) | 28,043 |
| + lnLExp wins (**V2 full**) | **17,657** |
| + `fracDeV_r = 0` (**V1 full**) | **9,120** |

Petrosian swap (`petroMag_r` in [12, 21], same other cuts): V2 **17,138**, V1 **8,668**
— slightly *fewer*, not more. Overlap (V2): both mags in range 17,096; model-only 561;
petro-only 42.

## What `fracDeV` and `lnL*` mean (SDSS photometry)

Per [SDSS DR16 magnitudes](https://www.sdss4.org/dr16/algorithms/magnitudes/):

1. Fit pure exponential → `lnLExp_r`
2. Fit pure de Vaucouleurs → `lnLDeV_r`
3. Linear mix of those two (ellipses fixed) → `fracDeV_r` = deV weight in cModel

So `lnLDeV < lnLExp` only asks which *single* profile wins. `fracDeV = 0` requires the
**composite** to put zero bulge light — the stronger pure-disk cut.

## Results

| Version | N (SQL) | N (ba > {Q0:g}) | median cos(i) | frac of V2 with fracDeV=0 |
|---------|--------:|----------------:|--------------:|--------------------------:|
| V1 fracDev0+lnL | {r1['n_sql']:,} | {r1['n_strict']:,} | {r1['median_cosi']} | — |
| V2 lnL only | {r2['n_sql']:,} | {r2['n_strict']:,} | {r2['median_cosi']} | {r2.get('frac_v2_with_fracDev0', float('nan'))} |

TOP truncated: V1={truncated.get('v1_fracDev0_and_lnL')}, V2={truncated.get('v2_lnL_exp')} (top_n={top_n}).

## Files

- `sql/v1_fracDev0_and_lnL.sql`, `sql/v2_lnL_exp.sql`
- `catalog/v1_*.csv`, `catalog/v2_*.csv`
- `plots/cdf_*.png`, `summary.csv`

```bash
python scripts/build_jimin_sdss_replication.py
python scripts/build_jimin_sdss_replication.py --mag-cut petro
```
"""
    (out_dir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument("--timeout", type=int, default=900)
    p.add_argument("--retries", type=int, default=3)
    p.add_argument(
        "--mag-cut",
        choices=sorted(MAG_CUTS.keys()),
        default="model",
        help="Magnitude column for WHERE [12,21]: model (=p.r) or petro (petroMag_r).",
    )
    p.add_argument(
        "--top",
        type=int,
        default=-1,
        help="SQL TOP N (historical SkyServer=10000). Omit/negative = no TOP (default).",
    )
    p.add_argument(
        "--skip-fetch",
        action="store_true",
        help="Reuse existing catalog CSVs (no SDSS query).",
    )
    args = p.parse_args()
    top_n: int | None
    if args.top is None or args.top < 0:
        top_n = DEFAULT_TOP_N
    else:
        top_n = int(args.top)

    mag_meta = MAG_CUTS[args.mag_cut]
    mag_sql = mag_meta["sql"]
    mag_label = mag_meta["label"]
    if args.out_dir is not None:
        out = args.out_dir
    elif mag_meta["subdir"]:
        out = OUT_ROOT / mag_meta["subdir"]
    else:
        out = OUT_ROOT

    sql_dir = out / "sql"
    cat_dir = out / "catalog"
    plot_dir = out / "plots"
    for d in (sql_dir, cat_dir, plot_dir):
        d.mkdir(parents=True, exist_ok=True)

    frames: dict[str, pd.DataFrame] = {}
    truncated: dict[str, bool] = {}

    for key, meta in VARIANTS.items():
        sql = build_sql(meta["morph_sql"], top_n=top_n, mag_sql=mag_sql)
        (sql_dir / f"{key}.sql").write_text(sql + "\n", encoding="utf-8")
        csv_path = cat_dir / f"{key}.csv"
        if args.skip_fetch and csv_path.exists():
            print(f"[*] Loading cached {csv_path.name}", flush=True)
            df = pd.read_csv(csv_path)
        else:
            df = query_variant(
                key,
                meta["morph_sql"],
                top_n=top_n,
                mag_sql=mag_sql,
                timeout=args.timeout,
                retries=args.retries,
            )
            df = normalize_columns(df)
            df.to_csv(csv_path, index=False)
        if "expAB_r" not in df.columns and "expab_r" in df.columns:
            df = normalize_columns(df)
        frames[key] = df
        truncated[key] = bool(top_n and top_n > 0 and len(df) == top_n)
        print(
            f"    {key}: N_sql={len(df):,}  mag={args.mag_cut}  truncated={truncated[key]}",
            flush=True,
        )

    # Cos(i) CDFs
    x = np.linspace(0, 1, 401)
    rows: list[dict] = []
    series_for_overlay: list[tuple[str, np.ndarray, str, str]] = []

    plot_names = {
        "v1_fracDev0_and_lnL": "cdf_fracDev0_and_lnL.png",
        "v2_lnL_exp": "cdf_lnL_exp.png",
    }

    for key, meta in VARIANTS.items():
        df = frames[key]
        cosi, strict = strict_cosi(df)
        med = plot_cdf(
            cosi,
            title=f"Jimin DR16 — {meta['label']} | {mag_label}",
            color=meta["color"],
            out=plot_dir / plot_names[key],
        )
        series_for_overlay.append((meta["short"], cosi, meta["color"], meta["label"]))

        frac_dev0 = float("nan")
        if key.startswith("v2") and "fracDeV_r" in df.columns:
            fd = pd.to_numeric(df["fracDeV_r"], errors="coerce")
            frac_dev0 = float((fd == 0).mean()) if len(fd) else float("nan")

        row = {
            "variant": key,
            "label": meta["label"],
            "mag_cut": args.mag_cut,
            "mag_sql": mag_sql,
            "n_sql": len(df),
            "n_strict": len(cosi),
            "median_cosi": round(med, 4) if len(cosi) else np.nan,
            "top_truncated": truncated[key],
            "top_n": top_n if top_n is not None else "",
            "q0": Q0,
        }
        if key.startswith("v2"):
            row["frac_v2_with_fracDev0"] = round(frac_dev0, 4)
        rows.append(row)
        print(
            f"    {key}: N_strict={len(cosi):,}  med={med:.3f}",
            flush=True,
        )

    # Overlay
    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    ax.plot((0, 1), (0, 1), "k--", lw=1.2, label="Uniform")
    for short, cosi, col, lab in series_for_overlay:
        med = float(np.median(cosi)) if len(cosi) else float("nan")
        ax.plot(
            x,
            empirical_cdf(cosi, x),
            color=col,
            lw=2.0,
            label=f"{short}  N={len(cosi):,}, med={med:.3f}",
        )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(rf"$\cos(i)$ (Hubble, $q_0={Q0:g}$)")
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(
        f"Jimin DR16 — V1 vs V2 | mag cut: {mag_label}\n"
        rf"expAB$_r$ > {Q0:g}; RA 148–152, Dec 0–4; {mag_sql}"
    )
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(plot_dir / "cdf_overlay.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(rows).to_csv(out / "summary.csv", index=False)
    write_readme(
        out,
        rows=rows,
        truncated=truncated,
        top_n=top_n,
        mag_key=args.mag_cut,
        mag_label=mag_label,
        mag_sql=mag_sql,
    )
    print(f"[*] Wrote outputs -> {out}", flush=True)


if __name__ == "__main__":
    main()
