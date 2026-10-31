"""
Unpack the tracked HSC control catalogues into catalog/ (git-ignored) and verify them.

The analysis reads catalog/HSC_kawinwanichakij_exp_analogue.csv (HSC disc pool) and
catalog/HSC_kawinwanichakij_sample_500k.csv (parent sample for the q0 / Sersic-window
appendix controls). Compressed copies live in data/catalog_snapshots/.

Run from repo root::

    python scripts/restore_catalogs.py
"""
from __future__ import annotations

import gzip
import hashlib
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SNAP = ROOT / "data" / "catalog_snapshots"
DEST = ROOT / "catalog"
SHA256 = {
    "HSC_kawinwanichakij_exp_analogue.csv": "d0838cc0d34e5e9c67e52bc5baa29abc7ef54af9b465e16264f3276f1d470786",
    "HSC_kawinwanichakij_sample_500k.csv": "9239cc03e6a9a576108638cb67569891f9745e17fcbd85abeaafddc5571b3559",
}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    DEST.mkdir(exist_ok=True)
    for name, digest in SHA256.items():
        out = DEST / name
        if not (out.exists() and sha256(out) == digest):
            with gzip.open(SNAP / (name + ".gz"), "rb") as a, open(out, "wb") as b:
                shutil.copyfileobj(a, b)
        ok = sha256(out) == digest
        print(f"{name}: {'OK' if ok else 'CHECKSUM MISMATCH'}")
        if not ok:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
