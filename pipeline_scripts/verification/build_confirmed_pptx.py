"""Build one-slide-per-panel PPTX for confirmed science-cut hosts.

Panel PNG files on disk are **never modified**. For slide embedding only,
RGBA is flattened to RGB in memory (PowerPoint renders RGBA poorly) and each
picture is placed at a fixed high DPI so the full pixel grid is stored and
drawn at the correct physical size (no shape rescale after insert).

See ``VISUAL_PANELS.md``.
"""

from __future__ import annotations

import io
import re
import sys
from pathlib import Path

import pandas as pd
from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE_TYPE
from pptx.enum.text import PP_ALIGN
from pptx.util import Emu, Pt

VER = Path(__file__).resolve().parent
if str(VER) not in sys.path:
    sys.path.insert(0, str(VER))

import vercommon as vc  # noqa: E402

PANELS = Path(vc.OUT_ROOT) / "panels"
HC = VER / "host_confirmation.csv"
OUT = Path(vc.OUT_ROOT) / "confirmed_fit_panels.pptx"

# 1749 px @ 200 DPI ≈ 8.745" — sharp on retina / export; full PNG in package.
EMBED_DPI = 200
TITLE_H_IN = 0.40
MARGIN_IN = 0.12
EMU_PER_INCH = 914400

_ALT = re.compile(
    r"outputs/panels/([A-Za-z0-9]+_(?:n1_sky|n1|sky|psf))\.png"
)


def panel_for(frb: str, notes: str) -> Path:
    m = _ALT.search(notes or "")
    if m:
        p = PANELS / f"{m.group(1)}.png"
        if p.is_file():
            return p
    return PANELS / f"{frb}.png"


def _embed_stream(path: Path) -> tuple[io.BytesIO, int, int]:
    """Load panel for embedding; flatten alpha → white RGB without touching disk."""
    with Image.open(path) as im:
        im.load()
        w, h = im.size
        if im.mode == "RGBA":
            bg = Image.new("RGB", im.size, (255, 255, 255))
            bg.paste(im, mask=im.split()[3])
            rgb = bg
        elif im.mode != "RGB":
            rgb = im.convert("RGB")
        else:
            rgb = im.copy()
    buf = io.BytesIO()
    rgb.save(buf, format="PNG", compress_level=3, optimize=False)
    buf.seek(0)
    return buf, w, h


def _px_to_emu(px: int, dpi: int = EMBED_DPI) -> int:
    return int(round(px / dpi * EMU_PER_INCH))


def _configure_slide_size(prs: Presentation, px_w: int, px_h: int) -> None:
    prs.slide_width = Emu(_px_to_emu(px_w) + int(2 * MARGIN_IN * EMU_PER_INCH))
    prs.slide_height = Emu(
        _px_to_emu(px_h) + int((TITLE_H_IN + 2 * MARGIN_IN) * EMU_PER_INCH)
    )


def add_image_slide(
    prs: Presentation,
    title: str,
    img: Path,
    *,
    title_color: RGBColor | None = None,
    subtitle: str = "",
) -> None:
    stream, px_w, px_h = _embed_stream(img)
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    sw = int(prs.slide_width)

    box = slide.shapes.add_textbox(
        Emu(int(MARGIN_IN * EMU_PER_INCH)),
        Emu(int(0.06 * EMU_PER_INCH)),
        Emu(sw - int(2 * MARGIN_IN * EMU_PER_INCH)),
        Emu(int(TITLE_H_IN * EMU_PER_INCH)),
    )
    tf = box.text_frame
    tf.clear()
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = title
    p.font.size = Pt(18)
    p.font.bold = True
    p.font.name = "Calibri"
    p.font.color.rgb = title_color or RGBColor(0x1A, 0x1A, 0x1A)
    if subtitle:
        p2 = tf.add_paragraph()
        p2.text = subtitle
        p2.font.size = Pt(11)
        p2.font.name = "Calibri"
        p2.font.color.rgb = RGBColor(0x55, 0x55, 0x55)

    pic_w = _px_to_emu(px_w)
    pic_h = _px_to_emu(px_h)
    left = (sw - pic_w) // 2
    top = int((MARGIN_IN + TITLE_H_IN) * EMU_PER_INCH)
    slide.shapes.add_picture(
        stream, Emu(left), Emu(top), width=Emu(pic_w), height=Emu(pic_h)
    )


def add_title_slide(prs: Presentation, n_confirmed: int) -> None:
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    sw, sh = int(prs.slide_width), int(prs.slide_height)
    box = slide.shapes.add_textbox(
        Emu(int(0.5 * EMU_PER_INCH)),
        Emu(sh // 2 - int(1.2 * EMU_PER_INCH)),
        Emu(sw - int(1.0 * EMU_PER_INCH)),
        Emu(int(2.4 * EMU_PER_INCH)),
    )
    tf = box.text_frame
    tf.clear()
    p = tf.paragraphs[0]
    p.text = "FRB host GALFIT verification panels"
    p.font.size = Pt(32)
    p.font.bold = True
    p.font.name = "Calibri"
    p.alignment = PP_ALIGN.CENTER
    p2 = tf.add_paragraph()
    p2.text = (
        f"{n_confirmed} confirmed hosts (science cut: mag ≤ 22, b/a > 0.2)\n"
        f"Full-resolution panels embedded at {EMBED_DPI} DPI"
    )
    p2.font.size = Pt(16)
    p2.font.name = "Calibri"
    p2.alignment = PP_ALIGN.CENTER
    p2.font.color.rgb = RGBColor(0x44, 0x44, 0x44)


def main() -> int:
    hc = pd.read_csv(HC)
    hc["confirmed"] = hc["confirmed"].astype(str).str.lower().eq("true")
    cohort = vc.cohort("all64")[["frb", "mag", "b_a", "in_53"]]
    m = hc.merge(cohort, on="frb", how="left")
    m["in_53"] = m["in_53"].astype(bool)

    confirmed = m[m["confirmed"] & m["in_53"]].sort_values("frb")
    rejected = m[(~m["confirmed"]) & m["in_53"]].sort_values("frb")

    # Resolve panel paths first; size the deck from a real panel.
    jobs: list[tuple[str, str, Path, str]] = []
    missing: list[str] = []
    for i, (_, row) in enumerate(confirmed.iterrows(), 1):
        frb = str(row["frb"])
        notes = "" if pd.isna(row.get("notes")) else str(row["notes"])
        img = panel_for(frb, notes)
        if not img.is_file():
            missing.append(frb)
            continue
        alt = _ALT.search(notes)
        leg = alt.group(1).split("_", 1)[1] if alt else "production"
        jobs.append(
            (f"{i}/{len(confirmed)}   {frb}   [{leg}]", "", img, "ok")
        )

    reject_jobs: list[tuple[str, str, Path]] = []
    for _, row in rejected.iterrows():
        frb = str(row["frb"])
        notes = "" if pd.isna(row.get("notes")) else str(row["notes"])
        img = PANELS / f"{frb}.png"
        if img.is_file():
            reject_jobs.append(
                (
                    f"{frb}   REJECTED — not in confirmed sample",
                    (notes[:160] + "…") if len(notes) > 160 else notes,
                    img,
                )
            )

    probe = jobs[0][2] if jobs else (reject_jobs[0][2] if reject_jobs else None)
    if probe is None:
        print("no panels to embed", file=sys.stderr)
        return 1
    with Image.open(probe) as im:
        px_w, px_h = im.size

    prs = Presentation()
    _configure_slide_size(prs, px_w, px_h)
    add_title_slide(prs, len(confirmed))

    for title, _sub, img, _ in jobs:
        add_image_slide(prs, title, img)

    for title, sub, img in reject_jobs:
        add_image_slide(
            prs,
            title,
            img,
            title_color=RGBColor(0xA0, 0x20, 0x20),
            subtitle=sub,
        )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    prs.save(str(OUT))

    # Re-open to verify blob sizes (full PNG in package).
    verify = Presentation(str(OUT))
    sizes = [
        len(shape.image.blob)
        for slide in verify.slides
        for shape in slide.shapes
        if shape.shape_type == MSO_SHAPE_TYPE.PICTURE
    ]
    src = probe.stat().st_size
    print(f"wrote {OUT}  ({OUT.stat().st_size // 1024} KB)")
    print(f"confirmed slides: {len(jobs)}  missing: {missing}")
    print(f"rejected slides: {len(reject_jobs)}")
    print(
        f"panel pixels {px_w}x{px_h} @ {EMBED_DPI} DPI -> "
        f'{px_w / EMBED_DPI:.2f}" x {px_h / EMBED_DPI:.2f}"'
    )
    if sizes:
        print(
            f"embedded blobs: min={min(sizes) // 1024}KB  "
            f"med={sorted(sizes)[len(sizes) // 2] // 1024}KB  "
            f"max={max(sizes) // 1024}KB  (source file {src // 1024}KB, untouched)"
        )
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
