"""Render FIT_VERIFICATION_CHECKS.md to a KaTeX HTML page.

Cursor's native Markdown Preview does not evaluate KaTeX. Open the generated
HTML in a browser (or Cursor Simple Browser) to see the equations.
"""
from __future__ import annotations

import re
import webbrowser
from pathlib import Path

import markdown

HERE = Path(__file__).resolve().parent
SRC = HERE / "FIT_VERIFICATION_CHECKS.md"
DST = HERE / "FIT_VERIFICATION_CHECKS.html"

# Protect math before markdown so underscores in \rm / \mathrm are not eaten
# as emphasis. Prefer $$ / $ in the source; \[ \] / \( \) are accepted too.
DISPLAY_DD = re.compile(r"\$\$.+?\$\$", re.DOTALL)
DISPLAY_BR = re.compile(r"\\\[.+?\\\]", re.DOTALL)
INLINE_DOLLAR = re.compile(r"(?<!\$)\$(?!\$).+?(?<!\$)\$(?!\$)", re.DOTALL)
INLINE_PAREN = re.compile(r"\\\(.+?\\\)", re.DOTALL)


def protect(text: str) -> tuple[str, list[str]]:
    slots: list[str] = []

    def stash(match: re.Match[str]) -> str:
        slots.append(match.group(0))
        return f"§§MATH{len(slots) - 1}§§"

    out = DISPLAY_DD.sub(stash, text)
    out = DISPLAY_BR.sub(stash, out)
    out = INLINE_DOLLAR.sub(stash, out)
    out = INLINE_PAREN.sub(stash, out)
    return out, slots


def restore(html: str, slots: list[str]) -> str:
    for i, tex in enumerate(slots):
        html = html.replace(f"§§MATH{i}§§", tex)
    return html


TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>GALFIT fit verification — method, implementation and results</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.22/dist/katex.min.css"/>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.22/dist/katex.min.js"></script>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.22/dist/contrib/auto-render.min.js"></script>
<style>
  :root {{
    color-scheme: light dark;
    --fg: #1a1a1a;
    --bg: #fafafa;
    --muted: #555;
    --border: #ddd;
    --code-bg: #f0f0f0;
  }}
  @media (prefers-color-scheme: dark) {{
    :root {{
      --fg: #e8e8e8;
      --bg: #1e1e1e;
      --muted: #aaa;
      --border: #333;
      --code-bg: #2a2a2a;
    }}
  }}
  body {{
    font-family: "Segoe UI", system-ui, sans-serif;
    line-height: 1.55;
    max-width: 52rem;
    margin: 2rem auto;
    padding: 0 1.25rem 4rem;
    color: var(--fg);
    background: var(--bg);
  }}
  h1, h2, h3, h4 {{ line-height: 1.25; }}
  h1 {{ font-size: 1.75rem; }}
  h2 {{ margin-top: 2.2rem; border-bottom: 1px solid var(--border); padding-bottom: 0.25rem; }}
  h3 {{ margin-top: 1.6rem; }}
  code, pre {{ font-family: Consolas, "Cascadia Mono", monospace; font-size: 0.92em; }}
  code {{ background: var(--code-bg); padding: 0.1em 0.35em; border-radius: 3px; }}
  pre {{ background: var(--code-bg); padding: 0.9rem 1rem; overflow-x: auto; border-radius: 6px; }}
  pre code {{ background: none; padding: 0; }}
  table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; font-size: 0.95em; }}
  th, td {{ border: 1px solid var(--border); padding: 0.4rem 0.55rem; text-align: left; vertical-align: top; }}
  th {{ background: var(--code-bg); }}
  blockquote {{
    margin: 1rem 0; padding: 0.5rem 1rem; border-left: 4px solid #5b8def;
    background: var(--code-bg); color: var(--muted);
  }}
  .banner {{
    background: #fff3cd; color: #664d03; border: 1px solid #ffecb5;
    padding: 0.75rem 1rem; border-radius: 6px; margin-bottom: 1.5rem;
  }}
  @media (prefers-color-scheme: dark) {{
    .banner {{ background: #3d3415; color: #ffe08a; border-color: #6b5a1e; }}
  }}
  .katex-display {{ margin: 1rem 0; overflow-x: auto; overflow-y: hidden; }}
</style>
</head>
<body>
<div class="banner">
  Rendered HTML with KaTeX. Cursor's native Markdown Preview does not evaluate
  math; use this page (or classic <code>Ctrl+Shift+V</code> preview) instead.
  Source: <code>FIT_VERIFICATION_CHECKS.md</code>.
</div>
{body}
<script>
document.addEventListener("DOMContentLoaded", function () {{
  renderMathInElement(document.body, {{
    delimiters: [
      {{left: "$$", right: "$$", display: true}},
      {{left: "$", right: "$", display: false}},
      {{left: "\\\\(", right: "\\\\)", display: false}},
      {{left: "\\\\[", right: "\\\\]", display: true}}
    ],
    throwOnError: false,
    ignoredTags: ["script", "noscript", "style", "textarea", "pre", "code"]
  }});
}});
</script>
</body>
</html>
"""


def main(open_browser: bool = True) -> Path:
    raw = SRC.read_text(encoding="utf-8")
    # Drop the Cursor-specific rendering callout; the HTML banner replaces it.
    raw = re.sub(
        r"> \*\*How to render the equations in Cursor\.\*\*.*?(?=\n\*\*Contents\*\*)",
        "",
        raw,
        count=1,
        flags=re.DOTALL,
    )
    protected, slots = protect(raw)
    body = markdown.markdown(
        protected,
        extensions=["tables", "fenced_code", "sane_lists", "toc"],
        extension_configs={"toc": {"permalink": False}},
    )
    body = restore(body, slots)
    DST.write_text(TEMPLATE.format(body=body), encoding="utf-8")
    print(f"Wrote {DST} ({DST.stat().st_size // 1024} KB, {len(slots)} math spans)")
    if open_browser:
        webbrowser.open(DST.as_uri())
    return DST


if __name__ == "__main__":
    import sys

    main(open_browser="--no-open" not in sys.argv)
