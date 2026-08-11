"""Lint the spec doc: KaTeX delimiters, table integrity, unescaped pipes.

Cursor's built-in preview only understands ``$...$`` / ``$$...$$``. It also
splits table rows on every raw ``|``, so a ``|x|`` inside math silently
destroys the table.
"""

import os
import re
import sys

DOC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "FIT_VERIFICATION_CHECKS.md")

text = open(DOC, encoding="utf-8").read()
lines = text.splitlines()
problems: list[str] = []

for i, line in enumerate(lines, 1):
    # Inline code may legitimately quote the delimiters that do not render.
    prose = re.sub(r"`[^`]*`", "", line)
    for bad in (r"\(", r"\)", r"\[", r"\]"):
        if bad in prose and not line.lstrip().startswith("```"):
            problems.append(f"{i}: unsupported delimiter {bad!r}: {line.strip()[:70]}")

# Fenced code blocks are not math and not tables.
in_fence = False
fence: list[bool] = []
for line in lines:
    if line.lstrip().startswith("```"):
        in_fence = not in_fence
        fence.append(True)
    else:
        fence.append(in_fence)

# $ must balance on every non-fenced line (display $$ live on their own lines).
for i, line in enumerate(lines, 1):
    if fence[i - 1] or line.strip() == "$$":
        continue
    stripped = re.sub(r"`[^`]*`", "", line)
    if stripped.count("$") % 2:
        problems.append(f"{i}: odd number of '$': {line.strip()[:70]}")

# A raw | inside math inside a table cell breaks the row.
for i, line in enumerate(lines, 1):
    if fence[i - 1] or not line.strip().startswith("|"):
        continue
    for m in re.finditer(r"\$([^$]*)\$", line):
        if "|" in m.group(1):
            problems.append(f"{i}: raw '|' in math inside a table: {m.group(0)[:50]}")

# Table rows must all have the same cell count within a block.
block: list[tuple[int, int]] = []
for i, line in enumerate(lines + [""], 1):
    if line.strip().startswith("|") and not fence[min(i, len(fence)) - 1]:
        block.append((i, line.strip().strip("|").count("|") + 1))
        continue
    if block:
        widths = {w for _, w in block}
        if len(widths) > 1:
            rows = ", ".join(f"L{ln}={w}" for ln, w in block)
            problems.append(f"table starting L{block[0][0]}: ragged widths ({rows})")
        block = []

katex_risky = re.findall(r"\\arctan_|\\mbox|\\bf\b|\\rm\{", text)
if katex_risky:
    problems.append(f"KaTeX-unsupported macros: {sorted(set(katex_risky))}")

n_inline = len(re.findall(r"(?<!\$)\$(?!\$)[^$\n]+\$(?!\$)", text))
n_display = text.count("$$") // 2
print(f"{len(lines)} lines, {n_inline} inline math spans, {n_display} display blocks")

if problems:
    print(f"\n{len(problems)} problems:")
    for p in problems:
        print("  " + p)
    sys.exit(1)
print("clean")
