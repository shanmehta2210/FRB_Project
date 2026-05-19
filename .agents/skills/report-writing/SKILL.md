---
name: report-writing
description: Ensures reports are concise, correctly formatted, and contain scientifically precise values. Use this skill whenever generating, formatting, or updating reports, PDFs, or analytical CSVs.
---

# Report Writing and Formatting

Follow these rules strictly whenever compiling data, writing reports, or generating summary documents for the user.

## Core Rules

1. **Be Concise**: Get straight to the point. Eliminate filler words and redundant summaries. Focus entirely on the data and the core message.
2. **Precision and Least Count**: Ensure that all numerical values being written are theoretically correct with respect to least count errors. Values should be calculated and rounded strictly up to those precision limits (e.g., significant figures relevant to the instrument or calculation).
3. **No Cutoffs**: Double-check the formatting layout. Ensure that no text, equations, or tables are cutting off the page boundaries or overflowing standard display constraints.
4. **Clean Layout (No Excess Whitespace)**: Strip out any excess whitespace, redundant empty lines, or uncontrolled padding. The report should look tightly formatted and professional.

## How to use it

- When writing table rows or summary data, apply the appropriate `round()` or mathematical limits to your variables before writing them to the file.
- Inspect your markdown or text outputs to ensure block sizes are appropriate and won't overflow when rendered.
- Avoid restating the prompt. Just deliver the requested analysis or table directly.
