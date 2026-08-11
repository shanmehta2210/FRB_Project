// Compile every math span in the spec doc with the real KaTeX, so a macro the
// preview cannot render is caught here rather than in the preview pane.
const fs = require('fs');
const katex = require(process.env.KATEX_PATH || 'katex');

const doc = fs.readFileSync(process.argv[2], 'utf8');
const lines = doc.split(/\r?\n/);

// Strip fenced code blocks: they are not math.
let inFence = false;
const prose = lines.map((l) => {
  if (l.trimStart().startsWith('```')) { inFence = !inFence; return ''; }
  return inFence ? '' : l;
}).join('\n');

const spans = [];
const display = prose.matchAll(/\$\$([\s\S]+?)\$\$/g);
for (const m of display) spans.push({ tex: m[1], display: true });
const noDisplay = prose.replace(/\$\$[\s\S]+?\$\$/g, '');
// Inline code can quote dollar signs that are not math.
const inline = noDisplay.replace(/`[^`]*`/g, '').matchAll(/\$([^$\n]+)\$/g);
for (const m of inline) spans.push({ tex: m[1], display: false });

let failed = 0;
for (const s of spans) {
  try {
    katex.renderToString(s.tex, { displayMode: s.display, throwOnError: true });
  } catch (e) {
    failed += 1;
    console.log(`FAIL ${s.display ? 'display' : 'inline'}: ${s.tex.slice(0, 70)}`);
    console.log(`     ${e.message.split('\n')[0]}`);
  }
}
console.log(`${spans.length} math spans compiled, ${failed} failed`);
process.exit(failed ? 1 : 0);
