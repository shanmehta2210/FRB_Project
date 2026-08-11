# Unterborn A1 math audit

## Formula (code = paper)

\[
\Delta m_r = 1.27\,(\log_{10} q)^2,\qquad m^f = m - \Delta m
\]

| \(q\) | \(\Delta m\) | \(m=21\to m^f\) |
|------:|-------------:|----------------:|
| 1.0 | 0.000 | 21.000 |
| 0.5 | 0.115 | 20.885 |
| 0.2 | 0.620 | 20.380 |
| 0.1 | 1.270 | 19.730 |

Identities checked: face-on \(m^f=m\); edge-on \(m^f<m\); every raw member (\(m\le\lim\)) satisfies \(m^f\le\lim\) (A1 is always a **superset**). Wrong sign \(m+\Delta m\) would **shrink** N — we observe N increase, so sign is correct.

## Why DES/HSC medians barely move (not a coding bug)

A1 does **not** reweight galaxies already inside the cut. It only **adds** objects with \(m>\lim\) but \(m^f\le\lim\). Shift size ≈ how many enterers × how edge-on they are.

Lim = 22, strict \(b/a>0.2\):

| | med \(b/a\) pool | med \(\Delta m\) | enter N / raw N | enter med \(b/a\) | enter med cos(i) | raw→A1 med cos(i) |
|--|----------------:|-----------------:|----------------:|------------------:|-----------------:|------------------:|
| LS | 0.447 | 0.156 | 116429 / 522951 = **22%** | 0.363 | 0.391 | 0.537→0.506 |
| HSC | 0.581 | 0.071 | 4218 / 36789 = **11%** | 0.385 | 0.336 | 0.576→0.545 |
| DES | 0.707 | 0.029 | 5092 / 60456 = **8%** | 0.377 | 0.327 | 0.628→0.605 |

Enterers are real edge-ons (med \(b/a\sim0.37\), med cos(i)\(\sim0.33\)) with med \(\Delta m\sim0.22\) — math is doing the right thing for those objects. DES/HSC simply have few of them near the limit (face-on-heavy pools → tiny \(\Delta m\) for most galaxies → thin entry strip).

Note: 0.537→0.506 is **LS** at mag 22, not HSC. HSC@22 is 0.576→0.545; DES@22 is 0.628→0.605.

## LS scaling vs dust

REX cos(i) scaling does not change \(q\) used in \(\Delta m(q)\). LS remains dust-active because native \(b/a\) is low.
