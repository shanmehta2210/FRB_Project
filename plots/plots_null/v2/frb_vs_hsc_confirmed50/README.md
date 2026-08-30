# Confirmed-50 vs HSC MC \(\cos(i)\) CDF

Paper sample: `in_53` \(\cap\) `confirmed` = **50**. Winning-leg \(q\) (13 re-fit, 37 production). Hubble \(q_0=0.2\).

HSC pool: Kawinwanichakij EXP analogue (`goodfits=1`, \(0.4<n<1.5\)), photometric \(r\le 22\), \(b/a>0.2\). **N = 24,450**. Each of 10,000 draws samples 50 galaxies without replacement.

| | |
|---|---|
| median \(q\) (FRB) | 0.630 |
| median \(\cos(i)\) (point) | 0.610 |
| median \(\sigma_{q,A}\) (\(10\times q_{\rm err}\)) | 0.100 |
| median \(\sigma_{q,B}\) (sky \(\oplus\) \(q_{\rm err}\) \(\oplus\) \(5^\circ\)) | 0.067 |
| hosts missing sky \(\pm\) | 0 |

## Protocol A — `frb_vs_hsc_inflate10.png`

\[\sigma_q = 10\,q_{\rm err},\qquad i = i_{\rm Hubble}(q)+\mathcal{N}(0,5^\circ)\]

## Protocol B — `frb_vs_hsc_sky_quad.png`

\[\sigma_{q,{\rm sky}}=\frac{|q_{+}-q_{-}|}{2},\qquad
\sigma_{q,5^\circ}=\tfrac12\bigl|q(i+5^\circ)-q(i-5^\circ)\bigr|\]

\[\sigma_q=\sqrt{q_{\rm err}^2+\sigma_{q,{\rm sky}}^2+\sigma_{q,5^\circ}^2}\]

No extra \(5^\circ\) after Hubble. \(q_{\rm err}\) is not \(\times 10\).

Table: `confirmed50_q.csv`. Tests: [`TEST_RESULTS.md`](TEST_RESULTS.md).

Literature on GALFIT **axis-ratio** errors (not \(n\)/\(R_e\)/mag): [`GALFIT_Q_ERROR_LITERATURE.md`](GALFIT_Q_ERROR_LITERATURE.md).
