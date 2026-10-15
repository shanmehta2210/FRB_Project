# Re-fits

Staging area for hosts that need constrained re-fits or manual sandboxes.
Do **not** overwrite production `Output/<FRB>_all/` from here.

Full docs: [`../REFIT_AND_REJECT_GRID.md`](../REFIT_AND_REJECT_GRID.md),
[`../SKY_PROTOCOL.md`](../SKY_PROTOCOL.md), [`../SANDBOX.md`](../SANDBOX.md).

## Reject-grid panels (per in-cut reject)

Built by `../run_reject_grid.py`:

| file | what |
|---|---|
| `panel_production.png` | **Byte copy** of `outputs/panels/<FRB>.png`. Never regenerated here. |
| `panel_n1.png` | **n = 1** fixed, sky free. Full verification. |
| `panel_sky.png` | Protocol sky fixed, n free. Full verification. |
| `panel_n1_sky.png` | **n = 1** + protocol sky. Full verification. |

```
Re-fits/<FRB>/
  sky_protocol.json          # independent sky (B/E/F)
  n1/   sky/   n1_sky/       # full GALFIT + verification suite
  psf/                       # optional: host Sérsic + PSF at unmasked star
  n1_psf/                    # same PSF companion with host n=1 fixed
  sersic/                    # optional: host + neighbor Sérsic at the same island
  sandbox/                   # optional hand-edited feedme (run_sandbox.py)
  reject_grid_summary.json
```

`q` is never held fixed in these legs — only `n` and/or sky.
