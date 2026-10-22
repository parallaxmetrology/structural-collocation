# How a Synthetic Control Manufactured Shared Error

Short technical note on the T1 non-mass-conserving blur (widen cells, keep peaks) and its
correction in T1-bis (mass-conserving Gaussian convolution). Local only; not published.

| file | what |
|---|---|
| `NOTE.html` | the note, self-contained (figures inlined as base64) |
| `minimal.py` | minimal reproduction: same random cells, both blurs; prints analytic checks, S1–S2 error correlation (full canvas and 15 × 15 window), triple-collocation noise fractions and rankings |
| `figures.py` | figure code, called by `minimal.py --figures` |
| `fig1_mechanism.png`, `fig2_shared_error.png`, `fig3_collocation.png` | figures used in the note |

## Run

From this folder, with the project's virtual environment (numpy, scipy; matplotlib only for figures):

```bash
../../.venv/bin/python minimal.py --figures
```

About 30 s; fixed seeds (20261015 fields, 7 bootstrap), 10,000 fields; requires numpy and scipy (matplotlib for figures). The source of `minimal.py` is also embedded in `NOTE.html`. Without `--figures` it prints numbers only.
The figures in `NOTE.html` are inlined copies of the PNGs; if they are regenerated, re-inline them.
