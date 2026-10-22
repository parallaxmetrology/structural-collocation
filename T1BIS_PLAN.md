# T1-bis: synthetic assumption audit (frozen 2026-10-15, before any T1-bis value was computed)

Why: an external review, checked against the code, found that T1's simulated blur widened
each Gaussian cell with its peak fixed (rain mass not conserved; small cells grow
proportionally more), added widths linearly rather than in quadrature, and read fields
inside 15 × 15 windows that clip smoother observers more. Each can create a shared
coordinate error by itself. T1-bis rebuilds the control with an analytical reference,
measures the assumption collocation needs (independent coordinate errors) directly on
the extracted coordinates, and re-tests T1's claims.

## Fields

- Canvas 64 × 64 cells (0.1° each). Truth: 1 + Poisson(1.5) Gaussian cells, centres
  uniform within ±5 cells of the canvas centre, widths U(0.8, 2.5), peaks
  lognormal(1.5, 0.8). Truth is never thresholded in T1-bis.
- 20,000 truth windows, seed 20261016.
- **Blur = convolution** of the rendered field with a normalised Gaussian kernel of width
  σ_b (`scipy.ndimage.gaussian_filter`, mode constant, truncate 4). This conserves mass
  and, on the full canvas, leaves the centroid unchanged and adds 2σ_b² to squared
  dispersion.
- Observers re-render the truth cells with per-cell errors, then blur:

| observer | position error sd | blur σ_b | intensity error (lognormal sd) |
|---|---|---|---|
| R (radar-like) | 0.5 | 0 | 0.3 |
| S1 (microwave-like) | 0.8 | 1.0 | 0.4 |
| S2 (infrared-like) | 1.5 | 2.0 | 0.6 |
| Gs (sharp independent) | 1.0 | 0 | 0.5 |
| Gb (blurred independent) | 1.0 | 1.5 | 0.5 |
| H (second sharp independent, benchmark only) | 0.7 | 0 | 0.4 |

## Reads

Centroid (cy, cx), dispersion d, squared dispersion d², each computed (i) on the **full
canvas** and (ii) in the **central 15 × 15 window**, plus (iii) windows shifted by
(+3, +3) and (−3, +4) cells to test boundary placement. A read is valid if the field's
mass in the read region is ≥ 1.

## Stages, predictions and pass criteria

**S-A analytical reference** (truth only, no errors, blur σ_b ∈ {1, 2}, full canvas).
Pass: median |centroid shift| < 0.01 cell and median relative error of
(d²_blurred − d²_truth − 2σ_b²) < 1%. If S-A fails, the rendering is wrong: stop and
report.

**S-B independence on the full canvas** (errors as in the table, no planted shared error).
Measure the true coordinate-error correlations (residual of each observer on truth,
pairwise). Prediction: for cy, cx and d², |corr(S1, S2)| < 0.05; d is reported (no
prediction).

**S-C cropping.** The same fields read in the central and shifted 15 × 15 windows.
Report each pairwise error correlation and its change from S-B. No directional
prediction; this measures the windowing contribution.

**S-D collocation accuracy and ranking** (triplet R, S2, S1; full canvas and central
window). Criteria per variable: ranking matches truth (R < S1 < S2 or whatever the truth
order is), and |estimate − truth| ≤ 0.05 (cy, cx) or ≤ 0.10 (d, d²). Prediction: on the
full canvas, cy, cx and d² meet both criteria; d is reported. Cropped results reported.

**S-E planted shared displacement** (one translation per window applied to S1 and S2,
σ_sh ∈ {0.5, 1.0}). Prediction: on the full canvas the centroid ranking inverts (estimate
S1 < R while truth R < S1) at σ_sh = 1.0, as in T1.

**S-F known-dependence benchmark** (five observers R, S1, S2, Gs, H; full canvas; cy).
Scenarios: none planted; S1–S2 shared displacement (σ_sh = 1.0); R–Gs shared displacement
(σ_sh = 1.0). Each scenario split into 10 blocks of 2,000 windows; `tetrad_audit`
(domain-qualification/independence_audit.py) run per block with 200 resamples.
Report per scenario: share of blocks flagging the planted pair (hit rate) and share
flagging any pair when nothing is planted (false-alarm rate).

## Reading

- S-B ≈ 0 and S-C > 0 for S1–S2: the T1 baseline shared error came from windowing;
  F6 is replaced by a windowing finding, and the real-data window control gains priority.
- S-B ≈ 0 and S-C ≈ 0: the baseline came from the non-conserving blur; F6 is withdrawn.
- S-B > 0: smoothing itself correlates coordinate errors; F6 is reinstated in corrected form.
- S-E inverts: T1's mechanism-class support for F3 is restored under the corrected control.
- S-D decides the dispersion variable (d vs d²) for any future use.
