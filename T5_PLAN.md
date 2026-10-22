# T5: shared smoothness vs shared physics (frozen 2026-10-15, before any T5 value was computed)

From T1 (`T1_RESULTS.md`): the real ranking microwave < radar < infrared is the
signature of an error shared by the two satellite observers. Two causes fit:
**shared physics** (both respond to ice aloft) and **shared smoothness** (both are
coarse relative to radar, and T1 showed smoothing alone correlates structural
errors). T5 makes radar as smooth as the microwave field and asks whether the
ranking survives.

## Observers and variables

The Probe A2 triplet (R, GOES IR, MW), the physically independent one. Variables:
**centroid y, centroid x, dispersion only** (T1 validated these; islands and wet
fraction are excluded). Same windows rule, estimators and bootstrap as Probe A2.

## Radar coarsening

- Radar 0.1° field smoothed with a Gaussian kernel of width σ cells
  (normalised convolution, so missing cells do not bleed zeros), **before**
  windowing and before the inclusion filter.
- **Ladder:** σ = 0, 0.5, 1.0, 1.5, 2.0, 3.0 cells (descriptive sweep,
  resolution-floor style).
- **Matched level σ\* (primary):** the ladder value (in steps of 0.25, from 0 to 3)
  whose radar median lag-1 spatial autocorrelation inside sample windows is
  closest to the microwave field's. Computed from field properties only (no
  collocation output), on every 4th hour of the Probe A2 sample.

## Synthetic calibration arm (T1 fields, before the real run)

Re-run T1 with radar-like R blurred (size + 0, + 1, + 2), at σ_sh = 0 and σ_sh = 1.0.
This shows what each cause predicts in a world where the answer is known:
- planted shared displacement (σ_sh = 1.0, a "physics" stand-in): the inversion
  S1 < R is expected to **persist** when R is blurred, because blurring R does not
  give it the shared displacement;
- smoothness only (σ_sh = 0): T1 showed no inversion; the arm reports whether
  blurring R creates or removes one.

## Predictions (real data, at σ*)

- **If shared smoothness drives the real ranking:** nf_R − nf_MW on centroid y and
  x falls to ≤ 0 (radar no longer noisier than microwave), i.e. the 95% CI of the
  difference includes or lies below 0.
- **If shared physics drives it:** nf_R − nf_MW stays clearly > 0 (CI above 0) at σ*.
- Either way, report the full ladder curve of nf_R, nf_MW, nf_GIR for the three
  variables.

## Interpretation rule

The two outcomes are read against the synthetic arm. Persistence at σ* does not
prove ice aloft; it removes smoothness as a sufficient explanation. Disappearance
at σ* makes smoothness sufficient and leaves physics unnecessary for this pattern.
