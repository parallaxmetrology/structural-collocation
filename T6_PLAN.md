# T6: do the two satellites share an error by keying on convective cores? (frozen 2026-10-15, before any T6 value was computed)

From T5: the real ranking (microwave < radar) needs an error shared by the two
satellite observers; smoothness is not sufficient and GOES parallax is not shared.
Candidate tested here, **H_core**: both satellite retrievals are pulled toward the
convective core of a rain system (cold tops, strong ice scattering), while radar
sees the full surface rain area, including light and stratiform rain. Relative to
surface rain, the satellites would then err in the same direction.

## Data and sample

Probe A2 triplet (radar, GOES IR, microwave), unblurred, same windows rule, same
cluster bootstrap (500 resamples, UTC day × 5° tile), centroid y and x only (T1:
the reliable variables). Dispersion reported descriptively.

## Definitions (fixed now)

- **Core** of an observer's window: cells at or above the 80th percentile of that
  observer's positive cells in the window (all positive cells if fewer than 5).
  Core field = field × core mask. Unit-free, so it applies to the cold-cloud field too.
- **Light-rain fraction** of a window: share of radar rain mass in cells below
  5 mm/h (0.1° hourly block mean).

## Predictions under H_core

- **P1 (core reading):** with every observer read on its own core, the inversion
  nf_R − nf_MW on centroid y and x shrinks: the upper end of its 95% CI falls below
  the full-area point estimate (+0.224 y, +0.193 x, from T5 σ = 0).
- **P2 (stratification):** the full-area inversion grows with light-rain fraction:
  (top tercile − bottom tercile) of nf_R − nf_MW > 0 with a joint cluster-bootstrap
  95% CI above 0, on centroid y and x.
- **P3 (geometry):** in windows where radar's core and full-area centroids differ by
  ≥ 1 cell, both satellites' centroids are closer to the radar **core** centroid than
  to the radar **full-area** centroid: median paired difference
  |S − R_full| − |S − R_core| > 0 with a bootstrap 95% CI above 0, for GOES IR and
  for microwave.

## Interpretation

All three hold: core keying is supported as the shared error. None holds: H_core is
not supported, and the shared error stays unidentified (ice aloft remains). Mixed:
reported as mixed, with the specific failures named.
