# Probe A results: structural collocation (2026-10-15)

Run as frozen in [PROBE_A_PLAN.md](PROBE_A_PLAN.md). Code `scripts/probe_a.py`;
numbers `data/probe_a_results.json`. One post-hoc check, labelled as such, in
`scripts/probe_a_posthoc.py` / `data/probe_a_posthoc.json`.

## Sample

25,468 window-hours (15 × 15 cells of 0.1°, July 2024) where radar (R), IMERG
infrared (IR), IMERG microwave (MW) and IMERG merged (M) all carried rain and MW
covered ≥ 80% of the window; 1,199 clusters (UTC day × 5° tile).

## Triple collocation on (R, IR, MW): noise fraction 1 − ρ² with the common signal

| variable | R | IR | MW | boot. share with a negative error variance | stable? |
|---|---|---|---|---|---|
| wet fraction | 0.349 [0.32, 0.38] | 0.611 [0.58, 0.64] | 0.007 [−0.03, 0.04] | 0.34 | no |
| centroid y | 0.254 [0.24, 0.27] | 0.538 [0.52, 0.56] | 0.070 [0.06, 0.08] | 0.00 | yes |
| centroid x | 0.290 [0.27, 0.31] | 0.567 [0.55, 0.58] | 0.122 [0.11, 0.14] | 0.00 | yes |
| dispersion | 0.570 [0.55, 0.59] | 0.778 [0.76, 0.79] | 0.278 [0.24, 0.31] | 0.00 | yes |
| islands | 0.677 [0.63, 0.71] | 0.932 [0.92, 0.95] | 0.628 [0.58, 0.67] | 0.00 | yes |
| log mean rain (value TC, incumbent target) | 0.408 [0.38, 0.43] | 0.609 [0.59, 0.63] | −0.028 [−0.05, −0.01] | 1.00 | no |

**Stop rule: not triggered** (4 of 5 structural variables stable). The value-based
target, the one TC is normally applied to, was unstable on the same windows.

## Quadruple check: error correlation of merged IMERG with its inputs

| variable | M–IR | M–MW |
|---|---|---|
| centroid y | 0.02 [−0.00, 0.05] | 0.26 [0.07, 0.39] |
| centroid x | −0.02 [−0.05, 0.01] | 0.36 [0.28, 0.44] |
| dispersion | −0.03 [−0.05, 0.00] | 0.31 [0.22, 0.39] |
| islands | 0.05 [0.03, 0.07] | 0.37 [0.32, 0.41] |
| wet fraction | −0.09 [−0.16, −0.03] | undefined (negative error variance) |

## Predictions

1. **Radar smallest noise fraction: not met (0 of 5).** MW is lowest on every
   structural variable.
2. **IR noise fraction > MW: met (5 of 5).**
3. **M flagged as error-correlated with its inputs, largest for M–MW: met for
   M–MW** (0.26–0.37, intervals exclude 0, on four variables). M–IR was predicted
   positive but smaller; it is ≈ 0.

## Post-hoc check (not pre-registered)

Hypothesis tested: MW looks best because R and IR share an hour-integrated
footprint that the MW snapshot lacks, biasing MW's error variance low. If so,
windows where MW covers both half-hours should show a *higher* MW noise fraction.

Observed: lower (3,487 both-half windows vs 21,981 one-half): centroid y −0.028
vs 0.085, centroid x 0.048 vs 0.132, dispersion 0.201 vs 0.285. **The footprint
explanation is not supported.** Alternative triplets give the same ranking (MW
noise 0.05–0.41 with instruments R and M).

## Reading

- Negative and near-zero MW error variances are not physical. In TC they arise
  when an error cross-covariance is present among the three observers. The
  pattern seen (R inflated, MW deflated) is what positive IR–MW error
  covariance produces: it inflates C(IR, MW), which lowers MW's (and IR's)
  estimated error and raises R's.
- A known construction route for that covariance, **verified in the IMERG V07
  ATBD** (`literature/IMERG_V07_ATBD_Huffman2023.pdf`): the IR estimate is now
  PDIR (PERSIANN Dynamic Infrared-Rain Rate, replacing PERSIANN-CCS), and the ATBD
  describes "the PMW calibration of PDIR" through a rotating 30-day accumulator
  (trailing for the Late run). IMERG IR is therefore not independent of IMERG MW.
  The calibration is distributional (a Tb-to-rain-rate mapping updated over 30
  days), so whether it alone accounts for the hour-level error covariance seen
  here is not established; it does break the independence premise.
- What held: structural collocation was numerically stable where value
  collocation was not, and it detected the dependence of M on MW that exists by
  construction.

## Operating window this maps

Structural collocation on (MRMS radar, IMERG IR, IMERG MW) at 0.1° / 1.5° windows /
hourly: stable estimates, but the per-observer ranking is not interpretable
because the IMERG infrared estimate is (very likely) calibrated to the microwave.

## Next step implied

Replace IMERG `IRprecipitation` with an infrared observer that has no microwave
calibration: GOES-16 ABI band 13 (10.3 µm) brightness temperature, public on AWS
(NOAA Open Data), with structural reads on the cold-cloud field (e.g. Tb < 235 K).
Then (R, GOES-IR, MW) is a physically independent triplet and predictions 1–3 can
be re-tested as written. Raw brightness temperature has no microwave calibration, so it removes the
dependence at its source.
