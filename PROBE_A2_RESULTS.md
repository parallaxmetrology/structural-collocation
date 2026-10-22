# Probe A2 results: structural collocation with raw GOES infrared (2026-10-15)

Run as frozen in [PROBE_A2_PLAN.md](PROBE_A2_PLAN.md). Code `scripts/fetch_goes.py`,
`scripts/probe_a2.py`; numbers `data/probe_a2_results.json`.

## Data and sample

- GOES-16 ABI band 13, CONUS sector, 4 scans per hour, 502 of 504 hours. Missing,
  per plan (< 3 scans): 2024-07-02 18:00 (0 scans) and 2024-07-10 14:00 (2 scans).
- 14,585 window-hours, 965 clusters (fewer than Probe A's 25,468 because GOES
  requires complete coverage in the window and the CONUS sector thins west of
  ~120°W; the 130–120°W band is 44% valid).
- Test hour alignment: radar vs cold-cloud field r = 0.39 at zero shift, best
  (0.398) one cell shifted, consistent with geostationary parallax of high cloud
  tops. Not corrected (a property of the observer).

## Triple collocation on (R, GIR, MW): noise fraction 1 − ρ²

| variable | R | GIR | MW | boot. share negative | stable? | raw r: R–GIR / R–MW / GIR–MW |
|---|---|---|---|---|---|---|
| wet fraction | 0.334 | 0.644 | −0.007 | 0.66 | no | 0.49 / 0.82 / 0.60 |
| centroid y | 0.257 [0.24, 0.28] | 0.505 [0.48, 0.53] | 0.036 [0.02, 0.05] | 0.00 | yes | 0.61 / 0.85 / 0.69 |
| centroid x | 0.282 [0.26, 0.30] | 0.524 [0.50, 0.55] | 0.089 [0.07, 0.11] | 0.00 | yes | 0.58 / 0.81 / 0.66 |
| dispersion | 0.534 [0.50, 0.57] | 0.841 [0.82, 0.86] | 0.244 [0.19, 0.29] | 0.00 | yes | 0.27 / 0.59 / 0.35 |
| islands | 0.659 [0.61, 0.71] | 0.924 [0.91, 0.94] | 0.585 [0.52, 0.64] | 0.00 | yes | 0.16 / 0.38 / 0.18 |
| log mean (value control) | 0.408 | 0.682 | −0.034 | 0.98 | no | 0.43 / 0.78 / 0.57 |

**Stop rule: not triggered** (4 of 5 structural variables stable).

## Quadruple check (M with its inputs)

| variable | M–GIR | M–MW |
|---|---|---|
| centroid y | 0.01 [−0.03, 0.05] | −0.34 [−1.69, 0.10] (unstable) |
| centroid x | −0.06 [−0.10, −0.01] | 0.20 [−0.05, 0.33] |
| dispersion | −0.05 [−0.09, −0.02] | 0.27 [0.10, 0.38] |
| islands | −0.01 [−0.03, 0.02] | 0.33 [0.25, 0.38] |

## Predictions

1. **Radar smallest noise fraction: not met (0 of 5)**, same as Probe A.
2. **GIR > MW: met (5 of 5).**
3. **M–MW largest, M–GIR positive but smaller: partly met.** M–MW is clearly
   positive on 2 of 4 stable variables (dispersion, islands); M–GIR is ≈ 0 or
   slightly negative.
4. (Descriptive, stated in the plan) Negative error variances did **not** fall:
   `wet` 0.34 → 0.66, `logmean` 1.00 → 0.98.

## What changed and what did not

Replacing IMERG's microwave-calibrated infrared with raw GOES brightness
temperature left the ranking unchanged (MW < R < IR on every structural
variable). **The IMERG-IR calibration route is therefore not the cause** of the
Probe A ranking.

The arithmetic of the exactly-identified triplet shows what drives it. For
centroid y, ρ²_MW = r(R,MW)·r(GIR,MW)/r(R,GIR) = 0.85 × 0.69 / 0.61 ≈ 0.96, while
ρ²_R = 0.85 × 0.61 / 0.69 ≈ 0.75. MW wins because the cloud-top infrared agrees
more with MW (0.69) than with radar (0.61).

## Reading (hypothesis, not tested here)

Two documented facts and one physical expectation fit the pattern:

- IMERG MW estimates are intercalibrated to the GPM combined radar–radiometer
  product (CORRA) by quantile matching on a rolling 45-day basis (IMERG V07 ATBD,
  verified in `literature/`). This is distributional and should affect amounts
  more than positions.
- Over land, passive-microwave precipitation retrievals lean on high-frequency
  ice-scattering signals (desk knowledge, not verified in this folder). Cold
  cloud tops are also an ice-aloft signal. If both satellites respond to ice
  aloft rather than surface rain, their errors relative to surface rain are
  positively correlated. That inflates C(GIR, MW), which lowers MW's estimated
  error and raises radar's: the pattern observed.

If this holds, the triplet is **instrument-independent but not
mechanism-independent**. Two of the three observers measure the same
intermediate (ice aloft).

## Operating window this maps

Structural collocation of radar, geostationary IR and passive-microwave rain at
0.1° / 1.5° windows / hourly gives numerically stable error estimates for
centroid, dispersion and island count. The per-observer ranking cannot be read as
accuracy, because the satellite observers plausibly share an ice-aloft error.
Value-based collocation (log mean rain) is unstable on the same windows in both
probes.

## Lesson for the rubric (proposed, not yet applied)

Gate 6 currently treats "different sensors" as provisional independence. This
probe gives a counter-case: radar, infrared and microwave are different sensors,
yet two of them plausibly respond to the same physical intermediate. Proposed
addition to Gate 6: independence is judged at the level of the **mechanism each
observer responds to**, not the instrument.

## Options from here

- A surface-rain observer with no ice sensitivity as the third leg: an areal
  gauge average where networks are dense (e.g. state mesonets), or disdrometers.
  Coverage is the limit.
- Close the precipitation line here and score Probe B (convective mode) with the
  rubric.
