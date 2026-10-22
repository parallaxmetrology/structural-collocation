# Feasibility plan (frozen 2026-10-15, before any predictor or error value was computed)

Exploratory feasibility, not the confirmatory test. Its job is to decide whether a
pre-registered test is worth running, and to fail fast if the structural read adds
nothing over plain baselines.

## Data

- Period: 504 hours, accumulations ending 2024-07-01 01:00 UTC through
  2024-07-22 00:00 UTC.
- **Discovery:** hours ending 1–10 July (incl. 11 July 00:00). **Validation:**
  the rest. Anything fitted is fitted on discovery only; validation is read once.
- Radar: MRMS `RadarOnly_QPE_01H` (not gauge-adjusted), top-of-hour files, 0.01°.
- Incumbent: MRMS `RadarAccumulationQualityIndex_01H` (RQI, 0–1).
- Satellite: IMERG Late V07 half-hourly; hourly = mean of the two half-hours
  (mm/h). Grid 0.1°, cropped to the MRMS domain.
- Gauges (label only): USCRN hourly `P_CALC` (mm), CONUS; ASOS hourly `p01i`
  (inches → mm) from IEM, routine reports at minute 50–59 assigned to the hour
  ending at the next top of hour.

## Grids

MRMS 0.01° is block-averaged 10 × 10 to 0.1° (aligned with IMERG cell centres);
a block needs ≥ 80% valid radar pixels, else missing. Radar-only and RQI are
also sampled at the gauge's own 0.01° pixel.

## Sample

A gauge-hour enters when: the gauge value exists; the radar pixel at the gauge is
valid; the 15 × 15 window of 0.1° cells centred on the gauge cell (±0.7°) has
≥ 90% valid radar blocks and no missing IMERG; and **both** sensors have rain in
the window (window sum ≥ 1 mm in each, so centroids are defined).

## Target

E = |radar-only at gauge pixel − gauge| (mm).

## Predictors

- Incumbent: `1 − RQI` at the gauge pixel.
- Naive 1: |R − S| at the gauge's 0.1° cell.
- Naive 2: mean |R − S| over the window.
- Magnitude controls: R at the gauge pixel, S at the gauge cell (known at
  prediction time).
- **Ω_loc (primary structural):** distance in cells between the intensity-weighted
  centroids of R and S in the window.
- **Ω_pat (secondary structural):** 1 − Pearson r between R and S over the window.

## Decision metric (validation block)

Rank-transform everything. On discovery fit OLS: E ~ baseline (incumbent + naive
1 + naive 2 + magnitude controls), and E ~ baseline + Ω_loc + Ω_pat. On
validation compute Spearman ρ(predicted, E) for each; Δρ = full − baseline. CI by
cluster bootstrap over (UTC day × 5° tile) blocks, 2,000 resamples.

- **Go** to a pre-registered full-season test: Δρ ≥ 0.02 and the 95% CI excludes 0.
- **Stop:** otherwise. Report the result as the operating window it maps.

Also reported (descriptive, not decision): Spearman of each predictor alone and
its partial Spearman controlling for the magnitude controls; separately for USCRN
and ASOS.

## Gate 6 (independence)

Signal check: share of sample windows where each sensor's window is non-constant.
Error independence: |Pearson r| between (R_cell − G) and (S_cell − G) at gauges;
PASS if < 0.30. `independence_audit.py` also run on the predictor set to show
redundancy among predictors.

## Known limits, stated in advance

Point gauge vs 1 km / 10 km cells adds noise to every predictor equally. 21 days
of one summer is one regime (warm-season convection). ASOS coordinates can be off
by a grid box. IMERG Late is smooth below ~200 km, so Ω_loc inside a ±0.7°
window may be noise-limited: that is part of what this checks.

---

## Amendment 1 (2026-10-15, before any target E or predictor value was computed)

Reason: Russell asked that the structural read not be so crude it misses the
signal, and pointed to ZTOYBOX (Saltimbanques) and Pathology_new (Core42, serial-
section persistence). Only data downloads and two grid/timing alignment checks
(radar vs satellite field correlation; no gauge or error values) had run.

Added structural reads (all R vs S unless noted; wet = ≥ 0.5 mm/h; dry = < 0.1 mm/h):

| read | borrowed from | definition |
|---|---|---|
| `om_iou` | Saltimbanques region overlap | 1 − Jaccard of the two wet masks |
| `om_islands` | Saltimbanques islands / n_groups | abs difference in number of wet connected components (≥ 2 cells, 8-connected) |
| `om_disp` | Canon-3 SDI / `field_dispersion` | abs difference in mass-weighted RMS distance from own centroid, / window diagonal |
| `om_void` | `field_void` | abs difference in dry fraction |
| `om_scale` | Core42 multiscale (g4_s1/s8/s32), persistence under resolution | 1 − Pearson r after 3 × 3 box smoothing |
| `om_scale_rec` | same | r(3 × 3) − r(1 × 1): how much disagreement resolves at coarser scale |
| `om_edge` | Saltimbanques Dependency (local load-bearing) | abs difference in distance from gauge cell to nearest wet cell (capped at 8) |
| `om_pat_persist` | HTAN serial-section persistence | mean of `om_pat` at t and t − 1 (causal) |
| `om_loc_ir`, `om_pat_ir` | Ω with a third observer | as `om_loc` / `om_pat` but radar vs IMERG infrared-only estimate (geostationary IR, physically distinct from radar and from microwave) |

Microwave-only (IMERG `MWprecipitation`) covers ~20% of cells per half-hour; it is
used only in a descriptive subset (windows with ≥ 80% MW coverage), not in the
decision model.

**Decision model (headline):** FULL_A1 = baseline + all 12 structural reads
(`om_loc`, `om_pat` and the 10 above), same rank-OLS, same discovery/validation
split, same Δρ ≥ 0.02 with cluster-bootstrap 95% CI > 0 rule. The original
two-read model is still computed and reported as "frozen v0". Both results are
reported whichever way they fall.

---

## Amendment 2 (2026-10-15, still before any target E or predictor value was computed)

Reason: Russell asked to also consult ZTOYBOX seismology, Material and
semiconductor. Three program lessons apply:

1. **Variance of disagreement is the durable unit**
   (`ZTOYBOX/Material/DETERMINATION_variance_of_disagreement.md`): in earlier work
   the mean of Ω did not carry signal while its spread did. Added:
   - `om_sigma_space`: std, over the 9 non-overlapping 5 × 5 sub-windows of the
     15 × 15 window, of local pattern disagreement (1 − r within each sub-window;
     sub-windows where either sensor is constant are skipped; needs ≥ 3 valid).
   - `om_sigma_time`: std of `om_loc` over hours t − 2, t − 1, t at the same gauge
     (causal; needs all three).
2. **Resolution floor** (`ZTOYBOX/Material/RESOLUTION_FLOOR_PROTOCOL.md`): predicted
   sign of every structural read vs E is **positive** (more disagreement → larger
   radar error). Descriptive sweep of `om_pat` and `om_loc` over window half-widths
   3, 5, 7, 10 cells (0.3° to 1.0°), reporting ρ(g) on the discovery block, with
   the sample fixed to gauge-hours valid at the largest window.
3. **Seismology guardrails** (`ZTOYBOX/seismology/GUARDRAILS.md`): (G1) report the
   share of near-zero Ω values and confirm they are not computation failures;
   (G4) the null must break spatial structure, not timing: **spatial-scramble
   null** = within each hour, reassign satellite windows to a different gauge's
   location (radar, gauge and RQI kept), recompute the structural reads, refit
   on discovery, report Δρ on validation. It should collapse toward 0.

Considered and not added: semiconductor NILS/LER (edge sharpness / roughness).
Satellite edges are smooth by resolution, so a sharpness disagreement would mostly
re-measure the resolution gap; noted for a radar-vs-radar extension.

**Headline decision model is now FULL_A2** = baseline + the 12 Amendment-1 reads
+ `om_sigma_space` + `om_sigma_time` (14 structural reads). Decision rule unchanged.
Also reported: frozen v0 (2 reads), and a variance-only model (baseline +
`om_sigma_space` + `om_sigma_time`), since that is the program's thesis.
