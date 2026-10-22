# Feasibility results (2026-10-15)

Run exactly as frozen in [FEASIBILITY_PLAN.md](FEASIBILITY_PLAN.md) (v0 + Amendments 1
and 2). Code: `scripts/fetch_*.py`, `scripts/analyze.py`. Numbers: `data/results.json`
(git-ignored with the rest of `data/`).

## Sample

- 288,690 gauge-hours where both radar and satellite carried rain in the 15 × 15
  window: 132,755 discovery (1–10 July 2024), 155,935 validation (11–21 July).
- 271,511 ASOS, 17,179 USCRN. 423 discovery and 533 validation clusters
  (UTC day × 5° tile).
- Signal check: every sample window was non-constant for both sensors. Near-zero
  `om_loc` (< 0.01 cells): 0.002% of rows (seismology guardrail G1 satisfied).
- Coverage: infrared reads 81%, `om_sigma_space` 81%, `om_sigma_time` 72%,
  microwave reads 39%.

## Decision (validation block, read once)

Baseline = MRMS radar quality index (1 − RQI) + |R − S| at the gauge cell + mean
|R − S| over the window + radar value at the gauge + satellite value at the gauge
cell. Validation Spearman ρ of the baseline prediction with radar error: **0.685**.

| Model | ρ (validation) | Δρ vs baseline | 95% CI (cluster bootstrap) | Rule (Δρ ≥ 0.02, CI > 0) |
|---|---|---|---|---|
| frozen v0 (+ `om_loc`, `om_pat`) | 0.6855 | +0.00002 | [−0.0003, +0.0003] | not met |
| variance only (+ `om_sigma_space`, `om_sigma_time`) | 0.6868 | +0.0014 | [+0.0009, +0.0019] | not met |
| **FULL_A2** (+ all 14 structural reads), headline | 0.6869 | **+0.0014** | [+0.0008, +0.0021] | **not met** |

**Spatial-scramble null** (satellite windows reassigned to another gauge's location
within the same hour):

| Model | Δρ under null | 95% CI |
|---|---|---|
| frozen v0 | +0.0003 | [−0.0005, +0.0011] |
| variance only | +0.0031 | [+0.0019, +0.0042] |
| FULL_A2 | +0.0031 | [+0.0019, +0.0043] |

The null increment is larger than the real one. The small positive Δρ in the
real run is therefore not carried by radar–satellite spatial agreement; it is
consistent with the temporal and missingness terms (lagged reads, coverage
indicators), which the scramble does not remove.

**Outcome: stop rule met.** No structural read adds held-out information about
hourly radar error at gauges beyond the baseline at this scale and label.

## What carries the error (discovery block, descriptive)

Partial Spearman with radar error after controlling for the radar and satellite
values at the gauge (magnitude):

| Read | partial ρ (ALL) | partial ρ (USCRN) |
|---|---|---|
| 1 − RQI (MRMS incumbent) | −0.019 | +0.030 |
| \|R − S\| at gauge cell | −0.207 | −0.209 |
| mean \|R − S\| in window | −0.033 | −0.043 |
| `om_loc` / `om_pat` | −0.013 / −0.015 | −0.005 / −0.002 |
| `om_iou` | −0.001 | +0.010 |
| `om_edge` | −0.041 | −0.027 |
| `om_sigma_space` | −0.004 | −0.020 |
| `om_sigma_time` | +0.002 | +0.032 |
| `om_loc_ir` / `om_pat_ir` | +0.003 / +0.004 | +0.018 / +0.033 |
| `om_loc_mw` / `om_pat_mw` | −0.009 / −0.013 | +0.010 / +0.015 |

All structural reads have |partial ρ| ≤ 0.04, and most are negative against the
pre-registered positive sign. MRMS's own quality index also carries no
information beyond magnitude here. Raw (unpartialled) correlations of the
structural reads are mostly negative (e.g. `om_iou` −0.37): larger disagreement
co-occurs with lighter, scattered rain, where absolute errors are small.

## Resolution-floor sweep (discovery, fixed sample n = 98,751, partial ρ)

| half-width | 0.3° | 0.5° | 0.7° | 1.0° |
|---|---|---|---|---|
| `om_pat` | −0.028 | −0.024 | −0.017 | −0.013 |
| `om_loc` | −0.027 | −0.018 | −0.014 | −0.009 |

No rung is correctly signed. No resolution floor exists in 0.3°–1.0° for this
label; the curve trends toward zero as the window grows.

## Gate 6 (independence)

- Radar error vs satellite error at gauges: Pearson r = 0.47, Spearman 0.40 →
  **PARTIAL** (0.30–0.60). Both errors are measured against the same gauge, so
  the gauge's own representativeness error is shared and inflates this r.
- Radar vs satellite at the gauge cell: r = 0.51.
- Predictor set: I = 0.196 (PASS) with one near-duplicate pair, `om_pat` ~
  `om_scale` r = 0.95 (by construction: the same pattern correlation before and
  after 3 × 3 smoothing), then `om_pat` ~ `om_pat_persist` 0.86.

## Operating window this maps

- Label: hourly absolute error of 1 km radar-only QPE at a point gauge, warm
  season, CONUS. This error is dominated by rain magnitude (baseline ρ 0.685,
  mostly from the radar value itself) and by point-vs-pixel noise. Neither the
  incumbent quality index nor any structural read separates beyond magnitude.
- Scale: 0.1° satellite windows of 0.3°–1.0° half-width.
- The rubric's Gate 4 for framing B was scored PARTIAL at the desk ("a point
  error is not an arrangement"); this run measures it: at this label the
  structural reads have nothing to track.

Not tested here: accumulations longer than one hour (daily or event totals,
where point noise averages down); gauge-network areal averages as the label;
cool-season or orographic regimes; radar-vs-radar or radar-vs-gauge-analysis
disagreement.
