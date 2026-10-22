# Audit: independent re-implementation of the July structural-collocation numbers

Date: 2026-10-15. Scope: Probe A2 (area triplet), T7 (point-support triplets, error correlation), T7-R (tetrads and
triplet table). Point estimates only; bootstraps were not re-run.

## Method

1. Wrote `reimpl.py` (sample construction) and `analyze.py` (estimators) from `PROBE_A_PLAN.md`, `PROBE_A2_PLAN.md`,
   `T7_PLAN.md`, `T7R_PLAN.md` and the raw `data/` files (`mrms`, `goes`, `imerg_x`, `imerg`, `gauges.csv`,
   `stations.csv`). The existing scripts and `*_windows.csv` / `t7_points.csv` were not opened until the numbers below
   existed. The grid origin (55.0°N, 130.0°W, 0.1° cells) was taken from the MRMS domain description in
   `FEASIBILITY_PLAN.md`.
2. Compared to `data/probe_a2_results.json`, `data/t7_results.json`, `data/t7r_results.json`
   (`compare.py`, full table in `match_table.md`, 101 quantities).
3. Only then read `scripts/probe_a.py`, `probe_a2.py`, `t7.py` and explained each difference.

`t5.py`, `t7_robust.py`, `t7s.py` were not read; `t7s` (the ρ_RG sweep) was out of the requested scope and was not
re-implemented. `t5.windows` / `t5.load_hour` (used by `t7.py`) were not opened either.

## Result

Final version of the re-implementation: sample sizes identical (Probe A2 n = 14,585 windows; T7 n = 2,402 window-hours,
1,201 weak and 1,201 strong; median radar peak 15.833), and **all 101 compared quantities agree with the stored
results to within 3.5e-9** (largest: T7 weak cy, RGM nf_MW, 0.06345299520 vs 0.06345299272).

| Item | Quantities compared | Max abs difference | Status |
|---|---|---|---|
| Probe A2 nf_R, nf_GIR, nf_MW and nf_R − nf_MW, six variables (wet, cy, cx, disp, isl, logmean) | 24 | < 1e-8 | match |
| T7 triplet (R, GIR, MW) and (R, MW, G), weak and strong, cy and cx: nf of each member, nf_R − nf_MW | 32 | < 1e-8 | match |
| Error correlation c(GIR, MW) with R and G as instruments | 4 | < 1e-8 | match |
| T7-R triplets (R, GIR, G) and (GIR, MW, G), all members | 24 | < 1e-8 | match |
| T7-R R2 difference nf_R(RIRG) − nf_R(RMG) | 4 | < 1e-8 | match |
| Tetrads t1, t2, t3 | 12 | < 1e-8 | match |
| Median radar peak | 1 | < 1e-8 | match |

Headline values reproduced (weak stratum, cy / cx): nf_R − nf_MW in (R, GIR, MW) = +0.195 / +0.203; in (R, MW, G)
= −0.196 / −0.259; c(GIR, MW) = +0.350 / +0.361; t3 = +0.014 / +0.012, t1 = +0.120 / +0.137, t2 = +0.134 / +0.148.

## What did not match at first, and why

The first pass (kept as `reimpl_first_pass.py.txt`) matched the specs as written but differed from the stored numbers
by up to 0.0049 on a noise fraction. Each difference traced to a detail the written plans leave open. None is a
mistake in the original code relative to the plans; all four are choices the plans did not state.

| # | Difference | Effect in first pass | Cause |
|---|---|---|---|
| 1 | Wet threshold: first pass used rain > 0.5 mm/h; original uses ≥ 0.5 (`probe_a.reads`) | Probe A2 `isl` nf_R +0.0049, nf_MW −0.0026; `wet` nf_MW −0.0007 | Plans say "cells ≥ 0.5 mm/h" for `wet` and leave the island threshold implicit. Reading plan literally gives ≥; first pass mis-coded it. Hourly MW means (two half-hours) fall exactly on 0.5 often enough to matter. |
| 2 | Float precision of the window-sum test (sum ≥ 1 mm): first pass in float64, original in float32 (merged product `M` is the mean of two float32 files) | 1 window (2024-07-15 15:00, window 20,33) excluded in first pass: M sum is 0.99999998 in float64 and exactly 1.0 in float32 | The plan states the ≥ 1 mm rule, not the precision. |
| 3 | T7 sample: first pass applied the Probe A2 conditions including M (merged IMERG complete and sum ≥ 1); original T7 does not use M | 1 point (2024-07-15 21:00, window 9,29; M sum 0.64) excluded; this shifted the median and moved the weak/strong boundary, giving weak-stratum differences up to 0.0016 (e.g. `cx` RGM nf_MW 0.0882 vs 0.0892; `cy` RGM nf_IR 0.5131 vs 0.5120) | T7_PLAN.md says "Probe A2 windows (radar, GOES IR, microwave; unblurred)", which can be read either way. The original omits M, which is sensible since M is not an observer in T7. |
| 4 | Stratum tie rule: original weak = radar peak ≤ median, strong = > median | none after #3, since no ties straddle the median here | Plan says "below vs above the median" without a tie rule. |

After corrections 1 to 3 the remaining differences are rounding only (< 4e-9; the MW hourly mean must be taken in
float32 as stored to reproduce exact threshold ties).

## Things confirmed by the re-implementation

- Triple-collocation formulas (error variance C_xx − C_xy C_xz / C_yz, noise fraction 1 − C_xy C_xz / (C_xx C_yz)) with
  sample covariance, as in the plans.
- The "point-centroid" in T7 is the gauge-value-weighted (observer-value-weighted) mean of integer cell offsets from the
  window centre (7, 7); a constant shift would not change any covariance in any case.
- c(GIR, MW) = [C(I,M) − C(I,R) C(M,G)/C(R,G)] / √(e_I e_M), with e_x = C_xx − C_xR C_xG / C_RG. This normalisation is
  not spelled out in T7_PLAN.md ("normalised"); this choice reproduces the stored value.
- Station-to-cell assignment by floor((55 − lat)/0.1), floor((lon + 130)/0.1) reproduces the exact T7 sample, which
  also supports the assumed grid origin.

## Not covered

- Cluster-bootstrap CIs, share of negative error variances, and the stop/stability rules (point estimates only).
- Quadruple check for M (Probe A) and the T7-S ρ_RG sweep.
- This is a re-implementation by a second reader of the same written plans and raw data. It checks that the code does
  what the plans say, and that the arithmetic and sample construction are consistent. It does not test the plans'
  assumptions (additive model, independence of radar and gauge errors, correctness of the raw-data fetch and the
  grid/timing alignment, hour assignment of gauge reports).
- Since the first-pass ambiguities (#1 to #4) moved results by up to 0.005, the August replication (T10) is more
  exactly specified if the wet threshold (≥ 0.5), the float32 handling, the T7 sample (no M), and the tie rule (≤ median
  is weak) are written into `T10_PLAN.md` before it runs. The existing scripts already do this; the plan text does not.

## Files

`reimpl.py` (samples), `analyze.py` (estimators), `compare.py` (comparison), `match_table.md` (all 101 values),
`reimpl_results.json`, `a2_windows.csv`, `t7_points.csv` (re-implementation outputs), `reimpl_first_pass.py.txt`
(superseded first version, kept for the record).
