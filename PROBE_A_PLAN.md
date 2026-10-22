# Probe A: structural collocation (frozen 2026-10-15, before any collocation value was computed)

From `RESEARCH_SYNTHESIS.md` §3. No gauges are used.

## Observers (four, hourly, 0.1° grid, 504 hours of July 2024 already in `data/`)

- **R**: MRMS radar-only QPE, 10 × 10 block mean (`data/mrms`).
- **IR**: IMERG Late `IRprecipitation`, mean of the two half-hours (`data/imerg_x`).
- **MW**: IMERG Late `MWprecipitation`, nan-mean of the two half-hours.
- **M**: IMERG Late merged `precipitation`, mean of the two half-hours (`data/imerg`).

R, IR and MW are physically different measurements (radar backscatter; cloud-top
infrared; microwave emission/scattering). M is built from IR and MW, so its errors
are expected to correlate with theirs.

## Sample

Non-overlapping 15 × 15-cell windows tiling the grid (rows 2–346, cols 2–696).
A window-hour enters when: R has ≥ 90% valid cells; IR and M are complete; MW has
≥ 80% valid cells (missing MW cells set to 0 only where valid share ≥ 80%); and
every observer's window sum ≥ 1 mm (so centroids exist).

## Structural variables (one scalar per observer per window-hour)

- `wet`: fraction of cells ≥ 0.5 mm/h
- `cy`, `cx`: intensity-weighted centroid (cells from window centre)
- `disp`: mass-weighted RMS distance from own centroid (cells)
- `isl`: number of wet connected components ≥ 2 cells (8-connected)
- `logmean` (value control, the incumbent TC target): log(window mean + 0.1)

## Estimators (additive model x_i = a_i + b_i T + e_i)

- **Triple collocation / ETC** (McColl et al. 2014) for triplet (x, y, z):
  error variance e_x = C_xx − C_xy C_xz / C_yz; squared correlation with truth
  ρ²_x = C_xy C_xz / (C_xx C_yz). Reported as **noise fraction** 1 − ρ² (scale-free).
- **Primary triplet:** (R, IR, MW), the only triplet of physically independent observers.
- **Quadruple check:** error cross-covariance of M with IR and with MW, using the
  other two independent observers as instruments:
  c(M, IR) = C(M, IR) − C(M, R) C(IR, MW) / C(R, MW), and likewise c(M, MW) with
  instruments R and IR; reported as error correlations.
- **Uncertainty:** 500 cluster-bootstrap resamples over (UTC day × 5° tile);
  report median and 95% interval, and the share of resamples where any error
  variance is negative.

## Predictions (stated before computing)

1. **Radar has the smallest noise fraction** of R, IR, MW, for most structural
   variables (radar is the highest-resolution observer).
2. **IR noise fraction > MW noise fraction** (IR sees cloud tops, not rain).
3. **The quadruple check flags M as error-correlated with its own inputs.**
   Clarification, written before computing: because the sample requires ≥ 80% MW
   coverage, M in these windows is dominated by its microwave input, so the
   largest flagged error correlation is expected for **M–MW**, with M–IR positive
   but smaller. If neither is clearly positive, the method is not seeing a
   dependence that is known by construction, and the probe stops.

## Stop rule

Stop and record the operating window if, for the majority of the five structural
variables, the (R, IR, MW) estimates are unstable: > 20% of bootstrap resamples
give a negative error variance for any observer, or the 95% interval of any noise
fraction spans more than 0.5.

## What a positive result would mean (and not mean)

Positive = stable structural error estimates without truth, with the known
dependence detected. It would establish structural collocation as a working
measurement, not a product; the next question would be whether structural noise
fractions carry information the value-based TC does not (descriptive here).
