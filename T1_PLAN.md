# T1: synthetic anchor for structural collocation (frozen 2026-10-15, before any synthetic value was computed)

House rule (ISR I-20): never quote a read without a machine-perfect control of the
same class. Probes A and A2 produced two claims that have no control yet:
F2 (structural collocation is stable where value collocation is not) and F3
(the persistent ranking MW < R is the signature of a shared satellite error).
T1 builds rain-like windows with **known truth** and simulated observers with
**planted** noise and **planted** error covariance, and asks whether the
collocation estimators recover what was planted.

## Truth fields

- 15 × 15-cell windows (as in Probes A/A2); 20,000 windows per scenario; seed 20261015.
- Each window: K = 1 + Poisson(1.5) Gaussian rain cells; centre uniform in
  [−5, 5]² cells; size s ~ U(0.8, 2.5) cells; peak intensity ~ lognormal(μ = 1.5, σ = 0.8)
  mm/h. Truth T = sum of cells; values < 0.1 set to 0.

## Simulated observers (each re-renders the same cells, so mass is not lost to blurring)

| observer | role (real analogue) | per-cell position error | size added (blur) | per-cell intensity error | shared window error |
|---|---|---|---|---|---|
| R | radar | N(0, 0.5) | 0 | ×lognormal(0, 0.3) | none |
| S1 | microwave-like | N(0, 0.8) | +1.0 | ×lognormal(0, 0.4) | + d_sh |
| S2 | infrared-like | N(0, 1.5) | +2.0 | ×lognormal(0, 0.6) | + d_sh |
| M | merged-like | 0.6·S1 + 0.4·S2 (same window) | | | inherits S1, S2 errors |
| G | independent surface-like (for P4 only) | N(0, 1.0) | +1.5 | ×lognormal(0, 0.5) | none |

d_sh = one displacement vector per window, N(0, σ_sh) in each axis, applied to
every cell of **both** S1 and S2: the planted shared error (a stand-in for "both
satellites see the same ice aloft, displaced from the surface rain").

Scenarios: σ_sh ∈ {0 (S0, all errors independent), 0.5, 1.0, 1.5} cells.

## Reads and estimators

Identical to Probe A (`scripts/probe_a.py`): `reads()` on every observer and on the
truth (wet ≥ 0.5, centroid y/x, dispersion, islands, log mean); same inclusion
filter (every observer window sum ≥ 1); TC/ETC noise fractions on (R, S2, S1) in
the (R, IR, MW) slots; quadruple check of M; 500 bootstrap resamples with
clusters of 20 consecutive windows.

**Ground truth for scoring:** true noise fraction of observer o on variable v =
1 − corr(x_o, x_T)²; true error correlation between o and p = corr of their OLS
residuals on x_T.

## Predictions

- **P1 (S0, independent errors):** structural TC recovers the true noise fractions
  within 0.05 (absolute) for centroid y and x, and within 0.10 for dispersion,
  islands and wet fraction, for all three observers.
- **P2 (S0):** value TC (log mean) is stable (< 20% negative bootstrap error
  variances) and within 0.05 of truth. If it fails here, value collocation fails
  on these fields even with every assumption met, and F2's contrast is a property
  of the variable type, not of the real-data dependence.
- **P3 (planted shared error):** as σ_sh grows, TC underestimates S1's noise and
  overestimates R's, and at σ_sh ≥ 1.0 the TC ranking puts S1 below R on centroid
  y and x although the true ranking is R < S1. (This is the pattern seen in A/A2.)
- **P4 (detection):** with the independent fifth observer G as an instrument,
  quadruple collocation c(S1, S2) = C(S1,S2) − C(S1,R)·C(S2,G)/C(R,G) recovers the
  planted S1–S2 error correlation within 0.10 on centroid y and x, and is within
  0.05 of 0 at σ_sh = 0.
- **P5 (lineage, as in A/A2):** the quadruple check reports clearly positive
  M–S1 and M–S2 error correlations in every scenario (M is built from both).

## Stop / interpretation rule

If P1 fails on centroids (the simplest case, all assumptions met), structural
collocation is not a valid estimator on windows of this kind, and F2/F3 are
withdrawn as findings. If P1 holds and P3 holds, the A/A2 ranking is reproducible
by a planted shared satellite error; this supports F3's mechanism class without
proving the ice-aloft cause.
