# T1 results: synthetic anchor for structural collocation (2026-10-15)

Run as frozen in [T1_PLAN.md](T1_PLAN.md). Code `scripts/t1_synthetic.py`
(confirmatory) and `scripts/t1_explore.py` (exploratory, after scoring). Numbers in
`data/t1_results.json`, `data/t1_explore.json`. 19,997 of 20,000 windows passed the
inclusion filter in every scenario.

## Predictions, scored

| | prediction | result |
|---|---|---|
| P1 | S0: TC within 0.05 of truth on centroids; within 0.10 on dispersion, islands, wet | **met for centroid y/x** (max error 0.034) **and dispersion** (0.079); **not met for islands** (radar 0.55 off, S1 0.47) **or wet fraction** (0.12–0.17 off; 100% negative bootstrap variances) |
| P2 | S0: value TC (log mean) stable and within 0.05 | **stable** (0% negative); accuracy **not met** (0.049–0.063 off) |
| P3 | shared error makes TC underrate S1 and overrate R; ranking flips at σ_sh ≥ 1.0 | **met**, and the flip already appears at σ_sh = 0.5 (centroid y: TC R 0.113 vs S1 0.082; truth R 0.045 vs S1 0.150) |
| P4 | with instrument G, recover planted S1–S2 error correlation within 0.10 on centroids; ≈ 0 at σ_sh = 0 | **partly**: within 0.10 at σ_sh 1.0 (0.43 vs 0.50) and 1.5 (0.63 vs 0.66); 0.15–0.21 low at σ_sh 0 and 0.5. The "≈ 0 at σ_sh = 0" premise was wrong: true error correlation is 0.16 there (see below) |
| P5 | lineage check positive for M–S1 and M–S2 in every scenario | **not met**: M–S1 undefined (negative error variance in the estimator), M–S2 out of bounds (0.79–1.14) for centroids and dispersion; defined and close on islands only |

Stop rule: P1 holds on centroids, so structural collocation is not withdrawn.

## What the anchor showed that the plan did not anticipate

**Smoothing is itself a shared error.** At σ_sh = 0 every per-cell error is
independent, yet the two blurred observers' structural errors correlate (centroid
0.16, dispersion 0.34, wet fraction 0.54). Blurring pulls both observers' reads
the same way relative to the sharp truth.

- **E1 (exploratory):** with a **sharp** instrument G (no added blur), quadruple
  collocation recovers that baseline almost exactly: centroid y 0.14 vs true 0.16 at
  σ_sh = 0; 0.50 vs 0.50 at σ_sh = 1.0. The P4 shortfall came from G sharing the
  blur.
- **E2 (exploratory):** giving the merged observer its own processing noise does not
  rescue the lineage check: M–S1 stays undefined because its instrument (S2)
  shares error with S1; M–S2 is overestimated (0.86 vs true 0.69 on centroid y).
  The lineage check registers dependence in direction but not in magnitude when its
  instruments are themselves dependent.

## Consequences for the real-data findings

- **F2 (structure stable, value unstable) is revised.** Value collocation is stable
  on synthetic fields, so its instability on real data is a property of the real
  observers' amount relationships, not of value collocation itself. Islands and wet
  fraction are not valid collocation variables (nonlinear in resolution); the real
  results on those two variables are withdrawn.
- **F3 (shared satellite error) is strengthened in mechanism class and widened in
  cause.** A planted shared error reproduces the real ranking (S1 < R) at modest
  strength. Two candidate causes now fit: shared physics (ice aloft) and
  **shared smoothness** (both satellites are coarse relative to radar).
- **The real lineage numbers (merged–MW 0.26–0.37)** are directional evidence only;
  their magnitude is not trustworthy when the instruments are dependent.

## Next test implied (T5)

Coarsen radar to the satellites' effective resolution before collocation. If the
real ranking is driven by shared smoothness, MW < R should disappear once radar is
equally smooth; if it persists, the shared-physics explanation gains.
