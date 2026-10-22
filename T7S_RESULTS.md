# T7-S results: sensitivity to radar–gauge error dependence (2026-10-15)

Run as frozen in [T7S_PLAN.md](T7S_PLAN.md). Code `scripts/t7s.py`; numbers
`data/t7s_results.json`. Sample: T7 point-support set (2,402 window-hours), weak and strong
halves, centroid y and x. 300 cluster-bootstrap resamples per grid point.

## The equivalence class

The maximum-likelihood fit is identical at every assumed ρ_RG (e.g. discrepancy 0.00185 at
0, 0.2 and 0.4, weak, centroid y): the data cannot say how correlated radar and gauge errors
are. The sweep shows what each assumption implies. At ρ_RG = 0 the model reproduces T7
(weak, centroid y: ρ_IM +0.361, radar 0.061, gauges 0.226, microwave 0.259).

## Tipping points

| case | ρ_IM at ρ_RG = 0 | ρ*₁ (satellite dependence no longer detected) | ρ*₂ (ranking flips) | admissible up to |
|---|---|---|---|---|
| weak, centroid y | +0.361 [0.294, 0.427] | 0.50 | 0.425 | ≥ 0.50 (grid end) |
| weak, centroid x | +0.368 [0.279, 0.447] | 0.50 | 0.475 | ≥ 0.50 |
| strong, centroid y | +0.431 [0.361, 0.498] | not reached (≤ 0.50) | 0.475 | ≥ 0.50 |
| strong, centroid x | +0.500 [0.442, 0.555] | not reached | 0.50 | ≥ 0.50 |

Path (weak, centroid y): ρ_RG 0.1 → ρ_IM +0.349, radar − MW −0.172; 0.2 → +0.333, −0.136;
0.3 → +0.307, −0.087; 0.4 → +0.263, −0.017; 0.5 → +0.172 [−0.087, 0.305], +0.090.

## Reading (per the plan's rule)

Tipping points lie between 0.10 and the end of the admissible range, so they are reported
as the dependence the conclusion tolerates: **the satellite dependence remains detected, and
radar remains ranked least noisy, for assumed radar–gauge error correlations up to about 0.4**
in every case. Only a radar–gauge error correlation of roughly 0.43–0.5 removes the inferred
satellite dependence or flips the ranking. Whether radar and gauge errors could be that
correlated is not something these data can decide; the grid stopped at 0.5, so the
admissible range is censored there.

## Correction after review (2026-10-15)

The reading above follows the plan's definitions: ρ*₂ is where the ranking's **point estimate** flips. By
interval, radar's lead (ν_R − ν_MW below zero with its 95% interval) holds only up to 0.275 (weak, centroid y),
0.375 (weak, x), 0.35 (strong, y) and 0.4 (strong, x); beyond those values the interval includes reversal.
The satellite dependence statement (interval excluding zero up to 0.475–0.5) is unaffected. Documents now
report both, by interval.
