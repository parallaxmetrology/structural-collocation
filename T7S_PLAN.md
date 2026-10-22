# T7-S: sensitivity to radar–gauge error dependence (frozen 2026-10-15, before any T7-S value was computed)

Why: T7's satellite attribution (GOES–microwave error correlation +0.35 to +0.51) and the
reading of collocation scores as true error rankings assume radar and gauge errors are
independent. Four observers cannot establish that. T7-S asks how much radar–gauge error
dependence would erase the inferred satellite dependence or change the ranking.

## Data

`data/t7_points.csv` (2,402 window-hours, point support), weak and strong strata split at the
median radar peak, centroid y and centroid x.

## Model

Four observers R, G, I (GOES), M (microwave); one common signal T with unit variance:

    Σ = λ λᵀ + Ψ,   Ψ diagonal ψ_i, plus ψ_RG = ρ_RG·√(ψ_R ψ_G) (fixed by assumption)
                    and ψ_IM free.

Parameters: 4 loadings, 4 error variances, ψ_IM (9) against 10 moments, so the model is
identified for each fixed ρ_RG. With ρ_RG free it is not (the four-observer ambiguity);
the sweep traces that equivalence class. Fitted by maximum likelihood (Gaussian
discrepancy) with ψ_i > 0; a solution is **admissible** if all ψ_i > 0, |ρ_IM| ≤ 1 and Ψ and Σ
are positive definite. At ρ_RG = 0 it should reproduce T7's estimates closely.

Sweep: ρ_RG ∈ {0, 0.025, 0.05, …, 0.50}. At each value report ρ_IM = ψ_IM/√(ψ_I ψ_M), the noise
fractions ν_i = ψ_i/(λ_i² + ψ_i), admissibility, and fit discrepancy; 300 cluster-bootstrap
resamples (UTC day × 5° tile), seed 13, for ρ_IM and ν_R − ν_M.

## Tipping points (the deliverable)

- **ρ*₁:** smallest ρ_RG at which the 95% CI of ρ_IM includes 0 (satellite dependence no
  longer detected).
- **ρ*₂:** smallest ρ_RG at which ν_R ≥ ν_M in the point estimate (ranking flips).
- **ρ_max:** largest ρ_RG with an admissible solution.

## Reading (fixed now)

- A tipping point at ρ_RG ≤ 0.10: a material limitation, stated in the abstract.
- A tipping point between 0.10 and ρ_max: reported as the dependence the conclusion tolerates.
- No tipping point within the admissible range: the conclusion survives every feasible
  radar–gauge dependence in this model.

## Decision gate

After T1-bis and T7-S, stop and write `GATE.md`: what the corrected controls support, and
which next paper they point to (satellite dependence, measurement-window effects, or the
limits of structural collocation). No further experiments until that is discussed.
