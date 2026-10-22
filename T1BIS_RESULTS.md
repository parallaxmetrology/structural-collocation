# T1-bis results: synthetic assumption audit (2026-10-15)

Run as frozen in [T1BIS_PLAN.md](T1BIS_PLAN.md). Code `scripts/t1bis.py`; numbers
`data/t1bis_results.json`. 20,000 windows; 19,825–19,998 valid per region.

## S-A analytical reference: PASS

Mass-conserving convolution on a 64 × 64 canvas: median centroid shift 5 × 10⁻¹⁶ cell
(max 4 × 10⁻¹⁵); squared dispersion grew by 2σ_b² within 0.007% (σ_b = 1) and 0.035%
(σ_b = 2). The rendering behaves as theory says.

## S-B / S-C: are the coordinate errors independent? (true residual correlations, S1–S2)

| region | centroid y | centroid x | d | d² |
|---|---|---|---|---|
| full canvas | −0.001 | +0.009 | +0.021 | +0.021 |
| central 15 × 15 window | +0.009 | +0.018 | +0.057 | +0.057 |
| window shifted (+3, +3) | +0.019 | +0.015 | +0.123 | +0.101 |
| window shifted (−3, +4) | +0.023 | +0.020 | +0.140 | +0.119 |

**S-B prediction met** (|r| < 0.05 on the full canvas for cy, cx and d²). With a blur that
conserves mass, independent per-cell errors give independent coordinate errors.
Windowing adds a small shared error for centroids (up to 0.023) and a moderate one for
dispersion (up to 0.14), growing as the storm sits further from the window centre.

## S-D: collocation accuracy and ranking (R, S2, S1)

Maximum |estimate − truth| and whether the ranking is right:

| region | cy | cx | d | d² |
|---|---|---|---|---|
| full | 0.002 ✓ | 0.002 ✓ | 0.005 ✓ | 0.007 ✓ |
| central window | 0.003 ✓ | 0.003 ✓ | 0.011 ✓ | 0.012 ✓ |
| shifted (+3, +3) | 0.006 ✓ | 0.004 ✓ | 0.024 ✓ | 0.024 ✓ |
| shifted (−3, +4) | 0.005 ✓ | 0.004 ✓ | 0.031 ✓ | 0.031 ✓ |

**Prediction met**, and more: every variable, including plain dispersion d, is recovered
with the correct ranking in every region. T1's dispersion mis-ranking came from the
non-conserving blur. d and d² perform alike here; d² has no advantage in this setting.

## S-E: planted shared displacement (S1 and S2), centroid y

| σ_sh | region | estimate R / S2 / S1 | truth R / S2 / S1 | inverted |
|---|---|---|---|---|
| 0.5 | full | 0.090 / 0.250 / 0.095 | 0.046 / 0.284 / 0.135 | no |
| 0.5 | window | 0.084 / 0.242 / 0.090 | 0.044 / 0.273 / 0.127 | no |
| 1.0 | full | 0.199 / 0.227 / 0.084 | 0.046 / 0.350 / 0.230 | **yes** |
| 1.0 | window | 0.181 / 0.223 / 0.083 | 0.044 / 0.334 / 0.214 | **yes** |

**Prediction met** at σ_sh = 1.0 (as specified). At 0.5 the estimates move toward inversion
(R overrated, S1 underrated) but do not cross; T1's earlier flip at 0.5 was helped by the
blur artefact.

## S-F: known-dependence benchmark (five observers, 10 blocks of 2,000 windows, cy)

| planted | hit rate | false alarms |
|---|---|---|
| nothing | — | 0 / 10 blocks |
| S1–S2 shared displacement | 10 / 10 | 0 / 10 |
| R–Gs shared displacement (the mirror case) | 10 / 10 | 0 / 10 |

With five observers the tetrad audit identifies the planted pair every time, in both
directions, and raises no false alarm when nothing is planted.

## Reading (per the plan's rules)

- **S-B ≈ 0 and S-C small:** the T1 baseline shared error (0.16) came from the blur not
  conserving mass. **F6 ("smoothing creates shared structural error") is withdrawn.**
  Windowing contributes a smaller, measurable shared error, strongest for dispersion and
  off-centre storms.
- **Dispersion reinstated:** with a correct blur, d is recovered and ranked correctly.
- **S-E inverts at σ_sh = 1.0:** the synthetic support for "a shared satellite error
  produces the observed inversion" is restored under the corrected control, with the
  threshold for inversion higher than T1 suggested.
- **S-F:** five observers separate a shared satellite error from the mirror radar–gauge case
  in simulation. This is the design the real data needs (a fifth observer whose own
  independence can be argued).
