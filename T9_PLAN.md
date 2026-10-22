# T9: real-data window control and fixed-sample smoothness test (frozen 2026-10-15, before any T9 value was computed)

From GATE.md, the cheap next step. Two checks on data already on disk; then stop.

## T9-W: does the window grid drive the inversion?

T1-bis found that 15 × 15 windows add a small shared coordinate error in simulation
(centroid ≤ 0.02). T9-W checks this on real data.

Triplet (radar, GOES IR, microwave), unblurred, same inclusion rules as Probe A2. Variants:

| variant | window | grid offset (rows, cols) |
|---|---|---|
| W0 baseline | 15 × 15 | (2, 2), as Probe A2 |
| W1 half-shift | 15 × 15 | (9, 9) |
| W2 large | 25 × 25 | (2, 2) |
| W3 centred subset | 15 × 15 | (2, 2), only window-hours whose radar centroid lies within ±2 cells of the window centre (descriptive; selection uses radar) |

Plus the point-support reversal (T7 triplet radar, microwave, gauges, ≥ 5 gauges) under W1 and W2.

Statistic: nf_R − nf_MW (centroid y and x), cluster bootstrap 300 resamples (UTC day × 5° tile), seed 21.

**Predictions** (from T1-bis: windowing effects are small for centroids):
- the radar–GOES–microwave gap stays > 0 with 95% CI above 0 in W1, W2 and W3, on both axes, and its point estimate in W1 and W2 is within ±50% of W0;
- the point-support reversal (radar, microwave, gauges: nf_R − nf_MW < 0, CI below 0) holds under W1 and W2.

If the gap falls by more than half or loses its CI in W1 or W2, windowing is reported as a material contributor.

## T9-S: smoothness test on a fixed sample with an exact match

T5 used an approximate match (radar 0.759 vs microwave 0.700) and a sample that changed along
the blur ladder.

- **Fixed sample:** window-hours (Probe A2 rules, W0 grid) that pass the inclusion rules at
  every blur level of the ladder σ ∈ {0, 0.1, 0.2, …, 1.0, 1.5, 2.0, 3.0}.
- **Exact match:** median lag-1 autocorrelation of radar vs microwave on this fixed sample (every
  4th hour), evaluated on the fine part of the ladder; σ\* = the linear interpolation at which
  radar equals microwave, rounded to the nearest 0.05, then run as an extra level.
- Statistic as above, on the fixed sample at every level including σ\*.

**Prediction:** at σ\*, nf_R − nf_MW > 0 with 95% CI above 0 on centroid y and x.

## After

Write `T9_RESULTS.md`, update the lab book and preprint, and stop.
