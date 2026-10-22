# T9 results: real-data window control and fixed-sample smoothness test (2026-10-15)

Run as frozen in [T9_PLAN.md](T9_PLAN.md). Code `scripts/t9.py`; numbers `data/t9_results.json`.
Statistic: nf_R − nf_MW (positive = radar rated noisier than microwave), cluster bootstrap
300 resamples, seed 21.

## T9-W: window grid

Triplet radar, GOES IR, microwave:

| variant | n | centroid y | centroid x |
|---|---|---|---|
| W0 baseline (15 × 15, offset 2) | 14,751 | +0.224 [0.198, 0.255] | +0.193 [0.165, 0.224] |
| W1 half-shifted grid (offset 9) | 14,691 | +0.199 [0.174, 0.227] | +0.175 [0.151, 0.199] |
| W2 large windows (25 × 25) | 8,435 | +0.181 [0.146, 0.211] | +0.191 [0.162, 0.224] |
| W3 centred storms only (descriptive) | 2,500 | +0.666 [0.533, 0.819] | +0.466 [0.326, 0.621] |

Point-support reversal (radar, gauges, microwave; ≥ 5 gauges):

| variant | n | centroid y | centroid x |
|---|---|---|---|
| W1 half-shifted grid | 3,061 | −0.231 [−0.269, −0.196] | −0.297 [−0.342, −0.246] |
| W2 large windows | 3,084 | −0.210 [−0.253, −0.174] | −0.248 [−0.293, −0.208] |

**Predictions met.** The gap stays positive with intervals above zero under every variant;
W1 and W2 are within ±50% of W0 (−11% and −19% on centroid y). The point-support reversal
holds under both. By the pre-registered rule the window grid is not a material contributor; this covers the tested
shift and size, not windowing or representativeness effects in general. In the descriptive
centred-storm subset the gap is larger, not smaller (the subset is selected using radar's
centroid, so it is not interpreted further).

## T9-S: fixed sample, matched smoothness (the plan's "exact match"; lag-1 within 0.013)

Fixed sample: 14,725 window-hours that pass the inclusion rules at every blur level.
Median lag-1 autocorrelation: microwave 0.700; radar 0.634 (σ = 0), 0.640 (0.3), 0.687 (0.4),
0.760 (0.5). Interpolated match σ\* = 0.42, rounded to **0.40** (radar 0.687 vs microwave 0.700).

| radar blur σ | centroid y | centroid x |
|---|---|---|
| 0 | +0.221 [0.196, 0.254] | +0.193 [0.165, 0.223] |
| **0.4 (σ\*)** | **+0.218 [0.192, 0.250]** | **+0.188 [0.160, 0.218]** |
| 1.0 | +0.196 [0.170, 0.225] | +0.160 [0.133, 0.191] |
| 2.0 | +0.163 [0.138, 0.195] | +0.124 [0.098, 0.155] |
| 3.0 | +0.145 [0.121, 0.176] | +0.104 [0.079, 0.135] |

**Prediction met.** With the sample held fixed and smoothness matched within 0.013 in lag-1
autocorrelation, the inversion is essentially unchanged (+0.218 vs +0.221 unblurred), and it
stays clearly positive even at σ = 3, far smoother than microwave. Blur levels below 0.3 cell
barely change a 0.1° field, so the first rungs coincide.

## Reading

Within the tested variants, neither the window grid nor the radar–microwave smoothness difference
accounts for the inversion. (Wording tightened after review, 2026-10-15; numbers unchanged.)
