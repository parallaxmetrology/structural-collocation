# T7-R: robustness check of T7 (frozen 2026-10-15, before any value was computed)

Robustness check, not a new probe. Same sample and code path as T7
(`data/t7_points.csv`: 2,402 window-hours with ≥ 5 gauges, point support; weak and
strong strata split at the median radar peak; centroid y and x; cluster bootstrap
over UTC day × 5° tile, 500 resamples, seed 11).

T7 replaced GOES with gauges. Two checks it left open:

## R1–R2: the mirror triplet (radar, GOES, gauges)

- **R1:** radar is the least noisy observer in (R, GIR, G): nf_R − nf_GIR < 0 and
  nf_R − nf_G < 0, 95% CIs below 0.
- **R2:** radar's noise fraction does not depend on which satellite is dropped:
  nf_R from (R, GIR, G) and from (R, MW, G) differ by < 0.05, or their difference CI
  includes 0.

## R3: four-observer consistency (tetrads)

With four observers and one common signal, independent errors imply three equal
products of correlations (vanishing tetrads). If the only breach is a shared
GOES–microwave error, it inflates r(GIR, MW) alone, which predicts:

- t3 = r(R,GIR)·r(G,MW) − r(R,MW)·r(G,GIR): **≈ 0** (CI includes 0);
- t1 = r(R,G)·r(GIR,MW) − r(R,GIR)·r(G,MW): **> 0** (CI above 0);
- t2 = r(R,G)·r(GIR,MW) − r(R,MW)·r(G,GIR): **> 0** (CI above 0).

Also reported: every observer's noise fraction from each of the three triplets that
contain it (a 4 × 3 table), to show which triplets agree.

## Reading

R1–R3 met: T7's conclusion is robust to the choice of satellite and the four-observer
structure matches a single shared GOES–microwave error. A failed R3 (t3 ≠ 0) would
mean a second breach exists (for example a radar–gauge or radar–satellite shared
error) and T7's ranking would need qualifying.
