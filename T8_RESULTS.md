# T8 results: anvil fraction and the shared satellite error (2026-10-15)

Run as frozen in [T8_PLAN.md](T8_PLAN.md). Code `scripts/t8.py`; exploratory
`scripts/t8_abs.py`. Numbers `data/t8_results.json`, `data/t8_abs.json`.
Sample: the T7 point-support set, 2,402 window-hours; anvil-fraction terciles with
means 0.12, 0.43, 0.75 (≈ 800 each). Anvil fraction correlates −0.34 with radar
window peak (anvil-heavy windows are lighter-rain windows).

## Predictions, scored (centroid y; x in brackets)

| | prediction | by tercile (low → high anvil) | top − bottom | verdict |
|---|---|---|---|---|
| P1 | GOES–MW error correlation grows with anvil | 0.375, 0.409, 0.389 (0.425, 0.434, 0.421) | +0.014 [−0.097, 0.111] (−0.004 [−0.125, 0.112]) | **not met** |
| P2 | (R, GIR, MW) inversion grows with anvil | 0.054, 0.178, 0.362 (0.099, 0.201, 0.448) | +0.308 [0.194, 0.455] (+0.350 [0.164, 0.529]) | **met** |
| P3 | radar least noisy in (R, MW, G) in every tercile | radar 0.028 / 0.049 / 0.081 vs MW 0.217 / 0.217 / 0.325, gauges 0.180 / 0.207 / 0.248 | | **holds** |

**Outcome: mixed.** The normalised coupling between the two satellites' errors is
constant at about 0.4 regardless of anvil extent; what grows with anvil fraction is
how badly the radar–satellite–satellite triplet is fooled.

## Exploratory: absolute error sizes by anvil tercile (cells², centroid y)

| tercile | shared GOES–MW error covariance | GOES error var. | MW error var. | radar error var. | gauge error var. |
|---|---|---|---|---|---|
| low anvil | 1.31 | 5.32 | 2.29 | 0.30 | 2.31 |
| middle | 1.63 | 6.15 | 2.57 | 0.59 | 3.24 |
| high anvil | 2.64 | 9.75 | 4.72 | 1.23 | 4.70 |

The shared covariance doubles from low to high anvil, but so does each satellite's
error variance, and the gauges' too. Anvil-heavy windows are noisier for every
observer (they are lighter-rain windows), so anvil fraction is **confounded** with
general light-rain noisiness and is not identified as a specific driver.

## Where this leaves the cause

- The two satellites share an error of roughly constant strength (normalised ≈ 0.4)
  relative to surface rain, present at all anvil levels.
- Its absolute size scales with conditions that make every observer noisier.
- H_aloft in its tested form (coupling grows with aloft-only cloud) is not
  supported; a constant coupling is compatible with both satellites observing the
  cloud system (cloud-top and ice structure) rather than surface rain everywhere,
  not only in anvils. Not established.
