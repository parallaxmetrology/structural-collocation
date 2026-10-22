# T5 results: shared smoothness vs shared physics (2026-10-15)

Run as frozen in [T5_PLAN.md](T5_PLAN.md). Code `scripts/t5.py`; exploratory
follow-up `scripts/t5_offsets.py`. Numbers `data/t5_results.json`,
`data/t5_offsets.json`.

## Smoothness match

Median lag-1 spatial autocorrelation inside sample windows: microwave 0.700; radar
0.633 unblurred, 0.759 at σ = 0.5 cell, 0.904 at σ = 1.0. **σ\* = 0.5 cell.** At the
0.1° grid, radar is already nearly as smooth as the gridded microwave field.

## Synthetic calibration arm (T1 fields; R blurred by + 0, + 1, + 2 cells)

| scenario | R blur | estimated nf_R − nf_S1, centroid y | truth |
|---|---|---|---|
| smoothness only (σ_sh = 0) | +0 | −0.011 [−0.019, −0.004] | −0.072 |
| | +1 | −0.052 [−0.059, −0.045] | −0.049 |
| | +2 | −0.057 [−0.064, −0.050] | −0.020 |
| shared displacement (σ_sh = 1.0) | +0 | **+0.139** [+0.129, +0.148] | −0.192 |
| | +1 | **+0.104** [+0.094, +0.113] | −0.169 |
| | +2 | **+0.098** [+0.090, +0.107] | −0.140 |

Smoothness alone never produces the inversion; a shared displacement produces it,
and blurring radar only shrinks it slowly. (Dispersion behaves less cleanly: with
smoothness only and no radar blur, its estimated difference has the wrong sign,
+0.063 vs truth −0.091. Centroids are the reliable variables.)

## Real data (R, GOES IR, MW), radar blur ladder

| radar blur σ (cells) | n | nf_R − nf_MW, centroid y | centroid x | dispersion |
|---|---|---|---|---|
| 0 | 14,751 | +0.224 [+0.196, +0.251] | +0.193 [+0.164, +0.223] | +0.288 [+0.222, +0.357] |
| **0.5 (σ\*)** | 14,819 | **+0.214 [+0.186, +0.240]** | **+0.181 [+0.151, +0.210]** | +0.304 [+0.240, +0.374] |
| 1.0 | 14,978 | +0.193 [+0.165, +0.223] | +0.161 [+0.134, +0.191] | +0.333 |
| 1.5 | 15,129 | +0.181 [+0.154, +0.212] | +0.144 [+0.116, +0.174] | +0.373 |
| 2.0 | 15,263 | +0.166 [+0.141, +0.192] | +0.128 [+0.100, +0.154] | +0.404 |
| 3.0 | 15,467 | +0.155 [+0.126, +0.186] | +0.111 [+0.084, +0.140] | +0.463 |

**Prediction outcome:** at σ\* the difference stays clearly above zero on both
centroids (and dispersion). **Shared smoothness is not a sufficient explanation**
of the real ranking. The ladder's shape (slow shrink, never crossing zero, even at
σ = 3, far smoother than microwave) matches the synthetic shared-displacement case.

## Exploratory follow-up: direction of the satellite offsets

Mean centroid offset from radar (cells of 0.1°; negative y = north):

| latitude band | n | GOES IR − radar (y) | MW − radar (y) |
|---|---|---|---|
| 20–30°N | 2,315 | +0.46 | +0.41 |
| 30–35°N | 3,978 | −0.04 | +0.08 |
| 35–40°N | 2,588 | −0.16 | −0.05 |
| 40–45°N | 1,987 | −0.47 | +0.01 |
| 45–55°N | 3,717 | **−0.90** | −0.12 |

- GOES cold cloud sits increasingly **north** of radar rain with latitude, ~0.9 cell
  (~10 km) at 45–55°N: the known geostationary parallax of high cloud tops (GOES-East
  at 75°W over the equator). The structural reads detect it without being told.
- Microwave shows little latitude trend, so **shared viewing-angle parallax is not
  the shared error**.
- In 20–30°N both satellites sit south of radar by a similar amount: an offset common
  to the two satellites relative to radar there (possibly radar coverage edges near
  the Gulf coast; not examined).
- Offset correlation corr(GIR − R, MW − R) = 0.47 (y), 0.45 (x); this includes the
  shared −R term by construction and is not diagnostic on its own.

## Where this leaves F3

The shared satellite error is real in the sense T1 defines (a pattern that only a
shared error reproduces), smoothness is excluded as sufficient, and parallax is
excluded as the main shared component. Ice aloft remains a candidate, now joined by
anything else the two satellite retrievals have in common relative to surface rain
(for example, both detecting the convective core rather than the full rain area).
