# T6 results: core vs full area (2026-10-15)

Run as frozen in [T6_PLAN.md](T6_PLAN.md). Code `scripts/t6.py`; exploratory
follow-ups `scripts/t6_intensity.py`, `scripts/t6_rqi.py`. Numbers
`data/t6_results.json`, `data/t6_intensity.json`, `data/t6_rqi.json`.
Sample: 14,751 window-hours (Probe A2 triplet, unblurred).

## Predictions, scored

| | prediction (under H_core) | result |
|---|---|---|
| P1 | core reading shrinks nf_R − nf_MW below the full-area point (CI upper end < +0.224 y / +0.193 x) | **not met**: core +0.207 [0.179, 0.235] (y), +0.180 [0.150, 0.213] (x) |
| P2 | inversion grows with light-rain fraction (top − bottom tercile CI > 0) | **met in form**: +0.146 [0.089, 0.211] (y), +0.122 [0.062, 0.185] (x). Caveat: 2/3 of windows are entirely light rain (fraction = 1.0), so the split is effectively "has a ≥ 5 mm/h core" (bottom tercile, 0.126) vs "all light" (0.268, 0.273) |
| P3 | satellites closer to radar's core centroid than to its full-area centroid | **reversed**: median paired difference −0.96 [−0.99, −0.92] cells (GOES IR), −1.10 [−1.12, −1.08] (MW); only 26% / 16% of windows closer to the core. Caveat: a core centroid (smaller area) is noisier, which biases this comparison toward the full area |

**Outcome: mixed; core keying is not supported as the shared error.** The one met
prediction points the other way from the hypothesis' picture: the inversion is
largest where there is **no** heavy core.

## Exploratory follow-up 1: inversion by radar peak intensity (quintiles)

| radar window peak (mm/h) | n | nf_R − nf_MW, centroid y | nf_R | nf_MW |
|---|---|---|---|---|
| ≤ 2.1 | 2,951 | +0.333 [0.257, 0.406] | 0.394 | 0.062 |
| ≤ 5.2 | 2,950 | +0.224 [0.170, 0.283] | 0.273 | 0.048 |
| ≤ 10.0 | 2,950 | +0.206 [0.149, 0.264] | 0.238 | 0.031 |
| ≤ 17.7 | 2,950 | +0.180 [0.143, 0.222] | 0.187 | 0.008 |
| ≤ 95.4 | 2,950 | +0.115 [0.086, 0.147] | 0.134 | 0.019 |

Monotone: the inversion triples from the strongest to the weakest windows, driven by
radar's estimated noise rising in light rain.

## Exploratory follow-up 2: is it radar's genuine beam-quality error?

Alternative explanation: in light, shallow rain far from a radar the beam
overshoots, so radar is genuinely the worse observer there and collocation is right.
MRMS's radar quality index (RQI) measures exactly that. Inversion by RQI tercile ×
intensity (below / above median peak):

| window RQI (mean) | weak rain | strong rain |
|---|---|---|
| low (0.25–0.34) | +0.262 [0.206, 0.329] | +0.165 [0.109, 0.222] |
| middle (0.82) | +0.259 [0.181, 0.335] | +0.143 [0.104, 0.180] |
| high (0.96) | **+0.317 [0.225, 0.427]** | +0.169 [0.134, 0.204] |

The inversion does not fall with better beam quality; it is largest in the
best-beam, weak-rain windows. Intensity, not RQI, carries it (corr RQI–peak 0.24).
Radar beam overshoot, as RQI measures it, is not the explanation. In the
high-RQI strong-rain cell MW's noise fraction is −0.002, an impossible value that
flags an assumption breach there.

## Where this leaves F3

The shared satellite error concentrates in **light rain**, independent of radar
beam quality, and is not explained by core keying, smoothness or parallax.
Candidates still open: (a) the two satellites share a light-rain blind spot or
bias (light rain is often shallow, warm or ice-poor, so cold tops and ice
scattering are weak or absent in the same places); (b) radar has a light-rain error
that RQI does not encode (e.g. drizzle drop-size effects, evaporation below beam).
The two predict different things for whether the satellites agree with **each
other** on light-rain structure (T7).
