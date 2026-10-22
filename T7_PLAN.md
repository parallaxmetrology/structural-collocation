# T7: whose error is it in light rain? (frozen 2026-10-15, before any T7 value was computed)

From T6: the inversion (collocation rates microwave below radar) concentrates in
light rain and does not depend on radar beam quality. Two explanations remain:

- **H_sat:** the two satellites share a light-rain error (a common blind spot or
  bias relative to surface rain).
- **H_radar:** radar has a light-rain error that RQI does not encode, and the
  collocation ranking is correct there.

Satellite-versus-satellite agreement alone cannot separate these: with three
observers, the inversion is algebraically the statement r(GIR, MW) > r(R, GIR). A
**fourth, independent, sharp observer** can (T1: a sharp independent instrument
recovers a shared error that blurred instruments cannot). That observer is the
gauge network: point measurements of surface rain, used by neither MRMS radar-only,
IMERG Late nor GOES.

## Sample and point support

- Probe A2 windows (radar, GOES IR, microwave; unblurred), restricted to
  window-hours with **≥ 5 gauges reporting** (ASOS + USCRN, same hour assignment as
  Experiment 1).
- Every observer is read **at the same points**: radar, GOES cold-cloud and
  microwave values in each gauge's 0.1° cell, and the gauge values themselves.
- Point centroid (y, x): weighted mean of the gauges' window-local positions, with
  each observer's values at those points as weights. A window-hour enters when every
  observer's point values sum > 0.
- **Strata:** radar window peak (block mean) below vs above the median of this
  sample ("weak", "strong").
- Cluster bootstrap 500 resamples (UTC day × 5° tile), seed 7.

## Predictions

- **P0 (point support reproduces the effect):** in the weak stratum, triple
  collocation on (R, GIR, MW) at the points gives nf_R − nf_MW > 0 (CI above 0) on
  centroid y and x. If not, point support changes the picture and P1–P2 are reported
  as uninterpretable.
- **P1 (replace one satellite with gauges):** weak stratum, triplet (R, MW, G).
  H_sat predicts the inversion disappears or shrinks (nf_R − nf_MW CI includes or is
  below 0, or its point estimate < half of P0's). H_radar predicts it stays > 0 and
  similar to P0.
- **P2 (quadruple check):** weak stratum, error correlation of GIR and MW with R and
  G as instruments: c(GIR, MW) = C(GIR, MW) − C(GIR, R)·C(MW, G)/C(R, G), normalised.
  H_sat predicts > 0 (CI above 0); H_radar predicts a CI including 0.
- **P3 (descriptive):** the same three quantities in the strong stratum.

## Stability rule

If more than 20% of bootstrap resamples give a negative error variance in a
weak-stratum triplet, that triplet's result is reported as unstable and not scored.
