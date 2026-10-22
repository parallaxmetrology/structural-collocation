# T8: is the shared satellite error driven by cloud and ice aloft beyond the surface rain? (frozen 2026-10-15, before any T8 value was computed)

From T7: with gauges as a fourth observer, GOES and microwave errors relative to
surface rain correlate at +0.35 to +0.51. Leading candidate cause, **H_aloft**: both
satellites respond to cold cloud and ice aloft (GOES: cold cloud tops; microwave over
land: ice scattering), which extends beyond where rain reaches the ground (anvils,
stratiform shields). Where that "aloft-only" area is large, both satellites should
be pulled the same way, away from surface rain.

## Measure (fixed now)

**Anvil fraction** of a window-hour = share of GOES cold-cloud cells (235 − Tb > 0)
where radar shows no surface rain (< 0.5 mm/h). Range 0–1. Computed from the full
window fields.

## Sample and method

The T7 point-support sample (`data/t7_points.csv`, 2,402 window-hours with ≥ 5
gauges). Split into anvil-fraction terciles. Within each tercile, as in T7: GOES–MW
error correlation with radar and gauges as instruments, and the (R, GIR, MW)
inversion nf_R − nf_MW. Joint cluster bootstrap (UTC day × 5° tile), 500 resamples,
seed 9, for each tercile and for the top − bottom difference. Centroid y and x.

## Predictions under H_aloft

- **P1:** GOES–MW error correlation is larger in the top anvil tercile than the bottom:
  top − bottom > 0 with 95% CI above 0, on centroid y and x.
- **P2:** the (R, GIR, MW) inversion is larger in the top tercile: top − bottom > 0,
  CI above 0.
- **P3 (descriptive):** in every tercile, the (R, MW, G) triplet keeps radar least
  noisy (the T7 ranking holds where the satellites are least disturbed too).

Also reported: correlation of anvil fraction with radar window peak (T6 located the
shared error in light rain; anvils and stratiform shields are typically light-rain
regions).

## Interpretation

P1 and P2 met: cloud and ice aloft beyond surface rain is supported as a driver of
the shared error. Neither met: anvil fraction does not modulate it and the cause
stays open. One met: reported as mixed.
