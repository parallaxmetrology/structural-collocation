# T7 results: whose error is it? Gauges as a fourth observer (2026-10-15)

Run as frozen in [T7_PLAN.md](T7_PLAN.md). Code `scripts/t7.py`; numbers
`data/t7_results.json`; point table `data/t7_points.csv`.

## Sample

2,402 window-hours with ≥ 5 reporting gauges and rain at the points for every
observer (median 9 gauges per window), split at the sample's median radar peak
(15.8 mm/h) into 1,201 "weak" and 1,201 "strong". Note: requiring rain at the gauge
points skews this sample wetter than T6's; "weak" here is not T6's lightest quintile.

## Results (point support, centroid y; centroid x in the table below)

| | weak | strong |
|---|---|---|
| **P0** (R, GIR, MW): nf_R − nf_MW | **+0.195 [0.116, 0.281]** (nf R 0.258, GIR 0.512, MW 0.063) | +0.185 [0.119, 0.255] |
| **P1** (R, MW, G): nf_R − nf_MW | **−0.196 [−0.260, −0.134]** (nf R 0.063, MW 0.259, G 0.225) | −0.216 [−0.270, −0.167] |
| **P2** error corr. GIR–MW (instruments R, G) | **+0.350 [0.287, 0.418]** | +0.450 [0.386, 0.516] |

| centroid x | weak | strong |
|---|---|---|
| P0 nf_R − nf_MW | +0.203 [0.104, 0.312] | +0.302 [0.227, 0.393] |
| P1 nf_R − nf_MW | −0.259 [−0.335, −0.202] | −0.365 [−0.444, −0.282] |
| P2 GIR–MW error corr. | +0.361 [0.284, 0.437] | +0.508 [0.458, 0.557] |

Negative-error-variance shares: 0.00–0.01 in all weak-stratum triplets (0.07 in one
strong-stratum case), so every result is scored.

Raw correlations at the points (weak, centroid y): radar–gauge 0.85, MW–gauge 0.76,
GIR–gauge 0.53; radar–MW 0.83, radar–GIR 0.60, GIR–MW 0.68.

## Predictions, scored

- **P0 met:** point support reproduces the inversion.
- **P1 met for H_sat, decisively:** replacing the second satellite with gauges
  reverses the ranking. Radar becomes the least noisy observer (0.06), gauges next
  (0.22, consistent with point-vs-pixel representativeness), microwave noisiest
  (0.26). H_radar predicted the inversion would stay; it did not.
- **P2 met for H_sat:** with radar and gauges as independent instruments, the GOES
  and microwave errors correlate at +0.35 (weak) and +0.45 (strong), intervals well
  above zero.
- **P3 (strong stratum):** the same picture, slightly stronger.

## Reading

The inversion seen in Probes A, A2, T5 and T6 was an artefact of collocation with
two observers that share an error. With a sharp, independent surface observer as the
third leg, collocation ranks the observers in the physically expected order (radar
best). This is the T1 lesson confirmed on real data: a shared error is invisible to
an estimator whose instruments share it, and visible to one with a sharp independent
instrument.

**Not established:** why the two satellites share an error relative to surface rain.
Candidates that survive: shared physics (both respond to ice and cloud aloft rather
than rain at the ground) or retrieval coupling not yet identified. The light-rain
concentration seen in T6 is not tested cleanly here (this sample is wetter).

Caveat: c(R, MW) (instruments GIR and G) is unstable in the weak stratum (CI
−6.0 to +0.1) and +0.35 [0.05, 0.55] in the strong stratum; it uses GIR as an
instrument, which T7 itself shows shares error with MW, so it is not interpreted.

## Robustness check (T7-R, pre-registered in T7R_PLAN.md; `scripts/t7_robust.py`, `data/t7r_results.json`)

| check | weak, centroid y | weak, x | strong, y | strong, x | prediction |
|---|---|---|---|---|---|
| R1 nf_R − nf_GIR in (R, GIR, G) | −0.592 [−0.694, −0.494] | −0.624 | −0.409 | −0.560 | < 0: **met** |
| R1 nf_R − nf_G in (R, GIR, G) | −0.215 [−0.313, −0.125] | −0.255 | −0.097 [−0.177, −0.026] | −0.247 | < 0: **met** |
| R2 nf_R(R,GIR,G) − nf_R(R,MW,G) | −0.029 [−0.072, +0.006] | −0.027 [−0.086, +0.024] | +0.036 [−0.002, +0.070] | +0.034 [−0.024, +0.088] | ≈ 0: **met** |
| R3 t3 (excludes the GOES–MW pair) | +0.014 [−0.003, +0.033] | +0.012 [−0.010, +0.033] | −0.021 [−0.041, +0.001] | −0.014 [−0.037, +0.009] | ≈ 0: **met** |
| R3 t1 | +0.120 [+0.089, +0.153] | +0.137 | +0.140 | +0.203 | > 0: **met** |
| R3 t2 | +0.134 [+0.099, +0.173] | +0.148 | +0.120 | +0.189 | > 0: **met** |

Every observer's noise fraction from each triplet that contains it (weak, centroid y):
triplets holding both satellites give microwave 0.03–0.06 and gauges 0.41; triplets
holding one satellite give microwave 0.26 and gauges 0.22–0.25, radar 0.03–0.06.
The four-observer structure is a single two-group split, {radar, gauges} vs
{GOES, MW}, with no other breach. T7's conclusion is robust to which satellite is dropped.

**Identification caveat** (from `tetrad_audit`, added afterwards): four observers cannot
distinguish a shared GOES–MW error relative to surface rain from a shared radar–gauge
error relative to a cloud-system signal. The satellite reading assumes surface rain is
the target; a fifth observer of a different physical mechanism would separate them.
