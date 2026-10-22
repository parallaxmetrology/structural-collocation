# Reviewer guide: structural collocation of radar, satellite and gauge rain fields

Prepared 2026-10-15 for outside technical review. Self-contained entry point; the full record is
`LABBOOK.html` (chronological, including dead ends) and the manuscript draft is `PAPER.html`.

## The claim, in its bounded form

1. **Empirical (replicated).** Triple collocation on the intensity-weighted rain centroid in 1.5° windows
   gives observer rankings that depend on which observers form the triplet. With MRMS radar-only, GOES-16
   10.3 µm and IMERG passive microwave, microwave is rated least noisy. At gauge locations, replacing either
   satellite with rain gauges reverses the ranking (radar least noisy). This held in July 2024 and, under a
   frozen pipeline, in August 2024.
2. **Interpretation (conditional).** Four-observer checks (two independent vanishing-tetrad constraints)
   show a pattern consistent with a single split, {radar, gauges} vs {GOES, microwave}. Read as an error
   shared by the two satellites relative to surface rain, its correlation is about +0.35 to +0.51 (July),
   +0.39 to +0.47 (August). Four observers cannot distinguish this from correlated radar and gauge errors.
   Sensitivity, by 95% interval: the satellite dependence stays detected for an assumed radar–gauge error
   correlation up to 0.475–0.5 (July) and through 0.5, the top of the tested range (August); radar's lead
   in the ranking is clear up to 0.275–0.4 (July) and 0.375–0.4 (August).
3. **Method (simulation).** On a synthetic control with known truth and a mass-conserving blur, structural
   collocation recovers centroid and dispersion noise fractions with correct rankings; a planted shared
   satellite displacement reproduces the inversion (at σ_sh = 1 cell); with five observers, a tetrad audit
   identifies a planted shared pair in either direction (10/10 blocks, no false alarms), under the
   simulation's assumptions.

**Not claimed:** the physical cause of the dependence; that radar–gauge errors are in fact weakly
correlated; that the replication bears on attribution (same instruments, processing and season, so a
persistent measurement or representativeness effect would replicate too); anything outside warm-season
CONUS, hourly, 0.1°. A pre-registered negative result (structural disagreement does not predict hourly
radar error at gauges) is reported separately (paper Appendix B).

## Claim ledger

| claim | support | boundary |
|---|---|---|
| Observer choice changes the estimated ranking: with gauges in place of either satellite, radar ranks least noisy. | Same point support; mirror triplet; frozen August replication. Robust to the tested window shifts and sizes, smoothing match and core-only reading. | Estimated rankings, not independently established accuracy rankings. |
| One common signal with mutually independent errors does not describe the four observers. | Tetrads t1, t2 > 0 (intervals excluding 0) in July and August; pattern compatible with the {radar, gauges} / {GOES, microwave} split. | A model failure; its physical cause is not identified. |
| Centroid and dispersion support collocation in the tested simulation. | Mass-conserving control recovers true rankings; planted shared displacement reproduces the inversion. | This simulation only, not every field or observing configuration. |
| Structural disagreement did not give the pre-specified improvement in predicting hourly radar error at gauges. | Validation Δρ +0.0014 against a +0.02 requirement; the spatial-scramble null gave +0.0031. | The tested hourly point-error task. |

**Conditional interpretation.** If radar and gauge errors carry sufficiently independent information about
surface rain, the results imply a shared GOES–microwave structural error (correlation about +0.35 to +0.51 in
July, assuming independence). The sensitivity analysis bounds the tolerated radar–gauge dependence; it does
not measure it.

**Not presently claimed:** a satellite defect; the physical mechanism; an advantage over existing
amount-based or representativeness analyses (not compared); a diagnostic that improves on available ones
(not compared); generality across seasons and regions.

**Open questions.** (a) Does geometry reveal dependence that amount-based collocation misses? (b) Do the two
groups respond to different physical quantities (surface rain vs aspects of the cloud system)? Processing,
resolution and sampling are alternative explanations. (c) Can the procedure serve as a warning that a
collocation ranking is not trustworthy under its assumptions, compared with existing diagnostics?

## Evidence chain

Every experiment had its plan frozen before any of its values were computed. Exploratory splits are labelled.

| step | question | plan | results | script | outcome |
|---|---|---|---|---|---|
| Probe A2 | ranking from (radar, GOES IR, microwave) | `PROBE_A2_PLAN.md` | `PROBE_A2_RESULTS.md` | `scripts/probe_a2.py` | microwave rated least noisy |
| T5 / T9 | does smoothness explain it? | `T5_PLAN.md`, `T9_PLAN.md` | `T5_RESULTS.md`, `T9_RESULTS.md` | `t5.py`, `t9.py` | not within the tested variants (close match on a fixed sample) |
| T9 | does the window grid explain it? | `T9_PLAN.md` | `T9_RESULTS.md` | `t9.py` | robust to tested shift and size |
| T6 | convective-core keying, beam quality | `T6_PLAN.md` | `T6_RESULTS.md` | `t6*.py` | not supported (RQI exploratory) |
| T7 | gauges as a fourth observer | `T7_PLAN.md` | `T7_RESULTS.md` | `t7.py` | reversal; GOES–MW error correlation |
| T7-R | mirror triplet, tetrads | `T7R_PLAN.md` | in `T7_RESULTS.md` | `t7_robust.py` | consistent with one two-group split |
| T8 | anvil fraction | `T8_PLAN.md` | `T8_RESULTS.md` | `t8.py` | no detectable change in coupling |
| T1-bis | corrected synthetic control | `T1BIS_PLAN.md` | `T1BIS_RESULTS.md` | `t1bis.py` | recovery, planted inversion, 5-observer identification |
| T7-S | radar–gauge sensitivity | `T7S_PLAN.md` | `T7S_RESULTS.md` (see its correction note) | `t7s.py` | bounds above |
| T10 | unseen-period replication | `T10_PLAN.md` | `T10_RESULTS.md` | `t10_area.py`, `t10_score.py` | replicated, 6/6 predictions |

Superseded: the first synthetic control (`T1_*`, non-mass-conserving blur) and the finding built on it (F6,
withdrawn). Both are kept in the lab book for the record. `GATE.md` holds the decision gate and its updates.

## Reproduction

Python 3 with numpy, scipy, pandas, h5py, eccodes (environment versions in `LABBOOK.html` §17). IMERG needs a
NASA Earthdata token in `secrets/earthdata_token.txt`; MRMS, GOES and gauges are open.

```
# July (from the project root): fetch, then analyse
.venv/bin/python scripts/fetch_gauges.py && .venv/bin/python scripts/make_stations.py
.venv/bin/python scripts/fetch_mrms.py && .venv/bin/python scripts/fetch_goes.py && .venv/bin/python scripts/fetch_imerg_both.py
.venv/bin/python scripts/probe_a2.py && .venv/bin/python scripts/t7.py && .venv/bin/python scripts/t7_robust.py && .venv/bin/python scripts/t7s.py

# August (T10): same scripts from runs/aug, where data -> ../../data_aug
cd runs/aug && export PRECIP_PERIOD=aug   # then the same fetch and analysis commands with ../../ paths, plus t10_area.py and t10_score.py
```

Reduced inputs are about 1.7 GB per month (digests in `LABBOOK.html` §17 and `T10_RESULTS.md`). Upstream
archives can be reprocessed, so fresh downloads may not reproduce the digests exactly.

## Checks on the computation

- **Independent re-implementation** (`audit/AUDIT.md`): written from the plan files alone, without reading
  the analysis code; reproduced 101 July quantities to 4 × 10⁻⁹ with identical samples.
- **Details fixed by code, not by plans** (found by that audit): ≥ (not >) at the 0.5 mm h⁻¹ wet threshold;
  float32 window sums for the ≥ 1 mm inclusion test; T7 does not require the merged IMERG product; the
  weak stratum is ≤ the median radar peak.
- **Period switch:** after adding it, all July result files reproduced byte-for-byte.
- **Deviations** from each plan are listed in its results file and in `LABBOOK.html` §17.

## Questions for a precipitation specialist

1. Does the centroid analysis reveal dependence information beyond established collocation and
   spatial-verification methods?
2. Is a radar–gauge error correlation below about 0.3–0.4, for hourly centroids at gauge points, physically
   defensible in this configuration (MRMS radar-only, not gauge-adjusted; ASOS and USCRN gauges)?
3. Is the contribution sufficiently distinct from existing work on representativeness error in
   collocation (e.g. Vogelzang and Stoffelen 2021)?
4. Are point-support centroids from about 8–9 gauges per 1.5° window a defensible observable, or could their
   coarseness itself produce the radar–gauge agreement?

## Contents of the frozen snapshot

`snapshots/precip-disagreement_2026-10-15.tar.gz` with `MANIFEST.sha256`: plans, results, scripts, builders,
`PAPER.html`, `LABBOOK.html`, figures, result JSONs and point tables for both months, and the audit. Raw and
reduced input arrays are not included (recreated by the fetch scripts; digests recorded).
