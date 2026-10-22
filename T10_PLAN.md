# T10: replication on an unseen period (frozen 2026-10-15, before any August 2024 data was fetched)

From GATE.md: the next substantial test. Every result so far is in-sample (July 2024). T10 reruns
the frozen July pipeline, unchanged, on a period none of the analysis has touched, and scores it
against rules fixed here.

## Period and data

- **Period:** hourly accumulations ending 2024-08-01 01:00 UTC to 2024-08-22 00:00 UTC (504 hours;
  same length and season as July, immediately after it).
- **Observers, products, domain, grid:** identical to July (MRMS RadarOnly_QPE_01H + RQI; GOES-16
  ABI L2 CMIPC band 13, 4 scans/h, field max(0, 235 − Tb); IMERG Late V07 MWprecipitation and merged;
  USCRN + ASOS hourly gauges; CONUS, 0.1° MRMS-aligned grid).
- **Storage:** reduced arrays in `data_aug/` (same layout as `data/`); raw files streamed and
  discarded as in July. Data digests recorded in `T10_RESULTS.md`.

## Code

The analysis is the July code as of this freeze (sha256, first 16 hex):

| script | hash |
|---|---|
| probe_a.py | 65ad569f802fc929 |
| probe_a2.py | 259d8714d7ebf466 |
| t5.py | 40d1dbe4fcfcc001 |
| t7.py | 84b440bf959a51b7 |
| t7_robust.py | 4adfcdcf0daafdb3 |
| t7s.py | 0ebabbecc99c7262 |
| fetch_mrms.py | 8a843af2ba9958be |
| fetch_goes.py | 88e5cb3b9f4d652e |
| fetch_imerg.py | 56c923ebf8c3e5e6 |
| fetch_imerg_extra.py | 07e26986dd66a73a |
| fetch_gauges.py | a840ae23cc493863 |

**Permitted change:** a period switch (an environment variable selecting the August hours, dates and
the `data_aug/` root), with July remaining the default. After the switch is added, the July outputs
must reproduce byte-for-byte in their numbers (Probe A2, T7, T7-R, T7-S results files) before any
August analysis runs. Any other code change is a deviation, logged in `T10_RESULTS.md` with its reason.
No parameter (windows, inclusion rules, thresholds, strata, seeds, resample counts) may change.

## Analyses (all exactly as in July)

- **Area triplet** (Probe A2): (R, GIR, MW), 15 × 15 windows, grid offset (2, 2), inclusion rules as
  Probe A2, cluster bootstrap 500 resamples (UTC day × 5° tile), seed 11.
- **Point support** (T7): window-hours with ≥ 5 reporting gauges, all observers read at gauge cells,
  weak/strong strata split at the median radar peak of the August sample, 500 resamples, seed 7.
- **Mirror and tetrads** (T7-R): same sample, 500 resamples, seed 11.
- **Sensitivity** (T7-S): ρ_RG grid 0 to 0.5 step 0.025, 300 resamples, seed 13.

## Primary predictions (confirmatory)

Notation: Δ = nf_R − nf_MW. "Both axes" = centroid y and centroid x.

- **H1 (area inversion):** (R, GIR, MW): Δ > 0, 95% CI above 0, both axes.
- **H2 (point inversion, T7 P0):** weak stratum, (R, GIR, MW) at points: Δ > 0, CI above 0, both axes.
- **H3 (point reversal, T7 P1):** (R, MW, G): Δ < 0, CI below 0, both axes, in **both** strata.
- **H4 (mirror, T7-R R1):** weak stratum, (R, GIR, G): nf_R − nf_GIR < 0 and nf_R − nf_G < 0, CIs below 0, both axes.
- **H5 (satellite error correlation, T7 P2):** weak stratum, c(GIR, MW) with R and G as instruments: > 0, CI above 0, both axes.
- **H6 (tetrad pattern, T7-R R3):** weak stratum, centroid y: t1 > 0 and t2 > 0 (CIs above 0); t3 CI
  includes 0. As in July, t3 ≈ 0 is compatibility with a single two-group split, not proof that no
  other dependence exists, and t3 = t2 − t1.

## Secondary (descriptive, not scored)

- Point estimates and CIs next to July's for every quantity above.
- T7-S: the ρ_RG at which the ranking flips and at which the satellite dependence ceases to be
  detected. Expectation (not scored): both ≥ 0.3.
- H3 to H5 in the strong stratum where not already required; H6 on centroid x.

## Scoring

- **Inconclusive** (scored before any prediction is read): fewer than 1,000 point-support window-hours
  (July: 2,402), fewer than 5,000 area window-hours (July: 14,751), or a primary triplet failing the
  July stability rule (> 20% of resamples with a negative error variance).
- **Replicated:** H1 to H6 all met.
- **Not replicated:** H1 or H3 fails in sign (point estimate of the wrong sign) on either axis, in any
  stratum where it is required.
- **Partial:** anything else; every prediction is reported individually, met or not.

No outcome changes the July results; T10 is reported alongside them. A failed or partial replication
is a measurement of the period dependence of the effect, recorded as such.

## After

Write `T10_RESULTS.md`, add T10 to the lab book and the paper (Section 4 and the abstract's last
sentence), and stop for discussion.
