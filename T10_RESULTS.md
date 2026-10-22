# T10 results: replication on an unseen period, August 2024 (2026-10-15)

Run as frozen in [T10_PLAN.md](T10_PLAN.md). August data were fetched only after the plan was frozen.
Numbers: `data_aug/t10_score.json` (August), `data/t10_score_july.json` (July, same scorer).
Δ = nf_R − nf_MW; 95% cluster-bootstrap intervals.

## Verdict: Replicated

Gates first: 1,933 point-support window-hours (≥ 1,000), 14,738 area window-hours (≥ 5,000), largest
negative-variance share 0.10 in a primary triplet (≤ 0.20). Not inconclusive.

All six primary predictions are met in August on both axes and in every required stratum.

| | quantity | July 2024 (in-sample) | August 2024 (T10) | met |
|---|---|---|---|---|
| H1 | area Δ cy | +0.221 [+0.194, +0.250] | +0.259 [+0.228, +0.289] | yes |
| H1 | area Δ cx | +0.193 [+0.164, +0.224] | +0.189 [+0.161, +0.214] | yes |
| H2 | point (R,GIR,MW) Δ weak cy | +0.195 [+0.116, +0.281] | +0.244 [+0.124, +0.366] | yes |
| H2 | point (R,GIR,MW) Δ weak cx | +0.203 [+0.104, +0.312] | +0.244 [+0.140, +0.362] | yes |
| H3 | point (R,MW,G) Δ weak cy | -0.196 [-0.260, -0.134] | -0.286 [-0.374, -0.181] | yes |
| H3 | point (R,MW,G) Δ weak cx | -0.259 [-0.335, -0.202] | -0.329 [-0.419, -0.249] | yes |
| H3 | point (R,MW,G) Δ strong cy | -0.216 [-0.270, -0.167] | -0.210 [-0.286, -0.138] | yes |
| H3 | point (R,MW,G) Δ strong cx | -0.365 [-0.444, -0.282] | -0.279 [-0.355, -0.215] | yes |
| H4 | mirror R1_R_minus_IR weak cy | -0.592 [-0.694, -0.494] | -0.572 [-0.675, -0.468] | yes |
| H4 | mirror R1_R_minus_G weak cy | -0.215 [-0.313, -0.125] | -0.174 [-0.276, -0.061] | yes |
| H4 | mirror R1_R_minus_IR weak cx | -0.624 [-0.730, -0.520] | -0.659 [-0.780, -0.541] | yes |
| H4 | mirror R1_R_minus_G weak cx | -0.255 [-0.372, -0.158] | -0.299 [-0.449, -0.154] | yes |
| H5 | c(GIR,MW) weak cy | +0.350 [+0.287, +0.418] | +0.425 [+0.317, +0.518] | yes |
| H5 | c(GIR,MW) weak cx | +0.361 [+0.284, +0.437] | +0.402 [+0.322, +0.487] | yes |
| H6 | tetrad t1 weak cy | +0.120 [+0.089, +0.153] | +0.164 [+0.108, +0.219] | yes |
| H6 | tetrad t2 weak cy | +0.134 [+0.099, +0.173] | +0.158 [+0.101, +0.215] | yes |
| H6 | tetrad t3 weak cy | +0.014 [-0.003, +0.033] | -0.006 [-0.025, +0.013] | yes |

## Secondary (descriptive)

- **Sample.** August: 1,933 point-support window-hours (966/967 per stratum; July 2,402), median radar
  peak 13.1 mm h⁻¹ (July 15.8), median 8 gauges per window (July 9); 2,295 gauge stations (July 2,305).
- **Magnitudes.** The area inversion is similar on centroid x (+0.189 vs +0.193) and larger on centroid y
  (+0.259 vs +0.221). The point-support reversal and the GOES–microwave error correlation are somewhat
  larger in August (weak stratum c(GIR, MW) +0.425 / +0.402 vs +0.350 / +0.361). t3 is −0.006
  [−0.025, +0.013] (July +0.014). As in July, t3 ≈ 0 is compatibility with a single two-group split, not
  proof that no other dependence exists.
- **Sensitivity (T7-S).** In July the ranking's point estimate flipped at an assumed radar–gauge error
  correlation of 0.425–0.5. In August it does not flip within the grid (0 to 0.5): at ρ_RG = 0.5 the
  point estimate of Δ is still just negative (−0.017 to −0.002), but its interval includes reversal from
  0.375–0.4 onward. The GOES–microwave error correlation stays detected through 0.5, the top of the tested
  range (weak centroid y +0.299 [+0.088, +0.446]); values above 0.5 are untested. Expectation
  (≥ 0.3) met. As before, this bounds how far the assumption can move; it does not measure the actual
  radar–gauge correlation.
- **Other Probe A2 variables (not scored).** Dispersion: nf R 0.558, GIR 0.852, MW 0.199 (microwave rated
  least noisy, as in July). Wet fraction and log mean unstable, as in July.

## Implementation and deviations

1. **Period switch.** `scripts/period.py` (`PRECIP_PERIOD=jul|aug`) supplies the hours and days; it replaced
   the hard-coded July dates in `probe_a.py` and the five fetch scripts. The data root is selected by the
   working directory: August runs from `runs/aug/`, where `data` links to `data_aug/`. (The plan said the
   variable would also select the root; using the working directory left every path string in the analysis
   code unchanged.) After the switch, July reproduced **byte-for-byte**: `probe_a2_results.json`,
   `t7_results.json`, `t7r_results.json`, `t7s_results.json`, `t7_points.csv`, `probe_a2_windows.csv`.
   The analysis scripts `probe_a2.py`, `t5.py`, `t7.py`, `t7_robust.py`, `t7s.py` have the hashes frozen in the plan;
   `probe_a.py` changed only in its HOURS line.
2. **Missing producer.** No saved script wrote `data/stations.csv`. It is the unique gauge stations sorted by
   (src, sid); `scripts/make_stations.py` reproduces the July file byte-for-byte and wrote August's.
3. **H1 interval.** The Probe A2 analysis stores per-observer intervals, not one for Δ. `scripts/t10_area.py`
   runs `probe_a.analyse` unchanged with Δ added to the quantities its existing `quad` hook returns, so Δ uses
   the same bootstrap draws (seed 11, 500). On July all stored values agree to 1e-14 (CSV round-trip).
4. **Scoring.** `scripts/t10_score.py` implements the plan's rules; on July it returns Replicated (in-sample check).
   The stability gate uses the recorded negative-variance shares (area triplet; T7 triplets in both strata);
   T7-R's mirror triplet records none, in July or August.
5. **Fetch efficiency.** `scripts/fetch_imerg_both.py` downloads each IMERG file once and writes both
   outputs with the unchanged crop functions; on a July file its arrays are identical to the stored ones.
   No partial files were left by the stopped single-product fetchers (every file loaded and checked).
6. **Area-sample reference.** The plan quoted July's area sample as 14,751 (the T9 W0 sample); the Probe A2
   sample, which also requires the merged product, is 14,585. The 5,000 gate is unaffected.
7. **Unstated details.** The independent audit (`audit/AUDIT.md`) found four details fixed by the code, not the
   plans (≥ at the wet threshold, float32 window sums, T7 not requiring the merged product, the median tie rule).
   They apply identically to both periods because the code is the same.

## Data (reduced inputs; digest = SHA-256 over sorted file names and contents, first 16)

| path | files | size | digest |
|---|---|---|---|
| data_aug/mrms/ | 504 | 229 MB | f77e0243e4895db7 |
| data_aug/imerg/ | 1008 | 988 MB | d24aecdedf7c3ea9 |
| data_aug/imerg_x/ | 1008 | 138 MB | bbdf2afc5a5c6241 |
| data_aug/goes/ | 504 | 320 MB | 9db03850a89f7f60 |
| data_aug/gauges.csv | 1 | 55.1 MB | 5f17e6d3f214b021 |
| data_aug/stations.csv | 1 | 0.1 MB | bfd5529bfa27ea04 |

All 504 GOES hours had ≥ 3 scans (none dropped; July dropped 2).

## Reading

The July results replicate on an untouched period under a frozen pipeline: the inversion with two
satellites as legs, its reversal when gauges replace a satellite, the mirror triplet, the GOES–microwave
error correlation and the tetrad pattern. The effect is therefore not specific to July 2024 (the two months are adjacent and both warm-season). It
does not change what the four-observer data can identify: the satellite attribution remains conditional on
radar and gauge errors being less correlated than about 0.4 to 0.5, and both months are warm-season CONUS.
Replication across adjacent months with the same instruments, processing and season shows repeatability, not
attribution: a persistent measurement or representativeness effect would replicate too.
