# Probe A2: structural collocation with raw GOES infrared (frozen 2026-10-15, before any GOES value was read)

Identical to [PROBE_A_PLAN.md](PROBE_A_PLAN.md) except the infrared observer.
Reason: Probe A's IMERG infrared estimate (PDIR) is calibrated to the microwave
estimates (IMERG V07 ATBD), so (R, IR, MW) was not an independent triplet
(`PROBE_A_RESULTS.md`).

## Changed observer

- **GIR** (replaces IR): GOES-16 ABI L2 Cloud and Moisture Imagery, CONUS sector,
  band 13 (10.3 µm) brightness temperature Tb, from `noaa-goes16` on AWS. No
  microwave or rain-rate calibration.
- Per hour (accumulation ending t): the four scans starting at minutes 01, 16, 31,
  46 of hour t − 1; per-pixel Tb averaged over the scans present (≥ 3 required).
- Regridded to the 0.1° MRMS/IMERG grid by averaging all 2 km ABI pixels whose
  centres fall in a cell; a cell needs ≥ 4 pixels, else missing.
- Converted to a non-negative cold-cloud field: **GIRc = max(0, 235 − Tb)** (K),
  the standard 235 K convective-cloud threshold. "Wet" for GIR = GIRc > 0. All
  structural reads use GIRc exactly as they used rain rate for the other observers.

R, MW and M are unchanged. The observer set is now physically independent for
(R, GIR, MW): radar backscatter, raw cloud-top infrared, microwave. M (merged
IMERG) is built from MW and from PDIR, which in turn is built from geostationary
IR, so M shares inputs with both GIR and MW.

## Sample, variables, estimators, bootstrap, stop rule

Exactly as in PROBE_A_PLAN.md, with GIR in place of IR (GIR must be complete in
the window; GIRc window sum ≥ 1 K-cell). The value-control variable `logmean`
for GIR is log(mean GIRc + 0.1).

## Predictions (unchanged wording, re-tested)

1. Radar has the smallest noise fraction of R, GIR, MW for most structural variables.
2. GIR noise fraction > MW noise fraction.
3. Quadruple check flags M as error-correlated with its own inputs; largest
   expected for M–MW (MW-covered windows), with M–GIR positive but smaller.

Added check (descriptive, stated now): the share of bootstrap resamples with a
negative error variance should fall relative to Probe A if the IMERG-IR
dependence was the cause of Probe A's instability on `wet` and `logmean`.
