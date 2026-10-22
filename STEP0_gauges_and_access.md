# Step 0: gauge independence and data access (2026-10-15)

Desk check plus link and file-listing checks. No bulk data downloaded.

## Access (all checked today)

| Source | Check | Result |
|---|---|---|
| IMERG Late V07 half-hourly (NASA GES DISC) | 16-byte range request with the Earthdata token | HTTP 206, valid HDF5 header: **token works** |
| MRMS on AWS (`noaa-mrms-pds`) | listed `RadarOnly_QPE_01H`, `RadarAccumulationQualityIndex_01H`, `GaugeInflIndex_01H_Pass2`, `MultiSensor_QPE_01H_Pass2` for 2024-07-01 | all present; no login |
| USCRN hourly (NCEI) | https://www.ncei.noaa.gov/pub/data/uscrn/products/hourly02/2024/ | HTTP 200; 159 station files (about 130 in CONUS after removing AK 24, HI 2, ON 1) |
| ASOS hourly precipitation (Iowa Environmental Mesonet) | `asos.py` request for DSM, 1–2 Jul 2024, field `p01i` with lat/lon | CSV returned |

## Which gauges are already inside the products

- MRMS gauge-based and gauge-corrected products ingest roughly 10,000 gauges an
  hour from **HADS and MADIS**, after automated QC (Qi et al. 2016; Martinaitis
  et al. 2021; NWS WDTD gauge QC page).
- Stage IV RFC analyses merge radar with HADS gauges (and METAR where timelier),
  with manual QC at the RFCs.
- **MRMS RadarOnly QPE is not gauge-adjusted.** IMERG Late uses no same-time
  gauges (climatological calibration only; from memory, not re-verified).

## Consequence for the design

Finding a gauge network that *nothing* ingests is not needed if the target is
re-stated:

- **Target:** error of MRMS **radar-only** QPE at a gauge, |radar-only − gauge|,
  hourly.
- **Predictors:** structural disagreement between radar-only and IMERG Late
  (the Ω read), against the incumbents: MRMS `RadarAccumulationQualityIndex`
  and a plain pixel difference |radar-only − IMERG|.
- None of the predictors uses same-time gauges, so **every** gauge is
  independent of them: ASOS (several hundred CONUS stations, public via IEM)
  and USCRN (~130 CONUS stations, triple-redundant weighing gauges, highest
  quality).
- What this gives up: it tests error in radar-only rainfall, not in the
  gauge-merged analyses that insurers and hydrologists usually use. Testing the
  merged products would need a truly withheld network; whether USCRN is
  ingested by MADIS was not established.

Known label issues: a gauge is a point and a pixel is 1 km (radar) or ~10 km
(IMERG), which adds noise to every predictor equally; ASOS site coordinates are
sometimes off by a grid box (an OU REU study found 79% of ASOS stations sat in a
different MRMS box after a location correction); heated tipping buckets undercatch
in wind.

## Download size (corrected)

The earlier ~0.2 MB/file figure came from a dry January hour. Summer hourly files
on 2024-07-01: radar-only 1.2 MB, radar quality index 3.5 MB, gauge influence
index 4.8 MB, multi-sensor 1.0 MB (top-of-hour files; radar-only is also
published every 2 minutes, which we skip).

| Pull | 3 weeks (504 h) | Full season (Jun–Aug, 2,208 h) |
|---|---|---|
| radar-only + radar quality index (needed) | ~2.4 GB | ~10 GB |
| + gauge influence + multi-sensor (optional) | +~2.9 GB | +~13 GB |
| IMERG Late cropped to CONUS (estimate) | <0.5 GB | ~1–2 GB |
| gauges | negligible | negligible |

Files can be cropped to small windows around gauges right after download and the
originals deleted, so disk use stays small; the cost is bandwidth.

## Go / no-go

**GO** for feasibility, with the re-stated target. Next: pull 3 weeks of
radar-only + quality index + IMERG Late + gauges, regrid to a common grid, run the
signal check and `independence_audit.py` on the radar vs satellite structural
reads, then a first look at Ω vs RQI vs plain difference against gauge error.

## Sources

- MRMS gauge inputs (HADS, MADIS) and QC: https://vlab.noaa.gov/web/wdtd/-/gauge-quality-control ; Qi et al. 2016: https://repository.library.noaa.gov/view/noaa/32147 ; MRMS overview: https://www.weather.gov/media/ctp/Advanced%20Radar%20Class%20Material/MRMS%20Deeper%20Dive.pdf
- Stage IV inputs: https://www.emc.ncep.noaa.gov/mmb/SREF/pcpanl/ ; HMT/DTC observation description: https://verif.rap.ucar.edu/eval/hmt/2010/obs.php
- ASOS location errors vs MRMS grid (OU REU): https://caps.ou.edu/reu/reu16/finalpapers/Harkema-finalpaper.pdf
- USCRN hourly: https://www.ncei.noaa.gov/pub/data/uscrn/products/hourly02/
- IEM ASOS service: https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py
