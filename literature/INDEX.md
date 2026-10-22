# Literature index (pulled 2026-10-15)

Open-access or author-hosted copies, for local research reading. PDFs and the
extracted `txt/` are git-ignored: several are author copies of copyrighted
articles and should not be redistributed. One paragraph per paper: what it is and
why it is here.

## Disagreement among independent observers (the formal maths)

- **TC_Alemohammad2015_multiplicative_triple_collocation_HESS.pdf** (CC BY 3.0).
  Triple collocation (TC) estimates each product's error variance from three
  collocated products with no ground truth. Applied to NEXRAD radar, TRMM, GPCP and
  GPI over the central US (biweekly, 2° × 2°), with a multiplicative (log) error
  model. Key cautions: TC assumes zero error cross-covariance between products and
  zero covariance between errors and truth; violating the first biases TC RMSE low.
  doi:10.5194/hess-19-3489-2015
- Not downloaded (link returned HTML): McColl et al. 2014, *Extended triple
  collocation* (GRL): adds each product's correlation with the unknown truth
  (squared correlation = unbiased signal-to-noise). Wei et al. 2023 (GRL) applies
  an ETC variant to precipitation merging. Quadruple collocation (Pierdicca et al.)
  relaxes the independent-error assumption when a fourth product exists.
- Not downloaded (paywalled): Tian & Peters-Lidard 2010 (GRL), global uncertainty
  map from the **spread among six satellite products**; spread understates error
  where products share a bias. doi:10.1029/2010GL046008

## Why point-gauge labels are noisy (explains the feasibility null)

- **GAUGE_Habib_Krajewski_radar_error_11.pdf**: Ciach, Habib & Krajewski 2003,
  *Zero-covariance hypothesis in the error variance separation method*. Tests the
  assumption that radar and gauge area-rainfall errors are uncorrelated; finds it
  is not rigorously met. The area–point (gauge representativeness) term is large;
  radar errors are comparable to or larger than single-gauge errors.
- **GAUGE_Habib_Krajewski_gauge_representativeness_9.pdf**: Habib, Ciach &
  Krajewski 2004, *conditional distribution transformation* to filter gauge
  representativeness errors out of radar–gauge verification samples.
- Search summary (not downloaded): Ciach & Krajewski 1999 found the area–point
  component dominates radar–gauge differences at short accumulations and remains
  significant even at 4 days (4 × 4 km grid).

## Scale: where radar and satellite agree

- **SPEC_Guilloteau2021_spacetime_spectra_satellite_vs_MRMS.pdf**: Guilloteau,
  Foufoula-Georgiou, Kirstetter, Tan & Huffman 2021 (J. Hydrometeor.). Wavenumber–
  frequency spectra and coherence of five satellite products vs GV-MRMS: all are
  too smooth below ~200 km and ~4 h, with low coherence there; IMERG-Final best.
  Defines "effective resolution" from coherence.
- **MULTIFRAC_hess2014_Belgian_radar_stable_scaling.pdf**: van den Berg, Delobbe &
  Verhoest 2014 (HESS, CC BY 3.0). Imperfect scaling of radar rainfall
  distributions across scales (α-stable modelling).
- Search summary: Lovejoy & Schertzer multifractal cascades (TRMM radar scaling
  over 5–20,000 km; codimension of mean field C1 ≈ 0.63, Lévy index α ≈ 1.5;
  the "zeros" support problem).

## Spatial verification incumbents (forecast vs observation)

- **OT_Francis2024_unbalanced_OT_Sinkhorn_verification_arXiv.pdf**: Francis,
  Cotter & Mittermaier (Met. Apps 2025). Unbalanced entropic optimal transport /
  Sinkhorn divergence; robust to double penalty and on average aligns with expert
  scores on ICP cases.
- **OT_SkokLledo2024_precip_attribution_distance_global_arXiv.pdf**: Skok & Lledó
  (QJRMS). Precipitation Attribution Distance, an upper bound on Wasserstein
  distance, on the sphere; 7 years of ECMWF IFS.
- **WAV_Buschow2019_wavelet_verification_GMD.pdf** (CC BY 4.0): wavelet-based
  structure verification tested on a stochastic rain model; extended by the SAD
  method (scale, anisotropy, direction; QJRMS 2021; R package `sad`).
- **OBSUNC_Allison2015_radar_ensemble_FSS_EMS_abstract.pdf** (CC BY 3.0): Met
  Office radar ensemble changes FSS by ~10% of the single-field score.

## Organization of convection (a structural property with no agreed metric)

- **ORG_MandorliStubenrauch2024_organization_indices_GMD.pdf** (CC BY 4.0): nine
  object-based organization indices compared; "no consensus on the best method";
  SCAI ≈ object count, ROME/COP ≈ object size, Iorg ≈ clustering.
- **ORG_Retsch2020_convective_organization.pdf**: ROME (Radar Organization
  Metric), defined from object size and pairwise proximity on Darwin radar.

## Convective mode: a human-assigned structural label

- **MODE_Smith2012_SPC_convective_mode_climatology_WAF.pdf**: SPC manual radar
  classification of 22,901 tornado and significant-severe grid-hour events
  (2003–11) into supercell / QLCS / disorganized and sub-modes. Section 3 on
  classification difficulty: unambiguous discrimination between clusters and lines,
  and between closely spaced cells and clusters, "was not always possible".
- **MODE_Thompson2012_convective_mode_environment_WAF.pdf**: environments of the
  same sample.
- **MCS_HaberlieAshley2018_partI_segmentation_JAMC.pdf**,
  **MCS_HaberlieAshley2018_partII_tracking_JAMC.pdf**: MCS segmentation,
  classification (hand-labelled training data) and tracking in CONUS radar mosaics.
- Dataset found (not a paper): **SVRIMG** (Haberlie et al.), radar reflectivity
  images centred on severe reports 1996–2017. Server live 2026-10-15: tornado
  index table 8.3 MB; classification table `sample_classifications_96-17.csv`
  58 KB; repo `ahaberlie/SVRIMG` (MIT) holds ~1,765 labelled tornado images:
  Cellular 669, QLCS 659, Tropical 269, Other 49, Missing 101, Noise 19. Whether
  these labels are human or model output still needs confirming from the SVRIMG
  dataset paper (BAMS/par.nsf.gov 10289754, download failed).
- Related: SVRIMG-DM (`jthielen/svrimg47`, Apache-2.0), detailed-morphology
  samples and a labelling app; PNNL FLEXTRKR MCS database (radar + IR, 2004–2017,
  no use constraints per ARM listing; pixel files on request).

## Topology of rain fields

- **TOPO_Biasotti2016_rainfall_maxima_topological_tracking.pdf** (encrypted PDF;
  text not extracted): topological extraction of prominent rainfall maxima and
  bottleneck matching to track them (Genoa radar and gauges).
- Search summary: no study found that uses persistent homology, Minkowski
  functionals or Euler characteristic as a formal precipitation verification score
  benchmarked against FSS/SAL. Theory for Gaussian-field Euler characteristic
  (Bobrowski & Borman) gives null baselines.

## Operational context

- **MRMS_DeeperDive_NWS_training.pdf**: MRMS product overview (2022), inputs,
  QPE passes, RQI.
- Search summary: RQI = RQI_blockage × RQI_height; static blockage term; height
  term uses freezing level and bright band; "does not fully capture radar QPE
  uncertainty" (Chen et al. 2013).
- Search summary (parametric insurance): Senegal index insurance, 65% of surveyed
  farmers in Kaolack reported gauge–satellite discrepancies; NZ insurer FMG paused
  a parametric product at 67% best trigger accuracy on 10 km grids.
- **IMERG_V07_ATBD_Huffman2023.pdf** (NASA): IMERG V07 algorithm description. IR estimate = PDIR, calibrated to passive-microwave estimates via a 30-day rotating accumulator; KF morphing combines PMW and IR.
- Search note (2026-10-15, extended): no structural/object-based triple collocation found; closest is binned TC for regime-dependent uncertainty, https://ieeexplore.ieee.org/document/9743452 .
