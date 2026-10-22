# Research synthesis: precipitation structure and multi-observer disagreement (2026-10-15)

Follows the feasibility null (`FEASIBILITY_RESULTS.md`). Sources and notes are in
`literature/INDEX.md`; the PDFs are in `literature/`. This file states what the
literature says, how it lines up with the program's established methods, and a
ranked list of probes. Each probe names its data, its incumbent and its stop rule.

## 1. What the literature says about the null

The feasibility label was hourly |radar − gauge| at a point. The radar-verification
literature has measured this directly:

- The area–point (gauge representativeness) term is a dominant part of radar–gauge
  differences at short accumulations and remains significant even at 4-day totals
  on a 4 km grid (Ciach & Krajewski 1999, via search; Ciach, Habib & Krajewski
  2003 in `literature/`).
- Radar and gauge errors are not rigorously uncorrelated (Ciach et al. 2003), so
  methods that separate the two error terms need care.
- Satellite products lose spectral power and coherence with MRMS below ~200 km and
  ~4 h (Guilloteau et al. 2021). The structural reads ran on 0.3°–1.0° windows,
  inside that regime.

Together these predict the measured pattern: the label is dominated by point
noise and magnitude, and satellite structure inside a ±0.7° window is near the
satellite's effective-resolution limit. The null is consistent with known
physics. It maps the operating window rather than closing the domain.

## 2. Mathematics available, mapped to program methods

| Literature tool | What it does | Program counterpart | Overlap / gap |
|---|---|---|---|
| **Triple / extended / quadruple collocation** (Stoffelen 1998; McColl 2014; Alemohammad 2015) | error variance of each of ≥ 3 observers **without truth**, from their pairwise covariances; ETC adds each observer's correlation with the unknown truth | Ω = disagreement among independent reads; `independence_audit.py` | TC's central assumption (zero error cross-covariance) is exactly what the independence gate checks. TC is applied to values (rain amounts); no study was found applying it to **structural coordinates** (centroid, coverage, organization). |
| **Product spread as uncertainty** (Tian & Peters-Lidard 2010) | spread among six satellite products as a global error map | mean Ω | This is mean-of-disagreement, which the program found collapses in pathology; the spread **of** disagreement (σΩ) is the program's distinct unit. |
| **Spectral coherence / effective resolution** (Guilloteau 2021) | the scales at which a product reproduces the reference | resolution-floor protocol | Same idea, different object: their coherence is between values; the floor protocol asks where a structural coordinate starts to carry signal. |
| **Optimal transport, Sinkhorn, PAD** (Francis 2025; Skok & Lledó 2024) | displacement-aware distance between rain fields | `om_loc` (centroid distance) is a crude transport cost | Mature and actively developed for forecast verification; competing here repeats the Biomechanical pattern. |
| **Wavelet SAD** (Buschow 2019, 2021) | scale, anisotropy, direction of structure | Core42 orientation and multiscale axes | Mature in verification; the anisotropy/direction reads are portable as observables. |
| **Organization indices** (Iorg, SCAI, COP, ROME; Mandorli & Stubenrauch 2024; Retsch 2020) | single numbers for how clustered convection is | Saltimbanques islands, void, dispersion, indispensability | Nine indices that "do not agree", each one read of one field. That is the transform-diversity situation the program's Gate 6 was built to diagnose. |
| **Convective-mode labels** (Smith et al. 2012; SVRIMG) | human classification of storm structure (cellular, QLCS, tropical) from radar | ISR reading of structure; the galaxy lesson on label type | A structural label assigned by eye, with documented ambiguity between lines and clusters. The incumbent is the trained eye plus organization indices. |
| **Persistent homology / Minkowski functionals** (Biasotti 2016; Bobrowski) | topological counts and persistence of rain structures | ISR islands and persistence | No benchmark against FSS/SAL was found; an open space, though without a named buyer. |
| **Multifractal cascades** (Lovejoy & Schertzer) | scale-invariant statistics of rain intensity | resolution floor, scale persistence | Explains why fine-scale structure is intermittent; mostly a null model here. |

## 3. Probes, ranked

Ranking weighs fit to the program's record (structural label, trained-eye
incumbent, physically independent observers), data in hand, cost, and how cleanly
a negative would read.

### Probe A: Structural collocation (data in hand, about a day)

- **Question:** treating radar, infrared, microwave and merged-satellite as four
  observers, can collocation recover each observer's **structural** error variance
  (centroid position, wet-area fraction, island count, dispersion) without any
  gauge, and does the zero-cross-covariance assumption hold for structural reads
  where it fails for values?
- **Why it fits:** removes the point-gauge label that drove the null. It is the
  formal version of the program's thesis (disagreement among independent observers
  as information), with an existing statistical theory to stand on.
- **Data:** the 3 weeks already reduced in `data/` (radar, IR, MW, merged IMERG).
  MW covers ~39% of sample windows.
- **Incumbent:** value-based TC (Alemohammad 2015); product spread (Tian 2010).
- **Checkable predictions:** (1) radar's structural error variance is smallest
  (highest resolution); (2) IR's error is larger than MW's (IR sees cloud tops, not
  rain); (3) quadruple collocation flags merged-IMERG ~ IR as error-correlated
  (IMERG is built from IR between overpasses). If (3) fails, the method is not
  seeing what it should and the probe stops.
- **Stop rule:** if structural TC gives negative or unstable error variances across
  bootstrap resamples (a known TC failure when assumptions break), record the
  operating window and stop.

### Probe B: Convective mode as a structural label (new domain, scorecard first)

- **Question:** do structural reads (islands, void, dispersion, elongation,
  indispensability) separate human-labelled Cellular vs QLCS vs Tropical storms,
  beyond the organization indices (Iorg, SCAI, COP, ROME) as incumbents? And do
  cases humans find hard (Smith 2012 §3: line vs cluster) read as structurally
  ambiguous?
- **Why it fits:** label is structural and human-assigned; the incumbent is a
  trained eye plus disagreeing indices; Saltimbanques' multi-mass machinery is
  built for exactly "how many organizations, and how do they relate".
- **Data:** SVRIMG, public: about 1,765 labelled tornado-centred radar images
  (1996–2017; Cellular 669, QLCS 659, Tropical 269), live server, MIT-licensed
  code. **Open questions before any build:** confirm the labels are human, not
  model output; check image licence/terms; a second physically independent
  observer would need GOES infrared (GridSat-B1, public) collocated by time and
  place.
- **Galaxy caution:** with radar alone, the Ω would be cross-construction-family
  (ISR style), not cross-sensor. Ambiguity between line and cluster may be a
  detection-at-resolution question (Gate 4 risk). Score it with RUBRIC.md before
  building.
- **Stop rule:** standard rubric gates; then increment over the four
  organization indices on a held-out year block.

### Probe C: Re-label at the scale the literature says is clean (cheap, mostly a check)

- Same reads, label = daily or event totals, and areal gauge averages where USCRN
  and ASOS are dense. The literature says the area–point term stays significant
  even at 4 days, so the expected gain is small. Worth one run only to confirm the
  operating-window map, not as a discovery bet.

### Probe D: Variance of disagreement across scale (exploratory, data in hand)

- Per event, wavelet coherence between radar and satellite across scales; read the
  spread of coherence across sub-areas (σ of disagreement, the program's durable
  unit) rather than its level. Incumbent: Guilloteau's spectral coherence. Value
  is a map of where and at what scale the products structurally agree; buyer
  unclear.

### Not recommended

- Forecast verification with transport, wavelet or neighbourhood scores: active,
  well-funded incumbents (Met Office, ECMWF, NCAR), with open software.
- Radar artifact geometry (blockage wedges, bright-band rings): real geometric
  signatures, but RQI and climatology-based Fourier detection already cover the
  main case; worth revisiting only with a radar-operator partner.

## 4. Recommendation

Run **Probe A** next: it uses data already on disk, needs no label, and tests the
program's core idea against the one body of mathematics built for the same
question. In parallel, desk-score **Probe B** with the rubric (Gates 2 and 4 are
the open questions: label provenance and detection-vs-structure).
