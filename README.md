# Structural collocation of precipitation observers

Triple collocation estimates each observer's error from the covariances of three observers, with no
ground truth, provided their errors are independent. This project applies it to the **geometry** of
the rain field (intensity-weighted centroid and dispersion in 1.5° windows) instead of rain amounts.
It uses four observers over the contiguous United States:

- MRMS radar
- GOES-16 infrared
- IMERG passive microwave
- rain gauges

The data are hourly, July 2024, replicated on August 2024.

- **Preprint:** [`PAPER.html`](PAPER.html) (also at
  [parallaxmetrology.com](https://www.parallaxmetrology.com/Precipitation_Structural_Collocation_PAPER.html))
- **Lab book:** [`LABBOOK.html`](LABBOOK.html). Every experiment in order, including dead ends.
- **Start here for review:** [`REVIEW.md`](REVIEW.md). Bounded claim, evidence chain, reproduction, open questions.

## The result, in its bounded form

1. **Which observer looks most accurate depends on which three you compare.** With radar, infrared
   and microwave, microwave is rated least noisy. At gauge locations, replacing either satellite
   with gauges reverses this, and radar comes out least noisy. Replicated on an unseen month under a
   frozen pipeline: 6 of 6 pre-registered predictions met.
2. **Four-observer checks (vanishing tetrads) fit a single split: {radar, gauges} vs {GOES, microwave}.**
   - Read as an error the two satellites share relative to surface rain, its correlation is about
     +0.35 to +0.51.
   - That dependence is physically plausible, since both satellites respond to ice and cloud aloft.
   - The contribution is a way to detect, localise and bound it without ground truth.
   - The reading is conditional. It holds while radar–gauge error correlation stays below about
     0.3–0.5 (sensitivity analysis T7-S).
3. **Synthetic control.** With known truth and a mass-conserving blur, structural collocation
   recovers centroid and dispersion rankings, and a planted shared error reproduces the inversion.
   With five observers, the shared pair is identified in simulation.
4. **A pre-registered negative.** Structural disagreement does not predict hourly radar error at
   gauges. This was the original question; see `FEASIBILITY_RESULTS.md`.

**Not claimed:**
- the physical cause of the dependence;
- that radar and gauge errors are in fact weakly correlated;
- anything beyond warm-season CONUS, hourly, 0.1°.

## How the work was done

- Each experiment has a plan (`*_PLAN.md`) frozen before its values were computed, and a results
  file (`*_RESULTS.md`) listing deviations. `GATE.md` records the decision points.
- **Independent re-implementation.** [`audit/AUDIT.md`](audit/AUDIT.md) was written from the plans
  alone, without reading the analysis code. It reproduced 101 July quantities to within 4 × 10⁻⁹.
  Its scripts keep the absolute paths of the machine it ran on.
- **Re-run check (2026-10-22).** Every analysis was re-run from the reduced inputs in a clean copy.
  All saved July and August outputs reproduced byte-for-byte, including T9, T1-bis and the T10
  scoring.
- The first synthetic control (T1) used a blur that did not conserve mass, and a finding built on it
  was withdrawn. Both are kept for the record. [`notes/synthetic-shared-error/`](notes/synthetic-shared-error/)
  explains the mechanism and gives a checklist.

## Layout

| path | contents |
|---|---|
| `*_PLAN.md`, `*_RESULTS.md` | frozen plans and results, one pair per experiment |
| `scripts/` | fetch, analysis and figure scripts; `independence_audit.py` (independence and tetrad audits, numpy only) |
| `results/july/`, `results/august/` | saved outputs: result JSONs, point and window tables |
| `figures/` | figures used in the paper and lab book |
| `build_paper.py`, `build_labbook.py` | build `PAPER.html` and `LABBOOK.html` from the result files |
| `audit/` | independent re-implementation and its comparison table |
| `notes/synthetic-shared-error/` | technical note with a minimal reproduction (`minimal.py`) |
| `literature/INDEX.md` | annotated reading list (PDFs not included) |
| `runs/aug/` | working folder for the August replication |

## Reproduce

Python 3.12 with the packages in `requirements.txt`:

```
python -m venv .venv && .venv/bin/pip install -r requirements.txt
```

**Data.** About 1.7 GB per month, not included. The fetch scripts stream the open archives and save
reduced arrays under `data/`:
- MRMS, GOES-16 and gauges (USCRN, ASOS) need no login.
- IMERG needs a free NASA Earthdata token. Put it in `secrets/earthdata_token.txt`, which is
  git-ignored. The scripts read it at run time and never print it.

```
# July 2024, from the repository root
.venv/bin/python scripts/fetch_gauges.py && .venv/bin/python scripts/make_stations.py
.venv/bin/python scripts/fetch_mrms.py && .venv/bin/python scripts/fetch_goes.py && .venv/bin/python scripts/fetch_imerg_both.py
.venv/bin/python scripts/probe_a2.py && .venv/bin/python scripts/t7.py && .venv/bin/python scripts/t7_robust.py && .venv/bin/python scripts/t7s.py
.venv/bin/python scripts/t9.py && .venv/bin/python scripts/t1bis.py
.venv/bin/python scripts/t10_area.py _july && .venv/bin/python scripts/t10_score.py _july

# August 2024 replication: same scripts from runs/aug, where data -> ../../data_aug
cd runs/aug && export PRECIP_PERIOD=aug
../../.venv/bin/python ../../scripts/fetch_gauges.py   # ...then the same fetch and analysis commands
../../.venv/bin/python ../../scripts/t10_area.py && ../../.venv/bin/python ../../scripts/t10_score.py
```

Compare your outputs with `results/`. The upstream archives can be reprocessed, so a fresh download
may not reproduce the input digests in `LABBOOK.html` exactly.

`scripts/t1bis.py` (synthetic control, about 20 minutes) needs no downloaded data.

To rebuild the paper from the saved results only, see [`results/README.md`](results/README.md).

## Cite

Parrish, R. (2026). *Structural collocation: estimating precipitation observer error from disagreement
in rain-field geometry, and detecting shared error without ground truth.* Preprint, Parallax Metrology.
https://www.parallaxmetrology.com/Precipitation_Structural_Collocation_PAPER.html

## License

Code and documents: MIT ([`LICENSE`](LICENSE)). The data come from NOAA, NASA, NCEI and the Iowa
Environmental Mesonet under their own terms.
