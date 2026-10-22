# Decision gate after T1-bis and T7-S (2026-10-15)

Required by T7S_PLAN.md: no further experiments until this is discussed.

## What the corrected controls support

1. **The empirical result stands unchanged:** collocation rankings depend on which observers
   form the triplet; replacing either satellite with gauges reverses the radar–microwave
   ranking at the same point support.
2. **The synthetic support is restored, with corrections (T1-bis).** With a mass-conserving
   blur, structural collocation recovers centroid and dispersion, with correct rankings, on
   the full canvas and in windows; a planted shared satellite error inverts the ranking (at
   σ_sh = 1.0 cell, not 0.5); five observers identify a planted shared pair, in either
   direction, 10/10 with no false alarms (under the simulation's assumptions). **F6 is withdrawn**: the earlier "smoothing
   creates shared error" came from the blur not conserving mass. (Centroid errors become uncorrelated on the
   full canvas; dispersion keeps a small residual, +0.021.) Windowing contributes a
   small shared error (centroid ≤ 0.02, dispersion up to 0.14).
3. **The attribution tolerates substantial radar–gauge dependence (T7-S).** The satellite
   dependence stays detected up to an assumed radar–gauge error correlation of 0.475–0.5; radar's
   ranking as least noisy is clear by interval up to 0.275–0.4, and its point estimate flips at
   0.425–0.5. (Corrected after review: the earlier "ranking survives to about 0.4" was point-based.)
   The data cannot fix that correlation. This bounds how far the assumption can move; it does not
   measure where the actual correlation lies.

## What is still conditional or open

Status as of the T10 addendum in brackets.

- The attribution still rests on radar–gauge errors being less correlated than ≈ 0.4. [Still conditional.
  August: dependence detected through 0.5, the top of the tested range; ranking clear to 0.375–0.4; no point
  flip by 0.5. Values above 0.5 untested.]
- Windowing effects in the real data are untested (synthetic says small for centroids). [T9: robust to the
  tested grid shift and window size; windowing in general not excluded.]
- The smoothness test (T5) used an approximate match and a changing sample. [T9: fixed sample, close match;
  gap unchanged.]
- One month, one region, warm season. No unseen-period replication yet. [T10: replicated on August 2024;
  still one region and one season.]
- The physical cause of the shared satellite error. [Open.]

## Which paper the controls point to

The evidence now supports **a methods paper on structural collocation with a dependence audit**,
with the precipitation case as its worked example: (a) structural collocation validated on
corrected synthetic controls; (b) a real case where a plausible ranking is driven by shared
error; (c) detection with a fourth observer, bounded by a sensitivity analysis, replicated on an unseen month (T10); (d) a
five-observer design that identifies the shared pair in simulation. The satellite-dependence
finding is its main empirical result, stated as conditional with a quantified tolerance.
The measurement-window story is a minor component (small effect), not the headline.

## Candidate next steps (as set at the gate; 1 and 2 since done)

1. Real-data window control and fixed-sample T5 (cheap; data on disk). [Done: T9.]
2. Freeze the corrected analysis; replicate on an unseen period (August 2024; new pull). [Done: T10, replicated.]
3. A real fifth observer of surface rain whose independence can be argued (disdrometers, or areal
   averages from dense surface networks; not lightning, which needs ice aloft and would likely side with the satellites),
   using the five-observer design S-F validated.
4. Lineage benchmark on real products with documented dependence.

## Addendum: T9 run (2026-10-15)

The cheap step was run as pre-registered (`T9_PLAN.md`, `T9_RESULTS.md`). Window grid: a half-
shifted grid and 25 × 25 windows leave the inversion (+0.199, +0.181 vs +0.224) and the
point-support reversal (−0.231, −0.210) intact. Smoothness on a fixed sample with a close
match (σ* = 0.4; radar 0.687 vs microwave 0.700 lag-1): gap +0.218 [0.192, 0.250], essentially
unchanged. Neither accounts for the inversion within the tested variants; windowing and
representativeness effects in general are not excluded (the 25 × 25 gap is 19% smaller). Remaining candidate steps: unseen-period replication,
a real fifth observer, lineage benchmark. Paused for discussion.

## Wording note (2026-10-15, after second external review)

"Exact match" is replaced by "close match" (lag-1 0.687 vs 0.700; σ* rounded from 0.42 to 0.40); "closed"
by "robust to the tested variants". Threshold counts are "unvalidated, excluded from interpretation" (not
re-tested on the corrected control). The T7-S bound and the five-observer result are stated as conditional on
their assumptions. No numbers changed. The plan files keep their frozen wording.

## Next step frozen (2026-10-15)

`T10_PLAN.md`: replication on August 2024 with the July code (hashes recorded), six primary
predictions and a scoring rule. Frozen before any August data was fetched. Fetch awaits go-ahead.

## Addendum: T10 run (2026-10-15)

Replicated. August 2024, frozen July pipeline (July reproduced byte-for-byte after the period switch; an
independent re-implementation also matched July). All six primary predictions met on both axes and required
strata; in August T7-S dependence stays detected through 0.5 (the top of the tested range) and the
ranking point estimate does not flip by 0.5, though its interval includes reversal from 0.375–0.4. Replication shows
repeatability, not attribution. Details and deviations: `T10_RESULTS.md`. Remaining candidate
steps: a real fifth observer, the lineage benchmark, another season or region. Paused for discussion.
