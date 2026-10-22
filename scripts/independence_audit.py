"""
independence_audit.py  —  Gate 6 (Independent) of the Domain Qualification Rubric.

Self-contained. Depends only on numpy. This folder brings over NOTHING from the
rest of the corpus: the correlation-audit logic below is reimplemented fresh from
the Parallax method (the same idea used in the Markets ten-room audit and
specified in the meta-brief's independence-audit section), so there is zero
coupling to any other folder's code.

Two things live here:

    independence_audit(observables) -> dict
        Given N named observable series (or flattened spatial maps), returns the
        independence score I = mean off-diagonal |Pearson r|, the max |r|, the
        flagged redundant pairs (|r| > 0.85), a drop-highest-pair sensitivity
        check, and a PASS / PARTIAL / FAIL gate verdict.

    omega(centroids) -> float
        Given N centroid positions, returns Omega = mean pairwise Euclidean
        distance. This is the disagreement scalar itself — computed only once a
        domain has cleared the gate. Included so the folder is complete, not so
        you run it on inputs that failed the audit.

    tetrad_audit(observables, clusters=None) -> dict          (added 2026-10-15)
        Shared-error check for four or more observers. independence_audit
        catches redundancy (high |r|); it cannot see two observers that share an
        ERROR relative to the common signal, because that leaves |r| moderate.
        Under one common signal with mutually independent errors, every tetrad
        difference of the correlation matrix (r_ab*r_cd - r_ac*r_bd, ...) is 0.
        A shared error between x and y inflates r_xy alone, so exactly the two
        tetrads that contain the pairing {x-y | z-w} move off zero and the third
        stays at 0. Bootstrap CIs (plain, or by cluster) decide "non-zero".

Run `python independence_audit.py` for a self-demo on synthetic data: one
genuinely-independent observable set (Gate 2 PASS) and one that is five
transforms of a single surface (Gate 2 FAIL).

Run `python independence_audit.py --tetrad` for the shared-error self-demo, and
`python independence_audit.py --csv FILE --cols A,B,C,D [--cluster COL]
[--split-median COL]` to run tetrad_audit on a CSV (stdlib csv, no pandas).
"""

from itertools import combinations
import numpy as np

# Verdict bands (from the meta-brief independence-audit spec).
I_PASS = 0.30          # I < 0.30  -> genuinely independent
I_FAIL = 0.60          # I >= 0.60 -> transform diversity, not field disagreement
PAIR_FLAG = 0.85       # |r| above this flags a pair as potentially redundant
SENSITIVITY_WARN = 0.20  # relative ΔI on drop-highest that should make you look


def _as_named_matrix(observables):
    """Normalize input into (names, X) where X is (N, M): N observables, M samples.

    Accepts a dict {name: array} or a sequence of arrays. Each observable is
    flattened to 1D; all must share the same length.
    """
    if isinstance(observables, dict):
        names = list(observables.keys())
        series = [np.asarray(observables[k], dtype=float).ravel() for k in names]
    else:
        series = [np.asarray(v, dtype=float).ravel() for v in observables]
        names = [f"obs{i}" for i in range(len(series))]

    if len(series) < 2:
        raise ValueError("Need at least 2 observables to audit independence.")
    lengths = {s.size for s in series}
    if len(lengths) != 1:
        raise ValueError(f"Observables have mismatched lengths: {sorted(lengths)}")
    # Fix 2026-10-15: a NaN or a constant (dead) series makes corrcoef NaN, which
    # used to be zeroed and counted as "perfectly independent", pulling I down.
    bad = [n for n, s in zip(names, series) if not np.all(np.isfinite(s))]
    if bad:
        raise ValueError(f"Non-finite values in observables {bad}; drop or impute first.")
    dead = [n for n, s in zip(names, series) if np.ptp(s) == 0]
    if dead:
        raise ValueError(f"Constant (dead) observables {dead} carry no signal; remove them.")
    return names, np.vstack(series)


def _abs_corr_matrix(X):
    """Absolute Pearson correlation matrix for rows of X, NaNs -> 0 (degenerate)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        C = np.corrcoef(X)
    C = np.abs(np.atleast_2d(C))
    C[~np.isfinite(C)] = 0.0
    return C


def independence_audit(observables):
    """Audit whether N observables are independent enough for Omega to mean
    'field disagreement' rather than 'transform diversity'.

    Returns a dict with: names, I, max_r, pairs (sorted high->low), flagged,
    sensitivity {dropped, I_after, delta_rel}, verdict, note.
    """
    names, X = _as_named_matrix(observables)
    n = len(names)
    C = _abs_corr_matrix(X)

    pairs = []
    for i, j in combinations(range(n), 2):
        pairs.append(((names[i], names[j]), float(C[i, j])))
    pairs.sort(key=lambda p: p[1], reverse=True)

    off_diag = [r for _, r in pairs]
    I = float(np.mean(off_diag)) if off_diag else float("nan")
    max_r = float(np.max(off_diag)) if off_diag else float("nan")
    flagged = [(pair, r) for pair, r in pairs if r > PAIR_FLAG]

    # Sensitivity: drop ONE observable from the highest-correlated pair, recompute I.
    sensitivity = {"dropped": None, "I_after": I, "delta_rel": 0.0}
    if n > 2 and pairs:
        (a, b), _ = pairs[0]
        keep = [k for k in names if k != b]           # drop the second of the top pair
        sub = {k: X[names.index(k)] for k in keep}
        _, Xs = _as_named_matrix(sub)
        Cs = _abs_corr_matrix(Xs)
        sub_off = [Cs[i, j] for i, j in combinations(range(len(keep)), 2)]
        I_after = float(np.mean(sub_off)) if sub_off else float("nan")
        delta_rel = abs(I_after - I) / I if I else 0.0
        sensitivity = {"dropped": b, "I_after": I_after, "delta_rel": float(delta_rel)}

    if I < I_PASS:
        verdict = "PASS"
        note = "Genuinely independent. Omega measures real field disagreement."
    elif I < I_FAIL:
        verdict = "PARTIAL"
        note = "Mixed. Flag it, and replace the most-correlated pair before trusting Omega."
    else:
        verdict = "FAIL"
        note = "Reads are one surface in disguise. Omega measures transform diversity, not the field."

    return {
        "names": names,
        "I": I,
        "max_r": max_r,
        "pairs": pairs,
        "flagged": flagged,
        "sensitivity": sensitivity,
        "verdict": verdict,
        "note": note,
    }


def omega(centroids):
    """Omega = mean pairwise Euclidean distance between N centroid positions.

    centroids: array-like of shape (N, D) (D is usually 2, in normalized coords).
    """
    pts = np.asarray(centroids, dtype=float)
    if pts.ndim != 2 or pts.shape[0] < 2:
        raise ValueError("Need at least 2 centroids, shape (N, D).")
    dists = [np.linalg.norm(pts[i] - pts[j])
             for i, j in combinations(range(pts.shape[0]), 2)]
    return float(np.mean(dists))


# ---------------------------------------------------------------------------
# Shared-error (tetrad) check, added 2026-10-15.
# Origin: precip-disagreement T7-R. GOES infrared and IMERG microwave are
# physically different sensors with moderate |r|, yet share an error relative
# to surface rain. independence_audit cannot see that; the tetrads can.
# ---------------------------------------------------------------------------

TETRAD_MIN_EFFECT = 0.03   # practical floor: |t| below this is not called non-zero
TETRAD_N_BOOT = 500


def _quad_tetrads(R, a, b, c, d):
    """The three tetrad differences for observers (a, b, c, d) on correlation R.

    Three ways to split four observers into two pairs (a "pairing"):
        P_ab|cd = r_ab*r_cd    P_ac|bd = r_ac*r_bd    P_ad|bc = r_ad*r_bc
    t1 = P_ad|bc - P_ab|cd   (excludes pairing ac|bd)
    t2 = P_ad|bc - P_ac|bd   (excludes pairing ab|cd)
    t3 = P_ab|cd - P_ac|bd   (excludes pairing ad|bc)
    With order (R, IR, MW, G) these are exactly T7-R's t1, t2, t3.
    """
    p_ab = R[a, b] * R[c, d]
    p_ac = R[a, c] * R[b, d]
    p_ad = R[a, d] * R[b, c]
    return np.array([p_ad - p_ab, p_ad - p_ac, p_ab - p_ac])


def _all_quad_tetrads(R, quads):
    return np.array([_quad_tetrads(R, *q) for q in quads])


def _heywood(R, x, y, z, w):
    """Implied squared loadings if the shared error sits in pair (x, y).

    Uses only correlations that hypothesis treats as clean (never r_xy).
    Any value outside (0, 1] makes that hypothesis inadmissible.
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        return {
            x: float(R[x, z] * R[x, w] / R[z, w]),
            y: float(R[y, z] * R[y, w] / R[z, w]),
            z: float(R[z, x] * R[z, w] / R[x, w]),
            w: float(R[w, x] * R[w, z] / R[x, z]),
        }


def tetrad_audit(observables, clusters=None, n_boot=TETRAD_N_BOOT,
                 min_effect=TETRAD_MIN_EFFECT, seed=0):
    """Shared-error audit for N >= 4 observers of one common signal.

    observables: dict {name: array} or sequence of arrays (same rules as
        independence_audit: equal length, finite, non-constant).
    clusters: optional array of cluster labels (one per sample); the bootstrap
        then resamples whole clusters. None -> plain bootstrap over samples.
    A tetrad is called non-zero when its 95% bootstrap CI excludes 0 AND
    |t| >= min_effect.

    Per quad, the pattern is read as:
        clean     all three tetrads ~0  (consistent with one signal + independent errors)
        pairing   exactly the two tetrads containing one pairing {x-y | z-w} are
                  non-zero, the one excluding it is ~0 -> one shared error in x-y
                  OR in z-w (four observers cannot separate complementary pairs)
        unresolved  any other pattern (several breaches, or a non-single-factor field)

    Across quads (N >= 5) a pair's support = share of the quads containing it that
    implicate it. A breach in x-y gives x-y support 1.0 and every other pair < 1.
    With N = 4 the complementary pair also scores 1.0; the Heywood check (implied
    squared loadings must lie in (0, 1]) is reported for both readings, and a
    reading that fails it is inadmissible.

    Returns dict: names, n, quads [{observers, t, ci, nonzero, pattern,
    implicated, heywood}], pair_support, identified, ambiguous, verdict, note.
    """
    names, X = _as_named_matrix(observables)
    n_obs, m = X.shape
    if n_obs < 4:
        raise ValueError("tetrad_audit needs at least 4 observables.")
    quads = list(combinations(range(n_obs), 4))
    R = np.corrcoef(X)
    T = _all_quad_tetrads(R, quads)                       # (Q, 3)

    rng = np.random.default_rng(seed)
    if clusters is None:
        groups = [np.array([i]) for i in range(m)]
    else:
        lab = np.asarray(clusters).ravel()
        if lab.size != m:
            raise ValueError("clusters must have one label per sample.")
        _, inv = np.unique(lab, return_inverse=True)
        groups = [np.where(inv == g)[0] for g in range(inv.max() + 1)]
    G = len(groups)
    B = np.empty((n_boot,) + T.shape)
    for k in range(n_boot):
        pick = rng.integers(0, G, G)
        idx = np.concatenate([groups[g] for g in pick]) if clusters is not None else pick
        with np.errstate(invalid="ignore", divide="ignore"):
            B[k] = _all_quad_tetrads(np.corrcoef(X[:, idx]), quads)
    lo = np.nanpercentile(B, 2.5, axis=0)
    hi = np.nanpercentile(B, 97.5, axis=0)
    nonzero = ((lo > 0) | (hi < 0)) & (np.abs(T) >= min_effect)

    # Which pairing (as index pairs within the quad) each tetrad excludes:
    # t1 excludes ac|bd, t2 excludes ab|cd, t3 excludes ad|bc.
    def pairings(a, b, c, d):
        return [((a, c), (b, d)), ((a, b), (c, d)), ((a, d), (b, c))]

    support = {pair: [0, 0] for pair in combinations(range(n_obs), 2)}
    quad_rows = []
    for qi, q in enumerate(quads):
        nz = nonzero[qi]
        excl = pairings(*q)
        if not nz.any():
            pattern, implicated, hey = "clean", None, None
        elif nz.sum() == 2:
            k = int(np.where(~nz)[0][0])               # the ~0 tetrad excludes pairing k
            pattern, implicated = "pairing", excl[k]
            (x, y), (z, w) = implicated
            hey = {f"{names[x]}-{names[y]}": {names[o]: v for o, v in _heywood(R, x, y, z, w).items()},
                   f"{names[z]}-{names[w]}": {names[o]: v for o, v in _heywood(R, z, w, x, y).items()}}
        else:
            pattern, implicated, hey = "unresolved", None, None
        for pair in combinations(q, 2):
            support[pair][1] += 1
            if implicated is not None and pair in implicated:
                support[pair][0] += 1
        quad_rows.append({
            "observers": [names[o] for o in q],
            "t": [float(v) for v in T[qi]],
            "ci": [[float(lo[qi, j]), float(hi[qi, j])] for j in range(3)],
            "nonzero": [bool(v) for v in nz],
            "pattern": pattern,
            "implicated": None if implicated is None else
                [f"{names[x]}-{names[y]}" for x, y in implicated],
            "heywood": hey,
        })

    pair_support = {f"{names[i]}-{names[j]}": (s / tot if tot else 0.0)
                    for (i, j), (s, tot) in support.items()}
    top = [p for p, v in pair_support.items() if v == 1.0]
    n_unres = sum(r["pattern"] == "unresolved" for r in quad_rows)
    n_clean = sum(r["pattern"] == "clean" for r in quad_rows)

    identified, ambiguous = [], []
    if n_obs == 4 and top:
        # Complementary pairs tie; keep the admissible reading(s).
        hey = quad_rows[0]["heywood"]
        ok = [p for p in top if all(0 < v <= 1 for v in hey[p].values())]
        if len(ok) == 1:
            identified = ok
        else:
            ambiguous = top
    else:
        identified = top

    if n_clean == len(quad_rows):
        verdict = "PASS"
        note = "All tetrads ~0: consistent with one signal seen through independent errors."
    elif n_unres:
        verdict = "FAIL"
        note = (f"{n_unres} quad(s) show no single-pair pattern: several shared errors, or the "
                "observers do not share one common signal. Independence is not established.")
    else:
        verdict = "FAIL"
        who = identified or ambiguous
        note = ("Shared error between " + " or ".join(who) +
                ". These observers are not independent relative to the signal; "
                "drop one, or replace it with a reader of a different physical mechanism.")
        if ambiguous:
            note += (" Four observers cannot separate complementary pairs and the Heywood "
                     "check does not decide; use physics or add a fifth observer.")

    return {"names": names, "n": int(m), "n_groups": G, "quads": quad_rows,
            "pair_support": pair_support, "identified": identified,
            "ambiguous": ambiguous, "verdict": verdict, "note": note}


def print_tetrad_report(report, max_quads=10):
    """Human-readable dump of a tetrad_audit result."""
    print(f"  samples: {report['n']}   bootstrap units: {report['n_groups']}   "
          "(* = 95% CI excludes 0 and |t| >= min_effect)")
    print("  quad (a,b,c,d): t1 = r_ad·r_bc − r_ab·r_cd, t2 = r_ad·r_bc − r_ac·r_bd, "
          "t3 = r_ab·r_cd − r_ac·r_bd")
    for row in report["quads"][:max_quads]:
        a, b, c, d = row["observers"]
        print(f"  quad ({a}, {b}, {c}, {d}):  pattern = {row['pattern']}"
              + (f"  ->  {' | '.join(row['implicated'])}" if row["implicated"] else ""))
        labels = [f"t1 (excludes {a}-{c} | {b}-{d})",
                  f"t2 (excludes {a}-{b} | {c}-{d})",
                  f"t3 (excludes {a}-{d} | {b}-{c})"]
        w = max(len(x) for x in labels)
        for lab, t, (l, h), nz in zip(labels, row["t"], row["ci"], row["nonzero"]):
            print(f"      {lab:<{w}s}  {t:+.3f}  [{l:+.3f}, {h:+.3f}]{'  *' if nz else ''}")
        if row["heywood"]:
            for hyp, l2 in row["heywood"].items():
                ok = all(0 < v <= 1 for v in l2.values())
                vals = ", ".join(f"{k}={v:.2f}" for k, v in l2.items())
                print(f"      if shared error in {hyp}: loading² {vals}"
                      f"  {'admissible' if ok else 'INADMISSIBLE (Heywood)'}")
    if len(report["quads"]) > max_quads:
        print(f"  ... {len(report['quads']) - max_quads} more quads")
    if len(report["names"]) > 4:
        hits = {p: v for p, v in report["pair_support"].items() if v > 0}
        if hits:
            print("  pair support (share of quads containing the pair that implicate it):")
            for p, v in sorted(hits.items(), key=lambda kv: -kv[1]):
                print(f"      {p:<20s} {v:.2f}")
    print(f"  verdict: {report['verdict']}  —  {report['note']}")


def _tetrad_demo():
    rng = np.random.default_rng(1)
    M = 3000
    signal = rng.normal(size=M)
    load = {"A": 0.8, "B": 0.7, "C": 0.75, "D": 0.6, "E": 0.65}
    errs = {k: rng.normal(size=M) for k in load}
    shared = rng.normal(size=M)

    def obs(keys, planted=()):
        out = {}
        for k in keys:
            e = np.sqrt(1 - load[k] ** 2) * errs[k]
            if k in planted:
                e = 0.6 * e + 0.5 * shared          # B and C share part of their error
            out[k] = load[k] * signal + e
        return out

    print("=" * 64)
    print("TETRAD DEMO A — four observers, independent errors  (expect PASS)")
    print("=" * 64)
    clean = obs("ABCD")
    print_tetrad_report(tetrad_audit(clean))
    print(f"  (independence_audit on the same set: I={independence_audit(clean)['I']:.3f})")

    print()
    print("=" * 64)
    print("TETRAD DEMO B — same, shared error planted in B–C  (expect FAIL, B–C | A–D)")
    print("=" * 64)
    planted = obs("ABCD", planted="BC")
    print_tetrad_report(tetrad_audit(planted))
    print(f"  (independence_audit on the same set: I={independence_audit(planted)['I']:.3f}, "
          f"verdict {independence_audit(planted)['verdict']}; I moves only "
          f"{independence_audit(clean)['I']:.3f} -> {independence_audit(planted)['I']:.3f}, "
          "nothing names B-C)")

    print()
    print("=" * 64)
    print("TETRAD DEMO C — five observers, shared error in B–C  (expect B–C identified)")
    print("=" * 64)
    print_tetrad_report(tetrad_audit(obs("ABCDE", planted="BC")), max_quads=0)


def _read_csv_columns(path, cols, extra=()):
    import csv
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    data = {c: np.array([float(r[c]) for r in rows]) for c in cols}
    more = {c: np.array([r[c] for r in rows]) for c in extra if c}
    return data, more


def _csv_main(argv):
    import argparse
    ap = argparse.ArgumentParser(description="tetrad_audit on CSV columns")
    ap.add_argument("--csv", required=True)
    ap.add_argument("--cols", required=True, help="comma-separated observer columns (>= 4)")
    ap.add_argument("--cluster", default=None, help="column of cluster labels for the bootstrap")
    ap.add_argument("--split-median", default=None,
                    help="numeric column: also report the <= median and > median halves")
    ap.add_argument("--n-boot", type=int, default=TETRAD_N_BOOT)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    cols = a.cols.split(",")
    data, more = _read_csv_columns(a.csv, cols, [a.cluster, a.split_median])
    subsets = [("all", np.ones(len(data[cols[0]]), bool))]
    if a.split_median:
        v = more[a.split_median].astype(float)
        med = float(np.median(v))
        subsets += [(f"{a.split_median} <= {med:.3g}", v <= med),
                    (f"{a.split_median} > {med:.3g}", v > med)]
    for label, mask in subsets:
        print("=" * 64)
        print(f"{a.csv}  [{label}]")
        print("=" * 64)
        obs = {c: data[c][mask] for c in cols}
        clu = more[a.cluster][mask] if a.cluster else None
        print_tetrad_report(tetrad_audit(obs, clusters=clu, n_boot=a.n_boot, seed=a.seed))
        print()


def print_report(report):
    """Human-readable dump of an independence_audit result."""
    print(f"  independence score I (mean|r|): {report['I']:.3f}   max|r|: {report['max_r']:.3f}")
    print(f"  verdict: {report['verdict']}  —  {report['note']}")
    s = report["sensitivity"]
    if s["dropped"] is not None:
        print(f"  sensitivity: drop '{s['dropped']}' -> I={s['I_after']:.3f} "
              f"(relative ΔI={s['delta_rel']*100:.0f}%"
              f"{'  ⚠ driving pair' if s['delta_rel'] > SENSITIVITY_WARN else ''})")
    if report["flagged"]:
        print("  flagged redundant pairs (|r|>0.85):")
        for (a, b), r in report["flagged"]:
            print(f"      {a:>14s} <-> {b:<14s}  |r|={r:.3f}")
    else:
        print("  no pair exceeds |r|=0.85")


def _demo():
    rng = np.random.default_rng(0)
    M = 2000

    print("=" * 64)
    print("DEMO A — five genuinely independent observables  (expect PASS)")
    print("=" * 64)
    independent = {
        "amplitude":    rng.normal(size=M),
        "onset_time":   rng.normal(size=M),
        "frequency":    rng.normal(size=M),
        "similarity":   rng.normal(size=M),
        "arrival_geom": rng.normal(size=M),
    }
    print_report(independence_audit(independent))

    print()
    print("=" * 64)
    print("DEMO B — five transforms of ONE surface  (expect FAIL)")
    print("=" * 64)
    base = rng.normal(size=M)
    noise = lambda s: base * s + 0.05 * rng.normal(size=M)  # tiny independent jitter
    degenerate = {
        "raw":       base,
        "scaled":    noise(2.0),
        "offset":    base + 3.0,
        "smoothed":  np.convolve(base, np.ones(5) / 5, mode="same"),
        "squared":   np.sign(base) * base ** 2,
    }
    print_report(independence_audit(degenerate))

    print()
    print("=" * 64)
    print("Omega sanity check")
    print("=" * 64)
    tight = [[0.01, 0.00], [-0.01, 0.02], [0.00, -0.01]]   # agreeing reads
    scattered = [[0.4, -0.3], [-0.35, 0.4], [0.05, -0.45]]  # disagreeing reads
    print(f"  Omega(agreeing reads)   = {omega(tight):.3f}  (near 0 -> coherent)")
    print(f"  Omega(scattered reads)  = {omega(scattered):.3f}  (large -> incoherent)")


if __name__ == "__main__":
    import sys
    if "--tetrad" in sys.argv[1:]:
        _tetrad_demo()
    elif "--csv" in sys.argv[1:]:
        _csv_main(sys.argv[1:])
    else:
        _demo()
