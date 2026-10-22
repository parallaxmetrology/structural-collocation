"""T1-bis (T1BIS_PLAN.md): synthetic assumption audit with mass-conserving blur, full canvas vs windows."""
import json, sys
import numpy as np, pandas as pd
from scipy import ndimage

sys.path.insert(0, "../domain-qualification")
from independence_audit import tetrad_audit  # read-only use

N, SEED, C = 20000, 20261016, 64
YY, XX = np.mgrid[0:C, 0:C] - (C - 1) / 2
OBS = {"R": (0.5, 0.0, 0.3), "S1": (0.8, 1.0, 0.4), "S2": (1.5, 2.0, 0.6),
       "Gs": (1.0, 0.0, 0.5), "Gb": (1.0, 1.5, 0.5), "H": (0.7, 0.0, 0.4)}
c0 = C // 2 - 7  # central 15x15 window rows/cols c0..c0+14 (centre at cell 31.5-0.5)
REGIONS = {"full": (slice(0, C), slice(0, C)),
           "win": (slice(c0, c0 + 15), slice(c0, c0 + 15)),
           "win_s1": (slice(c0 + 3, c0 + 18), slice(c0 + 3, c0 + 18)),
           "win_s2": (slice(c0 - 3, c0 + 12), slice(c0 + 4, c0 + 19))}
VARS = ["cy", "cx", "d", "d2"]


def render(cy, cx, s, I):
    f = np.zeros((C, C))
    for a, b, ss, ii in zip(cy, cx, s, I):
        f += ii * np.exp(-((YY - a) ** 2 + (XX - b) ** 2) / (2 * ss ** 2))
    return f


def blur(f, sb):
    return ndimage.gaussian_filter(f, sb, mode="constant", truncate=4.0) if sb > 0 else f


def reads(f, reg):
    sl = REGIONS[reg]
    A = f[sl]; Y = YY[sl]; X = XX[sl]
    yc, xc = Y.mean(), X.mean()
    m = A.sum()
    if m < 1:
        return None
    cy = (A * (Y - yc)).sum() / m; cx = (A * (X - xc)).sum() / m
    d2 = (A * ((Y - yc - cy) ** 2 + (X - xc - cx) ** 2)).sum() / m
    return {"cy": cy, "cx": cx, "d": np.sqrt(d2), "d2": d2}


def simulate():
    rng = np.random.default_rng(SEED)
    W = []
    for _ in range(N):
        K = 1 + rng.poisson(1.5)
        cells = dict(cy=rng.uniform(-5, 5, K), cx=rng.uniform(-5, 5, K), s=rng.uniform(0.8, 2.5, K), I=rng.lognormal(1.5, 0.8, K))
        err = {o: dict(dy=rng.normal(0, sd, K), dx=rng.normal(0, sd, K), m=rng.lognormal(0, li, K)) for o, (sd, _, li) in OBS.items()}
        W.append((cells, err, rng.normal(0, 1, 2), rng.normal(0, 1, 2)))
    return W


def build(W, shared_pair=None, sig=0.0, regions=("full", "win", "win_s1", "win_s2"), obs=("R", "S1", "S2", "Gs", "Gb")):
    rows = []
    for k, (cells, err, u1, u2) in enumerate(W):
        T = render(cells["cy"], cells["cx"], cells["s"], cells["I"])
        fields = {"T": T}
        for o in obs:
            sd, sb, li = OBS[o]; e = err[o]
            sh = np.zeros(2)
            if shared_pair and o in shared_pair:
                sh = (u1 if shared_pair == ("S1", "S2") else u2) * sig
            fields[o] = blur(render(cells["cy"] + e["dy"] + sh[0], cells["cx"] + e["dx"] + sh[1], cells["s"], cells["I"] * e["m"]), OBS[o][1])
        for reg in regions:
            rr = {n: reads(f, reg) for n, f in fields.items()}
            if any(v is None for v in rr.values()):
                continue
            row = {"k": k, "reg": reg}
            for n, v in rr.items():
                for var in VARS:
                    row[f"{n}_{var}"] = v[var]
            rows.append(row)
    return pd.DataFrame(rows)


def errcorr(d, a, b, var):
    T = d[f"T_{var}"].values
    res = {}
    for o in (a, b):
        x = d[f"{o}_{var}"].values; p = np.polyfit(T, x, 1); res[o] = x - np.polyval(p, T)
    return float(np.corrcoef(res[a], res[b])[0, 1])


def tc(d, var, trip=("R", "S2", "S1")):
    X = d[[f"{o}_{var}" for o in trip]].values; Cm = np.cov(X.T)
    nf = {}
    for i, o in enumerate(trip):
        j, k = [x for x in range(3) if x != i]
        nf[o] = 1 - Cm[i, j] * Cm[i, k] / (Cm[i, i] * Cm[j, k])
    T = d[f"T_{var}"].values
    true = {o: 1 - np.corrcoef(d[f"{o}_{var}"], T)[0, 1] ** 2 for o in trip}
    order = lambda dd: sorted(dd, key=dd.get)
    return {"est": {o: float(v) for o, v in nf.items()}, "true": {o: float(v) for o, v in true.items()},
            "max_abs_err": float(max(abs(nf[o] - true[o]) for o in trip)), "rank_ok": order(nf) == order(true),
            "est_order": order(nf), "true_order": order(true)}


if __name__ == "__main__":
    out = {}
    W = simulate()
    # S-A analytical reference
    shifts, rel = {1.0: [], 2.0: []}, {1.0: [], 2.0: []}
    for cells, *_ in W[:2000]:
        T = render(cells["cy"], cells["cx"], cells["s"], cells["I"]); r0 = reads(T, "full")
        for sb in (1.0, 2.0):
            r1 = reads(blur(T, sb), "full")
            shifts[sb].append(np.hypot(r1["cy"] - r0["cy"], r1["cx"] - r0["cx"]))
            rel[sb].append(abs(r1["d2"] - r0["d2"] - 2 * sb ** 2) / (2 * sb ** 2))
    out["S_A"] = {str(sb): {"median_centroid_shift": float(np.median(shifts[sb])), "max_centroid_shift": float(np.max(shifts[sb])),
                            "median_rel_err_d2": float(np.median(rel[sb]))} for sb in (1.0, 2.0)}
    out["S_A"]["pass"] = all(out["S_A"][s]["median_centroid_shift"] < 0.01 and out["S_A"][s]["median_rel_err_d2"] < 0.01 for s in ("1.0", "2.0"))
    print("S-A", json.dumps(out["S_A"]), flush=True)
    if not out["S_A"]["pass"]:
        json.dump(out, open("data/t1bis_results.json", "w"), indent=1); sys.exit("S-A failed: rendering wrong, stopping per plan")
    # S-B / S-C / S-D (no planted error)
    d = build(W)
    out["n"] = {r: int((d.reg == r).sum()) for r in REGIONS}
    out["S_BC"] = {}; out["S_D"] = {}
    for reg in REGIONS:
        s = d[d.reg == reg]
        out["S_BC"][reg] = {var: {f"{a}-{b}": errcorr(s, a, b, var) for a, b in [("S1", "S2"), ("R", "S1"), ("R", "S2"), ("R", "Gs"), ("S1", "Gb"), ("S2", "Gb"), ("Gs", "Gb")]} for var in VARS}
        out["S_D"][reg] = {var: tc(s, var) for var in VARS}
        print("S-B/C", reg, "n", len(s), {v: round(out["S_BC"][reg][v]["S1-S2"], 3) for v in VARS}, flush=True)
        print("  S-D", reg, {v: (round(out["S_D"][reg][v]["max_abs_err"], 3), out["S_D"][reg][v]["rank_ok"]) for v in VARS}, flush=True)
    # S-E planted shared displacement S1-S2
    out["S_E"] = {}
    for sig in (0.5, 1.0):
        de = build(W, ("S1", "S2"), sig, regions=("full", "win"), obs=("R", "S1", "S2"))
        out["S_E"][str(sig)] = {}
        for reg in ("full", "win"):
            r_ = tc(de[de.reg == reg], "cy")
            inv = r_["est"]["S1"] < r_["est"]["R"] and r_["true"]["R"] < r_["true"]["S1"]
            out["S_E"][str(sig)][reg] = {**r_, "inverted": inv}
            print("S-E", sig, reg, {k: round(v, 3) for k, v in r_["est"].items()}, "true", {k: round(v, 3) for k, v in r_["true"].items()}, "inverted", inv, flush=True)
    # S-F benchmark
    out["S_F"] = {}
    for name, pair in [("none", None), ("S1-S2", ("S1", "S2")), ("R-Gs", ("R", "Gs"))]:
        df = build(W, pair, 1.0 if pair else 0.0, regions=("full",), obs=("R", "S1", "S2", "Gs", "H"))
        flags = []
        for b in range(10):
            blk = df.iloc[b * 2000:(b + 1) * 2000]
            res = tetrad_audit({o: blk[f"{o}_cy"].values for o in ["R", "S1", "S2", "Gs", "H"]}, n_boot=200, seed=b)
            flags.append({"identified": res["identified"], "verdict": res["verdict"]})
        target = None if pair is None else f"{pair[0]}-{pair[1]}"
        hit = sum(target in f["identified"] for f in flags) / 10 if target else None
        fa = sum(f["verdict"] != "PASS" for f in flags) / 10 if target is None else sum(any(p != target for p in f["identified"]) for f in flags) / 10
        out["S_F"][name] = {"hit_rate": hit, "false_alarm_rate": fa, "blocks": flags}
        print("S-F", name, "hit", hit, "false alarm", fa, [f["identified"] for f in flags][:3], flush=True)
    json.dump(out, open("data/t1bis_results.json", "w"), indent=1, default=str)
