"""T9 (T9_PLAN.md): real-data window control (T9-W) and fixed-sample smoothness test (T9-S)."""
import json
import numpy as np, pandas as pd
import probe_a as P
import t5 as T5
import t7 as T7

SEED = 21


def reads(A):
    n = A.shape[0]; yy, xx = np.mgrid[0:n, 0:n] - (n - 1) / 2; s = A.sum()
    return (A * yy).sum() / s, (A * xx).sum() / s


def tiles(r0, c0, W):
    for i in range((348 - r0) // W):
        for j in range((698 - c0) // W):
            yield i, j, (slice(r0 + W * i, r0 + W * (i + 1)), slice(c0 + W * j, c0 + W * (j + 1)))


def ok(w):
    if np.isfinite(w["R"]).mean() < 0.9 or not np.isfinite(w["IR"]).all() or np.isfinite(w["MW"]).mean() < 0.8:
        return None
    w = {k: np.nan_to_num(v) for k, v in w.items()}
    return w if min(v.sum() for v in w.values()) >= 1 else None


def clu(t, r0, c0, W, i, j):
    lat = 55.0 - 0.1 * (r0 + W * i + W / 2); lon = -130.0 + 0.1 * (c0 + W * j + W / 2)
    return f"{t:%Y-%m-%d}|{int(np.floor(lat / 5))}_{int(np.floor(lon / 5))}"


def nf_boot(d, cols, n_boot=300):
    rng = np.random.default_rng(SEED)
    c, uc = pd.factorize(d.clu); idx = [np.where(c == k)[0] for k in range(len(uc))]
    out = {}
    for v in ["cy", "cx"]:
        X = d[[f"{o}_{v}" for o in cols]].values
        def nf(rows):
            C = np.cov(X[rows].T)
            n0 = 1 - C[0, 1] * C[0, 2] / (C[0, 0] * C[1, 2]); n2 = 1 - C[2, 0] * C[2, 1] / (C[2, 2] * C[0, 1])
            return n0, n2
        a, b = nf(np.arange(len(d)))
        B = np.array([np.subtract(*nf(np.concatenate([idx[k] for k in rng.integers(0, len(idx), len(idx))]))) for _ in range(n_boot)])
        out[v] = {"nf_first": float(a), "nf_third": float(b), "diff": float(a - b), "ci": [float(np.percentile(B, 2.5)), float(np.percentile(B, 97.5))]}
    return {"n": int(len(d)), "clusters": int(len(uc)), **out}


def build_area(r0, c0, W, centred=False, blur=0.0, keyset=None):
    rows = []
    for t in P.HOURS:
        h = T5.load_hour(t)
        if h is None:
            continue
        R, GIR, MW = h
        Rb = T5.blur_nan(R, blur)
        for i, j, sl in tiles(r0, c0, W):
            if keyset is not None and (t, i, j) not in keyset:
                continue
            w = ok({"R": Rb[sl], "IR": GIR[sl], "MW": MW[sl]})
            if w is None:
                continue
            rc = reads(w["R"])
            if centred and max(abs(rc[0]), abs(rc[1])) > 2:
                continue
            row = {"t": t, "i": i, "j": j, "clu": clu(t, r0, c0, W, i, j)}
            for k, v in w.items():
                cy, cx = reads(v); row[f"{k}_cy"] = cy; row[f"{k}_cx"] = cx
            rows.append(row)
    return pd.DataFrame(rows)


def build_points(r0, c0, W):
    rows = []
    for t in P.HOURS:
        h = T5.load_hour(t)
        if h is None or t not in T7.G_BY_T:
            continue
        R, GIR, MW = h
        gt = T7.G_BY_T[t]
        gr, gc, gv = T7.st.br.values[gt.si.values], T7.st.bc.values[gt.si.values], gt.p_mm.values
        for i, j, sl in tiles(r0, c0, W):
            w = ok({"R": R[sl], "IR": GIR[sl], "MW": MW[sl]})
            if w is None:
                continue
            a0, b0 = r0 + W * i, c0 + W * j
            m = (gr >= a0) & (gr < a0 + W) & (gc >= b0) & (gc < b0 + W)
            if m.sum() < 5:
                continue
            lr, lc = gr[m] - a0, gc[m] - b0; y = lr - (W - 1) / 2; x = lc - (W - 1) / 2
            vals = {"R": w["R"][lr, lc], "MW": w["MW"][lr, lc], "G": gv[m]}
            if min(v.sum() for v in vals.values()) <= 0:
                continue
            row = {"t": t, "clu": clu(t, r0, c0, W, i, j)}
            for k, v in vals.items():
                row[f"{k}_cy"] = float((v * y).sum() / v.sum()); row[f"{k}_cx"] = float((v * x).sum() / v.sum())
            rows.append(row)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    out = {"W": {}, "W_points": {}}
    for name, (r0, c0, W, cen) in {"W0": (2, 2, 15, False), "W1": (9, 9, 15, False), "W2": (2, 2, 25, False), "W3": (2, 2, 15, True)}.items():
        d = build_area(r0, c0, W, centred=cen)
        out["W"][name] = nf_boot(d, ["R", "IR", "MW"])
        r = out["W"][name]
        print(f"T9-W {name} n={r['n']}: cy R-MW {r['cy']['diff']:+.3f} {np.round(r['cy']['ci'],3)}  cx {r['cx']['diff']:+.3f} {np.round(r['cx']['ci'],3)}", flush=True)
    for name, (r0, c0, W) in {"W1": (9, 9, 15), "W2": (2, 2, 25)}.items():
        d = build_points(r0, c0, W)
        out["W_points"][name] = nf_boot(d, ["R", "G", "MW"])
        r = out["W_points"][name]
        print(f"T9-W points {name} n={r['n']}: cy R-MW {r['cy']['diff']:+.3f} {np.round(r['cy']['ci'],3)}  cx {r['cx']['diff']:+.3f} {np.round(r['cx']['ci'],3)}", flush=True)
    # T9-S fixed sample
    ladder = [round(0.1 * k, 1) for k in range(11)] + [1.5, 2.0, 3.0]
    keysets = []
    for s in ladder:
        ks = set()
        for t in P.HOURS:
            h = T5.load_hour(t)
            if h is None:
                continue
            R, GIR, MW = h; Rb = T5.blur_nan(R, s)
            for i, j, sl in tiles(2, 2, 15):
                if ok({"R": Rb[sl], "IR": GIR[sl], "MW": MW[sl]}) is not None:
                    ks.add((t, i, j))
        keysets.append(ks)
    fixed = set.intersection(*keysets)
    out["S_fixed_n"] = len(fixed)
    lag = {"MW": []}; lag.update({s: [] for s in ladder[:11]})
    for t in P.HOURS[::4]:
        h = T5.load_hour(t)
        if h is None:
            continue
        R, GIR, MW = h; blurred = {s: T5.blur_nan(R, s) for s in ladder[:11]}
        for i, j, sl in tiles(2, 2, 15):
            if (t, i, j) not in fixed:
                continue
            lag["MW"].append(T5.lag1(np.nan_to_num(MW[sl])))
            for s in ladder[:11]:
                lag[s].append(T5.lag1(np.nan_to_num(blurred[s][sl])))
    med = {str(k): float(np.nanmedian(v)) for k, v in lag.items()}
    xs = ladder[:11]; ys = [med[str(s)] for s in xs]; target = med["MW"]
    star = None
    for a, b, ya, yb in zip(xs[:-1], xs[1:], ys[:-1], ys[1:]):
        if (ya - target) * (yb - target) <= 0 and ya != yb:
            star = a + (target - ya) * (b - a) / (yb - ya); break
    star = round(round(star / 0.05) * 0.05, 2) if star is not None else None
    out["S_lag1"] = med; out["S_sigma_star"] = star
    print("T9-S fixed n", len(fixed), "lag1 MW", round(target, 3), {s: round(med[str(s)], 3) for s in xs}, "sigma*", star, flush=True)
    out["S"] = {}
    for s in sorted(set(ladder + ([star] if star is not None else []))):
        d = build_area(2, 2, 15, blur=s, keyset=fixed)
        out["S"][str(s)] = nf_boot(d, ["R", "IR", "MW"])
        r = out["S"][str(s)]
        print(f"T9-S sigma={s} n={r['n']}: cy R-MW {r['cy']['diff']:+.3f} {np.round(r['cy']['ci'],3)}  cx {r['cx']['diff']:+.3f} {np.round(r['cx']['ci'],3)}", flush=True)
    json.dump(out, open("data/t9_results.json", "w"), indent=1, default=str)
