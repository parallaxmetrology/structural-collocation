"""T6 (T6_PLAN.md): core-vs-full-area test of the shared satellite error."""
import json
import numpy as np, pandas as pd
import probe_a as P
import t5 as T5

VARS = ["cy", "cx", "disp"]


def core(A):
    pos = A[A > 0]
    if len(pos) < 5:
        return A
    return A * (A >= np.quantile(pos, 0.8))


def build():
    rows = []
    for t in P.HOURS:
        h = T5.load_hour(t)
        if h is None:
            continue
        R, GIR, MW = h
        for i, j, w in T5.windows(R, GIR, MW):
            row = {"t": t, "i": i, "j": j}
            for k, v in w.items():
                full = P.reads(v); cr = P.reads(core(v))
                for var in VARS:
                    row[f"{k}_{var}"] = full[var]; row[f"{k}c_{var}"] = cr[var]
            r = w["R"]
            row["light_frac"] = float(r[r < 5].sum() / r.sum())
            rows.append(row)
    d = pd.DataFrame(rows)
    lat = 55.0 - 0.1 * (P.R0 + P.WIN * d.i + P.WIN / 2); lon = -130.0 + 0.1 * (P.C0 + P.WIN * d.j + P.WIN / 2)
    d["clu"] = d.t.dt.floor("D").astype(str) + "|" + np.floor(lat / 5).astype(int).astype(str) + "_" + np.floor(lon / 5).astype(int).astype(str)
    return d


def core_frame(d):
    c = d[["t", "i", "j", "clu"]].copy()
    for k in ["R", "IR", "MW"]:
        for v in VARS:
            c[f"{k}_{v}"] = d[f"{k}c_{v}"]
    return c


def tercile_diff(d, n_boot=500, seed=6):
    """joint cluster bootstrap of (top - bottom tercile) of nf_R - nf_MW."""
    rng = np.random.default_rng(seed)
    d = d.copy(); d["ter"] = pd.qcut(d.light_frac.rank(method="first"), 3, labels=False)
    clus, uc = pd.factorize(d.clu); nclu = len(uc)
    out = {"tercile_edges": [float(x) for x in d.groupby("ter").light_frac.max()]}
    for v in VARS:
        stats = {}
        for ter in range(3):
            m = (d.ter == ter).values
            X = d.loc[m, [f"R_{v}", f"IR_{v}", f"MW_{v}"]].values; cl = clus[m]
            n_c = np.bincount(cl, minlength=nclu).astype(float)
            s_c = np.stack([np.bincount(cl, X[:, k], minlength=nclu) for k in range(3)], 1)
            ss = np.einsum("ri,rj->rij", X, X)
            ss_c = np.stack([np.stack([np.bincount(cl, ss[:, a, b], minlength=nclu) for b in range(3)], 1) for a in range(3)], 1)
            stats[ter] = (n_c, s_c, ss_c)
        def diff(w, ter):
            n_c, s_c, ss_c = stats[ter]
            n = (w * n_c).sum(); s = (w[:, None] * s_c).sum(0); q = (w[:, None, None] * ss_c).sum(0)
            C = q / n - np.outer(s / n, s / n)
            nfR = 1 - C[0, 1] * C[0, 2] / (C[0, 0] * C[1, 2]); nfM = 1 - C[2, 0] * C[2, 1] / (C[2, 2] * C[0, 1])
            return nfR - nfM
        ones = np.ones(nclu)
        pt = {ter: float(diff(ones, ter)) for ter in range(3)}
        B = []
        for _ in range(n_boot):
            w = rng.multinomial(nclu, ones / nclu).astype(float)
            B.append([diff(w, ter) for ter in range(3)])
        B = np.array(B); td = B[:, 2] - B[:, 0]
        out[v] = {"by_tercile": pt, "by_tercile_ci": {ter: [float(np.percentile(B[:, ter], 2.5)), float(np.percentile(B[:, ter], 97.5))] for ter in range(3)},
                  "top_minus_bottom": pt[2] - pt[0], "top_minus_bottom_ci": [float(np.percentile(td, 2.5)), float(np.percentile(td, 97.5))],
                  "n_per_tercile": [int((d.ter == ter).sum()) for ter in range(3)]}
    return out


def geometry(d, n_boot=1000, seed=8):
    rng = np.random.default_rng(seed)
    gap = np.hypot(d.R_cy - d.Rc_cy, d.R_cx - d.Rc_cx)
    s = d[gap >= 1].copy()
    out = {"n_windows": int(len(s)), "share_of_all": float(len(s) / len(d))}
    clus, uc = pd.factorize(s.clu); idx = [np.where(clus == c)[0] for c in range(len(uc))]
    for k, nm in [("IR", "GOES_IR"), ("MW", "MW")]:
        dfull = np.hypot(s[f"{k}_cy"] - s.R_cy, s[f"{k}_cx"] - s.R_cx)
        dcore = np.hypot(s[f"{k}_cy"] - s.Rc_cy, s[f"{k}_cx"] - s.Rc_cx)
        diff = (dfull - dcore).values
        bs = [np.median(diff[np.concatenate([idx[i] for i in rng.integers(0, len(idx), len(idx))])]) for _ in range(n_boot)]
        out[nm] = {"median_dist_to_R_full": float(np.median(dfull)), "median_dist_to_R_core": float(np.median(dcore)),
                   "median_paired_diff": float(np.median(diff)), "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                   "share_closer_to_core": float((diff > 0).mean())}
    return out


if __name__ == "__main__":
    d = build(); d.to_csv("data/t6_windows.csv", index=False)
    res = {"n": int(len(d))}
    res["full"] = T5.nf_diff_boot(d)
    res["core"] = T5.nf_diff_boot(core_frame(d))
    res["strata"] = tercile_diff(d)
    res["geometry"] = geometry(d)
    for v in VARS:
        print(f"{v}: full R-MW {res['full'][v]['diff_R_minus_MW']:+.3f} {res['full'][v]['diff_ci95']}  "
              f"core R-MW {res['core'][v]['diff_R_minus_MW']:+.3f} {res['core'][v]['diff_ci95']}  "
              f"terciles {[round(x, 3) for x in res['strata'][v]['by_tercile'].values()]} top-bottom "
              f"{res['strata'][v]['top_minus_bottom']:+.3f} {res['strata'][v]['top_minus_bottom_ci']}")
    print("core nf:", {v: (round(res['core'][v]['nf_R'], 3), round(res['core'][v]['nf_GIR'], 3), round(res['core'][v]['nf_MW'], 3)) for v in VARS})
    print("geometry:", json.dumps(res["geometry"], indent=1))
    print("tercile edges (light fraction):", res["strata"]["tercile_edges"])
    json.dump(res, open("data/t6_results.json", "w"), indent=1, default=str)
