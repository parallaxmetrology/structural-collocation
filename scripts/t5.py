"""T5 (T5_PLAN.md): coarsen radar to the microwave field's smoothness; does microwave < radar survive?"""
import json, os
import numpy as np, pandas as pd
from scipy import ndimage
import probe_a as P
import t1_synthetic as T1

VARS = ["cy", "cx", "disp"]
LADDER = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]
FINE = [x * 0.25 for x in range(13)]


def nf_diff_boot(d, n_boot=500, seed=5):
    """noise fractions (R, IR, MW) and the bootstrap CI of nf_R - nf_MW, per variable."""
    rng = np.random.default_rng(seed)
    clus, uc = pd.factorize(d.clu); nclu = len(uc)
    out = {}
    for v in VARS:
        X = d[[f"R_{v}", f"IR_{v}", f"MW_{v}"]].values
        n_c = np.bincount(clus, minlength=nclu).astype(float)
        s_c = np.stack([np.bincount(clus, X[:, k], minlength=nclu) for k in range(3)], 1)
        ss = np.einsum("ri,rj->rij", X, X)
        ss_c = np.stack([np.stack([np.bincount(clus, ss[:, a, b], minlength=nclu) for b in range(3)], 1) for a in range(3)], 1)
        def nfs(w):
            n = (w * n_c).sum(); s = (w[:, None] * s_c).sum(0); q = (w[:, None, None] * ss_c).sum(0)
            C = q / n - np.outer(s / n, s / n)
            return np.array([1 - C[0, 1] * C[0, 2] / (C[0, 0] * C[1, 2]),
                             1 - C[1, 0] * C[1, 2] / (C[1, 1] * C[0, 2]),
                             1 - C[2, 0] * C[2, 1] / (C[2, 2] * C[0, 1])])
        pt = nfs(np.ones(nclu))
        B = np.array([nfs(rng.multinomial(nclu, np.ones(nclu) / nclu).astype(float)) for _ in range(n_boot)])
        dd = B[:, 0] - B[:, 2]
        out[v] = {"nf_R": float(pt[0]), "nf_GIR": float(pt[1]), "nf_MW": float(pt[2]),
                  "nf_ci": {k: [float(np.percentile(B[:, i], 2.5)), float(np.percentile(B[:, i], 97.5))]
                            for i, k in enumerate(["R", "GIR", "MW"])},
                  "diff_R_minus_MW": float(pt[0] - pt[2]),
                  "diff_ci95": [float(np.percentile(dd, 2.5)), float(np.percentile(dd, 97.5))]}
    return out


def blur_nan(A, s):
    if s == 0:
        return A
    v = np.isfinite(A).astype(float)
    num = ndimage.gaussian_filter(np.nan_to_num(A), s, mode="constant")
    den = ndimage.gaussian_filter(v, s, mode="constant")
    with np.errstate(invalid="ignore", divide="ignore"):
        out = num / den
    out[v == 0] = np.nan
    return out


def lag1(A):
    a = [np.corrcoef(A[:, :-1].ravel(), A[:, 1:].ravel())[0, 1], np.corrcoef(A[:-1].ravel(), A[1:].ravel())[0, 1]]
    return np.nanmean(a)


def load_hour(t):
    if not os.path.exists(f"data/goes/{t:%Y%m%d%H}.npz"):
        return None
    R = np.load(f"data/mrms/{t:%Y%m%d%H}.npz")["R"]
    Tb = np.load(f"data/goes/{t:%Y%m%d%H}.npz")["Tb"]
    GIR = np.where(np.isfinite(Tb), np.maximum(0.0, 235.0 - Tb), np.nan)
    t0 = t - pd.Timedelta(hours=1)
    a, b = np.load(f"data/imerg_x/{t0:%Y%m%d%H}00.npz"), np.load(f"data/imerg_x/{t0:%Y%m%d%H}30.npz")
    with np.errstate(all="ignore"):
        MW = np.nanmean(np.stack([a["MW"], b["MW"]]), axis=0)
    return R, GIR, MW


def windows(R, GIR, MW):
    for i in range(P.NR):
        for j in range(P.NC):
            sl = (slice(P.R0 + P.WIN * i, P.R0 + P.WIN * (i + 1)), slice(P.C0 + P.WIN * j, P.C0 + P.WIN * (j + 1)))
            w = {"R": R[sl], "IR": GIR[sl], "MW": MW[sl]}
            if np.isfinite(w["R"]).mean() < 0.9 or not np.isfinite(w["IR"]).all() or np.isfinite(w["MW"]).mean() < 0.8:
                continue
            w = {k: np.nan_to_num(v) for k, v in w.items()}
            if min(v.sum() for v in w.values()) < 1:
                continue
            yield i, j, w


def smoothness(hours):
    """median lag-1 autocorrelation: MW, and radar at each fine blur level (field properties only)."""
    acc = {"MW": [], **{s: [] for s in FINE}}
    for t in hours:
        h = load_hour(t)
        if h is None:
            continue
        R, GIR, MW = h
        blurred = {s: blur_nan(R, s) for s in FINE}
        for i, j, w in windows(R, GIR, MW):
            acc["MW"].append(lag1(w["MW"]))
            sl = (slice(P.R0 + P.WIN * i, P.R0 + P.WIN * (i + 1)), slice(P.C0 + P.WIN * j, P.C0 + P.WIN * (j + 1)))
            for s in FINE:
                acc[s].append(lag1(np.nan_to_num(blurred[s][sl])))
    med = {k: float(np.nanmedian(v)) for k, v in acc.items()}
    star = min(FINE, key=lambda s: abs(med[s] - med["MW"]))
    return med, star


def build_real(sig):
    rows = []
    for t in P.HOURS:
        h = load_hour(t)
        if h is None:
            continue
        R, GIR, MW = h
        Rb = blur_nan(R, sig)
        for i, j, w in windows(Rb, GIR, MW):
            row = {"t": t, "i": i, "j": j}
            for k, v in w.items():
                rd = P.reads(v)
                for var in VARS:
                    row[f"{k}_{var}"] = rd[var]
            rows.append(row)
    d = pd.DataFrame(rows)
    lat = 55.0 - 0.1 * (P.R0 + P.WIN * d.i + P.WIN / 2); lon = -130.0 + 0.1 * (P.C0 + P.WIN * d.j + P.WIN / 2)
    d["clu"] = d.t.dt.floor("D").astype(str) + "|" + np.floor(lat / 5).astype(int).astype(str) + "_" + np.floor(lon / 5).astype(int).astype(str)
    return d


def synthetic_arm():
    wins = T1.simulate()
    out = {}
    base = T1.OBS["R"]
    for sig_sh in [0.0, 1.0]:
        for add in [0.0, 1.0, 2.0]:
            T1.OBS["R"] = (base[0], add, base[2])
            d = T1.build(wins, sig_sh)
            r = nf_diff_boot(d, n_boot=200)
            tr = T1.truth_scores(d)
            out[f"sh{sig_sh}_Rblur{add}"] = {v: {**r[v], "true_R": tr[v]["nf_true"]["R"], "true_MW": tr[v]["nf_true"]["MW"]} for v in VARS}
            print(f"synthetic sigma_sh={sig_sh} R blur +{add}: " + "  ".join(
                f"{v}: R-S1 {r[v]['diff_R_minus_MW']:+.3f} [{r[v]['diff_ci95'][0]:+.3f},{r[v]['diff_ci95'][1]:+.3f}] (true {tr[v]['nf_true']['R'] - tr[v]['nf_true']['MW']:+.3f})"
                for v in VARS), flush=True)
    T1.OBS["R"] = base
    return out


if __name__ == "__main__":
    res = {"synthetic": synthetic_arm()}
    med, star = smoothness(P.HOURS[::4])
    res["smoothness_median_lag1"] = {str(k): v for k, v in med.items()}; res["sigma_star"] = star
    print("lag-1 MW", round(med["MW"], 3), {s: round(med[s], 3) for s in FINE}, "sigma*", star, flush=True)
    levels = sorted(set(LADDER + [star]))
    res["real"] = {}
    for s in levels:
        d = build_real(s)
        r = nf_diff_boot(d)
        res["real"][str(s)] = {"n": int(len(d)), **r}
        print(f"real sigma={s} n={len(d)}: " + "  ".join(
            f"{v}: R {r[v]['nf_R']:.3f} GIR {r[v]['nf_GIR']:.3f} MW {r[v]['nf_MW']:.3f} R-MW {r[v]['diff_R_minus_MW']:+.3f} "
            f"[{r[v]['diff_ci95'][0]:+.3f},{r[v]['diff_ci95'][1]:+.3f}]" for v in VARS), flush=True)
    json.dump(res, open("data/t5_results.json", "w"), indent=1)
