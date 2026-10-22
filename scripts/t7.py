"""T7 (T7_PLAN.md): gauges as a fourth, sharp, independent observer at point support."""
import json
import numpy as np, pandas as pd
import probe_a as P
import t5 as T5

SEED = 7
st = pd.read_csv("data/stations.csv")
st["key"] = st.src + ":" + st.sid.astype(str)
st["br"] = np.clip(np.floor((55.0 - st.lat) / 0.1).astype(int), 0, 349)
st["bc"] = np.clip(np.floor((st.lon + 130.0) / 0.1).astype(int), 0, 699)
g = pd.read_csv("data/gauges.csv", parse_dates=["t_end"])
g["si"] = (g.src + ":" + g.sid.astype(str)).map({k: i for i, k in enumerate(st.key)})
G_BY_T = {t: grp for t, grp in g.groupby("t_end")}


def build():
    rows = []
    for t in P.HOURS:
        h = T5.load_hour(t)
        if h is None or t not in G_BY_T:
            continue
        R, GIR, MW = h
        gt = G_BY_T[t]
        gr, gc, gv = st.br.values[gt.si.values], st.bc.values[gt.si.values], gt.p_mm.values
        for i, j, w in T5.windows(R, GIR, MW):
            r0, c0 = P.R0 + P.WIN * i, P.C0 + P.WIN * j
            m = (gr >= r0) & (gr < r0 + P.WIN) & (gc >= c0) & (gc < c0 + P.WIN)
            if m.sum() < 5:
                continue
            lr, lc = gr[m] - r0, gc[m] - c0
            y = lr - (P.WIN - 1) / 2; x = lc - (P.WIN - 1) / 2
            vals = {"R": w["R"][lr, lc], "IR": w["IR"][lr, lc], "MW": w["MW"][lr, lc], "G": gv[m]}
            if min(v.sum() for v in vals.values()) <= 0:
                continue
            row = {"t": t, "i": i, "j": j, "n_g": int(m.sum()), "rpeak": float(w["R"].max())}
            for k, v in vals.items():
                row[f"{k}_cy"] = float((v * y).sum() / v.sum()); row[f"{k}_cx"] = float((v * x).sum() / v.sum())
            rows.append(row)
    d = pd.DataFrame(rows)
    lat = 55.0 - 0.1 * (P.R0 + P.WIN * d.i + P.WIN / 2); lon = -130.0 + 0.1 * (P.C0 + P.WIN * d.j + P.WIN / 2)
    d["clu"] = d.t.dt.floor("D").astype(str) + "|" + np.floor(lat / 5).astype(int).astype(str) + "_" + np.floor(lon / 5).astype(int).astype(str)
    return d


def tc_nf(C, a, b, c):
    """noise fraction of a in triplet (a, b, c); and its error variance."""
    return 1 - C[a, b] * C[a, c] / (C[a, a] * C[b, c]), C[a, a] - C[a, b] * C[a, c] / C[b, c]


def quad(C, p, q, u, v):
    c = C[p, q] - C[p, u] * C[q, v] / C[u, v]
    ep = C[p, p] - C[p, u] * C[p, v] / C[u, v]; eq = C[q, q] - C[q, u] * C[q, v] / C[u, v]
    return c / np.sqrt(ep * eq) if ep > 0 and eq > 0 else np.nan


IDX = {"R": 0, "IR": 1, "MW": 2, "G": 3}


def stats(X):
    C = np.cov(X.T)
    out = {}
    for name, (a, b, c) in {"RGM": ("R", "IR", "MW"), "RMG": ("R", "MW", "G")}.items():
        ia, ib, ic = IDX[a], IDX[b], IDX[c]
        nfa, eva = tc_nf(C, ia, ib, ic); nfb, evb = tc_nf(C, ib, ia, ic); nfc, evc = tc_nf(C, ic, ia, ib)
        out[f"{name}_nf_{a}"] = nfa; out[f"{name}_nf_{b}"] = nfb; out[f"{name}_nf_{c}"] = nfc
        out[f"{name}_neg"] = float(min(eva, evb, evc) < 0)
        out[f"{name}_diff_R_MW"] = out[f"{name}_nf_R"] - out[f"{name}_nf_MW"]
    out["c_IR_MW"] = quad(C, IDX["IR"], IDX["MW"], IDX["R"], IDX["G"])
    out["c_R_MW"] = quad(C, IDX["R"], IDX["MW"], IDX["IR"], IDX["G"])
    out["r"] = {f"{a}-{b}": float(np.corrcoef(X[:, IDX[a]], X[:, IDX[b]])[0, 1])
                for a, b in [("R", "IR"), ("R", "MW"), ("IR", "MW"), ("R", "G"), ("IR", "G"), ("MW", "G")]}
    return out


def analyse(d, var, n_boot=500):
    rng = np.random.default_rng(SEED)
    X = d[[f"{k}_{var}" for k in ["R", "IR", "MW", "G"]]].values
    clus, uc = pd.factorize(d.clu); idx = [np.where(clus == c)[0] for c in range(len(uc))]
    pt = stats(X)
    B = []
    for _ in range(n_boot):
        rows = np.concatenate([idx[k] for k in rng.integers(0, len(idx), len(idx))])
        s = stats(X[rows]); s.pop("r"); B.append(s)
    B = pd.DataFrame(B)
    ci = {k: [float(B[k].quantile(0.025)), float(B[k].quantile(0.975))] for k in B.columns if not k.endswith("_neg")}
    neg = {k: float(B[k].mean()) for k in B.columns if k.endswith("_neg")}
    return {"n": int(len(d)), "clusters": int(len(uc)), "point": {k: (float(v) if not isinstance(v, dict) else v) for k, v in pt.items()},
            "ci95": ci, "share_negative": neg}


if __name__ == "__main__":
    d = build(); d.to_csv("data/t7_points.csv", index=False)
    med = float(d.rpeak.median())
    res = {"n_total": int(len(d)), "median_peak": med, "gauges_per_window_median": float(d.n_g.median())}
    for stratum, s in [("weak", d[d.rpeak <= med]), ("strong", d[d.rpeak > med])]:
        res[stratum] = {v: analyse(s, v) for v in ["cy", "cx"]}
        for v in ["cy", "cx"]:
            a = res[stratum][v]; p = a["point"]; c = a["ci95"]
            print(f"{stratum} {v} n={a['n']}: P0 RGM R-MW {p['RGM_diff_R_MW']:+.3f} {np.round(c['RGM_diff_R_MW'],3)} "
                  f"(nf R {p['RGM_nf_R']:.3f} GIR {p['RGM_nf_IR']:.3f} MW {p['RGM_nf_MW']:.3f}) neg {a['share_negative']['RGM_neg']:.2f} | "
                  f"P1 RMG R-MW {p['RMG_diff_R_MW']:+.3f} {np.round(c['RMG_diff_R_MW'],3)} (nf R {p['RMG_nf_R']:.3f} MW {p['RMG_nf_MW']:.3f} G {p['RMG_nf_G']:.3f}) neg {a['share_negative']['RMG_neg']:.2f} | "
                  f"P2 c(GIR,MW) {p['c_IR_MW']:+.3f} {np.round(c['c_IR_MW'],3)} | c(R,MW) {p['c_R_MW']:+.3f} {np.round(c['c_R_MW'],3)} | r " + str({k: round(x, 2) for k, x in p["r"].items()}))
    print("total", res["n_total"], "median peak", round(med, 2), "median gauges/window", res["gauges_per_window_median"])
    json.dump(res, open("data/t7_results.json", "w"), indent=1)
