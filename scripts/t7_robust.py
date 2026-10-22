"""T7-R (T7R_PLAN.md): mirror triplet and four-observer (tetrad) consistency on the T7 sample."""
import json
import numpy as np, pandas as pd

OBS = ["R", "IR", "MW", "G"]
TRIPLETS = [("R", "IR", "MW"), ("R", "IR", "G"), ("R", "MW", "G"), ("IR", "MW", "G")]


def nf(C, a, b, c):
    i = {o: k for k, o in enumerate(OBS)}
    return 1 - C[i[a], i[b]] * C[i[a], i[c]] / (C[i[a], i[a]] * C[i[b], i[c]])


def stats(X):
    C = np.cov(X.T); r = np.corrcoef(X.T); i = {o: k for k, o in enumerate(OBS)}
    rr = lambda a, b: r[i[a], i[b]]
    out = {}
    for tri in TRIPLETS:
        for o in tri:
            others = [x for x in tri if x != o]
            out[f"nf_{o}|{''.join(tri)}"] = nf(C, o, *others)
    out["R1_R_minus_IR"] = out["nf_R|RIRG"] - out["nf_IR|RIRG"]
    out["R1_R_minus_G"] = out["nf_R|RIRG"] - out["nf_G|RIRG"]
    out["R2_nfR_RIRG_minus_RMWG"] = out["nf_R|RIRG"] - out["nf_R|RMWG"]
    out["t1"] = rr("R", "G") * rr("IR", "MW") - rr("R", "IR") * rr("G", "MW")
    out["t2"] = rr("R", "G") * rr("IR", "MW") - rr("R", "MW") * rr("G", "IR")
    out["t3"] = rr("R", "IR") * rr("G", "MW") - rr("R", "MW") * rr("G", "IR")
    return out


def analyse(s, var, n_boot=500, seed=11):
    rng = np.random.default_rng(seed)
    X = s[[f"{o}_{var}" for o in OBS]].values
    clus, uc = pd.factorize(s.clu); idx = [np.where(clus == c)[0] for c in range(len(uc))]
    pt = stats(X)
    B = pd.DataFrame([stats(X[np.concatenate([idx[k] for k in rng.integers(0, len(idx), len(idx))])]) for _ in range(n_boot)])
    return {"n": int(len(s)), "point": {k: float(v) for k, v in pt.items()},
            "ci95": {k: [float(B[k].quantile(0.025)), float(B[k].quantile(0.975))] for k in B.columns}}


if __name__ == "__main__":
    d = pd.read_csv("data/t7_points.csv", parse_dates=["t"])
    med = float(d.rpeak.median())
    res = {"median_peak": med}
    for stratum, s in [("weak", d[d.rpeak <= med]), ("strong", d[d.rpeak > med])]:
        res[stratum] = {v: analyse(s, v) for v in ["cy", "cx"]}
        for v in ["cy", "cx"]:
            p, c = res[stratum][v]["point"], res[stratum][v]["ci95"]
            f = lambda k: f"{p[k]:+.3f} [{c[k][0]:+.3f},{c[k][1]:+.3f}]"
            print(f"{stratum} {v}: R1 R-GIR {f('R1_R_minus_IR')} R-G {f('R1_R_minus_G')} | R2 {f('R2_nfR_RIRG_minus_RMWG')} | "
                  f"t1 {f('t1')} t2 {f('t2')} t3 {f('t3')}")
            print("   nf table:", {k: round(x, 3) for k, x in p.items() if k.startswith("nf_")})
    json.dump(res, open("data/t7r_results.json", "w"), indent=1)
