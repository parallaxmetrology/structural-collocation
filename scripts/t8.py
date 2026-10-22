"""T8 (T8_PLAN.md): anvil fraction as a driver of the shared GOES-microwave error."""
import json
import numpy as np, pandas as pd
import probe_a as P
import t5 as T5
import t7 as T7

d = pd.read_csv("data/t7_points.csv", parse_dates=["t"])
anv = {}
for t, grp in d.groupby("t"):
    h = T5.load_hour(t)
    R, GIR, MW = h
    for i, j in zip(grp.i, grp.j):
        sl = (slice(P.R0 + P.WIN * i, P.R0 + P.WIN * (i + 1)), slice(P.C0 + P.WIN * j, P.C0 + P.WIN * (j + 1)))
        r, g = np.nan_to_num(R[sl]), np.nan_to_num(GIR[sl])
        cold = g > 0
        anv[(t, i, j)] = float((cold & (r < 0.5)).sum() / cold.sum()) if cold.sum() else np.nan
d["anvil"] = [anv[k] for k in zip(d.t, d.i, d.j)]
d = d.dropna(subset=["anvil"]).copy()
d["ter"] = pd.qcut(d.anvil.rank(method="first"), 3, labels=False)

rng = np.random.default_rng(9)
clus, uc = pd.factorize(d.clu); idx = [np.where(clus == c)[0] for c in range(len(uc))]
res = {"n": int(len(d)), "anvil_tercile_edges": [float(x) for x in d.groupby("ter").anvil.max()],
       "anvil_mean_by_tercile": [float(x) for x in d.groupby("ter").anvil.mean()],
       "corr_anvil_rpeak": float(np.corrcoef(d.anvil, d.rpeak)[0, 1]),
       "n_by_tercile": [int(x) for x in d.ter.value_counts().sort_index()]}
KEYS = ["c_IR_MW", "RGM_diff_R_MW", "RMG_nf_R", "RMG_nf_MW", "RMG_nf_G"]
for var in ["cy", "cx"]:
    X = d[[f"{k}_{var}" for k in ["R", "IR", "MW", "G"]]].values; ter = d.ter.values
    def at(rows):
        out = {}
        for q in range(3):
            rr = rows[ter[rows] == q]
            s = T7.stats(X[rr])
            for k in KEYS:
                out[f"{k}_{q}"] = s[k]
        return out
    pt = at(np.arange(len(d)))
    B = pd.DataFrame([at(np.concatenate([idx[k] for k in rng.integers(0, len(idx), len(idx))])) for _ in range(500)])
    r = {"point": {k: float(v) for k, v in pt.items()}}
    r["ci95"] = {k: [float(B[k].quantile(0.025)), float(B[k].quantile(0.975))] for k in B.columns}
    for k in ["c_IR_MW", "RGM_diff_R_MW"]:
        dd = B[f"{k}_2"] - B[f"{k}_0"]
        r[f"{k}_top_minus_bottom"] = float(pt[f"{k}_2"] - pt[f"{k}_0"])
        r[f"{k}_top_minus_bottom_ci"] = [float(dd.quantile(0.025)), float(dd.quantile(0.975))]
    res[var] = r
    print(var, "c(GIR,MW) by tercile", [round(pt[f"c_IR_MW_{q}"], 3) for q in range(3)],
          "top-bottom", round(r["c_IR_MW_top_minus_bottom"], 3), np.round(r["c_IR_MW_top_minus_bottom_ci"], 3),
          "| inversion by tercile", [round(pt[f"RGM_diff_R_MW_{q}"], 3) for q in range(3)],
          "top-bottom", round(r["RGM_diff_R_MW_top_minus_bottom"], 3), np.round(r["RGM_diff_R_MW_top_minus_bottom_ci"], 3),
          "| RMG nf R/MW/G", [(round(pt[f"RMG_nf_R_{q}"], 3), round(pt[f"RMG_nf_MW_{q}"], 3), round(pt[f"RMG_nf_G_{q}"], 3)) for q in range(3)])
print("n", res["n"], "anvil tercile means", np.round(res["anvil_mean_by_tercile"], 3), "corr(anvil, rpeak)", round(res["corr_anvil_rpeak"], 3))
json.dump(res, open("data/t8_results.json", "w"), indent=1)
