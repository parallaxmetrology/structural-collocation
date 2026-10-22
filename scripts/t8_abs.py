"""T8 EXPLORATORY (not pre-registered): absolute (unnormalised) shared GOES-MW error covariance and each observer's
error variance by anvil tercile, instruments radar and gauges. Units: centroid cells^2."""
import json, numpy as np, pandas as pd
import t5 as T5, probe_a as P
d = pd.read_csv("data/t7_points.csv", parse_dates=["t"])
r8 = json.load(open("data/t8_results.json"))
anv = {}
for t, grp in d.groupby("t"):
    R, GIR, MW = T5.load_hour(t)
    for i, j in zip(grp.i, grp.j):
        sl = (slice(P.R0 + P.WIN * i, P.R0 + P.WIN * (i + 1)), slice(P.C0 + P.WIN * j, P.C0 + P.WIN * (j + 1)))
        r, g = np.nan_to_num(R[sl]), np.nan_to_num(GIR[sl]); cold = g > 0
        anv[(t, i, j)] = float((cold & (r < 0.5)).sum() / cold.sum()) if cold.sum() else np.nan
d["anvil"] = [anv[k] for k in zip(d.t, d.i, d.j)]; d = d.dropna(subset=["anvil"])
d["ter"] = pd.qcut(d.anvil.rank(method="first"), 3, labels=False)
out = {}
for var in ["cy", "cx"]:
    out[var] = {}
    for q in range(3):
        X = d[d.ter == q][[f"{k}_{var}" for k in ["R", "IR", "MW", "G"]]].values; C = np.cov(X.T)
        R_, I_, M_, G_ = 0, 1, 2, 3
        cov_IM = C[I_, M_] - C[I_, R_] * C[M_, G_] / C[R_, G_]
        ev = {"GIR": C[I_, I_] - C[I_, R_] * C[I_, G_] / C[R_, G_], "MW": C[M_, M_] - C[M_, R_] * C[M_, G_] / C[R_, G_],
              "R": C[R_, R_] - C[R_, M_] * C[R_, G_] / C[M_, G_], "G": C[G_, G_] - C[G_, R_] * C[G_, M_] / C[R_, M_]}
        out[var][q] = {"shared_cov_GIR_MW": float(cov_IM), **{f"errvar_{k}": float(v) for k, v in ev.items()},
                       "signal_var_R": float(C[R_, M_] * C[R_, G_] / C[M_, G_])}
        print(var, "tercile", q, {k: round(v, 3) for k, v in out[var][q].items()})
json.dump(out, open("data/t8_abs.json", "w"), indent=1)
