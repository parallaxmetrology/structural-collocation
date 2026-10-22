"""T6 EXPLORATORY (not pre-registered): is the inversion radar's genuine beam-quality error?
Inversion nf_R - nf_MW by tercile of window-mean MRMS RQI, and by RQI x intensity."""
import json, numpy as np, pandas as pd
import probe_a as P, t5 as T5
d = pd.read_csv("data/t6_windows.csv", parse_dates=["t"])
rq, pk = {}, {}
for t in P.HOURS:
    h = T5.load_hour(t)
    if h is None:
        continue
    R, GIR, MW = h
    Q = np.load(f"data/mrms/{t:%Y%m%d%H}.npz")["Q"]
    for i, j, w in T5.windows(R, GIR, MW):
        sl = (slice(P.R0 + P.WIN * i, P.R0 + P.WIN * (i + 1)), slice(P.C0 + P.WIN * j, P.C0 + P.WIN * (j + 1)))
        rq[(t, i, j)] = float(np.nanmean(Q[sl])); pk[(t, i, j)] = float(w["R"].max())
d["rqi"] = [rq[k] for k in zip(d.t, d.i, d.j)]; d["rpeak"] = [pk[k] for k in zip(d.t, d.i, d.j)]
d["rq3"] = pd.qcut(d.rqi.rank(method="first"), 3, labels=False)
d["pk2"] = (d.rpeak > d.rpeak.median()).astype(int)
out = {"rqi_edges": [float(x) for x in d.groupby("rq3").rqi.max()], "corr_rqi_peak": float(np.corrcoef(d.rqi, d.rpeak)[0, 1])}
for a in range(3):
    for b in range(2):
        s = d[(d.rq3 == a) & (d.pk2 == b)]
        r = T5.nf_diff_boot(s, n_boot=300)
        out[f"rqi{a}_peak{b}"] = {"n": int(len(s)), "rqi_mean": float(s.rqi.mean()),
                                  **{v: {"diff": r[v]["diff_R_minus_MW"], "ci": r[v]["diff_ci95"], "nf_R": r[v]["nf_R"], "nf_MW": r[v]["nf_MW"]} for v in ["cy", "cx"]}}
        x = out[f"rqi{a}_peak{b}"]
        print(f"RQI tercile {a} (mean {x['rqi_mean']:.2f}) {'strong' if b else 'weak'} n={x['n']}: "
              f"cy R-MW {x['cy']['diff']:+.3f} [{x['cy']['ci'][0]:+.3f},{x['cy']['ci'][1]:+.3f}] nf_R {x['cy']['nf_R']:.3f} nf_MW {x['cy']['nf_MW']:.3f} | "
              f"cx {x['cx']['diff']:+.3f}")
print("RQI edges", out["rqi_edges"], "corr(RQI, peak)", round(out["corr_rqi_peak"], 3))
json.dump(out, open("data/t6_rqi.json", "w"), indent=1)
