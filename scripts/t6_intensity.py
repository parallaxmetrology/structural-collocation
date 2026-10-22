"""T6 EXPLORATORY (not pre-registered): where does the shared satellite error live? Inversion nf_R - nf_MW by
quintile of radar window peak intensity (block-mean mm/h). Uses data/t6_windows.csv + radar peaks recomputed."""
import json, numpy as np, pandas as pd
import probe_a as P, t5 as T5, t6 as T6
d = pd.read_csv("data/t6_windows.csv", parse_dates=["t"])
peaks = {}
for t in P.HOURS:
    h = T5.load_hour(t)
    if h is None:
        continue
    R, GIR, MW = h
    for i, j, w in T5.windows(R, GIR, MW):
        peaks[(t, i, j)] = float(w["R"].max())
d["rpeak"] = [peaks[(t, i, j)] for t, i, j in zip(d.t, d.i, d.j)]
d["q"] = pd.qcut(d.rpeak, 5, labels=False)
out = {"edges_mm_h": [float(x) for x in d.groupby("q").rpeak.max()]}
for q in range(5):
    r = T5.nf_diff_boot(d[d.q == q], n_boot=300)
    out[str(q)] = {"n": int((d.q == q).sum()), **{v: {"diff": r[v]["diff_R_minus_MW"], "ci": r[v]["diff_ci95"],
                   "nf_R": r[v]["nf_R"], "nf_GIR": r[v]["nf_GIR"], "nf_MW": r[v]["nf_MW"]} for v in ["cy", "cx"]}}
    print(q, "peak<=", round(out["edges_mm_h"][q], 2), "n", out[str(q)]["n"],
          {v: (round(out[str(q)][v]["diff"], 3), [round(x, 3) for x in out[str(q)][v]["ci"]],
               round(out[str(q)][v]["nf_R"], 3), round(out[str(q)][v]["nf_MW"], 3)) for v in ["cy", "cx"]})
json.dump(out, open("data/t6_intensity.json", "w"), indent=1)
