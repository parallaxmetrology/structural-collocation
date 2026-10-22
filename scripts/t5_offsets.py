"""T5 EXPLORATORY follow-up (not pre-registered): are the satellite observers displaced from radar in a
shared direction? Uses Probe A2 window reads (data/probe_a2_windows.csv). Offsets in 0.1-degree cells;
centroid y is measured downward (row index), so a negative dy means the satellite centroid lies NORTH of radar."""
import json, numpy as np, pandas as pd
d = pd.read_csv("data/probe_a2_windows.csv", parse_dates=["t"])
out = {}
for ax in ["cy", "cx"]:
    g = d[f"IR_{ax}"] - d[f"R_{ax}"]; m = d[f"MW_{ax}"] - d[f"R_{ax}"]
    out[ax] = {"mean_GIR_minus_R": float(g.mean()), "mean_MW_minus_R": float(m.mean()),
               "se_GIR": float(g.std() / np.sqrt(len(g))), "se_MW": float(m.std() / np.sqrt(len(m))),
               "corr_offsets": float(np.corrcoef(g, m)[0, 1])}
lat = 55.0 - 0.1 * (2 + 15 * d.i + 7.5)
d["latband"] = pd.cut(lat, [20, 30, 35, 40, 45, 55])
by = d.assign(g=d.IR_cy - d.R_cy, m=d.MW_cy - d.R_cy).groupby("latband", observed=True).agg(
    n=("g", "size"), GIR_minus_R=("g", "mean"), MW_minus_R=("m", "mean"))
out["by_latitude_cy"] = {str(k): {kk: float(vv) for kk, vv in v.items()} for k, v in by.to_dict("index").items()}
json.dump(out, open("data/t5_offsets.json", "w"), indent=1)
print(json.dumps(out, indent=1))
