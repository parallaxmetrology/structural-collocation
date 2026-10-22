"""POST-HOC (not in PROBE_A_PLAN.md): test whether the MW noise fraction is biased low by a
shared hour-integrated footprint in R and IR. Hypothesis: windows where MW covers BOTH half-hours
(an hour-integrated MW) give a higher MW noise fraction than windows where MW covers one half-hour.
Also: TC with the alternative triplets (R, IR, M), (R, MW, M) for comparison."""
import json, numpy as np, pandas as pd
import probe_a as P

rows = []
for t in P.HOURS:
    t0 = t - pd.Timedelta(hours=1)
    a, b = np.load(f"data/imerg_x/{t0:%Y%m%d%H}00.npz"), np.load(f"data/imerg_x/{t0:%Y%m%d%H}30.npz")
    fa, fb = np.isfinite(a["MW"]), np.isfinite(b["MW"])
    for i in range(P.NR):
        for j in range(P.NC):
            sl = (slice(P.R0 + P.WIN * i, P.R0 + P.WIN * (i + 1)), slice(P.C0 + P.WIN * j, P.C0 + P.WIN * (j + 1)))
            rows.append((t, i, j, fa[sl].mean() >= 0.8 and fb[sl].mean() >= 0.8))
cov = pd.DataFrame(rows, columns=["t", "i", "j", "both"])
d = pd.read_csv("data/probe_a_windows.csv", parse_dates=["t"]).merge(cov, on=["t", "i", "j"])
out = {"n_both": int(d.both.sum()), "n_one": int((~d.both).sum())}
for name, sub in (("both_halves", d[d.both]), ("one_half", d[~d.both])):
    r = P.analyse(sub, n_boot=300)
    out[name] = {v: {o: round(r[v]["point"][f"nf_{o}"], 3) for o in ("R", "IR", "MW")} |
                 {"nf_MW_ci": [round(x, 3) for x in r[v]["ci95"]["nf_MW"]]} for v in P.VARS}
# alternative triplets on the full sample
def tri(C, x, y, z): return 1 - C[x, y] * C[x, z] / (C[x, x] * C[y, z])
alt = {}
for v in P.VARS:
    C = np.cov(d[[f"{o}_{v}" for o in P.OBS]].values.T)
    alt[v] = {"R|(IR,M)": round(tri(C, 0, 1, 3), 3), "R|(MW,M)": round(tri(C, 0, 2, 3), 3),
              "IR|(R,M)": round(tri(C, 1, 0, 3), 3), "MW|(R,M)": round(tri(C, 2, 0, 3), 3)}
out["alt_triplets_nf"] = alt
json.dump(out, open("data/probe_a_posthoc.json", "w"), indent=1)
print(json.dumps(out, indent=1))
