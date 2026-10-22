"""Probe A2 (PROBE_A2_PLAN.md): Probe A with GOES-16 band 13 cold-cloud field (GIRc = max(0, 235 - Tb))
replacing IMERG IRprecipitation. Analysis code is reused unchanged from probe_a.py; the GIR observer
occupies the 'IR' slot so the estimators and checks are identical."""
import json
import numpy as np, pandas as pd
import probe_a as P


def build():
    rows = []
    for t in P.HOURS:
        import os
        if not os.path.exists(f"data/goes/{t:%Y%m%d%H}.npz"):  # <3 GOES scans that hour (plan)
            continue
        R = np.load(f"data/mrms/{t:%Y%m%d%H}.npz")["R"]
        g = np.load(f"data/goes/{t:%Y%m%d%H}.npz")
        GIR = np.maximum(0.0, 235.0 - g["Tb"]); GIR[~np.isfinite(g["Tb"])] = np.nan
        t0 = t - pd.Timedelta(hours=1)
        M = (np.load(f"data/imerg/{t0:%Y%m%d%H}00.npy") + np.load(f"data/imerg/{t0:%Y%m%d%H}30.npy")) / 2
        a, b = np.load(f"data/imerg_x/{t0:%Y%m%d%H}00.npz"), np.load(f"data/imerg_x/{t0:%Y%m%d%H}30.npz")
        with np.errstate(all="ignore"):
            MW = np.nanmean(np.stack([a["MW"], b["MW"]]), axis=0)
        for i in range(P.NR):
            for j in range(P.NC):
                sl = (slice(P.R0 + P.WIN * i, P.R0 + P.WIN * (i + 1)), slice(P.C0 + P.WIN * j, P.C0 + P.WIN * (j + 1)))
                w = {"R": R[sl], "IR": GIR[sl], "MW": MW[sl], "M": M[sl]}
                if np.isfinite(w["R"]).mean() < 0.9 or not np.isfinite(w["IR"]).all() \
                        or not np.isfinite(w["M"]).all() or np.isfinite(w["MW"]).mean() < 0.8:
                    continue
                w = {k: np.nan_to_num(v) for k, v in w.items()}
                if min(v.sum() for v in w.values()) < 1:
                    continue
                row = {"t": t, "i": i, "j": j}
                for k, v in w.items():
                    rd = P.reads(v)
                    if k == "IR":  # GIR wet = GIRc > 0 (Tb < 235 K), per plan
                        rd["wet"] = float((v > 0).mean())
                        from scipy import ndimage
                        lab, n = ndimage.label(v > 0, structure=np.ones((3, 3)))
                        rd["isl"] = int((np.bincount(lab.ravel())[1:] >= 2).sum()) if n else 0
                    for name, val in rd.items():
                        row[f"{k}_{name}"] = val
                rows.append(row)
    d = pd.DataFrame(rows)
    lat = 55.0 - 0.1 * (P.R0 + P.WIN * d.i + P.WIN / 2); lon = -130.0 + 0.1 * (P.C0 + P.WIN * d.j + P.WIN / 2)
    d["clu"] = d.t.dt.floor("D").astype(str) + "|" + np.floor(lat / 5).astype(int).astype(str) + "_" + \
        np.floor(lon / 5).astype(int).astype(str)
    return d


if __name__ == "__main__":
    d = build(); d.to_csv("data/probe_a2_windows.csv", index=False)
    res = P.analyse(d)
    json.dump(res, open("data/probe_a2_results.json", "w"), indent=1)
    for v in P.VARS:
        p = res[v]["point"]; c = res[v]["ci95"]
        print(f"{v:8s} nf R={p['nf_R']:.3f} [{c['nf_R'][0]:.2f},{c['nf_R'][1]:.2f}] "
              f"GIR={p['nf_IR']:.3f} [{c['nf_IR'][0]:.2f},{c['nf_IR'][1]:.2f}] "
              f"MW={p['nf_MW']:.3f} [{c['nf_MW'][0]:.2f},{c['nf_MW'][1]:.2f}] | "
              f"errcorr M-GIR={p['M_IR']:.2f} [{c['M_IR'][0]:.2f},{c['M_IR'][1]:.2f}] "
              f"M-MW={p['M_MW']:.2f} [{c['M_MW'][0]:.2f},{c['M_MW'][1]:.2f}] | "
              f"neg={res[v]['share_boot_negative_ev']:.2f} unstable={res[v]['unstable']} | "
              f"r R-GIR={p['r_R_IR']:.2f} R-MW={p['r_R_MW']:.2f} GIR-MW={p['r_IR_MW']:.2f}")
    print("n", res["n"], "clusters", res["clusters"], "stop_rule", res["stop_rule_triggered"],
          "pred1", res["pred1_count"], "/5 pred2", res["pred2_count"], "/5")
