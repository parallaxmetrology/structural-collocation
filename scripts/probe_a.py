"""Probe A: structural collocation, exactly as frozen in PROBE_A_PLAN.md."""
import json
import numpy as np, pandas as pd
import period as PERIOD
from scipy import ndimage

HOURS = PERIOD.HOURS
OBS = ["R", "IR", "MW", "M"]
VARS = ["wet", "cy", "cx", "disp", "isl", "logmean"]
WIN, R0, C0 = 15, 2, 2
NR, NC = (347 - R0) // WIN, (697 - C0) // WIN
YY, XX = np.mgrid[0:WIN, 0:WIN] - (WIN - 1) / 2


def reads(A):
    s = A.sum()
    cy, cx = (A * YY).sum() / s, (A * XX).sum() / s
    lab, n = ndimage.label(A >= 0.5, structure=np.ones((3, 3)))
    isl = int((np.bincount(lab.ravel())[1:] >= 2).sum()) if n else 0
    return dict(wet=float((A >= 0.5).mean()), cy=cy, cx=cx,
                disp=float(np.sqrt((A * ((YY - cy) ** 2 + (XX - cx) ** 2)).sum() / s)),
                isl=isl, logmean=float(np.log(A.mean() + 0.1)))


def build():
    rows = []
    for t in HOURS:
        R = np.load(f"data/mrms/{t:%Y%m%d%H}.npz")["R"]
        t0 = t - pd.Timedelta(hours=1)
        M = (np.load(f"data/imerg/{t0:%Y%m%d%H}00.npy") + np.load(f"data/imerg/{t0:%Y%m%d%H}30.npy")) / 2
        a, b = np.load(f"data/imerg_x/{t0:%Y%m%d%H}00.npz"), np.load(f"data/imerg_x/{t0:%Y%m%d%H}30.npz")
        IR = (a["IR"] + b["IR"]) / 2
        with np.errstate(all="ignore"):
            MW = np.nanmean(np.stack([a["MW"], b["MW"]]), axis=0)
        for i in range(NR):
            for j in range(NC):
                sl = (slice(R0 + WIN * i, R0 + WIN * (i + 1)), slice(C0 + WIN * j, C0 + WIN * (j + 1)))
                w = {"R": R[sl], "IR": IR[sl], "MW": MW[sl], "M": M[sl]}
                if np.isfinite(w["R"]).mean() < 0.9 or not np.isfinite(w["IR"]).all() \
                        or not np.isfinite(w["M"]).all() or np.isfinite(w["MW"]).mean() < 0.8:
                    continue
                w = {k: np.nan_to_num(v) for k, v in w.items()}
                if min(v.sum() for v in w.values()) < 1:
                    continue
                row = {"t": t, "i": i, "j": j}
                for k, v in w.items():
                    for name, val in reads(v).items():
                        row[f"{k}_{name}"] = val
                rows.append(row)
    d = pd.DataFrame(rows)
    lat = 55.0 - 0.1 * (R0 + WIN * d.i + WIN / 2); lon = -130.0 + 0.1 * (C0 + WIN * d.j + WIN / 2)
    d["clu"] = d.t.dt.floor("D").astype(str) + "|" + np.floor(lat / 5).astype(int).astype(str) + "_" + \
        np.floor(lon / 5).astype(int).astype(str)
    return d


def cov_from_sums(n, s, ss):
    """s: (k,), ss: (k,k) sums over rows -> covariance matrix."""
    m = s / n
    return ss / n - np.outer(m, m)


def tc(C, x, y, z):
    ev = C[x, x] - C[x, y] * C[x, z] / C[y, z]
    rho2 = C[x, y] * C[x, z] / (C[x, x] * C[y, z])
    return ev, 1 - rho2


def quad(C):
    """indices: 0 R, 1 IR, 2 MW, 3 M. Error correlations of M with IR and with MW."""
    out = {}
    for name, (p, q, u, v) in {"M_IR": (3, 1, 0, 2), "M_MW": (3, 2, 0, 1)}.items():
        # c(p,q) = C(p,q) - C(p,u) C(q,v) / C(u,v)   (u, v independent instruments)
        c = C[p, q] - C[p, u] * C[q, v] / C[u, v]
        ep = C[p, p] - C[p, u] * C[p, v] / C[u, v]
        eq = C[q, q] - C[q, u] * C[q, v] / C[u, v]
        out[name] = c / np.sqrt(ep * eq) if ep > 0 and eq > 0 else np.nan
    return out


def analyse(d, n_boot=500, seed=11):
    rng = np.random.default_rng(seed)
    clus, uc = pd.factorize(d.clu)
    nclu = len(uc)
    res = {"n": int(len(d)), "clusters": int(nclu)}
    for var in VARS:
        X = d[[f"{o}_{var}" for o in OBS]].values.astype(float)
        # per-cluster sufficient statistics
        n_c = np.bincount(clus, minlength=nclu).astype(float)
        s_c = np.stack([np.bincount(clus, X[:, k], minlength=nclu) for k in range(4)], 1)
        ss_c = np.einsum("ri,rj->rij", X, X)
        ss_c = np.stack([np.stack([np.bincount(clus, ss_c[:, a, b], minlength=nclu) for b in range(4)], 1)
                         for a in range(4)], 1)
        def stats(wts):
            n = (wts * n_c).sum(); s = (wts[:, None] * s_c).sum(0); ss = (wts[:, None, None] * ss_c).sum(0)
            C = cov_from_sums(n, s, ss)
            out = {}
            for o, (x, y, z) in {"R": (0, 1, 2), "IR": (1, 0, 2), "MW": (2, 0, 1)}.items():
                ev, nf = tc(C, x, y, z)
                out[f"ev_{o}"] = ev; out[f"nf_{o}"] = nf
            out.update(quad(C))
            out["r_R_IR"] = C[0, 1] / np.sqrt(C[0, 0] * C[1, 1])
            out["r_R_MW"] = C[0, 2] / np.sqrt(C[0, 0] * C[2, 2])
            out["r_IR_MW"] = C[1, 2] / np.sqrt(C[1, 1] * C[2, 2])
            return out
        point = stats(np.ones(nclu))
        boots = [stats(rng.multinomial(nclu, np.ones(nclu) / nclu).astype(float)) for _ in range(n_boot)]
        B = pd.DataFrame(boots)
        neg = float(((B[["ev_R", "ev_IR", "ev_MW"]] < 0).any(axis=1)).mean())
        ci = {k: [float(B[k].quantile(0.025)), float(B[k].quantile(0.975))] for k in B.columns}
        widths = {o: ci[f"nf_{o}"][1] - ci[f"nf_{o}"][0] for o in ("R", "IR", "MW")}
        unstable = neg > 0.20 or max(widths.values()) > 0.5
        res[var] = {"point": {k: float(v) for k, v in point.items()}, "ci95": ci,
                    "share_boot_negative_ev": neg, "nf_ci_width": widths, "unstable": bool(unstable),
                    "pred1_radar_smallest_nf": bool(point["nf_R"] < min(point["nf_IR"], point["nf_MW"])),
                    "pred2_IR_gt_MW": bool(point["nf_IR"] > point["nf_MW"])}
    struct = [v for v in VARS if v != "logmean"]
    res["stop_rule_triggered"] = bool(sum(res[v]["unstable"] for v in struct) > len(struct) / 2)
    res["pred1_count"] = sum(res[v]["pred1_radar_smallest_nf"] for v in struct)
    res["pred2_count"] = sum(res[v]["pred2_IR_gt_MW"] for v in struct)
    return res


if __name__ == "__main__":
    d = build(); d.to_csv("data/probe_a_windows.csv", index=False)
    res = analyse(d)
    json.dump(res, open("data/probe_a_results.json", "w"), indent=1)
    for v in VARS:
        p = res[v]["point"]; c = res[v]["ci95"]
        print(f"{v:8s} nf R={p['nf_R']:.3f} [{c['nf_R'][0]:.2f},{c['nf_R'][1]:.2f}] "
              f"IR={p['nf_IR']:.3f} [{c['nf_IR'][0]:.2f},{c['nf_IR'][1]:.2f}] "
              f"MW={p['nf_MW']:.3f} [{c['nf_MW'][0]:.2f},{c['nf_MW'][1]:.2f}] | "
              f"errcorr M-IR={p['M_IR']:.2f} [{c['M_IR'][0]:.2f},{c['M_IR'][1]:.2f}] "
              f"M-MW={p['M_MW']:.2f} [{c['M_MW'][0]:.2f},{c['M_MW'][1]:.2f}] | "
              f"neg={res[v]['share_boot_negative_ev']:.2f} unstable={res[v]['unstable']}")
    print("n", res["n"], "clusters", res["clusters"], "stop_rule", res["stop_rule_triggered"],
          "pred1", res["pred1_count"], "/5 pred2", res["pred2_count"], "/5")
