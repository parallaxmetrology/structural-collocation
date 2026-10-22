"""T1 synthetic anchor (T1_PLAN.md). Known truth, planted observer errors, planted shared error.
Reads and estimators are reused unchanged from probe_a.py."""
import json
import numpy as np, pandas as pd
import probe_a as P

N, SEED = 20000, 20261015
SIGMAS = [0.0, 0.5, 1.0, 1.5]
YY, XX = P.YY, P.XX
OBS = {  # name: (pos sd, size add, log-intensity sd)
    "R": (0.5, 0.0, 0.3), "S1": (0.8, 1.0, 0.4), "S2": (1.5, 2.0, 0.6), "G": (1.0, 1.5, 0.5)}
SHARED = {"S1", "S2"}


def render(cy, cx, s, I):
    f = np.zeros((P.WIN, P.WIN))
    for a, b, ss, ii in zip(cy, cx, s, I):
        f += ii * np.exp(-((YY - a) ** 2 + (XX - b) ** 2) / (2 * ss ** 2))
    f[f < 0.1] = 0.0
    return f


def simulate():
    rng = np.random.default_rng(SEED)
    wins = []
    for _ in range(N):
        K = 1 + rng.poisson(1.5)
        cells = dict(cy=rng.uniform(-5, 5, K), cx=rng.uniform(-5, 5, K), s=rng.uniform(0.8, 2.5, K),
                     I=rng.lognormal(1.5, 0.8, K))
        err = {o: dict(dy=rng.normal(0, sd, K), dx=rng.normal(0, sd, K), m=rng.lognormal(0, li, K))
               for o, (sd, _, li) in OBS.items()}
        unit = rng.normal(0, 1, 2)  # shared window displacement direction/size, scaled by sigma_sh
        wins.append((cells, err, unit))
    return wins


def fields(win, sig):
    cells, err, unit = win
    out = {"T": render(cells["cy"], cells["cx"], cells["s"], cells["I"])}
    for o, (_, add, _) in OBS.items():
        e = err[o]; sh = unit * sig if o in SHARED else np.zeros(2)
        out[o] = render(cells["cy"] + e["dy"] + sh[0], cells["cx"] + e["dx"] + sh[1], cells["s"] + add, cells["I"] * e["m"])
    out["M"] = 0.6 * out["S1"] + 0.4 * out["S2"]
    return out


def build(wins, sig):
    rows = []
    for k, w in enumerate(wins):
        f = fields(w, sig)
        if min(v.sum() for v in f.values()) < 1:
            continue
        row = {"k": k}
        for name, slot in [("T", "T"), ("R", "R"), ("S2", "IR"), ("S1", "MW"), ("M", "M"), ("G", "G")]:
            for var, val in P.reads(f[name]).items():
                row[f"{slot}_{var}"] = val
        rows.append(row)
    d = pd.DataFrame(rows)
    d["clu"] = (d.k // 20).astype(str)
    return d


def truth_scores(d):
    out = {}
    for v in P.VARS:
        T = d[f"T_{v}"].values
        res, nf = {}, {}
        for o in ["R", "IR", "MW", "G", "M"]:
            x = d[f"{o}_{v}"].values
            r = np.corrcoef(x, T)[0, 1]; nf[o] = 1 - r ** 2
            b = np.polyfit(T, x, 1); res[o] = x - np.polyval(b, T)
        ec = lambda a, b: float(np.corrcoef(res[a], res[b])[0, 1])
        out[v] = {"nf_true": {o: float(nf[o]) for o in nf},
                  "errcorr_true": {"MW_IR": ec("MW", "IR"), "M_IR": ec("M", "IR"), "M_MW": ec("M", "MW"), "R_IR": ec("R", "IR")}}
    return out


def quad_G(d, n_boot=300, seed=3):
    """c(S1,S2) with independent instruments R and G (P4)."""
    rng = np.random.default_rng(seed)
    clus, uc = pd.factorize(d.clu); idx = [np.where(clus == c)[0] for c in range(len(uc))]
    out = {}
    for v in P.VARS:
        X = d[[f"MW_{v}", f"IR_{v}", f"R_{v}", f"G_{v}"]].values
        def est(rows):
            C = np.cov(X[rows].T)
            p, q, u, w = 0, 1, 2, 3
            c = C[p, q] - C[p, u] * C[q, w] / C[u, w]
            ep = C[p, p] - C[p, u] * C[p, w] / C[u, w]; eq = C[q, q] - C[q, u] * C[q, w] / C[u, w]
            return c / np.sqrt(ep * eq) if ep > 0 and eq > 0 else np.nan
        pt = est(np.arange(len(d)))
        bs = [est(np.concatenate([idx[i] for i in rng.integers(0, len(idx), len(idx))])) for _ in range(n_boot)]
        out[v] = {"point": float(pt), "ci95": [float(np.nanpercentile(bs, 2.5)), float(np.nanpercentile(bs, 97.5))]}
    return out


if __name__ == "__main__":
    wins = simulate()
    allres = {}
    for sig in SIGMAS:
        d = build(wins, sig)
        tc = P.analyse(d)
        allres[str(sig)] = {"n": int(len(d)), "tc": tc, "truth": truth_scores(d), "quadG": quad_G(d)}
        t = allres[str(sig)]
        print(f"\n=== sigma_sh={sig}  n={len(d)}")
        for v in P.VARS:
            p = tc[v]["point"]; tr = t["truth"][v]["nf_true"]
            print(f"{v:8s} TC R={p['nf_R']:.3f}(true {tr['R']:.3f}) S2={p['nf_IR']:.3f}(true {tr['IR']:.3f}) "
                  f"S1={p['nf_MW']:.3f}(true {tr['MW']:.3f}) neg={tc[v]['share_boot_negative_ev']:.2f} | "
                  f"S1-S2 err: quadG={t['quadG'][v]['point']:.2f} true={t['truth'][v]['errcorr_true']['MW_IR']:.2f} | "
                  f"M-S1 {p['M_MW']:.2f} (true {t['truth'][v]['errcorr_true']['M_MW']:.2f}) "
                  f"M-S2 {p['M_IR']:.2f} (true {t['truth'][v]['errcorr_true']['M_IR']:.2f})")
    json.dump(allres, open("data/t1_results.json", "w"), indent=1)
