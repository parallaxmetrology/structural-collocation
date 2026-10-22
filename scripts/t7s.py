"""T7-S (T7S_PLAN.md): sensitivity of the satellite attribution to assumed radar-gauge error correlation.
One-factor model for R, G, I (GOES), M (microwave) with fixed rho_RG and free rho_IM, fitted by ML."""
import json
import numpy as np, pandas as pd
from scipy.optimize import minimize

OBS = ["R", "G", "IR", "MW"]  # order: 0 R, 1 G, 2 I, 3 M
GRID = [round(0.025 * k, 3) for k in range(21)]
SEED = 13


def sigma(th, rho_rg):
    lam = th[:4]; psi = np.exp(th[4:8]); rho_im = np.tanh(th[8])
    Psi = np.diag(psi)
    Psi[0, 1] = Psi[1, 0] = rho_rg * np.sqrt(psi[0] * psi[1])
    Psi[2, 3] = Psi[3, 2] = rho_im * np.sqrt(psi[2] * psi[3])
    return np.outer(lam, lam) + Psi, lam, psi, rho_im, Psi


def fit(S, rho_rg, x0=None):
    ldS = np.linalg.slogdet(S)[1]
    def F(th):
        Sg = sigma(th, rho_rg)[0]
        sgn, ld = np.linalg.slogdet(Sg)
        if sgn <= 0:
            return 1e6
        return ld + np.trace(S @ np.linalg.inv(Sg)) - ldS - 4
    if x0 is None:
        x0 = np.r_[np.sqrt(np.clip(np.diag(S) * 0.7, 1e-3, None)), np.log(np.diag(S) * 0.3), 0.0]
    best = None
    for start in ([x0] if x0 is not None else []) + [np.r_[np.sqrt(np.diag(S) * 0.5), np.log(np.diag(S) * 0.5), 0.3]]:
        r = minimize(F, start, method="L-BFGS-B", options={"maxiter": 2000})
        if best is None or r.fun < best.fun:
            best = r
    Sg, lam, psi, rho_im, Psi = sigma(best.x, rho_rg)
    nf = psi / (lam ** 2 + psi)
    admissible = bool(np.all(psi > 1e-4 * np.diag(S)) and abs(rho_im) < 0.999
                      and np.all(np.linalg.eigvalsh(Psi) > 0) and np.all(np.linalg.eigvalsh(Sg) > 0))
    return {"rho_IM": float(rho_im), "nf": dict(zip(OBS, map(float, nf))), "F": float(best.fun), "admissible": admissible, "x": best.x}


def run(s, var, n_boot=300):
    X = s[[f"{o}_{var}" for o in OBS]].values
    X = (X - X.mean(0)) / X.std(0)
    clus, uc = pd.factorize(s.clu); idx = [np.where(clus == c)[0] for c in range(len(uc))]
    rng = np.random.default_rng(SEED)
    boots = [np.concatenate([idx[k] for k in rng.integers(0, len(idx), len(idx))]) for _ in range(n_boot)]
    res = {}
    for g in GRID:
        p = fit(np.cov(X.T), g)
        bi, bd = [], []
        for rows in boots:
            q = fit(np.cov(X[rows].T), g, x0=p["x"])
            bi.append(q["rho_IM"]); bd.append(q["nf"]["R"] - q["nf"]["MW"])
        res[str(g)] = {"rho_IM": p["rho_IM"], "rho_IM_ci": [float(np.percentile(bi, 2.5)), float(np.percentile(bi, 97.5))],
                       "nf": p["nf"], "diff_R_MW": p["nf"]["R"] - p["nf"]["MW"],
                       "diff_ci": [float(np.percentile(bd, 2.5)), float(np.percentile(bd, 97.5))],
                       "F": p["F"], "admissible": p["admissible"]}
    adm = [g for g in GRID if res[str(g)]["admissible"]]
    t1 = next((g for g in GRID if res[str(g)]["rho_IM_ci"][0] <= 0 <= res[str(g)]["rho_IM_ci"][1]), None)
    t2 = next((g for g in GRID if res[str(g)]["diff_R_MW"] >= 0), None)
    return {"grid": res, "rho_star_1": t1, "rho_star_2": t2, "rho_max": max(adm) if adm else None}


if __name__ == "__main__":
    d = pd.read_csv("data/t7_points.csv", parse_dates=["t"])
    med = d.rpeak.median()
    out = {}
    for stratum, s in [("weak", d[d.rpeak <= med]), ("strong", d[d.rpeak > med])]:
        for var in ["cy", "cx"]:
            r = run(s, var)
            out[f"{stratum}_{var}"] = r
            g0 = r["grid"]["0.0"]
            print(f"{stratum} {var}: at rho_RG=0 rho_IM {g0['rho_IM']:+.3f} {np.round(g0['rho_IM_ci'],3)} nf R {g0['nf']['R']:.3f} MW {g0['nf']['MW']:.3f} G {g0['nf']['G']:.3f} | "
                  f"rho*1 {r['rho_star_1']} rho*2 {r['rho_star_2']} rho_max {r['rho_max']}", flush=True)
            for g in ["0.1", "0.2", "0.3", "0.4", "0.5"]:
                x = r["grid"][g]
                print(f"   rho_RG={g}: rho_IM {x['rho_IM']:+.3f} {np.round(x['rho_IM_ci'],3)} R-MW {x['diff_R_MW']:+.3f} adm {x['admissible']}", flush=True)
    json.dump(out, open("data/t7s_results.json", "w"), indent=1)
