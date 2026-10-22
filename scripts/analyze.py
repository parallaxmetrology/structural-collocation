"""Feasibility analysis as frozen in FEASIBILITY_PLAN.md (v0 + Amendments 1 and 2).

Stage A: build gauge-hour samples and structural reads (real, and spatial-scramble null).
Stage B: Gate-6 checks, descriptive associations, resolution-floor sweep, and the
discovery -> validation increment test for each pre-stated model.
"""
import json, sys
import numpy as np, pandas as pd
from scipy import ndimage
from scipy.stats import spearmanr, rankdata, pearsonr

sys.path.insert(0, "../domain-qualification")
from independence_audit import independence_audit  # read-only use

HOURS = pd.date_range("2024-07-01 01:00", "2024-07-22 00:00", freq="h")
SPLIT = pd.Timestamp("2024-07-11 00:00")  # discovery: t_end <= SPLIT
W = 7                 # 15 x 15 window
SWEEP = (3, 5, 7, 10)  # resolution-floor half-widths
WET, DRY = 0.5, 0.1
rng = np.random.default_rng(20261015)

st = pd.read_csv("data/stations.csv")
st["key"] = st.src + ":" + st.sid.astype(str)
st["br"] = np.clip(np.floor((55.0 - st.lat) / 0.1).astype(int), 0, 349)
st["bc"] = np.clip(np.floor((st.lon + 130.0) / 0.1).astype(int), 0, 699)
g = pd.read_csv("data/gauges.csv", parse_dates=["t_end"])
g["si"] = (g.src + ":" + g.sid.astype(str)).map({k: i for i, k in enumerate(st.key)})
BR, BC = st.br.values, st.bc.values
YY, XX = np.mgrid[-W:W + 1, -W:W + 1]


def corr(a, b):
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a.ravel(), b.ravel())[0, 1])


def centroid(A, w):
    yy, xx = np.mgrid[-w:w + 1, -w:w + 1]
    s = A.sum()
    return np.array([(A * yy).sum(), (A * xx).sum()]) / s


def n_islands(A):
    lab, n = ndimage.label(A >= WET, structure=np.ones((3, 3)))
    if n == 0:
        return 0
    sizes = ndimage.sum(np.ones_like(A), lab, index=range(1, n + 1))
    return int((np.asarray(sizes) >= 2).sum())


def dispersion(A):
    c = centroid(A, W)
    d2 = (YY - c[0]) ** 2 + (XX - c[1]) ** 2
    return float(np.sqrt((A * d2).sum() / A.sum()) / np.hypot(2 * W + 1, 2 * W + 1))


def edge_dist(A):
    yx = np.argwhere(A >= WET)
    if len(yx) == 0:
        return 8.0
    return float(min(8.0, np.hypot(yx[:, 0] - W, yx[:, 1] - W).min()))


def sigma_space(A, B):
    vals = []
    for i in range(3):
        for j in range(3):
            a = A[5 * i:5 * i + 5, 5 * j:5 * j + 5]; b = B[5 * i:5 * i + 5, 5 * j:5 * j + 5]
            r = corr(a, b)
            if np.isfinite(r):
                vals.append(1 - r)
    return float(np.std(vals)) if len(vals) >= 3 else np.nan


def reads(Rw, Sw):
    """All structural disagreement reads between two 15x15 windows (both must carry mass)."""
    r1 = corr(Rw, Sw); r1 = 0.0 if not np.isfinite(r1) else r1
    Rs = ndimage.uniform_filter(Rw, 3, mode="nearest"); Ss = ndimage.uniform_filter(Sw, 3, mode="nearest")
    r3 = corr(Rs, Ss); r3 = 0.0 if not np.isfinite(r3) else r3
    mR, mS = Rw >= WET, Sw >= WET
    union = (mR | mS).sum()
    return dict(
        om_loc=float(np.linalg.norm(centroid(Rw, W) - centroid(Sw, W))),
        om_pat=1 - r1,
        om_iou=1 - ((mR & mS).sum() / union) if union else 0.0,
        om_islands=abs(n_islands(Rw) - n_islands(Sw)),
        om_disp=abs(dispersion(Rw) - dispersion(Sw)),
        om_void=abs(float((Rw < DRY).mean() - (Sw < DRY).mean())),
        om_scale=1 - r3,
        om_scale_rec=r3 - r1,
        om_edge=abs(edge_dist(Rw) - edge_dist(Sw)),
        om_sigma_space=sigma_space(Rw, Sw),
    )


def win(A, r0, c0, w=W):
    return A[r0 - w:r0 + w + 1, c0 - w:c0 + w + 1]


def build(scramble=False):
    rows = []
    for t in HOURS:
        m = np.load(f"data/mrms/{t:%Y%m%d%H}.npz")
        R, Rg, Qg = m["R"], m["Rg"], m["Qg"]
        t0 = t - pd.Timedelta(hours=1)
        S = (np.load(f"data/imerg/{t0:%Y%m%d%H}00.npy") + np.load(f"data/imerg/{t0:%Y%m%d%H}30.npy")) / 2
        xa, xb = np.load(f"data/imerg_x/{t0:%Y%m%d%H}00.npz"), np.load(f"data/imerg_x/{t0:%Y%m%d%H}30.npz")
        IR = (xa["IR"] + xb["IR"]) / 2
        MW = np.nanmean(np.stack([xa["MW"], xb["MW"]]), axis=0)
        # fast filter: window validity and mass for every cell
        R0 = np.nan_to_num(R); Rv = np.isfinite(R).astype(float)
        boxR = ndimage.uniform_filter(R0, 2 * W + 1, mode="constant") * (2 * W + 1) ** 2
        boxS = ndimage.uniform_filter(np.nan_to_num(S), 2 * W + 1, mode="constant") * (2 * W + 1) ** 2
        boxV = ndimage.uniform_filter(Rv, 2 * W + 1, mode="constant")
        Sbad = ndimage.uniform_filter((~np.isfinite(S)).astype(float), 2 * W + 1, mode="constant")
        gt = g[g.t_end == t]
        cand = []
        for si, gp in zip(gt.si.values, gt.p_mm.values):
            r0, c0 = BR[si], BC[si]
            if r0 < 10 or c0 < 10 or r0 >= 340 or c0 >= 690:
                continue
            if not (np.isfinite(Rg[si]) and np.isfinite(Qg[si]) and np.isfinite(R[r0, c0])):
                continue
            if boxV[r0, c0] < 0.9 or Sbad[r0, c0] > 0 or boxR[r0, c0] < 1 or boxS[r0, c0] < 1:
                continue
            cand.append((si, gp))
        if not cand:
            continue
        perm = np.arange(len(cand))
        if scramble and len(cand) > 1:
            perm = rng.permutation(len(cand))
        for k, (si, gp) in enumerate(cand):
            r0, c0 = BR[si], BC[si]
            sj = cand[perm[k]][0]; s0, d0 = BR[sj], BC[sj]   # satellite taken from (possibly) another gauge
            Rw = np.nan_to_num(win(R, r0, c0)); Sw = win(S, s0, d0); IRw = win(IR, s0, d0)
            if not np.isfinite(Sw).all() or Sw.sum() < 1:   # only possible under scramble
                continue
            row = dict(t=t, src=st.src.iat[si], si=si, lat=st.lat.iat[si], lon=st.lon.iat[si],
                       G=gp, Rpix=float(Rg[si]), Rcell=float(R[r0, c0]), Scell=float(S[s0, d0]),
                       RQI=float(Qg[si]), naive1=abs(float(R[r0, c0]) - float(S[s0, d0])),
                       naive2=float(np.abs(Rw - Sw).mean()),
                       sigR=bool(Rw.std() > 0), sigS=bool(Sw.std() > 0))
            row.update(reads(Rw, Sw))
            if np.isfinite(IRw).all() and IRw.sum() >= 1:
                ri = reads(Rw, IRw); row["om_loc_ir"] = ri["om_loc"]; row["om_pat_ir"] = ri["om_pat"]
            MWw = win(MW, s0, d0)
            if np.isfinite(MWw).mean() >= 0.8 and np.nansum(MWw) >= 1:
                rm = reads(Rw, np.nan_to_num(MWw)); row["om_loc_mw"] = rm["om_loc"]; row["om_pat_mw"] = rm["om_pat"]
            if not scramble:  # resolution-floor sweep (fixed sample: valid at w=10)
                R10 = win(R, r0, c0, 10); S10 = win(S, s0, d0, 10)
                if np.isfinite(R10).mean() >= 0.9 and np.isfinite(S10).all():
                    for w in SWEEP:
                        a = np.nan_to_num(win(R, r0, c0, w)); b = win(S, s0, d0, w)
                        if a.sum() > 0 and b.sum() > 0:
                            rr = corr(a, b)
                            row[f"pat_w{w}"] = 1 - (rr if np.isfinite(rr) else 0.0)
                            row[f"loc_w{w}"] = float(np.linalg.norm(centroid(a, w) - centroid(b, w)))
            rows.append(row)
    d = pd.DataFrame(rows)
    # causal temporal reads (same gauge, previous hours)
    d = d.sort_values(["si", "t"]).reset_index(drop=True)
    key = d.set_index(["si", "t"])
    def lag(col, h):
        idx = pd.MultiIndex.from_arrays([d.si, d.t - pd.Timedelta(hours=h)])
        return key[col].reindex(idx).values
    d["om_pat_persist"] = np.nanmean(np.vstack([d.om_pat.values, lag("om_pat", 1)]), axis=0)
    L = np.vstack([d.om_loc.values, lag("om_loc", 1), lag("om_loc", 2)])
    d["om_sigma_time"] = np.where(np.isfinite(L).all(axis=0), np.nanstd(L, axis=0), np.nan)
    d["E"] = (d.Rpix - d.G).abs()
    d["inc"] = 1 - d.RQI
    d["day"] = d.t.dt.floor("D")
    d["clu"] = d.day.astype(str) + "|" + (np.floor(d.lat / 5)).astype(int).astype(str) + "_" + \
        (np.floor(d.lon / 5)).astype(int).astype(str)
    d["disc"] = d.t <= SPLIT
    return d


BASE = ["inc", "naive1", "naive2", "Rpix", "Scell"]
V0 = ["om_loc", "om_pat"]
A1 = V0 + ["om_iou", "om_islands", "om_disp", "om_void", "om_scale", "om_scale_rec",
           "om_edge", "om_pat_persist", "om_loc_ir", "om_pat_ir"]
A2 = A1 + ["om_sigma_space", "om_sigma_time"]
VAR = ["om_sigma_space", "om_sigma_time"]
MODELS = {"frozen_v0": BASE + V0, "variance_only": BASE + VAR, "FULL_A2": BASE + A2}


def design(df, cols, ref):
    """Rank each column against the discovery distribution; missing -> discovery median rank
    plus a missingness indicator (fitted on discovery only)."""
    X = []
    for c in cols:
        refv = np.sort(ref[c].dropna().values)
        v = df[c].values.astype(float)
        r = np.searchsorted(refv, v, side="right") / len(refv)
        miss = ~np.isfinite(v)
        r[miss] = 0.5
        X.append(r)
        if ref[c].isna().any():
            X.append(miss.astype(float))
    return np.column_stack([np.ones(len(df))] + X)


def increment(d, cols_full, cols_base=BASE, n_boot=2000):
    D, V = d[d.disc], d[~d.disc]
    yD = rankdata(D.E) / len(D)
    pred = {}
    for name, cols in (("base", cols_base), ("full", cols_full)):
        beta, *_ = np.linalg.lstsq(design(D, cols, D), yD, rcond=None)
        pred[name] = design(V, cols, D) @ beta
    rb = spearmanr(pred["base"], V.E)[0]; rf = spearmanr(pred["full"], V.E)[0]
    clus = V.clu.values; uc = np.unique(clus); ix = {c: np.where(clus == c)[0] for c in uc}
    br = np.random.default_rng(7)
    boots = []
    for _ in range(n_boot):
        pick = np.concatenate([ix[c] for c in br.choice(uc, len(uc), replace=True)])
        e = V.E.values[pick]
        boots.append(spearmanr(pred["full"][pick], e)[0] - spearmanr(pred["base"][pick], e)[0])
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {"rho_base": float(rb), "rho_full": float(rf), "delta": float(rf - rb),
            "ci95": [float(lo), float(hi)], "go": bool(rf - rb >= 0.02 and lo > 0),
            "n_val": int(len(V)), "n_clusters_val": int(len(uc))}


def partial(df, x):
    s = df[[x, "E", "Rpix", "Scell"]].dropna()
    rk = lambda v: rankdata(v) / len(v)
    Z = np.column_stack([np.ones(len(s)), rk(s.Rpix), rk(s.Scell)])
    def resid(v):
        b, *_ = np.linalg.lstsq(Z, rk(v), rcond=None); return rk(v) - Z @ b
    return float(spearmanr(resid(s[x]), resid(s.E))[0]), int(len(s))


if __name__ == "__main__":
    d = build(); d.to_csv("data/samples.csv", index=False)
    dn = build(scramble=True); dn.to_csv("data/samples_null.csv", index=False)
    res = {"n_samples": int(len(d)), "n_disc": int(d.disc.sum()), "n_val": int((~d.disc).sum()),
           "by_src": d.src.value_counts().to_dict(),
           "clusters": {"disc": int(d[d.disc].clu.nunique()), "val": int(d[~d.disc].clu.nunique())},
           "signal_share": {"radar": float(d.sigR.mean()), "satellite": float(d.sigS.mean())},
           "near_zero_omloc_share": float((d.om_loc < 0.01).mean()),
           "coverage": {c: float(d[c].notna().mean()) for c in A2 + ["om_loc_mw"]}}
    # Gate 6
    eR, eS = (d.Rcell - d.G).values, (d.Scell - d.G).values
    res["gate6"] = {"error_r_pearson_radar_vs_sat": float(pearsonr(eR, eS)[0]),
                    "error_r_spearman_radar_vs_sat": float(spearmanr(eR, eS)[0]),
                    "field_r_cell": float(pearsonr(d.Rcell, d.Scell)[0])}
    full = d[BASE[:3] + A2].dropna()
    aud = independence_audit({k: full[k].values for k in BASE[:3] + A2})
    res["predictor_audit"] = {"n": int(len(full)), "I": aud["I"], "max_r": aud["max_r"],
                              "top_pairs": [(f"{a}~{b}", round(r, 3)) for (a, b), r in aud["pairs"][:8]]}
    # descriptive (discovery only, so validation stays unread until the decision step)
    D = d[d.disc]
    res["descriptive_discovery"] = {}
    for src in ("ALL", "USCRN", "ASOS"):
        s = D if src == "ALL" else D[D.src == src]
        out = {}
        for x in BASE[:3] + A2 + ["om_loc_mw", "om_pat_mw"]:
            ss = s[[x, "E"]].dropna()
            if len(ss) > 30:
                out[x] = {"rho": round(float(spearmanr(ss[x], ss.E)[0]), 4),
                          "partial": round(partial(s, x)[0], 4), "n": int(len(ss))}
        res["descriptive_discovery"][src] = out
    # resolution floor sweep (discovery, fixed sample valid at w=10)
    sw = D.dropna(subset=[f"pat_w{w}" for w in SWEEP] + [f"loc_w{w}" for w in SWEEP])
    res["resolution_floor_discovery"] = {"n": int(len(sw)), **{
        f"{k}_w{w}": round(partial(sw, f"{k}_w{w}")[0], 4) for k in ("pat", "loc") for w in SWEEP}}
    # decision (validation read once, here)
    res["decision"] = {name: increment(d, cols) for name, cols in MODELS.items()}
    res["null_spatial_scramble"] = {name: increment(dn, cols, n_boot=500) for name, cols in MODELS.items()}
    json.dump(res, open("data/results.json", "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))
