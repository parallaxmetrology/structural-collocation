"""Figures for LABBOOK.html. Reads only files already produced by the pipeline (data/*.json, *.csv,
reduced arrays). Writes figures/*.png. Deterministic."""
import json, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

os.makedirs("figures", exist_ok=True)
plt.rcParams.update({"font.size": 9, "font.family": "sans-serif", "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 130, "savefig.bbox": "tight"})
INK, RUST, MOSS, SAND, SLATE = "#1a1a1a", "#7a3b2e", "#3f6b3a", "#b59a5b", "#4a6280"
OBSCOL = {"R": SLATE, "IR": SAND, "GIR": SAND, "MW": MOSS, "M": RUST}
EXT = [-130, -60, 20, 55]

# ---------- F1: one hour through every observer ----------
t = pd.Timestamp("2024-07-21 23:00"); t0 = t - pd.Timedelta(hours=1)
R = np.load(f"data/mrms/{t:%Y%m%d%H}.npz")["R"]
M = (np.load(f"data/imerg/{t0:%Y%m%d%H}00.npy") + np.load(f"data/imerg/{t0:%Y%m%d%H}30.npy")) / 2
a, b = np.load(f"data/imerg_x/{t0:%Y%m%d%H}00.npz"), np.load(f"data/imerg_x/{t0:%Y%m%d%H}30.npz")
IR = (a["IR"] + b["IR"]) / 2
with np.errstate(all="ignore"):
    MW = np.nanmean(np.stack([a["MW"], b["MW"]]), axis=0)
Tb = np.load(f"data/goes/{t:%Y%m%d%H}.npz")["Tb"]
GIR = np.where(np.isfinite(Tb), np.maximum(0, 235 - Tb), np.nan)
panels = [("MRMS radar-only QPE (mm/h)", R, "rain"), ("IMERG merged (mm/h)", M, "rain"),
          ("IMERG infrared estimate, PDIR (mm/h)", IR, "rain"), ("IMERG microwave (mm/h; swath gaps blank)", MW, "rain"),
          ("GOES-16 band 13 cold cloud, 235 − Tb (K)", GIR, "cold")]
fig, axs = plt.subplots(3, 2, figsize=(10, 8.6))
for ax, (title, A, kind) in zip(axs.ravel(), panels):
    cm = plt.get_cmap("YlGnBu" if kind == "rain" else "magma_r").copy(); cm.set_bad("#e9e5dc")
    A2 = np.where(A > (0.1 if kind == "rain" else 0.5), A, np.nan) if True else A
    base = np.where(np.isfinite(A), 0.0, np.nan)
    ax.imshow(base, extent=EXT, cmap=matplotlib.colors.ListedColormap(["#f7f5ef"]), origin="upper")
    im = ax.imshow(A2, extent=EXT, cmap=cm, origin="upper",
                   norm=LogNorm(0.1, 50) if kind == "rain" else None, vmin=None if kind == "rain" else 0,
                   vmax=None if kind == "rain" else 40)
    ax.set_title(title, fontsize=9); ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(im, ax=ax, fraction=0.03, pad=0.01)
ax = axs.ravel()[5]
w = pd.read_csv("data/probe_a2_windows.csv", parse_dates=["t"])
wt = w[w.t == t]
ax.imshow(np.where(np.isfinite(R), 0.0, np.nan), extent=EXT, cmap=matplotlib.colors.ListedColormap(["#f7f5ef"]), origin="upper")
for _, r in wt.iterrows():
    y0 = 55 - 0.1 * (2 + 15 * r.i); x0 = -130 + 0.1 * (2 + 15 * r.j)
    ax.add_patch(plt.Rectangle((x0, y0 - 1.5), 1.5, 1.5, fill=True, fc=RUST, alpha=0.55, ec="none"))
ax.set_title(f"Probe A2 sample windows this hour (n={len(wt)})", fontsize=9); ax.set_xticks([]); ax.set_yticks([])
ax.set_xlim(EXT[:2]); ax.set_ylim(EXT[2:])
fig.suptitle(f"One hour through every observer: accumulation ending {t:%Y-%m-%d %H:%M} UTC", fontsize=11)
fig.savefig("figures/f1_observers.png"); plt.close(fig)

# ---------- F2: feasibility decision ----------
res = json.load(open("data/results.json"))
fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.8), gridspec_kw={"width_ratios": [1.1, 1.6, 1], "wspace": 0.45})
ax = axs[0]
names = ["frozen_v0", "variance_only", "FULL_A2"]; lab = ["frozen v0\n(2 reads)", "variance\nonly", "FULL A2\n(14 reads)"]
for k, (src, col, off) in enumerate([("decision", RUST, -0.17), ("null_spatial_scramble", SLATE, 0.17)]):
    d = [res[src][n]["delta"] for n in names]; lo = [res[src][n]["ci95"][0] for n in names]; hi = [res[src][n]["ci95"][1] for n in names]
    x = np.arange(3) + off
    ax.bar(x, d, 0.32, color=col, label="real" if k == 0 else "spatial-scramble null")
    ax.errorbar(x, d, yerr=[np.array(d) - lo, np.array(hi) - d], fmt="none", ecolor=INK, lw=0.8, capsize=2)
ax.axhline(0.02, color=RUST, ls="--", lw=1); ax.text(2.4, 0.0205, "go bar +0.02", color=RUST, fontsize=8, ha="right", va="bottom")
ax.set_xticks(range(3)); ax.set_xticklabels(lab, fontsize=8); ax.set_ylabel("Δρ vs baseline (validation)")
ax.set_ylim(-0.002, 0.024); ax.legend(fontsize=7, frameon=False, loc="center left"); ax.set_title("A · increment over baseline")
ax = axs[1]
dd = res["descriptive_discovery"]
keys = ["inc", "naive1", "naive2", "om_loc", "om_pat", "om_iou", "om_islands", "om_disp", "om_void", "om_scale",
        "om_scale_rec", "om_edge", "om_pat_persist", "om_loc_ir", "om_pat_ir", "om_sigma_space", "om_sigma_time"]
lbl = ["1−RQI", "|R−S| cell", "|R−S| window"] + [k.replace("om_", "") for k in keys[3:]]
y = np.arange(len(keys))
ax.barh(y + 0.2, [dd["ALL"][k]["partial"] for k in keys], 0.4, color=SLATE, label="all gauges")
ax.barh(y - 0.2, [dd["USCRN"][k]["partial"] for k in keys], 0.4, color=SAND, label="USCRN only")
ax.axvline(0, color=INK, lw=0.6); ax.set_yticks(y); ax.set_yticklabels(lbl, fontsize=7.5); ax.invert_yaxis()
ax.set_xlabel("partial Spearman with radar error | magnitude (discovery)"); ax.legend(fontsize=7, frameon=False)
ax.set_title("B · each read, rain amount controlled")
ax = axs[2]
rf = res["resolution_floor_discovery"]; ws = [3, 5, 7, 10]
ax.plot([w * 0.1 for w in ws], [rf[f"pat_w{w}"] for w in ws], "o-", color=RUST, label="om_pat")
ax.plot([w * 0.1 for w in ws], [rf[f"loc_w{w}"] for w in ws], "s-", color=SLATE, label="om_loc")
ax.axhline(0, color=INK, lw=0.6); ax.set_xlabel("window half-width (°)"); ax.set_ylabel("partial ρ (predicted sign: +)")
ax.legend(fontsize=7, frameon=False); ax.set_title(f"C · resolution sweep (n={rf['n']:,})")
fig.savefig("figures/f2_feasibility.png"); plt.close(fig)

# ---------- F3: anatomy of the label ----------
s = pd.read_csv("data/samples.csv", usecols=["Rpix", "G", "E", "src", "om_iou", "om_loc"])
fig, axs = plt.subplots(1, 3, figsize=(13, 3.6), gridspec_kw={"wspace": 0.35})
ax = axs[0]
hb = ax.hexbin(np.log10(s.Rpix + 0.1), np.log10(s.G + 0.1), gridsize=60, bins="log", cmap="Greys")
ax.plot([-1, 2], [-1, 2], color=RUST, lw=0.8); ax.set_xlabel("log10(radar at gauge pixel + 0.1)"); ax.set_ylabel("log10(gauge + 0.1)")
ax.set_title("A · radar vs gauge, hourly (n=%s)" % f"{len(s):,}")
ax = axs[1]
wet = s[s.Rpix > 0.1]
edges = np.logspace(-1, 2, 16); mid = []; mean = []; q3 = []; q1 = []
for lo_, hi_ in zip(edges[:-1], edges[1:]):
    e = wet.E[(wet.Rpix >= lo_) & (wet.Rpix < hi_)]
    if len(e) > 50:
        mid.append(np.sqrt(lo_ * hi_)); mean.append(e.mean()); q1.append(e.quantile(0.25)); q3.append(e.quantile(0.75))
ax.fill_between(mid, q1, q3, color=SAND, alpha=0.5, label="interquartile")
ax.plot(mid, mean, color=RUST, label="mean")
ax.set_xscale("log"); ax.set_yscale("log"); ax.legend(fontsize=7, frameon=False)
ax.set_xlabel("radar at gauge pixel (mm/h; > 0.1)"); ax.set_ylabel("|radar − gauge| (mm)")
ax.set_title(f"B · error scales with amount (n={len(wet):,})")
ax = axs[2]
q = pd.qcut(s.om_iou.rank(method="first"), 10, labels=False)
g_ = s.groupby(q).agg(iou=("om_iou", "mean"), E=("E", "mean"), R=("Rpix", "mean"))
ax.plot(g_.iou, g_.E, "o-", color=RUST, label="mean |radar − gauge|")
ax.plot(g_.iou, g_.R, "s--", color=SLATE, label="mean radar rain at gauge")
ax.set_yscale("log"); ax.set_xlabel("om_iou decile (1 − overlap of wet areas)"); ax.set_ylabel("mm"); ax.legend(fontsize=7, frameon=False)
ax.set_title("C · more disagreement, lighter rain")
fig.savefig("figures/f3_label.png"); plt.close(fig)

# ---------- F4: collocation noise fractions, A vs A2 ----------
A = json.load(open("data/probe_a_results.json")); A2 = json.load(open("data/probe_a2_results.json"))
VARS = ["wet", "cy", "cx", "disp", "isl", "logmean"]
VL = ["wet fraction", "centroid y", "centroid x", "dispersion", "islands", "log mean\n(value)"]
fig, axs = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True)
for ax, res_, irname, title in [(axs[0], A, "IR", "Probe A · (R, IMERG IR, MW)"), (axs[1], A2, "GIR", "Probe A2 · (R, GOES IR, MW)")]:
    x = np.arange(len(VARS))
    for k, (o, off) in enumerate([("R", -0.27), ("IR", 0), ("MW", 0.27)]):
        v = [res_[q]["point"][f"nf_{o}"] for q in VARS]
        lo = [res_[q]["ci95"][f"nf_{o}"][0] for q in VARS]; hi = [res_[q]["ci95"][f"nf_{o}"][1] for q in VARS]
        nm = irname if o == "IR" else o
        ax.bar(x + off, v, 0.26, color=OBSCOL[nm], label={"R": "radar", "IR": irname, "MW": "microwave"}[o])
        ax.errorbar(x + off, v, yerr=[np.clip(np.array(v) - lo, 0, None), np.clip(np.array(hi) - v, 0, None)],
                    fmt="none", ecolor=INK, lw=0.7, capsize=1.5)
    for q, xx in zip(VARS, x):
        if res_[q]["unstable"]:
            ax.text(xx, -0.105, "unstable", ha="center", fontsize=7, color=RUST)
    ax.axhline(0, color=INK, lw=0.6); ax.set_xticks(x); ax.set_xticklabels(VL, fontsize=8)
    ax.set_title(f"{title}  n={res_['n']:,}"); ax.legend(fontsize=7, frameon=False, loc="upper center", ncol=3)
axs[0].set_ylabel("noise fraction 1 − ρ² (lower = closer to common signal)"); axs[0].set_ylim(-0.13, 1.08)
fig.savefig("figures/f4_collocation.png"); plt.close(fig)

# ---------- F5: the arithmetic (raw correlations) + stability + quadruple check ----------
fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.8), gridspec_kw={"wspace": 0.42})
ax = axs[0]; x = np.arange(len(VARS))
for k, (key, col, lab_) in enumerate([("r_R_IR", SAND, "radar–GOES IR"), ("r_R_MW", SLATE, "radar–MW"), ("r_IR_MW", MOSS, "GOES IR–MW")]):
    ax.bar(x + (k - 1) * 0.27, [A2[q]["point"][key] for q in VARS], 0.26, color=col, label=lab_)
ax.set_xticks(x); ax.set_xticklabels([v.replace("\n", " ") for v in VL], fontsize=7, rotation=30, ha="right")
ax.set_ylabel("Pearson r (Probe A2)"); ax.legend(fontsize=7, frameon=False); ax.set_title("A · who agrees with whom")
ax = axs[1]
ax.bar(x - 0.18, [A[q]["share_boot_negative_ev"] for q in VARS], 0.34, color=SAND, label="Probe A")
ax.bar(x + 0.18, [A2[q]["share_boot_negative_ev"] for q in VARS], 0.34, color=SLATE, label="Probe A2")
ax.axhline(0.2, color=RUST, ls="--", lw=0.9); ax.text(3.0, 0.22, "stop threshold", fontsize=7, color=RUST, ha="right")
ax.set_xticks(x); ax.set_xticklabels([v.replace("\n", " ") for v in VL], fontsize=7, rotation=30, ha="right")
ax.set_ylabel("share of bootstrap draws with an\nimpossible (negative) error variance"); ax.legend(fontsize=7, frameon=False)
ax.set_title("B · structure stable, value not")
ax = axs[2]; SV = ["cy", "cx", "disp", "isl"]; xs = np.arange(len(SV))
for k, (res_, nm, col, off) in enumerate([(A, "A", SAND, -0.2), (A2, "A2", SLATE, 0.2)]):
    for key, mk in [("M_MW", "o"), ("M_IR", "^")]:
        v = [res_[q]["point"][key] for q in SV]; lo = [res_[q]["ci95"][key][0] for q in SV]; hi = [res_[q]["ci95"][key][1] for q in SV]
        lo = np.maximum(lo, -0.6)
        ax.errorbar(xs + off + (0.06 if key == "M_IR" else -0.06), v, yerr=[np.array(v) - lo, np.array(hi) - v], fmt=mk,
                    color=col, ms=4, lw=0.8, capsize=1.5, label=f"{nm}: merged–{'MW' if key == 'M_MW' else 'IR'}")
ax.axhline(0, color=INK, lw=0.6); ax.set_xticks(xs); ax.set_xticklabels(["centroid y", "centroid x", "dispersion", "islands"], fontsize=7)
ax.set_ylim(-0.6, 0.6); ax.set_ylabel("estimated error correlation"); ax.legend(fontsize=6.5, frameon=False, ncol=2)
ax.set_title("C · lineage detected without truth")
fig.savefig("figures/f5_arithmetic.png"); plt.close(fig)

# ---------- F6: post-hoc footprint check ----------
ph = json.load(open("data/probe_a_posthoc.json"))
fig, ax = plt.subplots(figsize=(6, 3))
SV = ["cy", "cx", "disp", "isl"]; xs = np.arange(len(SV))
for k, (key, col, lab_) in enumerate([("one_half", SAND, f"MW in one half-hour (n={ph['n_one']:,})"),
                                      ("both_halves", SLATE, f"MW in both half-hours (n={ph['n_both']:,})")]):
    v = [ph[key][q]["MW"] for q in SV]; ci = [ph[key][q]["nf_MW_ci"] for q in SV]
    ax.bar(xs + (k - 0.5) * 0.36, v, 0.34, color=col, label=lab_)
    ax.errorbar(xs + (k - 0.5) * 0.36, v, yerr=[[vv - c[0] for vv, c in zip(v, ci)], [c[1] - vv for vv, c in zip(v, ci)]],
                fmt="none", ecolor=INK, lw=0.7, capsize=1.5)
ax.axhline(0, color=INK, lw=0.6); ax.set_xticks(xs); ax.set_xticklabels(["centroid y", "centroid x", "dispersion", "islands"])
ax.set_ylabel("MW noise fraction"); ax.legend(fontsize=7, frameon=False)
ax.set_title("Post-hoc: footprint hypothesis predicted 'both' > 'one'. Observed the reverse.", fontsize=9)
fig.savefig("figures/f6_posthoc.png"); plt.close(fig)

# ---------- F7: where the evidence came from ----------
st = pd.read_csv("data/stations.csv")
fig, axs = plt.subplots(1, 2, figsize=(11, 3.4))
ax = axs[0]
ax.scatter(st[st.src == "ASOS"].lon, st[st.src == "ASOS"].lat, s=2, color=SLATE, label=f"ASOS ({(st.src=='ASOS').sum():,})")
ax.scatter(st[st.src == "USCRN"].lon, st[st.src == "USCRN"].lat, s=9, color=RUST, marker="^", label=f"USCRN ({(st.src=='USCRN').sum()})")
ax.set_xlim(-128, -65); ax.set_ylim(24, 50); ax.set_aspect(1.25); ax.legend(fontsize=7, frameon=False, loc="lower left")
ax.set_title("Gauges used as labels (feasibility test)")
ax = axs[1]
cnt = np.zeros((23, 46))
for (i, j), n in w.groupby(["i", "j"]).size().items():
    cnt[i, j] = n
im = ax.imshow(np.where(cnt > 0, cnt, np.nan), extent=[-129.8, -60.8, 20.3, 54.8], cmap="YlOrBr", origin="upper")
plt.colorbar(im, ax=ax, fraction=0.025, pad=0.01, label="window-hours")
ax.set_title(f"Probe A2 sample: 1.5° windows, July 2024 (n={len(w):,})")
fig.savefig("figures/f7_coverage.png"); plt.close(fig)
print("figures:", sorted(os.listdir("figures")))

# ---------- F8: T1 synthetic anchor ----------
if os.path.exists("data/t1_results.json"):
    T1 = json.load(open("data/t1_results.json")); TX = json.load(open("data/t1_explore.json"))
    sig = ["0.0", "0.5", "1.0", "1.5"]; sx = [float(s_) for s_ in sig]
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.9), gridspec_kw={"wspace": 0.35})
    ax = axs[0]
    mk = {"cy": "o", "cx": "o", "disp": "s", "isl": "^", "wet": "v", "logmean": "D"}
    vc = {"cy": SLATE, "cx": SLATE, "disp": MOSS, "isl": RUST, "wet": SAND, "logmean": INK}
    for v in ["cy", "cx", "disp", "isl", "wet", "logmean"]:
        for slot in ["R", "IR", "MW"]:
            ax.scatter(T1["0.0"]["truth"][v]["nf_true"][slot], T1["0.0"]["tc"][v]["point"][f"nf_{slot}"], marker=mk[v],
                       color=vc[v], s=28, label=VN_ if (VN_ := {"cy": "centroid", "disp": "dispersion", "isl": "islands",
                       "wet": "wet fraction", "logmean": "log mean (value)"}.get(v)) and slot == "R" else None)
    ax.plot([0, 1], [0, 1], color=INK, lw=0.7); ax.set_xlim(-0.05, 1.02); ax.set_ylim(-0.05, 1.02)
    ax.set_xlabel("true noise fraction (known truth)"); ax.set_ylabel("collocation estimate")
    ax.legend(fontsize=7, frameon=False, loc="upper left"); ax.set_title("A · S0, all assumptions met")
    ax = axs[1]
    for slot, col, nm in [("R", SLATE, "radar-like R"), ("MW", MOSS, "microwave-like S1"), ("IR", SAND, "infrared-like S2")]:
        ax.plot(sx, [T1[s_]["tc"]["cy"]["point"][f"nf_{slot}"] for s_ in sig], "o-", color=col, label=f"{nm}: estimate")
        ax.plot(sx, [T1[s_]["truth"]["cy"]["nf_true"][slot] for s_ in sig], "o--", color=col, alpha=0.6, label=f"{nm}: truth")
    ax.set_xlabel("planted shared S1–S2 displacement σ_sh (cells)"); ax.set_ylabel("noise fraction, centroid y")
    ax.legend(fontsize=6.5, frameon=False, ncol=1); ax.set_title("B · a shared error flips the ranking")
    ax = axs[2]
    ax.plot(sx, [T1[s_]["truth"]["cy"]["errcorr_true"]["MW_IR"] for s_ in sig], "o--", color=INK, label="true S1–S2 error corr.")
    ax.plot(sx, [T1[s_]["quadG"]["cy"]["point"] for s_ in sig], "s-", color=RUST, label="estimate, blurred instrument G")
    ax.plot([0.0, 1.0], [TX["E1_sig0.0"]["cy"]["quadG"]["point"], TX["E1_sig1.0"]["cy"]["quadG"]["point"]], "D", color=MOSS,
            ms=7, label="estimate, sharp instrument G (exploratory)")
    ax.axhline(0, color=INK, lw=0.5); ax.set_xlabel("σ_sh (cells)"); ax.set_ylabel("error correlation, centroid y")
    ax.legend(fontsize=6.5, frameon=False); ax.set_title("C · detecting the shared error")
    fig.savefig("figures/f8_t1.png"); plt.close(fig)
    print("f8 written")

# ---------- F9: T5 ladder, synthetic arm, parallax ----------
if os.path.exists("data/t5_results.json"):
    T5 = json.load(open("data/t5_results.json")); OFF = json.load(open("data/t5_offsets.json"))
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.9), gridspec_kw={"wspace": 0.35})
    ax = axs[0]
    lv = sorted(T5["real"], key=float); x = [float(l) for l in lv]
    for v, col, nm in [("cy", SLATE, "centroid y"), ("cx", MOSS, "centroid x"), ("disp", SAND, "dispersion")]:
        y = [T5["real"][l][v]["diff_R_minus_MW"] for l in lv]
        lo = [T5["real"][l][v]["diff_ci95"][0] for l in lv]; hi = [T5["real"][l][v]["diff_ci95"][1] for l in lv]
        ax.plot(x, y, "o-", color=col, label=nm); ax.fill_between(x, lo, hi, color=col, alpha=0.2)
    ax.axhline(0, color=INK, lw=0.7); ax.axvline(T5["sigma_star"], color=RUST, ls="--", lw=0.9)
    ax.text(T5["sigma_star"] + 0.05, 0.02, "σ* (matches MW\nsmoothness)", color=RUST, fontsize=7)
    ax.set_xlabel("radar blur σ (0.1° cells)"); ax.set_ylabel("nf_radar − nf_microwave")
    ax.legend(fontsize=7, frameon=False); ax.set_title("A · real data: the inversion survives blurring")
    ax = axs[1]
    for sh, col, nm in [("0.0", MOSS, "smoothness only"), ("1.0", RUST, "shared displacement")]:
        ys = [T5["synthetic"][f"sh{sh}_Rblur{a}"]["cy"]["diff_R_minus_MW"] for a in ["0.0", "1.0", "2.0"]]
        tr = [T5["synthetic"][f"sh{sh}_Rblur{a}"]["cy"]["true_R"] - T5["synthetic"][f"sh{sh}_Rblur{a}"]["cy"]["true_MW"] for a in ["0.0", "1.0", "2.0"]]
        ax.plot([0, 1, 2], ys, "o-", color=col, label=f"{nm}: estimate"); ax.plot([0, 1, 2], tr, "o--", color=col, alpha=0.5, label=f"{nm}: truth")
    ax.axhline(0, color=INK, lw=0.7); ax.set_xlabel("radar-like blur (cells)"); ax.set_ylabel("nf_R − nf_S1, centroid y")
    ax.legend(fontsize=6.5, frameon=False); ax.set_title("B · synthetic calibration")
    ax = axs[2]
    bands = list(OFF["by_latitude_cy"]); mids = [25, 32.5, 37.5, 42.5, 50]
    ax.plot(mids, [OFF["by_latitude_cy"][b]["GIR_minus_R"] for b in bands], "o-", color=SAND, label="GOES IR − radar")
    ax.plot(mids, [OFF["by_latitude_cy"][b]["MW_minus_R"] for b in bands], "s-", color=MOSS, label="microwave − radar")
    ax.axhline(0, color=INK, lw=0.7); ax.invert_yaxis(); ax.set_xlabel("latitude (°N)"); ax.set_ylabel("mean centroid offset (cells; up = north)")
    ax.legend(fontsize=7, frameon=False); ax.set_title("C · GOES parallax seen, not shared (exploratory)")
    fig.savefig("figures/f9_t5.png"); plt.close(fig)
    print("f9 written")

# ---------- F10: T6 core vs area, intensity, RQI ----------
if os.path.exists("data/t6_results.json"):
    T6 = json.load(open("data/t6_results.json")); TI = json.load(open("data/t6_intensity.json")); TQ = json.load(open("data/t6_rqi.json"))
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.9), gridspec_kw={"wspace": 0.35})
    ax = axs[0]; xs = np.arange(2)
    for k, (key, col, nm) in enumerate([("full", SLATE, "full area"), ("core", RUST, "core only (top 20%)")]):
        v = [T6[key][q]["diff_R_minus_MW"] for q in ["cy", "cx"]]; c = [T6[key][q]["diff_ci95"] for q in ["cy", "cx"]]
        ax.bar(xs + (k - 0.5) * 0.36, v, 0.34, color=col, label=nm)
        ax.errorbar(xs + (k - 0.5) * 0.36, v, yerr=[[a - b[0] for a, b in zip(v, c)], [b[1] - a for a, b in zip(v, c)]], fmt="none", ecolor=INK, lw=0.8, capsize=2)
    ax.set_xticks(xs); ax.set_xticklabels(["centroid y", "centroid x"]); ax.set_ylabel("nf_radar − nf_microwave")
    ax.legend(fontsize=7, frameon=False); ax.set_title("A · reading cores does not remove it")
    ax = axs[1]
    e = TI["edges_mm_h"]; mids = [str(round(x, 1)) for x in e]
    for v, col in [("cy", SLATE), ("cx", MOSS)]:
        y = [TI[str(q)][v]["diff"] for q in range(5)]; lo = [TI[str(q)][v]["ci"][0] for q in range(5)]; hi = [TI[str(q)][v]["ci"][1] for q in range(5)]
        ax.plot(range(5), y, "o-", color=col, label={"cy": "centroid y", "cx": "centroid x"}[v]); ax.fill_between(range(5), lo, hi, color=col, alpha=0.2)
    ax.set_xticks(range(5)); ax.set_xticklabels([f"≤{m}" for m in mids], fontsize=8); ax.set_xlabel("radar window peak quintile (mm/h)")
    ax.axhline(0, color=INK, lw=0.6); ax.set_ylabel("nf_radar − nf_microwave"); ax.legend(fontsize=7, frameon=False)
    ax.set_title("B · it lives in light rain (exploratory)")
    ax = axs[2]
    for b, col, nm in [(0, SAND, "weak rain"), (1, SLATE, "strong rain")]:
        y = [TQ[f"rqi{a}_peak{b}"]["cy"]["diff"] for a in range(3)]; c = [TQ[f"rqi{a}_peak{b}"]["cy"]["ci"] for a in range(3)]
        xx = np.arange(3) + (b - 0.5) * 0.3
        ax.errorbar(xx, y, yerr=[[a - cc[0] for a, cc in zip(y, c)], [cc[1] - a for a, cc in zip(y, c)]], fmt="o", color=col, capsize=2, label=nm)
    ax.set_xticks(range(3)); ax.set_xticklabels([f"RQI {TQ[f'rqi{a}_peak0']['rqi_mean']:.2f}" for a in range(3)])
    ax.axhline(0, color=INK, lw=0.6); ax.set_ylabel("nf_radar − nf_microwave (centroid y)"); ax.legend(fontsize=7, frameon=False)
    ax.set_title("C · not radar beam quality (exploratory)")
    fig.savefig("figures/f10_t6.png"); plt.close(fig)
    print("f10 written")

# ---------- F11: T7 gauges as fourth observer ----------
if os.path.exists("data/t7_results.json"):
    T7 = json.load(open("data/t7_results.json"))
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.9), gridspec_kw={"wspace": 0.35})
    for ax, stratum in [(axs[0], "weak"), (axs[1], "strong")]:
        a = T7[stratum]["cy"]["point"]
        groups = [("radar, GOES IR, MW", [("R", a["RGM_nf_R"]), ("GOES IR", a["RGM_nf_IR"]), ("MW", a["RGM_nf_MW"])]),
                  ("radar, MW, gauges", [("R", a["RMG_nf_R"]), ("MW", a["RMG_nf_MW"]), ("gauges", a["RMG_nf_G"])])]
        col = {"R": SLATE, "GOES IR": SAND, "MW": MOSS, "gauges": RUST}
        x0 = 0
        for gname, items in groups:
            for k, (nm, v) in enumerate(items):
                ax.bar(x0 + k, v, 0.8, color=col[nm]); ax.text(x0 + k, v + 0.01, nm, ha="center", fontsize=7, rotation=0)
            ax.text(x0 + 1, -0.07, gname, ha="center", fontsize=8)
            x0 += 4
        ax.set_xticks([]); ax.set_ylim(-0.1, 0.6); ax.axhline(0, color=INK, lw=0.6)
        ax.set_ylabel("noise fraction, centroid y (point support)")
        ax.set_title(f"{'A' if stratum == 'weak' else 'B'} · {stratum} stratum (n={T7[stratum]['cy']['n']:,})")
    ax = axs[2]
    xs = np.arange(2)
    for k, (stratum, colr) in enumerate([("weak", SAND), ("strong", SLATE)]):
        v = [T7[stratum][q]["point"]["c_IR_MW"] for q in ["cy", "cx"]]; c = [T7[stratum][q]["ci95"]["c_IR_MW"] for q in ["cy", "cx"]]
        ax.bar(xs + (k - 0.5) * 0.36, v, 0.34, color=colr, label=stratum)
        ax.errorbar(xs + (k - 0.5) * 0.36, v, yerr=[[a_ - b[0] for a_, b in zip(v, c)], [b[1] - a_ for a_, b in zip(v, c)]], fmt="none", ecolor=INK, lw=0.8, capsize=2)
    ax.set_xticks(xs); ax.set_xticklabels(["centroid y", "centroid x"]); ax.axhline(0, color=INK, lw=0.6)
    ax.set_ylabel("GOES IR–MW error correlation\n(instruments: radar, gauges)"); ax.legend(fontsize=7, frameon=False)
    ax.set_title("C · the shared satellite error, measured")
    fig.savefig("figures/f11_t7.png"); plt.close(fig)
    print("f11 written")

# ---------- F12: T8 anvil fraction ----------
if os.path.exists("data/t8_results.json"):
    T8 = json.load(open("data/t8_results.json")); T8A = json.load(open("data/t8_abs.json"))
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.9), gridspec_kw={"wspace": 0.35})
    am = T8["anvil_mean_by_tercile"]
    for ax, key, title, ylab in [(axs[0], "c_IR_MW", "A · coupling of satellite errors: flat", "GOES–MW error correlation"),
                                 (axs[1], "RGM_diff_R_MW", "B · how badly the triplet is fooled: grows", "nf_radar − nf_MW (R, GOES, MW)")]:
        for v, col in [("cy", SLATE), ("cx", MOSS)]:
            y = [T8[v]["point"][f"{key}_{q}"] for q in range(3)]; c = [T8[v]["ci95"][f"{key}_{q}"] for q in range(3)]
            ax.errorbar(am, y, yerr=[[a - b[0] for a, b in zip(y, c)], [b[1] - a for a, b in zip(y, c)]], fmt="o-", color=col, capsize=2,
                        label={"cy": "centroid y", "cx": "centroid x"}[v])
        ax.axhline(0, color=INK, lw=0.6); ax.set_xlabel("anvil fraction (tercile mean)"); ax.set_ylabel(ylab)
        ax.legend(fontsize=7, frameon=False); ax.set_title(title)
    ax = axs[2]
    for key, col, nm in [("shared_cov_GIR_MW", RUST, "shared GOES–MW cov."), ("errvar_GIR", SAND, "GOES error var."),
                         ("errvar_MW", MOSS, "MW error var."), ("errvar_G", INK, "gauge error var."), ("errvar_R", SLATE, "radar error var.")]:
        ax.plot(am, [T8A["cy"][str(q)][key] for q in range(3)], "o-", color=col, label=nm)
    ax.set_yscale("log"); ax.set_xlabel("anvil fraction (tercile mean)"); ax.set_ylabel("cells² (centroid y)")
    ax.legend(fontsize=6.5, frameon=False); ax.set_title("C · everything grows together (exploratory)")
    fig.savefig("figures/f12_t8.png"); plt.close(fig)
    print("f12 written")

# ---------- F13: the collocation thread in one picture ----------
if os.path.exists("data/t7_results.json"):
    A_ = json.load(open("data/probe_a_results.json")); A2_ = json.load(open("data/probe_a2_results.json"))
    T5_ = json.load(open("data/t5_results.json")); T6_ = json.load(open("data/t6_results.json")); T7_ = json.load(open("data/t7_results.json"))
    T1_ = json.load(open("data/t1_results.json"))
    def gap(r, v="cy"):
        return r[v]["point"]["nf_R"] - r[v]["point"]["nf_MW"]
    steps = [
        ("Probe A\nR, IMERG IR, MW", gap(A_), None, SAND),
        ("Probe A2\nR, GOES IR, MW", gap(A2_), None, SAND),
        ("T5\nradar blurred to\nMW smoothness", T5_["real"][str(T5_["sigma_star"])]["cy"]["diff_R_minus_MW"], T5_["real"][str(T5_["sigma_star"])]["cy"]["diff_ci95"], SAND),
        ("T6\ncores only", T6_["core"]["cy"]["diff_R_minus_MW"], T6_["core"]["cy"]["diff_ci95"], SAND),
        ("T7\nR, GOES, MW\nat gauge points", T7_["weak"]["cy"]["point"]["RGM_diff_R_MW"], T7_["weak"]["cy"]["ci95"]["RGM_diff_R_MW"], SAND),
        ("T7\nR, MW, gauges\n(sharp, independent)", T7_["weak"]["cy"]["point"]["RMG_diff_R_MW"], T7_["weak"]["cy"]["ci95"]["RMG_diff_R_MW"], MOSS),
    ]
    fig, axs = plt.subplots(1, 2, figsize=(13.5, 4.2), gridspec_kw={"width_ratios": [2.2, 1], "wspace": 0.3})
    ax = axs[0]
    for k, (lab_, v, c, col) in enumerate(steps):
        ax.bar(k, v, 0.65, color=col)
        if c:
            ax.errorbar(k, v, yerr=[[v - c[0]], [c[1] - v]], fmt="none", ecolor=INK, lw=0.8, capsize=2)
    ax.axhline(0, color=INK, lw=0.7)
    ax.set_xticks(range(len(steps))); ax.set_xticklabels([s[0] for s in steps], fontsize=7.5)
    ax.set_ylabel("nf_radar − nf_microwave (centroid y)\n> 0: collocation rates radar worse")
    ax.text(0.2, 0.26, "satellite pair as two legs: radar looks worse", fontsize=8, color=RUST)
    ax.text(4.55, -0.15, "gauges as a leg:\nradar best", fontsize=8, color=MOSS, ha="right")
    ax.set_ylim(-0.3, 0.3); ax.set_title("A · the inversion through the experiments")
    ax = axs[1]
    TBx = json.load(open("data/t1bis_results.json"))
    pts = [(0.0, TBx["S_D"]["full"]["cy"]), (0.5, TBx["S_E"]["0.5"]["full"]), (1.0, TBx["S_E"]["1.0"]["full"])]
    ax.plot([p_[0] for p_ in pts], [p_[1]["est"]["R"] - p_[1]["est"]["S1"] for p_ in pts], "o-", color=SAND, label="collocation estimate")
    ax.plot([p_[0] for p_ in pts], [p_[1]["true"]["R"] - p_[1]["true"]["S1"] for p_ in pts], "o--", color=MOSS, label="truth")
    ax.axhline(0, color=INK, lw=0.7); ax.set_xlabel("planted shared satellite error σ_sh (cells)")
    ax.set_ylabel("nf_R − nf_S1 (centroid y)"); ax.legend(fontsize=7, frameon=False)
    ax.set_title("B · known truth (corrected control): a shared error does this")
    fig.savefig("figures/f13_thread.png"); plt.close(fig)
    print("f13 written")

# ---------- F14: T1-bis and T7-S ----------
if os.path.exists("data/t1bis_results.json") and os.path.exists("data/t7s_results.json"):
    TB = json.load(open("data/t1bis_results.json")); TS = json.load(open("data/t7s_results.json"))
    fig, axs = plt.subplots(1, 3, figsize=(14.5, 4.1), gridspec_kw={"wspace": 0.4})
    ax = axs[0]
    regs = ["full", "win", "win_s1", "win_s2"]; rl = ["full canvas", "window", "window\nshift 1", "window\nshift 2"]
    for k, (v, col, nm) in enumerate([("cy", SLATE, "centroid y"), ("cx", MOSS, "centroid x"), ("d", SAND, "dispersion d")]):
        ax.bar(np.arange(4) + (k - 1) * 0.26, [TB["S_BC"][r][v]["S1-S2"] for r in regs], 0.25, color=col, label=nm)
    ax.axhline(0, color=INK, lw=0.6); ax.axhline(0.16, color=RUST, ls="--", lw=0.9)
    ax.text(1.5, 0.152, "T1 value (blur not mass-conserving)", color=RUST, fontsize=7, ha="center", va="top"); ax.set_ylim(-0.01, 0.18)
    ax.set_xticks(range(4)); ax.set_xticklabels(rl, fontsize=8); ax.set_ylabel("true S1–S2 coordinate-error corr.")
    ax.legend(fontsize=7, frameon=False, loc="center left"); ax.set_title("A · T1-bis: no planted error")
    ax = axs[1]
    for sig, col in [("0.5", SAND), ("1.0", RUST)]:
        e = TB["S_E"][sig]["full"]
        ax.plot(["R", "S1", "S2"], [e["est"][o] for o in ["R", "S1", "S2"]], "o-", color=col, label=f"σ_sh={sig} estimate")
        ax.plot(["R", "S1", "S2"], [e["true"][o] for o in ["R", "S1", "S2"]], "o--", color=col, alpha=0.5, label=f"σ_sh={sig} truth")
    ax.set_ylabel("noise fraction, centroid y (full canvas)"); ax.legend(fontsize=7, frameon=False)
    ax.set_title("B · T1-bis: planted shared error")
    ax = axs[2]
    g = [float(x) for x in TS["weak_cy"]["grid"]]
    for case, col, nm in [("weak_cy", SLATE, "weak, y"), ("weak_cx", MOSS, "weak, x"), ("strong_cy", SAND, "strong, y"), ("strong_cx", RUST, "strong, x")]:
        gr = TS[case]["grid"]; ks = sorted(gr, key=float)
        ax.plot([float(k) for k in ks], [gr[k]["rho_IM"] for k in ks], "-", color=col, label=f"{nm}: GOES–MW error corr.")
        ax.plot([float(k) for k in ks], [gr[k]["diff_R_MW"] for k in ks], ":", color=col)
    ax.axhline(0, color=INK, lw=0.6); ax.set_xlabel("assumed radar–gauge error correlation")
    ax.set_ylabel("solid: GOES–MW corr.; dotted: nf_R − nf_MW"); ax.legend(fontsize=6.5, frameon=False)
    ax.set_title("C · T7-S: radar–gauge sensitivity")
    fig.savefig("figures/f14_t1bis_t7s.png"); plt.close(fig)
    print("f14 written")
