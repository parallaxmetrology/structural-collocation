"""Build LABBOOK.html: self-contained (inline CSS, base64 figures), numbers read from data/*.json.
Run from precip-disagreement/:  .venv/bin/python scripts/make_figures.py && .venv/bin/python build_labbook.py
"""
import base64, hashlib, json, pathlib, platform, datetime
import numpy as np, scipy, pandas, h5py, matplotlib

H = pathlib.Path(".")
J = lambda p: json.load(open(H / p))
F = J("data/results.json"); A = J("data/probe_a_results.json"); A2 = J("data/probe_a2_results.json")
PH = J("data/probe_a_posthoc.json")
T1 = J("data/t1_results.json"); TX = J("data/t1_explore.json")
T5 = J("data/t5_results.json"); OFF = J("data/t5_offsets.json")
T6 = J("data/t6_results.json"); TI = J("data/t6_intensity.json"); TQ = J("data/t6_rqi.json")
T7 = J("data/t7_results.json"); T7R = J("data/t7r_results.json")
TB = J("data/t1bis_results.json"); TS = J("data/t7s_results.json"); T9 = J("data/t9_results.json")
X10 = J("data_aug/t10_score.json"); X10J = J("data/t10_score_july.json"); TSA = J("data_aug/t7s_results.json")
T8 = J("data/t8_results.json"); T8A = J("data/t8_abs.json")


def fig(name):
    return "data:image/png;base64," + base64.b64encode((H / "figures" / name).read_bytes()).decode()


def sha(p):
    return hashlib.sha256((H / p).read_bytes()).hexdigest()[:16]


def f3(x): return f"{x:.3f}"
def ci(c): return f"[{c[0]:.3f}, {c[1]:.3f}]"


VARS = ["wet", "cy", "cx", "disp", "isl", "logmean"]
VN = {"wet": "wet fraction", "cy": "centroid y", "cx": "centroid x", "disp": "dispersion", "isl": "islands",
      "logmean": "log mean rain (value)"}


def tc_rows(R, irname):
    out = []
    for v in VARS:
        p, c = R[v]["point"], R[v]["ci95"]
        out.append(f"<tr><td>{VN[v]}</td><td>{f3(p['nf_R'])} {ci(c['nf_R'])}</td><td>{f3(p['nf_IR'])} {ci(c['nf_IR'])}</td>"
                   f"<td>{f3(p['nf_MW'])} {ci(c['nf_MW'])}</td><td>{R[v]['share_boot_negative_ev']:.2f}</td>"
                   f"<td>{'no' if R[v]['unstable'] else 'yes'}</td></tr>")
    return "\n".join(out)


def quad_rows(R):
    out = []
    for v in ["cy", "cx", "disp", "isl"]:
        p, c = R[v]["point"], R[v]["ci95"]
        mw = "n/a" if p["M_MW"] != p["M_MW"] else f"{p['M_MW']:.2f} [{c['M_MW'][0]:.2f}, {c['M_MW'][1]:.2f}]"
        out.append(f"<tr><td>{VN[v]}</td><td>{p['M_IR']:.2f} [{c['M_IR'][0]:.2f}, {c['M_IR'][1]:.2f}]</td><td>{mw}</td></tr>")
    return "\n".join(out)


def t1_rows():
    s0 = T1["0.0"]; out = []
    for v in ["cy", "cx", "disp", "isl", "wet", "logmean"]:
        cells = []
        for slot in ["R", "MW", "IR"]:
            e = s0["tc"][v]["point"][f"nf_{slot}"]; tr = s0["truth"][v]["nf_true"][slot]
            cells.append(f"{e:.3f} / {tr:.3f}")
        out.append(f"<tr><td>{VN[v]}</td><td>{cells[0]}</td><td>{cells[1]}</td><td>{cells[2]}</td>"
                   f"<td>{s0['tc'][v]['share_boot_negative_ev']:.2f}</td></tr>")
    return "\n".join(out)


def t1_sweep_rows():
    out = []
    for s_ in ["0.0", "0.5", "1.0", "1.5"]:
        p = T1[s_]["tc"]["cy"]["point"]; tr = T1[s_]["truth"]["cy"]
        out.append(f"<tr><td>{s_}</td><td>{p['nf_R']:.3f} / {tr['nf_true']['R']:.3f}</td><td>{p['nf_MW']:.3f} / {tr['nf_true']['MW']:.3f}</td>"
                   f"<td>{p['nf_IR']:.3f} / {tr['nf_true']['IR']:.3f}</td><td>{T1[s_]['quadG']['cy']['point']:.2f} / {tr['errcorr_true']['MW_IR']:.2f}</td></tr>")
    return "\n".join(out)


def t5_rows():
    out = []
    for l in sorted(T5["real"], key=float):
        r = T5["real"][l]; star = " (σ*)" if float(l) == T5["sigma_star"] else ""
        cells = "".join(f"<td>{r[v]['diff_R_minus_MW']:+.3f} {ci(r[v]['diff_ci95'])}</td>" for v in ["cy", "cx", "disp"])
        out.append(f"<tr><td>{l}{star}</td><td>{r['n']:,}</td>{cells}</tr>")
    return "\n".join(out)


def off_rows():
    return "\n".join(f"<tr><td>{b}</td><td>{int(v['n']):,}</td><td>{v['GIR_minus_R']:+.2f}</td><td>{v['MW_minus_R']:+.2f}</td></tr>"
                     for b, v in OFF["by_latitude_cy"].items())


def t6_int_rows():
    return "\n".join(f"<tr><td>≤ {TI['edges_mm_h'][q]:.1f}</td><td>{TI[str(q)]['n']:,}</td><td>{TI[str(q)]['cy']['diff']:+.3f} {ci(TI[str(q)]['cy']['ci'])}</td>"
                     f"<td>{TI[str(q)]['cy']['nf_R']:.3f}</td><td>{TI[str(q)]['cy']['nf_MW']:.3f}</td></tr>" for q in range(5))


def t6_rqi_rows():
    out = []
    for a in range(3):
        w, s = TQ[f"rqi{a}_peak0"], TQ[f"rqi{a}_peak1"]
        out.append(f"<tr><td>{w['rqi_mean']:.2f}</td><td>{w['cy']['diff']:+.3f} {ci(w['cy']['ci'])}</td><td>{s['cy']['diff']:+.3f} {ci(s['cy']['ci'])}</td></tr>")
    return "\n".join(out)


def t7_rows():
    out = []
    for v, vn in [("cy", "centroid y"), ("cx", "centroid x")]:
        for s in ["weak", "strong"]:
            a = T7[s][v]; p = a["point"]; c = a["ci95"]
            out.append(f"<tr><td>{vn}</td><td>{s} (n={a['n']:,})</td><td>{p['RGM_diff_R_MW']:+.3f} {ci(c['RGM_diff_R_MW'])}</td>"
                       f"<td>{p['RMG_diff_R_MW']:+.3f} {ci(c['RMG_diff_R_MW'])}</td><td>{p['c_IR_MW']:+.3f} {ci(c['c_IR_MW'])}</td></tr>")
    return "\n".join(out)


def t8_rows():
    out = []
    for q in range(3):
        a = T8["cy"]["point"]; c = T8["cy"]["ci95"]; ab = T8A["cy"][str(q)]
        out.append(f"<tr><td>{T8['anvil_mean_by_tercile'][q]:.2f}</td><td>{a[f'c_IR_MW_{q}']:+.3f} {ci(c[f'c_IR_MW_{q}'])}</td>"
                   f"<td>{a[f'RGM_diff_R_MW_{q}']:+.3f} {ci(c[f'RGM_diff_R_MW_{q}'])}</td><td>{a[f'RMG_nf_R_{q}']:.3f} / {a[f'RMG_nf_MW_{q}']:.3f} / {a[f'RMG_nf_G_{q}']:.3f}</td>"
                   f"<td>{ab['shared_cov_GIR_MW']:.2f}</td><td>{ab['errvar_GIR']:.2f} / {ab['errvar_MW']:.2f} / {ab['errvar_G']:.2f}</td></tr>")
    return "\n".join(out)


def t7r_rows():
    out = []
    for key, lab, pred in [("R1_R_minus_IR", "R1: nf_R − nf_GOES in (R, GOES, gauges)", "&lt; 0"), ("R1_R_minus_G", "R1: nf_R − nf_gauges in (R, GOES, gauges)", "&lt; 0"),
                           ("R2_nfR_RIRG_minus_RMWG", "R2: nf_R with GOES dropped − with MW dropped", "≈ 0"),
                           ("t3", "R3: tetrad t3 (excludes the GOES–MW pair)", "≈ 0"), ("t1", "R3: tetrad t1", "&gt; 0"), ("t2", "R3: tetrad t2", "&gt; 0")]:
        cells = "".join(f"<td>{T7R[s][v]['point'][key]:+.3f} {ci(T7R[s][v]['ci95'][key])}</td>" for s in ["weak", "strong"] for v in ["cy", "cx"])
        out.append(f"<tr><td>{lab}</td>{cells}<td>{pred}</td></tr>")
    return "\n".join(out)


def t7r_nf_rows():
    p = T7R["weak"]["cy"]["point"]; names = {"R": "radar", "IR": "GOES IR", "MW": "microwave", "G": "gauges"}
    tris = [("RIRMW", "radar, GOES, MW"), ("RIRG", "radar, GOES, gauges"), ("RMWG", "radar, MW, gauges"), ("IRMWG", "GOES, MW, gauges")]
    out = []
    for tk, tl in tris:
        cells = "".join(f"<td>{p[f'nf_{o}|{tk}']:.3f}</td>" if f"nf_{o}|{tk}" in p else "<td>—</td>" for o in ["R", "IR", "MW", "G"])
        out.append(f"<tr><td>{tl}</td>{cells}</tr>")
    return "\n".join(out)


def tb_rows():
    out = []
    for r, rl in [("full", "full canvas"), ("win", "central window"), ("win_s1", "window shifted (+3,+3)"), ("win_s2", "window shifted (−3,+4)")]:
        bc = TB["S_BC"][r]; sd = TB["S_D"][r]
        cells = "".join(f"<td>{bc[v]['S1-S2']:+.3f}</td>" for v in ["cy", "cx", "d"])
        cells += "".join(f"<td>{sd[v]['max_abs_err']:.3f} {'✓' if sd[v]['rank_ok'] else '✗'}</td>" for v in ["cy", "d", "d2"])
        out.append(f"<tr><td>{rl}</td>{cells}</tr>")
    return "\n".join(out)


def ts_rows():
    out = []
    for c, cl in [("weak_cy", "weak, centroid y"), ("weak_cx", "weak, centroid x"), ("strong_cy", "strong, centroid y"), ("strong_cx", "strong, centroid x")]:
        r = TS[c]; g0 = r["grid"]["0.0"]
        out.append(f"<tr><td>{cl}</td><td>{g0['rho_IM']:+.3f} {ci(g0['rho_IM_ci'])}</td><td>{r['rho_star_1'] if r['rho_star_1'] is not None else 'not reached (≤ 0.5)'}</td>"
                   f"<td>{r['rho_star_2']}</td><td>≥ {r['rho_max']}</td></tr>")
    return "\n".join(out)


dec = F["decision"]; nul = F["null_spatial_scramble"]
DEC_ROWS = "\n".join(
    f"<tr><td>{lab}</td><td>{dec[k]['rho_full']:.4f}</td><td>{dec[k]['delta']:+.4f}</td><td>{ci(dec[k]['ci95'])}</td>"
    f"<td>{nul[k]['delta']:+.4f} {ci(nul[k]['ci95'])}</td><td>{'met' if dec[k]['go'] else 'not met'}</td></tr>"
    for k, lab in [("frozen_v0", "frozen v0 (om_loc, om_pat)"), ("variance_only", "variance only (σ space, σ time)"),
                   ("FULL_A2", "<b>FULL A2 (14 reads), headline</b>")])
DD = F["descriptive_discovery"]
DESC_ROWS = "\n".join(
    f"<tr><td><code>{k}</code></td><td>{DD['ALL'][k]['rho']:+.3f}</td><td>{DD['ALL'][k]['partial']:+.3f}</td>"
    f"<td>{DD['USCRN'][k]['partial']:+.3f}</td><td>{DD['ALL'][k]['n']:,}</td></tr>"
    for k in ["inc", "naive1", "naive2", "om_loc", "om_pat", "om_iou", "om_islands", "om_disp", "om_void", "om_scale",
              "om_scale_rec", "om_edge", "om_pat_persist", "om_loc_ir", "om_pat_ir", "om_sigma_space", "om_sigma_time",
              "om_loc_mw", "om_pat_mw"])
RF = F["resolution_floor_discovery"]
G6 = F["gate6"]; PA = F["predictor_audit"]

SCRIPTS = ["scripts/fetch_gauges.py", "scripts/fetch_mrms.py", "scripts/fetch_imerg.py", "scripts/fetch_imerg_extra.py",
           "scripts/fetch_goes.py", "scripts/analyze.py", "scripts/probe_a.py", "scripts/probe_a_posthoc.py",
           "scripts/probe_a2.py", "scripts/t1_synthetic.py", "scripts/t1_explore.py", "scripts/t5.py", "scripts/t5_offsets.py", "scripts/t6.py", "scripts/t6_intensity.py", "scripts/t6_rqi.py", "scripts/t7.py", "scripts/t7_robust.py", "scripts/t8.py", "scripts/t8_abs.py", "scripts/t1bis.py", "scripts/t7s.py", "scripts/t9.py", "scripts/period.py", "scripts/make_stations.py", "scripts/fetch_imerg_both.py", "scripts/t10_area.py", "scripts/t10_score.py", "scripts/make_figures.py", "build_labbook.py"]
DOCS = ["FEASIBILITY_PLAN.md", "FEASIBILITY_RESULTS.md", "RESEARCH_SYNTHESIS.md", "PROBE_A_PLAN.md",
        "PROBE_A_RESULTS.md", "PROBE_A2_PLAN.md", "PROBE_A2_RESULTS.md", "T1_PLAN.md", "T1_RESULTS.md", "T5_PLAN.md", "T5_RESULTS.md", "T6_PLAN.md", "T6_RESULTS.md", "T7_PLAN.md", "T7_RESULTS.md", "T7R_PLAN.md", "T8_PLAN.md", "T8_RESULTS.md", "T1BIS_PLAN.md", "T1BIS_RESULTS.md", "T7S_PLAN.md", "T7S_RESULTS.md", "T9_PLAN.md", "T9_RESULTS.md", "T10_PLAN.md", "T10_RESULTS.md", "GATE.md", "literature/INDEX.md"]
RESULTS = ["data/results.json", "data/probe_a_results.json", "data/probe_a_posthoc.json", "data/probe_a2_results.json", "data/t1_results.json", "data/t1_explore.json", "data/t5_results.json", "data/t5_offsets.json", "data/t6_results.json", "data/t6_intensity.json", "data/t6_rqi.json", "data/t7_results.json", "data/t7r_results.json", "data/t8_results.json", "data/t8_abs.json", "data/t1bis_results.json", "data/t7s_results.json", "data/t9_results.json", "data/t10_area_july.json", "data/t10_score_july.json", "data_aug/t10_area.json", "data_aug/t10_score.json", "data_aug/t7_results.json", "data_aug/t7r_results.json", "data_aug/t7s_results.json", "data_aug/probe_a2_results.json"]
def dir_digest(d):
    files = sorted((H / d).glob("*"))
    h = hashlib.sha256()
    for f in files:
        h.update(f.name.encode()); h.update(f.read_bytes())
    return len(files), sum(f.stat().st_size for f in files), h.hexdigest()[:16]


DATA_ROWS = "\n".join(f"<tr><td><code>{d}/</code></td><td>{n:,}</td><td>{b/1e6:,.0f} MB</td><td><code>{dg}</code></td></tr>"
                      for d in ["data/mrms", "data/imerg", "data/imerg_x", "data/goes", "data_aug/mrms", "data_aug/imerg", "data_aug/imerg_x", "data_aug/goes"] for n, b, dg in [dir_digest(d)])
DATA_ROWS += "\n" + "\n".join(f"<tr><td><code>{f}</code></td><td>1</td><td>{(H/f).stat().st_size/1e6:,.1f} MB</td><td><code>{sha(f)}</code></td></tr>"
                               for f in ["data/gauges.csv", "data/stations.csv", "data_aug/gauges.csv", "data_aug/stations.csv"])
HASH_ROWS = "\n".join(f"<tr><td><code>{p}</code></td><td><code>{sha(p)}</code></td></tr>" for p in SCRIPTS + DOCS + RESULTS)

CSS = """
body{max-width:940px;margin:40px auto;padding:0 22px;font-family:Georgia,'Times New Roman',serif;color:#1a1a1a;line-height:1.55;background:#fbfaf7}
h1{font-size:27px;margin-bottom:2px} h2{font-size:19px;margin-top:38px;border-bottom:2px solid #7a3b2e;padding-bottom:4px}
h3{font-size:15px;margin-top:22px;color:#7a3b2e}
.masthead{color:#6b6b6b;font-size:14px;margin-bottom:6px}
.find{background:#eef3ec;border-left:5px solid #3f6b3a;padding:12px 15px;font-size:14px;margin:14px 0}
.status{background:#f3efe6;border-left:4px solid #7a3b2e;padding:10px 14px;font-size:14px;margin:14px 0}
.kill{background:#f6ecea;border-left:5px solid #9c3b2a;padding:12px 15px;font-size:14px;margin:14px 0}
.open{background:#ecf0f5;border-left:5px solid #4a6280;padding:12px 15px;font-size:14px;margin:14px 0}
table{border-collapse:collapse;width:100%;margin:14px 0;font-size:13px;font-family:-apple-system,system-ui,sans-serif}
th,td{border:1px solid #e6e2da;padding:5px 8px;text-align:left;vertical-align:top} th{background:#f3efe6}
td{font-variant-numeric:tabular-nums}
img{max-width:100%;margin:10px 0;border:1px solid #e6e2da;background:#fff}
code{background:#f0ece3;padding:1px 4px;border-radius:3px;font-size:12.5px}
pre{background:#f0ece3;padding:10px 12px;font-size:12px;overflow-x:auto;border-radius:3px}
.note{font-size:13px;color:#6b6b6b} figcaption{font-size:12px;color:#6b6b6b;margin-bottom:16px;font-family:-apple-system,system-ui,sans-serif}
.tag{font-size:11px;font-weight:bold;padding:1px 5px;border-radius:3px;font-family:-apple-system,system-ui,sans-serif}
.conf{background:#dbe7d6;color:#2f5b2a} .expl{background:#e7e0cf;color:#7a3b2e} .hyp{background:#dfe5ee;color:#4a6280}
ol.trail li{margin-bottom:5px;font-size:14px} .toc{font-size:13.5px;font-family:-apple-system,system-ui,sans-serif}
"""

T1ROWS = t1_rows(); T1SWEEP = t1_sweep_rows(); FIG8 = fig("f8_t1.png")
T5ROWS = t5_rows(); OFFROWS = off_rows(); FIG9 = fig("f9_t5.png")
T6INT = t6_int_rows(); T6RQI = t6_rqi_rows(); FIG10 = fig("f10_t6.png")
T7ROWS = t7_rows(); FIG11 = fig("f11_t7.png"); T7RROWS = t7r_rows(); T7RNF = t7r_nf_rows()
T8ROWS = t8_rows(); FIG12 = fig("f12_t8.png"); FIG13 = fig("f13_thread.png")
TBROWS = tb_rows(); TSROWS = ts_rows(); FIG14 = fig("f14_t1bis_t7s.png")
HTML = f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Precipitation Disagreement: Lab Book</title>
<meta name="description" content="Lab book for the precipitation disagreement study: a negative domain-qualification result and the structural collocation thread it opened, replicated on unseen data.">
<link rel="canonical" href="https://www.parallaxmetrology.com/Precipitation_LABBOOK.html"><style>{CSS}</style></head><body>

<div class="masthead">Parallax Domain Qualification · Lab Book · 2026-10-15 ·
data: MRMS radar-only QPE, IMERG Late V07 (merged, IR, MW), GOES-16 ABI band 13, USCRN + ASOS gauges ·
July 2024 (replication: August 2024), CONUS · built {datetime.date.today().isoformat()}</div>
<h1>Precipitation: does structural disagreement among independent observers carry information?</h1>
<div class="note">An investigation that started as a domain-qualification test and followed the evidence.
Read each claim as <span class="tag conf">confirmatory</span> (pre-registered before measuring),
<span class="tag expl">exploratory</span> (found after) or <span class="tag hyp">hypothesis</span> (proposed, untested).
Every claim carries result + sample + where-it-fails.
Code, plans and saved results: <a href="https://github.com/parallaxmetrology/structural-collocation">github.com/parallaxmetrology/structural-collocation</a>.</div>

<div class="status"><b>Status (2026-10-15, paused after T10).</b> The original question is closed, negative. The thread it opened,
structural collocation, has one robust empirical result: the observer ranking depends on which observers form the triplet, and
replacing either satellite with gauges reverses it. <b>It replicated on an unseen period</b> (August 2024, frozen pipeline,
all six pre-registered predictions met; T10, §12). Reading the four-observer pattern as an error shared by the two satellites is
conditional on radar and gauge errors being independent, under assumed radar–gauge error correlations, satellite dependence detected up to 0.475–0.5 (July) and through the tested maximum of 0.5 (August); radar's ranking clear (interval excluding reversal) up to 0.275–0.4 (July), 0.375–0.4 (August) (T7-S); the actual value is not measured.
After an external review, the synthetic control was rebuilt with a mass-conserving blur (T1-bis): structural collocation recovers
centroid and dispersion with correct rankings, a planted shared satellite error reproduces the inversion, and five observers
identify a shared pair in simulation (under its assumptions); the earlier "smoothing creates shared error" (F6) is withdrawn. Decision gate reached
(<code>GATE.md</code>); T9 and the T10 replication run since; paused for discussion. Start with the summary; experiments follow in order.</div>

<div class="kill"><b>Review corrections (2026-10-15).</b> An external review (GPT, reading the lab book and preprint without the code),
checked here against the code, found: (1) the T1 dispersion estimates rank the observers wrongly with nothing planted
(radar-like {T1['0.0']['tc']['disp']['point']['nf_R']:.3f} estimated vs {T1['0.0']['truth']['disp']['nf_true']['R']:.3f} true; microwave-like {T1['0.0']['tc']['disp']['point']['nf_MW']:.3f} vs {T1['0.0']['truth']['disp']['nf_true']['MW']:.3f}), so dispersion is demoted; (2) "all assumptions met" was wrong
for T1's σ_sh = 0 case, since the extracted coordinates' errors were already correlated; (3) the simulated blur widened cells with
fixed peaks, which does not conserve rain mass, added widths linearly, and was read inside clipping windows, so F6 and the baseline
shared error may be rendering or windowing artefacts; (4) the three tetrads satisfy t2 − t1 = t3, so they are two independent
constraints, and an interval containing zero shows compatibility, not "no other breach"; (5) "constant" coupling across anvil
terciles should read "no detectable change"; (6) T5's smoothness match was approximate (radar 0.759 at the chosen step vs
microwave 0.700) and its sample changed along the ladder; (7) the satellite attribution and the +0.35 to +0.51 correlations are
conditional on radar–gauge error independence. The text below has been corrected; superseded claims are marked.</div>

<div class="toc"><b>Contents</b> · Summary · 1 Why precipitation · 2 Observers and data · 3 Experiment 1: radar error at gauges ·
4 Literature step · 5 Experiment 2: structural collocation (Probe A) · 6 Experiment 3: raw GOES infrared (Probe A2) ·
7 Experiment 4: T1 synthetic anchor · 8 Experiment 5: T5 smoothness vs physics · 9 Experiment 6: T6 core vs area · 10 Experiment 7: T7 gauges as fourth observer · 11 Experiment 8: T8 anvil fraction · 12 Corrected controls and replication: T1-bis, T7-S, T9, T10 · 13 What the evidence shows · 14 What it does not show · 15 Open threads · 16 Chronological trail · 17 Reproducibility</div>

<h2>Summary of the collocation thread</h2>
<p><b>What was asked.</b> Can disagreement among physically independent observers of one rain field be turned into a measurement
of each observer's error, with no ground truth? The mathematics exists for amounts (triple collocation: each observer's error
variance from three observers' covariances). This thread applied it to <i>structural</i> coordinates (where the rain sits, how
spread it is) and followed what came out.</p>
<p><b>What happened, in one line per step.</b></p>
<table><tr><th>step</th><th>question</th><th>result</th><th>what it changed</th></tr>
<tr><td>Probe A</td><td>collocate radar, IMERG infrared, microwave</td><td>stable on structure; ranks microwave above radar</td><td>an unexpected ranking</td></tr>
<tr><td>Probe A2</td><td>replace IMERG infrared (microwave-calibrated) with raw GOES</td><td>same ranking</td><td>calibration link not supported as the cause</td></tr>
<tr><td>T1</td><td>synthetic fields with known truth</td><td>centroid recovered; dispersion mis-ranked; a planted shared satellite error reproduces the ranking (superseded by T1-bis: blur not mass-conserving)</td><td>the ranking reread as a possible shared-error signature; islands and wet fraction dropped</td></tr>
<tr><td>T5</td><td>blur radar to microwave smoothness</td><td>ranking persists under the tested smoothing (approximate match, sample changes)</td><td>smoothing procedure not sufficient; GOES parallax detected, not shared</td></tr>
<tr><td>T6</td><td>read cores only; split by intensity and radar beam quality</td><td>not core keying; worst in light rain; not radar beam quality</td><td>where the error lives</td></tr>
<tr><td>T7</td><td>gauges as a sharp, independent fourth observer</td><td><b>with gauges as a leg, the ranking reverses (radar least noisy); GOES–microwave error correlation +0.35 to +0.51 if radar and gauge errors are independent</b></td><td>the reversal established; the satellite attribution conditional</td></tr>
<tr><td>T7-R</td><td>robustness: drop the other satellite; four-observer tetrads</td><td>radar still least noisy; tetrads show two groups, {{radar, gauges}} vs {{GOES, MW}}, compatible with one split (two independent tetrad constraints)</td><td>T7 robust to the choice of satellite</td></tr>
<tr><td>T8</td><td>does the coupling grow with anvil (cloud beyond rain)?</td><td>no detectable change in coupling (≈ 0.4); triplet bias grows where every observer is noisier</td><td>a simple "anvil" cause not supported</td></tr>
<tr><td>T1-bis</td><td>rebuild the synthetic control with a mass-conserving blur; measure coordinate-error independence directly</td><td>full canvas: centroid errors uncorrelated, dispersion small residual (+0.021); windows add small shared error; all variables recovered with correct ranking; planted shared error inverts at σ_sh = 1.0; five observers identify the pair 10/10, no false alarms</td><td>F6 withdrawn; dispersion reinstated; synthetic support restored</td></tr>
<tr><td>T9</td><td>real-data window grid; smoothness on a fixed sample with a close match (lag-1 within 0.013)</td><td>shifted and larger windows: inversion and point reversal intact; close smoothness match (σ* 0.4): gap +0.218, unchanged</td><td>robust to the tested grids and the smoothness match; windowing in general not excluded</td></tr>
<tr><td>T10</td><td>does it replicate on an unseen period? August 2024, frozen July code</td><td>all six predictions met: area inversion {X10['rows'][0]['point']:+.3f}, point reversal {X10['rows'][4]['point']:+.3f}, c(GOES, MW) {X10['rows'][12]['point']:+.3f}; T7-S does not tip by 0.5</td><td>effect not specific to July; identification limit unchanged</td></tr>
<tr><td>T7-S</td><td>how much radar–gauge error dependence would undo the attribution?</td><td>dependence detected to 0.475–0.5; ranking clear (by interval) to 0.275–0.4, point flip 0.425–0.5 (August: dependence through 0.5, ranking clear to 0.375–0.4, no point flip by 0.5)</td><td>the conditional attribution has a quantified tolerance</td></tr></table>
<img src="{FIG13}" alt="thread">
<figcaption>A: radar minus microwave noise fraction (centroid y) through the experiments. Whenever the two satellites are two of the
three legs, collocation rates radar worse; swapping one satellite for gauges reverses it. B: on synthetic fields with known truth, a
planted shared satellite error produces this reversal (original T1; confirmed at σ_sh = 1.0 by the corrected T1-bis).</figcaption>
<p><b>Where knowledge stands.</b></p>
<ul>
<li><b>Established</b> (pre-registered; replicated on an unseen month, T10): the ranking structural collocation gives depends on
which observers form the triplet; with radar, GOES and microwave, microwave is rated least noisy; replacing either satellite with gauges,
at the same point support, reverses it (radar 0.06 &lt; gauges 0.22 &lt; microwave 0.26, centroid y), and radar's score does not depend
on which satellite is replaced. August 2024: radar 0.04 &lt; gauges 0.24 &lt; microwave 0.33 (weak stratum).</li>
<li><b>Conditional</b>: reading these scores as true error rankings, and attributing the four-observer pattern ({{radar, gauges}} vs
{{GOES, microwave}}) to an error shared by the satellites (estimated +0.35 to +0.51), both assume radar and gauge errors are not strongly
correlated. Four observers cannot establish that; under assumed radar–gauge error correlations, satellite dependence detected up to 0.475–0.5 (July) and through the tested maximum of 0.5 (August); radar's ranking clear (interval excluding reversal) up to 0.275–0.4 (July), 0.375–0.4 (August) (T7-S); the actual value is not measured.</li>
<li><b>Supported by the corrected synthetic control</b> (T1-bis): centroid and dispersion recovered with correct rankings; a planted shared
satellite error inverts the ranking at σ_sh = 1.0 cell; five observers identify a shared pair without false alarms (under the simulation's assumptions). F6 withdrawn.</li>
<li><b>Not supported</b> as explanations of the inversion (tests did not detect them): IMERG's infrared-to-microwave calibration, the tested smoothing, GOES parallax,
convective-core keying, radar beam quality, and (in its tested form) aloft-only anvil cloud. Parallax and beam quality were
checked exploratorily; the others by pre-registered tests. "Not supported" means not detected within the tested variants.</li>
<li><b>Model failure</b> (both months): tetrads t1 and t2 are above zero in July and August, so one common signal with mutually
independent errors does not describe the four observers; the satellite attribution is one reading of that failure, not the only one.</li>
<li><b>Not presently claimed:</b> a satellite defect; the physical mechanism; an advantage over amount-based or representativeness
analyses (not compared); a diagnostic that improves on available ones (not compared); generality beyond warm-season CONUS.</li>
<li><b>Open:</b> what the two satellites share, if the attribution holds. A coupling with no detectable change across anvil extent fits both observing the cloud system rather than surface rain;
untested.</li></ul>
<p><b>Rules for structural collocation, learned here</b> (for the tested variables, simulation and observing configuration):</p>
<ol>
<li>Prefer geometric coordinates; with a correct (mass-conserving) blur, centroid and dispersion were recovered and ranked correctly;
threshold counts (islands, wet fraction) failed in the original simulation, were not re-tested, and are excluded from interpretation.</li>
<li>A shared error among the legs is invisible to the estimator; it shows up only as a plausible-looking but wrong ranking.</li>
<li>Check independence after feature extraction and windowing, not only from the hardware; simulation details (here, a blur that did not
conserve mass) can create shared errors of their own, and windows add a small one.</li>
<li>Independence is a property of what each observer responds to, not of the hardware (now in RUBRIC.md Gate 6).</li>
<li>The bias from a shared error is largest where all observers are noisiest (light rain, anvil-heavy windows).</li>
<li>Always anchor on synthetic fields with known truth before reading a real-data ranking (house rule I-20).</li></ol>
<p class="note"><b>Terms.</b> <i>Observer</i>: one independent measurement of the same rain field (radar, GOES infrared, microwave,
gauges). <i>Collocation</i>: estimating each observer's error from the covariances of three or more observers, with no truth.
<i>Leg</i>: one of the three observers in a triplet. <i>Instrument</i>: in the four-observer check, the two observers assumed
independent that let one estimate the error correlation of the other two. <i>Noise fraction</i>: 1 − ρ², the share of an
observer's variance not shared with the common signal (0 = perfect, 1 = pure noise). <i>Inversion</i>: collocation rating
microwave less noisy than radar. <i>Point support</i>: every observer read at the same gauge locations.</p>


<h2>1. Why precipitation</h2>
<p>The domain-qualification rubric (v2, <code>domain-qualification/RUBRIC.md</code>) scored four desk candidates. Weather
spatial verification, reframed as <i>observation disagreement</i> (radar vs satellite vs gauge), was the only one with public,
use-unrestricted data large enough for a powered test (scorecard <code>candidates/04_precipitation_spatial_verification.md</code>).
Its weaker gates were stated in advance: the label (a point error is not an arrangement, Gate 4 PARTIAL) and the incumbents
(SAL, MODE, FSS, MRMS quality layers). The first test was built to fail fast on exactly those.</p>

<h2>2. Observers and data</h2>
<table><tr><th>Observer</th><th>What it measures physically</th><th>Source (access)</th><th>Used in</th></tr>
<tr><td>R · MRMS RadarOnly_QPE_01H</td><td>radar backscatter → rain rate, not gauge-adjusted; 0.01° → 0.1° block mean</td><td>NOAA Open Data, AWS <code>noaa-mrms-pds</code> (open)</td><td>all</td></tr>
<tr><td>RQI · RadarAccumulationQualityIndex_01H</td><td>MRMS's own radar quality (beam height, blockage, freezing level)</td><td>same</td><td>Exp. 1 (incumbent)</td></tr>
<tr><td>M · IMERG Late merged</td><td>Kalman-morphed blend of microwave and infrared estimates</td><td>NASA GES DISC (Earthdata token)</td><td>all</td></tr>
<tr><td>IR · IMERG <code>IRprecipitation</code> (PDIR)</td><td>geostationary IR → rain, <b>calibrated to the microwave estimates</b> (ATBD)</td><td>same file</td><td>Exp. 1, Probe A</td></tr>
<tr><td>MW · IMERG <code>MWprecipitation</code></td><td>passive-microwave retrievals, quantile-matched to CORRA (spaceborne radar–radiometer)</td><td>same file</td><td>Probes A, A2</td></tr>
<tr><td>GIR · GOES-16 ABI band 13</td><td>raw 10.3 µm brightness temperature; cold-cloud field 235 − Tb</td><td>NOAA Open Data, AWS <code>noaa-goes16</code> (open)</td><td>Probe A2</td></tr>
<tr><td>G · gauges</td><td>USCRN (130 CONUS, triple weighing gauges) and ASOS (2,175) hourly precipitation</td><td>NCEI; Iowa Environmental Mesonet</td><td>Exp. 1 (label only)</td></tr></table>
<img src="{fig('f1_observers.png')}" alt="one hour through every observer">
<figcaption>One wet hour seen by each observer, on the common 0.1° grid. Radar and microwave place the rain bands similarly;
the infrared estimate and the raw cold-cloud field are broader and displaced (anvils, parallax). Bottom right: the 1.5°
windows that entered Probe A2 this hour.</figcaption>
<img src="{fig('f7_coverage.png')}" alt="coverage">
<figcaption>Where the evidence came from. Left: label gauges for Experiment 1. Right: Probe A2 windows over July 2024;
the GOES CONUS sector thins west of ~120°W, which removes western windows.</figcaption>

<h2>3. Experiment 1: does structural disagreement predict radar error at gauges? <span class="tag conf">confirmatory</span></h2>
<p><b>Design</b> (frozen in <code>FEASIBILITY_PLAN.md</code>, with two amendments written before any target or predictor
value existed). Target: E = |radar-only at the gauge pixel − gauge|, hourly. 504 hours, 1–21 July 2024; discovery 1–10 July,
validation 11–21 July, read once. Baseline: 1 − RQI, two plain radar–satellite difference maps, and the radar and satellite
values at the gauge. Structural reads: 14, borrowed from the program's instruments (Saltimbanques islands / void /
dispersion / edge dependency; Core42 multiscale; HTAN persistence; the variance-of-disagreement determination; an infrared
third observer). Rule: go if Δρ ≥ 0.02 with a cluster-bootstrap 95% CI above 0 (clusters = UTC day × 5° tile).
Seismology guardrails applied: near-zero Ω check, and a spatial-scramble null.</p>
<p><b>Sample:</b> {F['n_samples']:,} gauge-hours (ASOS {F['by_src']['ASOS']:,}; USCRN {F['by_src']['USCRN']:,});
{F['clusters']['disc']} discovery and {F['clusters']['val']} validation clusters. Near-zero om_loc share
{F['near_zero_omloc_share']*100:.4f}%.</p>
<table><tr><th>Model</th><th>ρ validation</th><th>Δρ vs baseline</th><th>95% CI</th><th>Δρ under scramble null</th><th>Rule</th></tr>
{DEC_ROWS}
<tr><td>baseline</td><td>{dec['FULL_A2']['rho_base']:.4f}</td><td colspan="4"></td></tr></table>
<div class="kill"><b>Result:</b> stop rule met. No structural read adds held-out information about hourly radar error at
gauges. The small positive Δρ is not spatial: the scramble null, which breaks the radar–satellite pairing, gives a
<i>larger</i> increment ({nul['FULL_A2']['delta']:+.4f}), consistent with the lag and missingness terms.
<b>Sample:</b> {F['n_val']:,} validation gauge-hours. <b>Where it fails:</b> one summer regime; the label is a point.</div>
<img src="{fig('f2_feasibility.png')}" alt="feasibility">
<figcaption>A: increment over baseline, real vs scramble null, with the +0.02 go bar. B: each read's partial Spearman with
radar error after controlling for rain amount (discovery block): every structural read is within ±0.04, most against the
predicted positive sign; MRMS's own quality index also adds nothing. C: resolution-floor sweep (0.3°–1.0°): no rung is
correctly signed.</figcaption>
<h3>Anatomy of the label</h3>
<img src="{fig('f3_label.png')}" alt="label">
<figcaption>A: hourly radar vs gauge at the gauge pixel. B: absolute error grows with rain amount. C: structural
disagreement is largest in light, scattered rain, where absolute errors are small, so raw disagreement correlates
<i>negatively</i> with error (om_iou raw ρ {DD['ALL']['om_iou']['rho']:+.2f}).</figcaption>
<table><tr><th>read</th><th>raw ρ</th><th>partial ρ (all)</th><th>partial ρ (USCRN)</th><th>n (discovery)</th></tr>
{DESC_ROWS}</table>
<p class="note">Resolution sweep (partial ρ, n={RF['n']:,}): om_pat {RF['pat_w3']:+.3f} / {RF['pat_w5']:+.3f} /
{RF['pat_w7']:+.3f} / {RF['pat_w10']:+.3f}; om_loc {RF['loc_w3']:+.3f} / {RF['loc_w5']:+.3f} / {RF['loc_w7']:+.3f} /
{RF['loc_w10']:+.3f} at 0.3° / 0.5° / 0.7° / 1.0°. Gate 6 here: radar-error vs satellite-error at gauges r =
{G6['error_r_pearson_radar_vs_sat']:.2f} (PARTIAL; both are measured against the same gauge, which inflates r).
Predictor audit I = {PA['I']:.3f}, top pair {PA['top_pairs'][0][0]} {PA['top_pairs'][0][1]:.2f} (same correlation
before and after 3 × 3 smoothing).</p>

<h2>4. Literature step: why the null, and what maths exists</h2>
<p>Nineteen papers were pulled (<code>literature/INDEX.md</code>). Two bodies of work mattered:</p>
<ul>
<li><b>The null was predictable.</b> The area–point (gauge representativeness) term dominates radar–gauge differences at short
accumulations and stays significant even at 4-day totals (Ciach & Krajewski 1999; Ciach, Habib & Krajewski 2003). Satellite
products lose power and coherence with MRMS below ~200 km and ~4 h (Guilloteau et al. 2021), which is the scale of the
windows used.</li>
<li><b>The formal maths for the program's idea already exists</b> as triple / extended / quadruple collocation (Stoffelen 1998;
McColl et al. 2014; Alemohammad et al. 2015): each observer's error variance from three observers' covariances, with no truth.
Its core assumption, zero error cross-covariance, is what the rubric's Gate 6 checks. Precipitation applications use amounts:
multiplicative TC over the United States (Alemohammad et al. 2015), MRMS, gauge-only and IMERG for extreme events (Li et al. 2020),
and "binned" TC conditioned on rain regime (Pellet 2022), akin to the T6 stratification. Extensions handle more than three
observers and errors correlated within product groups (Pan et al. 2015; Gruber et al. 2016), and displacement and representation
error within the framework (Jacques 2025). No application to <i>structural</i> coordinates was found. Known: violating zero
cross-covariance makes collocation underestimate errors, the effect seen here.</li>
<li><b>A satellite–satellite dependence is physically expected.</b> Cold cloud tops and passive microwave over land both respond to
ice aloft; land microwave retrievals rely largely on the ice-scattering signal (Spencer, Goodman & Hood 1989). The open question is
whether the dependence can be seen and bounded from the observers alone.</li></ul>
<p>Mapping to program methods, and the four probes ranked from it, are in <code>RESEARCH_SYNTHESIS.md</code>. Probe A was run first
because it needs no gauge label.</p>

<h2>5. Experiment 2: structural collocation, (R, IMERG IR, MW) <span class="tag conf">confirmatory</span></h2>
<p><b>Design</b> (<code>PROBE_A_PLAN.md</code>). Non-overlapping 15 × 15-cell (1.5°) windows, every hour; a window-hour enters
when all four observers carry rain and microwave covers ≥ 80% of it. Per observer: wet fraction, intensity-weighted centroid
(y, x), dispersion, island count, and log mean rain as the value control. Estimator: ETC noise fraction 1 − ρ² for the triplet
(R, IR, MW); quadruple check of merged IMERG's error correlation with IR and MW using the other two as instruments; 500
cluster-bootstrap resamples. Predictions: (1) radar least noisy; (2) IR noisier than MW; (3) merged most error-correlated with MW.
Stop if most structural variables are unstable.</p>
<p><b>Sample:</b> {A['n']:,} window-hours, {A['clusters']} clusters.</p>
<table><tr><th>variable</th><th>R</th><th>IR (IMERG)</th><th>MW</th><th>boot. share negative</th><th>stable</th></tr>
{tc_rows(A, 'IR')}</table>
<table><tr><th>variable</th><th>merged–IR error corr.</th><th>merged–MW error corr.</th></tr>
{quad_rows(A)}</table>
<div class="find"><b>Result:</b> stop rule not triggered (4 of 5 structural variables stable; the value control unstable in
100% of resamples). Prediction 1 not met (0/5; microwave least noisy). Prediction 2 met (5/5). Prediction 3 met for
merged–MW (0.26–0.37); merged–IR ≈ 0, not the smaller positive value predicted.</div>
<h3>Post-hoc check <span class="tag expl">exploratory</span></h3>
<p>Hypothesis: radar and IR share an hour-integrated footprint the microwave snapshot lacks, biasing MW's error low. If so, MW
covering both half-hours should raise its noise fraction.</p>
<img src="{fig('f6_posthoc.png')}" alt="posthoc">
<figcaption>Observed the reverse (n both = {PH['n_both']:,}, one = {PH['n_one']:,}). The footprint explanation is not supported.</figcaption>
<p>The IMERG V07 ATBD (in <code>literature/</code>) then confirmed that IMERG's infrared estimate is PDIR, calibrated to the
passive-microwave estimates through a rotating 30-day accumulator, so (R, IR, MW) was not an independent triplet. That
motivated Experiment 3.</p>

<h2>6. Experiment 3: raw GOES infrared in place of IMERG IR <span class="tag conf">confirmatory</span></h2>
<p><b>Design</b> (<code>PROBE_A2_PLAN.md</code>): identical, with GOES-16 band 13 brightness temperature (4 scans per hour,
binned to 0.1°, cold-cloud field 235 − Tb) replacing IMERG IR. 502 of 504 hours (2 hours had &lt; 3 scans and were dropped per
plan). Test-hour alignment: radar vs cold cloud r = 0.39 at zero shift, best one cell displaced (geostationary parallax),
not corrected.</p>
<p><b>Sample:</b> {A2['n']:,} window-hours, {A2['clusters']} clusters.</p>
<table><tr><th>variable</th><th>R</th><th>GIR (GOES)</th><th>MW</th><th>boot. share negative</th><th>stable</th></tr>
{tc_rows(A2, 'GIR')}</table>
<table><tr><th>variable</th><th>merged–GOES IR error corr.</th><th>merged–MW error corr.</th></tr>
{quad_rows(A2)}</table>
<img src="{fig('f4_collocation.png')}" alt="collocation">
<figcaption>Noise fraction per observer and variable, Probe A (left) and A2 (right), with 95% bootstrap intervals. Swapping in
raw GOES infrared changed nothing about the ranking: microwave &lt; radar &lt; infrared on every structural variable.</figcaption>
<img src="{fig('f5_arithmetic.png')}" alt="arithmetic">
<figcaption>A: raw correlations (A2). Radar and microwave agree most; the cold-cloud field agrees more with microwave than with
radar, and that asymmetry drives the ranking: for centroid y, ρ²_MW = 0.85 × 0.69 / 0.61 ≈ 0.96 and ρ²_R = 0.85 × 0.61 / 0.69
≈ 0.75. B: structural variables stable, value variables not, in both probes. C: estimated error correlations of merged IMERG
with its inputs, recovered without truth.</figcaption>
<div class="find"><b>Result:</b> stop rule not triggered. Prediction 1 not met (0/5); prediction 2 met (5/5); prediction 3 partly
(merged–MW clearly positive on dispersion and islands; merged–GOES IR ≈ 0). The IMERG-IR calibration route was <b>not</b> the
cause of Probe A's ranking.</div>

<h2>7. Experiment 4: T1 synthetic anchor <span class="tag conf">confirmatory</span></h2>
<p><b>Why:</b> house rule I-20, never quote a read without a machine-perfect control. Probes A/A2 produced claims with no
control. <b>Design</b> (<code>T1_PLAN.md</code>): 20,000 rain-like 15 × 15 windows with known truth (1 + Poisson(1.5) Gaussian
cells, lognormal intensities). Simulated observers re-render the same cells: radar-like R (position sd 0.5, no blur),
microwave-like S1 (sd 0.8, blur +1), infrared-like S2 (sd 1.5, blur +2), merged M = 0.6·S1 + 0.4·S2, and an independent
fifth observer G (blur +1.5). A <b>planted shared error</b> d_sh, one displacement per window applied to both S1 and S2, swept
over σ_sh = 0, 0.5, 1.0, 1.5 cells. Reads and estimators are the Probe A code, unchanged.</p>
<table><tr><th>variable (σ_sh = 0)</th><th>R estimate / truth</th><th>S1 estimate / truth</th><th>S2 estimate / truth</th><th>boot. share negative</th></tr>
{T1ROWS}</table>
<table><tr><th>σ_sh</th><th>R (centroid y) est. / truth</th><th>S1 est. / truth</th><th>S2 est. / truth</th><th>S1–S2 error corr. (instrument G) est. / truth</th></tr>
{T1SWEEP}</table>
<img src="{FIG8}" alt="T1">
<figcaption>A: with no planted shared error (not "every assumption met": coordinate errors were already correlated), centroid and dispersion estimates sit near the identity line in size; islands and wet
fraction do not. B: as the planted shared satellite error grows, the estimate underrates S1 and overrates R; the ranking
flips from σ_sh = 0.5, reproducing the real-data pattern. C: the estimated S1–S2 error correlation tracks the planted one
but reads low with a blurred instrument; a sharp instrument (exploratory) recovers it, including the 0.16 present with
nothing planted.</figcaption>
<div class="find"><b>Scored:</b> P1 met for centroid (max error 0.034) and, in absolute terms only, dispersion (0.079; ranking wrong: radar-like {T1['0.0']['tc']['disp']['point']['nf_R']:.3f} vs true {T1['0.0']['truth']['disp']['nf_true']['R']:.3f}, microwave-like {T1['0.0']['tc']['disp']['point']['nf_MW']:.3f} vs {T1['0.0']['truth']['disp']['nf_true']['MW']:.3f}), not for islands or wet
fraction. P2: value collocation stable (0% negative) but ~0.06 off. P3 met (flip from σ_sh = 0.5). P4 partly (within 0.10 at
σ_sh ≥ 1.0; 0.15–0.21 low below). P5 not met (merged check degenerate: undefined or out of bounds on centroids).
<b>Unplanned result:</b> with nothing planted, the two blurred observers' structural errors correlate (centroid 0.16); a sharp instrument
recovers it (0.14 vs 0.16; 0.50 vs 0.50 at σ_sh = 1.0), a blurred one cannot. <b>Superseded in part (review):</b> the blur widened cells with
fixed peaks (mass not conserved; widths added linearly) and fields were clipped by windows, so this correlation may be a rendering or
windowing artefact rather than an effect of smoothing; all T1 results were provisional until T1-bis, which now supersedes them (§12).</div>

<h2>8. Experiment 5: T5: shared smoothness or shared physics? <span class="tag conf">confirmatory</span></h2>
<p><b>Design</b> (<code>T5_PLAN.md</code>): blur radar with a Gaussian (normalised convolution) before windowing, over a ladder
σ = 0–3 cells, plus σ* chosen from field properties alone: the blur at which radar's median lag-1 autocorrelation matches the
microwave field's. Collocation on the independent triplet (radar, GOES IR, microwave), centroid y, centroid x and dispersion
only. A synthetic calibration arm (T1 fields, radar-like observer blurred) shows what each cause predicts.
<b>Prediction:</b> if smoothness drives the real ranking, nf_radar − nf_microwave falls to ≤ 0 at σ*; if shared physics
does, it stays above 0.</p>
<p>Smoothness match: microwave lag-1 {T5['smoothness_median_lag1']['MW']:.3f}; radar {T5['smoothness_median_lag1']['0.0']:.3f}
unblurred, {T5['smoothness_median_lag1']['0.5']:.3f} at 0.5 cell → σ* = {T5['sigma_star']} cell.</p>
<table><tr><th>radar blur σ</th><th>n</th><th>nf_R − nf_MW, centroid y</th><th>centroid x</th><th>dispersion</th></tr>
{T5ROWS}</table>
<img src="{FIG9}" alt="T5">
<figcaption>A: on real data the radar-minus-microwave gap shrinks slowly with blur but never crosses zero, even at σ = 3 (far
smoother than microwave). B: synthetic calibration: smoothness alone never creates the inversion; a shared displacement does,
and blurring radar only shrinks it, the same shape as A. C (exploratory): GOES cold cloud sits increasingly north of radar rain
with latitude (geostationary parallax, ~0.9 cell at 45–55°N); microwave does not follow, so parallax is not the shared error.</figcaption>
<table><tr><th>latitude band</th><th>n</th><th>GOES IR − radar (centroid y, cells; − = north)</th><th>MW − radar</th></tr>
{OFFROWS}</table>
<div class="find"><b>Result:</b> at σ* the gap stays above zero on both centroids (centroid y
{T5['real'][str(T5['sigma_star'])]['cy']['diff_R_minus_MW']:+.3f} {ci(T5['real'][str(T5['sigma_star'])]['cy']['diff_ci95'])}). (Superseded by T9's fixed-sample close match, which gives the same answer.) The tested
smoothing procedure is <b>not sufficient</b> to remove the inversion (match approximate: radar {T5['smoothness_median_lag1']['0.5']:.3f}
at the chosen step vs microwave {T5['smoothness_median_lag1']['MW']:.3f}; sample changes along the ladder). This does not exclude resolution
or representativeness effects in general. The real ladder has the shape of the synthetic shared-displacement case. Exploratory: the
structural reads detect GOES parallax unprompted; it is GOES-specific, not shared with microwave.</div>

<h2>9. Experiment 6: T6: do the satellites share an error by keying on convective cores? <span class="tag conf">confirmatory</span></h2>
<p><b>Design</b> (<code>T6_PLAN.md</code>): H_core says both satellites are pulled toward convective cores while radar sees the full
surface rain area. Tests on the (radar, GOES IR, microwave) triplet: P1 re-read every observer on its own core (cells ≥ its
window's 80th percentile); P2 split windows by radar light-rain mass fraction (&lt; 5 mm/h); P3 check whether the satellites
sit nearer radar's core centroid than its full-area centroid.</p>
<table><tr><th></th><th>prediction under H_core</th><th>result</th></tr>
<tr><td>P1</td><td>core reading shrinks nf_R − nf_MW</td><td><b>not met</b>: core {T6['core']['cy']['diff_R_minus_MW']:+.3f} {ci(T6['core']['cy']['diff_ci95'])} vs full {T6['full']['cy']['diff_R_minus_MW']:+.3f} (centroid y)</td></tr>
<tr><td>P2</td><td>inversion grows with light-rain fraction</td><td><b>met in form</b>: top − bottom tercile {T6['strata']['cy']['top_minus_bottom']:+.3f} {ci(T6['strata']['cy']['top_minus_bottom_ci'])}; but two-thirds of windows are all light rain, so the split is "has a heavy core" vs "all light"</td></tr>
<tr><td>P3</td><td>satellites nearer radar's core</td><td><b>reversed</b>: nearer the full area (median paired difference {T6['geometry']['GOES_IR']['median_paired_diff']:+.2f} cells GOES, {T6['geometry']['MW']['median_paired_diff']:+.2f} MW); caveat: core centroids are noisier</td></tr></table>
<p><b>Exploratory follow-ups:</b> where does the inversion live? By radar window peak intensity, and by MRMS radar quality index
(beam overshoot is the obvious alternative: radar could be genuinely worse in light, shallow rain far from the radar).</p>
<table><tr><th>radar peak (mm/h)</th><th>n</th><th>nf_R − nf_MW, centroid y</th><th>nf_R</th><th>nf_MW</th></tr>
{T6INT}</table>
<table><tr><th>window RQI</th><th>weak rain (below median peak)</th><th>strong rain</th></tr>
{T6RQI}</table>
<img src="{FIG10}" alt="T6">
<figcaption>A: reading cores leaves the inversion in place. B (exploratory): the inversion triples from the strongest to the weakest
windows. C (exploratory): it does not fall with better radar beam quality; it is largest in the best-beam, weak-rain windows.</figcaption>
<div class="find"><b>Result:</b> mixed against H_core; core keying is not supported. The shared satellite error concentrates in
light rain, and radar beam quality (RQI) does not explain it.</div>

<h2>10. Experiment 7: T7: whose error is it? Gauges as a fourth observer <span class="tag conf">confirmatory</span></h2>
<p><b>Why:</b> with three observers the inversion is algebraically r(GOES, MW) &gt; r(radar, GOES); it cannot say whose error it is.
T1 showed a sharp, independent fourth observer can. Rain gauges are that observer: point measurements of surface rain used by
neither MRMS radar-only, IMERG Late nor GOES. <b>Design</b> (<code>T7_PLAN.md</code>): Probe A2 windows with ≥ 5 reporting gauges;
every observer read at the gauge points (its value in each gauge's 0.1° cell); point centroids; split at the median radar peak.
P0: the (radar, GOES, MW) triplet reproduces the inversion at point support. P1: replace GOES with gauges; a shared satellite
error predicts the inversion vanishes, a radar light-rain error predicts it stays. P2: GOES–MW error correlation with radar and
gauges as instruments; shared satellite error predicts &gt; 0.</p>
<p><b>Sample:</b> {T7['n_total']:,} window-hours (median {T7['gauges_per_window_median']:.0f} gauges per window); split at
{T7['median_peak']:.1f} mm/h. This sample is wetter than T6's (rain is required at the points).</p>
<table><tr><th>variable</th><th>stratum</th><th>P0: (R, GOES, MW) nf_R − nf_MW</th><th>P1: (R, MW, gauges) nf_R − nf_MW</th><th>P2: GOES–MW error corr.</th></tr>
{T7ROWS}</table>
<img src="{FIG11}" alt="T7">
<figcaption>A, B: noise fractions at point support. With the two satellites as the other legs, radar looks noisier than microwave;
with gauges in place of GOES, radar is the least noisy observer and microwave the noisiest. C: the GOES–microwave error
correlation, measured with radar and gauges as instruments.</figcaption>
<div class="find"><b>Result:</b> P0, P1 and P2 met, in both strata. Replacing a satellite with gauges reverses the ranking
(weak, centroid y: radar 0.063, gauges 0.225, microwave 0.259). With radar and gauges as instruments, the GOES and microwave errors
correlate at +0.35 to +0.51; this estimate, and reading the inversion in Probes A, A2, T5 and T6 as an artefact of a shared satellite
error, are conditional on radar and gauge errors being independent. Not established: the cause.</div>
<h3>Robustness check of T7 <span class="tag conf">confirmatory</span> (pre-registered, <code>T7R_PLAN.md</code>)</h3>
<p>T7 replaced GOES with gauges. Two checks it left open: the mirror triplet (replace microwave instead), and four-observer
consistency. With four observers and one common signal, independent errors make three products of correlations equal
("vanishing tetrads"). A single shared GOES–microwave error inflates only r(GOES, MW), so the tetrad that excludes that pair (t3)
should stay ≈ 0 while the two that include it (t1, t2) turn positive. Same sample, bootstrap seed 11.</p>
<table><tr><th>check</th><th>weak, y</th><th>weak, x</th><th>strong, y</th><th>strong, x</th><th>predicted</th></tr>
{T7RROWS}</table>
<table><tr><th>triplet (weak, centroid y)</th><th>radar</th><th>GOES IR</th><th>microwave</th><th>gauges</th></tr>
{T7RNF}</table>
<div class="find"><b>Result:</b> R1, R2 and R3 met in both strata and on both axes. Radar is least noisy whichever satellite is
dropped, and its noise fraction does not depend on the choice. Triplets holding both satellites deflate microwave and inflate the
third leg; triplets holding one satellite agree with each other. The tetrads are compatible with a single two-group split, {{radar, gauges}} vs
{{GOES, microwave}}; since t2 − t1 = t3 they are two independent constraints, and an interval containing zero does not rule out other
dependence. <b>Identification caveat</b> (found afterwards by <code>tetrad_audit</code>): four observers
cannot tell a shared GOES–microwave error (relative to surface rain) from a shared radar–gauge error (relative to a cloud-system
signal). Reading it as a satellite error rests on taking surface rain as the target; a fifth observer, under explicit assumptions about its own signal and errors, could
separate the two. It should sense surface rain (disdrometers, dense-network areal averages); lightning needs ice aloft and would
likely side with the satellites.</div>

<h2>11. Experiment 8: T8: is the shared error driven by cloud and ice aloft beyond the rain? <span class="tag conf">confirmatory</span></h2>
<p><b>Design</b> (<code>T8_PLAN.md</code>): anvil fraction = share of GOES cold-cloud cells with no radar surface rain. On the T7
gauge-instrumented sample, split into anvil terciles (means {T8['anvil_mean_by_tercile'][0]:.2f}, {T8['anvil_mean_by_tercile'][1]:.2f},
{T8['anvil_mean_by_tercile'][2]:.2f}). P1: GOES–MW error correlation grows with anvil. P2: the (R, GOES, MW) inversion grows. P3: the
(R, MW, gauges) ranking keeps radar least noisy throughout. Anvil fraction correlates {T8['corr_anvil_rpeak']:+.2f} with radar peak.</p>
<table><tr><th>anvil (mean)</th><th>GOES–MW error corr.</th><th>(R, GOES, MW) nf_R − nf_MW</th><th>(R, MW, G) nf R / MW / G</th><th>shared cov. (cells²)</th><th>error var. GOES / MW / gauges</th></tr>
{T8ROWS}</table>
<img src="{FIG12}" alt="T8">
<figcaption>A: no detectable change in the normalised coupling of the two satellites' errors with anvil extent. B: how badly the triplet is
fooled grows threefold. C (exploratory): in absolute terms the shared covariance doubles, but every observer's error variance grows
with it, gauges included.</figcaption>
<div class="find"><b>Result:</b> mixed. P1 not met (top − bottom {T8['cy']['c_IR_MW_top_minus_bottom']:+.3f}
{ci(T8['cy']['c_IR_MW_top_minus_bottom_ci'])}); P2 met ({T8['cy']['RGM_diff_R_MW_top_minus_bottom']:+.3f}
{ci(T8['cy']['RGM_diff_R_MW_top_minus_bottom_ci'])}); P3 holds. No change in the satellites' error coupling was detected across terciles; anvil-heavy (lighter-rain) windows are noisier for every observer, so anvil fraction is confounded and not identified as a
specific driver.</div>

<h2>12. Corrected controls and replication: T1-bis, T7-S, T9, T10 <span class="tag conf">confirmatory</span></h2>
<p><b>Why:</b> the external review (see the corrections box above) found that T1's blur did not conserve mass and that the satellite
attribution assumed independent radar and gauge errors. Two pre-registered checks followed, with a decision gate after them.</p>
<h3>T1-bis: synthetic assumption audit (<code>T1BIS_PLAN.md</code>)</h3>
<p>Truth rendered on a 64 × 64 canvas; blur is a true, mass-conserving Gaussian convolution; reads on the full canvas and in 15 × 15
windows (central and shifted). Analytical reference passed: centroid shift ≈ 10⁻¹⁵ cell, squared dispersion grows by 2σ² within 0.04%.</p>
<table><tr><th>region</th><th>S1–S2 error corr. (cy)</th><th>(cx)</th><th>(d)</th><th>max |error|, ranking: cy</th><th>d</th><th>d²</th></tr>
{TBROWS}</table>
<p>Planted shared S1–S2 displacement, centroid y, full canvas: at σ_sh = 0.5 estimates R {TB['S_E']['0.5']['full']['est']['R']:.3f} /
S1 {TB['S_E']['0.5']['full']['est']['S1']:.3f} (truth {TB['S_E']['0.5']['full']['true']['R']:.3f} / {TB['S_E']['0.5']['full']['true']['S1']:.3f}), not inverted;
at σ_sh = 1.0 R {TB['S_E']['1.0']['full']['est']['R']:.3f} / S1 {TB['S_E']['1.0']['full']['est']['S1']:.3f} (truth {TB['S_E']['1.0']['full']['true']['R']:.3f} /
{TB['S_E']['1.0']['full']['true']['S1']:.3f}), <b>inverted</b>. Five-observer benchmark (tetrad audit, 10 blocks): planted S1–S2 found 10/10,
planted R–Gs (mirror) found 10/10, no false alarms with nothing planted.</p>
<h3>T7-S: radar–gauge sensitivity (<code>T7S_PLAN.md</code>)</h3>
<p>One-signal model for radar, gauges, GOES and microwave with an assumed radar–gauge error correlation ρ_RG (swept 0–0.5) and a free
GOES–microwave error correlation, fitted by maximum likelihood. The fit is identical at every ρ_RG (the data cannot fix it), so the sweep
traces the equivalence class.</p>
<table><tr><th>case</th><th>GOES–MW error corr. at ρ_RG = 0</th><th>ρ*₁ (no longer detected)</th><th>ρ*₂ (ranking flips)</th><th>admissible</th></tr>
{TSROWS}</table>
<img src="{FIG14}" alt="T1-bis and T7-S">
<figcaption>A: with a mass-conserving blur, the satellite-like observers' centroid errors are uncorrelated on the full canvas (dispersion keeps +0.021); windows add a
little (dashed: the T1 value from the non-conserving blur). B: a planted shared satellite error inverts the ranking at σ_sh = 1.0. C: the
inferred GOES–microwave error correlation (solid) and radar minus microwave noise fraction (dotted) as the assumed radar–gauge error
correlation grows.</figcaption>
<div class="find"><b>Result:</b> T1-bis met its predictions: centroid errors uncorrelated on the full canvas (dispersion +0.021, within the pre-registered |r| &lt; 0.05 but not zero), all variables recovered with
correct rankings, the planted-error inversion reproduced at σ_sh = 1.0, and the five-observer design identifies a shared pair in either
direction. F6 is withdrawn. T7-S: the satellite dependence stays detected for assumed radar–gauge error correlations up to 0.475–0.5; radar's ranking as least
noisy is clear by interval up to 0.275–0.4 and its point estimate flips at 0.425–0.5. (Corrected after the third review: the earlier
"ranking survives up to about 0.4" was read from point estimates.) See <code>GATE.md</code>.</div>
<h3>T9: real-data window grid and fixed-sample smoothness (<code>T9_PLAN.md</code>)</h3>
<table><tr><th>check</th><th>n</th><th>nf_R − nf_MW, centroid y</th><th>centroid x</th></tr>
<tr><td>baseline grid (15 × 15)</td><td>{T9['W']['W0']['n']:,}</td><td>{T9['W']['W0']['cy']['diff']:+.3f} {ci(T9['W']['W0']['cy']['ci'])}</td><td>{T9['W']['W0']['cx']['diff']:+.3f}</td></tr>
<tr><td>half-shifted grid</td><td>{T9['W']['W1']['n']:,}</td><td>{T9['W']['W1']['cy']['diff']:+.3f} {ci(T9['W']['W1']['cy']['ci'])}</td><td>{T9['W']['W1']['cx']['diff']:+.3f}</td></tr>
<tr><td>25 × 25 windows</td><td>{T9['W']['W2']['n']:,}</td><td>{T9['W']['W2']['cy']['diff']:+.3f} {ci(T9['W']['W2']['cy']['ci'])}</td><td>{T9['W']['W2']['cx']['diff']:+.3f}</td></tr>
<tr><td>point reversal (radar, gauges, MW), half-shifted grid</td><td>{T9['W_points']['W1']['n']:,}</td><td>{T9['W_points']['W1']['cy']['diff']:+.3f} {ci(T9['W_points']['W1']['cy']['ci'])}</td><td>{T9['W_points']['W1']['cx']['diff']:+.3f}</td></tr>
<tr><td>point reversal, 25 × 25 windows</td><td>{T9['W_points']['W2']['n']:,}</td><td>{T9['W_points']['W2']['cy']['diff']:+.3f} {ci(T9['W_points']['W2']['cy']['ci'])}</td><td>{T9['W_points']['W2']['cx']['diff']:+.3f}</td></tr>
<tr><td>fixed sample, radar blurred to matched smoothness σ* = {T9['S_sigma_star']}</td><td>{T9['S'][str(T9['S_sigma_star'])]['n']:,}</td><td>{T9['S'][str(T9['S_sigma_star'])]['cy']['diff']:+.3f} {ci(T9['S'][str(T9['S_sigma_star'])]['cy']['ci'])}</td><td>{T9['S'][str(T9['S_sigma_star'])]['cx']['diff']:+.3f}</td></tr>
<tr><td>fixed sample, σ = 3</td><td>{T9['S']['3.0']['n']:,}</td><td>{T9['S']['3.0']['cy']['diff']:+.3f} {ci(T9['S']['3.0']['cy']['ci'])}</td><td>{T9['S']['3.0']['cx']['diff']:+.3f}</td></tr></table>
<div class="find"><b>Result:</b> predictions met. The window grid does not drive the inversion (shifted and larger windows within
20% of baseline; point reversal intact), and on a fixed sample with radar matched to microwave smoothness (lag-1
{T9['S_lag1'][str(T9['S_sigma_star'])]:.3f} vs {T9['S_lag1']['MW']:.3f}) the gap is unchanged. The inversion is robust to the tested grid shift and window size and to the smoothness match; windowing and
representativeness effects in general are not excluded (the 25 × 25 gap is about a fifth smaller).</div>

<h3>T10: replication on an unseen period (<code>T10_PLAN.md</code>)</h3>
<p>The July pipeline, frozen with script hashes, run on 2024-08-01 01:00 to 08-22 00:00 UTC; August data were fetched only after
the plan was frozen. Before any August analysis the period switch was added and July was rerun: every results file reproduced
byte-for-byte. An independent re-implementation from the plan files alone (<code>audit/AUDIT.md</code>) also reproduced 101 July
quantities to 4 × 10⁻⁹. Gates: {X10['n_points']:,} point-support window-hours, {X10['n_area']:,} area window-hours, largest
negative-variance share {X10['max_negative_share']:.2f}; not inconclusive.</p>
<table><tr><th></th><th>quantity</th><th>July (in-sample)</th><th>August (T10)</th><th>met</th></tr>
{"".join(f"<tr><td>{a['h']}</td><td>{a['quantity']}</td><td>{j['point']:+.3f} {ci(j['ci95'])}</td><td>{a['point']:+.3f} {ci(a['ci95'])}</td><td>{'yes' if a['met'] else 'no'}</td></tr>" for j, a in zip(X10J['rows'], X10['rows']))}</table>
<div class="find"><b>Result: {X10['verdict']}.</b> H1 area inversion, H2 point-support inversion, H3 reversal with gauges (both strata),
H4 mirror triplet, H5 GOES–microwave error correlation and H6 tetrad pattern all hold in August on both axes. Magnitudes are similar
or somewhat larger (c(GOES, MW) weak stratum {X10['rows'][12]['point']:+.3f} / {X10['rows'][13]['point']:+.3f} vs July
{X10J['rows'][12]['point']:+.3f} / {X10J['rows'][13]['point']:+.3f}). Sensitivity: in August the satellite dependence stays detected
through 0.5, the top of the tested range; the ranking's point estimate has not reversed by 0.5 (gap {TSA['weak_cy']['grid']['0.5']['diff_R_MW']:+.3f}),
but its interval includes reversal from 0.375–0.4. In July the point estimate flipped at 0.425–0.5. Replication across adjacent months with the
same instruments, processing and season shows repeatability, not attribution. Deviations (all mechanical, listed in <code>T10_RESULTS.md</code>): the
data root is chosen by working directory; a missing stations-table producer was added (byte-identical on July); the H1 interval uses
the Probe A2 bootstrap draws with Δ added; IMERG was fetched with one download per file (arrays identical on a July test file).
What T10 does not change: four observers still cannot separate the satellite reading from correlated radar–gauge errors, and both
months are adjacent warm-season CONUS.</div>

<h2>13. What the evidence shows</h2>
<p class="note">Finding IDs are stable identifiers (the trail refers to them); the list is ordered by weight, not by ID.</p>
<div class="find"><b>F1: the original question is answered, negative <span class="tag conf">confirmatory</span>.</b>
Hourly point error of radar QPE is carried by rain amount and point-vs-pixel noise; no structural read, and not MRMS's own
quality index, separates beyond magnitude. <b>Sample:</b> {F['n_samples']:,} gauge-hours. <b>Where it fails:</b> one
warm-season month; longer accumulations and areal labels untested.</div>
<div class="find"><b>F2 (revised after T1-bis): collocation on structural coordinates recovers centroid and dispersion with correct
rankings on a corrected synthetic control <span class="tag conf">confirmatory</span>.</b> With a mass-conserving blur, maximum error
≤ 0.006 (centroid) and ≤ 0.031 (dispersion) across full canvas and windows, rankings correct. The original T1 dispersion mis-ranking
came from the non-conserving blur. Islands and wet fraction are unvalidated (they failed in the superseded simulation and were not re-tested)
and excluded from interpretation. The real-data contrast (value collocation unstable: {A['logmean']['share_boot_negative_ev']*100:.0f}%
and {A2['logmean']['share_boot_negative_ev']*100:.0f}% impossible variances) is <b>not</b> a property of value collocation, which is
stable on synthetic fields; it is a property of the real observers' amount relationships.</div>
<div class="find"><b>F3 (revised after review): the four-observer pattern is consistent with an error shared by the two satellites
relative to surface rain <span class="tag conf">confirmatory</span> (the reversal) / <span class="tag hyp">conditional</span>
(the attribution).</b> With gauges as a fourth observer the ranking reverses to radar &lt; gauges &lt; microwave (T7); if radar and gauge
errors are independent, the GOES–microwave error correlation is +0.35 to +0.51, and the inference survives assumed radar–gauge error
correlations, satellite dependence detected up to 0.475–0.5 (July) and through the tested maximum of 0.5 (August); radar's ranking clear (interval excluding reversal) up to 0.275–0.4 (July), 0.375–0.4 (August) (T7-S; the actual value is not measured). The reversal, the error correlation and the tetrad pattern
replicated on August 2024 (T10). A planted shared error reproduces the original inversion on the corrected control (T1-bis, σ_sh = 1.0); the tested smoothing, parallax, convective-core keying and radar beam quality did not
explain it (T5, T6; parallax and beam quality exploratory). The shared error looked
concentrated in light rain in T6; T7's wetter sample shows it in both halves. T8: no change in its normalised strength (~0.4)
was detected across anvil extent, so "aloft-only cloud" in that tested form is not supported. Cause open; a coupling without detectable
change fits both
satellites observing the cloud system rather than surface rain everywhere.</div>
<div class="find"><b>F8: collocation with a sharp independent observer recovers the physically expected ranking on real data
<span class="tag conf">confirmatory</span>.</b> Radar 0.06, gauges 0.22 (point-vs-pixel representativeness), microwave 0.26,
from three observers and no truth, at point support; the same with GOES in place of microwave (radar 0.03), and the four-observer
tetrads are compatible with a single two-group split (T7-R; two independent constraints; four observers cannot say which side is "the error" without fixing the target, and the scores are true error rankings only if radar and gauge errors are independent). The method works on real data once its independence premise is met, and it detects when
that premise fails.</div>
<div class="find"><b>F6: WITHDRAWN (was: smoothing is a shared error mechanism).</b> T1-bis showed that with a mass-conserving blur the
centroid errors of blurred observers are uncorrelated on the full canvas (|r| ≤ 0.01; dispersion keeps a small +0.021); the T1 value (0.16) came from the blur
not conserving mass. Windows add a small shared error (centroid ≤ 0.02, dispersion up to 0.14). What survives: check independence after
feature extraction and windowing, and audit the simulation itself.</div>
<div class="find"><b>F10: five observers identify a shared pair <span class="tag conf">confirmatory</span> (in simulation).</b> On the
corrected control, the tetrad audit with five observers flagged a planted shared satellite pair and the mirror radar–gauge pair 10/10 each,
with no false alarms when nothing was planted. This is the design that would remove the real-data ambiguity, given a fifth observer whose
independence can be argued (an observer of surface rain, not lightning, which shares the dependence on ice aloft); the 10/10 holds
under the simulation's assumptions.</div>
<div class="find"><b>F9: how badly a three-observer collocation is fooled depends on conditions even when the coupling does not
<span class="tag conf">confirmatory</span>.</b> The radar–GOES–MW inversion rises from {T8['cy']['point']['RGM_diff_R_MW_0']:+.2f} to
{T8['cy']['point']['RGM_diff_R_MW_2']:+.2f} across anvil terciles while no change in the GOES–MW error correlation (near 0.4) is detected:
larger absolute errors appear to amplify the coupling's effect. A practical rule for collocation users: the bias from a shared error is worst where
all observers are noisy.</div>
<div class="find"><b>F7: the structural reads detect geostationary parallax unprompted
<span class="tag expl">exploratory</span>.</b> GOES cold-cloud centroids sit north of radar rain by an amount that grows with
latitude, {OFF['by_latitude_cy']['(45, 55]']['GIR_minus_R']:+.2f} cells at 45–55°N (~10 km), the known high-cloud parallax for a
satellite over 75°W. A validation that the reads see real geometry.</div>
<div class="find"><b>F4: who agrees with whom <span class="tag expl">exploratory</span>.</b> On centroid position, radar–microwave
r ≈ {A2['cy']['point']['r_R_MW']:.2f}; radar–GOES IR ≈ {A2['cy']['point']['r_R_IR']:.2f}; GOES IR–microwave ≈
{A2['cy']['point']['r_IR_MW']:.2f}.</div>
<div class="find"><b>F5: a rubric lesson <span class="tag expl">exploratory</span>.</b> Independence must be judged by the
physical quantity each observer responds to, not by the instrument, and derived products can be linked by calibration chains.
Added to RUBRIC.md Gate 6 (v2.1).</div>

<h2>14. What this does not show</h2>
<ul>
<li>It does not show the instrument is wrong: every read ran, the guards fired, and the nulls behaved.</li>
<li>The observer ranking is interpretable only with the gauge leg (T7), and only at point support and for centroid position;
the area-support rankings in Probes A, A2, T5 and T6 are not accuracy rankings.</li>
<li>It does not establish what the two satellites share; the tested "aloft-only anvil cloud" form is not supported (T8), and
the cloud-system explanation is untested.</li>
<li>It is two adjacent months (July, and the August replication), one season (warm convection), one domain (CONUS), hourly, 0.1°.</li>
<li>Several splits were exploratory (T5 offsets; T6 by intensity and RQI; T8 absolute errors). No multiple-comparison correction
was applied to them; they locate effects and set next tests, they do not confirm.</li>
<li>Selection: every collocation sample requires microwave coverage of ≥ 80% of the window, which selects hours near a microwave
overpass; T7 also requires rain at the gauge points, which skews it wetter.</li>
<li>The T7 conclusions rest on gauge-dense windows (≥ 5 gauges); they replicated in August (T10) under the same design.</li></ul>

<h2>15. Open threads (investigator next steps)</h2>
<div class="open"><b>T1: synthetic anchor: done (§7).</b></div>
<div class="open"><b>T5: shared smoothness vs shared physics: done (§8).</b> The tested smoothing is not sufficient; a fixed-sample,
close-match version was run in T9 (§12).</div>
<div class="open"><b>T1-bis, T7-S and T9: done (§12).</b> Decision gate: <code>GATE.md</code>. On real data the inversion is robust to the tested window grids and the smoothness match.
T10 (August replication): done, replicated (§12). Remaining, not started: a real fifth observer of surface rain (disdrometers or dense-network
areal averages; not lightning); lineage benchmark; a different season or region.</div>
<div class="open"><b>T6: core vs full area: done (§9).</b> Core keying not supported; error lives in light rain.</div>
<div class="open"><b>T7: whose error is it: done (§10).</b> A dependence pattern consistent with a shared satellite error, conditional on radar–gauge independence (§10, §12); replicated in T10. T3 (surface observer) is answered by T7.</div>
<div class="open"><b>T8: anvil fraction: done (§11).</b> No detectable change in coupling; anvil confounded with light rain.</div>
<div class="open"><b>Cause of the shared error (candidates, not started).</b> (a) Joint anvil × intensity stratification to separate the
two; (b) cloud-system vs rain-field geometry: compare the satellites' centroids with the GOES cold-cloud centroid itself (if MW
tracks the cloud system, its error should point toward the cold-cloud centroid); (c) a different microwave product (GPROF swath,
not IMERG-gridded) to rule out IMERG gridding/processing as the shared step.</div>
<div class="open"><b>T2: blind lineage test.</b> Products with a documented family tree: IMERG Early / Late / Final (shared
inputs, different morphing and gauge calibration), MRMS radar-only vs multi-sensor, Stage IV. Does structural collocation
recover the known tree without being told it?</div>
<div class="open"><b>T3: the ice-aloft hypothesis.</b> Replace one satellite with a surface-rain observer that has no ice
sensitivity and is sharp (T1: a blurred instrument cannot see shared error), e.g. a dense-mesonet areal gauge field.
Prediction: radar's noise fraction falls below microwave's.</div>
<div class="open"><b>T4: longer accumulations.</b> Re-run Experiment 1 at daily or event totals and areal labels, to map how the
point-noise floor moves with accumulation.</div>

<h2>16. Chronological trail</h2>
<ol class="trail">
<li>Rubric v2 written; four candidates desk-scored; precipitation (observation disagreement) chosen for the only fully public,
powered data.</li>
<li>Step 0: confirmed open access to MRMS, IMERG (Earthdata token), USCRN, ASOS. Re-stated the target to radar-only error so
every gauge is independent of the predictors.</li>
<li>Feasibility plan frozen; Amendment 1 (Saltimbanques, Core42, HTAN reads; infrared third observer) and Amendment 2
(variance of disagreement, resolution floor, seismology guardrails) added before any result.</li>
<li>Data pulled and reduced in memory (MRMS 1,008 files, IMERG 2 × 1,008, gauges 1.1 M rows); grid and timing alignment checked
(correlation peaks at zero offset).</li>
<li>Experiment 1 read once: stop rule met; scramble null larger than real.</li>
<li>Literature pulled (19 papers): representativeness explains the null; triple collocation is the formal maths.</li>
<li>Probe A frozen and run: stable structure, ranking MW &lt; R &lt; IR; post-hoc footprint idea refuted; ATBD confirms IMERG IR
is microwave-calibrated.</li>
<li>Probe A2 frozen and run with raw GOES infrared: same ranking; calibration route excluded; ice-aloft hypothesis proposed;
rubric Gate 6 updated.</li>
<li>Original question closed; collocation thread kept open; lab book started.</li>
<li>T1 synthetic anchor frozen and run: structural collocation valid for centroid and dispersion; shared error reproduces the real ranking; smoothing found to be a shared error mechanism; F2 and F3 revised; T5 set as next.</li>
<li>T5 frozen and run: radar blurred to (and beyond) microwave smoothness, inversion persists; synthetic arm matches the shared-displacement case; exploratory check finds GOES parallax, not shared with microwave; F3 revised, F7 added.</li>
<li>T6 frozen and run: core keying not supported (P1 not met, P3 reversed, P2 met in form); exploratory splits put the shared error in light rain, independent of radar beam quality; F3 revised; T7 set.</li>
<li>T7 frozen and run: gauges as a sharp fourth observer; shared GOES–microwave error measured (+0.35 to +0.51); replacing a satellite with gauges reverses the ranking (radar least noisy); F3 resolved in existence and size; F8 added; T8 set.</li>
<li>T8 frozen and run: no detectable change in GOES–MW coupling across anvil terciles (P1 not met); the triplet inversion grows threefold (P2 met); exploratory absolute errors show all observers noisier in anvil-heavy windows; F3 updated, F9 added; cause candidates set (first labelled T9; the name was later given to the window and smoothness test).</li>
<li>Lab book consolidated (summary, glossary, data manifest, limitations). T7-R robustness check frozen and run: mirror triplet and tetrads all as predicted.</li>
<li>External review (GPT) checked against the code: blur not mass-conserving, dispersion mis-ranked, tetrads dependent, several
"not detected" read as "excluded", attribution overstated. Text corrected; F2, F3 revised; F6 under review; T1-bis and T7-S frozen.</li>
<li>T1-bis run: corrected blur passes the analytical reference; coordinate errors independent on the full canvas; all variables recovered
with correct rankings; planted inversion at σ_sh = 1.0; five-observer benchmark 10/10 with no false alarms. F6 withdrawn, F2 reinstated,
F10 added. T7-S run: attribution survives assumed radar–gauge error correlation up to ≈ 0.4. Decision gate (GATE.md).</li>
<li>T9 frozen and run: shifted and larger windows leave the inversion and point reversal intact; fixed-sample close smoothness match (lag-1 within 0.013) leaves
the gap unchanged. Paused at the gate for discussion.</li>
<li>Second external review: wording tightened (close match, not exact; robust to tested variants, not closed; threshold counts
unvalidated; T7-S bound and five-observer result stated as conditional on their assumptions). No numbers changed. Replication plan
T10 written and frozen before any August data was fetched.</li>
<li>Third external review (after T10): "beyond 0.5" overstated a grid that stops at 0.5; sensitivity now reported by interval in both
months (dependence to 0.475–0.5; ranking clear to 0.275–0.4 in July, 0.375–0.4 in August; the earlier "ranking survives to about 0.4" was
point-estimate based). "Excluded" softened to "tested, not supported within the tested variants". Replication stated as repeatability,
not attribution.</li>
<li>Technical note review: "coordinate errors independent on the full canvas" narrowed to centroids; dispersion keeps a small
residual (+0.021 in T1-bis; +0.031 [+0.005, +0.056] in the note's minimal example), origin unresolved. Note at
<code>notes/synthetic-shared-error/</code>.</li>
<li>Claim ledger added (REVIEW.md, this section, paper §5.1): the tetrad result stated first as a failure of the one-signal,
independent-error model; untested comparisons with existing methods listed as not claimed.</li>
<li>Independent re-implementation audit (separate session, from plan files only): 101 July quantities reproduced to 4 × 10⁻⁹; four
details fixed by code rather than plans recorded in <code>audit/AUDIT.md</code>.</li>
<li>T10 run: period switch added, July reproduced byte-for-byte; August fetched (504 MRMS and GOES hours, 1,008 IMERG half-hours,
2,295 gauge stations); all six pre-registered predictions met. Verdict: Replicated. Paused for discussion.</li></ol>

<h2>17. Reproducibility</h2>
<p>Everything is in <code>Documents/russell/precip-disagreement/</code>. Raw files were never stored: each was streamed, reduced
(block means, crops, gauge-pixel samples) and discarded. Reduced arrays live in <code>data/</code> (git-ignored, ~1.9 GB).
Credentials: an Earthdata token in <code>secrets/</code> (git-ignored, owner-only).</p>
<pre>cd precip-disagreement
uv venv --python 3.12 .venv &amp;&amp; uv pip install --python .venv/bin/python numpy scipy pandas h5py eccodes eccodeslib matplotlib pypdf
.venv/bin/python scripts/fetch_gauges.py          # USCRN + ASOS -> data/gauges.csv
.venv/bin/python scripts/make_stations.py        # unique stations -> data/stations.csv (reproduces the original inline step)
.venv/bin/python scripts/fetch_mrms.py            # MRMS radar-only + RQI, top-of-hour -> data/mrms/
.venv/bin/python scripts/fetch_imerg.py           # IMERG Late merged, cropped -> data/imerg/
.venv/bin/python scripts/fetch_imerg_extra.py     # IMERG IR, MW, randomError -> data/imerg_x/
.venv/bin/python scripts/analyze.py               # Experiment 1 -> data/results.json, samples*.csv
PYTHONPATH=scripts .venv/bin/python scripts/probe_a.py          # Probe A
PYTHONPATH=scripts .venv/bin/python scripts/probe_a_posthoc.py  # post-hoc (labelled)
.venv/bin/python scripts/fetch_goes.py            # GOES-16 band 13 -> data/goes/
PYTHONPATH=scripts .venv/bin/python scripts/probe_a2.py         # Probe A2
PYTHONPATH=scripts .venv/bin/python scripts/t1_synthetic.py     # T1 synthetic anchor (seed 20261015)
PYTHONPATH=scripts .venv/bin/python scripts/t1_explore.py       # T1 exploratory E1, E2
PYTHONPATH=scripts .venv/bin/python scripts/t5.py               # T5 (synthetic arm, smoothness match, ladder)
.venv/bin/python scripts/t5_offsets.py                          # T5 exploratory offsets
PYTHONPATH=scripts .venv/bin/python scripts/t6.py               # T6 core vs area
PYTHONPATH=scripts .venv/bin/python scripts/t6_intensity.py     # T6 exploratory: by intensity
PYTHONPATH=scripts .venv/bin/python scripts/t6_rqi.py           # T6 exploratory: by RQI
PYTHONPATH=scripts .venv/bin/python scripts/t7.py               # T7 gauges as fourth observer
.venv/bin/python scripts/t7_robust.py                           # T7-R robustness (mirror triplet, tetrads)
.venv/bin/python scripts/t1bis.py                                # T1-bis corrected synthetic audit (seed 20261016)
.venv/bin/python scripts/t7s.py                                  # T7-S radar–gauge sensitivity (seed 13)
PYTHONPATH=scripts .venv/bin/python scripts/t9.py               # T9 window grid + fixed-sample smoothness (seed 21)
PYTHONPATH=scripts .venv/bin/python scripts/t10_area.py _july    # T10 H1 interval on July (same draws as Probe A2)
.venv/bin/python scripts/t10_score.py _july                      # T10 scorer, July reference
# T10 August: cd runs/aug (data -> ../../data_aug), export PRECIP_PERIOD=aug, then
#   fetch_gauges, make_stations, fetch_mrms, fetch_goes, fetch_imerg_both; probe_a2, t10_area, t7, t7_robust, t7s, t10_score
PYTHONPATH=scripts .venv/bin/python scripts/t8.py               # T8 anvil fraction
PYTHONPATH=scripts .venv/bin/python scripts/t8_abs.py           # T8 exploratory absolute errors
.venv/bin/python scripts/make_figures.py &amp;&amp; .venv/bin/python build_labbook.py</pre>
<p><b>Environment at build:</b> Python {platform.python_version()} · numpy {np.__version__} · scipy {scipy.__version__} ·
pandas {pandas.__version__} · h5py {h5py.__version__} · matplotlib {matplotlib.__version__} · {platform.platform()}.
<b>Randomness:</b> bootstrap seeds fixed in code (scramble null 20261015; decision bootstrap 7; collocation 11; T1 fields 20261015, G bootstrap 3, E2 noise 99; T5 bootstrap 5; T6 bootstraps 5, 6, 8; T7 bootstrap 7; T8 bootstrap 9; T7-R bootstrap 11; T7-S 13; T9 21; T10 reuses each analysis's own seed).
<b>Known deviations from plan:</b> Experiment 1 excludes gauges within 1° of the domain edge (so the 1.0° sweep window fits);
Probe A2 dropped 2 GOES hours per its own rule; the stations table was built inline after the gauge fetch (now reproduced by <code>scripts/make_stations.py</code>, byte-identical).</p>
<table><tr><th>file</th><th>SHA-256 (first 16)</th></tr>
{HASH_ROWS}</table>
<p><b>Reduced input data</b> (what every analysis reads; digest = SHA-256 over sorted file names and contents, first 16). Upstream
archives can be reprocessed, so a rerun of the fetch scripts may not reproduce these digests exactly.</p>
<table><tr><th>directory / file</th><th>files</th><th>size</th><th>digest</th></tr>
{DATA_ROWS}</table>
<p class="note">Data sources: MRMS <code>s3://noaa-mrms-pds</code>; GOES-16 <code>s3://noaa-goes16</code> (NOAA Open Data:
open, attribution requested); IMERG <code>GPM_3IMERGHHL.07</code> at NASA GES DISC (Earthdata login); USCRN
<code>ncei.noaa.gov/pub/data/uscrn/products/hourly02</code>; ASOS via Iowa Environmental Mesonet <code>asos.py</code>.
Literature in <code>literature/</code> (git-ignored PDFs; index in <code>literature/INDEX.md</code>).</p>
</body></html>"""

(H / "LABBOOK.html").write_text(HTML)
print("wrote LABBOOK.html", round(len(HTML) / 1e6, 2), "MB")
