"""Build PAPER.html: arXiv-style preprint companion to LABBOOK.html. Self-contained (inline CSS, base64 figures).
Numbers are read from data/*.json so the paper and the lab book cannot drift.
Run from precip-disagreement/:  .venv/bin/python scripts/make_figures.py && .venv/bin/python build_paper.py
"""
import base64, json, pathlib, datetime

H = pathlib.Path(".")
J = lambda p: json.load(open(H / p))
F = J("data/results.json"); A = J("data/probe_a_results.json"); A2 = J("data/probe_a2_results.json")
T1 = J("data/t1_results.json"); TX = J("data/t1_explore.json"); T5 = J("data/t5_results.json"); OFF = J("data/t5_offsets.json")
T6 = J("data/t6_results.json"); TI = J("data/t6_intensity.json"); TQ = J("data/t6_rqi.json")
T7 = J("data/t7_results.json"); T7R = J("data/t7r_results.json"); T8 = J("data/t8_results.json"); T8A = J("data/t8_abs.json")
TB = J("data/t1bis_results.json"); TS = J("data/t7s_results.json"); T9 = J("data/t9_results.json")
X10 = J("data_aug/t10_score.json"); X10J = J("data/t10_score_july.json"); TSA = J("data_aug/t7s_results.json")


def fig(n):
    return "data:image/png;base64," + base64.b64encode((H / "figures" / n).read_bytes()).decode()


def ci(c, d=3):
    return f"[{c[0]:+.{d}f}, {c[1]:+.{d}f}]"


s0 = T1["0.0"]
star = str(T5["sigma_star"])
def _r(S, i): r = S["rows"][i]; return f"{r['point']:+.3f} [{r['ci95'][0]:+.3f}, {r['ci95'][1]:+.3f}]"
TAB2 = "".join(f"<tr><td>{a['h']}</td><td>{a['quantity']}</td><td>{_r(X10J, k)}</td><td>{_r(X10, k)}</td><td>{'yes' if a['met'] else 'no'}</td></tr>"
               for k, a in enumerate(X10["rows"]))
V = {
    "TAB2": TAB2, "X10_NP": f"{X10['n_points']:,}", "X10_NA": f"{X10['n_area']:,}", "X10_V": X10["verdict"],
    "X10_A": f"{X10['rows'][0]['point']:+.3f}", "X10_P": f"{X10['rows'][4]['point']:+.3f}", "X10_C": f"{X10['rows'][12]['point']:+.3f}",
    "X10_TS": f"{TSA['weak_cy']['grid']['0.5']['diff_R_MW']:+.3f}",
    # experiment 1
    "E1_N": f"{F['n_samples']:,}", "E1_D": f"{F['decision']['FULL_A2']['delta']:+.4f}", "E1_CI": ci(F['decision']['FULL_A2']['ci95'], 4),
    "E1_NULL": f"{F['null_spatial_scramble']['FULL_A2']['delta']:+.4f}", "E1_BASE": f"{F['decision']['FULL_A2']['rho_base']:.3f}",
    # T1
    "T1_CY_ERR": f"{max(abs(s0['tc'][v]['point'][f'nf_{k}'] - s0['truth'][v]['nf_true'][k]) for v in ['cy','cx'] for k in ['R','IR','MW']):.3f}",
    "T1_DISP_ERR": f"{max(abs(s0['tc']['disp']['point'][f'nf_{k}'] - s0['truth']['disp']['nf_true'][k]) for k in ['R','IR','MW']):.3f}",
    "T1_ISL_ERR": f"{max(abs(s0['tc']['isl']['point'][f'nf_{k}'] - s0['truth']['isl']['nf_true'][k]) for k in ['R','IR','MW']):.2f}",
    "T1_WET_NEG": f"{s0['tc']['wet']['share_boot_negative_ev']*100:.0f}",
    "T1_LOG_ERR": f"{max(abs(s0['tc']['logmean']['point'][f'nf_{k}'] - s0['truth']['logmean']['nf_true'][k]) for k in ['R','IR','MW']):.3f}",
    "T1_05_R": f"{T1['0.5']['tc']['cy']['point']['nf_R']:.3f}", "T1_05_S": f"{T1['0.5']['tc']['cy']['point']['nf_MW']:.3f}",
    "T1_05_TR": f"{T1['0.5']['truth']['cy']['nf_true']['R']:.3f}", "T1_05_TS": f"{T1['0.5']['truth']['cy']['nf_true']['MW']:.3f}",
    "T1_BASE_EC": f"{s0['truth']['cy']['errcorr_true']['MW_IR']:.2f}",
    "T1_G0": f"{T1['0.0']['quadG']['cy']['point']:+.2f}", "T1_G10": f"{T1['1.0']['quadG']['cy']['point']:.2f}",
    "T1_T10": f"{T1['1.0']['truth']['cy']['errcorr_true']['MW_IR']:.2f}",
    "E1S0": f"{TX['E1_sig0.0']['cy']['quadG']['point']:.2f}", "E1S0T": f"{TX['E1_sig0.0']['cy']['true']:.2f}",
    "E1S1": f"{TX['E1_sig1.0']['cy']['quadG']['point']:.2f}", "E1S1T": f"{TX['E1_sig1.0']['cy']['true']:.2f}",
    # A / A2
    "A_N": f"{A['n']:,}", "A2_N": f"{A2['n']:,}",
    "A_R": f"{A['cy']['point']['nf_R']:.3f}", "A_I": f"{A['cy']['point']['nf_IR']:.3f}", "A_M": f"{A['cy']['point']['nf_MW']:.3f}",
    "A2_R": f"{A2['cy']['point']['nf_R']:.3f}", "A2_I": f"{A2['cy']['point']['nf_IR']:.3f}", "A2_M": f"{A2['cy']['point']['nf_MW']:.3f}",
    "A_LOGNEG": f"{A['logmean']['share_boot_negative_ev']*100:.0f}", "A2_LOGNEG": f"{A2['logmean']['share_boot_negative_ev']*100:.0f}",
    "A2_RRI": f"{A2['cy']['point']['r_R_IR']:.2f}", "A2_RRM": f"{A2['cy']['point']['r_R_MW']:.2f}", "A2_RIM": f"{A2['cy']['point']['r_IR_MW']:.2f}",
    "A2_CLU": f"{A2['clusters']:,}", "A_CLU": f"{A['clusters']:,}",
    "A_MMW_LO": f"{min(A[v]['point']['M_MW'] for v in ['cy','cx','disp','isl'] if A[v]['point']['M_MW']==A[v]['point']['M_MW']):.2f}",
    "A_MMW_HI": f"{max(A[v]['point']['M_MW'] for v in ['cy','cx','disp','isl'] if A[v]['point']['M_MW']==A[v]['point']['M_MW']):.2f}",
    "T1D_RE": f"{s0['tc']['disp']['point']['nf_R']:.3f}", "T1D_RT": f"{s0['truth']['disp']['nf_true']['R']:.3f}",
    "T1D_SE": f"{s0['tc']['disp']['point']['nf_MW']:.3f}", "T1D_ST": f"{s0['truth']['disp']['nf_true']['MW']:.3f}",
    "T5_LSTAR": f"{T5['smoothness_median_lag1'][str(T5['sigma_star'])]:.2f}",
    "T5_N0": f"{T5['real']['0.0']['n']:,}", "T5_N3": f"{T5['real']['3.0']['n']:,}",
    "TB_CY": f"{max(TB['S_D'][r]['cy']['max_abs_err'] for r in TB['S_D']):.3f}", "TB_D": f"{max(TB['S_D'][r]['d']['max_abs_err'] for r in TB['S_D']):.3f}",
    "TB_BC_CY": f"{TB['S_BC']['full']['cy']['S1-S2']:+.3f}", "TB_WIN_CY": f"{max(TB['S_BC'][r]['cy']['S1-S2'] for r in ['win','win_s1','win_s2']):.3f}",
    "TB_WIN_D": f"{max(TB['S_BC'][r]['d']['S1-S2'] for r in ['win','win_s1','win_s2']):.2f}",
    "TB_05R": f"{TB['S_E']['0.5']['full']['est']['R']:.3f}", "TB_05S": f"{TB['S_E']['0.5']['full']['est']['S1']:.3f}",
    "TB_05TR": f"{TB['S_E']['0.5']['full']['true']['R']:.3f}", "TB_05TS": f"{TB['S_E']['0.5']['full']['true']['S1']:.3f}",
    "TB_E1R": f"{TB['S_E']['1.0']['full']['est']['R']:.3f}", "TB_E1S": f"{TB['S_E']['1.0']['full']['est']['S1']:.3f}",
    "TB_E1TR": f"{TB['S_E']['1.0']['full']['true']['R']:.3f}", "TB_E1TS": f"{TB['S_E']['1.0']['full']['true']['S1']:.3f}",
    "TS_R2MIN": f"{min(TS[c]['rho_star_2'] for c in TS):.3f}", "TS_R2MAX": f"{max(TS[c]['rho_star_2'] for c in TS):.2f}",
    "A2_GAP": f"{T9['W']['W0']['cy']['diff']:+.3f}",
    "T9_W1": f"{T9['W']['W1']['cy']['diff']:+.3f}", "T9_W2": f"{T9['W']['W2']['cy']['diff']:+.3f}",
    "T9_P1": f"{T9['W_points']['W1']['cy']['diff']:+.3f}", "T9_P2": f"{T9['W_points']['W2']['cy']['diff']:+.3f}",
    "T9_N": f"{T9['S_fixed_n']:,}", "T9_STAR": f"{T9['S_sigma_star']}", "T9_LR": f"{T9['S_lag1'][str(T9['S_sigma_star'])]:.3f}",
    "T9_LM": f"{T9['S_lag1']['MW']:.3f}", "T9_G": f"{T9['S'][str(T9['S_sigma_star'])]['cy']['diff']:+.3f}",
    "T9_GCI": ci(T9['S'][str(T9['S_sigma_star'])]['cy']['ci']), "T9_G3": f"{T9['S']['3.0']['cy']['diff']:+.3f}",
    # T5
    "T5_STAR": star, "T5_GAP": f"{T5['real'][star]['cy']['diff_R_minus_MW']:+.3f}", "T5_GAPCI": ci(T5['real'][star]['cy']['diff_ci95']),
    "T5_GAP3": f"{T5['real']['3.0']['cy']['diff_R_minus_MW']:+.3f}", "T5_LMW": f"{T5['smoothness_median_lag1']['MW']:.2f}",
    "T5_LR0": f"{T5['smoothness_median_lag1']['0.0']:.2f}",
    "PAR": f"{OFF['by_latitude_cy']['(45, 55]']['GIR_minus_R']:+.2f}", "PARMW": f"{OFF['by_latitude_cy']['(45, 55]']['MW_minus_R']:+.2f}",
    # T6
    "T6_CORE": f"{T6['core']['cy']['diff_R_minus_MW']:+.3f}", "T6_FULL": f"{T6['full']['cy']['diff_R_minus_MW']:+.3f}",
    "T6_Q0": f"{TI['0']['cy']['diff']:+.3f}", "T6_Q4": f"{TI['4']['cy']['diff']:+.3f}", "T6_RQIW": f"{TQ['rqi2_peak0']['cy']['diff']:+.3f}",
    "T6_RQIM": f"{TQ['rqi2_peak0']['rqi_mean']:.2f}",
    # T7
    "T7_N": f"{T7['n_total']:,}", "T7_G": f"{T7['gauges_per_window_median']:.0f}",
    "T7_P0": f"{T7['weak']['cy']['point']['RGM_diff_R_MW']:+.3f}", "T7_P0CI": ci(T7['weak']['cy']['ci95']['RGM_diff_R_MW']),
    "T7_P1": f"{T7['weak']['cy']['point']['RMG_diff_R_MW']:+.3f}", "T7_P1CI": ci(T7['weak']['cy']['ci95']['RMG_diff_R_MW']),
    "T7_NR": f"{T7['weak']['cy']['point']['RMG_nf_R']:.3f}", "T7_NM": f"{T7['weak']['cy']['point']['RMG_nf_MW']:.3f}",
    "T7_NG": f"{T7['weak']['cy']['point']['RMG_nf_G']:.3f}",
    "T7_CW": f"{T7['weak']['cy']['point']['c_IR_MW']:+.2f}", "T7_CWCI": ci(T7['weak']['cy']['ci95']['c_IR_MW'], 2),
    "T7_CS": f"{T7['strong']['cx']['point']['c_IR_MW']:+.2f}",
    "T7R_T1": f"{T7R['weak']['cy']['point']['t1']:+.3f}", "T7R_T1CI": ci(T7R['weak']['cy']['ci95']['t1']),
    "T7R_T2": f"{T7R['weak']['cy']['point']['t2']:+.3f}", "T7R_T2CI": ci(T7R['weak']['cy']['ci95']['t2']),
    "T7R_T3": f"{T7R['weak']['cy']['point']['t3']:+.3f}", "T7R_T3CI": ci(T7R['weak']['cy']['ci95']['t3']),
    "T7R_RG": f"{T7R['weak']['cy']['point']['nf_R|RIRG']:.3f}", "T7R_IG": f"{T7R['weak']['cy']['point']['nf_IR|RIRG']:.3f}",
    "T7R_R2": f"{T7R['weak']['cy']['point']['R2_nfR_RIRG_minus_RMWG']:+.3f}", "T7R_R2CI": ci(T7R['weak']['cy']['ci95']['R2_nfR_RIRG_minus_RMWG']),
    # T8
    "T8_C0": f"{T8['cy']['point']['c_IR_MW_0']:.2f}", "T8_C2": f"{T8['cy']['point']['c_IR_MW_2']:.2f}",
    "T8_TB": f"{T8['cy']['c_IR_MW_top_minus_bottom']:+.3f}", "T8_TBCI": ci(T8['cy']['c_IR_MW_top_minus_bottom_ci']),
    "T8_I0": f"{T8['cy']['point']['RGM_diff_R_MW_0']:+.3f}", "T8_I2": f"{T8['cy']['point']['RGM_diff_R_MW_2']:+.3f}",
    "T8_COV0": f"{T8A['cy']['0']['shared_cov_GIR_MW']:.2f}", "T8_COV2": f"{T8A['cy']['2']['shared_cov_GIR_MW']:.2f}",
    "T8_G0": f"{T8A['cy']['0']['errvar_G']:.2f}", "T8_G2": f"{T8A['cy']['2']['errvar_G']:.2f}", "T8_CORR": f"{T8['corr_anvil_rpeak']:+.2f}",
    "DATE": datetime.date.today().isoformat(),
}
for k, f_ in [("FIG_THREAD", "f13_thread.png"), ("FIG_OBS", "f1_observers.png"), ("FIG_A", "f4_collocation.png"),
              ("FIG_T5", "f9_t5.png"), ("FIG_T6", "f10_t6.png"), ("FIG_T7", "f11_t7.png"), ("FIG_T8", "f12_t8.png"), ("FIG_E1", "f2_feasibility.png"), ("FIG_TB", "f14_t1bis_t7s.png")]:
    V[k] = fig(f_)

CSS = """
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:Georgia,'Times New Roman',serif;font-size:11pt;line-height:1.62;color:#1a1a1a;background:#f4f4f0}
.page{max-width:840px;margin:0 auto;background:#fff;padding:64px 76px 80px;box-shadow:0 0 30px rgba(0,0,0,.12);min-height:100vh}
.rh{font-size:8.5pt;color:#666;border-bottom:1px solid #ccc;padding-bottom:6px;margin-bottom:36px;font-family:'Helvetica Neue',Arial,sans-serif}
h1{font-size:20pt;line-height:1.25;text-align:center;margin-bottom:12px}
.authors{text-align:center;font-size:10.5pt;color:#333;margin-bottom:4px}
.aff{text-align:center;font-size:9pt;color:#666;margin-bottom:22px;font-family:'Helvetica Neue',Arial,sans-serif}
.abstract{border-top:1px solid #ccc;border-bottom:1px solid #ccc;padding:14px 18px;margin:18px 0 26px;font-size:10pt;background:#fafaf7}
.abstract b{font-variant:small-caps;letter-spacing:.04em}
.kw{font-size:9pt;color:#555;margin-top:8px;font-family:'Helvetica Neue',Arial,sans-serif}
h2{font-size:13pt;margin:30px 0 8px} h3{font-size:11pt;margin:20px 0 6px;font-style:italic}
p{margin:8px 0;text-align:justify}
ul,ol{margin:8px 0 8px 22px} li{margin:3px 0}
figure{margin:18px 0} figure img{width:100%;border:1px solid #e3e3e3}
figcaption{font-size:9pt;color:#333;margin-top:6px;font-family:'Helvetica Neue',Arial,sans-serif;line-height:1.45}
figcaption b{font-weight:700}
table{border-collapse:collapse;width:100%;margin:12px 0;font-size:9pt;font-family:'Helvetica Neue',Arial,sans-serif}
th,td{border-top:1px solid #ccc;border-bottom:1px solid #ccc;padding:4px 6px;text-align:left;vertical-align:top}
th{background:#f3f3ee} td{font-variant-numeric:tabular-nums}
.tcap{font-size:9pt;font-family:'Helvetica Neue',Arial,sans-serif;margin-top:14px}
.eq{margin:10px 0 10px 30px;font-family:'Cambria Math','STIX Two Math','Times New Roman',serif;font-size:11pt}
.eq span.n{float:right;color:#555}
.note{font-size:9pt;color:#555}
code{font-family:Menlo,monospace;font-size:8.8pt;background:#f2f2ec;padding:0 3px}
.refs li{font-size:9.2pt;margin:4px 0}
@media (max-width:700px){.page{padding:28px 18px}}
"""

BODY = r"""
<div class="page">
<div class="rh">Structural collocation of precipitation observers &nbsp;|&nbsp; Parrish, R. &nbsp;|&nbsp; parallaxmetrology.com &nbsp;|&nbsp; Preprint draft @@DATE@@, not peer reviewed</div>

<h1>Structural Collocation: Estimating Precipitation Observer Error from Disagreement in Rain-Field Geometry, and Detecting Shared Error Without Ground Truth</h1>
<div class="authors">Russell Parrish</div>
<div class="aff">Parallax Metrology · www.parallaxmetrology.com · ORCID 0009-0008-9781-7995</div>

<div class="abstract"><b>Abstract.</b>
Triple collocation estimates each of three observers' error from their covariances alone, provided their errors are independent.
We apply it to <i>structural</i> coordinates of the rain field (intensity-weighted centroid and dispersion in 1.5° windows), using
hourly MRMS radar, GOES-16 10.3 µm brightness temperature, IMERG passive-microwave estimates and rain gauges over the contiguous
United States in July 2024. The rankings depend on which observers form the triplet. With radar, infrared and microwave, microwave
is rated least noisy (centroid y: @@A2_M@@ vs radar @@A2_R@@), robustly to a close smoothness match, two alternative window grids, core-only
reading and radar beam quality. At gauge locations, replacing either satellite with gauges reverses it (radar @@T7_NR@@, gauges
@@T7_NG@@, microwave @@T7_NM@@). Vanishing-tetrad checks give a pattern consistent with a single split, {radar, gauges}
versus {GOES, microwave}. Read as an error shared by the two satellites relative to surface rain, its estimated correlation is
@@T7_CW@@ to @@T7_CS@@. Such a dependence is physically plausible, since both satellites respond to ice and cloud aloft rather than to rain at
the surface; the contribution is a procedure that detects, localises and bounds it without ground truth. Four observers cannot exclude the mirror reading (correlated radar and gauge errors); a sensitivity analysis
shows the satellite dependence stays detected for an assumed radar–gauge error correlation up to 0.475–0.5, while radar's lead in the ranking is clear only up to 0.275–0.4; the actual value is not measured. A synthetic control with known truth and a
mass-conserving blur recovers centroid and dispersion with correct rankings, reproduces the reversal when a shared satellite error is
planted, and, with five observers, identifies the shared pair in either direction under the simulation's assumptions. Structural disagreement did not predict hourly radar
error at gauges (pre-registered negative result). Under a frozen pipeline, all six pre-registered predictions replicated on an
unseen month (August 2024). A real fifth observer of surface rain, and other seasons and regions, remain.
<div class="kw"><b>Keywords:</b> triple collocation; quadruple collocation; vanishing tetrads; precipitation verification; MRMS;
IMERG; GOES; observation error; independence; synthetic validation</div></div>

<figure><img src="@@FIG_THREAD@@" alt="summary">
<figcaption><b>Figure 1.</b> The result in one picture. (A) Radar minus microwave noise fraction (centroid y) across the
experiments: whenever the two satellites are two of the three observers, collocation rates radar worse; with gauges in place of a
satellite the sign reverses. Error bars, where shown, are 95% cluster-bootstrap intervals. (B) On synthetic fields with known truth (corrected
control), a planted error shared by the two satellite-like observers drives the estimate away from the truth; at
σ<sub>sh</sub> = 1 cell the estimated ranking inverts.</figcaption></figure>

<h2>1. Introduction</h2>
<p>Precipitation analyses from radar, geostationary infrared and passive microwave disagree, and where gauges are sparse there is
no direct way to say which is closer to the rain at the ground. Triple collocation (TC; Stoffelen 1998) resolves this in principle:
if three observers measure the same signal with mutually independent errors, each observer's error variance follows from the three
pairwise covariances. Extended TC adds each observer's correlation with the unknown truth (McColl et al. 2014). TC has been applied
to precipitation amounts over the United States with a multiplicative error model (Alemohammad et al. 2015), to MRMS, gauge-only and
IMERG estimates of extreme rainfall (Li et al. 2020), and in regime-binned form (Pellet 2022). Its central weakness is well documented:
when two observers share an error, TC underestimates their errors and overestimates the third. Extensions handle more than three
observers and errors correlated within groups of products (Pan et al. 2015; Gruber et al. 2016).</p>
<p>Two gaps motivate this work. First, TC has been applied to rain <i>amounts</i>. Much of what makes satellite and radar
precipitation fields disagree is <i>geometric</i>: where rain is placed, how it is spread, how many rain objects there are. Spatial
verification methods for forecasts (SAL, MODE, FSS, optimal transport; Wernli et al. 2008; Davis et al. 2006; Roberts and Lean
2008; Francis et al. 2025) measure such geometry against a reference, but they need a reference. We found no application of
collocation to structural coordinates; the nearest work treats displacement and representation error within the collocation
framework (Jacques 2025), not the geometry of the field as the collocated quantity. Second, the independence assumption is usually
argued from instrument physics ("different sensors") rather than tested. Here a dependence between the two satellites is physically
expected: both respond to ice and cloud aloft, and passive-microwave retrievals over land rely largely on the scattering signal of ice
above the freezing level (Spencer et al. 1989). What is not established is whether such a dependence can be seen, and bounded, from the
observers alone. We show a case where physically different sensors are not independent observers, how the error appears as a
plausible but wrong ranking, and how it can be detected and localised.</p>
<p><b>Contributions.</b> (i) Structural collocation on centroid and dispersion, validated on a synthetic control with a mass-conserving
blur (and a record of how a non-conserving blur created spurious shared errors); threshold counts unvalidated and excluded from interpretation. (ii) A real-data case in which a radar–infrared–microwave triplet produces a
confident inverted ranking; four candidate explanations were tested pre-registered and did not account for it within the tested variants (the IMERG calibration chain,
smoothness with a close match, the window grid (tested shift and size), convective-core keying) and two exploratory checks found the same (geostationary parallax, radar beam quality). (iii) A ranking reversal when
a satellite is replaced by rain gauges, robust to which satellite is replaced, and a four-observer dependence pattern consistent with
a shared satellite error, with its identification limit stated, replicated on an unseen month under a frozen pipeline. (iv) A sensitivity analysis bounding
the attribution by the radar–gauge error dependence it tolerates, and a five-observer design validated in simulation. (v) Six
practical rules, scoped to the tested configuration. A pre-registered negative result on predicting point radar error
from structural disagreement is reported in Appendix B.</p>

<h2>2. Data</h2>
<p>All data cover the contiguous United States (CONUS), hourly accumulations ending 2024-07-01 01:00 UTC to 2024-07-22 00:00 UTC
(504 hours; 502 for analyses using GOES, see Appendix A), and, for the replication in Section 4.7, hourly accumulations ending
2024-08-01 01:00 UTC to 2024-08-22 00:00 UTC (504 hours, all with GOES), fetched only after the replication plan was frozen. Every field is placed on the 0.1° MRMS-aligned grid. Table 1 lists the
observers.</p>
<div class="tcap"><b>Table 1.</b> Observers.</div>
<table><tr><th>Observer</th><th>Physical basis</th><th>Processing here</th><th>Source</th></tr>
<tr><td>R: MRMS radar-only QPE (1 h)</td><td>S-band reflectivity → rain rate; not gauge-adjusted (Zhang et al. 2016)</td><td>0.01° → 0.1° block mean (≥ 80% valid)</td><td>NOAA Open Data (AWS)</td></tr>
<tr><td>RQI: MRMS radar quality index</td><td>beam height, blockage, freezing level</td><td>window mean</td><td>same</td></tr>
<tr><td>GIR: GOES-16 ABI band 13</td><td>10.3 µm brightness temperature Tb (cloud top; Schmit et al. 2017)</td><td>4 scans/h averaged, binned to 0.1°; field 235 K − Tb, floored at 0</td><td>NOAA Open Data (AWS)</td></tr>
<tr><td>MW: IMERG Late V07 <code>MWprecipitation</code></td><td>passive microwave (emission/scattering)</td><td>mean of available half-hours</td><td>NASA GES DISC</td></tr>
<tr><td>IR: IMERG Late V07 <code>IRprecipitation</code></td><td>geostationary IR → rain (PDIR), calibrated to MW</td><td>mean of half-hours</td><td>same</td></tr>
<tr><td>M: IMERG Late V07 merged</td><td>Kalman-morphed MW + IR</td><td>mean of half-hours</td><td>same</td></tr>
<tr><td>G: gauges</td><td>surface catch: USCRN (130) and ASOS (2,175) hourly</td><td>label (App. B) and point observer (§4.4)</td><td>NCEI; Iowa Environmental Mesonet</td></tr></table>
<p>None of R, GIR, MW uses gauges at the time of observation; the IMERG infrared estimate is calibrated to the microwave estimate
through a 30-day rotating accumulator, and the microwave estimates are quantile-matched to the GPM combined radar–radiometer
product (Huffman et al. 2023). Figure 2 shows one hour through each observer.</p>
<figure><img src="@@FIG_OBS@@" alt="observers">
<figcaption><b>Figure 2.</b> One hour (ending 2024-07-21 23:00 UTC) through each observer on the common grid; bottom right, the
1.5° windows entering the (radar, GOES, microwave) sample that hour.</figcaption></figure>

<h2>3. Method</h2>
<h3>3.1 Windows and structural reads</h3>
<p>The domain is tiled with non-overlapping 15 × 15-cell (1.5°) windows. A window-hour enters when radar is ≥ 90% valid, the
infrared field is complete, the microwave field covers ≥ 80% of the window, and every observer's window sum is ≥ 1. For each
observer field A over window coordinates (y, x) centred on the window, we compute the intensity-weighted centroid and the
dispersion:</p>
<div class="eq">c = Σ A(y,x)·(y,x) / Σ A(y,x), &nbsp;&nbsp; d = [ Σ A·‖(y,x) − c‖² / Σ A ]<sup>½</sup> <span class="n">(1)</span></div>
<p>We also tested the wet-area fraction (cells ≥ 0.5 mm h⁻¹, or cold cloud for GIR), the number of rain objects (8-connected
components of ≥ 2 cells) and log mean rain as an amount-based control.</p>
<h3>3.2 Collocation estimators</h3>
<p>Under the additive model x<sub>i</sub> = a<sub>i</sub> + b<sub>i</sub>T + e<sub>i</sub> with errors independent of each other
and of T, the ETC estimates for a triplet (i, j, k) with covariance matrix C are</p>
<div class="eq">σ²<sub>e,i</sub> = C<sub>ii</sub> − C<sub>ij</sub>C<sub>ik</sub>/C<sub>jk</sub>, &nbsp;&nbsp;
ρ²<sub>i</sub> = C<sub>ij</sub>C<sub>ik</sub> / (C<sub>ii</sub>C<sub>jk</sub>), &nbsp;&nbsp; noise fraction ν<sub>i</sub> = 1 − ρ²<sub>i</sub> <span class="n">(2)</span></div>
<p>ν is written <i>nf</i> in the figures. Positions are in grid cells of 0.1° (about 10 km).</p>
<p>With a fourth observer, the error covariance of a pair (p, q) can be estimated using two observers (u, v) assumed independent of
everything else as instruments:</p>
<div class="eq">c<sub>pq</sub> = C<sub>pq</sub> − C<sub>pu</sub>C<sub>qv</sub>/C<sub>uv</sub>, normalised by (σ²<sub>e,p</sub>σ²<sub>e,q</sub>)<sup>½</sup> <span class="n">(3)</span></div>
<p>A one-factor model with independent errors also implies vanishing tetrads among four observers (Spearman and Holzinger 1924):</p>
<div class="eq">t<sub>1</sub> = r<sub>RG</sub>r<sub>IM</sub> − r<sub>RI</sub>r<sub>GM</sub>, &nbsp; t<sub>2</sub> = r<sub>RG</sub>r<sub>IM</sub> − r<sub>RM</sub>r<sub>GI</sub>, &nbsp; t<sub>3</sub> = r<sub>RI</sub>r<sub>GM</sub> − r<sub>RM</sub>r<sub>GI</sub> <span class="n">(4)</span></div>
<p>with R radar, I the infrared observer, M microwave and G gauges.</p>
<p>A single shared error between I and M inflates r<sub>IM</sub> only, so t<sub>3</sub> ≈ 0 while t<sub>1</sub>, t<sub>2</sub> &gt; 0.
The same pattern arises from a shared error between R and G; four observers cannot separate the two (Section 5.3). By
construction t<sub>2</sub> − t<sub>1</sub> = t<sub>3</sub>, so the three tetrads carry two independent constraints, not three.</p>
<h3>3.3 Uncertainty and pre-registration</h3>
<p>Intervals are 95% cluster-bootstrap percentiles, resampling (UTC day × 5° tile) blocks (300–500 resamples, fixed seeds).
Every experiment below was specified, with its predictions and decision rule, in a plan file frozen before any value of that
experiment was computed; the Appendix B plan received two amendments, both written before any target or predictor value existed.
Follow-up splits added after scoring are labelled exploratory.</p>
<h3>3.4 Synthetic anchor</h3>
<p>We generated 20,000 fields of known truth on a 64 × 64-cell canvas (1 + Poisson(1.5) Gaussian rain cells, widths U(0.8, 2.5),
lognormal peak intensities). Simulated observers re-render the same cells with independent per-cell errors: radar-like R (position
sd 0.5 cell, no blur, intensity noise lognormal 0.3), microwave-like S1 (sd 0.8, blur 1 cell, 0.4), infrared-like S2 (sd 1.5, blur 2,
0.6) and, for the five-observer benchmark, a gauge-like G<sub>s</sub> and a fifth observer. Blur is a mass-conserving Gaussian
convolution, which leaves the centroid unchanged and adds 2σ² to the squared dispersion. Reads are taken on the full canvas and in
central and shifted 15 × 15 windows. A shared displacement (one per field, N(0, σ<sub>sh</sub>) per axis) can be added to S1 and S2,
or to R and G<sub>s</sub> (the mirror case). True noise fractions and coordinate-error correlations are computed against the known
truth reads. A first version of this simulation widened cells with fixed peaks, which does not conserve rain mass; it produced
artefacts and is superseded (Section 4.1, Appendix A).</p>

<h2>4. Results</h2>
<h3>4.1 Synthetic validation</h3>
<p>The rendering behaved as theory requires (centroid unchanged to 10⁻¹⁵ cell; squared dispersion increased by 2σ² within 0.04%).
With nothing planted, the satellite-like observers' centroid errors were uncorrelated on the full canvas (centroid y
@@TB_BC_CY@@), while dispersion errors kept a small correlation (+0.021, within the pre-registered |r| &lt; 0.05 but not zero); windows added a small shared error (centroid up to @@TB_WIN_CY@@, dispersion up to @@TB_WIN_D@@), growing as the storm
sat further from the window centre (Figure 3A). Collocation recovered centroid and dispersion with the correct ranking on the full
canvas and in every window (maximum error @@TB_CY@@ centroid, @@TB_D@@ dispersion).</p>
<p>A planted shared satellite error moved the estimates toward inversion: at σ<sub>sh</sub> = 0.5 cell, R @@TB_05R@@ and S1 @@TB_05S@@
(truths @@TB_05TR@@, @@TB_05TS@@), not yet crossed; at σ<sub>sh</sub> = 1.0, R @@TB_E1R@@ and S1 @@TB_E1S@@ (truths @@TB_E1TR@@,
@@TB_E1TS@@), inverted (Figure 3B). With five observers, the tetrad audit flagged a planted shared satellite pair and the mirror
radar–gauge pair in 10 of 10 blocks each, with no false alarms when nothing was planted.</p>
<p><i>First simulation (superseded).</i> The non-conserving first version reported a baseline correlation between the two blurred
observers' errors (@@T1_BASE_EC@@ on centroid y), and a dispersion mis-ranking, both artefacts of the blur, and an inversion already at
σ<sub>sh</sub> = 0.5, which the artefact helped produce. Its findings for the variables the corrected control did not re-test still stand as provisional: the
object count failed (errors up to @@T1_ISL_ERR@@), wet fraction gave impossible negative error variances in @@T1_WET_NEG@@% of
resamples, and log mean rain was stable but biased by up to @@T1_LOG_ERR@@. Threshold counts are therefore not interpreted below.</p>
<figure><img src="@@FIG_TB@@" alt="synthetic control and sensitivity">
<figcaption><b>Figure 3.</b> (A) Corrected synthetic control: true coordinate-error correlation between the satellite-like observers
with nothing planted, full canvas and windows; dashed, the value from the superseded first simulation. (B) A planted shared satellite
error inverts the estimated ranking at σ<sub>sh</sub> = 1.0. (C) Sensitivity of the real-data attribution (Section 4.6): inferred
GOES–microwave error correlation (solid) and radar minus microwave noise fraction (dotted) against the assumed radar–gauge error
correlation.</figcaption></figure>

<h3>4.2 Real data: an inverted ranking</h3>
<p>On @@A_N@@ window-hours, the triplet (radar, IMERG IR, MW) gave centroid-y noise fractions R @@A_R@@, IR @@A_I@@, MW
@@A_M@@: microwave rated least noisy, radar second. Because the IMERG infrared estimate is calibrated to the microwave, we replaced it
with raw GOES-16 brightness temperature (@@A2_N@@ window-hours): R @@A2_R@@, GIR @@A2_I@@, MW @@A2_M@@ (Figure 4). The
ranking was unchanged on centroid and dispersion (the variables validated in §4.1; wet fraction and object count are shown in
Figure 4 for completeness but are not interpreted). The amount-based control was unstable in both triplets (impossible variances in
@@A_LOGNEG@@% and @@A2_LOGNEG@@% of resamples) although it was stable on synthetic fields, so its instability reflects the real
observers' amount relationships rather than the estimator. With three observers the ranking follows from the raw correlations
(centroid y: r<sub>R,GIR</sub> = @@A2_RRI@@, r<sub>R,MW</sub> = @@A2_RRM@@, r<sub>GIR,MW</sub> = @@A2_RIM@@): microwave is
rated best because the infrared field agrees more with microwave than with radar. Clusters: @@A_CLU@@ (A) and @@A2_CLU@@ (A2).
The merged IMERG product, used as a fourth observer in the IMERG-IR triplet, showed a positive error correlation with its microwave
input (@@A_MMW_LO@@ to @@A_MMW_HI@@ on the stable variables), a dependence that exists by construction; in the first
simulation the magnitude of such estimates was unreliable when the instruments were themselves dependent, and the corrected control did
not re-test it, so it is not interpreted further.</p>
<figure><img src="@@FIG_A@@" alt="A A2">
<figcaption><b>Figure 4.</b> Noise fractions per observer and variable for the (R, IMERG IR, MW) and (R, GOES IR, MW) triplets,
with 95% intervals; "unstable" marks variables failing the pre-registered stability rule.</figcaption></figure>

<h3>4.3 Alternative explanations tested</h3>
<p><i>Smoothness.</i> Radar was smoothed with a Gaussian before windowing, over a ladder up to σ = 3 cells. Microwave's median lag-1
autocorrelation was @@T5_LMW@@; radar's was @@T5_LR0@@ unsmoothed and @@T5_LSTAR@@ at the nearest ladder step, σ = @@T5_STAR@@ cell, so the
match was approximate. The gap ν<sub>R</sub> − ν<sub>MW</sub> remained @@T5_GAP@@ @@T5_GAPCI@@ at that step and @@T5_GAP3@@ at σ = 3
(Figure 5A); the sample also changed along the ladder (@@T5_N0@@ to @@T5_N3@@ window-hours). A pre-registered repeat on a fixed sample of @@T9_N@@ window-hours, with
radar blurred to match microwave's smoothness to within 0.013 in lag-1 autocorrelation (σ = @@T9_STAR@@; lag-1 @@T9_LR@@ vs microwave @@T9_LM@@), gave @@T9_G@@ @@T9_GCI@@, essentially
unchanged, and @@T9_G3@@ at σ = 3. The synthetic calibration in Figure 5B comes from the superseded first simulation and is shown for the record only.</p>
<p><i>Window grid.</i> Re-tiling with a half-shifted grid and with 25 × 25 windows left the gap at @@T9_W1@@ and
@@T9_W2@@ (baseline @@A2_GAP@@), and the point-support reversal at @@T9_P1@@ and @@T9_P2@@ (pre-registered). This shows robustness to the tested grid shift and window size; it does not exclude windowing or
representativeness effects in general (the 25 × 25 gap is about a fifth smaller).</p>
<p><i>Parallax.</i> GOES cold-cloud centroids lie increasingly north of radar rain with latitude (@@PAR@@ cells, about 10 km, at 45–55°N), the
expected geostationary parallax; microwave does not follow (@@PARMW@@), so parallax is not shared (Figure 5C, exploratory).</p>
<figure><img src="@@FIG_T5@@" alt="T5">
<figcaption><b>Figure 5.</b> (A) Radar minus microwave noise fraction as radar is smoothed; dashed line, the matched smoothness.
(B) Synthetic calibration from the superseded first simulation (record only). (C) Mean centroid offset from radar by latitude (exploratory).</figcaption></figure>
<p><i>Convective-core keying.</i> Reading every observer on its own core (cells at or above its 80th percentile) left the gap at
@@T6_CORE@@ (full area @@T6_FULL@@). The gap was largest in the weakest windows (@@T6_Q0@@ in the lowest peak-intensity quintile vs
@@T6_Q4@@ in the highest; exploratory).</p>
<p><i>Radar beam quality.</i> The gap did not decrease with radar quality index; in the best-beam windows (RQI ≈ @@T6_RQIM@@) with weak
rain it was @@T6_RQIW@@ (Figure 6, exploratory).</p>
<figure><img src="@@FIG_T6@@" alt="T6">
<figcaption><b>Figure 6.</b> (A) Full-area vs core-only reading. (B) The gap by radar peak-intensity quintile. (C) The gap by radar
quality index and intensity. B and C are exploratory.</figcaption></figure>

<h3>4.4 A sharp, independent fourth observer</h3>
<p>A three-observer ranking cannot reveal a shared error among its legs, so we added a fourth observer, rain gauges,
which none of the other observers uses. We restricted to window-hours with ≥ 5 reporting gauges (@@T7_N@@ window-hours, median @@T7_G@@ gauges) and read
every observer at the gauge locations (point support), splitting at the median radar peak. With (R, GIR, MW) at point support the
inversion reappeared (@@T7_P0@@ @@T7_P0CI@@, weak half). Replacing GOES with gauges reversed it: (R, MW, G) gave
@@T7_P1@@ @@T7_P1CI@@, with noise fractions radar @@T7_NR@@, gauges @@T7_NG@@, microwave @@T7_NM@@ (Figure 7). The GOES–microwave
error correlation with radar and gauges as instruments was @@T7_CW@@ @@T7_CWCI@@ (weak) and up to @@T7_CS@@ (strong, centroid x).</p>
<p>A pre-registered robustness check replaced the other satellite: in (R, GIR, G) radar remained least noisy (@@T7R_RG@@ vs GOES
@@T7R_IG@@), and radar's noise fraction did not depend on which satellite was dropped (difference @@T7R_R2@@ @@T7R_R2CI@@). The
tetrads (Eq. 4) gave t<sub>1</sub> = @@T7R_T1@@ @@T7R_T1CI@@, t<sub>2</sub> = @@T7R_T2@@ @@T7R_T2CI@@ and
t<sub>3</sub> = @@T7R_T3@@ @@T7R_T3CI@@ (weak, centroid y). Since t<sub>3</sub> = t<sub>2</sub> − t<sub>1</sub>,
this is two independent constraints, and it is compatible with a single two-group split; an interval containing zero for
t<sub>3</sub> does not establish that no other dependence exists. All pre-registered predictions held in both strata and both axes.
The estimated GOES–microwave correlations above use radar and gauges as instruments and are conditional on their errors being
independent (Section 5.3).</p>
<figure><img src="@@FIG_T7@@" alt="T7">
<figcaption><b>Figure 7.</b> Noise fractions at gauge points for the (R, GOES, MW) and (R, MW, gauges) triplets, weak and strong
halves; (C) GOES–microwave error correlation with radar and gauges as instruments.</figcaption></figure>

<h3>4.5 Where the coupling grows</h3>
<p>We tested whether the coupling grows with the share of cold cloud lying over no surface rain (anvil fraction). No change was detected: the
normalised GOES–microwave error correlation was @@T8_C0@@ in the lowest and @@T8_C2@@ in the highest anvil tercile (difference
@@T8_TB@@ @@T8_TBCI@@). The effect of the coupling on the three-observer ranking did grow, from @@T8_I0@@ to @@T8_I2@@. In absolute
terms the shared covariance roughly doubled (@@T8_COV0@@ → @@T8_COV2@@ cells²), but so did every observer's error variance,
gauges included (@@T8_G0@@ → @@T8_G2@@); anvil fraction correlates @@T8_CORR@@ with radar peak intensity (Figure 8).</p>
<figure><img src="@@FIG_T8@@" alt="T8">
<figcaption><b>Figure 8.</b> (A) Coupling of the satellite errors by anvil tercile. (B) Three-observer gap by anvil tercile.
(C) Absolute error variances and shared covariance (exploratory).</figcaption></figure>

<h3>4.6 Sensitivity to radar–gauge dependence</h3>
<p>We fitted a one-signal model to the four observers at gauge points with an assumed radar–gauge error correlation
ρ<sub>RG</sub>, swept from 0 to 0.5, and a free GOES–microwave error correlation, by maximum likelihood. The fit is identical
at every ρ<sub>RG</sub>, so the data cannot fix it; the sweep traces the equivalence class. We report two outcomes separately, by
interval. The inferred satellite dependence (GOES–microwave error correlation, interval excluding zero) stayed detected up to
ρ<sub>RG</sub> = 0.475–0.5 across strata and axes in July, and through 0.5, the top of the tested range, in August. Radar's ranking as
least noisy (interval of ν<sub>R</sub> − ν<sub>MW</sub> excluding zero) held up to 0.275–0.4 in July and 0.375–0.4 in August; beyond that
the interval includes reversal. The point estimate of the ranking reversed at @@TS_R2MIN@@–@@TS_R2MAX@@ in July and had not reversed by
0.5 in August (gap @@X10_TS@@ at 0.5, interval spanning zero) (Figure 3C). Neither month tests values above 0.5. In both months this
bounds how far the assumption can move; it does not measure the actual correlation.</p>

<h3>4.7 Replication on an unseen period</h3>
<p>Every real-data result above is in-sample (July). We froze the July pipeline (script hashes, inclusion rules, strata, seeds, six primary
predictions and a scoring rule) before fetching any data for August 2024, then ran it unchanged. Before the August run, the July
outputs were regenerated with the period switch in place and reproduced byte-for-byte; an independent re-implementation written from
the plan files alone also reproduced the July quantities to within 4 × 10⁻⁹. August gave @@X10_NA@@ area and @@X10_NP@@ point-support
window-hours. All six predictions held on both axes and in every required stratum (Table 2): the area inversion (@@X10_A@@ on centroid
y), the reversal when gauges replace GOES (@@X10_P@@, weak stratum), the mirror triplet, the GOES–microwave error correlation
(@@X10_C@@) and the tetrad pattern. Magnitudes were similar or somewhat larger than in July. The effect is therefore not specific
to July; the replication does not change what four observers can identify, and both months are adjacent and warm-season. Replication across adjacent months with the same instruments, processing and season shows repeatability, not attribution: a persistent measurement or representativeness effect would replicate too.</p>
<div class="tcap"><b>Table 2.</b> Pre-registered replication (95% cluster-bootstrap intervals). Δ = ν<sub>R</sub> − ν<sub>MW</sub>.
Verdict: @@X10_V@@.</div>
<table><tr><th></th><th>quantity</th><th>July 2024 (in-sample)</th><th>August 2024</th><th>met</th></tr>
@@TAB2@@</table>

<h2>5. Discussion</h2>
<h3>5.1 What is new and what is not</h3>
<p>Triple collocation, its sensitivity to shared error, and vanishing tetrads are established. Satellite precipitation products are
known not to be fully independent (e.g. Tian and Peters-Lidard 2010 use inter-product spread as uncertainty, acknowledging shared
features). The contributions here are narrower: collocation on structural coordinates, with a synthetic validation, on a mass-conserving control, of which
structural variables are admissible; a worked case in which physically different sensors fail the independence premise, where the
failure surfaces only as a plausible ranking; and a detection and localisation procedure (sharp independent observer, mirror
triplet, tetrads) that turns the failure into a measurement. Read most directly, the nonzero tetrads in both months reject a
single common signal with mutually independent errors for these four observers; the satellite attribution is one reading of that
failure. Two things were not tested: whether geometric coordinates reveal dependence that amount-based or representativeness analyses
miss, and whether the procedure, used as a warning that a collocation ranking is not trustworthy, improves on existing diagnostics.</p>
<h3>5.2 Implications for satellite-only evaluation</h3>
<p>Collocation of satellite products against each other is attractive exactly where ground observations are absent. In this sample,
the four-observer pattern is consistent with an error shared by an infrared and a microwave observer relative to surface rain, with an estimated normalised strength near 0.4 and no detectable change
across anvil extent; its effect on rankings was largest where all observers were noisiest. If that reading holds, collocation using two
satellites as two of three legs can flatter both in such conditions. We did not establish the physical source. Both satellites
respond to cloud and hydrometeors aloft rather than to rain at the surface (over land, passive-microwave retrievals lean on ice
scattering; Spencer et al. 1989); a coupling that does not vary detectably with anvil extent is compatible with
both observing the cloud system, and the simple anvil form was not supported. That the two groups respond to different physical
quantities (surface rain versus aspects of the cloud system) is a hypothesis; processing, resolution and sampling are alternative
explanations.</p>
<h3>5.3 Identification</h3>
<p>With four observers, a shared error between GOES and microwave relative to surface rain is observationally equivalent to a shared
error between radar and gauges relative to a cloud-system signal. Taking surface rain as the target gives the analysis its
purpose but does not show that radar and gauge errors are uncorrelated, so the satellite attribution and the estimated correlations
are conditional on that assumption. Vogelzang and Stoffelen (2021) make the same points for wind collocation: resolution
(representativeness) differences generate error covariance, and a covariance solution that fits does not by itself establish a correct
error model. A fifth observer would help only under explicit assumptions about its own signal and errors; different physics alone does
not guarantee identification. The sensitivity analysis (Section 4.6) quantifies this: the satellite dependence stays detected for radar–gauge error
correlations up to 0.475–0.5 (through the tested maximum of 0.5 in August), and radar's ranking is clear up to 0.275–0.4. It shows how far the assumption can move before the ranking changes, not where the actual correlation lies. In simulation, five
observers resolved the ambiguity in either direction (Section 4.1), under that simulation's assumptions about the fifth observer's signal
and errors; a real fifth observer would need those assumptions justified. The choice matters. Lightning (for example the GOES
Geostationary Lightning Mapper) is a poor candidate: charge separation requires ice aloft, so it would likely share the satellites'
dependence and leave the two readings unseparated. The natural candidate senses rain at the surface with errors independent of radar
and gauges, such as disdrometers or areal averages from dense surface networks, where they exist.</p>
<h3>5.4 Practical rules</h3>
<p>These are findings for the tested variables, simulation and observing configuration, not general laws.</p>
<ol>
<li>Prefer geometric coordinates. With a correct blur, centroid and dispersion were recovered and ranked correctly; threshold counts
failed in the superseded simulation, were not re-tested, and are excluded from interpretation.</li>
<li>A shared error among the legs is invisible to the estimator; it appears only as a plausible but wrong ranking.</li>
<li>Check independence after feature extraction and windowing, and audit the simulation itself: a blur that did not conserve mass
created shared centroid errors that a correct blur did not; a small dispersion correlation remained under the correct blur, and
windows add a small shared error. Validating the rendering and validating the measurements' error assumptions are separate checks.</li>
<li>Judge independence by the quantity each observer responds to, and check derived products for calibration chains.</li>
<li>The bias from a shared error is largest where all observers are noisy.</li>
<li>Validate on synthetic fields with known truth before reading a real-data ranking.</li></ol>
<h3>5.5 Limitations</h3>
<p>Two adjacent months (warm-season convection; July, and the August replication) over CONUS, hourly, 0.1°; the replication shows
repeatability under the same instruments, processing and season, not attribution. Several splits are exploratory and uncorrected for multiple
comparisons. Microwave coverage ≥ 80% selects hours near overpasses; the gauge sample requires rain at the gauge points and is wetter
than the full sample. The synthetic model is a simple cell model (corrected version in Section 4.1). Point-support centroids from about nine gauges per window are coarse. Collocation
assumes each observer is linear in a common signal; for the cold-cloud field, whose centroid tracks cloud rather than rain, this is an
approximation, and the four-observer checks are the evidence that it does not drive the results.</p>

<h2>6. Conclusions</h2>
<p>Structural collocation produces observer rankings that depend on which observers form the triplet. With radar, GOES infrared and
passive microwave it rated microwave least noisy; replacing either satellite with gauges reversed that, at the same point support. The
four-observer results show a dependence pattern consistent with an error shared by the two satellites relative to surface rain, under
the assumption that radar and gauge errors are not strongly correlated (satellite dependence detected up to 0.475–0.5 (July) and through 0.5, the top of the tested range (August); radar ranked least noisy with intervals excluding reversal up to 0.275–0.4 (July) and 0.375–0.4 (August); the actual value is not measured). The results replicated on an unseen month under a frozen pipeline. A corrected
synthetic control supports the variable choices and the mechanism, and shows that five observers can resolve the remaining ambiguity under its assumptions. The practical lesson is to test the
independence premise of a collocation, after feature extraction, before reading its ranking.</p>

<h2>7. Data and Resources</h2>
<p>MRMS: NOAA Open Data Dissemination, <code>s3://noaa-mrms-pds</code>. GOES-16 ABI L2 CMIPC band 13:
<code>s3://noaa-goes16</code>. IMERG Late V07 half-hourly (<code>GPM_3IMERGHHL.07</code>): NASA GES DISC (Earthdata login). USCRN
hourly: NCEI. ASOS hourly precipitation: Iowa Environmental Mesonet. Code, frozen plans, results and the full chronological lab book
(<code>LABBOOK.html</code>, including dead ends, amendments and file hashes) are in the project folder
<code>precip-disagreement/</code>; the four-observer tetrad audit is <code>tetrad_audit</code> in <code>independence_audit.py</code>.
These are available from the author on request; [a public repository link to be added before submission].</p>

<h2>Acknowledgements</h2>
<p>Analysis code, figures and a first draft of this text were produced with the assistance of an AI model (Claude, ChatGPT, Gemini)
working under the author's direction; the author is responsible for the content.</p>

<h2>Appendix A. Transparency and Reproducibility</h2>
<p>Each experiment has a plan file frozen before computation (<code>FEASIBILITY_PLAN.md</code>, <code>PROBE_A_PLAN.md</code>,
<code>PROBE_A2_PLAN.md</code>, <code>T1_PLAN.md</code>, <code>T5_PLAN.md</code>, <code>T6_PLAN.md</code>, <code>T7_PLAN.md</code>,
<code>T7R_PLAN.md</code>, <code>T8_PLAN.md</code>, <code>T1BIS_PLAN.md</code>, <code>T7S_PLAN.md</code>, <code>T9_PLAN.md</code>, <code>T10_PLAN.md</code>) and a results file. Deviations: the Appendix B experiment excluded gauges within 1° of the domain
edge; two GOES hours had fewer than three scans and were dropped per plan. One prediction premise was wrong and is reported as such
(the synthetic S1–S2 error correlation was not zero with nothing planted; the corrected control showed this was an artefact of the
non-conserving blur, not of smoothing, and the claim built on it was withdrawn). An interpretation was later
corrected: the tetrad pattern was first described as "one shared GOES–microwave error"; the identification limit in §5.3 was found
afterwards. An external review identified that the synthetic blur did not conserve rain mass, that the
synthetic dispersion ranking was wrong with nothing planted, and that the three tetrads are not independent; the text was corrected
accordingly and the corrected synthetic control (T1-bis) and a radar–gauge sensitivity analysis (T7-S) were pre-registered and run; their results
replace the provisional statements. Before the replication, an independent re-implementation from the plan files reproduced the July
quantities and recorded four details fixed by the code rather than the plans (≥ at the wet threshold, float32 window sums, T7 not
requiring the merged product, the median tie rule). The replication's deviations were mechanical (data root chosen by working
directory, a missing stations-table producer added and checked byte-identical on July, the Δ interval computed from the existing
bootstrap draws, a single-download IMERG fetcher checked identical on July) and are listed in <code>T10_RESULTS.md</code>. The first simulation's results are kept in the lab book for the record; the paper reports them only
where the corrected control did not re-test a variable (Section 4.1). All raw files were streamed and reduced in memory; reduced arrays and their digests are listed in the lab book.</p>

<h2>Appendix B. A negative result: structural disagreement does not predict point radar error</h2>
<p>The work began with a different question: does structural disagreement between radar and satellite predict where radar is wrong
at a rain gauge? On @@E1_N@@ gauge-hours, fourteen structural reads (location, pattern, overlap, object count, dispersion, void,
multiscale, edge, persistence, spatial and temporal variance of disagreement, and an infrared third observer) were added to a
baseline (MRMS quality index, plain differences, magnitudes). The validation increment was Δρ = @@E1_D@@ @@E1_CI@@ against a
pre-registered threshold of 0.02 (baseline ρ = @@E1_BASE@@), and a spatial-scramble null gave a larger increment (@@E1_NULL@@).
Hourly point error is dominated by rain amount and by the area–point representativeness term, as radar-verification studies predict
(Ciach and Krajewski 1999; Ciach et al. 2003).</p>
<figure><img src="@@FIG_E1@@" alt="E1">
<figcaption><b>Figure B1.</b> (A) Increment over baseline, real and scramble null. (B) Partial correlations of each read with
radar error after controlling for magnitude. (C) Resolution sweep.</figcaption></figure>

<h2>References</h2>
<ol class="refs">
<li>Alemohammad, S. H., McColl, K. A., Konings, A. G., Entekhabi, D., and Stoffelen, A. (2015). Characterization of precipitation product errors across the United States using multiplicative triple collocation. <i>Hydrol. Earth Syst. Sci.</i>, 19, 3489–3503.</li>
<li>Allison, L., Mittermaier, M., Scovell, R., Vernon, C., Norman, K., and Pierce, C. (2015). Quantifying the effect of radar observational uncertainty on the verification of kilometer-scale NWP precipitation forecasts using the Fractions Skill Score. EMS Annual Meeting Abstracts, EMS2015-75.</li>
<li>Ciach, G. J., and Krajewski, W. F. (1999). On the estimation of radar rainfall error variance. <i>Adv. Water Resour.</i>, 22(6).</li>
<li>Ciach, G. J., Habib, E., and Krajewski, W. F. (2003). Zero-covariance hypothesis in the error variance separation method of radar rainfall verification. <i>Adv. Water Resour.</i>, 26(5).</li>
<li>Davis, C., Brown, B., and Bullock, R. (2006). Object-based verification of precipitation forecasts. Part I. <i>Mon. Wea. Rev.</i>, 134.</li>
<li>Francis, J. J. M., Cotter, C. J., and Mittermaier, M. P. (2025). Examining entropic unbalanced optimal transport and Sinkhorn divergences for spatial forecast verification. <i>Meteorol. Appl.</i>, 32(4).</li>
<li>Gruber, A., Su, C.-H., Zwieback, S., Crow, W., Dorigo, W., and Wagner, W. (2016). Recent advances in (soil moisture) triple collocation analysis. <i>Int. J. Appl. Earth Obs. Geoinf.</i>, 45, 200–211.</li>
<li>Guilloteau, C., Foufoula-Georgiou, E., Kirstetter, P., Tan, J., and Huffman, G. J. (2021). How well do multisatellite products capture the space–time dynamics of precipitation? Part I. <i>J. Hydrometeor.</i>, 22.</li>
<li>Habib, E., Ciach, G. J., and Krajewski, W. F. (2004). A method for filtering out raingauge representativeness errors from the verification distributions of radar and raingauge rainfall. <i>Adv. Water Resour.</i>, 27.</li>
<li>Huffman, G. J., et al. (2023). <i>Algorithm Theoretical Basis Document (ATBD), NASA Global Precipitation Measurement (GPM) Integrated Multi-satellitE Retrievals for GPM (IMERG), Version 07</i>. NASA.</li>
<li>Jacques, D. (2025). Revisiting measurement and representation errors in geophysical sciences. <i>Q. J. R. Meteorol. Soc.</i>, 151(773). doi:10.1002/qj.5073.</li>
<li>Li, Z., Chen, M., Gao, S., Hong, Z., Tang, G., Wen, Y., Gourley, J. J., and Hong, Y. (2020). Cross-examination of similarity, difference and deficiency of gauge, radar and satellite precipitation measuring uncertainties for extreme events using conventional metrics and multiplicative triple collocation. <i>Remote Sens.</i>, 12(8), 1258.</li>
<li>McColl, K. A., Vogelzang, J., Konings, A. G., Entekhabi, D., Piles, M., and Stoffelen, A. (2014). Extended triple collocation: estimating errors and correlation coefficients with respect to an unknown target. <i>Geophys. Res. Lett.</i>, 41.</li>
<li>Pan, M., Fisher, C. K., Chaney, N. W., Zhan, W., Crow, W. T., Aires, F., Entekhabi, D., and Wood, E. F. (2015). Triple collocation: beyond three estimates and separation of structural/non-structural errors. <i>Remote Sens. Environ.</i>, 171, 299–310.</li>
<li>Pellet, V. (2022). A binned triple collocation for estimating regime-dependent uncertainties of precipitation. <i>IEEE Geosci. Remote Sens. Lett.</i>, 19. doi:10.1109/LGRS.2022.3162703.</li>
<li>Roberts, N. M., and Lean, H. W. (2008). Scale-selective verification of rainfall accumulations from high-resolution forecasts of convective events. <i>Mon. Wea. Rev.</i>, 136.</li>
<li>Schmit, T. J., Griffith, P., Gunshor, M. M., Daniels, J. M., Goodman, S. J., and Lebair, W. J. (2017). A closer look at the ABI on the GOES-R series. <i>Bull. Amer. Meteor. Soc.</i>, 98.</li>
<li>Spearman, C., and Holzinger, K. (1924). The sampling error in the theory of two factors. <i>Br. J. Psychol.</i>, 15.</li>
<li>Spencer, R. W., Goodman, H. M., and Hood, R. E. (1989). Precipitation retrieval over land and ocean with the SSM/I: identification and characteristics of the scattering signal. <i>J. Atmos. Oceanic Technol.</i>, 6, 254–273.</li>
<li>Stoffelen, A. (1998). Toward the true near-surface wind speed: error modeling and calibration using triple collocation. <i>J. Geophys. Res.</i>, 103(C4), 7755–7766.</li>
<li>Tian, Y., and Peters-Lidard, C. D. (2010). A global map of uncertainties in satellite-based precipitation measurements. <i>Geophys. Res. Lett.</i>, 37, L24407.</li>
<li>Vogelzang, J., and Stoffelen, A. (2021). Quadruple collocation analysis of in-situ, scatterometer, and NWP winds. <i>J. Geophys. Res. Oceans</i>, 126, e2021JC017189.</li>
<li>Wernli, H., Paulat, M., Hagen, M., and Frei, C. (2008). SAL: a novel quality measure for the verification of quantitative precipitation forecasts. <i>Mon. Wea. Rev.</i>, 136.</li>
<li>Zhang, J., et al. (2016). Multi-Radar Multi-Sensor (MRMS) quantitative precipitation estimation: initial operating capabilities. <i>Bull. Amer. Meteor. Soc.</i>, 97.</li>
</ol>
<p class="note">Reference details were compiled from the project's literature folder and searches on 2026-10-15; volume and page
numbers should be checked against the publishers before submission.</p>
</div>
"""

html = BODY
for k, v in V.items():
    html = html.replace(f"@@{k}@@", v)
assert "@@" not in html, [x for x in html.split("@@")[1::2]][:5]
out = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Structural Collocation Preprint</title>
<meta name="description" content="Preprint: triple collocation on rain-field geometry from radar, GOES-16 infrared, IMERG microwave and rain gauges, and detecting shared observer error without ground truth.">
<link rel="canonical" href="https://www.parallaxmetrology.com/Precipitation_Structural_Collocation_PAPER.html"><style>{CSS}</style></head><body>{html}</body></html>"""
(H / "PAPER.html").write_text(out)
print("wrote PAPER.html", round(len(out) / 1e6, 2), "MB")
