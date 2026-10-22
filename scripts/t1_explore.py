"""T1 EXPLORATORY follow-ups (not in T1_PLAN.md), run after the confirmatory scoring.
E1: instrument G with NO added blur: does quadruple collocation then recover the baseline S1-S2
    error correlation that smoothing alone creates (0.16 at sigma_sh = 0)?
E2: merged observer M rendered with its own processing noise (not an exact linear combination):
    does the lineage check (P5) stop being degenerate?
"""
import json
import numpy as np, pandas as pd
import probe_a as P
import t1_synthetic as T

out = {}
wins = T.simulate()

# E1: sharp instrument
T.OBS["G"] = (1.0, 0.0, 0.5)
for sig in [0.0, 1.0]:
    d = T.build(wins, sig)
    q = T.quad_G(d); tr = T.truth_scores(d)
    out[f"E1_sig{sig}"] = {v: {"quadG": q[v], "true": tr[v]["errcorr_true"]["MW_IR"]} for v in ["cy", "cx", "disp"]}
    print("E1", sig, {v: (round(q[v]["point"], 2), round(tr[v]["errcorr_true"]["MW_IR"], 2)) for v in ["cy", "cx", "disp"]})
T.OBS["G"] = (1.0, 1.5, 0.5)

# E2: merged with its own noise
rng2 = np.random.default_rng(99)
own = [dict(dy=rng2.normal(0, 0.5, len(w[0]["cy"])), dx=rng2.normal(0, 0.5, len(w[0]["cy"])),
            m=rng2.lognormal(0, 0.2, len(w[0]["cy"]))) for w in wins]
orig_fields = T.fields


def fields_e2(win, sig, _own={id(w): o for w, o in zip(wins, own)}):
    f = orig_fields(win, sig)
    cells, err, unit = win; o = _own[id(win)]; sh = unit * sig
    cy = 0.6 * (cells["cy"] + err["S1"]["dy"]) + 0.4 * (cells["cy"] + err["S2"]["dy"]) + sh[0] + o["dy"]
    cx = 0.6 * (cells["cx"] + err["S1"]["dx"]) + 0.4 * (cells["cx"] + err["S2"]["dx"]) + sh[1] + o["dx"]
    I = cells["I"] * err["S1"]["m"] ** 0.6 * err["S2"]["m"] ** 0.4 * o["m"]
    f["M"] = T.render(cy, cx, cells["s"] + 1.4, I)
    return f


T.fields = fields_e2
for sig in [0.0, 1.0]:
    d = T.build(wins, sig)
    tc = P.analyse(d, n_boot=200); tr = T.truth_scores(d)
    out[f"E2_sig{sig}"] = {v: {"M_MW": tc[v]["point"]["M_MW"], "M_IR": tc[v]["point"]["M_IR"],
                               "M_MW_ci": tc[v]["ci95"]["M_MW"], "M_IR_ci": tc[v]["ci95"]["M_IR"],
                               "true_M_MW": tr[v]["errcorr_true"]["M_MW"], "true_M_IR": tr[v]["errcorr_true"]["M_IR"]}
                           for v in ["cy", "cx", "disp", "isl"]}
    print("E2", sig, {v: (round(x["M_MW"], 2), round(x["true_M_MW"], 2), round(x["M_IR"], 2), round(x["true_M_IR"], 2))
                      for v, x in out[f"E2_sig{sig}"].items()})
json.dump(out, open("data/t1_explore.json", "w"), indent=1)
