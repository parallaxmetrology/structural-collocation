"""Score T10 (T10_PLAN.md) for the period whose results sit in ./data (run from the root for July, runs/aug for August)."""
import json, sys

SUF = sys.argv[1] if len(sys.argv) > 1 else ""
D = lambda f: json.load(open(f"data/{f}.json"))
A, T7, TR, TS = D("t10_area" + SUF), D("t7_results"), D("t7r_results"), D("t7s_results")
AX = ["cy", "cx"]
above = lambda ci: ci[0] > 0
below = lambda ci: ci[1] < 0
rows, primary = [], {}


def rec(h, label, val, ci, ok):
    rows.append((h, label, val, ci, ok)); primary.setdefault(h, []).append(ok)


for v in AX:
    rec("H1", f"area Δ {v}", A[v]["point"]["diff_R_MW"], A[v]["ci95"]["diff_R_MW"], above(A[v]["ci95"]["diff_R_MW"]))
for v in AX:
    w = T7["weak"][v]
    rec("H2", f"point (R,GIR,MW) Δ weak {v}", w["point"]["RGM_diff_R_MW"], w["ci95"]["RGM_diff_R_MW"], above(w["ci95"]["RGM_diff_R_MW"]))
for s in ["weak", "strong"]:
    for v in AX:
        w = T7[s][v]
        rec("H3", f"point (R,MW,G) Δ {s} {v}", w["point"]["RMG_diff_R_MW"], w["ci95"]["RMG_diff_R_MW"], below(w["ci95"]["RMG_diff_R_MW"]))
for v in AX:
    w = TR["weak"][v]
    for k in ["R1_R_minus_IR", "R1_R_minus_G"]:
        rec("H4", f"mirror {k} weak {v}", w["point"][k], w["ci95"][k], below(w["ci95"][k]))
for v in AX:
    w = T7["weak"][v]
    rec("H5", f"c(GIR,MW) weak {v}", w["point"]["c_IR_MW"], w["ci95"]["c_IR_MW"], above(w["ci95"]["c_IR_MW"]))
w = TR["weak"]["cy"]
for k, test in [("t1", above), ("t2", above), ("t3", lambda ci: ci[0] <= 0 <= ci[1])]:
    rec("H6", f"tetrad {k} weak cy", w["point"][k], w["ci95"][k], test(w["ci95"][k]))

# inconclusive gates (scored before predictions)
n_pts, n_area = T7["n_total"], A["n"]
neg = {"area " + v: A[v]["share_boot_negative_ev"] for v in AX}
neg.update({f"T7 {s} {v} {k}": T7[s][v]["share_negative"][k] for s in ["weak", "strong"] for v in AX for k in ["RGM_neg", "RMG_neg"]})
inconclusive = n_pts < 1000 or n_area < 5000 or max(neg.values()) > 0.20

met = {h: all(x) for h, x in primary.items()}
wrong_sign = any(r[2] <= 0 for r in rows if r[0] == "H1") or any(r[2] >= 0 for r in rows if r[0] == "H3")
verdict = "Inconclusive" if inconclusive else "Replicated" if all(met.values()) else "Not replicated" if wrong_sign else "Partial"

out = {"n_points": n_pts, "n_area": n_area, "max_negative_share": max(neg.values()), "negative_shares": neg,
       "rows": [{"h": h, "quantity": l, "point": p, "ci95": c, "met": o} for h, l, p, c, o in rows],
       "met": met, "verdict": verdict,
       "secondary_T7S": {k: {kk: TS[k][kk] for kk in ("rho_star_1", "rho_star_2", "rho_max")} for k in TS if "_" in k and isinstance(TS[k], dict) and "rho_star_1" in TS[k]}}
json.dump(out, open(f"data/t10_score{SUF}.json", "w"), indent=1)
print(f"n points {n_pts}, n area {n_area}, max negative-variance share {max(neg.values()):.3f}")
for h, l, p, c, o in rows:
    print(f"{h} {l:34s} {p:+.3f} [{c[0]:+.3f}, {c[1]:+.3f}] {'met' if o else 'NOT met'}")
print("T7-S", out["secondary_T7S"])
print("VERDICT", verdict)
