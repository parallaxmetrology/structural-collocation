"""T10 H1: the Probe A2 area analysis, with nf_R - nf_MW added to each bootstrap draw.
probe_a.analyse is run unchanged (same seed 11, 500 resamples, same draws); the quad() hook it already
calls is wrapped so the difference is computed from the same covariance matrix. Logged in T10_RESULTS.md."""
import json, sys
import numpy as np, pandas as pd
import probe_a as P

_quad = P.quad


def quad(C):
    out = _quad(C)
    nf_R = 1 - C[0, 1] * C[0, 2] / (C[0, 0] * C[1, 2])
    nf_MW = 1 - C[2, 0] * C[2, 1] / (C[2, 2] * C[0, 1])
    out["diff_R_MW"] = nf_R - nf_MW
    return out


P.quad = quad

if __name__ == "__main__":
    d = pd.read_csv("data/probe_a2_windows.csv", parse_dates=["t"])
    res = P.analyse(d)
    json.dump(res, open("data/t10_area" + (sys.argv[1] if len(sys.argv) > 1 else "") + ".json", "w"), indent=1)
    for v in ["cy", "cx"]:
        print(v, "diff", round(res[v]["point"]["diff_R_MW"], 4), np.round(res[v]["ci95"]["diff_R_MW"], 4))
