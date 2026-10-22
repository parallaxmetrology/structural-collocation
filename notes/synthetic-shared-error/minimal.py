"""Minimal reproduction: a non-mass-conserving blur manufactures shared error between simulated observers.

Same random rain cells, two ways to 'blur' a coarse observer:
  widen    : re-render each Gaussian cell with width w + b and the same peak (as in the original T1 control)
  convolve : render at width w, then convolve with a normalised Gaussian of sd b (mass-conserving, as in T1-bis)
Run:  python minimal.py            (numbers only)
      python minimal.py --figures  (also writes fig1-3 PNGs; needs matplotlib)
"""
import sys
import numpy as np
from scipy.ndimage import gaussian_filter

N, SEED, C = 10000, 20261015, 64
YY, XX = np.mgrid[0:C, 0:C] - (C - 1) / 2
OBS = {"R": (0.5, 0.0, 0.3), "S1": (0.8, 1.0, 0.4), "S2": (1.5, 2.0, 0.6)}  # position sd, blur b, log-intensity sd
WIN = (slice(C // 2 - 7, C // 2 + 8),) * 2  # central 15 x 15 window


def render(cy, cx, w, I):
    f = np.zeros((C, C))
    for a, b, s, i in zip(cy, cx, w, I):
        f += i * np.exp(-((YY - a) ** 2 + (XX - b) ** 2) / (2 * s ** 2))
    return f


def observe(cells, b, mode):
    cy, cx, w, I = cells
    if mode == "widen":
        return render(cy, cx, w + b, I)
    f = render(cy, cx, w, I)
    return gaussian_filter(f, b, mode="constant", truncate=4.0) if b > 0 else f


def reads(f, region=(slice(None), slice(None))):
    A, Y, X = f[region], YY[region], XX[region]
    m = A.sum()
    cy, cx = (A * Y).sum() / m, (A * X).sum() / m
    return m, cy, cx, np.sqrt((A * ((Y - cy) ** 2 + (X - cx) ** 2)).sum() / m)


def simulate():
    rng = np.random.default_rng(SEED)
    out = []
    for _ in range(N):
        K = 1 + rng.poisson(1.5)
        truth = (rng.uniform(-5, 5, K), rng.uniform(-5, 5, K), rng.uniform(0.8, 2.5, K), rng.lognormal(1.5, 0.8, K))
        obs = {o: (truth[0] + rng.normal(0, sd, K), truth[1] + rng.normal(0, sd, K), truth[2], truth[3] * rng.lognormal(0, li, K))
               for o, (sd, _, li) in OBS.items()}
        out.append((truth, obs))
    return out


def table(W, mode, region):
    """rows = windows; dict of observer -> array (N, 4) of (mass, cy, cx, d)."""
    T = np.array([reads(render(*t), region) for t, _ in W])
    return {"T": T, **{o: np.array([reads(observe(ob[o], OBS[o][1], mode), region) for _, ob in W]) for o in OBS}}


def errcorr(D, a, b, j):
    res = [D[o][:, j] - np.polyval(np.polyfit(D["T"][:, j], D[o][:, j], 1), D["T"][:, j]) for o in (a, b)]
    return np.corrcoef(*res)[0, 1]


def errcorr_ci(D, a, b, j, n_boot=1000, seed=7):
    """95% percentile interval for errcorr, resampling whole synthetic fields (rows) with replacement."""
    rng = np.random.default_rng(seed); n = len(D["T"])
    bs = []
    for _ in range(n_boot):
        k = rng.integers(0, n, n)
        bs.append(errcorr({o: v[k] for o, v in D.items()}, a, b, j))
    return np.percentile(bs, [2.5, 97.5])


def tc(D, j, trip=("R", "S1", "S2")):
    Cm = np.cov(np.array([D[o][:, j] for o in trip]))
    est = {o: 1 - Cm[i, (i + 1) % 3] * Cm[i, (i + 2) % 3] / (Cm[i, i] * Cm[(i + 1) % 3, (i + 2) % 3]) for i, o in enumerate(trip)}
    true = {o: 1 - np.corrcoef(D[o][:, j], D["T"][:, j])[0, 1] ** 2 for o in trip}
    return est, true


if __name__ == "__main__":
    W = simulate()
    print("1. Analytic checks on the truth field (first 1000 windows), blur b = 1 and 2")
    for mode in ("widen", "convolve"):
        for b in (1.0, 2.0):
            r = np.array([(*reads(observe(t, b, mode)), *reads(render(*t))) for t, _ in W[:1000]])
            mass = np.median(r[:, 0] / r[:, 4])
            shift = np.median(np.hypot(r[:, 1] - r[:, 5], r[:, 2] - r[:, 6]))
            d2 = np.median(np.abs(r[:, 3] ** 2 - r[:, 7] ** 2 - 2 * b ** 2) / (2 * b ** 2))
            ok = abs(mass - 1) < 1e-6 and shift < 0.01 and d2 < 0.01
            print(f"   {mode:8s} b={b:.0f}: mass ratio {mass:6.3f}  centroid shift {shift:.1e}  "
                  f"|d2 growth - 2b^2|/2b^2 {d2:6.2%}  -> {'PASS' if ok else 'FAIL'}")
    tabs = {(m, g): table(W, m, reg) for m in ("widen", "convolve") for g, reg in (("full", (slice(None),) * 2), ("15x15", WIN))}
    print("\n2. S1-S2 coordinate-error correlation with nothing shared planted (should be ~0)")
    print("   blur      region   centroid y  centroid x  dispersion   (95% field-bootstrap intervals, 1000 resamples)")
    for (m, g), D in tabs.items():
        cells = []
        for j in (1, 2, 3):
            lo, hi = errcorr_ci(D, "S1", "S2", j)
            cells.append(f"{errcorr(D, 'S1', 'S2', j):+.3f} [{lo:+.3f}, {hi:+.3f}]")
        print(f"   {m:8s}  {g:6s}   " + "  ".join(cells))
    print("\n3. Triple-collocation noise fractions (estimate / truth), full canvas")
    for m in ("widen", "convolve"):
        for j, v in ((1, "centroid y"), (3, "dispersion")):
            est, true = tc(tabs[(m, "full")], j)
            order = lambda d: "<".join(sorted(d, key=d.get))
            print(f"   {m:8s} {v:10s} " + "  ".join(f"{o} {est[o]:.3f}/{true[o]:.3f}" for o in OBS)
                  + f"   ranking est {order(est)}, true {order(true)}" + ("" if order(est) == order(true) else "  MIS-RANKED"))
    if "--figures" in sys.argv:
        import figures
        figures.make(W, tabs, observe, render, errcorr, tc, OBS)
