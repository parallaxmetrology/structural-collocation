"""Figures for NOTE.html. Called by `python minimal.py --figures`; reuses minimal.py's functions and results."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WIDEN, CONV, INK, MUTED, GRID = "#eb6834", "#2a78d6", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})


def make(W, tabs, observe, render, errcorr, tc, OBS):
    # Fig 1: the mechanism. (a) profiles of a narrow and a wide cell; (b) mass ratio vs cell width.
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7))
    x = np.linspace(-8, 8, 400); b = 2.0
    for w, peak, ls in ((0.8, 1.0, "-"), (2.5, 1.0, "--")):
        ax[0].plot(x, peak * np.exp(-x**2 / (2 * w**2)), color=INK, lw=1.2, ls=ls)
        ax[0].plot(x, peak * np.exp(-x**2 / (2 * (w + b)**2)), color=WIDEN, lw=2, ls=ls)
        s = np.sqrt(w**2 + b**2); ax[0].plot(x, peak * (w / s)**2 * np.exp(-x**2 / (2 * s**2)), color=CONV, lw=2, ls=ls)
    ax[0].text(0.3, 1.02, "widened (peak fixed)", color=WIDEN); ax[0].text(4.6, 0.1, "convolved", color=CONV)
    ax[0].text(-7.8, 0.85, "truth\nsolid w = 0.8\ndashed w = 2.5", color=INK, fontsize=8)
    ax[0].set(xlabel="distance from cell centre (cells)", ylabel="intensity (peak = 1)", title="(a) one cell, b = 2")
    w = np.linspace(0.8, 2.5, 100)
    for bb, ls in ((1.0, "--"), (2.0, "-")):
        ax[1].plot(w, ((w + bb) / w)**2, color=WIDEN, lw=2, ls=ls); ax[1].text(2.52, ((2.5 + bb) / 2.5)**2, f"widen, b = {bb:.0f}", color=WIDEN, va="center")
    ax[1].plot(w, np.ones_like(w), color=CONV, lw=2); ax[1].text(2.52, 1, "convolve, any b", color=CONV, va="center")
    ax[1].set(xlabel="true cell width w (cells)", ylabel="mass after / before", title="(b) mass multiplier ((w+b)/w)²", xlim=(0.8, 3.6))
    for a in ax: a.grid(color=GRID, lw=0.6)
    fig.tight_layout(); fig.savefig("fig1_mechanism.png"); plt.close(fig)

    # Fig 2: S1-S2 error correlation, nothing shared planted.
    labels = ["centroid y", "centroid x", "dispersion"]
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.6), sharey=True)
    for k, g in enumerate(("full", "15x15")):
        for i, (m, c) in enumerate((("widen", WIDEN), ("convolve", CONV))):
            v = [errcorr(tabs[(m, g)], "S1", "S2", j) for j in (1, 2, 3)]
            xs = np.arange(3) + (i - 0.5) * 0.36
            ax[k].bar(xs, v, 0.34, color=c, label=m)
            for xx, vv in zip(xs, v): ax[k].text(xx, max(vv, 0) + 0.01, f"{vv:+.2f}", ha="center", fontsize=7, color=INK)
        ax[k].axhline(0, color=MUTED, lw=0.8); ax[k].set_xticks(range(3), labels)
        ax[k].set_title("full 64 × 64 canvas" if g == "full" else "central 15 × 15 window"); ax[k].grid(axis="y", color=GRID, lw=0.6)
    ax[0].set_ylabel("corr(S1 error, S2 error)"); ax[0].legend(frameon=False, loc="upper left")
    fig.tight_layout(); fig.savefig("fig2_shared_error.png"); plt.close(fig)

    # Fig 3: TC estimate vs truth, full canvas.
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
    for k, (j, v) in enumerate(((1, "centroid y"), (3, "dispersion"))):
        ax[k].plot([0, 0.35], [0, 0.35], color=MUTED, lw=0.8, ls=":")
        for m, c in (("widen", WIDEN), ("convolve", CONV)):
            est, true = tc(tabs[(m, "full")], j)
            ax[k].scatter([true[o] for o in OBS], [est[o] for o in OBS], s=40, color=c, edgecolor="white", lw=1.5, label=m, zorder=3)
            for o in OBS: ax[k].annotate(o, (true[o], est[o]), xytext=(6, -8) if m == "widen" else (-12, 4), textcoords="offset points", fontsize=7, color=INK)
        ax[k].set(title=v, xlabel="true noise fraction", xlim=(0, 0.35), ylim=(0, 0.35)); ax[k].set_xticks(np.arange(0, 0.31, 0.1)); ax[k].grid(color=GRID, lw=0.6)
    ax[0].set_ylabel("triple-collocation estimate"); ax[0].legend(frameon=False, loc="upper left")
    fig.tight_layout(); fig.savefig("fig3_collocation.png"); plt.close(fig)
    print("\nwrote fig1_mechanism.png, fig2_shared_error.png, fig3_collocation.png")
