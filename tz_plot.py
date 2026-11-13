#!/usr/bin/env python3
"""
tz_plot.py -- look at the main halo's assembly history from tz_tree_catalog output.

    python tz_plot.py g205_z4_host_halo.csv --tag g205_z4

Writes g205_z4_mass.png, _metallicity.png, _sfr.png, _size.png and prints a
summary at the final snapshot. Add --show to open them interactively.

Every panel is an overlay. Quantities with different units go in separate
figures rather than onto a second y-axis, because a twin axis lets you imply
any crossing point you like by rescaling.

The stellar particle-type index is detected from the data, not assumed:
SubhaloStellarPhotometricsMassInRad is built from star particles, so whichever
SubhaloMassInRadType_k tracks it is the star index. This run has NTYPES=7.
"""

import argparse
import sys

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

# Validated categorical palette (CVD-safe in this order; do not reorder).
C = ["#4269d0", "#e8871a", "#a463f2", "#3ca951"]
INK, MUTED = "#1c1c1c", "#8a8a8a"


def style():
    matplotlib.rcParams.update({
        "figure.figsize": (7.2, 4.6), "figure.dpi": 140,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "axes.titlesize": 11,
        "axes.grid": True, "grid.color": "#e6e6e6", "grid.linewidth": 0.6,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "text.color": INK, "font.size": 10, "legend.frameon": False,
        "lines.linewidth": 2.0, "lines.markersize": 4,
    })


def star_index(d, ntypes=7):
    """Which PartType is stars, determined from the file."""
    ref = "Sub_SubhaloStellarPhotometricsMassInRad"
    cands = [k for k in range(ntypes) if "Sub_SubhaloMassInRadType_%d" % k in d]
    if ref not in d or not cands:
        return 4, "assumed (reference column absent)"
    r = d[ref].to_numpy(float)
    best, score = None, np.inf
    for k in cands:
        v = d["Sub_SubhaloMassInRadType_%d" % k].to_numpy(float)
        m = np.isfinite(r) & np.isfinite(v) & (r > 0)
        if m.sum() < 5:
            continue
        s = np.median(np.abs(v[m] - r[m]) / r[m])
        if s < score:
            best, score = k, s
    if best is None:
        return 4, "assumed (no usable overlap)"
    return best, "detected, median |diff| = %.2g" % score


def overlay(d, series, ylabel, title, path, logy=True, show=False):
    """One axes, N lines, legend plus a direct label at the last finite point."""
    series = [(lab, col) for lab, col in series if col in d]
    if not series:
        print("  skip %s: no columns present" % path)
        return
    fig, ax = plt.subplots()
    z = d["Redshift"].to_numpy(float)
    ends = []
    for i, (lab, col) in enumerate(series):
        y = d[col].to_numpy(float)
        m = np.isfinite(y) & (y > 0 if logy else np.isfinite(y))
        if not m.any():
            continue
        ax.plot(z[m], y[m], color=C[i % len(C)], label=lab, solid_capstyle="round")
        ends.append((z[m][-1], y[m][-1], lab, C[i % len(C)]))
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel("redshift")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", pad=10)
    ax.invert_xaxis()                      # time runs left to right
    ax.margins(x=0.04)
    if len(series) > 1:
        ax.legend(loc="upper left", fontsize=9)

    # Direct labels in a column outside the axes, pushed apart so converging
    # series (halo mass and total bound mass, typically) stay readable.
    fig.canvas.draw()
    to_axes = ax.transAxes.inverted()
    placed = []
    for x, y, lab, col in ends:
        f = float(to_axes.transform(ax.transData.transform((x, y)))[1])
        placed.append([f, lab, col])
    placed.sort(key=lambda r: r[0])
    gap = 0.058
    for j in range(1, len(placed)):
        if placed[j][0] - placed[j - 1][0] < gap:
            placed[j][0] = placed[j - 1][0] + gap
    over = placed[-1][0] - 1.0 if placed and placed[-1][0] > 1.0 else 0.0
    for f, lab, col in placed:
        ax.annotate(lab, xy=(1.015, f - over), xycoords="axes fraction",
                    color=col, fontsize=8.5, va="center", annotation_clip=False)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    print("  wrote %s" % path)
    if not show:
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--tag", default="host_halo")
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()

    d = pd.read_csv(args.csv)
    if "Redshift" not in d:
        sys.exit("no Redshift column; is this a tz_tree_catalog output?")
    style()

    k, how = star_index(d)
    print("stellar PartType index: %d  (%s)" % (k, how))
    S, G = "Sub_", "Group_"
    ms = lambda c: "phys_%s_Msun" % c
    kp = lambda c: "phys_%s_kpc" % c
    zs = lambda c: "phys_%s_Zsun" % c

    print("writing figures:")
    overlay(d, [
        ("halo $M_{200c}$",   ms(G + "Group_M_Crit200")),
        ("stellar, $<2r_{1/2}$", ms(S + "SubhaloMassInRadType_%d" % k)),
        ("gas, $<2r_{1/2}$",  ms(S + "SubhaloMassInRadType_0")),
        ("total bound",       ms(S + "SubhaloMass")),
    ], r"$M\ [\mathrm{M_\odot}]$", "Main halo mass assembly",
        "%s_mass.png" % args.tag, show=args.show)

    overlay(d, [
        ("gas",              zs(S + "SubhaloGasMetallicity")),
        ("gas, SFR-weighted", zs(S + "SubhaloGasMetallicitySfr")),
        ("stars",            zs(S + "SubhaloStarMetallicity")),
    ], r"$Z\ [Z_\odot]$", "Main galaxy metallicity",
        "%s_metallicity.png" % args.tag, show=args.show)

    overlay(d, [
        ("total",            S + "SubhaloSFR"),
        (r"within $2r_{1/2}$", S + "SubhaloSfrInRad"),
    ], r"SFR $[\mathrm{M_\odot\,yr^{-1}}]$", "Main galaxy star formation rate",
        "%s_sfr.png" % args.tag, show=args.show)

    overlay(d, [
        ("$R_{200c}$",          kp(G + "Group_R_Crit200")),
        ("stellar $r_{1/2}$",   kp(S + "SubhaloHalfmassRadType_%d" % k)),
        ("gas $r_{1/2}$",       kp(S + "SubhaloHalfmassRadType_0")),
    ], "physical kpc", "Main halo and galaxy sizes",
        "%s_size.png" % args.tag, show=args.show)

    # ---- summary at the final snapshot ----------------------------------
    last = d.dropna(subset=["Redshift"]).iloc[-1]
    print("\nfinal snapshot: snap %d, z = %.3f"
          % (int(last["snap"]), float(last["Redshift"])))
    for lab, col in [
        ("M_200c            [Msun]", ms(G + "Group_M_Crit200")),
        ("R_200c            [pkpc]", kp(G + "Group_R_Crit200")),
        ("M_star  (<2r1/2)  [Msun]", ms(S + "SubhaloMassInRadType_%d" % k)),
        ("M_gas   (<2r1/2)  [Msun]", ms(S + "SubhaloMassInRadType_0")),
        ("SFR               [Msun/yr]", S + "SubhaloSFR"),
        ("Z_gas             [Zsun]", zs(S + "SubhaloGasMetallicity")),
        ("Z_gas,SFR         [Zsun]", zs(S + "SubhaloGasMetallicitySfr")),
        ("Z_star            [Zsun]", zs(S + "SubhaloStarMetallicity")),
    ]:
        if col in d:
            print("  %-28s %.4g" % (lab, float(last[col])))

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()