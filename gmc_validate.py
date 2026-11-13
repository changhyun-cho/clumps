#!/usr/bin/env python3
"""
gmc_validate.py -- check the local GMC catalogs against the published numbers.

Every claim below is one that Wang et al. (2025, MNRAS 544, 2675) state in the
text, so agreement means the units and cuts in tz_gmc.py are right, and
disagreement points at a specific wrong assumption rather than a vague one.

    python gmc_validate.py            # numbers only
    python gmc_validate.py --figures  # also writes gmc_overview.png

The fiducial sample throughout is the 8x resolution level, resolved clouds only
(>= 30 elements), which is what the paper uses for cloud-scale analysis.
"""

import argparse

import numpy as np
import pandas as pd

import tz_gmc


def mle_slope(m, m_min):
    """MLE index of dN/dM ~ M^gamma above m_min, with its standard error."""
    x = np.asarray(m, float)
    x = x[np.isfinite(x) & (x >= m_min)]
    n = x.size
    if n < 30:
        return np.nan, np.nan, n
    s = np.sum(np.log(x / m_min))
    gamma = -1.0 - n / s
    return gamma, (gamma + 1.0) / np.sqrt(n), n


def fit_loglog(x, y):
    """Least-squares slope and normalisation of log10 y vs log10 x."""
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if m.sum() < 10:
        return np.nan, np.nan
    b, a = np.polyfit(np.log10(x[m]), np.log10(y[m]), 1)
    return b, a

def _ratio_at_fixed_sigma(a, b, lo=1.3, hi=2.6, nbin=6, nmin=10):
    """Median Sigma_SFR ratio b/a in matched Sigma_gas bins, so the comparison
    is not driven by the two samples living at different surface densities."""
    edges = np.logspace(lo, hi, nbin + 1)
    sa = a.sfr_yso / (np.pi * a.R_eff_pc ** 2)
    sb = b.sfr_yso / (np.pi * b.R_eff_pc ** 2)
    out = []
    for x, y in zip(edges[:-1], edges[1:]):
        ma = (a.Sigma_gas >= x) & (a.Sigma_gas < y)
        mb = (b.Sigma_gas >= x) & (b.Sigma_gas < y)
        if ma.sum() < nmin or mb.sum() < nmin:
            continue
        out.append(np.median(sb[mb]) / np.median(sa[ma]))
    if not out:
        return "-"
    return "%.2f  (%d bins)" % (np.median(out), len(out))


def report(label, got, expected, tol=None, unit="", ref=None):
    """expected is the paper's wording; ref/tol give a numeric tolerance check."""
    ok = ""
    if ref is not None and tol is not None and np.isfinite(got):
        ok = "  ok" if abs(got - ref) <= tol else "  <-- off"
    print("  %-46s %10.4g %-10s (paper: %s)%s"
          % (label, got, unit, expected, ok))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figures", action="store_true")
    args = ap.parse_args()

    pd.set_option("display.width", 200)
    d = tz_gmc.load(resolved_only=True)
    fid8 = d[(d.res == 8) & (d.variant == "fiducial")]
    no8 = d[(d.res == 8) & (d.variant == "noESF")]
    var8 = d[(d.res == 8) & (d.variant == "varSFE")]

    print("=" * 96)
    print("SAMPLE")
    print("=" * 96)
    print(d.groupby(["res", "variant"]).agg(
        N=("z", "size"), runs=("Run", "nunique"), haloes=("HaloID", "nunique"),
        z_hi=("z", "max"), z_lo=("z", "min")).to_string())

    print("\n" + "=" * 96)
    print("1. CLOUD SIZE AND MASS  (Sec. 4.2, Fig. 5)")
    print("=" * 96)
    report("median M_gas, 8x fiducial", np.median(fid8.M_gas), "~1e5", 5e4,
           "Msun", ref=1e5)
    report("median R_eff, 8x fiducial", np.median(fid8.R_eff_pc),
           "comparable to observed GMCs", None, "pc")
    report("median alpha, 8x fiducial", np.median(fid8.alpha), "~4", 1.5, ref=4)
    print("  fraction with alpha < 10, 8x fiducial: %.2f  (paper: 'majority')"
          % (fid8.alpha < 10).mean())

    # Mass function. The 30-element floor sets the completeness limit, so fit
    # a decade above it.
    print("\n  mass function  dN/dM ~ M^gamma   (paper: gamma = -2.49 +0.13 -0.07)")
    for lvl, sub in [(4, d[(d.res == 4) & (d.variant == "fiducial")]),
                     (8, fid8), (16, d[(d.res == 16)])]:
        m_min = 10 * np.median(sub.M_gas)
        g, err, n = mle_slope(sub.M_gas, m_min)
        print("    %2dx  M > %8.2e Msun   gamma = %6.2f +/- %.2f   (N = %d)"
              % (lvl, m_min, g, err, n))

    print("\n" + "=" * 96)
    print("2. SURFACE DENSITY  (Sec. 4.2, Fig. 8)")
    print("=" * 96)
    report("median Sigma_gas, 8x fiducial", np.median(fid8.Sigma_gas),
           "~70", 8, "Msun/pc^2", ref=70)
    report("median Sigma_gas, 8x noESF", np.median(no8.Sigma_gas),
           "~51", 8, "Msun/pc^2", ref=51)
    report("median Sigma_gas, 8x varSFE", np.median(var8.Sigma_gas),
           "~ fiducial", 12, "Msun/pc^2", ref=np.median(fid8.Sigma_gas))
    print("\n  redshift dependence of Sigma_gas, 8x fiducial "
          "(paper: no significant evolution)")
    zb = pd.cut(fid8.z, [3, 3.5, 4, 5, 6, 8, 10])
    t = fid8.groupby(zb, observed=True).agg(
        N=("Sigma_gas", "size"), Sigma=("Sigma_gas", "median"),
        R_pc=("R_eff_pc", "median"), Z=("Z_Zsun", "median"))
    print(t.to_string())

    print("\n" + "=" * 96)
    print("3. LARSON RELATION  sigma_v ~ R^b   (Sec. 4.2, Fig. 6)")
    print("=" * 96)
    print("  sigma_v is inverted from alpha under a uniform-sphere assumption,")
    print("  so the normalisation is approximate; the slope is the check.")
    for lab, sub in [("8x fiducial, z ~ 3-4", fid8[fid8.z.between(3, 4)]),
                     ("8x fiducial, z ~ 8-10", fid8[fid8.z.between(8, 10)]),
                     ("8x noESF,    z ~ 3-4", no8[no8.z.between(3, 4)])]:
        b, a = fit_loglog(sub.R_eff_pc.to_numpy(), sub.sigma_v.to_numpy())
        s20 = 10 ** (a + b * np.log10(20.0))
        print("    %-22s b = %5.2f   sigma_v(20 pc) = %4.2f km/s   (N = %d)"
              % (lab, b, s20, len(sub)))
    print("    %-22s b = %5.2f   sigma_v(20 pc) = %4.2f km/s"
          % ("Solomon+87 (MW)", 0.5, 0.72 * 20 ** 0.5))
    print("  paper: b ~ 0.6 across all systems, normalisation above local GMCs,")
    print("         and lower sigma_v at fixed R in noESF")

    print("\n" + "=" * 96)
    print("4. CLOUD-SCALE SFE  (Sec. 4.3, Figs. 9 and 10)")
    print("=" * 96)
    print("  SFR = M_young/3 Myr is exact from the catalogue, so t_dep = M_gas/SFR")
    print("  needs no assumption. eps_ff does: the paper's t_ff is the mass-weighted")
    print("  harmonic mean over cells, which is shorter than the mean-density t_ff")
    print("  used here, so our eps_ff is an upper bound.")
    for lab, sub in [("8x fiducial", fid8), ("8x noESF", no8),
                     ("8x varSFE", var8)]:
        sf = sub[sub.sfr_yso > 0]
        print("    %-14s star-forming fraction %.2f   median t_dep = %6.1f Myr   "
              "eps_ff <= %.3f   (mean-density t_ff = %.2f Myr)"
              % (lab, (sub.sfr_yso > 0).mean(),
                 np.median(sf.M_gas / sf.sfr_yso) / 1e6,
                 np.median(sf.eps_ff), np.median(sf.t_ff_Myr)))
    report("median t_dep, 8x fiducial",
           np.median(fid8[fid8.sfr_yso > 0].M_gas
                     / fid8[fid8.sfr_yso > 0].sfr_yso) / 1e6,
           "~100", 25, "Myr", ref=100)
    b, _ = fit_loglog(
        fid8[fid8.sfr_yso > 0].Sigma_gas.to_numpy(),
        (fid8[fid8.sfr_yso > 0].sfr_yso
         / (np.pi * fid8[fid8.sfr_yso > 0].R_eff_pc ** 2)).to_numpy())
    report("Sigma_SFR ~ Sigma_gas^b, 8x fiducial", b, "b = 1", 0.25, "", ref=1.0)

    # Inverting the published eps_ff gives the cell-level t_ff the catalogue
    # cannot supply, which is the cleanest way to state the systematic.
    sf = fid8[fid8.sfr_yso > 0]
    implied = 0.03 * 3.0 / np.median(sf.M_star / sf.M_gas)
    print("\n  Inverting the paper's median eps_ff = 0.03 with our exact SFR gives")
    print("  a cell-level <t_ff> of %.1f Myr, against %.1f Myr at the mean density."
          % (implied, np.median(sf.t_ff_Myr)))

    print("\n" + "=" * 96)
    print("4b. PHYSICS VARIANTS, HALO-MATCHED  (Sec. 4.2-4.3)")
    print("=" * 96)
    print("  The variant runs cover different haloes from the fiducial suite, so")
    print("  comparing whole catalogues confounds the model with the sample. Each")
    print("  row below is the same target halo at the same resolution.")
    pairs = [("g33206", 8, "g33206z8", "g33206noESF_z8", "noESF"),
             ("g10304", 4, "g10304z4", "g10304noESF_z4", "noESF"),
             ("g5760", 4, "g5760z4", "g5760noESF_z4", "noESF"),
             ("g33206", 8, "g33206z8", "g33206varEff_z8", "varSFE"),
             ("g5760", 4, "g5760z4", "g5760varEff_z4", "varSFE")]
    print("    %-8s %4s %-8s  %14s  %14s  %s"
          % ("halo", "res", "variant", "Sigma ratio", "t_dep ratio",
             "SFR ratio at fixed Sigma"))
    for halo, lvl, a, b_, var in pairs:
        A = d[d.Run == a]
        B = d[d.Run == b_]
        sA, sB = A[A.sfr_yso > 0], B[B.sfr_yso > 0]
        if len(sA) < 20 or len(sB) < 20:
            continue
        sig = np.median(B.Sigma_gas) / np.median(A.Sigma_gas)
        tdep = (np.median(sA.M_gas / sA.sfr_yso)
                / np.median(sB.M_gas / sB.sfr_yso))
        print("    %-8s %3dx %-8s  %14.2f  %14.2f  %s"
              % (halo, lvl, var, sig, tdep,
                 _ratio_at_fixed_sigma(sA, sB)))
    print("  paper: Sigma_gas 70 -> 51 without ESF (ratio 0.73), varSFE ~ fiducial;")
    print("         SFR at fixed Sigma is ~3x higher without ESF (Fig. 10, at 8x)")

    print("\n" + "=" * 96)
    print("5. METALLICITY  (not analysed in either paper; available here)")
    print("=" * 96)
    print("  Z is the sum of the seven tracked metals; H and He are normalised")
    print("  separately in these files so 1 - H - He is unusable.")
    print("  Z_sfr is NaN for clouds with no star-forming cell "
          "(%.0f per cent of 8x fiducial)."
          % (100 * fid8.Z_sfr.isna().mean()))
    print(fid8.groupby(pd.cut(fid8.z, [3, 4, 5, 6, 8, 10]), observed=True).agg(
        N=("Z_Zsun", "size"), Z_med=("Z_Zsun", "median"),
        Z16=("Z_Zsun", lambda x: x.quantile(0.16)),
        Z84=("Z_Zsun", lambda x: x.quantile(0.84))).to_string())

    if args.figures:
        make_figures(d, fid8, no8, var8)


def make_figures(d, fid8, no8, var8):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    C = ["#4269d0", "#e8871a", "#a463f2", "#3ca951"]
    INK, MUTED = "#1c1c1c", "#8a8a8a"
    matplotlib.rcParams.update({
        "figure.dpi": 140, "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "axes.titlesize": 10,
        "axes.grid": True, "grid.color": "#e9e9e9", "grid.linewidth": 0.6,
        "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
        "font.size": 9, "legend.frameon": False, "lines.linewidth": 1.8,
    })
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.4))

    # (a) mass function, one line per resolution level
    ax = axes[0, 0]
    for i, lvl in enumerate([4, 8, 16]):
        sub = d[(d.res == lvl) & (d.variant == "fiducial")]
        m = np.sort(sub.M_gas.to_numpy())
        ax.plot(m, np.arange(m.size, 0, -1) / m.size, color=C[i],
                label="%dx  (N=%d)" % (lvl, m.size))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$M_{\rm gas}\ [\mathrm{M_\odot}]$")
    ax.set_ylabel(r"$N(>M)/N$")
    ax.set_title("(a) cloud mass function by resolution", loc="left")
    ax.legend(fontsize=8)

    # (b) surface density vs redshift, fiducial vs variants
    ax = axes[0, 1]
    zb = [3, 3.5, 4, 4.5, 5, 6, 7, 8, 10]
    for i, (lab, sub) in enumerate([("fiducial", fid8), ("noESF", no8),
                                    ("varSFE", var8)]):
        g = sub.groupby(pd.cut(sub.z, zb), observed=True).Sigma_gas
        med, cnt = g.median(), g.size()
        keep = cnt >= 25          # otherwise single-snapshot bins dominate
        x = [iv.mid for iv in med.index[keep]]
        ax.plot(x, med[keep].to_numpy(), "o-", color=C[i], label=lab, ms=4)
    ax.axhline(70, color=MUTED, ls="--", lw=1)
    ax.annotate(r"$70\ \mathrm{M_\odot\,pc^{-2}}$", (9.3, 74), color=MUTED,
                fontsize=8, ha="right")
    ax.set_yscale("log"); ax.set_ylim(20, 200)
    ax.set_xlabel("redshift"); ax.invert_xaxis()
    ax.set_ylabel(r"$\Sigma_{\rm gas}\ [\mathrm{M_\odot\,pc^{-2}}]$")
    ax.set_title("(b) cloud surface density, 8x", loc="left")
    ax.legend(fontsize=8)

    # (c) Larson
    ax = axes[0, 2]
    sub = fid8.sample(min(4000, len(fid8)), random_state=0)
    ax.scatter(sub.R_eff_pc, sub.sigma_v, s=3, alpha=0.12, color=C[0],
               edgecolors="none", rasterized=True)
    r = np.logspace(np.log10(6), np.log10(120), 40)
    b, a = fit_loglog(fid8.R_eff_pc.to_numpy(), fid8.sigma_v.to_numpy())
    ax.plot(r, 10 ** (a + b * np.log10(r)), color=C[0],
            label=r"8x fid.  $b=%.2f$" % b)
    b2, a2 = fit_loglog(no8.R_eff_pc.to_numpy(), no8.sigma_v.to_numpy())
    ax.plot(r, 10 ** (a2 + b2 * np.log10(r)), color=C[1],
            label=r"noESF  $b=%.2f$" % b2)
    ax.plot(r, 0.72 * r ** 0.5, color=INK, ls="--", lw=1.2,
            label=r"Solomon+87  $b=0.5$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$R_{\rm eff}\ [\mathrm{pc}]$")
    ax.set_ylabel(r"$\sigma_v\ [\mathrm{km\,s^{-1}}]$")
    ax.set_title(r"(c) Larson relation ($\sigma_v$ from $\alpha$)", loc="left")
    ax.legend(fontsize=8, loc="upper left")

    # (d) eps_ff distribution
    ax = axes[1, 0]
    bins = np.logspace(-3.5, 0.3, 45)
    for i, (lab, sub) in enumerate([("fiducial", fid8), ("noESF", no8),
                                    ("varSFE", var8)]):
        v = sub.eps_ff[sub.sfr_yso > 0]
        ax.hist(v, bins=bins, histtype="step", density=True, color=C[i],
                label="%s  (med %.3f)" % (lab, np.median(v)))
    ax.set_xscale("log")
    ax.set_xlabel(r"$\epsilon_{\rm ff}^{\rm GMC}$ (mean-density $t_{\rm ff}$)")
    ax.set_ylabel("PDF")
    ax.set_title("(d) cloud-scale SFE, star-forming clouds", loc="left")
    ax.legend(fontsize=8)

    # (e) cloud KS relation, halo-matched so the noESF offset is not a sample effect
    ax = axes[1, 1]
    fidm = d[(d.Run == "g33206z8") & (d.sfr_yso > 0)]
    nom = d[(d.Run == "g33206noESF_z8") & (d.sfr_yso > 0)]
    edges = np.logspace(1.3, 2.6, 7)
    mids = np.sqrt(edges[:-1] * edges[1:])
    sub = fid8[fid8.sfr_yso > 0].sample(3000, random_state=0)
    ax.scatter(sub.Sigma_gas, sub.sfr_yso / (np.pi * sub.R_eff_pc ** 2),
               s=3, alpha=0.08, color=MUTED, edgecolors="none", rasterized=True)
    for i, (lab, sb) in enumerate([("g33206 fiducial", fidm),
                                   ("g33206 noESF", nom)]):
        ss = sb.sfr_yso / (np.pi * sb.R_eff_pc ** 2)
        med = [np.median(ss[(sb.Sigma_gas >= lo) & (sb.Sigma_gas < hi)])
               if ((sb.Sigma_gas >= lo) & (sb.Sigma_gas < hi)).sum() >= 10
               else np.nan for lo, hi in zip(edges[:-1], edges[1:])]
        ax.plot(mids, med, "o-", color=C[i], ms=4, label=lab)
    ax.plot(edges, edges / 1e8, color=INK, ls="--", lw=1.2,
            label=r"$t_{\rm dep}=100$ Myr")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(15, 500); ax.set_ylim(3e-8, 3e-5)
    ax.set_xlabel(r"$\Sigma_{\rm gas}\ [\mathrm{M_\odot\,pc^{-2}}]$")
    ax.set_ylabel(r"$\Sigma_{\rm SFR}\ [\mathrm{M_\odot\,yr^{-1}\,pc^{-2}}]$")
    ax.set_title("(e) cloud KS relation, same halo with and without ESF",
                 loc="left")
    ax.legend(fontsize=8, loc="upper left")

    # (f) metallicity vs redshift
    ax = axes[1, 2]
    zb = [3, 3.5, 4, 4.5, 5, 6, 7, 8, 10]
    g = fid8.groupby(pd.cut(fid8.z, zb), observed=True).Z_Zsun
    x = [iv.mid for iv in g.median().index]
    ax.plot(x, g.median().to_numpy(), "o-", color=C[0], ms=4, label="GMCs, 8x")
    ax.fill_between(x, g.quantile(0.16).to_numpy(), g.quantile(0.84).to_numpy(),
                    color=C[0], alpha=0.15, lw=0)
    host = pd.read_csv("g205_z4_host_halo.csv")
    hm = host.phys_Sub_SubhaloGasMetallicity_Zsun > 0
    ax.plot(host.Redshift[hm], host.phys_Sub_SubhaloGasMetallicity_Zsun[hm],
            color=C[1], lw=1.4, label="g205 host galaxy, gas")
    ax.set_yscale("log"); ax.set_xlim(10.5, 2.8); ax.set_ylim(3e-3, 1)
    ax.set_xlabel("redshift")
    ax.set_ylabel(r"$Z\ [Z_\odot]$")
    ax.set_title("(f) cloud metallicity vs host", loc="left")
    ax.legend(fontsize=8, loc="upper left")

    fig.tight_layout()
    fig.savefig("gmc_overview.png", bbox_inches="tight")
    print("\nwrote gmc_overview.png")


if __name__ == "__main__":
    main()
