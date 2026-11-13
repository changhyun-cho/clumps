#!/usr/bin/env python3
"""
clump_host_analysis.py -- independent version of the g205z4 host-vs-clump analysis.

    python clump_host_analysis.py

Writes four figures mirroring the student's set, so the two can be compared
panel by panel:

    indep_host_core.png       host metallicity, masses and SFR, z = 16 -> 3
    indep_option1_dZ.png      clump vs host gas metallicity, and their difference
    indep_option2_radial.png  clump metallicity vs distance from the host centre
    indep_option3_sfr.png     host SFR vs summed clump SFR and clump counts

Definitions, fixed before looking at results:

  clumps      CloudPhinder clouds of run g205z4 (4x) whose HaloID equals the
              main-branch Group_GroupID from the tree in that snapshot, and with
              Number_of_Particles >= 30 (readme resolution limit).
  Z_clump     sum of the C, N, O, Ne, Mg, Si, Fe mass fractions / 0.0127.
  Z_host      the same seven-species sum from Sub_SubhaloGasMetalFractions
              (all gas bound to the central subhalo), / 0.0127, so that the
              difference compares like with like.
  distance    |CenterPosition|, which is already relative to GroupPos, in
              physical kpc; also shown as d / R200c from the tree.
  host SFR    Sub_SubhaloSFR, instantaneous, Msun/yr.
  clump SFR   sum of SFR_inst (instantaneous, Msun/yr) over the clumps in a
              snapshot; the 3 Myr SFR_yso/1000 is shown as a second estimate.
  host M*     SubhaloMassInRadType_4 (stars within 2 r_1/2,*); M_gas likewise
              with type 0.
"""

import numpy as np
import pandas as pd

import tz_gmc

ZSUN = 0.0127
C = ["#4269d0", "#e8871a", "#a463f2", "#3ca951"]
INK, MUTED = "#1c1c1c", "#8a8a8a"


def load():
    h = pd.read_csv(tz_gmc.host_tree_path("g205_z4"))
    h["snap"] = h["snap"].astype(int)
    h["Zhost"] = sum(h["Sub_SubhaloGasMetalFractions_%d" % i]
                     for i in range(2, 9)) / ZSUN
    h["Zhost_gfm"] = h["phys_Sub_SubhaloGasMetallicity_Zsun"]
    h["Mstar"] = h["phys_Sub_SubhaloMassInRadType_4_Msun"]
    h["Mgas"] = h["phys_Sub_SubhaloMassInRadType_0_Msun"]
    h["SFR"] = h["Sub_SubhaloSFR"]
    h["R200"] = h["phys_Group_Group_R_Crit200_kpc"]

    d = tz_gmc.load(res=4, variant="fiducial")
    g = d[d.Run == "g205z4"].copy()
    main_gid = dict(zip(h.snap, h.Group_GroupID.astype(int)))
    g = g[(g.HaloID == g.Snapshot.map(main_gid))
          & (g.Number_of_Particles >= 30)].copy()
    for col in ["Zhost", "Zhost_gfm", "R200", "SFR"]:
        g[col] = g.Snapshot.map(dict(zip(h.snap, h[col])))
    g["dZ"] = g.Z_Zsun - g.Zhost
    g["x"] = g.d_center_kpc / g.R200
    return h, g


def style():
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({
        "figure.dpi": 140, "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "axes.titlesize": 10,
        "axes.grid": True, "grid.color": "#e9e9e9", "grid.linewidth": 0.6,
        "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
        "font.size": 9, "legend.frameon": False, "lines.linewidth": 1.8})


def fig_host(h):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(3, 1, figsize=(7.2, 8.6), sharex=True)
    ax[0].plot(h.Redshift, h.Zhost, color=C[0], label="species sum (C..Fe)")
    ax[0].plot(h.Redshift, h.Zhost_gfm, color=C[1], lw=1.2,
               label="SubhaloGasMetallicity (GFM)")
    ax[0].set_ylabel(r"host gas $Z\ [Z_\odot]$"); ax[0].set_yscale("log")
    ax[0].legend(fontsize=8, loc="upper left")
    ax[0].set_title("Host galaxy g205 (main branch)", loc="left")
    m = h.Mstar > 0
    ax[1].plot(h.Redshift[m], h.Mstar[m], color=C[0],
               label=r"$M_*$, $<2r_{1/2,*}$")
    ax[1].plot(h.Redshift, h.Mgas, color=C[1],
               label=r"$M_{\rm gas}$, $<2r_{1/2,*}$")
    ax[1].set_yscale("log"); ax[1].set_ylabel(r"$M\ [\mathrm{M_\odot}]$")
    ax[1].legend(fontsize=8, loc="upper left")
    m = h.SFR > 0
    ax[2].plot(h.Redshift[m], h.SFR[m], color=C[3])
    ax[2].set_yscale("log")
    ax[2].set_ylabel(r"SFR $[\mathrm{M_\odot\,yr^{-1}}]$")
    ax[2].set_xlabel("redshift"); ax[2].invert_xaxis()
    fig.tight_layout(); fig.savefig("indep_host_core.png", bbox_inches="tight")
    plt.close(fig)


def fig_option1(h, g):
    import matplotlib.pyplot as plt
    per = g.groupby("Snapshot").agg(z=("z", "first"), Zc=("Z_Zsun", "median"),
                                   dZ=("dZ", "median"), Zh=("Zhost", "first"))
    fig, ax = plt.subplots(2, 1, figsize=(7.2, 7.4), sharex=True)
    a = ax[0]
    a.scatter(g.z, g.Z_Zsun, s=3, alpha=.08, color=C[0], edgecolors="none",
              rasterized=True)
    a.plot(per.z, per.Zc, "o", ms=3.5, color=C[0], label="clumps, per-snapshot median")
    hm = h.Redshift <= g.z.max() + 0.3
    a.plot(h.Redshift[hm], h.Zhost[hm], color=C[1], label="host gas")
    a.set_yscale("log"); a.set_ylabel(r"$Z\ [Z_\odot]$")
    a.set_title("Option 1, independent: clump and host gas metallicity", loc="left")
    a.legend(fontsize=8, loc="upper left")

    a = ax[1]
    a.scatter(g.z, g.dZ, s=3, alpha=.08, color=C[2], edgecolors="none",
              rasterized=True, label="individual clumps")
    a.plot(per.z, per.dZ, "o", ms=3.5, color=C[2], label="per-snapshot median")
    a.axhline(0, color=INK, lw=1, ls="--")
    a.set_ylim(-1.0, 0.4)
    a.set_ylabel(r"$\Delta Z = Z_{\rm clump}-Z_{\rm host}\ [Z_\odot]$")
    a.set_xlabel("redshift"); a.invert_xaxis()
    leg = a.legend(fontsize=8, loc="lower left")
    for lh in leg.legend_handles:
        lh.set_alpha(1)
    fig.tight_layout(); fig.savefig("indep_option1_dZ.png", bbox_inches="tight")
    plt.close(fig)
    return per


def fig_option2(g):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    a = ax[0]
    sc = a.scatter(g.d_center_kpc, g.Z_Zsun, c=g.z, s=4, alpha=.35,
                   cmap="viridis_r", edgecolors="none", rasterized=True)
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel(r"distance from host centre, $|\mathrm{CenterPosition}|$ [kpc]")
    a.set_ylabel(r"$Z_{\rm clump}\ [Z_\odot]$")
    a.set_title("Option 2, independent: all clumps, coloured by redshift",
                loc="left")
    cb = fig.colorbar(sc, ax=a); cb.set_label("redshift")

    a = ax[1]
    edges = [0.005, 0.02, 0.05, 0.1, 0.2, 1.0]
    for i, (lo, hi) in enumerate([(3, 3.5), (3.5, 4), (4, 5), (5, 6)]):
        s = g[g.z.between(lo, hi)]
        gb = s.groupby(pd.cut(s.x, edges), observed=True).Z_Zsun
        med, n = gb.median(), gb.size()
        k = n >= 20
        if k.sum() < 2:
            continue
        a.plot([iv.mid for iv in med.index[k]], med[k].to_numpy(), "o-",
               color=C[i], ms=4, label=r"$z=%.1f$-$%.1f$" % (lo, hi))
    a.set_xscale("log")
    a.set_xlabel(r"$d\,/\,R_{\rm 200c}$")
    a.set_ylabel(r"median $Z_{\rm clump}\ [Z_\odot]$")
    a.set_title("at fixed redshift (bins with >= 20 clumps)", loc="left")
    a.legend(fontsize=8)
    fig.tight_layout(); fig.savefig("indep_option2_radial.png", bbox_inches="tight")
    plt.close(fig)


def fig_option3(h, g):
    import matplotlib.pyplot as plt
    agg = g.groupby("Snapshot").agg(
        sfr_inst=("sfr_inst", "sum"), sfr_yso=("sfr_yso", "sum"),
        n_all=("z", "size"), n_sf=("sfr_inst", lambda s: (s > 0).sum()))
    full = h.set_index("snap")[["Redshift", "SFR"]].join(agg)
    full[["n_all", "n_sf"]] = full[["n_all", "n_sf"]].fillna(0)

    fig, ax = plt.subplots(2, 1, figsize=(7.6, 7.4), sharex=True)
    a = ax[0]
    sel = full.Redshift <= 9.8
    f = full[sel]
    a.plot(f.Redshift, f.SFR.where(f.SFR > 0), color=C[3], lw=1.6,
           label="host SFR (all snapshots)")
    a.plot(f.Redshift, f.sfr_inst.where(f.sfr_inst > 0), "o", ms=3,
           color=C[0], label=r"$\Sigma$ clump SFR$_{\rm inst}$")
    a.plot(f.Redshift, f.sfr_yso.where(f.sfr_yso > 0), "^", ms=3,
           color=C[1], label=r"$\Sigma$ clump SFR$_{\rm yso}$ (3 Myr)")
    a.set_yscale("log"); a.set_ylabel(r"SFR $[\mathrm{M_\odot\,yr^{-1}}]$")
    a.set_title("Option 3, independent: host SFR vs clump activity", loc="left")
    a.legend(fontsize=8, loc="upper left")
    a = ax[1]
    a.plot(f.Redshift, f.n_all, color=C[0], lw=1.4, label="resolved clumps")
    a.plot(f.Redshift, f.n_sf, color=C[2], lw=1.4,
           label=r"star-forming (SFR$_{\rm inst}>0$)")
    a.set_ylabel("clumps per snapshot"); a.set_xlabel("redshift")
    a.invert_xaxis(); a.legend(fontsize=8, loc="upper left")
    fig.tight_layout(); fig.savefig("indep_option3_sfr.png", bbox_inches="tight")
    plt.close(fig)
    return full


def lag_correlation(full, zmax=6.0, maxlag=4):
    """Spearman correlation of log host SFR with log clump SFR shifted by k
    snapshots (~10 Myr each), over snapshots where both are positive."""
    from scipy.stats import spearmanr
    f = full[full.Redshift <= zmax].sort_index()
    out = []
    for k in range(-maxlag, maxlag + 1):
        cl = f.sfr_inst.shift(-k)
        m = (f.SFR > 0) & (cl > 0)
        if m.sum() < 10:
            continue
        r, _ = spearmanr(np.log10(f.SFR[m]), np.log10(cl[m]))
        out.append((k, r, int(m.sum())))
    return out


def main():
    pd.set_option("display.width", 200)
    h, g = load()
    style()
    print("clumps: %d resolved, main branch, %d snapshots, z = %.2f - %.2f"
          % (len(g), g.Snapshot.nunique(), g.z.min(), g.z.max()))

    fig_host(h)
    per = fig_option1(h, g)
    fig_option2(g)
    full = fig_option3(h, g)

    print("\nOPTION 1  clump vs host gas metallicity")
    t = per.groupby(pd.cut(per.z, [3, 3.5, 4, 5, 6, 8, 10]), observed=True).agg(
        snaps=("z", "size"), Z_clump=("Zc", "median"), Z_host=("Zh", "median"),
        dZ=("dZ", "median"))
    print(t.round(3).to_string())
    lo = g[g.z < 4]
    print("  individual clumps at z<4: dZ 16/50/84th pct = %s"
          % np.round(np.percentile(lo.dZ, [16, 50, 84]), 3))

    print("\nOPTION 2  median Z_clump vs d/R200c at fixed redshift")
    for zlo, zhi in [(3, 3.5), (4, 5)]:
        s = g[g.z.between(zlo, zhi)]
        tt = s.groupby(pd.cut(s.x, [0, 0.02, 0.05, 0.1, 0.2, 1.0]),
                       observed=True).agg(N=("x", "size"), Z=("Z_Zsun", "median"),
                                          d_kpc=("d_center_kpc", "median"))
        print("  z = %.1f-%.1f" % (zlo, zhi))
        print("   " + tt[tt.N >= 20].round(3).to_string().replace("\n", "\n   "))
    print("  distance range: %.2f - %.1f kpc (median %.2f)"
          % (g.d_center_kpc.min(), g.d_center_kpc.max(), g.d_center_kpc.median()))

    print("\nOPTION 3  host SFR vs summed clump SFR")
    f = full.dropna(subset=["sfr_inst"])
    ratio = (f.sfr_inst / f.SFR).replace([np.inf, -np.inf], np.nan).dropna()
    print("  snapshots with clumps: %d of %d" % (len(f), len(full)))
    print("  sum(clump SFR_inst) / host SFR: median %.2f (16-84%%: %.2f - %.2f)"
          % (ratio.median(), ratio.quantile(.16), ratio.quantile(.84)))
    print("  lagged rank correlation, z<6 (k>0: clumps lag the host by k snapshots)")
    for k, r, n in lag_correlation(full):
        print("    k = %+d   rho = %.2f   (N = %d)" % (k, r, n))

    print("\nwrote indep_host_core.png, indep_option1_dZ.png, "
          "indep_option2_radial.png, indep_option3_sfr.png")


if __name__ == "__main__":
    main()
