#!/usr/bin/env python3
"""
audit_checks.py -- test every derived quantity against something independent
of the assumption that produced it.

    python audit_checks.py

Each check prints PASS / WARN / FAIL with the numbers behind it. Nothing here
modifies data; it reads the catalogs, the g205 tree and inspect.output.
"""

import re

import h5py
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import tz_gmc

H0, OM, OL, HUB = 67.74, 0.3089, 0.6911, 0.6774      # Planck 2015, as in the run
G_KPC = 4.30091e-6                                     # kpc (km/s)^2 / Msun
METALS = ["C", "N", "O", "Ne", "Mg", "Si", "Fe"]
EL = ["H", "He"] + METALS


def verdict(ok, warn=False):
    return "PASS" if ok else ("WARN" if warn else "FAIL")


def section(t):
    print("\n" + "=" * 92 + "\n" + t + "\n" + "=" * 92)


tree = pd.read_csv("g205_z4_host_halo.csv")
tree["snap"] = tree["snap"].astype(int)
tree = tree.set_index("snap")

# ------------------------------------------------------------------ 1
section("1. Snapshot -> redshift map")
txt = open("inspect.output").read()
hdr = {}
for m in re.finditer(r"groups_(\d{3})/fof_subhalo_tab_\d{3}\.0\.hdf5.*?Redshift\s+=\s+([\d.]+)", txt, re.S):
    hdr[int(m.group(1))] = float(m.group(2))
for s, z in sorted(hdr.items()):
    dz = abs(tree.loc[s, "Redshift"] - z)
    print("  snap %3d  catalog header z = %.6f   tree z = %.6f   %s"
          % (s, z, tree.loc[s, "Redshift"], verdict(dz < 1e-6)))
zz = tree["Redshift"].to_numpy()
print("  monotonic decreasing: %s;  range %.3f -> %.3f over %d snapshots"
      % (verdict(np.all(np.diff(zz) < 0)), zz[0], zz[-1], len(zz)))
print("  NOTE: the same grid is assumed for every run; only g205's headers are available to check.")

# ------------------------------------------------------------------ 2
section("2. Tree unit conversions: does M200c match R200c for the run's cosmology?")
a = 1 / (1 + tree["Redshift"])
Hz = H0 * np.sqrt(OM * (1 + tree["Redshift"]) ** 3 + OL)                  # km/s/Mpc
rho_c = 3 * (Hz / 1e3) ** 2 / (8 * np.pi * G_KPC)                          # Msun/kpc^3
pred = 200 * rho_c * 4 / 3 * np.pi * tree["phys_Group_Group_R_Crit200_kpc"] ** 3
ratio = (tree["phys_Group_Group_M_Crit200_Msun"] / pred)
ratio = ratio[tree["phys_Group_Group_M_Crit200_Msun"] > 0]
print("  M200c / (200 rho_crit(z) 4/3 pi R200c^3): median %.4f, range %.4f - %.4f   %s"
      % (ratio.median(), ratio.min(), ratio.max(), verdict(abs(ratio - 1).max() < 0.01)))
print("  (checks mass x1e10/h, length x a/h and the redshift together)")

# ------------------------------------------------------------------ 3
section("3. Metallicity conventions")
with h5py.File("gmc_catalog/gmc_output.hdf5") as f:
    g = f["GMC_alpha20_n100"]
    X = {e: g[f"{e}_Fraction"][...] for e in EL}
    npart = g["Number_of_Particles"][...]
r = npart >= 30
tot = sum(X[e][r] for e in EL)
Zc = sum(X[e][r] for e in METALS)
excess = X["H"][r] + X["He"][r] - 1
print("  catalog: sum of 9 species: median %.5f  (should be ~1)" % np.median(tot))
print("  catalog: (H + He - 1) / Z_metals: median %.3f  16-84%% %.3f - %.3f"
      % tuple(np.percentile(excess / Zc, [50, 16, 84])))
t9 = sum(tree[f"Sub_SubhaloGasMetalFractions_{i}"] for i in range(9))
t9 = t9[t9 > 0]
print("  tree:    sum of 9 species: median %.5f, range %.5f - %.5f"
      % (t9.median(), t9.min(), t9.max()))
tsum7 = sum(tree[f"Sub_SubhaloGasMetalFractions_{i}"] for i in range(2, 9))
gfm = tree["Sub_SubhaloGasMetallicity"]
ok = (gfm > 0) & (tsum7 > 0)
print("  tree:    GFM metallicity / sum(7 metals): median %.3f  (untracked metals)" % (gfm[ok] / tsum7[ok]).median())
fO = X["O"][r] / Zc
print("  O / Z_metals in clouds: median %.3f   (solar, Asplund 2009: 0.43 -> O-enhanced by %+.2f dex)"
      % (np.median(fO), np.log10(np.median(fO) / 0.43)))
oh = 12 + np.log10((X["O"][r] / 15.999) / (X["H"][r] / 1.008))
dz = oh - (8.69 + np.log10(Zc / 0.0127))
print("  12+log(O/H) - [8.69 + log(Z/Zsun)]: median %+.3f dex  (the O-enhancement shows up here)"
      % np.median(dz))

# ------------------------------------------------------------------ 4
section("4. SFR units")
with h5py.File("gmc_catalog/gmc_output.hdf5") as f:
    g = f["GMC_alpha20_n100"]
    Ms, Sy, Si = g["StellarMass"][...], g["SFR_yso"][...], g["SFR_inst"][...]
m = Sy > 0
print("  StellarMass / SFR_yso = %s exactly for all %d non-zero rows   %s"
      % (np.unique(np.round(Ms[m] / Sy[m], 6)), m.sum(), verdict(len(np.unique(np.round(Ms[m] / Sy[m], 6))) == 1)))
print("  -> SFR_yso/1000 = M(<3 Myr)/3 Myr in Msun/yr. SFR_inst vs host SubhaloSFR: see notebook 1b (sum <= host).")
d4 = tz_gmc.load(res=4, variant="fiducial")
gid = tree["Group_GroupID"].astype(int)
mm = d4[(d4.Run == "g205z4") & (d4.HaloID == d4.Snapshot.map(gid))]
per = mm.groupby("Snapshot")[["sfr_inst", "sfr_yso"]].sum()
per["host"] = tree.loc[per.index, "Sub_SubhaloSFR"]
q = (per.sfr_yso / per.sfr_inst).replace([np.inf, -np.inf], np.nan).dropna()
print("  per snapshot, sum SFR_yso/1000 / sum SFR_inst: median %.2f, 16-84%% %.2f - %.2f   %s"
      % (q.median(), q.quantile(.16), q.quantile(.84), verdict(0.2 < q.median() < 5)))
print("  (3 Myr average vs instantaneous: same order of magnitude expected, not equal)")

# ------------------------------------------------------------------ 5
section("5. CenterPosition units: physical kpc or comoving kpc/h?")
# Clouds are disjoint sets of cells, so two cloud centres should rarely sit
# closer than the sum of their radii. Under the wrong unit assumption the
# separations shrink (comoving/h -> physical multiplies by a/h ~ 0.15-0.37)
# and overlaps become common.
res = []
d4r = d4[d4.Number_of_Particles >= 30]
for (run, snap, hid), c in d4r.groupby(["Run", "Snapshot", "HaloID"]):
    if len(c) < 5:
        continue
    pos = c[["CenterPosition_x", "CenterPosition_y", "CenterPosition_z"]].to_numpy()
    R = c["EffectiveRadius"].to_numpy()
    af = 1 / (1 + c["z"].iloc[0])
    for lab, scale in [("physical kpc (assumed)", 1.0), ("comoving kpc/h", af / HUB)]:
        p = pos * scale
        dist, idx = cKDTree(p).query(p, k=2)
        overlap = dist[:, 1] < (R + R[idx[:, 1]])
        res.append((lab, overlap.mean()))
res = pd.DataFrame(res, columns=["assume", "frac"])
for lab, s in res.groupby("assume").frac:
    print("  if CenterPosition is %-24s nearest neighbour inside R_eff,i + R_eff,j: %.1f%% of clouds"
          % (lab, 100 * s.mean()))
# CloudPhinder clouds are disjoint sets of cells, but R_eff is a sqrt(5/3)
# r_rms size and the clouds are filamentary, so some overlap is fine.
fr = res.groupby("assume").frac.mean()
print("  %s: the physical reading gives far fewer overlapping clouds"
      % verdict(fr["physical kpc (assumed)"] < 0.5 * fr["comoving kpc/h"]))
print("  NOTE: indirect evidence; worth confirming with Zihao.")

# ------------------------------------------------------------------ 6
section("6. Host apertures: are they ordered sensibly?")
for key, suf in [("all", ""), ("halfrad", "HalfRad"), ("maxrad", "MaxRad"), ("sf", "Sfr"), ("sfrw", "SfrWeighted")]:
    z7 = sum(tree[f"Sub_SubhaloGasMetalFractions{suf}_{i}"] for i in range(2, 9)) / 0.0127
    zz = z7[z7 > 0]
    print("  %-8s median Z over snapshots at z<4: %.3f Zsun" % (key, zz[tree.loc[zz.index, "Redshift"] < 4].median()))
print("  expected: inner apertures (halfrad) most enriched, 'all' (incl. CGM) lower")

# ------------------------------------------------------------------ 7
section("7. The 'host outside clumps' subtraction")
print("  see notebook section 3: defined in 72-74 of 76 snapshots; checked here for negatives")
mm = mm[mm.Number_of_Particles >= 30]
sf = mm[mm.sfr_inst > 0]
w = sf.sfr_inst
pc = pd.DataFrame({"S": w.groupby(sf.Snapshot).sum(),
                   "XO": (w * sf.O_Fraction_sfr).groupby(sf.Snapshot).sum()})
hs = tree.loc[pc.index]
rest = hs["Sub_SubhaloSFR"] - pc.S
xo = (hs["Sub_SubhaloSFR"] * hs["Sub_SubhaloGasMetalFractionsSfrWeighted_4"] - pc.XO) / rest
okm = (rest > 0.1 * hs["Sub_SubhaloSFR"])
print("  snapshots with clumps < 90%% of host SFR: %d;  of those, residual O fraction <= 0: %d   %s"
      % (okm.sum(), (xo[okm] <= 0).sum(), verdict((xo[okm] <= 0).sum() == 0)))
