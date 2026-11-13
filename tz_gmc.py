#!/usr/bin/env python3
"""
tz_gmc.py -- load the THESAN-ZOOM CloudPhinder GMC catalogs with physical units attached.

The catalogs in gmc_catalog/ are flat tables of clouds identified by CLOUDPHINDER
(Guszejnov et al. 2019; as refined in Wang et al. 2025, MNRAS 544, 2675). One row
is one cloud in one snapshot of one run. The HDF5 group name records the
identification setup: 'GMC_alpha20_n100' means alpha_max = 20 and n_H > 100 cm^-3.

    import tz_gmc
    d = tz_gmc.load()                 # everything, all runs, all resolutions
    d = tz_gmc.load(res=8, resolved_only=True)

UNITS, as stored in the files (readme.txt is right except where noted):

    CenterPosition      kpc, physical, relative to the host GroupPos
    EffectiveRadius     kpc, physical      (see R_eff_pc below)
    GasMass, StellarMass    Msun
    SFR_inst            Msun/yr
    SFR_yso             NOT Msun/yr.  Verified against the file: SFR_yso is
                        exactly StellarMass / 3000 for every row, i.e. the
                        3 Myr normalisation was done in kyr, not yr. The stored
                        number is 1000x the SFR in Msun/yr. Cross-check: the
                        median SFR_yso/SFR_inst ratio for clouds with both
                        non-zero is 825, i.e. ~1000. `sfr_yso` below applies the
                        /1000 correction; the raw column is kept as
                        `SFR_yso_raw`.
    X_Fraction          mass fraction of element X, mass-weighted over gas
    X_Fraction_sfr      the same, weighted by instantaneous SFR. NaN for the
                        ~85 per cent of clouds with no star-forming cell.
    VirialParameter     alpha = 2(E_kin + E_th)/|E_grav|, capped at 20 by
                        construction
    Number_of_Particles total gas cells + young star particles. Clouds with
                        < 30 are unresolved and should be cut (readme, and
                        Sec. 4.1 of the paper).

The nine element fractions sum to 1.0014, not 1: H and He are stored normalised
separately from the metals, so 1 - H - He is slightly negative and useless. Use
the explicit metal sum, which is what `Z` below does.

RUN NAMES: 'g205z4' is target halo g205 at the 4x resolution level. The trailing
z4/z8/z16 is the ZOOM FACTOR, not a redshift. The resolution sets the smallest
resolvable cloud, so the catalogues are not directly comparable across levels:
median R_eff is 29 / 21 / 13 pc and median M_gas 6.5e5 / 1.0e5 / 1.3e4 Msun at
4x / 8x / 16x. The paper adopts 8x as fiducial for cloud-scale work.
"""

import glob
import os

import h5py
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CATDIR = os.path.join(HERE, "gmc_catalog")
HOST_TREE_DIR = os.path.join(HERE, "host_trees")


def host_tree_path(tag):
    """Path of a host tree, '<target>_z<zoom>_host_halo.csv' (tag 'g205_z4' or the
    full file name), looked up in host_trees/ first and then the repo root."""
    name = tag if tag.endswith(".csv") else "%s_host_halo.csv" % tag
    for d in (HOST_TREE_DIR, HERE):
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    raise FileNotFoundError("%s not found in %s or %s" % (name, HOST_TREE_DIR, HERE))

ZSUN = 0.0127                     # Asplund+2009, as used by the Illustris/THESAN model
METALS = ["C", "N", "O", "Ne", "Mg", "Si", "Fe"]
ELEMENTS = ["H", "He"] + METALS   # GFM species order
G_CGS = 6.674e-8
MSUN_PER_PC3_TO_CGS = 6.770e-23   # 1 Msun/pc^3 in g/cm^3
G_PC_MSUN_KMS = 4.30091e-3        # G in pc Msun^-1 (km/s)^2
SEC_PER_MYR = 3.1557e13

# file -> (resolution level, physics variant)
FILES = {
    "gmc_output.hdf5":        (4,  "fiducial"),
    "gmc_output_z8.hdf5":     (8,  "fiducial"),
    "gmc_output_z16.hdf5":    (16, "fiducial"),
    "gmc_out_noESFz4.hdf5":   (4,  "noESF"),
    "gmc_out_noESFz8.hdf5":   (8,  "noESF"),
    "gmc_out_varSFEz4.hdf5":  (4,  "varSFE"),
    "gmc_out_varSFEz8.hdf5":  (8,  "varSFE"),
}


def snapshot_redshifts(path=None):
    """snap -> z. THESAN-ZOOM writes 189 snapshots from z=16 to z=3 on a common
    time grid, so the map taken from the g205 tree applies to every run."""
    t = pd.read_csv(path or host_tree_path("g205_z4"), usecols=["snap", "Redshift"])
    return dict(zip(t["snap"].astype(int), t["Redshift"].astype(float)))


def _read(path, group="GMC_alpha20_n100"):
    with h5py.File(path, "r") as h:
        if group not in h:
            raise KeyError("%s has no group %s (has %s)"
                           % (path, group, list(h.keys())))
        g = h[group]
        cols = {}
        for k in g:
            a = g[k][...]
            if a.ndim == 2:
                for i, c in enumerate("xyz"):
                    cols["%s_%s" % (k, c)] = a[:, i]
            elif a.dtype == object:
                cols[k] = a.astype(str)
            else:
                cols[k] = a
    return pd.DataFrame(cols)


def derive(d, zmap=None):
    """Attach redshift and the physical quantities the papers actually use."""
    if zmap is None:
        zmap = snapshot_redshifts()
    d = d.copy()
    d["z"] = d["Snapshot"].map(zmap)

    d["R_eff_pc"] = d["EffectiveRadius"] * 1e3
    d["M_gas"] = d["GasMass"]
    d["M_star"] = d["StellarMass"]
    d["M_tot"] = d["GasMass"] + d["StellarMass"]
    d["alpha"] = d["VirialParameter"]

    # SFR_yso is stored as StellarMass/3000; see the module docstring.
    d["SFR_yso_raw"] = d["SFR_yso"]
    d["sfr_yso"] = d["SFR_yso"] / 1e3                      # Msun/yr
    d["sfr_inst"] = d["SFR_inst"]                          # Msun/yr

    # Surface density from the effective radius, as in Sec. 4.2: Sigma = M/(pi R^2).
    # The clouds are filamentary, so this is a compactness measure, not a
    # line-of-sight column density.
    d["Sigma_gas"] = d["M_gas"] / (np.pi * d["R_eff_pc"] ** 2)     # Msun/pc^2

    # Free-fall time at the cloud mean gas density. The paper instead uses the
    # mass-weighted harmonic mean of the per-cell t_ff, which is shorter because
    # dense cells dominate it; this estimate is therefore an upper bound and the
    # eps_ff below a corresponding upper bound.
    rho = d["M_gas"] / (4.0 / 3.0 * np.pi * d["R_eff_pc"] ** 3)    # Msun/pc^3
    d["rho_gas"] = rho
    d["n_H"] = rho * MSUN_PER_PC3_TO_CGS * 0.76 / 1.6726e-24       # cm^-3, X=0.76
    t_ff_s = np.sqrt(3 * np.pi / (32 * G_CGS * rho * MSUN_PER_PC3_TO_CGS))
    d["t_ff_Myr"] = t_ff_s / SEC_PER_MYR

    # eps_ff = SFR * t_ff / M_gas  (Eq. 11)
    d["eps_ff"] = d["sfr_yso"] * (t_ff_s / 3.1557e7) / d["M_gas"]

    # Velocity dispersion inverted from the virial parameter under the uniform
    # sphere approximation alpha = 5 sigma^2 R / (G M). CloudPhinder computes
    # E_grav from the exact particle sum, so this is indicative, not exact.
    d["sigma_v"] = np.sqrt(d["alpha"] * G_PC_MSUN_KMS * d["M_tot"]
                           / (5 * d["R_eff_pc"]))                  # km/s

    # Metallicity. H and He are normalised separately from the metals in these
    # files, so sum the metals explicitly rather than taking 1 - H - He.
    d["Z"] = sum(d["%s_Fraction" % e] for e in METALS)
    d["Z_Zsun"] = d["Z"] / ZSUN
    d["Z_sfr"] = sum(d["%s_Fraction_sfr" % e] for e in METALS)
    d["Z_sfr_Zsun"] = d["Z_sfr"] / ZSUN
    d["OFe"] = d["O_Fraction"] / d["Fe_Fraction"].replace(0.0, np.nan)

    d["d_center_kpc"] = np.sqrt(sum(d["CenterPosition_%s" % c] ** 2
                                    for c in "xyz"))
    return d


def load(res=None, variant=None, files=None, resolved_only=False,
         min_particles=30, group="GMC_alpha20_n100"):
    """Concatenate catalogues into one frame.

    res      4, 8, 16 or None for all resolution levels
    variant  'fiducial', 'noESF', 'varSFE' or None for all
    """
    zmap = snapshot_redshifts()
    want = files or sorted(FILES)
    out = []
    for fn in want:
        if fn not in FILES:
            raise KeyError("unknown catalogue %s" % fn)
        lvl, var = FILES[fn]
        if res is not None and lvl != res:
            continue
        if variant is not None and var != variant:
            continue
        path = os.path.join(CATDIR, fn)
        if not os.path.exists(path):
            continue
        d = derive(_read(path, group), zmap)
        d["res"] = lvl
        d["variant"] = var
        d["catalog"] = fn
        out.append(d)
    if not out:
        raise RuntimeError("no catalogues matched res=%s variant=%s" % (res, variant))
    d = pd.concat(out, ignore_index=True)
    if resolved_only:
        d = d[d["Number_of_Particles"] >= min_particles].reset_index(drop=True)
    return d


def available():
    """What is on disk, with the group names each file carries."""
    rows = []
    for fn, (lvl, var) in sorted(FILES.items()):
        path = os.path.join(CATDIR, fn)
        if not os.path.exists(path):
            continue
        with h5py.File(path, "r") as h:
            groups = list(h.keys())
            n = h[groups[0]]["Run"].shape[0]
            runs = sorted(set(h[groups[0]]["Run"][...].astype(str)))
        rows.append({"file": fn, "res": lvl, "variant": var, "N": n,
                     "groups": ",".join(groups), "runs": ",".join(runs)})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    pd.set_option("display.width", 200)
    print(available().to_string(index=False))
    d = load(resolved_only=True)
    print("\n%d resolved clouds across %d catalogues" % (len(d), d.catalog.nunique()))
    print(d.groupby(["res", "variant"]).agg(
        N=("z", "size"), z_max=("z", "max"), z_min=("z", "min"),
        R_pc=("R_eff_pc", "median"), M_gas=("M_gas", "median"),
        Sigma=("Sigma_gas", "median"), alpha=("alpha", "median")).to_string())
