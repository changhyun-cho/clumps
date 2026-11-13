#!/usr/bin/env python3
"""
tz_tree_catalog.py -- host-halo catalog for a Thesan-Zoom run, straight from tree.hdf5.

tree.hdf5 in each run's output directory is not a generic merger tree. It is the
main-progenitor branch of the target halo, already extracted: every Group/ and
Subhalo/ dataset has leading dimension 189, one entry per snapshot, alongside
/Snapshots and /Redshifts. So the whole job is a flatten and a write.

    python tz_tree_catalog.py /orcd/data/mvogelsb/004/Thesan-Zooms/g205/z4/output \\
        -o g205_z4_host_halo.csv

GroupID and SubhaloID in the tree are the indices that object occupies in that
snapshot's fof_subhalo_tab catalog. They are NOT always 0. Use them to index the
full catalogs rather than assuming index 0.
"""

import argparse
import glob
import os
import sys

import h5py
import numpy as np

ZSUN = 0.0127  # Asplund+2009, the value used in the Illustris/TNG/Thesan model


def find_hubble(outdir, tree):
    """h lives in the Parameters group of the catalogs, not in the tree header."""
    for src in (tree.get("Header"), tree):
        if src is not None and "HubbleParam" in src.attrs:
            return float(src.attrs["HubbleParam"]), "tree"
    cand = sorted(glob.glob(os.path.join(outdir, "groups_*", "fof_subhalo_tab_*.0.hdf5")))
    for path in cand[:1]:
        with h5py.File(path, "r") as h:
            for grp in ("Parameters", "Header"):
                if grp in h and "HubbleParam" in h[grp].attrs:
                    return float(h[grp].attrs["HubbleParam"]), os.path.basename(path)
    return None, None


def flatten(prefix, node, n):
    cols = {}
    for key in node:
        item = node[key]
        if not isinstance(item, h5py.Dataset) or item.shape[0] != n:
            continue
        arr = np.asarray(item[...])
        if arr.ndim == 1:
            cols["%s%s" % (prefix, key)] = arr
        else:
            flat = arr.reshape(n, -1)
            for i in range(flat.shape[1]):
                cols["%s%s_%d" % (prefix, key, i)] = flat[:, i]
    return cols


# Comoving lengths, listed explicitly. A name-matching heuristic gets this wrong:
# SubhaloGasMetallicityHalfRad and friends carry "Rad" but are metallicities
# measured inside an aperture, not lengths.
LENGTH_FIELDS = {
    "GroupPos", "Group_R_Crit200", "Group_R_Crit500", "Group_R_Mean200",
    "Group_R_TopHat200", "SubhaloPos", "SubhaloCM", "SubhaloHalfmassRad",
    "SubhaloHalfmassRadType", "SubhaloVmaxRad", "SubhaloStellarPhotometricsRad",
}


def base_name(col):
    """Strip the Group_/Sub_ prefix and any _<i> component suffix flatten added."""
    for p in ("Group_", "Sub_"):
        if col.startswith(p):
            col = col[len(p):]
            break
    head, _, tail = col.rpartition("_")
    return head if head and tail.isdigit() else col


def add_physical(cols, hubble):
    a = 1.0 / (1.0 + cols["Redshift"])
    derived = {}
    for c, v in list(cols.items()):
        if c in ("snap", "Redshift"):
            continue
        if not np.issubdtype(v.dtype, np.number):
            continue
        b = base_name(c)
        is_mass = ("Mass" in b or b.startswith("Group_M_")) and "Metallicity" not in b
        is_len = b in LENGTH_FIELDS
        is_z = "Metallicity" in b and "Dust" not in b
        if is_mass:
            derived["phys_%s_Msun" % c] = v * 1e10 / hubble
        elif is_len:
            derived["phys_%s_kpc" % c] = v * a / hubble
        elif is_z:
            derived["phys_%s_Zsun" % c] = v / ZSUN
    return derived


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir", help="run output directory containing tree.hdf5")
    ap.add_argument("-o", "--output", default="host_halo.csv")
    ap.add_argument("--hubble", type=float, default=None,
                    help="override h if it cannot be found automatically")
    ap.add_argument("--no-derived", action="store_true")
    args = ap.parse_args()

    tree_path = (args.outdir if args.outdir.endswith(".hdf5")
                 else os.path.join(args.outdir, "tree.hdf5"))
    outdir = os.path.dirname(tree_path)
    if not os.path.exists(tree_path):
        sys.exit("no tree.hdf5 at %s" % tree_path)

    with h5py.File(tree_path, "r") as t:
        snaps = np.asarray(t["Snapshots"][...])
        n = snaps.shape[0]
        cols = {"snap": snaps, "Redshift": np.asarray(t["Redshifts"][...])}
        for name, prefix in (("Group", "Group_"), ("Subhalo", "Sub_")):
            if name in t:
                cols.update(flatten(prefix, t[name], n))
        hubble, src = (args.hubble, "--hubble") if args.hubble else find_hubble(outdir, t)

    derived = {}
    if not args.no_derived:
        if hubble is None:
            print("WARNING: could not determine h; skipping phys_ columns. "
                  "Pass --hubble to enable them.", file=sys.stderr)
        else:
            derived = add_physical(cols, hubble)
    cols.update(derived)

    names = list(cols.keys())
    with open(args.output, "w") as f:
        f.write(",".join(names) + "\n")
        for i in range(n):
            f.write(",".join(repr(float(cols[c][i])) for c in names) + "\n")

    with open(args.output + ".units.txt", "w") as f:
        f.write("Source: %s\n" % tree_path)
        f.write("h = %s (from %s)\n\n" % (hubble, src))
        f.write("Raw columns are AREPO code units:\n")
        f.write("  mass         1e10 Msun/h\n")
        f.write("  length       comoving kpc/h\n")
        f.write("  SFR          Msun/yr\n")
        f.write("  metallicity  M_Z/M_tot mass fraction, not solar-normalised\n")
        f.write("  *MetalFractions_i  mass fraction of element i, GFM order:\n")
        f.write("     0 H, 1 He, 2 C, 3 N, 4 O, 5 Ne, 6 Mg, 7 Si, 8 Fe\n")
        f.write("     (confirm against your run's GFM species list)\n\n")
        f.write("Derived:\n")
        f.write("  *_Msun = raw * 1e10 / h\n")
        f.write("  *_kpc  = raw * a / h,  a = 1/(1+z)\n")
        f.write("  *_Zsun = raw / %s\n\n" % ZSUN)
        for c in names:
            f.write("  %s%s\n" % (c, "   [derived]" if c in derived else ""))

    print("wrote %s  (%d snapshots x %d columns)" % (args.output, n, len(names)))
    print("wrote %s.units.txt" % args.output)

    # The point of the exercise: index 0 is not always the target.
    gid = cols.get("Group_GroupID")
    sid = cols.get("Sub_SubhaloID")
    rank = cols.get("Sub_SubhaloRankInGr")
    if gid is not None:
        off = np.where(gid != 0)[0]
        print("\nsnapshots where the target is NOT FoF group 0: %d / %d" % (len(off), n))
        for i in off:
            print("  snap %3d  z=%6.3f  GroupID=%d  SubhaloID=%d"
                  % (cols["snap"][i], cols["Redshift"][i], gid[i],
                     sid[i] if sid is not None else -1))
    if sid is not None:
        offs = np.where(sid != 0)[0]
        print("snapshots where the target is NOT subhalo 0: %d / %d" % (len(offs), n))
    if rank is not None:
        bad = np.where(rank != 0)[0]
        print("snapshots where the target is not the central of its group: %d / %d"
              % (len(bad), n))


if __name__ == "__main__":
    main()