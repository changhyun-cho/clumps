#!/usr/bin/env python3
"""
tz_inspect.py -- dump the structure of Thesan-Zoom group catalogs and the merger tree.

Read-only. Touches nothing. Run it before the extraction script so field names,
shapes, dtypes and unit attributes are confirmed rather than assumed.

    python tz_inspect.py /orcd/data/mvogelsb/004/Thesan-Zooms/g205/z4/output

Optionally give the snapshots to inspect (default: last, middle, first non-empty):

    python tz_inspect.py <output_dir> --snaps 188 100 010
"""

import argparse
import glob
import os
import sys

import h5py
import numpy as np


def chunk_files(outdir, snap):
    """All HDF5 chunks of one group catalog, in chunk order."""
    d = os.path.join(outdir, "groups_%03d" % snap)
    if not os.path.isdir(d):
        return []
    files = glob.glob(os.path.join(d, "*.hdf5"))

    def chunk_index(p):
        # .../fof_subhalo_tab_188.3.hdf5 -> 3 ; single-file catalogs -> 0
        base = os.path.basename(p)[:-5]
        tail = base.rsplit(".", 1)[-1]
        return int(tail) if tail.isdigit() else 0

    return sorted(files, key=chunk_index)


def describe_attrs(obj, indent):
    for k in obj.attrs:
        v = obj.attrs[k]
        if isinstance(v, np.ndarray):
            v = v.tolist()
        if isinstance(v, bytes):
            v = v.decode()
        print("%s%-28s = %s" % (indent, k, v))


def walk(node, indent="  ", max_depth=4, depth=0):
    """Print datasets with shape/dtype and any unit-ish attributes."""
    for key in node:
        item = node[key]
        path = item.name
        if isinstance(item, h5py.Dataset):
            units = {k: item.attrs[k] for k in item.attrs}
            unit_str = ""
            if units:
                bits = []
                for k, v in units.items():
                    if isinstance(v, np.ndarray):
                        v = v.tolist()
                    if isinstance(v, bytes):
                        v = v.decode()
                    bits.append("%s=%s" % (k, v))
                unit_str = "   [%s]" % ", ".join(bits)
            print("%s%-46s %-16s %s%s"
                  % (indent, path, str(item.shape), item.dtype, unit_str))
        else:
            print("%s%s/" % (indent, path))
            if depth + 1 < max_depth:
                walk(item, indent + "  ", max_depth, depth + 1)


def inspect_catalog(outdir, snap):
    files = chunk_files(outdir, snap)
    print("\n" + "=" * 90)
    print("GROUP CATALOG  snap %03d   (%d chunk file(s))" % (snap, len(files)))
    print("=" * 90)
    if not files:
        print("  no groups_%03d directory" % snap)
        return
    for f in files[:1] + ([files[-1]] if len(files) > 1 else []):
        print("\n-- %s" % f)
        with h5py.File(f, "r") as h:
            for grp in ("Header", "Parameters", "Config"):
                if grp in h:
                    print("\n  [%s attributes]" % grp)
                    describe_attrs(h[grp], "    ")
            print("\n  [datasets]")
            for top in h:
                if top in ("Header", "Parameters", "Config"):
                    continue
                print("  %s/" % top)
                walk(h[top], indent="    ")


def inspect_tree(outdir):
    path = os.path.join(outdir, "tree.hdf5")
    print("\n" + "=" * 90)
    print("MERGER TREE  %s" % path)
    print("=" * 90)
    if not os.path.exists(path):
        print("  not found")
        return
    with h5py.File(path, "r") as h:
        print("\n  [root attributes]")
        describe_attrs(h, "    ")
        print("\n  [structure]")
        walk(h, indent="    ", max_depth=4)

        # Show the first handful of values of every small 1-D integer dataset,
        # which is usually enough to recognise the tree convention in use
        # (LHaloTree / SubLink / rockstar-style descendant pointers).
        print("\n  [head of small integer datasets]")

        def peek(node):
            for key in node:
                item = node[key]
                if isinstance(item, h5py.Dataset):
                    if item.ndim == 1 and np.issubdtype(item.dtype, np.integer):
                        n = min(12, item.shape[0])
                        print("    %-46s %s ..." % (item.name, item[:n].tolist()))
                else:
                    peek(item)

        peek(h)


def first_nonempty(outdir):
    """Lowest snapshot number whose catalog contains at least one FoF group."""
    for snap in range(0, 189):
        files = chunk_files(outdir, snap)
        if not files:
            continue
        with h5py.File(files[0], "r") as h:
            ng = int(h["Header"].attrs.get("Ngroups_Total",
                     h["Header"].attrs.get("Ngroups_ThisFile", 0)))
        if ng > 0:
            return snap
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--snaps", nargs="*", type=int, default=None)
    args = ap.parse_args()

    if not os.path.isdir(args.outdir):
        sys.exit("no such directory: %s" % args.outdir)

    snaps = args.snaps
    if snaps is None:
        fne = first_nonempty(args.outdir)
        snaps = [s for s in (188, 100, fne) if s is not None]
        snaps = sorted(set(snaps), reverse=True)
        print("first snapshot with >=1 FoF group: %s" % fne)

    for s in snaps:
        inspect_catalog(args.outdir, s)

    inspect_tree(args.outdir)

    # Quick census of how Group 0's mass and position behave across time; a
    # discontinuity here is the cheap tell that index 0 is not one object.
    print("\n" + "=" * 90)
    print("GROUP 0 CENSUS (Ngroups, M_Crit200 code units, GroupPos)")
    print("=" * 90)
    print("%5s %6s %10s %14s %26s" % ("snap", "z", "Ngroups", "M200c", "GroupPos"))
    for snap in range(0, 189, 8):
        files = chunk_files(args.outdir, snap)
        if not files:
            continue
        with h5py.File(files[0], "r") as h:
            hd = h["Header"].attrs
            z = float(hd.get("Redshift", np.nan))
            ng = int(hd.get("Ngroups_Total", 0))
            if ng == 0 or "Group" not in h:
                print("%5d %6.2f %10d %14s %26s" % (snap, z, ng, "-", "-"))
                continue
            m200 = float(h["Group/Group_M_Crit200"][0]) if "Group_M_Crit200" in h["Group"] else np.nan
            pos = h["Group/GroupPos"][0]
            print("%5d %6.2f %10d %14.5e   [%9.3f %9.3f %9.3f]"
                  % (snap, z, ng, m200, pos[0], pos[1], pos[2]))


if __name__ == "__main__":
    main()
