# THESAN-Zoom clumps: progress log

Running record of the science question, what has been established, the
decisions behind the analysis, and what is still open. Newest log entries at the
bottom of each dated section; numbers are quoted as measured so they can be
checked against `clumps.ipynb`.

---

## 1. Science question

JWST spectroscopy reports that star-forming clumps in high-redshift galaxies are
metal-poor relative to their host, by about **−0.1 dex** in 12+log(O/H)
(Estrada-Carpenter et al. 2025, ApJ 991, 188; 2025, ApJ 979, 14).

**Does THESAN-Zoom reproduce that offset, and if so, which clumps are diluted
and why?**

Background papers (in `Thesan_Zoom_Clumps/`):

- **Wang et al. 2025, MNRAS 544, 2675** — GMC-scale to galaxy-scale star
  formation efficiency. Produced the CloudPhinder catalog we use. Clouds are
  near-universal: Σ_gas ≈ 70 M⊙ pc⁻², ε_ff ≈ 3%, t_dep ≈ 100 Myr, mass
  function slope −2.49. The global depletion time is set by the GMC mass
  fraction, not by cloud-scale SFE.
- **Wang et al. (2026 preprint), clumpiness** — power spectra of mass and
  light. Stellar mass is smooth; UV/Hα light is clumpy, and clumpiness tracks
  the burst cycle, so observed clumpy galaxies may be biased toward bursts.
  Never identifies individual clumps.

The two papers do not connect: one catalogues clouds, the other measures
appearance. Our question sits between them.

---

## 2. Data in hand

| item | what it is |
|---|---|
| `gmc_catalog/gmc_output.hdf5` | CloudPhinder clouds, **4×** resolution, 8 targets, z ≈ 15 → 3 |
| `gmc_catalog/gmc_output_z8.hdf5` | **8×** resolution, 4 targets, z ≈ 9 → 3 (paper's fiducial for cloud properties) |
| `gmc_catalog/gmc_output_z16.hdf5` | **16×** resolution, 2 targets, z = 3.05 → 3.46 only |
| `gmc_catalog/gmc_out_{noESF,varSFE}z{4,8}.hdf5` | physics variants (no early stellar feedback; density-dependent SFE) |
| `g205_z4_host_halo.csv` | main-progenitor branch of target g205 (4×) from `tree.hdf5`; the **only** host tree we have locally |
| `inspect.output` | structure dump of the g205 group catalogs and tree |

Full snapshots, group catalogs and trees live on MIT ORCD
(`/orcd/data/mvogelsb/004/Thesan-Zooms/<target>/<zoom>/output`); everything
here is derived summary data. Snapshots are
`output/snapdir_<snap>/snapshot_<snap>.<chunk>.hdf5` (68 chunks for g205/z4);
group catalogs are `output/groups_<snap>/`, the main branch is
`output/tree.hdf5`.

Resolution sets the smallest cloud: median R_eff = 29.5 / 21.6 / 13.0 pc and
median M_gas = 6.7×10⁵ / 1.0×10⁵ / 1.3×10⁴ M⊙ at 4× / 8× / 16×.
Catalogs at different resolution are not directly comparable.

---

## 3. Data conventions (verified against the files)

These are traps that have already produced wrong results once. All are
implemented in `tz_gmc.py`.

- **`z4`/`z8`/`z16` are zoom factors, not redshifts.** Redshift comes from the
  snapshot number (189 snapshots, z = 16 → 3, same grid for every run).
- **`HaloID` is the FoF group index in that snapshot**, sorted by mass, so it
  changes over time. g205's main branch is *not* group 0 in 36 of 189
  snapshots (z = 6.6–12.5). Match with the tree's `Group_GroupID`.
- **Only centrals are stored.** When the host is a satellite (g205 at snapshot
  29, a FoF bridging event) it has no clouds.
- **`CenterPosition` is already relative to `GroupPos`.** Distance = its norm.
  Physical kpc: read that way, 3.6% of clouds overlap their nearest neighbour
  (within the sum of effective radii); read as comoving kpc/h, 54% would.
- **`SFR_yso` = `StellarMass / 3000` exactly**, i.e. 1000 × M⊙/yr. Divide by
  1000. `SFR_inst` is M⊙/yr.
- **H and He are normalised separately from the metals**: H + He ≈ 1 + Z, so
  `1 − H − He` returns −Z. Build Z as the sum of C, N, O, Ne, Mg, Si, Fe.
  The tree's species fractions do the same (sum 1.0002–1.014), so this is
  how the simulation stores them, not a catalog error; the effect on Z and
  O/H is below 1%.
- **Oxygen is enhanced relative to solar**: O is 57% of the metal mass in
  clouds vs 43% in the Sun, so 12+log(O/H) sits +0.09 dex above
  8.69 + log(Z/Z⊙). Compare O/H with O/H; converting between Z and O/H with a
  solar pattern introduces this offset.
- **`X_Fraction_sfr` is NaN** for clouds with no star-forming cell.
- **Group vs subhalo.** Group-level quantities (FoF) are fragile: at snapshot
  29 `GroupMass` triples and `GroupPos` jumps 59 kpc for one snapshot.
  Subhalo-level quantities are robust. "Halo mass" in the papers is
  `Group_M_Crit200`.
- **Tree metallicity columns**: `SubhaloGasMetallicity` is GFM total Z and runs
  ~5% above the seven-species sum for the same gas. Build host Z from
  `SubhaloGasMetalFractions*` so clump and host match.
- **Particle types** (NTYPES = 7): 0 gas, 1 DM, 4 stars; 2, 3, 5, 6 empty for
  g205 (no low-res contamination, no black holes). Type 4 is stars only: the
  run has no `GFM_WINDS` (TNG-style wind particles) and no black holes;
  feedback is explicit (`FM_STAR_FEEDBACK`, `FM_DISCRETE_SN`, early stellar
  feedback, radiative transfer). `FM_AGB_WINDS` is AGB mass return, not a
  particle type. Confirmed in g205/z4 snapshot 188 on ORCD: 5,170,761
  type-4 particles, none with `GFM_StellarFormationTime <= 0`. The earliest
  star formed at a = 0.0422 (z ≈ 22.7), before the first snapshot.

---

## 4. Decisions

| date | decision | reason |
|---|---|---|
| 2026-09-10 | Load catalogs through `tz_gmc.py`, which attaches redshift and physical units | Units in the readme are wrong for `SFR_yso`, and Z needs care |
| 2026-09-10 | Resolved clouds only: `Number_of_Particles >= 30` | Readme and paper resolution limit; removes ~43% of 4× rows |
| 2026-11-08 | Statistics per snapshot, not per clump | Clumps in one snapshot share a gas reservoir; per-clump medians understate scatter |
| 2026-11-08 | Radial trends at fixed redshift, distance as d / R200c | Pooling redshifts turns time evolution into a fake radial gradient |
| 2026-11-13 | Host metallicity from the tree, never from the clumps | A mean over clumps makes ΔZ zero by construction |
| 2026-11-13 | Clump-host match: `HaloID == Group_GroupID` **and** `Sub_SubhaloRankInGr == 0` | Verified by the section 1b tests below |
| 2026-11-13 | Z⊙ = 0.0127 for clumps and host; O/H by number from O and H mass fractions | Model value and tree `phys_*_Zsun` convention; one value for both sides |
| 2026-11-13 | Report five host apertures; headline = SFR-weighted clumps vs host star-forming gas *outside* the clumps | Closest analogue to emission-line measurements; to be revisited once the observational host definition is known |
| 2026-11-13 | Host comparison done at 4× | g205 is the only local tree, and it is 4× |
| 2026-11-13 | `clumps.ipynb` is the working analysis file | Scripts superseded (see §8) |
| 2026-11-13 | Settle catalog questions from the data, not by asking the catalog's author | Treat `CenterPosition` as physical kpc (overlap test, §5.5); treat missing snapshots as unprocessed gaps; test the burst-window question on ORCD |
| 2026-11-13 | Project lives in its own repo: `~/GitKraken/clumps`, remote `changhyun-cho/clumps` | Too much unrelated material in the `astrophysics` repo; history restarted there |

---

## 5. Results so far

### 5.1 Catalog reproduces Wang et al. 2025 (8×, `gmc_validate.py`)

| quantity | ours | paper |
|---|---|---|
| median M_gas | 9.9×10⁴ M⊙ | ~10⁵ |
| median α | 4.5 | ~4 |
| median Σ_gas, fiducial / noESF | 68.1 / 52.4 M⊙ pc⁻² | ~70 / ~51 |
| median cloud t_dep | 105 Myr | ~100 |
| Σ_SFR ∝ Σ_gas^b | b = 0.95 | 1 |
| mass function slope (4× / 8× / 16×) | −2.39 / −2.67 / −2.59 | −2.49 |
| noESF SFR at fixed Σ, same halo, 8× | 2.6× | ~3× |

Limits: ε_ff from catalog data is an upper bound (0.056 vs 0.03) because only
the mean-density t_ff is available; the paper's value implies a cell-level
⟨t_ff⟩ ≈ 3.1 Myr. σ_v is inverted from α, so the Larson slope (0.41 vs 0.6) is
only indicative.

### 5.2 Clump catalog ↔ host tree matching is verified (g205, notebook §1b)

All clouds, unresolved included:

- Matched clouds never exceed the host's gas (max 1.6%) or the group's SFR.
  Their SFR exceeds the central subhalo's at 2 snapshots (98, 99; by 1% and
  10%), so clouds include a little group gas not bound to the central.
- Other HaloIDs exceed the host SFR in 35 of 231 cases, up to 47×: other
  galaxies.
- In the 36 reshuffled snapshots, the tree's GroupID passes (max 0.82 of host
  SFR, metallicity ratio 1.02); `HaloID == 0` fails (up to 47×).
- Matched clouds track the bursty host SFR, Spearman ρ = 0.97; others 0.11.
- No host clouds while the host is a satellite.

### 5.3 Catalog coverage is incomplete (g205)

- No clouds for any halo at snapshots 0–24, 46–47, 63–66, 147, 174–175.
  Snapshots 174–175 have host SFR 44–55 M⊙/yr, so they were not processed.
- Of 155 processed snapshots, g205 has clouds in 79. At the others the
  catalog holds clouds of other, small haloes only, so the host was not
  processed there.
- **The gaps are irregular.** Neither gas content nor burst timing predicts
  them: 72 of the 75 snapshots without host clouds have more star-forming
  cells than the poorest snapshot *with* clouds (28 cells); the catalog stores
  clouds down to 4 elements. Reproducing the paper's starburst windows (local
  SFR peak over one t_dyn = R200/V200 ≈ 100–300 Myr, > 2× the window median)
  puts 130 of 155 processed snapshots inside a window, and host clouds appear
  in 54% inside vs 36% outside. Snapshot 139 is itself a starburst peak
  (43 M⊙/yr, ~4,600 star-forming cells) with no host clouds.
- **Does it bias the result?** Snapshots with host clouds lean to higher SFR
  (median 4.3 vs 1.4 M⊙/yr). The headline offset (star-forming clumps vs host
  star-forming gas outside clumps) does not depend on burst phase (Spearman
  ρ = −0.09, p = 0.43; median +0.010 dex in the burstier half, +0.023 in the
  calmer half). The offset against *all* host gas does (ρ = +0.43,
  p < 0.001), so that aperture is sensitive to the uneven sampling.

### 5.4 Clump vs host metallicity (g205, 4×, notebook §3–7)

4,189 resolved main-branch clumps (3,251 star-forming) in 76 snapshots,
z = 3.0–9.6.

Δ12+log(O/H), clump − host, median over snapshots [16th, 84th]:

| host aperture | Δ (dex) |
|---|---|
| all subhalo gas | +0.002 [−0.054, +0.064] |
| within stellar r½ | −0.016 [−0.068, +0.019] |
| star-forming gas (mass-weighted) | +0.043 [−0.012, +0.078] |
| SFR-weighted | +0.033 [−0.009, +0.063] |
| **SFR-weighted, outside clumps** | **+0.012 [−0.018, +0.063]** |

- **No −0.1 dex dilution in the median.** But the aperture choice moves Δ by
  ~0.06 dex, comparable to the effect being tested.
- 11.2% of star-forming clumps sit ≥ 0.1 dex below the host's star-forming gas
  outside clumps (per snapshot: median 8.5%, 0–18%).
- Diluted clumps are not distinctive: slightly more gas and SFR, same radius.
  Within-snapshot rank correlations of clump O/H with mass, SFR, α, R_eff,
  Σ_gas, density and radius all have |ρ| < 0.1.
- Clumps in the same galaxy at the same time are chemically similar:
  16th–84th spread −0.07 to +0.06 dex.
- The low-metallicity tail grows with resolution: 13% / 17% / 20% of clumps
  sit > 0.1 dex below their snapshot median at 4× / 8× / 16×. Different
  targets and redshifts, so not yet a clean resolution test.
- Radial (fixed z): O/H flat inside ~0.1 R200c, falls in the outskirts
  (z = 4–5: Z 0.29 → 0.18 Z⊙ beyond 0.2 R200c). The virial parameter falls
  outward (α ≈ 10 → 6.4): inner clouds are less self-gravitating.
- Resolved clumps carry a median 47% of the host SFR, and their summed SFR
  follows host bursts with no lag beyond one snapshot (~10 Myr).

### 5.5 Audit of the calculations (2026-11-13, `audit_checks.py`)

| check | result |
|---|---|
| snapshot → z map vs group-catalog headers (snap 0, 100) | exact |
| tree M200c vs 200 ρ_crit(z) 4/3 π R200c³ (Planck 2015) | ratio 0.9998, all snapshots: mass, length, a and h conversions right |
| `SFR_yso` vs `StellarMass` | exactly /3000 in all 34,311 non-zero rows |
| Σ SFR_yso/1000 vs Σ SFR_inst per snapshot | ratio 0.55 (0.37–0.88): same order, as expected for 3 Myr vs instantaneous |
| `CenterPosition` units | physical kpc strongly favoured (above) |
| host "outside clumps" subtraction | never negative where defined |
| species sums, O enhancement | see §3 |

Not a bug but worth knowing: at z < 4 the host's star-forming gas
(Z ≈ 0.33–0.34 Z⊙ mass- or SFR-weighted) is *less* enriched than its gas
within r½ (0.42) or all its bound gas (0.39). The bound gas outside the
star-forming phase is the more metal-rich, plausibly metal-loaded outflows.
An emission-line host value would track the star-forming gas.

Not checkable locally: the snapshot grid for runs other than g205 (assumed
identical); exact meaning of the `HalfRad` aperture (TNG convention: within
the stellar half-mass radius).

### 5.6 Corrections (withdrawn statements)

Kept here so they are not repeated.

- "Clumps are slightly *more* enriched than the host at z = 4–6" — withdrawn.
  It came from per-clump medians against one aperture; per snapshot and
  across apertures the sign is not determined.
- "The cloud population only exists during starbursts" — withdrawn. May be
  the catalog's snapshot selection (§5.3), not physics.
- The exact cause of an earlier, inconsistent clump-host metallicity offset
  (values of −0.3 to −0.95) was never established from figures alone. What is
  established is the correct result in §5.4.

---

## 6. Open questions

**About the catalog** (to be settled from the data)
1. ~~Were clouds searched for only inside starburst windows?~~ No: the gaps
   are irregular and not explained by gas or by the paper's windows; they do
   not bias the headline offset (§5.3).
2. ~~Why are snapshots missing?~~ Treated as unprocessed gaps (the host forms
   up to 55 M⊙/yr at 174–175).
3. ~~Is `CenterPosition` physical kpc?~~ Yes, by the overlap test (§5.5).

**About the observations**
4. How do Estrada-Carpenter et al. define the host ("disk") metallicity:
   integrated spectrum including the clumps, inter-clump light, something
   else? Which line diagnostic, and at what spatial resolution? This decides
   the matching host aperture.
5. What are the sizes and masses of the observed clumps? CloudPhinder clouds
   at 4× are ~6×10⁵ M⊙ and ~30 pc; JWST clumps are likely 10⁷ M⊙-scale
   complexes of hundreds of pc. A fair comparison may need clouds aggregated
   into complexes.

**Science**
6. If THESAN-Zoom clumps are not diluted, is it because these galaxies have
   no settled disc and clouds condense from an already well-mixed turbulent
   ISM? Observed dilution would then point to clumps fed by metal-poor inflow.
   Test: are the diluted ~11% associated with recent inflow (needs snapshots).
7. Does the low-metallicity tail really grow with resolution? Needs
   same-target pairs (g33206, g5760, g10304 exist at 4× and 8×).

---

## 7. Next steps

- [x] Catalog coverage: gaps are irregular, headline offset unaffected (§5.3).
- [ ] Read the observational papers for the host definition and clump scale
      (questions 4–5); add them to `Thesan_Zoom_Clumps/`.
- [ ] On ORCD: run `tz_tree_catalog.py` for the other targets and save
      `<target>_z<zoom>_host_halo.csv` here; notebook §1b checks each one
      automatically.
- [ ] Clump complexes: link clouds within ~200–500 pc and repeat §3–5.
- [ ] Same-target resolution comparison of clump metallicity scatter.
- [ ] Fix the stellar-type detection in `tz_plot.py` (picks the empty type 2
      instead of 4, so the stellar curves silently vanish from
      `g205_z4_mass.png` and `g205_z4_size.png`).
- [x] Keep the catalog, papers and audit out of the public repo (§9).
- [x] Decide on the host CSV, `inspect.output` and the notebook backup (§9): kept private.

---

## 8. File map

| file | role |
|---|---|
| `clumps.ipynb` | **working analysis** (sections 1–8, host-matching tests in §1b) |
| `tz_gmc.py` | catalog loader with units, redshift and derived quantities |
| `gmc_validate.py`, `gmc_overview.png` | reproduction of Wang et al. 2025 numbers |
| `audit_checks.py` | independent checks of units, conventions and derived quantities (§5.5) |
| `clump_host_analysis.py`, `indep_*.png` | independent host-vs-clump analysis (superseded by the notebook) |
| `tz_tree_catalog.py`, `tz_inspect.py`, `tz_plot.py`, `sanity.py` | tree extraction and inspection (run on ORCD) |
| `g205_z4_host_halo.csv`, `.units.txt` | g205 main-branch tree, flattened |
| `figure 6*.py`, `correlation_*.py`, `reproduced_figure6.png` | early versions; host defined from the clumps themselves and redshift not separated. **Superseded** by notebook §4 and §7 |
| `student_audit.py`, `student_audit.png` | review of an earlier analysis |
| `clumps_backup_2026-11-13.ipynb` | notebook before the 2026-11-13 rebuild |

---

## 9. Repository note

Repo: `~/GitKraken/clumps`, remote `changhyun-cho/clumps` (**public**).
The clumps work started in `changhyun-cho/astrophysics` (`clumps/`, commits
b32d0b5, 292965b, eb538f5); those old versions stay in that repo's history.

**Ignored** (`.gitignore`): `gmc_catalog/` files (unreleased catalog), the two
papers in `Thesan_Zoom_Clumps/` (incl. the unpublished preprint),
`student_audit.*`, `g205_z4_host_halo.csv` and `.units.txt`, `inspect.output`,
`clumps_backup_2026-11-13.ipynb`, `.DS_Store`, the `tz_gmc` bytecode file.
The rules name files individually, so a new catalog file, paper or `.pyc`
must be added by hand. Checked 2026-11-13: no private file in any commit.

Cloning the repo does not give a runnable notebook: the catalog and host CSV
have to be copied in separately.

---

## Log

- **2026-09-10** — Read both papers; inventoried the catalogs; found the
  `SFR_yso` and metallicity traps; wrote `tz_gmc.py`; reproduced Wang et al.
  2025 at 8× (`gmc_validate.py`). First look at clump radial distribution in
  g205.
- **2026-11-04** — Reviewed an earlier g205 clump-host analysis; it gave
  negative metallicities and clump distances of 8–20 Mpc.
- **2026-11-08** — Independent re-analysis (`clump_host_analysis.py`): clump
  and host metallicity agree within ~±0.07 Z⊙; host SFR bursts are tracked by
  clump SFR.
- **2026-11-13** — Reviewed `clumps.ipynb` and the Figure 6 / correlation
  scripts; host had been built from the clumps. Rebuilt the notebook with the
  host from the tree, five apertures, per-snapshot statistics. Added the
  host-matching tests (§1b) and the rank-0 requirement; all pass for g205.
  Found the snapshot-coverage gaps. Moved the project into the git repo.
  Started this log. Set up `.gitignore` for the private material. Moved the
  project into its own public repo `changhyun-cho/clumps`; verified that no
  private file was committed and that the notebook runs there (all §1b tests
  pass). Audited all derived quantities (`audit_checks.py`, §5.5): no
  errors found; documented the oxygen enhancement and the host's
  star-forming gas being less enriched than its bound gas. Confirmed on ORCD
  that particle type 4 holds stars only (no wind particles). Decided to
  settle catalog questions from the data; the coverage gaps turned out to be
  irregular (not starburst windows) and do not bias the headline offset.
