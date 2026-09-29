# Quantum ESPRESSO decks — the DFT+U value oracle (`hp.x`)

`hp.x` (Timrov, Marzari, Cococcioni: DFPT linear-response Hubbard U) is the ORACLE for U values (doc/OpenWork.md
step 5, item 4).  An oracle is only as good as its match, so these decks run **our own pseudopotential and our own
Hubbard projector**: the UPF files are written by `CLIapps/gth2upf` from the GTH parameter database with
`PP_CHI` = the qchem GTH pseudo-atom's orbitals (the `HubbardU_Atomic` radial) — QE's `atomic` projector IS our χ.

```
build/Release/CLIapps/gth2upf --element Mn --q 7     # -> Mn.pz-gth-q7.UPF   (4s, 3d; E_atom -14.2442 Ha)
build/Release/CLIapps/gth2upf --element O  --q 6     # -> O.pz-gth-q6.UPF    (2s, 2p; E_atom -15.7484 Ha)
mpirun -np 1 ~/Code/q-e/bin/pw.x -in mn_atom.in      # UPF validation: spherical 3d5 4s2 in a 22-bohr box
mpirun -np 1 ~/Code/q-e/bin/pw.x -in mno.scf.in      # step 1: MnO AFM-II LDA (sla+vwn), k 2x2x2, smearing finds the AFM state
mpirun -np 1 ~/Code/q-e/bin/pw.x -in mno.scf2.in     # step 2: fixed occupations, tot_magnetization=0, from step 1 (hp.x's
                                                     #         2-step recipe for magnetic INSULATORS -- it refuses a gapped
                                                     #         system run as a metal; nbnd must equal step 1's)
mpirun -np 1 ~/Code/q-e/bin/hp.x -in mno.hp.in       # U(Mn 3d), U(O 2p) -> mno.Hubbard_parameters.dat
```
⚠ **nq = 1 is NOT a result**: with q = Γ only the perturbation repeats with the 4-atom cell and hp.x returned
U(Mn 3d) = −0.41 eV, U(O 2p) = 24.4 eV (2026-09-21).  The q-mesh is the supercell size of the linear-response
method; 2×2×2 is the first meaningful mesh (~1 h serial at 100 Ry).
⛔ Every QE binary here is an MPI build: ALWAYS `mpirun -np 1`, never the bare executable (it hangs).  The UPF
files are regenerated, not committed (they are a function of the database + the atom code).

**MnO ground state (pw.x, 2026-09-21):** E = −123.03661 Ry = −61.5183 Ha at k 2×2×2 (CP2K / qchem at Γ on
the VA span: −61.4706 / −61.4112), moments ±4.51 μB, LDA gap 1.32 eV, **atomic-projector 3d occupations
4.988↑ / 0.481↓** against qchem's atomic projector 4.93↑ / 0.45↓ at Γ — the two codes agree on what the
manifold holds.

**UPF validation (2026-09-21):** `mn_atom.in` with the spherical 3d5 4s2 occupations fixed gives
E = −28.48652 Ry = −14.2433 Ha against CP2K's ATOM code −14.2414 and our atom −14.2442 (2 mHa), eigenvalues
3d −6.838 / 4s −5.109 eV against CP2K −7.004 / −5.279: both shifted by +0.17 eV = the periodic box's G=0
alignment of the local potential; the 3d−4s splitting agrees to 4 meV.  With smearing instead of fixed
occupations the isolated atom breaks spherical symmetry (E 64 mHa lower) — compare like with like.
Geometry = `src/Calculation/Data/materials.json` `MnO_AFM2` and `IntegrationTests/CP2K/mno_afm2_gpw_va.inp`.

**hp.x, atomic projector, Mn 3d + O 2p, nq 2×2×2 (2026-09-21, 42 min serial):** U(Mn 3d) = **0.198 eV**, U(O 2p) = 26.56 eV.
The response matrices (`tmp/HP/mno.chi.dat`) say why: χ₀(Mn,Mn) = −0.0726 and χ(Mn,Mn) = −0.0712 — the SCF-screened
response of the Mn d occupation is almost the bare one, so χ₀⁻¹ − χ⁻¹ nearly cancels (a filled majority / empty
minority shell responds weakly; the "closed-shell problem" of linear-response U), while O behaves normally
(χ₀ −0.050 → χ −0.022 → 26 eV).  Non-orthogonalised `atomic` projectors are also the fragile choice in hp.x
(its documentation recommends ortho-atomic for insulators).  The controlled variant `mnoO.*` = **ortho-atomic,
Mn 3d only** (ortho-atomic 3d occupations 4.980↑ / 0.253↓).

**hp.x, ortho-atomic, Mn 3d only, nq 2×2×2, 100 Ry (2026-09-21, 24 min):** U(Mn 3d) = **1.008 eV** (χ₀ −0.0453 →
χ −0.0432: again ~95 % of the bare response survives the screening).  Projector-independent, so not the projector.
⚠ **The cutoff is NOT converged**: E(100 Ry) −123.037, E(140) −123.326, E(200) −123.409 Ry — GTH Mn q7 is hard
(projector radii ~0.2–0.3 bohr) and hp.x had warned "numerical instabilities due to too low cutoff for hard
pseudopotentials".  The occupations are stable (5.231 → 5.231) but a linear RESPONSE need not be; the hp.x
numbers above are therefore not yet the oracle value.  Next: the ground state at a converged cutoff (scan to
360 Ry), then hp.x there.

**Cutoff (2026-09-21):** E(Ry) at 100/140/200/280/360 Ry = −123.037 / −123.326 / −123.409 / −123.422 / −123.424:
**280 Ry** is converged to 1.3 mRy (`mnoO280.*`; the 100 Ry decks stay as the record of the mistake).

**hp.x, ortho-atomic, Mn 3d only, nq 2×2×2, 280 Ry (4 ranks, 49 min):** U(Mn 3d) = **0.958 eV** (χ₀ −0.0447 →
χ −0.0428).  Stable against the cutoff (1.008 at 100 Ry), so the small value belongs to this setup, not to numerics.

**The hp.x BUILD is sane:** QE's own `test-suite/hp_insulator_us_magn` NiO benchmark reproduces on this build,
**7.0521 vs 7.0514 eV** (PBEsol, US-PP, 25 Ry, χ₀ −0.223 → χ −0.086: a strongly screened response, 5× larger
bare response than our MnO's).  `NiOg.*` = the same NiO cell with OUR GTH Ni q10 / O q6 UPFs (LDA sla+vwn,
280 Ry, ortho-atomic): the discriminator between "our UPF route" and "LDA-GTH MnO".

⚠ **EVERY hp.x NUMBER BELOW IS CONDITIONED ON ITS STARTING U** (noticed 2026-09-22, when the ACBN0
comparison needed a matched ground state).  `hp.x` computes the LINEAR RESPONSE of the ground state it is
handed, so a deck's `HUBBARD` block is part of the answer.  The `NiOg.*` decks carry `U Ni-3d 3.0` — they
were derived from QE's own `test-suite/hp_insulator_us_magn` benchmark, which starts there — so the NiO
row is \f$U_{\rm LR}(U_{\rm in}=3\,{\rm eV})\f$, **not** a U=0 response and **not** a self-consistent U
(that needs iterating to \f$U_{\rm out}=U_{\rm in}\f$, which we have not run).  The `mno*.*` decks carry no
`HUBBARD` U, so the MnO rows ARE \f$U_{\rm in}=0\f$.  Quote the condition with the value, and match it
before comparing anything to it.

**NiO with OUR GTH UPFs (`NiOg.*`: Ni q10 + O q6, LDA sla+vwn, 280 Ry, ortho-atomic, k/q 2×2×2, U_in = 3 eV,
42 min on 4 ranks):** U(Ni 3d) = **5.267 eV** (χ₀ −0.113 → χ −0.0706; ortho-atomic 3d occupations 4.980↑ / 3.368↓, gap 2.86 eV) beside
QE's benchmark 7.05 (PBEsol, US-PP, atomic).  ⇒ **the gth2upf route is sound**: our pseudopotential and our
pseudo-atom projector give a normal linear-response U on NiO.  MnO's 0.96 eV is therefore what hp.x's linear
response gives for LDA-GTH MnO: a d⁵ high-spin shell, filled majority / empty minority, whose bare response is
2.5× weaker than NiO's (χ₀ −0.045 vs −0.113) and almost unscreened (χ/χ₀ 0.96 vs 0.62).

| system (this build, LDA, our GTH, 280 Ry, 2×2×2) | projector | U_in | χ₀ | χ | U(3d) |
|---|---|---|---|---|---|
| MnO AFM-II | atomic (+O 2p) | 0 | −0.0726 | −0.0712 | 0.20 eV (100 Ry) |
| MnO AFM-II | ortho-atomic | 0 | −0.0447 | −0.0428 | **0.96 eV** |
| NiO AFM-II | ortho-atomic | **3 eV** | −0.1130 | −0.0706 | **5.27 eV** |
| NiO, QE benchmark (PBEsol, US, 25 Ry) | atomic | **3 eV** | −0.223 | −0.086 | 7.05 eV (reproduced 7.052) |

## A2b (`doc/HubbardUPlan.md`): NiO O-2p, matched-PP LRT — 2026-09-25

The plan's ACBN0-vs-oracle Track A needed a value for U(O 2p) on NiO and Carta et al. (arXiv:2505.03698,
`doc/HubbardUPlan.md` §7) argue LRT stays well-behaved on entangled bands where cRPA does not — so the
right move was extending the ALREADY-TRUSTED matched-PP `NiOg.*` recipe to put O-2p in the Hubbard block
too, not another ABINIT cRPA sweep.  `NiOgO.*` = `NiOg.*` (same cell, same UPFs, same `U_in=3` on Ni)
+ `U O-2p 1.d-8` added to the `HUBBARD {ortho-atomic}` block, prefix changed so the original `NiOg` result
stays intact for comparison.  Same 3-step recipe (`scf.1` finds the AFM state with smearing, `scf.2` fixes
occupations from it, `hp.x` runs the 2×2×2 linear response) — 4 ranks, 280 Ry, ~1h28m wall (`hp.x` alone;
roughly double the Ni-only run's 42 min, sensible with twice the perturbed sites).

**Result:** U(Ni 3d) = **5.4343 eV**, U(O 2p) = **8.5139 eV** (gap 2.86 eV, unchanged from the Ni-only run
— HOMO 10.4685 / LUMO 13.3285 eV).  Ni's value moved only +0.17 eV (+3.2%) from the Ni-only 5.267 eV once
O joined the interacting subspace — a SMALL, sane shift, not the order-of-magnitude swing Carta et al.
report for cRPA under redefinition of the interacting/screening split.  That is itself evidence for their
central claim: **LRT is comparatively insensitive to where the D/R boundary is drawn; cRPA is not.**

⇒ **this is now the oracle for A3's gate-3 refutation test on BOTH orbitals**, replacing the literature's
unstated-convention "≳4 eV" O-2p bound with a real, matched-PP number.  It is also strikingly close to our
own raw ACBN0 estimate for NiO O-2p (7.27 eV, `doc/HubbardUPlan.md` §1 table) — a ratio of ~1.17, far
tighter than Ni-3d's ACBN0-vs-hp.x ratio of ~2.6× — suggesting ACBN0's basis-completeness overshoot
(§6.1's finding) is itself orbital-dependent: a diffuse O-2p likely already reaches near-complete
projector coverage in a way a compact Ni-3d does not.  Decks: `NiOgO.scf.1.in`, `NiOgO.scf.2.in`,
`NiOgO.hp.in`.  As with `NiOg.*`, the UPFs, wavefunctions and `HP/` scratch are regenerated, not committed.

## A6 (`doc/HubbardUPlan.md`): SrVO₃, matched-PP LRT — 2026-09-25

First material off A6's queue (broadening the matched-PP hp.x oracle set beyond MnO/NiO).  SrVO₃ was
picked first because it is the simplest structure in the set — cubic perovskite (no Jahn-Teller distortion
to source), a single-d-electron correlated METAL (no AFM ordering, so none of MnO/NiO's smearing-then-
fixed-occupation 2-step recipe is needed — one plain `scf` step with `occupations='smearing'` is enough,
matching QE's own `test-suite/hp_metal_us_magn` metal recipe, not the insulator one).  Geometry: ABINIT's
own validated `tests/tutoparal/Input/tucalc_crpa_1.abi` (SrVO3 cRPA tutorial) cell, `acell 3*7.2605` bohr,
cubic, V at (0,0,0), Sr at (½,½,½), O at (½,0,0)/(0,½,0)/(0,0,½) — an externally-sourced geometry, not a
guess.

⚠ **A real gth2upf limitation found and routed around, not fixed**: `gth2upf --element Sr --q 10` (the
semicore 4s²4p⁶5s² valence) builds and reports SUCCESS but silently integrates to only 2 electrons, not
10 — `PseudoAtom_EC`'s occupation model fills at most ONE shell per angular momentum `l` (`nv[l]` capped at
`2(2l+1)`), so it cannot represent Sr q10's TWO occupied s-shells (4s and 5s) at once; it silently drops
the second one.  This differs from Mn/Ni q7/q10, where "3d+4s" is two DIFFERENT l channels and the model
has no trouble.  Routed around by using **Sr q2** (5s² only, single s-shell, the light-valence choice —
consistent with this project's existing convention of using the light, non-`default` GTH variant for Mn/Ni
too).  Not fixed because it would need representing multiple radial shells per `l` in the pseudo-atom EC, a
real structural change, not a quick patch — worth a `doc/OpenWork.md` row if a future material needs a
genuine alkaline-earth/alkali semicore potential. V (q5, light valence, `3d³4s²`) and O (q6) built and
integrated correctly on the first attempt.

**Cutoff scan (LDA, `sla+vwn`, k 4×4×4, nonmagnetic, 2026-09-25):** E(Ry) at 80/120/160/200/250/300/350 Ry
= −109.5400 / −110.2918 / −110.4913 / −110.5463 / −110.5638 / −110.5678 / −110.5687 — converged to
**0.94 mRy at 300 Ry** (ecutrho 1200 Ry, dual 4); GTH V-q5 is a hard pseudopotential, same story as Mn q7.

**pw.x ground state (300 Ry, k 4×4×4, `HUBBARD {ortho-atomic}` `U V-3d 1.d-8`, 2m50s serial):**
E = −110.56776433 Ry, Fermi energy 7.6385 eV, converged in 19 iterations.  Projected (ortho-atomic) V-3d
occupation 3.51 electrons (out of 10) — well above the ionic d¹ picture, i.e. the same basis-completeness/
covalency inflation of the atomic-projector occupation already seen on MnO/NiO (§6.1), now on a THIRD
material and a genuinely different (metallic, non-magnetic) electronic structure — evidence the effect is
about the projector, not about magnetism or the insulating gap.

**hp.x, ortho-atomic, V-3d only, nq 2×2×2, U_in≈0 (2026-09-25, 55m21s serial):**
**U(V 3d) = 6.2502 eV.**  χ₀(V,V) = −1.7822 → χ(V,V) = −0.1436 (χ/χ₀ ≈ 0.081) — a MUCH more strongly
screened response than MnO (χ/χ₀ ≈ 0.96) or NiO (χ/χ₀ ≈ 0.62), exactly as expected: SrVO₃ is a real metal
with itinerant carriers available to screen, where MnO/NiO are gapped insulators.  Decks: `srvo3.scf.in`,
`srvo3.hp.in`.  UPFs, wavefunctions and `HP/` scratch regenerated (`gth2upf --element V --q 5`,
`--element Sr --q 2`, `--element O --q 6`), not committed, per this file's convention.

## A6 (`doc/HubbardUPlan.md`): LiCoO₂, matched-PP LRT — 2026-09-26

Second member of A6's LiMO₂ series (started with Co rather than V/Cr — the queue's own item 4 note flags
Cr as untested and likely needing the Cu-style anneal fix, so a well-characterized member was picked first
to shake out the RECIPE, not the element coverage).  Structure: R-3m (#166) layered rock-salt, from a real
single-crystal XRD refinement (Pinsard-Gaudart et al., *J. Crystal Growth* 334 (2011) 165–169, Table 2,
x=1 stoichiometric LiCoO₂) — hexagonal-setting `a = 2.81280(10) Å`, `c = 14.0272(9) Å`, O at Wyckoff 6c
`(0,0,z)` with `z = 1 − 0.7604(3) = 0.2396`.  Converted to the 4-atom rhombohedral PRIMITIVE cell
(QE `ibrav=5`) by hand (R-centering reduction, verified numerically with a small `numpy` script — not
guessed): `celldm(1) = 9.353622` bohr, `celldm(4) = cos α = 0.838532`, Co (½,½,½), Li (0,0,0), O
(0.2396,0.2396,0.2396) and (0.7604,0.7604,0.7604).  **Bond-length check against the reduction**: Co–O
1.9194 Å, Li–O 2.0895 Å — both match the literature CoO₆/LiO₆ octahedral distances (~1.92 Å / ~2.09 Å) to
better than 0.001 Å, confirming the hex→rhombohedral conversion is correct.

⚠ **A real `hp.x` requirement found the hard way**: it refuses to run unless the Hubbard-active atom(s)
are listed FIRST in `ATOMIC_POSITIONS` ("All Hubbard atoms must be listed first...").  The first attempt's
deck listed `Li, Co, O, O`, and `hp.x` stopped immediately with that error.  Reordering to put Co first
(species and positions both) fixed it.  In hindsight this was always
true of the working `NiOg`/`mnoO` decks (`Ni1,Ni2,O` / `Mn1,Mn2,O`) and of `srvo3.scf.in` (`V,Sr,O`) — none
of them happened to need Li/Sr first — but it was never stated as a rule until this deck violated it.
**Rule for every future deck in this series: Hubbard atom(s) first, always.**

**Cutoff scan (LDA, k 4×4×4, nonmagnetic, GTH Co-q9/Li-q1/O-q6):** E(Ry) at 80/120/160/200/250/300/350 =
−119.7478 / −121.0616 / −121.2397 / −121.2783 / −121.2900 / −121.2926 / −121.2933 — converged to 0.62 mRy
at 300 Ry, same working cutoff as SrVO₃.

⚠ **LiCoO₂ is a real band insulator, not a metal** (low-spin Co³⁺, d⁶, t₂g⁶eg⁰ — the paper's own
"spin-unpolarized" choice degenerates to a genuinely gapped nonmagnetic ground state here, unlike SrVO₃'s
metal).  The first `hp.x` attempt on the plain smeared `scf` ground state failed outright: "DOS(E_Fermi) is
too small... most likely the system has a gap, and hence it should NOT be treated as a metal."  Needed the
SAME 2-step recipe as MnO/NiO's magnetic insulators (`scf.1` smeared to find the ground state, `scf.2`
`occupations='fixed'` with `nbnd` pinned to `scf.1`'s count, from file) even though there is no magnetism
here at all — the 2-step recipe is about the GAP, not about AFM order specifically, a distinction this
project's own earlier MnO/NiO write-ups had conflated with "magnetic insulator."  LDA gap (fixed-occupation
run): 1.56 eV.

**hp.x, ortho-atomic, Co-3d only, nq 2×2×2, U_in≈0 (2026-09-26, 1h48m serial):**
**U(Co 3d) = 7.3070 eV.**  χ₀(Co,Co) = −0.3393 → χ(Co,Co) = −0.0971 (χ/χ₀ ≈ 0.286) — screened more than
MnO/NiO (0.96/0.62) but far less than the metallic SrVO₃ (0.081): a real gap, but a smaller one (1.56 eV
vs MnO/NiO's larger LDA gaps) with more covalent/polarizable Co–O bonding to screen with.  Projected
(ortho-atomic) Co-3d occupation 7.41 electrons (out of 10) — the same basis-completeness/covalency
inflation over the ionic d⁶ picture already seen on Mn/Ni/V, now confirmed on a FOURTH material and a
low-spin, fully-paired d-shell, ruling out "unpaired/open-shell character" as the cause.  Decks:
`licoo2.scf.1.in`, `licoo2.scf.2.in`, `licoo2.hp.in`.  UPFs, wavefunctions and `HP/` scratch regenerated
(`gth2upf --element Co --q 9`, `--element Li --q 1`, `--element O --q 6`), not committed.

## A6 (`doc/HubbardUPlan.md`): LiVO₂, matched-PP LRT — 2026-09-26

Second member of the LiMO₂ row.  R-3m structure from a real citation (user-supplied): **Mat. Res. Bull.
27, 555–562 (1992)**, `a = 2.8388(18) Å`, `c = 14.828(13) Å`, V at Wyckoff 3a `(0,0,0)`, Li at 3b
`(0,0,½)`, O at 6c `(0,0,z)` with `z = 0.25749(22)`.  ⚠ **A self-relaxed geometry was tried FIRST and then
discarded** — before this citation surfaced, `pw.x vc-relax` (our GTH-LDA pseudopotentials, QE's own BFGS
optimizer — QE relaxes, not `qchem`; this project's own code has no geometry-relaxation capability at all)
gave `a=2.753 Å, c=14.485 Å` (a sane ~3% LDA-contraction from experiment, bond lengths V–O 1.966 Å /
Li–O 2.026 Å checked sane), but the real citation superseded it the moment it existed — never guess when a
citable number is available, even a defensible self-consistent stand-in.  Converted hexagonal→rhombohedral
primitive by hand as in LiCoO₂; **bond-length check against the citation**: V–O = 1.988 Å, Li–O = 2.121 Å
(both sane for these ionic radii, both correctly larger than LiCoO₂'s 1.919/2.090 Å — V³⁺/Li⁺ vs the
smaller Co³⁺ environment).  `celldm(1) = 9.840415` bohr, `celldm(4) = 0.851403`; V listed first
(Hubbard-atom-first rule, learned the hard way on LiCoO₂).

⚠ **LiVO₂ in this idealized untrimerized cell is a METAL under nonmagnetic LDA** (V³⁺ d², partially-filled
t₂g, no Jahn-Teller/trimer distortion imposed) — confirmed by running `hp.x` directly on the plain smeared
`scf` ground state with no complaint, the SrVO₃ pattern, not the LiCoO₂/MnO/NiO 2-step.  This is the
*idealized* structure the cRPA-comparison paper's "isolated d-manifold LiMO₂" benchmark uses — the REAL
room-temperature LiVO₂ has V-trimer short-range order (Kojima et al., arXiv:1910.01337 and arXiv:2301.03833)
that would gap it; that physics is deliberately out of scope here, matching the reference paper's own
convention, not a mistake.

**hp.x, ortho-atomic, V-3d only, nq 2×2×2, U_in≈0 (2026-09-26, 2h39m serial):**
**U(V 3d) = 5.9526 eV.**  χ₀(V,V) = −2.6592 → χ(V,V) = −0.1530 (χ/χ₀ ≈ 0.0575) — even more strongly
screened than SrVO₃ (0.081), consistent with both being real metals.  Projected (ortho-atomic) V-3d
occupation 3.65 electrons (out of 10) — same covalency-inflation pattern over the formal d² picture seen
on every material so far.  Decks: `livo2.scf.in`, `livo2.hp.in`.  UPFs/wavefunctions/`HP/` regenerated,
not committed.

## A6 (`doc/HubbardUPlan.md`): LiCrO₂, matched-PP LRT — 2026-09-26

Third LiMO₂ member.  **The flagged aufbau risk did NOT materialize**: `gth2upf --element Cr --q 6` converged
on the FIRST attempt to the correct 3d⁵4s¹ ground state (`int rhoatom = 6 (expect 6)`), no kT-anneal
fallback needed — unlike Cu, Cr's near-degenerate 3d⁵4s¹/3d⁴4s² configurations did not trigger a limit
cycle here.  Structure: R-3m, `a = 2.8941(3) Å`, `c = 14.391(3) Å`, Cr at 3a `(0,0,0)`, Li at 3b `(0,0,½)`,
O at 6c `(0,0,z)` with `z = 0.7433(5)` (Garg et al., *Crystals* 9(1), 2 (2019), a real single-crystal-XRD
paper found via web search, browser-fetched since MDPI blocks plain `curl`/`WebFetch`).
**Bond-length check against the paper's own reported averages**: Cr–O = 2.0020 Å (paper: 2.003 Å), Li–O =
2.1144 Å (paper: 2.113 Å) — both agree to <0.1%.  (A units-double-conversion bug in the verification
script itself briefly produced a spurious ~1.06 Å "bond length" during this check — the geometry was right
all along; worth remembering that a bond-length sanity check is only as good as the script computing it,
so re-derive independently when a number looks wrong rather than trusting the first red flag.)
`celldm(1) = 9.599204` bohr, `celldm(4) = 0.837698`; Cr listed first.  Cutoff carried over at 300 Ry from
the established Mn/Ni/V/Co precedent for this PP family, not independently re-scanned.

**hp.x, ortho-atomic, Cr-3d only, nq 2×2×2, U_in≈0 (2026-09-26, 2h38m serial):** ran directly on the
smeared ground state with no complaint — **LiCrO₂ is METALLIC under nonmagnetic LDA** (Cr³⁺ d³, a
half-filled t₂g shell that is only a real Mott insulator once magnetic order or +U opens the gap; forcing
it nonmagnetic here, matching the reference paper's own convention, leaves a partially-filled degenerate
manifold).  **U(Cr 3d) = 5.8111 eV.**  χ₀(Cr,Cr) = −5.5041 → χ(Cr,Cr) = −0.1595 (χ/χ₀ ≈ 0.029) — the
LARGEST bare response of any material so far (three d-electrons available to respond) but also the most
strongly screened, netting a U comparable to V's.  Projected (ortho-atomic) Cr-3d occupation 4.72 electrons
(out of 10) — the covalency-inflation pattern over the formal-ionic count, now confirmed on a FIFTH
material.  Decks: `licro2.scf.in`, `licro2.hp.in`.  UPFs/wavefunctions/`HP/` regenerated, not committed.

## A6 (`doc/HubbardUPlan.md`): LiFeO₂, matched-PP LRT — 2026-09-26/27

Fourth LiMO₂ member.  **Real LiFeO₂ is NOT the layered R-3m phase** — its stable forms are cubic
disordered-rock-salt (α) or orthorhombic (β), since high-spin Fe³⁺ d⁵ has zero crystal-field stabilisation
energy pushing it toward the same cation-ordered layered structure Co³⁺/Ni³⁺/Cr³⁺/V³⁺ favour (user flagged
this from memory before sourcing began — correctly).  **Decision (user, 2026-09-26): idealize it into the
SAME untrimerized R-3m template as the rest of the row anyway**, matching the reference cRPA-comparison
paper's own apparent convention — checked from that paper's own computational-details table, which groups
LiFeO₂ with LiCrO₂/LiCoO₂ under the SAME isotropic `13×13×13` k-mesh (consistent with all three sharing the
same compact rhombohedral-primitive cell shape), unlike LiMnO₂'s oddly-shaped `6×13×6` mesh (consistent with
LiMnO₂ needing its own real, distorted cell in their work).  This is valid for **validating our own U
calculation methodology**, which is A6's actual goal — not a claim about LiFeO₂'s real ground state.

Geometry: Materials Project mp-19419 (GGA+U=5.3 eV relaxed, ICSD-backed: 78712/51759/51207), already given
in the PRIMITIVE rhombohedral form directly (no hex→rh conversion needed): `a = b = c = 5.052 Å`,
`α = β = γ = 33.159°`, Li (0,0,0), Fe (½,½,½), O (0.2404,0.2404,0.2404) and (0.7596,0.7596,0.7596).
Bond-length check: Fe–O = 1.971 Å, Li–O = 2.131 Å — both physically sane for high-spin Fe³⁺/Li⁺ octahedra
(a touch shorter than experiment would likely give, expected for a GGA+U-relaxed source geometry, the same
direction our own `vc-relax` shrank LiVO₂).  `celldm(1) = 9.546896` bohr, `celldm(4) = 0.837156`; Fe listed
first.  Cutoff carried over at 300 Ry.

**hp.x, ortho-atomic, Fe-3d only, nq 2×2×2, U_in≈0 (2026-09-27, 2h54m serial):** ran directly on the smeared
ground state with no complaint — **also METALLIC under nonmagnetic LDA** (high-spin Fe³⁺ d⁵ forced
spin-restricted has no majority/minority split to fill, so it is simply another partially-filled-manifold
metal like V/Cr, not a repeat of MnO's d⁵ "the inverse difference cancels" trap — that trap was specific to
MnO's AFM/spin-polarized response channel, which does not exist in a spin-restricted nspin=1 calculation).
**U(Fe 3d) = 7.5915 eV.**  χ₀(Fe,Fe) = −4.3821 → χ(Fe,Fe) = −0.1242 (χ/χ₀ ≈ 0.028) — as strongly screened as
Cr's.  Projected (ortho-atomic) Fe-3d occupation 6.46 electrons (out of 10) — covalency inflation over the
formal d⁵ picture, now confirmed on a SIXTH material.  Decks: `lifeo2.scf.in`, `lifeo2.hp.in`.
UPFs/wavefunctions/`HP/` regenerated, not committed.

## A6 (`doc/HubbardUPlan.md`): LiNiO₂, matched-PP LRT — 2026-09-27 — **LiMO₂ ROW COMPLETE**

Fifth and last LiMO₂ member.  Real LiNiO₂ is also Jahn-Teller-active (low-spin Ni³⁺ d⁷, singly-occupied
eg — confirmed directly: Materials Project's own DFT+U-relaxed entry (mp-25411) comes back with UNEQUAL
rhombohedral angles and two different oxygen z-components, i.e. genuinely monoclinic-distorted, not the
idealized R-3m this row uses).  Same decision as Fe: idealize into the row's untrimerized R-3m template.

Geometry from a real Rietveld refinement, Seo et al., *J. Electrochem. Soc.* 165 (2018) A2554
("Updating the Structure and Electrochemistry of Li$_x$NiO₂"): `a = 2.8751(1) Å`, `c = 14.2000(3) Å`,
Li 3a (0,0,0), Ni 3b (0,0,½), O 6c (0,0,z) with `z = 0.2424(1)`.  Bond-length check: Ni–O = 1.978 Å
(paper's own cation-mixing-corrected refinement reports ~1.95 Å average — same ballpark, the small gap is
the 1.81% Li/Ni antisite mixing this idealized 0%-mixing cell doesn't carry), Li–O = 2.103 Å (matches ~2.10
Å almost exactly).  `celldm(1) = 9.478789` bohr, `celldm(4) = 0.835726`; Ni listed first (reuses the same
Ni q10 pseudopotential already validated for rocksalt NiO — note LiNiO₂'s Ni is 3+/d⁷, a different ion
entirely from NiO's 2+/d⁸, so their U values are not expected to relate simply).

**hp.x, ortho-atomic, Ni-3d only, nq 2×2×2, U_in≈0 (2026-09-27, 2h48m serial):** ran directly on the smeared
ground state, no complaint — **METALLIC under nonmagnetic LDA** (the undistorted cell leaves eg¹ split
across two degenerate orbitals, a partially-filled manifold like every other member of this row bar Co).
**U(Ni 3d) = 9.1730 eV** — the largest of the whole LiMO₂ row.  χ₀(Ni,Ni) = −2.0136 → χ(Ni,Ni) = −0.1006
(χ/χ₀ ≈ 0.050), in the same strongly-screened-metal range as V/Cr/Fe.  Projected (ortho-atomic) Ni-3d
occupation 8.32 electrons (out of 10) — covalency inflation over the formal d⁷ picture, now confirmed on
EVERY material run this session (seven for seven).  Decks: `linio2.scf.in`, `linio2.hp.in`.

★ **First use of the new checkpoint practice** (`doc/OpenWork.md`'s SCF-restart feature row, 2026-09-27):
the converged U_in≈0 ground state's `.save` (72 MB — wavefunctions, ρ(G), `occup.txt`'s Hubbard occupation
matrix, full input/output metadata) was archived to `IntegrationTests/QE/checkpoints/linio2_U0.save/`
BEFORE `hp.x` touched it, rather than deleted after the run like every earlier material this session.  A
future self-consistent-U rerun on LiNiO₂ can `startingpot/startingwfc='file'` from this instead of
reconverging from an atomic guess.  (The five materials done earlier this session — SrVO₃, LiCoO₂, LiVO₂,
LiCrO₂, LiFeO₂ — do NOT have this checkpoint; their U=0 states were deleted before this practice started,
so their self-consistent-U reruns will need a fresh SCF, a modest one-time cost of a few minutes each.)

**LiMO₂ row summary (M = V, Cr, Fe, Co, Ni; all matched-PP `hp.x`, U_in ≈ 0, nspin=1):**

| M | U(M 3d) eV | χ/χ₀ | electronic character |
|---|---|---|---|
| V | 6.2502 (SrVO₃, perovskite) / 5.9526 (LiVO₂) | 0.081 / 0.058 | metal (both hosts) |
| Cr | 5.8111 | 0.029 | metal |
| Fe | 7.5915 | 0.028 | metal |
| Co | 7.3070 | 0.286 | insulator (LDA gap 1.56 eV) |
| Ni | 9.1730 | 0.050 | metal |

Co is the odd one out electronically (real gapped low-spin d⁶) as well as structurally (the only member
whose real ground state IS this same R-3m cell — V, Cr and Ni's real ground states are metallic/JT-active
in ways this idealized treatment does not capture, and Fe's real ground state isn't this phase at all).

## A6 (`doc/HubbardUPlan.md`): TiO₂ (rutile), matched-PP LRT — 2026-09-27

First of the ACBN0-paper's remaining benchmark set (TiO₂/ZnO/FeS₂).  A genuinely different chemistry from
the whole LiMO₂ row: Ti here is formally **Ti⁴⁺, d⁰** — no partially-occupied shell at all — so this tests
whether the +U machinery gives a sane number for a NOMINALLY empty correlated orbital (the effect DFT+U has
on a d⁰ semiconductor's conduction band, a real and common benchmark case, not a mistake).

Geometry: rutile, tetragonal `P4₂/mnm` (#136), from **QE's own `PP/examples/example08`** — a real,
already-built deck authored by Iurii Timrov himself (one of the `hp.x` method's own authors) —
`a = b = 4.5941 Å`, `c = 2.9589 Å`, Ti (0,0,0)/(½,½,½), O (0.3057,0.3057,0) and symmetric partners.
Independently cross-checked against CP2K's own `tests/QS/regtest-sym-2/c_17_rutile.inp` (a symmetry-only
regtest, but its cell is cited to Wyckoff's *Crystal Structures* Vol. I pp. 250–2): `a=b=4.59373 Å`,
`c=2.95812 Å` — agrees with the QE source to <0.01%.  QE `ibrav=6` (tetragonal P) used directly:
`celldm(1) = 8.681591` bohr, `celldm(3) = 0.644065`.  Ti listed first.  GTH Ti q4 (light valence, 3d²4s²,
matching the row's convention) converts cleanly.

**pw.x ground state (300 Ry, k 4×4×4, `occupations='fixed'`, 41 iterations):** E = −143.26265558 Ry.
No 2-step recipe needed — unlike every partially-filled-shell material this session, a genuine d⁰ gapped
semiconductor has no occupation ambiguity for `occupations='fixed'` to resolve.

**hp.x, ortho-atomic, Ti-3d only, nq 2×2×2, U_in≈0 (2026-09-27, 2h32m serial):**
**U(Ti 3d) = 4.6368 eV** (both symmetry-equivalent Ti sites agree exactly, as expected).  χ₀(Ti,Ti) =
−0.4313 → χ(Ti,Ti) = −0.1438 (χ/χ₀ ≈ 0.333).  Projected (ortho-atomic) Ti-3d occupation 4.665 electrons
(out of 10) — even a FORMALLY EMPTY d-shell picks up substantial weight from O-2p→Ti-3d covalent mixing
under the atomic projector, the same basis-completeness/covalency inflation seen on every other material
this session, now shown to have NOTHING to do with how many d-electrons are formally present.  Checkpoint
archived: `checkpoints/tio2_U0.save/`.  Decks: `tio2.scf.in`, `tio2.hp.in`.

## A6 (`doc/HubbardUPlan.md`): ZnO (wurtzite), matched-PP LRT — 2026-09-27 — **NOT A USABLE ORACLE**

Second of the remaining ACBN0-paper set.  Zn here is **Zn²⁺, d¹⁰ — a genuinely CLOSED shell**, not merely
spin-cancelled the way MnO's high-spin d⁵ is.  Geometry: wurtzite, hexagonal `P6₃mc` (#186), the standard
literature cell (multiple independent sources converge tightly): `a = 3.2495 Å`, `c = 5.2069 Å`, Zn at
(⅓,⅔,0)/(⅔,⅓,½), O at (⅓,⅔,u)/(⅔,⅓,u+½) with `u = 0.3825`.  Bond-length check: Zn–O = 1.973 Å (literature
~1.973–1.99 Å).  QE `ibrav=4` (hexagonal): `celldm(1) = 6.140665` bohr, `celldm(3) = 1.602370`.  Zn listed
first.  GTH Zn q12 (full 3d¹⁰4s² in the valence, needed since the whole point is putting +U ON the d-shell)
converts cleanly.  `pw.x` ground state (300 Ry, `occupations='fixed'`, no ambiguity for a real closed-shell
gapped insulator): E = −306.18412719 Ry, 18 iterations.

**hp.x, ortho-atomic, Zn-3d only, nq 2×2×2, U_in≈0 (2026-09-27, 50m44s serial):** ran without complaint and
returned **U(Zn 3d) = 35.3358 eV — an order of magnitude larger than every other material this session**,
and it is **NOT a trustworthy number**.  χ₀(Zn,Zn) = −0.00347 → χ(Zn,Zn) = −0.00309 (χ/χ₀ ≈ 0.89): both are
TINY compared to every partially-filled-shell material (χ₀ ranged −0.34 to −5.50 for Co/V/Cr/Fe/Ni) and
nearly equal to each other.  This is **exactly the "closed-shell problem of linear-response U" already named
in this file for MnO's d⁵ AFM case** (§ MnO entry above: "a filled majority/empty minority shell responds
weakly... χ₀ and χ nearly equal, so χ₀⁻¹−χ⁻¹ nearly cancels") — but reached by a DIFFERENT mechanism this
time: MnO's version comes from spin-cancellation (majority filled, minority empty, in a spin-polarized AFM
calculation); ZnO's comes from genuine electron-shell closure (d¹⁰, no available states for the Hubbard
perturbation to shift into AT ALL, so the bare response χ₀ itself is tiny, not merely screened away).  Same
mathematics — inverting the difference of two nearly-equal small numbers — two different physical routes to
it.  **Do not use 35.34 eV as an oracle value**; flag ZnO/Zn-3d alongside MnO/Mn-3d as a material where
same-site LRT U is close to ill-posed, not merely "large".  Checkpoint archived anyway (`checkpoints/
zno_U0.save/`, 104 MB) since the ground state itself is perfectly good — only the U extraction is the
problem.  Decks: `zno.scf.in`, `zno.hp.in`.

## A6 (`doc/HubbardUPlan.md`): KCuF₃, matched-PP LRT — 2026-09-28

**Dr. Carta replied with his actual input files** (`~/Code/reprints/materials_cloud_submission/`, KCuF₃/
Sr₂FeO₄/CrO₂/NiO all included) — real deck, no more guessing.  `KCuF3/LRT/dp/kcuf.1.scf.in`: cubic
`Pm-3̄m`, `a = 4.066704097 Å` (exact, from their relaxation), K (0,0,0), Cu (½,½,½), F (0,½,½)/(½,0,½)/
(½,½,0) — the SAME simple-cubic-perovskite template as SrVO₃, just K/Cu/F in place of Sr/V/O.  ⚠ **Their
own deck reveals real values that correct the SI's stated convention**: `ecutrho/ecutwfc = 672/84 = 8`, not
the SI's stated "four times" (Sr₂FeO₄'s own deck similarly runs a 10× dual, see below) — always read the
actual input file, a methods-section summary rounds off exactly the kind of detail that matters here.
Their `U_projection_type='ortho-atomic'` is the GROUND-STATE +U potential's projector (matches ours); their
actual Hubbard-PARAMETER determination is MLWF-based (Wannier90) per the paper's own stated method, a
different projector convention from our `hp.x` route — so their number and ours are not directly
comparable even now that the geometry is exact (same caveat flagged before the reply arrived).  We use
OUR OWN established recipe (GTH via `gth2upf`, `hp.x` ortho-atomic, U_in≈0) on their exact cell.

⚠ **GTH Cu-q11 needed the full cutoff-scan treatment, unlike the rest of this session's PPs.**  300 Ry (this
session's default carry-over) was NOT converged: E(Ry) at 300/350/400/450 = −242.14275/−242.15196/
−242.15491/−242.15588 — 9.2 mRy between 300→350, only settling to ~1 mRy at 450→400.  **Production cutoff:
450 Ry** (`ecutrho` 1800).  Cu q11's aufbau fallback fired exactly as recorded in `doc/HubbardUPlan.md`
(`gth2upf`'s kT-anneal-then-cold-MOM path, E_atom = −47.926315 Ha, bit-for-bit the same number logged when
the fallback was first built) — reproducible, not a fluke.  K used q1 (light valence, matching the
established Li/Na alkali-metal convention: not the semicore q9).

**pw.x ground state (450 Ry, k 4×4×4, `occupations='smearing'`, 12 iterations):** E = −242.15587589 Ry.
**hp.x, ortho-atomic, Cu-3d only, nq 2×2×2, U_in≈0 (2026-09-28, 2h20m serial):** ran directly on the smeared
state, no complaint — **METALLIC under nonmagnetic LDA** (Cu²⁺ d⁹ in the undistorted cubic cell has no
Jahn-Teller gap, a partially-filled eg shell like every other undistorted-cell member of this session).
**U(Cu 3d) = 8.1629 eV.**  χ₀(Cu,Cu) = −0.7349 → χ(Cu,Cu) = −0.0970 (χ/χ₀ ≈ 0.132).  Projected occupation
9.326 electrons (out of 10) — the MILDEST covalency inflation of any material this session (d⁹ has only
one hole's worth of headroom to inflate into, unlike d⁵–d⁸ elsewhere).  Checkpoint archived:
`checkpoints/kcuf3_U0.save/`.  Decks: `kcuf3.scf.in`, `kcuf3.hp.in`.

## A6 (`doc/HubbardUPlan.md`): Sr₂FeO₄, matched-PP LRT — 2026-09-28/29

Second of Dr. Carta's real decks.  `Sr2FeO4/LRT/dp/sfo.1.scf.in`: body-centered tetragonal `I4/mmm`
(K₂NiF₄-type), 7-atom primitive cell (Fe₁Sr₂O₄, one formula unit), given as EXACT Cartesian
`CELL_PARAMETERS {angstrom}`.  Fe (0,0,0); O at (0.842024,0.842024,0)/(0.157976,0.157976,0) (equatorial)
and (0.5,0,0.5)/(0,0.5,0.5) (apical); Sr at (0.642827,0.642827,0)/(0.357173,0.357173,0).  Bond-length
sanity check: Fe–O equatorial ×4 = 1.9547 Å, apical ×2 = 1.9835 Å (both sane for an FeO₆ octahedron, close
to regular); Sr–O 9-fold coordination, 2.50–2.77 Å (sane for Sr²⁺).  Fe is formally Fe⁴⁺ (d⁴) here — a
different oxidation state entirely from LiFeO₂'s Fe³⁺ (d⁵), a genuine second data point on the SAME
element.  Same caveat as KCuF₃: their own Hubbard-parameter route is MLWF-based, ours is `hp.x`
ortho-atomic — matched geometry, still not matched projector.

⚠ **`hp.x` crashed outright on Carta's raw relaxed cell**: `Error in routine d_matrix (9): D_S (l=2) for
this symmetry operation is not orthogonal`.  Diagnosis: their three cell vectors are equal in magnitude to
only ~8 significant figures (a real DFT relaxation's numerical tolerance, not a transcription error), and
QE's symmetry-finder detects a symmetry operation consistent with the IDEALIZED tetragonal lattice that
the actual (very slightly asymmetric) vectors don't EXACTLY satisfy — building the l=2 (d-orbital) Wigner
rotation matrix for that operation then fails an exact-orthogonality check.  **Fix: symmetrize the cell**
— average the three vectors' common |x|,|y| component and z component (agreement to 8 figures either way,
so this changes nothing physical) and rebuild `CELL_PARAMETERS` from the exactly-symmetric values, same
atomic fractional coordinates unchanged.  Re-ran: identical total energy to the last printed digit,
confirming the fix only removed noise.  **General lesson for any future real-DFT-relaxed cell fed to
`hp.x`**: expect this, and fix it by symmetrizing the geometry before debugging anything else — QE's own
symmetry-detection is more exacting about EXACT invariance than a relaxation's convergence threshold
guarantees.

**pw.x ground state (300 Ry, k 4×4×4, `occupations='smearing'`, 23 iterations):** E = −172.96422568 Ry.
**hp.x, ortho-atomic, Fe-3d only, nq 2×2×2, U_in≈0 (2026-09-28/29, 11h4m serial — the longest run this
session by a wide margin, the larger 7-atom lower-symmetry cell costing far more per q-point):** ran
without further complaint once the cell was symmetrized — METALLIC under nonmagnetic LDA.
**U(Fe 3d) = 8.0116 eV** — close to LiFeO₂'s 7.5915 eV despite the different oxidation state and host,
a sane cross-check.  χ₀(Fe,Fe) = −4.1515 → χ(Fe,Fe) = −0.1194 (χ/χ₀ ≈ 0.029), the same strongly-screened
range as Cr/Fe/Ni.  Projected occupation 6.32 electrons (out of 10) against the formal d⁴ picture — the
largest covalency inflation of any material this session, consistent with Sr₂FeO₄'s known negative-charge-
transfer character (Fe⁴⁺ is a strong enough oxidant that real electronic structure carries substantial
O-2p hole/ligand character — exactly the "entangled d-p" case Carta et al.'s own paper flags this material
as).  Checkpoint archived: `checkpoints/sr2feo4_U0.save/`.  Decks: `sr2feo4.scf.in`, `sr2feo4.hp.in`.

**Both of Dr. Carta's materials are now done — KCuF₃/Sr₂FeO₄ queue items closed.**
