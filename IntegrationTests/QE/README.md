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
