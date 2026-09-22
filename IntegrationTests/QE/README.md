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
