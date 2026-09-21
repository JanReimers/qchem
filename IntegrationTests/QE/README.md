# Quantum ESPRESSO decks — the DFT+U value oracle (`hp.x`)

`hp.x` (Timrov, Marzari, Cococcioni: DFPT linear-response Hubbard U) is the ORACLE for U values (doc/OpenWork.md
step 5, item 4).  An oracle is only as good as its match, so these decks run **our own pseudopotential and our own
Hubbard projector**: the UPF files are written by `CLIapps/gth2upf` from the GTH parameter database with
`PP_CHI` = the qchem GTH pseudo-atom's orbitals (the `HubbardU_Atomic` radial) — QE's `atomic` projector IS our χ.

```
build/Release/CLIapps/gth2upf --element Mn --q 7     # -> Mn.pz-gth-q7.UPF   (4s, 3d; E_atom -14.2442 Ha)
build/Release/CLIapps/gth2upf --element O  --q 6     # -> O.pz-gth-q6.UPF    (2s, 2p; E_atom -15.7484 Ha)
mpirun -np 1 ~/Code/q-e/bin/pw.x -in mn_atom.in      # UPF validation: spherical 3d5 4s2 in a 22-bohr box
mpirun -np 1 ~/Code/q-e/bin/pw.x -in mno.scf.in      # MnO AFM-II LDA (sla+vwn), k 2x2x2, U=1e-8 declares the manifolds
mpirun -np 1 ~/Code/q-e/bin/hp.x -in mno.hp.in       # U(Mn 3d), U(O 2p)
```
⛔ Every QE binary here is an MPI build: ALWAYS `mpirun -np 1`, never the bare executable (it hangs).  The UPF
files are regenerated, not committed (they are a function of the database + the atom code).

**UPF validation (2026-09-21):** `mn_atom.in` with the spherical 3d5 4s2 occupations fixed gives
E = −28.48652 Ry = −14.2433 Ha against CP2K's ATOM code −14.2414 and our atom −14.2442 (2 mHa), eigenvalues
3d −6.838 / 4s −5.109 eV against CP2K −7.004 / −5.279: both shifted by +0.17 eV = the periodic box's G=0
alignment of the local potential; the 3d−4s splitting agrees to 4 meV.  With smearing instead of fixed
occupations the isolated atom breaks spherical symmetry (E 64 mHa lower) — compare like with like.
Geometry = `src/Calculation/Data/materials.json` `MnO_AFM2` and `IntegrationTests/CP2K/mno_afm2_gpw_va.inp`.
