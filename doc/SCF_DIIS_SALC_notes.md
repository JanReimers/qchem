# DIIS with SALC (symmetry-adapted) bases — notes & open issue

Context: a `SymmetryAdaptedBasisSet` presents one `IrrepBasisSet` per point-group irrep, so the
SCF runs the per-irrep machinery (one `IrrepWF` + one `SCFIrrepAccelerator` per irrep) that was
built for atoms. Three things had to change to make DIIS behave for molecules.

> **Resolution (2026-06-19).** Fixing only issue 1 (the empty-irrep DIIS bail) turned out to be
> enough: with DIIS engaging cleanly, symmetric H₂O HF converges *identically* to the non-symmetric
> reference (DIIS at iter 7, converged ~iter 10) — and issues 2 and 3 **did not reproduce**. The
> occupation never flips and there is no crater. Best reading: on the current code (after the
> DBCache re-entrancy + geometry-key fixes) issue 2 was a *consequence* of issue 3's crater — the
> aufbau flipped on the garbage extrapolated Fock — and issue 3 itself no longer occurs. MOM (the
> overlap-based occupation tracker) was prototyped as the principled cure for issue 2; it is in the
> tree but **parked behind `constexpr bool EnableMOM=false`** (CompositeWF) since it is not
> load-bearing for these closed-shell cases. Keep it for excited-state / ΔSCF / hard open-shell.
> Per-section status is updated below.

## 1. DIIS never engaged — empty irreps blocked it  (FIXED)

**Fixed** with a process-level *seeded* discriminator in `SCFAcceleratorDIIS::CalculateProjections`:
once any irrep has shown a nonzero error we are past the zero-density start (`itsSeeded=true`), so
from then on an *exactly*-zero error means a **permanently-empty** irrep → skip it (it contributes
nothing) instead of bailing the whole extrapolation. Before seeding, a zero error still bails (the
not-yet-seeded first iterations atoms rely on). Atoms are doubly safe: they seed on iteration 1 and
have no empty DIIS channels anyway (`Atom_EC::GetN` gives empties a `Null` accelerator). The ~18
atomic tests the earlier naive "skip zero-error" attempt broke all pass.

### Original analysis (kept for context)

`SCFAcceleratorDIIS::CalculateProjections` sums the per-irrep commutator error and **bails the
whole extrapolation** the moment any irrep reports `GetError()==0`. For a molecule that includes
**unoccupied irreps** (e.g. A₂ for H₂O's ground state): `D'=0` there, so `[F',D']≡0` every
iteration → DIIS never runs (`Bail=Enk==0.0`), only plain linear mixing.

Why atoms are immune: `Atom_EC::GetN` returns the *per-channel* occupation, so empty channels get
`occ=0` → a `Null` accelerator that never reaches the loop. `Molecule_EC::GetN` returns the
**total** (the per-irrep split is an aufbau *output*, unknown when the accelerators are created),
so every irrep — empty ones included — gets a real DIIS accelerator.

Attempted fix (REVERTED): skip zero-error irreps instead of bailing.  **This broke ~18 atomic
HF/DHF/DFT tests** (e.g. Li HF off by 3.4%): the `Enk==0.0` bail *also* guards the zero-initial-
density first iterations, where every occupied channel legitimately has `[F',D']==0` until it has a
density.  Skipping it let DIIS extrapolate from a zero/garbage seed → bad early step → wrong
convergence.  So the bail is restored and DIIS again does NOT engage for symmetric molecules.

A correct fix must distinguish a **permanently-empty** irrep (skip it) from a **not-yet-seeded /
transiently-zero** occupied channel (still bail) — e.g. give the empty irrep a real `Null`
accelerator once the aufbau occupation is known (after iteration 1), rather than a per-iteration
`GetError()==0` heuristic.  And it must still solve issue 3 (the overshoot) to be worthwhile.

## 2. Occupation flipped under acceleration — global aufbau on the extrapolated Fock  (DID NOT REPRODUCE)

Once issue 1 was fixed, the documented flip (`A1=6,B1=2,B2=2` → `A1=8,B1=2,B2=0`) **did not occur**:
the Configuration trace shows `(1a₁)²(2a₁)²(1b₁)²(3a₁)²(1b₂)²` stable from iteration 2 through DIIS
engagement and convergence. Verified by running with MOM disabled *and* the old sticky-freeze
removed — still no flip. Conclusion: the flip was driven by issue 3's crater (a garbage extrapolated
Fock), not an independent defect; with no crater the aufbau stays put. The sticky-freeze was
therefore removed in favour of (parked) MOM. The per-iteration **Configuration** column in the SCF
trace — e.g. `(1a₁)²(2a₁)²(1b₁)²(3a₁)²(1b₂)²` — makes any occupation churn visible.

### Original analysis (kept for context)
With the molecular aufbau the per-irrep occupation is recomputed every iteration from the current
eigenvalues. Once DIIS engages, those eigenvalues come from the **DIIS-extrapolated** (non-physical)
Fock, and for H₂O the close B₂/A₁ valence levels reorder → the aufbau flips `A1=6,B1=2,B2=2` to
`A1=8,B1=2,B2=0` → wrong state.  (Original fix: a sticky per-irrep occupation freeze; superseded.)

## 3. The shared-coefficient multi-irrep DIIS step overshoots  (DID NOT REPRODUCE on current code)

With issue 1 fixed, symmetric H₂O HF no longer craters: DIIS drives Δ[F,D] down monotonically
(iters 7–10: 4.0e-2 → 1.4e-3 → 1.8e-4 → 2.0e-5) exactly like the non-symmetric reference — the
−70.17 crater and the ~57-iteration recovery are gone. The original analysis below may still apply
to harder/larger systems (genuinely disparate per-irrep convergence rates), in which case the
shared-coefficient safeguard is the right tool; for now it is not needed. **Watch for it** when the
single shared-coefficient extrapolation meets very stiff multi-block problems.

### Original analysis (kept for context)
Even with (1) and (2) fixed, the **first** multi-irrep DIIS step deterministically craters
(H₂O/DZVP: total energy → −70.1699437, the same value regardless of `EMax`), then recovers over
~12 iterations — twice. Net: ~57 iterations vs 26 with no DIIS, i.e. DIIS currently makes the
**symmetric** case *worse*.

It is **not** the occupation (frozen correct through the crater) and **not** the empty irrep
(frozen at 0, orbitals unused). Root cause: DIIS solves **one** coefficient set from the *summed* B
matrix and applies it to **all** irreps. When irreps converge at very different rates — A₁ (O 1s
core, fast) vs B₂ (valence, slow) — the shared coefficients over-extrapolate the slow block. The
single-IBS (Unit-symmetry) run never sees this because there is only one block, so its 2-point DIIS
is well behaved.

Note this is the *standard* full-error-vector DIIS (B = Σ_irrep ⟨eᵢ|eⱼ⟩, one coefficient set), so
the issue is the **lack of a robustness safeguard on the early step**, not the formulation per se.
Candidate fixes (SCF-convergence policy — TBD):
- reject / damp a DIIS step that raises the total energy (trust-region style), or
- cap the DIIS coefficient norm and fall back to the plain step when it is exceeded, or
- restart (purge) the DIIS history on a detected energy rise so it doesn't re-overshoot.

As of the resolution above, DIIS *engages* for SALC bases and converges as fast as the
non-symmetric run; the overshoot safeguard is held in reserve, not currently needed.

# Some notes from my research

Research into controlling electron configuration oscillations during Self-Consistent Field (SCF) iterations highlights several key algorithmic interventions. 

Orbital Mixing and Tracking The Method of Occupied-Moving (MOM) is a primary technique for oscillating systems, where occupied orbitals are selected by "shape" similarity to a previous stable iteration rather than by eigenvalue, preventing unwanted orbital swapping.  This is often combined with DIIS (Direct Inversion in the Iterative Subspace) but should be activated only once oscillatory behavior is detected to ensure stability. 

Damping and Level Shifting To dampen large fluctuations in early iterations, algorithms employ Fock matrix damping (static mixing of adjacent Fock matrices) or keyword-based damping protocols like SlowConv and VerySlowConv in ORCA, which increase mixing parameters to control charge sloshing.  Additionally, Level Shifting artificially raises the energy of virtual orbitals to widen the HOMO-LUMO gap, thereby stabilizing the iteration process for systems with small energy gaps or metallic character. 

Advanced Convergence Algorithms For pathological cases where standard methods fail, Second-Order SCF (SOSCF) provides quadratic convergence near the solution, while the Trust Region Augmented Hessian (TRAH) algorithm acts as a robust second-order fallback that automatically activates when DIIS struggles.  These methods ensure convergence to a true local minimum on the orbital rotation surface, even for open-shell transition metals and complex metal clusters.

## Maximum Overlap Method (MOM)
The foundational paper for preventing orbital swapping and maintaining specific electronic configurations is:

Gilbert, A. T. B.; Besley, N. A.; Gill, P. M. W. "Self-Consistent Field Calculations of Excited States Using the Maximum Overlap Method (MOM)." J. Phys.  Chem. A 2008, 112 (50), 13164–13171. DOI: 10.1021/jp801738f.
This work introduces the algorithm of maximizing overlap between occupied orbitals of successive iterations to avoid variational collapse and oscillation.
Barca, G. M. J.; Gilbert, A. T. B.; Gill, P. M. W. "Simple Algorithms for Excited-State SCF Calculations." J. Chem.  Theory Comput. 2018, 14, 1501–1509.  DOI: 10.1021/acs.jctc.7b01234.
Discusses the Initial Maximum Overlap Method (IMOM), which maximizes overlap with the initial guess rather than the previous cycle, offering improved stability for difficult cases. 

## Trust-Region Augmented Hessian (TRAH)
For robust convergence in open-shell and antiferromagnetic systems where DIIS fails:

Helmich-Paris, B. "A trust-region augmented Hessian implementation for restricted and unrestricted Hartree–Fock and Kohn–Sham methods." J. Chem.  Phys. 2021, 154, 164104. DOI: 10.1063/5.0040798.
Demonstrates that TRAH-SCF achieves convergence for notoriously difficult open-shell molecules where standard DIIS diverges, often finding lower-energy symmetry-broken solutions. 

## Damping, Level Shifting, and DIIS Variants
Key references for stabilization techniques involving mixing and energy gaps:

Pulay, P. "Convergence acceleration of iterative sequences.  The case of SCF iteration." Chem. Phys. Lett. 1980, 73 (2), 393–398.  DOI: 10.1016/0009-2614(80)80396-4.
The original formulation of DIIS, the standard acceleration method which can induce oscillations if applied too early. 
Kudin, K. N.; Scuseria, G. E.; Cancès, E. "A black-box self-consistent field convergence algorithm: One step closer." J. Chem.  Phys. 2002, 116, 8255. DOI: 10.1063/1.1470195.
Introduces EDIIS (Energy-DIIS) and discusses combining it with CDIIS to handle oscillatory behavior and ensure global convergence. 
Rabuck, D. D.; Scuseria, G. E. "Improving self-consistent field convergence by varying occupation numbers." J. Chem.  Phys. 1999, 110, 695. DOI: 10.1063/1.478176.
Discusses the use of Fermi-Dirac smearing (fractional occupations) and dynamic damping to stabilize convergence in systems with small HOMO-LUMO gaps. 

## General Convergence Reviews
Cancès, E.; Mennucci, B.; Tomasi, J. "A new integral equation formalism for the polarizable continuum model: Theoretical background and applications to isotropic and anisotropic dielectrics." J. Chem. Phys. 1997, 107, 3032. (Note: While focused on PCM, Cancès has extensive reviews on SCF math).
Klamroth, T. "On the convergence of the SCF iteration for the Kohn-Sham equations." ESAIM: M2AN 2007, 41 (2), 307–322.  DOI: 10.1051/m2an:2007016.
Provides a mathematical analysis of discontinuous changes in occupation numbers and how broadening techniques stabilize the iteration. 
