# Gaussian and Plane Waves (GPW) — Design Spec

Gaussian orbitals + plane-wave density/Coulomb. Strategic placement and first design. This is the
**convergence point** of the molecular Gaussian machinery and the PW Coulomb machinery — not a side
method. Depends on molecular pseudopotentials ([[doc/MolecularPseudopotentialPlan]]). Written for a
fresh session; altitude is "strategic + first milestone," because GPW is bigger and further off than
the PP work.

## 1. What GPW is, and why it's the unifier

Orbitals stay compact **Gaussians** (reuse the molecular integral engine). The density
`ρ = Σ_ab D_ab χ_a χ_b` is **collocated on a uniform grid**, FFT'd, and the Hartree potential solved
in G-space (`V_H(G) = 4π ρ(G)/G²`), then integrated back to `⟨χ_a|V_H|χ_b⟩`. The FFT Poisson solve is
O(N log N) — it sidesteps the 4-centre ERI entirely (this is how CP2K reaches thousands of atoms).

**The unification claim:** with GPW a molecule and a solid become the *same calculation* — identical
Gaussian orbital basis, identical integral engine — differing only in (a) the **Coulomb strategy** and
(b) boundary conditions. You don't build a separate orbital world for solids (as pure plane waves do);
you reuse the molecular one and swap the Coulomb term. That is the cleanest expression of the
"three axes" goal ({atom/molecule/Lattice} × {1E,DFT,HF} × {Pol,unPol}).

## 2. Hard dependency: molecular pseudopotentials

Plain GPW is clean **only when the density is smooth** enough to collocate on a modest grid.
All-electron Gaussian cores are too sharp → you would need **GAPW augmentation** (PAW-like, hard).
**Pseudopotentials give smooth valence density → plain GPW, no augmentation.** So molecular PPs are
GPW's enabler — do them first. GAPW (augmentation for hard cores / all-electron) is explicitly **out
of scope** for the first pass; revisit only if all-electron GPW is ever wanted.

## 3. Architectural placement: a Coulomb/Hartree *strategy* at the term level

GPW is **orthogonal to the orbital basis** — it's a third Hartree strategy beside the two you already
have:

| strategy | mechanism | home |
|---|---|---|
| exact 4-centre | `⟨ab|cd⟩` ERIs | HF (`Vee`/`Vxc`) |
| density fitting | fit ρ in aux Gaussians, 3-centre | current DFT (`FittedVee`) |
| **GPW** | collocate ρ → FFT → Poisson → integrate back | **new** |

Same `⟨χ|V_H|χ⟩` out; different internals. It slots beside the others — it does **not** restructure
the orbital basis or the SCF. (Cf. "the general matrix is the open design target" — Coulomb is a
swappable strategy.)

## 4. AO/FT clarification (important): GPW is all-`double`, FFT is hidden

GPW is **entirely the `double` world** — real orbitals, real density. The FFT is an *internal
Poisson-solve technique of the Hartree term*, **not** a representation of complex Bloch states. So GPW
is **not** "double orbitals + dcmplx density"; it's all-`double`, with a `FourierMap` living
**privately inside** the GPW Hartree term and never surfacing in any interface (the term exposes a
`double` matrix). GPW is therefore the **strongest example yet of "FourierMap is an implementation
detail"** — it reinforces the [[project_ao_ft_projection_cleanup]] conclusion rather than challenging
it. The collocation grid + FFT are term-private; nothing leaks.

## 5. The k-point seam (does the unification survive periodicity?)

- **Γ-point / molecule:** orbitals are real → the unification is exact and **all-`double`**. Molecule
  and Γ-solid differ only in the Coulomb term (+ periodic images in the Poisson solve). This is the
  clean first milestone.
- **General k:** Bloch orbitals `χ_k(r) = Σ_R e^{ik·R} χ(r−R)` become **complex** → the orbital side
  enters the `dcmplx`/Bloch world with a k-mesh. But that is the **same `dcmplx` machinery the pure-PW
  path already templates** (`Hamiltonian<dcmplx>`/`WaveFunction<dcmplx>`/`SCFIterator<dcmplx>`). So
  general-k GPW = a **Bloch-summed Gaussian IBS\<dcmplx\>** (lattice sum of real Gaussian integrals
  with `e^{ik·R}` phases) + the GPW Hartree, reusing the templated SCF stack.

**Verdict:** the unification holds, and it reveals that "Gaussian orbital basis" spans **both** `double`
(molecule/Γ) and `dcmplx` (general-k Bloch) via the templating you already have. GPW needs **no new SCF
machinery** — the Gaussian engine gains a Bloch-sum wrapper. Γ-point first (all-`double`), general-k
second (the periodic Gaussian IBS is the real new orbital-side piece there).

## 6. GPW vs the density fitting you already have (when does FFT win?)

They **coexist** as strategies; neither replaces the other:
- **Density fitting (aux Gaussians, 3-centre)** wins for **small–medium molecules** — compact aux
  basis, no grid overhead.
- **GPW (FFT-Poisson)** wins for **large and/or periodic** systems — O(N log N) Poisson, and periodic
  boundary conditions are *natural* for the FFT (periodic density fitting must wrestle the Coulomb
  divergence + lattice sums). For the **battery north-star** (periodic supercells), GPW's periodic-FFT
  Coulomb is the natural choice.

So: molecules → keep density fitting; solids / large supercells → GPW. The Hamiltonian factory picks
the Coulomb strategy from the structure type, same way it picks terms today.

## 7. The one genuinely new piece: collocation (and its adjoint)

Almost everything is reuse — Gaussian integrals ✓, PW FFT/Poisson ✓ (the `qchem.FFT` module + the PW
G-space Hartree), uniform grid ✓ (the PW grid; the Becke mesh is atom-centred for XC, distinct).
**The new code is collocation:**
- **Forward:** map `ρ = Σ D_ab χ_a χ_b` onto the uniform grid. Gaussian products are Gaussians
  (product theorem), each `ρ_ab` collocated within a cutoff radius.
- **Adjoint:** integrate `V_H(r)` back to `⟨χ_a|V_H|χ_b⟩` — the transpose of collocation.
- **Refinement (later):** multi-grid (map sharp vs smooth Gaussians to finer/coarser grids) — CP2K's
  efficiency trick; defer to a v2.

This is well-understood but fiddly, and it is the heart of GPW.

## 8. Phasing
- **Γ-point GPW (all-`double`).** Gaussian orbitals + PP (from the PP plan) + a `GPW_Hartree` term:
  collocate → FFT → Poisson → adjoint. *Gate: a Γ-point solid (or a molecule in a box) matches a
  density-fit DFT energy to grid-cutoff tolerance.*
- **General-k GPW.** Bloch-summed Gaussian IBS\<dcmplx\> + the templated `SCFIterator<dcmplx>` + the
  GPW Hartree. *Gate: a k-mesh band-structure / total energy vs a reference.*
- **(Deferred) GAPW** augmentation only if all-electron GPW is ever needed.

## 9. Reuse vs new
- **Reuse:** Gaussian integral engine (PG/libcint/MnD); `qchem.FFT`; PW G-space Hartree/Poisson;
  uniform grid; templated `Hamiltonian/WaveFunction/SCFIterator<T>`; molecular PPs (smooth density).
- **New:** collocation + adjoint (`GPW_Hartree` term internals); the Bloch-summed Gaussian IBS\<dcmplx\>
  (general-k only); Hamiltonian-factory wiring to pick the GPW Coulomb strategy by structure type.

## 10. Validation
- Γ-point GPW vs density-fit DFT on the *same* molecule/box → agree to grid-cutoff tolerance (the
  controlled approximation is the grid cutoff, not the physics).
- Energy must be variational-stable in the grid cutoff (converges as cutoff ↑).
- General-k: total energy / bands vs a PW or reference result.

## 11. Key file pointers (anchors for reuse)
- FFT/Poisson: `qchem.FFT` (`src/Common`); PW G-space Hartree in the PW DFT path
  (`UnitTests/PlaneWaveDFTUT.C` `HartreeVtilde`, and the `PW_Hartree` term).
- Templated SCF stack: `src/Hamiltonian` (`tHamiltonian<T>`), `src/WaveFunction`, `src/SCFIterator`
  (`tSCFIterator<T>`, already `double`|`dcmplx`).
- Gaussian integral engine: `src/BasisSet/Molecule` (PG / libcint / MnD).
- Coulomb-strategy precedent to mirror: `FittedVee` (`src/Hamiltonian/Internal/Imp/FittedVee.C`).
- Molecular PP enabler: [[doc/MolecularPseudopotentialPlan]].
