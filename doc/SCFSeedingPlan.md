# SCF Seeding — Design Spec

Improve how the SCF loop is seeded with an initial charge density. Two pain points, one
architecture. Written for a fresh session; everything below is grounded in the current code.

## 1. Current state (what actually happens today)

The seed enters as the **optional last argument** of the iterator:

```cpp
// src/SCFIterator/SCFIterator.C
tSCFIterator(const tbs_t<T>*, const ElectronConfiguration*, ham_t*, acc_t*, tDM_CD<T>* cd = 0);
```

`Initialize(cd)` (src/SCFIterator/Imp/SCFIterator.C:76) just forwards it:

```cpp
itsWaveFunction->DoSCFIteration(*itsHamiltonian, cd);   // cd == nullptr is legal
itsWaveFunction->FillOrbitals(0.0001);
itsCD = itsWaveFunction->GetChargeDensity();
```

So `cd == nullptr` ⇒ the Hamiltonian is built with **no density** ⇒ only the density-independent
terms (kinetic + external/nuclear) contribute ⇒ a **core-Hamiltonian guess**. There is **no SAD /
atomic-density / superposition machinery anywhere** (verified by grep).

Two callers, two behaviours:

- **Molecular** (`UnitTests/QchemTesterImp.C:72`): passes *no* `cd` ⇒ core guess. Clean, no
  boilerplate, works — but for ionic systems the core guess is far from the answer.
- **Plane-wave** (`UnitTests/PlaneWaveDFTUT.C:1060-1065`): hand-builds a uniform seed, every test:
  ```cpp
  hmat_t<dcmplx> D0 = blazem::zeroH<dcmplx>(n);
  for (size_t i=0;i<n;i++) D0(i,i)=double(Nelec)/double(n);   // D = (N/n) I  -> uniform rho
  auto* seed = new ChargeDensity::IrrepCD<dcmplx>(D0, pw, irr);
  cSCFIterator scf(bs.get(), &ec, ham, acc, seed);
  ```
  (The standalone `RunSCF` test-driver, line ~438, does the same as `R.rho[(0,0,0)] = Nelec/Omega`.)
  PW seeds uniform *on purpose* — the comment notes "Hartree+XC active from iteration 0, as real PW
  codes do." Core-guess (`cd=0`) would leave Hartree/XC off for iteration 0 (free-electron orbitals).

### The two problems
1. **UX / boilerplate.** The PW path makes every test/app hand-assemble a density matrix. The seed
   strategy should be a one-liner (or a sensible default), symmetric with the molecular path.
2. **Crude for charge-transfer.** Uniform (PW) and core-guess (molecular) are both far from an ionic
   ground state. NaF / CsI converge but **grind** — the SCF rebuilds the entire Na→F charge transfer
   from a structureless start.

## 2. Goals / non-goals

**Goals**
- A single `SeedStrategy` abstraction with a sensible **default per basis type**, so no caller
  hand-builds a density matrix.
- A **superposition-of-atomic-densities (SAD)** seed that captures atomic shell structure ⇒ far
  fewer SCF iterations, especially for ionic solids.
- Keep the converged answer **bit-identical**: the seed changes only the *path*, never the minimum.

**Non-goals (this work)**
- Not the qcMath refactor, not new accelerators, not forces.
- Not perfect oxidation-state inference — IonicSAD is an opt-in Phase 3 with a simple heuristic.

## 3. Architecture

The seed is conceptually one function: *given {Structure, basis, electron config}, produce an initial
`tDM_CD<T>` (or null for core-guess)*. Introduce a strategy + factory:

```cpp
enum class SeedStrategy { CoreGuess, Uniform, SAD, IonicSAD };

template<class T>
tDM_CD<T>* MakeSeedDensity(SeedStrategy, const tbs_t<T>* bs,
                           const Structure* st, const ElectronConfiguration* ec);
```

- `CoreGuess` → `nullptr` (today's `cd=0`).
- `Uniform`   → the centralized version of the current PW boilerplate (`D=(N/n)I`).
- `SAD`       → superpose neutral atomic densities (§3.2).
- `IonicSAD`  → superpose **ion** atomic densities (§3.3).

The iterator picks a **default** when none is supplied. Add the strategy to the ctor (or a small
overload) with a basis-resolved default:

```cpp
tSCFIterator(..., SeedStrategy seed = SeedStrategy::Default);
// Default resolves: molecular -> CoreGuess (today), PW -> Uniform (today),
//                   then BOTH -> SAD once §3.2/§3.3 land.
```

This deletes the PW test boilerplate (tests pass a strategy enum, or rely on the default).

### 3.1 The AO/FT axis returns — do NOT invent a third pattern

SAD is "sum atomic densities at the nuclei." The **concept** is shared; the **container bifurcates by
basis type**, exactly like the existing splits ([[project_functionfitter_isp_split]],
[[project_ao_ft_projection_cleanup]], [[project_dft_fitting_boundary_pin]]):

- **Molecular (`double`, AO):** block-diagonal superposition of atomic **density matrices**. Solve
  each unique element once (atomic LDA SCF), drop its converged DM into that atom's diagonal block of
  the molecular DM; off-diagonal blocks zero. Textbook SAD ⇒ an `IrrepCD<double>`.
- **Plane-wave (`dcmplx`, FT):** sum atomic density **form factors** with structure factors,
  `rho(G) = Σ_atoms ρ_atom(|G|) e^{-iG·R_atom}` ⇒ a Fourier density (`IrrepCD<dcmplx>`). This is the
  **same structure-factor assembly the pseudopotential already uses** (`MakeNuclear` / the local-PP
  form factor in `Integrals_Pseudo`/`PW_Pseudo`) — feed the *atomic valence density* radial function
  through that path instead of the PP.

So realize `MakeSeedDensity` as a **polymorphic capability split across the two basis libraries**
(molecule_BS provides the AO/block-DM face, lattice_BS provides the FT/form-factor face), with the
ChargeDensity/iterator layer agnostic — mirroring the cross-cast capability pattern already in the
tree. **Do not add a variant or a third mechanism.**

### 3.2 SAD — atomic densities

- Source: the existing atomic SCF (Atom basis + `tSCFIterator<double>`), solved **once per unique Z**
  (neutral, spherical, same LDA functional as the target), result cached like the GTH lookup.
- Molecular: cache the atomic DM in the atom's radial basis; map molecular basis fn → atom → block.
- PW: cache the atomic **valence** radial density; FT to ρ_atom(|G|) form factor.
- Spin: thread `Pol` — open-shell wants a spin-polarized atomic seed (atomic spin densities).

### 3.3 IonicSAD (Phase 3, opt-in)

Neutral SAD still leaves ~1 e⁻ of Na→F transfer for the SCF. IonicSAD seeds Na⁺ + F⁻ atomic
densities ⇒ near-converged start for NaF/CsI. Needs **oxidation states**: a simple
electronegativity/formal-charge heuristic, or user-supplied (could ride on the Structure alongside
the Zion/PP work — see [[project_hgh_params_table]]). Opt-in and user-overridable; the heuristic is
the only genuinely fuzzy piece in this whole plan, so isolate it.

## 4. Phased plan (each phase regression-safe)

- **Phase 0 — refactor, zero behaviour change.** Add `SeedStrategy` + `MakeSeedDensity`; move the PW
  uniform boilerplate into `MakeSeedDensity(Uniform,…)`; iterator default resolves PW→Uniform,
  molecular→CoreGuess. *Gate: every SCF test identical iteration count + energy.*
- **Phase 1 — SAD molecular** (block-diagonal atomic DM). *Gate: molecular DFT energies unchanged;
  iteration count ≤ today.*
- **Phase 2 — SAD plane-wave** (G-space form factors, reuse structure-factor assembly). Flip the PW
  default to SAD. *Gate: Si / NaF / CsI energies unchanged; iteration count strictly down.*
- **Phase 3 — IonicSAD** (opt-in, oxidation heuristic). *Gate: NaF/CsI iteration count down further;
  user-override path tested.*

## 5. Validation
- **Hard invariant:** converged total energy bit-identical across seed strategies (pin existing
  energy regressions; they must not move). The seed only changes the path.
- **Soft guards:** iteration-count assertions on NaF/CsI/Si with *generous* bounds (iteration count is
  sensitive to accelerator + tol; assert "≤ baseline", not an exact number).
- Keep `CoreGuess` reachable — some tests/edge cases may want the structureless start.

## 6. Risks / nuances
- **Total charge** must integrate to exactly `Nelec` (neutral SAD does by construction; ionic must
  balance the formal charges).
- **Representability:** PW seed populates only G-vectors present in the basis; molecular seed only the
  atom blocks present. Truncation must not lose charge (renormalize after projection).
- **Atomic-SCF cost** is negligible (≤ #species solves, cached) but must use the *same functional* as
  the target or the seed fights the SCF.
- **"Hartree active from iter 0":** SAD/Uniform give nonzero ρ at iter 0 (good); CoreGuess does not.
  SAD strictly dominates Uniform (it has shell structure), so SAD becomes the default once it lands.
- **Don't re-fork the AO/FT axis** — follow §3.1.

## 7. Open questions for the implementer
1. Atomic densities **live-solved** (always functional-consistent) vs **tabulated** (faster startup,
   like `saito.json`)? Recommend live-solved + cached; tabulate later only if startup cost bites.
2. Home for `MakeSeedDensity`: a tiny qcChargeDensity module for the dispatch, with the two faces in
   molecule_BS / lattice_BS (recommended, consistent with the tree) — confirm against the current
   library DAG before wiring.
3. IonicSAD oxidation source: heuristic vs Structure-carried formal charges vs explicit user arg.

## 8. Key file pointers
- Seed entry: `src/SCFIterator/Imp/SCFIterator.C:76` (`Initialize`), ctor `:61`.
- Molecular caller (clean default): `UnitTests/QchemTesterImp.C:72`.
- PW boilerplate to delete: `UnitTests/PlaneWaveDFTUT.C:1060-1065` (and `RunSCF` ~438).
- PW structure-factor assembly to reuse for SAD form factors: `PW_Pseudo` / `Integrals_Pseudo`
  (`src/Hamiltonian/Internal/Imp/PWTerms.C`, `src/BasisSet/Lattice_3D/Imp/PlaneWave_IBS.C`
  `MakeNuclear`/`PseudoG0Energy`).
- Density container: `src/ChargeDensity/Internal/Imp/IrrepCD.C` (`IrrepCD<T>`, the `tDM_CD<T>`).
