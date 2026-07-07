# Fitting / PW-DFT-fit Cleanup Plan

Follow-up cleanups surfaced while reviewing the plane-wave DFT-fit harmonization (see
`doc/MolecularPP_HarmonizationFindings.md`). None block current use; this is a menu to sequence
deliberately. Each item notes scope and whether it's bit-identical.

## Guiding principle — the `dynamic_cast` criterion

The deciding question for every cast is **why**, not what-type:

- **"I want more"** — broaden a narrow handle to a richer *capability* you intend to invoke
  (e.g. `rFIT_CD_ABS → FIT_CD_NonOrtho` to reach the Coulomb metric). This is fine; it's a request.
  Give it a good throw message.
- **"what are you"** — interrogate concrete identity to *branch* on it (or cross-cast between siblings
  that only coincidentally share a concrete class). This is the smell; the cast is standing in for a
  decision the object should own. Replace with a virtual predicate (`isOrtho()`), a redesign, or a
  guaranteed-by-contract requirement.

Abstract-vs-concrete target is only a *proxy* for this; intent is the real test. This is the criterion
for the system-wide cast survey (item C).

---

## Items

### A. `isOrtho()` + rename the fitter factories  — ✅ DONE
Landed: `isOrtho()` on the `FIT_CD_ABS`/`FIT_SF_ABS` cores (molecular `Fit_IBS`→false, `PlaneWaveFit_IBS`→true);
factory guards the down-cast by the contract; `Make{Scalar,Density}Fitter` → 4 `Fitting::Factory` overloads.
Bit-identical (170 non-A + PW/DFT anchors).

Replace the core→NonOrtho `dynamic_cast` inside `MakeDensityFitter`/`MakeScalarFitter` (the two `<double>`
overloads) with an explicit, mandatory predicate:
- Add `virtual bool isOrtho() const = 0;` to the `FIT_*_ABS` cores (every fit basis must declare its metric).
- Factory branches on it; the `else` still casts to `FIT_*_NonOrtho`, but now **guarded by contract**
  (`isOrtho()==false ⟹ IS-A _NonOrtho`) rather than assumed. Keep a good assert message.
- Rename `Make{Density,Scalar}Fitter` → overloaded `Fitting::Factory(...)` (4 arg-typed overloads; distinct
  return types are fine; matches `BasisSet::Factory`/`LASolver::Factory`). Fallback names `FactoryCD`/
  `FactoryScalar` are not needed — the overload works (only caveat: a caller holding the *union* `Fit_IBS`
  would be ambiguous, but callers hold the narrow face).
Supersedes the earlier "trait return face" idea (that baked in real↔non-ortho; `isOrtho()` keeps them
decoupled). Bit-identical.

### B. Lift `ScalarFunction<double>` to the fitter *cores* + give the ortho fitters real `op()(r)`  — ✅ DONE
`ScalarFunction<double>` moved from `FunctionFitter_{Scalar,Density}_NonOrtho` up to the cores; scalar
`_NonOrtho` collapsed entirely (`FunctionFitterImp` derives the core; density `_NonOrtho` keeps only
`FitGetSelfRepulsion`/`Integral`/`ReScale`). The ortho fitters implement `op()(r)`/`Gradient` by inverse-
transforming their `ΔG_Map` via a NEW **`G_FieldEvaluator`** abstract capability (qcBasisSet) — the SOLID-DIP
seam: `PlaneWaveFit_IBS` implements it off its own `GetGCartesian` (owns B), the ortho fitter cross-casts its
held `cFIT_*_ABS` → it ("I want more"), so no reciprocal lattice leaks into the fitting interface and no
concrete `PW_Evaluator` cast. New test `OrthoFitterRealSpaceField` (op(r) vs independent transform + RhoOnGrid
+ finite-diff gradient). Existing paths bit-identical (170 + 33 anchors); the new `op()(r)` is added capability.

The fitted *field* (`ρ_fit`, `v_xc,fit`) is **real and evaluatable** even for plane waves —
`ρ_fit(r)=Σ_G c_G e^{iG·r}` is real (Hermitian coeffs). Distinct from the fit *basis*, whose complex
`e^{iG·r}` functions force `VectorFunction<dcmplx>`. So:
- Move `: public virtual ScalarFunction<double>` from `FunctionFitter_*_NonOrtho` **up to the cores**
  `FunctionFitter_Density<T>` / `FunctionFitter_Scalar<T>`.
- Implement `operator()(r)`/`Gradient` on `OrthoFunctionFitter`/`OrthoScalarFitter` by evaluating
  `Σ_dm map[dm] e^{i(B·dm)·r}` (direct sum / inverse-FFT off `itsMap`). **Not a stub** — real capability.
  Get `B·dm → cartesian` by *delegating to the held fit basis* ("evaluate this `FourierMap` at r"), which
  avoids an abstract→concrete cast into `PW_Evaluator`.
- Consequence: the **scalar `_NonOrtho` collapses entirely** (nothing left once `ScalarFunction` leaves);
  the density `_NonOrtho` keeps only the genuinely-Coulomb/charge extras (`FitGetSelfRepulsion`, `Integral`,
  `ReScale`). Also fixes the `ReScale` asymmetry (scalar core vs density `_NonOrtho`) and aligns the factory
  returns (both `<double>` overloads return `_NonOrtho`; ortho overloads return the core).
**Driver:** the GUI needs to plot `ρ_fit`, `v_xc,fit`, and especially the delta `v_xc − v_xc,fit` — which
*is* the spatial fit-residual (see item G on measuring fit quality). One feature, two payoffs.

### C. `dynamic_cast` survey (existing CLAUDE.md item)
Sweep all casts; classify each by the A-principle. Known finds so far: the seed CD→SF cross-cast (item F),
the `FittedCD::DoFit` capability-probe (de-greyed by item E), the ortho fitters' `ProjectedDensity_G`/
`ProjectedScalar_G` casts (these are "I want more" — they pass). Give the survivors good throw messages.

### D. Drop `Band_DFT_IBS<dcmplx>` from `PlaneWave_IBS`  — ✅ DONE (Option 1, GPW-aware)
Decision (user Q on GPW): the abstract `Band_DFT_IBS<T>` **module/interface is KEPT** — a future GPW basis
(Gaussian orbitals, PW/FFT density) is its intended implementer (as `<double>`), per the file's own design
note. Only the dead `<dcmplx>` **base on `PlaneWave_IBS`** is dropped (production went fully G-space via
`Band_FT_IBS`; nothing casts to `Band_DFT_IBS<dcmplx>` — the sole prod reference was a stale PWTerms comment,
fixed). The 3 real-space methods (`Overlap`/`Repulsion`/`Integral`(SF)) stay as concrete **test-only oracles**
(independent analytic cross-checks of the shared FFT/Poisson machinery GPW will reuse). Bit-identical.

Production-dead: the G-space (`Band_FT_IBS`) route won; the real-space `Overlap(SF)`/`Repulsion(SF)`/
`Integral(SF)` methods survive only as independent test oracles in `PlaneWaveDFTUT.C`. Decide: delete the
base (and re-express/retire those tests) vs. keep them as independent oracles. Weigh their genuine
regression value against carrying the dead base.

### E. Seed owns its projection — move `ScalarSeedProjection_AO` out of `FittedCD`  — ✅ DONE
Landed per §E-spec (E now, F later): `NumericCD` IS-A `ProjectedDensity_AO`; `GetRepulsion3C` relocated onto it;
`FittedCD::DoFit` collapsed to one guaranteed cross-cast (seed overload deleted). Bit-identical (A_DFT/M_DFT/
M_HF_U_SAD/PlaneWave anchors). The CD→SF cross-cast now lives on `NumericCD::GetRepulsion3C` — item F.

`FittedCD` should not know seeding exists. Relocate the seed's overlap-fit projection onto the seed density
itself (`NumericCD`), making it a first-class `ProjectedDensity_AO` — mirroring how the PW seed
(`FourierSeedCD`) already owns its `GetFourierDensity`. Bit-identical. **Detailed spec below (§E-spec).**
NOTE: this relocation carries the CD→SF cross-cast with it (item F) — do F's design first if we want to
avoid relocating-then-reworking, or accept the relocated smell and fix in F.

### F. Fix the seed-fit Coulomb round-trip + the CD→SF cross-cast  — ✅ DONE
Design (user-agreed): "unconstrained fit metric" = a STRATEGY dispatched by polymorphism, not a runtime enum.
New `ProjectedDensity_AO::GetUnconstrainedFit(fbs)` = what the fitter calls; **default** = the Coulomb solve
`J⁻¹·GetRepulsion3C` (real densities inherit unchanged; composites still sum the RHS then one `J⁻¹`). The seed
**overrides** it with its overlap fit `S⁻¹⟨f|ρ⟩` directly — the fake `J⁻¹·(J·e)=e` round-trip AND the Coulomb-
face (`no`) cross-cast are GONE; only one honest overlap cross-cast (an "I want more" for the S metric) remains.
`GetRepulsion3C` kept for the composite RHS-summing (throwing default, so the seed needn't stub it); `J⁻¹` moved
from the fitter into the projection default; the fitter now owns only the Dunlap constraint. Not bit-identical
by design (drops `J·J⁻¹` fp-noise from the seed) but the SAD/DFT pins held — no re-pin needed. 171 + 34 anchors.

Root cause, confirmed in code: the seed's `GetRepulsion3C` returns `J·e` (`e=S⁻¹⟨f|ρ⟩`, the overlap-fit),
and `ConstrainedFF::DoFitUnconstrained` immediately computes `c0 = J⁻¹·(J·e) = e`. **The Coulomb metric `J`
cancels.** So the seed's fit is really *overlap-fit `e` + Dunlap charge constraint*, faked through the
Coulomb-metric `ProjectedDensity_AO` interface — and the `FIT_CD_ABS → FIT_SF_NonOrtho` cross-cast is the
symptom (the seed needs the `S` metric but is handed only the `J`/CD face, so it reaches sideways through
the `Fit_IBS` union). Proper fix: express the matrix-free seed fit as what it is — an overlap-metric fit +
charge constraint — using the fit basis's `FIT_SF_NonOrtho` face **directly** (no `J` round-trip, no cross-
branch cast). Interacts with E (changes what face the seed presents), so design F before/with E. Not a
mechanical change.

### G. `ProjectedScalar_AO` → `ProjectedScalar_R`  — ✅ DONE
Renamed the scalar projection to `_R` (carries only a real-space field; a Slater/BSpline basis uses it
identically). Density `ProjectedDensity_AO` left `_AO` (it genuinely exposes AO 3-centre `GetRepulsion3C`).
Pure rename, bit-identical.

The scalar one carries *only* a real-space field (`GetScalarFunction`) — nothing AO-specific; a Slater/
BSpline basis uses it identically. Rename to `_R` (real-space, paired with `_G`). Separately decide whether
to align the *density* `ProjectedDensity_AO` (which genuinely exposes AO 3-centre `GetRepulsion3C`) to `_R`
for symmetry, or leave it `_AO` as the more literal name. Cosmetic; bit-identical.

### H. Derive `nUniform` from a grid `E_cut` (Nyquist)
`MeshParams::nUniform` (default 20) is the last manual aliasing knob — the PW *FFT* grid already self-sizes
(`AutoGrid = 4·maxComp+1`), but the *real-space integration mesh* forces the user to guess `n`. Derive it:
`n ≳ a√(2E_cut)/π` (×2 for a density-bandwidth field). User sets a physical `E_cut`; the mesh follows.
The GPW-flavored item (see `MolecularPP_HarmonizationFindings.md` §6.4).

### I. GGA prerequisites (two seams to lay before gradient functionals)

**I.1 — `E_xc`/`V_xc` consistency (retire the ¾-virial special case).** Route `E_xc = ∫ε_xc·ρ` through the
fitter uniformly instead of the LDA-exchange `¾⟨ρ|v_x⟩` shortcut. The ¾ virial breaks for gradient
functionals. Existing TODO.

**I.2 — The grid-sizing seam (`relCutoff`).** For GGA the Vxc fit grid must be denser (the gradient
enhancement adds bandwidth to both `v_x` and `v_c`), and *how much* denser is a property of the **functional
type**, which only the Hamiltonian side knows. The creation point and the functionals are **co-located** in
`Ham_PW_DFT::BuildTerms` (a few lines apart), so this is a local seam, not deep threading.
- **Hard constraint — pass a *number*, never the functional.** `CreateVxcFitBasisSet` lives in `qcBasisSet`;
  `ExFunctional` lives in `qcHamiltonian`, which *depends on* `qcBasisSet`. Handing an `ExFunctional&` to a
  basis method would be a **library cycle** (linker-rejected). So the Hamiltonian distills the functional's
  appetite to a scalar and passes that; the basis stays functional-agnostic.
- **The seam:** `virtual double ExFunctional::GridCutoffFactor() const {return 1.0;}` (LDA→1; GGA→~1.5–2, the
  CP2K `REL_CUTOFF` idea); a `relCutoff` **field on `MeshParams`** (default 1.0 — so `CreateVxcFitBasisSet`'s
  signature grows a *field*, not a positional arg — dovetails with item H making `MeshParams` the grid-knob
  vehicle). `BuildTerms` builds the functionals first, takes `max(exchange, correlation)` (**shared grid →
  the denser of the two**, per the exchange-dominates analysis), and passes it. The fit basis scales its grid
  `E_cut` by `relCutoff` (it already holds `itsEcut`).
- **Scaffold it at `relCutoff = 1.0` now, wired and live** — introduce `GridCutoffFactor()` (default 1.0) +
  the `MeshParams.relCutoff` field + the `BuildTerms` threading with the value 1, alongside the denser-grid /
  `E_cut`-derived-grid work so that `relCutoff=1` **reproduces today's grid bit-identically**. The seam is
  then visible and exercised (not dead), and the GGA coder's whole job on this axis is to override
  `GridCutoffFactor()` — the grid follows automatically. (Until the `E_cut`-derived Vxc grid lands, the field
  is threaded but the grid stays the difference set; land them together so `relCutoff` is consumed from day
  one.)
- Acceptance for any grid change stays a **grid-convergence study of ρ / a property vs a fine reference**, not
  ΔE_total (the fits are non-variational — see the fit-quality note in `MolecularPP_HarmonizationFindings.md`).

### (folded) `Integrals_Overlap` placement
Already achieved by the Stage-C refinement — `Integrals_Overlap` is untouched and inherited only where the
overlap metric is genuinely used (`FIT_SF_NonOrtho`, `Orbital_1E_IBS`); the fit cores are clean
`IrrepBasisSet<T>` markers. Just verify no stray inheritance; no move.

---

## §E-spec — Seed owns its `ProjectedDensity_AO` (item E, the bit-identical relocation)

Findings: `NumericCD` (`src/ChargeDensity/NumericCD.C`) is the molecular SAD seed — a
`tChargeDensity<double>` that IS-A `ScalarFunction<double>` (superposed atomic ρ(r)) + a charge; **not** a
`ProjectedDensity_AO` today (its overlap-fit was pushed onto `FittedCD` by the stage-3 `NumericCD` refactor).
Real densities are `ProjectedDensity_AO` via `ProjectedDensityBase<T>=conditional_t<double,
ProjectedDensity_AO, NoProjectedDensity>`. The seed is created at `src/ChargeDensity/Imp/Seed.C:96`. The
`DoFit(ScalarFunction,charge)` overload is called from exactly one place — the `else` fallback in
`FittedCDImp::DoFit`.

Changes (bit-identical):
1. **`src/ChargeDensity/NumericCD.C`** — base list `+ public virtual Fitting::ProjectedDensity_AO`; add
   imports `qchem.Fitting.FunctionFitter` + `qchem.BasisSet.Fit_IBS`; declare
   `double FitGetConstraint() const override {return GetTotalCharge();}` and
   `rvec_t GetRepulsion3C(const BasisSet::rFIT_CD_ABS*) const override;`. Update the header comment (overlap-
   fit is back on the seed, as a genuine projection, not a stub — reversing the stage-3 note).
2. **`src/ChargeDensity/Imp/NumericCD.C`** — implement `GetRepulsion3C` = the old `ScalarSeedProjection_AO`
   body with `itsRho → *this` (⚠ still contains the CD→SF cross-cast — that's item F). Add `Fit_IBS`+`Blaze`
   imports.
3. **`src/ChargeDensity/FittedCD.C`** — delete `DoFit(const ScalarFunction<double>&, double)`; rewrite the
   `DoFit(const rChargeDensity&)` doc (every finite density presents its own `ProjectedDensity_AO`).
4. **`src/ChargeDensity/Internal/Imp/FittedCDImp.C`** — delete the anon `ScalarSeedProjection_AO` and the
   `DoFit(ScalarFunction,charge)` overload; collapse `DoFit(rChargeDensity)` to
   `auto* ao=dynamic_cast<const Fitting::ProjectedDensity_AO*>(&cd); assert(ao && "…finite density must
   present its AO projection"); itsFitter->DoFit(*ao);`. Prune now-unused imports (drop `Fit_IBS`; keep
   `ScalarFunction` — `FittedCD` IS-A one — and `FunctionFitter`).
5. **`src/ChargeDensity/Imp/Seed.C`** — no change.

Bit-identity: `sf->Overlap(*this)` samples the same `NumericCD` object the old fallback passed as a
`ScalarFunction`; `FitGetConstraint()==GetTotalCharge()` (scale included) matches the old
`cd.GetTotalCharge()`. SAD-seeded iterations are byte-for-byte identical.

Verify: `ninja UTMain` + full `./UnitTests/UTMain`; anchors `A_DFT.*`, `M_DFT_Water`, `M_HF_U_SAD`
(energies + iteration counts). `ninja allTests`.

Payoff: `FittedCD` loses all seed knowledge; the line-66 cast de-greys to a guaranteed capability
requirement (no `else`); the molecular seed matches `FourierSeedCD`'s "own your projection" pattern.

---

### K. Densify the Vxc/CD fit-{G} + swap `ProjectedScalar_G` `FourierMap`→`cvec_t` (the one real fit-quality knob)

The **only** genuine fit-accuracy lever in PW DFT. Two deferrals collapse into it:
- **`cvec_t` projection.** `ProjectedScalar_G` currently carries a `FourierMap` (kept for byte-identity while
  the fitter route landed). The projection `⟨c|v_xc⟩` is really a plain `cvec_t` (one coeff per fit-{G}, PW
  orthonormal ⇒ `InvOverlap=I` ⇒ `c = ⟨c|v_xc⟩`). Swap `FourierMap`→`cvec_t` here (and mirror in the density
  `ProjectedDensity_G` if it stays symmetric). Pure interface-honesty once the grid is free to differ.
- **Denser {G}, Vxc *and* CD together.** Do them as **one simultaneous upgrade** (user pin) — the CD-fit and
  Vxc-fit grids move in lockstep. This is what makes `relCutoff` (I.2) and the `E_cut`-derived Vxc grid (H)
  *do* something; before this they reproduce today's difference-set grid.
- **Acceptance is grid-convergence, never ΔE_total.** The fits are non-variational — a *better* fit can *raise*
  E_total (Dunlap RI-J is a lower bound). Measure convergence of **ρ / a physical property vs a fine
  reference** (ties to G's spatial fit-residual `‖ρ − ρ_fit‖` and `‖v_xc − v_xc,fit‖` diagnostics).

Depends on: H (grid self-sizes from `E_cut`) + I.2 (`relCutoff` seam live). Lands the deferred half of the
Item-2 XC work. Non-bit-identical *by design* (the grid changes) → guard with a grid-convergence test, not an
energy anchor.

### J. Rename `FourierMap` → `ΔG_Map` (it's just a data structure)  — ✅ DONE
Type `FourierMap` → `ΔG_Map` (Unicode Δ; toolchain verified g++15/C++20); module `qchem.FourierMap` →
`qchem.Math.GMap`; file `src/Math/FourierMap.C` → `GMap.C`; CMake + `qchem.Math` umbrella note updated; the
3 `export import` lines + ~30 type uses swept. Method names (`MakeFourierDensity` etc.) left as-is (honest
FFT provenance). ninja dyndep recovery applied to build/Release. Pure rename, bit-identical (170 + PW/DFT).


`FourierMap` (`using = std::map<ivec3_t, dcmplx, IVec3Less>` in `src/Math/FourierMap.C`) is a plain
G-space coefficient map keyed by a reciprocal-index **difference** `ΔG = Gᵢ−Gⱼ`. The name `Fourier`
implies provenance ("where it came from"); the type doesn't care. Rename to `ΔG_Map` to emphasise it's a
container, not a transform.

Scope (user-chosen): **type + module, ASCII module name.**
1. `src/Math/FourierMap.C`: `export using FourierMap` → `export using ΔG_Map`; rename the **module** to an
   ASCII, non-Fourier name (`qchem.Math.GMap` — folds into the `qchem.Math` sub-lib; NOT unicode, so the
   module name/filename stay ASCII-safe for ninja/dyndep) and rename the file `FourierMap.C → GMap.C`.
2. Update the 5 `import qchem.FourierMap;` lines → `import qchem.Math.GMap;` (Band_FT_IBS, FourierDensity,
   FunctionFitter, + the re-export in OrthoFunctionFitter's comment; also the `qchem.Math` umbrella note in
   `Math.C`) and the CMake `FILES` entry.
3. Sweep the ~30 type-name uses `FourierMap` → `ΔG_Map` (signatures in Band_FT_IBS, PlaneWave_IBS,
   FourierDensity/FourierSeedCD/CompositeCD/IrrepCD, ProjectedDensity_G/ProjectedScalar_G, OrthoFunctionFitter
   `itsMap`). Pure rename — no type change, bit-identical.
4. ninja dyndep recovery on the other build dirs (module file renamed — see the recovery memo).

Note: the method names `MakeFourierDensity`/`GetFourierDensity`/`ForwardGrid` are a **separate** question
(they describe the FFT that produces the map, so "Fourier" there is arguably honest) — leave them unless we
decide the density-carrier should also shed the word. Verify: full `UTMain` + `ninja allTests`; nothing
numeric moves.

---

## Suggested sequencing (rough)

1. **B** (ScalarFunction→cores + ortho `op()(r)`) — enables GUI plotting + fit diagnostics; symmetrizes the
   split; self-contained.
2. **A** (`isOrtho()` + `Factory()` rename) — small, high-legibility, bit-identical.
3. **F then E** (seed-fit redesign, then relocation) — do F's design first so E doesn't relocate a smell;
   or E now + F later if we want the separation-of-concerns win sooner.
4. **G**, **J**, **(folded)** — cosmetic renames/verify, cheap, standalone; do whenever (J is a pure
   rename, bit-identical, but touches a module file → dyndep recovery).
5. **D** (drop `Band_DFT_IBS<dcmplx>`) — independent; decide the test-oracle question.
6. **H** + **I.2's `relCutoff=1` scaffold** together — when the `E_cut`-derived Vxc grid lands, wire the
   `GridCutoffFactor()`/`MeshParams.relCutoff` seam at value 1 so LDA is bit-identical and the GGA hook is
   live. **I.1** (E_xc/`¾`-virial consistency) is a separable increment.
7. **K** (densify Vxc/CD {G} + `cvec_t`) — *builds on* H+I.2; the first deliberately non-bit-identical
   step. **I.1** should land before it so `E_xc` is functional-consistent under the new grid.
8. **C** (`dynamic_cast` survey) — **dead last.** It's a codebase-wide *hardening* pass (wrap every surviving
   cast in a throwing, info-rich guard). Every prior item churns the cast landscape — E relocates, F deletes
   (CD→SF cross-cast), B reshapes the fitter faces, K swaps the projection carriers — so surveying earlier
   audits casts about to move. The *principle* (capability vs identity, top of this doc) guides the refactors
   as we go; the *survey* pins the final state once it stops moving.

All numerics-affecting steps (H, I, K) use a **grid-convergence acceptance test** (NOT ΔE_total — see the
fit-quality note in the harmonization findings).

---

## Code-review findings (review of the A–K work, 2026-07-06)

Two-session multi-agent review of `d34abd02..HEAD`. **No active correctness bug** in the shipped
LDA/Γ path — the bit-identity discipline held; the removed-behavior and cross-file finders came back
empty. What follows is latent fragility, per-iteration waste, and cleanup. Each is tagged with the item
that should absorb it (not all are C — that was the reviewer's first, wrong, instinct).

**Done already:**
- ~~#4 — `FittedEpsXc::GetMatrix` re-ran the eps_xc `DoFit` on every energy query (2× per SCF iteration).~~
  Fixed `3bc17cb8`: guard the fit behind an `itsFitVersion` serial (mirrors `FittedVxc`'s `newCD`); the
  per-irrep `Overlap` still runs every call. 175 green, bit-identical.

**→ Fold into C (the `dynamic_cast` survey):**
- ~~**#1 `OrthoScalarFitter::Overlap`** — `dynamic_cast<G_FieldEvaluator*>(bs)` on the caller-supplied orbital
  basis, `assert`-only → release null-deref UB for a future non-PW complex orbital basis.~~ Done `1e1e9d83`:
  reference-cast `dynamic_cast<G_FieldEvaluator&>(*bs)` → throws `std::bad_cast` (the `Symmetry::Atom` pry-out
  convention). The systematic sweep of the rest of the codebase's casts is still C.
- **#7 two owners of the fit-grid cross-cast** — `PW_XC` caches its own `itsFitGrid` (a `G_FieldEvaluator`
  cross-cast) for `GetEnergy`, in parallel with the fitter that already holds the same engine. **Update
  (after #6, `1e1e9d83`):** the *correctness* hazard is GONE — `E_xc` and `v_xc` now both read the one shared
  `itsRhoGrid`, so they can't be integrated/fit on different grids. What REMAINS is the structural duality:
  `PW_XC::itsFitGrid` and the fitter's `itsFitBasis` are two cross-casts of the same `fb`. Push the energy
  quadrature *behind the fitter* so there is one owner of the grid. **✅ DONE:** new narrow
  `Fitting::GriddedScalarFitter` face (the PW scalar fitter's refinement) exposes `Grid()` → the
  `G_FieldEvaluator` it already holds; `Factory(cFIT_SF_ABS)` returns it (PW_XC is its only consumer). `PW_XC`
  dropped `itsFitGrid` and now borrows the ONE grid via `itsScalarFitter->Grid()` (RhoOnGrid / Integral /
  GridPoints) — no second cross-cast of `fb`. Pure structural refactor, bit-identical (173 + 36 PW/DFT).
- **Refuted-residual (compile-time-over-runtime)** — `ProjectedDensity_AO::GetRepulsion3C`'s default throws.
  **Correction:** the reviewer's "pure-virtual / `=delete`" fix does NOT apply — this is an *override-exactly-
  one-of-two* contract (real densities override `GetRepulsion3C` + use the default Coulomb `GetUnconstrainedFit`;
  the seed overrides `GetUnconstrainedFit` and never provides a Coulomb RHS). Neither method can be pure without
  forcing an unwanted override on the other camp, and a virtual can't be `=delete`d — so the throwing default is
  the CORRECT tool (C++ can't express "override one of two" at compile time). No action required. IF the
  compile-time-over-runtime tenet is still wanted: make `GetUnconstrainedFit` pure + a protected `CoulombFit(fbs)`
  helper (real densities override with a one-line `return CoulombFit(fbs);`), costing ~3 boilerplate overrides.
  C-survey tier.

**→ Fold into K (guards for the future k≠0 / densification path):**
- ~~**#2 `PWVxcField::operator()(points)`** ignores its `points` argument, validated only by `size()==size()` —
  a future diagnostic sampling a different same-cardinality point set pairs values to the wrong points.~~ Done
  `1e1e9d83`: the field is now explicitly GRID-BOUND — pointwise `op(r)` throws (nothing samples it pointwise),
  and the bulk op asserts point-set IDENTITY (not just size) against the cached `GridPoints()`.
- ~~**#3 unasserted alias-free invariant in `GridCoeff`** — wraps `Δm` mod N with no `|Δm|<N/2` check.~~ Done
  `1e1e9d83`: `assert(2|Δm|<N per axis, "densify the fit grid / raise relCutoff")` — the future k≠0 XC block
  now fails loudly instead of silently aliasing.

**→ New perf item (PW-DFT hot path, no correctness impact):**
- ~~**#5 PW grid geometry recomputed every call** — `GridPoints()` (N-point mesh + 3×3 inversion) and
  `AutoGrid()` (O(nG) scan, hit O(n²) times in assembly).~~ Done `1e1e9d83`: `AutoGrid`/`FFTGrid`/`GridPoints`
  cached lazily on `PW_Evaluator` (depend only on `itsG`); `GridPoints()` returns `const&` (no per-DoFit copy).
- ~~**#6 `PW_XC::GetEnergy` recomputes `RhoOnGrid`** (a 2nd inverse FFT of the same `ρ̃`); `PWVxcField` copies
  the `ΔG_Map` by value.~~ Done `1e1e9d83`: newCD-guarded `RefreshRhoGrid` shares ONE inverse FFT across
  `CalcMatrix`/`GetEnergy`; the field holds the precomputed grid by `const&` and just maps the functional over it.

**→ Cleanup tier (do whenever you touch the file — same tier as G/J):**
- ~~**#8** `ExFunctional::GetVxcs` batch API is DEAD (0 callers).~~ Done: DELETED (YAGNI — reintroduce a batch
  API when a 2nd caller / vectorized path exists; the one genuine batch site, `PWVxcField`, samples the field's
  own `op(rvec3vec_t)`, not this). Dropped the now-unused `qchem.Blaze` import from the impl. Bit-identical.
- ~~#9 `ProjectedScalar<T>` vestigial → collapse to standalone `ProjectedScalar_R`.~~ Done `daa19555`.
- ~~#10 stale `ProjectedScalar_G` / `Make*Fitter` comments (`FunctionFitter.C`, `PWTerms.C:14`).~~ Done `daa19555`.
