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

### A. `isOrtho()` + rename the fitter factories
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

### B. Lift `ScalarFunction<double>` to the fitter *cores* + give the ortho fitters real `op()(r)`
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

### D. Drop `Band_DFT_IBS<dcmplx>` from `PlaneWave_IBS`
Production-dead: the G-space (`Band_FT_IBS`) route won; the real-space `Overlap(SF)`/`Repulsion(SF)`/
`Integral(SF)` methods survive only as independent test oracles in `PlaneWaveDFTUT.C`. Decide: delete the
base (and re-express/retire those tests) vs. keep them as independent oracles. Weigh their genuine
regression value against carrying the dead base.

### E. Seed owns its projection — move `ScalarSeedProjection_AO` out of `FittedCD`
`FittedCD` should not know seeding exists. Relocate the seed's overlap-fit projection onto the seed density
itself (`NumericCD`), making it a first-class `ProjectedDensity_AO` — mirroring how the PW seed
(`FourierSeedCD`) already owns its `GetFourierDensity`. Bit-identical. **Detailed spec below (§E-spec).**
NOTE: this relocation carries the CD→SF cross-cast with it (item F) — do F's design first if we want to
avoid relocating-then-reworking, or accept the relocated smell and fix in F.

### F. Fix the seed-fit Coulomb round-trip + the CD→SF cross-cast  *(needs design)*
Root cause, confirmed in code: the seed's `GetRepulsion3C` returns `J·e` (`e=S⁻¹⟨f|ρ⟩`, the overlap-fit),
and `ConstrainedFF::DoFitUnconstrained` immediately computes `c0 = J⁻¹·(J·e) = e`. **The Coulomb metric `J`
cancels.** So the seed's fit is really *overlap-fit `e` + Dunlap charge constraint*, faked through the
Coulomb-metric `ProjectedDensity_AO` interface — and the `FIT_CD_ABS → FIT_SF_NonOrtho` cross-cast is the
symptom (the seed needs the `S` metric but is handed only the `J`/CD face, so it reaches sideways through
the `Fit_IBS` union). Proper fix: express the matrix-free seed fit as what it is — an overlap-metric fit +
charge constraint — using the fit basis's `FIT_SF_NonOrtho` face **directly** (no `J` round-trip, no cross-
branch cast). Interacts with E (changes what face the seed presents), so design F before/with E. Not a
mechanical change.

### G. `ProjectedScalar_AO` → `ProjectedScalar_R`
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

## Suggested sequencing (rough)

1. **B** (ScalarFunction→cores + ortho `op()(r)`) — enables GUI plotting + fit diagnostics; symmetrizes the
   split; self-contained.
2. **A** (`isOrtho()` + `Factory()` rename) — small, high-legibility, bit-identical.
3. **F then E** (seed-fit redesign, then relocation) — do F's design first so E doesn't relocate a smell;
   or E now + F later if we want the separation-of-concerns win sooner.
4. **G**, **(folded)** — cosmetic/verify, cheap, fold into whatever touches those files.
5. **C** (cast survey) — after A/E/F land (they remove the worst offenders).
6. **D** (drop `Band_DFT_IBS<dcmplx>`) — independent; decide the test-oracle question.
7. **H** + **I.2's `relCutoff=1` scaffold** together — when the `E_cut`-derived Vxc grid lands, wire the
   `GridCutoffFactor()`/`MeshParams.relCutoff` seam at value 1 so LDA is bit-identical and the GGA hook is
   live. **I.1** (E_xc/`¾`-virial consistency) is a separable increment. All numerics-affecting steps use a
   grid-convergence acceptance test (NOT ΔE_total — see the fit-quality note in the harmonization findings).
