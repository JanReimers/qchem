# AO/FT projection cleanup — spec (clean interfaces, zero NA-asserts)

> **STATUS: DONE (all 3 moves).** Full UTMain green (134) after each move; tree stays
> -Winconsistent-missing-override clean. Outcome:
> - **Move 1 (done):** `tDM_CD<T>` is now pure (ScalarFunction + DM contract); the AO face is a cross-cast
>   capability via a new `ProjectedDensityBase<T>` (=`ProjectedDensity_AO` for double, empty for dcmplx),
>   mirroring `FourierDensityBase<T>`. `IrrepCD`/`tComposite_CD` derive both capability bases;
>   `Polarized_CD` derives the AO one. `FittedCD::DoFit` now takes `const DM_CD&` and cross-casts to
>   `ProjectedDensity_AO` (assert non-null). The `IrrepCD<dcmplx>::GetRepulsion3C` NA-assert is gone (the
>   double-only 3-centre code is now behind `if constexpr`; the dcmplx body is an inert `return {}`).
>   `Composite`/`Polarized` `GetRepulsion3C` cross-cast their blocks to the AO face (sanctioned
>   abstract->abstract guard, like the existing FourierDensity cross-cast).
> - **Move 2 (done):** `FourierFunctionFitter` is now a standalone class (no abstract bases) with just
>   `DoFit(FourierMap)`/`Repulsion`/`Overlap`/`Write` — the ~12 NA-assert overrides are deleted.
>   `PW_Hartree`/`PW_XC` unchanged (already concrete). (Kept the named object; did not take the
>   delete-and-inline aggressive alternative.)
> - **Move 3 (done):** `DoFit(ProjectedDensity_FT)` removed from `FunctionFitter_Density<T>` and the
>   molecular `ConstrainedFF` NA override deleted. The `ProjectedDensity_FT = FourierMap` alias is kept.
> - **Assert state:** the Fitting layer now has ZERO `assert(false)`. Remaining CD asserts are OUT of the
>   three moves' scope and pre-existing: the symmetric `GetFourierDensity` double-path domain-guards
>   (IrrepCD/Composite) and the `AccumulateDirect/Exchange` HF-on-complex guards. The AO empty-side was made
>   fully inert (no assert) per Move 1; the FT empty-side keeps its domain-guard assert (untouched) — a minor
>   asymmetry, easily symmetrized later if desired.

Self-contained plan for a fresh session/agent. Goal: get the `<ρ|c>` projection's AO/FT distinction
OUT of the shared abstract interfaces (into the concrete `/Imp`-level types), and in doing so **delete the
`assert(false,"not implemented")` orphans** — each is an abstract forcing a method onto a concrete that has
no business with it.

## Not std::variant — the split aligns with the type system

`variant > polymorphic-base` was the right call *if you must unify two value-types at runtime*. But here you
don't: the AO/FT split already lines up with the existing capability structure — a **finite/molecular**
density carries an **AO** projection (`rvec_t`, Coulomb-metric solve); a **periodic/PW** density carries an
**FT** projection (`FourierMap`). So this is a **cross-cast to the capability you need**, exactly the
`PlaneWave_IBS : Band_FT_IBS` + `PW_External` pattern — distinct types, NO unboxing at all (cleaner than
variant: no `std::get`/visit, no boxing). The FT side ALREADY works this way (`PW_Hartree` does
`dynamic_cast<const FourierDensity*>(cd)`); the cleanup makes the AO side symmetric.

## Current asymmetry (the bug)

- `tDM_CD<T>` (`src/ChargeDensity/ChargeDensity.C:52`) `: public virtual Fitting::ProjectedDensity_AO` —
  so **every** density is forced to be AO-projectable, incl. the periodic one, which NA-asserts it:
  `IrrepCD<dcmplx>::GetRepulsion3C` (`Internal/Imp/IrrepCD.C:182`) `assert(false && "density fitting not
  used by the plane-wave path")`. **← NA-assert #1.**
- The FT side is already a capability, NOT forced on `tDM_CD`: `IrrepCD<T> : … , FourierDensityBase<T>`
  ("FourierDensity on the periodic (dcmplx) path; empty on the finite path", `Internal/IrrepCD.C:22`), and
  `PW_Hartree` cross-casts to it. This is the pattern to mirror.

## The three moves

### Move 1 — symmetrize the CD projection capabilities (kills NA-assert #1)

Make AO a cross-cast capability like FT, not a forced base of `tDM_CD`.

- `tDM_CD<T>`: **drop** the `: ProjectedDensity_AO` base and the `GetRepulsion3C` declaration
  (`ChargeDensity.C:54,95`). `tDM_CD` becomes pure (ScalarFunction + DM contract).
- Add a `ProjectedDensityBase<T>` mirroring `FourierDensityBase<T>`: primary (`<double>`, the finite path)
  `: ProjectedDensity_AO`; `<dcmplx>` specialization empty. `IrrepCD<T> : tDM_CD<T>, ProjectedDensityBase<T>,
  FourierDensityBase<T>` — so `IrrepCD<double>` has AO (real) + FT (empty); `IrrepCD<dcmplx>` has AO (empty)
  + FT (real). **Delete the `IrrepCD<dcmplx>::GetRepulsion3C` NA-assert specialization** (`Imp/IrrepCD.C:182-183`).
- Consumer: `FittedCD` currently up-casts the cd to `ProjectedDensity_AO` (`Internal/FittedCDImp.C:26`,
  `DoFit(const ProjectedDensity_AO&)`). Change `FittedCD::DoFit` to take the `tDM_CD`/`DM_CD` and **cross-cast**
  `dynamic_cast<const ProjectedDensity_AO*>(&cd)` (assert non-null — FittedCD is the finite/molecular path, so
  it always succeeds), then drive its fitter. (`FittedVee` already calls `DoFit(*cd)`.)
- Apply the same drop-forced-AO treatment to the other CD implementers that carry `GetRepulsion3C`:
  `Polarized_CD`, `tComposite_CD` (`CompositeCD.C`, `ChargeDensity.C` — the ones the FunctionFitter-split
  rippled through). Each gets the AO face via `ProjectedDensityBase<T>` (finite) instead of from `tDM_CD`.

### Move 2 — de-abstract `FourierFunctionFitter` (kills NA-asserts #2: the big cluster)

`FourierFunctionFitter` (`src/Fitting/FourierFunctionFitter.C`) derives BOTH `FunctionFitter_Scalar<dcmplx>`
and `FunctionFitter_Density<dcmplx>` (and `ScalarFunction<double>`), forcing it to implement ~12 molecular
methods it NA-asserts (`DoFit(ScalarFFClient)`, `DoFit(ProjectedDensity_AO)`, `ReScale`, `FitMixIn`×2,
`FitGetChangeFrom`×2, `FitGetSelfRepulsion`, `Integral`, `operator()`, `Gradient`). **It is used CONCRETELY,
never polymorphically** — `PW_Hartree`/`PW_XC` both declare `Fitting::FourierFunctionFitter fitter;` on the
stack (`PWTerms.C:116,152`), never a `FunctionFitter_Density<dcmplx>*`. So the abstract bases buy nothing but
the NA-asserts.

- **Drop the abstract bases.** Make `FourierFunctionFitter` a standalone class: `DoFit(const FourierMap&)`
  (store), `Repulsion(const obs_t<dcmplx>*)` (Hartree), `Overlap(const obs_t<dcmplx>*)` (XC), `Write`. Delete
  ALL the NA-assert overrides. `PW_Hartree`/`PW_XC` need no change (already concrete).
- The shared-abstract "harmonization" (FourierFunctionFitter IS-A FunctionFitter) is the thing that *forced*
  the NA-asserts; the call-pattern parallel (`DoFit(...)` then `Repulsion`/`Overlap`) remains in the source.
  Per the "minimal NA-asserts" rule this is the right trade.
- **Aggressive alternative:** the class is a thin map-holder + two delegations to `Band_FT_IBS`, and the
  terms already hold the basis. Could DELETE `FourierFunctionFitter` entirely and inline in `PW_Hartree`/
  `PW_XC`: `auto map=fd->GetFourierDensity(); return pw->Repulsion(map,Eh);` / `pw->Overlap(pw->ForwardGrid(vxc))`.
  Fewer types, zero NA-asserts. (Pick this if you don't value the named "PW fitter" object.)

### Move 3 — remove `DoFit(ProjectedDensity_FT)` from the abstract Density face (kills NA-assert #3)

`FunctionFitter_Density<T>` (`src/Fitting/FunctionFitter.C`) declares `DoFit(const ProjectedDensity_FT&)`,
which only the (now standalone) `FourierFunctionFitter` ever serviced — the molecular `ConstrainedFF`/
`IntegralConstrainedFF` NA-assert it, and **no one calls it through the abstract** (PW uses the concrete
fitter). Remove the `DoFit(ProjectedDensity_FT)` declaration from `FunctionFitter_Density<T>` and the
NA-assert override in the molecular Density impl. Keep the `ProjectedDensity_FT = FourierMap` alias (now used
only by the concrete `FourierFunctionFitter`).

## End state

- `FunctionFitter_Scalar<T>` / `FunctionFitter_Density<T>`: pure MOLECULAR fitters — `DoFit(ScalarFFClient)`+
  `Overlap`, and `DoFit(ProjectedDensity_AO)`+`Repulsion`/`Integral`/self-energy. No FT, no NA-asserts.
- `FourierFunctionFitter`: standalone PW G-space fitter (or gone/inlined). No NA-asserts.
- `tDM_CD<T>`: pure (no projection face). `IrrepCD<T>` carries the T-appropriate projection via the two
  symmetric capability bases (`ProjectedDensityBase<T>`, `FourierDensityBase<T>`); consumers cross-cast.
- **Zero `assert(false,"not implemented")` in the fitting/CD projection layer.** The AO/FT distinction lives
  in the concrete types + the `/Imp` files, not in the shared abstracts.

## Verify / green-gate

- Confirm (grep) no remaining polymorphic use of `FourierFunctionFitter` via `FunctionFitter_*<dcmplx>*`
  before Move 2 (currently none).
- Full `UTMain` after each move: molecular DFT (FittedVxc=Scalar, FittedCD/Vee=Density) + PW Si Γ −7.2273,
  NaF −20.3293, CsI −11.3868. Keep the tree -Winconsistent-missing-override clean.
- Sweep for any NEW `assert(false`/`"not implemented"` you might introduce — the whole point is to end with
  fewer, not more.
