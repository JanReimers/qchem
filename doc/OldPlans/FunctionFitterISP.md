# FunctionFitter ISP split — spec (Scalar vs Density)

> **STATUS: DONE.** Implemented and verified against the full UTMain regression (134 tests green;
> molecular DFT exercises both molecular faces, PW Si Γ −7.2273 / NaF −20.3293 / CsI −11.3868 exercise
> the FT both-faces path). Tree stays -Winconsistent-missing-override clean. Highlights vs the plan below:
> - `FunctionFitter<T>` replaced by `FunctionFitter_Scalar<T>` + `FunctionFitter_Density<T>`, each carrying
>   its own narrow `bs_t` (`FIT_SF_ABS` / `FIT_CD_ABS`). `FitFlavour` retired.
> - Shared coefficient/real-space machinery (operator()/Gradient/ReScale/FitMixIn/FitGetChangeFrom/Write)
>   lives in a new `FitImpBase<T,Face,FBS>` template (in `Internal/FunctionFitterImp.C`); the Scalar impl
>   (`FunctionFitterImp`) and the Density impls (`ConstrainedFF`/`IntegralConstrainedFF`) derive from it
>   with their own Face + FBS. **Both** `FitImpBase` instantiations are emitted in `Imp/FunctionFitterImp.C`
>   where the shared member bodies are visible (else the Density-face members link-fail).
> - `FourierFunctionFitter` now derives BOTH faces (one object, metric degenerate for PW).
> - Down-casts deleted: `AsFitIBS` (FittedVxc) and the `dynamic_pointer_cast<Fit_IBS>` in FittedVee.
> - `ProjectedDensity_AO::GetRepulsion3C` now takes `const FIT_CD_ABS*` directly (rippled through the
>   ChargeDensity implementers: IrrepCD, Polarized_CD, tComposite_CD); `FittedCD_Factory` takes FIT_CD_ABS.

Self-contained plan for a fresh session. This is **stage (b)** of the Fit_IBS ISP split (commit
`b380635c`): that split made `FIT_CD_ABS`/`FIT_SF_ABS` but left the `FunctionFitter` combined, so
`FittedVee`/`FittedVxc` still recover the full `Fit_IBS` via a contained `dynamic_pointer_cast`
(`AsFitIBS` in `FittedVxc.C`). Splitting the fitter removes that down-cast and gives "one metric at a
time." Do this BEFORE re-assessing the AO/FT data-structure scheme (held deliberately until after).

Tree is currently **-Winconsistent-missing-override clean** (the cache work cleaned it) — keep it so.

## The two orthogonal axes

`FunctionFitter<T>` (`src/Fitting/FunctionFitter.C`) conflates two independent axes:

| | **Scalar / overlap metric** | **Density / Coulomb metric** |
|---|---|---|
| fit | `DoFit(ScalarFFClient)` | `DoFit(ProjectedDensity_AO)`, `DoFit(ProjectedDensity_FT)` |
| contract | `Overlap(obs)` | `Repulsion(obs)`, `Integral`, `FitGetSelfRepulsion` |
| shared | `ReScale`, `FitMixIn`, `FitGetChangeFrom`, `Write`, `operator()(r)` (it IS-A `ScalarFunction<double>`) | ″ |

- **This split = the METRIC axis** (Scalar/overlap vs Density/Coulomb). Already half-named: the existing
  `FitFlavour::Unconstrained` ≡ Scalar, `ChargeConstrained` ≡ Density.
- **The AO/FT tension = the CONTAINER axis** (dense `ProjectedDensity_AO` + metric solve vs diagonal
  `ProjectedDensity_FT`=`FourierMap` that IS the fit). It lives ENTIRELY inside the Density face (the two
  density `DoFit` overloads). Leave it as two overloads here; do not unify (re-adds a Δm index to AO /
  pessimises — see `project_dft_fitting_boundary_pin`). Re-assess with the held AO/FT scheme afterward.

## Target interfaces

Two interfaces, each `: public virtual ScalarFunction<double>` (the fit is evaluable at r), each carrying
the shared ops typed to its own kind (`FitMixIn`/`FitGetChangeFrom` mix two fitters of the SAME face):

```cpp
template <class T> class FunctionFitter_Scalar  : public virtual ScalarFunction<double> {
  typedef std::shared_ptr<const BasisSet::FIT_SF_ABS> bs_t;   // overlap-metric aux basis (narrow!)
  virtual void      DoFit(const ScalarFFClient&) = 0;
  virtual hmat_t<T> Overlap(const obs_t<T>*) const = 0;       // Sum_a c_a <Oi|f_a|Oj>
  virtual void   ReScale(double); virtual void FitMixIn(const FunctionFitter_Scalar&,double);
  virtual double FitGetChangeFrom(const FunctionFitter_Scalar&) const; virtual std::ostream& Write(...) const;
};
template <class T> class FunctionFitter_Density : public virtual ScalarFunction<double> {
  typedef std::shared_ptr<const BasisSet::FIT_CD_ABS> bs_t;   // Coulomb-metric aux basis (narrow!)
  virtual void      DoFit(const ProjectedDensity_AO&) = 0;    // dense, metric solve (molecular RI)
  virtual void      DoFit(const ProjectedDensity_FT&) = 0;    // FourierMap, IS the fit (PW)
  virtual hmat_t<T> Repulsion(const obs_t<T>*) const = 0;     // Sum_a c_a <Oi|f_a/r12|Oj>
  virtual double    Integral() const = 0; virtual double FitGetSelfRepulsion() const = 0;
  virtual void   ReScale(double); virtual void FitMixIn(const FunctionFitter_Density&,double);
  virtual double FitGetChangeFrom(const FunctionFitter_Density&) const; virtual std::ostream& Write(...) const;
};
```

The narrow `bs_t` (FIT_SF_ABS / FIT_CD_ABS) is the payoff: it lets the down-cast go (below).

## THE KEY INSIGHT — `FourierFunctionFitter` implements BOTH faces

The metric distinction is real for AO (Gaussian RI: `InvOverlap·⟨f|f'⟩` vs `InvRepulsion·⟨ab|c⟩`+charge
constraint — two genuinely different solves), and **degenerate for FT**: plane waves are orthonormal/exact,
so projection IS the fit for both metrics. That's why `FourierFunctionFitter` already does Hartree
(`DoFit(ρ̃)→Repulsion`) AND XC (`DoFit(Ṽ_xc)→Overlap`) from the same `FourierMap`.

So: **`FourierFunctionFitter` implements both `FunctionFitter_Scalar<dcmplx>` and
`FunctionFitter_Density<dcmplx>`** (one object, both faces — PW collapses the metric), while the molecular
side has **two distinct impls**. This is the exact shape of the Fit_IBS split (the Gaussian aux basis
implements both `FIT_CD_ABS`+`FIT_SF_ABS`; the PW "fit basis" is implicit/degenerate). The AO/FT "tension"
is just the metric being real for Gaussians and degenerate for plane waves — the split makes that explicit
instead of buried in NA-asserts.

## Impl + Factory

- Impls (`src/Fitting/Internal/`): `FunctionFitterImp` (base, overlap metric) → the **Scalar** impl;
  `ConstrainedFF`/`IntegralConstrainedFF` (Coulomb metric + Dunlap charge constraint) → the **Density**
  impl. They already roughly cleave this way; re-cut so the Scalar impl carries `DoFit(ScalarFFClient)`+
  `Overlap` and the Density impl carries the two density `DoFit`s + `Repulsion`/`Integral`/self-energy.
- `FourierFunctionFitter` (`src/Fitting/FourierFunctionFitter.C`): change its base from
  `FunctionFitter<dcmplx>` to `: FunctionFitter_Scalar<dcmplx>, FunctionFitter_Density<dcmplx>`; it already
  has `DoFit(FourierMap)` + both `Overlap`/`Repulsion`. Its NA-asserted scalar/AO methods get pruned.
- Factory (`src/Fitting/Imp/FunctionFitter.C`): replace `MakeFunctionFitter(FitFlavour, fbs)` with
  `MakeScalarFitter(shared_ptr<const FIT_SF_ABS>&)` and `MakeDensityFitter(shared_ptr<const FIT_CD_ABS>&)`.
  The `FitFlavour` enum can retire (the type now IS the flavour); `IntegralConstrainedFF` vs `ConstrainedFF`
  is an internal Density-impl detail.
- `src/Fitting/Types.C`: keep `obs_t<T>=Orbital_1E_IBS<T>`; `fbs_t` splits into the two narrow faces (or
  drop it — the bs_t typedefs now name FIT_CD_ABS/FIT_SF_ABS directly).

## Consumer rewiring (the down-cast removal)

- **FittedVxc / FittedEpsXc / FittedVcorr** (`src/Hamiltonian/Internal/Imp/FittedVxc.C`): `itsFitter` becomes
  `FunctionFitter_Scalar<double>`; `MakeFunctionFitter(Unconstrained, AsFitIBS(bs))` → `MakeScalarFitter(bs)`
  with `bs` already the `FIT_SF_ABS` shared_ptr. **Delete the `AsFitIBS` helper + its down-cast.**
- **FittedVee** (`src/Hamiltonian/Internal/Imp/FittedVee.C`) → **FittedCD**
  (`src/ChargeDensity/Internal/Imp/FittedCDImp.C`): `itsFitter` becomes `FunctionFitter_Density<T>`;
  `MakeFunctionFitter(ChargeConstrained, …)` → `MakeDensityFitter(cdbasis)` with the `FIT_CD_ABS` shared_ptr
  threaded through. **Delete the FittedVee down-cast** (it currently `dynamic_pointer_cast`s FIT_CD_ABS→Fit_IBS).
- **PW_Hartree** (`Imp/PWTerms.C`): drives the Density face (`DoFit(ρ̃)`+`Repulsion`). **PW_XC**: drives the
  Scalar face (`DoFit(Ṽ_xc)`+`Overlap`). Both off the same `FourierFunctionFitter` (which is both faces).
- `FittedVxcPol` delegates to `FittedVxc` (unchanged shape).

## Green-gate anchors

Full `UTMain` (134) after each step; molecular DFT energy tests exercise both faces (FittedVxc=Scalar,
FittedCD/Vee=Density); PW: Si Γ −7.2273, NaF −20.3293, CsI −11.3868 exercise the FT both-faces path. Keep
the tree -Winconsistent-missing-override clean.

## After this

Re-assess the AO/FT container scheme (held). With the metric axis split out, the only remaining axis inside
`FunctionFitter_Density` is the AO/FT container (`ProjectedDensity_AO` dense vs `FourierMap` diagonal) — now
cleanly isolated, which is the point of doing this first.
