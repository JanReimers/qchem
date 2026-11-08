# V-ONE-OF-N sub-step 3 (D2), redesigned: the density projects, the Hartree term applies the physics (design note, lane S 2026-11-08)

RECORD, for the user's approval before any hub code.  It supersedes `doc/Records/VOneOfNDesign.md` §5–§6 for D2, and §6's
`ProjectOnto` line (sub-step 2's ruled option 1) is revisited in §2 below.  Read-only from lane-s-h-sprint at d09e98dc1;
nothing built or run for this note.

**Revised 2026-11-08 for the user's rulings on the first draft**, relayed by the controller:
- **(a)** The record is named **`Projections`**.  The "slightly wrong" tension came from one type serving both the forward's
  OUTPUT (projections) and the adjoint's INPUT (coefficients).
- **(b)** **Fit coefficients do NOT use the one-of-N record**: *"if the fitter imp data structure for the fit coefficients needs
  to be different for two different fit basis set types, then just make different fitter imp classes for each."*  Each fitter
  implementation holds its own coefficient type privately, and the metric solve sits inside it (new §2a).
- **(c)** Q1 is YES: the PW Hartree term becomes the first-order fitted term on a diagonal Poisson operator.  Dunlap stays
  Gaussian-only.
- **(e)** Q2 and Q4 RULED (user 2026-11-08):
  - *"getting clean code is more important than pinning out to old number.  The priority 1) get right numbers, is really get the
    right absolute numbers, not 'don't let the numbers drift 1e-10 when you refactor'."*  This is now CLAUDE.md priority 1 (main
    648ad4edd).  So the matrix-free densities take the CLEAN route (§2), and the gate is energies converged within ε, with any
    anchor that moves beyond roundoff re-judged against an independent route (pin 10) and listed (§4).
  - Renames: `Vee_DFT` → **`Vee_Fit`**, `Vee_DFT_Dunlap` → **`Vee_FitDunlap`**, `Vee_Hartree`'s response half →
    **`Vee_Response`**.  The exact-J `Vee` is unchanged.  The note uses the new names from here on.
- **Still open:** Q5 (A vs B, §1a), being re-explained by the controller.
- **(d)** A correction (user 2026-11-08): 1/r₁₂ plays TWO roles, and the note must keep them apart.  One is the fit METRIC, a
  fitter's choice; the other is the PHYSICS.  §1a is new, and §1's "k(G) in the fitter" is withdrawn.

**The user's ruling, 2026-11-08, relayed by the controller.**
1. The density answers ONE contraction of its own D with a three-index object the caller hands in.  It returns projections,
   index- or G-keyed.  No fit concept and no "Coefficients" appear in its names, and the record is renamed.
2. The periodic projection is ρ̃(G) only, over the overlap tensor.  The 4π/|G|² weight leaves the tensor and is applied by the
   Hartree term, the way the Gaussian route's repulsion lives with Vee and the fitter's metric.
3. Aim: a PW fit basis is an ordinary fit basis whose overlap and repulsion happen to be diagonal.  Inventory what retires and
   what stays PW-specific.
4. Say where numbers can move, and which gates prove 0 moves.
5. Moving `DensitySampler_*` / `FittedCD` out of qcChargeDensity is a separate row, not this step.

## 1. The finding that makes it small: the PW Hartree term already IS `Vee_Fit`, written by hand

The molecular "normal" Hartree term `Vee_Fit` (`src/Hamiltonian/Internal/Imp/Vee_DFT.C`) runs on `GridDensityFitter`
(`src/Fitting/Internal/GridDensityFitter.C`), and its shape is exactly what the user describes:

| step | `Vee_Fit` + `GridDensityFitter` (δ mesh, molecular) | `Vee_Hartree` (PW / GPW, today) |
|---|---|---|
| projection | the density contracts D into the fitter's forward: ρ(r_g) | the density contracts D into `Repulsion3C`: **k·ρ̃(G)**, kernel baked |
| "the fit" | orthogonal basis: the projection IS the fit | orthonormal {G}: the projection IS the fit (`OrthoFunctionFitter` says so) |
| the physics | `CoulombOperator::Potential(q) → v`, OWNED BY THE FITTER (a misplacement too, §1a) | inside the tensor's forward |
| matrix | the same pair's adjoint, Σ_g w_g v_g χ_iχ_j | `ContractAdjoint(Repulsion3C, V_H)` (overlap adjoint) |
| energy | `CoulombOperator::HartreeEnergy(trDV, q, v)` = ½⟨q,v⟩ | ½Ω Σ\|V_H\|²/k, a hand-written Parseval loop in the term |

So the cure is not a new PW design.  The PW route becomes this same row with a diagonal Poisson operator:

- **The density** returns ρ̃(G) through the OVERLAP tensor.  That is today's `GetFourierDensity`, which already exists and is
  already star-averaged.
- **The fit** is ρ̃ itself.  An orthonormal basis needs no metric solve, and no kernel appears anywhere in the fitter (§1a).
- **The physics, k(G) = 4π/|G|²** with k(0)=0, is the Hartree term's: the second-round repulsion of §1a.
- **The term** is the first-order fitted Hartree term (`Vee_Fit`).

The kernel then lives in ONE place, the physics round, instead of three:
- the PW `Repulsion3CTensor`'s `kernel` column array;
- the GPW collocator's `coulomb` flag;
- `SeedCD` / `FourierMixCD`'s hand-applied `k*rt`.

Dunlap stays Gaussian-only (RULED (c)).

## 1a. 1/r₁₂ has two roles: the fit METRIC and the PHYSICS (correction (d), user 2026-11-08)

In the user's words, Dunlap *"does everything with (6D) repulsion integrals: G=⟨f_i|1/r₁₂|f_j⟩ and p_i=⟨f_i|1/r₁₂|ρ⟩ and once
the fit coefficients are calculated from G⁻¹ + constraints, it then finally applies the actual repulsion physics with a second round
of 3 center repulsion integrals."*

| role | what it is | where it belongs | Gaussian (Dunlap) | δ mesh | PW / GPW |
|---|---|---|---|---|---|
| **1. FIT METRIC** | the norm the fit minimises: G = ⟨f|M|f'⟩, p = ⟨f|M|ρ⟩, c = G⁻¹p (+ constraint) | the FITTER, as an implementation choice | M = 1/r₁₂ (`ConstrainedFF`: `InvRepulsion()`, `Repulsion3C` forward) | orthogonal: c = p/w, metric irrelevant | orthonormal: c = k⁻¹(kρ̃) = ρ̃.  **Metric irrelevant, no kernel in the fitter** |
| **2. PHYSICS** | V_ab = Σ_f c_f ⟨ab|1/r|f⟩ and E = ½ cᵀ⟨f|1/r|f'⟩c, or the first-order ½Tr(DV) | the HARTREE TERM's question | `FittedCD::GetRepulsion` / `GetSelfRepulsion`: the second round of `Repulsion3C` and `Repulsion()` | `CoulombOperator` (direct sum, radial Poisson) | k(G)=4π/\|G\|², diagonal |

So the PW route has two misplacements today, and the mesh route has the same kind:
- **PW:** the kernel sits in the density's PROJECTION (`Repulsion3C`'s forward, `SeedCD`/`FourierMixCD`'s `k*rt`), which is
  role 2 inside role 1's input.  After: the density projects through the OVERLAP tensor, the fit is ρ̃, and k(G) appears only in
  role 2.
- **Mesh:** `GridDensityFitter` OWNS its `CoulombOperator` and runs `Potential(q)` inside `DoFit`, which puts role 2 inside the
  fitter.  After: the operator belongs to the Hartree term.

**How the term applies the physics while the coefficients stay fitter-private (ruling (b)).**  The term may not read the fitter's c
(ρ̃ for PW).  There are three shapes:
- **(A) The physics is the BASIS's Coulomb integrals, and the term asks the fitter to contract with them.**  ✔ RECOMMENDED.  This
  is exactly Dunlap: `Vee_FitDunlap::MakeMatrix` → `FittedCD::GetRepulsion(bs)`, contracting the private c with the fit basis's
  ⟨ab|1/r|f⟩.  The PW basis's Coulomb 3-centre integral is ⟨ab|1/r|G⟩ = k(G)·⟨ab|G⟩.  So `Repulsion3C` KEEPS the kernel, but
  only in the physics round: its adjoint applies k, and nobody calls its forward any more.  Today it is the other way round: the
  kernel is forward-only, and the forward is the projection.  The PW fit basis's Coulomb metric ⟨G|1/r|G'⟩ = Ω k δ gives the energy.
  So k(G) lives in the PW basis's REPULSION integrals, as 1/r₁₂ lives in the Gaussian basis's, and is reached only when the Hartree
  term asks for repulsion.  `GSpaceDensityFitter` holds ρ̃ and no kernel code: its `Repulsion(bs)` contracts ρ̃ with
  `orb.Repulsion3C(fit)`'s adjoint.  For the mesh, the δ basis has no analytic Coulomb integrals, so its role-2 object stays a
  `CoulombOperator`.  It moves out of the fitter (next bullet).
- **(B) The term owns the operator and lends it per call**: `fitter.Repulsion(bs, phys)`, `fitter.HartreeEnergy(trDV, phys)`.
  This is the literal "applied by the term".  But each fitter class needs its own physics type (ruling (b): `rvec_t` for the mesh,
  ΔG for PW), so the term's fitter face is no longer one face.  It is a template or a cast.
- **(C) The fitter hands ρ̃ to the term, which applies k.**  This contradicts ruling (b).  Listed only to rule it out.

**The mesh misplacement (`GridDensityFitter` owns `CoulombOperator`).**  Under (A) the cheapest correct home is the term, which
builds the operator (the factory already chooses it) and hands it to the fitter's two physics calls.  That is (B), on the one fitter
whose coefficient type is fixed (`rvec_t`), so (B)'s typing objection does not arise there.  It touches:
- `GridDensityFitter`: `DoFit` stops calling `Potential`; `Repulsion` and `HartreeEnergy` take the operator;
- `Vee_Fit`: it holds the operator;
- its factory.

That is three files and 0 moves (the same calls in the same order).  **I recommend it in the same step (3b)**, since 3b already makes
`Vee_Fit` the shared term.  If the user prefers 3b smaller, it becomes its own `OOD-SOLID-Cleanup` row.

## 2. The density's one operation, and the record's name

```
// on tChargeDensity<T> (ChargeDensity hub) -- replaces ProjectOnto(ScalarProjector), GetRepulsion3C(DensityProjector),
// GetFourierDensity, GetRepulsion3C(fit), GetRepulsion3C_Raw, StarAverage, GetRhoOnGrid
virtual Projections Project(const BasisSet::PairProjector<T>& p) const=0;   // <c|rho> for every function c of p's basis
```

**`PairProjector<T>`** is a BasisSet face, not a Fitting one.  It is "the three-index object the caller hands in", per orbital
block.  It merges today's `Fitting::ScalarProjector` and `Fitting::DensityProjector`, which differ only in the scalar axes they
are templated on.

```
template <class TFit> class PairProjector               // qcBasisSet; implemented by the FITTERS (they build it from their fit basis)
{
    virtual const qcMesh::MatrixForward<double>& Forward(const Orbital_DFT_IBS<double,TFit>&) const=0;   // a real block
    virtual const qcMesh::MatrixForward<dcmplx>& Forward(const Orbital_DFT_IBS<dcmplx,TFit>&) const=0;   // a complex block
    virtual Projections Project(const ScalarFunction<double>& f) const=0;   // a matrix-free density projects its FIELD
};
```

- **A DM density** (`IrrepCD`, the composite) contracts its D into `p.Forward(block)` and sums blocks with `+=`.  The periodic
  composite then star-averages, as `GetFourierDensity` does today.  The only realization-specific line is that one star-average
  call, and the composite already owns the crystal point ops.
- **A matrix-free density** (`NumericCD`, `SeedCD`, `FourierMixCD`) has no D.  It answers `p.Project(*this)`, the projection face
  every fit basis already has (`OrthogonalFit`, `Overlap(ScalarFunction)`).  **The clean route (RULED Q2)** has no
  seed-specific or mixer-specific projector entries.  A G-space projector asks the function for a GENERIC capability:
  - **`AnalyticFourier`**: a `ScalarFunction` that knows its own Fourier components f̃(G).  It is data-free, beside
    `ScalarFunction`, and asked with the usual "I want more" cross-cast between abstract faces.  The projector then evaluates f̃ on
    its own {G}; a function without the capability is sampled on the projector's raster and FFT'd.
  - `SeedCD` has it: its structure-factor sum Σ_species S(G)·f̃_atom(G) is that function's analytic transform.
  - `FourierMixCD` has it trivially: it IS a G-series.  A mixed density projected onto the SAME basis it was mixed on may also
    return its own record directly, a basis-identity check that throws on another basis.  Both readings give the same numbers.
  - `NumericCD` (molecular, real space) has no transform, so it takes `Overlap(ScalarFunction)` on the Gaussian fit basis, as today.
    The S⁻¹ solve moves from the density into `ConstrainedFF` (§3).

  Numbers may move here, and that is accepted: the seed's ρ̃ evaluated by the projector rather than by the seed, and NumericCD's
  solve in a different place.  §4 lists them, and the gate is energies within ε (ruling (e)).

**The record's name: `Projections` (RULED (a)).**  It carries ⟨c|ρ⟩, the forward's output, and nothing else.  The adjoint's input
is a fitter's own coefficients (§2a), so the record no longer serves two roles.  The rename is mechanical: one module
(`src/Math/FitCoefficients.C` → `Projections.C`) and about 25 sites.

## 2a. Fit coefficients are fitter-private (RULED (b)): every place one crosses a generic interface today

The rule: a fitter implementation holds its coefficients in its own type and does its own metric solve (S⁻¹, G⁻¹, Dunlap's
constraint, or the identity).  Where two fit-basis types need two coefficient types, they get two fitter classes.  What leaves a
fitter is a matrix (`Repulsion`, `Overlap`) or a number (`HartreeEnergy`, `Integral`, `SelfRepulsion`), never a coefficient vector.

| # | where a coefficient vector crosses an interface today (tree at d09e98dc1) | who holds it | after |
|---|---|---|---|
| 1 | `qcMesh::MatrixAdjoint<T>::Adjoint(const FitCoefficients&)` and `Integrate(const FitCoefficients&)`, sub-step 2's signatures | the caller's coefficients | **reverts to `rvec_t`** (the vector realizations: mesh, raster, δ, Gaussian, occupations).  The coefficient-side `.Vector()` reads sub-step 2 added go away: `DenseMatrixIntegrator`, `ScreenedMatrixIntegrator`, `DensePairProjectionIntegrator`, `LowdinProjector` |
| 2 | `GSpacePairProjectionIntegrator::Adjoint(const FitCoefficients&)` (`.GMap()` read) | `Vee_Hartree`'s V_H | it **stops implementing `MatrixAdjoint`**.  It keeps `MatrixForward` (its output is `Projections` holding ρ̃) and gains a typed `Adjoint(const ΔG_Map&)`, called only by the PW fitter that owns it.  Nothing generic calls a G-space adjoint today: every fitter holds its integrator CONCRETELY (`ConstrainedFF`, `GridDensityFitter`, `DeltaFunctionFitter`, `FunctionFitterImp`, `Hubbard_U`), and the one face-level caller, `DensitySampler_Pair::BlockAdjoint`, is raster-only |
| 3 | `ConstrainedFF::Repulsion`, `FunctionFitterImp::Overlap` → `Adjoint(FitCoefficients(itsFitCoeff))` | `itsFitCoeff` (rvec_t) | `Adjoint(itsFitCoeff)`.  It was already private; only sub-step 2's wrap goes |
| 4 | `DeltaFunctionFitter::Contract` → `Adjoint(itsC)`; sub-step 2 made `itsC` a `FitCoefficients` | `itsC` | back to `rvec_t itsC`: the δ fitter's own type |
| 5 | `GridDensityFitter::Repulsion` → `Adjoint(FitCoefficients(itsV))`; `CoulombOperator::Potential(q)→v`, `SelfRepulsion(q,v)` | `itsC`, `itsV` (rvec_t) | `rvec_t`, private.  The `CoulombOperator` (role 2, physics) moves to the Hartree term and is lent to `Repulsion` / `HartreeEnergy` (§1a), so `itsV` is computed at the physics call, not in `DoFit` |
| 6 | `Hubbard_U` → `LowdinProjector::Adjoint(FitCoefficients(ch.W))`, `Integrate` | `ch.W` (its potential on occupation axes) | `Adjoint(ch.W)`.  Hubbard_U is its own paired implementation |
| 7 | `Vee_Hartree`'s field memo `itsField` (sub-step 2: `FitCoefficients`), its adjoint, and the `.GMap()` Parseval energy | the TERM holds V_H today, which is the leak | **a new PW density fitter** (`GSpaceDensityFitter`, replacing `OrthoFunctionFitter`) holds ρ̃, its fit, as a private `ΔG_Map`.  It has **no k(G) and no V_H**: the physics is the PW basis's Coulomb integrals (§1a (A)), reached when the term asks `Repulsion(bs)` (the `Repulsion3C` adjoint, k applied there) and `HartreeEnergy(trDV)` (½ ρ̃ᵀ⟨G|1/r|G'⟩ρ̃ = ½Ω Σ k\|ρ̃\|²).  The term holds no field |
| 8 | `DensitySampler::Matrix(bs, const rvec_t& v)`: v_xc on the δ points / raster, computed by the XC TERM pointwise | the term | **not this step.**  It is the XC side of the same leak (the term computes the potential coefficients and hands them to the "fitter", the sampler).  It is V-XC-ONE-TERM's job: the sampler takes the functional and keeps v |
| 9 | `FittedVxc`'s v and ε fits | the fitters (private) | unchanged |
| 10 | the density mixers (ρ̃ maps, raster values) | the mixers | these are PROJECTIONS (the density's ⟨c|ρ⟩), not coefficients, so they ride `Projections` (sub-step 4) |
| 11 | `ConstrainedFF::FitGetRepulsion(const ConstrainedFF*)`: another fitter's `itsFitCoeff` | same class | unchanged.  It is one class reading its sibling instance |

So after D2 the record sits only on the FORWARD side: `MatrixForward::Forward → Projections`, and `tChargeDensity::Project →
Projections`.  Its typed reads are done only by the fitters that consume their own projections:
- `ConstrainedFF` reads `.Vector()` before its J⁻¹ solve;
- the δ samplers read `.Vector()` as values at their points;
- `Hubbard_U` reads `.Vector()` as occupations;
- `GSpaceDensityFitter` reads `.GMap()` as ρ̃.

**Sub-step 2 (committed, unlanded) is reworked by this, not discarded.**  Its record, `GSpacePairProjectionIntegrator`, the
integrators' new signatures and every paired `.Vector()` / `.GMap()` read stand.  What changes is the rename, plus
`ProjectOnto(ScalarProjector) → Project(PairProjector)` on the hub.  So the ruling on sub-step 2's option 1 carries over, with the
new names.

## 3. Inventory: what retires, what stays PW-specific

**Retires.**

| what | where | why it can go |
|---|---|---|
| `FourierDensity` face, `FourierDensityBase<T>`, `NoFourierDensity`, the `qchem.ChargeDensity.FourierDensity` module | `src/ChargeDensity/FourierDensity.C` | `Project(PairProjector)` answers all of it |
| the 11 casts to it | `PWTerms_Hartree.C`, `DensitySampler_Pair.C`, `DensityMixer.C` ×4, `PolarizedDensityMixer.C` ×2, `FieldMixer.C`, `CompositeCD.C`, `IntegrationTests/PW/Boxes.C` | they ask the one hub method |
| `GetRepulsion3C(cFIT_CD_ABS)` on every density (IrrepCD_Fourier, Composite, SeedCD, PolarizedSeedCD, FourierMixCD, PolarizedMixCD) | ChargeDensity | the kernel moves to the operator (§1) |
| `GetRepulsion3C_Raw` and `StarAverage` on the face | `CompositeCD.C` | the average happens once, inside the composite's `Project` (it already does for ρ̃) |
| `GetRepulsion3C(DensityProjector)` + `CoulombMetric_ProjectedDensity`, `DensityProjector`, `ScalarProjector` | Fitting hub | merged into `PairProjector` (BasisSet) |
| `OverlapMetric_ProjectedDensity::GetUnconstrainedFit` (the S⁻¹ solve inside `NumericCD`) | Fitting hub + `NumericCD` | the seed projects (`p.Project(*this)`), and `ConstrainedFF` does its own solve |
| `ProjectedDensity_G`, `OrthoFunctionFitter`'s cast to it | Fitting hub | the fitter receives `Projections` |
| `OrthoFunctionFitter::Repulsion`, a throwing "obsolete" stub | Fitting | it becomes the real Hartree matrix (§1) |
| `Vee_Hartree`'s ground-state half: field memo, Parseval energy loop, `Volume()` | `PWTerms_Hartree.C` | it becomes `Vee_Fit` (see the next table for the response half) |
| the kernel in the PW tensor's FORWARD (`Contract`'s `kernel` columns) and GPW's `MakeCollocator(coulomb=true)` | `PairProjection.C`, `PW Evaluator`, `GPW_Evaluator.C` | nothing asks for a kernel-baked PROJECTION any more.  Under §1a (A), `Repulsion3C` survives with the kernel moved to its ADJOINT (the physics round) |
| `OrthoFunctionFitter` | Fitting | replaced by `GSpaceDensityFitter` (§2a row 7), the PW fitter class that owns its G-keyed coefficients |

**Stays PW-specific, because it genuinely is.**

| what | why |
|---|---|
| the diagonal Poisson kernel (`ReciprocalLattice::CoulombKernel`, k(0)=0) | the periodic 1/r: role 2 only (§1a), in the PW / GPW basis's Coulomb integrals.  `Repulsion3C`'s ADJOINT applies it (the physics round), and the Coulomb metric is Ω·k.  Never in a projection, never in fitter code |
| the IBZ star-average (`SymmetrizeGMap`) inside the periodic composite's `Project` | crystal point ops belong to the k-composite |
| the response's δV_H(G+q) (`TransitionRepulsion`, `CoulombKernel(m,q)`, `ΔGq_Map`) | a Bloch-q field.  It stays on `Vee_Hartree`'s response half until the response route is redone.  ⚠ This means `Vee_Hartree` does not fully retire in this step: its response methods move to a small `Vee_Hartree_Response`, or stay with a ground-state body that delegates |
| `GetRhoOnGrid` / raster ρ_DM(r) | sub-step 4 (R1.0n/o), unchanged here |
| `OrthoNormalScalarFitter` (the XC raster fit, wrap-around lookup) | not a Hartree object; V-XC-ONE-TERM's |

**What has to become scalar-generic for this.**  `Vee_Fit` is a molecular term (`rsmat_t MakeMatrix`).  The PW term has both
block corners (`cobs_t` and `robs_t`).  So the merged term is `tVee_Fit<T>`, or `Vee_Fit` gains the two periodic `MakeMatrix`
overloads.  It holds a fitter through the same narrow face it uses now (`DoFit` / `Repulsion` / `HartreeEnergy` / `Integral`),
whichever class implements it.  `GridDensityFitter` / `CoulombOperator` stay `rvec_t` (ruling (b)), and the energy's inner
product lives inside each fitter, so no new record dispatcher is needed for the Hartree term.  The mixers' inner product is still
sub-step 4's.

## 4. Where numbers can move, and the gates that say they did not

Under ruling (e), numbers MAY move.  The gate is the converged energy within the run's ε, and an anchor that moves beyond roundoff
is re-judged against an independent route (pin 10: CP2K / QE / ABINIT for PW and GPW, PySCF for molecular) and LISTED, never
refreshed silently.  Rows 1–5 below are reorderings of the same arithmetic, so roundoff (≲1e-15 relative) is expected there.  Rows 6
and 8 are real changes of route.

| # | change | today | after | expected |
|---|---|---|---|---|
| 1 | the kernel vs the k-sum | Σ_k (k·ρ̃_k) per block, then merged | k·(Σ_k ρ̃_k) | roundoff |
| 2 | the kernel vs the star-average | avg(k·x) over a star | k·avg(x).  ⚠ k is computed per member from G=B·m, so it is not bit-constant on a star | roundoff |
| 3 | G=0 | k(0)=0 drops it | same `CoulombKernel`, and ρ̃(0)=N/Ω is carried and multiplied by 0 | exact |
| 4 | the energy | ½Ω Σ \|kρ̃\|²/k | ½Ω Σ k\|ρ̃\|² (or ½⟨ρ̃, kρ̃⟩) | roundoff; it was already "not bit-identical" with the trace form (`eh_trace`) |
| 5 | GPW forward | a Coulomb closure: collocate, FFT, combine, ×k | the overlap closure (collocate / memo replay, FFT, combine), then ×k in the physics round (`Repulsion3C`'s adjoint, §1a (A)) | ×k on the same combined value is the same IEEE op; the CollocMemo replays one closure instead of two |
| 6 | seeds / mixed density | `SeedCD`/`FourierMixCD` hand-apply k·rt to their own ρ̃ | the projector evaluates their `AnalyticFourier` f̃(G) on its {G}, and the physics round applies k | the same analytic values should give ~roundoff.  Any difference lands only on the iteration-0 seed and the mixed iterates, so the CONVERGED energy should not move beyond ε |
| 8 | `NumericCD`'s S⁻¹ solve | inside the density (`GetUnconstrainedFit`) | inside `ConstrainedFF`, on `Overlap(ScalarFunction)` | same linear algebra in another place: roundoff, at the molecular seed only |
| 7 | `Vee_Fit` charge diagnostic | absent on PW | `HartreeChargeError` = Ω ρ̃(0) − N, a new diagnostic line | a report change, not a number change |

**Gates, in order.**
1. **Before anything is deleted:** a unit gate, in `UTLattice_BS`, the PW + GPW toys, checking that k·`Contract(Overlap3C, D)`
   equals `Contract(Repulsion3C, D)` to 1e-15 relative, per G.  It proves #1, #2 and #5 on their own.
2. `eh_trace` A/B on the PW and GPW anchors: G-space energy vs Tr(D V_H), the pin-32 shape check.  Before vs after, side by side.
3. The short tier, the long tier, and `-L "panel|periodic"` (the mixers move in this step), all green.  `scripts/anchordiff`
   is the INSTRUMENT that lists every move, not the pass criterion.  Each anchor that moves beyond roundoff gets a line in the
   commit and in `doc/Records`: its old value, new value, the run's ε, and its independent-route value (pin 10).  An iteration-count
   change is listed the same way.  The PW/GPW anchors that carry the claim are listed below; the `WILL_FAIL` panel cases are listed
   apart, each with its failure reason checked.
   - `GPW_NaF.Γ_Imp_eqExactResume` / `eqColdStart`
   - the MnO AFM-II anchor
   - PW Si, NaF and CsI
   - `PW/Boxes`, including the polarized seeds
   - `GPW.RawXCConsistencyFD`
   - the TiO₂ LR CP2K/QE gates, which are the response half, unchanged
4. RelChecked, vet-decks, and the hub audit, regenerated citing the ruling.
5. **One timed MnO run**, before and after, quoting the buckets that change:
   - "scf: V_H per-block ΔG_Map merge"
   - "V_H IBZ star-average"
   - "rho ball closure"
   - "E_H G-space pairing"

   The Coulomb closure call goes away, and the merge and the average now ride ρ̃'s.  So I expect a small saving, but that is a
   claim to measure, not to assume.

## 5. Order, if approved

- **3a.**
  - The rename (`FitCoefficients` → `Projections`).
  - `MatrixAdjoint::Adjoint` / `Integrate` revert to `rvec_t`, and the coefficient-side wraps go (§2a rows 1–6).
  - `GSpacePairProjectionIntegrator` gets its typed adjoint.
  - `PairProjector` in BasisSet; `tChargeDensity::Project`, implemented by every density.  The old faces stay alive one commit,
    forwarding to it.

  Expected roundoff; NumericCD's solve moves (§4 row 8).
- **3b.**
  - `GSpaceDensityFitter` (ρ̃ private, no kernel), replacing `OrthoFunctionFitter`.
  - PW / GPW `Repulsion3C` with the kernel moved into its adjoint, plus the PW fit basis's Coulomb metric (§1a (A)).
  - `Vee_Fit` scalar-generic, and holding the mesh `CoulombOperator` (§1a, the mesh misplacement).
  - The PW / GPW Hartree term becomes `Vee_Fit` on it, with gate 1 first, and `Vee_Hartree`'s response half becomes `Vee_Response`.

  Judged by §4's gates: ε, with moves listed and re-judged.
- **3c.** Delete everything §3 retires, the `FourierDensity` face and the 11 casts with it.  Regenerate the hub snapshot.
- Then sub-step 4 (R1.0n/o) as planned.  It is smaller now, since the mixers already ask `Project`.

Hubs touched: ChargeDensity (`Project`; the `FourierDensity` face deleted), Fitting.FunctionFitter (`ScalarProjector`,
`DensityProjector`, the two metric faces and `ProjectedDensity_G` leave), BasisSet (`PairProjector` added; `Orbital_DFT_IBS`'s
periodic `Repulsion3C` loses its last periodic client, while the molecular one stays: it is the real Coulomb ⟨ab|c⟩).

## Questions for the user

1. ✅ RULED (c): the PW Hartree term becomes the first-order fitted term on a diagonal Poisson operator; Dunlap is Gaussian-only.
2. ✅ RULED (e): the clean route.  Matrix-free densities project through a generic `AnalyticFourier` capability, or are sampled;
   numbers may move, judged against ε and re-judged against an independent route (§2, §4).
3. ✅ RULED (a): `Projections`.
5. (NEW, correction (d)) How the term applies the physics while the coefficients stay private: §1a's (A), physics = the basis's
   Coulomb integrals, as in Dunlap (my pick), or (B), the term lends its operator per call?  And does the mesh `CoulombOperator`
   move to the term in 3b (my pick), or as its own row?
4. ✅ RULED (e): `Vee_Fit` → `Vee_Fit`, `Vee_FitDunlap` → `Vee_FitDunlap`, and `Vee_Hartree`'s response half → `Vee_Response`
   (δV_H(G+q), until the response route is redone).
