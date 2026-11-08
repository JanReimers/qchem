# V-ONE-OF-N: the coefficient record (design note, lane S 2026-11-07)

RECORD, cited by `doc/HamiltonianSprintPlan.md` step 5 and `doc/OOD-SOLID-Cleanup.md` V-ONE-OF-N.  **RULED (user 2026-11-07, "Yes" to all four
questions at the end):** the §6 surface including `MatrixForward`'s return type; ONE NON-template `FitCoefficients` (an empty
member is cheap; §2's `template <class T>` sketch is superseded); `Projector3` → `PairProjection<T>`; `Overlap`'s identity /
diagonal realizations wait for a sparse PW basis (§7).  Written before step 4 landed, so line numbers are as of 06e976a3.
**RULED (user 2026-11-08, sub-step 2, option (B)):** *"In general there are multiple levels of clients. I was thinking of the
client for the union data structure, so that means the Fitter derived implementations."*  So the record also has
`const rvec_t& Vector() const` and `const ΔG_Map& GMap() const`, each THROWING on a mismatch (§3a).  CLAUDE.md "Two levels of
client" (main d374041dd) is the convention.

## 3a. Two levels of client (ruled 2026-11-08)

The paired implementations, the code that built or was handed the producer and so already knows the realization, read it
through the typed accessors `Vector()` / `GMap()`: the `FunctionFitter` derived classes (`ConstrainedFF`, `DeltaFunctionFitter`,
`GridDensityFitter`, the ortho fitter), their samplers (`DensitySampler_*` on a δ fit), `Hubbard_U`'s own integrator, and the
fitters' own adjoints.  None of them branches on which member is filled; a wrong call throws.  Everything above them (the SCF,
the mixers, the Hamiltonian terms) uses the dispatchers only (`+=`, `*=`, `MaxAbs`, `Dot`, `Contract` / `ContractAdjoint`).
A grep for `.Vector()` and `.GMap()` lists exactly the realization-aware code.  The alternatives the user did not take: (A) an
operation on the record per numeric consumer (a solve, a pointwise map, a positivity guard), most of them meaningless for a
ΔG map; (C) `MatrixForward` keeps returning `rvec_t`, which would have changed §6 and closed R1.0q by declaring it out of scope.

Read-only, from the tree at lane-s-h-sprint 06e976a3 (step 4 landed on it, not yet on main).  Nothing built or run for this note.
It touches three tier-1 hubs (`qchem.ChargeDensity`, `qchem.Fitting.FunctionFitter`, possibly `qchem.BasisSet` through
`Projector3`), so the user rules on the surface in §6 before code.  The code is sprint step 5, as the controller ruled.

## 1. The problem, in one sentence

One operation, "this density's (or this field's) coefficients in that fit basis", has two realizations chosen by the
PRODUCER: an `rvec_t` (a Gaussian fit's coefficients; a δ fit's values at its points) or a `ΔG_Map` (a plane-wave fit's
ρ̃(G)).  Today each realization has its own METHOD and its own FACE: `MatrixForward::Forward → rvec_t`,
`FourierDensity::GetFourierDensity → ΔG_Map`, `GetRepulsion3C → rvec_t | ΔG_Map` by overload, plus
`GetRhoOnGrid → rvec_t` (empty means "none").  So a client must know which basis it has, and it learns that by
cross-casting (`dynamic_cast<const FourierDensity*>`: 10 sites in src and 1 in a test, listed in §5).

## 2. The record

```
//! ρ's (or a field's) coefficients in ONE fit basis: exactly one realization is populated, chosen by the producer.
template <class T> class FitCoefficients
{
public:
    explicit FitCoefficients(rvec_t a);      //!< a Gaussian fit's c_a, or a δ fit's values at its points
    explicit FitCoefficients(ΔG_Map g);      //!< a plane-wave fit's coefficients ρ̃(G)
    // -- dispatching operations only; no member is public --
    FitCoefficients& operator+=(const FitCoefficients&);   // throws on a realization mismatch
    FitCoefficients& operator*=(double);
    double  Dot(const FitCoefficients&, const Metric&) const; // the mixers' inner products (Pulay, the residual norm)
    double  MaxAbs() const;
    size_t  NumCoefficients() const;
    ...
private:
    rvec_t itsA;    // one of these two is empty
    ΔG_Map itsG;
    friend ... the two dispatchers in §3 (Contract / ContractAdjoint), nothing else
};
```

- **Name.** `FitCoefficients` names what the caller gets.  I considered `FitField` and `Coefficients`, but the first suggests
  a field on points (exactly the conflation R1.0n/o rules out) and the second is too generic to grep.
- **Members:** `rvec_t` and `ΔG_Map` only.  `ΔGq_Map` (the response's q-shifted map) stays out: it has one producer and one
  consumer, both in the response route, and folding it in would put q on every client.  A third realization is added only
  when a producer needs it.
- **T:** the run scalar, as `ProjectedDensity<T>`, so `FitCoefficients<double>` on a molecular run never carries an empty
  `ΔG_Map` by type.  Open question for the user: whether one non-template record (an empty map is cheap, the user's
  V-BLOCKED-MATRIX ruling) is preferred.  I lean to non-template, for the same reason.

## 3. Clients reach it only through dispatching functions

- `Contract(P, D) → FitCoefficients` and `ContractAdjoint(P, const FitCoefficients&) → hmat_t<T>` are the two dispatchers,
  living beside `Projector3` (§4).  They pick the realization from `P` and fill or read the matching member; every other
  client uses the record's operations (§2) or hands it on.
- `MatrixForward<T>::Forward(D)` returns `FitCoefficients`, and `MatrixAdjoint<T>::Adjoint` takes it.  `qchem.Mesh.Integrator`
  is not on the HUBS list, but every fitter and density hub uses it, so I list it in §6 anyway.  With that, `Vee_Hartree`
  joins the pair: its V_H(G) is a `FitCoefficients` holding a `ΔG_Map`, and `R1.0q` closes on its ORIGINAL target.
- **Today's direct readers of a realization**, all of which go:
  - `Projector3::dense` is read outside its dispatchers by `fnorm` / `relative_fnorm` (`src/BasisSet/Imp/Projector3.C`, test
    metrics) and by its producers (`PG_LibCint/Evaluator.C`, `SymmetryAdapted_IBS.C`, the Gaussian and Radial IBS).  Producers
    may keep writing it; `DenseProjector3Integrator` is a dispatcher.  The row's "IrrepCD reads `dense`" is already gone
    (`FiniteIrrepCD::GetRepulsion3C` goes through `DensityProjector` since 6f71e1a6); the comment at `IrrepCD.C:156` says so.
  - `ProjectedDensity_G::Map()` returns the `ΔG_Map`; the ortho fitter reads it.  It becomes `ProjectedDensity<T>` carrying a
    `FitCoefficients` (and `ProjectedDensity_G` goes).

## 4. A mismatch throws

`Contract`'s `assert(g.dense.empty() && …)` and `ContractAdjoint`'s `assert(g.applyAdjoint && …)` become throws (an assert
is dead in Release, CLAUDE.md).  The record's operations throw on a realization mismatch (`rvec_t += ΔG_Map`).  The checks sit
once per call, not per element, so they cost nothing measurable.

**Projector3's rename.**  It is the map from an orbital-pair product to a fit basis's coefficients, 3-, 2- or 1-index by
realization (dense ⟨ab|c⟩, PW columns plus kernel, or closures).  I propose **`PairProjection<T>`**: what it projects (a pair
product), and what it does.  The alternative `FitTensor<T>` is shorter but says "tensor", which the closure realization is not.
The rename rides with this change (one module, `src/BasisSet/Projector3.C`, and its includes), not as a separate sweep.

## 5. D2 and R1.0n/o close inside it (user rulings 2026-11-07)

**D2: `FourierDensity` becomes a `tChargeDensity`, with Fourier an implementation detail.**  Its methods map onto the record:

| today (`FourierDensity`, conditional base `FourierDensityBase<T>`) | becomes (on `tChargeDensity<T>`) |
|---|---|
| `GetFourierDensity(const cFIT_SF_ABS&) → ΔG_Map` | `Coefficients(const FIT_SF_ABS<T>&) → FitCoefficients` |
| `GetRepulsion3C(const cFIT_CD_ABS&) → ΔG_Map` | `RepulsionCoefficients(const FIT_CD_ABS<T>&) → FitCoefficients` (the finite `GetRepulsion3C(DensityProjector)` folds in) |
| `GetRepulsion3C_Raw` | goes: it is the composite's own un-star-averaged sum, an implementation detail of `Composite_Fourier` |
| `StarAverage(ΔG_Map&)` | goes off the face: its only callers are inside `CompositeCD.C` |
| `GetRhoOnGrid(fit) → rvec_t` (empty = none) | goes (R1.0n/o, below) |

So `FourierDensityBase<T>` and its `std::conditional_t` go, and every `dynamic_cast` to `FourierDensity` goes with them:
`PWTerms_Hartree.C:71`, `DensitySampler_Pair.C:100`, `DensityMixer.C:54,75,76,108`, `PolarizedDensityMixer.C:136,313`,
`FieldMixer.C:119`, `CompositeCD.C:453`, plus `IntegrationTests/PW/Boxes.C:102`.

**R1.0n/o: the `GetRhoOnGrid` question is wrong.**  Its four clients ask for the record instead:
- `DensityMixer` and `FieldMixer` build `GField{GetFourierDensity, GetRhoOnGrid}`, two views of one density, because the
  ρ̃ mixer wants coefficients and the real-space residual wants values.  They take `Coefficients(fit)`.  A client that needs
  ρ at INTEGRATION points asks the mesh (ρ[D], pin 28) through `DensitySampler`, which already owns that question.
- `PolarizedDensityMixer` does the same per channel (lines 265-266, 305-306).
- `DensitySampler_Pair` (`DensitySampler_Pair.C:102`) is the one client that wants values on a δ fit's own points, because the
  δ basis IS that point set.  For a δ fit those values ARE its coefficients, so `Coefficients(δfit)` answers it with no new
  name.  The user's distinction stays visible: the δ fit basis's points are the BASIS; the integration mesh is a separate
  object the sampler holds.  Not conflated, and not routed through each other.

## 6. The hub surface, for the user (nothing edited)

- **`qchem.ChargeDensity`:**
  - `tChargeDensity<T>` gains `Coefficients(const FIT_SF_ABS<T>&) const` and `RepulsionCoefficients(const FIT_CD_ABS<T>&) const`,
    both returning `FitCoefficients`.
  - The `FourierDensity` face and `FourierDensityBase<T>` are deleted, and the re-export `qchem.ChargeDensity.FourierDensity` goes.
- **`qchem.Fitting.FunctionFitter`:**
  - `ProjectedDensity<T>` carries a `FitCoefficients`, and `ProjectedDensity_G` is deleted.
  - `DensityProjector::Forward` / `ScalarProjector::Forward` keep their signatures: they return a `MatrixForward`, whose `Forward`
    changes (below).
  - The fitters' coefficients as VALUES (V-XC-ONE-TERM: `FittedVxc`'s v and ε fits handed back, not held) is a later,
    separate member each, `Coefficients() const → FitCoefficients`.  I list it here so the user sees where it leads; I would
    not add it until V-XC-ONE-TERM needs it.
- **`qchem.Mesh.Integrator` (not a listed hub, but every hub leans on it):** `MatrixForward<T>::Forward(D)` and
  `Forward(L) → FitCoefficients`; `MatrixAdjoint<T>::Adjoint(const FitCoefficients&)`.
- **`qchem.BasisSet` (only if `Projector3` is reached from it):** `Projector3<T>` → `PairProjection<T>`, and the
  `Overlap3C` / `Repulsion3C` return types follow the rename.

## 7. Overlap's dense / diagonal / identity: wait for sparse PW

`BasisSet::Overlap` of an orthonormal PW basis is the identity, so a dense `n×n` matrix of it is waste.  But no client pays
for it today: the PW SCF builds it once per block, and `LASolver` wants a dense S anyway.  The record's shape leaves room for
an `identity` / `diagonal` member when a sparse PW basis exists to need it.  **Recommendation: wait.**  Adding it now is a
third hub (`Orbital_1E_IBS`) for no measured gain.

## 8. Pin 32: each fitted term this touches is gradient-consistent

Each fitted term's H must be the gradient of its own E: one object on both pairs.

| term | shape | same object on both pairs? |
|---|---|---|
| `Vee_DFT` (normal) | Oᵀ K O | yes: one `Overlap3C` (forward `Contract`, adjoint `ContractAdjoint`), `K` the `CoulombOperator` |
| `Vee_DFT_Dunlap` | Cᵀ G⁻¹ C | yes: `ConstrainedFF`'s one `DenseProjector3Integrator` per block (forward via `DensityProjector`, adjoint in `Repulsion`) |
| `Vee_Hartree` (periodic) | Oᵀ k O, k = 4π/\|G\|² diagonal | yes: one `Repulsion3C` (`Contract` / `ContractAdjoint`); joins the `MatrixForward` pair here |
| `FittedVxc` | v fit for H, ε fit for E | **NOT gradient-consistent, and this change does not fix it**: H = Oᵀ S⁻¹⟨f\|v_xc(ρ)⟩ is not ∂/∂D of E = Tr(D ⟨ε_fit⟩).  That belongs to step 5's V-XC-ONE-TERM / VxcProductSpace, which owns the XC fit.  I flag it here so the record does not appear to settle it. |

The record changes none of these shapes: it changes only the type the coefficients travel in.

## 9. Order inside step 5, and what it moves

1. `FitCoefficients` and the two dispatchers with throws, plus `PairProjection` (the rename only).  0 moves.
2. `MatrixForward` / `MatrixAdjoint` on the record; `Vee_Hartree` joins the pair; R1.0q closes.  0 moves (same arithmetic).
3. D2: `Coefficients` / `RepulsionCoefficients` on `tChargeDensity`; the `FourierDensity` face, its conditional base and the 11
   casts go.  0 moves.
4. R1.0n/o: the four `GetRhoOnGrid` clients move to `Coefficients`, and `GetRhoOnGrid` goes.  0 moves, because a δ fit's
   coefficients are the values it returns today.  The panel runs here, since the mixers are inside the SCF loop.
5. Then V-XC-ONE-TERM, and V-NONLINEAR-MESH, as planned.

## Questions for the user

1. Is the surface in §6 approved, including the `MatrixForward` return type?
2. `FitCoefficients<T>` (templated on the run scalar) or one non-template record?
3. `PairProjection<T>` as Projector3's new name?
4. Is "Overlap waits for sparse PW" (§7) right?
