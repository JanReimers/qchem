# Real vs complex — what decides the scalar type, and how it flows

**Status: PLAN, no code beyond Step 0.**  Written 2026-08-16 out of the GPW runtime work
(doc/GPWPlan1.md "THE RUNTIME GAP, MEASURED"), which surfaced the question but does not settle it.
The decision to write this instead of starting: the change is large, and several open OOD items sit
on the very faces it would restructure (see "Prerequisites").

---

## 1. What the physics dictates

A block can be real **iff its basis functions can be chosen real AND every Hamiltonian term
restricted to it is real-symmetric.**  A conjunction — the two halves fail for different reasons and
are answered by different objects.

| case | what forces complex | decided by | TR survives | structural side-effect |
|---|---|---|---|---|
| Bloch k ∉ TRIM | \f$\chi^k=\sum_R e^{ikR}\chi\f$ | **irrep** (`BlochQN`) | yes (pairs ±k) | none |
| spin–orbit \f$\xi\,L\cdot S\f$ | \f$L=-i\,r\times\nabla\f$ | **term** | yes → Kramers | **large**: spin channels merge into 2-spinors, `ms` stops being a good label, double groups |
| magnetic field / vector potential (A·p, GIAO, current-DFT) | \f$p\to p+A/c\f$ | **term** | **no** | none, but no Kramers to lean on |
| molecular 4-/2-component Dirac | \f$\alpha\cdot p\f$ coupling | **term** | yes | large/small structure |
| ATOMIC Dirac (our `A_DHF`) | — stays REAL | convention | — | the \f$i\f$ absorbs into the small component; works by spherical symmetry |
| electric field | — stays REAL | — | yes | none.  "External fields" is NOT uniform: E real, B complex |
| complex harmonics / complex point-group irreps (E of C₃, C₄…) | — convention only | — | yes | real combinations always exist for a real H; the SALC machinery already takes them |
| non-collinear spin / spin spirals / GHF | rotating quantization axis ⇒ spinor | **ansatz** | depends | 2-component |
| complex SCF instability (real H, complex minimum) | the nonlinear functional, not the operator | **ansatz policy** | yes | none |
| CAP / resonances | non-Hermitian complex-SYMMETRIC H | term | no | different linear algebra entirely — not `LASolver`'s problem |

Two things worth stating because they are easy to get wrong:

- **"Some open-shell systems need complex" is false for the linear problem.**  A real symmetric H
  always admits real eigenvectors.  Complex can only help the NONLINEAR problem: a real HF/KS
  stationary point unstable toward a lower complex one (the classic complex RHF/UHF instability).
  That is an explicit ansatz policy, and it can only ever DOWNGRADE a block real→complex.
- **SOC is not "T becomes complex".**  It merges the spin channels into 2-component spinors, so the
  `ms` half of `Irrep` stops being a good quantum number.  A scalar-type rule cannot express that;
  it is a block-STRUCTURE change.  Anything designed here must not hard-wire "2 independent spin
  channels" either.

### The rule

```
block is real  ⇔  irrep.IsReal() ∧ (every term.PreservesReal())    [∧ not overridden by ansatz policy]
```

Each participant answers about ITSELF, which is what makes it derivable rather than declared:
`BlochQN::IsReal()` from \f$2\,ik\equiv0\ (\mathrm{mod}\ N)\f$ — **exact integer arithmetic, no
tolerance**; `Spherical` / molecular point-group irreps trivially true; an SO or vector-potential
term always false; kinetic / Coulomb / XC / local-PP / KB always true.  Adding SOC to a run then
flips every block complex automatically, with no special-casing anywhere.

## 2. TWO types, not one — the distinction that must not be lost

```
basis type    = irrep.IsReal()                          // "are my functions (and S, T, V) real?"
working type  = basis.IsReal() ∧ ham.PreservesReal()    // "are H, C, D real?"        working <= basis
```

With SOC the basis functions are still real — real Gaussians, real SALCs — so S, T, V_loc, the KB
projectors and **the collocation streams** stay real; only H, C, D go complex.  Since the geometry
layer is what dominates GPW, conflating the two types would needlessly double the cost of the
expensive half the moment SOC arrives.

**A THIRD axis, found by V1.1(i) (2026-08-16) — the FIT basis.**  `Orbital_DFT_IBS` is now
`template <class T, class TFit = T>`, because the fit basis is a property of the RUN, not of a block:
every call site builds it once from the whole `tBasisSet` and all irrep/k blocks share it.  So when
`BlochQN::IsReal()` makes TRIM blocks real, those real blocks still share the run's COMPLEX G-space fit
basis — `Orbital_DFT_IBS<double,dcmplx>`, a combination that had no spelling before.  Three of four
combos are meaningful (`<double,double>`, `<dcmplx,dcmplx>`, `<double,dcmplx>`); `<dcmplx,double>` is
never instantiated.  Note this is NOT the same split as basis-vs-working type below — it is orbital-vs-fit
— so the full picture is three axes, and any two of them can differ.

**Today they are welded:** `tSCFIterator<T>` takes `tbs_t<T>*` — one T for basis and everything
downstream.  For the Bloch-only case that is harmless (basis real ⇒ working real), so the first
increment need not decouple them.  **The standing instruction is only: do not harden it further.**
When the `dcmplx`-pinned solid classes are touched, keep basis type and working type as separate
parameters even while they are always equal.

## 3. Where the decision is computed

**Not in the basis, and it cannot be:** `qcHamiltonian` links `qcBasisSet`, and every driver builds
basis → Hamiltonian → iterator (`Hamiltonian::Factory(..., bs.get(), ...)` takes the basis).  A
basis that consulted H would be a dependency cycle the linker rejects.

- basis / irrep answers about itself → `Symmetry::IsReal()` (virtual, default true; `BlochQN`
  overrides), forwarded by `IrrepBasisSet::IsReal()`.
- each term answers about itself → `PreservesReal()`, constant per term type; `tHamiltonian<T>`
  folds its terms with ∧ (a one-liner).
- **the composition root computes the AND** — the layer that already constructs both and then builds
  the WF, CD and accelerator (`RunGpw` / `Hamiltonian::Factory` / `qchem::Calculation`).  It is the
  only place with both facts, and it is already the assembly point.

## 4. The containers that inherit it

One decision, carried by the type system, in dependency order:

1. `BasisSetImp<T>` → its `IrrepBasisSet` children — **the decision point**
2. Hamiltonian `tHT_Common<T>::CacheMap` (`std::map<Irrep,hmat_t<T>>`) → per-irrep H blocks
3. `CompositeWF` → `tIrrepWF` children (`itsIWFs`/`itsQNWFs`) — where the real eigensolve lives
4. `tComposite_CD` → `tDM_CD` children (`itsCDs`)
5. **`tSCFAccelerator<T>` — the open one.**  See §6.

**Child slot** = the element type of the container: today `tDM_CD<T>*` / `tIrrepWF<T>*` (one T per
run, so all children necessarily agree); under this plan
`std::variant<IrrepCD<double>, IrrepCD<dcmplx>>` or equivalent, so children of one composite may
differ.  `std::visit` with a GENERIC LAMBDA keeps every algorithm single-source — that, not the
container, is the answer to the dual-maintenance concern.

**Keep the composites templated at their FACE.**  A genuinely non-template `tComposite_CD` would
still have to derive from `tDM_CD<dcmplx>` to serve `tHamiltonian<dcmplx>`, so molecules would need
a SECOND composite class — two copies of the aggregation logic, i.e. exactly the dual maintenance
being avoided.  Face stays `tComposite_CD<T>`; only the child slot varies.  Molecules then hold a
variant that never takes the complex alternative: a negligible tag, one implementation.

## 5. Staging

- **Step 0 — DONE `1b8b9a83`.**  Exact ±1 Bloch phases at TRIM, so realness is a BITWISE fact, not a
  tolerance.  This is what lets a later narrowing of H be an ASSERTED narrowing.  Gate
  `GPW.TRIM_BlochMatricesAreExactlyReal`.
- **Step 1 — DONE 2026-08-17.**  `Symmetry::IsReal()` (default true; `BlochQN` overrides with EXACT
  integer arithmetic in the ctor: real ⇔ \f$N_i\,|\,2(ik_i+\mathrm{shift}_i)\f$ per component, a
  non-half-integer shift fails exactness outright) + `IrrepBasisSet::IsReal()` (pure forward) +
  `PreservesReal()` on all three term families (default true) with `tHamiltonianImp` AND-folding it
  in `Add()` beside `IsVirialValid`.  Pure queries, no behaviour change, no types moved.  Both
  `basis.perIrrep` emitters (molecular `MakeIrrepWFs`, GPW `VetGpwConditioning`) report a `real`
  field per block.  The GPW evaluator's float-k `IsTRIM()` is DELETED: the ctor takes `kIsReal`
  beside `kFrac` (the same pried-out-fact seam as `Getk`) and `GPW_IBS` passes `irrep->IsReal()`.
  Unit test `SymQNTests.IsReal` (Γ / zone boundary / negative index / MP-shift); full sweep 760/760.
- **Step 2 — the composite child slot** (`tComposite_CD` first: fewer methods, ~8 of 18 are T-typed
  and all are aggregation points; `IrrepCD<T>` is already templated so both alternatives exist).
  Prerequisite increment: the §6 accelerator face change (scalar-agnostic manager, typed `Create`
  overloads) — **LANDED 2026-08-17, see §6**.
  Buys the WF/D half — diagonalization and C/D storage.  Dominant in the many-PW/ultrasoft regime;
  small in GPW.  Increments:
  - **2a DONE 2026-08-17** — `tComposite_CD`'s children are now
    `std::variant<unique_ptr<tDM_CD<double>>, unique_ptr<tDM_CD<dcmplx>>>` (the §4 child slot),
    with two typed `Insert` overloads.  Scalar-independent aggregation (charge, ρ(r), rescale,
    the Fourier trio — `FourierDensity` is a non-template face) is single-source generic lambdas;
    the T-typed-argument operations (contract clients, block/Φ maps, MixIn pairing, HF sweep)
    forward through a `SameT<T>` view whose cross-scalar arm THROWS until Step 3 makes a mixed
    child reachable.  Pure structure, no behaviour change; full sweep green.
  - **2b DONE 2026-08-17** — the `CompositeWF` child slot: `itsIWFs`/`itsQNWFs`/`itsSpinWFs` are
    variant slots over `tIrrepWF<double>`/`<dcmplx>` (owning + non-owning-ref mirrors).  The
    scalar-signature calls (GetIrrep, GetOrbitals — the `Orbitals` base is non-template,
    DoSCFIteration, ComputeStep, MoveOrbitals, basis populations) are single-source generic visits
    (`RefOf`/`IrrepOf`/`OrbitalsOf` helpers); the T-typed calls (CalculateH's Hamiltonian, the
    fills' OccupationPolicy, the MOM OrbitalView cast) forward through `SameT<T>` with the same
    throwing cross arm as 2a.  `GetChargeDensity` is the first MIXED-READY seam: each child builds
    a density typed by its own scalar and 2a's typed `Insert` overloads take either.  Pure
    structure, no behaviour change; full sweep green.
  - **2c** — the NARROWING that makes a real child real: when `irrep.IsReal() ∧ ham.PreservesReal()`
    (the fact must be threaded into the WF build — `WaveFunction::Factory` has the ham), the block's
    complex H at TRIM is ASSERT-narrowed (bitwise, by Step 0) to real, solved with `LASolver<double>`,
    and C/D stored real — its `GetChargeDensity()` then Inserts the `<double>` alternative of 2a.
    Acceptance: a Γ-only GPW run bit-comparing the narrowed path against the complex one.
- **Step 3 — un-pin the basis** (`GPW_IBS` from `IrrepBasisSet<dcmplx>`), so terms produce real H
  directly and Φ + the quadrature GEMMs go real.  Buys the Hamiltonian/quadrature half — dominant in
  the GPW regime.  The composite variant then falls out as a CONSEQUENCE of blocks having different
  basis types, which is the coherent order.  **ORDER DECIDED 2026-08-17 (user): Step 3 before 2c** —
  no throwaway narrowing-child class; the real WF child falls out of the real basis.  Increments:
  - **3a DONE 2026-08-17** — `tGPW_IBS<T>` exists with BOTH instantiations (`using GPW_IBS =
    tGPW_IBS<dcmplx>` keeps every spelling).  The EPW mixins are `<E, T=dcmplx>`
    (PlaneWave/fit bases untouched); a `<double>` block is `IrrepBasisSet<double>` +
    `Orbital_DFT_IBS<double,dcmplx>` (the V1.1 combination, now explicitly instantiated) +
    `Integrals_Pseudo<double>`, with the evaluator's exactly-real results ASSERT-narrowed by
    `ToScalar` (bitwise, per Step 0) at the mixin/PP seams — and cached REAL in `theCache<double>`,
    so a real block never stores widened matrices.  The DFT tensor tier needs no narrowing at all
    (tensors follow TFit).  The evaluator itself stays complex — making its streams real is a later,
    purely performance increment.  Gate: `GPW.TRIM_RealBlockMatchesComplexBitwise` (S/T/Ven/V_NL/
    V_loc + orbital values, `EXPECT_EQ` on doubles, Γ and the zone corner).
  - **3b DONE 2026-08-17** — the `BasisSetImp` child slot (§4 item 1: the DECISION POINT):
    children are `std::variant<unique_ptr<Orbital_1E_IBS<double>>, unique_ptr<Orbital_1E_IBS<dcmplx>>>`
    with typed `Insert` overloads; the scalar-independent whole-set faces (`GetNumFunctions`,
    `GetIrreps`, `Write`) are generic visits; `GetIBS` (the face `Iterate<D>` casts from) is the
    SAME-scalar view, throwing loudly on a cross child; `GetChild(i)` is the typed view 3c's
    mixed-aware consumers use.  NO factory wiring yet (deliberate: enabling
    `IsReal() ? tGPW_IBS<double> : ...` in `GPW_BasisSet` before 3c would break every
    `Iterate<Complex_OIBS>` consumer) — the decision wires in WITH 3c.  Unit tests
    `RealComplexBasisSlot.*` build the first genuinely MIXED set (real Γ + complex ¼-k in one
    `BasisSetImp<dcmplx>`) and pin the aggregate/typed-view/throw contract.
  - **3c — re-thread the consumers per block.**  Increments:
    - **3c-1 DONE 2026-08-17 — the TERM STACK serves the real block.**  New capability faces
      `Static_HT_RealBlock` / `Dynamic_HT_RealBlock` (the V1.6/V1.7 cross-cast idiom: only terms
      that can serve a real block carry them; the assembly cross-casts and fails loudly otherwise),
      with caching Imp mixins mirroring `tStatic/tDynamic_HT_Imp`'s Irrep-keyed discipline over an
      own real cache (the §4-item-2 container, realized as a typed pair — disjoint irrep key sets).
      ALL ten periodic terms implement them via ONE scalar-generic assembly body each
      (`MakeMatrixT<U>`: cross-cast to `Integrals_Pseudo<U>` / `Orbital_DFT_IBS<U,dcmplx>`, narrow
      exactly at the end): Kinetic and IonIon (conditional `StaticRealBlockBase<T>` on the dcmplx
      instantiation), the three PP terms, Vee_Hartree, PWFittedVxc (raw route; the legacy ball-fit
      route throws — nothing real-block reaches it), and the three Becke terms.  `XC_GridEngine`
      gained the typed `Matrix(robs_t)` + a real Φ-table cache — a real block's quadrature GEMM runs
      blaze's REAL kernel, the first realized Step-3 win.  Gates `RealComplexTerms.*`:
      statics + Hartree BITWISE vs the complex block; the Becke GEMM machine-equal (~3e-15 — the
      real kernel's different summation order, i.e. the win itself; everything elementwise stays
      `EXPECT_EQ`).
    - **3c-2 DONE 2026-08-17 (Fock half) — the ASSEMBLY + SCF wiring.**
      `tBasisSet` face: public `GetNumIBS` + `GetRealIBS(i)` (the cross-scalar block view, default
      null; `BasisSetImp` serves the `<double>` alternative of a complex-faced set).  `Ham_RealBlock`
      face conditionally on `tHamiltonian<dcmplx>`; `tHamiltonianImp` folds the terms'
      Static/Dynamic_HT_RealBlock capabilities (HF terms throw — no periodic exact exchange).
      `tIrrepWF<double>` gained the run-typed `CalculateH(Ham_RealBlock&, cChargeDensity*, cbs_t*)`;
      `tCompositeWF`'s Fock-build loops dispatch per child (`CalcH<T>` — 2b's reserved cross arm,
      live), and `MakeIrrepWFs` is the MIXED walk (`GetRealIBS` first → `MakeOneIrrepWF<double>`,
      single-source member template; the real child's accelerator comes from the §6 typed `Create`
      and its density Inserts 2a's `<double>` arm).  The GPW preflight (`VetGpwConditioning`,
      `EmitGpwGrids`) walks the typed child slot generically.  Gate:
      `RealComplexTerms.HamiltonianAssemblyServesTheRealBlockBitwise` — the `Ham_RealBlock` fold
      (kinetic + 3 PP + Hartree) equals the native complex assembly BITWISE, i.e. exactly the matrix
      a real WF child receives.
      **Deferred to 3c-2b (the ENERGY/DENSITY half, prerequisite for 3c-3):** the composite
      density's cross-scalar contract arms (static reuses `tStatic_CC<double>` on the term mixins;
      dynamic needs a run-typed `GetEMatrixR` client face), the mixed `DM_RhoAtPoints` pointwise
      arm, and the REAL child's ρ̃ contribution — `IrrepCD`'s Fourier trio is conditioned on
      T==dcmplx, conflating scalar with lineage; Step 3 breaks that identification, so the
      conditional must move from the scalar to a basis-capability probe.
    - **3c-3** — wire the factory decision (`irrep->IsReal() ∧ PreservesReal`, threaded from the
      composition root into `GPW_BasisSet`); 2c falls out; then Step 4's mixed-mesh acceptance.
- **Step 4 — the accelerator** (§6), and only then a mixed-mesh run as the acceptance test.

An all-TRIM mesh (Γ-only, or Γ-centred 2×2×2 = `MNO_KMESH=2`) makes an entire run real with no
heterogeneity at all.  Mixed meshes (3×3×3, shifted MP) are the only case that truly needs per-block
variation — worth remembering when scoping.

## 6. Open questions

- **The accelerator — SETTLED 2026-08-17 (code survey of all four: DIIS, GDM, Ladder, Null).  The
  feared redesign does not exist: the history is ALREADY per block.**  Every block-typed object —
  the DIIS \f$F'\f$/\f$E=[F',D']\f$ deques, GDM's orbitals, Fock, CG tangents and geodesic factors —
  lives in the PER-IRREP accelerator (`tSCFIrrepAcceleratorDIIS<T>` / `...GDM<T>`), created one per
  block by `Create(lasb,qns,occ)` in `MakeIrrepWFs`.  The run-level `tSCFAccelerator<T>` holds NO
  T-typed state at all: DIIS's cross-block coupling is the REAL B matrix
  \f$B_{ij}=\sum_\mathrm{blocks}\mathrm{Re\,tr}(E_i^\dagger E_j)\f$ and REAL coefficients \f$c\f$
  (each child reads them back through a `const rvec_t&`); GDM's manager is explicitly "no global
  coupling"; the Ladder holds rungs and a rung index.  Every manager↔child interaction is scalar:
  `GetError(i,j)→double`, `GetNproj`, `Append1`/`Purge1` (the lockstep depth), `Ready`, `Engageable`,
  `RejectStep`.  So in a mixed-irrep run a real TRIM block gets a `<double>` irrep accelerator and
  its history is stored real — the memory win in the many-PW regime realizes AUTOMATICALLY, and no
  widening point exists anywhere, because `UseFD` (the only face that takes block-typed matrices)
  sits on the per-block object.
  **The residual work is the FACE, and it is small:** the run-level `T` appears only in `Create`'s
  signature and the child-pointer element type.  Decision: make the run-level accelerator
  SCALAR-AGNOSTIC — one non-template manager per algorithm with two typed `Create` overloads
  (`LASolver<double>*→tSCFIrrepAccelerator<double>*`, `dcmplx` likewise; both alternatives of the
  templated per-irrep classes already exist), holding its children through a small non-template
  aggregation sub-face carrying exactly the scalar methods listed above.  Real extrapolation
  coefficients over complex \f$F'\f$ are mathematically unchanged (Pulay's \f$c\f$ is real by
  construction of the metric), so ONE \f$c\f$ still serves all blocks.  This face change is a
  PREREQUISITE INCREMENT for Step 2 (mixed children need a `Create` per block type) and touches the
  `tSCFIterator<T>`/`tCompositeWF<T>` accelerator-pointer seams — do it as its own commit, keeping
  basis and working type separate parameters per the §2 standing instruction.  **Gate open for
  Step 2.**
  **The face increment LANDED 2026-08-17:** `SCFAccelerator` is the ONE non-template manager (typed
  `Create` overloads; DIIS/GDM coordinate their children through the scalar `DIIS_Block`/`GDM_Block`
  sub-faces, the Ladder's rungs are themselves scalar-agnostic, both Factory doors return the same
  type, and the friend-based manager↔child couplings are gone).  The `tSCFIterator<T>` /
  `tCompositeWF<T>` / WF-factory seams hold plain `SCFAccelerator*`.  Behavioural anchor: NEW
  `M_Calculation.BoronUHF{LadderGDM,PureGDM}` — boron's open-2p¹ doublet through the facade's
  Ladder and standalone-GDM options (the molecular GDM path's first gtest coverage; both land on
  the same UHF/dzvp energy −24.525932 to ~1e-12 Ha) — plus the full sweep.  **Step 2 may start.**
- **2-component structure.**  Non-collinear spin (likely to arrive before SOC-for-anisotropy, given
  the frustrated-magnet/battery direction) needs spinor blocks, not just complex scalars.  Nothing
  here should hard-wire "exactly two independent spin channels".
- **Ansatz policy.**  Where does "allow complex orbitals" (complex instability, GHF) live?  It is a
  run-level policy that downgrades blocks; probably `SCFParams`, but it interacts with the
  composition root's AND and has no home yet.

## 7. Prerequisites (doc/CleanupCandidates.md)

**STATE 2026-08-17: THE GATE IS OPEN — every prerequisite is done.**  V1.1 ✅, V1.5 ✅, V1.6 ✅
(`2d0f6982` + `bd6c648a`, which completed it at the composite/polarized level), V1.7 ✅, V1.8 ✅,
V1.10 ✅, V1.11 ✅ (`43bbebad`..`2398dd07`, five increments), R2.16 ✅ (done 2026-08-07; only its header
marker was missing).  Note V1.5 was never actually a GATE for this work — `G_FieldEvaluator` sits off
the T axis (see below) — but it landed anyway, and its four client-named faces make the re-typing
cheaper, so the distinction is now academic.

Load-bearing, because they sit on the faces this restructures — finishing them SHRINKS this work,
and none of them is invalidated by it (the type change alters child containers, not the method lists
the ISP splits produce):

- **V1.8** `IrrepCD`↔`IrrepCD` concrete same-class casts in the hot path — they target the very class
  that becomes a variant alternative.
- **V1.6** `tDM_CD::Accumulate*` face split — the composite's T-typed methods ARE the widening points.
- **V1.7** the periodic trio's 9 asserting overrides — each multiplies against a new alternative.
- **V1.11 ✅ DONE 2026-08-17 (`43bbebad`..`2398dd07`, five increments)** — the occupation seam is the
  SCFIterator's `OccupationPolicy` slot; the WF child carries NO occupation state, so the variant-child /
  narrowing work lands on a clean face.  **ALL §7 prerequisites are now DONE.**
- **V1.1 ✅ DONE 2026-08-16 (`d49db261`+`4b37221d`+`2bfb83b4`)** — the basis-face merge landed:
  `Projector3<T>` unifies the 3C tensor, `Band_FT_IBS` is deleted, and `Orbital_DFT_IBS<double,dcmplx>`
  (a real TRIM block on the run's complex fit basis — exactly this plan's case) is a live spelling whose
  tensor caches under `theCache<TFit>()`.  **V1.5 ✅ DONE 2026-08-16 (`f18a6ee9`+`9ebaebdb`)** — the
  grid-engine union is four client-named faces (`G_FieldEvaluator`/`G_Quadrature`/`G_StructureFactor`/
  `G_SpectralFilter`), so re-typing touches one face per consumer.  **V1.10 ✅.**  The basis was the
  decision point and the head of the type flow — **all basis-side prerequisites are now DONE; only
  V1.11 (the occupation seam) remains.**
- **R2.16** construction-time facts re-asked at run time — the same principle this rests on.

NOT prerequisites (do not gate on them): R1.0b (shared-radial reader), R2.20 (oracle helpers in a
test module), R2.5's remaining `exit(-1)` sites.

**V1.1's question — *"can `Band_FT_IBS` be `Orbital_DFT_IBS<dcmplx,dcmplx>`?"* — was ANSWERED on
2026-08-16, and the answer was better than the question:** `Band_FT_IBS` is DELETED, `Projector3<T>`
unifies the 3C tensor, and `Orbital_DFT_IBS<double,dcmplx>` is a live spelling.  That combination — a
real TRIM block sharing the run's complex fit basis — is exactly this plan's case, so Step 3's basis
un-pinning now starts from a merged DFT face rather than two classes.
