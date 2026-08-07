# Cleanup candidates (running list for the SOLID/OOD/CleanCode session)

Things noticed in passing while adding features — flagged here instead of fixed inline, so the
refactoring session can batch them.  (User keeps the master list; merge freely.)

Reorganized 2026-08-05 (Claude review pass): everything below the SOLID section is grouped by
status — **READY / VERIFY / DECIDED-ELSEWHERE / DONE** — with per-item verdicts from a 4-sweep
code verification.  User prose retained (typos fixed only); user replies of 2026-08-05 folded in.

## High-level goals (user)

- I am mostly concerned about abstract interfaces.  These are what the client code for any module
  is supposed to see.  This project is an attempt (probably never been done before) to capture one
  set of high-level interfaces that support SCF electronic structure calculations for any
  structure: atoms, molecules and 3D lattices (also 1D polymers, 2D graphene), using any basis
  set.  In other words 95% of the high-level abstract interfaces (charge density, Hamiltonian,
  orbitals, wavefunction, SCF iterator, accelerators) should be structure neutral and basis-set
  neutral.  Even the basis set interfaces (qcBasisSet) are (were) structure neutral.
- I have so far identified 4 libraries that have a structure-neutral face with specific structure
  specializations:
  1. qcStructure: Structure has derived classes for Atom, Molecule, UnitCell.  Lattice_3D and
     ReciprocalLattice have no inheritance relation with Structure.
  2. qcSymmetry: separate folders for Atom/Molecule/Lattice_3D symmetry types.  Spin degrees of
     freedom are orthogonal to spatial symmetry ... but that only works for non-relativistic irreps.
  3. ElectronConfigurations: we need separate classes for non-aufbau (specific # of electrons per
     irrep) filling.  We also have separate classes for Atom, Molecule and Crystal aufbau filling
     ... can these be combined into one general aufbau filler?
  4. qcBasisSet -> qcAtom_BS, qcMolecule_BS, qcLattice_3D_BS
- With the introduction of lattice calculations some structure-specific classes have been creeping
  into the qcChargeDensity and qcHamiltonian libraries.  Also possibly into qcWaveFunction and
  qcOrbitals.  If we can find ways of getting these refactored back to more generic status, that
  would be a big win.

## SOLID principles

1. SRP: Single Responsibility Principle
   - This principle is more about concrete classes than it is about interfaces.
   - We have a group of classes called Evaluators that do the low-level work of evaluating all
     required integrals, op(r)/grad(r) and a couple of other minor things, for a particular basis
     function type.  These have multiple responsibilities but they are mostly built up through a
     network of mixins in order to achieve the end result.  The objective with the Evaluators is to
     have a well-defined list of functions you need to implement in order to have the framework
     just work for an SCF calculation.  So we are not trying to satisfy SOLID at this level of code.
   - Mixins are a very powerful technique for building up a final concrete class.  We can apply the
     SRP to the mixin components, but obviously it does not make sense for the final class.  It is
     by definition multi-responsibility.
2. OCP: Open Closed Principle
   - We are still in the R&D stages for this project so abstract interfaces will be modified.  We
     just need a good reason to do so.
   - For example if we extend a structure-neutral interface in order to get some Lattice_3D feature
     working, we need to ask the question: what does this change mean for atoms and molecules?
3. LSP: Liskov Substitution Principle
   - Virtual functions that default to some sort of "not implemented" behaviour violate the LSP.
4. ISP: Interface Segregation Principle
   - There is always an urge to add getters and setters.  For setters, always ask: is this
     something I can set at construction time?  For getters always ask: what are we going to do
     with the getter data?  Why not ask the owning class to do that task instead (maintain
     encapsulation).
   - I do not claim that all the abstract interfaces in qchem are truly segregated.  A good example
     is the irrep basis set interfaces.  They do three things:
     a. Deliver integral tables (interface mixin from multiple).
     b. Evaluate op(r) and grad(r) (interface mixin from VectorFunction<T>).
     c. Expose the irrep symmetry, as an abstract interface pointer.
   - a & b are built up from many mixin interfaces.  ISP?  I think so.
   - Any time our abstract interfaces get augmented, those new functions should be part of an
     existing responsibility, not introducing a new one.
5. DIP: Dependency Inversion Principle
   - This is possibly the most powerful concept, and often enables adherence to the previous 4
     principles.
   - This applies equally well for C++ classes, C++ modules (DAG enforced by compiler) and
     libraries (DAG enforced by linker).
   - Example: lib A depends on lib B.  B needs a way to send info to A.  Create an abstract
     interface AI (ha ha) in B; classes in A derive from B::AI; when instances of the A classes are
     passed to B (as an AI*), B can then call back into A.
   - src/ChargeDensity/ChargeDensity.C tStatic_CC, tDynamic_CC are a working example that the
     Hamiltonian library classes derive from and pass back into the ChargeDensity library.
   - As well as being a dependency inversion, this probably has a "Gang of Four" pattern name.
     (Claude: the closest GoF name is **Observer** — a callback interface owned by the lower layer;
     in modern terms simply a "callback interface".)

There are 6 other principles related to package cohesion and coupling.  My experience is that the
first 5 are the most commonly misunderstood and misapplied.

---

# Worklist (status-organized, 2026-08-05)

Legend:
- **READY** — claim verified against the tree; fix direction clear; executable without a design
  session.  (Correctness-adjacent items front-loaded.)
- **VERIFY** — needs a design discussion, a measurement, or a repro before acting.
- **DECIDED-ELSEWHERE** — the governing decision/record lives in another doc or pin; engage it
  there rather than relitigating here.
- **DONE** — kept for the record.

Per-item status markers (added 2026-08-06), so status is visible where you READ the item rather than
only in the LANDED section below:
- **✅ DONE `<sha>`** — landed and green on branch `solid-cleanup`.
- **⏳** — written and building, but not yet fully verified.
- **🔶** — the DESIGN is settled (usually by a user ruling) but no code has been written yet.
- unmarked — untouched.

Execution venue (user, 2026-08-05): a branch in the second working tree **~/Code/qchem1**, so the
MnO campaign proceeds undisturbed in qchem6.

## LANDED on branch `solid-cleanup` (qchem1, 2026-08-05) — 665/665 ctest green

- `06e23f5d` — **R1.1, R1.2, R1.6, R1.8, R2.1.**  The `te.Exc=0.0` clobber deleted; `=`→`+=` sweep
  (Kinetic, DiracKinetic, RestMass, Ven, PW_Pseudo Een+E_alphaZ, PW_Kinetic, Vee); FittedVee
  re-expressed through locals (its `te.Eee` is a Dunlap combination of its OWN two pieces, so a
  bare `+=` would have read back another term's accumulation); Vee's two dead zero-stores deleted;
  `DM_ContractBlocks` pure virtual + its false "periodic path asserts out" doc fixed; two `Write()`
  pointer-streams dereferenced; FittedVee's missing cast assert added.
  **REFINEMENT TO R1.2 (found while doing it):** `GridChargeLost` must KEEP `=` — Delta_XC_Pol and
  Delta_VcorrPol both write the same value, so `+=` would double a health diagnostic.  A blanket
  sweep breaks it.  That it needs the OPPOSITE rule from every neighbouring field is independent
  evidence for V1.12 (get the diagnostic out of the energy value object).
- `80fc2ae8` — **V1.4** (Phi key → `map<Irrep,...>`, both sides on `Spin::None`).  Verified before
  landing: `Irrep::operator<` compares `SequenceIndex()` by VALUE (not handle identity) and
  `BlochQN::SequenceIndex()` uniquely encodes the k index, so periodic k-blocks keep separate
  tables instead of collapsing into one — the one way this change could have silently mixed one
  block's basis table into another's density.  All GPW/PW anchors unmoved.

- **R2.2** — `Kinetic` + `PW_Kinetic` collapsed to `Kinetic<T>` (new module
  `qchem.Hamiltonian.Internal.Kinetic`, inline, following the `IonIon<T>` recipe); `Internal/Imp/
  Kinetic.C` deleted, `PW_Kinetic` deleted from PWTerms.  Confirmed while doing it that PW_Kinetic's
  stated justification was FALSE: `Integrals_Kinetic<dcmplx>` IS instantiated, and the file that
  instantiates it says in so many words that the periodic dcmplx bases use the cached accessors.  So
  the periodic path now gets the `<p^2>` CACHE for free (it was calling the uncached `MakeKinetic()`
  every build).  Same value either way — the cache key carries k/Ecut/nG — and all PW/GPW anchors are
  unmoved.
  **GOTCHA WORTH KNOWING (now in CLAUDE.md):** moving the body into a template broke the build on
  `0.5*bs->Kinetic()` — *"operator* neither visible in the template definition nor found by ADL"*.
  A module that DEFINES a template using Blaze operators must `import qchem.Blaze` ITSELF; ADL at
  instantiation consults the template's DEFINITION context, not the instantiating TU's imports.
  Verified empirically: the import alone is sufficient and NO `<blaze/Math.h>` include is needed
  (both variants were built and compared).  Non-template code in an `Imp/` TU never hits this, so it
  will recur on every future T-templating collapse (`DiracKinetic`, `RestMass` if they go periodic).

- **R2.4 + V1.9 + R1.3** — the `Structure`→`UnitCell` cast sites are GONE, replaced by the user's
  free pry-out helpers (the `Symmetry::Atom::Getl` pattern), added to `qchem.ReciprocalLattice`
  (which already re-exports `qchem.UnitCell`, and sits in qcStructure below every consumer):
  `isPeriodicCell` (non-throwing probe), `GetUnitCell`, `GetReciprocalCell`, `GetReciprocalLattice`
  (the `Get` forms throw `std::bad_cast` via the reference cast).  All four sites converted:
  `SeedCD`'s anon `ReciprocalOf` DELETED (it was this helper, privately), SeedCD's flip-group cast,
  `MakeDensityMixer`, and `SCFIterator`.
  - **Design point that fell out:** the mixer factory needs a GRACEFUL fallback (it degrades to
    linear D-mixing and warns), so a purely throwing pry-out could not serve it — hence the
    `isPeriodicCell` probe alongside the throwing accessors.  Worth keeping in mind for the other
    capability faces: "can you?" and "give me it" are two different questions.
  - **R1.3 (the slicing copy) is FIXED as a side effect, not stopgapped**: SCFIterator now does
    `st->Clone()`, which is polymorphic (no slice of a derived cell) and returns the
    `shared_ptr<Structure>` the member already held.  V1.10b can still delete the member entirely.
  - Small ISP win: `MakeDensityMixer` was handing the concrete `UnitCell*` to
    `CreateVxcFitBasisSet`, which takes a `Structure*` — it now passes the neutral pointer.
  - Stale-comment batch: `Band_DFT_IBS.C`'s header no longer claims `PlaneWave_IBS` implements it
    (it documents the D1 re-decision instead), `PlaneWaveDFTUT`'s comment now names the real cast
    (`Integrals_Pseudo<dcmplx>`), and FOUR unused imports left SCFIterator (FourierDensity,
    FourierMixCD, Band_FT_IBS, ReciprocalLattice — each had exactly one occurrence in the file:
    its own import line).

- **V1.10b** — mixer creation moved OFF the structure-neutral iterator.  `MakeDensityMixer` (one factory
  running a three-way capability probe with a `cerr` fallback to linear D-mixing) is split into
  `MakeLinearMixer<T>` and `MakePeriodicMixer`; the choice is now a protected virtual
  `tSCFIterator<T>::CreateMixer` (base = linear) overridden by `SolidSCFIterator`.  The class that KNOWS
  it is periodic chooses, so the periodic factory takes Band_FT_IBS + UnitCell + FourierDensity as
  PRECONDITIONS and THROWS — user ruling 2026-08-06: "silent problems always end up consuming a lot of
  time... in the R&D phase failing loudly is encouraged" (a whole multi-hour run could previously
  converge on the wrong mixer behind one `cerr` line).  Everything the override needs is passed in, so
  iterator state stays private.  Verified first that all six Kerker/Pulay call sites are genuine GPW runs
  and that every `tSCFIterator<dcmplx>` in the tree IS a `SolidSCFIterator` (the lone `cSCFIterator` is a
  base-pointer holder), so the override is always reached.
  - **The `isPeriodicCell` prediction was HALF right (worth recording honestly).**  The probe did NOT
    disappear; it survives in exactly two places, and both are now CONTRACT CHECKS rather than
    behavioural branches: the `MakePeriodicMixer` precondition throw, and an `assert` in the base ctor.
    The base ctor keeps the cell snapshot because the `Structure` dangles by `Iterate` time and
    `SolidSCFIterator` inherits its constructors — but the runtime branch there became an assert, since
    the `if constexpr (is_same_v<T,dcmplx>)` above it already answers periodic-vs-molecular at compile
    time.  That doubled question (compile-time gate + runtime re-check) was the actual smell.
  - **Generalizable rule from this:** a periodic-vs-molecular `if` is a decision at the wrong altitude;
    the same probe reads fine as a *precondition*.  "Can you?" and "give me it" stay separate questions
    (V1.9), but "which am I?" should be answered by the type, not re-asked at run time.

- `72fecf8d` — **R1.4 + R1.5 + V1.3 (mechanism)** — three items, one session (2026-08-07); 667/667 ctest green.
  - **R1.4** the two silent-zero `Gradient()` overrides (FourierMixCD, `IrrepCD<dcmplx>`) now THROW.
    The plan said "assert"; asserts are compiled OUT of the build we test (`build/Release` is `-DNDEBUG`
    unless `QCHEM_RELCHECKED=ON`), so an asserting stub would still have returned the silent zero under
    `ctest`.  Throwing is loud in both configurations — the V1.10b ruling applied one level down.
  - **R1.5** `tChargeDensity::EvalBatch` DELETED; the fast batch is now an override of the inherited
    `ScalarFunction::operator()(rvec3vec_t)`.  The fork turned out to be LATENT (the one density-batching
    caller happened to use the `EvalBatch` spelling), not the live perf trap the item claimed.
  - **V1.3** both ε-adapter classes (`FittedEpsXc`, `FittedEpsCPol`) DELETED — but by naming
    `tDynamic_CC`'s method `GetEMatrix`, not by the planned `DM_ContractBlocks` reuse, which turned out
    to be blocked (a `tDynamic_HT` never receives `wholeBasis`, so it cannot enumerate the irrep blocks a
    block map needs — only `tDynamic_HF_HT` can).  The user's own V/E face IS the mechanism: one spelling
    was the entire collision.  Details + the blocked-route writeup under V1.3.
  - **Generalizable rule from R1.4:** "interim: assert" is only a real diagnostic where asserts are LIVE.
    Check the build's NDEBUG status before choosing assert-vs-throw for a *wrong-value* (as opposed to
    crash-soon) failure mode; the same NDEBUG hazard is already on record under V1.6.

**Process note for the next session:** build **`allTests`**, not just `ITMain`.  A first `ctest` in
qchem1 reported 160/590 failures that were ENTIRELY stale per-library `UT*` binaries (undefined
symbols + segfaults from the tree jumping ~30 commits while only ITMain was relinked) — zero real
failures.  After `ninja allTests` the same tree is 665/665.

**R1.3 is now SUPERSEDED-IF-V1.10b-LANDS and was deliberately NOT done** — the mixer-layering fix
deletes `itsKerkerCell` outright, so the one-line `Clone()` stopgap is wasted work if V1.10b lands
in the same session.

## READY

### R1 — correctness-adjacent (do these first)

- **R1.0b SHARED-RADIAL (SP/"L"-shell) support in the Gaussian94 reader — USER IDEA 2026-08-06.**
  A Cartesian d shell carries an s-type contaminant r²e^{−αr²}, which is what made the Mn GPW basis
  rank-deficient (λmin 1.15e-07, cond 8.2e7 — see `GPW_SCF.MnAtomInBoxDChannel`).  The cure that landed
  is empirical (drop the s window to 2 functions and let the contaminants span the rest).  The user's
  structural alternative: generate bases at the **valgen** stage with the SAME RADIAL SET SHARED ACROSS
  l — the standard trick Pople calls an **"SP shell"** (a.k.a. **"L shell"**, e.g. 6-31G), generalized
  here to s+d.  CP2K's format encodes it natively (one set with `lmin..lmax`, e.g. `2 0 1 4 4 4`) and the
  MOLOPT families share one exponent set across s, p AND d.
  **BLOCKED ON:** the flagged `PG_Cart::IrrepBasisSet` bug — the Gaussian94 reader MERGES same-exponent
  shells across l, which is why valgen carries the standing "KEEP EXPONENTS DISJOINT ACROSS l" rule.
  Multi-l shells must be representable before shared radials can be emitted.
  **MEASURE, DON'T ASSUME:** sharing exponents does not by itself remove the contaminant redundancy —
  the contaminant is r²e^{−αr²} (n=2) while an s at the same α is e^{−αr²} (n=0), i.e. independent
  functions; the measured cause was the NUMBER of s functions whose span mimics r²e^{−αr²}.  Shared
  radials buy compactness + CP2K-comparable structure; whether they also buy conditioning is an
  experiment (`GPW_SCF.MnAtomInBoxDChannel` + the vet's λmin/cond readout is the instrument).
  **ALSO WORTH FIXING WHILE THERE:** our GPW path is Cartesian-only — SPHERICAL bases throw
  ("not a molecular Gaussian basis (no `Molecule::LatticeSum1E`)", the parked S3b work).  Spherical d
  has no contaminant at all and is what CP2K actually solves in (its log for our own shell list: 55
  Cartesian vs 47 spherical functions), so wiring the spherical lineage into the lattice sums would
  retire this whole class of problem AND make the CP2K comparisons apples-to-apples.

- **R1.0 The "FAKE RADIAL" `op(r)` on atomic bases — USER DIRECTION 2026-08-06, own session.**
  An ATOM irrep block's 3-D face `operator()(rvec3_t)` returns the purely RADIAL chi_i(|r|) with the
  irrep's Y_lm silently omitted.  Any consumer that quadratures it as a real 3-D function gets nonsense
  when the integrand carries angular structure — this produced the occupied-d KB defect (l=0 projector
  leaking into EVERY l block; every l>=1 projector integrating to ~1e-33), invisible for the whole life
  of the PP code because MnO is the first system with OCCUPIED d projectors (fix c2d86ec9;
  doc/SymmetryUpgradePlan.md §7 step 7).
  *Interim fix that landed:* `BasisSet::ImplicitAngular_IBS` (`ImplicitL()`/`RadialValues(r)`) makes the
  fakeness explicit in the type system and `PP_NonLocal` cross-casts to it for a per-l radial assembly —
  that CONTAINS the trap, it does not remove it.
  *The real fix (user), two steps:* **(1)** consumers stop touching `op(r)` directly — code that only
  wants an INTEGRAL asks the basis for it, so the basis owns how its own functions are represented (the
  CLAUDE.md "prefer classes to do/answer high-level operations" bias; kills this whole bug class).
  NAMING (user, 2026-08-06): NOT a new `Integrate` verb — these are **`MakeOverlap` / `MakeOverlap3C`
  OVERLOADS**, joining the existing `Integrals_Overlap` family, with the arity following the call site:
    * \f$\langle\chi_i|f\rangle\f$ → a VECTOR → a `MakeOverlap(f, mesh)` overload
      (`PP_NonLocal`'s projection vector b_i).
    * \f$\langle\chi_i|f|\chi_j\rangle\f$ → a MATRIX → a `MakeOverlap3C(f, mesh)` overload
      (`PP_Local`'s V_loc block — f in the operator/C slot, exactly the existing 3C shape).
  Call sites on the raw-op(r) mesh route today: `PP_Local::CalculateMatrix`
  (`WeightedOverlap(mesh,*bs,VlocField)`), `PP_NonLocal::CalculateMatrix` explicit-angular branch
  (`Overlap(mesh,*bs,BetaYlmField)`), `Fit_IBS` (`Overlap(itsMesh,*this,f)`), plus the XC/fitting paths
  sharing the WeightedOverlap shape.
  **(2)** then make `op(r)` HONEST (real Y_lm linear combinations), nothing depending on the fake any
  more.  Design question to settle FIRST in that session: an atomic block carries l plus an m-LIST and
  the spherical solver keeps m-degeneracy in the OCCUPATIONS, so decide what `op(r)` returns per (i,m)
  before touching the density/XC paths, which currently rely on the radial-only form.
  *On completion:* RETIRE `ImplicitAngular_IBS` and collapse `PP_NonLocal::CalculateMatrixRadial` back
  into the single 3-D route.  Keep `A_PP.PerLKleinmanBylanderOracle` (s/p/d/f + cross-l zeros) and the
  CP2K atom oracles as the refactor invariant — they pin the physics independently of which route
  computes it.

- **R1.1 ✅ DONE `06e23f5d`. `FittedVxcPol::GetEnergy` clobbers `te.Exc`** — `te.Exc = 0.0;` before delegating
  (Imp/FittedVxcPol.C:83).  Correct today only because exchange precedes correlation in
  Ham_DFTcorr_P's term list; a silent order-dependent zeroing of any prior XC contributor.
  Fix: remove the zeroing (struct is zero-initialized), let the children `+=`.
- **R1.2 ✅ DONE `06e23f5d` (with one CORRECTION, below). `=` vs `+=` on `EnergyBreakdown`.**  Assigners:
  Kinetic, DiracKinetic, RestMass, Ven, Vee, IonIon, PW_Kinetic, PW_Pseudo (`te.Een=`);
  accumulators: PP_Local/NonLocal, Vxc, FittedVxc, PW_Hartree (`te.Een+=`), Delta_*.  PW_Pseudo(=)
  and PW_Hartree(+=) share `Een` in one Hamiltonian, correct only because the static list is
  iterated before the dynamic list (Imp/HamiltonianImp.C:79-81).  Sweep to `+=` everywhere.
- **R1.3 ✅ DONE `38a1ebd6` — fixed via `Clone()`, not stopgapped. `UnitCell` SLICING copy** — SCFIterator.C:163
  `make_shared<const UnitCell>(*cell)` slices any derived cell (e.g. FCCUnitCell) to a plain
  UnitCell; `Structure::Clone()` exists.  (The cast itself → V1.9.)
  **What the Kerker code does with the cell (user Q, answered 2026-08-05):** the snapshot is only
  passed to `MakeDensityMixer` (SCFIterator.C:244), which uses it for exactly two things
  (DensityMixer.C:327-328): (1) as the plain `Structure*` argument to `CreateVxcFitBasisSet` (no
  concrete-UnitCell need at all), and (2) `MakeReciprocalCell()` → the `ReciprocalLattice` for the
  q²/(q²+q₀²) Kerker filter.  So the entire concrete requirement collapses to "give me your
  reciprocal lattice" — exactly the V1.9 periodic-geometry face.  Slicing severity: BENIGN today
  (FCCUnitCell's ctor just forwards a special cell matrix — no extra state/behavior, so the sliced
  copy keeps itsA intact); latent, becomes live the moment UnitCell grows derived state.  The deep
  copy exists only because `Lattice_3D::GetStructure` returns a temporary that dangles by Iterate
  time (the in-code comment says so) — `Clone()` fixes the slice now; the lifetime fix upstream
  would remove the snapshot entirely.  **SUPERSEDED-IF-V1.10b-LANDS**: the mixer-layering fix
  deletes `itsKerkerCell` outright.  Do the one-line `Clone()` only as a stopgap if V1.10b is not
  in the same session.
- **R1.4 ✅ DONE `72fecf8d` (THROW, not assert — deviation explained). Silent zero `Gradient()` overrides** — FourierMixCD.C:75 and IrrepCD<dcmplx>
  (Internal/Imp/IrrepCD.C:333) return `rvec3_t(0,0,0)` with no assert; a future GGA/plotting
  consumer through the neutral ScalarFunction face gets silently wrong ∇ρ.  ~~Interim: assert.~~
  Real fix (a `DifferentiableField` capability face) → V1.6-adjacent design.
  **Both now `throw std::logic_error` naming the missing implementation.**  The planned `assert` would
  have been a NO-OP in the build we actually test: `build/Release` is `-DNDEBUG` (`QCHEM_RELCHECKED=OFF`,
  CMakeLists.txt:77-93), so an asserting stub still returns the silent zero under `ctest`.  A throw is loud
  in BOTH configurations, and matches the V1.10b "fail loudly in the R&D phase" ruling.  Verified no live
  caller: the density `Gradient` consumers are all `<double>` (SlaterExchange, the fitters, the composite/
  polarized/spin forwarders) — nothing on the periodic path asks for ∇ρ, and 665/665 stays green.
- **R1.5 ✅ DONE `72fecf8d`. `tChargeDensity::EvalBatch` duplicates `ScalarFunction::operator()(rvec3vec_t)`.**
  Identical signature after alias expansion.  `EvalBatch` DELETED (the inherited `ScalarFunction` batch
  op carries the same pointwise-loop default); FourierMixCD's fast factorized-phase path now overrides
  `operator()(rvec3vec_t)`; the four call sites in `XC_GridEngine` (PWTerms.C:375,403,404,411) spell it
  `(*cd)(mesh->Points())`.
  **CORRECTION to the claim above (verified while doing it): the fork was LATENT, not live.**
  OrthoFunctionFitter.C:101 batches `ps.GetScalarFunction()` — a `ProjectedScalar_R`'s field (v_xc, ε_xc),
  never a density — so no FourierMixCD reaches it today.  The ONLY site that batched a density was
  `XC_GridEngine`, and it happened to use the `EvalBatch` spelling, i.e. it got the fast path BY LUCK OF
  SPELLING.  The trap was one neutral-face caller away, not already sprung.
- **R1.6 ✅ DONE `06e23f5d`. `Write()` streams raw POINTERS (hex addresses)** — Imp/FittedVxc.C:111 (`os << itsLDAVxc`)
  and Imp/FittedVxcPol.C:93 (`os << itsUpVxc << itsDownVxc`); `qchem::op<<` binds only
  `const Streamable&`, so these hit `void const*`.  Need `*`.
- **R1.7 🔶 DESIGN SETTLED (user ruling), not yet implemented. `SymmetryAdapted_IBS::MakeDirect/
  MakeExchange` return empty `ERI4{}` silently**
  (SymmetryAdapted_IBS.C:80-81, comment "never called") — fails QUIETLY: a future generic-ERI4
  caller gets a zero Fock contribution.  Interim: assert.
  **USER RULING (2026-08-05): SymmetryAdapted_IBS should not be inheriting the (ERI4 part of the)
  `Orbital_HF_IBS<double>` interface — direction now settled.**  The interface's own doc comment
  already points there: "Client code rarely need ERIs directly, they only need the contraction
  over a density matrix."  Split `Orbital_HF_IBS` into (a) the CONTRACTION face (`Accumulate*` —
  what the HF terms consume) and (b) the ERI4-SUBSTRATE face (`MakeDirect/MakeExchange` +
  `Direct/Exchange` caches — the implementation detail behind the default contraction path).
  Concrete AO bases implement both; the SALC decorator implements only (a) (it builds in the AO
  basis and slices to its irrep block); terms consume only (a); the dummy `ERI4{}` bodies die.
- **R1.8 ✅ DONE `06e23f5d`. `FittedVee` casts `bs` and dereferences with NO assert** (Imp/FittedVee.C:41-42) — the
  sibling sites at least assert.  (The odftbs_t casts themselves are sanctioned abstract→abstract.)

### R2 — mechanical hygiene

- **R2.1 ✅ DONE `06e23f5d`. `tDM_CD::DM_ContractBlocks` → pure virtual.**  The asserting default is DEAD — all three
  concrete families override it (IrrepCD both T, Composite, Polarized); pure-virtual is free today.
  Also fix its false doxygen ("the periodic path asserts out" — the periodic path CALLS it,
  PWTerms.C:158).
- **R2.2 ✅ DONE `48e25b74`. Collapse `Kinetic` + `PW_Kinetic` → `Kinetic<T>`.**  Both are 0.5×(kinetic matrix);
  PW_Kinetic's stated justification (PWTerms.C:59-60, "symmetric cache bypassed for complex") is
  FALSE — `Integrals_Kinetic<dcmplx>` is instantiated (Imp/Orbital_1E_IBS.C:30).  Follow the
  `IonIon<T>` collapse recipe (Internal/IonIon.C documents it); PW_Kinetic dies.
- **R2.3 ⛔ WITHDRAWN — NOT a free dedup; re-filed as part of V1.1 (verified 2026-08-07).**
  ~~Dedup the literal 4-line Imp mirror — Imp/Band_FT_IBS.C:15-25 == Imp/Orbital_DFT_IBS.C:10-20 (only
  `theCache<T>` differs).  Mechanically removable independent of the V1.1 merge.~~
  The claim "only `theCache<T>` differs" is FALSE.  Read side by side, the two bodies differ in THREE ways:
  (i) `theCache<dcmplx>` vs `theCache<T>`; (ii) the RETURN type — `const G_ERI3&` vs `const ERI3<T>&`;
  (iii) the ARGUMENT types — `cFIT_SF_ABS`/`cFIT_CD_ABS` vs `rFIT_SF_ABS`/`rFIT_CD_ABS`.  (ii) and (iii)
  are precisely V1.1's two listed blockers (the `ERI3<T>` vs `G_ERI3` return-type question and the
  hard-wired-REAL fit-face arguments).  So there is no dedup to do here that is not the V1.1 merge itself —
  a shared body cannot be written until those two are settled.  **Do it inside V1.1, not before it.**
- **R2.4 ✅ DONE `38a1ebd6`. Stale-comment/import batch**: Band_DFT_IBS.C header claims PlaneWave_IBS implements it
  (false); PlaneWaveDFTUT.C:1079 claims PW_Pseudo routes through Band_DFT_IBS (it casts
  `Integrals_Pseudo<dcmplx>`, PWTerms.C:46); SCFIterator.C:30-31 imports of
  FourierDensity/FourierMixCD are stale (not consumers); PWTerms.C:59-60 false cache comment dies
  with R2.2.
- **R2.5 ⚗️ HAMILTONIAN HALF DONE (3 of 5 sites); the 2 ChargeDensity sites remain. `exit(-1)` in library code → throw** (5 sites): ~~Imp/LDAVxc.C:33-43 (dies with R2.6),
  Imp/FittedVxcPol.C:49-53, Imp/VxcPol.C:41~~, and `tPolarized_CD::MixIn`/`GetChangeFrom`
  (Imp/ChargeDensity.C:149-167 — also an LSP narrowing: accepts any tDM_CD, requires Polarized).
  Contrast the correct pattern at Imp/HF_HT.C:28-30.  These kill the pybind GUI / test runner.
  Fold into the D7 cast-survey custom-exception work.
  **DONE 2026-08-07:** LDAVxc's two died with the class (R2.6); FittedVxcPol/VxcPol now
  `throw std::runtime_error` naming the Spin::None-on-a-polarized-term mistake.  `qcHamiltonian` is now
  `exit()`-free.  **STILL OPEN: the two `tPolarized_CD` sites in qcChargeDensity** — left deliberately,
  because their fix is not just a throw: the LSP narrowing (the signature accepts any `tDM_CD` but the
  body requires a Polarized) is the actual defect, and that is V1.6/V1.8's seam.
- **R2.6 ✅ DONE 2026-08-07. The `LDAVxc` bundle** — a "Hamiltonian term" whose `CalcMatrix`/`GetEnergy` call
  `exit(-1)`; its only real job is `GetScalarFunction()` (the fit callback).  `FittedVxc` holds it
  as a raw OWNING pointer to the CONCRETE class (Terms.C:291, DIP violation);
  `tFittablePotential` (HamiltonianTerm.C:113-119) exists solely for it; `FittedVxc::
  UseChargeDensity` is dead (overrides nothing, no caller — `newCD()` already triggers the refit).
  One fix kills all four: collapse LDAVxc to a plain `Fitting::ProjectedScalar_R` adapter, hold as
  `unique_ptr<ProjectedScalar_R>`, delete `tFittablePotential`.
  **All four verified before the cut, and all four are gone.**  The module
  `qchem.Hamiltonian.Internal.LDAVxc` is DELETED outright (both TUs + the two CMake entries); LDAVxc had
  exactly one user in the whole tree, `FittedVxc`.
  - **Went one step further than the plan, for free.**  Rather than an adapter that hands the fitter the
    `ExFunctional` (what LDAVxc did), the adapter SELF-EVALUATES: `VxcDensity{ex,cd}` returning
    `ex->GetVxc((*cd)(r))`, an exact mirror of the `EpsXcDensity` already sitting beside it (and of
    FittedVcorrPol's `PolVcDensity`/`PolEpsCDensity` pair).  Verified byte-identical first: ALL THREE
    `ExFunctional` subclasses define `operator()(r)` as literally `GetVxc((*itsChargeDensity)(r))`
    (SlaterExchange, VWN_Correlation, Libxc_LDA).  The V and E fields are now visibly the same shape,
    differing in one method call.
  - **Consequence worth acting on: this makes V1.13 a DELETION rather than a redesign.**
    `ExFunctional::InsertChargeDensity` now has ZERO callers tree-wide (it had exactly one, in LDAVxc),
    so `itsChargeDensity` is never set and the whole FIELD face of `ExFunctional`
    (`operator()`, `Gradient`, the member) is dead code.  The fitter reaches the field through the
    adapters instead, which take the density as a ctor argument — the hidden-init landmine V1.13
    describes cannot happen on this path any more.  NOT deleted here: removing the field face means
    `ExFunctional` stops being a `ScalarFunction<double>`, which IS V1.13's value-face/field-face split,
    and it also retires `SlaterExchange::Gradient`.  Do it as V1.13, now cheap and compiler-verifiable.
  - Not touched (still V1.13): `SetPolarized`/`isPolarized` (still no caller, so `isPolarized` is still
    permanently true).
- **R2.7 `FittedCD::Clone()` — delete.**  Pure virtual (FittedCD.C:28) whose SOLE implementation
  asserts false and returns nullptr (Imp/FittedCDImp.C:58-64).  Dead contract clause; restore when
  the polarized-from-unpolarized path exists (real blocker per the assert message: a cloneable
  FunctionFitter).
- **R2.8 `InsertStandardTerms<dcmplx>` = assert(false)** (Imp/HamiltonianImp.C:49-53) — a
  molecular-only convenience on the T-generic base; move it down to the real lineage.
- **R2.9 Small Hamiltonian hardening**: (i) `XC_GridEngine` constness laundering —
  Rho/RhoPol/Matrix/Phi non-const, called from const term methods via a non-const `shared_ptr`;
  every other cache in the module is `mutable`+const-method — align it (and note: no
  cross-invalidation between its two rho caches, scalar + Up/Dn, which can be live simultaneously).
  (ii) `tDynamic_HT_Imp_NoCache::GetMatrix` returns a reference to shared scratch
  (HamiltonianTerm.C:104-107) — safe today only by immediate consumption; document or return by
  value.  (iii) `Dynamic_HF_HT_Imp::itsWholeBasis` latches on first use (Imp/HF_HT.C:31) — assert
  on change.
- **R2.10 `Fit_IBS::SetMesh` → ctor parameter.**  Two-phase construction; the construction-time
  principle was already settled for the XC quadrature — build the mesh first, hand it in.
- **R2.11 `DB_Cache_RAM.C`** — a screenful of `-Winconsistent-missing-override` warnings on every
  qcBasisSet build (`Get`/`Register`/`GetCache*`).  Mechanical `override` sweep.
- **R2.12 `UnmatchedCounts`/fold `tol` defaults** — 1e-8 fractional as a literal in three places
  (GPW `CreateXCQuadrature`, SymmetrizeMesh overloads); name it once.
- **R2.13 Becke strings/labels rename in `Delta_*`/`XC_GridEngine`.**  Verified: the classes are
  already behaviorally mesh-neutral — ZERO branching on Becke; only 3 `Write()` strings + 3
  profiler labels (the one factory branch lives in Imp/Hamiltonians.C:201 where policy belongs).
  Rename "becke" → "XC mesh"/"XC quadrature".  (The REAL non-neutrality — the
  `Symmetry::Lattice_3D::Fold` + dcmplx dependency — is V1.5's problem.)
- **R2.14 Hamiltonian term-naming sweep** (user conventions):
  - Enn/Vnn nuclear-nuclear repulsion; Een/Ven electron-nuclear attraction; Eee/Vee
    electron-electron repulsion; Eex/Vex exchange; Ecorr/Vcorr correlation; Exc/Vxc
    exchange-correlation.
  - For solids these don't work perfectly: at G=0 only certain combinations of Enn+Een+Eee are
    finite.  Also PPs are a combination of Ven+Vee + KB projectors; use Ven for PPs, where n is
    understood to be a shielded nucleus.
  - Renames wanted: `PP_Local`→`Ven_Local`, `PP_NonLocal`→`Ven_NonLocal`, `PW_Pseudo`→Ven-flavored
    name (user suggested Ven_PP), `PW_Hartree`→Vee-flavored name (user suggested Vee_PP) — final
    names = user's call; PW_Kinetic needs no rename (dies in R2.2).
  - **USER 2026-08-07: `Ven_PP` / `Vee_PP` both approved in principle.  Also asked whether PW_Pseudo's
    short/long components should be renamed too, assuming PW_Hartree is strictly a PP term and strictly
    ee repulsion.  BOTH ASSUMPTIONS ARE WRONG — verified against PWTerms.C before renaming anything:**
    - **`PW_Pseudo` has NO long-range component.**  Its matrix is `MakeLocalPotentialShort` + (optional)
      `MakeSeparablePotential` (Imp/PWTerms.C:50-53) and its energy is `te.Een` + the SHORT G=0
      alignment.  The LONG-range local part was moved OUT to `PW_Hartree` — the CP2K local-PP split
      (doc/GPWPlan.md 0e-PP): the deep-well erf potential folds into the ONE G-space Poisson solve
      instead of a per-orbital-pair sharp-field sweep.  So there is no short/long pair to rename here;
      the term is entirely short+nonlocal.  `Ven_PP` alone would read as "the whole PP", so prefer
      `Ven_PP_Short` (or `Ven_PP` + a doc line saying the long part lives in the electrostatics term).
    - **`PW_Hartree` is NOT strictly a PP term.**  `itsLocal` may be NULL — a pure all-electron / no-PP
      run gets plain Hartree with no core-charge fold (PWTerms.C:75-89 and the ctor doc).  A `_PP`
      suffix would be false for that configuration.
    - **`PW_Hartree` is NOT strictly ee repulsion.**  It contributes to BOTH energies (Imp/PWTerms.C:
      135-141): `te.Eee += 0.5*(total-eLong)` (electron-electron Hartree, WITH the ½ double-counting)
      and `te.Een += eLong` (electron-ion long-range, NO ½), plus the LONG G=0 alignment into
      `te.E_alphaZ`.  One Poisson solve, two physical energies — exactly the user's own "for solids
      these don't work perfectly" caveat made concrete.
    - **Naming consequence:** `Vee_PP` is misleading twice over.  Candidates that survive the facts:
      `Vee_VenLong` (says both energies, keeps the Vee/Ven vocabulary), or a what-it-IS name like
      `PW_Electrostatics` / `Poisson_PW` (one Poisson solve; the energy SPLIT is then a documented
      property rather than a promise the name makes).  Recommend `Vee_VenLong`, since the term
      genuinely owns two contributions and the name should not hide the second.  **Final call: user's.**
  - `Delta_XC` → e.g. `DeltaFittedVxc` (it IS a FittedVxc with a δ-function fit basis); rename
    along with the V2.1 decision.
- **R2.15 `nAngular` → degree-typed angular interface.**  `nAngular` is a COUNT for Lebedev but a
  DEGREE for GL/EM (and the imposed site-adapted builder consumes it as the degree) — the dual
  semantics BLOCKS flipping the free-run Becke default to the measured-equal Leb-302 (67% of
  GL-29's directions).  Fix: `angularDegree` + per-scheme count resolution; the default flip rides
  along.  (The warn/auto-resolve ergonomics half of this item stays in D5.)

## VERIFY

### V1 — interface-design questions

- **V1.1 `Orbital_DFT_IBS` ⇄ `Band_FT_IBS` merge.**  User (2026-08-05): Orbital_DFT_IBS simply
  specifies what integrals an IBS must supply to support DFT; "Band" and "FT" have no place in that
  specification.  History: Band_FT_IBS once had ~12 extra Fourier/Grid members, twiddled down over
  many sessions to essentially ONE — the main conceptual gap was reluctance to acknowledge that the
  PW representation of ρ is a fit, just a trivial one whose expansion coefficients come in one
  step.  Verified current state: (i) the scalar-type mismatch (Orbital_DFT_IBS<T> templated but
  fit-face args hard-wired REAL, rFIT_*_ABS in Fit_IBS.C:45-46,87-88, vs Band_FT_IBS hard-wired
  dcmplx) — user: trivial to fix, and the FIT_*_ABS<T> templating is already pinned with the fitter
  work; (ii) `ERI3<T>` vs `G_ERI3` return types (G_ERI3 was built as the harmonizing data-structure
  spec); (iii) extras: the using-decl is noise; `CreateXCQuadrature` (new, 2026-08 W1) is
  HOISTABLE to the neutral DFT face (molecules implement it returning their Becke quadrature — a
  unification, not a divergence); leaving exactly ONE genuinely Fourier member,
  `MakeOverlap(f(G))` (Band_FT_IBS.C:77) — itself re-expressible through the fit abstraction
  (fitted-potential coefficients → contraction), removing "FT" from the face.  Engage the pinned
  FACTOR-not-FUSE analysis (fitting-boundary pin + doc/FittingCleanupPlan.md) — the argument-type
  question is settled there; what remains is execution sequencing with the fitter templating.
  **Absorbed from the withdrawn R2.3 (2026-08-07):** the 4-line `Overlap3C`/`Repulsion3C` cache-lookup
  bodies in Imp/Band_FT_IBS.C:15-25 and Imp/Orbital_DFT_IBS.C:10-20 are the SAME code modulo exactly the
  two blockers listed above (return type `G_ERI3` vs `ERI3<T>`, argument `cFIT_*` vs `rFIT_*`) plus
  `theCache<dcmplx>` vs `theCache<T>`.  They are the smallest concrete instance of the merge, so they make
  a good FIRST target once those are decided — and a good litmus test that the decision actually works.
  **USER (2026-08-05): merge would be fantastic; THE one big remaining issue = for molecules we do
  a Dunlap fit for ρ (Coulomb metric, Repulsion integrals, charge-constrained) while for solids we
  don't (orthonormal projection ⇒ metric degenerate) — requires discussion.**  Groundwork for that
  discussion already exists in the fitting-boundary pin: the metric axis is REAL for AO (two
  different solves) and DEGENERATE for FT (projection IS the fit for both metrics), which is why
  the PW fitter implements both metric faces at once — the merge discussion is "how does the
  merged face express the metric choice without naming it", not "which metric wins".
- **V1.1b 🔶 ANALYSIS DONE, awaiting the user's re-read of the paper. The `Eee = 2·EeeFit − EeeFitFit` expression is DUNLAP-SPECIFIC — it is part of V1.1's
  metric discussion, not a free-standing formula (user, 2026-08-05; user wants to re-read the
  paper).**  Verified conventions: `GetSelfRepulsion()`=½⟨ρ̃|ρ̃⟩ (Imp/FittedCDImp.C:55) and
  `eeeFit`=½⟨ρ|ρ̃⟩, so the expression is exactly Dunlap's ROBUST form
  \f$E_J\approx\langle\rho\tilde\rho\rangle-\tfrac12\langle\tilde\rho\tilde\rho\rangle\f$, whose
  error is SECOND order in (ρ−ρ̃) **only because the fit is done in the COULOMB metric** — the
  stationarity condition ⟨ρ−ρ̃|c⟩_Coulomb=0 is what annihilates the first-order term.  The ctor
  says so out loud: "Charge-CONSTRAINED Coulomb-metric density fit (Dunlap-Connolly-Sabin 1979)"
  (Imp/FittedCDImp.C:21-23).
  - **Under an OVERLAP-metric fit the robustness justification evaporates** (stationarity is then
    ⟨ρ−ρ̃|c⟩_overlap=0, which does NOT kill the Coulomb first-order error) — so a "straight rho fit"
    needs a different energy expression.  User's instinct confirmed.
  - **Under an EXACT/orthonormal fit it degenerates CORRECTLY**: ρ̃=ρ ⇒ ⟨ρρ̃⟩=⟨ρ̃ρ̃⟩=⟨ρρ⟩ ⇒
    2A−B = ½⟨ρρ⟩ = the exact Hartree.  So the PW/GPW path is safe as-is; it is an overlap-metric
    GAUSSIAN fit that would break.  **This is exactly V1.1's "molecules Dunlap-fit, solids don't"
    issue seen from the ENERGY side** — whatever the merged face does about the metric, this
    expression moves with it.
  - **An overlap-metric fit path ALREADY EXISTS**: `NumericCD` (the SAD seed) overrides
    `GetUnconstrainedFit` with c₀=S⁻¹⟨f|ρ⟩ (NumericCD.C:48-51).  Today that is harmless and the
    invariant holds — but only BY ACCIDENT OF TYPING: `FittedVee::GetEnergy` and
    `tSCFIterator::TotalEnergy` take `tDM_CD*`, and the matrix-free seeds are `tChargeDensity` but
    NOT `tDM_CD`, so a seed can never reach the energy expression (it only feeds `CalcMatrix`,
    where the metric does not matter — J(ρ̃) is just the fitted density's potential).  Nothing
    states or enforces this; widen the energy path to `tChargeDensity`, or make a matrix-free
    density a `tDM_CD`, and the formula silently loses its second-order property.
  - Action: name the invariant where it lives (the fitter/energy pair), and fold the "which metric
    ⇒ which energy expression" rule into V1.16's explicit metric-strategy face.  Also ties the
    third leg of the CD taxonomy (D2's seeds) to a real correctness boundary.
- **V1.2 🔶 FEASIBILITY PROBE DONE + user APPROVED attempting it (appendix below). `Orbital_PP_IBS` — a structure-neutral PP-integral face (dependency INVERSION).**  User
  framing (2026-08-05): PPs require certain NEW TYPES of integrals from the IBS; the question is
  whether there is a structure-neutral way to ask for them without spilling PP details — if yes, we
  can break the qcBasisSet(qcLattice_BS)→qcPseudopotential dependence (today PlaneWave_IBS/GPW_IBS
  implement `Pseudopotential::Integrals_Pseudo<dcmplx>` whose args are PP types).  Candidate
  neutral primitives: (1) ⟨i|V|j⟩ for a species-attached local radial field (lattice-summed; the
  long/short split as a range parameter, not a PP concept); (2) projector brackets ⟨i|β_lm⟩ for
  species-attached radial×Y_lm functions.  qcPseudopotential then calls the face from ABOVE —
  the dependency edge inverts.  Passes the pseudo-wall pin (these ARE new integral types).
  Molecular PP (term-side quadrature today, by design) could optionally adopt the face later.
  Related dead surface: `Integrals_Pseudo::MakeLocalPotential` (unsplit matrix) is documented
  unit-test-only.  **USER (2026-08-05): approved to attempt — try it and see if we hit any
  roadblocks.**  (Feasibility probe of the actual `LocalPotential`/`SeparablePotential` payloads
  and the DAG: see the probe notes appended below when available.)
- **V1.3 ✅ MECHANISM DONE `72fecf8d` (both ε-adapters deleted; via `GetEMatrix`, NOT `DM_ContractBlocks` — see "WHAT
  LANDED INSTEAD" below).  The quadrature-term face (second list) is still open.  `FittedEpsXc`/`FittedVxc` simplification.**  Physics answer (user asked 2026-08-05):
  yes, genuinely different matrices — v_xc = δE_xc/δρ = ε_xc + ρ·∂ε_xc/∂ρ, so
  D·⟨i|v_xc|j⟩ = ∫v_xc·ρ ≠ ∫ε_xc·ρ = E_xc (Slater exchange: factor 4/3 — the retired ¾-virial
  shortcut, FittingCleanupPlan §I.1).  A class merge is impossible: both need `GetMatrix` with the
  SAME signature returning DIFFERENT matrices, and `DM_Contract` dispatches through that one
  signature.  The over-complication is real: FittedEpsXc is a full Dynamic_CC with its own fitter,
  version stamp, and cached matrix, to deliver ONE scalar per cycle — and it has already been
  cloned (`FittedEpsCPol`, Imp/FittedVcorrPol.C:78-98).
  **USER DESIGN DIRECTION (2026-08-05): a framework RE-ALIGNMENT for everything ex/corr/xc, in the
  Hamiltonian library, structure-neutral (not PW/GPW-specific).**  All other terms assume
  E_zz = D·V_zz — energy computed from the (possibly cached) matrix V_zz.  The ex/corr/xc terms
  break that assumption (E_xc = D·⟨ε̃⟩ ≠ D·⟨v_xc⟩).  So: a new interface (extension of
  `tDynamic_HT`?) that forces **`GetVMatrix()` and `GetEMatrix()`** and sorts out which gets used
  where (V → Fock assembly; E-matrix → the energy contraction).  This dissolves the DM_Contract
  signature collision (no second Dynamic_CC object needed), subsumes FittedEpsXc AND its
  FittedEpsCPol clone into one seam, and gives GGA/corr a home from day one.  The quadrature
  routes (Delta_XC/PW_XC compute E by grid integral, NO E-matrix) are handled by the user's
  companion design (2026-08-05): **a second variant/alternate of the `tDynamic_HT` face tailored
  for quadrature-evaluated terms; the Hamiltonian keeps SEPARATE LISTS of each term type and
  loops through them — perfect LSP in action** (no term is ever forced to fake a matrix it
  doesn't have; new term types are easy to add as new lists).

  **IMPLEMENTATION PLAN (grounded 2026-08-06, compute-free pass over the term hierarchy).**  Two
  facts found while reading it change the shape of this item:
  1. **The "separate lists per term type" design is ALREADY IN PRODUCTION.**  `tHamiltonianImp` keeps
     THREE lists — `itsSHTs` / `itsDHTs` / `itsHF_HTs` — behind three faces (`tStatic_HT`,
     `tDynamic_HT`, `tDynamic_HF_HT`) with three different `GetMatrix` arities, and both the matrix
     and energy loops simply walk each list in turn (Imp/HamiltonianImp.C:63-82).  So adding a fourth
     kind is idiomatic here, not novel — it is the user's own pattern, already load-bearing.
  2. **A term whose energy is NOT D·V already exists, and it is the precedent to copy.**
     `tDynamic_HF_HT`'s own doc: it "is deliberately NOT a `tDynamic_CC`: its energy comes from its
     OWN cached blocks (`DM_ContractBlocks`), not a per-irrep `GetMatrix` round-trip."  That is
     exactly the E≠D·V escape the xc family needs, and it requires NO new machinery.
  - ~~**Mechanism (small, and it deletes code):** `FittedVxc` builds its ⟨i|ε̃_xc|j⟩ blocks into a map
    and takes E_xc from `cd->DM_ContractBlocks(...)`, exactly as `Vee`/`Vxc` already do.~~
    **THIS ROUTE IS BLOCKED — found on contact, 2026-08-07.**  `Vee`/`Vxc` can build a whole-system block
    map only because they are `tDynamic_HF_HT`, whose `GetMatrix` RECEIVES `wholeBasis` (the composite
    basis — the cross-irrep view).  `FittedVxc` is a `tDynamic_HT`: its 3-arg `GetMatrix` gets ONE irrep's
    basis and never the composite, so it has no way to ENUMERATE the irrep blocks a `DM_ContractBlocks`
    map must contain.  Faking it (latch a `map<Irrep,const odftbs_t*>` from successive `CalcMatrix` calls,
    à la `Dynamic_HF_HT_Imp::itsWholeBasis`) would add exactly the raw-pointer latching R2.9(iii) already
    flags, AND is unsound in timing: `GetEnergy` runs on the density AFTER the step (SCFIterator.C:279
    contracts ρ_out while the Fock was built from ρ_in), so cached ε blocks would be fit against the wrong
    density.  Plus the string→Irrep key change below.  Enumerating the irrep blocks is precisely the job
    `DM_Contract` already does — that callback IS the density telling the term its own leaf bases.
  - **WHAT LANDED INSTEAD — the user's own V/E face, used as the mechanism.**  The collision was never
    `DM_Contract`; it was one SPELLING.  `tDynamic_HT` already IS-A `tDynamic_CC` (Hamiltonian.C:56), and
    both faces named the method `GetMatrix`, so one term could expose only ONE matrix — hence a second
    object to carry the other.  So `tDynamic_CC`'s method is now **`GetEMatrix`** (ChargeDensity.C), with
    `tDynamic_HT::GetEMatrix` defaulting to `GetMatrix` (E = D·V, true for every term but the xc family)
    and `IrrepCD::DM_Contract` calling `GetEMatrix`.  `FittedVxc` and `FittedVcorrPol` now override it
    with their ε fit — so **`FittedEpsXc` AND its clone `FittedEpsCPol` are both DELETED**, and the two
    ε fitters became plain members beside the v fitters (same fit basis, so the 3-centre setup is still
    shared).  Zero new state, no key change, no basis latching, and no numbers move: `DM_Contract` calls
    exactly the same fit+`Overlap` on exactly the same `(bs, spin, cd)` as the adapters did.
    - The V half deliberately KEPT the name `GetMatrix` rather than becoming `GetVMatrix`: that rename is
      a large mechanical sweep (`tStatic_HT`/`tDynamic_HT`/`tDynamic_HF_HT`/`tHamiltonian` + every term
      and call site) with no behavioural content.  Do it as a standalone naming commit if wanted — the
      V/E DISTINCTION is already in the type system and documented on both faces.
  - **Still outstanding (the declarative half, unchanged):** the quadrature terms (`Delta_XC`, `PW_XC`)
    own a mesh, integrate E directly, and have NO E-matrix — they must never be forced to fabricate one.
    That is the SECOND term list / separate face, and it is untouched by this commit.
  - ~~**Snag to fix en route:** `DM_ContractBlocks` is still keyed by `std::string` BasisSetID~~ — moot for
    V1.3 now (nothing here touches `DM_ContractBlocks`), but the observation stands on its own and is worth
    keeping: V1.4 moved `DM_RhoAtPoints` to `Irrep`, and `tHT_Common`'s term cache was ALREADY
    `std::map<Irrep,hmat_t<T>>` (HamiltonianTerm.C:23).  So Irrep is the established key on the term
    side and BasisSetID is the odd one out; extend V1.4 to `DM_ContractBlocks` in its own pass.
    (This also retro-justifies V1.4: the term caches had been Irrep-keyed all along.)
- **V1.4 ✅ DONE `80fc2ae8`. `DM_RhoAtPoints` Phi key → Irrep (USER RULING 2026-08-05).**
  User model: a CompositeCD is a list of IrrepCDs — not BasisSetCDs; one irrep points to one IBS.
  `Irrep` has meaning in the real world outside of code; `BasisSetID` is purely a code construct,
  whose job is the DB cache's CROSS-RUN caching (same irrep, different radial functions) — caching
  grids inside one SCF run is a different concern and should not borrow its key.  Verified
  enablers: `Irrep` is designed as a std::map key (`operator<`, src/Symmetry/Irrep.C:25) and every
  IBS exposes it (`IrrepBasisSet::GetIrrep(const Spin&)`, IrrepBasisSet.C:62).  The spin detail
  is resolved (user, 2026-08-05): **Spin=None is a valid spin state** — key on
  `GetIrrep(Spin::None)`, the spatial irrep.  → **READY**: change `XC_GridEngine::itsPhi` +
  `DM_RhoAtPoints` signature from `map<std::string,...>` to `map<Irrep,...>`.
- **V1.5 `G_FieldEvaluator` 3-face ISP split.**  Verified: 11 public methods with near-1:1
  method↔consumer mapping — `EvalField*`/`ForwardFFT`/`GridCoeff`/`FieldCoeffs` →
  OrthoFunctionFitter only; `RhoOnGrid`/`Integral`/`EmitGridReport` → PWTerms only;
  `ApplySpectralFilter` → DensityMixer only (Kerker); `MakeFourierDensity` → SeedCD only.  Three
  disjoint responsibilities (evaluate-a-ΔG_Map, FFT quadrature, one-off analytic services).
  Caveat: deliberately promoted to the "full PW grid-engine DIP seam" in FittingCleanupPlan §K —
  the split must engage that decision.  Correction to the old bullet: PWVxcField is not a 4th
  consumer (helper inside PWTerms; the Hamiltonian reaches the engine via `FunctionFitter::Grid()`,
  FunctionFitter.C:139).  Related: XC_GridEngine's `Lattice_3D::Fold` + dcmplx dependency bars any
  molecular reuse of the quadrature engine — same design session.
- **V1.6 `tDM_CD::Accumulate*` — face split, NOT pure-virtual.**  Verified override matrix: every
  concrete family relies on the default for exactly 2 of the 4 (IrrepCD lacks `*All`;
  Composite/Polarized lack `*Both`) — pure-virtual just forces 6 new asserting stubs.  The `*Both`
  pair is an internal Composite↔leaf collaboration protocol (only called from
  Imp/CompositeCD.C:48,59) — no business on the public face.  Shape: a whole-system face and a
  pair-partner leaf face (the `tSpinResolved_CD` cross-cast idiom).  NDEBUG hazard on record:
  these void assert-only bodies are silent NO-OPS in Release — a bare IrrepCD through
  `Vee::AccumulateAll` yields a zeroed J and a silently wrong Fock.
- **V1.7 The periodic trio (`GetFourierDensity`/`GetRhoOnGrid`/`GetRepulsion3C`) — 9 asserting
  stubs, the largest LSP block in qcChargeDensity.**  Re-declared + NA-asserted on
  Polarized/Composite/IrrepCD for BOTH T (Imp/ChargeDensity.C:97,120,138; Imp/CompositeCD.C:197,
  226,249; Imp/IrrepCD.C:267,286,303).  The correct mechanism ALREADY EXISTS in the same file —
  `FourierDensityBase<T>` gives dcmplx the capability and double an empty base — but the derived
  classes re-declare outside it and assert.  Fix: declarations live only on the dcmplx side
  (if-constexpr-guarded definitions or a `tPeriodic_CD` mixin).
- **V1.8 `IrrepCD`↔`IrrepCD` concrete same-class casts in the hot path** (Imp/IrrepCD.C:84,98,
  218,227: `Accumulate*Both`/`MixIn`/`GetChangeFrom` take abstract `tDM_CD&` and narrow to the
  concrete leaf to touch `itsDensityMatrix`; the in-file comment names "the IrrepCD↔IrrepCD
  idiom").  Abstract→concrete, the pattern the project rule forbids; also makes MixIn
  unimplementable for any future leaf.  Wants a double-dispatch primitive or an abstract
  density-block face.  (Design with V1.6 — same seam.)
- **V1.9 ✅ DONE `38a1ebd6`. `Structure`→concrete-`UnitCell` down-casts in 4 libraries**
  (SCFIterator.C:163 [+ R1.3 slicing], DensityMixer.C:319, Imp/SeedCD.C:26,94).  UnitCell.C:38-41
  itself documents the pattern as the thing to avoid.
  **USER FIX SHAPE (2026-08-05): free pry-out HELPERS, exactly like the atom-symmetry precedent.**
  `Symmetry::Atom::Getl(const Symmetry&)` (src/Symmetry/Atom/Spherical.C:68-73) exists because l is
  needed literally *everywhere*; the free helpers downcast the abstract base to the atomic concrete
  and throw `std::bad_cast` on mismatch — "casting and error handling in ONE place" and readable
  client code.  Do the same for Lattice_3D: `GetReciprocalCell(const Structure*)` /
  `GetReciprocalLattice(const Structure*)` (+ others as they show up — cell volume, frac↔cart).
  Semantically identical to the cast, but the 4 sites become one-liners naming what they WANT.
  Verified demand (the R1.3 answer): every consumer of the Kerker cell snapshot wants exactly
  "give me your reciprocal lattice" — `CreateVxcFitBasisSet` takes the plain `Structure*`, and the
  only concrete use is `MakeReciprocalCell()` → `ReciprocalLattice` (DensityMixer.C:327-328).
  Open (minor): helper-in-a-lib placement — the helpers need `UnitCell` visible, so they live
  wherever `Getl`'s analogue would (qcStructure or a Lattice_3D helper module), NOT in the
  consuming libraries.  Note this is the pragmatic sibling of the abstract periodic-geometry
  capability face; the helpers can land FIRST and the face later behind them (the helper body is
  the only thing that changes).
- **V1.10 Two abstract→CONCRETE basis casts in src/** — Imp/SymmetryAdapted_IBS.C:109,118
  (Orbital_HF_IBS* → concrete SymmetryAdapted_IBS, solely to reach `itsO`) and
  Internal/Imp/Orbital_DHF_IBS.C:89,109 (Orbital_HF_IBS& → Orbital_RKB_HF_IBS_Imp&).  Both are
  "give me your private state" reaches — promote the needed answer to an abstract question on the
  face.  (Unit-test exemption does not apply; these are src/.)
- **V1.10b ✅ DONE (see LANDED). Mixer
  LAYERING: `SCFIterator` should only ever see `tDensityMixer` (USER DESIGN,
  2026-08-05).**  Today the iterator knows about Kerker specifically: it holds a Kerker-only
  `itsKerkerCell` snapshot (SCFIterator.C:182, built by the R1.3 cast at :163) and calls
  `MakeDensityMixer(..., itsBS, itsKerkerCell.get(), itsCD.get())` at :244 — i.e. mixer CREATION
  (and its periodic-basis/cell/seed feasibility probe, DensityMixer.C:319-328, which `cerr`s
  "[Mixer] DISABLED" when the probe fails) lives inside the structure-neutral iterator.  Correct
  placement: build the mixer ABOVE the iterator (CalculateSolid / RunGPW facade level) and inject
  the `tDensityMixer*`, or resolve it in a `SolidSCFIterator` via virtual dispatch.  **This
  SUBSUMES R1.3 and one of V1.9's four casts**: with the mixer arriving pre-built, `itsKerkerCell`
  and its slicing deep-copy disappear entirely (the copy exists only to dodge the dangling
  `Lattice_3D::GetStructure` temporary — a problem the facade doesn't have).  Also relocates the
  DISABLED-fallback decision to a layer that can report it properly.
- **V1.11 Occupation seam: two-phase SCF-WF construction + the `TakeElectrons*` family.**
  (i) `SCFWaveFunction::Init/SetMOM/SetSmearing/AdoptMOMReference/ReleaseMOMReference`
  (SCFWaveFunction.C:39-76): run-config known at construction time delivered by post-ctor setters;
  every concrete WF defends against the un-configured state.  (ii) Five `TakeElectrons*` virtuals
  on `TOrbitals<T>` (Orbitals.C:128-152) with a dual-meaning tuple slot (the doc comment itself
  warns "here the double is −TS, NOT leftover electrons"); every new occupation policy so far added
  a virtual to the abstract face (OCP).  One `Fill(const OccupationPolicy&)` returning a named
  struct; the policy object also feeds the WF ctor.  Cross-ref: this IS SCFStrategyPlan's
  "occupation seam" — design it there.
- **V1.12 `EnergyBreakdown` — 13 public data members (OCP+SRP).**  Every new term family edits the
  struct + totals + `op+=` + Display; `GridChargeLost` is a GPW health DIAGNOSTIC ("not an energy",
  its own comment) and `MinusTS` is WF-side entropy — both riding the energy value object; also
  lattice-only `E_alphaZ` in the neutral struct.  Candidate: keyed contributions + a small fixed
  set of roles for the totals; move GridChargeLost to the run report/IterationTrace (which already
  carries it).
- **V1.13 `ExFunctional` — data + setters on an abstract face, with a hidden-init landmine.**
  Carries `itsChargeDensity*` + `isPolarized` (violates the data-free-interface convention).  Its
  ScalarFunction face dereferences `itsChargeDensity`, set ONLY by `LDAVxc::UseChargeDensity` — on
  the PW/GPW path it is NULL and `op()(r)` would segfault.  `SetPolarized` has NO caller ⇒
  `isPolarized` is permanently true and Gradient's unpol branch is dead; meanwhile GetVxc branches
  on a DIFFERENT flag (SlaterExchange::itsSpin).  Split the value face (GetVxc/GetEpsXc(ρ)) from
  the field face; polarized = a type, not a bool (the `SpinCorrelation` face below it is the
  correct data-free shape).  ~~(Coordinate with R2.6 — LDAVxc is the only setter caller.)~~
  **PROMOTED after R2.6 landed (2026-08-07): this is now a DELETION, not a redesign.**  R2.6 removed the
  ONLY caller of `InsertChargeDensity` (it was `LDAVxc::UseChargeDensity`), so `itsChargeDensity` is never
  set and the entire FIELD face — `ExFunctional::operator()`, `Gradient`, the member — is dead tree-wide.
  The fitter now reaches the field through the `VxcDensity`/`EpsXcDensity` adapters, which take the density
  as a CTOR ARGUMENT, so the hidden init cannot recur.  Executing it = delete the member + the two
  virtuals + `InsertChargeDensity`, drop `ExFunctional`'s `ScalarFunction<double>` base, and retire
  `SlaterExchange::Gradient` (its only reason to exist was that face).  The COMPILER verifies the claim:
  anything still needing the field face fails to build.  `SetPolarized`/`isPolarized` (still callerless,
  hence permanently true) rides along in the same pass.
- **V1.14 Report-emission creep on neutral faces + global bool toggles** — `EmitBasisUsage`
  (WaveFunction.C:48, defaulted no-op), `EmitRadialReport` (IrrepBasisSet.C:68), `EmitGridReport`
  (G_FieldEvaluator.C:60, PURE — forces every implementor), plus function-local-static
  `bool& ReportBandGap()`/`ReportGridCharge()` process-globals that leak state between tests (the
  SCFIterator comment admits it).  Fix: a reporter/visitor that PULLS; toggles on SCFParams.
- **V1.15 `tBasisSet<T>::Create*FitBasisSet` defaults** — the generic body hard-codes
  `Orbital_DFT_IBS<double>` regardless of T (only the explicit dcmplx specializations save it) and
  derefs an unguarded iterator (a 1E/HF-only basis ⇒ null ⇒ UB in Release).  Also
  `CreateXCQuadrature`'s default body is byte-identical in two places (Imp/BasisSet.C:42 ==
  Band_FT_IBS.C:53).  Decide: pure, or a documented "this basis does not fit" contract; hoist the
  shared default.  (Interacts with the V1.1 CreateXCQuadrature hoist.)
- **V1.16 `ProjectedDensity_AO::GetRepulsion3C` asserting default** (Fitting/Imp/FunctionFitter.C:
  28-35) — the overlap-vs-Coulomb metric choice is encoded as "override the OTHER method and this
  one becomes poison", an implicit unenforceable pairing (NumericCD must just remember).  Make the
  metric an explicit refinement face.
- **V1.17 `tWaveFunction::GetSpinDensity()` returns null as the unpolarized answer**
  (WaveFunction.C:39) — a capability half the hierarchy lacks, on the base, every client
  null-checking a raw pointer.  Correct idiom one library over: `tSpinResolved_CD` as a cross-cast
  face.  Move to a `SpinResolvedWF` face — also aligns with the spin-native-is-primary bias (the
  polarized WF is the primary type, not a special case bolted on via nullable getter).
- **V1.18 `FourierMixCD` tell-don't-ask + `MakeDensityMixer` ISP.**  `RhoTilde()` hands out the
  raw ΔG_Map and PulayMixer runs the whole DIIS algebra outside the density
  (DensityMixer.C:183-233); `SetRawRho` + external `RasterKerker` is a get/compute/set straddle —
  the mixing algebra wants to be density-side ops or a dedicated mixable-density face.  Also
  `MakeDensityMixer` takes `const tDM_CD*` but uses only GetTotalCharge + a FourierDensity cast
  (DensityMixer.C:312-320) — excludes the matrix-free seeds from seeding the mixer BY TYPE, not
  intent.
- **V1.19 Seed assembly: give `Structure` the question.**  Seed code reads concrete `Atom` public
  fields (Imp/Seed.C:102-104 `a->itsZ`,`a->itsR`; Imp/SeedCD.C:91 `itsSpinFlip`) and clones the
  UnitCell into (unflipped, flipped) groups because `MakeFourierDensity(st, formFactor(Z,g2))` is
  species-keyed — the flip-group sub-cell duplication is the SYMPTOM; the neutral face yielding
  raw atoms is the CAUSE.  A `ForEachSite(fn(Z,R,flip))` / per-atom form-factor visitor (or the
  per-atom-index `G_FieldEvaluator` overload — weigh against the pseudo-wall pin) removes the
  field access, the UnitCell casts, AND the sub-cell duplication.  Also `MakeSeedDensity`'s closed
  enum switch has assert(false)+nullptr arms (Imp/Seed.C:152,158) — in Release an unsupported
  (strategy × T) cell returns null and SCFIterator silently runs a core guess;
  compile-time-over-runtime says make unsupported combos unconstructible or throw.
- **V1.20 `SymmetryAdapted_IBS` → Internal: decide the library-family rule first.**  Verified: all
  consumers inside the src/BasisSet/ tree BUT in a different CMake target (qcMolecule_BS).  Moving
  to `.Internal.` creates a cross-LIBRARY internal import under the CLAUDE.md rule.  Precedent the
  rule is already bent: `qchem.BasisSet.Internal.GMap` imported by src/Fitting + src/ChargeDensity.
  Decide: is the qcBasisSet* family "one library" for Internal purposes?
- **V1.21 `BandStructure.C`: promote or demote.**  Confirmed test-only (sole import =
  tests/BandStructureUT.C:13); worse, tests/PlaneWaveUT.C:59 defines its OWN local `SolveBands`
  instead of importing.  Either promote (band plots are on the viz roadmap) or demote into the
  test tree; either way kill the duplicate.
- **V1.22 `MakePeriodicBeckeMesh` ε-tail drops vs orbit consistency (W2c find).**  The builder's
  borderline drop decisions (`<eps` screens + `w>0` keep) are per-point and bit-sensitive, so the
  site-adapted caller post-filters orbit-incomplete points (`CreateSiteAdaptedBeckeMesh`).
  Cleaner: make the drop decision ONCE per representative (angular dir × radial shell) and apply
  it to the whole atom orbit inside the builder — removes the second fold pass + the filter.
- **V1.23 `Symmetry::Lattice_3D::DirectOf`** — currently unused after the CreateXCQuadrature move
  (GPW uses its native direct ops).  Keep (documents the U=Wᵀ convention; T3 will want it) or fold
  its doc into `ReciprocalOp` and drop — decide at the refactor session.
- **V1.24 `GDMParams::FDMax` naming + the fallback commit (2026-08-03 imposed×GDM investigation).**
  (i) `FDMax` reads as a step size but is the ENGAGEMENT gate on ‖[F′,D′]‖ (the geodesic step size
  is the quadratic-model `itsStdef` capped by `Trust` radians); rename to something like
  `EngageBelowFD` and consider making the norm intensive (per-element or per-electron) so one value
  transfers across basis sizes.  (ii) `DirectMinStep`'s 12-backtrack line search COMMITS a tiny
  non-descent step on fallback (SCFIterator.C ~L424) — the measured uphill leak on imposed NaF-SR
  (+23–56 mHa over 100 iterations along projector-curved diffuse directions).  Fix: hold position
  on fallback (or accept `best` only within a noise floor); pair with soft-direction
  preconditioning — the 1/(ε_a−ε_i) diagonal Hessian blows up the step exactly along the
  near-degenerate diffuse modes.  Reproducers: DISABLED_ImposedGDMProbe_SiDiamondIBZ (healthy),
  DISABLED_NaFImposedGDMSmearProbe (pathological, NAFGDM_* knobs); GPW_GDMTRACE=1 shows
  DESCENT/FALLBACK per step.
- **V1.25 Minor CD-interface trims** — non-const `Polarized_CD::GetChargeDensity(Spin)` overload
  has no external consumer (removable); `tSpinDensity` holds two raw non-owning `tDM_CD*` with
  unmanaged lifetime.

### V2 — measurements / sweeps

- **V2.1 Spin-native collapse: `Delta_XC` ×2 vs `Delta_XC_Pol`+`Delta_VcorrPol`.**  Per the
  spin-native-is-primary bias the unpol pair should become the ζ=0 collapse: the polarized pair at
  ρ↑=ρ↓=ρ/2 IS the unpolarized answer (singlet-collapse gate proves it to 0.12 mHa).  Collapsing
  halves the class count at ~2× XC pointwise work for closed shells (the Ham_PP trade).  Decide
  with a perf measure.  Falls out for free if it lands: XC_GridEngine's second (scalar) rho cache
  dies (pair only; total = up+dn).  Supporting duplication catalog (verified near-verbatim):
  Delta_XC vs Delta_XC_Pol CalcMatrix/GetEnergy differ only in Rho vs RhoPol; VxcPol vs
  FittedVxcPol are both "own two channel terms, cast to Polarized_CD, dispatch on spin" — a
  generic pair-dispatcher should span the HF pair too (the user's original "generic pol Vxc" ask).
- **V2.2 GPW default seed policy.**  `GpwOptions.seed` defaults to `Uniform`, which the Na-doublet
  campaign showed has a STABLE wrong basin for electron-sparse systems (lone-electron doublet
  converged 72 mHa high with every health metric green).  The molecular facade already defaults
  DFT to SAD.  Candidate: default GPW to `IonicSAD` (SAD-family), Uniform = explicit opt-in —
  needs a suite sweep since every pinned GPW anchor re-seeds.
- **V2.3 Polarized PLANE-WAVE Vxc fit route** — `Ham_PW_DFT` polarized currently THROWS for
  `VxcFit::PlaneWave`: per-channel PW_XC needs per-spin rho-grid caches (PW_XC's `itsRhoGrid` is
  keyed on `cd->Version()` alone, which a polarized density aliases across channels — the trap the
  engine's RhoPol pair-cache fixes).  Design note in the throw message.

### V3 — repro / campaign bugs (Spin-SAD, 2026-08-04)

- **V3.1 Polarized ATOMIC solver fails on an empty minority channel** — `AtomCalculation(11, 0,
  {.pseudopotential=true, .valence=1, .pol=Polarized})` (Na q1 doublet: nUp=1, nDown=0) dies with
  "Invalid setup of symmetric matrix" before the first Fock.  Same failure class as the tier-4b
  finding (`cSCFAcceleratorDIIS::GetNProj` segfault on an empty EC).  Worked around for the Na
  library entry (1 valence electron ⇒ the pair is EXACT without an SCF, hand-constructed in
  atomic_valence_densities.json), but any future fully-polarized one-electron species hits it
  again — the valgen `--spin` path should either fix the empty-channel atom SCF or synthesize the
  exact pair itself when nDown==0.
- **V3.2 Unoccupied-l shells break the atom SCF** — the same Na valgen run ALSO failed when the
  recipe carried the (unoccupied) p window from valence_lowq_sr (`--shell 1:2:0.09:0.3`): "Invalid
  setup of symmetric matrix" even unpolarized-adjacent.  GenerateValenceBasis documents "higher-l
  polarization shells ride along un-validated (as intended)" — but at least the polarized
  GenerateSeedDensity path chokes on them.  Reproduce + fix or document the restriction.

### V4 — watch triggers (act when the trigger appears)

- **V4.1 `CollocMemo` dual duty** — replay memo (exact-D level densities) + adjoint D-screen in
  one struct with different lifetimes (last-D vs union-Dscr).  Fine today; split when a THIRD
  consumer appears.
- **V4.2 `SolveSPD`/NNLS in `SymmetrizeMesh.C`** — module-private dense Cholesky + Lawson-Hanson
  NNLS beside the mesh code; promote to qcMath when a SECOND consumer of small dense LS/NNLS
  appears (Blaze has no NNLS).

## DECIDED-ELSEWHERE

- **D1 `Band_DFT_IBS`** — verified fully dead in code (zero imports/derivations/casts; only the
  CMake entry).  Governing record: FittingCleanupPlan §D deliberately KEPT the module as the
  intended future-GPW `<double>` interface.  NEW FACT flagged for re-decision: GPW has since
  landed (`GPW_IBS.C`) and does NOT implement it.  (Stale comments → R2.4.)
- **D2 "Why is a FourierCD different than a tChargeDensity?"** — governed by the
  pw-fitting-uniform-interface pin (no "Fourier" in abstract faces; PW fitting looks identical to
  molecular fitting) — same campaign, don't re-derive here.  User context preserved: at a high
  level Fourier/PW is an implementation detail; `G_ERI3` (src/BasisSet/Internal/GMap.C) is the
  multi-purpose data structure harmonizing the return-type differences — it ends up being a spec
  for the needed data-structure flexibility; high-level code should not have to worry about G or
  no G, PW, Fourier.  Verified corrections to the old bullet: the class is `FourierDensity` — a
  standalone cross-cast MIXIN, not a tChargeDensity; SCFIterator.C is NOT a consumer (stale
  imports → R2.4); missed consumers = DensityMixer (3 cast sites), tPolarized_CD ↑+↓ forwarding,
  tComposite_CD k-sum forwarding, + 3 direct implementors (SeedCD, PolarizedSeedCD, FourierMixCD).
  Taxonomy note: the seeds (NumericCD/SeedCD/PolarizedSeedCD) are tChargeDensity but neither
  DM_CD nor FittedCD — a third leg beside the doc's two primary variants (DM_CD = ERI-capable
  first stop after any SCF iteration; FittedCD = transient DFT fit: Create/DoFit/GetRepulsion),
  orthogonal aspects: Polarized/Unpolarized, Composite (multiple irreps — is this really
  different than Polarized, since irreps can hold spin?), version tracking, mixing algorithms,
  spin density (plotting).
- **D3 `FIT_SF_ABS::SymmetrizeRaster`** — (SRP) symmetry op on a fit-basis interface.  Removable
  once a τ-acting DIRECT-raster fold variant of `FoldGrid` exists (T3 groundwork): precompute the
  raster fold, ctor-inject into `tComposite_CD` beside `itsPointOps`, delete the virtual (uniform
  route's voxel-shift moves out of `PlaneWaveFit_IBS`).
- **D4 The basis-as-policy-carrier virtual family** — `GetReciprocalPointOps` /
  `GetDetectedReciprocalOps` on `tBasisSet` + `SetSymmetryOps` on the GPW evaluator +
  `GPWParams.imposeSymmetry` (renamed from `reduceBZ` 2026-08-02): all exist because the run's §3
  symmetry policy rides the basis.  Verified sole consumers: CompositeWF.C:169 (handing ops to
  tComposite_CD's stream fold; empty for atoms/molecules — lattice policy riding the neutral
  base) and GPW_SCF_UT.C:310 (src/-dead diagnostic getter).  Governing record: the §3
  `SymmetryPolicy` object (SymmetryUpgradePlan / Lattice_3D-owned space group) — when plumbed,
  these collapse to consumers reading ONE context and the per-factory bool dies.
- **D5 Becke `MeshParams` ergonomics under `imposeSymmetry`** — `mp.nAngular` IS honored but
  `mp.angular` (the GL/Lebedev scheme) is silently REPLACED by the site-adapted rule (announced,
  not warned).  Wanted: `nAngular<=0` → auto-resolve the calibrated default; console warn on
  scheme override; a real error (not a bare assert) when the requested L is unachievable for a
  low-symmetry site (C1/Cs seed-pool exhaustion).  Lands with the `SymmetryPolicy`/facade pass.
  (The degree-typed `angularDegree` interface half is executable now → R2.15.)
- **D6 `BeckeXCParams()` lives in the TEST file + `ResolveXCMesh` (test driver)** — the de-facto
  PRODUCTION Becke recipe and the run-policy resolution (Auto grid × imposeSymmetry interplay)
  living in the integration-test harness; both belong with the facade/driver once the policy
  object exists (library beside `MeshParams`; tests read it, not the other way round).
- **D7 The `dynamic_cast` survey = FittingCleanupPlan §C** (the one surviving item there; the
  "I want more" vs "what are you" criterion is written there).  Run §C as part of THIS session —
  the cast findings above (V1.8, V1.9, V1.10) are its seed list; give survivors the custom
  exceptions (subsumes R2.5's throw work).
- **D8 Standing pin governing every fit-touching item here**: fit quality is measured by
  grid-convergence of ρ/property vs a fine reference — NEVER ΔE_total (fits are non-variational).

## DONE (record)

- **BZ creep on the neutral `Symmetry` base (ISP) — DONE 2026-08-04 (same day): the ATOM SHELL
  CONVENTION landed.**  `BlochQN::GetDegeneracy()` = star (ctor converts the k-mesh layer's w_k;
  asserted integer), `GetWeight()` = uniform 1/N_mesh, `StarSize()` + the `EnergyLevel` report
  scaling DELETED (the 8/8 wedge display is now the physical occupation), `Crystal_EC::GetN` =
  star×Nval per block in insulator mode / plain Nval for the global-μ total.  No factory/KBlock/
  call-site changes (BlochFactory keeps taking w_k); free meshes (star=1) bit-identical; g=w·degen
  invariance verified (Al metal); imposed Si diamond IBZ reproduces E to all digits; 656/656.
  The base is down to TWO Bloch-flavored members (GetWeight, MergeAcrossIrreps) — the remaining
  capability-face split is now small enough to fold into any future Symmetry-base touch.
  *(Original analysis, for the record:)*  The base carried three Bloch-flavored defaulted virtuals
  (`GetWeight()` pre-existing; `MergeAcrossIrreps()` + `StarSize()` from the 2026-08-03
  MergeTol/IBZ-table fixes).  Rather than split a capability face, the user's design was adopted:
  a wedge k-block is a SHELL exactly like an atom's l-block (one stored representative,
  degeneracy = the symmetry copies), so **`BlochQN::GetDegeneracy()` = star size** (spatial; spin
  stays layered by `Irrep`) and `StarSize()` is DELETED.  A coordinated convention switch — the
  invariant is w_k × (per-block quantity), so the star factor must leave the weight in the same
  commit: `GetWeight()` becomes plain 1/N_mesh (most of the weight plumbing through `ReduceToIBZ`
  → `KBlock` → `BlochQN` evaporates); `Crystal_EC::GetN` per wedge block becomes star×Nelec (each
  band holds star×2 natively; the level table needs NO report-layer scaling); every w_k consumer
  rebalanced (density sums, `tIrrepWF::GetEntropyTerm`; `FillOrbitalsGlobalFermi`'s g = w·degen is
  INVARIANT — sanity check).  `MergeAcrossIrreps()` stays (unfolded meshes still carry star
  partners as separate blocks).

## Answered questions (from the pre-reorg sections)

- **"XC_GridEngine — I don't know what this is."**  It is the shared XC quadrature engine: holds
  the mesh, the fold, the geometry-fixed Φ tables, and the ρ caches, SHARED via one `shared_ptr`
  across Delta_XC / Delta_XC_Pol / Delta_VcorrPol (Imp/Hamiltonians.C:217) so the Φ tables are
  built once.  "Engine" earns its keep through the sharing — merging it into Delta_XC would
  triple the Φ tables.  Rename candidate: `XC_Quadrature`.  The "Becke route is pure QUADRATURE"
  phrasing = fit basis is (pseudo) delta functions, per the user's note.
- **"PW_XC: what distinguishes it from FittedVxc?"**  Genuinely different mechanics, not just
  PWs: raster ρ by inverse FFT per density serial, the RAW-adjoint vs BALL-fit fork, energy by
  direct grid quadrature (no ε fit, no DM_Contract).  And yes, GPW uses it: there is no Ham_GPW —
  `Ham_PW_DFT` + GPW basis with `UnitCellKind::Uniform` instantiates PW_XC
  (Imp/Hamiltonians.C:246-249).
- **"Band_DFT_IBS — is this no longer used?"**  Dead in code; deliberately kept — see D1.
- **"BandStructure.C only consumers seem to be unit tests."**  Confirmed — see V1.21.

---

# Appendix — V1.2 Orbital_PP_IBS feasibility probe (2026-08-05, read-only)

**Verdict: NO hard blocker.  The inversion works; four real roadblocks (all with mitigations)
and two minor ones.**

What the basis side ACTUALLY consumes from `LocalPotential`/`SeparablePotential` (verified — no
implementation reaches for rloc/C1..C4/h_ij/Zion):
- Local: `FormFactor{,Long,Short}(Z,G2)`, `FormFactorG0Short(Z)`, `Vloc(Z,r)`, and the optional
  `ShortRangeGaussian(Z) -> {c,n,alpha}` cross-cast face.
- Separable: `NumProjectors/AngularMomentum/Coefficient` (Coefficient is already a DIAGONALIZED
  per-projector scalar — no D-matrix reaches the basis), `Projector(Z,p,q)`, `BetaR(Z,p,r)`, and
  the optional `BetaGaussian` face.
- `Zion` is consumed only term-side (IonIon).

Roadblocks, ranked:
1. **The analytic-FT fast path is load-bearing, not optional.**  PW has NO real-space local path
  at all, and GPW's G-space local route is a BOX-INDEPENDENCE requirement (Evaluator.C:836-841),
  not a speed hack.  An opaque-radial-function argument is a non-starter.  Mitigation: the
  neutral field argument is DUAL-SPECTRAL — `ValueR(Z,r,range)` AND `ValueQ(Z,q2,range)` — which
  is what `LocalPotential_Q` already is, renamed (the TYPES are already neutral; only the NAMES
  are PP-flavored).
2. **Three optional capability faces select the fast paths via dynamic_cast**
  (`LocalPotential_Gaussian`, `SeparablePotential_R`, `SeparablePotential_Gaussian`) — re-express
  as neutral optional faces; preserve the `MultiSpecies_*` forwarding.  Bulk of the mechanical work.
3. **The G=0 alignment leaks into the basis** — `FormFactorG0Short` is read basis-side
  (Evaluator.C:1039) to keep the analytic lattice sum consistent with the drop-ΔG=0 convention.
  Neutral face needs a `CellMeanQ(Z,range)` scalar or the two GPW local paths silently disagree.
4. **The β sharpness hint is a parameter leak** — `ShortRangeGaussian(Z)[0].alpha` pulled purely
  to size grid levels (Evaluator.C:901-906, the NaF diving-ghost note).  Neutral face needs an
  explicit `EffectiveExponent(Z,range)` (0 = unknown → old behavior, already the fallback).
5. (minor) **Projector BRACKETS are the wrong primitive for PW** — `MakeSeparablePotential` uses
  the (2l+1)P_l(cosγ) addition theorem to avoid the m-loop and never forms ⟨i|β_lm⟩ vectors.
  Keep the projector primitive at MATRIX level (`MakeProjectorMatrix(st, projectorSet)`, D
  inside); GPW's bracket loop stays an implementation detail.
6. (minor) `Integrals_Pseudo::MakeLocalPotential` (unsplit) is unit-test-only — drop, don't port
  (re-plumb its 7 test call sites).

DAG findings (better than hoped):
- qcPseudopotential imports NOTHING from qcBasisSet (only qcMath/qcCommon/qcStructure) — the
  inversion edge qcPseudopotential→qcBasisSet is cycle-free.  BUT the adapter that converts
  `LocalPotential` → the neutral field argument can live in **qcHamiltonian** (which already
  links both and owns the terms) — then qcPseudopotential stays a pure LEAF with zero new edges.
  Strictly better DAG.
- The face + argument interfaces → qcBasisSet; the pure-data `RadialGaussianTerm{c,n,alpha}`
  struct → qcMath (already hosts Monomial/CartTerm, links nothing).  Do NOT reuse
  `Molecule::LatticeSum1E::GaussianFunction` as the public type (qcMolecule_BS links UP to
  qcBasisSet → cycle).
- After the retarget, `PWTerms.C`'s two casts go to `BasisSet::Orbital_PP_IBS<dcmplx>` and
  qcPseudopotential drops off qcLattice_BS's link line (modulo the two test TUs importing
  GTH_Potentials directly).

Structure-neutrality check: the primitives are NOT inherently periodic — `PP_Local`/`PP_NonLocal`
compute exactly ⟨i|V_species|j⟩ and the D|β⟩⟨β| loop by mesh quadrature, so a molecular basis can
implement the same face later.  Three PERIODIC CONVENTIONS must be parameters, not baked in:
(1) dropping ΔG=0; (2) the cell-mean subtraction (already gated on `isFinite()`); (3) the
long/short RANGE SPLIT itself (exists to route V_long into the Hartree Poisson — a molecular
implementor answers `Full` and ignores it, exactly what `FormFactorShort`'s default-0 already does).

Minimal neutral signature sketch (2 virtuals instead of 4 — Long/Short/Full collapse into a
`FieldRange` argument):

```cpp
// qcMath — pure data, no deps
struct RadialGaussianTerm { double c; int n; double alpha; };   // c r^{l+2n} e^{-alpha r^2}

// qcBasisSet
enum class FieldRange { Full, Long, Short };                    // a range split, NOT a PP concept

class SpeciesRadialField {
  virtual double ValueR   (int Z, double r,  FieldRange) const =0;  // <- Vloc
  virtual double ValueQ   (int Z, double q2, FieldRange) const =0;  // <- FormFactorLong/Short (fast path kept)
  virtual double CellMeanQ(int Z,            FieldRange) const {return 0;}  // <- FormFactorG0*
  virtual double EffectiveExponent(int Z,    FieldRange) const {return 0;}  // <- the beta grid hint
};
class SpeciesRadialField_Gaussian   // optional capability, cross-cast
  { virtual std::vector<Math::RadialGaussianTerm> AsGaussians(int Z, FieldRange) const =0; };

class SpeciesProjectorSet {
  virtual size_t Count (int Z) const =0;
  virtual double Weight(int Z, size_t p) const =0;              // <- Coefficient (scalar; no D-matrix)
  virtual int    L     (int Z, size_t p) const =0;
  virtual double RadialQ(int Z, size_t p, double q) const =0;   // <- Projector(q)
};
class SpeciesProjectorSet_R        { virtual double RadialR(int Z, size_t p, double r) const =0; };
class SpeciesProjectorSet_Gaussian { virtual std::vector<Math::RadialGaussianTerm> AsGaussians(int Z, size_t p) const =0; };

template<class T> class Orbital_PP_IBS {                        // the face, in qcBasisSet
  virtual hmat_t<T> MakeSpeciesFieldMatrix(const Structure*, const SpeciesRadialField&, FieldRange) const =0;
  virtual hmat_t<T> MakeProjectorMatrix   (const Structure*, const SpeciesProjectorSet&)            const =0;
};
```

`HGH_LocalPotential`/`HGH_SeparablePotential` then simply ALSO implement the neutral faces
(rename-level: `FormFactorLong(Z,G2)` → `ValueQ(Z,G2,Long)`).
