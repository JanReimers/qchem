# Open Work — single tracker for the in-flight threads

One place to see everything dangling, so it isn't spread across sessions. Each item: status, the
next concrete step, and where its plan/code lives. **Recommended discipline: close these before
opening new threads.** Recommended order at the bottom.

Branch: `main` (local, ahead of origin). Anchor build/test: `cd build/Release && ninja UTMain && ./UnitTests/UTMain` (146 green as of this writing).

---

## A. Spherical SALC  ·  PARTIAL (clean checkpoint)  ·  plan: `doc/SphericalSALCPlan.md`

Extend point-group SALC adaptation from the Cartesian PG basis to the two spherical bases.

| Stage | What | Status |
|---|---|---|
| S1 | `SphericalShellRep` (real-spherical operation rep, `qcSymmetry`) | ✅ committed `194f9971`, 4 tests |
| S2 | generalize `AoShell`/`BuildOperationRep` to dispatch Cartesian vs spherical | ⬜ next |
| S3 | the two convention-matched extractors (PG_Spherical in-house; libcint-spherical foreign) | ⬜ the bulk / bug-prone |
| S4 | dispatch in `PG::SymmetryAdapt` by orbital-IBS type | ⬜ |
| S5 | end-to-end tests: spherical-adapted SCF == un-adapted spherical SCF | ⬜ |

**Next step:** S2 — add a spherical angular descriptor to `AoShell` and let `BuildOperationRep` call
`SphericalShellRep` (built S1) for spherical shells; keep the Cartesian path byte-identical.
**Lifts:** the facade's `{.symmetry=true}` Cartesian-only guard, once S3a (PG_Spherical) lands.

---

## B. Spin-native DFT (was "D2 polarized")  ·  IN PROGRESS (B1+B2 done)  ·  plan: `doc/SpinNativeDFTPlan.md` + tenet `feedback_spin_polarized_primary`

Reframed per the design tenet: **spin-polarized is the native formulation; unpolarized is the
ζ=0 efficiency collapse** — not "add the polarized special case." Four pieces, all spin-first:
1. spin-native VWN5 correlation (ζ-dependent) — *replace* the paramagnetic-only `VWN_Correlation`.
2. `Ham_DFTcorr` carrying both channels (U = collapse) — today only `Ham_DFTcorr_U` exists.
3. open-shell molecular occupation `(n↑,n↓)` — `Molecule_EC(Ne)` is closed-shell aufbau only.
4. facade multiplicity on `CalcOptions`.

**Plan doc DONE** (`doc/SpinNativeDFTPlan.md`): scopes the four pieces into staged B1–B4, grounded in the
real types. Key finding surfaced: **correlation does NOT separate by spin channel** (exchange does) — so
`FittedVxcPol`'s two-independent-channel split can't carry correlation; v_c^σ(ρ↑,ρ↓) couples both channels
through r_s and ζ, needing a two-channel functional face + a `FittedVcorrPol` that fits against the full
`Polarized_CD`. Full VWN5 spin params (ferro + spin-stiffness) tabulated in the doc.

**B1 DONE** (uncommitted): spin-native VWN5 in `VWN_Correlation.C` — para/ferro/spin-stiffness branches via
a shared `Gval`/`dGdx`, `EvalRZ(ρ,ζ)` returning ε_c + (r_s,ζ) partials, public two-channel face
`GetEpsC(rUp,rDn)`/`GetVc(rUp,rDn,Spin)`. Scalar face = ζ=0 collapse, byte-identical (anchors unmoved).
`LDA_XC_UT` extended: vs libxc `LDA_C_VWN` polarized on an (r_s,ζ) grid to 1e-9 + a collapse-consistency
check. **150/150 UTMain green.**

**B2 DONE** (uncommitted): abstract `SpinCorrelation` face + `VWN_Correlation` implements it; `FittedVcorrPol`
fits v_c^σ against the full `Polarized_CD` (+ SAD-seed ρ/2 fallback), energy via `FittedEpsCPol` Dynamic_CC
(polarized DM_Contract sums channels ⇒ ∫ε_c·ρ); `Ham_DFTcorr_P` mirrors `_U`; both Factory throws un-gated.
Anchor `M_DFT.WaterPolarizedLDA` collapses to the unpolarized LDA anchor to 1e-6. **151/151 green.**

**Next step:** **B3** — open-shell molecular occupation `Molecule_EC(nUp,nDown)`; closed-shell `Molecule_EC(Ne)`
becomes the collapse (delegating ctor, existing call sites untouched). `GetN(Irrep)` returns per-spin counts;
the aufbau-per-spin SCF loop already handles asymmetric channels (no SCF change). Then B4 (facade multiplicity).
**Done already (D1):** closed-shell molecular HF+DFT via the facade — unified `Model` enum + Factory
resolver, LDA + Xα, `8b8df1d0`.
**Fixed (was the B entry bug):** polarized Xα + **SAD seed** segfault — `FittedVxcPol::CalcMatrix`
null-derefed on the spin-agnostic seed (`cd85d13c`); polarized Xα water now converges to the unpolarized
anchor to ~1e-11, and `M_DFT.WaterPolarizedSAD` guards it. So polarized Xα itself is sound; this track is
about the *correlation* side (spin-native VWN5 + `Ham_DFTcorr` two-channel) + open-shell occupation.

---

## C. Namespace unification (review item 4)  ·  ✅ DONE `108ced3b`  ·  plan: `doc/APIErgonomicsReview.md` §4

Whole-tree sweep — **full unification** (user-chosen): not just the four review-named symbols but
*everything* under `qchem::`. Moved: `ScalarFunction`, `Spin`, the `Vector3D` geometry family
(`Vector3D`/`Matrix3D`/`Vector2D`/`Matrix2D`), the `BasisSet::` tree, the stray library namespaces
(`Symmetry`/`qcMesh`/`Pseudopotential`/evaluator-helpers/`SpecialFunctions`/`blazem`), and the ~39
global-scope class-definition files (`Structure`/`Atom`/`Molecule`, `UnitCell`, `Lattice_3D`,
`ElectronConfiguration`+`*_EC`, `Irrep`/`MolecularIrrep`/`Orbital_QNs`, `SCFParams`, `Real_BS`/
`Real_OIBS`, `Cache2/3/4`, `ERI*`, `FourierMap`, …) + their `Imp/` units.

Left global by design (intentional vocabulary, per CLAUDE.md "dcmplx is global"): the lowercase `_t`
aliases (`rvec_t`/`rmat_t`/`rvec3_t`/`dcmplx`/…) in `Types.C`, the cmath/IntPow passthroughs, display
tables (`Strings`), `std`-container `op<<` (`stl_io` — ADL needs std reach), `PeriodicTable` lookups,
the `qchem` umbrella, and internal `detail` namespaces.

Gotchas captured for posterity: (a) blaze + stl_io operator re-exports had to be *mirrored* into
`namespace qchem` or they get shadowed by the in-qchem `Vector3D` operators for qchem-scope code
(`Blaze.C`/`stl_io.C`); (b) test/harness `using namespace qchem;` collides with the
namespace-vs-class pattern (`qchem::Hamiltonian::Hamiltonian`) → qualify those class uses; (c)
`forward.H` test-friend decls qualified to `::XxxTests`.

Verified: UTMain 146/146 (matches pre-sweep anchor); `allTests`+`scfrun` build; sampled per-module
suites green. pybind (OFF by default) updated by-name, not compiler-checked.

---

## D. Test → Facade migration + slim QchemTester  ·  ✅ DONE  ·  plan: `doc/TestFacadeMigrationPlan.md`

`qchem::Calculation` is now THE molecular recipe — every molecular fixture drives the public facade and
`TestMolecule` is deleted.

| Stage | What | Commit |
|---|---|---|
| 1 | `CalcOptions` `Engine{MnD,LibCint}` / `Angular{Cartesian,Spherical}` (+ `.seed`) | `585086bb`,`8d936edc` |
| 2 | free `RelativeError(E,Eref)` helper (`UnitTests/TestUtils.C`) | `585086bb` |
| 3 | migrate `M_HF_U` (6), `M_DFT` (3), `M_Sym`→`{.symmetry=true}` (7), `A_DFT` `A_PG` (oracle) | `8d936edc`,`c96188c2`,`b0e85305` |
| 4 | delete `TestMolecule`; slim `QchemTester` to the atom/Dirac harness | `3c0becdd` |

Anchors byte-for-byte (the facade's default `SCFParams` == the old `scf` literal; its auto DFT recipe
== the hand-set SAD+DIIS-from-start; default mesh == `TestMolecule::GetMeshParams()`).  Stayed on the
scaffold by design: `TestAtom`/`TestDiracAtom` + the Z-keyed NIST/Dirac oracle asserts; `scfrun`.

Two additive facade extensions were needed to migrate cleanly (both general-purpose):
`CalcOptions.seed` (for the SAD-seed HF bootstrap test) and `AcceleratorOptions.eMax` already existed
(used by `A_PG` to pin the Z-scaled DIIS gate).  `A_PG`'s NIST checks keep the scaffold's SIGNED
relative-error bound (bounds over-binding only) so pass/fail is identical.

⚠️ **Bug surfaced:** polarized Xα + SAD seed **segfaults** through the facade (`M_Sym` polarized DFT hit
it; worked around with `seed=CoreGuess`).  Latent in the not-yet-built spin-native DFT path — belongs to
track **B**.

---

## E. Retire QchemTester — atom tests onto qchem::AtomCalculation  ·  ✅ DONE  ·  (follow-on to D)

Goal (user, "otherwise new tests appear on the old system") EXCEEDED: not just `QchemTester::Init` retired
— the **whole `QchemTester`/`TestAtom`/`TestDiracAtom`/`TestMolecule` scaffold is deleted** (`97963b63`).
A and B did not block it.

Built a SIBLING `qchem::AtomCalculation` (NOT an extension — the GUI team is actively on `Calculation`, so
its surface stays untouched).  AtomCalcOptions: atomic exponent-pool basis (AtomType + BasisSetAccuracy or
{N,emin,emax}); model/pol/xalpha; the `xc` exchange-functional selector; `pseudopotential`/`valence`;
nAngular=1 mesh; seed=CoreGuess; an `accelerator` json escape hatch.  Model-driven EC, ctor converges once
with SCFParams, GetIrreps/GetOrbital/GetStructure accessors.  Free Z-keyed oracle helpers
RelativeHF/DFT/DHFError(E,Z) in UnitTests/TestUtils.C.

New **public Hamiltonian-library API** (the long-wanted exchange-functional selector):
`qchem::Hamiltonian::XCFunctional{kind,alpha,libxcId}` + `enum XC{SlaterXalpha,DiracVWN,LibXC}` +
`Factory(Pol,st,XCFunctional,mesh,bs)`, and a public PP `Factory(st,element,valence,mesh,bs)`.  Internals
(ExFunctional / Ham_PP_U) no longer leak.

| File | Status |
|---|---|
| A_HF_U / A_HF_P | ✅ `bf36a405`,`8f72981c` |
| A_DHF (DE1 ions, DHF oracle, Phir, fine structure) | ✅ `0e70b05f` |
| A_DFT atomic (Slater-Xα + LSDA) + A_PG | ✅ `54dfcd51` |
| A_DFT_U (libxc via `{.xc=XCFunctional{.kind=XC::LibXC,.libxcId=7}}`) | ✅ `45f88cf9` |
| A_PP (Si PP via `{.pseudopotential=true}`) | ✅ `45f88cf9` |
| scfrun (CLI driver; all models + `--model PP` ions + accelerator tuning) | ✅ `97963b63` |
| delete the QchemTester scaffold | ✅ `97963b63` |

146/146 UTMain green throughout; allTests + scfrun build; anchors byte-identical.  Future: `Calculation`
(molecular) can adopt the same `XCFunctional` selector when the GUI team is ready.

## Bigger milestone on the horizon (not dangling, just flagged)

- **PBE / GGA** — the highest-value functional for the battery north-star, but a real library
  increment (density-gradient machinery on the mesh), not an enum value. The unified `Model` enum is
  ready to list it with a "not wired" throw. See `doc/FacadeDFTPlan.md`.

## Deferred / descoped (recorded so they're not re-litigated)

- Fold `QchemTester` + the pybind bridge onto the facade (test-harness/binding cleanup; the facade
  makes most of `QchemTester` redundant — a good litmus test, but not lib-surface work).
- Container utils to `src/` (`sample_scalar`/`sample_gradient`, `Structure::BoundingBox`) — binding convenience.
- `SCFParams` ASCII rename — DROPPED; solved instead by C++20 designated initializers (`34ccf302`),
  keeping the compact unicode names.
- `MolecularSym_EC` → `FixedIrrepOcc_EC` rename — belongs with the queued symmetry-naming cleanup.

## Recently closed (for context)

- API ergonomics additive layer: facade (`a7ef787b`), umbrella `import qchem;` (`53206930`),
  `{.symmetry=true}` Cartesian-guarded (`ff5d0a39`), DFT D1 (`8b8df1d0`).
- HF→DFT cross-test energy contamination (fit-basis `Norm` cache not mesh-keyed) — FIXED `c570b20b`.

---

## Recommended order (user-chosen 2026-06-30)

1. ~~**C (namespace sweep)**~~ — ✅ DONE `108ced3b` (full unification: whole tree under `qchem::`).
2. ~~**D (test → facade migration + slim QchemTester)**~~ — ✅ DONE (`585086bb`…`3c0becdd`).
3. **A (finish Spherical SALC S2–S5)** — NEXT. Parked at a clean S1 checkpoint; resume when convenient.
4. **B (spin-native DFT)** — plan first, build spin-first, sequence toward PBE+U / batteries.
   (The polarized-Xα + SAD-seed segfault that D surfaced is already fixed, `cd85d13c`.)

Note: one library session at a time, GUI on its own branch → land each session at a clean commit and the
whole-tree sweep (C) has nothing to collide with. Hold the line on not opening new threads.
