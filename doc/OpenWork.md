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
| S2 | generalize `AoShell`/`BuildOperationRep` to dispatch Cartesian vs spherical | ✅ new `OperationRep` module (`f5abf1d2`) |
| S3a | in-house `PG_Spherical` extractor `ExtractAoShells(SphData)` | ✅ reads basis's own c2s (zero convention risk) |
| S3b | libcint-spherical extractor (foreign real-harmonic order/norm) | ⬜ only remaining piece |
| S4 | dispatch in `PG::SymmetryAdapt` by orbital-IBS type | ✅ PGData/SphData; libcint-sph still guarded |
| S5 | end-to-end tests: spherical-adapted SCF == un-adapted spherical SCF | ✅ `M_Sym.water_HF_spherical_*` |

**S2 DONE:** `AoShell`+`BuildOperationRep` moved into a NEW module `qchem.Symmetry.OperationRep` (imports
both `CartesianRep` + `SphericalRep`), resolving the `CartesianRep`↔`SphericalRep` cycle — the dispatcher
can't live in either per-shell module. `AoShell` gained a `c2s` (HarmonicC2S) field; `IsSpherical()`/
`nComponents()`; `BuildOperationRep` picks `SphericalShellRep` vs `CartesianShellRep` per shell. Cartesian
path byte-identical (M_CartesianRep/M_SALC/M_Sym green); new `OperationRep.spherical_d_shell_dispatch` test.
Only `SALC.C` + `M_CartesianRep.C` needed import edits (everything else gets `AoShell` transitively via SALC).

**S3a+S4+S5 DONE** (in-house MnD-spherical SALC works end to end): new module
`qchem.BasisSet.Molecule.PG_Spherical.Symmetry` — `ExtractAoShells(SphData)` reads the basis's OWN c2s
(`comps[].terms`) into `AoShell.c2s`, so the rep is built in the basis's exact m-order/coeffs (no
convention to match — the S3a risk evaporated). `SymmetryAdapt` dispatches PGData→Cartesian, SphData→
spherical, else throws (libcint-spherical still unsupported). Facade guard relaxed: `{.symmetry=true}` +
`angular=Spherical` now allowed for `engine=MnD`, still blocked for `engine=LibCint`. Also fixed a latent
S2 gap: `BuildSALCs` computed `nAO` from `monomials.size()` (0 for spherical) — now `nComponents()`.
`M_Sym.water_HF_spherical_{unpolarized,polarized}` (water/dzvp O d-shell). **157 UTMain, 14 UTSymmetry green.**

**Next step (only remaining piece):** **S3b** — the libcint-spherical extractor. The bug-prone one: it must
match **libcint's** real-harmonic ordering + normalization (a foreign convention), and libcint-spherical
presents AS a PGData with spherical components (a trap). This is genuinely separable; the in-house spherical
SALC is fully shippable without it.

---

## B. Spin-native DFT (was "D2 polarized")  ·  ✅ DONE (B1–B4)  ·  plan: `doc/SpinNativeDFTPlan.md` + tenet `feedback_spin_polarized_primary`

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

**ALL FOUR PIECES DONE** (B1 `51157449`, B2 `4189af69`, Factory cleanup `7132d645`, B3 `02b68b24`, B4 next commit):
- **B1** spin-native VWN5 (`VWN_Correlation` para/ferro/stiffness + two-channel `GetEpsC`/`GetVc`; vs libxc polarized).
- **B2** `SpinCorrelation` face + `FittedVcorrPol` (fits v_c^σ against full `Polarized_CD`, energy via `FittedEpsCPol`) + `Ham_DFTcorr_P` + Factory un-gate.
- **B3** `Molecule_EC(nUp,nDown)` spin-native; `Molecule_EC(Ne)` = minimal-spin collapse (byte-identical).
- **B4** `CalcOptions.multiplicity` → (nUp,nDown) + auto-promote to polarized + parity-validate.

Anchors: `M_DFT.OxygenTripletLDA` (-149.2562876393, real ζ≠0), `OxygenTripletBelowSinglet` (Hund), `BadMultiplicityThrows`,
`WaterPolarizedLDA` (closed-shell collapse). **154/154 UTMain green**, unpolarized anchors byte-identical throughout.
Also did the **Factory interface cleanup** the user asked for (one DFT build site, dropped the public raw-`ExFunctional*`
overload, honest LibXC-polarized throw).

**Remaining for the north-star (separate increments, NOT this track):** spin-native GGA (PBE — gradient machinery),
LibXC-polarized (the libxc wrapper needs two-channel `xc_lda_vxc`), and +U. The LDA *interface* is now spin-native end to end.
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
- **Automated BASIS-FUNCTION trimming for ill-conditioned spans (USER, 2026-08-14).**  The MnO
  campaign's verdict: every "clever" cure that discards DIRECTIONS in coefficient space failed the
  user's earlier trials, while removing WHOLE diffuse AO functions (manually, then via
  `CholeskyPivoted` + `orthoTol`) is what actually works — run 58 killed the variant-B collapse by
  dropping 6 named functions at the door.  Pivoted Cholesky already identifies WHICH functions are
  redundant (greedy residual-pivot selection, indices reported).  The future feature for
  inexperienced users: promote this from an ortho-time knob to a first-class VET-stage policy —
  auto-select the kept sub-basis (auto gap tolerance already exists, `LASolverLapack.C
  detect_null_gap`), report it as a basis decision (species/shell/exponent, not bare indices), and
  perhaps regenerate the trimmed basis via `qchem.ValenceBasisGen` so the user sees a BASIS, not a
  filter.  Direction-space (eigen/SVD canonical-ortho) trimming stays a numerical fallback, never
  the user-facing policy.
  **NOT display-only (user 2026-08-14): the trim must happen BEFORE anything downstream is built.**
  Everything falls out of the surviving function list — grid-ladder depth (Max/MinExponent), the
  collocation streams and their caches, KB projections, the whole per-pair machinery — so building
  all of it for the full span and then filtering at ortho time does the dropped functions' work for
  nothing AND leaves their pairs in every stream/cache.  Vet decides the kept sub-basis; the basis
  the run constructs IS the kept sub-basis.  Corollary (same user note, "wrong logic — harmless for
  now"): the rank decision is a property of S, i.e. of the BASIS, made ONCE — today each spin
  channel's LASolver independently re-derives it (the doubled `[ortho]` console line), which is
  only coherent because S is channel-independent; the vet-stage placement makes it structurally
  single-decision.
  **And it must be SYMMETRY-EQUIVARIANT (user catch, 2026-08-15): drop whole ORBITS under the
  (magnetic) space group, never individual AOs.**  The greedy per-function pivot resolves
  symmetry-TIED pivots by numerical noise: runs 58–60 dropped O₁'s p(0.18) but O₂'s s(0.15), and
  DIFFERENT d(0.18) m-component pairs on the two Mn — the kept span breaks the O-site equivalence
  and the Mn flip-mirror at the ~sub-1% level (pivot < 1e-2 ⇒ ≥99% reproducible), the same order
  as run 59's site-moment asymmetry (0.6756/−0.6804, m_net −0.005) and plausibly feeding its
  8%-class Shubnikov flip-defect audit reading.  Shell-level variant spans (VA/VB, the CP2K
  pairing) restore exact site symmetry by construction — another reason the exact-span pairs are
  the trustworthy comparison rows.

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
3. ~~**B (spin-native DFT)**~~ — ✅ DONE (B1–B4, see above). Done before A this session. The LDA interface is
   now spin-native end to end; PBE-GGA / LibXC-polarized / +U are separate increments toward the battery north-star.
4. **A (finish Spherical SALC S2–S5)** — NEXT remaining thread. Parked at a clean S1 checkpoint; resume when convenient.

Note: one library session at a time, GUI on its own branch → land each session at a clean commit and the
whole-tree sweep (C) has nothing to collide with. Hold the line on not opening new threads.

## Recommended order (USER re-prioritization, 2026-08-15 — supersedes the 2026-06-30 list above)

Context: the MnO campaign's {AFM,FM}×{qchem,CP2K} matrix (doc/SphericalLatticePlan.md) exposed the
COST gap as the binding constraint — CP2K runs the same cell in minutes/<1 GB where qchem takes
hours/12 GB — and the accuracy comparison itself is stalled behind code-health items.  The
symmetric VA matrix (runs 61/62 + CP2K va/vb decks) completes overnight; after that, NO further
long MnO runs until 1–3 and 5 land.

1. **Close the RUNTIME gap vs CP2K.**  ✅ **ROUND 1 DONE 2026-08-15 (MnO Γ benchmark 171 s → 101 s;
   per-iteration ≈30 s → ≈9 s, physics bit-for-bit unmoved) — record + next steps in
   doc/GPWPlan1.md "THE RUNTIME GAP, MEASURED".**  The charter's fast-recompute campaign was aimed
   at the wrong bucket for the runs we actually pay for: **with the streams cached the pair loops
   are 4% of the run**; the cost was the atom-centred XC mesh (Φ tables 55.6 s, ρ sampling 23.8 s,
   H_xc quadrature inside the 21.8 s iterate residue) — all SERIAL, while the pair loops were the
   one place with OMP.  Landed: `qchem.Parallel::WorkerThreads()` (the one `GPW_OMP_THREADS` knob),
   OMP at the four XC-mesh sites (all partitioned by OUTPUT ⇒ bit-identical at any thread count),
   a triangular H_xc quadrature, the shell-radial hoist in the pointwise basis sweep, a threaded
   batched `SeedCD`, plus the measurement apparatus itself (`GPW_REPORT=1` + per-iteration timing
   buckets) — that is what turned "the runtime gap" into a ranked list.
   **Round 2 (same day):** the `openblas_set_num_threads(1)` pin — lost in the machine migration —
   restored as `qchem::PinBlasToOneThread()` (one home, four mains, loud link failure if the BLAS is
   swapped), and `QCHEM_BLAZE_BLAS` (default ON) routing Blaze's dense products to OpenBLAS: 2.55 →
   41.7 GFlop/s on the XC-mesh shape.  **Durable finding: Blaze dispatches to BLAS only when the
   DESTINATION is a plain matrix** — a `submatrix` destination (or a Hermitian ADAPTOR operand)
   silently falls back at 1.87 GFlop/s, so round 1's hand-blocked triangular quadrature was losing
   13× to save 2×; under BLAS mode it is one whole-matrix product with no OpenMP at all.  ρ sampling
   4.17 → 0.18 s, per-iteration ≈ 8.7 → ≈ 7.6 s, `ctest -j8` 716/716 (sweep 615 → 540 s).
   **Round 2 leftovers, now resolved or re-scoped:** the real-Γ/TRIM path is DONE (doc/RealComplexPlan.md
   — `GpwOptions::realTRIMBlocks` defaults true since `46feb84a`).
   **Round 3 DONE 2026-08-19 (per-iteration SCF 1.42× on the MnO magnetic cell; record in
   doc/GPWPlan1.md "Round 3").**  It closed the flip's last named increment and REFUTED its premise:
   the ~150 s attributed to "the complex-internal collocation streams" is real (collocate 87.5 +
   integrate 41.1 + stream build 23.5) but is DRAM-BANDWIDTH bound, not complex-arithmetic bound — a
   `perf` annotation puts 49% of the scatter on the value + index loads, ~1% on arithmetic, and nothing
   measurable on the Bloch phase or the complex D contract (they are per (pair, offset), amortised over
   runs of 20–40 points; the streams were already `vector<double>`/`<float>`).  So the fix was fewer
   BYTES PER POINT: run-length stream geometry (`runBase`/`runLen`, no per-point index — bit-identical
   replay; 1.44× collocate, 1.41× integrate, and MnO stream RAM 5.78 → 3.70 GB at unchanged coverage,
   which is item 2's lever too), plus the genuinely complex-bound neighbour `FourierMixCD`'s batched
   inverse FT (half-space fold + hoisted (x,y) phase, 2.39× — it was the single largest per-iteration
   bucket at 87.6 s).  `ctest -j8` 747/747, 0 failed.
   **Round 4 STAGE A DONE 2026-08-19 (doc/GPWPlan1.md "Round 4"): the SHELL-BLOCKED box walk, on the
   stream build.**  Reading the kernel showed the charter under-sold it: not only the `exp`s but the
   pair→level assignment, the offset list, the box centre/reach/half-widths, the ellipsoid pre-screen,
   the incremental r walk and the modulo wrap are ALL shell properties (they read `radials[i]` only) —
   components differ only in `pols[i]`/`ns[i]`, and `PGData::Init` already lays a shell out contiguously.
   New `ForShellPairBox` walks the box once per shell pair and hands back the two per-component FACTOR
   arrays, so a d×d shell pair evaluates 2 radials + 12 polynomials per point instead of 72 + 72;
   `ForPairBox` is now its 1×1 case (one box walk in the codebase).  The build also collapsed to ONE
   path: shell-major walking forced the two-pass shape everywhere (the tiering order is row-major over
   COMPONENT pairs and is part of the result), retiring the serial one-pass path's transient-bound abort
   and streaming fp32 demotion.  **Stream build 28.73 → 17.01 s (1.69×); sweep 212 → 203 s (the serial
   path did not regress); stream cache byte-for-byte identical (same pts/offsets/runs/meanRun), physics
   unmoved; 747/747.**
   **Round 4 STAGE B DONE 2026-08-19 — and it RE-SCOPED the charter.**  Same blocking on the uncached
   collocate/integrate arms (cached pairs still replay row-major, so fully-cached runs stay
   bit-identical); the D-aware box resolves as the UNION + per-component screen, with `PairPrefactorExp`
   exposed so each component pair keeps its own WHOLE-TERM prefactor kill (without that the first cut
   measured only 1.03–1.08×).  **Static sharp-field sweeps 2.13× (`LocalPPKappaSelfConverged` 15.9 →
   7.5 s — uniform ε, no D-aware boxes); over-budget SCF only 1.17× (Si) / 1.07× (Mn).**  Energies
   identical to every printed digit, cached vs all-on-the-fly; 747/747.
   **★ THE FINDING THAT MATTERS: the over-budget regime is NOT geometry-bound.**  perf on the
   all-on-the-fly Si run: 61% is the irreducible per-component-pair EMIT (multiply + magnitude screen +
   the scattered `r[idx]+=`), 22% polynomials, 10% `exp` — so only ~32% is what blocking can remove.
   **Run 64's 4318 s was never going to fall by a large factor to a faster kernel; that regime needs
   fewer TERMS, not cheaper ones** (looser `GPW_DENSITY_EPS`, the T3 stream fold, or a smaller span).
   Charter item 2's "fast-recompute kernel unblocks the CP2K-span cell" is hereby narrowed: it delivered
   on the SETUP sweeps (build 1.69×, static field 2.13×) and not on the per-iteration scatter/gather.
   **Still open here:** Φ-table SCREENING (the O(N)-XC increment: batch the mesh, keep a per-batch
   significant-function list — the win grows with cell size) and the Becke mesh build (`BeckeCutoff`
   alone 11.5%).
2. **Close the RAM gap vs CP2K.**  Same lever: recompute fast enough → the stream cache tier
   shrinks → RAM falls with it (CP2K caches nothing; its kernels are just fast).
3. **Understand why CP2K holds the FULL 136-function span and qchem cannot.**  ⛔ **THE SCREEN-DISCIPLINE
   HYPOTHESIS IS REFUTED (run 64, 2026-08-17, doc/logs/mno_probe_run64_fullspan_tighteps.log).**  THE
   EXPERIMENT BELOW WAS RUN: 132-of-136 span (tol 1e-3 dropped 4 AOs), eps 1e-12 both screens, 4 iters
   — and it dove anyway, to E = −82.19, i.e. 20.7 Ha BELOW CP2K's variational −61.47 on a SUBSET of
   CP2K's own span (unphysical, same argument as run 58's −67.28).  Tighter eps is NOT the cure, so the
   next hypothesis has to come from somewhere else: candidates are the SVD/eigen-consistent F/S
   filtering the run-58 note parked, the symmetry-INEQUIVARIANT AO drop (run 64 dropped indices
   11/9/115/1 individually — see the vet-stage item above), and the possibility that the near-null
   directions need projecting out of F as well as S rather than merely being screened around.
   Re-run note: set `GPW_MNO_VERBOSE=1` (GPW_REPORT=1 gives the ledger but NOT the per-iteration table,
   so run 64 cannot show whether it dove or stalled).  Analysis in doc/GPWPlan1.md "Run 64".
   *(original hypothesis, kept for the record)* Hypothesis (banked,
   testable): screen discipline — CP2K's 1e-14 eps keeps the F/S inconsistency below what the
   λ~2e-5 near-null modes can amplify (CP2K itself collapsed 3.5 Ha at loose eps — the retracted
   SR oracle).  THE EXPERIMENT: v2 span, MNO_ORTHO_TOL=1e-3 (zero drops), GPW_SCREEN_EPS=
   GPW_DENSITY_EPS=1e-12, GPW_MNO_NMAX=4 — dive gone ⇒ confirmed, and qchem gains full-136
   capability priced in runtime (= item 1's business).
3b. **REAL vs COMPLEX — the scalar-type plan (NEW 2026-08-16, doc/RealComplexPlan.md).**  Fell out of
   item 1: a Bloch block at a TRIM k (2k ≡ 0 — Γ AND every zone-boundary point, so a Γ-centred
   2×2×2 mesh is TRIM throughout) is exactly real, and so are its H, C, D.  The plan derives the
   rule from the physics rather than declaring it — `block is real ⇔ irrep.IsReal() ∧ (every
   term.PreservesReal())`, each participant answering about itself — covers SOC / B-fields /
   non-collinear / Dirac in one table, and separates BASIS type from WORKING type (with SOC the
   basis stays real and only H/C/D go complex, so the expensive geometry layer must not be dragged
   complex with them).  Step 0 is DONE (`1b8b9a83`, exact ±1 phases ⇒ realness is bitwise).  The
   REST IS DELIBERATELY NOT STARTED: it is gated on the cleanup items that sit on the same faces
   (V1.6, V1.7, V1.8, V1.11, V1.1, V1.5, V1.10, R2.16 — listed in the plan), and one open question
   (the SCF accelerator's per-block history) could negate the memory win if settled late.

4. **Code cleanup batch.**  doc/CleanupCandidates.md D1–D13 + the vet-stage symmetric basis trim
   (this file, above) + Δρ/N convergence gate (doc/SCFStrategyPlan.md) + GDM fallback-diagonalize
   breadcrumb (run 59's silent +302 mHa hop) + per-channel ortho duplication + the fingerprint's
   overconfident "raise NMaxIter" advice.
5. **Finish the symmetry upgrade (doc/SymmetryUpgradePlan.md) — the FOLDS, armed and REPORTED.**
   Inventory 2026-08-15: the {G}-star fold (SymmetrizeGMap/EvaluateSymmetricGMap, cubic-star +
   non-symmorphic unit gates) is wired at exactly TWO static sites (the local-PP sweeps, imposed
   runs only) and reports nothing; the per-iteration G-space consumers (ρ̃, Poisson multiply, V_xc
   gathers, G_ERI3 columns, seed structure factors) are UNFOLDED — 12–24× on the MnO magnetic
   group, 48× cubic, left on the table.  The REAL-space counterpart (T3 pair-stream orbit fold,
   71× measured on the cache) is BUILT but still opt-in (GPW_STREAM_FOLD=1; the open-shell slosh
   retraction).  Work: default-arm T3 where safe + auto-arm at multi-k (T3.4b), extend the ball
   fold to the per-iteration sites, and print a fold-factor line in every grid/stream report
   (the user's standing complaint: the folding is invisible on cout).  Caveats: the FFT itself
   does not fold trivially; the dominant per-iteration cost is real-space, so T3 is the big win.
6. **Then** return to the MnO ground-state accuracy campaign (doc/SphericalLatticePlan.md arm-2
   verdict + the deep-moment basin + MNO_KMESH=2 — k-convergence moves the ordering MORE than the
   physical 6J₁+12J₂ ≈ 4 mHa scale, so it needs items 1–3 first to be affordable).
