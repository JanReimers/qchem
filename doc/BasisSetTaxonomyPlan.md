# BasisSet Taxonomy Plan — V1.33

**Status: RECORD — executed in full 2026-09-13 (eleven commits `c2cb79a3`..`d5ddb1a5`, 856/856).**  The rulings in §1 stand as
the design; the per-step notes below are the execution log.  Closed row: `CleanupHistory.md` "LANDED 2026-09-13 — V1.33".  Executes `CleanupCandidates.md` V1.33 ("the BasisSet
taxonomy is on the wrong axis").  Read §1 once; it is the ruling.  §4 is the running order.

Two proposals preceded this plan and they were NOT in conflict — they were the two axes:

- 2026-08-20 (tracker row): classify by *basis kind* — `Radial / Polarized{Cartesian|Spherical} / GPW / PW`.
  That is the **FAMILY** axis.
- 2026-09-13 (this session): classify by *symmetry* — `Translational / Spherical / Point / PointLattice`.
  That is the **GROUP** axis.

The resolution: **library boundaries follow the family (the integral engine); module names carry the
group.**  Each axis gets exactly one level of the hierarchy and neither leaks into the other.

---

## 1. The ruling (pins)

### 1.1 An IrrepBasisSet is a carrier space for ONE irrep of G

`BasisSet = ⊕_irreps IrrepBasisSet`, where the irrep is a label of **G, the symmetry group of the
Hamiltonian**.  "Categorise by symmetry" therefore means "by the group whose irreps label the blocks".
Two consequences that dissolve most of the historical confusion:

1. **Carriers are not invariants.**  \f$Y_{lm}\f$ is not spherically symmetric, \f$e^{i\mathbf{k}\cdot\mathbf{r}}\f$
   is not translation-invariant, a SALC is not point-group-invariant.  Each *transforms as* an irrep.
   Invariance was never the criterion.
2. **k is an irrep label of the translation group T**, exactly as \f$l\f$ labels O(3) irreps or \f$A_1\f$ labels
   \f$C_{2v}\f$ irreps.  The code already says so: `tGPW_IBS(const UnitCell&, const sym_t& irrep, ...)` — "the
   Bloch symmetry IS the k-label".  k is neither "in the basis" nor "in the eigenfunctions"; it is the
   block label, and the eigenfunctions inherit it because H is block-diagonal in it.  Bloch's theorem is the
   statement that the projection operator \f$P_k=\sum_R e^{-i\mathbf{k}\cdot\mathbf{R}}T_R\f$ block-diagonalises H.

### 1.2 Three orthogonal axes and a role

| axis | values today | lives in |
|---|---|---|
| **G** — the block labels | \f$1\f$, P (point group), O(3), O(3)\* (double group, Dirac), T (translations); later T⋊P (space group), P\*, T^d (d=1,2) | `qcSymmetry/{Atom, Molecule, Lattice_3D}` — already laid out by G |
| **Family** — the analytic seed | Gaussian {primitive, contracted}×{Cartesian, spherical}; Slater; BSpline; exponential (PW); delta; composites (APW/LAPW) | the evaluators (engines) |
| **Construction** — seed → irrep carrier | DERIVED from (G, family); takes only two values (§1.3) | `qcSymmetry` (SALC, Bloch phases, Fold) |
| **Role** | Orbital_{1E, HF, DFT, DHF, PP}; Fit_{CD, SF} | the `qcBasisSet` core faces — already orthogonal |

Because the axes are orthogonal, **a family belongs to no G-category**.  `Atom/Evaluators/Gaussian` and
`Molecule/PG_Spherical` are the SAME family with two engines; the atomic engine exists because O(3)
collapses every integral to a radial one.  "Gaussians are inefficient for atoms via the molecular engine"
is an engine property, not a family property.

### 1.3 The construction axis is not free

- **Subduction (A):** the family is already irrep-adapted to a bigger group \f$G_{big}\supseteq G\f$; block it by the
  branching rules \f$G_{big}\downarrow G\f$.  \f$Y_{lm}\f$ are O(3) carriers → in a molecule subduce O(3)↓P — that is
  `Symmetry::Molecule::{Cartesian,Spherical}ShellRep`.  \f$e^{i\mathbf{q}\cdot\mathbf{r}}\f$ are \f$\mathbb{R}^3\f$
  carriers → in a crystal subduce \f$\mathbb{R}^3\downarrow T\f$ — that is the \f$\mathbf{k}+\mathbf{G}\f$ labelling.
- **Induction (B):** the seed is site-local; apply G's projection operators to its orbit.  AO shells on
  equivalent atoms → SALCs (site group ↑ P).  An AO at a cell atom → the Bloch sum
  \f$\sum_R e^{i\mathbf{k}\cdot\mathbf{R}}\phi(\mathbf{r}-\mathbf{R})\f$ (\f$1\uparrow T\f$).  Deltas on a uniform grid → Bloch
  sums of deltas = plane waves restricted to the grid: **the uniform delta basis is the PW basis's adjoint**
  because both carry T, one induced, one subduced, and the DFT is the unitary between them.

Molecular Gaussians are B∘A (∘ is composition, read right to left: FIRST subduce the shell to the site group,
THEN induce the site up to P).  GPW is B_T∘(B_P∘A): Bloch-sum the molecular basis; for the full space group,
induce once more over the little group.  (Notation fixed to the standard right-to-left reading 2026-09-14,
when the Doxygen page `\ref basisset_taxonomy` in `src/BasisSet/BasisSet.C` defined it.)

### 1.4 Spin is a factor of G, until it is not

`Symmetry::Irrep = (sym_t spatial) × (Spin ms)` encodes \f$G = G_{spatial}\times SU(2)\f$: the non-relativistic
case, and the REASON Pol/UnPol is a block structure at all — the two spin blocks are the SU(2) factor's
irreps.  Dirac's \f$\Omega_{\kappa m_j}\f$ is not an odd man out; it is the row where G is a **double group** and the
product does not factor.  Consistent with `Pins.md` spin-native: spin-native = "G contains SU(2)"; UnPol =
the ζ=0 collapse of that factor; SOC/non-collinear = the factor dissolving into G.  This PREDICTS two rows
we do not have: P\* (4-component molecular: RKB Gaussians × double point group) and T⋊P\* (crystals with SOC).

### 1.5 "No symmetry" and non-uniform meshes

\f$G=1\f$ is a group with one irrep — the trivial row, not a category (non-symmetric molecules already sit
there).  A Becke grid on a repeating cell is \f$G=T\f$, delta family, construction B.  "Translational in real
space but not in reciprocal space" is a FAMILY statement (a non-uniform mesh has no G-ball dual), not a
symmetry one.

### 1.6 The two placement rules

1. **A LIBRARY is an engine.**  A basis library's mass is its integral engine, and an engine factorises
   over the family, not over G.  `UnitCell` inside the Gaussian engine is therefore legitimate: a lattice
   sum is a Gaussian-engine operation.  The GPW seam ("GPW is a new evaluator, not a new IBS"; the image
   sum pushed into the molecular seam so a transformed basis transforms once per point — the 4.8 s/iter
   sharing pin) STAYS inside the Gaussian engine.  A library cut on the G axis is ruled out.
2. **A MODULE carries the G.**  `qchem.BasisSet.Gaussian.Point.*` vs `qchem.BasisSet.Gaussian.Lattice.*`.
   The V1.33 evidence ("`qchem.UnitCell` imported inside `Molecule/`") becomes a compiler-adjacent
   invariant: *no `.Point.` module imports `qchem.UnitCell` or `qchem.Symmetry.Lattice_3D.*`* — one grep,
   run by ctest (§4 step 3).

Corollary — **the container layer.**  Each G has a ⊕_irreps container + factory
(`Atom/BasisSet+Factory`, `Molecule/SymmetryAdaptedBasisSet+Factory+Readers`, `Lattice_3D/BasisSet+
BandStructure+Factory`).  It is family-agnostic in principle.  **It becomes its own library only when it must
see ≥2 engine libraries** — true today only for T (PW and Gaussian engines, and APW/LAPW will add the
radial one).  For O(3) and P it stays a module tier inside the single engine library.

### 1.7 The four tiers, and what each one is tagged with

The atom library already shows the finished shape (`Atom/IrrepBasisSet.C`): every IBS flavour is a thin class
templated on an injected evaluator, mixins combine the radial / angular / non-rel / relativistic aspects, and
C++ **concepts** are the brief spec — *"meet this spec and everything just works"*.  Read through the axes,
that is four tiers, each tagged by exactly one axis:

| tier | tagged by | what it is | today |
|---|---|---|---|
| **role faces** | role | `Orbital_1E_IBS<T>`, `Orbital_DFT_IBS<T,TFit>`, `Fit_IBS`, … | `qcBasisSet` core |
| **spec** | G | the concepts an engine must meet to drive one G's blocks, plus the evaluator-injected mixins that turn a conforming engine into the role faces | atom: `Atom/Evaluators/Evaluator.C` + `Atom/IrrepBasisSet.C`; lattice: `Lattice_3D/IrrepBasisSet.C` + the concepts in `Evaluators/PW/Evaluator.C` |
| **engine** | family | the integral machinery; knows nothing of the spec, merely satisfies it | the `Evaluators/` trees |
| **IBS** | (G, family) | one thin class per pair: `Lattice::Orbital_1E_IBS<GPW_Evaluator, double>` reads as (G, role, family, scalar) — the whole taxonomy in one type name | `GPW_IBS`, `PlaneWave_IBS`, the atom `*_IBS<E>` |

Rules that follow:
- **The spec is engine-free, so it lives BELOW every engine — in the core.**  Measured: the lattice mixins
  import only core faces + the concept, and the lattice concepts are purely structural (no
  `derived_from`), which is exactly why `GPW_Evaluator` satisfies `isPW_1E_Evaluator` without inheriting
  from `PW_Evaluator`.  The atom concepts DO say `std::derived_from<E, Evaluator>` — that is what keeps the
  radial spec engine-bound today, and it is the one thing to loosen when the pattern is carried to molecules.
- **A concept is named for the G it serves, never for a family.**  `isPW_1E_Evaluator` is satisfied by a
  Gaussian engine; its content (complex Bloch blocks, `chmat_t`, `NuclearMatrix(cl)`) is a *lattice* spec.
  Rename: `isLattice_1E_Evaluator` / `isLattice_DFT_Evaluator`; mixins `Lattice::Orbital_1E_IBS<E,T>` etc.,
  mirroring `Atom::Orbital_1E_IBS<E>`.
- **A family assumption inside a spec is a defect.**  `isPW_DFT_Evaluator` demands
  `OverlapMatrix(std::function<dcmplx(const ivec3_t&)>)` — a potential→matrix bridge by G-vector LOOKUP.
  That is the PW *fit-family* pairing leaking into the G spec (the same frozen pairing as §5); it belongs to
  the fit side of `Orbital_DFT_IBS<T,TFit>`, and the audit at step 3 should flag any other such term.

---

## 2. Placement table — everything we have, and what we can anticipate

| basis | G | family | construction | role |
|---|---|---|---|---|
| Slater / Gaussian / BSpline × \f$Y_l\f$, \f$Y_{lm}\f$ | O(3) | radial | A | Orbital, Fit |
| RKB{Slater,Gaussian,BSpline} × \f$\Omega_{\kappa m_j}\f$ | O(3)\* | radial (RKB) | A | Orbital (DHF) |
| PG_Cart / PG_Spherical / PG_LibCint | P (incl. \f$1\f$) | Gaussian | B∘A | Orbital, Fit |
| GPW | T (T⋊P later) | Gaussian | B_T∘B_P∘A | Orbital |
| PW | T | exponential | A | Orbital |
| PlaneWaveFit | T | exponential | A | Fit |
| uniform delta (Vxc fit) | T | delta | B (= PW's adjoint) | Fit |
| Becke delta, molecule | \f$1\f$ or P (site-adapted invariant mesh) | delta | B | Fit |
| Becke delta, periodic | T | delta | B | Fit |
| **anticipated** | | | | |
| numerical atomic orbitals (SIESTA / FHI-aims) | any | radial {…, Numeric} × \f$Y_{lm}\f$ | as Gaussians | Orbital |
| APW / LAPW | T | **composite**: PW outside ⊗ radial×\f$Y_{lm}\f$ inside | A + sphere matching | Orbital |
| molecule in a PW box (CP2K style) | T (large cell) | exponential | A | Orbital |
| 2D slabs / 1D polymers | T², T¹ (+P) | any | unchanged | — |
| 4c molecular / SOC crystals | P\*, T⋊P\* | RKB × \f$\Omega\f$ | B∘A | Orbital |
| Wannier functions | T (inverse: k→R) | derived | B⁻¹ | — |
| multiwavelets / FEM | \f$1\f$ or T | delta-like | B | Fit, Orbital |
| DFT+U / KB projectors | — | radial family, subduced | — | a ROLE (`Orbital_PP_IBS`), not a basis |

Nothing needs a fourth axis.  The only genuinely new KIND of entry is the composite family (APW), and it
is what forces the T container layer to be its own library (rule 1.6 corollary).

### 2.1 Acceptance, on paper (step 4, 2026-09-13) — where each anticipated row LANDS in the finished tree

| row | library / module | tier | verdict |
|---|---|---|---|
| **APW / LAPW** | `qcLattice_BS` (`Lattice/{APW,LAPW}_IBS.C`, already there) — the composite is (PW outside) ⊗ (radial×$Y_{lm}$ inside), so the IBS lives in the T CONTAINER and pulls BOTH engines: `qcLattice_BS → qcPlaneWave_BS` (today) **and `→ qcRadial_BS`** (the day the sphere interior stops being hand-rolled: `LAPW_IBS.C` currently owns its own radial solver; when it asks `Radial.*` for $u_l(r)$ the edge appears, and NO new library does).  The sphere-matching is construction A + a boundary condition — a property of the composite family, not a new axis. | IBS (T, composite) | fits; forces the one edge the plan predicted |
| **numerical atomic orbitals** (SIESTA / FHI-aims) | a NEW engine library `qcNumeric_BS` (family = numeric radial × $Y_{lm}$ on MULTIPLE centres), with `Numeric.Point.*` / `Numeric.Lattice.*` tags exactly like `Gaussian`.  **NOT** a fourth `Radial.*` sub-directory reused by `Gaussian.Point`: `qcRadial_BS` is the O(3) engine — one centre, every integral radial — and a multi-centre numeric basis has none of that machinery (its integrals are two-centre tables + grids).  The single-atom NAO GENERATOR is a `Radial.Numeric.*` row (like `valgen`); the basis it emits is consumed by the multi-centre engine — two libraries, the same way `Radial.Gaussian` and `Gaussian.*` are two today. | engine (family) | fits; the point of the exercise |
| **double groups** (P\*, T⋊P\*, O(3)\* beyond the atom) | `qcSymmetry` rows (`Symmetry::{Molecule,Lattice_3D}::DoubleGroup…`), consumed through the SAME `Irrep` currency the containers already take.  Basis side: the RKB family gains `Gaussian.Point.RKB_*` / `Gaussian.Lattice.RKB_*` IBS classes (family = RKB Gaussian, G = P\*/T⋊P\*), no new library — the engine is still Gaussian.  Spin dissolving into G (§1.4) changes the BLOCK LABEL, not the axis count. | G (symmetry) | fits; a Symmetry row, not a BasisSet row |

None of the three needs a fourth axis; the plan stands and step 5 may run.

---

## 3. Library map — current → target

Dependency direction today (measured 2026-09-13): the PW half of `Lattice_3D` imports nothing from
`Molecule`; the GPW half imports `Molecule.LatticeSum1E` AND the PW half (`Evaluators.PW`,
`PlaneWaveFit_IBS`, `PeriodicGridEvaluator`) — GPW = Gaussian orbitals × PW density fit.  So the layering
below has no cycle: PlaneWave is a leaf over the core, Gaussian depends on PlaneWave, Lattice on both.

| target library | target module prefix | contents (from today's tree) | links |
|---|---|---|---|
| `qcBasisSet` (+ the lattice spec) | `qchem.BasisSet.*` | the role faces, caches, `DeltaFit_IBS`, `Projector3`, `SpeciesField`, `GMap`; **NEW: `qchem.BasisSet.Lattice_IBS`** — the G-tagged spec tier (§1.7): `isLattice_{1E,DFT}_Evaluator` + the `Lattice::Orbital*_IBS<E,T>` mixins (today `Lattice_3D/IrrepBasisSet.C` + the two concepts in `Evaluators/PW/Evaluator.C`) | qcElConfig qcStructure qcMesh |
| `qcRadial_BS` (= `qcAtom_BS`) | `qchem.BasisSet.Radial.*` | ALL of `Atom/` — the O(3) engine; the one place the group IS the engine.  `Radial.{Slater,Gaussian,BSpline}.*`, the shared angular integrals, the container + factory tier | qcBasisSet |
| `qcPlaneWave_BS` (new, split out of `qcLattice_BS`) | `qchem.BasisSet.PlaneWave.*` | `Evaluators/PW`, `Evaluators/Imp/PeriodicGridEvaluator`, `PlaneWave_IBS`, `PlaneWaveFit_IBS`, `Internal/{GVectors,KPlusG}` — the engine only; the spec it satisfies moves to the core (§1.7) | qcBasisSet qcSymmetry |
| `qcGaussian_BS` (= `qcMolecule_BS` + the GPW half) | `qchem.BasisSet.Gaussian.*` with the G as a sub-tag: `Gaussian.Evaluators.*` (engine, no G — serves both), `Gaussian.Point.*` (today's `Molecule.{IBS, PG_Cart, PG_Spherical, PG_LibCint, SymmetryAdaptedBasisSet, Factory, Reader(s), BasisFiles}`), `Gaussian.Lattice.*` (`LatticeSum1E`, `LatticeScreener`, `PG_Spherical.LatticeView`, `Evaluators/GPW`, `GPW_IBS`) | qcBasisSet qcPlaneWave_BS qcSymmetry cint |
| `qcLattice_BS` (thin; the T container layer) | `qchem.BasisSet.Lattice.*` | `Lattice_3D/BasisSet` (the ⊕_k container + `Factory`/`GPWFactory`/`GPWParams`), `BandStructure`, `APW_IBS`, `LAPW_IBS` | qcPlaneWave_BS qcGaussian_BS (qcRadial_BS when APW is real) qcLASolver |

Naming notes:
- `Radial` rather than `Spherical`: the engine is the radial reduction, and `Spherical` collides with
  `PG_Spherical` (spherical-harmonic Gaussians), which would be a standing confusion.
- `Lattice` rather than `Lattice_3D`: T^d is a row in §2.  `qcSymmetry/Lattice_3D` keeps its name for now —
  V1.33's scope precision stands: `Symmetry::Molecule` vs `Symmetry::Lattice_3D` is a REAL distinction, and
  the symmetry library's rename is a separate decision.
- The C++ namespaces (`BasisSet::Atom`, `::Molecule`, `::Lattice_3D`) follow the module names in the same
  step; the external blast radius is small (6 / 11 / 8 files outside `src/BasisSet`).

Blast radius outside `src/BasisSet` (measured): 14 `.C` files import the three sub-library module
prefixes; 38 `CMakeLists.txt` name the three library targets (mostly test links); `pybind/` — check at
step 2, flag to the binding owner, do not edit.

---

## 4. Running order

Every step lands green on `ctest -j8` (full suite), bit-identical energies on the GPW/PW/atom anchors, and
is a commit on its own.  Steps 1–3 are mechanical; nothing in them changes a number.

**Step 0 — pin the ruling.**  This file.  Add the §1.6 rules to `CleanupCandidates.md` R1.0 (the durable
design rulings).  ✅ 2026-09-13 (this file written; R1.0 entry pending the user's read).

**Step 1 — split `qcLattice_BS` along the engine seam** (no renames yet).
- 1a0. Move the spec to the core: `Lattice_3D/IrrepBasisSet.C` + the two concepts → `qchem.BasisSet.Lattice_IBS`
  (namespace `qchem::BasisSet::Lattice`), concepts renamed `isLattice_{1E,DFT}_Evaluator`, mixins renamed
  `Lattice::Orbital_1E_IBS<E,T>` / `Lattice::Orbital_DFT_IBS<E,T>` / `Lattice::Irrep_IBS<E,T>` (mirroring the
  atom names).  Dependency-free by measurement (§1.7).  `static_assert(isLattice_1E_Evaluator<GPW_Evaluator>)`
  at each concrete IBS — the spec is checked where engine meets spec, never inside the engine.
  ✅ 2026-09-13.  `src/BasisSet/Lattice_IBS.C` (git-mv'd from `Lattice_3D/IrrepBasisSet.C`, so `--follow`
  survives), `qcBasisSet` FILE_SET; the two concepts left `Evaluators/PW/Evaluator.C` and the two
  `static_assert`s left `Evaluators/GPW/Evaluator.C` for `GPW_IBS.C` (+ new ones beside `PlaneWave_IBS`).
  `isLattice_DFT_Evaluator` does NOT carry the old `OverlapMatrix(std::function<dcmplx(const ivec3_t&)>)`
  term — the §1.7 family leak; no mixin ever consumed it (the PW-term `applyAdjoint` closures reach the
  method directly), so dropping it changed no call.  Both engines still satisfy the spec structurally; the
  GPW engine module no longer names any concept.  Full `ctest -j8` green, no number touched (no numeric
  code moved — a pure relocation + rename).
- 1a. New library `qcPlaneWave_BS` from the PW half listed in §3.  `qcLattice_BS` links it.  Gate: build +
  `PlaneWaveDFTUT`, `GPW_UT`, `GPW_SCF_UT` unchanged.
  ✅ 2026-09-13.  `src/BasisSet/PlaneWave/{Evaluators/{Evaluator,PeriodicGridEvaluator}, PlaneWave_IBS,
  PlaneWaveFit_IBS, Internal/{GVectors,KPlusG}}` (all `git mv`), target `qcPlaneWave_BS` linking only
  `qcBasisSet qcStructure qcSymmetry`; `qcLattice_BS` links it PUBLIC so every client reaching those modules
  through `qcLattice_BS` is unchanged.  Module names/namespaces keep the `Lattice_3D` spelling until step 2.
  `ldd libqcPlaneWave_BS.so` already shows the 1c property: no `qcMolecule_BS`, no `qcLattice_BS`.
  ctest -j8 855/855.  Left where they were: `Lattice_3D/tests/{PlaneWaveUT,GMapUT,...}` — the PW-only
  unit tests still ride `UTLattice_3D_BS`; splitting a `UTPlaneWave_BS` off is the TE (test-suite axes)
  programme step's business, not a library-split step's.
- 1b. Move `Lattice_3D/Evaluators/GPW` + `Lattice_3D/GPW_IBS.C` into `qcMolecule_BS` (physically under
  `Molecule/` for now; they become `Gaussian/Lattice/` at step 2).  `qcMolecule_BS` links `qcPlaneWave_BS`.
  Gate: same + `A_PP`, `L_PP`, `RealComplexTermsUT`.
  ✅ 2026-09-13.  `src/BasisSet/Molecule/Lattice/{GPW_Evaluator,GPW_IBS}.C` (+ `Imp/`), all `git mv`
  (`Evaluators/GPW/Evaluator.C` → `GPW_Evaluator.C` so the name survives without its directory);
  `qcMolecule_BS` links `qcPlaneWave_BS` PUBLIC.  `Lattice_3D/Evaluators/` is gone.  Module names and
  the `BasisSet::Lattice_3D` namespace unchanged until step 2.  `ldd`: PlaneWave → {core}; Molecule →
  {core, PlaneWave}; Lattice → {core, PlaneWave, Molecule, LASolver} — exactly §3's layering.
  ctest -j8 855/855 (A_PP 1, L_PP 3, RealComplexTerms 10, PlaneWaveDFT 29, GPW 29, GPW_SCF 36).
- 1c. What is left in `Lattice_3D/` is the container tier (`BasisSet`, `BandStructure`, `APW_IBS`,
  `LAPW_IBS`); `qcLattice_BS` links both engines.  Confirm with `ldd`/link order that `qcPlaneWave_BS`
  pulls in nothing from `qcMolecule_BS`.
  ✅ 2026-09-13.  `Lattice_3D/` = `{BasisSet, BandStructure, APW_IBS, LAPW_IBS}.C` + `Imp/` + `tests/`.
  `ldd libqcPlaneWave_BS.so` → core only; `nm -DC --undefined-only` has no `Molecule::`, `GPW`, or
  container-tier symbol.  Step 1 complete — three commits (c2cb79a3, 8104e5fa, e4cc3370).

**Step 2 — rename onto the two axes.**  Directories, library targets, module names, C++ namespaces, in
one commit per library so `git log --follow` survives:
- `Atom/` → `Radial/`, `qcAtom_BS` → `qcRadial_BS`, `qchem.BasisSet.Atom.*` → `qchem.BasisSet.Radial.*`.
- `Molecule/` → `Gaussian/` with the `Point/` and `Lattice/` sub-directories of §3, `qcMolecule_BS` →
  `qcGaussian_BS`, `qchem.BasisSet.Molecule.X` → `qchem.BasisSet.Gaussian.Point.X` or `.Lattice.X`;
  `Molecule.Evaluators.*` → `Gaussian.Evaluators.*` (no G tag — the engine serves both).
- `Lattice_3D/` → `Lattice/`, `qchem.BasisSet.Lattice_3D.*` → `qchem.BasisSet.Lattice.*`.
- `.vscode/settings.json` TestMate exe names and the root `allTests` DEPENDS list are re-checked (the
  `UTSCFAccelerator` lesson: ctest's N must not drop).
- `pybind/`: report any breakage to the binding owner; do not edit.
  ✅ 2026-09-13, four commits (one per library, all `git mv` so `--follow` survives):
  - **Radial** f5461b8f: `Atom/`→`Radial/`, `qcRadial_BS`/`UTRadial_BS`, `qchem.BasisSet.Radial.*`,
    `BasisSet::Radial`.  (`Symmetry::Atom` and the `qchem::Atom` class untouched.)
  - **PlaneWave** 1f602b2f: modules `qchem.BasisSet.PlaneWave.{Evaluators, Evaluators.PeriodicGrid,
    PlaneWave_IBS, PlaneWaveFit_IBS, Internal.*}`, namespace `BasisSet::PlaneWave`.  `GPW_Evaluator`
    re-exports `RasterPolicy`/`RasterFields` (`using PlaneWave::…`) so a GPW client names the knobs where it
    names the block.
  - **Gaussian** b3a7b7c7: `Molecule/`→`Gaussian/` with `Evaluators/` (no G), `Point/` (IBS mixins, PG_*,
    SALC container, factory, readers) and `Lattice/` (LatticeSum1E, LatticeScreener,
    `SphericalLatticeView` ← `PG_Spherical.LatticeView`, GPW); `qcGaussian_BS`/`UTGaussian_BS`; namespace
    `BasisSet::Gaussian` for the whole library (the G is the MODULE tag; a `Gaussian::Lattice` sub-namespace
    would shadow the core `BasisSet::Lattice` spec).  CMake: the module file set is declared at the library
    root with `BASE_DIRS` so the three sub-directories contribute.
  - **Lattice** (this commit): `Lattice_3D/`→`Lattice/`, `qchem.BasisSet.Lattice.*`, `UTLattice_BS`;
    namespace `BasisSet::Lattice_3D`→`BasisSet::Lattice`, which MERGES the container tier with the 1a0 spec
    tier — the same arrangement as `BasisSet::Radial` (mixins + concrete IBS in one G namespace).  No
    name-hiding surfaced (the container never names a core face unqualified).
  - `pybind/qchem_bridge.cpp` imports `qchem.BasisSet.Molecule.Factory` / `BasisSet::Molecule::Factory` ⇒
    WILL break under `-DQCHEM_PYBIND=ON` — **flagged for the binding owner**, not edited.
  - Naming wart noted, not fixed (not mechanical): the radial engines' `class Radial` now reads
    `BasisSet::Radial::Evaluators::{Gaussian,Slater,BSpline}::Radial`; and `BasisSet::Radial::Evaluators::Gaussian`
    already coexisted with the math namespace `qchem::Gaussian` (the code writes `::qchem::Gaussian::`), now
    joined by `BasisSet::Gaussian`.  Compiles cleanly; a reader's wart only.
  - ctest N = 855 through all four (no exe dropped; 0 NOT_BUILT).

**Step 3 — make the G tag a checked invariant.**  `scripts/audit-basisset-gtags` in the style of
`scripts/audit-internal-reexports`, run by ctest: a module named `qchem.BasisSet.Gaussian.Point.*` must not
import `qchem.UnitCell`, `qchem.Lattice_3D`, or `qchem.Symmetry.Lattice_3D.*`; a module named
`qchem.BasisSet.Radial.*` must not import either of those nor `qchem.Symmetry.Molecule.*`.  This is
V1.33's original evidence turned into a test.  Any `.Point.` module that fails the audit at this step is a
REAL misplacement to be moved into `.Lattice.` — expected candidates: `PG_Cart/BasisSet.C`'s
`CollocateDensity` grid↔cell map and `PG_Cart/Imp/IrrepBasisSet.C`'s "PERIODIC caller" path.  Decide each
on its merits (move the capability to a `.Lattice.` face, or accept it as engine-level and untag it).
  ✅ 2026-09-13.  `scripts/audit-basisset-gtags` (interface AND implementation units; `Point.*` may not
  import `qchem.{UnitCell, Lattice_3D, ReciprocalLattice, Symmetry.Lattice_3D.*, BasisSet.{Lattice,
  PlaneWave, Gaussian.Lattice}.*}`; `Radial.*` additionally not `qchem.Symmetry.Molecule.*` nor
  `qchem.BasisSet.Gaussian.*`), ctest `BasisSetGTagAudit` (N 855→856), negative-tested on synthetic
  modules.  It fired on exactly the two predicted sites, which are ONE thing: `PG_Cart::Orbital_IBS`
  IS-A `Periodic_Gaussian_IBS` and forwards `LatticeSum1E` to its `PG_Cart_MnD` engine base.
  **Ruling: untag it** — `Point/PG_Cart/` → `Gaussian/PG_Cart/`, `qchem.BasisSet.Gaussian.PG_Cart.*`.
  Not a compromise: the PG_Cart block is the family's raw AO block, the **G=1 SEED** (§1.5: G=1 is a row)
  that both constructions act on from OUTSIDE it — `GetAoShells` feeds the P induction in
  `Point/SymmetryAdaptedBasisSet`, `LatticeSum1E` feeds the T induction in `Lattice/GPW`.  Both are
  questions about the seed.  Option (a), a `.Lattice.` subclass, would need the Point factory to construct
  it (a construction-path change, not mechanical) and is exactly what the §5 sequel's per-(G, engine) thin
  classes do properly.  Open question for that sequel: `PG_Spherical` / `PG_LibCint` are seeds by the same
  argument (they pass the audit only because their lattice ability rides `SphericalLatticeView`); untag
  them when the sequel splits the tiers, not before.  `PG_Spherical`, `PG_LibCint`, the IBS mixins, the
  SALC container, factory and readers stay `Point.*`; nothing in `Radial.*` fired.  ctest -j8 856/856.

**Step 4 — acceptance, on paper.**  Append to §2 a one-line placement for each of: APW/LAPW (the composite
forces `qcLattice_BS → qcRadial_BS`), a numerical-radial family (a fourth `Radial.*` sub-directory, then
reused by `Gaussian.Point`? — NO: NAOs are their own multi-centre engine; that is the point of the exercise),
the double groups (a `Symmetry` row, not a `BasisSet` row).  If any of them needs a fourth axis, the plan is
wrong and comes back here before step 5.
  ✅ 2026-09-13 — §2.1 above.  All three fit; the plan stands.

**Step 5 — tracker hygiene.**  V1.33 ✅ → `CleanupHistory.md` the day it closes; one-line stub here and in
`CleanupCandidates.md`; `README.md` moves this file to RECORD.
  ✅ 2026-09-13.

---

## 5. Out of scope, noted so they are not lost

- **The GPW pairing is frozen in the engine.**  `GPW_Evaluator` owns its PW density-fit grid, so
  `qcGaussian_BS → qcPlaneWave_BS` is a hard-coded (orbital family × fit family) pairing.  `Pins.md`
  "everything is a fit — pairings are POLICY, never hard-coded" says this should eventually be assembled
  at the factory (the container layer, where pairing belongs).  Not this plan; it is the R1.0h per-iteration
  scope's neighbour.
- **Carrying the evaluator-injection pattern to molecules** (`Point_IBS` spec: `isPoint_{1E,DFT,HF}_Evaluator`
  + `Point::Orbital*_IBS<E>` in the core, with `PG_Cart / PG_Spherical / PG_LibCint` collapsing to one thin
  class per engine) is the user's stated goal and the natural sequel; it also loosens the atom concepts'
  `derived_from<E, Evaluator>` so the radial spec can move to the core too.  Own plan, after this one.
- **`qcSymmetry` mirrors the G axis already** (`Atom / Molecule / Lattice_3D`).  Whether it should be
  renamed to the group names (`O3 / Point / Lattice`) is a separate, cheaper decision.
- **T⋊P (space-group irreps as block labels)** is the next G row; it lands in `qcSymmetry.Lattice_3D`
  + the `Gaussian.Lattice` container, and touches nothing in this plan's library boundaries.
