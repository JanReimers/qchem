# D-BUILD-BLAST: the import-graph evidence (RECORD; the open work is CleanCode row D-BUILD-BLAST)

Tools: `scripts/import-blast` (blast radius per module, `--module`, `--why`, `--diff A..B`), `scripts/stale-imports` (text heuristic; confirm with `scripts/importprobe`).

## Item 1: `scripts/import-blast` (lane B 2026-11-18)
**What it reads:**
- CMake's module scan writes one P1689 `.ddi` per TU (provides / requires). A TU that provides a module interface or partition writes its BMI.
- The compile rules carry NO `restat`. So a recompiled interface unit rewrites its BMI and recompiles every importer. The rebuild a change causes is the closure over ALL imports (`import` and `export import` alike), not only re-exports.
- Headers (forward.H, a few quoted includes) come from `ninja -t deps`.
- Orphan `.ddi` files of deleted or renamed sources are skipped. 26 such orphans had inflated the first count: 833 TUs / 367 modules became 811 / 355.

**Validation on a real build:** the merge of main into lane-b, feb17ccc6 (150 files, 17 interface units touched, one new module, a CMakeLists change).
- `--diff abacd64c9..feb17ccc6` predicted **313** TU compiles. ninja compiled **313**, the SAME set: every compiled TU was predicted, and none were predicted that were not compiled.
- With the orphans still counted, the first prediction read 318. The 6 over-predicted TUs were all orphans.
- A CMake re-run does not disturb the count, which a `ninja -n` dry run cannot say: its dry run plans only the regeneration step.
- **Caveat:** `--diff` uses the build tree's CURRENT graph. A commit that adds a module or changes imports is predicted on the graph before it. The new files are listed apart, and a changed import list moves the answer by its own few edges.

**Top of the table** (blast = TUs a change to the interface recompiles; 811 TUs):
- The leaf utilities are at 88%: Factorials, IntPow, Constants, CMath, Vector3D, Matrix3D, Types.
- Streamable 546.
- **FFT 490, NUFFT 488, PointSampling 486**, all through Mesh.
- Blaze 478, Structure 448.

## Item 3 (1): FFT under Mesh
- `qchem.Mesh` re-exports `qchem.PointSampling` legitimately: a Mesh IS a PointSampling.
- But `src/Math/PointSampling.C` held its implementation in the interface unit: `CachedBound`, the concrete `NUFFTSampler` / `RasterSampler`, and the `Make*` bodies. Those are the only users of `import qchem.NUFFT; import qchem.FFT;`.
- Moving that tail to `Imp/` cuts the two edges. Simulated on the graph (orphans still counted): FFT 514 → 5, NUFFT 512 → 1.
- **DONE (lane-b f04412464):** the tail is in `src/Math/Imp/PointSampling.C`, and the interface imports only `qchem.Types`. Measured on the real graph: **FFT 490 → 6, NUFFT 488 → 2**. `import-blast --diff WORKTREE` predicted 488 compiles for the move, and ninja built 488.

## Item 2 (the stale-import heuristic, row item 3): `scripts/stale-imports`
**Method:**
- For each `import X;` it checks whether any name X makes visible (its declarations at namespace / class scope, with function bodies skipped, plus what X re-exports) appears in the TU's code.
- `--verify SNAPSHOT` checks each suspect with the compiler: one removal, `-fsyntax-only`, importprobe's command, in a frozen `cp -al` copy of the build. It cannot be under /tmp, which is a tmpfs.

**Measured on lane-b feb17ccc6:**
- The first pass found 91 suspects, of which the compiler confirmed 86. All 5 misses were blind spots:
  - 4 × `qchem.Strings`, whose arrays and non-ASCII names the patterns missed. Fixed.
  - 1 × `qchem.Blaze`, an ADL operator.
- Now: **87 suspects, 86 confirmed** (precision 99%).
- Modules that export operators carry a ⚠ in the report, because an ADL use is invisible to a name match.
- RECALL is not measured: a stale import that shares a name with something else is kept.
- The 86 confirmed are listed by `scripts/stale-imports`. The widest:
  - 15 × `qchem.Structure` (blast 448);
  - 11 × `qchem.Matrix3D`;
  - 7 each of `qchem.CMath` and `qchem.Types`.
- None was removed. Removal needs the RelChecked build, because an import used only inside `assert()` looks stale to NDEBUG. The heuristic sees asserts as uses, but `--verify` compiles under NDEBUG.

## Item 3 (2): the BasisSet `.Internal.` cache fan-in, the `Cache()` INVENTORY (no code changed; awaiting the user's ruling)
**The carrier:**
- `qchem.BasisSet.IrrepBasisSet` imports `qchem.BasisSet.BasisCache` (blast 381 → 382) only because `IrrepBasisSet_IDs` declares `protected: virtual const BasisCache& Cache() const=0;`.
- A class owned by another module cannot be forward-declared, so naming it forces the import.
- BasisCache must import ERI4, Cache4 and PrimitiveStore (→ Cache2), because its members and signatures name them.
- **Graph simulation with that one edge cut:** BasisCache 382 → 169, ERI4 400 → 197, Cache4 411 → 200, PrimitiveStore 407 → 207, Cache2 415 → 215.

**Implementers of `Cache()`: ONE.**
- `src/BasisSet/Internal/IrrepBasisSetImp.C:38`: `IrrepBasisSetImp::Cache() override {return *itsCache;}`.
- Every family (PG_Cart, PG_Spherical, PG_LibCint, Radial, GPW, PW, ...) inherits it from there.

**Callers** (32 call sites in 13 files; the three `template <class U> ... Cache()` hits in Hubbard.C, GSpaceDensityProjector.C and DeltaFunctionFitter.C are unrelated members of the same name).

| where | through what | calls |
|---|---|---|
| `src/BasisSet/Imp/IrrepBasisSet.C:15` | **THE FACE**: `Integrals_Overlap<T>::Overlap()`, a mixin of the `IrrepBasisSet` face | `this->Cache().Get<T>(Own::Overlap, …)` |
| `src/BasisSet/Imp/Orbital_1E_IBS.C:16` | **THE FACE**: `Integrals_Kinetic<T>::Kinetic()`, a mixin of the tier-1 hub `Orbital_1E_IBS` | `Get<T>(Own::Kinetic, …)` |
| `src/BasisSet/Imp/Orbital_DFT_IBS.C:16,22` | **THE FACE**: the hub `Orbital_DFT_IBS` 3-centre accessors (⟨ab\|c⟩ Overlap / Repulsion) | `Get<TFit>(Three::…, …)` |
| `src/BasisSet/Imp/Orbital_DHF_IBS.C:13` | **THE FACE**: `qchem.BasisSet.Orbital_DHF_IBS` RestMass | `Get<T>(Own::RestMass, …)` |
| `src/BasisSet/Internal/Imp/Orbital_DHF_IBS.C:17` | Internal `Orbital_RKBL_IBS<T>::Kinetic(rkbs)` | `GetRKBCross<T>` |
| `src/BasisSet/Internal/Imp/Orbital_ERI4_IBS.C:38` | **ON A PARTNER**: `RequireSameParent(c)` compares `&c.Cache()` with `&this->Cache()`, with c an `Orbital_ERI4_IBS<T>` (an Internal type) | identity check |
| `src/BasisSet/Internal/Imp/Orbital_ERI4_IBS.C:50,56` | Internal `Orbital_ERI4_IBS` Direct / Exchange (the partner gives only `CacheIndex()`) | `Get(Pair::…)` |
| `src/BasisSet/Radial/IrrepBasisSet.C:417,462` | the radial family (this) | `GetCache4` |
| `src/BasisSet/Radial/BasisSet.C:107…326` (8) | the radial family's block ctors (this) | `Register4` |
| `src/BasisSet/Gaussian/{PG_Cart,Point/PG_Spherical,Point/PG_LibCint}/Imp/IrrepBasisSet.C` (8) | the Gaussian families' block ctors (this) | `Primitives()` |

**What the inventory shows (it changes the candidate fix in the item 3 verdict):**
- The CACHED ACCESSORS themselves live on the FACES. `Overlap()`, `Kinetic()`, the ⟨ab|c⟩ pair and `RestMass()` are template mixins of the faces (`Integrals_Overlap<T>`, `Integrals_Kinetic<T>`, `Orbital_DFT_IBS`, `Orbital_DHF_IBS`), with their bodies in the faces' own implementation units. Those bodies read `Cache()` virtually, through the face.
- So taking `Cache()` off `IrrepBasisSet_IDs` is not just "declare it on the Imp base". Those 5 face-level accessor bodies would lose their way to the cache. Three ways out, each a design question:
  - **(a) Move the accessor bodies** to the Internal implementation side (IrrepBasisSetImp or a mixin there), so the face keeps only `virtual Overlap()=0` and its kin. That changes the hubs' declarations (non-virtual → pure virtual), which is a hub change.
  - **(b) A cross-cast** from the face to a small Internal abstract `CachedBlock {Cache(); CacheIndex();}` inside each accessor. Abstract to abstract is legal under the cast rule, but it costs a dynamic_cast per cached read.
  - **(c) Leave it.** Accept 382 as the cost of the face carrying its cached accessors.
- **The partner call** (`RequireSameParent`) is on the Internal `Orbital_ERI4_IBS<T>` type, not on the face. Any of (a)-(c) keeps it working if the Internal type keeps `Cache()`.
- `CacheIndex()` (a size_t) needs no import and can stay on the face in every option.

## Item 4: function bodies to `Imp/`, by rank (lane B 2026-11-18)
**Scope (controller):** Structure, PositionJet, Mesh, PairProjection, Iterators and SlotMap. OccupationPolicy, ChargeDensity, Hamiltonian and Calculation/SCFIterator/setups wait for lane S's P1-P5.

**What was there to move:**
- **Outcome, Iterators, SlotMap, PositionJet and PairProjection are class templates.** Their "non-template bodies" in the ranking were template members, which must stay in the interface. Nothing to move.
- **Structure** (lane-b ceebb8f28):
  - Moved to `Imp/{Structure,Atom,Molecule}.C`: the SiteRange / iterator bodies, `~Structure`, `G0Background`, `Admits`, `Sites`, and the Atom / Molecule overrides (`Clone`, `GetNumAtoms`, `Admits`, `GetSite` ×2, `Molecule()`).
- **Mesh** (ceebb8f28, new `src/Mesh/Imp/Mesh.C`):
  - Moved: the constructors, `Points`, the Fourier-sampler install, `FourierTarget`, `SiteIntegrals`, `Integrate`, the imposed-symmetry names, and the non-exported helpers.
  - **INLINE BY DESIGN**, with the reason written at each:
    - `Weights`, `size`, `NSites`, `SiteBegin`, `SiteEnd`: they are read inside point loops, and this build has no LTO, so an out-of-line body would cost a call per point.
    - `BeckeCutoff`: an Enzyme derivative leaf.
  - `ExactAngularDegree` became `constexpr`.

**Measured:**
- Build: the Mesh+Structure move compiled 486 TUs, which is Mesh's blast of 485 plus the new TU.
- Results: no number moved; `anchordiff` showed 0 movers on the short tier (1703 compared) and the long tier (18).
- Checks: HubInterfaceAudit stayed green, and RelChecked `allTests` builds.

**The payoff is modest, and why:**
- Over the two months before the move, `Structure.C` had 9 commits and `Mesh.C` had 16.
- Most of those changed DECLARATIONS, which still rebuild every importer. Only 1 + 2 were body-only edits, which now rebuild a single TU.
- Estimate: about 1,400 TU compiles saved over those two months.
- The convention earns its keep where bodies churn while the declarations hold still. It does not help an interface that is still growing.
