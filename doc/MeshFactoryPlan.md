# V-MESH-FACTORY — the plan (lane H, 2026-11-21; UNDER EXECUTION)

**APPROVED as a whole: USER RULING 2026-11-21, "Please proceed with the plan"** (§8's steps, §4.4 (ii), §5's hub diffs as written;
a hub diff that departs from §5 goes back to the user first).  Born as lane H's design note in the session notes; moved here when its
execution began, per the doc/ rule that a plan lives at the top of doc/ only while it is executed.  Progress: §8 steps 1 and 2 landed; step 3 READY (landing); step 4 built (4a, 4b), sweeps in flight.

This note is design only: no code. It surveys main at 8952f826c (lane-h-group3b, which contains main 41581b896; landed as main 360fee3ff, so the line numbers hold there).
The row is the controller's, carrying the user's framing ("big compile time when we extend (OCP), as a hint there is a design flaw hiding
somewhere"; "in your factory pattern it would run the other way").

## 0. Verdict in five lines

1. **The flaw is real.** A RECIPE type (`qcMesh::MeshParams`, plus `RadialKind`, `AngularKind`, `CellKind`, and the Becke knobs `beckeOrder`
   and `beckeEps`) sits in the signatures of three tier-1 hubs. It also shares one module interface (`qchem.Mesh`) with the BUILT object
   (`qcMesh::Mesh`), which every hub needs anyway.
2. **The row's evidence needs one correction.** Group 3a's 757-TU rebuild came from `Symmetry/Group.C` and the Fit_IBS hub, not from
   `MeshParams`: 3a's diff does not touch `src/Mesh/Mesh.C`. §6 counts the import graph: `qchem.Mesh` alone reaches
   507 of 830 TUs, so the case stands without 3a.
3. **The target works and keeps what "the structure owns the strategy" was protecting** (§2). It needs three things together:
   - (a) a module split, the recipe moving out of `qchem.Mesh`;
   - (b) a mesh factory, one overload per concrete structure kind;
   - (c) library code receiving a built `Mesh` (or a `MeshMaker` where the library must choose the recipe itself).

   (a) without (c) buys nothing: a hub that still NAMES the recipe still imports it, and a module's importers rebuild when anything they
   import transitively changes.
4. **Hub diffs:**
   - `Structure` loses `CreateIntegrationMesh` (×3).
   - `BasisSet` and `Orbital_DFT_IBS` swap `MeshParams` in `CreateVxcFitBasisSet` for what the fits actually read.
   - **`MakeXCQuadrature` leaves the `qchem.BasisSet` hub.** This is a fourth hub signature the row did not list: it is in the snapshot.

   §5 has every hub diff exactly.
5. **The one decision is RULED (user, 2026-11-21): §4.4 takes (ii).** The concrete GPW basis receives a `MeshMaker` through its
   constructor and builds its PP-fallback mesh itself. No mesh argument goes on any `Orbital_1E_IBS` / `Orbital_PP_IBS` entry point: *"I
   think interfaces like Orbital_1E_IBS should not have mesh related entities in them."* §7 assesses the hub faces that would still
   carry mesh entities after the plan.

## 1. What names the recipe today

Of 45 module interfaces that name `MeshParams` (none under `Imp/` or `tests/`), by layer:

| layer | interfaces | after the target |
|---|---|---|
| **tier-1 hubs** | `Structure.C` (Structure, Atom, Molecule), `BasisSet.C` (`CreateVxcFitBasisSet`, **`MakeXCQuadrature`**), `Orbital_DFT_IBS.C` | **none** (§5) |
| **hub re-exports** | `Fit_IBS.C` and `Internal/Fit_IBS.C` (`export import qchem.Mesh` "for MeshParams"), and `Hamiltonian/Types.C` (`export import qchem.Mesh; // qcMesh::MeshParams`) | they re-export the built-mesh module only |
| **mesh library** | `Mesh.C`, `XCPolicy.C`, `RadialMesh.C`, `AngularMesh.C`, `BeckeClimb.C`, `Internal/RadialFactory.C`, `Internal/AngularFactory.C`, `SiteMultipoles.C` | unchanged: these ARE the factory's internals; they import `qchem.Mesh.Recipe` |
| **policy layer** | `Calculation/*` (RunOptions, Setup, Solid/Molecule/AtomSetup, the three facades, Deck, SetupDecisions) | unchanged: recipes are policy |
| **fit factories** | `PG_Cart/BasisSet.C`, `PG_Spherical/BasisSet.C`, `Radial/BasisSet.C`, `Radial/IrrepBasisSet.C` (EFit ctor), `PlaneWave_IBS.C`, `GPW_IBS.C`, `SymmetryAdapted_IBS.C` | take what they read (§4.2) |
| **library self-builders** | `Hamiltonian/Factory.C`, `Internal/Hamiltonians.C`, `Seed.C`, `SCFIterator.C`, `GPW_Evaluator.C`, `DensitySampler.C` (2), `UnitCell.C`, `Lattice_3D/Imposition.C` | §4 |

The library calls that BUILD a mesh from a recipe through the abstract `Structure&`, outside tests and facades:
- `Fit_IBS::FitMesh` (the A1 fit's private projection mesh);
- `MakeXCQuadrature`;
- `Hamiltonians.C:379` (`ppQuadrature`);
- `Seed.C:298` (the SAD seed's own mesh);
- `GPW_Evaluator.C:1529` (the PP nonlocal mesh at `PPMeshParams(cl,sep)`);
- `Imposition.C` (×3, which holds the concrete `UnitCell`).

## 2. Against "the structure owns the mesh; it owns the strategy"

The ruling is at `Structure/Imp/UnitCell.C:488`. It arrived with R2.17 (599f8eee1, `CreateSiteAdaptedBeckeMesh` → the ops-taking
overload), in the move that took the Uniform-vs-Becke switch OUT of GPW_IBS's XC-quadrature factory. Read for what it FORBADE:

> "a BASIS was switching on cellKind to decide how a STRUCTURE builds its own mesh."

That is a ruling about **who may switch on a recipe kind**: never a library that is not the structure's. The factory keeps it:
- the switch lives in `MakeMesh(const UnitCell&, …)`, keyed by the concrete structure type;
- no basis, term or density sees a kind.

**What changes is the second half, "it owns the strategy".** Today ownership is a VIRTUAL on the abstract face. That is exactly what
lets any library holding a `const Structure*` build a mesh from a recipe, so every holder of a Structure imports the recipe. The
`Structure.C:128` comment, "Pure virtual -- each geometry owns its most efficient mesh (no central dispatch)", survives as follows:
- **no central dispatch:** overload resolution on the concrete type. There is no runtime switch and no cast.
- **each geometry owns its most efficient mesh:** the builders are today's bodies, moved verbatim (`MakeMolecularMesh`,
  `MakePeriodicBeckeMesh`, the uniform raster, the site-adapted build). They stay in qcStructure, beside the geometry they read.

**What is honestly lost:**
- **(i) Building a mesh from a `Structure&` alone.** The callers that did so are exactly the library self-builders of §4, and this note
  argues they should not.
- **(ii) A second hierarchy keyed by the same kinds.** The setups need a `MeshMaker` per kind. One class template,
  `StructureMeshMaker<S>` (§3), closes over `MakeMesh(const S&, …)`, so the parallel hierarchy is ONE template and not hand-written.
- **(iii) A new STRUCTURE kind** (a slab) adds an overload instead of overriding a virtual. That is OCP-compatible (adding, not
  modifying), in a non-hub module.

It also matches the house DIP rule (CLAUDE.md, "Constructors pass concrete objects; factories concentrate the concrete associations"): a
setup that knows its kind should DEMAND the concrete structure. Today `Setup::structure` is a `shared_ptr<const Structure>` (Setup.C:50).

**My recommendation: go to the factory.** The previous move was right about the switch, and that part is kept.

## 3. The target, concretely

- **`qchem.Mesh`** (qcMesh) keeps the BUILT object only:
  - `Mesh`, `MeshBuilder`, `SiteIntegrals`, `Integrate`, `MeshWeightResponse`;
  - the `Forces()` stamp (a `GroupPtr`, so it keeps `import qchem.Symmetry.Group`).

  Every hub imports this, and it changes rarely.
- **`qchem.Mesh.Recipe`** (new, qcMesh) holds:
  - `MeshParams`, `RadialKind`, `AngularKind`, `UnitCellKind` / `CellKind`, kept as ENUMS, with every switch exhaustive and no
    `default:` (-Wswitch is the compile-time check, as the user prefers for a segregated factory interface);
  - the pin-30 descriptor `qcMesh::Forces(const MeshParams&, const Frame&)`, moved from `Mesh.C` (3b put it there; it reads a recipe);
  - ~~the ABSTRACT seam `MeshMaker`~~ **AS BUILT (step 2): `MeshMaker` has its own module, `qchem.Mesh.Maker`** (above both
    `qchem.Mesh` and `qchem.Mesh.Recipe`): `Make` returns a `qcMesh::Mesh`, and `qchem.Mesh` imports the recipe, so the seam cannot sit
    in the recipe module without a cycle.  It still names only a recipe and a built mesh: a library that holds one (the GPW basis, §4.4)
    imports `qchem.Mesh.Maker` but not the structure factory.  Dependency inversion: the library depends on the abstraction, and the
    concrete `StructureMeshMaker<S>` stays at the edge.

  `qchem.Mesh.XCPolicy` and `BeckeClimb` already sit above it.
- **`qchem.Structure.MeshFactory`** (new, qcStructure, the policy-facing wiring) holds:
  ```
  qcMesh::Mesh MakeMesh(const Atom&,     const qcMesh::MeshParams&);
  qcMesh::Mesh MakeMesh(const Molecule&, const qcMesh::MeshParams&);
  qcMesh::Mesh MakeMesh(const UnitCell&, const qcMesh::MeshParams&);
  qcMesh::Mesh MakeMesh(const UnitCell&, const qcMesh::MeshParams&, const std::vector<Symmetry::SymOp>& ops);   // invariant (§3a imposed)
  //! the concrete realization of qcMesh::MeshMaker (declared in qchem.Mesh.Recipe), handed to a library that must choose a recipe itself (§4)
  template <class S> class StructureMeshMaker : public qcMesh::MeshMaker { /* holds shared_ptr<const S>; Make = MakeMesh(*s, mp) */ };
  ```
  - The bodies of `CellMeshForces` and `AsUniform` (3b) move with the UnitCell overloads, so the pin-30 refusal before any point is built
    stays at the build.
  - ~~Each setup builds its `StructureMeshMaker<S>` and holds it next to the structure.~~ **AS BUILT (step 2; controller ruling
    2026-11-22, option (A)): A MAKER BELONGS TO A GEOMETRY POINT, NOT TO THE SETUP.**  A maker over the setup's REFERENCE structure would
    build at the reference geometry, and a geometry step moves the atoms.  Where `S` is known:
    - **Molecule.** `MoleculeSetup::Resolve` / `Admit` took `const Structure&`, and every point was a polymorphic `Clone()` + `MoveSites`,
      so `S` was lost at each point.  The doors are now TYPED (`const Molecule&`), and `Resolve` stores a typed POINT BUILDER that copies the
      reference as an `S`, moves its sites, adds the PP core charges, and pairs it with `StructureMeshMaker<S>` over that point.
      `MoleculeSetup::MakePoint()` / `MakePoint(R)` hand out a `Point {structure, meshMaker}`; `MoleculeCalculation` is constructed from
      one and holds `itsMeshMaker` beside `itsStructure`.  No cast, no hub change.  A caller holding only a `Structure&` spells out its
      concrete type (no forwarder); a door on `const UnitCell&` is added only if a real (non-test) path sends a cell here.
    - **Solid.** Every solid point is already a concrete `UnitCell` (`Lattice_3D::GetUnitCell`, `FollowCell`'s `FollowedCell::cell`), so the
      step-3 consumer builds `StructureMeshMaker<UnitCell>` from the point's own cell; the setup holds none.
    - **Atom.** `AtomSetup` builds its `Atom` itself: the same, at its step-3 consumers.
  - A one-atom MOLECULE (lane G's Q1 route) is a `Molecule`, so it gets the Becke build, as today.

## 4. The library self-builders, one by one

1. **`MakeXCQuadrature`** (BasisSet hub, `BasisSet/Imp/BasisSet.C:98`).
   - Today it cross-casts `SymmetryImposition` → `XCMeshSource` and otherwise calls `cl.CreateIntegrationMesh(mp)`.
   - It moves to the factory layer: the Hamiltonian factory's caller (the setup) builds the `FitQuadrature` and INJECTS it.
   - **AS BUILT (step 3; controller + user rulings 2026-11-22): the setup injects the MAKER and a SYMMETRY SOURCE, not a built quadrature.**
     The Hamiltonian builds its quadratures lazily and by route: `XCMeshShare` on its first ask, none on the raster route ("pays for no
     mesh"), a second at the molecular Hartree grid's own `nRadial`.  Injecting BUILT quadratures would build meshes nobody asked for, or
     copy those route conditions into every setup.  A quadrature depends only on geometry and recipe, never on the density, so building
     it lazily is allowed.
     - **One interface, not two class networks** (the user: "why do we need two parallel class networks for raw and folded meshes? And
       where does the Symmetry class enter for folding?").  `qcMesh::MeshMaker` has two overloads: `Mesh Make(mp)` and
       `FoldedMesh Make(mp, const SymmetrySource&)`.  The second is a non-virtual TEMPLATE METHOD: no ops → the plain mesh; ops but no
       density averaging (a free run's symmetric numerics) → the invariant mesh, unfolded; averaging → invariant and orbit-folded.  Only
       the two steps that need GEOMETRY are the realization's protected hooks, `MakeInvariant` and `FoldOrbits` (the fold needs the cell
       matrix and the Shubnikov flip-fixed test), realized in `StructureMeshMaker<S>` through `MakeMesh(S, mp, ops)` / `FoldOrbits(S, mesh,
       ops)`: the cell's are the moved body of `CrystalImposition::MakeXCMesh`; an atom's invariant mesh is its radial-only set (an angular
       recipe throws: a finite direction set claims a symmetry it lacks) with the trivial fold; a molecule's throws by name (OOD-SOLID
       V-FOLD-NEUTRAL: a point-group mesh is not built yet; until then the molecule setup's policy hands it no ops, and the empty op set
       returns before either hook).
     - **`qcMesh::SymmetrySource`** (data-free: `MeshOps()`, `AveragesDensityOverOrbits()`, named for what the maker does with the answer)
       supplies OPS ONLY; the maker owns the geometry.  So a source resolved at one cell can never move a mesh built over another (the
       user's "does the grid builder hold state information?"; gated by `MeshFactory.AStrainedCellsMakerBuildsOverItsOwnCellUnderTheReferencesOps`).
       Realized by `CrystalImposition` (listed after its other public faces, the cross-cast trap) and `NoSymmetry` (the plane-wave family,
       the HF bootstrap's sibling, an atom, a C1 molecule).  **A molecule HAS a point group** (the user, 2026-11-22: "molecules have a point
       group"): handing it `NoSymmetry` silently would be V-FOLD-NEUTRAL's defect under a new name.  So the library stays truthful (the
       molecular hooks refuse non-trivial ops by name) and the POLICY layer chooses: `MoleculeSetup` detects the group and, until the fold
       exists, picks the unfolded mesh (`MoleculeSetup::meshSymmetry`) and files policy row `MeshFoldNotBuilt` ("point group C2v (4 ops)
       detected; the XC mesh is not folded", citing V-FOLD-NEUTRAL) on every run that builds an XC mesh.  A molecular point-group source
       realizes the face later (V-FOLD-NEUTRAL, lane H's next task) with no library change.  *(Done 2026-11-23: V-FOLD-NEUTRAL step 2 `9fb57e3e7`.  `MolecularPointGroupSource` realizes it, with no library face change, and policy row `MeshSymmetry` replaces `MeshFoldNotBuilt`.)*
     - **No name carries the client's intended use** (the user: "MakeXCMesh should not have XC in the name"): the Hamiltonian's helper is
       `BuildMeshQuadrature`.  `CrystalImposition` no longer builds meshes: its `MakeXCMesh` is a one-line transitional forwarder for
       `XCMeshSource`'s last client, `MakeXCQuadrature` (a hub declaration with ~20 test callers); step 4 deletes the four together
       (`MakeXCQuadrature`, `XCMeshSource`, `MakeXCMesh`, `CrystalImposition::cell`) and moves those tests to `maker.Make(mp, source)`.
     - `Fold` / `FoldPoints` stay in `Symmetry::Lattice_3D` though they are geometry-free bookkeeping (noted for V-SYM-ONE-SYMOP).
   - The imposed path asks `CrystalImposition::MakeXCMesh`, which holds the concrete cell and calls
     `MakeMesh(cell, mp, ops)` / `MakeMesh(cell, mp)`.
   - The cross-cast leaves the library, as D-POINTSET-CAST prefers.
   - **V-FOLD-NEUTRAL synergy:** that row moves `XCMeshSource` out of `Lattice_3D` under the nonlinear-mesh name. The factory layer is
     where a point-group imposition would realize it, and then `MakeMesh(const Molecule&, mp, ops)` is one added overload with no hub
     edit. (Today the UnitCell.C:145 note says it would be "hoisted onto Structure", a hub edit.)
2. **The fit factories (`CreateVxcFitBasisSet`).** What each realization actually READS from `MeshParams`:
   - **PW and GPW:** `relCutoff` ONLY. It is the fit-grid multiplier (CP2K REL_CUTOFF), set from the functional's `GridCutoffFactor()`. It
     is not a mesh at all. Its home in `MeshParams` is itself a defect of the recipe.
   - **A1 Gaussian (PG_Cart, PG_Spherical) and the radial EFit:** a projection MESH, built privately via `Fit_IBS(st, mp)` →
     `st.CreateIntegrationMesh(mp)`.
     - On the molecular Auxiliary path, `Hamiltonians.C:199` builds the XC quadrature at the SAME `mp`, then the A1 fit builds a second
       mesh from that recipe.
     - **Handing the fit `quad->GetMesh()` removes a duplicate build.** The points and order are identical, so it is bit-identical. Verify
       with anchordiff: the console loses the `auxiliaryFit` mesh scope, so a molecular console byte test may need its expected text.
   - **δ, multipole:** the quadrature, already injected.
   - So the blocks take `(const FitQuadrature*, double fitCutoffFactor)`; §5 gives the exact form.
   - An A1 block handed no quadrature THROWS by name (it is a quadrature fit), the same shape as δ's existing "none was handed in".
3. **`ppQuadrature` (`Hamiltonians.C:379`) and the SAD seed (`Seed.C:298`).**
   - Both build at the run's XC recipe through `st->CreateIntegrationMesh(mp)`, and the PP term is V-NONLINEAR-MESH's `nextClient`.
   - The Hamiltonian factory receives a `MeshMaker` (or the built meshes) beside the recipe: `Hamiltonian/Factory.C` is a factory, so a
     recipe in its signature is allowed. `MakeSeedDensity` receives the built `shared_ptr<const Mesh>`.
   - **`SCFIterator`'s HF SAD bootstrap** builds an LDA sibling Hamiltonian from `optional<MeshParams> seedMesh` INSIDE the library. It
     takes the `MeshMaker` too.
   - **Flagged, not fixed here:** a library building a sibling Hamiltonian is a policy act that may belong to the facade.
4. **`GPW_Evaluator.C:1529`: RULED (ii), user 2026-11-21.** The GPW basis computes `PPMeshParams(cl, sep)` from its own projector reach
   (the V2.5 floor: the density cutoff lifted to C·α_max + α_β) and builds through `cl->CreateIntegrationMesh`. Here the recipe is a fact
   only the library knows. The mesh is the FALLBACK of `MakeSeparablePP` / `MakeSeparablePPByL` for a projector model with no
   closed-Gaussian face: a `SpeciesProjectorSet_Gaussian` (GTH, HGH) takes the analytic branch above it (`GPW_Evaluator.C:1495-1525`)
   and builds no mesh.
   - **The shape.** `tGPW_IBS<T>`'s constructors (`GPW_IBS.C:27,37,47`) DEMAND a `std::shared_ptr<const qcMesh::MeshMaker>`: no default
     and no null, per the concrete-constructor rule.
     - It is stored in `GPW_Evaluator` beside `itsCell`.
     - Line 1529 becomes `itsMeshMaker->Make(PPMeshParams(cl, sep))` inside the same `report::Scope("ppNonLocal")`, so the run report is
       unchanged.
     - The `cl` argument still supplies the PP centres and Z, as today.
   - **Who supplies it.** The whole-set factory (`BasisSet/Lattice/Imp/BasisSet.C:201`, the per-k-block `new tGPW_IBS<U>`) passes the
     maker it is given through `GPWParams` (non-hub), one maker shared by every k-block like `mol` and `ladder`. The solid setup builds
     it as `StructureMeshMaker<UnitCell>` (§3), the same maker its own meshes come from.
     - The unit tests that build a `tGPW_IBS` by hand (`GPW_UT.C`, `VxcProductSpace.C`) build one from the cell they already hold.
   - **AS BUILT (step 3):** the maker is a DEMANDED positional argument of `GPWFactory(lat, mol, meshMaker, p)` (and of `GPW_BasisSet`),
     not a `GPWParams` field: a struct field left out of a designated initializer compiles to a silent null, an argument cannot (CLAUDE.md,
     compile-time over runtime).  The two library passes that build a GPW basis for its analytic overlaps alone, `VetStageTrim` and
     `ShellConditioning`, take it too.  `MakeMeshMaker(const S&)` (`qchem.Structure.MeshFactory`) makes one over a copy of a held
     geometry, which every hand-built block in the ~100 test and facade sites uses.
   - **What no face gains.** `Orbital_PP_IBS::MakeProjectorMatrix` / `MakeProjectorMatrixByL` / `MakeProjectorPullback` keep their
     signatures. No `Orbital_1E_IBS`, `Orbital_PP_IBS` or `BasisSet` entry point names a mesh, a recipe or a maker. Today's
     `Orbital_1E_IBS` snapshot (lines 113-138) carries no mesh or quadrature entity, and the plan keeps it that way.
   - **Why not (i), the earlier recommendation.** It put a mesh argument on the PP entry point, which is exactly the hub-face mesh entity
     the user rules out. My case against (ii) was that "the basis then holds the ability to make meshes for its whole life". The user
     rejected it ("It's just a pointer no?"), and I withdraw it. The capability is one immutable pointer to a geometry the basis already
     holds by value (`itsCell`), so holding it costs nothing, and it is what the concrete-constructor rule asks for.
   - **What (ii) gives up from (i):** the PP mesh is not built in the factory's sequence. The `ppNonLocal` scope still nests where the
     PP term's build calls in, as today, so the run report still shows when it is built.
5. **`DensitySampler.C` (×2)** names `MeshParams` in comments and in the sampler's ctor documentation only. It integrates on the mesh it
   is handed: mechanical.
6. **`Imposition.C`** already holds the concrete `UnitCell`: it calls the factory overloads directly.

The **ε climbs** (pin 31: `ClimbBeckeLadder`, the mesh pilot, `xcMeshAudit`) live in the facades and move RECIPES. They rebuild a mesh
through the setup's maker, so the factory adds no new path.

## 5. Exact hub diffs (each needs a user ruling cited in the snapshot commit)

**`qchem.Structure`** (src/Structure/Structure.C):
```
- import qchem.Mesh;        // qcMesh::Mesh / MeshParams -- the real-space integration mesh a geometry produces
- virtual qcMesh::Mesh CreateIntegrationMesh(const qcMesh::MeshParams&) const=0;            // Structure
- virtual qcMesh::Mesh CreateIntegrationMesh(const qcMesh::MeshParams&) const override;     // Atom
- virtual qcMesh::Mesh CreateIntegrationMesh(const qcMesh::MeshParams&) const override;     // Molecule
- qcMesh::Mesh MakeMolecularMesh(const Structure&, const qcMesh::MeshParams&);              // module-internal; moves to MeshFactory's Imp
```
- **Not a hub:** `UnitCell.C`'s two `CreateIntegrationMesh` overloads leave the same way.
- After this, Structure names no mesh at all, and the `qchem.Structure` interface stops importing `qchem.Mesh`.

**`qchem.BasisSet`** (src/BasisSet/BasisSet.C):
```
- virtual FIT_SF_ABS<T>* CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams&, VxcFit fit=VxcFit::Auto, const FitQuadrature* quad=nullptr, const PhiChooser& phiTables={}) const;
+ virtual FIT_SF_ABS<T>* CreateVxcFitBasisSet(const Structure* cl, VxcFit fit, const FitQuadrature* quad, double fitCutoffFactor=1.0, const PhiChooser& phiTables={}) const;
- FitQuadrature MakeXCQuadrature(const Structure& cl, const qcMesh::MeshParams& mp, const Symmetry::SymmetryImposition* imp);
```
- `fitCutoffFactor` is today's `MeshParams::relCutoff`, under the name of what it is.
- `quad` stays a pointer: a plane-wave run has no real-space quadrature.

**`qchem.BasisSet.Orbital_DFT_IBS`** (src/BasisSet/Orbital_DFT_IBS.C):
```
- virtual FIT_SF_ABS<TFit>* CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams&) const=0;
+ virtual FIT_SF_ABS<TFit>* CreateVxcFitBasisSet(const Structure*, const FitQuadrature* quad, double fitCutoffFactor) const=0;
```
- **Alternative to weigh:** two faces, a quadrature fit and a G-grid fit. Each realization reads only one of the two arguments. I keep
  one face because the whole-set factory (`tBasisSet::CreateVxcFitBasisSet`) already chooses between realizations, and splitting the face
  would move that choice into every caller.

**`Orbital_DFT_IBS:Fit_IBS` partition:** `- export import qchem.Mesh;`. The face names no mesh in any signature (snapshot lines 63-112):
the re-export existed "for MeshParams" only.
- A re-export trim needs no design discussion (user, 2026-11-05). Importers that leaned on it import `qchem.Mesh` themselves.
- The `Internal/Fit_IBS` ctor `Fit_IBS(const Structure&, const MeshParams&)` becomes `Fit_IBS(std::shared_ptr<const qcMesh::Mesh>)`.
  It is module-internal and not in the snapshot.

**No change:** Orbital_1E_IBS (§4.4: the GPW PP mesh enters through the CONCRETE `tGPW_IBS` constructor, which is not a hub, and
`Orbital_PP_IBS`'s projector faces keep their signatures), Hamiltonian, ChargeDensity, FunctionFitter, Representation, WaveFunction,
LASolver and G_FieldEvaluator.

## 6. Evidence, measured from the module import graph

**Method.** For each module, count every TU that imports it TRANSITIVELY:
- an importer of an importer counts;
- an implementation unit counts as importing its own interface;
- partitions are resolved to their primary module.

The script is at `/tmp/claude-1000/-home-janr-Code-qchem6/1a9dcc67-8090-4d84-9156-cf625cf34b2f/scratchpad/blast.py`. It covers 830 TUs in
src/, IntegrationTests/ and CLIapps/.

**Why not a dry run.** `ninja -n` is NOT a measurement here. It cannot read the dyndep files, so touching the 4-importer leaf
`qchem.Mesh.BeckeClimb` lists the same 678 compiles as touching `qchem.Mesh` itself.

| module touched | TUs rebuilt (transitive importers) | why |
|---|---|---|
| `qchem.Mesh` (MeshParams + enums + the built Mesh, today) | **507** | Structure, BasisSet's Fit_IBS (`export import`) and Hamiltonian/Types (`export import`) all pull it in |
| `qchem.Symmetry.Group` | 510 | Structure, Mesh and every fit import it.  It drove 3a's rebuild |
| `qchem.Structure` | 469 | |
| `qchem.BasisSet.Orbital_DFT_IBS:Fit_IBS` | 372 | |
| `qchem.BasisSet.Fit_Types` | 290 | |
| `qchem.BasisSet` | 274 | |
| `qchem.Hamiltonian.Types` | 173 | |
| `qchem.Mesh.XCPolicy` | **96** | already a recipe-side module above qchem.Mesh: the proxy for `qchem.Mesh.Recipe`'s reach after the target |
| `qchem.Mesh.BeckeClimb` | 4 | |

- **3a and 3b.** The union of 3a's changed files reaches 511 TUs by the graph, and 3b's union also reaches 511. The controller's observed
  757 of 874 compiles exceeds the graph closure. So the build system recompiles more than the import graph strictly requires (CMake's
  collated per-target module maps are my guess, not verified). That overhead is a separate lever (ModuleToolchainPlan), and it scales
  every number above.
- **The acceptance number.** Today, adding a `RadialKind` enumerator (DE2) or a partition kind edits `qchem.Mesh`: 507 TUs or more.
  After the target, the same edit touches `qchem.Mesh.Recipe`, which should land near XCPolicy's ~100 (the mesh library, the factory,
  the policy layer, the Hamiltonian factory, the PP/GPW sites). §8 step 5 re-measures it with the same script.

## 7. The cross-cutting contracts

- **The symmetry-imposed path** is §4.1 and §4.6.
  - `MakeMesh(cell, mp, ops)` keeps W2b's site-adapted invariant build and the Uniform `MakeInvariant` branch. These are today's
    `UnitCell::CreateIntegrationMesh(mp, ops)` body, verbatim.
  - The imposition keeps owning the ops ("no op set or cell matrix leaves the imposition", Imposition.C:49).
  - The D5 warning (a stated angular rule REPLACED under an imposed group) moves with the body.
  - **Recommended alongside, not required:** it is a `static bool warned` printed with `std::cout`, which sits outside the policy scope
    of the run report. It is a policy in a library, so it should report under `policy`.
- **The per-site block contract** (`MeshBuilder::BeginSite`, `RequireSiteBlocks`, `qcMesh::SiteIntegrals`).
  - `SiteIntegrals` is a query on the BUILT mesh. Its consumers (`DensitySampler`'s site moments, the δ fit, `SiteMultipoles`) see only
    `Mesh` and are untouched.
  - The builders keep calling `BeginSite` per centre.
  - **One gain:** today `RequireSiteBlocks` runs on the imposed Becke path only (`UnitCell.C:583`). The factory is one place every
    atom-centred mesh passes through, so it can require site blocks on EVERY atom-centred overload (Atom, Molecule, periodic Becke, free or
    imposed) at O(sites) cost, with a throw that survives NDEBUG.
- **The Becke partition classes** (lane S: Becke / SSF / modified Becke behind a bulk virtual face). I have not seen lane S's class design;
  this is the contract the factory offers it.
  - The recipe carries the partition KIND as an enum beside `beckeOrder` / `beckeEps`, which are already in `MeshParams`.
  - The factory's exhaustive switch builds the partition object and hands it to `MakeMolecularMesh` / `MakePeriodicBeckeMesh`.
  - The built `Mesh` holds it as its `MeshWeightResponse`. `BeckePartition.C`'s header already says "whoever owns the mesh then, the
    weight's derivative is a property of the mesh", so the forces' grid response differentiates the same object the build used.
  - **Adding SSF then means** one class (in qcMesh, beside `BeckePartition`), one enumerator in `qchem.Mesh.Recipe`, one `case` in the
    factory, and zero hub edits.
  - **Today the same addition** edits `Mesh.C` (`MeshParams`), which `Structure.C` imports, so it rebuilds §6's count. It also edits both
    builders.
  - **The bulk face** (`Weights(points, owner, out)` and its gradient) suits the factory: the builder loops points per owner, and the
    periodic build's image gathering (`PeriodicBeckeWeight`) is a per-point competitor list that the partition consumes.
- **Hub faces that still carry a mesh entity after the plan.** **RULED: (a) and (b) stay: the user, 2026-11-21, "a and b were exactly
  the right places for grid stuff."** No follow-on removal rows. The assessment as written: After §5, no hub names a
  recipe (`MeshParams`), a structure-built mesh call or `MakeXCQuadrature`. What remains:
  - **`BasisSet::CreateVxcFitBasisSet` and `Orbital_DFT_IBS::CreateVxcFitBasisSet` take a `const FitQuadrature*`** (= `qcMesh::FoldedMesh`
    plus the Shubnikov sigmas, via `qchem.BasisSet.Fit_Types`).
    - **Needed today?** Yes, as an argument somewhere. The δ fit IS one function per quadrature point, and the A1 fit projects on those
      points (§4.2), so the points are the fit's DATA, not a policy. It is a built, stable type, so it does not cause the rebuild churn this
      row is about: recipes churn, built meshes rarely do.
    - **Can it go the same way?** Yes, at a price. `CreateVxcFitBasisSet` is a factory method on the orbital basis: the family decides which
      fit realizations it can pair with (a PW basis takes a PW fit). Taking the mesh off the hubs means moving that choice to a fit factory
      at the policy edge, which receives the CONCRETE orbital basis (the whole-set factories already know it) and passes the quadrature to
      the concrete fit's constructor, as §4.4 does for the GPW PP mesh. Both `CreateVxcFitBasisSet` virtuals would then leave the hubs.
    - **Cost:** one dispatch per basis family moves out of the libraries (about 7 realizations, §1's fit-factory row), and the fit factory
      is the dirty part by design.
    - **My read (superseded by the ruling above):** worth it as a follow-on row after this plan, not inside it. This plan already shrinks the argument to what the fits read.
  - **`FunctionFitter` re-exports `qchem.Mesh.Integrator`.** That module holds `MatrixForward` / `MatrixAdjoint`, the D ↔ operator maps: no
    mesh, recipe or point set appears in its types. It is a mesh entity by its LIBRARY name only (it sits in qcMesh). Needed: yes (every
    density fit is a forward/adjoint pair). If the user wants a hub to show no `qchem.Mesh.*` import at all, the fix is a module rename
    or move, not a design change.
  - **`Fit_IBS`'s `export import qchem.Mesh`** goes in §5.
- **DE2 is not sequenced behind this** (user: "we soldier on with a big build").
  - It lands as a `RadialKind` value in today's `Mesh.C`.
  - The split later moves it with the rest of the recipe, with no extra work.

## 8. Sequencing (each step green and anchordiff 0 at 1e-10; every step is a pure relocation, so the same points arrive in the same order)

1. **Split `qchem.Mesh` → `qchem.Mesh` + `qchem.Mesh.Recipe`.**
   - No hub signature changes. The hubs import both modules for now, so this step alone buys nothing; it separates the concerns.
   - The `Forces(mp, frame)` descriptor and the abstract `MeshMaker` move.
   - **First, a gate that a GTH run never builds the `ppNonLocal` mesh** (controller, 2026-11-21). The analytic branch must not silently
     fall back to the mesh, which is what §4.4's rewiring of line 1529 could break unseen. The test is
     `CalculationFacade.AGTHPseudopotentialOnGPWBuildsNoPPNonLocalMesh` in `src/Calculation/tests/CalculationFacade.C`, beside its
     molecular twin `APseudopotentialOnAGaussianBasisBuildsNoPPQuadrature`:
     - It opens a `report::Run` and runs GPW Si Γ with gth-lda (the cheapest existing GPW Si recipe).
     - At iteration 1 it searches `report::CurrentRunReport()` RECURSIVELY and requires no `ppNonLocal` key at any depth. The search is
       recursive because the scope nests wherever the PP term's build calls in, and a path-anchored check would rot with the run sequence.
     - It also requires `hamiltonian` to be present, so a report that never got there cannot pass vacuously.
     - Its counterpart, the fallback still working, is already gated: `GPW_UT.C` `AnalyticSeparablePPMatchesMesh` (and `_DChannel`) runs both arms.
     - It lands green on today's code, before anything moves.
2. **`qchem.Structure.MeshFactory`** with the overloads and `MeshMaker` / `StructureMeshMaker<S>`.
   - The `Structure` virtuals temporarily DELEGATE to it. An Imp unit may import a module that imports its own interface: the graph is
     still acyclic.
   - The setups build and hold their maker. **As built:** the molecule setup's typed point builder hands each point its maker (§3).
   - The factory gets the every-atom-centred `RequireSiteBlocks` check.
3. **The self-builders in §4, one commit each.**
   - **AS SEQUENCED (controller ruling 2026-11-22, option (ii)): the fit factories moved to step 4.**  The A1 block builds its mesh
     behind the `Orbital_DFT_IBS` face `CreateVxcFitBasisSet(const Structure*, const MeshParams&)`; handing it the quadrature IS §5's
     hub line, so it lands with the other §5 diffs and step 3 stays hub-free.  Step 4 rebases onto lane S's δ-default landing and keeps
     its placeholder comments in the Gaussian `CreateCDFitBasisSet` / `CreateVxcFitBasisSet` bodies.
   - Order: Seed and ppQuadrature, then the fit factories (they get the duplicate-mesh saving), then `MakeXCQuadrature` into the factory
     layer, then the GPW PP mesh (§4.4 (ii): the `tGPW_IBS` constructors demand the maker, `GPWParams` carries it, and step 1's gate must
     stay green), then the SCFIterator bootstrap.
4. **The hub diffs of §5**, with their ruling. The snapshot is regenerated, and the diff is deletions (with the `Fit_IBS` re-export trim) plus the two
   `CreateVxcFitBasisSet` lines.  `Orbital_1E_IBS`'s snapshot section does not move.
   - **AS BUILT (two commits, 2026-11-22).**  4a: the fit factories take `(quad, fitCutoffFactor)`; the A1 / radial fits SHARE the
     handed-in mesh (`Fit_IBS(shared_ptr<const Mesh>)`, throws on null), so the duplicate molecular A1 mesh and its `auxiliaryFit` scope
     are gone; `MakeXCQuadrature`, `XCMeshSource`, `CrystalImposition::MakeXCMesh` and `CrystalImposition::cell` are deleted, and
     `EmitGpwGrids` takes the point's cell for the lattice-sum report instead of reading it off the imposition.  `Fit_Types` is a PLAIN
     import of `Orbital_DFT_IBS` (no new hub re-export).  4b: `CreateIntegrationMesh` leaves `Structure` / `Atom` / `Molecule` /
     `UnitCell`, and `qchem.Structure` stops importing `qchem.Mesh`.
   - **A guard the plan did not foresee:** `Molecule` has an IMPLICIT `Molecule(const Structure&)`, so `MakeMesh(*abstract, mp)` compiled
     as a molecule mesh (16 test calls on the first port, cells among them).  `MakeMesh(const Structure&, ...)` was DELETED in the
     factory as a guard, so an abstract argument was a compile error; a caller holding only the face takes a `MeshMaker`.  Step 5a
     (USER RULING 2026-11-23, b59addbf8) made the constructor `explicit` and removed the guard: an abstract argument now finds no
     viable overload at all.
5. **`Hamiltonian/Types.C` stops re-exporting `qchem.Mesh` for `MeshParams`.**
   - Re-measure §6: touching `qchem.Mesh.Recipe` should rebuild only the mesh library, the factory, the policy layer, the Hamiltonian
     factory and the PP/GPW sites.
   - That count is the acceptance number for "extending a scheme no longer modifies hubs".
   - **AS BUILT (lane H 2026-11-23):**
     - 5a (b59addbf8, USER RULING 2026-11-23): `Molecule(const Structure&)` is `explicit`; step 4b's deleted `MakeMesh(const Structure&)`
       guards are removed (an abstract argument finds no viable overload).  Three test sites had used the conversion; no library code.
     - aedca3550 (lane G's clang 23 finding): `qchem.SymmetrizeMesh`'s bodies move to `Lattice_3D/Imp/SymmetrizeMesh.C`; its exported
       inlines had exposed the TU-local helpers `Wrap01`, `Fractional`, `OpsOf` (`-WTU-local-entity-exposure`, the only ones in lane G's
       full build).
     - 5b (70fc2d4cb): `Hamiltonian/Types.C` stops re-exporting `qchem.Mesh`, and `qchem.Mesh` ends its step-1 re-export of
       `qchem.Mesh.Recipe`.  Two moves the plan did not foresee were needed: `RadialMesh` / `AngularMesh` (built value types) move into
       `qchem.Mesh`, because `qchem.Mesh.Quadrature` imported the radial FACTORY module for the type alone (446 TUs on that one import);
       and `SiteMultipoles` loses its `MeshParams` convenience ctor (it sees built meshes only; its two callers build the rules).
     - **The measured number: a `qchem.Mesh.Recipe` edit reaches 214 TUs (of 844), down from 507.**  The "~100" above took XCPolicy's
       reach as the proxy; the floor is the union of the modules that hold a recipe or a `MeshMaker` BY DESIGN: the Hamiltonian factory
       122, `GPW_Evaluator` 120, `GPW_IBS` 118 (the §4.4 ruling's `MeshMaker`), `Imposition` 113, `Seed` 99, XCPolicy 99, `SCFIterator`
       78, the structure factory and the setups.  Below it would mean taking the `MeshMaker` back out of the GPW basis.

Rough size: about 60 files, most of them mechanical.
- **Behaviour changes:**
  - the duplicate A1 mesh build disappears (bit-identical);
  - the site-block check widens (a throw only on a broken builder).
- **Run-report effect:** the `auxiliaryFit` mesh scope leaves the molecular A1 runs. A console byte test may need its expected text, and
  that is a run-sequence improvement, not a refresh.
