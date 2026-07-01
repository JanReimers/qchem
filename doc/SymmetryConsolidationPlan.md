# qcSymmetry consolidation — fewer modules/files for the quantum-number layer

Plan for simplifying `src/Symmetry` (qcSymmetry). Prompted by the observation that the atom-symmetry
classes "just store quantum numbers, no functionality," spread across many tiny module/.C files. Living
document; Phase 1 shipped, the rest staged for review because the higher-value merge touches public API.

## Diagnosis (survey)

qcSymmetry is ~37 `.C` files but only ~800 lines of real logic. Two populations:

- **Algorithmic core — KEEP SEPARATE.** `PointGroup` (geometry + detection, ~350 L impl), `CharacterTable`,
  `AbelianGroup`, `SALC`, and the rep files `ShellRep`/`CartesianRep`/`SphericalRep`/`OperationRep`. These do
  real work (rotation algebra, character projection, Gram-Schmidt, the c2s projection). Leave them.
- **Quantum-number storage — OVER-FRAGMENTED.** ~12 small classes that hold 1–3 fields and implement only
  `SequenceIndex()` (an ordering key) + `Write()`: `Spin`, `Irrep`, `Orbital_QNs`, `MolecularIrrep`,
  `UnitQN`, `BlochQN`, and the four spherical/spinor containers `Yl`/`Ylm`/`Ωκ`/`Ωκmj`. Each tends to get its
  own module and/or `.C` file. This is the fragmentation to collapse.

The QN classes form the intended hierarchy *by added quantum number* (see `project_symmetry_naming_cleanup`):
`Irrep` (spin + spatial `sym_t`) ⊂ `Orbital_QNs` (+ principal n). `sym_t = shared_ptr<const Symmetry>` is the
polymorphic spatial-symmetry handle whose concretes are `MolecularIrrep` / `Yl`/`Ylm` / `Ωκ`/`Ωκmj` / `BlochQN`
/ `UnitQN`. `SequenceIndex()` is purely an ordering/caching key — no physics.

## Staged consolidation

Ordered by value/risk. Internal-only merges first; the public-API one last (it's also entangled with the
queued naming cleanup).

### Phase 1 — merge the 4 spherical/spinor QN impls into one file  ·  ✅ DONE  ·  zero API change
`Yl.C`, `Ylm.C`, `Okmj.C` are all implementation units of the SAME module
(`qchem.Symmetry.Internal.Spherical`) — three `.C` files for one module. Merge their bodies into a single
`Internal/Imp/Spherical.C`; the interface `Internal/Spherical.C` (which exports `Yl`/`Ylm`/`Ωκ`/`Ωκmj`) is
unchanged, so nothing outside sees a difference. CMake: 3 impl lines → 1.

### Phase 2 — fold `UnitQN` + `BlochQN` into one module  ·  low risk (internal, via Factory only)
Unlike Phase 1, these are separate *modules* (`qchem.Symmetry.Unit`, `qchem.Symmetry.BlochQN`), each a
concrete `Symmetry::Symmetry` used **only** through `Factory` (`UnitFactory`/`BlochFactory`). Merge into one
`qchem.Symmetry.Internal.BasicQNs` (both classes) + one impl file. Only `Imp/Factory.C` imports them, so the
blast radius is one file. Module-name change, hence "low" not "zero."

### Phase 3 — merge `Irrep` + `Orbital_QNs` → one `qchem.Symmetry.OrbitalQNs` module  ·  MODERATE (public API)
They're a linear hierarchy (`Orbital_QNs : Irrep`) that can't exist apart, so one module is natural. BUT
`Irrep` is imported widely across BasisSet (`IrrepBasisSet`, `BasisSet`, the atom evaluators, tests), so this
is an import-rename ripple (`import qchem.Symmetry.Irrep;` → `…OrbitalQNs;`) across ~20 files. Mechanical but
wide. **Do this together with the queued `project_symmetry_naming_cleanup`** (rename `Symmetry::Symmetry`→the
`Irrep ⊂ SpinIrrep ⊂ OrbitalQNs` naming + the `CarriesSpin()` predicate fix) so the public surface churns
once, not twice. Left for review — it's the one with real blast radius and a naming decision attached.

### Not merging (by design)
`PointGroup`, `SALC`, `CharacterTable`, `AbelianGroup`, and the `*Rep` files — real algorithms with clear
boundaries. `Factory` + `Spherical` (abstractions + `Getl()/Getκ()` cast helpers) stay separate: `Spherical`
is the "what can I query from a symmetry?" contract used across the atom basis, and folding it into `Factory`
would muddle construction vs. dispatch.

## Doxygen pass (browse-enabling)

The newer rep files (`ShellRep`/`CartesianRep`/`SphericalRep`/`OperationRep`) already use `//!`. The QN
interfaces (`Spin`, `Symmetry`, `Irrep`, `Orbital_QNs`, `Spherical`, `MolecularIrrep`) and the algorithmic
impls (`PointGroup`, `CharacterTable`, `AbelianGroup`, `SALC`) are mostly plain `//`. Convert the QN
interfaces to `//!`/`\brief` as each module is touched, so the hierarchy renders in Doxygen. (Done for the
Phase-1 file + the central QN interfaces.)

## Related parked item — an "unimplemented feature" exception

Separate but adjacent (GUI-team-facing): the facade's feature combinatorics — `{Cartesian,Spherical}` ×
`{MnD,libCint}` × `{symmetry,no}` × model — has holes (spherical-libcint SALC = S3b, deliberately not built).
Today those throw `std::runtime_error`. Worth a small dedicated `NotImplemented`/`UnsupportedCombination`
exception the GUI can catch, carrying the specific combination (engine, angular, symmetry, model) + a
suggested working alternative, so an unbuilt corner is a clean caught error, not a surprise. Spec TBD.
