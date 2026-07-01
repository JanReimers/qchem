# qcSymmetry structure — SALC pipeline + QN-hierarchy notes

Working notes on simplifying the recently-built molecular-symmetry code. **Companion to
`project_symmetry_naming_cleanup` (memory breadcrumb)** — the QN-hierarchy rename lives there; this doc
covers the SALC/ShellRep pipeline structure and coordinates the one overlapping change (the
`Irrep`+`Orbital_QNs` module merge, which belongs with that rename).

## The `Structure + RawBasisSet → SALC_IBS` call flow is clean

Traced `PG::SymmetryAdapt` end to end. It's a shallow, linear pipeline — the ShellRep DIP refactor already
removed the branching, so there is no tangle to unwind:

![Structure + RawBasisSet to SALC_IBS dataflow](diagrams/salc_call_flow.svg)

```
Structure ─► StructureToSymPoints ─► BuildAbelianGroup ─► AbelianGroup {char table, ops}
RawBasisSet ─► ExtractAoShells ─────────────────────────► AoShell[] {geom, norm, ShellRep}
                          │
   (AbelianGroup, AoShell[]) ─► BuildSALCs  ─► SALCs {O-matrix, irrep labels}
        BuildSALCs internally:  ∀ g: BuildOperationRep(shells, R) → shell.rep→Rep(R)   (no angular branch)
   (rawIBS, SALCs) ─► SymmetryAdaptedBasisSet ─► per-irrep SALC_IBS blocks
```

Each step is one well-named function with a clear output type. The algorithmic modules
(`PointGroup`/`AbelianGroup`/`CharacterTable`, `SALC`, the `*Rep` files) are appropriately separate and do
real work. **No structural rework needed here.**

## The one real smell — the `dynamic_cast` dispatch in `SymmetryAdapt`

`SymmetryAdapt` picks the shell extractor by casting the abstract `Real_OIBS` to the concrete `PGData` /
`SphData`:

```cpp
if (auto* pg  = dynamic_cast<const PGData*>(ibs))  shells = ExtractAoShells(*pg);
if (auto* sph = dynamic_cast<const SphData*>(ibs)) shells = PG_Spherical::ExtractAoShells(*sph);
```

This is (a) the abstract→concrete cast CLAUDE.md flags, and (b) the **combinatorics landmine**: every new
`{engine}×{angular}` delivery adds another `if (cast)`, and libcint-spherical is a *trap* precisely because
its cast to `PGData` silently succeeds with the wrong (Cartesian) components.

**Fix (DIP, same move as ShellRep): a virtual `GetAoShells()` on the orbital IBS.** Push the extraction into
the type that owns the data. Then:
```cpp
shells = ibs->GetAoShells();          // no cast, no PGData/SphData knowledge in the orchestrator
```
- Kills the abstract→concrete cast; `SymmetryAdapt` becomes basis-agnostic.
- A new basis (libcint-spherical S3b, or anything) becomes pluggable by implementing one method — it can't be
  a silent trap because there's no cast to succeed wrongly.
- Layering is fine: the orbital-IBS interface already depends on qcSymmetry (via `Irrep`), so returning
  `vector<Symmetry::AoShell>` adds no new edge.
- The two current `ExtractAoShells` (`PG_Cart`, `PG_Spherical`) share a `ShellTypeId` + block-grouping
  skeleton; a small shared helper dedupes it as they become the two method bodies.

**✅ DONE.** New abstract `qchem.BasisSet.Molecule.AoShellSource` (pure `GetAoShells()`); `PG_Cart::Orbital_IBS`
and `PG_Spherical::Orbital_IBS` implement it (delegating to their existing `ExtractAoShells` on their own
`PGData`/`SphData` base). **Follow-up (LSP): `AoShellSource` now REFINES `Real_OIBS`** (`: public virtual
Real_OIBS`), so `SymmetryAdapt` iterates `Iterate<const AoShellSource>` — the iterator does the encapsulated
cast, and since an `AoShellSource` is-a `Real_OIBS` the raw IBS falls out by plain upcast, eliminating the
client-side `dynamic_cast`-filter entirely (the "type-check instead of polymorphism" smell). Kept off the
templated `Orbital_1E_IBS<T>` base on purpose (AO shells are real-space, meaningless for the dcmplx
plane-wave instantiation). `SymmetryAdapt` now does one abstract→abstract cast and has zero PGData/SphData
knowledge; the evaluators stay symmetry-free (AoShellSource lives in qcMolecule_BS, which already links
qcSymmetry — no new edge). (Also fixed a latent `AoShell.monomials` breakage in M_PGSymmetry left by the
ShellRep refactor — that target wasn't in the build loop.) Pairs naturally with the parked
`NotImplemented`/`UnsupportedCombination` exception idea.

**Follow-up regression fix (`8254c42c`):** the refactor initially added `AoShellSource` to only the *MnD*
bases, silently regressing **`engine=LibCint` + Cartesian + symmetry** (the facade allows it; the old code
matched libcint via the `PGData` cast). No test caught it — every symmetry test used `engine=MnD`. Fixed:
`PG_LibCint::Orbital_IBS` now implements `AoShellSource` (Cartesian delegates to `PG_Cart::ExtractAoShells`;
spherical *throws* — its real-harmonic convention isn't matched, S3b, and its `PGData` base holds the
Cartesian layout, so reading it as Cartesian would be the silent trap). Added the missing
`{engine}×{symmetry}` guard test `M_Calculation.WaterSymmetryLibCint`. **158 UTMain green.** Lesson: the
`{engine}×{angular}×{symmetry}` matrix needs coverage — this is the exact landmine the parked exception idea
is about.

## Coordination with `project_symmetry_naming_cleanup`

That breadcrumb owns the QN-hierarchy rename: `Symmetry::Symmetry` → the `Irrep ⊂ SpinIrrep ⊂ OrbitalQNs`
naming (name by the added quantum number) + the `CarriesSpin()` predicate fix for the spinor double-group.
One structural change belongs with it, not here:

- **Merge `Irrep` + `Orbital_QNs` into one `qchem.Symmetry.OrbitalQNs` module.** They're a linear hierarchy
  (`Orbital_QNs : Irrep`) that can't exist apart. BUT `Irrep` is imported across ~20 BasisSet files, so the
  import-rename ripple should happen **once**, together with the naming rename — do it as part of
  `project_symmetry_naming_cleanup`, not as a standalone merge.

The core QN interfaces (`Spin`/`Symmetry`/`Irrep`/`Orbital_QNs`) are now doxygen'd (LaTeX formulas via
`\f$…\f$`) so the hierarchy browses while that rename is designed.

## Minor / low-interest (not prioritized)

Pure `.C`-file consolidation inside `Internal/Imp/` (merging tiny same-module impl units) — low value.
Phase 1 (the three spherical-QN impls `Yl`/`Ylm`/`Okmj` → one `Internal/Imp/Spherical.C`) is done because it
was zero-risk/zero-API; the rest (`Unit`+`BlochQN`) isn't worth the churn on its own and is not planned.
