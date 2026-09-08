# Spherical SALC Plan — symmetry adaptation for the spherical-harmonic bases

Design notes for extending point-group SALC adaptation from the Cartesian PolarizedGaussian
basis (`PG_Cart`, done) to the two **spherical** deliveries: `PG_Spherical` (in-house MnD
spherical) and libcint-spherical. Living document. Companion to `doc/MolecularSymmetryPlan.md`
(the group-theory engine, which this reuses wholesale).

## Why

Today `PG::SymmetryAdapt` works **only** for the Cartesian `PG_Cart_MnD::PGData` orbital IBS
(`src/BasisSet/Molecule/PG_Cart/Imp/SymmetryAdapt.C` — it `dynamic_cast`s each IBS to `PGData`
and asserts otherwise). The spherical bases are unsupported:

| `Molecule::Factory` delivery | orbital IBS type | `SymmetryAdapt` today |
|---|---|---|
| default (MnD, Cartesian) `PG_Cart` | IS-A `PGData` | ✅ works (`M_Sym`) |
| `engine=libcint`, Cartesian | IS-A `PGData` | ✅ works (same Cartesian layout) |
| `angular=spherical` `PG_Spherical` | IS-A `SphData` (NOT `PGData`) | ❌ asserts "no PolarizedGaussian orbital IBS" |
| `engine=libcint`, spherical | IS-A `PGData`, but 2l+1 spherical components | ⚠️ **trap** — cast succeeds, `ExtractAoShells` reads the wrong (Cartesian) shells; not caught |

Spherical is the physically cleaner basis (each shell expands into its 2l+1 real solid harmonics,
dropping the Cartesian d/f s-contaminant), so it is arguably the *more valuable* symmetry target —
yet `M_Sym` only proves the Cartesian path.

## The one thing that is angular-specific

The SALC pipeline is angular-representation-specific in exactly **one** place:

> `CartesianShellRep(R, exps)` (`src/Symmetry/CartesianRep.C`) — how a shell's Cartesian
> monomials `p_a(R⁻¹u)` map back onto the shell's monomials, giving the (#comp × #comp) rep matrix.

Everything downstream is angular-agnostic and **reused unchanged**: the center permutation,
`BuildOperationRep`, `BuildSALCs` (irrep projection), and the `SymmetryAdapted_IBS` decorator
(it just applies the SALC transform `O_Γᵀ·M·O_Γ` to matrices). So the task reduces to: produce a
correct per-shell rep matrix in the spherical components, in the basis's own m-ordering.

## Key reuse: a cart→spherical transform already exists

`src/BasisSet/Molecule/Evaluators/PG_Spherical_MnD/SolidHarmonics.C` already expresses each real
solid harmonic χ_{l,m} as a Cartesian expansion (`SphericalShell`, `CartTerm`), for l = 0..3
(s,p,d,f). Call its (2l+1 × n_cart) coefficient matrix `C_l`. Then the spherical operation rep is

    D_sph(R) = C_l · CartesianShellRep(R) · C_l⁺          ( C_l⁺ = right pseudo-inverse )

This is **exact** (not an approximation): the spherical subspace is rotation-invariant, so the
Cartesian contaminant is projected out cleanly. → the new math is cheap and leans on already-tested
code, instead of a from-scratch real-Wigner-D implementation.

## Work breakdown (dependency order)

### Stage S1 — spherical shell rep  ·  `qcSymmetry`  ·  low risk
- Add `SphericalShellRep(R, l)` → (2l+1)×(2l+1) real-Wigner-D, via `C_l · CartesianShellRep · C_l⁺`.
- Decide where `C_l` lives: it currently sits in the `PG_Spherical_MnD` evaluator. Either (a) lift the
  pure-math c2s coefficients down into `qcSymmetry` (keeps `qcSymmetry` self-contained, LAPACK-free —
  the pseudo-inverse is a fixed small matrix, can be precomputed/orthonormal so `C_l⁺ = C_lᵀ`), or
  (b) pass `C_l` into the rep builder from the basis side. **Lean (a)** to keep the convention with
  the math.
- Test: faithfulness `D(R1)D(R2)=D(R1R2)`; p ⇒ `D = R`; d/f against reference D-matrices. l≤3 only
  (matches `SolidHarmonics`; l≥4 ⇒ extend that table, it asserts past f today).

### Stage S2 — generalize `AoShell` / `BuildOperationRep`  ·  `qcSymmetry`  ·  ✅ DONE
- `AoShell` + `BuildOperationRep` MOVED from `CartesianRep` into a new module
  `qchem.Symmetry.OperationRep` (`export import`s both `CartesianRep` + `SphericalRep`).  Reason: the
  dispatcher needs both per-shell reps, but `SphericalRep` already imports `CartesianRep`, so it cannot
  live in either — it must sit ABOVE both.  `CartesianRep`/`SphericalRep` are now pure per-shell math.
- `AoShell` gained `HarmonicC2S c2s` (empty ⇒ Cartesian, set ⇒ spherical) + `IsSpherical()`/`nComponents()`;
  `BuildOperationRep` dispatches `SphericalShellRep(R,c2s)` vs `CartesianShellRep(R,monomials)` per shell.
  `c2s` is LAST in the struct so the existing positional aggregate initializers stay valid.
- `BuildSALCs` unchanged (consumes rep matrices).  Cartesian path byte-for-byte (M_CartesianRep/M_SALC/
  M_Sym green); new `OperationRep.spherical_d_shell_dispatch` test.  Import churn minimal: only `SALC.C`
  (re-export) + `M_CartesianRep.C` — every other consumer gets `AoShell`/`BuildOperationRep` transitively
  through `SALC`.

### Stage S3 — the two extractors  ·  basis tree  ·  THE REAL COST (two conventions)
The rep's component order must match the basis's m-order. The two spherical deliveries differ:
- **S3a `PG_Spherical`**  ·  ✅ DONE.  New module `qchem.BasisSet.Molecule.PG_Spherical.Symmetry` with
  `ExtractAoShells(const SphData&)`.  KEY SIMPLIFICATION realised while building: the extractor reads the
  basis's OWN Cartesian expansion straight out of `SphData::comps[].terms` (with `ns[]` as the per-harmonic
  norm) — it does NOT re-derive via `SphericalShell(l)`.  So the c2s handed to `SphericalShellRep` is, by
  construction, in the basis's exact m-order with the basis's exact coefficients — the "self-consistent, so
  lowest risk" note became "zero convention risk."  The predicted S3a hazard evaporated.
- **S3b libcint-spherical**  ·  ⬜ only remaining piece.  Extractor matching **libcint's** real-harmonic
  ordering *and* normalization (a foreign convention).  Bug-prone (sign/order); needs its own verification.
  Trap: libcint-spherical presents AS a `PGData` with spherical components, so the Cartesian extractor would
  silently misread it — the facade guards it out (`engine=MnD` only for symmetry+spherical) until S3b lands.

### Stage S4 — dispatch in `SymmetryAdapt`  ·  basis tree  ·  ✅ DONE
- `SymmetryAdapt` now `dynamic_cast`s each orbital IBS to `PGData` (Cartesian → existing extractor) then
  `SphData` (spherical → S3a); neither ⇒ a clear thrown error (libcint-spherical unsupported).  The
  `SymmetryAdaptedBasisSet` decorator is unchanged.  (These casts are abstract→concrete data structs — the
  established pattern for the data-carrying IBS types; flagged for the CLAUDE.md cast survey.)
- Facade guard relaxed: `{.symmetry=true}` + `angular=Spherical` allowed for `engine=MnD`, still blocked for
  `engine=LibCint` (S3b).
- Latent S2 gap fixed: `BuildSALCs` computed `nAO` from `monomials.size()` (0 for spherical) — now
  `nComponents()`, the angular-agnostic count.

### Stage S5 — tests  ·  ✅ DONE (mirror `M_Sym`)
- `M_Sym.water_HF_spherical_{unpolarized,polarized}`: spherical-adapted SCF total energy == the un-adapted
  spherical SCF, water/dzvp (the O d-shell exercises the harmonic path).  157 UTMain / 14 UTSymmetry green.

## Staging / shippability

Independent of the facade `{.symmetry=true}` flag (#2), which ships now **guarded to Cartesian**.
Each spherical mode is independently shippable:

1. **S1 + S2** (engine generalization) — no behavior change, pure capability.
2. **S3a + S4 + S5** for `PG_Spherical` — first user-visible spherical SALC (own conventions).
3. **S3b** for libcint-spherical — match the foreign ordering; lifts the last guard.

When S3a lands, the facade guard relaxes from "Cartesian only" to "Cartesian or MnD-spherical"; when
S3b lands, the guard disappears.

## Design note — `ShellRep` polymorphism (post-S5 refactor, DONE)

The S2 `AoShell` carried a discriminated union (`monomials` xor `c2s`, told apart by `IsSpherical()`), and
`BuildOperationRep` branched on it.  Resolved by Dependency Inversion: a new `qchem.Symmetry.ShellRep`
abstraction (`nComponents()` + `Rep(R)`), with `CartesianShellRep` / `SphericalShellRep` as its
implementations (they became classes holding `exps` / `c2s`).  `AoShell` now holds a
`shared_ptr<const ShellRep>` — no monomials/c2s, no angular flag, nothing to discriminate — and
`BuildOperationRep` just calls `rep->Rep(R)` (no branch).  Crucially the dependency inverts:
`OperationRep`/`SALC` depend ONLY on the `ShellRep` abstraction, not on the Cartesian/spherical concretes;
the per-basis extractors construct the concrete rep.  Diagnosis of the original flaw: `AoShell` was
conflating two jobs — a node in the permutation structure AND the source data for the angular rep; the
union was the symptom.  (Considered `std::variant` — rejected as a type-safe union that keeps the second
job in `AoShell` and still forces a `visit` branch, rather than removing the responsibility.)

## Risk notes

- Low-risk core (S1/S2 reuse `SolidHarmonics` + the tested Cartesian rep). The **m-ordering /
  normalization conventions** (S3) are the only real hazard — and there are two of them, which is
  why "both spherical versions" ≈ 2× the extractor+test work, not 1×.
- Keep `qcSymmetry` LAPACK-free (the existing invariant): `C_l⁺` is a fixed small matrix; precompute
  it orthonormal so no solve is needed.
- Higher angular momentum (l ≥ 4 / g) is out of scope until `SolidHarmonics` is extended.
