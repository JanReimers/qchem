# Molecular basis set upgrade

## Reading order

This document is split into two parts:

- **Near-term, decided work** — Stages 1–3 (PG refactor → cache_4 → cache_3/DFT). This is what we
  execute next, in order.
- **Longer-term direction** — Goals A–F (generic BasisSet-level code, new basis types,
  BasisSetSource). This is the *big picture* the near-term work is clearing a path toward. It is
  deliberately less pinned down; treat it as intent, not commitment.

Do the stages first. Everything in the "Longer-term direction" section assumes the PG tree has
already been simplified and proven.

## Optimize for clarity, not speed (whole plan)

The driving goal of this work is to get the OOD/SOLID framework right and demonstrate that new
integral algorithms and basis-set types can be added with small amounts of shared code. Where a
choice trades a little runtime for a lot of clarity, take the clarity. Any perf concern raised below
is flagged so we don't *mistake the refactor for a regression* later — not because we should
pre-optimize.

## Regression gate (applies to every stage)

No stage deletes or renames the old tree until these stay green:

- Build target is **`UTMain`** (other targets — `UTSym`, `UTBasisSet` — have pre-existing failures;
  do not gate on them).
- Molecular PG tests must reproduce the old energies to tolerance **before** the old tree is
  scrubbed:
    - `src/BasisSet/Molecule/tests/M_PGSymmetry.C` (`M_PG_Sym.*`, incl. `water_DFT_{un,}polarized`)
    - `src/BasisSet/Molecule/tests/M_PG_MeshIntegrals.C`
    - `src/BasisSet/Molecule/tests/libCint.C`
- HF **and** DFT paths both pass (see Stage 1 note — DFT is not optional for the scrub gate, because
  DFT+SALC consumers exist today).

## Dependency note

Stages 2 and 3 build directly on the just-landed **DBCache** work (defects A/B/C, client-owned
`BasisSetID()` geometry key, re-entrant `Get`). Any new cache (`Cache2`, `PG_cache4`) must obey that
identity/invalidation model — it is already decided. See `doc/DBCache_issues.md` and the
`project_integral_cache_geometry_key` memory. In particular: the cache key is the client's
`BasisSetID()`, the cache no longer assembles `(RadialID, AngularID)` itself, and `Get(type, ids…,
make)` is self-contained / re-entrant.

---

# Near-term work

## Stage 1: PG refactor — collapse the multiple dispatch

### Problem (verified)
The PG code emulates multiple dispatch over combinations of radial-function types in the integral
routines. C++ has no native multimethods, so this is hand-rolled as chained virtual re-dispatch:
each `Integrate` overload peels off one argument, swaps `this` into the now-resolved slot, and
re-calls — then a `dynamic_cast` recovers the concrete type once all arguments are resolved.
`GaussianRF.C` alone has **17 `dynamic_cast`s** and 8 `Integrate` overloads. With 4-index integrals
and 2 concrete radial types this is a 2⁴ combinatorial fan-out. This was over-application of the
dependency-inversion principle — too clever by half.

### Insight (verified)
There are only **two** concrete radial types: `GaussianRF` (primitive) and `ContractedGaussianRF`.
Collapse to **one** radial type that always carries contraction coefficients; a primitive is simply
a length-1 coefficient vector (`{1.0}` or the normalization constant). This removes the type
combinatorics — and therefore the entire dispatch/`dynamic_cast` machinery — at the root.

### Approach: parallel tree (full copy, not shared)
Build a new tree `src/BasisSet/Molecule/PolarizedGaussian1`
(namespace `BasisSet::Molecule::PolarizedGaussian1`) and grow the new single-radial mechanism from
scratch. **Copy everything else over wholesale** — `MnD/`, `Readers/`, `Reader.C`, `Symmetry.C`,
polarizations — into PG1 with the new namespace, rather than sharing it with the old tree.

Copying (not sharing) is deliberate: it means the two trees cannot collide on module names, the old
tree's `dynamic_cast`-laden code can't accidentally pull in the new radial type (or vice-versa), and
— most importantly — we are free to **tweak the MnD interface in PG1** without touching or risking the
old tree. The duplication is temporary: the old tree is scrubbed at the Stage 1 exit. (Trying to
share one MnD across both radial designs would be the same too-clever coupling we're removing.)

### The seam (designed 2026-06-19)
The old code has **two parallel primitive-vs-contracted type splits**, both encoding the same idea:
- Radial: `GaussianRF` (one exponent) vs `ContractedGaussianRF` (holds `gs` = vector of primitive
  `RadialFunction*`, `cs` = coeffs).
- Hermite: `GaussianH3` vs `ContractedGaussianH3` (holds vector of primitive `Hermite3*` + coeff ref;
  its `operator()` already *is* "sum of primitive blocks × contraction coeffs").

The entire `Integrate` peel-off dispatch + `dynamic_cast`s exist only to resolve which concrete type
each argument is, so as to reach the **primitive** M&D kernel — which already exists as standalone
statics taking concrete pointers: `GaussianRF::Integrate3C(...)` and `Integrate4C(...)` (the latter is
the `lambda`/RNLM/Hermite math).

**Single-type design (keep only the contracted form):**
- One radial type holding a list of `(exponent, coeff)` primitives; an uncontracted function is a
  length-1 list with `coeff = normalization`.
- `Integrate` collapses from six virtual overloads to **one non-virtual function**: loop over each
  function's primitives, accumulate `c_a·c_b·…·primitiveKernel(α_a, α_b, …)` by calling the existing
  `Integrate3C`/`Integrate4C` math per primitive-tuple. No peel-off, no re-dispatch, no `dynamic_cast`.
- `ContractedGaussianH3` stops being a type — contraction becomes the accumulation loop. Primitive
  `Hermite1`/`GaussianH3` (per primitive-pair) survive unchanged; they are what Stage 2 caches as Ω_ab.
- `GetH3` / the `Hermite3` base collapse to the primitive form only.

**Net:** the primitive math is untouched; we delete the dispatch scaffolding and replace it with a
contraction loop. Structure the loop so "compute/lookup Ω for primitive pair (i,j)" is a clean inner
call — that keeps Stage 2 (cache relocation) a localized change rather than another rewrite. Since the
copied MnD lives only in PG1, reshape it freely to fit.

### Performance note
Folding primitive + contracted into one type means the hot inner kernel iterates a length-1
coefficient vector with the attendant indirection where the old primitive path had none. M&D's cost
lives in the Hermite recursion, so this is expected to be negligible — but it is a clarity-over-speed
choice made on purpose (see top-of-doc note). Do not later read it as a regression without measuring.

### DFT in Stage 1? — Decided: **yes, carry it from the start (option a)**
DFT is not an afterthought. The scrub gate requires DFT regressions to pass, DFT+SALC consumers exist
today, and maintaining two divergent trees across three stages is exactly the kind of cleverness this
refactor is trying to retire. So PG1 carries the 3C/DFT path from the start, and "proven" at the
Stage 1 exit means HF **and** DFT both reproduce old-tree energies. (Stage 3 then becomes "move the
3C path onto DB_Cache," not "add DFT.")

### Verification
PG1 is verified two ways, which should agree:
- **Against old PG** — H2O and N2, HF **and** DFT, reproduce old-tree energies within tolerance. This
  is the primary gate (it's the apples-to-apples comparison through the same code paths).
- **Against libcint / libintx** — the `src/BasisSet/Molecule/tests/libCint.C` test already provides an
  independent integral reference. Cross-checking PG1 integrals against libcint (and/or libintx)
  catches errors that *both* PG and PG1 might share (e.g. a long-standing M&D convention bug), which
  the PG-vs-PG1 comparison cannot. Use it as the independent oracle.

### Exit criteria
Both verifications pass and the regression gate is green. Only then: scrub the old PG, rename
`PolarizedGaussian1` → `PolarizedGaussian`, move to Stage 2.

## Stage 2: New PG — relocate Ω caching into DB_Cache

Following M&D, the code uses a charge distribution Ω_ab (`GaussianCD` in
`src/BasisSet/Molecule/PolarizedGaussian/Internal/CDCache.C`) to store data for a primitive pair
ab. **The Ω_ab caching algorithm already exists** in `CDCache` — this stage *relocates ownership* of
that storage into the `src/BasisSet/Internal/DB_Cache.C` framework. It is a migration, not a new
caching scheme.

- Add `src/BasisSet/Internal/Cache2.C` for 2-index caching. The open design question is "**and maybe
  looping**": decide whether `Cache2` mirrors `Cacheable4`, or whether 2-index Ω_ab storage is just a
  map/array keyed by primitive pair. `CDCache` is the existing reference for what's actually needed.
- New `class PG_cache4 : public virtual Cacheable` holds pointers to Ω_ab and Ω_cd plus whatever else
  the 4C kernel needs. Determine that "whatever else" by reading
  `src/BasisSet/Molecule/PolarizedGaussian/Radial/Imp/GaussianRF.C` →
  `GaussianRF::Integrate4C`. (`double lambda = 2*Pi52/(ab.AlphaP*cd.AlphaP*sqrt(ab.AlphaP+cd.AlphaP));`
  // M&D 3.31 is one candidate per-4-index quantity — confirm against the rest of the kernel.)
- **Cache key & lifetime must follow the DBCache model** (see Dependency note): the 4-index key is
  derived from the clients' `BasisSetID()`, and entries live under the same geometry-aware
  invalidation already in place. PG_cache4 does not invent its own identity.

## Stage 3: New PG — use cache_3 for DFT

Largely mechanical once Stage 2's caching shape is settled: bring the 3C (DFT fit) path onto the same
DB_Cache framework. Reference `SymmetryAdapted_IBS`'s `MakeOverlap3C` / `MakeRepulsion3C` (the raw 3C
is computed once and shared across irreps — that arrangement must be preserved).

**After Stage 3, pause and reassess** before starting the longer-term work below.

---

# Longer-term direction (big picture — intent, not commitment)

Once the PG tree is clean and as much code as possible lives at the generic BasisSet level (out of
the `Molecule/` and `Atom/` subtrees), the following become tractable. These are ordered roughly, but
expect to re-order as the stages above teach us what the right abstractions are.

**A. Use `Cache4` (`src/BasisSet/Internal/Cache4.C`) for molecular integrals.** Store all data for a
four-index basis-function combination by deriving from `Cacheable4`. Used in the 4-index HF
Direct/Exchange loops.

**B. Lift the Evaluator idea to PG basis sets.** Move
`src/BasisSet/Atom/Evaluators/Evaluator.C` (minus `Getl()` and the `Angular` class) out to
`src/BasisSet/Internal/Evaluator.C` so PG can reuse it.

**C. Hoist the 1E integral loops to be basis-set-agnostic.** Move `Integrals_Overlap`,
`Integrals_Kinetic`, `Integrals_Nuclear` (from `src/BasisSet/Atom/IrrepBasisSet.C`) out to
`src/BasisSet/Orbital_1E_IBS.C`. The i,j matrix-building loops should be basis-set agnostic.
  - **⚠ Spike before committing:** the "even for solids" claim is unproven. The Lattice tree has its
    own `IBS_Evaluator` (plane waves, complex-valued, k-point phase factors). Confirm a single 1E
    loop signature genuinely unifies Atom / Molecule / Lattice — k-point phases may break the
    "just loop i,j" assumption — *before* this becomes committed work.

**D. Hoist the ERI loops** in `Orbital_HF_IBS<E>::MakeDirect` / `MakeExchange` out to
`src/BasisSet/Orbital_HF_IBS.C`.
  - **⚠ Spike (highest-risk item in the back half):** the atomic exchange `Ak` calculation may not
    fit a generic ERI loop. If it doesn't, this either doesn't happen or forces a new abstraction —
    and that abstraction would constrain the Evaluator interface from **B**, so it must be designed
    *before* **B** is finalized, not discovered during **D**. Do this spike early.

**E. Add more molecular basis-set types.** If A–D are done right, each new type is mostly just a new
Evaluator:
  - `SphericalGaussian` via M&D
  - `PolarizedGaussian` / `SphericalGaussian` via [libcint](https://github.com/sunqm/libcint)
    and/or [libintx](https://github.com/ValeevGroup/libintx) instead of M&D
  - Room for more modern algorithms (e.g. ACE — Yanai et al., *Int. J. Quantum Chem.* **76**,
    396–406 (2000)). Many integral algorithms exist (PH, HGP, OS, DRK, M&D-PRISM, HGP-PRISM); the
    point is a framework where each lands in a small amount of shared code, not implementing them all.

**F. Define a `BasisSetSource` interface**, with implementations to:
  - read basis sets from the local filesystem (start with GAUSSIAN94 —
    `src/BasisSet/Molecule/PolarizedGaussian/Reader.C`)
  - pull from [Basis Set Exchange](https://www.basissetexchange.org/) /
    [basis_set_exchange](https://molssi-bse.github.io/basis_set_exchange/)
