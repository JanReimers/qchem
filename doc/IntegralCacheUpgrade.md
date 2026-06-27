# Integral Cache Upgrade — spec (3 stages)

> **STATUS: DONE (all 3 stages).** Implemented and verified against the full UTMain regression (134
> tests green after each stage). See the per-stage **✅ DONE** notes below. Highlights:
> - Stage 1 accessor lives in a new same-module impl unit `src/BasisSet/Internal/Imp/DB_Cache.C`
>   (declared in the `DB_Cache` interface) so callers keep importing only the lightweight interface —
>   no import churn, single emitted instance.
> - Stage 2 needed the `<dcmplx>` explicit instantiation in **both** Imp TUs (`DB_Cache_RAM.C` and
>   `DB_Cache_RAM_Report.C` — report/GC members live in the latter; the split-instantiation pattern
>   already used for `<double>`).
> - Stage 3's `override` on `CacheDim()` surfaced `-Winconsistent-missing-override` warnings; fixed
>   `IrrepBasisSet<T>::GetVectorSize` and 10 pre-existing missing-`override`s in `PlaneWave_IBS.C`.

Self-contained plan for a fresh session. Motivated by two bugs the multi-species plane-wave work
exposed (both already FIXED in passing, see commits `606a54ff`, `11aefdc3`), which revealed the cache
is overdue for a proper lifecycle + a complex sibling + a self-check.

## Background / current state

- `theGlobalCache` is a raw global pointer `IntegralsCache<double>* theGlobalCache = 0;`
  (`src/BasisSet/Internal/DB_Cache.C:91`), **double-only**, that every `main()` must `new` by hand:
  `BasisSet::theGlobalCache = new BasisSet::IntegralsCache_RAM<double>(true);`
  (in `UnitTests/gtestmain.C`, `src/BasisSet/Atom/tests/gtestmain.C`, `.../A_Pool.C`,
  `src/BasisSet/Molecule/bench/MnDBench.C`, `UnitTests/scfrun.C`; and **commented out** in
  `src/Common/tests/gtestmain.C` — a latent null-deref). It also leaks (never deleted).
- The cache keys by `DBCacheClient::BasisSetID()` (a geometry string), NOT the `this` pointer
  (`DB_Cache_RAM.C:213+`). Good — distinct geometries don't collide.
- The **complex (plane-wave) path has no cache**: `theGlobalCache` is `IntegralsCache<double>` and the
  cache stores `smat_t<T>` (symmetric), but the PW overlap is **Hermitian** (`hmat_t<dcmplx>≡chmat_t ≠
  smat_t<dcmplx>`). So `Integrals_Overlap<dcmplx>::Overlap()` (`src/BasisSet/Imp/IrrepBasisSet.C`, the
  `template <>` specialization) lazily buffers `MakeOverlap()` in a `static map` — was keyed by the
  `this` pointer (stale-cache bug `11aefdc3`, now keyed by `BasisSetID()` as a stopgap).

Type facts: `hmat_t<double> ≡ smat_t<double> ≡ rsmat_t`; `hmat_t<dcmplx> ≡ chmat_t`; `smat_t<dcmplx>` is
complex-SYMMETRIC (≠ chmat_t). The `Make…()` makers (`MakeOverlap`/`MakeKinetic`/…) already return
`hmat_t<T>`.

---

## Stage 1 — Lifecycle (construct-on-first-use). Lowest risk; the "both of us got bitten" fix.  ✅ DONE

Replace the manual global with a **templated Meyers accessor** in `DB_Cache_RAM.C` (where
`IntegralsCache_RAM` is complete), exported:

```cpp
template <class T> IntegralsCache<T>& theCache()      // built on first use, destroyed at exit
{
    static IntegralsCache_RAM<T> instance(/*makelog=*/false);
    return instance;
}
```

- Fixes: no `main()` boilerplate (impossible to forget), no leak (function-static destructs at exit and
  still prints the RAM-usage report from its dtor), thread-safe init (C++11 magic statics).
- Migrate callers `theGlobalCache->Get(…)` → `theCache<T>().Get(…)`:
  `Imp/IrrepBasisSet.C` (Overlap), `Imp/Orbital_1E_IBS.C` (Kinetic, Nuclear), `Imp/Fit_IBS.C`
  (Charge/Repulsion/InvOverlap/InvRepulsion/Norm), `Imp/Orbital_DFT_IBS.C` (Overlap3C/Repulsion3C),
  `Internal/Imp/Orbital_DHF_IBS.C`.
- Delete the 5 `theGlobalCache = new …` lines from the mains; delete the global (`DB_Cache.C:91`).
- **Gotcha:** stage 1 keeps DOUBLE-only behaviour — only `theCache<double>()` is instantiated. Do NOT
  reference `theCache<dcmplx>()` yet (the `smat_t<dcmplx>` storage won't compile until stage 2). The
  dcmplx `Overlap()` stopgap stays. Confirm the templated callers (Kinetic/Nuclear) are not
  ODR-instantiated for `dcmplx` (the PW path calls `MakeKinetic()` directly and uses the Overlap stopgap).
- The `(true)` ctor arg = `makelog` (verbose per-make log); default it `false`.

## Stage 2 — Complex cache (`smat_t<T>` → `hmat_t<T>`).  ✅ DONE

Migrate the **symmetric** caches (`I2C`, `I2n`) from `smat_t<T>` to `hmat_t<T>`. Byte-identical for
double (`hmat_t<double>==smat_t<double>`, and the makers already return `hmat_t<T>`); Hermitian-correct
for dcmplx (`chmat_t`).

- `DB_Cache.C` + `DB_Cache_RAM.C`: change `Get(I2C,…)` and `Get(I2n,…)` return type + their
  `std::function<smat_t<T>()>` make-param to `hmat_t<T>`; change the internal storage maps
  (`itsSMats`, `itsNMats`) to `hmat_t<T>`. Leave `I2x` (cross, non-symmetric → `mat_t<T>`), `I1C`
  (vectors), `I3C`/`I4C` (ERI) unchanged.
- Add `template class IntegralsCache_RAM<dcmplx>;` and make the dcmplx instantiation compile (check the
  `I2x`/`ERI3`/`ERI4` members for the dcmplx path; PW only needs `I2C`).
- Route the dcmplx `Overlap()` through `theCache<dcmplx>().Get(I2C::Overlap, this, [this]{return
  MakeOverlap();})` and **delete the static-buf stopgap** in `Imp/IrrepBasisSet.C`. Then the general
  `Integrals_Overlap<T>::Overlap()` template serves both — consider deleting the `template <>` dcmplx
  specialization entirely.

## Stage 3 — Self-checking cache (the `CacheDim()` guard). Independent of 1/2; arguably do FIRST as insurance.  ✅ DONE

On a cache **hit**, verify the cached matrix's size matches what this client expects; on mismatch, dump
`{operator, BasisSetID(), cached dims, expected dim}` and throw. Would have caught NaF(113)/CsI(251)
instantly at the cache boundary instead of a Cholesky segfault three layers down. Catches any
`BasisSetID()` that isn't specific enough *whenever the colliding geometries differ in size* (the common
case; same-size-different-content still slips — rarer, separate problem).

- `DBCacheClient` (`DB_Cache.C`): add `virtual size_t CacheDim() const = 0;`
  **NAME IT `CacheDim()`, NOT `size()`** — renaming `GetNumFunctions`/`GetVectorSize` → `size()` is a
  known compiler-noise minefield in this `VectorFunction`/diamond hierarchy (history: it was attempted
  and abandoned). A distinct name sidesteps it entirely.
- `IrrepBasisSet<T>` (`src/BasisSet/IrrepBasisSet.C:58`): add the **single bridge**
  `virtual size_t CacheDim() const override { return GetNumFunctions(); }`. This is the final overrider
  for EVERY concrete cache client — they are all `IrrepBasisSet<T>` (the molecular geometry-aware
  `BasisSetID` is delegated *through* the basis, so `PGData` etc. are not standalone clients). The
  abstract mixins (`Integrals_Overlap<T>`, …) stay `CacheDim()`-pure, which is fine (abstract). Net
  noise: one bridge, possibly a one-liner or two the compiler points out — NOT a flood. (No `HasSize`
  super-base — it'd add a third parallel diamond root for no gain; `IrrepBasisSet<T>` is the natural
  single bridge.)
- It MUST live on `DBCacheClient` (not a cache-side `dynamic_cast` to `IrrepBasisSet`): the cache can't
  see `IrrepBasisSet` — that module imports the cache, so importing it back is a module cycle.
- In `IntegralsCache_RAM::Get` (`DB_Cache_RAM.C`): on hit verify `cached.rows()==bs->CacheDim()` (I2C/I2n);
  for the cross `I2x`, `rows==a->CacheDim() && cols==b->CacheDim()`. Also assert on insert.

## Suggested order

3 (cheap insurance) → 1 (lifecycle) → 2 (complex). Or 1 → 2 → 3 per the original framing. Each
green-gates independently against the full UTMain regression (Si Γ −7.2273, NaF −20.3293, CsI −11.3868).
