# Integral cache (DB_Cache) — issues to fix

Consolidated from the molecular-symmetry work, where the cache was hit repeatedly.  The cache is
`IntegralsCache_RAM` (`src/BasisSet/Internal/Imp/DB_Cache_RAM.C`), interface `IntegralsCache_Base`
(`src/BasisSet/Internal/DB_Cache.C`).  `theGlobalCache` is a process-global singleton.

## Interface improvement  (DONE)
Replaced the `Has()`/`GetXXX()`/`Set()` protocol with a single self-contained
`Get(type, ids…, make)` call, where `make` is a `std::function<Result()>` lambda that
captures any extra args (cluster, mesh, partner IBS).  Client code is now e.g.
```c++
return cache->Get(I2C::Kinetic, {RadialID(),AngularID()}, [this]{ return MakeKinetic(); });
return cache->Get(I2n::Nuclear, {RadialID(),AngularID()}, cl->ID(),
                  [this,cl]{ return MakeNuclear(cl); });
```
All 16 call sites migrated; the old `Has`/`Set`/`GetVec`/`GetSMat`/`GetMat`/`GetERI3`/`GetERI4`/
`SetDirect`/`SetExchange` API and the `Ix1`/`Ix2` variants are deleted.

### Original notes (kept for context)
Typical client code looked like this
```c++

return cache->Has(IntegralsCache_Base::I2C::Overlap,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()))
    ? cache->GetSMat() : cache->Set(MakeOverlap());

```
or
```c++

 return cache->Has(IntegralsCache_Base::I2n::Nuclear,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()),cl->ID())
        ? cache->GetSMat() : cache->Set(MakeNuclear(cl));

```
where cl is Cluster&.  This puts a burden on the client to make the if statement.  The reason for this is that I don't how to
pass the MakeXXX(sometimes 0 sometimes 1 argument) function into the cache.  If there is a c++20 (or earlier) way to do this, then  we can reduce the burden on the client to something like
```c++
 return cache->Get(IntegralType,UniqueIDForBasisSet,CallThisFunctionIfYouDontAlreadyHaveIt);
or
 return cache->Get(IntegralType,UniqueIDForBasisSet,CallThisFunctionIfYouDontAlreadyHaveIt,ClusterOrSOmeOtherArgument);

```
If you can think of a way to this that would be great.  I would like establish this first as it may impact our descision on how
to solve problems A&B below.

## Scope: How much stuff is worth caching?
Another thing to think about: Perhaps caching 1-electron integrals is a waste of effort, only cache 3 and 4 center ERI integrals.  One
option to deal with this is at construction the DBCache gets a list of bools telling which type to cache and which ones are direct methods (always calls the MakeXXX function) so to speak.  The other is just rip out all the 1E function calls.  Interested in your opinion.  Maybe you know the big commercial programs do.
Keep in mind we want go beyond atoms and molecules, and use this cache for solids as well. 

There were three distinct problems (B was the serious one); **all three are now fixed** — see each
section.  This doc is kept as the record of what they were and how they were resolved.

## A. The key is atom-biased — missing geometry  (FIXED)

**Fixed:** the cache no longer assembles the per-basis key from `(RadialID, AngularID)`.  Identity
is now the **client's** job, via a one-method interface owned by the cache:

```c++
struct DBCacheClient { virtual std::string BasisSetID() const = 0; };
```

`Get(...)` takes `const DBCacheClient*` (two of them for cross-IBS quantities) and uses
`client->BasisSetID()` verbatim as the per-basis key axis — it knows nothing about atoms /
molecules / solids.  ID assembly lives where the geometry knowledge is:
- `IrrepBasisSet_IDs` (inherits `DBCacheClient`) provides the default `RadialID()+"|"+AngularID()`
  — complete for an **atom** (centre pinned at the nucleus).  `RadialID`/`AngularID` are retained
  as the atom's building blocks; the radial/angular split is no longer the cache key.
- `PGData::BasisSetID()` (molecular) folds radial @ centre : pol per function — geometry-aware.
  `PGData::RadialID()`/`AngularID()` are restored to clean basis-only strings (the geometry hack is
  gone).  The molecular IBS overrides `BasisSetID()`; otherwise it would fall back to the
  geometry-free default and collide again.
- `SymmetryAdapted_IBS::BasisSetID()` = `raw->BasisSetID() + "[label]"` (per-irrep), replacing the
  old `AngularID()`-with-label hack; `RadialID`/`AngularID` delegate to the raw basis unchanged.

`IrrepBasisSet_IDs::GetID()` (used as basis identity by `IrrepCD`/`FittedFunction`) now delegates to
`BasisSetID()`, so there is a single identity source.

Original symptom (now covered by tests): two same-basis molecules at different orientations
collided — a rotated water reused a canonical water's ERIs and gave a non-variational -146 Ha
instead of -76.02.  `M_PG_Sym` rotated/translated water + DFT all pass.

## B. Has()/Set() is not re-entrant — shared stateful "last key"  (FIXED)

**Fixed** by the `Get(type, ids…, make)` redesign above.  Each `Get` does its own `find`,
computes `make()` into a local (which may itself perform nested cached `Get`s), then `insert`s —
holding no find-iterator across `make()`.  No `itsLastKeyX`/iterator members remain, so the
clobbering scenario below can no longer occur.  The 4-centre GC's "don't evict what I just made"
protection now takes the freshly-inserted key as a parameter (`RunGarbageCollector(protect)` /
`Purge(…, protect)`) instead of reading shared `itsLastKey4a/b`.

### Original analysis (kept for context)
`Has(...)` did the map lookup AND stashed the looked-up key and iterator in **mutable members** of
the cache, e.g. `itsLastKey2`, `itsLastKeyx`/`its2xIterator`, `itsLastKey3`/`its3CIterator`,
`itsLastKeyn` (DB_Cache_RAM.C:140-233).  `Set(value)` then relies on those stashed members to file
the value in the slot `Has()` found.

This breaks under **any nested cached access** between an outer `Has()` and its matching `Set()`:

```
outer.Has(K_outer)   // miss; stashes itsLastKeyX = K_outer, iterator = end()
  ... compute the value, which itself does:
  inner.Has(K_inner) // stashes itsLastKeyX = K_inner   <-- CLOBBERS the outer stash
  inner.Set(v_inner)
outer.Set(v_outer)   // files v_outer under K_inner / the wrong iterator  <-- CORRUPTION
```

Consequences seen this session:
- The SALC decorator's `MakeOverlap()` (it transforms the raw overlap) tripped it; worked around by
  calling the raw **compute** `raw->MakeOverlap()` instead of the cached `raw->Overlap()`.
- The **fitted-DFT 3-centre path** tripped it hard: symmetric-vs-non-symmetric water DFT gave
  *non-deterministic* totals (-75.86 / -75.56 / -75.66) depending on access/run order.  Worked
  around in the decorator by building the transformed 3C from the raw **compute** `MakeOverlap3C`/
  `MakeRepulsion3C` (so no nested cached `Has/Set`) — but that path is not committed; it is the
  reason DFT+SALC is parked until this is fixed.

**Fix:** make `Has`/`Set` self-contained — no shared `itsLastKeyX`/iterator state.  Options:
- have `Has()` return the slot (iterator / optional reference) and pass it to `Set()`, or fold
  lookup-or-insert into a single `GetOrCompute(key, lambda)` call (the cleanest — eliminates the
  Has-then-Set protocol entirely), or
- at minimum, have `Set(value, key)` take the key and do its own single `insert`/`operator[]`.

Once (B) is fixed, the decorator can use the **cached** raw 3C (fast) instead of recompute, and
multiple geometries / sym-vs-non-sym runs can share a process safely.

## C. Order / cross-run pollution  (FIXED via B)

Should be gone now that (B) is fixed — confirm by restoring DFT+SALC (next section) and
re-running the sym/non-sym and atomic-DFT/HF tests in one process.

### Original analysis (kept for context)
Because of (B), running several computations with the same basis in one process can corrupt later
ones (the DFT sym run after the non-sym run, or unrelated atomic DFT tests after the SALC tests).
The atomic DFT/HF tests in `UTMain` are sensitive to this.  Fixing (B) should remove it.

## After the DBCache fix — return to the symmetry work

- **DFT + SALC — DONE.**  `SymmetryAdapted_IBS` now derives from `Orbital_DFT_IBS<double>` and
  overrides only the compute hooks `MakeOverlap3C`/`MakeRepulsion3C` (+ fit-basis delegation to the
  raw DFT basis); the cached `Overlap3C`/`Repulsion3C` accessors are inherited unchanged.  The
  compute hooks transform the raw basis's **cached** 3C (`itsRawDFT->Overlap3C(c)`), now safe under
  the re-entrant `Get` — so the raw 3C is computed once and shared by every irrep, and the
  transformed block is cached per irrep under its own AngularID.  The earlier raw-COMPUTE
  work-around (and the per-fit-basis local memoization) is gone; the 1-e path was likewise switched
  from `itsRaw->MakeOverlap()` back to the cached `itsRaw->Overlap()`.  Tests:
  `M_PG_Sym.water_DFT_unpolarized` (matches to <1e-5) and `water_DFT_polarized` (≈2e-4 residual from
  the SALC DIIS/occupation topic below, not the transform).  (The reverted WIP lived in
  `git stash@{0}`; reconstructed cleaner here, so that stash is now superseded.)
- The molecular-symmetry **DIIS** is a *separate* topic (see doc/SCF_DIIS_SALC_notes.md): empty-
  irrep handling + the shared-coefficient overshoot, independent of the cache.
