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

There are three distinct problems; (B) is the serious one.

## A. The key is atom-biased — missing geometry  (worked around)

Every quantity is keyed by `IBS_ID_t(RadialID, AngularID)`.  For an **atom** that is complete (the
centre is pinned at the nucleus).  For a **molecule** it is not: overlap, kinetic and the 2-/3-
centre electron integrals are all orientation-dependent, but `RadialID`/`AngularID` historically
encoded only exponents and angular momenta — no atom centres.  Only `Nuclear` was geometry-aware
(it also keys on `cl->ID()`, and `Atom::ID()` includes `R=<pos>`).

Symptom: two same-basis molecules at different orientations collided — a rotated water reused a
canonical water's ERIs and gave a non-variational -146 Ha instead of -76.02.

Work-around in place: `PGData::RadialID()` now folds radial+centre+polarization per function into
the key and `PGData::AngularID()` is empty (commit "make the molecular-basis integral cache key
geometry-aware").  **Principled fix:** an explicit geometry/centre component in the key, not
centres stuffed into "RadialID"; and retire the radial/angular split for the molecular path (it
only made sense for atoms).

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

- Re-do **DFT + SALC** properly: the decorator (3C transform by O, fit-basis delegation, derive
  `SymmetryAdapted_IBS` from `Orbital_DFT_IBS`) was correct and gave matching sym/non-sym energies;
  it just needs the cache to be re-entrant so it can cache (not recompute) the raw 3C.  The diff
  was reverted out of `src/BasisSet/SymmetryAdapted_IBS.C` + `.../Imp/...` and `Orbital_DFT_IBS.C`.
- The molecular-symmetry **DIIS** is a *separate* topic (see doc/SCF_DIIS_SALC_notes.md): empty-
  irrep handling + the shared-coefficient overshoot, independent of the cache.
