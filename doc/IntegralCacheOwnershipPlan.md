# Integral-cache ownership: one caching framework, owned locally

**Status: LIVE PLAN (user rulings 2026-11-12, §8), to be executed by lane H after PADE, stage by stage (S0-S6).**
Each stage is its own READY, with "no number moves" checked by `scripts/anchordiff` at 1e-10 on short and long, plus panel4t at S4
(it touches the Hamiltonian terms).  M2 is measured before S4.  Written 2026-11-12 as a design note.  It was revised the same day for
the user's three follow-ups, relayed by the controller: the 3C owner, lifetime, threading, instance identity.

**Asked by the controller (2026-11-12), relaying the user's direction.**  Build ONE caching framework, owned locally: a cache is an
object with an owner, and there are no keyed global singletons.  BasisSet and Hamiltonian use it the same way.

The user's comments (verbatim, via the controller):
- On "reuse between runs at the same geometry": *"when I conceived of this I was thinking about runs that did all atoms in the
  periodic table 1->118.  This is not a very important capability."*
- *"'it's what exists' I don't give that any weight ... we want to get to the right place."*
- On sharing between two equal, independently built bases: *"I think that capability should be very low priority."*
- On the split along the geometry line: *"the split seems to have a negligible gain ... We need a really good reason to justify
  the inconsistency and added complexity of supporting two caching modes."*

This note finds NO such reason.  §4 measures the one candidate, the process-wide primitive caches (Cache2/3), and finds no gain.
So everything goes local.

Background: the exact-key fix (da67bfd27, `doc/Records/OpenWork_History5.md` "A cache key is exact") closed today's collision.  This
note removes the KEY itself: an owned cache needs no string identity, so the collision class goes away by construction.

---

## 1. Every cached quantity, what it depends on, and its proposed owner

Today all of the string-keyed entries below live in ONE process-wide object, `IntegralsCache_RAM<T>` (`theCache<T>()`, T = double
or dcmplx, `src/BasisSet/Internal/DB_Cache*.C`).  The primitive caches (Cache2/3/4) are registered in the same object.

| entry (today's slot) | who fills it today | depends on | proposed OWNER |
|---|---|---|---|
| `I2C::Overlap` | every IBS (`Imp/IrrepBasisSet.C`) | one block | **the block itself** (its parent's store, keyed by the block) |
| `I2C::Kinetic` (⟨p²⟩) | `Orbital_1E_IBS` | one block | the block.  ⚠ The Kinetic TERM also keeps ½⟨p²⟩ in its own static cache, so T is stored TWICE today; §7 |
| `I2C::RestMass` | `Orbital_DHF_IBS` | one block | the block |
| `I2x::Kinetic` (RKB L/S cross) | `Orbital_DHF_IBS` | the large block and its own small-component partner | the DHF block (it owns the partner) |
| `I2n::Nuclear` | `Orbital_1E_IBS::Nuclear(cl)` | block × the nuclear potential of `cl` | **the Ven TERM.**  It already caches its matrix per irrep (`tStatic_HT_Imp`); the basis-side copy is a second copy |
| `I2n::LocalPP` / `LocalPPLong` / `LocalPPShort` / `SeparablePP` | `GPW_IBS` (k-block) | k-block × PP species × structure | **the Ven_PP terms**, which also cache their blocks already.  `LocalPPShort` is keyed WITHOUT the grid on purpose, so a rebuilt basis (a grid ladder stage) reuses it.  Owned locally, that reuse becomes an explicit hand-over, if a measurement says it is worth one (§7, M2) |
| `I1C::Charge`, `I1C::Normalization`, `I2C::Repulsion`, `I2C::InvOverlap`, `I2C::InvRepulsion` | `Fit_IBS` | the fit basis alone | **the fit basis** (members of `Fit_IBS`, an Internal class) |
| `I2x::Repulsion` (fit × fit) | `Fit_IBS::Repulsion(b)`; today only `ConstrainedFF` with b = itself | two fit bases | the fit basis (the cross with itself is a one-block entry) |
| `I3C::Overlap`, `I3C::Repulsion` (⟨ab\|c⟩) | `Orbital_DFT_IBS::Overlap3C/Repulsion3C(fit)`, cached in `theCache<TFit>` | orbital block × fit basis | **the orbital parent**, keyed by (block, the fit's INSTANCE ID).  ⟨ij\|f⟩ = ⟨f\|ij⟩ is stored once (§2, the user's rule) |
| `I4C::Direct`, `I4C::Exchange` (ERI4, canonical block pair only) | `Orbital_ERI4_IBS` | two blocks of ONE parent basis | **the parent BasisSet** (§2) |
| `I1C` + mesh, `I2x` + mesh (`itsmVecs`, `itsmMats`) | nobody | -- | **delete**: no caller since `Fit_IBS::Norm` went analytic (2026-11-11) |
| `Cache4` (atom Rk Slater integrals; exponent grouper inside) | radial atom bases, by `RadialType` | the radial exponents, shared by every l block | **the atom BasisSet** (its l blocks share it) |
| `Cache2` "PG.Omega" (Ω of a primitive pair), "PG.RNLM", `Cache3` "PG.H3" | MnD kernels via `findΩ` / `findRNLM` / `findH3`, through a static pointer memo in `GaussianRF.C` | primitive pairs/triples (exponent, centre, l); content-interned ids | **the parent molecular BasisSet** (its blocks and the fills it runs share them; §4: nothing is lost) |
| primitive registry (the content-interned `PrimGaussian` index) | process-wide, never freed | (exponent, centre, l) | the same owner as Cache2 (ids local to it) |
| `GeometryCacheBudget` (RAII size-0 budget on the global Cache2/3 for lattice sums) | lattice drivers | -- | **goes**: a lattice driver owns a scratch primitive cache for its image loop and drops it at scope end.  The budget becomes an object's lifetime |
| ERI4 RAM cap (4 GB, oldest-first GC with timestamps) | the global cache | -- | **goes**: resident memory is bounded by the LIVE bases.  The cap becomes a setup POLICY: predicted ERI4 bytes against RAM, WARN AND RUN or REFUSE ON RESOURCES (CLAUDE.md policy rules) |

The Hamiltonian side already caches locally.  Today it has two shapes:
- `tHT_Common::itsCache`: a `map<Irrep, hmat>` with pre-created slots, filled on first read for static terms, and rebuilt at
  `RefreshForDensity` for dynamic ones.
- Some terms key by `bs->BasisSetID()` STRINGS: `HF_HT`'s J/K blocks and `PWTerms_PP`.

The fitters hold their own memo layers, keyed by the block's `sym_t` and locked:
- `DeltaScalarFitter::Handle`;
- `DensitySampler_OverlapProjected::Integrator`;
- `ConstrainedFF`.

These are the "same framework" the user asks for, written four ways.

## 2. Who owns each slot, and what is left of its key

**The user's rule** (2026-11-12, via the controller): *"If the orbital basis caches all the 3C integrals then so f's need an ID but
ij do not."*
- An owner never keys ITSELF.  A slot's key is whatever is left once the owner is implicit.
- ⟨ij|f⟩ and ⟨f|ij⟩ are one tensor, stored once, on the ORBITAL basis.  This reverses my first draft, which put the 3C on the fit
  basis.
- The fitters' memo layers keyed by `sym_t` then become reads of that one slot (§7).

| slot | owner | residual key | needs geometry? |
|---|---|---|---|
| Overlap, Kinetic ⟨p²⟩, RestMass | the parent's cache | the block index | no |
| RKB L/S cross kinetic | the parent's cache | the block index (the small partner belongs to the block) | no |
| ERI4 Direct / Exchange | the parent's cache | an UNORDERED pair of block indices (canonical i ≤ j): no string | no |
| ⟨ab\|c⟩ Overlap3C / Repulsion3C | the parent's cache | (block index, the fit basis's INSTANCE ID) | no |
| fit × fit `I2x::Repulsion` | the fit basis's own cache | the other fit basis's instance ID; today it is always itself, so nothing is left | no |
| fit basis's own Charge, Norm, Repulsion, InvOverlap, InvRepulsion | the fit basis | nothing | no |
| Nuclear, LocalPP (all three ranges), SeparablePP | the Ven / Ven_PP TERM | the block's `Irrep`, as the term caches today | no |
| the δ fit's Φ tables, the mesh-dependent XC data | the XC term / its fitter (`DeltaFit_IBS` already holds them per block) | the block | no |
| Cache4 Rk | the atom parent's cache | exponent indices from ITS grouper | no |
| Ω, RNLM, H3 (primitive pairs and triples) | the molecular parent's cache | primitive indices from ITS registry | no |

**No residual key needs geometry.**  Every geometry dependence is carried by the owner: a new geometry is a new BasisSet, so it is
a new cache.  The one entry that WANTS to cross owners is `LocalPPShort` across a grid ladder.  That becomes an explicit hand-over,
if measurement M2 (§7) says it is worth one, never a key.

The δ and multipole fits carry no geometry in their `BasisSetID` today.  Under instance IDs that stops mattering.  We do NOT put a
δ fit's point set in any key: the user, *"we don't want that, the key would be huge!!"*.  Two δ fits paired with one orbital basis
get two different instance IDs.

The element type of a 3C slot follows the FIT axis.  A real TRIM block paired with the run's complex G-space fit stores a
`PairProjection<dcmplx>`, so the parent's cache holds 3C slots of both scalars.  It does this the way `BasisSetImp` holds blocks of
both scalars: one typed map per scalar, never a cast.

## 3. Construction, lifetime, threading

**Lifetime (SETTLED by the user, 2026-11-12).**
- *"The cache becomes a constructor argument of the irrep basis set.  The actual BasisSet will make a cache."*
- The parent BasisSet's factory constructs the cache and passes it to every block it builds.  There is no default and no null
  (CLAUDE.md: constructors pass concrete objects; factories concentrate the wiring).
- A unit test holds one cache in its harness, a fixture member, and passes it in.  No block exists without one.
- The handle is a `std::shared_ptr<BasisCache>` (name open).  The parent and the blocks co-own it, so no lifetime order can
  dangle.
- A fit basis made by `CreateCDFitBasisSet` / `CreateVxcFitBasisSet` gets its OWN small cache from that factory, for its
  one-basis entries.  Its 3C rows live in the orbital parent's cache.

**Threading.**  The user asked: *"is it possible to lock small parts of the caches as opposed to the whole thing?"*  Yes.  There
are three options, each matched to what an entry kind needs.

| option | what it is | fits |
|---|---|---|
| (a) per-entry once-init | A SHORT lock on the slot map to find or insert the entry, then the fill OUTSIDE any lock.  Other askers of THAT entry wait on it alone (a per-entry once-flag or future); askers of other entries never wait.  It keeps today's re-entrancy: a fill may do nested cached gets, because no lock is held while computing | **S, T, ERI4, 3C, the fit metrics.**  Few, big entries.  Their fills are threaded INSIDE already (3c18237f4), so contention on the map lock is a few acquisitions per run |
| (b) sharded maps | N maps, each with its own lock, chosen by key hash | nothing here.  It wins only with many concurrent askers of many small entries, which is exactly the case (c) serves without a lock |
| (c) resolve serially, then read immutably | Before a parallel region, resolve every entry the workers need into an immutable table; workers take no lock | **Ω / RNLM / H3 in the hot loops.**  This is what 3c18237f4 does with `PrimitivePairTable`.  A lock per primitive lookup is a check in the hot loop, against CLAUDE.md's "check where it is cheap" |

So the cache takes (a) for its slots.  The primitive caches stay single-thread by contract: only a fill's serial prepass touches
them.  That contract is checked by an owner-thread assert (live in RelChecked) and a TSan gate.  This also closes the open "Cache2/3
thread safety" row by construction.  A `const` accessor on a `mutable` cache is honest: every slot is a pure function of fixed
inputs (CLAUDE.md's lazy-caching rule).

## 3b. The partner's identity: an instance token, not a pointer and not a string

Only one kind of residual key names another OBJECT: the fit basis in a 3C slot (and a second fit in `I2x`).  It is an INSTANCE
identity.

**The token.**
- Every basis that can be a partner owns `std::shared_ptr<const InstanceToken>`, stamped at construction from a process-wide
  `std::atomic<uint64_t>` serial and never reused.
- The cache stores, beside each such key, the serial and a `weak_ptr<const InstanceToken>`.

**A serial, not a GUID.**  A GUID buys uniqueness ACROSS processes.  No cache key outlives its process (the store is RAM only), so a
64-bit serial is enough: it is cheaper, readable in a report, and ordered.  If a key ever has to be persisted, a GUID is the
upgrade, and that is a different feature (a content hash, as in the June note below).

**Why not the pointer.**  An address is reused after the object is freed, so a new partner at the old address gets the dead one's
entries: a stale hit.  It has happened here.
- In June 2026 the Hamiltonian's dynamic-term cache keyed on the density's POINTER.
- A freed density's address was reused, and `M_DFT_Water` diverged at iteration 5, nondeterministically (allocator-dependent; clean
  under valgrind).
- d5a23e36 stopped the reuse (shared_ptr densities).  582f1d05 fixed the stale `GetMatrix`.  f354f890 replaced the pointer with a
  logical-clock serial stamped at construction (`NextDensityVersion`).
- In the user's words, it *"took a loooong time"* to find.  This token is the same cure, applied before the bug instead of after.

**Copies.**
- A copy constructs a NEW token, so it gets a new serial.  The token member's copy constructor mints; it does not share.
- Blocks are not copied today: there is no IBS-level `Clone`.  The only `Clone`s are MnD Hermite/Block internals.
- So the alternative, "an immutable exact copy may share the serial", buys nothing.  It would also re-open the question an ID
  exists to close.

**A dead partner: EVICT and report.**
- A serial is never reused, so an entry keyed by a dead partner can NEVER be hit again.  It is pure waste, and the cache can tell
  on its own (`weak_ptr::expired()`): no destructor coupling, no back-pointer, no deregistration call to forget.
- The cache sweeps its partner-keyed entries (a handful: one or two fit bases per orbital basis) on every insert, and at the
  owner's report point.
- It writes the count and bytes freed, with each dead partner's name and serial, to the run report under the caller's scope
  (`report::Set`, one level in), and one line to the console.
- A warning without eviction would keep memory that can never be read.  So evict.
- In a normal run the fits die WITH their orbital basis, so a non-zero count is itself informative: something rebuilt a fit
  basis mid-run.

**Gate (S3).**
1. Build a fit basis, fill its 3C on an orbital block, destroy the fit.
2. The orbital cache reports and frees EXACTLY those entries (the count and the bytes).
3. Build a second fit basis, at the same address if the allocator obliges.  Placement-new into the freed storage forces it.
4. It gets a different serial, MISSES, and fills fresh.

**Where content keys survive (pin 34's domain).**  Instance tokens cover every residual key that names another object.  Content
identity is left only where content IS the question: deduplication INSIDE one owner.
- the primitive registry (equal (exponent, centre, l) across a basis's blocks share Ω);
- the atom grouper (equal exponents across l blocks share Rk);
- `GaussianRF::TypeID` (SALC grouping of equal radials);
- the PG_Cart PP memo (`KeyOf`).

Those stay exact (`KeyString`, %.17g): an approximate content key there would merge two primitives.  Everything else in pin 34's
audit table goes away with its global key (`Structure::ID` as a cache key, `MeshParams::ID` as one, the `BasisSetID` writers).

## 4. Cache2/3 per owner: what is lost (MEASURED)

A probe (UTGaussian_BS, NOT committed; it lived in a scratch edit of `M_ThreeCBlocks.C`):
- built every integral a molecular DFT/HF point needs, serially: S, both Dunlap 3C tensors, ERI4 J and K;
- ran the H₂O₂ geometries below in one process, one after another;
- before each build, cleared every string-keyed entry, so all builds did the SAME work;
- varied only how warm the global Ω cache was.

`RNLM` and `H3` were never registered.  The grouped 3C and ERI4 fills keep their RNLM and H3 tables local, so only Ω is a live
global cache on the molecular hot path today.

| H₂O₂, serial | dzvp: build s | new Ω | resident Ω | dzvp-molopt: build s | new Ω | resident Ω |
|---|---|---|---|---|---|---|
| A, reference (Ω cold) | 0.331 / 0.355 | 932 | | 20.4 / 33.5 | 2695 | |
| A again (Ω warm, same work) | 0.328 / 0.358 | 0 | | 20.9 / 31.7 | 0 | |
| FD step, one O moved 1e-4 | 0.332 / 0.355 | 537 (42 % reused) | | 37.1 / 31.0 | 1323 (51 % reused) | |
| FD step, one H moved 1e-4 | 0.333 / 0.330 | 204 (78 % reused) | | 41.2 / 26.4 | 931 (65 % reused) | |
| relaxation step, every atom moved | 0.372 / 0.379 | 932 (0 % reused) | | 40.2 / 30.6 | 2695 (0 % reused) | |
| after 12 geometries | | | **5.8 MB, never freed** | | | **28.8 MB, never freed** |

Two repetitions are shown, "first / second".  The MOLOPT timings swing 2× because the box was shared with other lanes' runs.  The
dzvp ones are steady.

**The verdict.**
- Warm against cold Ω at IDENTICAL work: 0.331 against 0.328 s on dzvp, and within noise on MOLOPT.  Building Ω is under 1 % of
  a point's integrals.  The cost is the contraction, which no cross-geometry cache can share.
- A finite-difference step reuses 42-78 % of the Ω ENTRIES and saves no measurable TIME.
- A relaxation step reuses nothing: every atom moves.
- Multi-k needs no cross-geometry cache: the lattice paths already run Cache2/3 at a size-0 budget (`GeometryCacheBudget`), so
  they keep nothing today.
- Meanwhile the global Ω grows without bound (0.6 MB per dzvp geometry, 3 MB per MOLOPT one), and "the registry caches are NOT
  clearable" (`DB_Cache_RAM.C`).

So there is no "really good reason" to keep Cache2/3 global.  Per owner, the within-geometry reuse stays (that is where Ω's 99.9 %
hit rate comes from: across irreps and quartets of one point), and the dead geometries go.

⚠ This CORRECTS my proposal in `OpenWork_History5.md` ("Cache2/3 ... real reuse across one-atom FD steps ... stay global").  The
reuse is real in COUNT and worth nothing in TIME.

Not measured: `Cache4` (atom Rk) across atom runs.  Its only cross-run beneficiary is the 1→118 sweep, which the user rates as
"not a very important capability".

## 5. The hub signatures that change, staged so each step lands green

Every stage computes the same integrals with the same arithmetic.  So the gate at each stage is:
- `scripts/anchordiff` at 1e-10 shows **0 moves**;
- the short and long sweeps pass, and RelChecked passes.

No number may move.  A moved number is a bug in the stage, not roundoff.

| stage | what | hub? |
|---|---|---|
| **S0** | Delete the dead mesh-keyed overloads.  Key `HF_HT`'s J/K blocks and `PWTerms_PP`'s by `Irrep`, not the `BasisSetID` string.  Measure M2 (§7) | no |
| **S1** | The framework: `SlotMap<Key,V>` (Common, beside `LazyLock`).  It holds find-or-create under a lock, a per-slot once-fill, a byte count for the report, and a REFRESH discipline for density-dependent slots (`Clear` + rebuild in `RefreshForDensity`; the block read stays pure, V1.41).  Unit tests: a threaded fill is the serial one, nested fills do not deadlock, a refresh slot is never filled lazily.  Then `tHT_Common::itsCache` moves onto it | no (Hamiltonian internals) |
| **S2** | `Fit_IBS`'s own entries (Charge, Norm, Repulsion, InvOverlap, InvRepulsion, self fit × fit) become members of the fit basis | no (Internal) |
| **S3** | **Hub.**  `IrrepBasisSet_IDs : DBCacheClient` becomes a cache-client face: `BasisSetID()` and `CacheDim()` retire, replaced by the block's index in its parent.  Every concrete block constructor DEMANDS a `shared_ptr<BasisCache>`, which the parent's factory makes and a test harness holds (§3).  A partner basis (a fit) owns an `InstanceToken` (§3b).  The bodies of `Overlap()`, `Kinetic()`, `RestMass`, the RKB cross, `Direct`/`Exchange` and `Overlap3C`/`Repulsion3C` switch to the cache.  **The accessors keep their signatures**: only the base face changes.  Gates: the §3b dead-partner and address-reuse gate, and the threaded once-fill.  `HubInterfaceAudit` snapshot in a commit citing the ruling | **YES**: `qchem.BasisSet` (IrrepBasisSet), `Orbital_DFT_IBS` |
| **S4** | **Hub.**  Nuclear and PP leave the basis: `Orbital_1E_IBS::Nuclear(cl)` (the cached accessor) and the cached PP face of `GPW_IBS` go.  The Ven/Ven_PP terms call the uncached `Make*`, and their existing static cache is the one copy.  This also fixes a naming defect: `tGPW_IBS::MakeSpeciesFieldMatrix` / `MakeProjectorMatrix` are `Make`s that CACHE (CLAUDE.md: Get caches, Make computes).  If M2 says the ladder reuse matters, an EXPLICIT hand-over of the short-PP blocks from one stage's term to the next, a named phase (pin 11) | **YES**: `Orbital_1E_IBS` |
| **S5** | Cache2/3/4 and the primitive registry move into the parent store.  The MnD kernels take the cache from the fill that calls them, instead of from `GaussianRF.C`'s static memo.  A lattice driver's image loop owns a scratch primitive cache, and `GeometryCacheBudget` is deleted | no (BasisSet internals) |
| **S6** | Delete `IntegralsCache_RAM`, `theCache<T>`, `qchem.BasisSet.Internal.DB_Cache` and the ERI4 GC.  Each store reports its bytes under the CALLER's report scope (the basis section): the callee writes one level in, replacing the global `EmitIntegralsCacheReport` snapshot.  The ERI4 RAM cap becomes a setup policy (predicted bytes; WARN AND RUN / REFUSE ON RESOURCES, reported under `policy`) | the report face only |

The Hamiltonian hub (`qchem.Hamiltonian`) needs NO face change.  Its terms already own their caches, and S1 changes only what
type those caches are.  `pybind/` does not read `BasisSetID` (grep, 2026-11-12).  Any other compile break there at S3 goes to lane G;
do not edit `pybind/`.

## 6. How a relaxation frees the old geometry

- Each geometry point of `MoleculeCalculation` builds its own BasisSet and Hamiltonian.
- A point's `Carry` holds only matrices by VALUE (D, and the previous point's D and S).  It holds no basis.
- When the optimizer drops the previous point, its BasisSet, fit bases and terms die, and their stores go with them: ERI4, 3C, Ω,
  the PP blocks.  No eviction code and no GC: the memory belongs to the point.
- Gate (S6): an N-step relaxation whose resident integral bytes after every step are at most two points' worth (the live point and
  the one being built).  Today that number grows by one point per step.

The solid path does the same.  `SolidCalculation`'s geometry session builds each point's molecular basis fresh:
`point->mol=recipe.Make(*point->cell, …)` (`Imp/SolidCalculation.C` ~3225).  It hands that basis to `AtGeometry` through a
`shared_ptr` that the point owns, so it dies with the point.

## 7. Found along the way

- **T and V_nuc are stored twice.**  The basis caches ⟨p²⟩ and V_nuc, and the Kinetic and Ven terms each keep their own copy in
  their static cache.  S4 removes V_nuc's basis copy.  T's basis copy stays, because GPW's home-only path and
  `SphericalLatticeView` read `Kinetic()` too.
- **Dead code**: the mesh-keyed `Get` overloads (S0).
- **A cache keyed by a string inside a term**: `HF_HT` (`RefreshedBlocks(...).at(bs->BasisSetID())`) and `PWTerms_PP` (S0).
- **`Make` that caches**: `tGPW_IBS::MakeSpeciesFieldMatrix` / `MakeProjectorMatrix` (S4).
- **ERI4 canonical order by string**: `min(BasisSetID)`.  It becomes the block index (S3).
- **M2 (to measure at S0): the grid-ladder reuse of `LocalPPShort`.**  It is keyed without the grid so a later stage's rebuilt basis
  hits it.  Cost of the short PP against a stage's setup on NaF and MnO.  If it is small, S4 drops it; if it is large, S4 adds the
  explicit hand-over.
- **The never-freed primitive registry** and the never-clearable Cache2/3: the root of the "trail of dead geometries", gone at S5.

## 8. The user's rulings (2026-11-12, relayed verbatim by the controller)

These are the user's selected answers to the six questions this section posed.
1. **One local framework, nothing process-global** (Cache2/3/4 per owner): *"Approve all three (Recommended)"* (covers 1, 3, 4).
2. **The handle**: *"shared_ptr (Recommended)"*.  The parent BasisSet makes the cache; every IBS constructor takes it; tests hold
   one in the harness.
3. **Partner identity**: an instance serial in a shared token, held by the cache as a `weak_ptr`; dead-partner entries EVICTED and
   reported (§3b).  Approved under 1.
4. **Threading**: once-init for slots, serial-resolve-then-immutable for the primitive hot loops, no sharding (§3).  Approved under 1.
5. **The S3 and S4 hub changes**: *"Approve S3 and S4 (Recommended)"*.  This is the design-discussion ruling for the `BasisSet`
   hub, to be cited in the commit that regenerates `scripts/hub-interfaces.snapshot`.
6. **The ERI4 cap**: *"Setup policy (Recommended)"*: WARN AND RUN / REFUSE ON RESOURCES at setup, not a silent GC.

The §7 findings fall into the stages:
- S0: the dead mesh-keyed `Get`s, and the `HF_HT` / `PWTerms_PP` string keys.
- S4: V_nuc's second copy, and the caching `GPW_IBS` `Make*`.
- T's second copy (½⟨p²⟩ in the Kinetic term beside ⟨p²⟩ on the basis) stays: one N×N matrix per block, accepted.
