# ERI4 Rework: a custom ERI container that owns its contraction and bra–ket symmetry

Status: Stages 1+2 DONE (committed ffcd7f18 — `ERI4::MatMul` member + `ERI4::ScatterBoth` + unit test).
Stage 3 RESCOPED 2026-07-02 after a survey (see §5.4): the two-target driver is NOT in `ChargeDensity`, it
needs the cross-irrep view (the whole composite `BasisSet`) threaded from `CompositeWF` through the
Hamiltonian terms.
**Stage 3a DONE** (inert plumbing — threads `const tbs_t<T>*`, bit-identical — see §8.3a). **Stage 3b DONE (Coulomb/J + Exchange/K)**
(`Vee`/`Vxc`/`VxcPol` consume the context; whole-system `ScatterBoth`; `Jac` AND `Kab` RAM 521→369 MB —
see §8.3b/3b-K). Only the belt-and-suspenders canonical cache key (3c) remains. Author: design session
2026-07-02b; rescope + 3a + 3b(J+K) 2026-07-02.

## 1. Problem

`ERI4` is `SymmetricMatrix<SymmetricMatrix<double>>` (`ERI4T<double,smat_t>`,
[src/BasisSet/Internal/ERI4.C](../src/BasisSet/Internal/ERI4.C),
[src/BasisSet/Internal/ERI4T.C](../src/BasisSet/Internal/ERI4T.C)). It exploits **4 of the 8-fold** ERI
symmetry — `a↔b` (outer symmetric matrix) and `c↔d` (inner symmetric matrix) — but **not the bra–ket
swap** `(ab)↔(cd)`, i.e. `J(i,j) = J(j,i)ᵀ`. So for every pair of distinct irreps `i≠j` the code
computes and stores **both** `J(i,j)` and `J(j,i)`, doubling ERI storage and doubling the (expensive)
integral build. The `Jac`/`Kab` block cache is the dominant consumer of cache RAM.

The container is also "too clever by half": `operator()(a,b)` hands back a reference to a *whole inner
symmetric matrix* (`const smat_t&`). That exposed-storage interface is what blocks a lazy transpose view
(the transposed inner block is scattered across all stored inner blocks) — but it turns out production
doesn't need that interface at all (§3).

## 2. The symmetry we're leaving on the table, and the 4× trap

`J(i,j) = J(j,i)ᵀ` for both Direct and Exchange — already asserted by the suite
(`fnorm(J1, J1ba.Transpose())≈0` and the K analogue in A_Cache4 / DBCache). So storing only the canonical
`i≤j` block and deriving the other orientation is mathematically safe.

The trap: the contraction `MatMul` (§3) is **localized** in the `(cd)` direction —
`S_ab += Σ_cd J(a,b)⊙D_cd` reads the contiguous inner block `J(a,b)`. The transposed contraction
(`S_cd += Σ_ab J(a,b)·D_ab`, fixing `cd` and gathering `J(a,b)(cd)` across all `ab`) is a **scattered**
Schur product — empirically ~**4× slower**. So a "store canonical + serve the flip via a transpose
flag/view" design would tax every flipped block 4× **every SCF iteration** — a bad trade that gives back
in recurring CPU what it saved once in RAM.

## 3. Survey: the actual external contract (why this is deep-but-narrow)

Production touches an ERI4 through a **single choke-point**, not scattered indexing:

- **Contraction** — `MatMul(rsmat_t& Sab, const ERI4& J, const rsmat_t& Dcd)`
  ([ERI4.C:10-18](../src/BasisSet/Internal/Imp/ERI4.C)):
  ```cpp
  for (auto ia:iv_t(0,Nab))
      for (auto ib:iv_t(ia,Nab))
          Sab(ia,ib) += blazem::sum( gabcd(ia,ib) % Scd );   // S_ab += Σ_cd J(a,b)⊙D_cd
  ```
  Reached only via `Orbital_HF_IBS::AccumulateDirect/AccumulateExchange`
  ([Orbital_HF_IBS.C:24-38](../src/BasisSet/Imp/Orbital_HF_IBS.C)), which is what
  `SymmetryAdapted_IBS` and the `ChargeDensity` Fock build call. **No raw indexing in the SALC path.**
- **Builders** — `MakeDirect` / `MakeExchange` return an `ERI4` by value
  ([src/BasisSet/Atom/IrrepBasisSet.C](../src/BasisSet/Atom/IrrepBasisSet.C) ~277-365; molecular MnD path
  builds the analogue).
- **Dimensions** — `Nab()`, `Ncd()`, `size()`.
- **Debug/test only** — `operator()(a,b)(c,d)` raw read (the `MnDBench` sink at
  [MnDBench.C:56](../src/BasisSet/Molecule/bench/MnDBench.C), plus `operator==` / `fnorm` /
  `relative_fnorm` / `Transpose` in [ERI4.C](../src/BasisSet/Internal/Imp/ERI4.C)). `Transpose()` has **no
  production caller**.

So a replacement type must own: **build**, **contract** (`MatMul`, and the new fused op below),
**dimensions**, and a **debug element read**. That's a small surface — the blast radius is the container +
its contraction + two builders + the cache key + the Fock driver loop, not hundreds of indexing sites.

The J/K cache: `IntegralsCache_RAM<T>::Get(I4C, a, b, make)`
([DB_Cache_RAM.C](../src/BasisSet/Internal/Imp/DB_Cache_RAM.C) ~317-338) keys `map4_t = map<IBS_ID, map<IBS_ID, ERI4>>`
by `(a->BasisSetID(), b->BasisSetID())` with **no canonicalization** — `(A,B)` and `(B,A)` are distinct
entries today.

The Fock driver: `ChargeDensity` composites iterate their per-cd-irrep densities and each accumulates into
the full `Sab` via `bs_ab->AccumulateDirect(Sab, D_cd, bs_cd)`
([IrrepCD.C:58-69](../src/ChargeDensity/Internal/Imp/IrrepCD.C),
[CompositeCD.C:31-38](../src/ChargeDensity/Imp/CompositeCD.C)). The `(i,j)` and `(j,i)` blocks are
accumulated **independently** — the bra–ket symmetry is not exploited here.

## 4. The key idea: fuse the two contractions (no scattered penalty)

Do **not** store-canonical-then-flip. Instead, restructure the Fock build so one pass over a *single*
canonical block `J(i,j)` scatters into **both** Fock sub-blocks. Because `MatMul` already loads the
contiguous inner block `Jab = J(i,j)(a,b)`, the second (transposed) contribution is a **scaled whole-block
add**, never a transposed gather:

```
// contributes the (i,j) canonical block to BOTH S_i and S_j
for ab in irrep i (ia≤ib, w = (ia==ib?1:2)):
    const smat& Jab = J(a,b);
    S_i(ab) += sum( Jab % D_j );        // localized  → scalar into S_i(ia,ib)
    S_j     += (w * D_i(ia,ib)) * Jab;  // localized  scaled-block accumulate into S_j
```

Cost accounting per distinct pair `(i,j)`:

|                | RAM (blocks) | integral build | contraction / SCF iter |
|----------------|:---:|:---:|:---:|
| today          | 2   | 2   | 2 localized passes |
| fused canonical| **1** | **1** | **1 fused pass ≈ 2 localized** |

So the contraction is **cost-neutral**, while RAM and integral-build **halve**. The 4× scatter never
happens because we accumulate the whole block scaled by a scalar instead of gathering its transpose. This
is the "flat packed array walked once, scattering into all symmetry targets" pattern. `i==j` blocks stay a
plain `MatMul` (self-transpose, no partner).

## 5. Design

### 5.1 Custom `ERI4` (interface-first, storage-second)
**Design invariant: ALL symmetry/index bookkeeping lives inside `ERI4`.** The contraction loops, the
`ScatterBoth` weights, canonicalization/transpose orientation, and any packed-storage indexing are private
to `ERI4`; the Fock driver and consumers only ever say "contract yourself" / "scatter yourself into these
two blocks" and never see an index. (The integral *fill* in `MakeDirect`/`MakeExchange` stays with the
builders — that's physics + Rk-cache loops, not symmetry packing — but it writes through a controlled
`ERI4` seam, not by reaching into `smat<smat>`.)

Replace the `smat<smat>` type with a class whose **public** surface is operations, not storage:
- `MatMul(S_ab, D_cd)` — the localized contraction (diagonal blocks and general use).
- `ScatterBoth(S_i, S_j, D_i, D_j, w?)` — the fused §4 op (off-diagonal canonical pairs). Owns the
  symmetry bookkeeping.
- `Nab()`, `Ncd()`, `size()`.
- A **debug** element read (`at(a,b,c,d)` or keep `operator()`), used only by tests/bench and `==`/`fnorm`.
- `Transpose()` may be dropped (no production caller) or kept as a test helper.

Internal **storage is a separable decision** (do 5.1 with the current `smat<smat>` inside first, then tune):
- Option A: keep `smat<smat>` internally — minimal change, still gets the §4 win.
- Option B: flat packed array over the unique `(a≤b, c≤d)` quadruples with explicit symmetry weights —
  best locality for `ScatterBoth`, and the natural home for eventually packing the *diagonal* `i==j`
  block's internal bra–ket symmetry too. Defer unless profiling asks for it.

### 5.2 Canonical J/K cache key
In `Get(I4C, a, b, make)`: canonicalize to `(min,max)` by `BasisSetID` (or an explicit irrep order),
store once, and return `{block, needsTransposeBit}` so the driver knows the orientation. Only ever one
entry per unordered pair.

### 5.3 Fock driver restructure  — RESCOPED (survey 2026-07-02, see §5.4)
Original plan: change the `ChargeDensity` accumulation from "each cd-density accumulates into the full
`Sab` independently" to "iterate **canonical** irrep pairs `i≤j`, fetch the one block, call `ScatterBoth`
into the `i` and `j` slices of the **full Fock matrix**" (`i==j` → `MatMul`).

**That mental model was wrong about where the driver lives.** A read-only survey (§5.4) showed the Fock is
**not** a single "full matrix" the `ChargeDensity` layer slices — it is built **one irrep block at a time,
independently**, and no layer at or below `ChargeDensity` ever holds two irrep Fock targets at once. So the
restructure is not a `ChargeDensity` edit; it needs a new **cross-irrep seam** threaded from `CompositeWF`
down through the Hamiltonian terms (§5.4). This is the whole reason Stage 3 was paused for a rescope.

### 5.4 The real driver, and the cross-irrep context seam (survey finding)

**Call chain (atomic per-`l`-irrep ERI4 path):**
```
SCFIterator
  → CompositeWF::DoSCFIteration            [OUTER loop over irreps — CompositeWF.C:74]
      for each IrrepWF:
        IrrepWF::CalculateH                 [IrrepWF.C:41]
          → Hamiltonian::GetMatrix(bs = ONE irrep, S, cd)     [HamiltonianImp.C:48]
              Σ terms:  term->GetMatrix(bs, S[, cd])
                Vee::CalcMatrix(bs = ONE irrep, cd)            [Vee.C:27]
                  → dm->AccumulateDirect(single Jab target, hf_bs)
                      → CompositeCD iterates cd-irreps → IrrepCD
                          → Orbital_HF_IBS::AccumulateDirect
                              → MatMul(Jab, Direct(i,j), D_j)
```
The `Jab` reaching the `ChargeDensity` accumulation is **one** irrep's block; the composite loops the `cd`
irreps *into that single target*. `Jab(i)` and `Jab(j)` are produced in **two different** `IrrepWF`
`CalculateH` calls and never coexist. `ScatterBoth(Si,Sj,Di,Dj)` needs all four at once → **nothing below
`CompositeWF` can call it.**

**Scope facts:**
- Off-diagonal `i≠j` ERI4 pairs occur **only in the atomic path today** — the molecular MnD and SALC paths
  call `Get(I4C, raw, raw)` with the *same* whole-AO basis on both sides
  ([SymmetryAdapted_IBS.C:32-33](../src/BasisSet/Imp/SymmetryAdapted_IBS.C) passes `raw` twice), so those
  blocks are diagonal / self-transpose and canonicalization is a no-op for them. **But this is not a
  permanent property**: periodic HF for **solids** will carry **k-vector irreps _and_ SALC simultaneously**,
  giving genuine off-diagonal cross-irrep blocks. So the seam must be **general block algebra**, never
  atom-special-cased (mirrors the §9 "solids inherit it" note).
- **Bit-identity is unreachable** and should not be the gate: today `J(j,i)` is built by an *independent* Rk
  quadrature; any symmetry-exploiting version derives it from `J(i,j)ᵀ`, which differs at ~1e-13 (reordered
  sums). Pin a **tight regression anchor** (~1e-10…1e-13), not literal bit-identity. (The *`ScatterBoth`
  unit test* vs two `MatMul`s on the **same** block stays `EXPECT_EQ` — that arithmetic really is identical.)
- §5.2 (canonical key) and §5.3 are **coupled**: canonical-key *alone*, serving the flip via a transpose
  view, taxes every flipped block with the ~4× scattered contraction (§2) every SCF iteration. So Stage 3
  cannot ship the cache-key change without the fused scatter — they land together or not at all.

**The seam (design direction, OOD):** `CompositeWF` must **not** learn 4-index ERI details, but it *is* the
natural owner of the fact that a **global, cross-irrep view exists**. So:

1. `CompositeWF` builds a **context struct** — the list of irreps and their per-irrep data (bases,
   densities, Fock accumulators, and whatever else a term might need a global view of) — and **threads it
   through** the `Hamiltonian::GetMatrix` → `term->GetMatrix`/`CalcMatrix` call chain as an extra parameter.
2. **Every Hamiltonian _term_ may use or ignore it.** 1-e terms (Kinetic, Nuclear, Overlap) and the
   per-irrep DFT `Vxc` ignore it and stay exactly as they are. The **Coulomb/exchange (`Vee`) term** is the
   one consumer: given the global view it exploits the cross-irrep `J(i,j)=J(j,i)ᵀ` symmetry.
3. **Start as a plain struct; architect a polymorphic type later** (a pointer that *answers abstract
   questions* — "give me the canonical irrep pairs", "the Fock target / density for irrep `k`" — so atom vs
   molecule vs solid all plug in without the term branching) **only if it still makes sense** after the
   struct version exists.

**REALIZED (post-3b):** once the consumer was written it was clear the struct only ever held
`vector<const tobs_t*> irrepBases`, which is *exactly* what the whole (composite) `BasisSet` already is —
`Iterate<tobs_t>()` over it yields the per-irrep blocks. So `tHamiltonianContext` was **deleted** and the
seam threads `const tbs_t<T>*` (the composite basis, `= BasisSet::BasisSet<T>`, nullptr = "no view") through
`GetMatrix`.  `CompositeWF` just passes its own `itsBS`; `MakeContext` is gone.  The `BasisSet` *is* the
polymorphic "answers abstract questions" view the point above anticipated — no wrapper needed.

**Two candidate mechanisms for how `Vee` uses the view** (decide during implementation planning):
- **(a) Whole-system precompute, per-irrep pull unchanged (preferred — least invasive).** On its first call
  in an SCF iteration, `Vee`, handed the context, walks **canonical irrep pairs once**, and `ScatterBoth`s
  each fetched block into an internal per-irrep-keyed `J`/`K` map (`i==j` → `MatMul`). Each subsequent
  per-irrep `GetMatrix` just **slices its block out** of that map. `GetMatrix`'s one-matrix-per-irrep return
  contract is unchanged; the two-target scatter is hidden inside the term. This is the "flat packed array
  walked once, scattering into all targets" pattern, and it generalizes the existing
  [`SymFockCache`](../src/BasisSet/Imp/SymmetryAdapted_IBS.C) memoize-once idea.
- **(b) Term writes all targets directly.** The struct owns mutable per-irrep Fock accumulators; a new
  whole-system `term->Assemble(context)` scatters into them and the per-irrep `GetMatrix` skips the Coulomb
  term. More invasive (changes the assembly protocol + energy path); defer unless (a) proves awkward.

## 6. Generic sizing via a worst-case `Irrep` (the atomic LMax track — separate layer)

This is a *coordinated but distinct* track from §5, and it answers "is `LMax` generic enough for
molecules/solids?" — **no**, but `Irrep` is.

`LMax` (the R^k truncation over a shared exponent pool) is an **atomic-radial** concept; molecules/solids
compute whole shell-quartets and have no "grow the entry" analogue, so `LMax` must **not** climb into the
generic ERI4/basis layer. Instead the generic carrier is bare `sym_t` — the polymorphic **spatial** symmetry handle (the `sym`
field of `Irrep` in [src/Symmetry/Irrep.C](../src/Symmetry/Irrep.C); **no spin** — `LMax` is purely
spatial, so we don't want the `{Spin,sym}` bundle here):

- The basis answers a generic query — "your worst-case `sym_t` for ERI evaluation" — declaratively at
  registration (it knows all its shells up front). Atom/molecule/solid all speak `sym_t`.
- Each backend **projects** it: the atomic `Cache4` does `LMax = Getl(sym)`
  ([Symmetry/Spherical.C](../src/Symmetry/Spherical.C) exports `Getl(const sym_t&)`); the molecular/solid
  engines ignore it (nothing to size).

This retires the fragile part of the current atomic Rk cache (see [doc note] and
`project_atom_eri_caching_perf`): the incremental per-shell `maxls` discovery + the
eviction-on-`Register` + `Rk::isSupported`. Two clean options, both keyed off the declared worst-case
`Irrep`:
- **Declare-ceiling**: basis declares worst-case irrep → `Cache4` sizes every Rk to that `LMax` up front;
  cross-basis, track the max declared ceiling and grow.
- **Demand-grow**: the requester passes its `Irrep`(s); each `Rk` grows `Rabcd_k` monotonically to cover
  the requested `LMax`, never evicts. (Preferred — see the earlier assessment; `Register()` collapses to
  grouper-index assignment only.)

Keep §5 (generic, bra–ket, high value) and §6 (atomic, LMax→Irrep) as separate PRs at their own layers.

## 7. Correctness guards (fail loudly — this is the whole point)

- **Fock regression anchor** (NOT literal bit-identity — see §5.4): pin a converged HF energy AND the
  assembled Fock matrix before/after the rework (the A_HF_dfPin idiom + an M_* molecular case) to a **tight
  tolerance** (~1e-10…1e-13). The fused `ScatterBoth` reproduces the independent double-contraction only to
  ~1e-13 because today's `J(j,i)` is an independent Rk build, not `J(i,j)ᵀ`; demanding bit-identity would
  fail on reordered sums. (The *`ScatterBoth`-vs-two-`MatMul`s* unit test on the **same** block stays
  `EXPECT_EQ`.)
- **Symmetry invariant**: keep/extend the `fnorm(J(i,j), J(j,i)ᵀ)≈0` assertions for Direct and Exchange.
- **Dedup/reuse guard**: mirror `Cache4Tests.HF2_S{L,G}_CrossElementReuse` — assert that after touching
  pair `(i,j)` the cache holds ONE entry, and a request for `(j,i)` adds ZERO new blocks (regression here
  = someone re-materialized both orientations, silently doubling RAM).
- **`ScatterBoth` == two `MatMul`s** unit test: random `J`, `D_i`, `D_j` → the fused op equals the
  independent `MatMul(J,D_j)→S_i` + `MatMul(Jᵀ,D_i)→S_j`.

## 8. Staging

1. Introduce the custom `ERI4` with `smat<smat>` inside (Option A) + `MatMul` unchanged; drop the exposed
   `const smat&` return from the public API (route raw reads through a debug `at()`), fix bench/tests.
   Bit-identical, no behaviour change. **Gate: full suite green, Fock pinned.**
2. Add `ScatterBoth` + its unit test (§7). No driver change yet.
3. **RESCOPED (§5.4)** — the cross-irrep context seam + canonical cache key + fused scatter, as sub-stages:
   - **3a. DONE** (inert plumbing, 167 UTMain + 109 UTAtom_BS green, all pins bit-identical). Threads the
     whole (composite) basis `const tbs_t<T>*` via `IrrepWF::CalculateH` →
     `tHamiltonian::GetMatrix(bs,S,cd,wholeBasis)` → `t->GetMatrix(bs,S,cd,wholeBasis)` on the **dynamic**
     terms only (static terms keep the 2-arg form); `CompositeWF` passes its own `itsBS`, a consumer does
     `wholeBasis->Iterate<tobs_t>()`.  Least-churn seam: `tDynamic_HT::GetMatrix(…,wholeBasis)` is a NON-pure
     overload defaulting to the existing 3-arg (so every term ignores it for free — the ~11 `CalcMatrix`
     overrides are untouched); `tHamiltonian` keeps its 3-arg (null-basis) overload for the stand-alone test
     callers (PlaneWaveDFTUT, EigenSolverUT).  (Originally threaded a `tHamiltonianContext<T>` wrapper struct;
     deleted post-3b once it was clear the composite `BasisSet` already IS that view — see §5.4 REALIZED.)
   - **3b. DONE (Coulomb/J only)** — `Vee` consumes the context via mechanism (a).  New seam:
     `Orbital_HF_IBS::AccumulateDirectBoth(Ji,Jj,Di,Dj,cd)` = `Direct(*cd).ScatterBoth(...)` (fetches ONLY
     the canonical block); `tDM_CD::AccumulateDirectAll(Jall,abBases)` (whole-system) — `Polarized_CD` sums
     both spins, `tComposite_CD` runs the canonical-pair loop (diagonal → existing `AccumulateDirect`
     MatMul; off-diagonal → `IrrepCD::AccumulateDirectBoth`, which sibling-casts its partner as
     `MixIn`/`GetChangeFrom` do and keeps densities encapsulated).  `Vee` overrides BOTH `GetMatrix`
     overloads: the 4-arg (Fock, has ctx) stashes the run-stable irrep-basis list + builds the whole-system
     Coulomb once per density into `itsJ` (keyed by BasisSetID) + slices; the 3-arg (energy `DM_Contract`,
     no ctx) reuses the stashed list so the post-diagonalization density gets the SAME banked build (else
     the energy path would re-materialize `J(j,i)` every iteration).  **SALC caveat:** `SymmetryAdapted_IBS`
     has no per-irrep-pair ERI4 (empty `MakeDirect` → builds the whole-AO Fock and slices), so it OVERRIDES
     `AccumulateDirectBoth` to fall back to two independent AO slices (the AO build already banks the full
     8-fold symmetry).  **Gate:** 167 UTMain green; `Jac` RAM 521→369 MB (off-diagonal blocks now single-
     orientation; diagonal blocks are self-transpose and unhalvable, hence <50%).
   - **3b-K. DONE (Exchange/K)** — the mechanical mirror of 3b for `K(i,j)=K(j,i)ᵀ`: `AccumulateExchangeBoth`
     / `AccumulateExchangeAll` alongside the Direct versions.  Spin twist: exchange is SAME-SPIN, so the
     whole-system build runs per single-spin density.  `Vxc` (RHF) sums both spins (`Polarized_CD::
     AccumulateExchangeAll`) and stores blocks scaled −½; `VxcPol` (UHF) builds PER SPIN from
     `Polarized_CD::GetChargeDensity(s)` and stores scaled −1, keyed by (spin,irrep).  Both cache the scaled
     blocks so `GetMatrix` returns a reference; `SymmetryAdapted_IBS` gets the two-AO-slice fallback.
     **Gate:** 167 UTMain + 109 UTAtom_BS green; `Kab` RAM 521→369 MB (now symmetric with `Jac`).
   - **3c. DONE (rigid guard, not the transpose-serving key).** Instead of §5.2's "canonicalize + serve the
     transpose orientation" (which needs a `{block, needsTranspose}` return contract), we enforce
     **canonical-only** and make the partner a hard error: `IntegralsCache_RAM::Get(I4C,a,b,…)` **throws**
     (a `std::runtime_error`, live in Release — not a debug assert) when `a>b` by BasisSetID.  Production
     fetches are made canonical at the source: `Orbital_HF_IBS::AccumulateDirect/ExchangeBoth` request the
     min-BasisSetID block and swap the two `ScatterBoth` target/density pairs when the partner sorts first
     (identical result, `J(b,a)=J(a,b)ᵀ`).  So the bra-ket partner can never be built or stored — the
     invariant is structural, not a discipline every caller must remember.  Bonus: canonicalizing the fetch
     also de-duplicates *cross-basis* storage (two bases whose iteration orders differ now agree on the
     canonical key).  Integral-level unit tests that deliberately built both orientations
     (`A_Cache4`/`DBCache`/`BasisSet_Atom`) now fetch only the canonical block from the cache, verify
     symmetry against the **uncached** `MakeDirect`/`MakeExchange`, and `EXPECT_ANY_THROW` on the
     non-canonical request (turning the old symmetry check into a guard test).  **Gate:** 167 UTMain + 109
     UTAtom_BS green; zero guard throws in production.
4. (Optional) flat packed storage (Option B) if profiling wants it.
5. (Separate track) §6 atomic Rk `LMax`→`Irrep`: worst-case-irrep query + demand-grow Rk; delete
   eviction/`isSupported`; reuse guard stays green.

## 9. Risks / open questions

- **Driver restructure is the real cost — and it is NOT in `ChargeDensity` (§5.4).** The survey found the
  Fock is built one irrep at a time (outer loop in `CompositeWF::DoSCFIteration`), so the two-target scatter
  needs a cross-irrep **context struct** threaded from `CompositeWF` through the Hamiltonian terms, not a
  `ChargeDensity`/`CompositeCD` edit. The `ChargeDensity` layer never holds two irrep targets. Scope 3a
  (inert struct plumbing) as its own bit-identical step before touching the cache key.
- **`ScatterBoth` symmetry weights** (the `w=2` off-diagonal factor, `i==j` and `a==b`/`c==d` edges) are
  exactly the bookkeeping `smat<smat>` was hiding — unit-test them hard (§7).
- **Complex `T` / plane waves**: `IrrepCD<dcmplx>::AccumulateDirect` already asserts HF-not-applicable, so
  the rework is real-`T` only; keep the complex path asserting out.
- **Solids**: no HF-for-solids yet; design the driver so the canonical-pair scatter is not atom/molecule
  specific (it's pure block algebra), so a future periodic HF inherits it.
