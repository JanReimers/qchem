# ERI4 Rework: a custom ERI container that owns its contraction and bra–ket symmetry

Status: PLAN (survey done, not started). Author: design session 2026-07-02b.

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

### 5.3 Fock driver restructure
Change the `ChargeDensity` accumulation from "each cd-density accumulates into the full `Sab`
independently" to "iterate **canonical** irrep pairs `i≤j`, fetch the one block, call `ScatterBoth` into
the `i` and `j` **slices** of the full Fock matrix" (`i==j` → `MatMul`). This is the most invasive part:
the driver must hand the op **two** target sub-blocks at once, which the current one-target
`AccumulateDirect(Sab, D_cd, bs_cd)` signature does not express.

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

- **Fock bit-identical regression**: pin a converged HF energy AND the assembled Fock matrix before/after
  the rework (the A_HF_dfPin idiom + an M_* molecular case). The fused `ScatterBoth` must reproduce the
  independent double-contraction to ~1e-13.
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
3. Canonical J/K cache key (§5.2) + Fock driver restructure (§5.3) to use `ScatterBoth`. **Gate:** RAM
   report shows one block per unordered pair; Fock bit-identical; reuse guard green.
4. (Optional) flat packed storage (Option B) if profiling wants it.
5. (Separate track) §6 atomic Rk `LMax`→`Irrep`: worst-case-irrep query + demand-grow Rk; delete
   eviction/`isSupported`; reuse guard stays green.

## 9. Risks / open questions

- **Driver restructure (§5.3) is the real cost.** The current one-target `AccumulateDirect` must become a
  two-target scatter over canonical pairs; the `ChargeDensity`/`CompositeCD` iteration owns this. Scope it
  before committing to the cache-key change.
- **`ScatterBoth` symmetry weights** (the `w=2` off-diagonal factor, `i==j` and `a==b`/`c==d` edges) are
  exactly the bookkeeping `smat<smat>` was hiding — unit-test them hard (§7).
- **Complex `T` / plane waves**: `IrrepCD<dcmplx>::AccumulateDirect` already asserts HF-not-applicable, so
  the rework is real-`T` only; keep the complex path asserting out.
- **Solids**: no HF-for-solids yet; design the driver so the canonical-pair scatter is not atom/molecule
  specific (it's pure block algebra), so a future periodic HF inherits it.
