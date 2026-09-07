# Lattice-gas / cluster-expansion for Li ordering — the design, before the code (cut 2026-09-07)

**Status: SPECCED, NOT BUILT — and deliberately deferred** (user, 2026-09-07: *"We don't need multiple
configurations right now"*).  This file exists so the design is not re-derived.  The one piece that IS
built is `Supercell` (`src/Structure/UnitCell.C`), because `doc/ParallelAndOraclePlan.md` 2.1 needed
replication for the size question; its declaration carries the forward-compatibility reasoning.

**The north star this serves**: Li/Na cathode voltage curves — total energies for many Li configurations,
fitted to a lattice-gas model, Monte-Carlo'd.

---

## 1. THE SYSTEM, AND THE THING THAT SURPRISED US

The target is λ-Mn₂O₄ (delithiated spinel) with Li placed on a subset of the 8a tetrahedral sites, giving
Li\_{n/8}Mn₂O₄.

★ **THE CONVENTIONAL CUBIC CELL IS ALREADY THE 8-SITE LATTICE — the first campaign needs NO supercell.**
Spinel AB₂O₄ is Fd-3m with Z=8: A on 8a (tetrahedral), B on 16d, O on 32e.  So

| | |
|---|---|
| LiMn₂O₄ conventional | Li₈Mn₁₆O₃₂ — 56 atoms |
| λ-Mn₂O₄ (framework) | **Mn₁₆O₃₂ — 48 atoms** |
| 8a sites in that cell | **8** |

⇒ \f$x=n/8\f$ for \f$n=0\ldots8\f$, \f$2^8=256\f$ decorations, ~20 inequivalent.  The user's "2×2×2 …
Li in one of 8 cells" intuition was right about the SITE COUNT; the conventional cell delivers it directly.
(⚠ Confirm against the reference structure actually used before building on this.)

**Where a supercell then earns its place** — not the first campaign, the second:
- finer \f$x\f$ resolution (1/16, 3/16, …);
- ★ **cluster-expansion RANGE convergence**: a CE fitted in one cell cannot see pair interactions longer
  than that cell, so a larger one is needed to test that the fitted model predicts held-out energies.

## 2. THE GUEST-LATTICE MODEL

A configuration is an occupation bit-vector over an ORDERED site list — the lattice-gas variable itself —
so that is what the API should hand back, not a pile of decorated cells:

```
host    : UnitCell             the replicated framework (Mn2O4), NO guests
sites   : std::vector<rvec3_t> candidate guest positions, fractional in the host cell
```

Decorating is then the existing `AddAtom(3, sites[i])` for each occupied `i`; Li₀.₁₂₅ is `10000000`.
⚠ **The site list is SUPPLIED, never hardcoded to a Wyckoff set.**  8a only for now (it supports ordering
at half filling); 16c octahedral is *"something we should have the ability to try"* (user) — its energetics
and packing feasibility are open, and it may not even be seen in neutron diffraction.
⚠ For a combined 8a+16c list the sublattices share faces and cannot be co-occupied at short range, so that
model needs an EXCLUSION constraint in the enumerator.  Another reason the site set is a parameter.

## 3. SYMMETRY DEDUP — THE ORDER-OF-MAGNITUDE WIN, AND WE OWN MOST OF THE MACHINERY

256 decorations reduce to ~20 inequivalent, and the orbit SIZE is the multiplicity the CE and the MC both
need.  This is `FoldPointsPeriodic` one level up: instead of folding POINTS under ops, fold OCCUPATION
VECTORS under the permutation each op induces on the site list — and that permutation is exactly what
`FoldPointsPeriodic`'s `apply` already computes.

⚠ **THE ONE GAP: `SpaceGroup::Detect` documents "Assumes a primitive cell (one τ coset per W)"**, so
pointing it at a supercell MISSES the internal pure translations that carry most of the dedup power.  The
fix is better than detection: we BUILD the supercell, so those translations are known by construction —
compose `Detect(primitive)` with the replica translations and hand the group over.  Same discipline as the
imposed-symmetry path, where ops are ctor-injected rather than rediscovered.

## 4. NEIGHBOUR SHELLS — A GENUINE GAP IN `qcStructure`

The MC app is *"a completely separate app … but it should use qcStructure library to store the lattice and
locate neighbour shells"* (user).  There is **no neighbour or coordination-shell facility** anywhere in
`src/Structure` or `src/Symmetry` (verified 2026-09-07).  `Ewald.C` enumerates images, but for a convergent
sum, not as an addressable shell structure.  So: new code, in `qcStructure`, serving two consumers.

★ **AND IT DOES NOT VIOLATE THE "NO CUT IN LATTICE SUMS" PIN — name the difference in the interface.**
That rule is about never truncating a CONVERGENT PHYSICAL SERIES.  A cluster expansion is the opposite
case: it *is* a finite basis of cluster functions by construction, so a shell index or radius is a **MODEL
parameter**, not a numerical truncation.  Name it after the model (`clusterShells`, not `cutoff`) so the
two can never be confused by a later reader.

## 5. WHERE THINGS LIVE

`qcStructure` already depends on `qcSymmetry` (`SymmetrizeMesh.C` imports it), so the fold machinery is
reachable and all four pieces can live there:

| piece | status |
|---|---|
| `Supercell(prim, n)` | ✅ **BUILT** 2026-09-07 — free function, no interface change |
| guest lattice (host + ordered site list) | specced above |
| configuration enumeration + orbit multiplicities | specced above; sits on `Fold` |
| neighbour shells | specced above; genuine gap |

Consumers: the SCF driver decorates → `SolidCalculation` → energy; the MC app is a new `CLIapps/` binary
reading the same lattice.  **Spin is NOT the helper's business** (user, 2026-09-07: *"right now I think
spin structure is specified in the seed charge density … leave that out of the supercell helper for now"*).

## 6. ⚠ THE SIZING QUESTION, UNANSWERED

A decorated 48-atom cell has LOW symmetry — our 12–48× folds largely vanish — and that is a different
regime from the 4-atom MnO at Γ with 24 ops that every Phase 1 number was measured on.  **Size ONE
configuration before committing to a campaign of 20–200.**  If a single run is hours, the design changes:
the framework's setup (basis, Becke mesh, Φ tables ≈ 15 s/call on MnO) is identical across configurations
and would want sharing — the same "eager shared prologue" shape as `doc/OpenWork.md` item **KP**.
