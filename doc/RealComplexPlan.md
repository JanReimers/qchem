# Real vs complex — what decides the scalar type, and how it flows

**Status: PLAN, no code beyond Step 0.**  Written 2026-08-16 out of the GPW runtime work
(doc/GPWPlan1.md "THE RUNTIME GAP, MEASURED"), which surfaced the question but does not settle it.
The decision to write this instead of starting: the change is large, and several open OOD items sit
on the very faces it would restructure (see "Prerequisites").

---

## 1. What the physics dictates

A block can be real **iff its basis functions can be chosen real AND every Hamiltonian term
restricted to it is real-symmetric.**  A conjunction — the two halves fail for different reasons and
are answered by different objects.

| case | what forces complex | decided by | TR survives | structural side-effect |
|---|---|---|---|---|
| Bloch k ∉ TRIM | \f$\chi^k=\sum_R e^{ikR}\chi\f$ | **irrep** (`BlochQN`) | yes (pairs ±k) | none |
| spin–orbit \f$\xi\,L\cdot S\f$ | \f$L=-i\,r\times\nabla\f$ | **term** | yes → Kramers | **large**: spin channels merge into 2-spinors, `ms` stops being a good label, double groups |
| magnetic field / vector potential (A·p, GIAO, current-DFT) | \f$p\to p+A/c\f$ | **term** | **no** | none, but no Kramers to lean on |
| molecular 4-/2-component Dirac | \f$\alpha\cdot p\f$ coupling | **term** | yes | large/small structure |
| ATOMIC Dirac (our `A_DHF`) | — stays REAL | convention | — | the \f$i\f$ absorbs into the small component; works by spherical symmetry |
| electric field | — stays REAL | — | yes | none.  "External fields" is NOT uniform: E real, B complex |
| complex harmonics / complex point-group irreps (E of C₃, C₄…) | — convention only | — | yes | real combinations always exist for a real H; the SALC machinery already takes them |
| non-collinear spin / spin spirals / GHF | rotating quantization axis ⇒ spinor | **ansatz** | depends | 2-component |
| complex SCF instability (real H, complex minimum) | the nonlinear functional, not the operator | **ansatz policy** | yes | none |
| CAP / resonances | non-Hermitian complex-SYMMETRIC H | term | no | different linear algebra entirely — not `LASolver`'s problem |

Two things worth stating because they are easy to get wrong:

- **"Some open-shell systems need complex" is false for the linear problem.**  A real symmetric H
  always admits real eigenvectors.  Complex can only help the NONLINEAR problem: a real HF/KS
  stationary point unstable toward a lower complex one (the classic complex RHF/UHF instability).
  That is an explicit ansatz policy, and it can only ever DOWNGRADE a block real→complex.
- **SOC is not "T becomes complex".**  It merges the spin channels into 2-component spinors, so the
  `ms` half of `Irrep` stops being a good quantum number.  A scalar-type rule cannot express that;
  it is a block-STRUCTURE change.  Anything designed here must not hard-wire "2 independent spin
  channels" either.

### The rule

```
block is real  ⇔  irrep.IsReal() ∧ (every term.PreservesReal())    [∧ not overridden by ansatz policy]
```

Each participant answers about ITSELF, which is what makes it derivable rather than declared:
`BlochQN::IsReal()` from \f$2\,ik\equiv0\ (\mathrm{mod}\ N)\f$ — **exact integer arithmetic, no
tolerance**; `Spherical` / molecular point-group irreps trivially true; an SO or vector-potential
term always false; kinetic / Coulomb / XC / local-PP / KB always true.  Adding SOC to a run then
flips every block complex automatically, with no special-casing anywhere.

## 2. TWO types, not one — the distinction that must not be lost

```
basis type    = irrep.IsReal()                          // "are my functions (and S, T, V) real?"
working type  = basis.IsReal() ∧ ham.PreservesReal()    // "are H, C, D real?"        working <= basis
```

With SOC the basis functions are still real — real Gaussians, real SALCs — so S, T, V_loc, the KB
projectors and **the collocation streams** stay real; only H, C, D go complex.  Since the geometry
layer is what dominates GPW, conflating the two types would needlessly double the cost of the
expensive half the moment SOC arrives.

**Today they are welded:** `tSCFIterator<T>` takes `tbs_t<T>*` — one T for basis and everything
downstream.  For the Bloch-only case that is harmless (basis real ⇒ working real), so the first
increment need not decouple them.  **The standing instruction is only: do not harden it further.**
When the `dcmplx`-pinned solid classes are touched, keep basis type and working type as separate
parameters even while they are always equal.

## 3. Where the decision is computed

**Not in the basis, and it cannot be:** `qcHamiltonian` links `qcBasisSet`, and every driver builds
basis → Hamiltonian → iterator (`Hamiltonian::Factory(..., bs.get(), ...)` takes the basis).  A
basis that consulted H would be a dependency cycle the linker rejects.

- basis / irrep answers about itself → `Symmetry::IsReal()` (virtual, default true; `BlochQN`
  overrides), forwarded by `IrrepBasisSet::IsReal()`.
- each term answers about itself → `PreservesReal()`, constant per term type; `tHamiltonian<T>`
  folds its terms with ∧ (a one-liner).
- **the composition root computes the AND** — the layer that already constructs both and then builds
  the WF, CD and accelerator (`RunGpw` / `Hamiltonian::Factory` / `qchem::Calculation`).  It is the
  only place with both facts, and it is already the assembly point.

## 4. The containers that inherit it

One decision, carried by the type system, in dependency order:

1. `BasisSetImp<T>` → its `IrrepBasisSet` children — **the decision point**
2. Hamiltonian `tHT_Common<T>::CacheMap` (`std::map<Irrep,hmat_t<T>>`) → per-irrep H blocks
3. `CompositeWF` → `tIrrepWF` children (`itsIWFs`/`itsQNWFs`) — where the real eigensolve lives
4. `tComposite_CD` → `tDM_CD` children (`itsCDs`)
5. **`tSCFAccelerator<T>` — the open one.**  See §6.

**Child slot** = the element type of the container: today `tDM_CD<T>*` / `tIrrepWF<T>*` (one T per
run, so all children necessarily agree); under this plan
`std::variant<IrrepCD<double>, IrrepCD<dcmplx>>` or equivalent, so children of one composite may
differ.  `std::visit` with a GENERIC LAMBDA keeps every algorithm single-source — that, not the
container, is the answer to the dual-maintenance concern.

**Keep the composites templated at their FACE.**  A genuinely non-template `tComposite_CD` would
still have to derive from `tDM_CD<dcmplx>` to serve `tHamiltonian<dcmplx>`, so molecules would need
a SECOND composite class — two copies of the aggregation logic, i.e. exactly the dual maintenance
being avoided.  Face stays `tComposite_CD<T>`; only the child slot varies.  Molecules then hold a
variant that never takes the complex alternative: a negligible tag, one implementation.

## 5. Staging

- **Step 0 — DONE `1b8b9a83`.**  Exact ±1 Bloch phases at TRIM, so realness is a BITWISE fact, not a
  tolerance.  This is what lets a later narrowing of H be an ASSERTED narrowing.  Gate
  `GPW.TRIM_BlochMatricesAreExactlyReal`.
- **Step 1 — `Symmetry::IsReal()` + `IrrepBasisSet::IsReal()` + `Term::PreservesReal()`.**  Pure
  queries, no behaviour change, no types moved.  Report the answer per block in the run report
  (nobody should have to infer it).  `IsTRIM` in the GPW evaluator should then DEFER to the irrep
  rather than re-deriving realness from a float k.
- **Step 2 — the composite child slot** (`tComposite_CD` first: fewer methods, ~8 of 18 are T-typed
  and all are aggregation points; `IrrepCD<T>` is already templated so both alternatives exist).
  Buys the WF/D half — diagonalization and C/D storage.  Dominant in the many-PW/ultrasoft regime;
  small in GPW.
- **Step 3 — un-pin the basis** (`GPW_IBS` from `IrrepBasisSet<dcmplx>`), so terms produce real H
  directly and Φ + the quadrature GEMMs go real.  Buys the Hamiltonian/quadrature half — dominant in
  the GPW regime.  The composite variant then falls out as a CONSEQUENCE of blocks having different
  basis types, which is the coherent order.
- **Step 4 — the accelerator** (§6), and only then a mixed-mesh run as the acceptance test.

An all-TRIM mesh (Γ-only, or Γ-centred 2×2×2 = `MNO_KMESH=2`) makes an entire run real with no
heterogeneity at all.  Mixed meshes (3×3×3, shifted MP) are the only case that truly needs per-block
variation — worth remembering when scoping.

## 6. Open questions

- **The accelerator.**  `tSCFAccelerator<T>` is run-level.  In GPW (F is 122²) widening a real block
  at the face costs nothing.  In the many-PW regime DIIS history is over the COEFFICIENT arrays
  (N_pw × n_bands), which are exactly the block-typed objects one must not widen — the history is
  the memory hog.  So a mixed-irrep design eventually forces the accelerator to be block-aware, or
  to hold its history per block.  **This is the piece that could quietly negate the memory win in
  the regime that motivates the whole exercise; settle it before Step 2 ships.**
- **2-component structure.**  Non-collinear spin (likely to arrive before SOC-for-anisotropy, given
  the frustrated-magnet/battery direction) needs spinor blocks, not just complex scalars.  Nothing
  here should hard-wire "exactly two independent spin channels".
- **Ansatz policy.**  Where does "allow complex orbitals" (complex instability, GHF) live?  It is a
  run-level policy that downgrades blocks; probably `SCFParams`, but it interacts with the
  composition root's AND and has no home yet.

## 7. Prerequisites (doc/CleanupCandidates.md)

Load-bearing, because they sit on the faces this restructures — finishing them SHRINKS this work,
and none of them is invalidated by it (the type change alters child containers, not the method lists
the ISP splits produce):

- **V1.8** `IrrepCD`↔`IrrepCD` concrete same-class casts in the hot path — they target the very class
  that becomes a variant alternative.
- **V1.6** `tDM_CD::Accumulate*` face split — the composite's T-typed methods ARE the widening points.
- **V1.7** the periodic trio's 9 asserting overrides — each multiplies against a new alternative.
- **V1.11 ✅ DONE 2026-08-17 (`43bbebad`..`2398dd07`, five increments)** — the occupation seam is the
  SCFIterator's `OccupationPolicy` slot; the WF child carries NO occupation state, so the variant-child /
  narrowing work lands on a clean face.  **ALL §7 prerequisites are now DONE.**
- **V1.1 ✅ DONE 2026-08-16 (`d49db261`+`4b37221d`+`2bfb83b4`)** — the basis-face merge landed:
  `Projector3<T>` unifies the 3C tensor, `Band_FT_IBS` is deleted, and `Orbital_DFT_IBS<double,dcmplx>`
  (a real TRIM block on the run's complex fit basis — exactly this plan's case) is a live spelling whose
  tensor caches under `theCache<TFit>()`.  **V1.5 ✅ DONE 2026-08-16 (`f18a6ee9`+`9ebaebdb`)** — the
  grid-engine union is four client-named faces (`G_FieldEvaluator`/`G_Quadrature`/`G_StructureFactor`/
  `G_SpectralFilter`), so re-typing touches one face per consumer.  **V1.10 ✅.**  The basis was the
  decision point and the head of the type flow — **all basis-side prerequisites are now DONE; only
  V1.11 (the occupation seam) remains.**
- **R2.16** construction-time facts re-asked at run time — the same principle this rests on.

NOT prerequisites (do not gate on them): R1.0b (shared-radial reader), R2.20 (oracle helpers in a
test module), R2.5's remaining `exit(-1)` sites.
