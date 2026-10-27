# The k axis: block threading (KP), one q per star (Q9), the little group of q (D5)

**Status: DESIGN, 2026-10-27 — for the user's review; no code yet.**  Born at the top of `doc/` because it is about to be
executed.  It merges three open rows that all act on the SAME loop — the per-(k, σ)-block loop of an SCF iteration and of a
response kernel application: OpenWork §2 "k-point parallelism (KP)", "Response on the little group of q: revisit D5", and the
q-star half of "LR efficiency ledger" (`Records/LinearResponsePlan.md` Q9, D5).

## 1. What the time is, measured (no code)

**The response is two per-block loops.**  NiO k222 q222 (Oct 4, n = 128, `~/Code/qchem6-runs/NiO/nio_qmesh222_step5.log`, 8109 s):

| bucket | calls | s | share |
|---|---|---|---|
| response: integrate transition potential (pair GATHER) | 4256 | 3537 | 44 % |
| response: collocate transition density (pair SCATTER) | 4256 | 3046 | 38 % |
| response: XC-mesh transition adjoint + forward | 2×4256 | 1031 | 13 % |
| response: transition V_H (both halves) | 2×4256 | 258 | 3 % |
| everything else (SCF, setup) | | ~240 | 3 % |

4256 = 266 kernel applications × **16 blocks** (8 k × 2 σ).  Every bucket above is **per block**, and the blocks are visited
one after another.  Within a block the scatter/gather pair loops are OpenMP-threaded (`PairThreads`), but only when
`QCHEM_OPENMP_THREADS` is set — the 2026-10-26 overnight NiO run did not set it and runs on ~1.5 cores (~20 h).

**Where the redundancy is.**  (a) Symmetry-equivalent q are recomputed in full: tonight NiO (0,½,0) reproduced (0,0,½) to 7
digits after ~2.4 h serial.  (b) Inside one q, the k sum runs over the full mesh although the little group of q (that also
fixes the perturbed site) maps most k onto each other.  (c) The 16 blocks are independent until their δρ is summed.

**The SCF has the same shape** — `CompositeWF` (`src/WaveFunction/Internal/Imp/CompositeWF.C`): Fock assembly (`CalcH` per irrep),
the diagonalization, the GDM step and the density insert are plain serial `for (auto& w:itsIWFs)` loops.  The 2026-09-07
"eager refresh" hoisted every k-independent density-derived state out of the loop precisely so the blocks could run
concurrently; the threading itself was deferred because the benchmarks were Γ (width 2), with the note *"in the real world
multi-k is the norm"* (`Records/ParallelAndOraclePlan.md` §1.3a(iii)).

## 2. The three levers and what each saves

| lever | saves | on NiO k222 q222 (estimate) | needs |
|---|---|---|---|
| **KP** thread the block loop | wall time ÷ min(blocks, cores) | 16 blocks on 8 cores: up to ~8× over serial; vs today's inner pair threading, measure first (P0) | thread-safe const paths; a deterministic reduction |
| **Q9** one q per star | q count | 8 → 4 (Γ, X-type ×3 → 1, the (0,½,½)-type ×3 → 1, L) — ~2× | the op mapping q → q′ and site I → I′ (+ Bloch phase) |
| **D5** k on the little group of (q, I) | blocks per kernel application (time); with D5-b also ground-state RAM | ×|G_{q,I}|, several-fold on rutile/spinel (row); on AFM NiO small | orbital ROTATION (D_AO), phased symmetrizers |

They multiply.  KP is independent of symmetry and pays on every multi-k run (SCF too), so it goes first.

## 3. Design

### 3.1 KP — the block loop, threaded, bit-identical

- **Where**: the response kernel application (the AO transition density's pair loop, `src/ChargeDensity/Internal/TransitionDensity.C`,
  driving `CollocateTransition` / `IntegrateTransition`, `PG_Cart_MnD/Evaluator.C` ~2600–2760) and the SCF's `CompositeWF` loops.
- **Determinism (pin the anchors)**: GATHER is per block → per-block output, no reduction, threads freely.  SCATTER sums blocks
  into one δρ: each block scatters into ITS OWN buffer, then the buffers are summed **in block order**, serially.  Same bits at
  any thread count — so the sweep must show 0 moved, not "within tolerance".
- **One level of threads at a time.**  OpenMP's default (`max-active-levels` 1) already makes an inner `#pragma omp parallel`
  run single-threaded inside an outer region.  Rule: outer (blocks) when the width ≥ the thread count, else inner (pairs) as today.
  Γ runs (width 2) keep today's behaviour exactly.
- **Memory is a run choice** (user 2026-10-22): per-thread scatter buffers cost (threads × the multigrid δρ); print it on the
  `threads:` banner line, cap the outer width if it would exceed a stated budget.
- **Audit first**: every `mutable` memo on the const scatter/gather/XC-sampling paths (Cache2 "nested cached gets", the XC
  sampler's density-serial caches, the f_xc kernel version).  The eager refresh covered the SCF's k-independent state; the
  response's `RefreshForDensity(D0, δ)` and its δ-version cache are the response's equivalent and must be hoisted the same way.

### 3.2 Q9 — one q per star

For an op g = {W|τ, σ} of the ground state's group with W q ≡ q′ (mod G): χ_{I′J′}(q′) = e^{iφ} χ_{IJ}(q) with I′ = g I (+ lattice
R_I, which supplies the Bloch phase φ = q′·(R_{J′} − R_{I′})).  With `HubbardResolution::Manifold` (one scalar channel per
manifold, hp.x's) that is a **site permutation plus phases — no orbital rotation**.  Slot/Eigen resolution needs the manifold's
rotation (D^l, §3.3), so Q9 serves Manifold first and refuses the others until D^l exists.  The real-space χ_IJ(R) assembly
(`ToRealSpace`) then runs over the full q-mesh from the star representatives.  Time reversal (q′ = −q) is a complex conjugate.

### 3.3 The one new primitive: D_AO(g, k) — an op acting on the orbital basis

Both halves of D5 (and CK-3a) need the representation of g on the Bloch orbital basis, mapping block k to block gk:
**C_{gk} = D_AO(g,k) C_k** — an atom permutation (A → gA + R), the Bloch phase e^{i gk·R}, and per shell the rotation of its
angular part.  On our **Cartesian** shells that rotation is exact polynomial substitution x_i → W_ij x_j (no Wigner recursion);
the spherical view composes it with its fixed `Tc`.  σ swaps the spin blocks.  **Unit gates that need no SCF to be right**:
D_AO is unitary in the S metric (D† S_{gk} D = S_k) on any cell; on an IMPOSED converged state D† F_{gk} D = F_k to the SCF
tolerance.  This primitive is small, testable alone, and is where the risk concentrates — so it is its own phase.

### 3.4 D5-a — reduced k per (q, I) on a FULL-mesh ground state (time, not RAM)

1. **The group**: G_{q,I} = { g in the run's group : W q ≡ q (mod G), g I = I + R } (hp.x: the small group of q that does not
   move the perturbed Hubbard atom), plus time reversal combined with an op taking q → −q where the run has it.  On AFM NiO the
   sublattice-swapping ops carry σ and move I, so they drop out — expect a small group there; rutile/spinel are where it pays.
2. **The blocks**: K_{q,I} = the IBZ of the full mesh under G_{q,I} (`ReduceToIBZ` with that op subset), weights |star|/N.  The
   Krylov vector lives on K_{q,I} only; the partners k+q are stored (full mesh), so `Reference::Partners` needs no change.
3. **The kernel**: scatter the K blocks, then apply the PHASED symmetrizer P_q = (1/|G|) Σ_g e^{iφ_g} ĝ to δρ_q before it is
   consumed: on the raster levels (the existing τ-aware raster star-average, plus the q phase), on the XC quadrature (the
   folded Becke mesh, plus the phase; σ swaps δρ_↑/δρ_↓), and on the Hartree G-ball (`ΔGq_Map`: (G+q) → W(G+q), phase
   e^{−i(G+q)·τ}).  δρ_q is NOT lattice-periodic, hence the phases; `TransitionDensity.C` today says "never symmetrized"
   (§3d finding 5: a perturbation breaks the group) — that stays true of the FULL group; P_q uses only the subgroup that
   leaves this perturbation invariant.
4. **The Hubbard side** (the +U kernel and the measured χ rows): δn_J needs the full-BZ sum; rotate δD_k onto its star with
   D_AO (n × n matrix products — cheap next to a scatter) and contract as today.  Measured rows I′ ≠ I come out averaged over
   their G_{q,I} orbit, which is exact.
5. **Gate**: TiO₂ k2/q2 χ(q) reduced vs full mesh (both banked, row D5) to ~1e-8 relative; NiO imposed vs tonight's free run.

### 3.5 D5-b — an IBZ-reduced ground state (RAM)

Store the ground state on the IBZ of the full group (an imposed run already can); reconstruct each needed block k (and k+q) on
demand as D_AO(g,k_rep) C_{k_rep}, eigenvalues and occupations copied, cached per q.  RAM ÷ |G| for the ground state; the
response's per-q partner cache is the stated RAM-for-speed run choice.  The same reconstruction IS CK-3a's unfolding (restart a
free full-mesh run from an imposed saved state), so CK-3a lands here for free.  Gate: D5-b χ(q) == D5-a χ(q) bit-for-bit modulo
the rotation's rounding (1e-12).

## 4. Phases, in order, each with its gate

| phase | what | gate |
|---|---|---|
| **P0** | measure, no code: one NiO k222 kernel application serial / inner-8 / (after P1) outer-8; per-thread buffer size | numbers in Benchmark §10 |
| **P1** | KP in the response kernel application (audit, hoist, per-block buffers, ordered sum) | `GPW_Si.k211_*Chi*` χ BIT-identical at 1 and 8 threads; NiO P0 speedup |
| **P2** | KP in `CompositeWF` (Fock assembly, diagonalize, GDM step, density insert) | short + long sweep: 0 moved at 1 and 8 threads |
| **P3** | Q9 (Manifold resolution) | tonight's NiO per-q χ from 4 of 8 q; the Si supercell-equivalence gate |
| **P4** | D_AO(g, k) primitive | S-unitarity on TiO₂/NiO cells; F-invariance on an imposed converged state |
| **P5** | D5-a | TiO₂ k2/q2 reduced == full; NiO imposed == free |
| **P6** | D5-b (+ CK-3a) | D5-b == D5-a; CK-3a: a free run restarted from an imposed state starts at the answer |

## 5. Questions for the user (each changes the design)

1. **Assertion policy for Q9/D5 on a FREE run.**  Reducing by a group asserts the ground state has it.  On an imposed run that
   holds by construction.  On a free run (tonight's NiO) the options are: (a) never reduce — impose the group to get the
   speed; (b) an explicit deck opt-in (`response.assumeSymmetry`), audited by computing TWO members of one q-star and comparing
   (the §5b measurement that found 1e-5, made routine).  Recommendation: (b) — it keeps "imposition is opt-in" and makes the
   audit automatic.
2. **Threading level rule** (§3.1): outer-or-inner per loop (recommended: simple, no nested OMP), or a split budget
   (outer × inner ≈ cores, needs `max-active-levels` 2 and is harder to keep deterministic)?
3. **Order.**  P1–P2 (KP) first because every multi-k run gains and it is independent of symmetry; Q9 before D5 because it is
   small and tonight's NiO is its gate.  Or is the rutile/spinel RAM (D5-b) the priority?
