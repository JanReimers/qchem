# The k axis: block threading (KP), one q per star (Q9), the little group of q (D5)

**Status: DESIGN RULED 2026-10-27** (user: §5 recommendations agreed; §6 interface design "looks very nice").  No code yet.  Born at the top of `doc/` because it is about to be
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

## 4a. Execution log

**P0 (2026-10-27) found the build, not the loop.**  qchem1's `build/Release` had NO OpenMP: the /opt clang has no `omp.h`, so
`find_package(OpenMP)` failed quietly and every `#pragma omp` compiled serial — `QCHEM_OPENMP_THREADS` was a no-op (the "5.4×"
numbers in `Records/ParallelAndOraclePlan.md` were measured in `~/Code/qchem6/build/Release`, configured with apt's libomp).
Reconfigured with the flags in the root `CMakeLists.txt` comment; CLAUDE.md now says how to check (`ldd rundeck | grep omp`).

**P1 (2026-10-27): the response kernel application threads over its block pairs.**
- `qchem.Parallel`: `OuterFor(n, f)` (items concurrently when `OuterThreads(n)` > 1, each into its own slot; first worker
  exception rethrown), `InnerThreads()` (1 inside an outer item — every inner region now asks it), `OuterLoopScope`
  (thread_local marker), `LazyLock`, and the TEST-ONLY `ScopedWorkerThreads`.
- Threaded: the transition density's two pair sums (`GetTransitionRepulsion`, `ProjectOnto`) and `InducedFock`'s coupled-pair
  gather — slots summed in the serial order.
- Lazy caches a worker can fill, now `LazyLock`ed: `GPW_Evaluator::Ladder`, PG_Cart `Shells()` / `BoxTasks()`, the pair
  sampler's `itsAdj`.  The concurrency contract is written on `tResponse_HT::GetMatrix(…, δ)`.
- `report::Timed`: per-thread nesting stack, mutex on the ledger.  Under an outer loop a per-block bucket is THREAD-seconds
  (summed over workers) and the master's enclosing bucket charges the loop's wall.
- **The bit-identity bound**: identical to serial when the thread count ≤ the NARROWEST outer loop's width.  The polarized XC
  refresh projects each spin channel separately, so on a polarized run that width is the per-spin pair count; above it a loop
  hands its threads to the inner regions (§5.2), which sum in arrival order (measured 2e-15, run to run).  Gate:
  `ResponseKernel.GPW_Si_k311_q13_KernelIsBitIdenticalAcrossThreads_{Raster_Pol,Becke_UnPol}` (== on every element of δF).
- **For P2, found during the P1 audit**: `SinglesDensitySampler::MatrixT` (the SCF Fock's XC matrix) REFITS on a change of v and
  caches `itsFittedV` — the two spins' blocks ask with different v, so concurrent Fock blocks would race on the fit.  Hoist the
  per-spin fit into the eager refresh before P2 threads `CalcH`.

**P2 audit (2026-10-27) — what `CompositeWF`'s Fock loop WRITES, so it cannot run concurrently as is.**
- **The per-irrep Fock cache resets on first use.**  `tDynamic_HT_Imp::GetMatrix` clears its whole `itsCache` map and records
  the density version when it first sees a new `cd`.  The R1.0h pre-create phase (`PrepareSlots`) clears the map and makes
  every irrep's node, but does NOT record the version, so the first block of every iteration still wipes the map: concurrent
  blocks would be filling nodes another block is erasing.
- **Terms refit their shared state on first use** through `newCD(cd)` inside `MakeMatrix`: `Hubbard_U` (the occupation
  matrices, `Hubbard.C` ~1017), `Vee_Hartree` (`PWTerms_Hartree.C` ~49, ~135), `FittedVee` (molecular).  And
  `SinglesDensitySampler::MatrixT` refits the XC potential on a change of v (`itsFittedV`) — ↑ and ↓ blocks ask with different v.
- What is already right: `Vxc_Quadrature` warms its ρ in the eager refresh (`RefreshForDensity` → `WarmForDensity`), and
  `CompositeWF` calls `ham.RefreshForDensity(itsBS, cd)` before the block loop — the hook exists.
- **RULED + BUILT (user 2026-10-27: "Require RefreshForDensity before any GetMatrix everywhere; GetMatrix only asserts the
  stamped version").**  `HT_SlotOwner::PrepareSlots(bs, cd)` stamps the density (a dynamic cache drops its values on a NEW one;
  idempotent for the same); `tHamiltonian::RefreshForDensity(bs, cd)` = PrepareSlots + every term's refit; every `GetMatrix`
  (dynamic, no-cache, real-block) THROWS on an unstamped density; `GetTotalEnergy` refreshes first; the in-loop `newCD` writes
  in Hubbard/Hartree `MakeMatrix` are gone.  The sweep found the callers: 28 failures, all TESTS driving terms directly (plus a
  null-basis segfault in `Hubbard_U::PrepareSlots`) -- no production path skipped the refresh.  Sweep 1102 pass, 0 of 1100 moved.
  The stateful design itself is OOD debt **V1.41** (`doc/OOD-SOLID-Cleanup.md`): the cure is an immutable `H.At(cd)` value.
- *(the question as it was put)*: make the eager phase the ONLY place a term looks at a new density.
  `RefreshForDensity(bs, cd)` refits (Hubbard occupations, V_H, the XC fit per spin) and STAMPS the cache version; `GetMatrix`
  then only fills its own pre-created slot, and asserts the version instead of resetting.  That turns R1.0h's "a pre-warm,
  not a replacement" into a contract, so a caller outside the SCF loop (energy evaluation, unit tests) must refresh first — the
  question is whether to keep the self-correcting fallback for those callers (serial only) or require the refresh everywhere.

**P2 threading — what remains after the refresh contract (2026-10-27), two hazards one level down, each needing design:**
- **The XC refit**: `SinglesDensitySampler::MatrixT` holds ONE fitter and refits whenever v changes; ↑ and ↓ blocks ask with
  different v, so concurrent blocks would refit each other's fitter mid-contraction.  Fix: per-spin fits prepared in
  `Vxc_Quadrature::RefreshForDensity`, `MatrixT` a read (a small sampler face, or a v-keyed fitter set).
- **The gather memo**: PG_Cart `itsIntegrateMemos` (LRU of B_ij(n) per exact field) is written during the SCF gathers, and an
  insert ERASES the oldest -- freeing a memo another thread may be reading.  It exists so the k-blocks share one field's
  work, so disabling it under threads costs ×n_k.  Fix: build it in the refresh phase (each term's refresh warms its field's
  memo), the blocks only read -- plumbing from the term refresh to the evaluator.
Payoff on today's workload is small (NiO: the SCF was ~3% of a q-mesh run; the response is P1-threaded).  Both fixes land more
naturally inside V1.41's immutable `H.At(cd)` value, which owns exactly this density-derived state.  **PARKED 2026-10-27 (user): the SCF threading rides V1.41**; the plan continues with P3 (Q9).

**P2 LANDED 2026-10-28 on V1.41 stage 1** (`tHamiltonian::At(bs, cd)` -> `tHamiltonianAt`, OOD-SOLID-Cleanup V1.41).
`tCompositeWF::AssembleFock` runs ONE TASK PER BASIS BLOCK (its spin children together: a spin-independent term keys one slot per
block), the FIRST block alone (it builds the shared memos with the inner threads), the rest through `OuterFor`; and the wave
function's FIRST assembly is wholly serial -- found by the 2-thread sweep: concurrent `MakeKinetic -> findΩ -> Cache2` segfaulted
on every multi-k test, because the integral cache DB (Cache2/3/4, keep-newest eviction on the lattice) is not thread-safe and every
geometry-fixed lazy (statics, 3-centre tensors, projectors, Phi tables) is filled on that first build.  The two hazards above,
fixed: the XC refit -> a field-fit POOL (`SinglesDensitySampler::FitFor`, one fitter per v, recycled only when unheld; the Pair
sampler's ball fallback locked); the gather memo -> `IntegrateMemo` entries are locked `shared_ptr`s, a miss builds privately and
publishes when complete (so no reader is freed under, and a throw leaves nothing half-built).  Also locked: the Hubbard projector
lazy (the seed iteration's blocks fill it), XC `itsFockSource`, the KB per-l diagnostic bank, `DumpH`.
- Gates: UTResponse `ResponseKernel.GPW_Si_k311_FockIsBitIdenticalAcrossThreads_{Raster_Pol,Becke_UnPol}` (the loop, serial vs
  threaded, `==` on every element; from a CONVERGED state, so it cannot see first-build lazies -- the threaded sweep does).
  Serial short sweep 1112/1112, 0 moved (vs `q9_short` less `e280baec`'s 3; banked `~/Code/qchem6-runs/anchors/v141a_short.log`).
  2-thread sweep: no crash; moves <= 3e-10 and a few restart iteration counts, Γ tests included (width-0 outer loop) -- the inner
  regions' arrival-order sums; `Deck.ARunFromADeckEqualsTheSameRunBuiltByHand_Si_Gamma` (EXPECT_DOUBLE_EQ of two runs) varies run to
  run at 2 threads for the same reason: it holds serially only.
- Measured, TiO₂ k222 free (8 blocks), 6 iterations, Etot identical to every digit: 1 thread 906 s; 6 threads with the outer loop
  459 s vs 481 s inner-only (A/B in one binary) -- −4 %, the Fock assembly being a few s of a ~65 s iteration here (diagonalize +
  density dominate).  Inside an outer item the `h ball closure` (per-level `ΔG_Map` restriction + inverse FFT) costs 6.6× per call
  (0.052 -> 0.34 s): contention (allocation-heavy maps), not oversubscription (pocketfft is single-threaded) -- a lever if the
  outer loop is to pay on a gather-bound run.
- LONG tier (2026-10-28): serial 7/7, 0 moved of 5 comparable vs `dpulay2_long` (banked `v141b_long.log`); 8 threads 7/7 pass, wall
  686 -> 265 s, 3 moved -- Si k222 q211 −6.2e-9 (iters 13 -> 15), k222 q222 Becke −6e-10 (13 -> 14), k211 χ FD +2 iters -- none with
  the outer Fock loop engaged (2..8 blocks at 8 threads: width < threads), so the inner regions' arrival-order sums again.  The P2
  gate as written ("0 moved at 8 threads") cannot be met by ANY threaded run while the inner regions reduce in arrival order; what
  holds is: serial 0 moved, the outer loop bit-identical to serial (the UT gate), threaded runs within the SCF tolerance.

**P3 (Q9) DONE 2026-10-27: one q per star.**  `qchem.Response.QStars` (`FoldQMesh`: the q-mesh under the ground state's ops;
each member = a representative + a channel permutation + per-channel lattice shifts) and `LinearResponse(Nq, …, const QStars*)`:
solve the representatives and ONE audit member, generate the rest by X_{q'}(gI,gJ) = e^{−2πi q'·(t_I−t_J)} X_q(I,J),
q' = W^{−T}q, and REFUSE (Why::Configuration) if the audit disagrees beyond 1e-6.  Facade: `HubbardLinearResponse(…, qStars)`,
deck `postSCF.hubbardLinearResponse.qStars: true` (the opt-in of ruling 1); Manifold resolution only.  The ground state's group
is the Shubnikov group of the decoration the run was BUILT with (all-zero spins = the space group).  Gate, Si k222 q222 (both
sites perturbed): 3 of 8 q solved + 1 audit; audit 1.5e-15–5.5e-15; U 2.07451 / 2.07444 eV identical to the full mesh; wall
287 → 189 s.  UT `ResponseQStars.SiQ222FoldsIntoThreeStars`.  **Also: the facade now HOLDS its lattice** (user: "SolidCalculation
lacks the lattice -- this seems very strange"): it kept only fragments (`cellMatrix` for the fingerprint, then the fold's op list
and site positions); now `std::optional<Lattice_3D>` + the resolved decoration, and every former fragment is asked of it.

**P4 DONE 2026-10-27: `OrbitalImage(g)` — the D_AO primitive.**  `qchem.Symmetry.Lattice_3D.BlochOperationRep`
(`BuildBlochOperationRep`, `ImageOfK`): the lattice sibling of the molecular `BuildOperationRep` — shell permutation × each
shell's own angular rep (Cartesian or spherical, so the view needs no `Tc` composition of its own) × norms × the Bloch phase
e^{−2πi k′·t} of the image centre's lattice wrap, k′ = W^{−T}k.  The face `BasisSet::OrbitalImageSource` (one method, per
BLOCK: the block knows its own k and cell, so the ruled `OrbitalImage(g,k)` became `OrbitalImage(g)` on block k); `tGPW_IBS`
answers it.  σ is NOT applied (it swaps spin blocks, the caller's).  Gates:
- S-unitarity M^H S_{k′} M = S_k, every op, general k (so every target block differs), plus the representation property
  M(g₂;k′)M(g₁;k) = M(g₂g₁;k) and a refused non-symmetry: diamond Si (48 ops) 2e-15; rutile β-MnO₂ and TiO₂
  (`valence_lowq_tm.bsd`), 16 ops, Cartesian AND spherical, 5e-15 / 4e-15 (`UTLattice_BS GPW.OrbitalImage_*`).  KB and
  V_short carry over to 1e-14 (`GPW.OrbitalImage_SiDiamond_PPMatricesCarryOver`).
- F-invariance on an imposed converged state: Si k222, Becke mesh, 76 (block, op) pairs (Γ/X/L little groups 48 + 16 + 12):
  F 2.9e-9, D 7.5e-8 (`UTResponse OrbitalImage.GPW_Si_k222_ImposedBecke_FockAndDensityCarryOver`).
- **Two findings, both OpenWork §3 rows:** (a) on the DEFAULT uniform-raster XC route the same imposed state fails F by
  2.8e-4 / 7.4e-3 Ha (ops off the primitive-axis permutations / the glides) — a raster-XC equivariance defect, not M;
  (b) a FREE run's XC mesh is not invariant (Becke F 5e-4, D 2e-3), so ruling 1's measurement as written would reject every
  op on every free run — a decision before P5.

**P5 STAGE 1 DONE 2026-10-27: the little group of (q, J) on the GATHER side.**  Deck `hubbardLinearResponse.littleGroupK: true`
(`SolidCalculation::HubbardLinearResponse(…, littleGroupK)`), opt-in until the real-cell gates below pass.
- `Response::MakeOrbitalImages` (ruling 1, measured): every op of the run's Shubnikov set, on every block, as the orbital-basis
  image U = C_{k'}† S_{k'} M_g(k) C_k; KEPT only if U maps orbitals onto orbitals of the same ε and f and is unitary, to 1e-6.
  Console: kept / tried and every rejected op's residual.
- `Response::LittleGroupReduction` (one per (q, J)): of those, the ops that keep q (pair → pair) and the perturbation up to ONE
  phase χ_g, MEASURED (residual ~1e-14), never derived; the orbits; the reduced packing — √orbit per representative (an
  ISOMETRY of the full symmetric subspace, so GMRES takes the full solve's steps) and the stabilizer projection on every pack.
- The solver's reduced operator: unfold → ToAO (the WHOLE δρ is still scattered) → `ChargeDensity::RestrictCoupling` (the §6
  decorator; this stage only narrows the coupling, so the kernel's pair loop gathers δF on the representatives) → ToMO on the
  representatives → unfold δF → R0 (its q = 0 δμ sums every block) → pack.  Measurement on the unfolded δD (exact for every
  row).  The kernel, the terms and `InducedFock` are untouched.
- **Gates** (`UTResponse ResponseLittleGroup.*`): Si k222 free Becke, Si-p on site 0 — q211 χ 1.7e-10 (short tier), q222 χ
  6.7e-9, χ₀ exact, iterations within 1 (`_Long`; wall 106 s full → 30 s reduced at 4 threads; little groups 24 at Γ, 6–8
  elsewhere; 3–4 of 8 blocks solved).  **The safety net**: on the uniform-raster route the measurement keeps only the 6
  primitive-axis permutations (the 42 others rejected at 4e-4–9e-3 — the §3 raster defect), and χ still equals full to 3e-15.
- ⚠ COLLINEAR: a σ op relabels ↑/↓ blocks with the spatial U (comments at `MakeOrbitalImages`, `LittleGroup.C`, P4's face).
- **TiO₂ GATE PASSED 2026-10-28** (`decks/campaigns/TiO2/tio2_k222_q222_littleGroupGate.json`, full and reduced on ONE ground state;
  log `~/Code/qchem6-runs/TiO2/tio2_k222_q222_littleGroupGate_20261028.log`, record r007): all 8 χ₀(q), χ(q) equal to every printed
  digit (8 sig.), U 4.09898 eV both, kernel applications 97 → 83.  ⚠ BUT the reduction was weak: the measurement kept only 8 of the 16
  P4₂/mnm ops — the other 8 missed by 1.37e-6 against tol 1e-6, while the reference's OWN eigenpair residual is 5.4e-6 Ha (SCF stopped at
  Δρ 1e-6).  So 12–16 of 16 blocks were solved.  The tolerance is tighter than the reference it measures: OPEN (OpenWork handoff).
  **DIAGNOSED 2026-10-28** on the k221 ladder (same deck at k 2×2×1, q Γ; `~/Code/qchem6-runs/TiO2/k221_symmetryLadder/`), worst kept /
  ops kept: BallOnly Δρ 1e-6 → 7.0e-7, 4 kept (rejects all ≈1.0e-6 on the M-block eigenvalue, τ and non-τ ops alike); Pulay-on-D
  (no Kerker/cusp fit) → 8 kept, same Etot; manifold on both Ti → identical (Hubbard projector cleared); AliasFree Δρ 1e-6 → 16, 3.5e-7;
  **AliasFree Δρ 1e-9 → 16, 9.0e-9** (Si-level); **BallOnly Δρ 1e-9 → 16, 9.96e-7**.  Two causes: the SCF stopping level, and a ~1e-6
  equivariance floor from BallOnly aliasing (also 0.54 mHa in Etot; the OpenWork ⛔ raster row).  Nothing wrong in the image machinery.
- NOT YET: time reversal (q → −q); the NiO free AFM gate (the first σ-op exercise); STAGE 2 (scatter the representatives only and symmetrize each result inside `RestrictCoupling` — the |G|-fold
  scatter saving; stage 1 is its oracle).

## 5. Rulings (user, 2026-10-27)

1. **Assertion policy.**  A free run is reduced only by the subgroup its ground state is MEASURED to have: once D_AO exists (P4),
   check ‖D_{gk} − D_AO D_k D_AO†‖ per op per block and reduce by the ops that pass (tolerance ~ the SCF's, 1e-6; the console
   names the ops kept and their residuals; none pass → full mesh, said so).  Before P4 (i.e. for Q9) an explicit deck opt-in,
   audited by computing two members of one q-star.  An imposed run passes by construction.  This is MEASUREMENT, not imposition:
   nothing projects the SCF.
   **1a (user, 2026-10-27, after P4's finding (b)): a free run gets symmetric NUMERICS, not a symmetric answer.**  Its XC mesh is
   built invariant under the DETECTED group (the same `UnitCell::CreateIntegrationMesh(mp, ops)` an imposed run uses) but with
   NO fold, so nothing averages the density; a measurement failure is then physics, never the integrator.  Rejected: loosening
   the tolerance to the grid's asymmetry (it would hide the small breaking the test exists to catch); reducing imposed runs only
   (free runs second-class).  Landed: `GPW_IBS::itsMeshOps` (`UTResponse OrbitalImage.GPW_Si_k222_Becke_FreeRunCarriesOver`,
   384 pairs, F 2.6e-9 / D 3.4e-8 — was 5.1e-4 / 2.3e-3).  The uniform-raster route still fails (OpenWork §3, its own defect).
2. **Threading level: one level per loop.**  Outer (blocks) when the width ≥ the thread count, else inner (pairs) as today; no
   nested OpenMP.  Outer threading with the fixed-order block sum is BIT-IDENTICAL to serial (today's inner scatter is not:
   its per-thread buffers combine in arrival order).
3. **Order: P0 → P1–P2 (KP) → P3 (Q9) → P4–P6 (D5)**, re-checked against the 2026-10-27 SrVO₃ k444 run's peak RAM: if it is
   RAM-bound, P4 → P5 → P6 move ahead of P3 (only D5-b reduces RAM).

## 6. Interfaces (ruled 2026-10-27): symmetry enters as VALUES and DECORATORS, never as accessors

| phase | refused (implementation details on an interface) | instead | interface change |
|---|---|---|---|
| P1–P2 | `IsThreadSafe()` on terms; `GetBlock(i)` on `CompositeWF` | thread the loops where they already are | none; the CONTRACT "after `RefreshForDensity` a term's const queries are concurrency-safe" is written down; a `mutable` memo that breaks it is fixed inside its class |
| P3 | `GetSymmetryOps()` on `Reference`/WF; site positions on `HubbardChannels` | a q-star fold VALUE (members as op + channel permutation + Bloch phases), built by the `Response.C` adapter | none on a face; the q-mesh response takes the fold as an argument |
| P4 | `GetShells()`/`GetCenters()`/Cartesian powers on the orbital basis | the basis ANSWERS: a new cross-cast capability face, one method, "the matrix carrying block k's orbitals to gk under g" (Gaussian lattice basis; the spherical view composes `Tc`; a PW basis would permute G) | one new face |
| P5 | `Symmetrize(ops, q)` on every term; `GetBlocks()` on `TransitionDensity`; `GetFold()` on `TransitionProjector` | a DECORATOR implementing the existing transition-density faces, built where δ is built (a sibling of `AO_TransitionDensity_Factory`): scatters the reduced blocks, symmetrizes each RESULT (Hartree ΔGq_Map: free G algebra; XC values: the sampler's face below; Hubbard rows: δD rotated with P4).  The kernel, the terms and `ResponseKernel::InducedFock` never learn about symmetry | one new face on the XC sampler, one method: "average this sampled q-field over these ops, with these phases" (it owns its `FoldedMesh`) |
| P6 | `GetIrreducibleBlocks()`/`GetStarOf(k)` on the WF | a DECORATOR over the const WF face: an UNFOLDED view answering `GetQNs()`/`GetOrbitals()` for every mesh point (rebuilt with P4, cached); `MakeReference`, the response and a CK-3a restart walk it like a full-mesh run | none |

Already true, so nothing to build: the +U response is frozen-only (`Hubbard_U::RefreshForDensity` throws otherwise), so δF has no
+U part to symmetrize; only the measured χ rows are orbit-averaged (reusing P3's channel permutation).  `TransitionDensity.C`'s
"NEVER symmetrized" stays true of the FULL group; it gains the qualification (the decorator uses the subgroup that leaves THIS
perturbation invariant).

**LSP watch-points.**  (a) `TransitionDensity::Channel(spin)` of the decorator under a σ-op (δρ↑ ↔ δρ↓): a symmetrized channel
depends on BOTH inner channels — the views must symmetrize across them, never wrap one.  (b) The unfolded WF view is read-only
and post-SCF: it refuses to be iterated by an SCF rather than pretending to be updatable.

**Naming** (the client-consumes rule): `OrbitalImage(g, k)` for the basis face, `AverageOverOps(field, action)` for the sampler
face (not `RotationMatrix` / `Symmetrize`, which name mechanisms).  The little-group ACTION value lives in
`qchem.Symmetry.Lattice_3D` beside `BZReduction` and `Fold` (values about groups acting on meshes and sites).

**A temporary SOLID violation is allowed when stuck** (user, 2026-10-27) — but it is LOGGED the same day in
`doc/OOD-SOLID-Cleanup.md` (what, where, why, the intended shape), to be judged later in a wider context.
