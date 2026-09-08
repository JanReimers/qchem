# Open Work — the live tracker (v2, cut 2026-08-19)

**READ THIS AT SESSION START — and read `▶ WHAT IS OPEN` below first: it names the ONE next action, then
lists the open work, each row with its next concrete action and the one section to point a session at.
CLOSED items are not in that table (they keep their sections further down).**  Everything
after it is a mix of open items and the evidence that produced them; the evidence is kept deliberately
(three open items exist because a measurement refuted the obvious answer) but it interleaves, so the index
is the map.  Closed threads are retired to two history files:
`doc/OpenWork_History1.md` (threads A–E, the 2026-06-30 and 2026-08-15 orderings, runtime rounds 1–4) and
`doc/OpenWork_History2.md` (cut 2026-08-25 — the Vxc repair thread and the three-part fit-basis interface
refactor of 2026-08-21 → 08-24).  Durable design rulings live in `doc/CleanupCandidates.md` R1.0, not in
either history.

**Cut again 2026-09-08 → `doc/OpenWork_History3.md`** (user: *"done and todo items interleaved and lots of
history so it is hard for me to read and assess"*).  Eight sections — 1316 of 2859 lines, every one CLOSED,
ACTED ON, or self-labelled *no action here* — moved out whole, each leaving a stub that says what it
concluded.  Verified line-for-line: nothing lost.  ⚠ **What is STILL interleaved, deliberately**: `The plan,
in order` Steps 2, 3, 5 and 6 mix landed and open items inside one section each; splitting them needs a
judgement call per bullet rather than per section, and the table above is the map in the meantime.

▶ **THE STANDING RULE, same as `doc/CleanupCandidates.md`'s**: a section's ✅ verdict goes to the history
file the day it is written.  A closed section left in the tracker is indistinguishable from open work.

---

## ▶ WHAT IS OPEN — START HERE

> **▶ NEXT SESSION — ONE NEXT ACTION: SIZE THE BECKE GRID (put the uniform probe on the existing Becke
> ladder).**  Bin 1's like-for-like gap is one gather and that gather is behind N4; bin 2 is large but is
> gated on the grid RECIPE, not on the build's code — see item 1 below for both.
>
> ★★★ **THE PARITY GAP IS CALL COUNT, NOT KERNEL — MEASURED 2026-09-05 (§5f).**  Per call our gather is
> **0.88×** CP2K's and our collocation **0.93×**, and both codes spend ~99% of the run in those two
> routines.  But a fixed-point iteration issued **4 gathers + 2 collocations** where CP2K issues **2 and 2**
> (read off its own `T I M I N G` block: `integrate_v_rspace` 88/44 steps, `calculate_rho_elec` 90/44).
>
> ✅ **LEVER A LANDED (`9f4f4ae2`) — the Hartree ENERGY stops building a matrix.**  Two of our four gathers
> were the SAME TERM at two densities (the Fock's \f$V_H[\rho_{mix}]\f$ and `Vee_Hartree::GetEnergy`'s
> \f$\tfrac12\mathrm{Tr}(DV_H[\rho_{new}])\f$).  Parseval kills the second one:
> \f$E_H=\tfrac12\Omega\sum|V_H|^2/k\f$, no matrix.  **Parity probe 341.4 → 287.5 s CPU (−15.8%), gathers
> 95 → 66, 2.05× → 1.72×**; Si 2×2×2 rows −19.5% / −11.4% (\f$V_H\f$ is now memoized across irrep blocks
> too).  806/806, every iteration count and \f$E_{tot}\f$ reproduced.
>
> ⛔ **LEVER B IS REFUTED — do not build it.**  \f$V_H\f$ is a BALL field and \f$v_{xc}\f$ a RAW RASTER
> field, on purpose; routing \f$V_H\f$ through the raw adjoint changes the Hartree block by **6e-5
> relative** (measured on Si Γ), because the two adjoints truncate differently (per-level \f$\{G\}\f$ ball
> vs spectral box) — it breaks the very adjointness lever A rests on.  ⇒ **B becomes available only if XC
> gives up the raw \f$\rho_{DM}\ge0\f$ feed, i.e. only if N4 lands.**  Filed under N4, not as a bin-1 item.
>
> ★★★ **AND THE LIKE-FOR-LIKE NUMBER IS 1.13× (measured 2026-09-06).**  Our recipe anneals `Ladder,GDM`;
> CP2K's benchmarked decks DIAGONALISE AND MIX (`&DIAGONALIZATION` + Broyden/DIIS — checked in the decks and
> in every log's update-method column; **none of them run `&OT`**).  So only our FIXED-POINT stage has a
> counterpart there.  Differenced on its own (`GPW_MNO_NMAX` 6 against 2, single stage): **3 gathers + 2
> collocations per iteration against CP2K's 2 + 2 = 9.35 s against 8.29 = 1.13×**, and the whole residual is
> the third gather — lever B.  The two-stage probe's **1.72×** additionally carries the GDM line search's
> trial densities, which have no counterpart in a run that is mixing rather than minimising.
>
> ⏸ **LEVER C (the trial densities) IS A TODO UNTIL WE HAVE OT** (user, 2026-09-06: *"CP2K uses the OT
> method instead of GDM so we can't do proper parity timings against CP2K anyway"*).  It becomes a real
> question when OT exists and can be timed against GDM, minimiser to minimiser.  Do not optimise it against
> a comparison that does not exist.
>
> ⇒ **BIN 1's remaining like-for-like gap is ONE gather, and it is behind N4.**  So the live work is bin 2
> (below), and the row to quote is `CP2K_COMPAT=1`: **1.13× on the comparable stage, 1.72× for the
> two-stage probe.**
> **THE TABLE THAT SAYS SO IS FILLED**: `doc/Benchmark.md` **§5a**, now its own section holding ONE table
> (user, 2026-09-05), nine rows, BOTH codes pinned serial, and setup split out of the per-iteration figure.
>
> **Per SCF iteration we are ahead of CP2K on seven of nine rows** — MnO defaults **0.90×**, MnO FM
> **0.83×**, `QCHEM_BECKE_XC=0` **0.53×**, Si **0.13× / 0.84× / 0.68×**, NaF full-SR **0.19×**.  Two rows
> are behind: NaF SR2 **1.47×** (a Becke cost — really bin 2) and **`CP2K_COMPAT=1` 2.05×**.
> ⇒ **Single-thread parity now IS that one row**, because it is the only row with OUR accelerations off, and
> nothing else on this list moves it.  Its setup is 1.73 s (1%), so bin 2 cannot help it; its cost is the
> box walk (§5e: ~95% of that run, 74% of it the gather).
> ⚠ It is a `GPW_MNO_NMAX=10` capped probe with a different stage mix (rule 3d), so **step one is an
> uncapped or equal-cap comparison** before optimising against it.
> ⚠ `QCHEM_BECKE_XC=0` is ONE of §2's seven deviations, **not** parity (user, 2026-09-05); and
> `CP2K_COMPAT=1` is our best KNOWN parity, not proven parity — the list has grown 4 → 7 every time anyone
> looked, so **finding deviation #8 is part of this item**, not a distraction from it.
>
> **THEN, in order:**
>
> **1. ⛔ BIN 2 — THE BECKE MESH BUILD IS GATED ON ITS GRID SIZE, NOT ON ITS CODE** (user, 2026-09-06:
> *"I am reluctant to work on that until identify the proper becke grid size (Nradial=40, Nangular=29 is
> big) that yields the same accuracy some metric (‖Vxc−Vxcfit‖?) as the uniform 20³ mesh.  Becke is always
> more expensive to setup, no way around that.  It is parallel which helps."*)
>
> The cost is real — MnO's setup is **184.3 s = 47% of the DEFAULT run** against CP2K's 8.1 s (**23×**), of
> which **136.6 s is TWO Becke mesh builds** (68.3 s each) + 43.1 s of XC-mesh Φ tables — but the first
> question is how many points the recipe actually needs.  `nRadial=40, degree=29` sets the entire Becke side
> of every cost comparison in this tracker and the build scales with it directly, so a recipe that is 2×
> oversized is a 2× saving with no engineering at all.  ⇒ **CALIBRATE FIRST, AND POSSIBLY FIND THERE IS
> NOTHING LEFT TO OPTIMISE.**  ⚠ It does NOT touch the parity row either way (setup 1.76 s there, BEATING
> CP2K's 8.1 s).
>
> ✅ **MEASURED 2026-09-06 — the ladder now carries the uniform route as its yardstick** (`BeckeLadder()`
> in `IntegrationTests/GPW_SCF_UT.C`, V2.8 block; one frozen density per system, reference nR=100 GL-41,
> scored by \f$\Delta E_{xc}\f$ and \f$\max|\Delta V_{xc}(i,j)|\f$ — the error in the operator that is
> actually diagonalised).  Both routes are FITS of the same \f$v_{xc}\f$, and a fit is **(integration
> grid) × (fit basis)**, two INDEPENDENT axes: the uniform route is {uniform raster}×{plane-wave
> \f$\{G\}\f$}, the Becke route is {Becke mesh}×{delta basis}, both metrics orthonormal — which is what
> makes one scoreboard legitimate.  ⚠ The axes must stay orthogonal IN THE CODE (user, 2026-09-06): a grid
> does not imply a basis, and which pairings are worth using is high-level POLICY, never hard-coded:
>
> | system | UNIFORM: pts, max\|dVxc\| | cheapest Becke rung ≤ that | production nR=40 GL-29 |
> |---|---|---|---|
> | Si covalent | 15625, **3.54e-4** | nR=40 GL-15, 6968 pts (**3.5×** fewer) | 24472 pts, 1.75e-5 |
> | NaF ionic | 15625, **1.617e-1** | nR=40 GL-5, 920 pts (**26×**) | 23656 pts, 7.92e-5 |
> | Al metal | 512, **5.75e-3** | nR=40 GL-15, 3448 pts (**3.5×**) | 12074 pts, 2.58e-4 |
> | Mn atom-in-box | 262144, **1.303e-1** | nR=40 GL-5, 556 pts (**25×**) | 14104 pts, 1.31e-6 |
>
> ⇒ **(1) THE RECIPE IS OVER-GENEROUS ON EVERY SYSTEM** — 3.5× on the two hard ones, ~25× on the sharp
> ones — and the build scales with the point count, so that factor comes straight off the 136.6 s.
> ⇒ **(2) THE ENERGY AND THE MATRIX DISAGREE ABOUT WHICH ROUTE IS BETTER.**  The uniform route's
> \f$\Delta E_{xc}\f$ is tiny everywhere (3e-7 on Si — better than every Becke rung below nR=60) while its
> matrix error is 20× to 10⁵× worse.  Error that cancels in the integral does not cancel in the operator.
> **Score the matrix.**
> ⇒ **(3) SO "MATCH THE UNIFORM MESH" IS A WEAK TARGET EXACTLY WHERE BECKE EXISTS TO HELP** — on NaF and
> Mn it would license GL-5, which nobody should ship.  ▶ **The threshold is a POLICY CALL and needs the
> user**: at an ABSOLUTE \f$\max|\Delta V_{xc}|\le\f$ 1e-4 the answer is **nR=40, GL-17…21** — still
> 2–3× cheaper than production.
> ⚠ Both V2.6a rules survive: Al is NON-MONOTONIC on both axes here too (GL-9 beats GL-11, nR=25 beats
> nR=30), and a frozen ladder UNDERSTATES the self-consistent shift on a metal ⇒ **no default may be
> flipped on ladder evidence alone; it needs a converged A/B on Al.**  ⚠ And MnO — the system whose
> 136.6 s started this — has no ladder of its own; the Mn sextet ATOM is a proxy for its sharpness, not
> for its partition.
>
> ⇒ **TWO ITEMS SURVIVE THE GATE and are worth having whatever the calibration says**: (a) **why TWO
> builds** for the same cell — the anneal's two stages each build one, and that redundancy is independent of
> grid size; (b) the build THREADS (`src/Structure/Imp/UnitCell.C:273` — 68.3 s serial against the 16.7 s
> the 09-04 threaded ledger read, ~4×), which is the user's *"it is parallel which helps"* and which still
> collides with the NaF Amdahl reading below (`doc/Benchmark.md` §5e).
> ⚠ The NaF threading finding that promoted bin 2 still stands as written: on NaF the Amdahl-inferred
> serial time (9 s) matches the measured setup buckets (9.7 s, of which the Becke mesh build alone is 7.0 s
> on a TWO-ATOM cell) to ~7%; the SCF threads essentially perfectly.
>
> **2. ✅ DONE 2026-09-06 — THE THREADED TABLE IS FILLED, BOTH SIDES (`doc/Benchmark.md` §7), and it
> INVERTS the question.**  CP2K at 12 OMP threads is **0.82–1.09×** its own serial on these decks (its
> banner says 12 threads; it measures **196% CPU**), and 12 MPI ranks buy **1.44× for 8.3× the CPU**.
> ★ Its own timing block localises it: `grid_collocate_task_list` **1.09×** and `grid_integrate_task_list`
> **1.12×** — the two routines that are 98% of the run.  ⚠ Whether that is "the GPW route's OMP was never
> a priority" (their centre is hundreds of small molecules, where MPI over molecules is the axis) or "a
> 4-atom cell has too few tasks to spread" is NOT separated by this measurement; the discriminator is a
> supercell deck we do not have.  Ours:
> **2.10–4.67×** (best on the pure box-walk parity route).  Cross-code at 12 cores on MnO: **4.14 s per SCF
> step against CP2K's best 5.90 s (0.70×), on 596 s of CPU against 3107 s (0.19×)**.  ⇒ **We are not
> missing parallel opportunities CP2K exploits — on THIS class of system** (4-atom, high-symmetry: the
> regime §2 says most favours us; a 100-water box would invert it, and none of this is a many-node claim).
>
> ▶ **SEQUENCED IN `doc/ParallelAndOraclePlan.md` (2026-09-06) — start at its 1.1.**  The opportunities we
> ARE missing are our own (§7c), in priority order: (1) **~59 s of UNBUCKETED
> work** in a 128 s threaded MnO run — diagonalise/ortho/mix/fit-solve, none of it timed, so the first move
> is an INSTRUMENT not an optimisation; (2) the **XC-mesh quadrature GEMM at 1.21×** — ⚠ a
> DELIBERATE trade, not an oversight (user): our parallelism lives ABOVE the linear algebra (per
> k/irrep/spin) with BLAS pinned to one thread to avoid OMP nesting, and one dispatched whole-matrix
> `zgemm` measured 34.1 GFlop/s against 1.87 for any blocked form.  The thing to question is the WIDTH of
> the level above — at Γ with 2 spins it is 2-way, so ten cores idle in that bucket by construction.
> ▶ **Cheap test worth running**: `QCHEM_BLAZE_BLAS=ON` was correctly refused on 2026-08-15 against NETLIB
> BLAS; `libblas`/`liblapack` now resolve to **openblas-pthread**, so the 34 GFlop/s single-threaded path
> is available again and composes with the pin.  One rebuild, one run; (3) the **\f$V_H\f$ field
> build at 0.89×**, i.e. `SymmetrizeGMap`'s serial 48-op star-average, now 6.6% of the threaded wall.
> ✅ And §5e's open question is CLOSED: **the Becke mesh build threads at 8.2×** on MnO.  ⚠ NaF is still the
> anomaly (serial fraction = its setup) — likeliest SIZE, not structure; settle it with §7c's bucket table
> taken on NaF.
>
> **3. The cheap gaps**: ✅ MnO FM re-taken 2026-09-05 (411 s CPU / 474 MB against the stale 2321 s /
> 4947 MB; 0.83× CP2K per SCF iteration) · ✅ NaF's CP2K iteration counts read off its logs (SR2 16 steps,
> full-SR 27) · ⬜ an attribution A/B against the parent commit on the VA recipe.
>
> **4. Only after the serial fraction is down**: the residual 1.33–1.45× CPU inflation at 12 threads
> (load imbalance / memory bandwidth).  Not worth chasing while setup dominates.
>
> ✅ **CLOSED 2026-09-04, do not re-open**: the busy-wait barrier (`StopOmpThreadsBusyWaiting` — libomp
> spun 200 ms after every parallel region, 65% of billed CPU; §7 is filled because of it) · the gather's
> D-screen (removed: it was a self-fulfilling truncation of the Fock and broke the fold's
> orbit-invariance) · the LatticeScreener seam · the one-gather XC term · `-march=native` as the default.
>
> ⏸ **PARKED — the ‖V_xc − V_xc_fit‖ fit-quality study** (user: *"defocusing"*).  Still true that every
> Becke-vs-uniform cost number in the tree is taken at unknown-equal accuracy, so "Becke is a negative
> acceleration" is a COST statement and not a verdict — but it is not what bins 1 and 2 need.
> ⚠ Note it now COLLIDES with the bin-2 finding: the Becke mesh build is both the biggest setup cost and
> the threading ceiling, so whoever attacks setup will be standing next to this question anyway.
>
> ⚠ `doc/OldPlans/ScreeningPlan.md` is CLOSED (2026-09-04) — read it for WHY, do not take work from it.

Everything after this index is a MIX of open work and the evidence that produced it.  The evidence stays on
purpose — it is what stops items being re-litigated, and several items exist because a measurement refuted
the obvious answer — but it interleaves, and you cannot see the shape of the work by scrolling.
**This table is the whole of the open work.**  One row per OPEN item, in priority order, each with the next
CONCRETE action and the ONE section to point a session at.  If it is not in this table, it does not need doing.

⚠ **CLOSED items are NOT listed here any more** (2026-08-27, user: *"START HERE is now unclear because it
mostly lists done stuff"*).  An index whose rows are mostly ✅ cannot be read for what to do next, which is
the one job it has.  The closed threads keep their SECTIONS below — the evidence is the point — they just
stop competing for the reader's attention here:

| closed | what it settled | its section |
|---|---|---|
| **A — the on-the-fly box walk** | 2026-08-26/27: **14.5× on MnO** (570 → 39 s per iteration), CP2K standing 67× → ~5×.  Four bit-identical edits, then the separable-contraction kernel (`GPW_CONTRACT_CUBE`, built + gated OFF pending the anchor re-bank). | *"THE ON-THE-FLY BOX WALK"* + `doc/CollocationRewritePlan.md` |
| **N5** | `CP2K_COMPAT` + the self-describing banner.  ⚠ Remainder: `raster`/`cutoffFactor` are TYPED options the policy does not reach. | *"✅ T5 / N5"* |
| **N1** | T1–T5 all built.  ⚠ Remainder: the COVERAGE GAP (`RunGpw`/`RunGpwAnnealed` bypass the facade). ★ The detector RULES are now unit-tested (`src/Calculation/tests/RunDiagnostics.C`). | *"★★★ N1"* |
| **N2** | The ρ<0 lobes are BAND-LIMITING, not aliasing — which eliminated every cheap alternative to N4. | *"★ N2"* |
| **1** | `GPW_XC_DM_SOURCE` does not earn the default; the measurement convicts the MIXER, not the exact ρ. | *"✅ ITEM 1 MEASURED"* |
| **3** | The imposed XC mesh keeps its site blocks. | (below) |
| **C — the collocation rewrite** | 2026-08-27, steps 7–8: the 3.9 GB pair-stream VALUE cache is deleted and the (shell pair, offset) TASK LIST replaces it (~0.2 MB); `GPW_CONTRACT_CUBE` defaults ON. MnO peak RSS **3915 → 155 MB** with the box-walk buckets ~1.1× slower on the whole run — the trade the cache existed for had evaporated. 792/792 green on BOTH kernel settings. Two latent defects fell out: the integrate-back's `Re[D·conj(phase)]` screen (the shifted-MP defect, 4.1 Ha, previously reached only by over-budget pairs) and the walk's per-component `|v|` screen (broke the collocate/integrate adjoint at 1.2e-8 once the shared frozen stream was gone). | `doc/CollocationRewritePlan.md` step 7 |

| # | open item | the next concrete action | point a session at |
|---|---|---|---|
| **SCR** | ✅ **THE SEAM IS BUILT (2026-09-04).**  `LatticeScreener` + `GeometryOnlyScreener`/`DAwareScreener` in `src/BasisSet/Molecule/LatticeScreener.C`; both collocation faces take a `const LatticeScreener&`; the `GPW_DAWARE_SCREEN` bool is gone from the box walk and survives only as `RunPolicy::DAwareScreen`, a declared CP2K deviation.  D-aware stays the default; suite unchanged.  ⚠ The `DensityHandle` proxy this row anticipated was NOT built and should not be: the screener is **stateless** — the walk already computes each term's weight and hands it in, so no density, no reseat, no staleness (`ScreeningPlan.md` §4). | ⛔ **§5 IS CLOSED — REFUTED ON MEASUREMENT (2026-09-04), do not build it.**  `M_PG_BoxWalk.WhatTheGeometryHoistWouldBuy` prices the hoist CEILING at **13.2% of the kernel** (chord share 10–17% across box sizes, not the ~40% claimed) against a **+14.5% wall** price for the geometry-only screener that makes it legal, plus **132 MB** on a run whose peak RSS is 110 MB.  Best case is a wash.  The ~40% was the per-LINE work, most of which is the \f$e_2\f$ fold — which reads the density-weighted coefficients and is not hoistable under any screener.  ▶ **The next lever is NOT in the kernel**: see the per-step field count below. | `doc/OldPlans/ScreeningPlan.md` |
| **S** | ★★ **THE ANCHOR-MOVING SPRINT — A1 and A7 ARE DONE (2026-08-27), A2–A6 remain.** Five items that each move banked numbers, to be done in ONE re-bank so they do not mask each other (user, 2026-08-27). | Pick the sprint window. A5 (the `IonicSAD` seed default) re-seeds every GPW anchor, so it goes first or last. A4 (the Δρ/N gate) is now doubly motivated — see the Na2 note in the sprint section. | *"THE ANCHOR-MOVING SPRINT"* |
| **N4** | ★★★ **THE RIGHT TREE: MAKE EVERYTHING ELSE ROBUST WITH \f$V_{xc}[\rho\ge0]\f$** (user, 2026-08-25).  ⚡ **AND IT NOW CARRIES A BIN-1 PRIZE**: `doc/Benchmark.md` §5f lever B (one gather per spin, CP2K's `sum_up_and_integrate`, worth ~1 of our 3 gathers per iteration) is blocked ONLY by XC needing the raw \f$\rho_{DM}\ge0\f$ feed — \f$V_H\f$ is a ball field, \f$v_{xc}\f$ a raw raster field, and routing \f$V_H\f$ through the raw adjoint moves the Hartree block by 6e-5 relative (measured).  If N4 makes the ball XC route safe, B becomes exact and free. *"ρ̃_mix is not exactly garbage … but it is still pretty junky for Vxc"*, and improving the junk (N2) is barking up the wrong tree. ⇒ **"the flag does not earn the default" was the wrong headline for the right measurement**: what failed is the MIXER, not feeding \f$V_{xc}\f$ the exact ρ. | Build the **CUSP-DEFICIT** form \f$\rho_{XC}=\rho_{mix}+(\rho[D]_{exact}-\rho[D]_{BL})\f$ — XC keeps Hartree's OWN mixed array, so there is **no \f$\alpha_{eff}\f$ to choose** and the measured failure cannot occur. Plus **N3** (charge/spin channels) and **N1/T1-T3** (so a future collapse cannot masquerade as an answer). | *"★★★ N4 — THE RIGHT TREE"* |
| **N3** | ★★ **CHARGE AND SPIN NEED SEPARATE PRECONDITIONING — ⚠ HALF-BUILT ALREADY (corrected 2026-08-25): `QCHEM_MIX_RHO_M=1` in `MakePeriodicMixer` ALREADY selects the (ρ,m) basis with "Kerker on ρ, PLAIN LINEAR on m", carrying the same *"m has none"* argument. So this needs a MEASUREMENT and a promotion, not a build.** — Kerker is applied per spin channel, so by linearity it damps the SPIN channel too, and the spin channel has **no 4π/G² divergence to justify it** (user). It is charge medicine taken by the magnetisation; cf. VASP's independent `AMIX_MAG`/`BMIX_MAG`. | Split the mixing policy into charge + spin channels. ⚠ Do this KNOWING that today's AFM basin is propped up by the current behaviour (see ITEM 1 MEASURED) — so it needs the N1 detectors landed first, or it will look like a regression. | *"★★ N3 — THE MIXING POLICY"* |
| **2** | **BENCHMARK PROTOCOL — no timing table is comparable until this holds** (user, 2026-08-25). Two defects today: no table states its THREAD state per row, and qchem runs accelerations CP2K does not — the factored/low-rank ρ is **ON BY DEFAULT** (`QCHEM_DM_LOWRANK`), so every row since `07d13bf6` has it | (a) build the self-describing BANNER `doc/Benchmark.md` already asks for — thread counts + the qchem-only feature flags — so rows describe themselves instead of relying on discipline; (b) re-take the rows under the two-phase rule: **single-thread parity FIRST**, then N=8/16 for OMP-shaped gaps. | `doc/Benchmark.md` → *"BENCHMARK PROTOCOL"*, and Step 0 (instruments) |
| **4** | **Step 5 — MnO accuracy, name the operator**: the sharpest coordinate on the list, with a banked oracle, and its first move is cheap | ⚠ **PIN `GPW_XC_DM_SOURCE` first** — individual terms move ~100 mHa with it, so the term-by-term CP2K breakdown means nothing until item 1 is settled. Then the cheap first move. | Step 5 |
| **5** | ✅ **DONE 2026-09-06 — Step 0c, "the instruments report WHAT, not WHEN"**.  `report::RunElapsed()` (the run's own steady clock, zeroed at `Begin`) stamps EVERY emitted item: the console heading carries `[t=12.34 s]`, a scoped `Section` carries its whole SPAN `[t=3.21→12.34 s]` (the honest reading — it renders at scope CLOSE), `Log` and `[fold]` lines are stamped too, and the run record grows a chronological `timeline` array beside the nested document.  ⇒ **The design question — chronological stream vs nested render — is answered by keeping BOTH**, which costs one array.  5 unit tests in `UTCommon`. | ⚡ **It paid on its first run**: on the MnO ledger the gaps between stamps named a 35.8 s silent block (the stage-2 Hamiltonian rebuild) with no bucket added — see `doc/ParallelAndOraclePlan.md` 1.1(a). | Step 0 |
| **6** | **`FIT_SF_Ortho` — separate the metric axis into faces**: specced 2026-08-23, not built. `OverlapDiagonal` sits on the metric-NEUTRAL face, so `Fit_IBS` invents an answer in the wrong normalisation | Move it to `FIT_SF_Ortho<T>` (both fit faces in one increment) and `Fit_IBS` simply loses it — delete the landmine, do not correct it. ⚠ Acceptance criterion: must NOT become a `dynamic_cast` type switch. | *"★ SPECCED, NOT BUILT"* (near the end) |
| **OT** | ⏸ **ORBITAL TRANSFORMATION (OT) — the minimiser CP2K actually ships, and the ONLY way to time a minimiser against theirs** (user, 2026-09-06: *"CP2K uses the OT method instead of GDM so we can't do proper parity timings against CP2K anyway"*).  Two things hang off it: a like-for-like STAGE-2 comparison (today only our fixed-point stage has a counterpart — CP2K's benchmarked decks diagonalise and mix), and `doc/Benchmark.md` §5f **lever C**, the GDM line search's trial densities (42 of the parity probe's 82 collocations), which cannot be judged without one. | Not scheduled.  When it is: build OT beside GDM under the existing accelerator/loop-driver seam (`doc/SCFStrategyPlan.md`'s role seams already anticipate another direct minimiser), then time OT-vs-OT and re-open lever C. | `doc/Benchmark.md` §5f + §5a footnote ᵇ |
| **KP** | ★★ **K-POINT PARALLELISM — THE ONE PARALLEL AXIS EVERY OTHER CODE HAS AND WE DO NOT** (user, 2026-09-07: *"all these Γ runs are mostly for development, in the real world multi-k is the norm … if all the other codes do this and we don't because of internal caching, well then that is embarrassing"*).  ✅ **CONFIRMED IN CP2K'S OWN SOURCE**, not from memory: `PARALLEL_GROUP_SIZE` (*"Number of processors to be used for a single kpoint … the number of groups must divide the total number of kpoints"*), `kpoint%kp_range` assigning each group its k slice from `kp_dist`, and `para_env_inter_kp%sum(...)` reducing across groups (`src/input_cp2k_kpoints.F`, `src/kpoint_methods.F`).  ⚠ Believed true of QE (`-nk` pools), VASP (`KPAR`) and ABINIT (`npkpt`) as well — **NOT verified here**, no local sources; verify before quoting. | ⛔ **THE "WE CAN'T, BECAUSE OF OUR MEMOS" ANSWER DOES NOT SURVIVE CONTACT WITH THEIR DESIGN.**  Their groups are separate MPI ADDRESS SPACES, so each simply holds its own copy of the density-derived state — the duplication our objection treats as disqualifying, they accept, and it is cheap beside the per-k work.  ⇒ Our coupling is a CHOICE, not a constraint. ★ **AND THE FIX IS SMALLER THAN THE OBJECTION SUGGESTS**: every one of the 20 `mutable` memos on the term stack (V_H ΔG_Map, ρ rasters, XC-mix rasters, fitted v) is **k-INDEPENDENT** — in the pool model these are exactly the objects computed ONCE and read by every group.  The obstacle is not that they exist, it is that they are **lazily filled on FIRST BLOCK ACCESS**, which turns a read-only shared resource into a write-on-first-touch race.  ✅ **BUILT 2026-09-08 — `tHamiltonian::RefreshForDensity`.**  `tDynamic_HT` gains a `RefreshForDensity(cd)` hook (default no-op); `tHamiltonianImp` folds it over the DYNAMIC terms only (a static term is density-independent, so the phase must never reach it); `Vee_Hartree` warms \f$V_H[\rho]\f$ and the three periodic XC terms warm \f$\rho\f$ (or the \f${\uparrow,\downarrow}\f$ pair) through one `XC_Quadrature::WarmForDensity(cd, polarized)` entry point — expressed as ONE engine call rather than the term reaching in by name, so the warming AND the shape-exclusivity rule stay with the engine that owns the caches.  `tCompositeWF::DoSCFIteration` and `BuildFockAndComputeSteps` drive it in its own report bucket immediately before their block loops.  Gates: `EagerRefresh.*` in `UTHamiltonian` (3 unit tests, no SCF — the fold is pure plumbing over the term lists, so it is pinned by CALL COUNTS, not timings).  830/830.
> ⚠ **IT IS A PRE-WARM, NOT A REPLACEMENT, AND THE DOCS SAY SO.**  Every memo keeps its density-serial guard and those guards remain the correctness mechanism — energy evaluation and the unit tests drive terms outside any prologue, so an assert of the form *"the phase must have run first"* would be false. A test pins that non-contract deliberately (`AssemblyWithoutAPriorRefreshIsLegal`).
> ⛔ **AND THE BLOCK LOOP IS NOT YET READ-ONLY — one write remains, and it is a DIFFERENT problem.** `tDynamic_HT_Imp::GetMatrix` stores its result in `mutable CacheMap itsCache` keyed by `Irrep`, i.e. **one map entry written per block, inside the loop**.  That memo is k-DEPENDENT, so the eager phase cannot warm it by construction; it needs per-block storage (or no cache under a parallel loop), which is a separate increment.  ▶ Do not read this item as "the loop is now threadable" — it is one of two obstacles, and the one that was blocking the DESIGN.  ▶ The original argument stands: an explicit refresh phase — and is a better design regardless, since it makes the phase structure visible instead of implicit in call order.  Then: shared prologue → read-only parallel k loop → density reduction, i.e. CP2K's own decomposition.  ⚠ `itsByL`/`itsByLSeen` ("irrep blocks already decomposed") is the one genuinely per-block accumulator and needs separate handling. | ⏸ **NOT NOW, BY AGREEMENT, BUT IT STAYS ON THE LIST** (user: *"if it doesn't make sense to do it now that is fine … but it should stay on the list until we are suitably embarrassed that we have to do it"*).  Reasons to wait, not to drop: every row we own is Γ (width 2), so the payoff cannot even be MEASURED until the multi-k rows of **2.1/2.2** exist; and Phase 1's remaining lever (1.3b) is a bigger win on the cells we actually run.  ▶ **Re-open it at 2.1**, where the width is 16 and the question stops being hypothetical. | `doc/ParallelAndOraclePlan.md` (1.3 + Phase 2) |
| **BM** | ★ **THE BECKE MESH — FOUR FINDINGS; THE "OPEN BUG" WAS NOT ONE (2026-09-07/08, (4) CLOSED 2026-09-08).**  Came out of the supercell work but is INDEPENDENT of it; the supercell symmetry fix stands on its own evidence.  ⚠ **THE CAVEAT THAT TURNED OUT TO BE THE ANSWER** (user, 2026-09-08): *"The polyhedra truncations might make these counts difficult to interpret."*  Right instinct, and the mechanism was even simpler than truncation — the points are WRAPPED (see (4)) — but the ruling stands as a rule: **a raw count of distinct radii is a LEAD, never proof of a defect.**  It is the reason the 199-vs-49 was never acted on as a bug. | **(1) ⛔ A LIVE TRAP, already filed in `doc/CleanupCandidates.md`**: setting `MeshParams::cellKind=Becke` ALONE leaves the rest at the struct's own defaults — nR=30, α=1, **L=5** — where `BeckeXCParams`' are nR=40, α=2, **L=29**.  A degree-5 XC mesh is not a Becke run.  It cost a bogus 40 mHa "imposed vs free" discrepancy that read exactly like a symmetry bug.  **Ask for the RECIPE (`BeckeXCParams(-1,-1,-1)`), never the kind alone.**  ⚠ Related: `GPW_BECKE_L/NR/ALPHA` are consulted ONLY for arguments passed `<0`, so they silently do nothing against a caller-supplied degree — a sweep over `GPW_BECKE_L` produced three identical runs before that was noticed. **(2) ✅ AT THE PRODUCTION RECIPE IMPOSED AND FREE AGREE**: Si 1×1×1 Becke L=29, imposed \f$-7.11493826\f$ vs free \f$-7.114983942\f$ = **46 µHa**, and both sit ~0.1 mHa from the uniform-mesh anchor \f$-7.115067844\f$ (itself matching the banked \f$-7.11506\f$).  So the two XC routes and the two symmetry arms all agree; there is no energy-level defect. **(3) ★ A 2.9× SITTING UNCLAIMED**: `InvariantAngularMesh.StockLebedevIsAlreadyInvariantUnderSiTdSiteGroup` shows **Lebedev-29 (302 dirs) is ALREADY exactly invariant under Si's \f$T_d\f$ site group — 0 unmatched of 7248** rotated directions, because Lebedev rules are built from OCTAHEDRAL orbits and any axis-aligned cubic site group is a subgroup of \f$O_h\f$.  W2b nonetheless builds a site-adapted rule at **886 dirs/atom** (48128 mesh points against 16392), buying an invariance it already had.  ▶ Test the stock rule for site invariance FIRST and reuse it when it passes; keep the adapted construction as the fallback.  ⚠ Cell-dependent (big on high-symmetry cubic cells, vanishing as site symmetry drops — measure on MnO), and the invariance test must stay EXACT or a false positive silently reintroduces the bug W2b exists to prevent. **(4) ✅ CLOSED 2026-09-08 — THERE IS NO CORNER-ATOM DEFECT, AND THE SUPERCELL GRID IS THE PRIMITIVE GRID REPLICATED.**  Root cause of the whole lead: `MakePeriodicBeckeMesh` emits every point **WRAPPED INTO THE HOME CELL** (`kpt = r - A*n0`, UnitCell.C), so the stored coordinate is NOT \f$R_a+v\f$ and \f$\|p-R_a\|\f$ IS NOT THE OFFSET.  That single fact produced every number in the lead: 199 distinct "radii" on the corner atom, 49 on the interior one, and the 17.75 that is not a radial node.  The corner atom read worse for the obvious reason once named — an atom at (0,0,0) has its whole grid straddling three cell faces, so nearly every point wraps, while (¼,¼,¼) keeps more of its inner shells intact.  ▶ Recover the offset MODULO THE LATTICE (\f$v=p-R_a+An\f$, one \f$n\f$ in a bounded box) and the structure is exact: **both** Si sites, corner and interior, put **480/480 points on a radial node, 7 distinct radii, rMax = 10.8889 (a node), zero off-direction, zero ambiguous** — free AND imposed (868/868 there).  **The atom LABELS are right too, in BOTH cells**: each primitive block decomposes about its own atom 480/480 and about the other atom **0/480**, and on the 2×2×2 all 16 blocks take 868/868 of their own atom's nodes with a best WRONG-atom match of **0 points**.  And the user's real question, answered point by point in both settings: folding the 2×2×2 mesh back into the primitive cell, every one of its 16 site blocks matches its primitive partner **bijectively, max \f$|\Delta r|=6\times10^{-15}\f$, zero unmatched**, max \f$|\Delta w|=1.0\times10^{-7}\f$ ABSOLUTE and each site's Sum(w) equal to \f$1.3\times10^{-8}\f$ relative.  ⚠ **THE WEIGHT METRIC MUST BE ABSOLUTE, NOT RELATIVE** — the partition is an eps-converged (1e-6) image series gathered in Chebyshev CELL shells, and a supercell shell is 8 primitive cells with twice the interplanar floor, so the two settings truncate the same convergent series at different places.  A per-point RELATIVE comparison is noise in the tail and says so loudly: the worst relative deviation is 9.6% — **on a point whose weight is 3.8e-82**.  ⛔ **RETRACTED with the rest**: "the corner atom is markedly worse (199 vs 49)".  There was never a cell-imaging bug; the measurement was reading a wrapped coordinate. | ▶ **BM(4) IS CLOSED; (1) AND (3) REMAIN.**  (1) is the live `cellKind=Becke` trap, filed in `doc/CleanupCandidates.md`.  (3) is the unclaimed 2.9×: test the STOCK Lebedev rule for site invariance first and reuse it when it passes, keeping W2b's adapted construction as the fallback — measure the size of the win on MnO before spending anything, since it vanishes as site symmetry drops.  ▶ The gates that closed (4) are `BeckeMesh.*` in **`src/Structure/tests/BeckeMeshUT.C`** (6 tests, ~14 s in UTStructure, no SCF): the wrapped-product decomposition, the atom-label discriminator in both cells, the imposed radial decomposition, and the free + imposed supercell replication.  Anything that touches the Becke build's coordinates, wrapping, site blocks or partition should run them first. | `doc/SymmetryUpgradePlan.md` "SUPERCELLS" + this row |
| **PAR** | ★★★ **THE SEQUENCED PLAN FOR PARALLEL WORK + THE SECOND ORACLE — `doc/ParallelAndOraclePlan.md`** (cut 2026-09-06).  Phase 1 our own OMP gap (instrument the ~59 s unbucketed → BLAS-mode serial arm → the 2×6 nesting), Phase 2 the size question (our Si supercell curve; then CP2K's 32-atom MnO), **Phase 2.5 the SOLID/OOD cleanup campaign** (user: implementation detail keeps creeping into the abstract faces — and a new capability must land AFTER it), Phase 3 DFT+U oracle-first (✅ CP2K HAS `&DFT_PLUS_U`, so the oracle is already installed and validated), Phase 4 a second code only when Step 5 needs one. | ✅ **1.1 IS DONE (`22b7f0b9`), and it REFUTED the plan's own guess**: the SCF's linear algebra is ~5 s of a 393 s serial MnO run — the DIAGONALISATION is **0.038 s total**.  What the instrument found instead is **`setup: hamiltonian ctor` = 23.9 s**, exclusive of the Becke mesh build and the Φ tables, running ONCE and not threading = 20% of the 12-thread wall.  Unbucketed 49 s → 25 s. ✅ **1.1(a) IS DONE (2026-09-06) AND IT REFUTED ITS OWN PREMISE TOO**: the residual is NOT in `Iterate` — that bucket charges **0.33 s**, `Converge` 0.0009 s.  The 25 s is the **SECOND Hamiltonian**: `SolidCalculation::BuildStage` rebuilds the whole thing per ANNEAL STAGE and had no bucket, so 1.1's "it runs once" was wrong.  Ledger now `[×2, 23.18 s/call]`, **Hamiltonian construction = 68.4 s = 56% of a 121 s threaded run** (121.0 s = 39% of the 313.9 s unset arm), setup 69.5 s vs SCF 51.5 s, and **unbucketed 25 s → 0.03 s** (39 buckets sum to 120.98 of 121.01 s; 313.71 of 313.76 unset) — the ledger is a PARTITION now, in both arms.  ★ **The ctor's EXCLUSIVE half is the part that does not thread (1.11×)** — its children do (Φ tables 8.7×), which is why it hid. ✅ **1.1(b) IS DONE TOO (2026-09-06) — AND IT MADE THE RUN 1.45× FASTER.**  Six buckets took the 23.75 s/call ctor apart, and it was **two calls to one bad index**: the run folds the same ~97k mesh points TWICE per Hamiltonian (orbit-consistency filter in `UnitCell`, then `FoldMesh` in `GPW_IBS`), at 9.6 + 9.5 s = **38 s of a 121 s run**. ★★★ But the fix was not de-duplication — `TorusIndex` bucketed every mesh on a **constant 64³ grid** (its ctor started at 64 and only ever SHRANK, and the shrink condition never fires at tol 1e-8).  Average occupancy 0.37/bucket looked perfect and was the wrong statistic: an atom-centred RADIAL mesh is clustered, so every query near a nucleus scanned thousands of candidates.  Grid now as fine as the tolerance allows + centre-bucket-first probe ⇒ **both folds 50× (9.5 → 0.19 s/call)**, Hamiltonian ctor 34.4 → 15.5 s/call, **MnO 12-thread wall 120.8 → 83.2 s**, `Etot` bit-identical, 813/813, pinned by a grid-independence unit test.  ⚠ The duplicate fold is now 0.19 s/call — **no longer worth removing**. ✅ **THE BECKE-THREADING FLAG IS SETTLED, AGAINST MY ARM**: it threads **8.4×** (68.74 s/call at `GPW_OMP_THREADS=1` vs 8.15 at 12), §7c stands, §5e closed — and the real finding is a PROTOCOL DEFECT: **`GPW_OMP_THREADS` unset is NOT a serial arm** (the Becke build is parallel by DEFAULT and reads the variable only as a thread CAP), so a serial row must set `=1` and check for 99% CPU. ✅ **R2.22 IS DONE (`0210cfb9`, 2026-09-06): MnO 12-thread wall 83.2 → 67.6 s (1.23×), `Etot` bit-identical, 814/814.**  The iterator no longer deletes what it was handed, so the facade keeps ONE Hamiltonian for the whole schedule and the ledger's `hamiltonian ctor` / `becke mesh build` run once (the surviving `[x2]` buckets are the per-stage SCF residues, correctly).  The accelerator is still rebuilt per stage — stale Pulay/DIIS + the type changes — which was never the cost.  A pre-existing LEAK on RunGpw's non-annealed path fell out of it.  Then 1.2 (BLAS-mode serial arm) and 1.3 (the 2×6 nesting). ⚠ Source builds are the default for comparison codes; flang-21 at `/opt/LLVM-21.1.6-Linux-X64/bin`, gfortran + OpenMPI already present.  ⛔ **1.2 IS CLOSED AS A NO-OP**: `QCHEM_BLAZE_BLAS` has been **default ON since 2026-08-15** (cache ON, `BLAZE_BLAS_MODE=1`, `libopenblas.so.0` linked, `PinBlasToOneThread` in `gtestmain`) — the plan had read a STALE `CMakeLists.txt` comment sitting fifteen lines above the `option(... ON)` that contradicted it (comment fixed).  ⇒ Every campaign number was already through a pinned OpenBLAS zgemm, so §7c's 1.21× H_xc is a **WIDTH** problem, not a kernel one, and 1.3 stops being conditional. ▶ **NEXT = 1.3 (the 2×6) TOGETHER WITH \f$V_H\f$.**  MnO Γ is now **4.30×** (290.1 s serial / 67.5 s at 12 threads) — two thirds of the 3.08→5 criterion closed by REMOVING serial work, not by threading it.  ⚠ 1.3 alone reaches only **4.82×** (perfect 2×6 on the 11.15 s H_xc bucket), so the \f$V_H\f$ `SymmetrizeGMap` walk (8.4 s, threads at **1.00×**) must ride along or the increment will look like it missed the target.  Cheapest thing available (but not the next thing): the angular sets waste one NNLS restart at production L (`pool=82 FAIL → pool=164 ok`), ~30% of 3.7 s. | `doc/ParallelAndOraclePlan.md` |
| **7** | **Continuous cleanup** — ⚠ **now scheduled as PHASE 2.5 of `doc/ParallelAndOraclePlan.md`**, not a rhythm: the user wants a campaign before DFT+U lands on these faces | `doc/CleanupCandidates.md` R1/R2 + **V1 (the interface-design questions)**, item 6's `FIT_SF_Ortho` metric split, the `IsPolarized()`/`IsRelativistic()` identity smell, the NEW grid×fit-basis audit, the `dynamic_cast` survey, plus **V1.32** (de-template the finite `IrrepCD` leaf). | `doc/ParallelAndOraclePlan.md` Phase 2.5 |
| **TE** | ★★ **THE TEST SUITE — ORGANIZATION FIRST, THEN COST (user, 2026-09-08).**  Two complaints, and the ORGANIZATION one is primary: *"for the test review I am also concerned about organization."*  ▶ **THE SHAPE THE USER ASKED FOR — an SCF test is a POINT IN A PRODUCT SPACE, so name and file it as one:** `{basis: PW, GPW, LAPW, …} × {material: Si, NaF, MnO, Na, Al, …} × {real-space grid: Uniform, Becke} × {k: Γ, multi-k} × {symmetry imposed: yes, no} × {kT: 0, anneal}` (etc. — the axis list is the user's, and it is open-ended by design).  Then: pick the FILE BREAKDOWN (by basis set?), lay the chosen permutations out in a CONSISTENT ORDER, and give them a CONSISTENT NAMING CONVENTION so a reader can see which cells of the product are covered and which are holes.  ⚠ Today the naming is ad-hoc (`SR_2x2x2ShiftedMP_vs_CP2K`, `RealTRIMBlocksWithMOMMatchComplex_SiMixedMesh`, `BeckeXC_IBZ_SiDiamond`) — each name is individually reasonable and the SET is unreadable, which is exactly why nobody can answer "do we need them all?".  ★ The coverage question is a CONSEQUENCE of the layout, not a prerequisite for it: once the permutations are in a grid, the duplicates and the holes are both visible. **AND THE COST, MEASURED 2026-09-08** (`ctest -j8`, 870 tests, 937 CPU-s, 207 s wall): **`GPW_SCF` is 611 s = 65% of all test CPU from 35 enabled tests**, and **`GPW_SCF.PolarizedRunKeepsItsSpin` alone is 251 s = 27% of the suite — it sets the -j8 wall floor by itself** (Mn atom, 16-bohr box, 12 SCF iterations on a Becke mesh).  Also: **27 of `GPW_SCF`'s 62 tests are DISABLED** — hand-run ladders, sweeps and probes, i.e. INSTRUMENTS, not tests, and 4800 lines of `GPW_SCF_UT.C` that ctest never touches. | ▶ **THREE ACTIONS, in order.**  (a) **Define the axes and the naming convention on paper first**, then re-file — the layout is the deliverable, not the deletions.  (b) **`PolarizedRunKeepsItsSpin`: its claim is a MIXER property** — ρ̃ (Kerker) mixing on a polarized density must not collapse the spin — *not* an SCF property. It belongs in `src/ChargeDensity/tests` driven by a hand-built polarized density, with no SCF at all; that is the user rule (unit for the dev loop, integration for acceptance) and it roughly halves the suite wall time. ⚠ While there: its order probe samples m(r) at 0.7 bohr off the nucleus, which is a spin DENSITY, not a moment — the Becke site blocks now provably support a real INTEGRATED site moment (`BeckeMesh.*`, 2026-09-08), so the probe should become an integrated one and gets cheaper doing it.  (c) **TEST libcint-SPHERICAL (S3b)** — the one remainder of the retired `doc/OldPlans/SphericalSALCPlan.md`, guarded out today; the in-house spherical SALC is shippable and this is the last arm (user, 2026-09-08: *"SphericalSALCPlan.md can [be] retired, just add 'test libCint' into stage C"*).  (d) **Decide the 27 disabled tests as a CLASS**: promote to a `CLIapps/` probe binary, or delete the ones whose verdict is already banked in `doc/`. | this row + `doc/ParallelAndOraclePlan.md` PHASE 2.5 |
| **8** | **Step 6 — the 136-function span**: a capability gap, no longer a blocker | Time-boxed research only. Do not let it grow into a track. | Step 6 |

### ★★★ THE GAP-CLOSE PRIORITY ORDER (user, 2026-08-28) — and what goes in which bin

> *"Right now I think these gap close priorities make sense: 1) Per iteration CPU time, 2) Init (pre
> iteration) time, 3) top RAM usage (mostly solved I think), 4) very roughly match the # of iterations.
> I think items that fall into bin 4 can just be documented in OpenWork.md."*

| bin | the axis | where it stands (MnO AFM-II VA, 2026-08-28) |
|---|---|---|
| **1** | **per-iteration CPU** | ★★★ **THE TABLE IS FILLED — `doc/Benchmark.md` §5a, 2026-09-05**, nine rows, BOTH codes pinned serial and measured at 97–99% CPU, with setup split out so the per-iteration figure is SCF-only.  ★ Re-taken after §5f's lever A.  **Ahead of CP2K on seven of nine**: MnO defaults **0.82×**, MnO FM **0.77×**, `QCHEM_BECKE_XC=0` **0.45×**, Si Γ **0.14×**, Si 8 k **0.66× / 0.56×**, NaF full-SR **0.18×**.  Losses: NaF SR2 **1.45×** (a Becke cost) and `CP2K_COMPAT=1` **1.72×** (was 2.05×).  ⇒ **bin 1 is NOT closed, and it is now ONE row**: the parity row, which is also the only one with our accelerations off.  ⚠ Earlier cuts of this line (2.14×/1.47×/1.99×, and a RETRACTED 1.00×/1.07×) divided TOTAL CPU by iterations, i.e. they charged setup to bin 1 |
| **2** | **init / pre-iteration time** | ⛔ **PROMOTED — MEASURED SERIALLY 2026-09-05: MnO's setup is 184.3 s = 47% of the default run against CP2K's 8.1 s (23×)**, of which **136.6 s is TWO Becke mesh builds** (68.3 s each) + 43.1 s XC-mesh Φ tables.  ✅ With `QCHEM_BECKE_XC=0` our setup is **1.76 s and BEATS CP2K's 8.1 s**.  ⇒ bin 2 is a Becke-mesh question, exclusively — and on the default MnO route it is now worth more than everything left in bin 1.  ⚠ The old "16.7 s" mesh-build figure was a THREADED ledger bucket; the build is `#pragma omp parallel for` (UnitCell.C:273) |
| **3** | **peak RAM** | ✅ **solved, and we WIN**: 1323 → 476 MB default; **113–132 MB on the parity routes against CP2K's 217 MB**, re-confirmed 09-05 |
| **4** | **iteration count** | 31 (default) / 93-and-capped (parity) against CP2K's 44 — ⇒ DOCUMENT, do not chase.  ⚠ Removing the gather's D-screen (`a7561e92`) moved the SMALL rows' counts (Si Γ 11 → 17, Si 2×2×2 Γ-centred 7 → 16, shifted MP 16 → 14, NaF SR2 29 → 23) at unchanged \f$E_{tot}\f$; the MnO counts did not move.  The findings are below |

### ★★★ BIN 1's REMAINING GAP IS IN THE HAMILTONIAN, NOT THE KERNEL (2026-09-04)

The 09-04 ledger shows **~8.6 KS-field integrations per SCF iteration** on a 2-channel system where the
physics needs **3** (one \f$V_H\f$ + one \f$V_{xc}\f$ per spin).  `GPW_INTEGRATE_CENSUS=1` reports every
miss as a genuinely NEW field, so this is **not** redundancy a memo can remove — the code really is
integrating that many distinct potentials.  Traced to two independent causes, both structural:

**(1) THE HARTREE MATRIX IS BUILT ONCE PER SPIN AND IS SPIN-INDEPENDENT.**  `Vee_Hartree::MakeMatrixT`
takes `const Spin&` — *unnamed*, i.e. provably unused — so \f$\langle i|V_H|j\rangle\f$ is identical for up
and down.  But `tDynamic_HT_Imp::GetMatrix` keys its cache on `Irrep qns(bs->GetIrrep(s))`, which DOES vary
with spin, so the two channels miss each other and the same matrix is gathered twice.
⇒ ~2.1 of the 4.2 "h ball" calls per iteration are an exact duplicate.
✅ **A fix here is BIT-IDENTICAL** — the same matrix, computed once instead of twice.

**(2) EXCHANGE AND CORRELATION ARE SEPARATE TERMS, EACH PAYING ITS OWN GATHER.**
`Vxc_QuadraturePol::MakeMatrixT` builds \f$v_x\f$ and calls `itsQuad->Matrix(bs,v)`;
`Vcorr_QuadraturePol::MakeMatrixT` builds \f$v_c\f$ and calls it again — same basis, same spin, same grid.
The gather is LINEAR in the field, so
\f$\langle i|v_x|j\rangle+\langle i|v_c|j\rangle=\langle i|(v_x+v_c)|j\rangle\f$: summing the two
POTENTIALS pointwise and gathering once is mathematically identical and halves the XC gathers.
⇒ 4.4 "h raw" calls per iteration where 2.2 would do (matches the prediction 2 terms × 2 spins = 4).
⚠ **A fix here is NOT bit-identical** — it changes the summation order from (matrix + matrix) to
(matrix of the field sum), so it needs an anchor re-bank.

**WHAT IT WAS WORTH — BOTH FIXED, MEASURED 2026-09-04.**  Wall **2:11.2 → 1:31.4 (−30.6%)**; gather
misses 65 → 43, its bucket 92.4 → 56.8 s.  ⚠ The "1.54× → 1.07× CP2K" this line originally carried is
WITHDRAWN — that run used the default basis and one SCF stage, not the VA/annealed recipe CP2K's 8.5 s/step
belongs to.  The −30.6% is a same-config A/B and stands.
\f$E_{tot}\f$ **identical** to all 10 printed figures and the full suite green — no anchor re-banked.
⚠ **THE SPLIT WAS NOT EVEN**: (2) carried essentially all of it; (1) was worth ~0.3%, not the ~23% first
estimated, because its duplicates were ALREADY GatherMemo hits.  ⇒ **Read gather MISSES, not closure
calls** — the ledger prints both and only the first is work.  A cost estimate taken off call counts alone
is wrong whenever a memo sits underneath.

⚠ **THE DESIGN QUESTION THIS RAISES**, because the term architecture is otherwise good (open/closed: add a
term, change nothing else): the inefficiency is that **the SUM OF MATRICES could be the MATRIX OF THE SUM**
whenever terms share a quadrature and a grid.  Expressing that without destroying the term seam is the
actual work — a term would contribute its FIELD to a shared quadrature rather than its finished matrix, and
the quadrature would gather once.  ⇒ Do not paper over it with another cache; the caches are already
correct and are catching everything catchable.

### ★★★ THE k-SCALING GAP, AND THE CROSS-k GATHER MEMO THAT WOULD CLOSE IT (2026-09-04)

**THE GAP.** Per-ITERATION cost, Si SR, serial, against CP2K's own per-step:

| | qchem s/iter | CP2K s/iter | ratio |
|---|---|---|---|
| Γ | **0.063** | 0.417 | **0.15×** (we are 6.7× FASTER) |
| 2×2×2 Γ-centred (8 k) | 0.434 | 0.431 | 1.01× |
| 2×2×2 shifted MP | 0.411 | 0.429 | 0.96× |

⇒ **our per-iteration cost rises 6.9× from Γ to 8 k-points; CP2K's rises 1.03×.**  We start 6.7× ahead and
spend the whole lead on k-points.  (CP2K also folds **8 → 4 by TIME REVERSAL** on the shifted mesh —
confirmed in `bench_Si_222shift_cp2k.log`, with point-group symmetrization OFF — while we run all 8.  That
is a separate 2× on that row.)

**THE CAUSE, and it is a memo bypass.**  The per-offset reductions \f$B_{ij}(n)\f$ are k-INDEPENDENT (the
Bloch phase enters only the final contraction), and `IntegrateMemo` exists to share them.  But
`memoize = (screenD==nullptr && !sf)` bypasses it whenever the caller passes a density screen — which the
per-iteration KS path always does.  So every k-block redoes the whole real-space sweep.
`GPW_INTEGRATE_CENSUS=1` on Si 2×2×2: **32 of 51 gathers are the SAME FIELD** blocked only by the screen.
⚠ The census labels that case *"screen widened"*, but it compares screen HASHES — it establishes
*different*, not *wider*.  Across k-blocks the screens are merely different (each block's `Dscr` is the
union over ITS channels, built from `D(k)`), so neither covers the other.

**AN ATTEMPT THAT WORKED AND WAS STILL REVERTED (measured, then backed out).**  Withholding `screenD`
whenever the screener ignores weights (which `CP2K_COMPAT=1` guarantees, since it selects
`GeometryOnlyScreener`) unlocks the memo with no new machinery:

| Si 2×2×2, `CP2K_COMPAT=1` | gather/call | iterations | total CPU |
|---|---|---|---|
| Γ-centred, before | 0.04057 s | 7 | **4.08 s** |
| Γ-centred, after | **0.01467 s (2.77×)** | **16** | 6.59 s ⛔ |
| shifted MP, before | 0.0429 s | 16 | 7.59 s |
| shifted MP, after | **0.02187 s (1.96×)** | **14** | **5.30 s ✅** |

\f$E_{tot}\f$ identical on both rows.  ⇒ **The per-call win is robust (1.96–2.77×); the NET is a coin flip**,
because withholding the screen also changes the SEED Fock: the SAD seed density is DIAGONAL, so with a
screen the first Fock is diagonal-only and without one it is a full sweep.  That is the same coupling the
2026-08-18 cross-run pollution note describes.  Perturbing trajectories for a reason unrelated to the
optimisation is not acceptable, so it was reverted.

### ⛔ ATTEMPT 2 (2026-09-04): FIX THE DIAGONAL-SEED DEFECT FIRST — REFUTED BY THE STREAM FOLD

**THE DEFECT IS REAL, and naming it is the durable part** (user: *"the Screener system needs to properly
handle diagonal seed densities"*).  \f$h_{ij}=\langle\chi_i|V|\chi_j\rangle\f$ is NOT weighted by
\f$D_{ij}\f$.  Dropping a term because \f$D_{ij}=0\f$ is sound for the ENERGY alone — \f$\mathrm{Tr}(Dh)\f$
is blind to it — but \f$h\f$ is DIAGONALIZED to make the next density, so zeroing \f$h_{ij}\f$ wherever the
density vanishes is a SELF-FULFILLING truncation: a pair with no density can never acquire any.  With the
SAD seed, whose \f$D\f$ is DIAGONAL, that is every off-diagonal element.  The tree had recorded the symptom
(the cross-run note's *"a fresh process's ... is diagonal-only"*) without recognising it as a defect.
⚠ ROOT CAUSE is an overload in the screener interface: `cij==0` means BOTH "structurally absent"
(fold-dead) AND "zero density".  The first must be excluded; the second is NO INFORMATION and deserves the
floor.

**THE FIX WORKED ON BIN 1 AND WAS STILL REVERTED.**  Substituting a unit weight for a vanishing one, plus
withholding `screenD` when the rule ignores weights (which is then bit-identical), gave:

| Si 2×2×2, `CP2K_COMPAT=1` | gather/call | **CPU/iteration** | vs CP2K 0.43 s |
|---|---|---|---|
| Γ-centred | 0.0407 → **0.0146 (2.79×)** | 0.587 → **0.414 s** | 1.36× → **0.96×** |
| shifted MP | 0.0427 → **0.0220 (1.94×)** | 0.473 → **0.379 s** | 1.10× → **0.88×** |

\f$E_{tot}\f$ identical on both; the iteration-count movement (7→16, 16→14) is bin 4 and was accepted.

⛔ **BUT IT BREAKS THE STREAM FOLD** — `GPW.StreamFoldReducedMatchesFull_{DimerInBox,SiDiamond_HalfK}`,
gate `dHS` at **0.14 against a 2.4e-8 tolerance**.  The reduced arm screens on `FoldScreenMax`, the ORBIT
MAX, which exists precisely so every orbit member truncates IDENTICALLY; the full arm uses each member's
own \f$|D_{ij}|\f$.  A per-term substitution therefore gives a member whose own value is 0 the FLOOR in the
full arm while the reduced arm derives it from the representative at the orbit-max tolerance — different
boxes, different \f$h\f$.  ⇒ **"no information → floor" destroys the orbit-invariance the fold rests on.**
⚠ And hoisting the substitution BEFORE the fold does not rescue it: \f$\max(1.0,\,|D|<1)=1.0\f$ would
dominate every orbit and collapse D-aware screening to the floor everywhere.

▶ **WHERE THIS POINTS.**  The two symptoms have one source: **we D-screen the gather at all.**  CP2K does
not; `GeometryOnlyScreener` does not; the \f$\varepsilon/|c_{ij}|\f$ widening is an ENERGY-accuracy
argument and \f$h\f$ has a second consumer.  Dropping the density screen from the GATHER (keeping it on the
collocation, where the weight really is the scatter weight) would fix the diagonal seed, restore
orbit-invariance for free, and unlock the cross-k memo — at a cost that has never been measured on the
gather alone.  ⇒ That measurement is the next step, and it is a `RunPolicy`-level decision, not a local one.

▶ **THE ALTERNATIVE THAT PRESERVES TRAJECTORIES EXACTLY.**  \f$B_{ij}(n)\f$ does not depend on the screen at all —
the screen only decides which terms are COMPUTED.  So: memoize \f$B\f$ over the UNION of the active sets
seen, and have each caller contract only ITS OWN active set out of `PairB::nb` (it can test its own screen
per (pair, offset) in O(1) — that is just the existing pre-filter).  Replay is then valid whenever the
memo's set COVERS the caller's, each k-block gets exactly the terms it would have computed, and no
trajectory moves.  Cost: the stored sweep is the union, wider than any single block's.
⚠ And raise `kMaxIntegrateMemos` (currently **4**) with it: 8 k-blocks × 2–3 distinct fields per iteration
will thrash a depth-4 cache even once the key is right — the same shape as the `CollocMemo` depth-1 → 5 fix.

### ⚠ BIN 4 — THE ITERATION COUNTS ARE NOT COMPARABLE YET, AND HERE IS EXACTLY WHY (2026-08-28)

⛔ **THE TWO CODES DO NOT RUN THE SAME ρ-MIXING ALGORITHM** (checked on the user's instruction against
`IntegrationTests/CP2K/mno_afm2_gpw_va.inp`, which is the deck every MnO row is compared to):

| | CP2K's MnO deck | qchem's MnO recipe |
|---|---|---|
| density mixing | `BROYDEN_MIXING`, `ALPHA 0.2`, **`BETA 1.5`**, `NBUFFER 8` | `Kerker(G0=1.0)` α=0.45, **no history** |
| orbital step | `&DIAGONALIZATION STANDARD` | `Ladder` → **`GDM`** (direct minimisation) |
| smearing | Fermi-Dirac, kT 5.0e-3 au | Fermi-Dirac, kT 5e-3 ✅ **matches** |
| convergence | `EPS_SCF 1.0E-6` | `MinΔρ 1e-5` |
| cap | `MAX_SCF 200` + `IGNORE_CONVERGENCE_FAILURE` | `NMaxIter 80` |
| extra virtuals | `ADDED_MOS 20` | the full basis (N=118) |

Three things follow, and none of them is "qchem converges worse":
1. **CP2K composes Kerker WITH Broyden** — `BETA` inside `&MIXING` IS the Kerker damping denominator, so
   its recipe is *Kerker-damped 8-vector quasi-Newton*.  ✅ **AND SO IS OURS** — `MakeGSpaceMixer` passes
   `kerkerG0` straight into `PulayMixer`, and the run confirms it: `[Pulay] ENABLED (↑): G0=1`.  So
   `MNO_PULAY=8 MNO_ALPHA=0.2` IS the same CLASS of mixer as CP2K's deck (Kerker-damped history), and the
   remaining differences are the quasi-Newton flavour (Pulay/DIIS vs Broyden) and the Kerker parameter
   (our `G0=1` against their `BETA=1.5`).
   ⚠ **The SCF banner is misleading and that is worth fixing**: it prints
   `mixer: Pulay(depth 8, start 5) alpha=0.2` and never mentions Kerker, so a reader of a benchmark row
   would conclude the preconditioner had been swapped out.  The `[Pulay] ENABLED … G0=` line has the truth;
   the banner — which is the line `doc/Benchmark.md` tells people to copy — does not.
2. **The orbital steps are different classes.**  CP2K diagonalises and mixes the density; qchem's stage 2
   is GDM, a direct minimiser that does NOT mix.  An iteration means a different amount of work and a
   different amount of progress on the two sides — which is bin 4's whole point.
3. **The convergence measures are different quantities**, so "44 steps vs 31" compares two thresholds as
   much as two trajectories.

### ★★★ CP2K DOES **NOT** DO D-AWARE SCREENING — checked in the source, 2026-08-28

> *"The D-aware tolerance (ε/|c_ij|) has caused a number of problems, including very subtle bugs that were
> hard to track down.  Does CP2K also do D-aware screening?"* (user)

**No.**  From `/home/janr/Code/cp2k/src/task_list_methods.F`, the production task-list path:

```fortran
CALL compute_pgf_properties(cube_center, lb_cube, ub_cube, radius, ...,
                            ra, rab, rab2, dft_control%qs_control%eps_rho_rspace)
...
SUBROUTINE compute_pgf_properties(..., la_max, zeta, la_min, lb_max, zetb, lb_min, ra, rab, rab2, eps)
   cutoff    = 1.0_dp
   prefactor = EXP(-zeta*f*rab2)                     ! the GAUSSIAN PRODUCT prefactor E_ij -- geometry
   radius    = exp_radius_very_extended(..., eps=eps, prefactor=prefactor, cutoff=cutoff)
```

Three things settle it:
1. **The subroutine takes no density argument at all** — its whole parameter list is angular momenta,
   exponents, centres and one `eps`.
2. **`eps` is `eps_rho_rspace`, a single GLOBAL threshold.**  Not \f$\varepsilon/|c_{ij}|\f$.
3. **`radius_list[ntasks]` is fixed data on the task list**, built once, while `pab_blocks` is passed
   SEPARATELY to `grid_collocate_task_list` per call.  CP2K keeps the GEOMETRY and the DENSITY strictly
   apart; we mix them.

What CP2K screens on the density side is DBCSR **block sparsity** — which blocks exist at all — a
structural question, not a per-task radius that breathes with \f$|D|\f$.

⇒ **THIS IS THE MOST INTERESTING PARITY ITEM LEFT, because three separate threads converge on it:**
- **the bug record** — the D-aware screen is behind the shifted-MP `Re[]` defect (4.1 Ha), the union-vs-
  per-component box design, `FoldScreenMax` (a 2.3 Ha collapse when the orbit reduction used the signed
  average), and the reduced-vs-full tier split in the T3 gates.  It is the single most defect-dense idea
  in the collocation path;
- **CP2K parity** — it is a deviation we never declared, and it is not on the `CP2K_COMPAT` table;
- **BIN 1** — and this is the part that was not obvious.  \f$\varepsilon_{eff}\f$ is the ONLY
  \f$D\f$-dependent input to `BoxGeom` (user, same day: the chord quadratic *"does not depend on rho(r)"*
  — correct except for \f$R_q\f$, which comes from \f$\varepsilon_{eff}\f$).  **Remove it and the ENTIRE
  box geometry becomes iteration-invariant** — `hw`, `Rq`, the chord bounds, all of it — so the per-line
  `sqrt` that was just measured at ~40% of the kernel and called "close to irreducible" moves OUT of the
  SCF loop and into the task list, where it is paid once.

⚠ **AND IT IS NOT FREE — the trade needs measuring, not assuming.**  Dropping it widens every box
(\f$\varepsilon_{eff}=\varepsilon/|c|\ge\varepsilon\f$, so the D-aware box is always the smaller one)
and, more importantly, retires the per-(pair, offset) KILL `pf >= -log(e)`, which removes whole terms
rather than shrinking them.  Reach scales as \f$\sqrt{\ln}\f$ and work as its cube: a pair at
\f$|c|=10^{-2}\f$ loses ~1.4× on volume, at \f$10^{-4}\f$ ~2.1×.  ⇒ **The experiment is one knob
(`eps_eff := kDensityEps` unconditionally) and one bin-1 probe run**, and it is worth doing before any
more per-line micro-optimisation — it either buys the sqrt back for free or prices the D-aware screen
honestly for the first time.

### ⛔ AND MATCHING THE MIXER DID NOT HELP — MEASURED, 2026-08-28

Since our Pulay composes with Kerker exactly as CP2K's Broyden does, `MNO_ALPHA=0.2 MNO_PULAY=8` under
`CP2K_COMPAT=1` puts qchem in the same mixer CLASS as the deck.  Run it:

| MnO AFM-II VA, `CP2K_COMPAT=1` | qchem default mixer (α=0.45, no history) | **CP2K-matched (α=0.2, Pulay 8)** |
|---|---|---|
| stage 1 (kT=5e-3) | UNSETTLED at **13** iters, −61.41070717 | **OSCILLATING at 80**, −61.40861546, Eamp 0.53 |
| stage 2 (kT=0) | FIT-FLOOR STALL at 80, **−61.39789688**, Δρ 1.69e-5 | FIT-FLOOR STALL at 80, **−61.39789768**, Δρ 8.97e-5 |
| the AFM order | survived | survived |
| CPU / wall / RSS | 2736 s / 45m40s / 112 MB | 3654 s / 1h01m / 242 MB |

**Three findings, and the second is the useful one:**
1. ⛔ **The CP2K-matched mixer is WORSE here, not better.**  α=0.2 with an 8-vector history OSCILLATES for
   80 iterations in stage 1 where plain α=0.45 Kerker settles in 13.  Matching CP2K's mixer parameters is
   therefore NOT a route to matching its iteration count on this system.
2. ★★★ **BOTH MIXERS FLOOR AT THE SAME PLACE: −61.39789688 and −61.39789768 — 8e-7 Ha apart, from two
   different mixing algorithms.**  So the stall is a property of the MAP, not of either mixer, and no
   amount of mixer tuning will move it.
   ⛔ **BUT DO NOT READ THAT AS "THE UNIFORM GRID CANNOT CONVERGE" — THE CONTROL SAYS OTHERWISE.**  The
   detector prints "functional/grid" and an earlier cut of this section took it at face value.  It is
   wrong, and the run that refutes it was taken the same day: **`QCHEM_BECKE_XC=0` with the imposition
   KEPT — same uniform XC mesh, same pair route — CONVERGED in 25 iterations at Δρ 3.47e-6** (−61.40358773).
   The uniform mesh is therefore perfectly capable of converging this system.  What the two stalled runs
   have that the converged one does not is **the imposition switched OFF**, which is exactly the banked
   2026-08-26 finding: *"the imposed star-average was buying CONVERGENCE, not accuracy and not the magnetic
   basin"*.  ⇒ The parity stall belongs to the FREE run, not to the grid, and the 5.1 mHa is how far a
   stalled run sits from an answer it never reached — NOT a physical grid bias.
   ⚠ The grid's genuine ACCURACY cost is a separate, also-measured number: **6.1e-4 Ha** between the
   converged Becke and converged uniform imposed runs (−61.40297551 vs −61.40358773).
3. ⚠ **And the two runs FLOOR at Δρ values 5× apart (1.69e-5 vs 8.97e-5) while agreeing on E to 8e-7.**
   Δρ is measuring the mixer's own step, not the distance to the answer — the sharpest single argument for
   A4 that this session has produced.

★ **AND A4 HAS A CONCRETE, SMALL FINDING SITTING IN IT: `MinΔρ` IS COMPARED AGAINST THREE DIFFERENT
QUANTITIES.**  One threshold, three scales:

| site | what it returns |
|---|---|
| `LinearMixer::Mix` | `GetChangeFrom(old)/GetTotalCharge()` — **normalised**, "relative MaxAbs change" |
| `LinearMixer::ReDampMix` | `GetChangeFrom(old)` — **un-normalised**, and the code says so ("matches the legacy re-damp") |
| `DirectMinDriver::Step` (GDM) | `GetChangeFrom/GetTotalCharge` — normalised |
| `KerkerMixer::MixField` | \f$\lVert\tilde\rho_{out}-\tilde\rho_{in}\rVert_\infty\f$ over **G-space coefficients** — a different norm of a different object |

⇒ A fixed `MinΔρ` therefore means something different per RECIPE and per SYSTEM SIZE, which is exactly
what A4 (the Δρ/N gate, `doc/SCFStrategyPlan.md`) exists to fix — and it is now the **fourth** thing this
session has pointed at A4, after the Na2 gate's chaotic α-dependence, the density-degenerate Na2 state, and
the parity run's FIT-FLOOR STALL.  ⚠ Note the MnO parity stall is on the GDM path, which IS normalised, so
this finding does not explain that one; it is a separate defect found while looking.

### ✅ ITEM 1 MEASURED, 2026-08-25 — AND THE MECHANISM IS BIGGER THAN THE ITEM

**All rows: MnO AFM-II, VA (N=118), `MNO_IMPOSE=1`, SINGLE THREAD (`OMP_NUM_THREADS=1 GPW_OMP_THREADS=1`),
same binary as the banked row.**  The baseline reproduces `doc/Benchmark.md` to all 10 s.f.
(−61.40297621, 14+17 iterations) — **there was never a regression on the default path**, and every
collapse below is an armed arm.

| arm | Etot | **Eee** | m_stag | verdict |
|---|---|---|---|---|
| **baseline** | **−61.40297621** | **13.480** | 0.646 → 0.667 survived | conv, PASSED |
| `GPW_XC_DM_SOURCE=1` | −45.52875429 | 28.995 | died it 12 | FAILED |
| + `GPW_XC_DM_BOOST=2` | −45.36 (st1) | — | died it 10 | FAILED |
| + `GPW_XC_DM_BOOST=0.5` | −45.52847324 | — | died it 19 | FAILED |
| + `QCHEM_DM_LOWRANK=0` | −45.52875429 | — | died it 12 | FAILED — *identical to lowrank ON* |
| + `MNO_PULAY=8` | −46.32003443 | 30.060 | st1 SURVIVES (0.540), st2 collapses | FAILED |
| baseline + `MNO_PULAY=8` | **−61.40297621** | 13.480 | survived | conv — *so Pulay×GDM is fine* |
| `MNO_KERKER_G0=0.01` (flat filter), **NO FLAG** | **−38.45368688** | **35.098** | **died it 7** | FAILED |
| `MNO_KERKER_G0=0` (→ linear D-mix), **NO FLAG** | −56.38756586 | 13.609 | died it 13 | FAILED |

**★★★ THE MECHANISM — A MONOTONE DOSE-RESPONSE IN Eee.**  The AFM staggering's shortest mode is
\f$|G|=1.24\f$ (\f$f=0.61\f$ at \f$G_0=1\f$); the lowest CHARGE mode is \f$|G|=0.65\f$ (\f$f=0.30\f$) —
the arithmetic already in `GPW_SCF_UT.C`, which notes Kerker "already favours the magnetism 2:1".
Arming the flag replaces \f$\alpha f(G)\f$ on the XC channel with a FLAT \f$\alpha_{eff}=0.33\f$
(\f$\alpha=0.45\f$): the AFM mode goes 0.275 → 0.33 (**1.2× less damped**, barely touched) while the
charge mode goes 0.135 → 0.33 (**2.4× less damped**).  It destroys the selectivity, not the step size.

| low-G charge-mode damping | Eee | m_stag dies | Etot |
|---|---|---|---|
| 0.135 (baseline, f=0.30) | 13.480 | never | −61.403 |
| 0.33 (armed — XC channel only) | 28.995 | it 12 | −45.529 |
| 0.45 (flat Kerker — BOTH channels) | **35.098** | **it 7** | **−38.454** |

Monotone in every column.  ⇒ **KERKER'S LOW-G CHARGE-SLOSH DAMPING IS WHAT HOLDS THIS AFM BASIN.**  The
moment death is a CONSEQUENCE of the charge runaway, not a spin-channel effect.  ⚠ The linear-D-mix arm
is a **different** failure (Eee stays 13.6 — no slosh, moment still dies): same symptom, other mechanism,
**do not conflate them**.

**★ Eee IS NOW A VALIDATED CHARGE-SLOSH DETECTOR.**  Already computed, zero cost, and it ordered all four
failures correctly without having been designed for any of them.  See the grad-student item below.

**THREE CORRECTIONS TO THIS ROW'S OWN NUMBERS, all measured today:**
1. ⛔ The ρ̃ bucket is **15% of this row's CPU (76.2 s of 500.4 s)**, NOT the 51% recorded 2026-08-21 — the
   Becke-ε and Φ-table work shrank the row around it.  **The perf ceiling is ~1.17× on the row, not ~20×**,
   and 20× was never available on a row whose stage 2 is GDM (already D-backed).
2. ✅ The ACCURACY claim is **fully upheld**: matrix-free ρ̃ gives a mean **15.2%** of 97160 points with ρ<0
   (max 18.2%, min ρ −0.154); the DM route gives **0 negatives over 154 samplings**.  Per-sampling the DM
   route is **130×** cheaper (5.077 s → 0.039 s), against the banked ~124×.  *The item's GOAL is right; its
   ROUTE is what fails.*
3. Single-thread CPU is **500 s** against the table's threaded 663 s — the busy-wait inflation the
   benchmark protocol warns about, measured.

**⚠ AND IT BEARS ON STEP 5 — UNTESTED, BUT CHECK BEFORE BUILDING ON THE OFFSET.**  If Kerker's charge-slosh
damping is load-bearing for the AFM basin, the AFM and FM arms are not necessarily converged under equally
safe conditions — and Step 5's ~37 mHa *configuration-selective* FM-favouring bias is measured BETWEEN
those two arms.  Not asserted; worth one control.

**LOW-RANK FACTORING — EXONERATED, AND ONE CLAIM QUALIFIED.**  `QCHEM_DM_LOWRANK=0` vs ON give
**−45.52875429 identically** (10 digits, both stages) across a 160-iteration pathological trajectory: the
factoring is exact and is not the cause.  ⚠ But the code's justification — *"the rank is the same from tol
1e-6 to 1e-12, so the occupied block is cleanly separated … the cut is not a tuning decision"* — is a
**kT=0 statement** (user, 2026-08-25: at kT>0 the gap fills with weak modes and the cutoff gets hard to
pick).  MEASURED pivot spectrum at kT=5e-3: 14 occupied modes at O(0.1–5), a **thermal tail at ~7e-6**,
then roundoff ~2e-13 — a 4–5 decade gap, not 12, and the doc's own "safe" tol=1e-6 sits within ~8× of the
thermal modes.  It is SAFE here only because pstrf uses LAPACK's roundoff floor (not that tol) and because
MnO's gap swamps this kT: rank stayed **13 of 118 in 124 of 154 calls**, and only **4 calls** showed a tail
at all.  A smaller gap or hotter stage would populate it persistently.  Qualify the claim; do not delete it.

---

**CLOSED — reference, not work.**  Step 1 (the head-to-head table STANDS → `doc/Benchmark.md`) · Step 2
(the T3 pair-stream orbit fold is ARMED BY DEFAULT) · Step 4 (RAM largely answered by Step 2 — *do not open
it as a track*) · **the ρ-GEMM low-rank D** (BUILT, guarded, exact, and default-ON via `QCHEM_DM_LOWRANK`;
its own bullet carries the correction that closed it — the bucket was MISIDENTIFIED, the DM GEMM measured
**1.70 s** against the matrix-free ρ̃ sampling's **35.0 s**, which is index item 1, so the 7–8× rank win had
almost nothing to bite on) · and both history files.

**EVIDENCE DOSSIERS — no action inside them.**  *"NEW IDEA — fast evaluation of DM-ρ(r) by FACTORING D"*
and its Q1/Q2/Q3, the tier-0 results, the spectrum finding and the LSP design ruling (~440 lines) are the
worked evidence for **item 2**, not a separate thread.  Likewise the ✅/⛔ bullets inside Step 3 are the
record of what was taken or refuted — five of its twelve items are closed and stay for the reasoning.

---

## ★ ONLY THE CELL IS MATERIAL-SPECIFIC (user, 2026-08-25) — what the SECOND magnetic material needs

Reading the inlined MnO run, the user's observation: it does (1) make the lattice, (2) set
`SolidCalcOptions`, (3) set `SCFParams`, (4) define the order parameter, (5) define the anneal schedule,
(6) define the accelerator ladder — **and only (1) is about MnO.**

Sharpened: it is the SHAPE of 2–6 that is generic; a handful of VALUES are MnO's (the cell, `Nelec` /
`species` / `multiplicity`, the two Mn probe positions, and the tuned numbers — α=0.45 because MnO sloshes
worse than NaF, `orthoTol=1e-3` because cond(S)~7e8, kT=5e-3 for the open d manifold).  Everything else a
second magnetic material would repeat verbatim.

**★ THE HIGHEST-VALUE GENERALISATION IS (4), BECAUSE IT CAN NOW BE DERIVED.**  `m_stag` hardcodes MnO's two
Mn positions and a 0.7 bohr offset.  But the run now holds the magnetic decoration
(`SolidCalcOptions::siteSpins`, \f$\sigma_A=\pm1\f$) AND the integrated site moments (T2's
`Hamiltonian::SiteMoments`, \f$\mu_A\f$) in the same object, so for ANY collinear ordering
\f[ m_{order}=\frac{\sum_A \sigma_A\,\mu_A}{\sum_A|\sigma_A|} \f]
— no probe geometry, no guessed offset, no per-material closure, and it is the INTEGRATED observable rather
than the point sample that misled this campaign for months (Step 0a).  Deriving it would delete step (4)
for every future material and leave the schedule as the only thing a new cell must state.
⇒ Do it when the SECOND magnetic material arrives, not before: one cell is not enough to see which of the
tuned values are really material-specific and which were MnO's accidents.

## ✅ CLOSED — the N1 tiers, T1 through T5  →  `doc/OpenWork_History3.md`

T1 (a non-converged energy is unreturnable), T2 (the positive path is exercised), T3 (the Eee charge-slosh
detector and its measured threshold), T4 (the detectors moved into the library) and T5/N5
(self-description + `CP2K_COMPAT`) all landed 2026-08-26.  The evidence is kept in full in the history
file; the standing summary is index item **N1**, whose ONE remaining piece is the coverage gap.

## ✅ THE ON-THE-FLY BOX WALK — 2.21× ON MnO, DONE 2026-08-26  →  `doc/OpenWork_History3.md`

Power tables, incremental wrap, key/nn hoist, reach-sphere screen, then the chord/interval skip (1.40×)
= 3.11× on MnO's box walk; uncached MnO 573 s/iter and RAM 166 MB, which BEATS CP2K.  The exp recurrence
was TRIED AND REJECTED (1.03×, anisotropic, flips a degenerate SCF basin; branch
`exp-recurrence-experiment`).  Say "box walk", not "kernel".  Full record + the four edits in the
history file.

## ✅ WHY CP2K WAS ~22× FASTER — READ FROM ITS SOURCE  →  `doc/OpenWork_History3.md`

Product-centre re-expansion ⇒ separable exp tables + a tensor contraction, zero transcendentals per point.
**Acted on**: the live plan is `doc/CollocationRewritePlan.md` (COMPLETE — cache deleted, contract kernel,
`template<int LP>`, spin-native XC pair route).  Kept in history as the source reading, not as an action.

## ★★★ THE ANCHOR-MOVING SPRINT — the roster (assembled 2026-08-27)

> **USER, 2026-08-27:** *"the only thing pulling me in the other direction is that we actually have a queue
> of (other) anchor moving improvements.  And the thinking was to tackle them all in one sprint."*

**WHY BATCHING IS RIGHT, stated once so it does not have to be re-argued.**  Every anchor re-bank costs a
full MnO acceptance cycle plus a doc pass, and — the part that matters more — doing them ONE AT A TIME
means each re-bank partially MASKS the next one's effect, because the reference it is measured against has
just moved.  Batch them, re-bank once, and every delta is attributable.

**THE PROBLEM THIS SECTION FIXES:** the concept already existed in the docs ("the bit-moving batch", "the
V1.22/§K class") but the roster did not — it was scattered across `doc/OpenWork.md` and
`doc/CleanupCandidates.md` with no single place to read it off.  That is the thing most likely to make the
sprint cost more than it should.

| # | item | where it is described | state | delta |
|---|---|---|---|---|
| A1 | **the collocation contraction kernel** (`GPW_CONTRACT_CUBE`) | `doc/CollocationRewritePlan.md` steps 5–6 | ✅ **LANDED 2026-08-27 — DEFAULT ON**, with A7 | re-banked, below |
| A2 | **V1.22** — Becke per-representative partition | this file, *Continuous — CLEANUP* | not built | unmeasured; imposed runs only |
| A3 | **§K** | `doc/CleanupCandidates.md` (deferred, user) | not built | unmeasured |
| A4 | **the Δρ/N convergence gate** | `doc/SCFStrategyPlan.md` | not built | unmeasured |
| A5 | **GPW default seed → `IonicSAD`** | `doc/CleanupCandidates.md` ("every pinned GPW anchor re-seeds") | not built | unmeasured; re-seeds EVERY GPW anchor |
| A6 | **`SCFParams::XCCuspDeficit`** — the N4 XC feed | N4 above | flag exists, off | a TRAJECTORY change by its own description |
| A7 | **dropping the pair-stream cache** | `doc/CollocationRewritePlan.md` step 7 | ✅ **LANDED 2026-08-27**, with A1 | re-banked, below |

✅ **A1 + A7 LANDED TOGETHER, 2026-08-27, exactly as this section required.**  The re-measurement that
opened step 7 said the 3.9 GB was buying ~1.1× on the whole MnO run (2.91× on the two box-walk buckets, and
those are only ~58% of it), against 25× the RAM — so the cache went, the (shell pair, offset) task list took
its place, and the contraction kernel became the default because the cache-less walk is ~18× slower and the
two were therefore never separable.  ⇒ **A2–A6 are what is left of the sprint.**  Full record:
`doc/CollocationRewritePlan.md` step 7; re-taken rows: `doc/Benchmark.md`.

⚠ **AND A4 PICKED UP A SECOND MOTIVATION ON THE WAY.**  `ImposedOrderLostIsAPostconditionFailure_Na2Box`
went red, and sweeping it found NO stable mixing step: α = 0.65 ✔, 0.7 ✔, 0.75 ✘, 0.8 ✔ on the contracted
kernel.  In every arm — converged or capped — the moment dies and E lands on the same −0.332045 state; what
oscillates is one slow, nearly-degenerate DENSITY mode that E cannot see, and the run's own fingerprint
already calls it *"DENSITY-DEGENERATE (E settled, ρ rotates — benign)"* while `DidConverge()` calls it a
failure.  **A convergence criterion that disagrees with the run's own diagnosis is the thing A4 fixes.**
The gate now sweeps its mixing step instead of betting on one draw (four measured-good values, first that
converges wins), which is a stopgap, not the fix.

⚠ **A5 is the one to sequence FIRST or LAST, not in the middle**: it re-seeds every GPW anchor, so anything
measured against a pre-A5 reference has to be re-measured after it.

### A1 — the contraction kernel: its deltas, measured BEFORE it landed
Kept as the record of what the re-bank was told to expect (the re-taken rows are in `doc/Benchmark.md`):

| system | before | after | shift |
|---|---|---|---|
| Si Γ | −7.115067662 | −7.115067844 | −1.8e-7 |
| NaF Γ | −24.5468825477 | −24.5468834873 | −9.4e-7 |
| MnO AFM-II Γ | −59.69580383 | −59.69580301 | +8.2e-7 |

★ **AND THE SHIFT IS THE WALK'S ERROR, NOT THE KERNEL'S.**  Against a naive exact reference the contraction
tracks to 2.4e-14 relative and the walk to ~3e-14; and the walk additionally applies a per-component screen
the contracted cube has no analogue for, whose integrated cost measures 1.5e-7 absolute (5.8e-5 relative on
a pair whose weights span four decades).  ⇒ Re-banking to the contracted values is moving the anchors
TOWARD the truth, not away from it.  Speed: **4.98× MnO / 3.96× Si / 3.61× NaF** on the box walk.

⚠ **AND ONE THING TO DO WHEN IT LANDS** (the standing provenance rule): the run banner must state WHICH
collocation kernel produced a row, or a future `doc/Benchmark.md` entry is not reproducible.  It is NOT a
`CP2K_COMPAT` deviation — CP2K collocates exactly this way — so it does not belong on that table.

### ⚠ WHAT ROTS NOW THAT A1 IS THE DEFAULT — the polarity has REVERSED
- **The KERNEL cannot rot**: `src/BasisSet/Molecule/tests/M_PG_BoxWalk.C` calls `MakePairPoly` /
  `ContractCube` / `GatherCube` DIRECTLY, and it is now also the production route.
- **THE REFERENCE WALK is the one that can rot**: `GPW_CONTRACT_CUBE=0` is an investigation opt-out and
  nothing exercises it by default.  ⇒ Run `GPW_CONTRACT_CUBE=0 ctest -j8` alongside the plain sweep at any
  breakpoint that touches the collocation path.  Measured 2026-08-27: **792/792 both ways** (and getting the
  walk arm green is what surfaced its per-component `|v|` screen breaking the collocate/integrate adjoint —
  see the plan's step 7).

### ⚠ THE COVERAGE GAP — now the ONLY thing open under N1
`RunMnO` and the new Na2 gate go through `SolidCalculation`, so those are covered.  Other GPW tests still
construct `SolidSCFIterator` directly (`RunGpw`, `RunGpwAnnealed`), so **none of T1–T5 reaches them** — no
`Outcome`, no detectors, no banner.  `ImposedShubnikovHoldsAFMThroughSCF_Mn2Box`, which is the negative
control for the T2 gate, is one of them: the pair currently proves the detector discriminates only because
a human read both outputs.  Retiring both drivers onto the facade is the rest of
`doc/TestFacadeMigrationPlan.md`, and it is now a mechanical follow-through rather than a design problem.

★ **AND THE GAP RUNS THE OTHER WAY TOO (found 2026-08-26).**  The facade covers the OUTCOME side but not
the REPORTING side: `SolidCalculation` opens no `report::Begin/End` run, so a test that went through it
LOST the timing ledger, the `grids` section and the PEAK RSS line that `GpwReport` gives every other GPW
driver.  MnO is the case that mattered — it is the benchmark's most expensive row and it was the one with
no instrument (`e8339cf2` patches it at the call site).  So the migration has to carry the reporting
bracket ACROSS, not just the `Outcome`: the natural home is `SolidCalculation` bracketing its own run (it
already prints its banner unconditionally, which is the same argument), with `Depth()` deciding whether it
nests inside a driver's existing bracket.  Until then, every new `SolidCalculation` caller silently opts
out of the measurement discipline `doc/Benchmark.md` depends on.

### ⚠ TWO SMALLER LOOSE ENDS LEFT BY THIS SESSION
- **Na2's polarized singlet will not converge from an AFM seed at α=0.3** — Δρ parks at ~1e-2 and
  oscillates forever while E is flat to 1e-9, on DIIS, GDM, Kerker, Pulay and smearing alike; α=0.5
  converges in 66.  The same cell run UNPOLARIZED converges in 30 at α=0.3.  Not chased (the gate only
  needed one converging recipe), but "the two-channel run of a ζ=0 system is much harder to converge than
  its unpolarized twin" is a claim worth either explaining or fixing — it is the sort of thing that will
  cost a magnetic campaign later.
- **`VALENCE_LOWQ_VA` under CARTESIAN d is rank-deficient by 10** and the pivot filter auto-drops to
  exactly SR's 122-function span (same min pivot to six digits).  That is the free controlled experiment
  for the user's hand-trim-beats-auto-trim rule — see the vet-stage trim item under *Continuous — CLEANUP*.

## ★★★ N4 — THE RIGHT TREE: MAKE EVERYTHING ELSE ROBUST WITH \f$V_{xc}[\rho\ge0]\f$

> **USER, 2026-08-25:** *"'GPW_XC_DM_SOURCE does not earn the default' bothers me because it means we are
> forced to feed ρ̃_mix into Vxc … Granted ρ̃_mix is **not** exactly garbage … but I think it is still pretty
> junky for Vxc.  We can improve it with bigger G balls (the N2 sweeps) … but I still think that is barking
> up the wrong tree.  I think the right tree is figure out how to make everything else robust with
> Vxc[non junky ρ>0]."*

**⇒ THE HEADLINE WAS WRONG FOR THE MEASUREMENT.**  *"The flag does not earn the default"* is true of the
flag AS IMPLEMENTED, but it reads as "keep feeding ρ̃_mix to \f$V_{xc}\f$", which is NOT what the data say.
What the dose-response actually convicts is **the MIXER** — a flat \f$\alpha_{eff}\f$ destroying Kerker's
mode selectivity — **not the act of feeding \f$V_{xc}\f$ the exact ρ**.  The ρ≥0 goal is untouched by the
result and remains right (0 negatives in 154 samplings vs 15.2%).

**THE FLAG IS ONE OF TWO ROUTES, AND IT IS THE WORSE ONE.**  `GPW_XC_DM_SOURCE` implements the
**WHOLESALE REPLACEMENT** (route b): discard \f$\tilde\rho_{mix}\f$, hand XC a separately-damped
\f$\rho[D]\f$.  That is what forces a damping choice to exist at all — and the damping choice is precisely
what breaks.  The plan's own algebra had already derived the better form, then shelved it on cost:
\f[ \rho_{XC}(r)=\underbrace{\rho_{mix}(r)}_{\text{HARTREE'S OWN ARRAY}}+\underbrace{\big(\rho[D](r)_{exact}-\rho[D](r)_{BL}\big)}_{\text{the CUSP DEFICIT}} \f]
- **The mode-selectivity failure cannot occur here, by construction.**  The band-limited content XC sees is
  *literally the array Hartree sees* — same Kerker filter, same \f$f(G)\f$, same 2:1 selectivity.  There is
  **no \f$\alpha_{eff}\f$ to pick**, so the entire flat-vs-shaped mismatch is not solved but VACUOUS.
- **The added term needs little or no damping**: the cusp deficit is dominated by core density, which barely
  moves between iterations — accurate well before convergence, exact at it.
- ⛔ **THE COST OBJECTION WAS THE WRONG TIEBREAKER.**  The plan rejected this as *"a net LOSS unless the
  correction's {G} can be truncated"*, comparing GEMM counts.  Measurement changes the comparison: the
  wholesale route is not costlier-or-cheaper, it is **UNSTABLE** (−61.40 → −45.53).  Robustness outranks a
  GEMM.

**WHAT "EVERYTHING ELSE" HAS TO BECOME ROBUST — the concrete list:**
1. **The mixer must precondition CHANNELS, not densities** (→ **N3**): charge gets Kerker, spin gets its own
   policy, and the XC feed stops being a second, independently-damped copy of the density.
2. **The XC feed must not be a separate mixing decision at all** — the cusp-deficit form above.
3. **A collapse must not be able to masquerade as an answer** (→ **N1**, T1–T3), so that landing 1 and 2
   cannot silently trade one failure for another.

★ **N2 NOW SUPPLIES THE EVIDENCE THIS ARGUMENT WAS MISSING (2026-08-25).**  The ρ<0 lobes survive an
EXACT raster unchanged (8× the real-space points at identical {G} moves them 1%), so they are band-limiting
— content the ball cannot represent — and no raster geometry or mixing choice can remove them.  Widening the
ball does remove them but cannot scale (an 8× supercell overflows this box at today's C=2).  **That
eliminates every cheaper route and leaves the cusp deficit as the only one**, which is a stronger position
than the design argument alone.

⚠ **STILL NOT MEASURED — the FORM below is a design argument, not a result.**  Two things need testing before it is
believed: (a) whether \f$\rho_{mix}+\Delta_{cusp}\f$ is actually pointwise non-negative (the negative lobes
sit AT the cusps and \f$\Delta_{cusp}\f$ is exactly the positive sharp content missing there, so large
cancellation is expected but not guaranteed) — the existing `GPW_RHO_NEGATIVE` census answers it directly;
and (b) whether \f$\Delta_{cusp}\f$ is as iteration-static as the core-electron argument claims.

## ★★ N5 — `CP2K_COMPAT`: ONE SWITCH FOR EVERY DEVIATION (user, 2026-08-25)

> *"there will eventually be a number of CP2K deviations so we will need one env bool (like CP2K_COMPAT)
> to trigger off."*

**This is the SAME WORK as index item 2's self-describing banner, and it should be built once.**  A run has
to be able to (a) STATE which deviations it is running and (b) turn them ALL off with one switch.  Today
neither exists, and the list is already longer than anyone tracks by hand:

| deviation | today's default | knob |
|---|---|---|
| `QCHEM_DM_LOWRANK` — factored/low-rank ρ (a singles route CP2K does not run) | **ON** | `=0` |
| `SCFParams::XCCuspDeficit` — the N4 XC feed | off | factory policy |
| the T3 pair-stream orbit fold | **ARMED** | per-row declaration |
| `QCHEM_MIX_RHO_M` — (ρ,m) channel basis instead of (ρ↑,ρ↓) | off | env |
| `RasterPolicy` — `BallOnly` IS CP2K's bet (✅ vindicated, N2) | BallOnly | `GPW_RASTER_POLICY` |
| `cutoffFactor` — C=2 vs CP2K's own grid policy | 2 | `MNO_CUTOFF_FACTOR` |

**THE DESIGN, so it does not become a scattered set of `if (getenv(...))`:** `CP2K_COMPAT=1` should resolve
to a POLICY OBJECT consulted where the choices are MADE (the factories), exactly as `XCCuspDeficit` now is
— never a flag read deep inside a kernel.  That is what makes "CP2K parity is a property of what was
BUILT" true rather than aspirational, and it is what lets the banner print the resolved state instead of
guessing at it.
⇒ **Ordering:** build it WITH the banner (item 2), not before — the banner is what makes the switch
verifiable, and a switch nobody can verify is worse than no switch.
⚠ And a new accelerator is not finished until it is on the table above (`doc/Benchmark.md`'s standing rule).

## ★★★ N1 — PROTECTING THE USER FROM A PLAUSIBLE WRONG NUMBER

> **THE ASK (user, 2026-08-25):** *"how do we protect the user from getting wrong results by selecting a
> bad combo … I like to think of a grad student on their first day trying to get some results out of our
> code."*

Today's session produced **four** different collapsed states (−45.5, −46.3, −56.4, −38.5), every one of
them a plausible-looking Hartree number from a run that did not converge.  `SolidCalculation` documents
itself as *"A converged (or ATTEMPTED) periodic SCF"* and then offers `double Energy() const` with **no
precondition**.  So the trap does not even require a bad combination: it requires only not calling
`Converged()`.

**THE ORGANISING PRINCIPLE: you cannot enumerate bad COMBOS, but you can detect bad OUTCOMES.**  T1–T3
below catch combinations nobody has thought of yet; T4–T5 only catch the ones we have.  That asymmetry is
why the ordering is what it is.

- **T1 — MAKE A WRONG ANSWER UNRETURNABLE (types, not checks).  The single highest-value change found
  today.**  `Converge()` hands back a `ConvergedCalculation` (or `std::optional`); `Energy()`,
  `EnergyTerms()`, `Density()` exist ONLY on that type.  Then no one can *write* code that reads an
  unconverged energy — not by discipline, by construction.  This is exactly the project's standing
  *prefer build-failure to runtime crash / give capabilities only to types that have them* bias
  (`feedback_compile_time_over_runtime`), applied to the facade.
- **T2 — A POSTCONDITION ON THE IMPOSITION.**  Imposing an AFM Shubnikov group and converging to
  \f$m_A=m_B=0\f$ is a run that **contradicted its own constraint**; that is a postcondition, not a
  diagnostic, and it should fail hard.  ★ Newly possible: item 3 restored the site blocks, so this can use
  the **integrated** site moment instead of the point probe that misled this campaign for months.  The SSB
  release-audit is the existing symmetry analogue.
- **T3 — ALWAYS-ON OUTCOME DETECTORS.**  (a) **Eee doubling** — validated today (13.48 → 29.0 → 35.1),
  already computed, zero cost, and it ranked four failures correctly without being designed for any of
  them.  (b) **Moment collapse** — the `** DIED at iteration N` logic EXISTS but lives in the MnO *test*,
  not the library, so no user gets it.
- **T4 — STOP EXPOSING KNOBS THAT INTERACT.**  `KerkerG0`, `PulayDepth` and the XC ρ source are not three
  independent physics choices; today's evidence is that getting the combination wrong silently costs the
  magnetic basin.  Per *prefer classes to answer high-level questions*: derive a **MixingPolicy** from
  (ordering, cell, functional).  The user says "AFM MnO", not `G0=1.0, Pulay=0`.  Expert overrides stay,
  behind a named policy, so a bad combo becomes a deliberate act.
- **T5 — SELF-DESCRIPTION (index item 2, widened).**  A run must state mixer, \f$G_0\f$, Pulay depth, XC ρ
  source, thread state and accelerations.  ⚠ Measured today: with `MNO_KERKER_G0=0` the run prints **no
  mixer line at all** — the fallback to linear D-mixing is entirely silent.

## ✅ N2 — THE DENSITY G BALL: RESOLVED (the lobes are BAND-LIMITING)  →  `doc/OpenWork_History3.md`

The ball-vs-raw question was answered and lever B REFUTED on measurement (6e-5); filed under N4.

## ★★ N3 — THE MIXING POLICY: CHARGE AND SPIN ARE DIFFERENT CHANNELS

Kerker is applied PER SPIN CHANNEL (`[Kerker] ENABLED (↑)` / `(↓)`), and by linearity that damps the
CHARGE and the SPIN channel alike.  **But \f$v_H\f$ depends only on \f$\rho_\uparrow+\rho_\downarrow\f$, so
\f$\partial v_H/\partial m\equiv 0\f$: there is no \f$4\pi/G^2\f$ anywhere in the spin channel, for any
magnetic order.**  Kerker's \f$G^2/(G^2+G_0^2)\f$ shape is derived from screening a Coulomb divergence the
spin channel does not have.  cf. VASP's independent `AMIX_MAG`/`BMIX_MAG`.

⚠ **KERNEL ≠ RESPONSE (user's question, "only true for AFMs?").**  The KERNEL statement is universal.  The
RESPONSE is not: \f$\chi_m(q\to0)\f$ genuinely diverges approaching a Stoner instability, so the spin block
of the SCF Jacobian CAN be ill-conditioned at small q for a near-critical ferromagnet — with no Coulomb
divergence in sight.  "Finite kernel" therefore does NOT license "no spin preconditioner"; it licenses "not
*Kerker's*".
✅ **The \f$G=0\f$ carve-out is already correct, for a stated CHARGE reason.**  `FourierMixCD.C:36` exempts
\f$G=0\f$ (`f=1`, "not frozen") because ρ̃'s (0,0,0) coefficient is a shape-dependent fit projection and
*"freezing it strands the XC's mean density at the seed value."*  That same line is what keeps the TOTAL
MOMENT — the FM order parameter — mixable, and nothing records it.  The FM arm is protected by accident of
a charge-motivated decision: fragile documentation, not fragile code.

⚠⚠ **DO NOT LAND N3 BEFORE N1's DETECTORS.**  Today's AFM basin is propped up by the very behaviour N3
removes (ITEM 1 MEASURED, the dose-response table), so a correct fix will *look* like a regression on MnO.

## DEFECTS FOUND EN ROUTE, 2026-08-25 (recorded, not fixed)

- **SEGFAULT (rc=139), pre-existing and flag-INDEPENDENT:** `MNO_MOM=1 MNO_MOM_START=1 MNO_MOM_SEED=1`
  with `MNO_ANNEAL` + `MNO_ANNEAL_PENALTY` crashes at SCF start.  Armed AND control both die, so it is not
  `GPW_XC_DM_SOURCE`.  Not diagnosed.
- **STALE COMMENT** `GPW_SCF_UT.C` *"this KerkerG0 is currently INERT on the AFM arm"* (2026-08-07) — the
  spin-resolved ρ̃ mixer it was waiting on has landed; the baseline prints `[Kerker] ENABLED (↑)/(↓)`.
- **THE LINEAR-MIXER FALLBACK IS SILENT** — see N1/T5.
- ⚠ **A TRAP FOR THE NEXT SESSION:** `MNO_MOM=0` sets `UseMOM=false`, so `MNO_ANNEAL_PENALTY` has no MOM
  reference and silently does nothing.  Cost me one arm before I noticed.

## Where we are, in one paragraph

The **real-TRIM scalar-type track is COMPLETE** (`doc/RealComplexPlan.md`; TRIM blocks build real by
default).  The **runtime campaign has had four rounds** (`doc/GPWPlan1.md`): threading + BLAS (round 1–2),
the run-length collocation streams + the ρ̃ half-space fold (round 3), and the shell-blocked box walk
(round 4).  MnO free-run per-iteration SCF is down ~1.4× and stream RAM 5.78 → 3.70 GB.  Two charter
premises were **narrowed by measurement** in the process — "the collocation streams are complex-bound"
(they are DRAM-bandwidth bound) and "a CP2K-class recompute kernel unblocks the CP2K-span cell" (the
over-budget regime is per-term EMIT bound, so it needs fewer terms, not cheaper ones).  That is the
standing lesson: **cost attribution before optimisation.**  On accuracy, the sharpest coordinate we have
is the **VA (N=118) exact-span code-vs-code table** (`doc/SphericalLatticePlan.md`), which has been
waiting on a term-by-term breakdown to name its operator.

---

## The plan, in order

### ✅ Step 0 — FIX THE INSTRUMENTS — DONE  →  `doc/OpenWork_History3.md`

0a (the INTEGRATED site moment replacing the point probe), 0b, and 0c (`report::RunElapsed()` stamping
every emitted item — index item 5, done 2026-09-06) all landed.  ★ The durable ruling that came out of
0a is a standing rule, not history: **report an INTEGRATED measurable quantity, never a point sample of
a field.**

### ✅ Step 1 — THE HEAD-TO-HEAD TABLE — IT STANDS  →  `doc/Benchmark.md` (live) + `doc/OpenWork_History3.md`

The table is a standing instrument in its own file now; `doc/BenchmarkHistory.md` holds its record.
⚠ Read `doc/Benchmark.md` §5a and COPY the run command — never reconstruct it.

### Step 2 — ARM THE SYMMETRY FOLDS  ·  plan: `doc/SymmetryUpgradePlan.md`

**The biggest single runtime multiplier left, by an order of magnitude** — this is Step 3's work done
properly, not a separate track.  Inventory (2026-08-15, re-verified 2026-08-19):

- the {G}-star fold is wired at exactly TWO static sites (the local-PP sweeps, imposed runs only);
- the **per-iteration G-space consumers are UNFOLDED** — ρ̃, the Poisson multiply, the V_xc gathers, the
  G_ERI3 columns, the seed structure factors: **12–24× on the MnO magnetic group, 48× cubic**, unclaimed;
- ~~the REAL-space T3 pair-stream orbit fold is BUILT but opt-in~~ → **✅ ARMED BY DEFAULT 2026-08-19**
  on every imposed Γ run (below).

**★ T3 IS ARMED — and the "auto-arm criterion" the retraction asked for turned out not to exist, because
the problem was never *which runs are safe*.**  The 2026-08-03 retraction stood on this: the fold imposed
STRICTLY MORE than `imposeSymmetry` itself, so a degenerate open shell (imposed Si p²-in-a-box) flipped
into charge-transfer sloshing ~0.26 Ha off.  The actual defect was one line of the replay: reading the
representative's own \f$D_{ij}\f$ **SAMPLES** the pair orbit, and sampling equals projecting only if D is
already symmetric.  The replay now reads the **orbit-projected** D (`FoldProjectedD`), and since
collocation is linear and equivariant, \f$P\rho_{\rm red}[D]=\rho[PD]=P\rho_{\rm full}[D]\f$ for **any**
iterate — folded and unfolded imposed runs solve the same equations, so arming is a pure cost decision.
Measured on the exact retraction cell: **ΔE(fold on/off) = 1.3e-8 Ha where it was 0.26 Ha**
(`StreamFoldOpenShellMatchesUnfolded_SiAtomInBox`), and the closed-shell Si-diamond A/B is unchanged at
8e-7.  `GPW_STREAM_FOLD=0` is now the opt-OUT.

**★ MEASURED ON THE BENCHMARK ROW (MnO AFM-II Γ, VA span, `MNO_IMPOSE=1`, identical command, same box):
20m05s → 13m25s wall, 2240 → 1809 s CPU, and 4947 → 1349 MB peak RSS — 1.24× the CPU and 3.7× THE RAM**,
at an energy identical to nine significant figures (−61.402976200 → −61.40297623) with `m_stag` still
±0.6667.  Per bucket: pair scatter 263.4 → 41.0 s (6.4×), pair gather 167.6 → 23.1 s (7.3×), stream build
110.1 → 28.2 s (3.9×) — better than the 4.60× rep-pair reduction, because the pairs the fold drops are the
expensive ones.  **Two consequences for this plan.**  (1) **Step 4 (RAM) is largely ANSWERED** — the RAM
half of the MnO gap was mostly the streams, and it fell from 23× CP2K to 6.2× without a campaign.
(2) **The Φ-table build is now HALF the run** (379 s of 805, untouched by folding), so Step 3's Φ-screening
item is the whole of the remaining gap on this row.  Caveat worth keeping: 4.60× is this cell's orbit
factor under the σ=None subgroup (12 of its 24 magnetic ops — a flip op relates D↑ to D↓ and may not fold a
per-channel stream), not the 71× the high-symmetry diamond gate showed.  An orbit factor belongs to a cell.

**And the trap it walked into, worth remembering:** the first cut also projected the integrate-back's
`screenD` — which is a matrix of |D| MAGNITUDES, not a density matrix.  Signed averaging cancelled
mixed-σ orbits to ~0, the D-aware screen dropped live terms, and the imposed O2 triplet collapsed 2.3 Ha.
The suite caught it in one sweep.  Fix: a screen is reduced by the orbit **MAX** (`FoldScreenMax`), never
by the signed average.  **PIN: ask what a matrix MEANS before symmetrizing it.**  New coverage: the
screened integrate-back arm in `StreamFoldGate`, and a DIMER-IN-A-BOX cell — the single-atom Si cells miss
this bug entirely (8.9e-16), the dimer catches it (34% of scale).

Remaining work: **T3.4b** (multi-k per-block arming — union-of-reps stream caches or the star-summed joint
scatter; Γ-only is done) and extending the ball fold to the per-iteration G-space sites.  Caveats: the FFT
itself does not fold trivially.  For scale: rounds 3–4 bought 1.4× by hand-tuning kernels.

**What `MNO_IMPOSE=0` is for — ANSWERED by the user 2026-08-20, and it is a new track, not a knob.**  Two
routes stay legitimate: a FREE run is the DEFAULT and first-class ("some user just wants to run with no
symmetry and see what happens") — correct answer, honest report of the symmetry found, paying only time,
now visible as `NONE` at all three `[fold]` sites and measured at 1.5× wall / 3.7× RAM on MnO.  The
METHODICAL route is the **SSB DESCENT** (`doc/SymmetryUpgradePlan.md` §3b): converge imposed → save the CD
→ a symmetry-ANALYSIS run releases the imposition and ranks candidate SUBGROUPS with weights → re-impose
each and let the energies decide.  It removes the guess from today's release-check, which needs the
symmetry-broken seed (i.e. the answer) handed to it — MnO's AFM-II was assumed, never derived.  §3b carries
the design, the inventory of what exists (`SymmetryDefects`, the ops chokepoint) versus what does not
(`Impose::Subgroup`, subgroup closure, crystal irreps, **CD persistence — nothing at all**), and the one
repair the design needs: step 3 cannot be a single free iteration, because the symmetric solution is a
stationary point of the free map and SSB is second-order — it must measure GROWTH or CURVATURE.

### Step 3 — RUNTIME, CONTINUED  ·  plan: `doc/GPWPlan1.md`

Measure against Step 1, after Step 2 (folding changes what is hot).

- **Φ-table BUILD — ✅ 4.6× TAKEN 2026-08-20, and NOT by the mechanism this item predicted.**  The bucket
  went **379.2 → 83.3 s** on the MnO benchmark row (190 → 36.8 per anneal stage), taking the whole row to
  **8m36s wall / 1554 s CPU** — with Step 2 that is **2.34× wall and 1.44× CPU off the banked cut**, at an
  energy identical to nine significant figures.  Cost attribution first, as always, and it found three
  things — of which the item's own hypothesis ("Φ is stored dense, batch the mesh") was the SMALLEST:
  1. **The dense cart→spherical transform was inside the image loop (the big one, 113 → 36.8 s).**
     `SphericalLatticeView::operator()` applies `T^T v` per call, and the periodic caller
     (`GPW_Evaluator::Eval`) calls it ONCE PER LATTICE IMAGE per mesh point — so a 122×118 dense mat-vec
     ran ~150× more often than the physics needs, at **10.6% of the entire run's cycles**.  T is
     block-diagonal per shell and IDENTITY for every s/p shell, so it is a 1-to-6-term sum per output:
     evaluating only its nonzeros is bit-identical (the skipped terms contribute a hard 0.0) and ~50×
     cheaper.  **The deeper fix is its own item below — promoted out of this sub-bullet 2026-08-20 because
     a still-open fix filed under a ✅ heading is a fix that gets lost.**
  2. **The pointwise sweep had no magnitude screen (190 → 120 s).**  `PG_Cart::IrrepBasisSet::operator()`
     evaluated every contracted radial at every point — including an α=36 Mn d shell at 20+ bohr, where its
     value is ~e⁻¹⁴⁰⁰⁰ — because the image list necessarily reaches as far as the most DIFFUSE function.
     Now screened on a cached per-function radius (`PGData::Reaches`, the same ε-magnitude discipline as
     `BetaSupportRadius`; ε=1e-10 sits far below the ~1e-4 quadrature error).
  3. **Column-major fill: REFUTED, and the question is CLOSED.**  `mat_t` is column-major while the sweep
     produces rows, so the fill strides by npts per element — which looks like the classic cache disaster
     and is not.  Element (g,i) sits at `i*npts+g`, so consecutive POINTS write ADJACENT elements of the
     SAME line: 122 interleaved SEQUENTIAL streams, ~1 miss per 8 points per function.  **Timed in
     isolation at the real shape (99370×122): 13 ms column-major vs 8 ms row-major — 1.1 ns/element, the
     memory-BANDWIDTH floor for writing a 97 MB table at all — and the row-major route then owes a 21 ms
     transpose, so it is a net LOSS.**  The entire fill is 0.04% of the bucket.  An in-run A/B had
     suggested 7 s; that was single-sample noise beside a threaded Becke mesh, and the isolated
     measurement is what settled it.  **Lesson, again: measure the SUSPECT alone before believing a
     difference of two whole-run timings.**
- **⚠ THE LEDGER MEASURES WALL, `perf` MEASURES CYCLES — do not read one against the other.**  The single
  biggest CPU consumer in the whole code is `BeckeCutoff` at **~50% of cycles**, but that loop is the one
  that is THREADED BY DEFAULT, so on 16 cores it is only ~35 s of wall and reads as a small ledger bucket.
  The reverse held for Φ: 15.5% of cycles, 47% of wall, because it is serial.  Both numbers were right and
  reading the flat profile against the ledger cost a wrong conclusion before the arithmetic closed.
- **Becke partition — ✅ 4.44× TAKEN 2026-08-20, and BOTH named suspects were WRONG.**  The item said to
  choose between the shell-convergence retest and the final P-set double loop.  A per-point call census
  (`GPW_BECKE_COUNT=1`, now permanent and env-gated) answered it by COUNTING rather than by profiling a
  threaded loop — and then refuted the fix that the answer implied:

  | | MnO AFM-II VA |
  |---|---|
  | competitor images / live point | **3183** (max 27436, shells ≤10) |
  | P-set members / point | 100.5 |
  | `BeckeCutoff`, one mesh build | **2.53e10** |
  | share: final double loop / shell retest | **97.3% / 2.7%** |

  So the retest is not the problem despite growing as s³.  But nothing that REMOVES work from the double
  loop is big either, all measured on MnO: a κ-distance factor screen **1.17×**, PERFECT screening
  **1.71×** (an unreachable ceiling — 58% of the factors genuinely differ from 1), and skipping the **69%**
  of P-set members whose product underflows to exactly 0 is worth **1.02×**, because the existing `P>0`
  early exit already disposes of them in 2.3% of the work.  *The 31.7:1 ratio between the image list and
  the P-set is not waste.*
  **The lever was the multiplier neither suspect named: ε.**  ε fixes |im|, and |im| multiplies the whole
  O(|P-set|·|im|) loop — so ε scales the dominant cost instead of shaving it.  It had only ever been probed
  TIGHTER.  Loosening 1e-8 → **1e-6** takes the MnO becke build **36.95 → 8.33 s threaded (294 → ~66 s
  serial)** at Etot identical to 9 s.f. (−52.48387019 → −52.48387023) and site moments identical to 6
  digits.  **PIN: the binding gate is `BeckeEquivalentSitesOwnEqualShares`** (site shares equal to 1e-8
  RELATIVE — a PARTITION property, far sharper than any integration test); it survives 1e-6 and breaks at
  1e-5, so the default keeps a decade of margin.  ⚠ Unlike the Φ-build 4.6× this is a **TOLERANCE trade,
  not a bit-identical restructuring** — the weights move at ~1e-6 relative.  `GPW_BECKE_EPS` overrides.
  **★ ON THE BENCHMARK ROW: 8m36s → 6m56s wall, 1554 → 663 s CPU (2.34×), 1350 → 1323 MB**, at
  −61.40297622 → −61.40297621.  Cumulatively the row is now **2.90× wall / 3.38× CPU** off the banked cut,
  and the CP2K CPU gap is 6.0× → **1.8×**.
  **⚠ AND IT CAUGHT A FLAW IN THE BENCHMARK'S CPU COLUMN, which is the reverse of the wall/cycles trap
  above.**  The becke build is 294 s SERIAL but bills ~590 s of CPU threaded 16-way (36.95 s × 16), because
  the OpenMP threads BUSY-WAIT at the barrier and CPU time counts spinning as work.  Two anneal stages of
  that was ~1180 s of the row's 1554 s CPU (76% — the earlier "~73% of CPU" estimate was RIGHT), which is
  how 4.44× off a bucket that looked like 38% cut the row by 2.34×.  Reading the serial 294 s against a
  threaded CPU total is what produced the wrong "38%".  **The pin "compare CPU, not wall" assumes CPU
  tracks WORK; where qchem threads and CP2K does not, it does not.**  Serial is the honest comparison.
- **★★★ THE BECKE PARTITION IS COMPUTED PER POINT ON A MESH THAT HAS ONLY 4290 ORBITS — and the fix is
  V1.22, already filed as a CLEANUP item (2026-08-21).**  The imposed run prints
  `invariant mesh 98816 points in 4290 orbits` (avg orbit 23.0).  The Becke weight is SYMMETRY-INVARIANT:
  under an op the distance multiset \f$\{|r-R_b|\}\f$ is preserved with the atoms permuted, so
  \f$w_{g(a)}(gr)=w_a(r)\f$ EXACTLY — and the site-adapted construction already guarantees a partner point
  lies on the partner ATOM's grid.  So the partition is evaluated **98816 times where 4290 would do: up to
  ~23×** on the largest CPU bucket in the code (75 s serial ⇒ ~3 s), on top of the 4.44× ε already taken.
  **This is `doc/CleanupCandidates.md` V1.22** — *"make the drop decision ONCE per representative (angular
  dir × radial shell) and apply it to the whole atom orbit inside the builder — removes the second fold
  pass + the filter"* — filed for CORRECTNESS (the per-point `<eps`/`w>0` decisions are bit-sensitive and
  break orbit consistency, which is why the caller post-filters orbit-incomplete points).  **Its
  PERFORMANCE content was never noticed because it sat under cleanup.**  The same per-representative loop
  that fixes the drop asymmetry also removes the redundant partitions.
  ⚠ Scope: IMPOSED runs only — a free run has no orbits and keeps the full cost.  The benchmark row is
  imposed; the free production run is not.  And V1.22 is on the deliberately-anchor-moving list, so it
  belongs in the batched re-pin with §K and the Δρ/N gate.
- **Becke partition, what is LEFT in that loop** — ⚠ **RE-RANKED 2026-08-21, and the two interact:**
  the vectorization item is worth doing, the `norm()` table probably is not, and BOTH are small beside
  V1.22 above.
  0. **VERDICT (asked directly, 2026-08-21).**  `-march=native` alone measured **1.13×** on this bucket with
     the loop still SCALAR, so item 2 (vectorization) has real headroom — the μ computation, the 4-iteration
     polynomial and the product-reduction all vectorize once the data-dependent exit goes, and that exit
     saves only 2.3% of work.  Item 1 (the `norm()` table) is DOUBTFUL: natom²·(4s+1)³ ≈ 250k entries ≈
     2 MB reached through 3-D index arithmetic, against a HARDWARE sqrt at ~15 pipelined cycles — the same
     shape as the column-major fill that measured out a net LOSS.  **And the two point OPPOSITE ways:**
     vectorizing turns the sqrt into a 4-wide `sqrtpd` but turns the table into a GATHER.  So: do the
     vectorization first; the `norm()` item may then be moot or actively harmful, and should be re-measured
     rather than assumed.
  1. **The `norm()` is redundant across points.**  `R_ij = |R_i − R_j|` is recomputed inside the per-point
     loop — 2.46e10 square roots — but it is **point-INDEPENDENT**: image index t is always the same
     (atom, cell-offset) tuple, so `R_i − R_j` depends only on the offset DIFFERENCE, not on the point's
     own cell.  The distinct values number ~natom²·(2s+1)³ ≈ 148k (~1 MB), tabulable once per mesh build
     and shared across threads.  Worth up to ~2× of the pair cost, and it is exact.
  2. **The inner loop cannot vectorize.**  `BeckeImage` is array-of-structs, and the per-element `P>0`
     break is a data-dependent exit — so the hottest loop in the code is scalar and latency-bound (each
     `BeckeCutoff` is 4 serially-dependent iterations).  The census says those early exits save only 2.3%,
     which is a poor price for blocking SIMD; chunking the exit + an SoA layout would unlock it.
  3. **The build flags undercut both.**  `qcStructure` compiles with `-O3 … -O2 -g` — clang takes the LAST
     `-O`, so it builds at **-O2** — and there is no `-march=native`, i.e. baseline SSE2 at 2 doubles per
     vector.  (Whole-tree change, so it interacts with the BLAS pin; not a Becke-local decision.)
- **★ THE Φ BLOCH POINT-SUM SEAM — ✅ DONE 2026-08-20 (`55de8578`), and it is the seam Φ-sparsity extends.**  `GPW_Evaluator::Eval`
  (`Evaluator.C:801`) Bloch-sums the orbitals by calling `(*itsOrb)(r-R_k)` **once per lattice image per
  mesh point** — and when `GPW_SPHERICAL=1` that call IS the cart→spherical transform.  The transform is
  LINEAR, so \f$\sum_k \phi_k T^{\!\top} v_k = T^{\!\top}\sum_k \phi_k v_k\f$: summing the images in
  CARTESIAN and transforming ONCE per point is the honest structure.  `Eval` cannot do that, because it
  sees only an abstract `Real_BS` and has no idea a transform is hiding inside — so the Bloch POINT-sum has
  to move into the periodic seam (`Molecule::LatticeSum1E`), beside the per-shell-pair enumeration that
  already lives there for the 1E matrices.  **That is an interface question — it needs a new face method
  (a Bloch point-value), so it does not get taken unilaterally.**
  **Three reasons it goes BEFORE the sparsity item, none of them its size:**
  1. **Sparsity rewrites this exact caller.**  Φ is filled by `PhiAt` → `Eval(R[k])` per point
     (`Evaluator.C:1215`); sparsity batches that fill per atom / radial shell with a significant-function
     list.  Seam first ⇒ sparsity is built on it.  Sparsity first ⇒ the deeper fix becomes a SECOND rewrite
     of code just written, threading a significance list through an image loop that should not hold the
     transform at all.
  2. **It retires the code's own last exception.**  `Evaluator.C:243` states: *"THERE IS NO CUT … The ONE
     remaining explicit image list is the INTERNAL Bloch-orbital set for Eval/EvalGradient."*  It is
     eps-DERIVED and justified as a superset of what the screen keeps — not a violation — but it is the
     last place GPW enumerates images itself instead of letting the molecular seam do it.  This change
     closes the exception and the redundancy together.
  3. **The structures compose.**  T is block-diagonal per shell, which is the granularity sparsity wants to
     batch on; a once-per-point transform composes with that, one buried in the image loop does not.
  **MEASURED FIRST, and the measurement DEMOTED it before it was built** — `perf` put the transform layer
  (`SphericalView_IBS::operator()`) at **2.9% of run cycles**, so this was booked as a SEAM CLEANUP that
  sparsity builds on, NOT as a speedup.  Built on that basis and it delivered the predicted size:
  **Φ bucket 26.36 → 20.39 s serial (1.29×)**, Etot unchanged to every printed digit.
  **THE SHAPE (user, 2026-08-20):** `VectorFunction<T>` gains the point-set `operator()` **MIRRORING the one
  `ScalarFunction<T>` already had**, same defaulted pointwise-loop body — *"any derived class with a faster
  overload is welcome to re-implement."*  That also retired my own objection: the "ISP sin" that
  `VectorFunction.C`'s header records was a **Mesh** dependency, not a bulk overload, and the sibling face
  had carried exactly this signature all along.  Return is a MATRIX where `ScalarFunction`'s is a vector,
  because the pointwise form returns a vector where `ScalarFunction`'s returns a scalar.
  Chunked in `MakePhi` so its threading survives (qcMath is a leaf, cannot host OpenMP, and molecular bases
  reach `MakePhi` too).  The enumeration reproduces the caller's former rule exactly, so the untransformed
  result is bit-identical and only the transform reordering moves (~1e-16).
  ⚠ **`Eval`/`EvalGradient` still run their own image loop** — the per-point callers (KB quadrature) were
  NOT migrated, because the seam re-derives its offset list per call and pointwise delegation would pay
  that per point.  So the *"ONE remaining explicit image list"* is narrowed, not retired; retiring it needs
  the offsets cached on the seam side, which is its own increment.
- **★ AND THE MEASUREMENT FOUND SOMETHING BIGGER NEXT DOOR — ✅ `87da4854`, 1.35× for one line.**
  `IrrepBasisSet::size()` is VIRTUAL and cross-`.so`, so as the CONDITION of
  `for (size_t i=0;i<size();i++)` in the pointwise basis sweep it was re-dispatched through the PLT every
  iteration — billed as three symbols totalling **2.73% of run cycles**, i.e. MORE than the entire seam
  fix it was measured in service of.  Hoisted (bit-identical): **Φ bucket 35.69 → 26.36 s**.  The same
  `i<size()` pattern at four other sites is deliberately NOT touched — they are one-time setup or memoized
  and absent from the profile.  **Cumulative on the Φ bucket today: 35.69 → 20.39 s, 1.75×.**
- **Φ-table SPARSITY — ⛔ PREMISE REFUTED ON THIS CELL, 2026-08-20.  Do not build it for MnO.**  The item
  said Φ is stored dense (npts × n) while the true object is SPARSE, so batching the mesh (per atom / per
  radial shell) with a per-batch significant-function list makes every Φ-shaped cost — the ρ GEMM and the
  H_xc GEMM — O(npts·n_sig²).  The win is QUADRATIC in n_sig, so n_sig was worth measuring before building
  the machinery.  `GPW_PHI_SPARSITY=1` (new, permanent) reports it:

  | batching | imposed (97160 pts) | free (31664 pts) |
  |---|---|---|
  | **per-atom SITE blocks** (the item's own proposal) | *mesh has none — see below* | **1.05×** (n_sig mean 115/118) |
  | 4096-pt batches | 1.71× | 1.59× |
  | 1024-pt batches | 1.89× | 2.30× |
  | 256-pt batches | 1.97× | 2.89× |
  | 64-pt batches | 2.08× (n_sig max 110/118) | 3.25× (n_sig max 108/118) |
  | **Φ nonzero fraction** | **48.1%** | **46.5%** |

  **Φ is HALF DENSE, not sparse**, and the per-ATOM batching the item proposed is worth **1.05×** — with 4
  atoms in a small cell essentially every function reaches every atom's grid.  Even 64-point batches only
  take n_sig from 118 to 110.  And these are CEILINGS: they assume zero batching overhead against a GEMM
  whose efficiency small batches would wreck.  **The item's own caveat — "the win grows with cell size;
  MnO's 4 atoms understate it" — was exactly right, and it is the whole story.**  So Φ-sparsity is NOT the
  lever for the benchmark row; it stays a real idea for LARGE cells (the battery north-star's supercells),
  where it should be re-measured with this instrument before anyone builds it.  **The ρ GEMM is still the
  top wall bucket, and it now needs a different idea.**
- **★★★ Vxc MUST BE FED THE DM ρ(r) — AND THE PROJECT ALREADY DECIDED THAT.  The Becke XC path does not
  honour it, which is both an ACCURACY regression and the largest per-iteration cost (user, 2026-08-20).**
  `doc/GPWPlan.md:286` records the original insight verbatim — *"feed XC the DM-ρ, pointwise NON-NEGATIVE by
  construction (PSD D ⇒ φᵀDφ ≥ 0)"* — and §0.5(f2), *"the DM-ρ raw XC feed"*, was **BUILT + ACCEPTED
  2026-07-23** on measured evidence: the collapse basin removed, NaF 1.3 → **0.2 mHa** vs CP2K, and *"the raw
  feed removed the ball-Gibbs noise from the XC residual"* so the coarse stage converged in **45 iterations
  instead of 515**.
  **What happens now instead.**  `XC_SinglesQuadrature::RhoPol` takes the DM route only when the density it is
  handed carries a D (`cPolarized_CD` → `DM_RhoAtPoints`).  On the ρ̃-mixing recipes (Kerker/Pulay — the MnO
  production recipe) the mixed density is a `PolarizedMixCD over FourierMixCD` with **no D**, so XC falls to
  the matrix-free branch: a batched inverse FT over the whole {G} at every mesh point, every iteration.
  That is:
  - **an ACCURACY regression** — a BAND-LIMITED ρ, with Gibbs ringing, locally-negative lobes and the ρ>0
    guard, fed to the ATOM-CENTRED Becke mesh whose entire purpose is to resolve the sharp core features a
    band-limited grid cannot.  High-accuracy quadrature on a low-accuracy integrand.  The DM-ρ is exact
    pointwise AND non-negative by construction, so the guard becomes unnecessary rather than tuned;
  - **the largest per-iteration cost**: 35.0 s / 6 iterations serial against the DM GEMM's 1.70 s (84.1 s on
    the benchmark row) — so honouring (f2) here is worth roughly **20×** on that bucket, before the low-rank
    Cholesky route makes the GEMM cheaper still.

  **The design question to settle (do NOT assume the answer).**  The mixer exists to damp the SCF, and the
  Hamiltonian is normally built from the MIXED density; Kerker in particular is inherently a G-space
  preconditioner and has no DM-space form.  So feeding XC from D while Hartree keeps ρ̃_mix means the two
  terms see different densities off-convergence.  **Fixed-point argument in favour:** at self-consistency
  \f$\rho[D]=\rho_{mix}\f$, so the converged answer is UNCHANGED — only the SCF trajectory differs.  That
  makes it a convergence-behaviour question, not a correctness one, and it wants measuring on the recipes
  that currently lean on Kerker (NaF's low-G slosh, MnO's AFM basin) rather than deciding on principle.
  Options: (a) mix in DM space (`LinearMixer` is convex and PSD-preserving, but loses Kerker's low-G
  preconditioning); (b) Hartree from ρ̃_mix, XC from the DM — what (f2) appears to have intended;
  (c) reconstruct a DM-backed density from the mixed ρ̃ — not generally possible.
  **THE SHAPE — corrected twice by the user, and the algebra decides most of it.**
  \f$\tilde\rho[D]\f$ is the COLLOCATED (band-limited) density, so \f$\rho_{mix}-\rho[D]_{exact}\f$ and
  \f$\mathrm{IFT}[\tilde\rho_{mix}-\tilde\rho[D]]\f$ are DIFFERENT objects — they differ by the cusp content,
  the largest thing in the problem.  Writing both as "Δ" hid that.  The computable form is:
  \f[ \rho_{XC}(r)=\underbrace{\rho_{mix}(r)}_{\text{band-limited}}+\underbrace{\big(\rho[D](r)_{exact}-\rho[D](r)_{BL}\big)}_{\text{the CUSP DEFICIT of band-limiting}} \f]
  i.e. the mixed density with the sharp core content — the part only the DM can supply, and the part the
  atom-centred mesh exists for — restored.  The cusp deficit is nearly ITERATION-INVARIANT (core electrons
  barely move), so borrowing it from the current D is accurate well before convergence, and exact at it.
  ⚠ **Retraction:** this does NOT give Hartree and XC the identical density, as first claimed.  The honest
  defence is that the Poisson operator is LINEAR and diagonal in G (band-limiting converges fast there)
  while XC is a NONLINEAR POINTWISE functional needing the real-space cusp — each term gets the
  representation its operator requires.  That is weaker than "same array", so it no longer by itself
  outranks simply feeding XC from D.

  **★ IS THERE A PERFORMANCE WIN?  NOT FROM THE DM ROUTE ITSELF — and the algebra says where it can come
  from.**  `FourierMixCD.C:65` samples \f$\rho(r)=\sum_G\tilde\rho(G)e^{iG\cdot r}\f$ by DIRECT SUMMATION at
  every point, so the 35 s is O(npts·|{G}|) — 97k points × |{G}| per iteration.  \f$\rho[D]_{exact}\f$ ADDS
  the GEMM; it removes no sampling.  Subtracting in G-space first leaves ONE sampling, so the scheme costs
  today's sampling PLUS a GEMM — **a net LOSS unless the correction's {G} can be truncated.**
  With Kerker (\f$\tilde\rho_{mix}=\tilde\rho_{in}+\alpha f\tilde\delta\f$, \f$f=G^2/(G^2+G_0^2)\f$,
  \f$\tilde\delta=\tilde\rho_{out}-\tilde\rho_{in}\f$, and \f$\tilde\rho[D]=\tilde\rho_{out}\f$):
  \f[ \tilde\rho_{mix}-\tilde\rho[D]=\big(\alpha f(G)-1\big)\,\tilde\delta \f]
  ⚠ so at HIGH G, \f$f\to1\f$ and this tends to \f$-(1-\alpha)\tilde\delta\f$ — **75% of the residual's
  high-G content at α=0.25.  It is NOT low-G**, and the earlier "Kerker weights it toward low G" claim was
  WRONG (Kerker suppresses the INCREMENT at low G, which makes the correction ≈ \f$-\tilde\delta\f$ there
  and leaves it nearly untouched at high G).
  **But it is proportional to the SCF RESIDUAL**, which → 0.  So the truncation error is bounded by
  \f$(1-\alpha)|\tilde\delta_{highG}|\f$ and vanishes as the run converges: an **ADAPTIVE G-ball keyed to
  \f$|\tilde\delta|\f$** is cheap late, accurate early, and the CONVERGED answer is exactly
  \f$\rho[D]_{exact}\f$ either way.  Second route, better if the correction is smooth: FFT it to a coarse
  UNIFORM grid and INTERPOLATE to the mesh points — O(N log N + npts) instead of O(npts·|{G}|) — legitimate
  only because the correction has no cusps, which is precisely why the full ρ cannot be treated that way.

  **SO THE TWO HALVES SEPARATE, AND SHOULD BE JUDGED SEPARATELY:**
  - **ACCURACY — certain, and the priority (user).**  Exact cusps, ρ≥0 by construction, the ρ>0 guard made
    near-vacuous, and the (f2) property recovered.  Costs a GEMM.
  - **PERFORMANCE — contingent.**  ~20× is available from feeding XC the DM ρ ALONE (drop the correction
    entirely: cost becomes just the GEMM), but then XC sees \f$\rho_{out}\f$, not \f$\rho_{mix}\f$ — the
    SCF-trajectory question, which at the fixed point changes nothing.  With the correction kept, the win
    exists only via the adaptive ball or the interpolate route.
  **GATING INSTRUMENT (build first):** dump the radial spectrum of \f$(\alpha f-1)\tilde\delta\f$ vs |G| per
  iteration.  It is computable ENTIRELY INSIDE THE MIXER — \f$\tilde\delta\f$ is already formed there for
  `ApplySpectralFilter` — so no new plumbing between the mixer and XC is needed to answer it.

  **★★★ AND THE FACTORED FORM CHANGES THE OBJECT FROM PAIRS TO SINGLES (user, 2026-08-21) — which is a
  bigger structural win than the XC feed it came from.**
  \f$\rho=\sum_{ij}D_{ij}\chi_i\chi_j\f$ is a sum over PAIRS; \f$\rho=\sum_m|\psi_m|^2\f$ with
  \f$\psi_m=\sum_i L_{im}\chi_i\f$ is a sum over SINGLES.  On the MnO benchmark cell that is:

  | object | count |
  |---|---|
  | pairs (raw) | **8778** |
  | pairs after the T3 fold | 1909 |
  | **single basis functions** | **118** |
  | **orbitals after the low-rank contraction** | **14–17** |

  74× fewer objects than raw pairs, 16× fewer than folded representatives.  **And this is not confined to
  the XC feed:** the Φ table already IS a singles stream — which is exactly why `DM_RhoAtPoints` costs 1.7 s
  — while the COLLOCATION path still pays pairs, and it is the pair streams that drove Step 4's RAM (5.78 →
  3.70 → 1.35 GB) and the scatter/gather buckets (263 → 41 s under the fold).  A factored density would
  replace the pair scatter with a singles sweep + GEMM + row-norm, and shrink the stream caches by the same
  ratio.
  ⚠ **TWO CORRECTIONS TO THE ABOVE, both making the comparison harder than the table suggests.**
  1. **"Singles are more diffuse" is nearly a NON-argument (user, 2026-08-21): "only by a factor of 2 — a
     diffuse orbital paired with itself is still diffuse."**  The most diffuse PAIR is
     diffuse×itself at \f$2\alpha_{\min}\f$ against the single's \f$\alpha_{\min}\f$: a factor 2 in
     exponent, \f$\sqrt2\f$ in radius, ~2.8 in box VOLUME.  And collocation cost is dominated by the
     DIFFUSE objects — exactly the pairs that are barely more compact than singles.  So per-object cost is
     comparable and the count ratio translates fairly directly into a cost ratio.  My "singles cover more
     points" counterweight was much weaker than stated.
  2. **★ BUT THE 118-vs-8778 COMPARISON IS APPLES-TO-ORANGES, and that cuts the other way.**  8778 counts
     **(i, j, R)** terms INCLUDING lattice offsets; 118 counts functions with NO offsets.  The singles
     enumeration is **(i, R)** — 118 × the images each function reaches (the KB path quotes ~133 images
     elsewhere), so the comparable figure could be ~15k singles against 8778 pairs, REVERSING the sign.
     The mechanism behind that is real and is the strong counterweight: **pairs SCREEN far harder than
     singles.**  The Gaussian product theorem gives a pair the factor \f$e^{-\mu|R_{ij}|^2}\f$ with
     \f$\mu=\alpha_i\alpha_j/(\alpha_i+\alpha_j)\f$, so distant pairs die exponentially in the SEPARATION,
     while a single is screened only by its own reach.  That is why CP2K collocates pairs.
  **THE INSTRUMENT THIS DEMANDS, before any of it is built:** count \f$(i,R)\f$ SINGLES against
  \f$(i,j,R)\f$ PAIRS *at the same ε*, and weight each by its box volume — the object counts alone decide
  nothing.  The pair side already reports itself (`[fold] collocation streams … 8778 → 1909`); the singles
  side needs the same census over the orbital reach.  **Only that comparison settles the crossover**, and
  it is system-dependent: n=118 here, while the battery-north-star supercells grow n and make pair
  screening bite harder, so the answer may differ between the two targets and BOTH routes may be worth
  keeping.
  ⚠ Note the interaction with the T3 stream fold: the fold's 4.60× is a reduction on PAIRS.  Whether an
  equivalent fold exists on singles (orbitals are not symmetry-adapted individually) is an open question,
  and if not, the honest comparison is 118 singles against 1909 folded pairs, not 8778 raw.

  **PIN: this is the first thing to look at in Step 3, ahead of anything else** — it is the only item that is
  simultaneously an accuracy repair, a documented-decision regression, and a ~20× on the top bucket.

  **⇒ STATE 2026-08-25: THE CURE IS BUILT AND SITTING BEHIND AN ENV FLAG.**  The 2026-08-21/22 session wired
  it (record now in `doc/OpenWork_History2.md`): `ExactSourceOf` cross-casts the mixed density to the
  `cDM_Sourced_CD` the mixer already retains, so the exact ρ[D] is one cast away, and `DampXCChannel`
  applies the SAME α the field mix used (the half-damped map was the first wiring's defect — NaF DIIS 34 →
  100+ iterations until α_eff was read off the mix instead of set).  It is **OPT-IN**: `GPW_XC_DM_SOURCE=1`,
  on the stated ground that *"a trajectory change must earn its place against banked recipes."*  Production
  has therefore never taken it.  Nothing needs designing; what is missing is the measurement that promotes
  it.

  **TWO LIVE CONSEQUENCES, lifted out of the archived record because they still bind:**
  1. ⛔ *"at the fixed point ρ[D]=ρ_mix, so the converged answer is unchanged"* is **WRONG**.  ρ̃ is a
     band-limited FIT PROJECTION, so ρ̃_mix(r) ≠ ρ[D](r) pointwise even at convergence — they have DIFFERENT
     fixed points.  MEASURED: NaF **139 μHa**; MnO individual terms **~100–150 mHa** (Ekin +110, Een −148)
     while the TOTAL moves 8 μHa — the variational signature.  So route (b) is **MORE ACCURATE, not
     neutral**, which changes what the measurement below is even asking.
  2. ⚠ **Step 5 must PIN this flag before its term-by-term CP2K breakdown means anything** — the individual
     terms move by ~100 mHa with it.  (It does NOT explain the ~100 mHa offset Step 5 is chasing: the total
     barely moves.)

  **⇒ ✅ MEASURED 2026-08-25 — THE ANSWER IS NO.  See *"✅ ITEM 1 MEASURED"* above the plan for the full
  evidence.**  The flag takes MnO AFM-II from converged −61.40297621 to collapsed −45.5, in BOTH mixer
  configurations (with and without Pulay history) and with the low-rank factoring both ON and OFF.  The
  cause is NOT ρ[D]: it is that a flat \f$\alpha_{eff}\f$ destroys **Kerker's mode selectivity**, un-damping
  the low-G CHARGE mode 2.4× while barely touching the AFM mode — confirmed by a monotone dose-response in
  Eee (13.5 → 29.0 → 35.1).  ⚠ The three cost/accuracy numbers quoted in THIS bullet are superseded there:
  the ρ̃ bucket is **15% of the row, not 51%**, so the ceiling is **~1.17×, not ~20×**; the ρ≥0 accuracy
  claim is UPHELD (0 negatives in 154 samplings vs 15.2%).  The GOAL survives; the ROUTE does not — its
  successor is the **alias-free G ball** (index row N2).

  **THE ORIGINAL NEXT ACTION, kept for the reasoning (now executed)** (the standing lesson: cost attribution before
  optimisation).  Arm `GPW_XC_DM_SOURCE=1` on the MnO benchmark row and read three numbers off one run:
  per-iteration SCF time (the claim is the matrix-free inverse-FT sweep disappears), the ρ<0 point fraction
  and negative mass (`GPW_RHO_NEGATIVE=1` — the accuracy half; ρ[D] is pointwise non-negative for a PSD D by
  construction), and the iteration count.  Then the flag either earns the default or does not, on evidence
  rather than on the argument above.  ⚠ Compare **CPU time, not wall** (qchem threads here, CP2K does not),
  and note the class names in this item: the engine is `XC_SinglesQuadrature` / `XC_PairQuadrature` since
  the 2026-08-22 absorption.
- **★★ THE ρ GEMM — LOW-RANK D.  ✅ BUILT AND DEFAULT-ON (`QCHEM_DM_LOWRANK=0` is the A/B valve), and then
  ⛔ LARGELY NEUTERED BY ITS OWN MEASUREMENT — read the correction below before quoting the 7–8×.**  It was
  written as *"the successor to the refuted Φ-sparsity item, and the only remaining plan for the largest
  per-iteration bucket"*; the bucket turned out to be something else.  ⚠ It is also the headline example of
  an acceleration CP2K does not run, so it must be OFF for any head-to-head row (index item 2).  `IrrepCD_Core::DM_RhoAtPoints` forms
  \f$\rho(r_g)=\Phi_g^\dagger D\,\Phi_g\f$ as `(P·D)` row-dotted back into `P` — **O(npts·n²)**, i.e. on MnO
  97160 × 118² ≈ 1.35e9 MAC per spin per iteration, and it is a real BLAS `gemm` so there is no waste to
  screen out.  But **D is a density matrix: its rank is the OCCUPIED count, not n.**  With \f$D=LL^\dagger\f$,
  \f$\rho_g=\lVert L^\dagger\Phi_g\rVert^2\f$, so the GEMM becomes npts × n × **r** and the win is **n/r**.
  MnO: n=118, ~13 occupied per spin ⇒ potentially **6–9×**.  Four steps (user's decomposition):
  1. + 2. **ONE pivoted-Cholesky call** — it is rank-REVEALING, so the rank and the factor come together;
     there is no separate rank pass.  `CholeskyPivoted` + `detect_null_gap` already exist in `LASolver`
     (the ortho step drives them).  O(n³) ≈ 5e5 flops against a ~1e8 GEMM — free.
  3. \f$\Psi = P\,L\f$ — (npts×n)·(n×r), the same GEMM, thin instead of square.
  4. \f$\rho_g=\sum_m|\Psi_{gm}|^2\f$ — a row norm, O(npts·r).

  > **★ PIN (user, 2026-08-20), AND ITS LIMIT — the objection to trimmed eigen is INVERSION-driven, and we
  > do not invert here.**  First reading was "always pivoted Cholesky, never a trimmed eigendecomposition",
  > because near-zero eigenvalues arrive CLUSTERED, eigenvectors inside a clustered subspace are
  > ill-conditioned (free to rotate arbitrarily within it), and that noise rides into the kept factor — with
  > the observed consequence being damage to LATE GDM CONVERGENCE.  **The user then narrowed it: that
  > experience comes from ORBITAL work where the factor must be INVERTED** (orthogonalisation, \f$S^{-1/2}\f$),
  > and inversion amplifies near-null noise by \f$1/\lambda\f$.  **Here the factor is only ever MULTIPLIED**
  > (\f$\rho_g=\lVert L^\dagger\Phi_g\rVert^2\f$): the error is \f$\lVert\Delta L\rVert\,\lVert\Phi_g\rVert\f$,
  > bounded, with no small denominator anywhere.  So the historical failure does NOT transfer, and truncated
  > eigen is probably fine.  **Standing preference: pivoted Cholesky where D is PSD** (greedy on the
  > diagonal, truncation bounded by the trailing diagonal, no rotational noise) — but it is a preference,
  > not a prohibition, and if D is genuinely non-PSD then Cholesky is not "less preferred", it is
  > INAPPLICABLE and the eigen split is the route.  **MAY BE MATERIAL-DEPENDENT** (user): a system that
  > appears to need the eigen route is a finding to record, not a silent switch.

  **IS THE D REACHING `DM_RhoAtPoints` EVER NON-PSD?  The maths is one-directional, and the tree's evidence
  points to PSD.**  \f$\rho(r)=\Phi(r)^\dagger D\,\Phi(r)\f$ is a quadratic form in \f$v=\Phi(r)\f$, so
  \f$\rho(r)<0\f$ ANYWHERE proves D is not PSD — but NOT conversely: \f$\{\Phi(r)\}\f$ traces a
  3-dimensional manifold in \f$\mathbb{C}^n\f$, so a negative eigenvalue whose eigenvector is never realised
  as some \f$\Phi(r)\f$ stays invisible, and \f$\rho\ge0\f$ everywhere is compatible with a non-PSD D.
  **The user recalls ρ going negative in PP-DFT work — but WHICH ρ decides whether it bears on D:**
  - the DIRECT form \f$\Phi^\dagger D\Phi\f$ (what `DM_RhoAtPoints` computes) going negative ⇒ D is NOT PSD;
  - a FITTED or FOURIER-TRUNCATED ρ going negative ⇒ says NOTHING about D (Gibbs ringing / aliasing).

  **★ THE USER'S PP-DFT SURPRISE WAS THE FIRST KIND** — the DIRECT form going negative, which is why it was
  a surprise: the fitted/Fourier-truncated case was already well understood and accepted.
  **★★ AND THE CULPRIT IS A PP SUBTLETY THAT BREAKS THE IMPLICATION'S HIDDEN PREMISE (user).**  "Direct form
  negative ⇒ D not PSD" silently assumes \f$\rho(r)\f$ IS \f$\Phi^\dagger D\Phi\f$.  For ULTRASOFT PPs and
  PAW it is not: \f$\rho=\sum_m f_m|\tilde\psi_m|^2+\sum_{ij}\rho_{ij}Q_{ij}(r)\f$, and the AUGMENTATION
  functions \f$Q_{ij}\f$ are NOT positive-definite — localized, sign-changing objects that restore the true
  norm in the core.  So a USPP/PAW density dips negative near the nuclei with a perfectly PSD D: the
  negative ρ convicts the AUGMENTATION, not the density matrix.
  **Does not apply here:** this tree is GTH/HGH, i.e. NORM-CONSERVING (`Pseudopotential/GTH_Potentials.C`),
  which has no augmentation charge — \f$\rho=\sum_m f_m|\tilde\psi_m|^2\f$ exactly, so the premise HOLDS on
  the GPW path.  (The `augmentation` machinery in the tree is `APW_IBS`/`LAPW_IBS`, the all-electron
  augmented-plane-wave family, which does not feed `DM_RhoAtPoints`.)  ⚠ It would come back if USPP/PAW
  were ever added, or if the APW/LAPW route were pointed at this seam.
  **The governing property is CONVEXITY, and it explains the user's other observation — "for atoms and
  molecules doing HF this is impossible":** there D = C f C† with INTEGER occupations (PSD), and linear
  mixing \f$(1-\alpha)\rho_{in}+\alpha\rho_{out}\f$ with \f$\alpha\in[0,1]\f$ is a CONVEX combination —
  non-negative weights — so it cannot manufacture a negative eigenvalue.  PSD is preserved and ρ≥0 is
  guaranteed.  The property dies the moment a scheme EXTRAPOLATES instead of interpolating: Pulay/DIIS
  coefficients sum to 1 but individual \f$c_k\f$ may be NEGATIVE (affine, not convex, so the result leaves
  the PSD cone); likewise α>1, cross-geometry density extrapolation, or Methfessel-Paxton / cold smearing,
  which produce genuinely negative OCCUPATIONS (Fermi smearing cannot — f∈(0,1)).

  **IN THIS TREE, D IS PSD TODAY — verified by reading the mixers, not assumed:**
  - `LinearMixer` is the ONLY thing that touches D, and its `itsRelax` is clamped to `itsRelMax`, which
    starts at 1.0 and only ever DECREASES (to 0.5 near convergence) — so α∈[0,1], convex, PSD-preserving.
  - the EXTRAPOLATING mixers never touch D: `DensityMixer.C:159` — *"A mixer whose subject is a `GField` —
    i.e. one that works purely in G space (Kerker, Pulay; NOT the …)"*, and `:193` names that the "HISTORY
    face: a field mixer whose step is an EXTRAPOLATION over past iterations (Pulay/Broyden)".  The one
    operation that can leave the convex hull acts on \f$\tilde\rho\f$, not on D.
  - which is also why the two negative-ρ mechanisms stay separated here: the G-space extrapolation can drive
    \f$\tilde\rho\f$ negative (the documented ringing — `LDA_XC.C:163`, `PWTerms.C:392`, `GPW_IBS.C:53`),
    while D itself stays in the PSD cone.  `DM_RhoAtPoints` is pointwise with no band limit and does not ring.

  ⚠ **THE GUARD THIS IMPLIES.**  Cholesky is applicable *because* nothing currently extrapolates D — a
  property of today's mixer set, not a theorem.  It would break the day a DENSITY-SPACE Pulay lands (the
  `"Pul" (density-DIIS/Pulay)` tag at `DensityMixer.C:75` suggests one is contemplated), or an α>1 grow, or
  MP/cold smearing.  **The canary is free: a pivoted Cholesky that FAILS is exactly the signal that D left
  the cone** — so make that failure loud and route it to the eigen split rather than letting it degrade
  silently.  (And per the narrowed pin above, the eigen split is fine HERE because the factor is only
  multiplied, never inverted.)
  **★★ MEASURED 2026-08-20 (`GPW_DM_RANK=1`, new + permanent) — BOTH QUESTIONS ANSWERED, AND THE WIN IS REAL.**
  On the MnO benchmark recipe, per spin channel, over four iterations:

  | | measured |
  |---|---|
  | n | 118 |
  | **numerical rank of D** | **14–17** |
  | **⇒ ρ-GEMM speedup n/r** | **7.0–8.4×** |
  | rank stability over tol 1e-6 … 1e-12 | 14–17 (a CLEAN GAP, not a judgement call) |
  | \f$\lambda_{\min}/\lambda_{\max}\f$ | **−1e-16 ⇒ PSD to roundoff** |

  So **D is PSD** — the convexity argument above confirmed empirically, and pivoted Cholesky applies — and
  the rank sits just above the 13 occupied per spin, exactly as Fermi smearing should give (a few
  fractionally-occupied states on top).  The tolerance-independence is the important half: the occupied
  block is cleanly separated from the null space, so the truncation is unambiguous.
  **⛔ BUT THE BUCKET WAS MISIDENTIFIED — BY ME, FOLLOWING `doc/Benchmark.md` (corrected 2026-08-20).**
  IMPLEMENTED (guarded, exact, A/B at identical energy) and then measured on the production recipe:

  | bucket, 6 iterations serial | s |
  |---|---|
  | `scf: XC-mesh ρ sampling (matrix-free density)` | **35.0** |
  | `scf: XC-mesh quadrature H_xc` | 4.69 |
  | `scf: XC-mesh ρ sampling (all iterations)` — **THE DM GEMM** | **1.70** |

  `doc/Benchmark.md` called the 84.1 s top bucket "the Φ-shaped GEMM, i.e. the Φ-SPARSITY item".  It is
  NEITHER.  ***Matrix-free* means "carries no density matrix", so it cannot be the DM GEMM** — per
  `PWTerms.C:698` it is the **ρ̃-MIXED density sampled on the XC mesh, a batched inverse FT over the whole
  {G}, on EVERY Kerker/Pulay iteration**.  On these recipes the mixer hands XC a ρ̃-backed density from
  iteration 1 on, so **the DM GEMM is nearly BYPASSED**.  The code had already split those two buckets for
  exactly this reason (*"lumping it into the GEMM hid the fact that the mixed-density sampling, not the
  GEMM, was the iteration's largest XC cost"*) and the benchmark doc re-lumped them in prose.
  **Consequences.**  (1) The low-rank route is kept — it is exact, guarded, 7–8× on the GEMM, and that GEMM
  IS the hot path for DM-backed routes (molecular/atomic DFT, non-ρ̃-mixed recipes), growing with n — but it
  is **NOT a benchmark-row win and must not be quoted as one**.  (2) The Φ-SPARSITY refutation above stands
  on its own evidence (Φ really is ~48% dense) but its "targets the top bucket" framing was wrong for the
  same reason.  (3) **★ THE ACTUAL PER-ITERATION LEVER IS THE ρ̃-MIXED SAMPLING** — an inverse FT over {G}
  evaluated at ~97k mesh points every iteration — and nothing has looked at it.  That is the next item, and
  the same discipline applies: attribute the cost inside that bucket before optimising it.
  ⚠ Instrument note: the PSD test must be on a RELATIVE floor.  The first cut tested \f$\lambda_{\min}<0\f$
  and screamed "NOT PSD" on every run at \f$\lambda_{\min}\approx-1.8\times10^{-15}\f$ against
  \f$\lambda_{\max}=13.2\f$ — one ulp.  An eigensolver ALWAYS returns O(eps·λmax) negatives for a
  numerically-PSD matrix; a real negative eigenvalue from an extrapolated D would sit orders of magnitude
  above that floor.
  *Not a clever trick:* evaluating ρ in ORBITAL space rather than DM space is what most codes do; the only
  observation here is that this one took the DM route, which is a reasonable choice that happens to cost n/r.
- **★ ✅ FIXED 2026-08-25 — THE IMPOSED XC MESH LOST ITS SITE BLOCKS, and the Step 0a site-moment
  instrument was SILENTLY DEAD on every imposed run.**  *Cause:* the orbit-consistency filter in
  `UnitCell::CreateIntegrationMesh(imposed, Becke)` rebuilt the mesh by walking ORBITS, which
  interleaves the atoms' points and destroys the contiguous per-centre blocks; it now walks the
  ORIGINAL site order with a keep-mask (point ORDER is free — the caller rebuilds its fold FROM the
  finished mesh — the BLOCKS are not).  `RequireSiteBlocks()` now **throws** on both Becke arms so the
  invariant survives NDEBUG, and `EmitSiteMoments` announces its own silence instead of returning
  quietly.  *Verified:* Etot −61.40297621 to 10 s.f., identical stage counts (14+17), identical mesh
  (97160 pts / 4220 orbits / 414 dropped) — and **29 `[site moments]` reports where there were zero**,
  reading **Mn ±4.78 e** against CP2K's Mulliken ±4.654, O at 4e-08, net 3.5e-08.  The original
  description follows.  `MakePeriodicBeckeMesh` records
  one block per atom (`BeginSite`), which is what makes \f$\int w_A m\f$ a real integral instead of a point
  sample.  The imposed/invariant mesh path (symmetrise + `FoldPointsPeriodic`) drops them: the engine's
  mesh reports `NSites()==0`.  `qcMesh::SiteIntegrals` guards this with an `assert`, which is COMPILED OUT
  under `NDEBUG` — so a Release imposed run does not fail, it just prints no Becke site moment at all
  (verified: `QCHEM_SITE_MOMENTS=1 MNO_IMPOSE=1` emits only the legacy point probe).  This matters twice
  over: **every benchmark row is `MNO_IMPOSE=1`**, and Step 5's "re-read the campaign's moment conclusions
  against the integrated observable" would have silently produced nothing on exactly the runs it targets.
  Fix = carry the site blocks through the invariant-mesh construction (each orbit rep keeps its owning
  site), and turn that `assert` into something that survives Release.
- ~~**Becke mesh build** (~12–32 s): `BeckeCutoff` alone was 11.5% of the round-3 profile.~~ → covered by
  the two Becke items above.
- **NOT worth re-attacking without new information:** the per-iteration scatter/gather.  Round 4 measured
  it at 61% irreducible per-(pair,point) emit; it needs fewer TERMS (looser `GPW_DENSITY_EPS`, the T3
  fold, a smaller span), not a faster kernel.

### ✅ Step 4 — RAM — LARGELY ANSWERED BY STEP 2  →  `doc/OpenWork_History3.md`

Arming the T3 fold took MnO AFM peak RSS 4947 → 1349 MB (23× CP2K → 6.2×).  Do not open it as a track;
re-read the number off `doc/Benchmark.md` and reopen only if it still binds.

### Step 5 — MnO ACCURACY: NAME THE OPERATOR  ·  plan: `doc/SphericalLatticePlan.md`

**This has the sharpest coordinate on the list and a banked oracle — and its first move is CHEAP.**  The
VA (N=118) exact-span table, both codes at full rank, zero span/symmetry/ensemble excuses:

| | AFM | FM | Δ(AFM−FM) |
|---|---|---|---|
| qchem VA | −61.40298 | −61.44158 | +38.6 mHa (FM wins) |
| CP2K VA | −61.30333 | −61.30478 | +1.46 mHa (FM wins) |
| offset | **−99.7 mHa** | **−136.8 mHa** | **config-selective −37.2 mHa** |

So qchem carries (a) a ~100 mHa **configuration-BLIND** offset — and it sits *below* a variational
reference on an identical span, which convicts an operator or convention, not a basis (suspects: G=0 /
alignment conventions, XC quadrature, V_loc bookkeeping); and (b) a ~37 mHa **configuration-SELECTIVE,
SPAN-INDEPENDENT** FM-favouring bias (v2 span: ~45 mHa; the I0 d-selective signature).  Run 61 also
converged in 14 cold iterations at `m_stag ±0.667` — **the moment question and the energy question are the
same investigation on the same run.**

**★ RE-READ THE CAMPAIGN'S MOMENT CONCLUSIONS FIRST (new, 2026-08-19).**  Step 0a measured the integrated
Mn moment at ~3.6–4.8 e against the point probe's ~0.3, so every "weak-moment basin" / "the moment died"
conclusion in `doc/SymmetryUpgradePlan.md`, `doc/SphericalLatticePlan.md` and
`doc/SymmetryUpgradeHistory.md` needs re-reading against the honest observable before more physics is
built on it.  Collapse-to-zero findings survive (zero is zero); MAGNITUDE and site-ASYMMETRY findings do
not automatically.  Run 61's `m_stag ±0.667` in particular is a point-probe number and is NOT evidence
that VA sits in a weak-moment basin.  This is cheap — rerun the banked recipes with
`QCHEM_SITE_MOMENTS=1` — and it may re-scope the whole moment half of this step.

**Next concrete move on the ENERGY half (cheap, does NOT wait on Steps 1–4):** the term-by-term breakdown —
Ekin/Eee/Exc/E_loc/E_NL against CP2K's energy blocks, with `Een` split into V_loc/V_nl — run first on the
**Mn ATOM** (seconds, and the Mn-atom oracle is already banked at −14.2440 restricted / −14.658 atomic
polarized), then on the crystal.  A configuration-blind ~100 mHa offset should be visible on a single atom.
Then `MNO_KMESH=2` (k-convergence moves the ordering MORE than the physical 6J₁+12J₂ ≈ 4 mHa scale) and +U.

### Step 6 — THE 136-FUNCTION SPAN: a capability gap, no longer a blocker  ·  time-boxed research

Why CP2K holds the full diffuse span and qchem must strip it.  **Demoted from "blocks the comparison" to
"capability gap"**, because VA/VB give exact-span comparisons today (Step 1).

⛔ The **screen-discipline hypothesis is REFUTED** (run 64, `doc/logs/mno_probe_run64_fullspan_tighteps.log`):
132-of-136 span at eps 1e-12 both screens still dove to −82.19, i.e. 20.7 Ha BELOW CP2K's variational
−61.47 on a SUBSET of CP2K's own span.  Tighter eps is not the cure.  Surviving candidates:

1. SVD/eigen-consistent F/S filtering (the run-58 note, parked);
2. the **symmetry-INEQUIVARIANT AO drop** — run 64 dropped indices 11/9/115/1 individually; see the
   vet-stage item under *Continuous* below, which is a near-prerequisite;
3. projecting the near-null directions out of **F** as well as **S**, rather than screening around them.

Re-run note: set `GPW_MNO_VERBOSE=1` (`GPW_REPORT=1` gives the ledger but NOT the per-iteration table, so
run 64 cannot show whether it dove or stalled).

### Continuous — CLEANUP (a rhythm between the steps, not a phase)

- `doc/CleanupCandidates.md` D1–D13, plus **V1.32** (de-template the finite `IrrepCD` leaf).
- **★ PULLED FORWARD out of the cleanup bucket — the VET-STAGE symmetric basis trim**, because it is a
  live candidate for Step 6.  Promote basis trimming from an ortho-time knob to a first-class VET-stage
  policy: auto-select the kept sub-basis (pivoted Cholesky already reports which functions are redundant;
  auto gap tolerance exists in `LASolverLapack.C detect_null_gap`), report it as a BASIS decision
  (species/shell/exponent, not bare indices), perhaps regenerating it via `qchem.ValenceBasisGen` so the
  user sees a basis, not a filter.  Three pins, all from the user:
  - **NOT display-only** (2026-08-14): the trim must happen BEFORE anything downstream is built — grid-ladder
    depth, collocation streams and their caches, KB projections, the whole per-pair machinery all fall out
    of the surviving function list, so filtering at ortho time does the dropped functions' work for nothing
    AND leaves their pairs in every stream and cache.
  - The rank decision is a property of **S**, i.e. of the BASIS, made ONCE — today each spin channel's
    LASolver re-derives it independently (the doubled `[ortho]` line), coherent only because S is
    channel-independent.
  - **SYMMETRY-EQUIVARIANT** (2026-08-15): drop whole ORBITS under the (magnetic) space group, never
    individual AOs.  Greedy per-function pivoting resolves symmetry-TIED pivots by numerical noise (runs
    58–60 dropped O₁'s p(0.18) but O₂'s s(0.15), and different d(0.18) m-components on the two Mn),
    breaking site equivalence at the same ~sub-1% order as run 59's site-moment asymmetry.
    Direction-space (eigen/SVD canonical-ortho) trimming stays a numerical fallback, never user-facing policy.
  - **★ AND THE SYMMETRY PIN IS ALSO A CONVERGENCE PIN** (user, 2026-08-26): pivot-filtering *"would
    sometimes remove only 3/4"* of a symmetry-equivalent set, and hand-removing functions BEFORE S vetting
    has proved more reliable for convergence than letting the filter choose.  That is the SAME mechanism
    as the equivariance pin above, not a separate argument -- a partial orbit is a symmetry-broken basis,
    and what it costs is not only site-equivalence in the reported numbers but the run's ability to
    converge at all.  ⇒ The vet-stage trim is not a tidier way to do what ortho-time pivoting does; it is
    a DIFFERENT and better-behaved thing, and ortho-time filtering should end up as the fallback that
    fires when the vetted basis was still not good enough.
    ⚠ **The partial-orbit case is the one to reproduce** -- when the filter lands on a whole orbit there
    is by construction nothing to see.  (Measured 2026-08-26 in passing: on the MnO cell `VALENCE_LOWQ_SR`
    is 122 hand-trimmed at full rank and `VALENCE_LOWQ_VA` under Cartesian d is 132 auto-dropped to 122,
    both at min kept pivot 0.0236681 -- the same kept set two ways, so THAT pair is a null control, not a
    test.  Runs 58-60 above are the real evidence: O1's p(0.18) dropped against O2's s(0.15).)
- Δρ/N convergence gate (`doc/SCFStrategyPlan.md`); GDM fallback-diagonalize breadcrumb (run 59's silent
  +302 mHa hop); the per-channel ortho duplication above; the fingerprint's overconfident
  "raise NMaxIter" advice.

---

---

## EVIDENCE DOSSIER (no action here) — fast ρ by FACTORING D  →  `doc/OpenWork_History3.md`

The worked evidence for INDEX ITEM 2 (Step 3's low-rank-D ρ GEMM): Q1/Q2/Q3, the tier-0 results, the
spectrum finding, the LSP design ruling.  **There is nothing to start here.**  Read it in the history
file when you build that item, or when you are about to re-propose something it already refuted.

## ✅ CLOSED 2026-08-25 → `doc/OpenWork_History2.md`

The 2026-08-21 → 08-24 arc moved out whole: **the Vxc repair thread** (routes a/b/c, α_eff, the Kerker
residual spectrum, the GDM canary), **the XC separation-of-concerns design item** (four steps + the
absorption), **the fitting boundary**, and **one fit-basis interface**.  All done and verified; the durable
rulings live in `doc/CleanupCandidates.md` R1.0.  Two live consequences were LIFTED out of that record
rather than filed with it — they are in Step 3's runtime item above: route (b) is MORE accurate than ρ̃
(not neutral), and Step 5 must pin the flag before its term-by-term CP2K comparison means anything.

One thing from that arc is **specced and NOT built**, so it stays here:

### ★ SPECCED, NOT BUILT (user, 2026-08-23) — SEPARATE THE METRIC AXIS INTO FACES

`OverlapDiagonal()` sits on the metric-NEUTRAL fit face, so a basis that has no diagonal metric must invent
an answer — and `Fit_IBS`'s invented one is in a different normalisation from every other member of its own
face.  Move it to a new `FIT_SF_Ortho<T>` (mirror `FIT_CD_*`; **both sides in the same increment**, user)
and `Fit_IBS` simply loses it: the landmine is deleted, not corrected.  No orthonormal marker face —
orthonormal is orthogonal with a unit diagonal, and a memberless face is the null-object pattern already
rejected twice.

⚠ **ACCEPTANCE CRITERION, and it is the point of the item (user):** removing `isOrtho()` must NOT be
replaced by `if (dynamic_cast<FIT_SF_NonOrtho*>(fbs)) … else …` — a type switch wearing a cast is worse
than the bool.  Measured and reassuring: **all eight `isOrtho()` call sites today are `assert`s, zero live
branches**, so there is nothing to replace; narrowing a PARAMETER type (`OrthogonalFit(const
FIT_SF_Ortho<T>&)`) is the sanctioned substitute.  The one real branch in the tree — `Factory`'s δ-vs-PW
`dynamic_pointer_cast` — is a REPRESENTATION branch at a creation boundary, not a metric one, and must be
left alone rather than laundered into a metric test.  Full spec, including the open question that would
delete even that branch: `doc/CleanupCandidates.md` R1.0.

---



## Parked threads (real work, deliberately not in the plan above)

- **A. Spherical SALC — S3b, the libcint-spherical extractor**  ·  `doc/OldPlans/SphericalSALCPlan.md`.  The ONLY
  remaining piece; S1–S5 are done and the in-house spherical SALC is fully shippable without it.  The
  bug-prone one: it must match **libcint's** real-harmonic ordering + normalization (a foreign convention),
  and libcint-spherical presents AS a `PGData` with spherical components (a trap).  Genuinely separable.
- **PBE / GGA**  ·  `doc/FacadeDFTPlan.md`.  The highest-value functional for the battery north-star, but a
  real library increment (density-gradient machinery on the mesh), not an enum value.  The unified `Model`
  enum is ready to list it with a "not wired" throw.
- **LibXC-polarized** (the wrapper needs two-channel `xc_lda_vxc`) and **+U**.  The LDA *interface* is
  spin-native end to end already (`doc/SpinNativeDFTPlan.md`, closed).

## Deferred & descoped — recorded so they are not re-litigated

- Fold `QchemTester` + the pybind bridge onto the facade — test-harness/binding cleanup, not lib surface.
- Container utils to `src/` (`sample_scalar`/`sample_gradient`, `Structure::BoundingBox`) — binding convenience.
- `SCFParams` ASCII rename — DROPPED; solved by C++20 designated initializers (`34ccf302`).
- `MolecularSym_EC` → `FixedIrrepOcc_EC` rename — belongs with the queued symmetry-naming cleanup.

---

## Why this order (the reasoning, so it can be argued with)

1. **Instruments before measurements, measurements before optimisation.**  Rounds 3 and 4 each narrowed a
   charter premise that had been plausible for months, and both failed the same way — cost was attributed
   to the thing that was easiest to name.  Step 1 exists so Steps 2–4 have a finish line.
2. **Step 0 is cheap and blocks Step 1.**  A moment instrument in the wrong units and invisible fold
   factors would both be baked into the benchmark table.
3. **Step 2 before Step 3** because folding is 12–48× where kernel tuning was 1.4×, and because it changes
   what the profile says is hot.
4. **Step 4 is mostly a consequence**, not a track: RAM falls out of folding + re-tiering.
5. **Step 5 does not have to wait.**  Its first move is a term-by-term energy comparison on a single Mn
   atom against a banked oracle — minutes, not overnight.  It was parked for a cost that its cheapest
   experiment does not have.
6. **Step 6 is demoted** because VA/VB spans make the code-vs-code comparison possible today; it is a
   capability gap, and its most testable remaining candidate (the inequivariant AO drop) is the vet-stage
   cleanup item.


