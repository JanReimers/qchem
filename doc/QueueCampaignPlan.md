# Queue campaign — clear the critical path before forces (DRAFT v2 2026-10-29, user comments on v1 folded in)

Born at the top of `doc/` because it is about to be executed; retires to `doc/Records/` when the forces design note starts.
Scope: `doc/OpenWork.md` §2–§3, `doc/CleanCode.md`, `doc/OOD-SOLID-Cleanup.md`, `doc/ModuleToolchainPlan.md`.
The HubbardUPlan / KAxisPlan sequence continues in its own lane (H) and is not re-planned here.

## 1. The finish line: the CRITICAL PATH is clear, not the queues empty

Most new rows are born from campaign work (β-MnO₂, NiO, the PW route): every measurement finds a corner, so a zero-row goal is
unreachable while campaigns run.  The finish line is instead: **no open row tagged F, V or S.**

| tag | meaning | counts? |
|---|---|---|
| **F** | blocks or corrupts forces / geometry optimization (grid aliasing, SCF robustness, restart after a geometry step, basis that changes with geometry) | yes |
| **V** | blocks the spinel voltage curve (Li / PP / basis plumbing, conditioning of Li cells, GGA, the functional recorded truthfully) | yes |
| **S** | structural debt IN THE CODE THE FEATURES WILL TOUCH: Hamiltonian terms (V1.41, V1.35, R1.0q, V1.34, §I.1, V1.39), SCF iterator / accelerators / mixers (D9, V-MIXER-GBALL, V1.27, V1.24(i), D-DIIS-FDMIN), the Lattice basis + GPW evaluator (D-GPWPAIR, D3, D4, D13, D-ENV-TOL), the facade (D-ENV typed fields, D-PW-FACADE, and **V-CALCNET — see below**), plus D-WARNINGS (it finds bugs) | yes |
| **H** | owned by the Hubbard / KAxis lane | no |
| **X** | hygiene far from the feature path (R2.24, R2.23, R2.21, D10, R1.0f, D-RK-ONCE, D-CONTRACTED-SHELL, D-POINTPROBE, D-TIMING-GATE, D-XCROUTES-GROWTH) | no — but see the quota |
| **P** | parked / trigger-only / a feature, not debt | no |

**The compromise on clean code (user: "the cleaner the code, the easier anything big").**  v1 put all OOD/clean rows on the back
burner; that was wrong.  Two rules instead:
1. **Clean before you build**: every S row in the files a feature will touch closes before that feature starts.  That is why S counts.
2. **A standing hygiene quota**: lane B spends at least one chunk in three on X rows, so the far-away debt also shrinks, at a known rate.

**V-CALCNET, re-opened (user 2026-10-29: "a design review any time you deem appropriate")**: a geometry optimizer is a NEW CLIENT
that drives a Calculation many times at moving positions — `MoleculeCalculation` AND `SolidCalculation` (user), so it needs one face
over both — keeping the basis, the PP set, the grids and the previous density.  That is the network's placement question, and so is
C4's "where is PP data read" (D-STRUCTDATA step 2).  So the review runs FIRST in lane B, before C4 (§5).

The controller reports the **open F+V+S count** at every merge.

## 2. What the zoom-out found

- **25 ✅ rows still sit in `OpenWork.md`** (4 in §2, ~21 in §3) against the same-day rule; with their REMAINDERs about a third of the
  tracker's lines.  A remainder that is not F/V/S becomes a one-line CleanCode row or is dropped.
- **§2 rows that are not features**: PW-1 remainder (= OOD `D-PW-FACADE`), "Dense PW cost" + "LR efficiency ledger" + "LR on REAL
  TRIM blocks" (→ §4), "CP2K finite-α LRT" automation (→ CleanCode), "Second oracle" (a rule, in CLAUDE.md → delete), "Diffuse-basis
  ACTUATOR" (a defect with a ruling, pin 22 → §3).
- **Duplicates**: §3 "Seed library keyed by q" ⊂ D-STRUCTDATA step 2; §3 "`u_bands` after a GDM stage" is very likely the
  GDM-virtuals defect fixed in `7fb46357` (re-measure, then close).
- **Already done since v1 (checked in git 2026-10-29)**: [F,D] has one home + the SCF cap guard (`1e0b2069`, `2b7fc2ca`, merged
  `48f2e02f`); the converge-once `Calculation` constructor (D-CALC-DOUBLE-SCF's core); the linear mixer's irrational grow ratio +
  `MaxRelaxRo` (part of "Linear D-mixing never re-damps").  V1.41 stage 1 is on main (`78ff59de`); stages 2–3 are open.
- **Forces concerns not yet written down** (go into the forces design note): (a) pin 22's vet trim depends on λ_min(S(R)) — along a
  relaxation the basis must be FROZEN or the PES is discontinuous; (b) GPW aliasing is an egg-box force error; (c) a geometry step is a
  restart from a moved density (the Ladder row's unmeasured remainder; CK-3); (d) `H.At(cd)` is the natural carrier of dE/dR at a
  converged D; (e) R1.0r (§7) — a density symmetrized under more symmetry than the k-sample has is not variational, so its forces are not
  the derivative of its energy.
- **GGA goes before forces** (agreed): forces design note → GGA → forces → relax → voltage.

## 3. Convergence — stop and design (user 2026-10-29)

**The evidence.**  Most of the last month's lost time is SCF convergence, and the causes are a short list: pin 28 violations (XC fed a
fitted ρ); the BallOnly raster where AliasFree was needed; restarts (stale history, refill, the GDM rung); imposed symmetry that leaked
(the Shubnikov flip ops projected only in XC); GDM (stale virtuals, false convergence on a changed Hamiltonian, the crawl); and a ρ
mixer and a Fock accelerator each keeping a HISTORY of a different object.  And the β-MnO₂ discrepancy hunt is now stalled on convergence too.

**What the list has in common.**  None was a missing algorithm.  Each was an ASSUMPTION one component made about another that nothing
checked, found late (often hundreds of iterations, or days, after it went wrong) and diagnosed by elimination.  So the cure is policy +
instrumentation, not another mixer.  Five candidate policies (C1's design note decides):

| # | policy | today | what it would have saved |
|---|---|---|---|
| **P1** | **One history per run.**  Exactly one object extrapolates (DIIS on F, Pulay on ρ or D, or GDM's own); a ladder HANDS OFF (re-seat, drop the other's history) and says so.  Illegal pairs are unrepresentable at the facade, not warned about.  This is `Records/SCFStrategyPlan.md`'s "one shared extrapolator" — designed, not finished (OOD D9) | Kerker/Pulay + DIIS can both run | the competing-histories episodes, GDM's stale-mixer fallback |
| **P2** | **Invariants are measured every iteration and fail at iteration 1, not 300.**  A small fixed set, cheap, silent unless violated: the imposed group's residual ‖P_G D − D‖ (and on m), the XC ρ source (pin 28 — already a throw), the reference eigenpair residual, charge, aufbau holes.  A violation names the object that broke it | some exist as columns or one-off throws | the Shubnikov leak (found at iteration 300), pin 28, stale virtuals |
| **P3** | **Precision contracts.**  A consumer DECLARES the accuracy it needs (littleGroupK: orbital equivariance 1e-6; a response: reference Δρ 1e-9; a force: …); a producer REPORTS what it delivers (the raster's alias estimate — `[grid audit]` already measures it; the SCF stop level; the Becke mesh).  Pre-flight compares and refuses with the cause, or escalates the producer (BallOnly → AliasFree) on its own.  The response-reference refusal (`a2e664f1`) is the first instance; generalize it | ad hoc, one case | **the AliasFree day**: speed stays the default, the one consumer that needs AliasFree gets it, and nobody spends hours finding out |
| **P4** | **Named recipes, gated by a convergence PANEL.**  Campaign decks pick a named recipe (`Auto`, `MagneticOxide`, `Restart`, …), not free knob combinations.  A panel of the known-hard cases (MnO AFM-II free, β-MnO₂ A1AFM at U = 1, the Mn–O dimer, Na₂ AFM, NiO, O₂ box, the Si +U loop, exact and coarse→fine restarts, a TiO₂ little-group reference) runs as a ctest label with ITERATION BUDGETS: a change that helps one case and hurts another goes red before a campaign finds it | hard cases live in campaign logs | the "fix one, break another" cycle |
| **P4a** | **The regime matrix** (user 2026-10-29: "best practices by kT = 0 / kT > 0, metal / insulator").  One table: rows = regime (insulator kT = 0; magnetic insulator with a near-degenerate d frontier; metal / smeared kT > 0; restart, same Hamiltonian; carried stage, changed Hamiltonian (+U loop, geometry step)), columns = accelerator, mixer, preconditioner (Kerker for metals' long-wavelength sloshing, none for insulators), occupation policy, stop gate.  Each cell cites its panel case, or says UNMEASURED.  Seeds already banked: pin 15 (GDM is fixed-occupation: kT = 0 only), `Auto` = Pulay on D for seeded stages and Ladder@GDM for restores (2026-10-26), Kerker + Pulay for free rocksalt TMOs (CLAUDE.md), plain Kerker holds β-MnO₂'s U = 1 basin where Kerker + Pulay overshoots (e9765e33's row).  The matrix IS the `Auto` policy's spec (the code chooses, the table documents); a cell that was earned by a wrong number becomes a pin | scattered over CLAUDE.md, pin 15, OpenWork rows | choosing a recipe by trial and error per campaign deck |
| **P5** | **The stop rule for campaign lanes.**  A campaign session that spends more than ~2 h on one run's convergence STOPS: it files a row with the deck and the trace, and the run joins the panel.  Convergence is debugged in the convergence lane, on the panel, never by knob-tuning one campaign deck | — | β-MnO₂ stalls; keeps H on physics |

## 4. The chunks

Each chunk = one branch, one merge, one sweep + `scripts/anchordiff` against the main baseline.  "Touches" is for conflict planning.

| # | chunk | rows | tag | touches | anchors move? | model · effort |
|---|---|---|---|---|---|---|
| **C0** | Tracker hygiene | 25 ✅ rows → History5; §2 misfits moved; duplicates merged; every open row tagged F/V/S/H/X/P | — | `doc/` | no | controller |
| **C1** | **SCF robustness: design, then build** (replaces v1's C1+C2) | step 1 = a design note `doc/SCFRobustnessPlan.md`: classify every convergence incident in History5 / HubbardUHistory since 2026-09 by root cause, rule on P1–P5.  Then build in the order ruled: the panel (P4) FIRST so every later step is measured; P2 invariants; P1 (D9, V-MIXER-GBALL); then the open rows on the panel: "GDM crawls", "GDM false-converges", "Linear D-mixing" remainder, "Na₂ AFM seed", "Polarized Γ Si jumps"; D-DIIS-FDMIN, V1.24(i), V1.27, D-CALC-DOUBLE-SCF migration | F/S | `src/SCFAccelerator`, `src/SCFIterator`, mixers, `Calculation` | iteration counts | Opus · high (design); Opus · medium (build) |
| **C2** | Grid anchor window (re-bank ONCE) | design note first (§7): one accuracy target → every grid tolerance derived; the raster N sized from ε (`BallOnly`/`AliasFree` retired); V_loc long analytic or folded into Hartree; §K's one ball for ρ and v_xc + the XC-raster audit.  Then: AUTO density cutoff folds the PP sharpness; uniform-raster XC not equivariant; free-run Becke mesh cost lever; R1.0r (§7, ruled) | F | `GPW_Evaluator`, `Lattice`, mesh code | every GPW anchor — judged vs CP2K (pin 10) | Opus · high; holds the heavy slot |
| **C3** | Stateless Hamiltonian | V1.41 stages 2–3, then the V1.35 plan; §I.1 + V1.39 with GGA; exit deliverable = **the forces design note** | F/S | `src/Hamiltonian` terms (incl. `Hubbard_U` — coordinate with H) | no (bit-identical, staged) | Opus · high |
| **C4** | Inputs & provenance | D-STRUCTDATA step 2 + seed keyed by `q` + D-BASIS-PP + D-BASIS-REGISTRY + "valgen --seed-into"; HubbardUPlan B4 (Li q3, with the basis-vs-PP check); restart fingerprint (functional hard-coded "LDA" → recorded, WARM; consequences printed); PW `state.save` refused up front; D-ENV step 5 typed fields + D-ENV-TOL | V/S | `Calculation/Imp`, seeds, `BasisSet` factory, `Deck.C`, `valgen` | no (expect 0 moved) | Sonnet · high |
| **C5** | Basis conditioning | "Diffuse-basis ACTUATOR" (`Prune`, pin 22); "vet trim deletes the valence core of a densified ladder"; "valgen is optimistic about diffuse functions" (`--env`) — Li is the customer | V/F | `BasisSet/Lattice` vet, `valgen` | possibly | Opus · high |
| **C6** | Warnings + hygiene quota | D-WARNINGS first (count, then `[[clang::lifetimebound]]`); then X rows, one chunk in three | S/X | scattered | no | Sonnet · medium |
| **C7** | Toolchain, the safe half | ModuleToolchain Step 0 (Clang 22) + Step 1 spike + Step 5 gate (1) (does Boost.UT fail on an empty selection?) | X | CMake, one lib | no | Sonnet · high |

**Not in this campaign (P)**: ModuleToolchain Steps 2–5 (whole-tree: a FREEZE window with no other lane open, between this campaign
and forces, if at all before forces); V1.24(iii) (needs OT); D-ENV-RERUN stage 2 (H, last); D-CLOUD (your account); NiO KB ghost,
MnO span / ordering, charged-cell convention (H or P); every §2 feature row.  The Blaze heevx fix waits for Step 4 as ruled.

## 5. Sequencing

- **Lane A (critical path, Opus)**: C1 → C2 → C3 → forces design note.  C1's design step comes first because P3 decides HOW C2's raster
  default lands, and P4's panel is the instrument for everything after it.
- **Lane B (V + S + quota)**: V-CALCNET design review (doc only, Opus · high) → C4 → C6 → C5 → C7, X rows between.
- **Lane H** continues HubbardUPlan / KAxisPlan.  Under P5 it hands any convergence stall to lane A's panel instead of tuning.

H, A and B run concurrently; the controller, not the user, tracks them (user 2026-10-29: "2" was what a human can follow, not a
session cap).  The real cap is the BOX: two module builds at once plus one heavy job fills 14 GB, so a fourth lane is worth it only
if it mostly edits rather than builds.  The user talks to the controller; lanes talk to the controller.  C3 starts only after
KAxisPlan P5 lands, since H owns V1.41 so far.

**Lane G (GUI, ruled 2026-10-30, starts when lane B's steps 6–8 land):** worktree `~/Code/qchem-lane-g` on `lane-g` + a clone of `github.com/JanReimers/qchem-GUI` at `~/Code/qchem-GUI`; Sonnet.  Owns `pybind/` (the binding owner the lib lanes flag breaks to).  Step 1: a FRESH bridge on main (not a merge of the stale `origin/viz`, 1512 behind): bind only `MoleculeSetup::Resolve` + `Calculation` (later `MoleculeSession`), an opaque handle type (no `void*`), keep the module-unit / plain-C-header split; then make qchem-GUI's `QChemBackend` run on it, checking first whether `stream_scf` needs a per-iteration hook on the facade.  Off the critical-path count; build under the build lock.  **Requirements added (user 2026-10-30):** (1) an AI agent inside the app that helps build the structure and choose Hamiltonian, basis, U (incl. computing it), convergence criteria and methods; (2) a front end for SOLID², CP2K, QE, ABINIT and VASP -- **option 2.2 (user's pick): the user states WHAT they want, and the back ends that can do it are enabled** (not 2.1, features toggled by the chosen back end).  Consequence: the GUI's seam becomes a back-end-NEUTRAL run spec (structure, physics, task, ONE accuracy target) + a `Backend` capability/translate/run/parse interface; the agent's product is a spec, so order = spec → back ends (SOLID² + CP2K first) → agent.  Step 0 of lane G is that design note, in the qchem-GUI repo, on Opus; plumbing on Sonnet.

## 6. Controller protocol (this session's job)

**Working copies.**  The repo of record is `~/Code/qchem6`, and lane H works IN it (2026-10-29: committing on `main`).  Each new lane
gets a `git worktree` of qchem6 named for its branch — `~/Code/qchem-lane-a` (`lane-a`), `~/Code/qchem-lane-b` (`lane-b`), and
`~/Code/qchem-main` (`main`, the controller's: docs + fast-forwards, no build).  NB a worktree with submodules cannot be `git worktree
move`d: `mv` it, `git worktree repair`, then fix `core.worktree` in `.git/worktrees/<name>/modules/*/config` — so every
branch is visible in one repo and merging needs no fetch.  Each gets its own `build/Release`, configured WITH OpenMP (CLAUDE.md flags;
check `ldd … | grep omp`).  Disk is no issue (317 G free; a build tree is 4–9 G).
**H is on its own branch** (done 2026-10-29), so `main` is no one's working branch.  `main` then advances ONLY by fast-forward to a lane's swept merge commit (below), which needs no
build of `main` and no checkout of it in a working tree H is using.
**qchem1** (`openwork-sweep`): its work is merged into main (`48f2e02f`); it holds only `TODO` edits and this draft.  Retire it when convenient.

**Merging.**  At chunk START a lane merges `main`.  At chunk END: build `allTests`, take the heavy slot, `ctest -j4 -LE long` and
`-L long -LE "panel|periodic"` (serial), plus `-L "panel|periodic"` when the chunk touches the SCF path (CLAUDE.md, ruled 2026-10-30), keep `LastTest.log`, `scripts/anchordiff` against the controller's main baseline.  Because the lane merged main FIRST, the
swept commit IS the merge result: the controller fast-forwards main to it (`git merge --ff-only lane-a`), stores that sweep as the new
baseline `~/Code/qchem6-runs/anchors/main-<sha>.log`, and tells the other lanes to merge main.  If main moved in between, the lane
re-merges (a docs-only move needs no re-sweep — the controller's call).  Nothing is pushed anywhere.  H lands the same way, so every
baseline step has one named cause.

**Docs conflicts.**  Workers ADD a new row at the END of its table, close a row by moving it to the history file the same day, and
never reformat a table.  The controller resolves doc conflicts at merge time.

**Intake of new rows.**  At each merge the controller diffs the three queues (`git diff main...lane -- doc/OpenWork.md
doc/CleanCode.md doc/OOD-SOLID-Cleanup.md`), tags each new row, slots F/V/S rows into a pending chunk (or opens one), and reports the
open F+V+S count.  A new CONVERGENCE row also gets its case added to the panel (P4).

**Box resources (i7-10700: 8 cores / 16 threads, 14 GB + 32 G swap).**  RAM is the constraint, not cores.
- **One HEAVY slot**: a full ctest sweep, the long tier, the panel, or any calculation over ~2 GB or ~30 min, taken with
  `flock ~/Code/qchem6-runs/.heavy.lock scripts/memsafe <cmd>` — a second heavy job simply waits.  Overnight campaign runs hold it.
- **Builds**: ONE lane compiles at a time: `flock ~/Code/qchem6-runs/.build.lock scripts/memsafe -H 5G ninja -C build/Release -j4 <targets>`
  (2026-10-29: two lanes compiling at once beside a 5.7 GB PW run put the box at 12/14 GB with 8 GB swapped; a module TU peaks ~1.8 GB).
- **Threads**: a lane's targeted runs use `QCHEM_OPENMP_THREADS=4`; the heavy-slot holder gets 8 when alone; never two jobs at 8.
  `QCHEM_BLAS_THREADS` stays 1.  Every batch script exports its count.

**Models.**  Controller: Opus · medium.  Design and physics chunks (C1 design, C2, C3, C5): Opus · high.  C1 build: Opus · medium.
Plumbing and mechanical chunks (C4, C6, C7): Sonnet · high / medium.  No Haiku.  A Sonnet chunk that hits a physics surprise (an
anchor moves unexpectedly) is handed to Opus.

**No backward-compatibility doors** (user 2026-10-30): an API change migrates every caller and deletes the old entry point in the same chunk -- never "the old ctor delegates, so no caller changes".  `pybind/` is the one exception: flag a break for the binding owner.

**Sweeps that are compared run SERIAL** (`QCHEM_OPENMP_THREADS` unset), like every main baseline: several tests are not thread-count invariant (lane A 2026-10-29: 57 anchor moves up to 2.8e-6 at 4 threads; OpenWork row).  Threads are for panel and campaign runs, never a landing sweep.

**Landing sweep = short serial (`-LE long`) + routine long serial (`-L long -LE panel`, ~6 min); the convergence panel (`-L panel`, hours) runs only when a lane touches SCFAccelerator / SCFIterator / mixers / the Calculation SCF path, and is judged against its iteration budgets, not by anchordiff** (user 2026-10-30).  The threaded panel pass runs the LONG panel only (`ctest -L panel -L long`, both labels required); the fast `_Panel` cases are judged from the SERIAL short tier, because two of them (NaF exact resume 5 vs 3, O2Box anchor 7e-7) fail at 4 threads on one run and pass on another -- the thread-invariance row (2026-10-30, lane B gate).

## 7. The three grid-and-symmetry decisions, with context

**Raster default (agreed: AliasFree — with your caveat).**  BallOnly is much faster and was right for weeks; the loss was the HOURS
spent before anyone suspected the raster.  So the decision has two halves: (a) C2 MEASURES AliasFree's wall cost on the
`doc/Benchmark.md` §5a rows; (b) whatever the default, P3 makes the choice self-announcing: a consumer that needs equivariance or
force accuracy below the audit's alias estimate either escalates to AliasFree itself or refuses with the reason.  If (a) shows ≲ 1.3×
wall, AliasFree is the default and BallOnly a stated speed choice; if it is much dearer, BallOnly stays the default under (b).
Forces will be such a consumer either way.

**R1.0r — which group the density is star-averaged under.  RULED 2026-10-29 (user): the INTERSECTION.**  (C0 records the ruling on the
OOD row; C2 builds it.)  The density
is symmetrized under the FULL crystal group even when the k-MESH has less symmetry (Si 2×1×1: the mesh keeps 2 ops, ρ is averaged under
48).  A BZ sum over that mesh is only invariant under the mesh's ops, so today's fixed point is not the one that sampling defines: it
is arguably a nicer answer (the sampling error is symmetrized away), which is why it was left as a judgement.  VASP/QE use the
INTERSECTION (crystal ∩ mesh).  Why it matters now: with forces, the energy must be stationary in the density at the given sampling,
or the Hellmann–Feynman force is not dE/dR; and a relaxation that lowers the symmetry (Jahn–Teller in Li₂Mn₂O₄) needs the group to
follow the sampling and the structure, not a fixed big group.  Cost: small (`FoldGrid` already computes the mesh test; expose it,
intersect in `DetectPointOps`); Γ and isotropic meshes are unaffected, so nothing banked moves.  Only anisotropic meshes change.

**§K — the "perfectly tuned" fit balls.  RULED 2026-10-29 (user): yes, if a perfect XC fit can be defined like the AliasFree
ball — it can, to a stated tolerance, and it is the SAME ball.**  The reasoning, for C2 to verify against the code:
- *Two separate things set a density fit's error*: the BALL (which G are kept) and the RASTER (how finely r is sampled to get those
  G).  AliasFree fixes the raster: every G inside the ball comes out exact.  It does not size the ball.  On a PW orbital basis the
  ball can be exact (pin 28: Ecut_fit ≥ 4 Ecut, Parseval).  On Gaussians it cannot: a pair φ_iφ_j with exponent sum p has the
  spectrum e^{−G²/4p}, never exactly zero.  So "perfect" = **the ball is sized from a tolerance ε and the tightest pair**,
  e^{−G_max²/4p_max} ≤ ε (per multigrid level), not from `cutoffFactor · α_max` alone.  That is also the "AUTO density cutoff misses
  O q6" row's cure, with the PP's own sharpness as one more "pair".
- *For XC the ball needed is the SAME ball.*  The Fock matrix needs ⟨i|v_xc|j⟩ = Σ_G (φ_iφ_j)*(G) v_xc(G), and the pair spectrum is
  ≤ ε outside the density ball, so v_xc's components outside it do not reach H_ij.  A v_xc fit on the density ball is therefore
  perfect for the Fock matrix (to ε) — which is why §K always said "Vxc and CD in lockstep".
- *The one real knob left is the XC RASTER.*  v_xc(ρ(r)) is NOT band-limited (ρ^{1/3} where ρ is small), so sampling it on a finite
  raster folds its high-G content into the ball, and E_xc = ∫ε_xc(ρ) is a quadrature with the same error.  That is a convergence
  parameter, not a fit — QE has the same one (its default XC grid was 107 μHa off on β-MnO₂).  Make it self-reporting like the grid
  audit: evaluate on N and 2N rasters once, print the difference, escalate (P3) when a consumer needs better.  The Becke route stays
  the alternative quadrature.
- So C2's §K step = one tolerance-sized ball for ρ AND v_xc + an XC-raster audit, in the same re-bank window as the raster default.
  Acceptance is grid-convergence of ρ and of a property, never ΔE (pin 3).  (Seen in passing: `GPW_IBS`'s constructor defaults
  `raster = BallOnly` while `GPW_Evaluator`'s member initializer says `AliasFree` — C2 makes one default.)

**Refinement after reading `Records/GridSharpTerms.md` (2026-10-29) — the ε is on the TERM'S ERROR, not on the raw spectral tail.**
- Sizing the ball from the raw pair tail e^{−G²/4p} ≤ 1e-10 would cost Ecut ≈ 46 p — a ~250× larger raster.  What a term needs is
  its MATRIX ELEMENTS right to ε, and each kernel converts that differently: Hartree's 4π/G² suppresses the tail (its ball is
  forgiving); V_loc long uses the harmonic rule p β/(p+β) — **already ε-sized today** (`vlocEps` 1e-5); the local-PP sweeps already
  use a κ rule (`localPPRelCutoff` 30 Ha = e^{−15}, CP2K's REL_CUTOFF).  So the ε method already exists for two terms.
- The measured failure (GridSharpTerms §3) was the RASTER, not the ball: alias e^{−G_N²/p}.  `BallOnly` at the auto floor C = 2 gives
  ±1.5–2.7 Ha on a tight pair's V_long diagonal; Ecut/p ≥ 5.5 (alias ≈ e^{−11}) gives ≤ 1e-5 Ha.  So **the raster N is what ε must
  size**: G_N² ≥ p · ln(1/ε) for the tightest pair on each multigrid level (ε = 1e-6: G_N² ≥ 13.8 p; `BallOnly` gives 2p at C = 2,
  `AliasFree` 8p).  Per level, a pair-to-level routing rule in the REL_CUTOFF style, derived from ε.
- **The names go.**  With N derived from ε, `RasterPolicy {BallOnly, AliasFree}` has no job: the raster is whatever meets ε, and
  the audit becomes an assert that it did.  The explicit CP2K-parity overrides (`mgridEcuts`, `relCutoff`) stay as stated choices.
- **Is it the lattice-sum ε?  Same PRINCIPLE, not the same NUMBER.**  `screenEps` (1e-10) drops r-space terms whose Gaussian
  prefactor is below it; the reach grows like √ln(1/ε), so a tight value is nearly free.  A grid's cost grows like (ln 1/ε)^{3/2} in
  points, so the grid ε must be the loosest that meets the accuracy target.  The code's own comment already says they are "not converged
  by one tolerance".  **RULED 2026-10-29 (user): ONE user-facing accuracy target** (Ha per matrix element), from which every mechanism
  derives its internal tolerance by a stated formula (screenEps, densityEps, vlocEps, κ, the raster N), printed together on the
  banner.  Two reasons: (1) simple for the user; (2) **CONSISTENT error across terms** -- "there is no point in having one term
  accurate to 1e-10 and another only good to 1e-4".  So each formula is judged by the error it DELIVERS on its term, measured, and a
  term far tighter than the target is wasted cost just as one looser is a defect.  Candidate for a pin once built.
- **V_loc long analytic?  Yes, it can be, and there are two routes.**  −Z erf(√β r)/r is the potential of a Gaussian ion charge, so
  ⟨χ_i|V|χ_j⟩ is a smeared-nucleus attraction integral (a Boys function; libcint has the erf-attenuated 1/r).  The lattice sum is
  conditionally convergent, so it needs an Ewald split: erf(√β r)/r − erf(√η r)/r is short-ranged (an analytic real-space lattice sum,
  like the GTH short part already is), and erf(√η r)/r with a small η is so smooth that its G-sum Σ_G ρ_ij(G) V(G) needs only a tiny
  ball, with analytic pair transforms — no raster at all.  The other route is CP2K's: add the ions' Gaussian charges to ρ̃ and solve
  ONE Hartree problem for the neutral total, so V_loc long is no longer a separate raster-integrated term, though it still uses Hartree's
  raster.  Either removes the term with the biggest measured alias error, and an analytic route gives exact force derivatives
  (no egg-box from that term).  Hartree and XC stay on the raster either way.  C2 decides by measuring cost, first a design note.

## 8. Decisions waiting on you

1. ~~Does lane H count toward your two sessions?~~  Settled: no session cap; the box is the limit (§5).
2. ~~AliasFree default~~  Agreed, with the measure-then-decide caveat (§7).
3. ~~R1.0r~~  Ruled: the intersection (§7).
4. ~~GGA before forces~~  Agreed.
5. ~~§K~~  Ruled: one tolerance-sized ball for ρ and v_xc + an XC-raster audit, in C2 (§7).
6. ~~P1–P5~~  Accepted ("very sensible"), plus P4a, the regime matrix (§3).
7. ~~V-CALCNET~~  Review it first in lane B, covering Molecule and Solid (§1, §5).

## 9. Tag index (C0, 2026-10-29) — every open row, by tag and chunk

Row names abbreviated; the row itself lives in its queue.  The controller updates this index at each merge, and the open F+V+S count
is read from it.

**F (19)** — C1: "O₂ triplet does not converge on the molecular facade" (lane G S4 2026-10-30; a panel case), "Linear D-mixing diverges" (remainder), "Na₂ AFM seed", "GDM crawls",
"Kerker+Pulay overshoots β-MnO₂'s U = 1 basin" (a panel case), "Polarized Γ Si jumps", "GDM: should an exhausted line search END the
stage?", §2 "`MixingPolicy`" (becomes P4a's regime matrix); "the Kerker field mixer is not exactly group-equivariant", "the low-rank D factor can split a symmetry orbit" (B2 R2.1, 2026-10-30: a density that leaves its group is not variational, so its forces are not dE/dR).  (Closed by lane A, 218956da: "A threaded NaF exact resume STALLS GDM" -- the re-seat now seats the DM source, `7a016d0b`.)  Lane B: OpenWork "PES / run-setup remainders" (landed ce639536; gate (c) = E continuous along a geometry path, the NaF Γ path test, the Löwdin C carry -- the relaxation's own prerequisites; it now holds all of OOD V-CALCNET's open work, so V-CALCNET is counted here, not under S).  "NaF's F site is a MAXIMUM at the path test's level" (lane B 2026-10-30; suspect: the q1 Na PP lacks the semicore's Pauli repulsion; CP2K q1/q9 k222 arms next).  C2: "BallOnly raster aliases tight pairs", "AUTO density cutoff misses O
q6", "uniform-raster XC not equivariant", "a FREE run's XC mesh" (cost remainder), OOD R1.0r, CleanCode §K.  Forces note: §2 "SCF
checkpoint/restart" (CK-3 = restart across a geometry step).  (Closed on main: "`u_bands` after a GDM stage" 07af0909, "GDM false-converges" ab1c206a.)

**V (5)** — C4: "restart fingerprint does not record the run POLICY" (the remainder), "spinel charge order chosen by
atom index".  (Closed by lane B, 827d0e47: "valgen --seed-into", D-BASIS-REGISTRY, D-STATE-SCHEMA.)  (Closed by lane B, ce639536: "PW states cannot be saved", D-BASIS-PP, D-STRUCTDATA step 2.)  C5: "valgen is
optimistic about diffuse functions", "vet trim deletes the valence core", "Diffuse-basis ACTUATOR".  Feature on the path: §2 "PBE / GGA".

**S (26)** — C1: CleanCode D-CALC-DOUBLE-SCF (migration) (V1.24(i), D-DIIS-FDMIN closed by B3, 2a3852c9); OOD V1.27, D9; OpenWork "`ρ_lost/N` dead
column" (rides V1.27).  C2: OOD V-MIXER-GBALL (= R1.5, moved from C1 2026-10-30), R1.0 (metric axis), V1.1b, BM(1), D2, D3, D4, D13, D-GPWPAIR, R1.0n/o.  C3: OOD V1.41, V1.35, R1.0q,
V1.34.  With GGA: OOD §I.1, V1.39.  C4: CleanCode D-ENV, D-ENV-TOL;
Lane B design review: OOD V-CALCNET -- now counted under F as "PES / run-setup remainders".  OpenWork "The SCF observer has two homes" (controller 2026-10-30, lane B).  C1: OpenWork "R2.1 as a throw needs a measure that holds on a REDUCED evaluation".  C6: CleanCode D-WARNINGS.  Forces prerequisite (lane A C1/C3, 2026-10-29): OpenWork "short tier not green at `QCHEM_OPENMP_THREADS=4`" (thread-count invariance: forces and a relaxation need energies independent of the thread count).

**H (lane H)** — §1; §2 DFT+U, A7 estimators, D5 little group, the PW umbrella + PW-1/2/4, KP; §3 NiO ghost, MnO ordering, the campaign
re-run, NiO/SrVO₃ vs hp.x, the PW 4× XC ball re-run, the three β-MnO₂ rows; §4 LR ledger, dense PW cost, LR on real TRIM; OOD V1.40;
CleanCode D-ENV-RERUN, D-CKALPHA-AUTO.  OpenWork "CP2K oracle: is `ADDED_MOS` converged on the smeared TMO decks?" (lane G 2026-10-30).

**X (21, the quota)** — OOD V-XC-ONE-TERM (lane B, after the δ-on-the-real-path chunk; user 2026-10-31); OpenWork "The SCF iteration row has two homes" (lane B after `deck::Build`; the GUI's live plots wait on it); OpenWork "molecular facade default XC mesh 0.5 Ha/atom off" (lane B, FIRST in its quota: a wrong default every GUI molecule hits), "molecular PP route 4 mHa/atom off CP2K on F₂" (suspect-operator gate); OpenWork "one XC mesh for standalone runs and sessions?" (lane B 2026-10-30); OpenWork "`A_HF_dfPin.SelfPin` parameter names are raw bytes" (lane G 2026-10-30); OpenWork "`CalcOptions`'s default basis `sto-3g` is not registered" (lane B 2026-10-30, user's call: register STO-3G or default to `dzvp`); OpenWork "R2.1 cannot measure the raster XC route" (rides C2); OpenWork "every `_Long` test re-proves its claim at the cheapest resolution" (lane A, right after C1 B2 lands; it shares the vendored-state mechanism with panel #19); CleanCode D10, R2.24, R2.23, R2.21, R1.0f, D-POINTPROBE, D-TIMING-GATE, D-XCROUTES-GROWTH, D-RK-ONCE,
D-CONTRACTED-SHELL, D-FIT-MERGE; OpenWork "Blaze heevx/syevx" (rides ModuleToolchain Step 4).

**P (parked / features / triggers)** — §2 B-spline basis, hybrids, GW, transport/MD, USPP/PAW, OT (and with it CleanCode V1.24(iii)),
space-group irreps, SSB descent, smearing finished, molecular SAD seed, spherical SALC, valgen --auto, run report, lattice gas; §3 the
136-function span, ‖V_xc − V_xc_fit‖, charged cells, the anchor-batch pointer; CleanCode R1.0b, V1.29, D-MODTOOL (Steps 2–5),
D-CLOUD (yours), D-CP2K-LIBXC; OOD R1.0 step 2, V1.38, V4.1/V4.2, R2.18, V1.19.

**Count at C0: F 16 · V 11 · S 27 = 54 on the critical path.  Then (18ab4db5): F 16 · V 11 · S 28 = 55; b8a36fd3 +1 F = 56.  Then (ce639536): F 18 · V 8 · S 27 = 53; +1 S (observer's two homes) = 54.  Then (2a3852c9, lane A B3): F 18 · V 8 · S 26 = 52.  Then (827d0e47, lane B): F 19 · V 5 · S 26 = 50; +1 F (O₂ triplet) = 51.  Now (218956da, lane A): F 19 · V 5 · S 26 = 50** -- B2 closed one S and its R2.1 measurement found three critical-path defects (2 F, 1 S): the instrument working, not slippage.  Lane B steps 6–8 closed three V; V-CALCNET's remainder re-tagged S → F (the relaxation's path gate).

## 10. Landing log (controller)

| main | from | sweep | anchors moved | baseline kept |
|---|---|---|---|---|
| `a3f4b50f` (2026-10-29) | lane-b C4 step 1 (`a2439643`: `PseudopotentialSet`, `pseudopotentials.json`) | serial short 1121 + long 7/7 | 0 / 1114 short, 0 / 8 long vs `capguard_*` | `~/Code/qchem6-runs/anchors/main-a3f4b50f_{short,long}.log` |
| `e9955ebf` (2026-10-29) | lane-b C4 step 2 (inject the `PseudopotentialSet`; clean replace of the molecular PP fields) | serial short 1123 + long 8 | 0 / 1119 short, 0 / 8 long vs `main-a3f4b50f` | `~/Code/qchem6-runs/anchors/main-e9955ebf_{short,long}.log` |
| `87493dd9` (2026-10-29) | lane-a C1 B0+B1 (the convergence panel, 19 cases, RED as WILL_FAIL; u_bands + GDM false-converges closed; thread-invariance row) | serial short 1126; long 13/13 at 28fe4a31 (4 threads, panel incl.) | 0 / 1112 short vs `main-e9955ebf` | `~/Code/qchem6-runs/anchors/main-87493dd9_{short,long}.log` (long = 4-thread, re-baseline serial at the next long sweep) |
| `20873b21` (2026-10-30) | lane-b C4 steps 3–5 (`84d5569f` seed library keyed by PP variant q; `47f0816a` basis provenance PP-PROVENANCE + CheckPPProvenance) | serial short 1130 (merged) + serial long 7/8 pre-merge | 0 / 1124 short vs `main-87493dd9`, 0 / 8 long | `~/Code/qchem6-runs/anchors/main-20873b21_{short,long}.log` (long = SERIAL baseline) |
| `18ab4db5` (2026-10-30) | lane-a C1 B2 (`f6bb5996` + `77a99402`: R2.1 group residual REPORTED, R2.4, R2.5 siteSpins required, R2.7, panel #10 RED; panel out of the routine long tier) -- a merge commit, main's delta since the lane's base was docs only | serial short 1133; serial long `-L long -LE panel` 6/6; panel 20/20 at 4 threads | 0 / 1128 short, 0 / 7 long vs `main-20873b21` | `~/Code/qchem6-runs/anchors/main-18ab4db5_{short,long,panel4t}.log` |
| `ce639536` (2026-10-30) | lane-b C4 steps 6–8 (`SolidSetup`/`MoleculeSetup::Resolve`, `GaussianBasisRecipe`, `PotentialEnergySurface::ConvergeAt(cell, R)` + both sessions), PW state save, clean replace (old ctors deleted, 148 callers migrated) -- a merge commit, main's delta docs only | serial short 1147/1147; serial long `-L long -LE panel` 6/6; long panel 7/7 at 4 threads inside budget (MnO₂ K-P 4.4e-7 Ha / TiO₂ 76→84 it vs the threaded baseline: thread variance, the invariance row) | 0 / 1131 short, 0 / 7 long vs `main-18ab4db5` | `~/Code/qchem6-runs/anchors/main-ce639536_{short,long,panel4t}.log` |
| `2a3852c9` (2026-10-30) | lane-a C1 B3 (R1.1 joint Pulay hand-off, R1.3, R1.4 `EngageBelowFD`, FDMin deleted + the positional-DIIS-init fix `7d6222b0`; vendored-state twins; deterministic threaded collocation `5ec26b06`; GridEpsilonPlan doc) -- a merge commit, main's delta docs only | serial short 1149/1149; serial long `-L long -LE "panel|periodic"` 5/5; periodic 3/3; long panel 8/8 at 4 threads inside budget (MnO₂ ×3, TiO₂ moves: the now-deterministic threaded partition vs main's random draw) | 0 / 1144 short, 0 / 4 long, 0 periodic vs `main-ce639536` | `~/Code/qchem6-runs/anchors/main-2a3852c9_{short,long,panel4t,periodic}.log` |
| `b9650a5c` (2026-10-30) | lane-g step 1 (`4fb54e66`: the fresh nanobind bridge on `MoleculeSetup::Resolve` + `Calculation`; facade: constructor observer, `Levels`, `GetSpinDensity`) -- a fast-forward | serial short 1152/1152 (no SCF path changed: no long tier or panel) | 0 / 1146 short vs `main-2a3852c9` | `~/Code/qchem6-runs/anchors/main-b9650a5c_short.log` (long/panel/periodic: symlinks to `main-2a3852c9`) |
| `827d0e47` (2026-10-30) | lane-b (`0ec790b1` NaF Γ path test + a session freezes the plain atom-centred XC mesh, `xcMeshSymmetric`; Q3 UT; fingerprint gaps 2–4 / D-STATE-SCHEMA 1.1; D-BASIS-REGISTRY; valgen --seed-into; Löwdin carry, AO stays default) -- a merge commit, main's delta docs only | serial short 1160/1160; serial long `-L long -LE panel` 8/8 (no SCF path changed: no panel) | 0 / 1146 short, 0 / 5 long vs `main-2a3852c9` | `~/Code/qchem6-runs/anchors/main-827d0e47_long.log` (short: lane B's log) |
| `218956da` (2026-10-31) | lane-a (`7a016d0b` the fallback re-seat seats the DM source -- the threaded NaF resume; the rounding-level line-search accept REVERTED `8650f2d2` after the sweep caught GDM limit cycles; C2 M1 `3c3f5ee1`: AliasFree 2.2–6.2× wall, BallOnly stays the default) -- a fast-forward, main merged in at `37136132` | serial short 1160/1160; serial long 5/5 on the merge; `7a016d0b`: panel 4t + periodic 0 moved apart from the known MnO₂ threading moves | 0 / 1159 short, 0 / 6 long vs `main-827d0e47` | `~/Code/qchem6-runs/anchors/main-218956da_{short,long,panel4t,periodic}.log` (panel/periodic from `7a016d0b`) |
