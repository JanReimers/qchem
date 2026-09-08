# Open work — history 3 (cut 2026-09-08)

The third closed-record split out of `doc/OpenWork.md`, after `OpenWork_History1.md` (threads A–E, the
2026-06-30 and 2026-08-15 orderings, runtime rounds 1–4) and `OpenWork_History2.md` (the Vxc repair thread
and the fit-basis interface refactor).

**Why this cut** (user, 2026-09-08): *"OpenWork.md and doc/CleanupCandidates.md have done and todo items
interleaved and lots of history so it is hard for me to read and assess."*  The tracker had grown to 2859
lines of which these eight sections — every one of them CLOSED, ACTED ON, or explicitly marked *no action
here* — were 1316.  Each leaves a stub in the tracker naming what it concluded and pointing here.

NOTHING IS TRIMMED.  The same reasoning as `doc/CleanupHistory.md`'s preamble applies and has been earned
here too: the refutations are worth more than the landings.  Three of these sections exist ONLY because a
measurement killed the obvious answer — the exp recurrence, lever B, and the diagonal-seed fix — and each
would otherwise be re-proposed.


---

## The N1 outcome-detector tiers T1–T5 — harvested 2026-09-08 (was OpenWork.md L593-777)

*(verbatim)*

## ✅ CLOSED — the N1 tiers, T1 through T5 (kept for the evidence, not as a next action)

**T1 and T2 are BUILT (2026-08-25).**  What follows is the remainder, in order, with the state each one
starts from.  Everything they need now exists: `Outcome<T,E>`, `SCFFailure::Why`, the Hamiltonian's
`SiteMoments` aggregate, and the MnO Gamma test running through `SolidCalculation`.

### ✅ T1 — DONE.  A non-converged energy is UNRETURNABLE
`Converge()`/`Result()` return `Outcome<Converged, SCFFailure>`; `Energy`/`EnergyTerms`/`TotalCharge`/
`Density`/`SpinDensity`/`DensityMatrix` live ONLY on `Converged`.  `LastIterateTerms()`/
`LastIterateCharge()` keep deliberately-unconverged numbers reachable under a name that cannot lie.

### ✅ T2 — DONE, AND ITS POSITIVE PATH IS NOW EXERCISED (2026-08-26)
An imposed run that loses its magnetic order fails as `Why::OrderLost`, judged on the INTEGRATED Becke
site moment.  Self-calibrating (imposed ∧ the run carried order > 0.05 e ∧ collapse below 1% of its peak),
so a FREE run finding m=0 — physics — is never touched.

**★ THE POSITIVE PATH: `GPW_SCF.ImposedOrderLostIsAPostconditionFailure_Na2Box`.**  Na2 in a box at its
bond length: two neutral Na (so the SAD seed plants the library's spin pair, ±1 e), AFM flip on the
second, imposed, Becke mesh PINNED — and a ground state, the closed-shell σ², that has no order to keep.
It converges in 66 iterations and fails as `OrderLost` with the moment dead from iteration 9.  Its mirror
is `ImposedShubnikovHoldsAFMThroughSCF_Mn2Box` (same shape, real magnet, order SURVIVES), so the pair says
the detector *discriminates* rather than always firing.

**⛔ AND BUILDING IT FOUND A DEFECT IN T2's OWN BASELINE, now fixed.**  T2 took its yardstick from
`GetWaveFunction()->GetChargeDensity()` at construction — which is NOT the seed.  `Initialize` builds a
Fock from the seed, diagonalizes and FILLS, so the earliest reachable density is already one aufbau fill
downstream, and that fill is exactly where a fragile order dies.  Measured: a Na2 seed staggered at
**±1 e** reads **±0.07 e** one fill later — under any honest floor, so the postcondition SKIPPED the very
case it exists for.  ⇒ `SolidCalculation` now BUILDS THE SEED ITSELF (the same `MakeSeedDensity` call the
iterator would have made) and hands it to the explicit-seed ctor, measuring the raw seed's moments first:
**0.4686 e**, and the check has teeth.  Two ends, two ways to be misled — the yardstick is now
`max(raw seed, any iterate)`, because a moment that GROWS before it dies (MnO: 0.0046 → 0.106 → 7e-5)
needs the high-water mark just as much.

### ✅ T3 — DONE.  THE Eee CHARGE-SLOSH DETECTOR, AND ITS THRESHOLD IS MEASURED
`RunDiagnostics::ChargeSloshed()` on `SolidCalculation`, feeding `Why::ChargeSlosh`.  Eee is now carried on
`SCFProgress` (the whole `EnergyBreakdown`, not the one term today's detector reads), the facade records it
every iteration, and the rule is:

**Eee at the END, over the run's OWN floor, and only on a run that DID NOT CONVERGE.**  Three decisions,
each of which a measurement forced:

1. **END, not peak.**  Measured on the healthy MnO baseline 2026-08-26: Eee runs `14.17, 15.23, 12.51`
   and then settles at `13.18`.  A peak-based ratio reads **1.22** on a run that did nothing wrong —
   every SCF passes through a transient leaving the seed, and one that overshoots and comes back is
   a healthy run doing its job.
2. **The run's own floor, never an absolute number.**  Eee is extensive; 13.5 Ha is healthy for this MnO
   cell and would be a catastrophe for Si.  (This is the shape T2 uses, and the row below is why it had
   to be self-calibrating.)
3. **⛔ NON-CONVERGED RUNS ONLY — and the predicate says so ITSELF, rather than trusting the caller to
   remember.**  This is the sharpest thing the calibration turned up.  A converged density is stationary,
   so "it sloshed" cannot be true of it; meanwhile a run that legitimately RESTRUCTURES on its way to the
   answer can move Eee a long way — **Na2 measures 1.71×, converged and correct**.  Without the gate the
   detector fails a good run, trading a missed diagnosis for a WRONG one.

| run (2026-08-26 unless noted) | Eee floor → end | end/floor | converged | verdict |
|---|---|---|---|---|
| MnO AFM-II imposed, baseline | 12.51 → 13.18 | **1.05** | yes | healthy |
| MnO AFM-II imposed, `GPW_XC_DM_SOURCE=1`, MOM **on** | 10.99 → 13.15 | **1.20** | yes | healthy — see the recipe note below |
| MnO VA-spherical imposed, `GPW_XC_DM_SOURCE=1`, MOM **on** | 11.00 → 13.18 | **1.20** | yes | healthy — same |
| Na2-in-a-box, AFM seed → singlet | 0.114 → 0.195 | **1.71** | yes | healthy (a real restructuring) |
| **MnO AFM-II imposed, `MNO_KERKER_G0=0.01`** (flat filter) | **14.17 → 35.10** | **2.48** | **no** | **SLOSH — detector FIRES** |
| **MnO, the BANKED recipe + `GPW_XC_DM_SOURCE=1`** | **14.42 → 29.00** (peak 65.01) | **2.01** | **no** | **SLOSH — detector FIRES** |

⇒ **`kSloshFactor = 1.5`**, above every healthy run measured and below every collapse; and because the
detector is only consulted on a run that ALREADY failed, being wrong costs a LABEL, never a good answer.

**★ THE POSITIVE PATH IS EXERCISED TOO** (the same "an unexercised guard is a guess" rule T2 got held to).
The flat-Kerker arm re-run 2026-08-26 reproduces its 2026-08-25 collapse exactly and the facade reports it
as `Why::ChargeSlosh`, not as `NotConverged`:

```
[MnO AFM-II Gamma] NO ANSWER: the Hartree term ran away and the SCF never converged: Eee ended at
2.476827x the lowest value this run reached, which is the signature of an UNDAMPED low-G charge mode
[order(integrated site moment): seed 4.781 e, peak 4.781 e, final 4.3e-15 e -- DIED at step 7;
 Eee 14.17 -> 35.1 Ha (peak 37.99, end/floor 2.477)]
```

E = −38.510, 80 iterations, `OSCILLATING`.  And the Eee series is worth reading — after iteration 9 it is a
clean PERIOD-2 LIMIT CYCLE, `35.210, 35.098, 35.210, 35.098, …` for seventy iterations.  That is charge
sloshing in the most literal sense available, and it is the reason the term is the right instrument.
⚠ Note both detectors are live on this run and the ranking picks the mechanism: the order died at step 7
(4.781 e → 4e-15 e, the integrated moment) and rides in `details`, while `why` names the charge runaway
that caused it.  It is a demonstration of the ordering rule, not an accident.

✅ **AND THE ONE THING THAT FIRST LOOKED LIKE A NON-REPRODUCTION IS RESOLVED — IT WAS THE RECIPE.**
Re-running `GPW_XC_DM_SOURCE=1` on 2026-08-26 it first CONVERGED to −61.41457 with the order intact, which
did not match the −45.53 of 2026-08-25.  Two hypotheses were tested and both REFUTED before the real cause
turned up, and both are worth keeping because they are the obvious guesses:

- **NOT the diffuse end of the basis.**  `SR` and `VA` have the SAME \f$\alpha_{\min}=0.1\f$ on this cell.
  (`SR` is diffuse-trimmed relative to `valence_lowq`, α_min 0.03 → 0.1; `VA` descends from `SPH` which
  descends from `SR`, so it inherits that trim.)  What they differ in is CONDITIONING, not reach:
  `SR` = 122 hand-trimmed at full rank, min kept pivot **0.0236681**; `VA` under Cartesian d = 132 that the
  pivot filter auto-drops to 122, min kept pivot **0.0236681** (the same set, two routes — see the vet-stage
  trim item under *Continuous — CLEANUP*); `VA` spherical = 118 hand-trimmed at full rank, min kept pivot
  **0.00381699**, i.e. ~6× worse conditioned than SR while being full rank and equally diffuse.
- **NOT the basis at all.**  The armed flag on the banked `VA`+spherical basis STILL converges
  (−61.41127, 48 iterations, order survives, end/floor 1.20).

**IT IS THE OCCUPATION MACHINERY.**  The banked command is
`MNO_ANNEAL="5e-3,0" MNO_ACC="Ladder,GDM" MNO_MOM=0 MNO_ORTHO_TOL=1e-3 MNO_SHARED_MU=1` — **MOM OFF
(aufbau) and a SHARED Fermi level**, so the moment is free to relax and nothing holds the occupation.  This
test's own STATUS block already says it: *the order and ~9 Ha of binding are BOTH bought by holding the
occupation*.  Run verbatim with the flag, 2026-08-26, it reproduces the banked collapse **to all ten
significant figures**: **−45.52875429**, Eee 29.00 against the banked 28.995, the integrated site moment
4.781 e → 1.1e-14 e dead from step 4.
⇒ Item 1's row STANDS, and now with its scope stated: the flag collapses **the fragile recipe** (aufbau +
shared μ), and leaves the MOM-held recipe alone.  ⚠ Do not quote it as "the flag always collapses", and do
not quote today's converging rows as a retraction — they are a different recipe.

### ✅ T4 — DONE.  THE DETECTORS ARE IN THE LIBRARY
The `** DIED at iteration N` post-mortem is now `RunDiagnostics` on `SolidCalculation`, reachable through
`Diagnostics()`, so every caller gets it — the MnO test only PRINTS it now.  Two hooks, composed by the
facade so the caller's telemetry still works: the **order probe** is the only place this iteration's
density is in hand (where the integrated site moment is free — the XC term has already rastered that
density serial), the **observer** is the only thing that fires exactly once per iteration.  So the probe
measures and the observer files.  A caller that sets no probe of its own gets the integrated moment in
the trace column under `m_site` — an INTEGRATED observable, not a point sample of m(r).

### ✅ T5 / N5 — DONE.  SELF-DESCRIPTION AND `CP2K_COMPAT`, BUILT TOGETHER
`qchem.RunPolicy` (a qcCommon leaf) resolves the deviation set ONCE, and the four env flags that were read
at their point of use now come from it: `QCHEM_DM_LOWRANK` (ChargeDensity factory), `GPW_STREAM_FOLD`
(GPW basis factory), `QCHEM_MIX_RHO_M` (`MakePeriodicMixer`), `GPW_XC_DM_SOURCE` (the XC term).
`CP2K_COMPAT=1` turns them all off; **an explicitly-set individual knob still WINS** (saying
`CP2K_COMPAT=1 GPW_STREAM_FOLD=1` is a deliberate act, and the banner marks it `(stated)` so a run that
thinks it has parity and does not says so out loud).  `ReresolveRunPolicy()` is the named escape hatch for
the two acceptance gates that A/B the fold in one process — named, rather than making "resolved once" a
comment instead of a property.

`SolidCalculation` prints the banner UNCONDITIONALLY, four lines at construction plus one per `Converge`:

```
[MnO AFM-II Gamma run] system: 4 atoms, 26 valence e, multiplicity 1 (POLARIZED), seed=IonicSAD
[MnO AFM-II Gamma run] grids: densityEcut=auto C=2 raster=BallOnly xcMesh=Becke (nR=40 L=29)
[MnO AFM-II Gamma run] symmetry: IMPOSED (Shubnikov from the decoration);  threads: OMP_NUM_THREADS=1 GPW_OMP_THREADS=1 (BLAS pinned to 1)
[MnO AFM-II Gamma run] CP2K_COMPAT=0 -> DEVIATING;  QCHEM_DM_LOWRANK=on*  GPW_STREAM_FOLD=on*  QCHEM_MIX_RHO_M=off  GPW_XC_DM_SOURCE=off   [* = differs from CP2K]
[MnO AFM-II Gamma scf] mixer: Kerker(G0=1.000000) alpha=0.45;  XC rho source: rho_mix;  accel: Ladder;  kT=0.005 MOM=on NMaxIter=80
```

That closes the measured defect this item names — with `MNO_KERKER_G0=0` the fall back to linear D-mixing
was ENTIRELY silent (the mixer identity appeared only in the Verbose per-iteration column) — and it is
what makes a `doc/Benchmark.md` row self-describing instead of relying on discipline.
⚠ **STILL OPEN under N5 — and the table is SHORTER THAN THE TRUTH.**  Four deviations are wired; at least
four more are known and not, so a `CP2K_COMPAT=1` row today is closer to parity than the default but is
NOT parity.  Measured/checked 2026-08-26:

| missing deviation | why it is not in the table | measured effect |
|---|---|---|
| ~~**the pair-stream CACHE**~~ (gap 2, user) | ✅ **DELETED 2026-08-27** (`doc/CollocationRewritePlan.md` step 7) — there is no deviation left to declare: qchem re-evaluates every iteration exactly as CP2K does, off a ~0.2–0.4 MB task list | it WAS the single biggest RAM term: MnO peak RSS 3915 → 155 MB on the free probe, 1323 → 463 MB on the imposed benchmark row |
| ~~**`imposeSymmetry` ITSELF**~~ | ✅ **WIRED 2026-08-26** (user: *"CP2K_COMPAT should do (imply) imposeSymmetry=0"*) — the fifth declared deviation, knob `QCHEM_IMPOSE_SYMMETRY` | CP2K does **NO** symmetry work in these decks (see below); our imposed row folds the BZ, star-averages ρ, uses the site-adapted invariant XC mesh (~2×) and folds the streams (5.2× on MnO pairs) |
| `raster` (`BallOnly`) | typed option | BallOnly IS CP2K's bet (N2) — a deviation in mechanism only |
| `cutoffFactor` (C=2) | typed option | ~0.15 mHa of grid error at C=2 (N2) |

★ **CP2K DOES NO SYMMETRY BLOCKING AT ALL — verified locally, not assumed.**  The 1129-line
`IntegrationTests/CP2K/bench_MnO_AFM2_VA_cp2k.log` contains **zero** occurrences of "irrep", "symmetry" or
"point group": QuickStep keeps K and P as DBCSR sparse **atom-block** matrices over the full AO basis and
diagonalizes the whole thing — its blocking is atom-pair SPARSITY, not irrep.  No SALC blocking, no
k-block splitting by irrep.  The one symmetry knob that exists is BZ-side and is OFF in our own deck:
`BRILLOUIN| K-Point point group symmetrization  OFF`, with all 8 k-points of the 2×2×2 mesh listed.
⇒ **doc/Benchmark.md's MnO rows compared qchem-WITH-symmetry against CP2K-WITHOUT**, which flattered
qchem on exactly the axis the table measures.

✅ **FIXED 2026-08-26: `CP2K_COMPAT=1` now IMPLIES `imposeSymmetry=0`.**  It is the fifth entry in the
deviation table (`QCHEM_IMPOSE_SYMMETRY`), and it is the ONE that OVERRULES THE CALLER rather than merely
supplying a default — because every banked recipe sets `MNO_IMPOSE=1`, so a switch that let the recipe win
would need the recipe edited too, and then it would not be one switch.  The facade ANDs the caller's flag
with the policy's permission ONCE (`const bool imposed = opts.imposeSymmetry && ...`) and nothing below
reads the raw option again.  A vetoed imposition is NEVER silent: the banner prints
`symmetry: FREE  [asked for, VETOED by CP2K_COMPAT]`, the mirror of the hazard that made `imposeSymmetry`
opt-in in the first place.  `QCHEM_IMPOSE_SYMMETRY=1` is the stated escape hatch.
★ **MEASURED SAME DAY, AND THE EXPECTATION WAS HALF WRONG — IN THE INFORMATIVE HALF.**  The user
predicted this "may break MnO AFM convergence".  It breaks the CONVERGENCE and NOT the ORDER:
- stage 1 caps at 80 iterations at −60.431, stage 2 at 80 more at **−57.620** — 3.8 Ha short of the
  −61.40297618 the imposed compat run reaches in 24 iterations;
- but **m_stag 0.66 / 0.59 and the integrated site moment 4.781 → 4.222 e: the AFM order SURVIVED.**
⇒ **The imposed star-average was buying CONVERGENCE, not the magnetic basin.**  That is the opposite of
the standing assumption (S3/S4 read as "the imposition is what holds the order") and it moves the open
question from SYMMETRY to the MIXER — where the CP2K deck already points: 44 steps at Broyden α=0.2 /
NBUFFER 8 / MAX_SCF 200, against our α=0.45 / PulayDepth 0 / 80.
★ **AND IT VALIDATED T3's END-VS-PEAK RULE ON A RUN IT WAS NOT CALIBRATED AGAINST.**  This run's Eee
PEAKED at **39.98** and came back to 15.64 (end/floor 1.233).  A peak-based detector would have read 3.15×
and convicted it of charge slosh; the end-based one correctly reports `NotConverged` — which is exactly
what is wrong with it — and stays silent on both other channels.


---

## The on-the-fly box walk — harvested 2026-09-08 (was OpenWork.md L778-956)

*(verbatim)*

## ★★★ THE ON-THE-FLY BOX WALK — **2.21× ON MnO, DONE 2026-08-26** (was "~100× off CP2K")

> **USER, 2026-08-26:** *"I think it makes sense to chase the 100x gap in on the fly basis function
> evaluation.  Even if we just achieve 2 or 3x, then we can resolve other issues (MnO mixing) more
> rapidly."*

**RESULT: 2.21× on the MnO acceptance probe, 2.63× on Si, and the anchors did not move.**  Four edits to
the shell-blocked box walk (`NR_Evaluator::ForShellPairBox` and the two lambdas that consume it,
`scatterShell` / `integrateShell`), each one found by PROFILE and not by guessing.  Commits `5c0aca8d`,
`2cbdfaca`, `e8339cf2`.

⚠ **"KERNEL" WAS AN OVERLOADED WORD AND IS NOW RETIRED** (user, 2026-08-26).  Everything below says **BOX
WALK** and means one thing: the code that, for one shell pair at one cross-cell offset, walks the grid
points of the product's compact box and evaluates \f$\chi_i(r)\chi_j(r-R)\f$ there.  Every number is the sum
of exactly two timing-ledger buckets — `scf: integrate-back (pair gather)` + `scf: collocate density (pair
scatter)` — which are EXCLUSIVE and DISJOINT, so no setup subtraction appears anywhere in this section.

### ★ STEP 0 IS DONE, AND IT WAS AN INSTRUMENT PROBLEM, NOT A MEASUREMENT PROBLEM

The "853 s/iteration, read it as order 10²" caveat existed because `RunMnO` drives `SolidCalculation`
directly and **`SolidCalculation` opens no report run** — so the MnO arm was the one campaign run with no
timing ledger, and its cost had to be reconstructed as *3-iteration CPU minus an estimated ~182 s setup*.
`e8339cf2` gives the arm the same `GpwReport` bracket every other GPW driver holds.  The ledger then sums
to the wall clock (1201.3 s of 1202.1 s), so nothing is estimated:

⇒ **the true before-figure is 573 s/iteration, not 853 — the estimate was 1.5× pessimistic.**

### THE MEASUREMENT — MnO AFM-II Γ, uncached, A/B back-to-back on the same box

`MNO_SKIP_FM=1 GPW_MNO_NMAX=2 GPW_REPORT=1 GPW_STREAM_BUDGET_PTS=0 GPW_STREAM_BUDGET_PTS_F32=0`, the two
binaries differing ONLY in the box-walk diff (`git apply -R` of the src patch, same test source both sides).

```
[MnO AFM-II Gamma run] system: 4 atoms, 26 valence e, multiplicity 1 (POLARIZED), seed=IonicSAD
[MnO AFM-II Gamma run] grids: densityEcut=auto C=2 raster=BallOnly xcMesh=Becke (nR=40 L=29)
[MnO AFM-II Gamma run] symmetry: FREE;  threads: OMP_NUM_THREADS=unset GPW_OMP_THREADS=1 (BLAS pinned to 1)
[MnO AFM-II Gamma run] CP2K_COMPAT=0 -> DEVIATING;  QCHEM_DM_LOWRANK=on*  GPW_STREAM_FOLD=on*  ...
[fold] collocation streams (T3 pairs): NONE (7503 items evaluated in full)  [free/multi-k run]
```
★ **NO FOLD IS ACTIVE** on this row (free run) — the walk did the full unreduced pair set both sides, which
is what makes it a clean algorithm-to-algorithm A/B.  Measured 103% CPU, i.e. serial.

| bucket | before | after | ratio |
|---|---|---|---|
| `scf: collocate density (pair scatter)` | 753.78 s (41.88 s/call ×18) | **341.48 s** (18.97 ×18) | **2.21×** |
| `scf: integrate-back (pair gather)` | 386.06 s (38.61 s/call ×10) | **173.17 s** (17.32 ×10) | **2.23×** |
| **BOX WALK, total** | **1139.84 s** | **514.65 s** | **2.21×** |
| box walk PER ITERATION | 569.9 s | **257.3 s** | |
| `setup: collocation stream build` (same walk) | 30.83 s | 20.32 s | 1.52× |
| `setup: local-PP LONG` (same walk) | 15.46 s | 7.87 s | 1.96× |
| wall / CPU | 1202.1 / 1243.6 s | **558.9 / 596.3 s** | 2.15× / 2.09× |
| PEAK RSS | 165.8 MB | 166.0 MB | — |

**THE ANCHORS DID NOT MOVE.**  Identical trajectory both sides: `iters=2`, `lastΔρ=0.00854547`,
`m_stag 0.4121 → 0.3180`, `Eee 14.159 → 15.203`, integrated site moment `4.781 → 3.631 e`.  `Efinal`
differs by **2e-8 Ha** on a −59.7 Ha number (−59.69580385 → −59.69580383), entirely from the one edit that
is not bit-identical (below).

⇒ **Against CP2K's 8.5 s/iteration the standing on the uncached path goes from ~67× to ~31×.**  Note this
recalibrates the charter claim: it was never 100×.

### THE FOUR EDITS, in the order the profile produced them

| # | edit | Si box walk | bit-identical? |
|---|---|---|---|
| — | baseline | 4.818 s | — |
| 1+2 | **monomial power tables**, and `uintpow` split so it is not self-recursive | 3.803 s | ✅ |
| 4 | **incremental modulo wrap** | 2.541 s | ✅ |
| 3 | **`key/nn` decode hoisted out of the point loop** | 2.529 s | ✅ |
| 5 | **screen (3) = the reach SPHERE, not its rectangular hull** | **1.829 s** | sub-ε (see below) |

1. **THE HANDOFF'S CANDIDATES 1 AND 2 WERE ONE DEFECT.**  `pols[i](d)` per component was three
   `uintpow` calls, and BOTH `Polarization::operator()` and `uintpow` were out-of-line CROSS-DSO calls from
   the walk (12.1% + 1.7%@plt and 6.1% + 1.2%@plt ≈ 21%).  A shell's components are monomials over the SAME
   displacement, so the per-axis powers \f$x^0..x^{L_x}\f$ are shared: build them once per point and index.
   `uintpow` would not inline because clang will not inline a self-recursive function — the n>4 binary-
   powering tail moved to `uintpow_tail`, associations unchanged.  ★ **Grid-agnostic**, as the handoff said.
2. **THE MODULO WRAP** cost SIX integer divisions per point although the grid index advances by one per
   step.  `mx`/`my` and the row offset hoist to their own loops; the innermost keeps one compare.  A debug
   assert pins the incremental residue against the modulo.
3. **`fI[key/nn - si.begin]`** put a hardware integer divide on the innermost path — `divl` at **13.9%**,
   the single largest instruction in the profile.  The component-local slots are a property of the
   (pair, offset) pre-filter, not of the point.  Most of that 13.9% turned out to be the local-PP sweep,
   which is why the SCF buckets barely moved while `setup: local-PP LONG` fell 1.21 → 0.92 s.
4. ⚠ **THE ONE THAT IS NOT BIT-IDENTICAL — and the biggest single win (1.39×).**  Screen (3)'s bound was a
   fixed `lnE+12`, and for a DIFFUSE pair that exceeds the bounding box entirely (pMin=0.3 gives
   r_screen=10.8 au against reach=9.76), so **the screen never fired** and the CORNERS of the box — ~48% of
   its points, every one with an envelope already below `epsEff` — paid the full exp/poly evaluation, after
   which the consumer's own `|val|<eps` test threw the results away.  The bound is now
   \f$\min(\ln\varepsilon^{-1}+12,\ p_{\min}\mathrm{reach}^2+\mathrm{pf})\f$, i.e. the reach sphere the box
   is already the bounding box OF.  **The margin is not lost** — that expression is reach's own +1 a.u.
   polynomial margin restated in logs (+5.5 at pMin=0.3) — and the rectangular hull was an artifact of the
   walk order, not of any tolerance.  Evidence it sits beneath the noise: Si unchanged in every printed
   digit on both the cached and uncached path; NaF `Etot` identical to 12 s.f. against a cached-vs-uncached
   spread of 1.25e-6; MnO 2e-8 Ha; 771/771 twice.  `GPW_SPHERE_SCREEN=0` restores the hull for A/B.

### THE THREE-TIER HARNESS — CONFIRMED, with one correction to the handoff

- **Si** `GPW_SCF.SiliconGammaConverges` — the EDIT-MEASURE loop.  Box walk 4.818 → 1.829 s (**2.63×**);
  wall 7.13 → 4.00 s.  ⚠ Its XC mesh resolves to **Uniform**, so a Si profile is silent about Becke.
- **NaF** `DISABLED_NaFRocksaltGamma` — the PROFILING case.  ✅ It resolves to **BECKE**, so it does not
  share Si's blind spot.  Box walk 241.2 → 177.4 s for edit 4 alone; total run 271 → 206 s.
  ★ **AND IT SETTLES THE BECKE QUESTION**: the whole Becke side is Φ tables 1.03 s + ρ sampling 6.9 s +
  H_xc 0.8 s ≈ **8.8 s against the box walk's 241 s**.  The Φ tables are built once and cached, so that
  path is SETUP, not per-iteration.
  ★★ **AND WHY, from the user (2026-08-26): the Becke XC route uses a factored density MATRIX**
  (\f$D=LL^\dagger \Rightarrow \rho_g=\lVert L^\dagger\Phi_g\rVert^2\f$), so it needs **SINGLES, not pairs**
  — Φ per function plus a GEMM.  The two paths are disjoint by construction, which is why nothing here
  touches it and why the Si profile's Becke blind spot was harmless.
- **MnO** — the ACCEPTANCE run.  `MNO_SKIP_FM=1 GPW_MNO_NMAX=2` makes it a 9-minute probe rather than a
  45-minute one, and the ledger's per-call numbers make a bounded probe as good as a full run for cost.

### ⛔ THE exp RECURRENCE — TRIED, MEASURED, REJECTED (2026-08-26).  Branch `exp-recurrence-experiment`

> **USER:** *"We have to try it, and yes we need to be able to turn it on and off."*

Built behind `GPW_EXP_RECURRENCE` (=1 on, =2 audits against the direct \c exp, plus a debug assert on
every point), measured on all three tiers, and **parked on a branch rather than merged.**

**IT WORKS — the transcendental really is eliminated.**  \c exp falls from ~21% of the Si profile to
**1.39%**.  Numerically it is sound: max \f$|\Delta_{rel}|\f$ **4.8e-14**, max \f$|\Delta_{abs}|\f$
**1.5e-14**, and zero cases where the direct evaluation underflowed while the recurrence did not.

**AND IT IS STILL NOT WORTH IT.**  The box walk gains **1.13× (Si), 1.18× (NaF), 1.03× (MnO)** — and MnO
is the acceptance case.  The reason is structural and was not visible before measuring: the payoff scales
with LINE LENGTH, and the recurrence must be advanced on **every walked point** — including the ~48% the
sphere screen now rejects, where the direct form paid no \c exp at all.  MnO's lines are short, so a
transcendental on half the points was traded for four multiplies on all of them: a wash exactly where it
was needed.  ⇒ **Removing 20% of the profile bought 3%.**  Worth remembering as a general caution: a hot
symbol's share is an upper bound on the win, not an estimate of it.

⛔ **AND THE HAZARD IS WORTH MORE THAN THE SPEED — THE RECURRENCE IS ANISOTROPIC BY CONSTRUCTION.**  It
runs along z, so z is computed by a different arithmetic path than x and y.
`GPW_SCF.NaPseudoAtomInBoxDoublet` (one electron, 48 imposed ops, degenerate) converges to a DIFFERENT
BASIN, **0.97 mHa high**, in 9 iterations instead of 23 — while ISOTROPIC perturbations FOUR ORDERS
LARGER (`GPW_DENSITY_EPS` ×2 and ×0.5) all land in the correct basin to 9 digits.  Si and NaF, both
gapped, do not move at all.  **A 1e-14 axis asymmetry selecting an SCF basin is a different class of risk
from drift**, and it lands on exactly the degenerate open shells this campaign already fights (the
"Becke × degenerate open shell oscillates — fixed-axis angular grid vs rotating ρ" entry is the same
disease).
⚠ **AND MERELY HAVING THE CODE COSTS 5–8%**: the per-point branch put the default-OFF Si box walk at
1.92–2.00 s against **1.83 s** with the code absent.  A hot loop does not carry an off switch for free.

★ **A METHOD NOTE.**  The audit first reported a clean 4.8e-14 on the very run that was 1 mHa wrong,
because it returned early when the direct value underflowed to zero — silently excluding the single most
dangerous failure mode (direct 0, recurrence not).  The blind spot turned out to be empty here, but the
instrument was only trustworthy AFTER it was widened to the absolute deviation and made to count the
excluded cases.  **An audit that cannot see its own worst case is not evidence.**

**TO RESURRECT:** pair it with the interval skip below (which removes the wasted advances) and settle the
anisotropy first.  Neither is likely to change the MnO verdict.

### ✅ THE CHORD — bit-identical, and the acceptance case is the biggest winner (2026-08-27)

Screen (3) is a sphere about the product centre, so by the Gaussian product identity
\f$a_I|r-R_i|^2 + a_J|r-R_j|^2 \equiv p_{\min}|r-P|^2 + \mathrm{pf}\f$ the test is exactly
\f$|r-P|^2\le R_q\f$ — a CONVEX QUADRATIC in \f$t\f$ along a z-line, so the accepted set is an INTERVAL
with a closed form.  A line that misses the sphere is now skipped without visiting a single point of it.

| | box walk before | after | |
|---|---|---|---|
| Si | 1.829 s | 1.508 s | 1.21× |
| NaF | 177.4 s | 150.0 s | 1.18× |
| **MnO** | **514.6 s** | **366.9 s** | **1.40×** |

★ **MnO gains MOST, and the reason is the cell.**  A sphere spans more FRACTIONAL width in a skewed cell
than in a cubic one (the same sqrt(3)/sqrt(2) effect screen (2) was fixed for), so a rhombohedral cell's
bounding box wastes more of itself on corners.  The prediction was 1.15–1.25× from the Si-shaped estimate;
the acceptance case beat it.

**BIT-IDENTICAL BY CONSTRUCTION, not by measurement.**  \f$R_q\f$ is widened by a relative epsilon and the
solved interval by one point at each end, so the bracket is a strict SUPERSET of what the per-point test
accepts — that test stays the sole authority on what is kept.  And \c r is still ACCUMULATED through the
skipped head, never jumped to \f$r_y+t_0 s_z\f$: jumping would round differently AND make the walk
anisotropic, which is the defect that sank the exp recurrence.

**⇒ SESSION TOTAL on the MnO box walk: 1139.8 → 366.9 s = 3.11×** (per iteration 569.9 → 183.5 s), so the
standing against CP2K's 8.5 s/iteration goes **67× → ~22×**.


---

## Why CP2K was ~22× faster, read from its source — harvested 2026-09-08 (was OpenWork.md L957-1080)

*(verbatim)*

## ★★★ WHY CP2K WAS ~22× FASTER — READ FROM ITS SOURCE (and now acted on: see the rewrite plan)

> **USER, 2026-08-27:** *"I think it is time to investigate more deeply into what CP2K is actually doing."*

Source read at `/home/janr/Code/cp2k` (2026.1, `757bb76`), backend `src/grid/cpu/`.  **The gap is
ARCHITECTURAL, not micro** — every line below is a structural difference, and together they explain an
order of magnitude in a way no further tuning of our loop can.

**THE ROOT MOVE: CP2K NEVER CARRIES TWO GAUSSIANS TO A GRID POINT.**  `cab_to_cxyz`
(`grid_cpu_collint.h:1047`) binomially re-expands \f$(x-a)^{l_a}(x-b)^{l_b}=\sum_s\alpha_s(x-p)^s\f$ ONCE
per (shell pair, offset), collapsing the pair into **one** Gaussian \f$e^{-\zeta_p|r-P|^2}\f$ times **one**
polynomial in \f$(r-P)\f$.  We instead carry two Gaussians, two Cartesian monomials AND a loop over live
component pairs to every point.  Everything below follows from this one difference.

★★ **AND WE ALREADY HAVE THE COLLAPSE — IT IS `Ω` (user, 2026-08-27).**
`src/BasisSet/Molecule/Evaluators/PG_Cart_MnD/GaussianRF.C`, `struct Ω : public Cacheable2`, interned in
the process-global `Cache2` keyed on the primitive pair, carries exactly CP2K's four quantities:
\f$\alpha_p=a+b\f$ (their `zetp`), the product centre \f$P\f$ (their `rp`), the prefactor `Eij` (their
`prefactor`), and `H2` — the M&D expansion coefficients.  And `Hermite2` stores them as `d`,`e`,`f`
indexed \f$(N,n_a,n_b)\f$, \f$(L,l_a,l_b)\f$, \f$(M,m_a,m_b)\f$: **already PER-DIRECTION SEPARABLE**,
which is precisely the property that makes \f$\Lambda_{NLM}(r-P)\f$ factor into three 1-D functions.
⇒ **`cab_to_cxyz` is not something to write.**  Our starting point is better than CP2K's, not worse: they
re-derive the collapse per task, we have it cached as geometry.

⛔ **AND THE "FOLD THE DENSITY MATRIX IN UP FRONT" CAVEAT THIS SECTION FIRST CARRIED WAS WRONG — do NOT
copy that half.**  (It also said "D-block", which in a discussion of Mn *d shells* is an unforgivable
overload: read **density-matrix sub-block**, \f$D_{ij}\f$ over \f$i\in\f$ shell I, \f$j\in\f$ shell J,
never a block of d-type Gaussians.)  Three corrections:
1. **\f$D\f$ is the ONLY thing in the whole walk that changes between SCF iterations.**  \f$\Omega\f$,
   \f$\alpha_p\f$, \f$P\f$, `Eij`, `H2`, the boxes and the chords are all GEOMETRY.  CP2K folds `pab`
   in early because it re-derives the cube per task anyway; we do not, so the right split is: separable
   tables from \f$\Omega\f$ alone (cacheable across iterations AND across k-blocks), then contract the
   COEFFICIENT TENSOR with \f$D\f$ per iteration — \f$O(n_I n_J n_{NLM})\f$ per pair, **no grid points
   involved**.  ⇒ The no-cut discipline and the D-aware screen need not move at all, because \f$D\f$
   never enters the geometry object.
2. **It is not \f$D_{ij}\f$ that multiplies the pair product** but
   \f$\mathrm{fold}\cdot\mathrm{Re}[D_{ij}\overline{e^{ikR_n}}]\f$ — a REAL scalar per (pair,
   **offset**), since the Bloch phase rides on the image.  Any weighting is per-offset.
3. **The integrate direction has no \f$D\f$ at all.**  `IntegratePotential` PRODUCES \f$h_{ij}\f$;
   its `screenD` is a screening magnitude, not a weight.  The tensor structure still applies (grid →
   coefficients → block, CP2K's mirror), but "fold D in" was never a statement about that direction.
⚠ One distinction to keep: \f$D\f$ is over CONTRACTED, normalized components (`ns[i]`), \f$\Omega\f$
over PRIMITIVE pairs.  Every current basis is uncontracted so they coincide TODAY — `gi[p]` exists for a
reason and the mapping must not be assumed.

1. **ONE ISOTROPIC GAUSSIAN ⇒ THE EXPONENTIAL FACTORISES, so there are NO transcendentals per point.**
   - *Orthorhombic*: three 1-D tables `pol[dir][power][ig]` = \f$(x-x_p)^{l}e^{-\zeta_p(x-x_p)^2}\f$
     (`grid_cpu_collint.h:409`) — **O(n) exps for the whole cube**, and the polynomial power rides in the
     same table.
   - *General / triclinic — OUR CASE (FCC, rhombohedral)*: **"Mathieu's trick"**
     (`grid_cpu_collint.h:532, 752`).  The quadratic form factors into THREE 2-D tables,
     \f$e^{-\zeta_p Q(i,j,k)}=T_{ij}\,T_{jk}\,T_{ki}\f$ with
     \f$T_{ij}=e^{-\zeta_p(d_i^2h_{ii}+2d_id_jh_{ij})}\f$ — the cross terms distribute exactly.  Per grid
     point the 3-D Gaussian is then **three table lookups and two multiplies**.
2. **THE TABLES ARE BUILT BY THE MULTIPLICATIVE RECURRENCE WE JUST REJECTED** — the comment at
   `grid_cpu_collint.h:414` is literally *"Reuse the result from the previous gridpoint to avoid to many
   exps"*.  ⇒ **Our recurrence experiment was the right identity at the wrong LEVEL.**  They apply it to an
   O(n) (or O(n²)) TABLE built once per cube; we applied it inside an O(n³) point loop, so we still paid
   O(n²) seedings and still evaluated everything else per point.  They also seed *symmetrically outward
   from the cube centre* (`general_fill_exp_table`) — the numerically safe start we identified and skipped.
   And because their tables are per-DIRECTION-PAIR rather than along z only, the anisotropy that flipped
   `NaPseudoAtomInBoxDoublet` does not arise.
3. **THE CUBE IS A TENSOR CONTRACTION, NOT A POINT SWEEP.**  cxyz → cxy → cx → grid as nested 1-D passes.
   The innermost loop (`grid_cpu_collint.h:38`) is `lp+1` fused multiply-adds against the table, **four
   grid points at a time, under `#pragma omp simd`** — no exp, no monomial, no component-pair loop.  There
   is also a whole alternative backend (`src/grid/dgemm/`) that recasts the same contraction as BLAS calls.
4. **THE CHORD BOUNDS ARE CACHED, NOT RE-SOLVED.**  `grid_sphere_cache_lookup`
   (`common/grid_sphere_cache.h:43`) memoises the per-line sphere bounds by discretized radius and cell.
   We now compute the same interval (above) but solve a quadratic per line, every line, every pair.

**WHAT THIS MEANS FOR US, stated as a cost law rather than a wish.**  Per (pair, offset) cube of edge n:

| | transcendentals | per-point work |
|---|---|---|
| qchem today | \f$2n^3\f$ | 2 exps + 6 power tables + \f$(n_I{+}n_J)\f$ products + a loop over live pairs |
| CP2K (general cell) | \f$O(n^2)\f$ table entries, \f$O(n)\f$ exps | \f$l_p{+}1\f$ FMAs against a table |

⇒ **THE PLAN IS `doc/CollocationRewritePlan.md`** (2026-08-27), steps 0–8 with the gate first.
⇒ **There is no remaining 2× inside our current loop shape** — the exp experiment demonstrated that
directly (removing 20% of the profile bought 3%).  Closing the rest means adopting the SHAPE: the
product-centre re-expansion first, then separable tables, then contraction.  That is a real piece of work,
but it is a well-defined one, and CP2K's own file layout is a serviceable specification.

★ **AND IT ATTACKS BOTH BENCHMARK AXES WITH ONE CHANGE.**  The stream cache costs **4218 MB** because it
stores per-point VALUES — \f$O(n^3)\f$ per (pair, offset).  Separable tables are \f$O(n)\f$ (orthorhombic)
or \f$O(n^2)\f$ (triclinic).  So this route is the first one that can deliver CP2K's memory profile AND
its CPU, which is exactly what doc/Benchmark.md said was needed and what caching never gave.

⚠ **THE REAL DESIGN QUESTION IS THE SCREEN, NOT THE DENSITY MATRIX.**  Today the D-aware tolerance is PER
COMPONENT PAIR (`epsHere[k]`, applied to \f$|val|\f$ at each point) — and pre-filtering on it is what
made the shell hoist pay at all (1.03–1.08× without it against 2.13× with).  In a contracted cube the
per-component identity is GONE: one cube serves the whole shell pair.  So the screen has to move onto the
COEFFICIENT TENSOR, and whether that can be done without either loosening the no-cut discipline or
re-inflating the boxes is the question to settle BEFORE the kernel is written.  The T3 orbit fold (keyed
per (pair, offset), reading the orbit-projected \f$D\f$) needs the same re-derivation.

### WHAT WAS LEFT, and it is now CLOSED

**`exp` IS THE WHOLE REMAINING STRUCTURE: 29% of the NaF profile** (`__ieee754_exp_fma` 17.2% +
`exp@@GLIBC_2.29` 9.8% + 2.3% of plt stubs — note ~⅓ of that is the glibc wrapper, not the evaluation).
Everything else is flat: after the four edits no single instruction in the walk exceeds 3%.

⛔ **SUPERSEDED BY MEASUREMENT — see the rejected-experiment section above.  What follows was the
PREDICTION; the measured answer is 1.03× on MnO, and the remaining item is now the INTERVAL SKIP:**
under the sphere screen the accepted set on a line is a CONTIGUOUS dz interval (the screen is a convex
quadratic in t), so it can be solved analytically and the loop bounded to it.  That removes the wasted
`ri2`/screen work on the ~48% of points the screen rejects, it is **BIT-IDENTICAL** provided `r` is still
accumulated through the skipped region, and it needs no flag.  **Not yet built or measured.**

**The candidate was the handoff's #3**: along a run of collinear grid points
\f$e^{-\alpha(x+nd)^2}\f$ obeys a two-term MULTIPLICATIVE RECURRENCE, which is what CP2K's cube-mapped
collocate exploits.  Two things to weigh before starting:
- ⚠ **It cannot be bit-identical** (accumulated rounding along the line), and unlike the sphere screen it is
  not a sub-ε truncation — it is a different arithmetic for the SAME term.  That needs the user's call, and
  a drift bound measured against the direct form rather than asserted.
- ⚠ **It is uniform-ladder-only** (user, 2026-08-26) — the Becke Φ path has no collinear runs to walk.
  Given the point above, that costs nothing here: the Becke path is singles + GEMM and is already ~3%.
- Ceiling: if exp became free the walk would go 2.53 → ~1.8 s on Si, i.e. ~1.4×.  **Worth doing, but it is
  the last 1.4×, not another 2×** — the cheap structural fat is gone.

**Two smaller things the profile still shows** (both grid-agnostic, both ~2%): `IVec3Less` at 1.9% (the
handoff's candidate 4 — a `std::map<Vector3D<int>,…>` comparator, in `Projector3`/`SymmetrizeGMap`, NOT in
the box walk), and `qcMesh::BeckeCutoff` at 3.2% on NaF.


---

## N2 — the density G ball — harvested 2026-09-08 (was OpenWork.md L1298-1393)

*(verbatim)*

## ✅ N2 — THE DENSITY G BALL: RESOLVED (the lobes are BAND-LIMITING)

\f$\rho=\sum_{ij}D_{ij}\chi_i\chi_j\f$ is a product of Gaussians, so its G content runs to
\f$2\alpha_{\max}\f$ and \f$C=2\f$ is the NAIVE Nyquist floor — which is what MnO has always run
(`densityEcut 72 Ha` against `alpha_max 36`).  But a ball that only just reaches the product's own exponent
still ALIASES the cusp tail, and `GPW_IBS.C:53` already records the consequence: *"under-resolving it
aliases rho into large spurious negative lobes → the XC collapse … F needed ~8*alpha_max for negCharge
−9 e → −0.03 e."*  Measured here at C=2: **15.2% of the 97160 Becke points carry ρ<0** (min −0.154), and
`SlaterExchange::GetVxc` guards `ro>0`, so a sixth of the atom-centred quadrature contributes NOTHING to
\f$E_{xc}\f$.

**⇒ ✅ RESOLVED 2026-08-25 — AND THE MECHANISM IS THE OPPOSITE OF WHAT I FIRST WROTE.**

**(a) THE CONTROLLED EXPERIMENT: `RasterPolicy`, WITH {G} HELD FIXED.**  `RasterPolicy::AliasFree` chooses
divisions so the product of ANY two ball waves is sampled exactly; `BallOnly` (CP2K's bet, and what GPW
actually runs — `BasisSet.C:326`) resolves the ball alone and accepts the fold-back.  Crucially the BALL IS
IDENTICAL in both — nG 8623, 870 G-stars — only the raster differs (level-0 40³ vs 80³, 8× the points).
That separates ALIASING from BAND-LIMITING, which a `cutoffFactor` sweep cannot, because C moves both.

| matched end-of-stage-1 | BallOnly (prod) | AliasFree | ratio |
|---|---|---|---|
| points ρ<0 | 15464 (15.92%) | 15380 (15.83%) | 0.99× |
| negative mass | −1.812e-3 e | −1.816e-3 e | 1.00× |
| min ρ | −0.1460 | −0.1468 | 1.01× |
| Etot | −61.40297621 | −61.40297359 | **2.6 µHa** |
| CPU / RSS | 500 s / 1316 MB | 1155 s / 4652 MB | **2.3× / 3.5×** |

**8× the raster points moves the lobes by 1%.**  ⇒ they are **GIBBS RINGING from band-limiting**, exactly
what a truncated Fourier series does at a cusp — NOT fold-back.  ✅ **By-product: this closes the
calibration the code itself flags as never taken** (*"the A/B is the measurement"*, `Evaluator.C:250`):
**`BallOnly` is vindicated as the default.**

**(b) THE `cutoffFactor` SWEEP — real, but my mechanism reading was BACKWARDS.**  Matched end-of-stage-1:

| C | Ecut | pts ρ<0 | negative mass (e) | min ρ | Etot | st2 it | CPU / RSS |
|---|---|---|---|---|---|---|---|
| 2 | 72 | 15.9% | −1.81e-3 | −0.146 | −61.40297621 | 17 | 500 s / 1316 MB |
| 3 | 108 | 14.2% | −1.88e-4 | −0.0366 | −61.40282017 | **12** | 619 s / 1973 MB |
| 4 | 144 | 6.6% | −1.65e-5 | −0.0071 | −61.40283283 | 13 | 814 s / 2638 MB |
| 6 | 216 | 5.6% | **−3.02e-7** | −5.49e-4 | (stopped after st 1) | — | — |

Negative MASS falls ~geometrically (**6000× by C=6**) while the point COUNT flattens at 5–6% — ringing
amplitude decaying with its support unchanged.  C works because it **ENLARGES {G}** (less band-limiting),
NOT because it reduces aliasing; (a) proves that directly.
★ **Two findings worth acting on independently of N4:** production **C=2 carries ~0.15 mHa of grid error**
(C=2→3 moves 0.15 mHa, C=3→4 only 1.3e-5 Ha), and **C=3 CONVERGES FASTER** (12 vs 17 stage-2 iterations)
for 1.24× CPU.  The C=2 default deserves a review on its own merits.

**(c) ⛔ BUT C DOES NOT SCALE — AND THAT IS THE DECIDING ARGUMENT (user, 2026-08-25).**  *"Imagine doing
Li_0.125MnO2 in a 2x2x2 supercell of LiMnO2, with C=4."*  Raster cost is per unit VOLUME, so an 8× supercell
gives ~10.3 GB (C=2) / 15.4 GB (C=3) / 20.6 GB (C=4) of raster ALONE on a 14 GB box — **even today's C=2
does not fit**, before the pair streams.  Raising C pays a **GLOBAL** price (the whole cell volume) to fix a
**LOCAL** defect (cusps within a fraction of a bohr of each nucleus).
*(Caveat: the measured RSS exponent is C^1.00, not the ideal raster's C^1.5, so the C-scaling is the soft
part of that estimate; the 8× volume factor is solid and C=2 alone already overflows.)*

**⇒ NET: EVERY CHEAP ALTERNATIVE TO N4 IS ELIMINATED.**  Raster geometry does nothing; a bigger ball works
but cannot scale.  The only remaining cure is to give \f$V_{xc}\f$ content the ball CANNOT represent — the
cusp deficit, whose support is LOCAL (O(n_atoms × a small ball), linear in system size).
⛔ **RETRACTED, 2026-08-25 (user asked the right question: "is this about ρ̃ or the DM ρ?").**  An earlier
cut of this section said the grid ladder "already knows which content is sharp" and that the XC feed
"throws that away", implying N4 might be an ASSEMBLY of existing machinery.  **That is wrong.**  The levels
are FFT'd individually and \f$\tilde\rho\f$ is *"combined NESTED in G-space"* (`GPW/Evaluator.C:302`), so
coarse content is embedded into the finer ball and the result lives on the **RUNG's** ball — the sharpest
level there is.  Verified: the Kerker line's G-count equals the rung's nG on every run (C=2: 24263 at
ecut 144; C=3: 44873 at 216; C=4: 69261 at 288; AliasFree: 24263, same ball, finer raster).
**Nothing is discarded in the assembly.**  The band limit is that ball's radius, and the cusp content above
it was never captured at ANY level — **absent, not thrown away**.
⇒ There is no "stop collapsing the ladder" shortcut.  The missing content can only come from something with
NO grid at all, i.e. the analytic \f$\rho[D]\f$ — which **strengthens** N4 and makes it a BUILD, not an
assembly.
⚠ Note also: the production row's EFFECTIVE density cutoff is **144 Ha (4·α_max)**, not the nominal 72 —
the top rung extends it — so the C sweep's cutoffs are 144/216/288.  (Trend unaffected.)

★ **WHERE THE PAIRS ACTUALLY ARE** (user, 2026-08-25: *"I have lost track of where the pairs are used if
not in evaluating DM ρ(r) for Vxc"*).  `MakeCollocator(coulomb, grid)` has two outputs: `coulomb=false` →
\f$\tilde\rho(G)\f$ (what Kerker MIXES, and what the matrix-free XC route inverse-FTs), and `coulomb=true`
→ \f$V_H\f$ with \f$4\pi/G^2\f$ folded in (the POISSON solve).  Its adjoint `MakeIntegrator` turns
grid-borne potentials into \f$\langle\chi_i|V|\chi_j\rangle\f$ and backs the Overlap3C/Repulsion3C tensors.
**So pairs are the HARTREE/POISSON path and the 3-centre tensor path — XC never uses them** (it is Φ-table
SINGLES on the Becke mesh already).  Cost on the production row: scatter 41.3 s + gather 22.8 s + stream
build 27.3 s = **91 s of 500 s (18%)**, which bounds what any pairs→singles rework can win *on this cell*.
⇒ **THE OPEN QUESTION IS NARROWER THAN "DO WE NEED PAIRS":** it is whether the RASTER/Poisson path can go
singles too, given \f$D=LL^\dagger\f$.  ✅ NEW EVIDENCE: the Cholesky orbitals ARE localised — measured
**IPR = 3.49 effective basis functions of 118** (range 3.23–3.95 over 154 factorisations), which satisfies
the code's own stated criterion for the singles route being viable.  ⚠ Three counterweights: (a) the WEIGHT
is localised but the SUPPORT is not (73 of 118 functions carry |L|>1e-3, 90 carry >1e-5, and collocation
accuracy is set by `GPW_DENSITY_EPS=1e-10`, i.e. by the tails); (b) the enumeration is (i,R) vs (i,j,R) and
pairs screen EXPONENTIALLY in separation via the Gaussian product theorem while orbitals do not; (c) ★ **an
orbital MIXES ALL EXPONENTS, so it has no single bandwidth and cannot be assigned to a ladder level** — it
needs the FINEST grid everywhere, which dissolves the multigrid's entire economy (levels 1–4 are 20³/9³/5³/3³
against level 0's 40³).  ⇒ it would trade 1909 folded pairs *on cheap levels* for ~13 orbitals *all on the
fine grid*.  **The census the plan already demands** — (i,R) vs (i,j,R) at matched ε, weighted by box
volume — is what settles it, and the answer likely FLIPS at battery-supercell sizes, where screened pairs
grow ~N² while occupied orbitals grow ~N.


---

## Step 0 — fix the instruments — harvested 2026-09-08 (was OpenWork.md L1444-1598)

*(verbatim)*

### Step 0 — FIX THE INSTRUMENTS  ·  ~½ day  ·  DO THIS FIRST

Two readouts, both cheap, both prerequisites for Step 1 (a benchmark table that bakes in a wrong
instrument is a table we redo) — **both now DONE**, plus **0c**, added 2026-08-20 and still OPEN:
the console reports WHAT but not WHEN, which is a prerequisite for the pre-SCF work rather than for Step 1.

**0a — REPLACE the order parameter with an INTEGRATED site moment.**

> **STANDING RULE (user, 2026-08-19): always report a proper INTEGRATED, real-world-measurable quantity.
> Never a point sample of a field.**  It applies beyond this one probe — anywhere we are tempted to
> characterise a solution by evaluating a field somewhere.

What we report today is `m_stag = ½[m(Mn1) − m(Mn2)]` with `m(r)=ρ↑−ρ↓` evaluated at **one point, 0.7 bohr
off the nucleus, along +x** (`IntegrationTests/GPW_SCF_UT.C:3624`, `off(0.7,0,0)`).  Three things are wrong
with it:

1. **It is a spin DENSITY (e/bohr³), not a moment.**  CP2K's 4.654 is a *Mulliken site moment in μB*
   (`doc/CP2Kresults.md:153`, `doc/SymmetryUpgradeHistory.md:252`) — integrated and basis-partitioned.  So
   "qchem 0.67 vs CP2K 4.65" was never a like-for-like statement.  The same run already prints
   `|m̃(q_AFM)|·Ω/2 = 3.126 e⁻`, ~4.8× the point probe and much nearer CP2K's scale; earlier sessions did
   compare *that* to 4.65 (`doc/SphericalLatticePlan.md:124`).
2. **0.7 was never derived.**  The only justification in the tree is the parenthetical "(the d-shell peak —
   the d density VANISHES at the nucleus)".  The motivation is sound; the number is asserted.
3. **★ It samples ONE DIRECTION, which for a d shell is the confounder itself.**  A cubic-split d spin
   density is not spherically symmetric — e_g and t_2g have very different amplitude along [100] vs [111] —
   so this probe responds to the ORBITAL OCCUPATION PATTERN, not just the moment.  In a campaign whose
   central difficulty is occupation hopping among near-degenerate d configurations, the instrument moves
   when the orbitals rotate even if the moment does not.

**What it can and cannot support.**  As a COLLAPSE DETECTOR it is fine and most of its historical use is
safe — "m_stag 0.366 → 0.0046 at iteration 1" is a real finding in any units.  It does NOT support
(a) comparison to CP2K, (b) magnitude claims like "the weak-moment basin" / "m_stag ±0.667", or
(c) site-asymmetry at the few-percent level (`Mn1=+0.3058, Mn2=−0.5849`, `SymmetryUpgradeHistory.md:855`,
which fed run 59's flip-defect reading) — two sites can differ in ORBITAL occupation at identical moment.

> **THE PROPER DEFINITION IS BADER, AND IT IS A FUTURE FEATURE (user, 2026-08-19).**  The physically
> correct atomic basin is R. F. W. Bader's QTAIM basin — the region bounded by the ZERO-FLUX surface
> \f$\nabla\rho(r)\cdot n(r)=0\f$ — which is a property of the density itself, not of a chosen partition
> function, and is what makes an "atomic" moment or charge well defined rather than conventional.
> Implementing it (gradient-path/zero-flux basin assignment on the mesh, then \f$\int_{\Omega_A} m\f$) is
> its own increment, filed here as a wanted feature.  **For now the Becke fuzzy basins are good enough** —
> they are integrated, they carry units, and they are a vast improvement on a point sample.  Whatever we
> report must say WHICH partition it used, because until Bader lands the number is partition-dependent.
> (Note the run's `|m̃(q_AFM)|·Ω/2` is partition-FREE and is the closest thing we print to a
> neutron-diffraction observable — worth keeping alongside for exactly that reason.)

**The work:** a **Becke-partitioned site moment** \f$\int w_A(r)\,m(r)\,d^3r\f$ becomes THE order parameter
(we already build a site-adapted Becke mesh with per-atom partition weights for XC, so it is nearly free
per iteration), plus a **Mulliken site moment** for the exact CP2K-comparable number.  Any surviving point
probe is renamed so it cannot be read as a moment and carries its units.  The actual effort is PLUMBING:
`GpwOptions::orderProbe` receives only a `cDM_CD&` and does two point evaluations, so a partitioned
integral needs the mesh and its per-atom weights to reach that seam.
**★★ MEASURED 2026-08-19 — THE MOMENT GAP WAS THE INSTRUMENT.**  With `XC_SinglesQuadrature` now reporting the
Becke-partitioned \f$\int w_A(\rho_\uparrow-\rho_\downarrow)\f$ (`QCHEM_SITE_MOMENTS=1`), the MnO magnetic
cell reads, in ELECTRONS, per iteration:

| | Mn1 | Mn2 | O | O |
|---|---|---|---|---|
| seed | **+4.781** | −4.781 | 7e-10 | 4e-10 |
| iter 1 | **+4.663** | −4.653 | −0.006 | −0.003 |
| iter 2 | +1.871 | −1.763 | −0.087 | −0.022 |
| iter 3 | +4.208 | −4.172 | −0.024 | −0.011 |
| iter 4 | +3.631 | −3.621 | −0.005 | −0.005 |

…while the point probe on the SAME run reported `m_stag = 0.318`.

**The Mn moment is ~3.6–4.8 e, not ~0.3.  CP2K's Mulliken is ±4.654.**  The seed is a clean high-spin d⁵
(+4.78, O ≈ 0, net ≈ 3e-10 under the fixed n↑=n↓), which is exactly what a Mn²⁺ should be — so the
instrument is behaving.  **The "≈7× moment discrepancy" was essentially all units**, and the
"weak-moment basin" that shaped a large part of the MnO campaign is at minimum badly overstated by an
instrument reading ~0.3 where the observable is ~4.

Honest caveats: (a) this is a short unconverged run and the trajectory is still bouncing
(4.78 → 4.66 → 1.87 → 4.21 → 3.63), so it does NOT establish the converged moment — only its SCALE;
(b) Becke ≠ Mulliken, so 4.66 vs 4.654 agreeing to three digits is partly luck — what is solid is that
both are ≈4.7 for a high-spin d⁵.  **Historical `m_stag` numbers in the plan docs are collapse/survival
evidence only; they are not moments and must not be quoted as such.**  Re-reading the campaign's
"weak-moment" conclusions against the integrated observable is now its own item under Step 5.

**0b — make the FOLDS visible.  ✅ DONE 2026-08-19.**  One shared `qchem::report::EmitFold(site, nOps, raw,
reps, note)` — `[fold] <site>: N ops, raw -> reps = F×` on cout, plus a `folds.<site>` report entry — wired
at the three sites that were completely dark.  Run-scoped dedup lives in the reporter (several sites fire
once per k-block; eight identical lines train the reader to skip them), and ARMED-NESS is read from the
REDUCTION, not from an op count, because a `Fold` carries orbits and not the op set that made them.
Measured, first time any of this was on the console:

| site | free run | imposed |
|---|---|---|
| XC mesh (Becke star-average) | NONE | **21.6×** (Si IBZ) |
| `V_loc` {G}-star | NONE | **26.1×** (48 ops) |
| collocation streams (T3 pairs) | NONE `[opt-in: GPW_STREAM_FOLD=1 — armed by default since T3.5, Step 2]` | **12.5×** (48 ops) |

**And the message that matters for Step 2, now unmissable:** the production MnO magnetic run prints
`NONE` on all three — it folds NOTHING, on a cell whose magnetic group has 12–24 ops.  That is the
12–48× sitting unclaimed, said out loud on every run instead of inferred from a plan doc.

**0c — MAKE THE ORDER LEGIBLE.  THE CONSOLE CANNOT TELL YOU *WHEN* (user, 2026-08-20).  ✅ BUILT 2026-09-06.**

> ✅ **WHAT LANDED.**  `report::RunElapsed()` — the run's own `steady_clock`, zeroed at `Begin`, inert
> outside a run — stamps every emitted item.  On the console: `meta  [t=0.06 s]`,
> `grids ▸ xcQuadrature  [t=0.06 s]`, and a scoped `Section` prints its whole SPAN
> `assembled  [t=3.21→12.34 s]`, because that render fires at scope CLOSE and the span is the only honest
> thing to print there.  `Log`'s heartbeat and every `[fold]` line carry the stamp too.  In the record: a
> top-level `timeline` array of `{item, [t0,] t}`, run-scoped (a nested SAD bootstrap keeps its own) and
> dedup'd through `EmitAt`'s idempotence, so a provider that re-announces the same grid per k-block does
> not flood it.  **The design question this item posed — chronological stream vs nested render — is
> answered by KEEPING BOTH**: the nested document is unchanged, the chronological one is an index beside
> it, and that costs one array.  5 unit tests in `UTCommon` (`src/Common/tests/Reporting.C`).
>
> ⚡ **AND IT EARNED ITS KEEP ON THE FIRST RUN IT WAS USED ON.**  The MnO ledger hunt
> (`doc/ParallelAndOraclePlan.md` 1.1(a)) was looking for 25 s of unbucketed time.  The stamps found it
> before any bucket was read — 60.17 s at the end of anneal stage 1, 95.93 s at the first line of stage 2,
> **35.8 seconds in which the run printed nothing at all**, which is exactly the stage-2 Hamiltonian
> rebuild.  The item's own claim ("the GAPS BETWEEN SECTIONS are the unbucketed time") is now measured,
> not argued.
>
> ⏸ **NOT built, and deliberately**: the *"constructed X"* trace at real construction points (the sibling
> idea below).  The stamps made the pre-SCF sequence legible enough to find the block that mattered; if a
> future hunt needs finer resolution inside one silent gap, that is the increment to build then.

> **The user's view into a run is the console output.**  It says WHAT happened; it does not say WHEN, and
> today it can actively mislead about it.

**The evidence, from the session that raised it.**  A run printed `grids ▸ xcQuadrature kind Becke` and,
further down — *after the iteration header* — `grids ▸ xcQuadrature kind Uniform`.  Chasing the second one
cost the better part of an hour, and the sharpest wrong turn was this: **a report section RENDERS when its
ENCLOSING section closes, so a block's POSITION in the console is not its construction time.**  Reading
"early" and "late" off the log produced a lazy-construction story that was simply false — the object was
built where the log implied it was not.  A runtime BACKTRACE settled it in one shot.  (The key collision
itself is fixed; see the `vxcFitGrid` commit.  This item is the general defect it exposed.)

**The proposal (user):** give **each report item a TIMESTAMP**, and let the renderer OPTIONALLY guarantee
that items stream out in the true order — so ordering is *read*, not inferred.  Two things fall out of it
that are worth stating separately, because they are separable increments:
- a **monotonic timestamp per item** is the minimum, and it is nearly free: `Timed` already reads
  `steady_clock` (`Common/Imp/Reporting.C`), so the same clock can stamp every `EmitAt`/`Emit`;
- an **order-preserving render mode** is the part with a design question: sections currently nest, and a
  strictly chronological stream and a nested-section render are two different documents.  Options are to
  stream chronologically with the section as a FIELD on each item, or to keep the nesting and print each
  item's stamp so a reader can reconstruct order.  Deciding that is the increment.

**Sibling idea, complementary not alternative:** a one-line *"constructed X"* trace at the real
construction points would make the **pre-SCF sequence** legible directly — which matters because there are
known problems in that sequence, and today the only way to establish what runs when is to instrument and
dump a stack.  An instrument that makes ORDER visible is a prerequisite for hunting them, exactly as 0a/0b
were prerequisites for Step 1.

- **0d — THE SELF-DESCRIBING BENCHMARK BANNER (added 2026-08-25, user).**  qchem has no equivalent of
  CP2K's `GLOBAL| Number of threads` line, so every row in `doc/Benchmark.md` asserts its thread state by
  hand — and `doc/Benchmark.md` has now asked for this twice.  Print, at run start: the OpenMP thread count
  (`GPW_OMP_THREADS`), the BLAS thread state, and the **qchem-only accelerator flags** (`QCHEM_DM_LOWRANK`
  — ON by default — and `GPW_XC_DM_SOURCE`, plus the T3 fold state).  That makes a row self-describing
  instead of relying on the person taking it, which is the discipline that keeps failing.  See
  `doc/Benchmark.md` → **BENCHMARK PROTOCOL** for the three rules this instrument serves.


---

## Step 1 — the head-to-head table — harvested 2026-09-08 (was OpenWork.md L1599-1662)

*(verbatim)*

### Step 1 — THE HEAD-TO-HEAD TABLE, built as a standing benchmark  ·  **table: `doc/Benchmark.md`**  ·  ✅ THE TABLE STANDS

**Both columns are now MEASURED, on this box, through one wrapper** (`scripts/bench`, 2026-08-19, 1 thread
each).  The CP2K half stopped being banked prose the moment the packaged CP2K 2025.2 was validated against
five banked 2026.1 decks — all five reproduce to the printed digits (`doc/CP2KBuild.md`).  Seven rows carry
energy + wall + peak RAM on both sides; three cells remain (below).

| | Δ(E) | CPU q/c (wall) | RAM q/c |
|---|---|---|---|
| Si Γ | −10.6 µHa | 1.5× (6.1 s / 5.2 s) | 267 / 148 MB |
| Si 2×2×2 Γ-centred | −14.9 µHa | 1.7× (8.3 s / 5.8 s) | 269 / 153 MB |
| NaF Γ (SR2) | +0.877 mHa | **13.2×** (39.5 s / 7.4 s) | 577 / 173 MB |
| NaF Γ (full SR) | +1.349 mHa | 2.1× (2m44 / 1m42) | 3090 / 186 MB |
| **MnO AFM-II Γ (VA)** | **−99.65 mHa** | **6.0×** (20m05 / 6m14) | **4947 / 217 MB** |
| **MnO FM Γ (VA)** | **−136.80 mHa** | **12.1×** (21m45 / 3m13) | **4947 / 217 MB** |

**★ THE HEADLINE THE TABLE EXISTS TO PRODUCE: on MnO we are 6–12× the CPU time and 23× THE RAM** — 4947 MB
against CP2K's 217 MB on a function-for-function identical 118-basis.  That is the number Steps 2–4 are
now measured against, and it says the RAM half is not a rounding detail.
**⚠ COMPARE CPU TIME, NOT WALL** (user, 2026-08-19): CP2K is genuinely serial here, qchem is not — blaze
threads the BLAS whatever `GPW_OMP_THREADS` says (measured 115–239% CPU), so the wall column flatters us by
the cores taken and the true ratio is about 2× worse than a wall reading.  A knob is not a measurement.
**Loose end worth a look:** NaF SR2 is 13.2× while the LARGER full-SR span is 2.1× — same cell, same Γ, same
path, and the smaller basis is the worse ratio.
**And the profile named its own next lever:** the largest bucket in the 20-minute AFM run is the **XC-mesh
Φ-table build at 370 s (31%)** — bigger than either pair loop, and exactly Step 3's Φ-screening item.  The
scatter/gather that rounds 3–4 were spent on is 22% + 14%.

**Reproducibility repairs made along the way** (each was silently corrupting rows):
- **The `[fold]` line left `cout` at precision 2** for the rest of the run — `std::defaultfloat` restores the
  format flag, not the precision.  THAT is why energies printed as `Etot=-7.1`; the earlier table recorded it
  as a "detail level".  The fingerprint's `Efinal` had the same disease from the other direction: it printed
  12 digits in run 61 only because a verbose table had set `fixed(10)` upstream.  Both now state and restore
  their own precision, and the run summary prints Etot at 10 s.f.
- **The annealed driver reported no energy summary, no ledger and no PEAK RSS** — and every MnO row runs
  through it, so "every GPW run reports PEAK RSS" was false for exactly the runs whose RAM the table needs.
- **The VA span came from a script OVERWRITING a committed basis file** (`bisect_valence_sph.py` →
  `valence_lowq_sph.bsd`), so no MnO row was reproducible after the file was restored.  VA/VB are now
  committed basis sets selected by `GPW_BASIS_SPAN=va|vb|sph|sr`; the run prints its own `nFunctions`, and
  the re-run reproduced run 61 exactly (118 functions, λ_min 1.29e-3, E = −61.4029762 vs −61.4029762007).
- **`DISABLED_NaFRocksaltGamma` ran a 2×2×2 mesh against a Γ-era anchor** (and against Γ-only CP2K decks);
  it simply failed.  `NAF_KMESH` now selects the mesh and the anchor follows it (`NAF_SPAN=sr|sr2` likewise).

**The threaded repeat is deferred ON PURPOSE** (user, 2026-08-19): this cut is the serial baseline, the whole
table gets re-run at 12 threads **after** Steps 2–3 bring the qchem times down — re-measuring a 20-minute row
that is about to change by an order of magnitude buys a number with a short shelf life.  Keep both cuts when
it happens: serial = the algorithmic comparison, threaded = the user-facing time, ratio = parallel
efficiency, which nothing measures today.

**★ AND IT IMMEDIATELY PAID FOR ITSELF: the Si shifted-MP row was BROKEN, and the bug is now FIXED.**
`SR_2x2x2ShiftedMP_vs_CP2K` had rotted to −3.7351 against its −7.86744 anchor while sitting DISABLED — it is
the suite's ONLY fractional-k SCF coverage, so nothing caught it.  Root cause: the D-aware integrate-back
screen tested `|Re(D_ij·conj(phase))|` **as if a real part were a magnitude**.  At a quarter-integer k the
Bloch phase is purely imaginary on every odd offset, so that test discarded every odd-offset term and the
Hartree/XC matrix came out EXACTLY REAL (`maxIm(dV)=0` at k=¼ vs 0.067 next door); an H missing its
imaginary part has the wrong spectrum, hence 2.5 Ha.  Fixed by screening on the true magnitude `|D_ij|` —
which is what this project's own **"the magnitude screen is the only truncation"** rule always meant.
The row now reads −7.868473428 vs CP2K −7.867436530 (**1.04 mHa**), the test is ENABLED at ~14 s, and Si Γ /
Si 2×2×2 Γ-centred / NaF Γ are unchanged to every digit (TRIM k has real phases, where old and new agree).
Full detail in `doc/Benchmark.md` footnote ¹.

**Still open on this table:** the `MNO_KMESH=2` multi-k MnO row (cost unmeasured; CP2K needs a matching
`&KPOINTS` deck) and a k-point CP2K deck for NaF.


---

## Step 4 — RAM — harvested 2026-09-08 (was OpenWork.md L2215-2223)

*(verbatim)*

### Step 4 — RAM  ·  read it off Step 1's table, then decide

**✅ LARGELY ANSWERED BY STEP 2, as predicted — do not open it as a track.**  Arming the T3 fold took the
MnO AFM row's peak RSS **4947 → 1349 MB** (23× CP2K → 6.2×) on top of round 3's 5.78 → 3.70 GB, because the
RAM was the streams and folding cuts their *demand* by the orbit factor.  What remains is the Φ tables
(Step 3's screening item, which is a RAM lever as much as a time one) and the fact that a FREE run folds
nothing and therefore still pays the full 4947 MB.  Re-read the number off the table after Step 3; reopen
this only if it still binds then.


---

## Evidence dossier — factoring D — harvested 2026-09-08 (was OpenWork.md L2316-2786)

*(verbatim)*

## EVIDENCE DOSSIER (no action here) — fast evaluation of DM-ρ(r) by FACTORING D

> **This section and everything under it — Q1/Q2/Q3, the tier-0 results, the spectrum finding, the LSP
> design ruling — is the worked evidence for INDEX ITEM 2 (Step 3's low-rank-D ρ GEMM).  It is not a
> separate thread and there is nothing to start here.**  Read it when you build that item, or when you are
> about to re-propose something it already refuted.


**The idea.**  \f$\rho(r)=\Phi(r)^\dagger D\,\Phi(r)\f$ is evaluated everywhere in this code as a quadratic
form in D — a sum over PAIRS of basis functions.  D is a density matrix, so it is PSD with rank = the
occupied count, not n.  Factor it once, \f$D=LL^\dagger\f$ with L an (n × r) pivoted-Cholesky factor, and
\f[ \rho(r)=\lVert L^\dagger\Phi(r)\rVert^2 \]
which is a sum over SINGLES.  Two things change at once: the cost drops from O(npts·n²) to O(npts·n·r), and
the fundamental object stops being a pair.

**What is already MEASURED (2026-08-20/21), so this is not speculation:**

| | |
|---|---|
| numerical rank of the real mixed D, MnO, per spin | **14–17** against n=118 |
| ⇒ arithmetic ratio n/r | **7.0–8.4×** |
| rank stability, tol 1e-6 … 1e-12 | 14–17 — a CLEAN GAP, so the cut is read off, not tuned |
| \f$\lambda_{\min}/\lambda_{\max}\f$ | −1e-16 ⇒ **PSD to roundoff** |
| exactness | LAPACK's own pivot floor ⇒ residual at ROUNDOFF; A/B at IDENTICAL energy |
| implemented | `LowRankFactor` + the thin GEMM in `IrrepCD_Core::DM_RhoAtPoints` (`QCHEM_DM_LOWRANK=0` opts out) |
| instruments | `GPW_DM_RANK=1` (rank + PSD census) |

**Why D is PSD here, since the whole thing rests on it:** only `LinearMixer` touches D and its α is clamped
≤ 1, so it is a CONVEX combination of PSD matrices; the extrapolating mixers (Kerker, Pulay, Broyden) are
FIELD mixers and never touch D.  The guard is not assumed — `LowRankFactor` checks conserved mass
(\f$\lVert L\rVert_F^2\f$ vs Tr D) and returns false, keeping the caller on the full path, because LAPACK
`pstrf` does NOT error on an indefinite matrix; it stops early and reports a rank, which would silently
truncate ρ.

### Q1 — can we use this for Vxc IN COMBINATION WITH MIXING?

**Partial answer: yes, and there are three routes, but the cheap one is not obviously the accurate one.**
The Hamiltonian is built from the MIXED density, and on Kerker/Pulay recipes that density is a G-space
field with no D — which is exactly why XC currently falls back to sampling \f$\tilde\rho_{mix}\f$ by direct
summation at every mesh point (35 s / 6 iterations; the largest per-iteration cost).

- **(a) Mix D itself.**  `LinearMixer` is convex ⇒ PSD ⇒ everything stays factorable and exact.  Loses
  Kerker's G-dependent preconditioning, which NaF's low-G slosh and MnO's AFM basin lean on.  One-line A/B.
- **(b) Feed XC \f$\rho[D]\f$ alone.**  Cost collapses to the GEMM (~20×).  XC then sees \f$\rho_{out}\f$,
  not \f$\rho_{mix}\f$ — a TRAJECTORY change; at the fixed point \f$\rho[D]=\rho_{mix}\f$, so the converged
  answer is unchanged.  Cheapest, and the least justified.
- **(c) Cusp restoration.**  \f$\rho_{XC}=\rho_{mix}^{BL}+(\rho[D]_{exact}-\rho[D]_{BL})\f$ — the mixed
  density plus the sharp content band-limiting destroys, which is the content the atom-centred mesh exists
  to integrate.  The deficit is nearly ITERATION-INVARIANT (core electrons barely move) and → 0 at
  convergence.  **Cost is NOT automatically better**: it still needs one G-space sampling, so it is today's
  cost PLUS a GEMM unless the correction's {G} truncates.

  **The algebra says it does not truncate trivially.**  With Kerker,
  \f$\tilde\rho_{mix}-\tilde\rho[D]=(\alpha f(G)-1)\tilde\delta\f$, which at high G tends to
  \f$-(1-\alpha)\tilde\delta\f$ — 75% of the residual's high-G content at α=0.25.  What rescues it is that
  the whole correction is ∝ the SCF RESIDUAL, so an **adaptive G-ball keyed to \f$|\tilde\delta|\f$** is
  cheap late and accurate early, and the converged answer is \f$\rho[D]_{exact}\f$ either way.
- **GATING INSTRUMENT (build first):** the radial spectrum of \f$(\alpha f-1)\tilde\delta\f$ vs |G| per
  iteration — computable ENTIRELY INSIDE THE MIXER, since \f$\tilde\delta\f$ is already formed there for
  `ApplySpectralFilter`.  No mixer↔XC plumbing needed to answer it.
- ⚠ Honest note: (c) does NOT give Hartree and XC the identical array, as first claimed here.  The defence
  is that Poisson is LINEAR and diagonal in G (band-limiting converges fast) while XC is a NONLINEAR
  POINTWISE functional needing the cusp — each term gets the representation its operator requires.

### ✅ TIER-0 PRECURSOR RESULTS (2026-08-21)

**(1) IS L LOCALIZED?  MEASURED ON MnO — peaked, but NOT compactly supported.**  `GPW_DM_RANK=1` now also
reports the factor's structure.  Per orbital, over 4 iterations: **IPR (effective basis functions carrying
it) = 3.2–3.9** of n=118 — the WEIGHT is on 3–4 functions — but the coefficient decay is slow:

| \f$|L|>10^{-t}\cdot\max\f$ | t=1 | t=2 | t=3 | t=4 | t=5 |
|---|---|---|---|---|---|
| mean # functions | 8.6 | 38.3 | 64.7 | 82.6 | 89.6 |
| ≈ atoms (29 fns/atom) | <1 | 1.3 | 2.2 | 2.8 | 3.1 |

**★★ AND LOCALITY WAS NEVER THE GATE (user, 2026-08-21): "13 delocalized orbitals is a big win."**  Correct,
and it demotes everything below.  With r≈17 the object count alone carries the idea: **17 orbitals against
8778 (i,j,R) pair terms**, a dense **GEMM** replacing a SCATTER-bound emit loop (round 4 measured the pair
path at *"61% irreducible per-(pair,point) emit"*; a GEMM has no scatter), and a cached table of **O(grid·n) ≈ 60 MB against
1.35 GB of pair streams**.
⚠ **CORRECTION on WHAT IS CACHED (user, 2026-08-21): the ORBITALS cannot be cached, Φ can.**  D changes
every iteration, so the factor L (or \f$U\sqrt\lambda\f$) changes with it and the r orbitals are NOT
geometry-fixed.  What is cacheable is the **npts × n basis table Φ**, built once; each iteration CONTRACTS
it, \f$\Psi=\Phi L\f$ (npts × r), and row-norms.  So storage is O(grid·n) ≈ 60 MB for a 40³ grid — not the
O(grid·r) ≈ 8.7 MB first written here.  **This makes the proposal exactly "run the collocation grid the way
the XC mesh already runs":** Φ-cache + per-iteration GEMM is the existing `XC_SinglesQuadrature` pattern, so the
machinery is not new.  A delocalised orbital covering the whole cell costs one grid sweep, and 17 grid
sweeps is nothing beside 8778 boxes.  **Locality is a BONUS (compact boxes), not a precondition** — I had
elevated it to a gate, which was wrong.

**THE EIGEN/CHOLESKY SIDE-BY-SIDE (both now reported by `GPW_DM_RANK=1`).**  Eigen modes surviving
\f$\lambda>\text{tol}\cdot\lambda_{\max}\f$:

| tol | 1e-4 | 1e-6 | 1e-8 | 1e-10 | 1e-12 | 1e-14 |
|---|---|---|---|---|---|---|
| modes | 14 | 15 | **17** | **17** | **17** | 19 |

**The plateau at 17 across four decades IS the gap** — modes 18–19 appear only at 1e-14, i.e. roundoff.
Pivoted Cholesky at LAPACK's default floor gives 19, so eigen is marginally leaner (it is the minimal-rank
factorisation).  Locality of the two factors:

| factor | columns | IPR mean | IPR range |
|---|---|---|---|
| eigen (natural orbitals) | 17 | 6.8 | 1.9–12.2 |
| pivoted Cholesky | 19 | **3.4** | 1.4–5.7 |

Cholesky is ~2× more localized — matching the literature — but **both are tiny against n=118**, so the
natural orbitals are NOT the delocalised canonical picture: MnO's occupied manifold is atomic-like
(Mn 3s3p3d, O 2s2p).  **★★ AND THE PIVOTS ARE THE CHOLESKY ANALOGUE OF λ (user, 2026-08-21) — SO CHOLESKY IS SELF-SUFFICIENT.**
Split \f$U=D_p\bar U\f$ with \f$D_p=\mathrm{diag}(\text{pivots})\f$ and \f$\bar U\f$ unit-upper-triangular:
\f$D=(P\bar U^H)D_p^2(\bar U P^T)\f$, so \f$\rho=\sum_k d_k^2|\psi_k|^2\f$ — the same per-mode-density form,
with \f$d_k^2\f$ in λ's place.  And `pstrf` pivots on the largest remaining diagonal residual, so the
sequence is **MONOTONE BY CONSTRUCTION**: it IS the rank-revealing criterion.  MEASURED:

| | λ | pivot² |
|---|---|---|
| **kT=0** | 12.32 … 0.371, **13, hard stop** | 4.44 … 0.143, **13, hard stop** |
| **kT=5e-3** | 13.24 … 0.367 │ 2.6e-5, 1.2e-5, 1.2e-5 | 4.58 … 0.143 │ 8.0e-6, 7.4e-6, 7.1e-6 │ **2.8e-13, 2.1e-13** |

Same 13 at kT=0, same three-mode THERMAL tail at kT=5e-3 — **the two factorisations agree on where the
physics stops.**  (\f$d_k^2\neq\lambda_k\f$ — residual self-overlaps, not eigenvalues — but they order the
same content and terminate together, which is the property that matters.)
**★ And the pivots resolve ONE TIER MORE, which retires an earlier mis-reading:** they show 14 physical +
3 thermal + **2 at ~2e-13, i.e. ROUNDOFF**.  That is the whole of the "Cholesky rank 19 vs eigen 17"
discrepancy reported above — 19 = 14+3+2, with LAPACK's default floor admitting two roundoff modes.  **The
factorisations never disagreed about quality; I was comparing different TOLERANCES.**
**CONSEQUENCE: the eigendecomposition is not needed to see the spectrum.**  The pivot sequence carries the
same truncation information and falls out of the factorisation already being done, at O(n³/3) against
eigen's O(n³) with a larger constant.  Cholesky delivers factor + spectrum + rank in ONE call.
**KEEP BOTH ON THE TABLE (user, 2026-08-21).**  **So the choice is a mild trade, not a fork:** the ρ GEMM is O(npts·n·r) and wants the
smallest r (eigen, 17 vs 19); collocation wants compact boxes (Cholesky, IPR 3.4 vs 6.8).  Both work.
⚠ The user's narrowed pin still applies: the historical objection to trimmed eigen was INVERSION-driven and
we never invert here — but if D ever leaves the PSD cone, Cholesky is inapplicable and eigen is the route.

**Support is what sets a collocation box, not weight** — and since \f$\chi_i\sim e^{-\alpha r^2}\f$, a
\f$10^{-t}\f$ prefactor shrinks its reach only LOGARITHMICALLY.  So on MnO the orbital boxes are large
(~3 of 4 atoms), and the naive "Cholesky orbitals are localized ⇒ compact boxes ⇒ singles win" hope is NOT
supported *on this cell*.  ⚠ **But see (3): the literature's sparsity claim is ASYMPTOTIC, and a 4-atom
cell cannot exhibit it** — every orbital necessarily touches most of a 4-atom cell.  So this measurement is
**negative for MnO and INCONCLUSIVE for the sizes the idea actually targets.**  Re-measure on a supercell
before concluding.

**(2) DOES r STAY SMALL ON A METAL?  NOT MEASURED — and the attempt found something else.**
`DM_RhoAtPoints` never fires on the Si or Al tests, including the Becke-Al one: those runs do not reach the
DM route at all.  **The DM-ρ route is exercised on a NARROW set of configurations** (among those tried, only
the MnO recipe), which is itself an argument for the Vxc repair — that change is what would widen it.  The
metal-rank question stands open and needs a Becke-XC + DM-backed Al configuration to answer.

### Q2 — is there any literature on this?  ✅ **YES — it is established, and it has a name.**

**"Cholesky-decomposed density" (CDD) / "Cholesky orbitals".**  A pivoted Cholesky decomposition of the AO
density matrix is a known technique that *"preserves sparsity while reducing rank, with the rank at most
equal to the number of active occupied or virtual orbitals"*, and it *"can be considered as generating
localized molecular orbitals"* — obtained directly from the density matrix, **non-iteratively, with no
initial orbitals and no optimization**, and it is numerically stable and can be made linear-scaling for
matrices with a linear-scaling number of non-zeros.  Used in CDD-MP2 (including a relativistic variant),
DMRG-NEVPT2, and Edmiston–Ruedenberg localization over Cholesky-decomposed integrals.
**What this buys us, concretely:** the rank bound matches our measurement (r=14–19 against 13 occupied per
spin); the locality claim is REAL but ASYMPTOTIC, which reframes (1) above as a small-cell artifact rather
than a refutation; and "non-iterative, no initial guess" means the factor is safe to take per density
serial without a convergence story.
⚠ Distinguish from **Cholesky decomposition of the ERI tensor** (Beebe–Linderberg), which is the far
better-known use of the same word and a different object.
Sources: [Cholesky decomposition techniques in electronic structure theory (chapter)](https://www.diva-portal.org/smash/get/diva2:396223/FULLTEXT01.pdf) ·
[Relativistic Cholesky-decomposed density matrix MP2](https://www.sciencedirect.com/science/article/abs/pii/S0301010418311388) ·
[Multireference PT with Cholesky decomposition for DMRG](https://pubs.acs.org/doi/10.1021/acs.jctc.6b00778) ·
[Occupied and virtual Edmiston–Ruedenberg orbitals using Cholesky-decomposed integrals](https://pubs.acs.org/doi/10.1021/acs.jctc.2c00261)

**Superseded note** — what follows was the pre-search guess, kept only because its distinction still holds:
- I found two review articels in /home/janr/Documents/Reprints1/Qchem/Cholesky .  They are mostly focused on two different usages of CD: 1) Factoring ERI integrals, 2) RI Density fitting, using CD to find unbiased fit basis sets.

- **"Cholesky orbitals" / Cholesky decomposition of the DENSITY MATRIX** (Aquilante, Koch, Pedersen,
  Sánchez de Merás and co-workers, in the Cholesky-techniques line of work).  Pivoted Cholesky of D is, I
  believe, an established way to generate **localized occupied orbitals** — which matters far beyond
  citation etiquette; see Q3.  ⚠ Distinguish from **Cholesky decomposition of the ERI tensor** (Beebe &
  Linderberg), which is a different and much better-known use of the same word.
- **Plane-wave / PAW codes evaluate ρ from ORBITALS by construction** (\f$\rho=\sum_m f_m|\psi_m|^2\f$).
  So the factored form is the NORM outside the Gaussian world; what is unusual is the local-orbital
  convention of going through D.  The interesting literature is therefore on the CROSSOVER, not the idea.
- **Density-matrix vs orbital grid evaluation in DFT implementations**, and the screening-driven choice
  behind CP2K/Quickstep's pair collocation.
- Adjacent: density-matrix **purification** and linear-scaling methods, where PSD-ness and idempotency of D
  are load-bearing in the same way they are here.
- **Action:** a literature check should precede building — if Cholesky orbitals are standard, their known
  properties (locality, stability, behaviour under near-degeneracy) are free knowledge, and the failure
  modes are already documented by someone else.

### Q3 — can this ultimately kill the RAM-hungry PAIR STREAMS?

**Separate the two halves; they have different answers.**

- **RAM: plausibly YES, and for a structural reason.**  The pair streams' footprint is
  \f$\sum_{\rm pairs}({\rm box\ points})\f$ — a sum over variable-size boxes, which is why it reached
  **5.78 GB** and needed round 3 (→3.70) and the T3 fold (→1.35).  The factored route's footprint is a
  DENSE TABLE of known size: O(grid · n), or with the contraction applied first only **O(grid · r)** —
  the cached **Φ table, O(grid·n)** — 64000 × 118 doubles ≈ **60 MB** (the ORBITALS cannot be cached: D moves
  every iteration, so only Φ is geometry-fixed — user, 2026-08-21).  Bounded and predictable versus
  unbounded and data-dependent.
- **CPU: GENUINELY OPEN — this is the crossover, and my first framing of it was wrong twice.**
  - "118 singles vs 8778 pairs" was apples-to-oranges: 8778 counts **(i,j,R)** terms INCLUDING lattice
    offsets, 118 counts bare functions.  Singles enumerate **(i,R)**, so with ~133 images per function the
    comparable figure could be ~15k — possibly MORE objects than pairs.
  - **Pairs screen far harder**, and this is the real counterweight: the Gaussian product theorem gives
    \f$e^{-\mu|R_{ij}|^2}\f$ decay in the SEPARATION (\f$\mu=\alpha_i\alpha_j/(\alpha_i+\alpha_j)\f$),
    while a single is screened only by its own reach.  That is why CP2K collocates pairs.
  - "Singles are more diffuse" is NOT a real objection (user): the most diffuse pair is diffuse×itself at
    \f$2\alpha_{\min}\f$ — factor 2 in exponent, \f$\sqrt2\f$ in radius, ~2.8 in box volume — and cost is
    dominated by exactly those diffuse pairs.
  - **★ THE LOCALIZATION QUESTION MAY DECIDE IT.**  If pivoted Cholesky really does yield LOCALIZED
    orbitals (Q2), then \f$\psi_m\f$ has a COMPACT box rather than a whole-cell one, and the r=14–17
    orbitals collocate cheaply — which would make the singles route win outright.  If instead the columns
    of L are delocalized, each \f$\psi_m\f$ costs the whole grid and the route is much weaker.
    **MEASURE THE SPATIAL EXTENT OF L's COLUMNS.**  This is the single most decisive unmeasured quantity in
    the whole idea, and it is cheap: the factor already exists in `LowRankFactor`.
- **REQUIRED CENSUS before building:** count \f$(i,R)\f$ singles against \f$(i,j,R)\f$ pairs **at the same
  ε**, each weighted by BOX VOLUME.  Object counts alone decide nothing.  The pair side already
  self-reports (`[fold] collocation streams … 8778 → 1909`); the singles side needs the same census.
- ⚠⚠ **THE STRONGEST STRUCTURAL OBJECTION, and it is not locality: THE FACTORED FORM LOSES THE MULTIGRID.**
  GPW's multi-level collocation works because \f$\rho=\sum_{ij}D_{ij}\chi_i\chi_j\f$ is **LINEAR in the pair
  products** — each pair is assigned to the coarsest level that resolves \f$\alpha_i+\alpha_j\f$ and the
  per-level densities are simply SUMMED (`CollocateDensity` returns one grid density per level).
  \f$\rho=\sum_m|\Psi_m|^2\f$ is **QUADRATIC in Ψ**, so Ψ must be assembled at ONE resolution before it can
  be squared — level contributions cannot be summed after squaring.  The factored route therefore forces
  every function onto the FINE grid, including the DIFFUSE ones that currently live on cheap coarse levels,
  which is exactly the saving the ladder exists to capture.  (Ironically a single→level assignment would
  otherwise be SIMPLER than the pair one — one exponent \f$\alpha_i\f$, not a sum.)
  **Measurable before building:** read the ladder's level occupancies and per-level point counts, and price
  collocating everything on the fine grid against today's distribution.
  Workarounds if it bites: assemble Ψ per level and INTERPOLATE to fine before squaring (interpolation
  error on a SMOOTH Ψ, unlike on ρ), or keep the factored route for the XC MESH only — atom-centred, no
  ladder — and leave uniform-grid collocation on pairs.
- ⚠ **System-dependent, and it inverts.**  n=118 here; the battery north-star's supercells grow n and make
  pair screening bite harder.  Both routes may deserve to survive rather than one replacing the other.

### ★★ THE SPECTRUM ITSELF: r IS 14, λ DOES NOT PREDICT LOCALITY, AND {|φ_m|²} IS A FIT BASIS

**(a) ⛔ RETRACTED — "THE REAL RANK IS 14 AND THE REST IS NUMERICAL DROSS" WAS WRONG.  THE SMALL MODES ARE
THERMAL OCCUPATION.**  The user's challenge — *"I will buy into this iff you can confirm kT=0.000000000
for this run"* — was decisive.  **kT = 0.005 Ha on that run, not zero**, and λ≈2.6e-5 under Fermi smearing
means \f$(\varepsilon-\mu)/kT=\ln(1/f)\approx10.6\f$, i.e. **≈0.053 Ha ≈ 1.4 eV above μ** — an entirely
plausible thermally-occupied conduction state.  The experiment settles it:

| | rank at tol 1e-4 … 1e-14 | smallest λ | tail |
|---|---|---|---|
| **kT = 5e-3** | 14 → 15 → 17 → 17 → 17 → 19 | 0.367 | **3 modes at ~2.6e-5** |
| **kT = 0** | **13, 13, 13, 13, 13, 13** | 0.371 | **NONE** |

**Rank 13 at EVERY tolerance across ten decades, and the 2.6e-5 modes vanish with the smearing.**  13 is
exactly the electron-pair count (26/2) — a free correctness check on the whole census.
**Consequences, and they matter:**
- **Truncating at "r=14" would DROP PHYSICAL DENSITY** (~3×2.6e-5 ≈ 8e-5 e of thermal tail) while calling
  the result exact.  The IMPLEMENTATION is safe — LAPACK's default floor keeps rank 19 — but the
  interpretation was the error, and it is the kind that would have shipped as an "exact" optimisation.
- **THE RANK IS kT-DEPENDENT**: 13 at kT=0, 17–19 at kT=5e-3.  Exactness at kT=0 is a hard spectral gap at
  machine precision; at kT>0 "exact" becomes "accurate to whatever thermal tail you choose to drop".
- **★ This ANSWERS open question 4 (does r stay small for metals / larger smearing?): NO, not in general.**
  A fatter smearing tail means more thermally occupied states and a higher numerical rank, so the n/r win
  SHRINKS with kT.  MnO at kT=5e-3 *with a gap* is a favourable case; a metal at larger kT is the adverse
  one, and that is now an argued expectation rather than an open guess.
- Technical caveat: λ(D) are not occupation numbers, since the AO basis is non-orthogonal (occupations are
  eigenvalues of \f$S^{1/2}DS^{1/2}\f$).  The RANK is basis-independent, and rank=13=N/2 at kT=0 confirms
  the reading.
**Lesson: a four-decade gap in a spectrum is not automatically a numerical one — ask what physics could
put something there before calling it noise.**

**(b) DOES LOCALITY CORRELATE WITH λ (user's idea: if so, assign grid levels from λ)?  MEASURED: NO.**
Descending in λ, spin ↑:

| λ | 13.24 | 1.33 | 1.09 | 0.567 | 0.528×2 | 0.463×2 | 0.404 | 0.389×2 | 0.377×2 | 0.367 |
|---|---|---|---|---|---|---|---|---|---|---|
| IPR | 4.7 | **2.0** | 6.6 | **10.9** | 5.1 | 6.5 | **11.6** | 8.5 | 8.3 | 6.3 |

The largest λ is moderately localized, the SECOND largest is the most localized of all, and the most
DELOCALIZED modes sit mid-spectrum.  So a level assignment cannot be read off λ *via locality*.
(Repeated λ are symmetry-degenerate pairs; their IPRs match exactly, as they must — a free correctness
check on the whole census.)

**CAN IPR ASSIGN GRID LEVELS? (user, 2026-08-21)  It is FREE to compute but it is the WRONG QUANTITY —
and the right one is equally free.**
- **Cost: negligible.**  IPR is O(n·r) ≈ 118×14 ≈ 1650 ops against the O(n³) ≈ 1.6e6 eigendecomposition
  that produced U, i.e. 0.1% of the factorisation and invisible beside the ~1e8-MAC GEMM.
- **But a ladder level resolves BANDWIDTH, not EXTENT.**  IPR counts how many basis functions carry a
  mode (box size); a level is chosen by how SHARP the mode is (how high in G).  Independent properties: a
  mode can be localized-and-smooth (one atom, diffuse ⇒ coarse level, small box) or delocalized-and-sharp
  (tight functions on several atoms ⇒ fine level, big box).
- **The right quantity is the exact analogue of the EXISTING pair rule.**  `CollocateDensity` already puts
  each pair on the coarsest level resolving \f$\alpha_i+\alpha_j\f$.  The mode analogue is
  \f[ \alpha_{\rm eff}(m)=\max\{\alpha_i:\ |U_{im}|\ \text{above a magnitude screen}\} \f]
  giving \f$|\phi_m|^2\f$ a bandwidth \f$2\alpha_{\rm eff}(m)\f$ — also O(n) per mode.  The screen on
  \f$|U_{im}|\f$ is the project's standard ε discipline; without it a trace admixture of one tight function
  would drag a whole mode onto the fine grid.
- **The blocker is placement, not cost:** per-function \f$\alpha_i\f$ does not reach the `IrrepCD` seam
  (only whole-basis `MaxExponent`/`MinExponent` do).  It DOES live inside the molecular seam, where the
  pair→level rule already runs with primitives encapsulated — **so the mode→level assignment belongs
  there, handed the FACTOR instead of D.**  Same place, same rule, different object.
- **★ BUT ASK WHETHER THE LADDER IS NEEDED AT ALL HERE.**  The multigrid exists to keep ~8778 diffuse
  PAIRS off the fine grid; its whole economic case is object count.  With **14 modes**, putting EVERY mode
  on the fine grid costs 14 sweeps — which is the GEMM already being done.  **The factored route may
  simply not need a ladder**, which would delete the assignment problem rather than solve it.  Measure
  before building either: price 14 fine-grid modes against today's per-level pair distribution.

**(c) ⛔ RETRACTION — THE MULTIGRID IS NOT LOST.**  Keeping λ SEPARATE (user) makes it obvious:
\f$\rho=\sum_m\lambda_m|\phi_m|^2\f$ is a **SUM OF PER-MODE DENSITIES**, so each mode can be collocated on
its own ladder level and the level densities summed — EXACTLY as pairs are.  What the previous entry
actually showed was narrower: computing \f$\Psi=\Phi L\f$ as ONE DENSE GEMM forces one resolution.  That is
an implementation choice, not a constraint.  **The real trade is single-GEMM efficiency vs per-mode
multigrid**, and both are available.

**(d) ★★ AND THAT IS A FIT BASIS EXPANSION (user, 2026-08-21).**  \f$\rho=\sum_m\lambda_m|\phi_m(r)|^2\f$ is
an expansion of the density in the basis \f$\{|\phi_m(r)|^2\}\f$ with coefficients \f$\lambda_m\f$ — and it
is **EXACT, MINIMAL (14 functions), ADAPTIVE, and FREE**: no fitting equation, no auxiliary basis to
choose, no Dunlap constraint, no J-matrix solve.  The eigensolver hands over the coefficients.
This slots into the project's existing framing — `DeltaFittedVxc` is a *fitted* Vxc whose fit basis is
DELTA FUNCTIONS, i.e. quadrature is the δ-basis special case of fitting (user).  So:

| fit basis | size | character |
|---|---|---|
| δ functions | n_pts | quadrature — exact, expensive |
| auxiliary Gaussians (`VALENCE`, `A1_exch`) | hundreds | approximate, cheap ANALYTIC integrals |
| **\f$\{|\phi_m|^2\}\f$** | **14** | **exact, adaptive, GRID-evaluable** |

**Where it fits and where it does not:** \f$|\phi_m|^2\f$ is cheap to EVALUATE (GEMM then square) but its
analytic integrals against other Gaussians are not obviously cheap, so it serves GRID/QUADRATURE consumers
(XC, collocation) and not analytic-integral ones (the 3-centre Coulomb path).  And it is PER-ITERATION
adaptive, which suits grids and is awkward for anything cached across iterations.
**Open:** is there a consumer that wants an exact 14-function density expansion?  The Hartree route needs
one FFT of the total ρ, not 14 Poisson solves, so it is not obviously that one — the honest first use is
the density itself on the grids, which is where this whole section started.

### ★★★ THE OCCUPIED SUBSPACE FREEZES AFTER ~3 ITERATIONS — and that is the ONE lever the rank cannot give

**Question (user, 2026-08-21):** does the 13-mode subspace CHANGE across iterations, or do the iterations
merely ROTATE modes within a fixed subspace — and if the latter, can a projector be worked out before the
SCF starts?  **MEASURED** (kT=0, r=13, `GPW_DM_RANK=1` now reports it), as
\f$\lVert U_{prev}^\dagger U_{now}\rVert_F^2/r=\sum_k\cos^2\theta_k/r\f$:

| iteration | overlap | **dimensions that MOVED** |
|---|---|---|
| 1 → 2 | 0.613 | **5.03** |
| 2 → 3 | 0.874 | 1.63 |
| 3 → 4 | 0.9994 | **0.008** |
| 4 → 5 | 0.99995 | 0.0007 |
| 7 → 8 | 0.999999 | **0.00001** |

**BOTH readings are right, separated in time.**  It must move — the SCF is solving FOR it, and *"if it
didn't move you'd be converged at iteration 1"* (user: "right!!!").  But it moves 5 of 13 dimensions in the
FIRST step, 1.6 in the second, and from **iteration 3 onward is FROZEN to 4–6 digits**: the later
iterations really are just rotating modes inside a fixed subspace.
**So a projector cannot be formed BEFORE the SCF (5 dimensions would be wrong) — but it can after ~3
iterations**, and the union over the whole run is **~20 dimensions, not 118**.

**★ WHY THIS IS THE INTERESTING ONE: it hits the cost the RANK CANNOT.**  \f$h_{ij}=(\Phi^\dagger W\Phi)_{ij}\f$
is O(npts·n²) with NO D in it, so the low-rank factor does nothing for the H_xc GEMM (recorded above as a
correction to an over-claim).  Projecting Φ to (npts × d) ONCE makes the per-iteration H build
**O(npts·d²)** — with d≈20 against n=118 that is **~35×** on an otherwise immovable bucket.
**Safety:** freezing CONSTRAINS the variational problem, but the table bounds that constraint at <1e-4
after iteration 3, and one FULL-SPACE step at the end verifies stationarity.  That is textbook subspace
iteration with a convergence check, not a new gamble.  The safer variant is Davidson-style: keep the
ACCUMULATED subspace (~20 dims) and expand only when the residual demands.
**Bonus finding: the SEED already gets 8 of 13 dimensions right** (only 5 move).  So seed quality directly
sets how much must move — a sharper seed shrinks the transient AND the accumulated subspace.
**Open:** all of this is kT=0.  With smearing the "occupied subspace" is the thermally-occupied one and its
dimension is kT-dependent (13 → 17–19 here), so the freeze-out should be re-measured at kT>0 before use.

### DESIGN RULING — the factor is a POLICY, not a CD type (LSP, user 2026-08-21)

**Proposal considered:** two `FittedCD` flavours (eigen, Cholesky) taking D + Φ + a λ cutoff instead of a
fit basis, with a `CompositeFittedCD` holding a map of per-irrep ones.  **Rejected on three grounds, the
last of which is decisive.**

1. **`CompositeCD` ALREADY IS that composite** — `CompositeCD.C` is *"any array of Irrep DM_CDs"*, holding
   `variant<unique_ptr<tDM_CD<double>>, unique_ptr<tDM_CD<dcmplx>>>` children with `DM_RhoAtPoints`
   documented as "sum over irrep blocks".  A second composite would duplicate it, re-solving the mixed
   real/complex child problem the variant already solves.  **The multi-irrep question dissolves: the
   factorisation is per-block and block aggregation is done.**
2. **`FittedCD` is the wrong base.**  It exists for the ANALYTIC COULOMB path (`GetRepulsion`,
   `GetSelfRepulsion`) — which \f$\{|\phi_m|^2\}\f$ cannot serve cheaply — and its `DoFit` implies a fit
   with a residual and a quality measure (the standing pin: *fit quality = grid-convergence of ρ*).  The
   factorisation is EXACT, not approximate; naming it a fit imports a contract it need not honour.
3. **★ LSP KILLS IT (user's criterion, and it is the right one).**  `tDM_CD<T>` inherits
   `tMixableDensity<T>`, whose contract is `MixIn(other, c)` = *this = (1−c)·this + c·that*.
   **Low rank is NOT CLOSED under mixing:** \f$(1-c)L_1L_1^\dagger+cL_2L_2^\dagger\f$ has no rank-r factor;
   the honest result is \f$[\sqrt{1-c}L_1,\sqrt{c}L_2]\f$ of rank \f$r_1+r_2\f$, so the rank DOUBLES every
   iteration (13→26→52…).  Bounding it means reconstruct + refactor O(n³) per mix — storing D with extra
   steps.  And the contractions get WORSE:

   | `tDM_CD` operation | with D | with L |
   |---|---|---|
   | `DM_RhoAtPoints` | O(npts·n²) | **O(npts·n·r)** ✅ |
   | `DM_Contract` = Tr(DV) | O(n²) | Tr(L†VL) = O(n²r) ❌ |
   | `DM_ContractBlocks` | O(n²)/block | O(n²r) ❌ |
   | `MixIn` | O(n²) | **unbounded rank** ❌ |

   The factored form wins **one of four** operations.  A `FactoredCD` would satisfy every SIGNATURE while
   breaking the BEHAVIOURAL contract its callers rely on — and the failure is SILENT (rank creep, not a
   crash), which is the worst kind.

**RULING: D stays the truth; the factor is a DERIVED, CACHED representation used only by the operation that
benefits.**  ⚠ My first shape for that — a policy MEMBER on `IrrepCD_Core` — is SUPERSEDED by the user's
(better) one: it would have put factorisation state on EVERY density including those that never factor.

**★ THE DESIGN (user, 2026-08-21): a derived leaf that overrides ONLY `DM_RhoAtPoints`, selected by the
existing Factory.**  `DM_RhoAtPoints` is already `virtual` on `IrrepCD_Core<T>`, and D + `MixIn` +
the contractions are all INHERITED unchanged — so LSP holds by construction: same contract, same values
(exact to roundoff), different cost.
Note the existing mixins (`IrrepCD_Fourier<PeriodicIrrepCD<T>>`, `IrrepHF_PairBase<T,Leaf>`) ADD faces;
they do not override core virtuals, and overriding across a virtual-inheritance diamond invites dominance
ambiguity — so LINEAR derivation is the right mechanism here.
**One template, so the two orthogonal axes (leaf × factorisation) COMPOSE instead of multiplying:**
```cpp
//! Same density, same value, cheaper route: D stays the truth on the Leaf; only rho(r) is factored.
template <class Leaf, class Fact> class FactoredRho : public Leaf
{
public:
    using Leaf::Leaf;                                  // as the leaves already do from the core
    virtual rvec_t DM_RhoAtPoints(const rvec3vec_t&,
                                  const std::map<Irrep,mat_t<T>>&) const override;
private:
    mutable mat_t<T> itsL;  mutable size_t itsRank=0;
    mutable size_t   itsFactorVersion=size_t(-1);      // keyed on IrrepCD_Core::itsVersion
};
```
- **The invalidation key already exists**: `IrrepCD_Core` holds `size_t itsVersion //!< TRANSIENT freshness
  serial`.  `DM_RhoAtPoints` is `const`, so the memo is `mutable` — a normal derived cache.
- **The PSD fallback becomes inheritance, not a branch**: `LowRankFactor` returning false ⇒ call
  `Leaf::DM_RhoAtPoints(...)`.  That IS "D stays the truth", expressed structurally.
- **`CompositeCD` needs nothing**: these are `tDM_CD<T>`, so a factored block can sit beside an unfactored
  one in the same composite — useful, since two spin channels or a real TRIM block may reasonably differ.

**★★ CLIENTS GET INSTANCES FROM THE FACTORY; THE CONCRETE CLASSES STAY INTERNAL (user).**  The seam is
already there and already returns the ABSTRACT face:
`template <class T> tDM_CD<T>* IrrepCD_Factory(const hmat_t<T>&, const tobs_t<T>*, Irrep)`, with the leaves
in `qchem.ChargeDensity.Imp.IrrepCD` (off the public surface).  So this EXTENDS an existing seam:
```cpp
enum class RhoRoute { Direct, PivotedCholesky, EigenTrim };     // the ONLY thing that crosses the boundary
template <class T> tDM_CD<T>* IrrepCD_Factory(const hmat_t<T>&, const tobs_t<T>*, Irrep,
                                              RhoRoute = RhoRoute::Direct);
```
**ONE enum, not two — the LEAF axis must NOT become one.**  The factory already picks the leaf by
cross-casting the BASIS (`Orbital_DFT_IBS<T,dcmplx>` ⇒ periodic, else the finite leaf, with a
`if constexpr` guard that a finite COMPLEX density cannot exist).  Deriving it from the argument is
strictly better than an enum: a caller cannot request a leaf inconsistent with the basis it passed.
`using Leaf::Leaf` makes all routes constructible identically, so the switch is uniform; the default
`Direct` leaves every existing call site untouched; and a fourth route later is one enum value plus one
case, with no client change and nothing new exported.
**FILED, NOT ACTED ON — the ISP observation underneath:** `tDM_CD` BUNDLES four capabilities and the
factored form serves one.  If a factored TYPE is ever wanted as a first-class density (rather than a leaf
override), the prerequisite is splitting `tDM_CD` into narrower faces (grid-evaluable vs mixable vs
contractable).  Same shape as the `LatticeSum1E` ISP item; it belongs in that deferred session, not forced
now by a performance idea.

### Further questions (mine)

1. **Does the INTEGRATE-BACK factor too?**  Partly answered: \f$h_{ij}=\langle\chi_i|V|\chi_j\rangle
   =(\Phi^\dagger W\Phi)_{ij}\f$ is already a Φ-shaped GEMM in `XC_SinglesQuadrature` — so the XC side is ALREADY
   pair-free.  It is the HARTREE/density-collocation path on the uniform grid that still pays pairs, so Q3
   is really about that path alone.
2. **The rank does NOT help \f$h\f$.**  \f$\Phi^\dagger W\Phi\f$ is O(npts·n²) with no D in it, so the
   H_xc GEMM stays quadratic in n however small r is.  Any "everything becomes O(n·r)" claim is wrong.
3. **Multi-k:** each k-block has its own \f$D^k\f$ and factors independently, so the rank is per-block and
   the scheme should carry over — worth confirming that the per-block ranks stay small at MNO_KMESH=2.
4. **Does r stay small where it matters?**  r=14–17 was measured under Fermi smearing at kT=5e-3 with a
   gap.  METALS (Al) and larger smearing put more states in the fractional tail.  Re-run `GPW_DM_RANK=1`
   on Al and on a metallic recipe before assuming the ratio generalises.
5. **Does the factorisation cost stay negligible at scale?**  O(n³) per density serial is nothing at
   n=118 (~5e5 flops vs a ~1e8 GEMM), but it grows faster than the GEMM it saves.  Find the n where they
   cross; it is probably far away, but it should be a known number rather than an assumption.
6. **Interaction with the T3 stream fold.**  The fold's 4.60× is a reduction on PAIRS.  Orbitals are not
   individually symmetry-adapted, so an equivalent fold on singles may not exist — in which case the
   honest comparison is singles against **1909 folded** pairs, not 8778 raw.
