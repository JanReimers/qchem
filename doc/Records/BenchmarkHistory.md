# Benchmark HISTORY — how the qchem ↔ CP2K rows got to where they are

**This file is the RECORD, not the instrument.**  `doc/Benchmark.md` carries the live protocol, the live
rows and the standing rules; everything here is the reasoning, the wrong turns and the superseded numbers
that produced them.  Split out 2026-09-04, when Benchmark.md hit 809 lines and stopped being readable
(user: *"there is a lot of history in there, about wrong turns or parity mistakes"*).

⚠ **NOTHING HERE IS A CURRENT NUMBER.**  Every table below has been superseded at least once.  Quote
`doc/Benchmark.md`; come here only to find out WHY a row reads the way it does, or to avoid re-walking a
path that has already been walked.

---

## 1. The 2026-09-04 retraction — a re-take that measured a different system

⛔⛔ **A 2026-09-04 ATTEMPT TO RE-TAKE THESE ROWS WAS RETRACTED — READ THIS BEFORE QUOTING ANY 09-04 NUMBER.**
The re-take ran `MNO_SKIP_FM=1 GPW_REPORT=1` and NOTHING ELSE, i.e. it silently used the DEFAULT basis
(`VALENCE_LOWQ_SR`, Cartesian) and a SINGLE SCF stage — where every row in this table is
`GPW_SPHERICAL=1 GPW_BASIS_SPAN=va` (VA, 118 functions, spherical) with the two-stage
`MNO_ANNEAL="5e-3,0" MNO_ACC="Ladder,GDM" MNO_MOM=0 …` recipe printed in doc/Benchmark.md §4.
⇒ It was a DIFFERENT SYSTEM, so "1.35× / 1.54× / 1.07× / 1.00× vs CP2K" and the claim that these rows were
"stale by ~1.45×" are ALL WITHDRAWN.  CP2K's 8.5 s/step is a VA-deck number and only a VA run may be
divided by it.
★ **THE TELL, recorded so the next person catches it in one minute instead of an afternoon**: this probe is
documented as *"~6 minutes"* and the retracted runs finished in **90 seconds**.  A 4× discrepancy against
the written cost of the SAME probe is a CONFIGURATION difference until proven otherwise -- it is not your
speedup.  ⇒ **A row is comparable only if it names its recipe.  Copy the command, do not reconstruct it.**


**The corrected sequel.**  Re-run with the printed command, every row reproduced its banked \f$E_{tot}\f$
and iteration count exactly, and the honest deltas against 08-28 were −24.3% (ALL DEFAULTS), −30.9%
(`QCHEM_BECKE_XC=0`) and −7.3% (the parity probe).  Those are in Benchmark.md's live ladder.

### 1b. The A/B deltas from that session, which DO stand

⚠ **THE 09-04 A/B DELTAS ARE VALID; THEIR ABSOLUTE COLUMN IS NOT.**  Every pair below ran the SAME
(non-VA, single-stage) configuration on both sides in one session on one binary, so the RATIOS stand and
are what the changes are worth.  There is NO CP2K column, deliberately: that configuration has no CP2K
counterpart (see the retraction above).  ⇒ Read these as "what the change did", never as a scoreboard.

| MnO, `CP2K_COMPAT=1`, serial, NMAX=10 (⚠ default basis, 1 stage — NOT a table row) | collocate/call | gather/call | wall |
|---|---|---|---|
| D-aware screen, before the XC fix | 1.266 s | 1.254 s | 1:54.6 |
| geometry-only screen, before the XC fix | 1.520 s | 1.421 s | 2:11.2 |
| D-aware screen, AFTER the one-gather XC term | 1.256 s | 1.241 s | **1:25.1** |
| geometry-only screen, AFTER the one-gather XC term | 1.393 s | 1.321 s | **1:31.4** |

⇒ the D-aware screen is worth **+14.5% wall**; the one-gather XC term is worth **−30.6% wall** (gather
misses 65 → 43, its bucket 92.4 → 56.8 s) with \f$E_{tot}\f$ IDENTICAL to all 10 printed figures on both
arms.  **Whether either lands us near CP2K is UNMEASURED** — it needs the VA recipe.

Both arms: 65 gathers / 22 collocations, \f$E_{tot}\f$ agreeing to 2.3e-7 — identical trajectories, so the
per-call column is a clean A/B.  Peak RSS 110 MB either way.  Threads: `OMP_NUM_THREADS=1
GPW_OMP_THREADS=1`, BLAS pinned to 1.  ⚠ The `[ FAILED ]` on this probe is the `NMAX=10` cap, not a defect.

¹⁰ **`CP2K_COMPAT=1` NOW IMPLIES THE GEOMETRY-ONLY SCREENER** (2026-09-04, doc/OldPlans/ScreeningPlan.md §7).  CP2K
does not screen the collocation on the density, and this tree had been taking that deviation silently — so
the honest parity row is the second one, **1.54×**, and it is the first parity row that actually deserves
the name.  The D-aware arm above is the qchem-vs-qchem delta: our screen buys **+14.5% wall**.


### 1c. The superseded per-call table, and the footnote that was withdrawn with it

⁸ ⚠ **THE 93-ITERATION FIGURE IS PRE-FIX AND ITS STAGE MIX DIFFERS — do not read 29.4 → 19.3 as a 1.5×.**
The two runs make different numbers of collocations PER iteration (13.4 against 10.0) because a longer
stage 2 changes the GDM line-search mix, so per-ITERATION is not comparable across caps.  ⇒ **The
cap-independent measure is PER CALL**, straight off the ledger, and that is what the bin-1 pass moved:

| per call, MnO | before the 2026-08-28 per-line pass | after |
|---|---|---|
| collocate / gather, ALL DEFAULTS | 0.444 / 0.509 s | **0.402 / 0.419 s** |
| collocate / gather, `CP2K_COMPAT=1` | 2.03 / 2.24 s | 1.84 / 1.88 s ⁹ |

⁹ ⛔ **WITHDRAWN 2026-09-04.**  This footnote claimed the row above was "stale by ~1.45×"; the run behind
that claim used the DEFAULT basis and a single SCF stage, not the VA/annealed recipe these rows are taken
with, so it measured a different system and said nothing about this row.  The 08-28 numbers stand until
someone re-takes them WITH THE PRINTED COMMAND.  What the 09-04 session did establish, on its own
configuration and as A/B deltas only:


---

## 2. What `CP2K_COMPAT=1` cost us to discover — the 2026-08-26 reckoning

The distilled, current answer lives in Benchmark.md §2.  This is how each item was found.

### ★ THE `CP2K_COMPAT=1` ROW — WHAT THE DEFAULT ROW ABOVE IT WAS HIDING (2026-08-26)

**First, the reassuring half: the deviations are pure ACCELERATIONS, not physics.**  Turning all four off
moves the MnO total by **3e-8 Ha** (−61.40297621 → −61.40297618, agreeing to 10 s.f.) on a different
trajectory (10+14 iterations against 14+17).  That is the property `CP2K_COMPAT` most needed to demonstrate
about itself, and it is now measured rather than asserted.

**Then the uncomfortable half.**  Stripped of our own accelerations the standing against CP2K is
**2.5× CPU and 23× RAM**, not the 1.8×/6.1× the default row reports — and 1.79× per ITERATION (38.4 s
against 21.4 s), since the compat run happened to need fewer of them.  The default row is not wrong, but it
is qchem-with-accelerations against CP2K-plain, and it should never be quoted as an algorithm-to-algorithm
comparison.

**⚠ AND `CP2K_COMPAT=1` IS STILL NOT PARITY** — the switch covers four deviations and at least two more
matter (doc/OpenWork.md N5 carries the full list):
- ~~**The pair-stream CACHE**~~ ✅ **DELETED 2026-08-27** (`doc/OldPlans/CollocationRewritePlan.md` step 7), so this
  deviation no longer exists: like CP2K, qchem now re-evaluates the orbital pairs every iteration and keeps
  only a ~0.2–0.4 MB (shell pair, offset) TASK LIST.  The history below is kept because it is what turned
  the campaign toward making the on-the-fly evaluation fast instead of caching harder — and because its
  last row is the one that finally made the cache indefensible.  **What ACTUALLY closed it**: after the
  contraction kernel, the cache bought 2.91× on the two box-walk buckets and ~1.1–1.5× on a whole run,
  against 25× the RAM on the unfolded probe.  ⚠ Read the MnO row above before quoting this as a pure win.
  Original measurement, zeroing the budgets (`GPW_STREAM_BUDGET_PTS=0`, all 8778 pairs on-the-fly) on the
  MnO 3-iteration probe:

  | MnO, per SCF iteration | CPU/iter | peak RSS |
  |---|---|---|
  | CP2K (44 steps, 373 s CPU; its own log prints 8.4–8.5 s/step) | **8.5 s** | **217 MB** |
  | qchem, cache ON  | ~31 s (3.6×) | 4218 MB (19×) |
  | qchem, cache OFF — **as measured 2026-08-25, ESTIMATED** | ~853 s (100×) | 353 MB (1.6×) |
  | qchem, cache OFF — **MEASURED from the ledger, before the box-walk work** | **573 s (67×)** | **166 MB (0.8×)** |
  | qchem, cache OFF — **MEASURED, after it (2026-08-26)** | **260 s (31×)** | **166 MB (0.8×)** |
  | qchem, cache OFF + contracted kernel — **MEASURED 2026-08-27, and this is now the DEFAULT** | **35.5 s (4.2×)** | **155 MB (0.7×)** |

  ⇒ **The cache is not an advantage we hold over CP2K; it is a workaround for an on-the-fly box walk that
  was 67× off theirs, bought with 4 GB.**  Without it our RAM is BETTER than CP2K's (166 vs 217 MB), so
  the target was never "cache more cleverly" or "trade RAM for CPU" — it was **make the on-the-fly pair
  evaluation fast**.  That is now half done: **2.21× measured on this very probe** (doc/OpenWork.md, the
  box-walk section), which is the "even 2 or 3× would pay for itself" the user asked for.
  ⚠ **THE 853 FIGURE WAS AN INSTRUMENT ARTEFACT, NOT A MEASUREMENT.**  It was (2740 s CPU for 3
  iterations) minus an ESTIMATED ~182 s setup, because `RunMnO` drives `SolidCalculation`, which opens no
  report run — so the benchmark's most expensive row was the ONE campaign run with no timing ledger.
  `e8339cf2` gives the arm the same `GpwReport` bracket every other driver holds; the ledger's exclusive
  buckets then sum to the wall clock (1201.3 s of 1202.1 s) and nothing is subtracted.  **The estimate was
  1.5× pessimistic** — hence "67×", not "100×".  There is still no policy hook for the cache, only the raw
  budget knob.
  ⚠ Provenance for both measured rows: `MNO_SKIP_FM=1 GPW_MNO_NMAX=2 GPW_REPORT=1
  GPW_STREAM_BUDGET_PTS=0 GPW_STREAM_BUDGET_PTS_F32=0`, AFM arm, **symmetry FREE and NO fold active**
  (`[fold] collocation streams (T3 pairs): NONE`), `GPW_OMP_THREADS=1`, BLAS pinned to 1, measured 103%
  CPU — i.e. serial, unfolded, and the two binaries differed ONLY in the box-walk diff.  Trajectory
  identical both sides (iters, lastΔρ, m_stag, Eee, site moment); `Efinal` moves 2e-8 Ha.
- ~~**`imposeSymmetry` ITSELF.**~~  ✅ **WIRED 2026-08-26** — `CP2K_COMPAT=1` now implies
  `imposeSymmetry=0` (knob `QCHEM_IMPOSE_SYMMETRY`).

  ⛔ **AND THE RE-TAKEN ROW DOES NOT EXIST: AT TRUE PARITY OUR MnO RECIPE DOES NOT CONVERGE.**  Measured
  2026-08-26, the banked recipe under `CP2K_COMPAT=1` with the imposition vetoed: stage 1 caps at 80
  iterations at −60.431, stage 2 caps at 80 more at **−57.620** — 3.8 Ha short of the −61.40297618 the
  imposed compat run reaches in 24 iterations (2555 s CPU, 42m57s, 4.76 GB).  **So there is currently no
  honest MnO row to put in this table**, and that absence is the finding: the imposed star-average was
  buying CONVERGENCE, not accuracy and not the magnetic basin.
  ★ **The AFM ORDER SURVIVED the free run** (m_stag 0.66/0.59, integrated site moment 4.781 → 4.222 e), so
  the imposition was NOT what held the basin — which is the opposite of what was expected, and it moves
  the question from symmetry to the MIXER.  Note where CP2K stands on the same cell: 44 steps at 8.5 s
  each with Broyden α=0.2 / NBUFFER 8 / MAX_SCF 200, against our α=0.45 / PulayDepth 0 / 80.  ⇒ A fair
  MnO row needs a CP2K-like mixing recipe first; taking one before that would be comparing their converged
  answer against our iteration cap.
  Verified locally 2026-08-26: CP2K does NO symmetry work in these decks.
  The 1129-line `bench_MnO_AFM2_VA_cp2k.log` contains **zero** occurrences of "irrep", "symmetry" or
  "point group" — QuickStep keeps K and P as DBCSR sparse ATOM-BLOCK matrices over the full AO basis and
  diagonalizes the lot (blocking by atom-pair SPARSITY, not by irrep; no SALC blocking, no k-block
  splitting).  The one symmetry knob that exists is BZ-side and defaults OFF: our own Si 2×2×2 log reads
  `BRILLOUIN| K-Point point group symmetrization  OFF` and lists all 8 k-points.
  ⇒ Our imposed MnO row folds the BZ, star-averages ρ every iteration, uses the site-adapted invariant XC
  mesh (~2×) and folds the collocation streams (5.2× on pairs).  **The CP2K row does none of it.**  These
  rows compare a SYMMETRY-exploiting code against a SPARSITY-exploiting one on a small, high-symmetry
  cell — the regime that most favours us, and a 100-water box would invert it.  (CP2K's design centre is
  large disordered systems, where every one of {G}, {k}, {r} has a group of order 1, so the folding payoff
  is exactly 1× for what they build for.  That reading of WHY is inference; the WHAT above is from the logs.)

Δ(AFM−FM) on the VA span: **qchem +38.61 mHa, CP2K +1.46 mHa** — both order FM first, and the
**configuration-SELECTIVE part of the offset is −37.15 mHa** (`OpenWork` Step 5).  Every one of these four
energies reproduces its banked value (runs 61/62 and the CP2K VA pair) to the digits those were recorded at.

¹ **SOLVED 2026-08-19 — and it was a SCREEN, not a phase.**  This test had rotted to −3.7351 while DISABLED
(it is the suite's ONLY fractional-k SCF coverage; every other k is TRIM, where the defect is structurally
invisible).  It is now ENABLED at ~14 s and reads −7.868473428, converged in 16 iterations at Δρ 1.0e-9.

**The bug** (`PG_Cart_MnD/Evaluator.C`, the D-aware integrate-back screen): the term was dropped when
`|Re(D_ij · conj(phase))| · maxv < eps` — a REAL PART used as if it were a magnitude.  `Re[D e^{-ikR}]` is
the right coefficient on the COLLOCATION side, where it multiplies a real pair product and a zero means a
genuinely zero contribution to ρ; the integrate-back's term is `phase·b`, whose size is `|b|` however the
phase is oriented.  **At a quarter-integer k, \f$e^{2\pi ikn}=i^n\f$ is purely imaginary for every ODD
offset**, so for real-ish D the screen discarded every odd-offset term and the Hartree/XC matrix came out
EXACTLY REAL — measured `maxIm(dV) = 0` at k=¼ against 0.067 at k=0.25001.  An H missing its imaginary part
has the wrong spectrum, so the SCF converged 2.5 Ha high.  Fix: screen on the true magnitude `|D_ij|`
(= `|D_ij conj(phase)|`, since `|phase|=1`) — strictly more conservative, and what the project's own
"the magnitude screen is the only truncation" rule always meant.

**How it was found**, since the sequence is the reusable part: a single-k sweep showed E(k) smooth except at
exactly ¼ and ¾; k=¼+1e-9 was fine, so it was an exact-value branch and not physics; every operator and the
1E spectrum were proven element-wise continuous (three gates, still enabled); the symptom was an extra
singlet with the Λ₃ doublet straddling E_F; a Fock-matrix fingerprint then showed the density-dependent
potential losing its imaginary part **only** at ¼; and `GPW_DENSITY_EPS=1e-30` recovered the right answer,
naming the screen.  Verification: E(k) is now smooth (0.249 → −7.563844, ¼ → −7.565529, 0.251 → −7.567208)
and **k=¾ equals k=¼ exactly**, as time reversal requires.  Si Γ, Si 2×2×2 Γ-centred and NaF Γ are unchanged
to every digit — at TRIM k the phase is real and the old test agreed with the new one.

² CP2K's NaF decks carry no `&KPOINTS` section, i.e. they are Γ.  The qchem test's own default is a **2×2×2**
mesh (8 k → 3 irreducible) — 116 mHa of band dispersion below its Γ value — so the row that compares to CP2K
is the `NAF_KMESH=1` one, and the 2×2×2 row is a qchem-only cost/size datapoint until a k-point CP2K deck
exists.  This mismatch was live in the test until 2026-08-19: it carried a Γ-era anchor (−24.4304) against a
2×2×2 configuration and simply failed.
³ the full-SR span runs FULL RANK at Γ ("kept 32 of 32", λ_min 4.35e-4) — the historical near-null trouble
was the multi-k cell.  The retracted −27.93128 "oracle" (a 3.5 Ha `EPS_PGF_ORB` screening artifact) is gone
from this table for good.

⁴ **RE-MEASURED TWICE, 2026-08-19/20** — same command, same box.  All three cuts agree on the ENERGY to
nine significant figures (−61.402976200 → −61.40297623 → −61.40297622, a spread of ~1e-8 Ha) with `m_stag`
±0.6667 over 17 iterations: both changes below are cost changes, not physics.

| cut | wall | CPU | peak RSS |
|---|---|---|---|
| banked (no folds, dense Φ) | 20m05s | 2240 s | 4947 MB |
| + **Step 2**: T3 stream fold ARMED (plan T3.5) | 13m25s | 1809 s | **1349 MB** |
| + **Step 3**: Φ-table build (screen + sparse spherical transform) | 8m36s | 1554 s | 1350 MB |
| + **Step 3**: Becke partition ε 1e-8 → 1e-6 | **6m56s** | **663 s** | 1323 MB |
| | **2.90×** | **3.38×** | **3.7×** |

**Step 2** (seconds, banked → armed): pair **scatter 263.4 → 41.0** (6.4×), pair **gather 167.6 → 23.1**
(7.3×), **stream build 110.1 → 28.2** (3.9×) — a larger factor than the 4.60× rep-pair reduction, because
the pairs the fold drops are the expensive ones.  **THE RAM CAME FROM HERE**, and it was the streams.
**Step 3**: the **Φ-table build 379.2 → 83.3 s** (4.6×; 190 → 36.8 per anneal stage).  Its dominant cause
was NOT Φ's density but a dense 122×118 cart→spherical transform applied INSIDE the lattice-image loop —
~150 mat-vecs per mesh point where one would do, 10.6% of all cycles — plus a missing magnitude screen on
the pointwise sweep.  Details and the rejected third hypothesis in `doc/OpenWork.md` Step 3.
**Becke ε** (2026-08-20): the partition's ε-converged competitor series ran at ε=1e-8, which had only ever
been probed TIGHTER.  ε fixes |im| (3183 competitor images per live point on this cell) and the partition
costs O(|P-set|·|im|) per point, so ε scales the dominant loop rather than shaving it: **1e-8 → 1e-6 takes
the becke build 36.95 → 8.33 s threaded (294 → ~66 s SERIAL), 4.44×**, at an energy identical to nine
significant figures (−61.40297622 → −61.40297621).  Margin: the binding gate is
`BeckeEquivalentSitesOwnEqualShares` (site shares equal to 1e-8 RELATIVE), which survives 1e-6 and fails at
1e-5, while the quadrature error itself is ~1e-4.  ⚠ This is a **TOLERANCE trade, not a bit-identical
restructuring** like the Φ work — the weights move at ~1e-6 relative.  `GPW_BECKE_EPS` overrides for A/B.
Against CP2K the CPU gap narrows 6.0× → 4.2× → **1.8×** and RAM 23× → **6.1×**; on WALL this row is 1.11×.

> **⚠ AND IT EXPOSED A FLAW IN THIS TABLE'S OWN CPU COLUMN.**  The becke build is 294 s SERIAL but bills
> ~590 s of CPU when threaded 16-way (36.95 s wall × 16) — because the OpenMP threads BUSY-WAIT at the
> barrier, so CPU time counts spinning as work.  Two anneal stages of that was ~1180 s of the banked row's
> 1554 s CPU (76%), which is why removing 4.44× of a "38%" bucket cut the row by 2.34×.  The user's pin
> "compare CPU, not wall" assumes CPU tracks work done; against a serial CP2K it does not, wherever qchem
> threads.  **The serial column is the honest algorithmic comparison** — which is exactly why the whole
> table was cut at one thread.  Anyone reading a threaded CPU number here should divide by the parallel
> efficiency, which still nothing measures.

**What is hot now:** the top bucket is `scf: XC-mesh ρ sampling (matrix-free)` at 84.1 s.
⚠ **THIS ROW PREVIOUSLY MIS-NAMED THAT BUCKET "the Φ-shaped GEMM, i.e. the Φ-SPARSITY item" — IT IS
NEITHER** (corrected 2026-08-20).  *Matrix-free* means exactly "carries no density matrix", so it cannot be
the DM GEMM.  Per `PWTerms.C:698`, it is the **ρ̃-MIXED density sampled on the XC mesh — a batched inverse
FT over the whole {G}, on EVERY Kerker/Pulay iteration** (plus the iteration-0 seed).  The real DM GEMM is
the separate `scf: XC-mesh ρ sampling (all iterations)` line, and on this recipe the mixer hands XC a
ρ̃-backed density from iteration 1 on, so the GEMM is nearly BYPASSED: measured 1.70 s against the
matrix-free bucket's 35.0 s on a 6-iteration cut.  The code had already split these two buckets for this
exact reason — *"lumping it into the GEMM hid the fact that the mixed-density sampling, not the GEMM, was
the iteration's largest XC cost"* — and this table re-lumped them in prose.  **The per-iteration lever is
the ρ̃-mixed sampling, not anything Φ- or D-shaped.**
⁵ the FM row is still the PRE-Step-2 measurement (its energy is unaffected, its cost is not) — re-run it
with the AFM command when the threaded cut of this whole table is taken.  **The Si and NaF rows likewise
predate Step 3**: their ENERGIES are unchanged (verified — Si Γ still reads −7.115067665 to every digit),
but any row whose XC runs on the Becke mesh may now be faster than its cost columns say.  Only the MnO AFM
row has been re-measured end to end.  Nothing in the table is stale in the Δ column, which is the column it
exists for.


---

## 3. The pair-stream value cache, and what deleting it did

### ★★★ WHAT DELETING THE CACHE DID TO THIS TABLE — read the RAM column and the MnO row together

| system | CPU before → after | peak RSS before → after |
|---|---|---|
| Si Γ | 7.4 → **2.4 s** (3.1× faster) | 267 → **28 MB** (9.5×) |
| Si 2×2×2 Γ-centred | 9.6 → 8.9 s | 269 → **30 MB** (9.0×) |
| Si 2×2×2 shifted MP | 15.6 → 17.5 s | 269 → **31 MB** (8.7×) |
| NaF SR2 Γ | 94.5 → **39.7 s** (2.4×) | 577 → **54 MB** (10.7×) |
| NaF SR2 2×2×2 | 112 → **84.5 s** | 590 → **68 MB** (8.7×) |
| NaF full-SR Γ | 219 → **43.4 s** (5.0×) | 3090 → **58 MB** (**53×**) |
| **MnO AFM-II Γ (imposed, VA)** | 663 → 976 → 837 → 678 → 620 → **584 s** (**0.88× — FASTER**) | 1323 → **491 MB** (2.7×) |

**Two rows now BEAT CP2K on CPU outright** — Si Γ at 0.48× and NaF full-SR at 0.43× — and **every** qchem
row is now well under CP2K's RAM (28–68 MB against 148–186 MB on the small cells), which is the first time
that has been true.  The NaF full-SR row is the headline: 3090 MB and 219 s CPU became 58 MB and 43 s, and
it was the row whose 3 GB used to force `scripts/memsafe`.

⛔ **AND THE MnO ROW GOT SLOWER — 1.26×, and it is not noise.**  That row is a long IMPOSED run (14+17
iterations) where the two box-walk buckets were **477 s of 976 s CPU**, so the cache's measured 2.91× on
those buckets translated almost exactly into the 1.47× the row lost when it went.  ⇒ **The cache was still buying real
CPU on the one row that matters most**, and the case for deleting it rests on the RAM axis and on the two
latent defects it was hiding, not on a free lunch.  ⚠ The 663 s "before" is the 2026-08-19 banked row on an
older binary, so treat the MnO delta as indicative; the directly-measured, same-binary A/B is the
2.91×-on-the-buckets figure in `doc/OldPlans/CollocationRewritePlan.md` step 7.
✅ **AND THE ROW HAS SINCE GONE PAST WHERE THE CACHE LEFT IT — 976 → 837 → 678 → 620 → 584 s, against
the 663 s the 4 GB cache used to buy, on 2.7× less RAM.**  Four changes did it, none needing a re-bank —
the box-walk buckets went 477 → 344 (the `template<int LP>` dispatch) → 192 (the **collocation memo depth
fix**) → 123 (the **gather memo**) → **91 s** (the **exp-table recurrence**).  The first three are
bit-identical by construction; the fourth is not, and moved nothing anybody pins (below).  ⇒ Against CP2K
this row now stands at **1.57× CPU** (was 2.24×), **2.22× CPU per ITERATION** (was 3.2×) — and **0.88× on
WALL, i.e. faster in wall clock.**

---

## 4. The collocation memo had depth 1

### ★★★ THE COLLOCATION MEMO HAD DEPTH 1, AND A POLARIZED RUN ALTERNATES TWO DENSITIES (2026-08-28)

Found by a CALL CENSUS — bucketing the four closure sites so the ledger reports a per-site call count.  On
the MnO row: **64 closure calls, 4 memo hits — 6%.**  `sameD` compared against the LAST collocation only,
so \f$D_\uparrow\f$ and \f$D_\downarrow\f$ evicted each other every single call.  The unpolarized runs the
memo was written against never showed it (Si Γ collocates 1.45×/iteration; MnO was collocating ~20×).

Keeping four older (D, ρ) pairs — same EXACT match rule, so a replay is bit-identical — on the 3-iteration
probe:

| | collocations | memo hits | bucket |
|---|---|---|---|
| depth 1 (as shipped) | 60 | 4 (6%) | 36.7 s |
| **depth 5** | **16** | **48 (75%)** | **9.8 s (3.74×)** |

16 misses is the true number of DISTINCT densities, so depth 5 catches every repeat.  On the full row:
collocate calls **368 → 120**, its bucket **221 → 71 s**, `Etot` unmoved at −61.40297551, +10 MB of RSS.
⚠ This is not an acceleration CP2K lacks — it is removing OUR OWN redundancy; CP2K collocates ρ once per
step.  Knob `GPW_COLLOC_MEMO` (0 restores depth 1).


---

## 5. The gather census — 53% of the integrates were exact duplicates

### ★★★ AND THE SAME CENSUS ON THE GATHER SIDE: 53% OF THE INTEGRATES WERE EXACT DUPLICATES (2026-08-28)

`GPW_INTEGRATE_CENSUS=1` classifies each gather against a short history of (field, screen) hashes.  On the
MnO probe: **30 gathers, 14 distinct, 16 with V AND screen bit-identical** — and not one "same V, widened
screen", so the repeats were pure redundancy, not a screening artefact.

The molecular `IntegrateMemo` could not catch them, for two *correct* reasons: screened calls bypass it (its
key is \f$V\f$ alone while the active set moves with \f$D\f$), and FOLDED calls bypass it (its `nb` records
\f$b\f$ without the orbit multiplicity and the replay never runs `fillImages`).  An imposed polarized run is
both, so it never memoized at all.  ⇒ The fix caches the FINISHED \f$h\f$ on \f$(V_L,\ \text{screen})\f$
in `GPW_Evaluator` — per k-block by construction, which is what makes it exact whatever route produced
\f$h\f$.  Gathers **184 → 74**, bucket **121 → 50 s**, `Etot` unmoved.

⚠ **AND THE PROFILE HAS RE-ORDERED — the box walk is no longer the biggest block on this row.**  ⚠ "This
row" means the DEFAULT row (`CP2K_COMPAT=0`), which runs Becke; the parity row is a different measurement
and is dealt with below.

| block | s (of 367 s wall) | share |
|---|---|---|
| **Becke XC mesh** (ρ sampling 74 + Φ tables 42 + H_xc 23 + mesh build 17 + 3) | **158** | **43%** |
| collocate + integrate-back | 123 | 33% |
| unaccounted (LA, mixing, FFT) | ~75 | 20% |
| local-PP, 1E sums, closures | ~11 | 3% |

⛔ **AND MEASURING IT INVERTED WHAT THAT SHARE SEEMS TO IMPLY — BECKE IS A 3× NET WIN, NOT AN ADVERTISEMENT.**
"43% of the row" does NOT mean "43% to be had by turning it off".  The uniform grid this system needs is
**571787 points against Becke's 48320** — 11.8× more — which is exactly why the `Auto` selector picks Becke
here.  Measured 2026-08-28 with `QCHEM_BECKE_XC=0` and everything else at default (the imposition KEPT, so
the run still converges — full `CP2K_COMPAT` does not, see below):

| MnO AFM-II Γ, VA, imposed | Becke (default) | **uniform XC** |
|---|---|---|
| Etot | −61.40297551 | −61.40295935 (1.6e-5 Ha — the quadrature difference) |
| iterations | 17 | 13 |
| **CPU** | **584 s** | **1805 s (3.09× SLOWER)** |
| wall | 5m28.1s | 30m10.9s |
| **peak RSS** | **491 MB** | **4494 MB (9.2×)** |
| the XC buckets | 158 s | **1592 s** (ρ sampling 874 + Φ tables 488 + H_xc 230) |

⇒ **The problem was not that Becke is an unfair advantage; it was that our UNIFORM route — the parity
route — was ~10× dearer than Becke on this cell.**  ✅ **FIXED 2026-08-28, and it was a CONFLATION, not an
algorithm.**  `XC_PairQuadrature` — ρ via the GPW collocation, \f$v_{xc}\f$ pointwise, \f$H_{xc}\f$ via the
exact transpose, i.e. **CP2K's algorithm, no Φ table anywhere** — already existed and Si Γ already used it.
But `VxcFit::Auto` read `becke || polarized`, and that `polarized` half forced EVERY polarized run onto the
Φ table whatever its grid, because `XC_PairQuadrature::RhoPol` threw.  Nothing spin-specific was ever in the
way: `applyRaw` takes a \f$D\f$, so a channel is one call with that channel's density, and the adjoint
needed no change at all.

| MnO AFM-II Γ, VA, imposed, `QCHEM_BECKE_XC=0` | Φ-table route | **collocation (pair) route** |
|---|---|---|
| **CPU** | 1805 s | **246.5 s — 7.3× faster** |
| wall | 30m10.9s | **4m05.4s** |
| **peak RSS** | 4494 MB | **105 MB — 43× smaller** |
| the XC buckets | 1592 s | **5.2 s** |
| Etot / iterations | −61.40295935 / 13 | −61.40358773 / 25 |

★ **And on that route qchem BEATS CP2K on both axes — 246 s against 373 s CPU, 105 MB against 217 MB.**
⚠ It is still not a parity ROW: the imposition, the low-rank ρ and the stream fold are all still on.  What
it says is that the parity ROUTE is no longer the thing standing in the way.
⚠ The two uniform-XC energies differ by 6.3e-4 Ha because they are different discretisations of the same
quadrature; the pair route is the variational one (\f$H_{xc}=\partial E_{xc}/\partial D\f$ to machine
precision, gate `GPW.RawXCConsistencyFD`).

⇒ The atom-centred XC quadrature is the largest single cost of the DEFAULT row, and it is a **qchem-only
algorithm** — CP2K runs XC on the uniform grid.  ✅ **DECLARED 2026-08-28** as `QCHEM_BECKE_XC` (`RunPolicy::BeckeXC`), so
it is on the deviation table above and `CP2K_COMPAT=1` routes XC to a **basis-sized** uniform grid — the
sizing stays in `qcMesh::ResolveXCMesh`, because handing back a bare `cellKind` would leave `nUniform`'s
basis-blind default of 20 in charge, which is the under-resolution that selector exists to prevent.
⚠ The DEFAULT path is untouched by the change (with `BeckeXC()` true the new overload delegates to the
unchanged two-argument resolver), and that is measured, not asserted: Si Γ reproduces −7.115067844 to all
10 s.f. and `ctest -j8` is 793/793.  What DOES move is every `CP2K_COMPAT=1` row — of which this table has
exactly one, already marked ⚠ STALE.
⇒ `raster` and `cutoffFactor` remain the last two typed options outside the policy (N5).


---

## 6. The exp-table recurrence — the one place the kernel did not follow CP2K

### ✅ THE EXP-TABLE RECURRENCE — the one place the kernel did not follow CP2K (2026-08-28)

The table build was \f$3n^2\f$ scalar `std::exp` calls where CP2K uses a recurrence.  The exponent is
LINEAR in the inner index of each table (only the cross term \f$2e_ae_bh_{ab}\f$ moves), so a row is one
seed and \f$n\f$ multiplies — **2 exps per row instead of \f$n\f$**.  Seeded at the LARGEST entry and
walked downward, which is the 2026-08-26 underflow rule: a recurrence seeded in the tail can start below
the underflow floor and stay zero through entries that matter.

| | kernel at 32³ | box walk, Si Γ | NaF SR2 Γ | MnO row |
|---|---|---|---|---|
| direct `exp` | 111.8 µs | 0.290 s | 1.71 s | 123 s |
| **recurrence** | **72.9 µs (1.53×)** | **0.204 s** | **1.23 s** | **91 s** |

⚠ It is **NOT bit-identical** — a product of \f$n\f$ rounded factors is not the rounded product — so it was
built default-OFF and defaulted ON only on evidence: against a naive exact reference the contraction goes
1e-15 → **7e-15 relative, flat in box size** (it is \f$n\varepsilon\f$), still **4× better than the WALK's
3e-14**; `ctest -j8` is **793/793 on both settings**; and all three anchors above are unchanged **to all 10
printed s.f.**  ⇒ Anchor-moving in principle, moved nothing in practice.  `GPW_EXP_RECURRENCE=0` is the A/B.

⚠ **AND IT ONLY HALVED THE TABLE TERM, not eliminated it** (89 → 45 µs): with the exps gone the build is
bound by writing \f$n^2\f$ entries, which no algorithm removes — the table has to exist.  The kernel's
\f$O(N^2)\f$ share is 80% → **62%**.

⇒ **And the next KERNEL lever is still NOT the batching.**  Before the recurrence, **84% of the contraction
kernel was the three 2-D Mathieu `exp` tables** (\f$3n^2\f$ scalar `std::exp` calls), worth up to ~35% of this
row.  CP2K builds those tables by RECURRENCE where we call `exp`, and doing the same is the one route that
keeps this table apples-to-apples — a vectorised `exp` would be a qchem-only acceleration under rule 3
above, declared on the deviation line and switched OFF for every head-to-head row, i.e. speed we could not
quote here (user, 2026-08-28).

**Energies moved by ≤ 1e-6 Ha everywhere** — Si Γ −6.6e-7, Si 2×2×2 +1.9e-8, shifted MP 0 (10 s.f.),
NaF SR2 Γ +6.5e-9, NaF SR2 2×2×2 −3.1e-7, NaF full-SR +3.2e-6, MnO −7.0e-7 — which is the anchor re-bank
A1+A7 predicted (doc/OpenWork.md, the anchor-moving sprint) and it does not change any verdict in the Δ
column.


---

## 7. Fold state and cost profile of the MnO rows

### Fold state and cost profile of the MnO rows (their provenance, and Steps 2–3's target)

The MnO rows above run `MNO_IMPOSE=1`, so unlike the FREE production run they DO fold — printed by the run
itself (Step 0b):

| site | this row | free production run |
|---|---|---|
| XC mesh (Becke star-average) | **23.03×** (24 ops, magnetic/Shubnikov) | NONE |
| `V_loc`-long {G}-star | **10.37×** (12 ops) | NONE |
| collocation streams (T3 pairs) | **4.60×** (12 ops, 8778 → 1909 rep pairs) ⁴ | NONE |

So the imposed row now claims all three folds.  The stream fold is **12 ops, not 24** — the S3 pin: a
magnetic imposition may fold the PER-CHANNEL streams only under the σ=None (sublattice-preserving) subgroup,
since a flip op relates D↑ to D↓.  And 4.60× on a 4-atom cell is not the 71× the diamond gate cell showed:
the orbit factor is a property of the cell's symmetry, so the headline number belongs to a cell, never to
the feature.  **The FREE production run still folds NOTHING at any of the three sites** — that is a
`MNO_IMPOSE` decision, not a plumbing gap.

Where the time goes (AFM arm, `GPW_REPORT=1` ledger) — BEFORE the fold was armed (of 1205 s) and AFTER
(of 805 s):

| bucket | s (banked) | s (+Step 2 fold) | s (+Step 3 Φ) |
|---|---|---|---|
| setup: XC-mesh **Φ tables** | 370.0 (31%) | 379.2 (47%) | **83.3 (16%)** |
| scf: collocate density (pair scatter) | 263.4 (22%) | 41.0 (5%) | 42.1 (8%) |
| scf: integrate-back (pair gather) | 167.6 (14%) | 23.1 (3%) | 23.3 (5%) |
| setup: collocation stream build | 110.1 (9%) | 28.2 (4%) | 29.0 (6%) |
| setup: Becke mesh build | 70.6 (6%) | 69.4 (9%) | 71.5 (14%) |
| **scf: XC-mesh ρ sampling (matrix-free)** | 57.1 (5%) | 76.8 (10%) | **84.1 (16%)** |
| scf: XC-mesh H_xc quadrature | 33.8 (3%) | 44.1 (5%) | 28.3 (5%) |

**The profile has flattened and the head has moved.**  The pair loops (Step 2) and the Φ BUILD (Step 3) are
each down to single-digit-to-16% shares, and the largest bucket is now the **Φ-shaped ρ GEMM** — the
Φ-SPARSITY item, which is what `OpenWork` Step 3's "Φ-table screening" was always aimed at and which the
build's cost used to hide.  Second is the **Becke mesh build**, whose 71 s of wall is ~50% of the run's
CYCLES (it is the one loop threaded by default) — so on the CPU column, which is the honest one, the Becke
partition is now the single biggest item in the code.  ⁶ the XC-mesh buckets are not comparable one-to-one
across cuts (iteration counts differ per cut); per-iteration cost is what matters there.
⁶ the XC-mesh buckets are not comparable one-to-one across the two runs (the folded run converged its second
stage in 17 iterations; per-iteration cost is what matters there, and it did not change).


---

## 8. Superseded whole-run reading of the three MnO rows (2026-08-28)

★ **THE THREE `MnO AFM-II` ROWS ARE ONE SYSTEM AND ONE RECIPE, with progressively more of OUR deviations
switched off.**  Read them as a ladder, not as three experiments — the only thing changing is which of the
six declared deviations are on:

| the row | what is off | CPU | RSS |
|---|---|---|---|
| **ALL DEFAULTS** | nothing — every deviation at its qchem default (Becke XC, imposition, low-rank ρ, stream fold) | 584 s | 491 MB |
| **`QCHEM_BECKE_XC=0`** | the Becke XC mesh only — so XC runs CP2K's way, everything else still ours | **246 s** | **105 MB** |
| **`CP2K_COMPAT=1`** | ALL SIX — the honest algorithm-to-algorithm row | 2736 s | 112 MB |

⚠ The middle row is the fastest because Becke's atom-centred quadrature is expensive machinery, and the
bottom row is the slowest because parity also removes the STREAM FOLD (5.2× on MnO's pair count) and the
low-rank ρ.  ⇒ **Always say WHICH of the three** — "the MnO row" has meant all three of these in the space
of one day, and they differ by 11× in CPU.


⁷ **THE PARITY ROW — IT EXISTS AGAIN, AND THE OLD VERDICT IS RETRACTED (re-measured 2026-08-28).**
`doc/OpenWork.md` and this file have carried *"AT TRUE PARITY OUR MnO RECIPE DOES NOT CONVERGE ... stage 2
caps at −57.620, 3.8 Ha short"* since 2026-08-26.  Re-run on today's tree, same recipe, `CP2K_COMPAT=1`
(so the imposition, the low-rank ρ, the stream fold AND the Becke mesh are all off):

| | 2026-08-26 | **2026-08-28** |
|---|---|---|
| stage 1 | caps at 80, −60.431 | UNSETTLED at 13, −61.41070717 |
| stage 2 | caps at 80, **−57.620** | **FIT-FLOOR STALL at 80, −61.39789688** |
| against the imposed answer (−61.40297551) | **3.8 Ha short** | **5.1 mHa short** |
| the AFM order | — | **SURVIVED both stages** (m_stag 0.636, 0.610) |
| peak RSS | 5034 MB | **112 MB** |

⇒ **It still hits the cap, but for a completely different and far more benign reason.**  It is no longer
collapsing: the energy is settled (ΔE amplitude 2.7e-8), the magnetic order holds without any imposition,
and Δρ has FLOORED at 1.69e-5 against a 1e-5 target — the run's own detector calls it a *"FIT-FLOOR STALL
(Δρ floored, ΔE tiny -- functional/grid)"*, not an oscillation.  That is **A4 territory** (the Δρ/N
convergence gate, `doc/Records/SCFStrategyPlan.md`), which is the third independent thing this session has pointed
at A4.

⚠ **THE HONEST PARITY STANDING, then**: 2736 s against 373 s CPU is **7.3×**, or **3.5× PER ITERATION**
(29.4 s against 8.5 s) since we took 93 iterations to CP2K's 44 — and **0.52× on RAM**.  The per-iteration
gap is bigger than the default row's 2.22× because parity also strips the stream fold, so each collocation
covers ~5× more pairs (2.2 s/call against 0.5 s).  ⇒ **The fold is now the largest single thing the
comparison removes**, which is exactly the qchem-vs-qchem delta rule 3 asks to be reported separately.

⁶ **THE XC-PARITY ROW, and the first MnO row on which qchem beats CP2K on BOTH axes** — 246 s against
373 s CPU and 105 MB against 217 MB.  It is the default recipe with ONE deviation removed: the Becke mesh,
so XC runs the way CP2K runs it (`XC_PairQuadrature` — ρ via the GPW collocation, \f$v_{xc}\f$ pointwise,
\f$H_{xc}\f$ via the exact transpose).  ⚠ **It is NOT the parity row**: the imposition, the low-rank ρ and
the stream fold are all still on.  What it establishes is that the parity ROUTE is no longer what stands in
the way — before 2026-08-28 the same configuration cost 1805 s and 4.5 GB, because a polarized run could not
reach that route at all.


---

## 9. The superseded 2026-09-04 per-iteration cut

Replaced 2026-09-05 by `doc/Benchmark.md` §5a's single bin-1 table, which re-took every row on one binary
at a MEASURED one thread and split setup out of the per-iteration figure.  Kept here because the 09-04
numbers are quoted in commit messages and in `doc/OpenWork.md`, and because two of the rows moved for a
reason worth remembering: **the 09-04 rows were taken with the gather's density screen still in**, and
removing it (`a7561e92`) changed both the iteration counts and the per-iteration cost.

★★★ **THE SMALL ROWS, DECOMPOSED (2026-09-04)** — the table this replaced:

| row | CPU banked | **CPU now** | qchem iters | CP2K iters | **qchem s/iter** | CP2K s/iter | **per-iter ×** | old whole-run × |
|---|---|---|---|---|---|---|---|---|
| Si Γ | 2.2 s | **0.69 s** | 11 | 12 | **0.063** | 0.417 | **0.15×** | 0.44× |
| Si 2×2×2 Γ-centred | 8.9 s | **3.04 s** | 7 | 13 | **0.434** | 0.431 | **1.01×** | 1.6× |
| Si 2×2×2 shifted MP | 17.5 s | **6.57 s** | 16 | 14 | **0.411** | 0.429 | **0.96×** | 2.9× |
| NaF SR2 Γ | 37.9 s | **26.3 s** | 29 | — | 0.906 | — | — | 5.3× |

★ **THE MnO CODEGEN ARMS (2026-09-04)** — `-O3` against `-O3 -march=native`, the measurement that made
native the Release default.  Its `×` column is per-iteration TOTAL CPU against CP2K's 8.5 s/step, i.e. the
setup-contaminated figure §5a now separates:

| the row | iterations | banked 08-28 | **O3** | **NATIVE** | \f$E_{tot}\f$ |
|---|---|---|---|---|---|
| **ALL DEFAULTS** | 14+17 = **31** | 18.2 s (2.14×) | **13.78 s (1.62×)** | **12.45 s (1.47×)** | −61.40297551 = banked |
| **`QCHEM_BECKE_XC=0`** | 14+25 = **39** | 6.31 s (0.74×) ᵃ | **4.36 s (0.51×)** | **3.99 s (0.47×)** | −61.40358773 = banked |
| **`CP2K_COMPAT=1`**, `GPW_MNO_NMAX=10` probe | 10+10 = **20** | 19.3 s (2.3×) | **17.89 s (2.11×)** | **16.94 s (1.99×)** | (capped, both arms equal) |
| ⚠ `CP2K_COMPAT=1`, the full table row | 13+80 = **93**, both CAPPED | 29.4 s | not re-taken (45 min) | — | |
| CP2K (its own log) | 44 | 8.5 s | — | — | |

ᵃ the banked `QCHEM_BECKE_XC=0` row recorded 246 s CPU but only stage 2's iteration count; stage 1 is
assumed 14, so 6.31 s/iter is derived, not banked.

⇒ **WHAT MOVED against 08-28**: ALL DEFAULTS −24.3% CPU, `BECKE_XC=0` −30.9%, the parity probe −7.3%.
⚠ That was "current vs banked", NOT an attribution — other work landed between 08-28 and 09-04 and no
parent-commit A/B was run on this recipe.


---

# Benchmark.md verbatim as of 2026-10-01

Archived in full on 2026-10-01 when `doc/Benchmark.md` was trimmed to the live instrument (what is left, how to
run it, the lessons).  Section numbers are the ones other docs cite.  Nothing here is a current number.  The
copy sits inside one `~~~~` fence so it is byte-for-byte the old file (its own ``` blocks included); heading
anchors for it are the `## N.` lines inside the fence.

~~~~markdown
# The qchem ↔ CP2K head-to-head — runtime, RAM, energy

**This is an INSTRUMENT, not a report.**  A row here is a claim that two codes did the same work on the
same hardware; the sections below are what makes that claim checkable.
📖 **The reasoning, the wrong turns and every superseded number are in `doc/Records/BenchmarkHistory.md`** (split
out 2026-09-04).  Come back here for what is TRUE NOW; go there for why.

**Read in this order:** §1 the process → §2 what parity actually means → §3 the rules → §4 the commands →
**§5a the BIN-1 TABLE** → §5 the whole-run rows → **§5f why the parity row is 2.05×**.
★ **§5a is the one table to point at for per-iteration CPU** — it holds nothing else, by request
(2026-09-05); the deltas and open questions that used to crowd it are in §5e.

---

## 1. THE SYSTEMATIC PROCESS — single thread first, then threads; and the four bins

### 1a. Single-thread parity FIRST, then threads — in that order (user)

The systematic procedure, and the reason for it: a threaded comparison against a serial code measures two
different things at once — the algorithm and the parallel efficiency — and if the single-thread gap is
unknown you cannot tell which one you are looking at.  So:

1. **Get qchem's SINGLE-THREAD time roughly in line with CP2K.**  `OMP_NUM_THREADS=1 GPW_OMP_THREADS=1`,
   which pins both our own OpenMP regions and the BLAS underneath blaze.  This is the honest
   algorithm-to-algorithm number and the one that should drive optimisation.
2. **THEN turn on N=8 or 16 and look for OMP-SHAPED GAPS** — poor scaling, barrier waits, regions that do
   not thread at all.  Those are a different class of defect with different fixes, and they only read
   cleanly once step 1 is settled.

⚠ Related and already measured: our OpenMP threads **BUSY-WAIT at the barrier**, so a threaded qchem run
bills far more CPU than it uses (a 294 s serial build billed ~590 s CPU at 16 threads).  That is why the
CPU column overstates qchem wherever it threads — and another reason the single-thread row is the clean one.


### 1b. The four gap-close bins (user, 2026-08-28) — priority order

> *"1) Per iteration CPU time, 2) Init (pre iteration) time, 3) top RAM usage (mostly solved I think),
> 4) very roughly match the # of iterations.  I think items that fall into bin 4 can just be documented."*

| bin | the axis | where it stands (MnO AFM-II VA, 2026-09-04) |
|---|---|---|
| **1** | **per-iteration CPU** | ★ **PER SCF ITERATION WE ARE AHEAD ON 7 OF 9 ROWS** (§5a, 2026-09-05, both codes pinned serial): MnO ALL DEFAULTS **0.82×**, MnO FM **0.77×**, `QCHEM_BECKE_XC=0` **0.45×**, Si **0.14× / 0.66× / 0.56×**, NaF full-SR **0.18×**.  The two losses are NaF SR2 **1.45×** (a Becke cost) and — the one that counts — `CP2K_COMPAT=1` **1.72×** (was 2.05× before §5f's lever A).  ⇒ **NOT closed** — but on the SAME ALGORITHM (our fixed-point stage against CP2K's diagonalise-and-mix decks) it is **1.13×**, three gathers against their two, and the third gather is §5f lever B, behind N4.  The 1.72× two-stage figure adds a GDM line search CP2K's runs have no counterpart for |
| **2** | **init / pre-iteration** | ⛔ **PROMOTED, and now measured serially: MnO's setup is 184 s = 47% of the default run against CP2K's 8.1 s (23×)**, of which **136.6 s is two Becke mesh builds** (§5a).  The earlier "~57 s of 328 s" was a THREADED ledger bucket — the build's partition loop is `omp parallel for`.  ✅ With `QCHEM_BECKE_XC=0` our setup is **1.76 s and beats CP2K's 8.1 s** — so bin 2 is a Becke-mesh question, exclusively |
| **3** | **peak RAM** | ✅ **solved, and we win**: ~476 MB defaults, **113–132 MB on the parity routes against CP2K's 217 MB** |
| **4** | **iteration count** | ⚠ **RE-JUDGED 2026-09-20 under rule 3f** — the old "31 (defaults) / capped (parity) against CP2K's 44" compared a 1e-5 G-space residual with CP2K's 1e-6 max\|ΔP\|, i.e. two different questions.  On the SAME measure, threshold and loop shape (one density-side history, no Fock DIIS, no MOM): **imposed 22 / 24, free (`CP2K_COMPAT=1`) 33 / 37 against CP2K's 44 / 104** (U=0 / U=4 eV, §5 rows ⁸).  ⇒ CLOSED in our favour; what remains is bin 1's third gather |

The live tracker for these is `doc/OpenWork.md` (accuracy/features) and §10 below (perf levers); this file holds the measurements behind them.

---

## 2. WHAT `CP2K_COMPAT=1` ENCOMPASSES — ⚠ AN EMERGING LIST, NOT A FINISHED ONE

**Every deviation found so far is an ACCELERATION, not physics**: turning them all off moves the MnO total
by **3e-8 Ha** (agreeing to 10 s.f.).  That is the property the switch most needed to demonstrate about
itself, and it is measured rather than asserted (history §2).

⚠ **THIS LIST HAS GROWN EVERY TIME SOMEONE LOOKED.**  It started at four items; it is eight (the eighth, +U's
form, is the first that is physics rather than an acceleration — inert unless a Hubbard manifold is on the run,
so every banked row below is untouched by it).  Assume it is still incomplete — a parity row is only as honest as the last thing we noticed we were doing differently.

| # | knob | what qchem does that CP2K does not | found |
|---|---|---|---|
| 1 | `QCHEM_DM_LOWRANK` | factored/low-rank ρ (\f$D=LL^\dagger\f$) — a SINGLES route; CP2K collocates PAIRS | 08-25 |
| 2 | `GPW_STREAM_FOLD` | orbit fold on the collocation pair streams (5.2× on MnO's pair count) | 08-25 |
| 3 | `QCHEM_MIX_RHO_M` | (ρ,m) mixing channels instead of (up,dn) | 08-25 |
| 4 | `GPW_XC_DM_SOURCE` | \f$V_{xc}\f$ fed ρ[D] wholesale instead of ρ_mix | 08-25 |
| 5 | `QCHEM_IMPOSE_SYMMETRY` | space-group imposition: BZ fold + ρ star-average + site-adapted XC mesh.  ⚠ **OVERRULES the caller** | 08-26 |
| 6 | `QCHEM_BECKE_XC` | atom-centred (Becke) XC quadrature instead of the uniform grid.  ⚠ **OVERRULES the caller**; it was **43% of the MnO row** | 08-28 |
| 7 | `GPW_DAWARE_SCREEN` | D-aware collocation box tolerance \f$\varepsilon/|c_{ij}|\f$ instead of flat \f$\varepsilon\f$ | 09-04 |
| 8 | `QCHEM_U_EIGEN` | **the first PHYSICS deviation, not an acceleration** (user ruled it onto this list 2026-09-20): DFT+U evaluated on the Löwdin block's EIGENVALUES (Dudarev, rotationally invariant) where CP2K keeps only its DIAGONAL POPULATIONS (`dft_plus_u.F`, `IF (isgf == jsgf)`).  Inert on a run without a Hubbard manifold; on MnO VA at U=4 eV it is 0.08 vs 0.61 Ha of E_U | 09-20 |

**Verified locally** that CP2K does none of the symmetry work: the 1129-line `bench_MnO_AFM2_VA_cp2k.log`
contains **zero** occurrences of "irrep", "symmetry" or "point group"; QuickStep keeps K and P as DBCSR
sparse ATOM-BLOCK matrices and blocks by atom-pair SPARSITY, not by irrep.  The one symmetry knob that
exists is BZ-side and defaults OFF (`BRILLOUIN| K-Point point group symmetrization  OFF`).

⇒ **These rows compare a SYMMETRY-exploiting code against a SPARSITY-exploiting one on a small,
high-symmetry cell — the regime that most favours us.**  A 100-water box would invert it; CP2K's design
centre is large disordered systems where every group has order 1.

**THE OVERRIDE RULE**: an explicitly-set knob WINS over `CP2K_COMPAT=1`, and the banner marks it `(stated)`
so a run that thinks it has parity and does not says so out loud.  ⇒ **A new accelerator is NOT FINISHED
until it appears on that deviation line** (`qchem::RunPolicy`).

---

## 3. THE RULES a comparable row must satisfy

### 3a. COPY the command from §4 — never reconstruct it

⛔ The MnO rows need **eight** env vars.  Drop them and you silently get a different basis
(`VALENCE_LOWQ_SR`, Cartesian, not VA/spherical) and a single SCF stage instead of the two-stage anneal —
a different system, on which CP2K's 8.5 s/step means nothing.  This cost a full retraction on 2026-09-04
(history §1).
★ **THE TELL**: that probe is documented as *"~6 minutes"*; the bad runs took **90 seconds**.  **A large
discrepancy against the DOCUMENTED cost of the SAME probe is a CONFIGURATION difference until proven
otherwise — it is not your speedup.**
✅ **The check that catches it in one run**: a correct recipe reproduces the banked \f$E_{tot}\f$ AND
iteration count to all digits.  If it does not, you are not running the banked system.

### 3b. EVERY timing table states its thread state — per table, for BOTH codes

Not in prose somewhere above it: in the caption or a column, so a row cannot be quoted out of context.
Today's tables carry a prose warning and it is not enough — the warning is right (CP2K is genuinely serial
here at 97–99% CPU; qchem is not, at 115–239%) but a reader lifting one row will not carry it.

✅ **THE BANNER EXISTS AS OF 2026-08-26** (doc/OpenWork.md T5/N5).  `qchem::SolidCalculation` prints it
UNCONDITIONALLY — no verbose flag, no opt-in — four lines at construction plus one per `Converge`:

```
[<label> run] system: 4 atoms, 26 valence e, multiplicity 1 (POLARIZED), seed=IonicSAD
[<label> run] grids: densityEcut=auto C=2 raster=BallOnly xcMesh=Becke (nR=40 L=29)
[<label> run] symmetry: IMPOSED (Shubnikov from the decoration);  threads: OMP_NUM_THREADS=1 GPW_OMP_THREADS=1 (BLAS pinned to 1)
[<label> run] CP2K_COMPAT=0 -> DEVIATING;  QCHEM_DM_LOWRANK=on*  GPW_STREAM_FOLD=on*  QCHEM_MIX_RHO_M=off  GPW_XC_DM_SOURCE=off   [* = differs from CP2K]
[<label> scf] mixer: Kerker(G0=1.000000) alpha=0.45;  XC rho source: rho_mix;  accel: Ladder;  kT=0.005 MOM=on NMaxIter=80
```

**So a row is now taken by COPYING those lines beside it, not by remembering what was set.**  The
deviation line is generated from the SAME policy object the factories consult (`qchem.RunPolicy`), so it
cannot drift from what was actually built, and `CP2K_COMPAT=1` turns the whole set off in one place.
⚠ The standing rule is unchanged and now has somewhere to point: **a new accelerator is not finished until
it appears on that deviation line.**


### 3c. qchem-ONLY accelerations are OFF for a head-to-head row (user)

If qchem runs an algorithm CP2K does not, the row stops being a comparison and becomes an advertisement.
Turn them off, take the comparable number, and report the acceleration SEPARATELY as a qchem-vs-qchem
delta — which is the more useful statement anyway.

| flag | default | for a head-to-head row |

The full list is §2; the mechanism is `CP2K_COMPAT=1`.  ⇒ Report the acceleration SEPARATELY as a
qchem-vs-qchem delta — which is the more useful statement anyway (§6).

### 3d. The `CPU` column is a WHOLE-RUN TOTAL — divide by the iterations

Two runs in the table differ by 3× in iteration count, so a run that needs twice the steps looks twice as
slow on a per-ITERATION basis it may actually be winning.  ⇒ Quote **CPU/iteration** for bin 1.
⚠ And per-iteration is not comparable across different ITERATION CAPS either (the stage mix changes the
GDM line-search); the cap-independent measure is **per CALL**, straight off the ledger.

### 3e. Read gather MISSES, not closure calls

The ledger prints both and **only misses are work**.  A cost estimate taken off call counts is wrong
whenever a memo sits underneath — that is how a 2026-09-04 fix was over-estimated **70×** (history §1b).

### 3f. An ITERATION COUNT is comparable only on the SAME CONVERGENCE MEASURE at the SAME THRESHOLD (user, 2026-09-20)

*"Anytime I look at top and CP2K finishes many minutes before ITMain, my antenna goes up."*  It should: the two
codes were converging DIFFERENT QUANTITIES.  CP2K's `EPS_SCF` is **max|P_out − P_in| over the AO density-matrix
elements**, per spin, un-normalised (`qs_scf_loop_utils.F`, `self_consistency_check`), typically 1e-6 (MnO deck)
or 1e-7.  Our `MinΔρ` on a Kerker/Pulay recipe is the MIXER's residual **max|ρ̃_out(G) − ρ̃_in(G)|** in G-space
(the largest coefficient is N/Ω ≈ 0.08 on MnO, so 1e-5 there is ~1e-4 relative — **~100× looser** than the deck),
and on a linear-D recipe it is a third thing (Σ_blocks‖ΔD‖_F/N_e).  Measured on Si Γ (2026-09-20): the Kerker gate
"converges" in **5** iterations at 1e-3 on the residual (E = −7.114894); CP2K's measure at 1e-6 needs **27**
with Pulay(8) and lands on −7.115068 — 0.17 mHa lower, the true anchor.  On MnO +U the honest count was not
"55 vs 104".  ⇒ `SCFParams::Measure::MaxΔD` puts CP2K's measure on our loop (successive D_out, max over blocks
and spins); **a row that quotes an iteration count, a wall, or a "converges in" sets `Δρmeasure=MaxΔD`,
`MinΔρ=EPS_SCF` and `NMaxIter=MAX_SCF` off the deck** — and states it.  ⚠ Every iteration-count column in §5
written before this rule compares a looser qchem criterion with a tighter CP2K one; re-take before quoting.

---

## 4. HOW TO PRODUCE A ROW — the same wrapper for both codes

⚠ **THE `GPW_SCF.*` FILTERS IN THIS SECTION NAME RETIRED TESTS** (the 2026-09-15 test-suite renaming, rule 3a's
own failure mode).  The live rows and their current filters are in **`scripts/retake5a`** — `GPW_Si.Γ_CP2K`,
`GPW_Si.k222_CP2K`, `GPW_Si.k222s_Imp_CP2K`, `GPW_NaF.Γ_Imp_Anchor`, and `gpwprobe mno` for MnO — which also
sets the rule-3f knobs (`GPW_MEASURE=maxdd GPW_EPS GPW_NMAX GPW_ACC=null GPW_MOM=0 GPW_KERKER_G0=1 GPW_PULAY=8`,
all read by the harness's `EnvOverrides`).  Run `scripts/retake5a` to reproduce the rule-3f table in §5a.

```bash
scripts/bench "Si Gamma qchem" -- build/Release/IntegrationTests/ITMain --gtest_filter=GPW_SCF.SiliconGammaConverges
scripts/bench "Si Gamma cp2k"  -- cp2k -i IntegrationTests/CP2K/si_fcc_gpw.inp
```

**Peak RAM is measured from OUTSIDE, identically for both codes** (`/usr/bin/time -v` → "Maximum resident
set size" = the kernel's `VmHWM` high-water mark).  That is the whole reason the RAM column is now
obtainable without touching CP2K: no in-program hook was ever needed, and measuring both sides through one
wrapper is what makes the two columns comparable.  On a qchem row the wrapper also prints qchem's OWN
internal `VmHWM` as a cross-check — they should agree (measured: 266 vs 266.7 MB on Si Γ).

For the qchem detail (per-bucket ledger, folds, task-list geometry) add `GPW_REPORT=1`:

```bash
GPW_REPORT=1 build/Release/IntegrationTests/ITMain --gtest_filter=<test> --gtest_also_run_disabled_tests
```

Every GPW run now prints, without extra flags:

| what | where it comes from |
|---|---|
| `Etot` at **10 s.f.**, iterations, convergence verdict | the run fingerprint + summary lines |
| wall time, and the per-bucket ledger | the `timing` section |
| **PEAK RSS (MB, process high-water)** | `timing`, from Linux `VmHWM` (added 2026-08-19) |
| `[fold] <site>: … = F×` for every fold site | `EmitFold`, Step 0b |
| `[collocation] kernel=… ;  task list: … tasks, … MB` | the kernel + task-list readout (unconditional) |
| `[site moments] … [e]` (polarized) | `QCHEM_SITE_MOMENTS=1`, Step 0a |
| **`[t=12.34 s]` on every section heading, fold and log line** | the run clock, `report::RunElapsed()` (added 2026-09-06, `doc/OpenWork.md` Step 0c) |

⚠ **A STAMP IS A MOMENT, NOT A DURATION — AND A GAP BELONGS TO WHAT RAN *BEFORE* IT, NOT TO THE LINE THAT
CARRIES IT.**  The stamped line is the END of the silent stretch above it.  So in
`scf ▸ siteMoments  [t=95.93 s]` after a stage boundary at 60.17 s, the 35.8 s is the stage rebuild that
FINISHED at 95.93 — it says nothing about the cost of the site moments, which are **1.7 ms per call, 0.053 s
over the whole run** (`scf: order probe`, bucketed 2026-09-06 for exactly this reason).  Attributing a gap
to the item that closes it is the one wrong inference this instrument makes easy; the LEDGER is what prices
an item, the stamps only localise WHEN.

**READ THE STAMPS BEFORE THE LEDGER.**  The `timing` table says WHAT the run spent; the stamps say WHEN,
and the GAPS BETWEEN CONSECUTIVE STAMPS are exactly the time no bucket is charging.  That is how the last
25 s of unbucketed MnO time was located (a 35.8 s silent stretch between two anneal stages) without adding
a bucket first — cheaper than guessing where to put the next probe, and it works on a run you have already
taken.

…**including the ANNEALED driver**, which had none of the summary half until 2026-08-19 — and every MnO row
runs through it, so the "every GPW run reports PEAK RSS" claim above was false for precisely the runs whose
RAM this table most needs.  Two precision leaks were repaired at the same time: the `[fold]` line left
`cout` at 2 s.f. for the rest of the run (`std::defaultfloat` restores the format flag, not the precision),
and the fingerprint's `Efinal` inherited whatever a verbose table had left behind.  **A number's precision
must not depend on which other diagnostics were switched on** — every line above now states and restores
its own.

**Read the provenance, not just the number.**  `PEAK RSS` is the process high-water mark, so it is only a
clean per-config figure when the process runs ONE config — a `--gtest_filter` naming several tests reports
the watermark of the whole process.  `GPW_OMP_THREADS` governs the GPW **pair loops only** and says nothing
about the BLAS: with it unset, blaze still ran these rows at 115–239% CPU.  BLAS routing is
`QCHEM_BLAZE_BLAS` (default ON).  **qchem has no equivalent of CP2K's `GLOBAL| Number of threads` banner
line — a run cannot state how parallel it actually was**, which is why `scripts/bench` reports measured CPU%
and CPU seconds rather than trusting a knob.  Worth closing: one line at run start stating the effective
BLAS/OpenMP thread counts would make every future row self-describing.



```bash
# qchem  (GPW_REPORT=1 for the ledger; the energy line now prints 10 s.f.)
GPW_REPORT=1 scripts/bench "Si Gamma qchem"  -- build/Release/IntegrationTests/ITMain --gtest_filter=GPW_SCF.SiliconGammaConverges
# NaF: NAF_KMESH picks the mesh (1 = the CP2K-comparable Γ), NAF_SPAN the basis (sr2 default, sr = full)
GPW_REPORT=1 NAF_KMESH=1 NAF_SPAN=sr2 scripts/bench "NaF SR2 Gamma qchem" -- build/Release/IntegrationTests/ITMain \
    --gtest_filter=GPW_SCF.DISABLED_NaFRocksaltGamma --gtest_also_run_disabled_tests
# MnO, VA span, the runs-61/62 recipe -- now selected BY NAME instead of by overwriting a committed file.
# ONE ARM PER INVOCATION (MNO_SKIP_FM / MNO_SKIP_AFM): peak RSS is a PROCESS watermark, so a run that does
# both arms reports one number for the pair.
GPW_SPHERICAL=1 GPW_BASIS_SPAN=va MNO_ANNEAL="5e-3,0" MNO_ACC="Ladder,GDM" MNO_MOM=0 \
MNO_ORTHO_TOL=1e-3 MNO_SHARED_MU=1 MNO_IMPOSE=1 MNO_SKIP_FM=1 GPW_REPORT=1 \
    scripts/memsafe scripts/bench "MnO AFM2 VA qchem" -- build/Release/IntegrationTests/ITMain \
    --gtest_filter=GPW_SCF.DISABLED_MnO_AFM2_RhombohedralGamma --gtest_also_run_disabled_tests

# Si k-MESH ROWS -- ⚠ THESE WERE MISSING FROM THIS BLOCK UNTIL 2026-09-04, i.e. two table rows had no
# printed recipe at all (rule 3a's own failure mode).  Recovered by matching Etot to the banked value.
GPW_REPORT=1 scripts/bench "Si 222g qchem" -- build/Release/IntegrationTests/ITMain \
    --gtest_filter=GPW_SCF.DISABLED_SR_2x2x2GammaCentred_vs_CP2K --gtest_also_run_disabled_tests
GPW_REPORT=1 scripts/bench "Si 222shift qchem" -- build/Release/IntegrationTests/ITMain \
    --gtest_filter=GPW_SCF.SR_2x2x2ShiftedMP_vs_CP2K

# CP2K -- from the DECK'S directory (relative BASIS_SET_FILE_NAME) and at one thread
cd IntegrationTests/CP2K
OMP_NUM_THREADS=1 ../../scripts/bench "Si Gamma cp2k"    -- mpirun -np 1 cp2k.psmp -i si_fcc_gpw.inp
OMP_NUM_THREADS=1 ../../scripts/bench "MnO AFM2 VA cp2k" -- mpirun -np 1 cp2k.psmp -i mno_afm2_gpw_va.inp
```

Verify each side reproduces its own history before reading a Δ: the CP2K decks against `doc/Records/CP2Kresults.md`
(all five re-validated 2026-08-19, `doc/Records/CP2KBuild.md`), and the qchem runs against the tests' own anchors.


---

## 5a. ★★★ THE BIN-1 TABLE — per-ITERATION CPU, qchem vs CP2K, one thread each

**ONE TABLE, ON PURPOSE** (user, 2026-09-05: *"the important comparison table was the only table in
section 5a, but now the section is flooded with other tables so I cannot easily point to it"*).  Every A/B
delta, ledger split and open question that used to sit here is in **§5e**.  This is the bin-1 instrument
and the only table to quote for bin 1.

**THREAD STATE — SERIAL ON BOTH SIDES, AND MEASURED SO** (rule 3b): qchem `OMP_NUM_THREADS=1
GPW_OMP_THREADS=1`, measured **99% CPU** on every row; CP2K `OMP_NUM_THREADS=1`, measured 97–99%.
**Taken 2026-09-05**, one box (14 GB, 16 cores), qchem at `9f4f4ae2` built `-O3 -march=native`, CP2K 2025.2,
commands copied from §4.  ★ **Every row RE-TAKEN after §5f's lever A** (the Hartree energy stopped building
a matrix): the parity row fell 2.05× → 1.72×, the Si 8-k rows 0.84×/0.68× → 0.66×/0.56×, MnO defaults
0.90× → 0.82×, `BECKE_XC=0` 0.53× → 0.45×.  ⚠ **The qchem `setup` column is much larger than the 09-04 entry's ~57 s for MnO
because the Becke mesh build THREADS** (`src/Structure/Imp/UnitCell.C:273` — the partition loop is
`#pragma omp parallel for` over quadrature points, and `GPW_OMP_THREADS=1` pins it): 68.3 s serial per
build against the 16.7 s that ledger read with threads free.  The serial figure is the one that belongs
beside a serial CP2K row (`doc/Records/BenchmarkHistory.md` §9).

**HOW TO READ IT.**  `CPU` is whole-run user+sys.  `setup` is that run's own pre-SCF work — qchem: the sum
of the ledger's `setup:` buckets; CP2K: total CPU minus the sum of its printed per-step times.  So

- **`s/it (SCF)` = (CPU − setup) / iterations is the BIN-1 number**, and the `×` column is its ratio;
- **`setup` is the BIN-2 number**, in the same row, on the same run.

⇒ **Compare the SCF columns.**  The total column is kept only because it is what the whole-run table (§5)
divides — and the two disagree by 1.8× on MnO precisely because 44% of that run is setup.

⚠ **THE `q iters` / `c steps` COLUMNS BELOW ARE RULE-3f-STALE (2026-09-20)**: every qchem count was taken at
1e-3–1e-5 on the MIXER RESIDUAL (a G-space quantity), every CP2K count at `EPS_SCF` on max|ΔP| — not the same
measure, and on MnO ~100× looser on our side.  The **s/it** columns survive (a per-call cost does not care when
the loop stops); the counts and the `total ×` column do not — read them as "how many steps that recipe took to
ITS OWN criterion", never as a convergence-rate comparison.  The re-taken, like-for-like counts are the two
rule-3f MnO rows in §5 (footnote ⁸: 33 / 37 vs 44 / 104) and the Si example in rule 3f itself (5 vs 27 to reach
the same anchor).  A full re-take of this table on `Measure::MaxΔD` is queued in §10 below.

| row | span / k | q iters ⚠3f | q CPU | q setup | **q s/it (SCF)** | c steps ⚠3f | c CPU | c setup | **c s/it (SCF)** | **BIN 1 ×** | q s/it (total) | total × ⚠3f |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Si Γ ᵇ✓ | SIPP_SR, 1 k | 17 | 1.10 s | 0.16 s | **0.055** | 12 | 5.03 s | 0.4 s | 0.383 | **0.14×** ✅ | 0.065 | 0.15× |
| Si 2×2×2 Γ-centred ᵇ✓ | SIPP_SR, 8 k | 16 | 4.46 s | 0.35 s | **0.253** | 13 | 5.60 s | 0.6 s | 0.385 | **0.66×** ✅ | 0.275 | 0.64× |
| Si 2×2×2 shifted MP ᵇ✓ | SIPP_SR, 8 k | 14 | 4.26 s | 1.15 s | **0.218** | 14 | 5.91 s | 0.5 s | 0.386 | **0.56×** ✅ | 0.300 | 0.71× |
| NaF SR2 Γ ᵇ✓ | LOWQ_SR2, 1 k | 23 | 23.0 s | **9.82 s** | **0.573** | 16 | 7.18 s | 0.9 s | 0.394 | **1.45×** ⛔ | 1.000 | 2.23× |
| NaF full-SR Γ ᵇ✓ | LOWQ_SR, 1 k | 30 | 30.0 s | **11.17 s** | **0.628** | 27 | 101.85 s | 6.9 s | 3.519 | **0.18×** ✅ | 1.000 | 0.27× |
| **MnO AFM-II — ALL DEFAULTS** (5 of §2's 7 deviations active ᵃ) ᵇ⚠ | VA, 1 k | 14+17 = **31** | 395.6 s | **184.3 s** | **6.817** | 44 | 372.9 s | 8.1 s | 8.291 | **0.82×** ✅ | 12.76 | 1.51× |
| MnO AFM-II, `QCHEM_BECKE_XC=0` (4 of 7 — ONE rung down, NOT parity ᵃ) ᵇ⚠ | VA, 1 k | 14+25 = **39** | 147.0 s | 1.81 s | **3.723** | 44 | 372.9 s | 8.1 s | 8.291 | **0.45×** ✅ | 3.769 | 0.44× |
| ⚠ MnO AFM-II, `CP2K_COMPAT=1` probe (0 of 7 — `AT PARITY` as far as is KNOWN ᵃ) ᵇ⚠ | VA, 1 k | 10+10 = **20**, CAPPED | 287.5 s | 1.89 s | **14.28** | 44 | 372.9 s | 8.1 s | 8.291 | **1.72×** ⛔ | 14.38 | 1.70× |
| ★ **…the same row's FIXED-POINT stage alone — THE LIKE-FOR-LIKE NUMBER** ᵇ✓ | VA, 1 k | marginal | 37.4 s / 4 it | cancels | **9.35** | 44 | 372.9 s | 8.1 s | 8.291 | **1.13×** | 9.35 | 1.13× |
| MnO **FM** — ALL DEFAULTS ᵇ⚠ | VA, 1 k | 18+15 = **33** | 398.0 s | **184.8 s** | **6.460** | 22 | 192.4 s | 8.7 s | 8.350 | **0.77×** ✅ | 12.06 | 1.38× |

**★ THE RULE-3f RE-TAKE (2026-09-20, `scripts/retake5a` + the two `gpwprobe mno` rows of footnote ⁸) — the SAME
convergence measure (max|ΔD| / max|ΔP|) at the DECK'S threshold, the deck's loop shape (one density-side history:
Kerker G0=1 + Pulay 8; nothing on the Fock side; no MOM), `CP2K_COMPAT=1` (no Becke mesh, no imposition), both
codes serial.**  Whole-run wall and RSS here are like-for-like BECAUSE the loops now stop at the same place;
CP2K's numbers are the banked ones (its counts were already at `EPS_SCF`).  The one stated difference: the Si
decks mix P directly (α 0.4) where our stable no-Fock-side route is Kerker + Pulay (§10 row "Linear D-mixing
DIVERGES").

| row | EPS / MAX_SCF | **q iters** | q wall | q RSS | **c steps** | c wall | c RSS | iters q/c | qchem E | CP2K E |
|---|---|---|---|---|---|---|---|---|---|---|
| Si Γ (free) | 1e-7 / 60 | **16** | 2.7 s | 39 MB | 12 | 5.0 s | 148 MB | 1.33× | −7.115067447 | −7.115057882 |
| Si 2×2×2 Γ-centred (free) | 1e-7 / 60 | **15** | 8.6 s | 54 MB | 13 | 5.6 s | 153 MB | 1.15× | −7.778472674 | −7.778457865 |
| Si 2×2×2 shifted MP | 1e-7 / 60 | **14** | 10.6 s | 55 MB | 14 | 5.9 s | 153 MB | 1.00× | −7.867452508 | −7.867436530 |
| NaF SR2 Γ | 1e-6 / 200 | **16** | 7.8 s | 84 MB | 16 | 7.2 s | 173 MB | 1.00× | −24.43039482 | −24.431213375 |
| MnO AFM-II (free) | 1e-6 / 200 | **33** | 5m42s | 262 MB | 44 | 6m14s | 217 MB | 0.75× | −61.41154311 | −61.303325178 |
| MnO AFM-II +U 4 eV (free) | 1e-6 / 200 | **37** | 6m23s | 263 MB | 104 | 15m46s | 217 MB | 0.36× | −60.81349836 | −60.68597088 |

⇒ **Iteration counts are AT PARITY on Si/NaF (1.0–1.3×) and in our favour on MnO (0.75× / 0.36×)**; the Si Γ
1.33× is 16 vs 12 on a 3-second run.  ⚠ **The Si 2×2×2 shifted-MP energy is −7.867452508 here, 16 µHa from
CP2K's −7.867436530, where the table above has −7.868473429 (−1.04 mHa).  Checked the same day: it is the
IMPOSITION, not the criterion or the grid** — with `QCHEM_IMPOSE_SYMMETRY=1` every criterion and both XC grids
give −7.868473; free, −7.867452.  The symmetry fold of the shifted (k=±¼, non-TRIM) mesh lowers the energy by
1.02 mHa: an OPEN DEFECT, `doc/OpenWork.md` §3 (footnote ¹ is now history).  The other energies moved by
< 1 µHa; on MnO the 108 mHa absolute offset (§4a) is unchanged.

ᵇ **✓ = BOTH CODES RUN THE SAME KIND OF SCF STEP ON THIS ROW; ⚠ = THEY DO NOT.**  A per-iteration ratio is
only a comparison when the iteration is the same thing on both sides, so this was CHECKED per row, in each
run's own trace, rather than assumed:

| row | qchem's step | CP2K's step | |
|---|---|---|---|
| Si Γ, Si 2×2×2 (both) | DIIS + linear mixing on a diagonalise | `DIIS/Diag.` + `P_Mix` | ✓ |
| NaF SR2 Γ, NaF full-SR Γ | DIIS + Kerker on a diagonalise — the ladder's GDM rung is **NOT ENGAGEABLE** (Fermi-smeared occupations sit outside the integer-occupation manifold GDM rotates, and it says so) | `Broy./Diag.` | ✓ |
| the four MnO rows | **two stages**: `Ladder` (fixed point) then `accel: GDM` — DIRECT MINIMISATION with a geodesic line search | `Broy./Diag.` | ⚠ |

⇒ **Only the MnO rows are affected, and only through their second stage.**  CP2K's counterpart for a direct
minimiser is OT — and (checked in the decks and in every log's update-method column) **none of the
benchmarked decks run `&OT` either**: they all diagonalise and mix.  The one deck with an `&OT` section is
the unbenchmarked `naf_gpw.inp`, and its comment says why (diagonalisation diverged on that diffuse basis).
So on these rows the mismatch is not GDM-vs-OT, it is **minimiser vs mixer**, and stage 2's line-search
trial densities land in the probe's per-iteration average with nothing on the other side to match them.

The ★ row isolates the comparable half by DIFFERENCING `GPW_MNO_NMAX` (a single-stage `MNO_ANNEAL=5e-3` run
at N=6 against N=2, so setup and seed cancel and what is left is the marginal cost of one iteration):
**3 gathers + 2 collocations per iteration against CP2K's 2 + 2**, and at our per-call rate that is 9.35 s
against 8.29.  ⇒ **On the same algorithm we are 1.13×, and the whole residual is the third gather** —
\f$V_H\f$ gathered separately from \f$v_{xc}^\sigma\f$ (§5f lever B).

ᵃ **counted off each run's own banner**, not from memory (rule 3b).  Defaults:
`QCHEM_DM_LOWRANK=on* GPW_STREAM_FOLD=on* QCHEM_MIX_RHO_M=off GPW_XC_DM_SOURCE=off
QCHEM_IMPOSE_SYMMETRY=on* QCHEM_BECKE_XC=on* GPW_DAWARE_SCREEN=on*` (`*` = differs from CP2K, five of
them); the middle row is the same with `QCHEM_BECKE_XC=off(stated)`; the bottom row prints
`CP2K_COMPAT=1 -> AT PARITY` with every flag off.

⇒ **PER SCF ITERATION WE ARE AHEAD OF CP2K ON SEVEN OF THE NINE ROWS**, including MnO on its default
route (0.82×) and its FM arm (0.77×).  Two rows are behind, and they are behind for two different reasons:

⛔ **NaF SR2 Γ — 1.45×, and it is a Becke cost, not a GPW one.**  43% of that row's CPU is setup.  §5e.
✅ Its DRIVER is comparable (footnote ᵇ): both codes diagonalise and mix there — the ladder's GDM rung is
not engageable under Fermi smearing — so unlike the MnO rows, this 1.45× is a like-for-like step cost and
nothing about it is waiting on OT.

⛔ **THE PARITY ROW, AND WHAT IT ACTUALLY MEASURES.**  `CP2K_COMPAT=1` removes the stream fold (5.2× on
MnO's pair count) and the low-rank ρ as well, so this is the honest algorithm-to-algorithm row and the
other MnO rows are not.  It reads **1.72×** — but the two-stage probe runs a MINIMISER (GDM) that CP2K's
decks do not, so read the two rows separately (footnote ᵇ):
- **on the same algorithm — fixed point, diagonalise and mix, both codes — we are 1.13×**, and the entire
  residual is the third gather (§5f lever B, blocked behind N4);
- the **1.72×** adds our stage-2 GDM line search, whose trial densities have no counterpart in a CP2K run
  that is mixing rather than minimising.  ⇒ It is not a like-for-like step cost, and chasing it as one
  would be optimising against a comparison that does not exist.
⚠ Two further caveats, both real: `GPW_MNO_NMAX=10` CAPS the probe at 20 iterations with a different stage
mix, so rule 3d applies; and **`CP2K_COMPAT=1` is our best KNOWN parity, not proven parity** — §2's
deviation list has grown every time anyone has looked (4 → 7 items).
⇒ **BIN 1 IS NOT CLOSED, but it is now one gather wide on the comparable algorithm.**

★★★ **AND THE TABLE SIZES BIN 2, MEASURED.**  MnO's serial setup is **184.3 s — 47% of the default run —
against CP2K's 8.1 s (23×)**, and **136.6 s of it is TWO Becke mesh builds** (68.3 s each, one per anneal
stage; the rest is 43.1 s of XC-mesh Φ tables).  Both vanish with `QCHEM_BECKE_XC=0`, where our setup is
**1.76 s and BEATS CP2K's 8.1 s**.  ⇒ Bin 2 is a Becke-mesh question, exclusively.

⛔ **BUT IT IS A GRID-SIZE QUESTION BEFORE IT IS A CODE QUESTION** (user, 2026-09-06): *"I am reluctant to
work on that until identify the proper becke grid size (Nradial=40, Nangular=29 is big) that yields the
same accuracy … as the uniform 20³ mesh.  Becke is always more expensive to setup, no way around that.  It
is parallel which helps."*  The build scales with the point count, and `nR=40, degree=29` sets the entire
Becke side of every cost figure on this page — so the recipe, not the loop, is what to attack first.
⇒ **The calibration is the gate**, and the harness for it exists (`BeckeLadder` + `UniformXCProbe` in
`IntegrationTests/GPW_SCF_UT.C`, scoring \f$\Delta E_{xc}\f$ and \f$\max|\Delta V_{xc}|\f$ on one frozen
density); what is missing is putting the uniform probe on the ladder so the crossover can be read off.
`doc/OpenWork.md`'s next action carries it, including the two standing rules the V2.6a work earned (a
frozen density understates the self-consistent shift on a METAL; Al's angular convergence is
non-monotonic).
⚠ **Two things survive the gate whatever it says**: the build is `omp parallel for`
(`src/Structure/Imp/UnitCell.C:273`, ~4× — the user's *"it is parallel which helps"*), and **there are two
builds for the same cell**, one per anneal stage, which is a redundancy independent of grid size.

⚠ **WHAT MOVED, AND WHY — rule 3a's check, honestly reported.**  **Iteration counts**: removing the
gather's density screen (`a7561e92`) is not bit-identical, so the small rows took a different SCF path —
Si Γ 11 → 17, Si 2×2×2 Γ-centred 7 → 16, shifted MP 16 → 14, NaF SR2 29 → 23.  **Every MnO row still
reproduces its banked iteration count exactly** (31 / 39 / 20 / 18+15).
**Energies**: they now differ in their last digits from the banked values, by two roundoff-scale effects
stacked — the unscreened gather, and lever A's G-space \f$E_H\f$ (§5f), which sums the same integral in a
different order.  Si Γ reproduces to all 10 printed digits; the other Si rows to 2–7 nHa; NaF to 2.7–6.3
µHa; **MnO to 2.2–2.6e-7 Ha** (−61.40297529 against the banked −61.40297551), i.e. 4e-9 relative on a 61 Ha
total and four orders below the 99.65 mHa the row is actually measuring.  ⚠ Anchors pinned tighter than
1e-7 on a periodic total energy will need re-banking; none in the suite are (806/806 green).
The MnO **FM** row is a different case: it had not been re-taken since 08-19 (footnote ⁵), so this is its
first measurement on the current code — 398 s CPU and 481 MB against the stale row's 2321 s and 4947 MB.

---

## 5. THE ROWS — single thread

⚠ **THE WHOLE-RUN TABLE BELOW IS NOT THE BIN-1 INSTRUMENT — §5a ABOVE IS** (user, 2026-09-04: *"the
previous table in section 5 is not very informative"*).  Its `wall` and `CPU` columns are whole-run totals
over runs whose iteration counts differ by 3×, and they charge each run's setup to its per-step cost.
**Read it for two things: the ENERGY column (Δ vs CP2K, the accuracy claim) and peak RSS (bin 3).**
For runtime go to §5a — same runs, decomposed.

Energies in Ha.  **Both columns measured on this box (14 GB, 16 cores) through `scripts/bench`** — the CP2K
side is no longer banked prose: `apt`'s CP2K 2025.2 reproduces every banked 2026.1 deck value to the printed
digits (`doc/Records/CP2KBuild.md`), so both codes are measured under one wrapper.  Provenance per row is in
*How each row was produced*.

> ✅ **BOTH SIDES ARE NOW SERIAL AND MEASURED SO** (2026-09-05): qchem `OMP_NUM_THREADS=1
> GPW_OMP_THREADS=1` at 99% CPU on every re-taken row, CP2K `OMP_NUM_THREADS=1` at 97–99%.  Wall and CPU
> therefore agree to ~1% and the `CPU ×` column is finally a like-for-like ratio.
> ⚠ **The history this replaces, because it will be quoted**: rows before 09-05 had `GPW_OMP_THREADS` unset
> while blaze threaded the BLAS regardless — measured 115–239% CPU — so their wall column FLATTERED qchem
> by the number of cores it happened to take.  **A thread-count knob is not a measurement**; an earlier cut
> of this table printed "1 thr" on the qchem rows on the strength of the unset knob.

**★ EVERY ROW BUT TWO RE-TAKEN 2026-09-05, SERIAL, ON ONE BINARY** (`a3045163`, `-march=native`) — the same
runs §5a decomposes, so the two tables cannot disagree.  The `taken` column says which rows are not from
that cut.

| row | last taken | why it may have moved since |
|---|---|---|
| Si Γ, Si 2×2×2 (both), NaF SR2 Γ, NaF full-SR Γ, MnO AFM-II defaults, MnO `BECKE_XC=0`, **MnO FM** | **09-05, serial, current** | — |
| NaF SR2 2×2×2 | 08-28 | no CP2K counterpart (footnote ²), so it is not a head-to-head row; ⚠ cheap to re-take |
| **MnO `CP2K_COMPAT=1`** (the uncapped 93-step row) | **08-28** | ~45 min to re-take; the 09-05 cut used the documented `GPW_MNO_NMAX=10` probe instead (§5a, §5b).  The old row's *"does not converge"* verdict is RETRACTED — see ⁷ |

⇒ **The `CP2K_COMPAT=1` row is now real rather than aspirational**, because the XC pair route made the
parity ROUTE affordable — see footnote ⁷ (§5d).  CP2K column untouched throughout.

| system | k-mesh | span | qchem Etot | CP2K Etot | Δ (qchem−CP2K) | wall q / c | **CPU q / c** | **CPU ×** | peak RSS q / c |
|---|---|---|---|---|---|---|---|---|---|
| Si (FCC) | Γ | SIPP_SR | −7.115067844 | −7.115057882 | **−10.0 µHa** | **1.1 s** / 5.2 s | **1.1** / 5.0 s | **0.22×** | **30** / 148 MB |
| Si (FCC) | 2×2×2 Γ-centred | SIPP_SR | −7.778472833 | −7.778457865 | **−15.0 µHa** | 4.5 s / 5.8 s | **4.4** / 5.6 s | **0.79×** | **35** / 153 MB |
| Si (FCC) | 2×2×2 shifted MP | SIPP_SR | −7.868473429 ¹ | −7.867436530 | **−1.04 mHa** | 4.3 s / 6.1 s | **4.2** / 6.0 s | **0.71×** | **34** / 153 MB |
| NaF (rocksalt) | Γ | LOWQ_SR2 (both) | −24.4303428161 | −24.431213375 | **+0.871 mHa** | 23.2 s / 7.4 s | **23.0** / 7.2 s | **3.2×** | **55** / 173 MB |
| NaF (rocksalt) | 2×2×2 Γ-centred | LOWQ_SR2 | −24.5468834873 | — ² | | 1m07.8s / — | 84.5 s / — | — | **68** / — MB |
| NaF (rocksalt) | Γ | LOWQ_SR (full) | −24.4309445921 | −24.432293467 | **+1.349 mHa** | **30.2 s** / 1m42s | **30.0** / 102 s | **0.29×** | **58** / 186 MB |
| **MnO AFM-II — ALL DEFAULTS** ⁴ | Γ | **VA (N=118)** | −61.40297529 ⁴ | −61.303325178 | **−99.65 mHa** | **6m38s** / 6m14s | **396** / 373 s | **1.06×** | **476** / 217 MB |
| **MnO AFM-II, `QCHEM_BECKE_XC=0`** ⁶ | Γ | **VA (N=118)** | −61.40358753 | −61.303325178 | −100.26 mHa | **2m28s** / 6m14s | **147** / 373 s | **0.39×** | **113** / 217 MB |
| **MnO FM — ALL DEFAULTS** ⁵ | Γ | **VA (N=118)** | −61.44158219 ⁵ | −61.304782531 | **−136.80 mHa** | **6m40s** / 3m13s | **398** / 192 s | **2.07×** | **481** / 217 MB |
| **MnO AFM-II, `CP2K_COMPAT=1`** ⁷ | Γ | **VA (N=118)** | −61.39789688 ⁷ | −61.303325178 | −94.57 mHa | 45m40s / 6m14s | **2736** / 373 s | **7.3×** | **112** / 217 MB |
| **MnO AFM-II, `CP2K_COMPAT=1`, CP2K's LOOP + MEASURE** ⁸ | Γ | **VA (N=118)** | −61.41154311 (free, 33 it to max\|ΔD\|<1e-6) | −61.303325178 (44 it) | −108.2 mHa | **5m42s** / 6m14s | **10.3 s/it** / 8.29 s/it | **1.24×** | 262 / 217 MB |
| **MnO AFM-II +U (4 eV, Mn d), `CP2K_COMPAT=1`, CP2K's LOOP + MEASURE** ⁸ | Γ | **VA (N=118)** | −60.81349836 (free, 37 it) | −60.68597088 (104 it) | −127.5 mHa (ΔE(U) +0.5980 vs +0.6174) | **6m23s** / 15m46s | **10.3 s/it** / 9.05 s/it | **1.14×** | 263 / 217 MB |
| MnO AFM-II | 2×2×2 (`MNO_KMESH=2`) | VA | ❓ | ❓ | | ❓ | ❓ | | ❓ |


### 5d. Footnotes to the table

Compact here; the full stories are in `doc/Records/BenchmarkHistory.md` at the section named after each.

- **¹** Si 2×2×2 shifted MP (−1.04 mHa) — ⚠ **LOCALISED 2026-09-20 by the rule-3f re-take: it is the IMPOSITION.**
  The FREE run is −7.867452508, 16 µHa from CP2K; imposing the point group on the shifted (k=±¼) mesh gives
  −7.868473429 under every criterion and both XC grids.  Open defect, `doc/OpenWork.md` §3.  What follows is
  the history.  The residual after a **D-aware integrate-back SCREEN defect** was
  fixed 2026-08-19.  It is the suite's ONLY fractional-k SCF coverage (every other k is TRIM, where the
  defect is structurally invisible), which is why it had rotted to −3.7351 while DISABLED.  *(history: the
  screen was reading \f$\mathrm{Re}[D\overline{e^{ikR}}]\f$ where it needed \f$|D|\f$.)*
- **²** No CP2K counterpart: CP2K's NaF decks carry no `&KPOINTS` section, i.e. they are Γ, while the qchem
  test defaults to 2×2×2.
- **⁴** MnO ALL DEFAULTS — **re-measured four times** (08-19/20, 09-04, 09-05), same command, same box;
  every cut agrees on the energy to 4e-8 and on the iteration count (14+17 = 31) exactly.
- **⁵** ✅ **RE-TAKEN 2026-09-05** (it had been the only pre-Step-2 measurement left): serial, VA recipe,
  18+15 = 33 iterations, \f$E_{tot}\f$ −61.44158232 against the 08-19 row's −61.441583060 (7e-7, the
  gather-screen removal).  Its cost columns fell 5.6× in CPU and 10× in RAM; per SCF iteration it is
  **0.83× CP2K** (§5a).
- **⁶** ⚠ **ONE of §2's seven deviations off, NOT parity** (user, 2026-09-05) — it makes XC run CP2K's way
  and leaves everything else ours.  It is the first MnO row on which qchem beat CP2K on both axes and is
  still the standout: **0.53× per SCF iteration**, 107 MB against 217 MB, and a setup of 1.76 s against
  CP2K's 8.1 s.
- **⁸** The two CP2K-LOOP rows (2026-09-20, rule 3f's first rows): `gpwprobe mno` under `CP2K_COMPAT=1 GPW_REPORT=1 GPW_SPHERICAL=1 GPW_BASIS_SPAN=va MNO_MOM=0 MNO_ACC=Null MNO_PULAY=8 MNO_PULAY_START=5 MNO_MEASURE=maxdd MNO_EPS=1e-6 MNO_SKIP_FM=1` (+ `MNO_U=4 QCHEM_U_EIGEN=0` for the +U row; logs `probe_u0_compat.log`, `probe_u4_compat.log`) — the deck's loop shape (ONE density-side history, diagonalise, no Fock extrapolation, no MOM) on the deck's measure at the deck's threshold, both codes serial.  Both FREE runs CONVERGE and hold AFM-II (m̃(q)Ω/2 = 3.13 / 3.18 e).  The per-iteration residual is the banked one: **3 gathers + 2 collocations against CP2K's 2 + 2** (102 gathers / 33 it at 1.89 s, 68 collocations at 1.93 s = 9.7 of the 10.3 s), i.e. §5f lever B, V_H gathered separately from v_xc.  Setup 1.9 s vs 8.1 s.  The +U refresh costs ~15 ms/iter (`scf: eager refresh`), CP2K's full-matrix S½PS½ ~0.75 s/iter (8.29 → 9.05).  ⚠ The earlier attempt on the anchor's recipe (Ladder DIIS + MOM + `PulayDepth=8`, log `mno_u5.log`) sat at max|ΔD| 3e-3 (U=0) / a period-2 cycle (U=4) for 200 iterations — two histories extrapolating each other and MOM rotating a smeared degenerate frontier; it is what rule 3f's "5 vs 27 on Si" looked like on MnO.  The energies: the free run is 0.3 mHa below the imposed anchor (−61.41154 vs −61.41124: the star-averaged density is a constrained problem), and the 108 mHa absolute offset to CP2K is the open §4a row, unchanged by any of this.
- **⁷** The parity row: an earlier *"at true parity our recipe does not converge"* verdict was **RETRACTED**
  (08-28).  It still hits the iteration cap, but for a far more benign reason — 5.1 mHa short, not 3.8 Ha,
  with the AFM order surviving both stages.
- **⁸** The 93-iteration parity figure is pre-fix and its stage mix differs — **do not** read 29.4 → 19.3 as
  a 1.5×; per-iteration is not comparable across iteration caps (rule 3d).

### 5b. The three MnO rows are ONE system and ONE recipe

★ **THE THREE `MnO AFM-II` ROWS ARE ONE SYSTEM AND ONE RECIPE, with progressively more of OUR deviations
switched off.**  Read them as a ladder, not as three experiments — the only thing changing is which of §2's
**seven** declared deviations are on.  ⚠ `QCHEM_BECKE_XC=0` is ONE of the seven: it is a rung on the way to
parity, **not** parity (user, 2026-09-05: *"there is much more to parity than just BECKE_XC=0"*).

| the row | what is off | CPU (serial, 09-05) | RSS |
|---|---|---|---|
| **ALL DEFAULTS** | nothing — every deviation at its qchem default (Becke XC, imposition, low-rank ρ, stream fold) | 396 s (31 it) | 476 MB |
| **`QCHEM_BECKE_XC=0`** | ONE of the seven: the Becke XC mesh — so XC runs CP2K's way, everything else still ours | **147 s (39 it)** | **113 MB** |
| **`CP2K_COMPAT=1`** | ALL SEVEN KNOWN — the honest algorithm-to-algorithm row | 288 s (20 it, CAPPED probe) | 132 MB |

⚠ The middle row is the fastest because Becke's atom-centred quadrature is expensive machinery — and on the
default route it is 44% of the run before SCF even starts (§5a).  The bottom row costs the most per step
because parity also removes the STREAM FOLD (5.2× on MnO's pair count) and the low-rank ρ; its 341 s is a
`GPW_MNO_NMAX=10` probe, not a converged run (the uncapped row was 2736 s over 93 capped steps).
⇒ **Always say WHICH of the three** — "the MnO row" has meant all three of these in the space of one day.
⚠ And **`CP2K_COMPAT=1` is our best KNOWN parity, not proven parity**: §2's list started at four items and
is at seven, having grown every time anyone looked (user, 2026-09-05: *"CP2K_COMPAT=1 is supposed to be
parity, but we keep finding new differences"*).  A new find moves the bottom row, in either direction.


### 5c. Like-for-like: what "same span" costs


CP2K can be forced down to the spherical spans qchem holds at FULL RANK — `IntegrationTests/CP2K/`
`mno_{afm2,fm}_gpw_v{a,b}.inp` with the `VALENCE-LOWQ-V{A,B}` entries in `VALENCE-LOWQ-BASIS`; **VA = N 118**
(full rank in both codes), **VB = N 128**.  So the MnO comparison does NOT wait on the 136-span question
(`OpenWork` Step 6).

On the qchem side that span used to be produced by `doc/scripts/bisect_valence_sph.py` **overwriting the
committed `BasisSetData/valence_lowq_sph.bsd` in the working tree** — so a run could not state which span it
had used, and no row over it was reproducible after the file was restored.  VA and VB are now committed
basis sets (`valence_lowq_v{a,b}.bsd`, `BasisSetData::VALENCE_LOWQ_V{A,B}`) selected by
**`GPW_BASIS_SPAN=va|vb|sph|sr`**, and the run prints its own `nFunctions` — reproducing run 61's basis block
exactly (118 functions, λ_min 1.29e-3, cond 4.41e3).  Si and NaF already share one span per material.


---

## 5e. Bin-1 working notes — the deltas and the open questions BEHIND §5a's table

⚠ **§5a holds ONE table on purpose.**  Everything that used to sit beside it and made it hard to point at
lives here.  Nothing in this section is the bin-1 instrument; §5a is.

### What removing the gather density screen cost and bought (2026-09-04)

★★★ **THE UNSCREENED GATHER (2026-09-04) — A REAL SPLIT BY k-MESH, NOT A UNIFORM WIN.**  Removing the
density screen from the integrate-back (see §6, and §10/OpenWork_History4 for why it is also a CORRECTNESS fix)
unlocks the k-independent memo, so one real-space sweep serves every k-block.  Where there are no k-blocks
to share it across, it is pure added width:

| GATHER seconds PER ITERATION (ledger bucket ÷ iterations) | screened | **unscreened** | |
|---|---|---|---|
| Si Γ (1 k) | 0.00693 | 0.00784 | **+13%** ⛔ |
| Si 2×2×2 Γ-centred (**8 k**) | 0.196 | **0.112** | **−43%** ✅ |
| Si 2×2×2 shifted MP (**8 k**) | 0.187 | **0.111** | **−41%** ✅ |

and on the MnO rows — all Γ, and all with IDENTICAL iteration counts before and after, so total CPU is a
fair comparison there:

| MnO, VA recipe, serial, NATIVE | iters | screened | **unscreened** | |
|---|---|---|---|---|
| ALL DEFAULTS | 31 | 386.1 s | 408.6 s | **+5.8%** ⛔ |
| `QCHEM_BECKE_XC=0` | 39 | 155.7 s | 168.2 s | **+8.1%** ⛔ |
| `CP2K_COMPAT=1` probe | 20 | 338.9 s | 346.1 s | **+2.1%** ⛔ |

⇒ **IT SAVES ~40% OF GATHER TIME ON AN 8-k RUN AND COSTS 2–8% AT Γ.**  Every MnO row is Γ, so the flagship
numbers get slightly worse while the multi-k rows get materially better.  **KEPT**: it is first a
CORRECTNESS fix — the diagonal-seed defect and the stream fold's orbit-invariance, both in
§10 below — and the Γ cost is the price of not having a self-fulfilling truncation in the Fock.
⚠ It is NOT a `CP2K_COMPAT` matter: the screen is gone on every route, not just the parity one.

### ⚠ THE BECKE MESH BUILD THREADS — which collides with §7a's reading (open, 2026-09-05)

§5a's setup column is SERIAL, and 68.3 s of it per build is the Becke mesh.  That loop is
`#pragma omp parallel for` over quadrature points (`src/Structure/Imp/UnitCell.C:273`), so it is NOT
structurally serial — yet §7a's NaF threading run inferred a ~9 s serial residual and matched it to the
setup buckets, of which 7.0 s IS the mesh build.  Both cannot be the whole story.
⇒ **Resolve it while profiling bin 2** (§10 below, named next action): measure the mesh build's own
speed-up curve on MnO and on NaF separately.  Candidate explanation — NaF is a 2-atom cell whose point
count is too small to amortise the region, so it threads on MnO and does not on NaF.

★ **A note on units, since the two columns are subtracted**: the ledger's buckets are WALL seconds and
`CPU` is user+sys.  Every §5a row measured 97–99% CPU, so the subtraction is sound there; it would NOT be
on a threaded row, which is one more reason §5a is a serial-only table.

### The k-SCALING gap — the one structural bin-1 defect the table shows

Same system, same span, only the k-mesh changing (SCF-only seconds per iteration, §5a's numbers):

| | qchem s/iter | CP2K s/iter |
|---|---|---|
| Γ | **0.051** | 0.383 |
| 8 k-points (Γ-centred) | **0.323** | 0.385 |
| **rise** | **6.3×** | **1.01×** |

We start **7.5× ahead** of CP2K at Γ and keep only **1.2×** of it at 8 k — still a win, but the whole
shape of the small-cell story is in that rise.  Cause and the two refuted fixes:
§10 k-scaling row (the gather memo is bypassed whenever a density screen is passed, and 32
of 51 gathers on the 8-k run are the SAME FIELD).  ⚠ Part of their flatness is not an optimisation to copy
but a fold we do not do: CP2K folds a shifted MP mesh 8→4 by time reversal (§6b item 1).

### NaF's loss is a BECKE cost, not a GPW one

`NaF SR2 Γ` is the one small row where we are behind per iteration, and its ledger says why —
Becke mesh build **7.00 s (30% of the run)** + XC-mesh ρ sampling **6.55 s (28%)** + hamiltonian ctor
**2.02 s (9%)**, against a box walk (gather + collocate) of **0.80 s (3%)**.
⇒ On the system where Becke should be at its best, **Becke IS the runtime**, and 42% of that row's CPU is
pre-SCF setup — i.e. it is mostly a bin-2 row wearing a bin-1 label.  That makes §6's parked
\f$\lVert V_{xc}-V_{xc}^{fit}\rVert\f$ study the deciding measurement for this row too.

### The box walk is ~95% of an MnO parity run

109 s of 114.6 s, of which the GATHER is 74% (81.5 s over 65 calls, against 27.9 s over 22).
⇒ Amdahl leaves nothing outside the walk worth touching on that route, and any per-call win is worth ~3×
more on the gather side than on the collocate side.

### STANDING PROBE for bin 1

(user, 2026-08-28: *"we just cut off at ~10 or so iterations, just to get a decent average"*):
`CP2K_COMPAT=1 GPW_MNO_NMAX=10` — ~6 minutes, and quote per-call beside per-iteration.

### ✅ THE PER-STEP COMPARISON — ANSWERED 2026-09-05, see §5f (kept for the question it asked)

★ **RESOLVED BY COUNTING**: a fixed-point iteration issues **4 gathers + 2 collocations**, and both of the
readings below turned out to be true of a different part of it — the term assembly asks twice for Hartree
(§5f lever A) *and* the GDM line search does build extra densities once it engages (lever C, and it is
bin-4 currency).  CP2K's own counters say it does 2 and 2.  The original entry, which framed the question:

The 09-04 ledger shows **~9 distinct KS-field integrations per SCF iteration** on a 2-channel system where
the physics needs ~3 (one \f$V_H\f$ gather + one \f$V_{xc}\f$ per channel); `GPW_INTEGRATE_CENSUS=1` says
all of them are genuinely NEW fields, so it is NOT redundancy a memo could remove.  Two readings, and they
call for opposite responses:
- **the accelerator is buying bin 4 with bin 1 work** — if the Ladder line-search builds \f$H\f$ several
  times per "iteration", then our iteration is not CP2K's step and per-iteration is the WRONG metric; the
  honest one is total \f$H\f$-builds (or wall) to convergence;
- **or the term assembly asks more often than it needs to**, in which case it is a ~3× on the dominant
  bucket and would put the run BELOW CP2K.

⇒ **Distinguish them before optimising either way**: count \f$H\f$-builds per SCF step directly, and read
CP2K's own per-step \f$H\f$ count out of its log.  Cheap, and it decides whether bin 1 is finished.

---

## 5f. ★★★ WHERE THE 2.05× IS: **CALL COUNT, NOT KERNEL** (2026-09-05)

**Our kernel is already FASTER than CP2K's, per call.**  Read the two codes' own counters side by side —
CP2K's `T I M I N G` block (44 SCF steps) against our ledger (the parity probe, 20 iterations):

| | calls per SCF step | s per call | who wins the call |
|---|---|---|---|
| CP2K `integrate_v_rspace` | **2.00** (88/44) | 2.133 | |
| qchem gather | **4.0** (marginal, fixed-point) | **1.882** | **qchem, 0.88×** |
| CP2K `calculate_rho_elec` | **2.05** (90/44) | 2.029 | |
| qchem collocate | **2.0** (marginal, fixed-point) | **1.889** | **qchem, 0.93×** |

Both codes spend ~99% of the run in these two routines (CP2K 370 s of 374; qchem 334 s of 341).
⇒ **The whole 2.05× is that we call the gather TWICE as often as CP2K.**  Nothing in the kernel is behind.

**HOW THE COUNTS WERE TAKEN — differencing, not instrumentation** (`GPW_MNO_NMAX` = n against n+4, so the
setup/seed offset cancels and what is left is the marginal cost of one iteration):

| stage | gathers / iteration | collocations / iteration |
|---|---|---|
| Ladder (fixed point) | **4.0** | **2.0** |
| GDM (direct min) | 5.0 | 2.0, plus the LINE SEARCH's trials once it engages |

★ **AND THE LEDGER NAMES THE FOUR.**  A fixed-point iteration closes 2 `h ball` fields and 2 `h raw` fields.
On this route the raw adjoint is XC's (one field per spin) and the ball adjoint is Hartree's, so:

> **THE HARTREE MATRIX IS BUILT TWICE PER ITERATION** — once for the FOCK at \f$\rho_{mix}\f$
> (`Vee_Hartree::MakeMatrixT`) and once for the ENERGY at \f$\rho_{new}\f$
> (`Vee_Hartree::GetEnergy` → `0.5*cd->DM_Contract(this,cd)`).  Two different densities, so the gather memo
> cannot catch it.  ✅ Corroborated: with a pass-through mixer (\f$\alpha=1\f$, no Kerker) the two densities
> coincide on some iterations and the count falls 4.0 → 3.5.

**⇒ THREE LEVERS WERE ON THE TABLE.  ONE LANDED, ONE IS REFUTED, ONE REMAINS** (2026-09-05).  The
per-iteration budget was 4×1.882 + 2×1.889 = **11.3 s** against CP2K's 2×2.133 + 2×2.029 = **8.3 s**:

| | lever | verdict |
|---|---|---|
| **A** | **\f$E_H\f$ WITHOUT A MATRIX** — \f$E_H=\tfrac12\mathrm{Tr}(DV_H)=\tfrac12\Omega\sum_{\Delta G}\overline{\tilde\rho}V_H\f$, a G-space pairing on the fit ball.  Not an approximation: the gather is the exact adjoint of the collocation, so Parseval makes the two expressions equal by construction.  ★ **AND ONE MAP IS ENOUGH** — \f$V_H=k\tilde\rho\f$ with \f$k\f$ real, so \f$\overline{\tilde\rho}V_H=|V_H|^2/k\f$ and \f$\tilde\rho\f$ is never fetched (a second fetch costs a second IBZ star-average of the whole map: 16 ms/call on Si Γ) | ✅ **LANDED `9f4f4ae2`.  Parity probe 341.4 → 287.5 s CPU (−15.8%), gathers 95 → 66, 2.05× → 1.72×.**  Si 2×2×2 rows −19.5% / −11.4% (the same fix memoizes \f$V_H\f$ across irrep blocks).  806/806, every iteration count and \f$E_{tot}\f$ reproduced |
| **B** | **ONE GATHER PER SPIN: \f$\langle i|V_H+v_{xc}^\sigma|j\rangle\f$** — CP2K's `sum_up_and_integrate`, and the `CompositeExFunctional` argument one level up | ⛔ **REFUTED ON MEASUREMENT — the two fields are not on the same discretization.**  See below |
| **C** | the GDM LINE SEARCH's trial densities — **42 of the probe's 82 collocations**, i.e. ~28% of the run after A | ⏸ **PARKED AS A TODO UNTIL WE HAVE OT** (user, 2026-09-06).  See below |

⛔ **WHY B IS REFUTED (2026-09-05).**  Summing the two potentials needs them on ONE representation, and they
are deliberately on two: **\f$V_H\f$ is a BALL field** (band-limited, Poisson in G, gathered by the ball
adjoint's per-level \f$\{G\}\f$ restriction) while **\f$v_{xc}\f$ is a RAW RASTER field** (pointwise
nonlinear, gathered by the raw adjoint's per-level spectral BOX truncation).  Measured directly — build the
Hartree block both ways on Si Γ — they differ by **6e-5 relative**, so routing \f$V_H\f$ through the raw
adjoint is a different operator, not a re-association: it breaks the adjointness with the collocation that
made \f$\tilde\rho\f$, which is exactly the property lever A's Parseval identity rests on.
⇒ **B is available only if BOTH fields ride the ball route** — i.e. only if XC gives up the raw
\f$\rho_{DM}\ge0\f$ feed, which is the collapse-basin fix (`doc/GPWPlan` 0.5(f2)).  ⇒ **It becomes
available if N4 lands** (`doc/OpenWork.md`: make everything robust with \f$V_{xc}[\rho\ge0]\f$), and not
before.  Filed there, not here.

⏸ **WHY C IS PARKED, NOT PURSUED (user, 2026-09-06).**  *"CP2K uses the OT method instead of GDM so we
can't do proper parity timings against CP2K anyway."*  The line search is DIRECT MINIMISATION; CP2K's
counterpart for that is OT, which we do not have — and (checked in the decks and logs, §5a footnote ᵇ) the
benchmarked CP2K runs are not even using OT: they diagonalise and mix.  ⇒ There is no comparable step to
optimise C against.  **Trial-density cost becomes a real question when OT exists and can be timed against
GDM, minimiser to minimiser; until then it is a TODO, not a lever.**

⇒ **WHERE THAT LEAVES BIN 1.**  Measure the stage CP2K's decks actually match — our FIXED-POINT stage —
and it is **3 gathers + 2 collocations per iteration against CP2K's 2 + 2, i.e. 9.35 s against 8.29 =
1.13×** (§5a's fixed-point row, differenced 2026-09-06).  The term assembly is at its floor: 1 Hartree +
2 XC, and CP2K reaches 2 only by summing \f$V_H\f$ into \f$v_{xc}^\sigma\f$ before integrating.
⇒ **The whole remaining like-for-like gap is lever B, and lever B is behind N4.**

---

## 6. qchem ACCELERATIONS NOT IN CP2K — the qchem-vs-qchem deltas

These are what §2's knobs BUY.  They are excluded from a head-to-head row by rule 3c and reported here
instead, which is the more useful statement.

| acceleration | measured worth | provenance |
|---|---|---|
| Becke XC mesh (`QCHEM_BECKE_XC`) | ⚠ **NEGATIVE on MnO, and worse than it looked**: 1.8× on the SCF iteration (6.82 vs 3.72 s/iter) and **102× on setup** (184.3 vs 1.81 s — two mesh builds at 68.3 s each).  It buys atom-centred accuracy, not speed | §5a |
| symmetry imposition | buys CONVERGENCE, not accuracy and not the magnetic basin — the AFM order survives a free run | history §2 |
| stream fold (`GPW_STREAM_FOLD`) | 5.2× on MnO's pair COUNT | history §7 |
| D-aware screen (`GPW_DAWARE_SCREEN`) | +14.5% wall on the collocation route | doc/OldPlans/ScreeningPlan.md §6 |
| collocation memo depth 5 | 3.74× on its bucket (60 → 16 collocations) | history §4 |
| gather memo on (V, screen) | catches the exact duplicates; the remainder are genuinely distinct fields | history §5 |
| `-march=native` (now the Release DEFAULT) | −9.6% / −8.4% / −5.3% CPU on the three MnO rows, \f$E_{tot}\f$ bit-identical | CMakeLists.txt |

⚠ **NOT an acceleration and NOT a deviation**: `GPW_CONTRACT_CUBE` (the separable-contraction collocation
kernel).  CP2K collocates exactly this way (`grid_cpu_collint.h`, Mathieu's three 2-D tables).  It is on
this page only so a row can state which kernel produced it — every run prints `[collocation] kernel=…`.

---

## 6b. ACCELERATIONS **CP2K** HAS THAT **WE** DO NOT — ⚠ ALSO AN EMERGING LIST

§2 and §6 are one direction; this is the other, and it had no home until 2026-09-04.  A row is only fair
if BOTH lists are known, and this one started at zero because nobody had looked.

| # | what CP2K does | what it costs us | found |
|---|---|---|---|
| 1 | **TIME-REVERSAL k FOLDING** — a Monkhorst-Pack mesh is folded \f$8\to4\f$ (`BRILLOUIN\| List of Kpoints ... 4`, weights 0.25, with `K-Point point group symmetrization OFF`).  We run all 8. | up to **2×** on any non-TRIM mesh.  ⚠ It folds NOTHING on a Γ-centred 2×2×2, where every k is its own inverse — which is why the two Si k-rows behave differently | 09-04 |
| 3 | **`sum_up_and_integrate` — ONE integrate per spin.**  CP2K adds \f$V_H\f$ and \f$v_{xc}^\sigma\f$ into one real-space potential and integrates it ONCE per spin: 88 `integrate_v_rspace` calls over 44 steps.  We call the gather 4× per iteration (2 Hartree + 2 XC).  ⚠ This is the whole of the parity row's 2.05×, since our per-call gather is 0.88× theirs | the parity row (§5f) | 09-05 |
| 2 | **k-INDEPENDENT per-step cost** — their per-iteration time is flat from Γ to 8 k (0.383 → 0.385 s SCF-only) where ours rises **6.3×** (0.051 → 0.323).  Not a "feature" so much as a consequence of collocating the k-summed density once per step | most of the Γ lead (§5a, §5e) | 09-04 |

⚠ Item 1 is a REAL algorithmic advantage they hold, not a deviation to switch off; item 2 is a gap of ours
with a diagnosed cause.  ⇒ **Neither belongs in `CP2K_COMPAT`** — that switch turns OUR accelerations off,
and turning theirs off is not available to us.

⚠ CP2K also runs `Wavefunction type COMPLEX` on BOTH Si k-rows, including the all-TRIM Γ-centred one where
we use REAL blocks.  So on that row we hold an advantage they decline to take — and are still only at
parity per iteration.

---

## 7. THE ROWS — threaded (12 cores), AND THE OMP GAP ANALYSIS

✅ **FILLED 2026-09-06, BOTH SIDES.**  The cross-code half was the thing missing since 09-04; it is here
now, and it inverts the question that motivated it.  Same binaries and decks as §5a, 12 cores each,
\f$E_{tot}\f$ agreeing to all printed digits with the serial runs on every row.

### 7a. Our side: 12 threads against our own serial

| row | serial wall | **12-thread wall** | **speedup** | 12t CPU / serial | RSS ser / 12t |
|---|---|---|---|---|---|
| NaF SR2 Γ | 23.0 s | **10.96 s** | **2.10×** | 1.41× | 55 / 73 MB |
| MnO `QCHEM_BECKE_XC=0` | 147.0 s | **49.5 s** | **2.97×** | 1.50× | 113 / 197 MB |
| MnO `CP2K_COMPAT=1` probe | 287.5 s | **61.6 s** | **4.67×** | 1.58× | 132 / 210 MB |
| MnO ALL DEFAULTS | 395.6 s | **128.4 s** | **3.08×** | 1.51× | 476 / 563 MB |

★ **The parity route threads BEST (4.67×)** — it is nothing but the box walk, which is our OpenMP region.
The Becke routes carry more serial residue, which 7c breaks down.

### 7b. ★★★ CP2K's side — AND THE HEADLINE IS THAT ITS PARALLELISM DOES NOT HELP HERE

| CP2K row | serial wall | 12 OMP threads | 12 MPI ranks | best speedup |
|---|---|---|---|---|
| Si Γ | 5.18 s | 5.49 s (**0.94×**) | — | **none** |
| NaF SR2 Γ | 7.37 s | 8.95 s (**0.82×**) | — | **none** |
| MnO AFM-II VA | 374.2 s | 342.7 s (**1.09×**) | **259.6 s (1.44×)** | **1.44×** |

⚠ **These are CP2K's own thread counts, not a knob we set and hoped**: its banner reports
`GLOBAL| Number of threads for this process 12`, and the run still measured **196% CPU** — ~2 cores of
work spread over 12 threads.  On the two small cells it is *slower* than serial while billing 2× the CPU.
The MPI arm is the honest one for CP2K (its design centre), and 12 ranks buy **1.44× for 8.3× the CPU**
(3107 s against 374 s serial).

★ **AND ITS OWN TIMING BLOCK SAYS WHERE**: on MnO, the two GPW hot routines — which are **98% of the run**
— barely move between 1 and 12 threads:

| CP2K routine (SELF time) | serial | 12 threads | speedup |
|---|---|---|---|
| `grid_collocate_task_list` | 182.6 s | 167.9 s | **1.09×** |
| `grid_integrate_task_list` | 187.7 s | 167.2 s | **1.12×** |

⇒ **It is precisely the GPW collocate/integrate route whose OMP axis does not scale here** — the same two
routines whose per-call cost we beat by 0.88×/0.93× serially (§5f).  Two readings, and this measurement
does not separate them: (a) the route's OpenMP was never a priority — CP2K's research centre is hundreds
of small molecules in a box, where the parallel axis that matters is MPI over molecules/atoms, not OMP
inside a solid's task list (user, 2026-09-06); (b) a 4-atom cell simply has too few tasks per level to
spread.  ▶ **The discriminator is one deck we do not have**: the same run on a supercell (say 2×2×2 MnO,
32 atoms).  If its grid routines start scaling, it was size; if they do not, it was the code.  ⚠ Until
then, quote the measurement, not either explanation.

⇒ **CROSS-CODE, AT 12 CORES, ON MnO** — per SCF step, since the iteration counts differ (rule 3d):
**qchem 4.14 s/iteration against CP2K's best 5.90 s (0.70×)**, and we get there on **596 s of CPU against
their 3107 s (0.19×)**.  ⚠ Rule 3c applies: that qchem row runs our accelerations.  ⚠ And the standing §2
caveat is the whole explanation — this is a 4-atom, high-symmetry cell, the regime that most favours a
symmetry-exploiting code and least favours CP2K's distribute-over-atoms design.  **On a 100-water box the
verdict would invert**, and none of this says anything about many-node scaling.

⇒ **SO: "ARE WE MISSING PARALLEL OPPORTUNITIES CP2K EXPLOITS?"  NOT ON THIS CLASS OF SYSTEM.**  We
out-scale CP2K on both of its axes here.  The opportunities we are missing are our OWN, and 7c names them.

### 7c. Where OUR threaded time actually goes — bucket by bucket

MnO ALL DEFAULTS, the same run in both columns (serial 395.6 s → 12 threads 128.4 s):

| bucket | serial | 12 threads | speedup | |
|---|---|---|---|---|
| setup: Becke mesh build | 138.2 s | 16.9 s | **8.2×** | ✅ |
| scf: XC-mesh ρ sampling (matrix-free) | 74.9 s | 11.3 s | **6.6×** | ✅ |
| scf: collocate (pair scatter) | 44.8 s | 8.3 s | **5.4×** | ✅ |
| setup: XC-mesh Φ tables | 44.3 s | 6.8 s | **6.5×** | ✅ |
| scf: integrate-back (pair gather) | 16.4 s | 2.7 s | **6.2×** | ✅ |
| **scf: XC-mesh quadrature H_xc** | 11.2 s | **9.2 s** | **1.21×** | ⛔ |
| **scf: E_H V_H field build** | 7.6 s | **8.5 s** | **0.89×** | ⛔ |
| scf/setup: the FFT closures, local-PP short | ~3 s | ~3.4 s | ~0.9× | ⛔ |
| **everything not in a bucket** | ~52 s | **~59 s** | **~0.9×** | ⛔ |

⚠ **THE LAST ROW IS SUPERSEDED — IT IS NOW 0.03 s** (2026-09-06, `doc/Records/ParallelAndOraclePlan.md` 1.1 + 1.1(a)).
Instrumenting it did not find "diagonalisation, orthogonalisation, mixing, the fit solves": the SCF's whole
linear-algebra side is ~5 s of a serial run and the diagonalisation is **0.038 s**.  It found **the
Hamiltonian being built once per ANNEAL STAGE** — `SolidCalculation::BuildStage` re-`Factory`s it and had no
bucket at all.  On the same 12-thread MnO row (121.0 s wall, 39 buckets summing to 120.98 s):

| bucket | GPW_OMP unset | 12 threads | speedup | |
|---|---|---|---|---|
| `setup: hamiltonian ctor` (exclusive of the two below) | 51.4 s `[×2, 25.7/call]` | **46.4 s** `[×2, 23.2/call]` | **1.11×** | ⛔ |
| ⤷ `setup: XC-mesh Φ tables` | 54.1 s `[×2]` | 6.2 s `[×2]` | **8.7×** | ✅ |
| ⤷ `setup: becke mesh build` | 15.4 s `[×2]` | 15.8 s `[×2]` | **1.0×** | ⚠ see below |
| **Hamiltonian construction, all in** | **121.0 s = 39%** | **68.4 s = 56%** | 1.77× | ⛔ |
| the three residue buckets together | 0.36 s | 0.36 s | — | — |
| **everything not in a bucket** | **0.05 s** | **0.03 s** | — | ✅ |

Same build, same recipe, `Etot=-61.40297529` on both arms to all printed digits; 313.9 s wall / 532 s CPU /
**169%** unset, 121.0 s / 569.5 s CPU / 469% at 12 threads.  39 buckets summing to 313.71 of 313.76 s and
120.98 of 121.01 s respectively.

★ **SO THE SHAPE OF THIS ROW IS SETUP, NOT SCF — 69.5 s of setup against 51.5 s of SCF at 12 threads.**
That is bin 2 ("Init (pre-iteration) time") in the gap-close priority order, and it is now the biggest
single lever on the MnO wall — bigger than either pair loop.  **The half that does not thread is the ctor's
own 1.11×**; its children scale fine, which is precisely why it hid.  ▶ Next: bucket the fit-basis half of
that 23.2 s/call, and establish whether an anneal stage must rebuild the Hamiltonian at all (~34 s if not).

★★★ **AND THEN THE CTOR WAS OPENED UP, AND THE ROW MOVED — 2026-09-06, `ParallelAndOraclePlan.md` 1.1(b).**
Six more buckets found that the 23.75 s/call was **two calls to one bad index**: the run folds the same
~97k mesh points twice per Hamiltonian (the orbit-consistency filter in `UnitCell::CreateIntegrationMesh`,
then `FoldMesh` in `GPW_IBS::CreateXCQuadrature`), 9.6 + 9.5 s/call.  `TorusIndex` bucketed EVERY mesh on a
constant 64³ grid — 0.37 points per bucket on average, which is the wrong statistic for a clustered
atom-centred radial mesh.  Grid made as fine as the tolerance allows + centre-bucket-first probe:

| | GPW_OMP unset | 12 threads |
|---|---|---|
| orbit-consistency fold | 9.61 → **0.163** s/call | 9.61 → **0.193** s/call |
| `FoldMesh` | 9.47 → **0.151** s/call | 9.47 → **0.190** s/call |
| Hamiltonian construction, all in | — | 34.4 → **15.5** s/call |
| **WHOLE RUN** | 313.9 → **237.2 s** (**1.32×**) | 120.8 → **83.2 s** (**1.45×**) |

⚠ The "GPW_OMP unset" column is a MIXED arm, not a serial one (see the protocol note below) — both sides
of its 1.32× were taken the same way, so the ratio is sound, but do not read it as a serial figure.
| setup / SCF split | 71.8 / 165.2 s | 32.1 / 51.1 s |
| peak RSS | 496 → 506 MB | 562 → 561 MB |

`Etot = -61.40297529` on every arm before and after, to all printed digits — it is a data structure, not a
numerical method.  814/814 green.  The largest setup bucket is now the Becke mesh build (8.15 s/call), then
the site-adapted angular sets (4.03 s/call, never measured before).

⚠ **THE MnO ROWS IN §5a AND §7a ARE STALE, AND ARE NOT RE-TAKEN PER INCREMENT** (user, 2026-09-06).  A row
is a CROSS-CODE artefact: it costs a CP2K arm as well as ours, and re-taking one after every speedup both
burns the box and lets each increment mask the last.  Same discipline as the anchor-moving sprint
(`doc/OpenWork.md` item **S**): **let the wins accumulate, then re-bank the rows ONCE**, at the end of
Phase 1 — after 1.2 (BLAS-mode arm) and 1.3 (the 2×6), which move them again.  Until then read the ratios
in §5a/§7a as a floor on our side and quote 1.1(b)'s numbers for anything that turns on the MnO wall.

✅ **THE BECKE MESH BUILD "COLLISION" IS RESOLVED, AND THE BANKED ROW WAS RIGHT — MY ARM WAS WRONG**
(2026-09-06).  The flag raised here said the mesh build measured the same in both arms (15.4 vs 15.8 s)
where §7c banks 138.2 → 16.9 s (8.2×).  One run at `GPW_OMP_THREADS=1` settles it:

| arm | CPU% | `setup: becke mesh build` |
|---|---|---|
| `GPW_OMP_THREADS=1` | **99%** | **68.74 s/call** |
| `GPW_OMP_THREADS` unset | 188% | 7.10 s/call |
| `GPW_OMP_THREADS=12` | 469% | 8.15 s/call |

⇒ The mesh build threads at **8.4×**, §7c's 8.2× stands, and §5e's open question is closed.

★★★ **THE REAL FINDING IS A PROTOCOL DEFECT: `GPW_OMP_THREADS` UNSET IS NOT A SERIAL ARM.**  The Becke
build is *"parallel by DEFAULT under QCHEM_OPENMP; `GPW_OMP_THREADS`, when set, is honoured as the thread
CAP"* (`UnitCell.C`) — so UNSET means **all cores** for that loop while the GPW pair loops stay serial.
An "unset" row is therefore a MIXED arm, not a serial one, and its 169–188% CPU says so out loud.  This is
§4's own *"a knob is not a measurement"* rule biting the person who wrote it down.
▶ **A serial row MUST set `GPW_OMP_THREADS=1`** (99% CPU is the check).  Any row in this file whose serial
arm reads much above 100% CPU was taken the mixed way and should be re-taken at the Phase-1 re-bank.

**MnO ALL DEFAULTS, post-1.1(b), `Etot=-61.40297529` on all three arms:**

| arm | wall | vs serial |
|---|---|---|
| `GPW_OMP_THREADS=1` (true serial) | **361.2 s** | 1.00× |
| unset (mixed) | 237.2 s | 1.52× |
| `=12` | **83.2 s** | **4.34×** |

★ **THE BECKE MESH BUILD THREADS AT 8.2× — §5e's open question is CLOSED, and the answer is "it threads
fine on MnO".**  Whatever holds NaF back (7d) is not the loop being serial.

⛔ **THE THREE THINGS THAT DO NOT SCALE, in priority order:**
1. **~59 s of UNBUCKETED work** — the largest single block of a threaded run, and it is invisible because
   nothing times it: diagonalisation, orthogonalisation, mixing, the fit solves, the SCF bookkeeping.
   ▶ **The first action is an instrument, not an optimisation**: bucket the SCF's non-GPW work.
2. **The XC-mesh quadrature \f$H_{xc}=\Phi^\dagger\,\mathrm{diag}(w\,v)\,\Phi\f$ at 1.21×** — and this
   one is a DELIBERATE TRADE, not an oversight (user, 2026-09-06; the reasoning is written into
   `DeltaFit_IBS::AdjointT`).  **Our parallelism lives ABOVE the linear algebra** — per k-block / irrep /
   spin — with BLAS pinned to one thread (`qchem::PinBlasToOneThread`), precisely to avoid OMP nesting.
   The measurement behind it: one dispatched whole-matrix `zgemm` runs **34.1 GFlop/s against 1.87 for
   ANY blocked or viewed form**, so hand-blocking to spread over threads loses 13× to gain 8×.
   ⇒ **The bucket's 1.21× is the PRICE of that policy, and the thing to question is the WIDTH OF THE LEVEL
   ABOVE, not the policy.**  On MnO at Γ that level is 2 spins × 1 k-block = **2-way**, so ten of twelve
   cores have nothing to do in this bucket by construction.  On the 8-k Si rows the same policy has 16-way
   width and nothing is left on the table.
   ▶ **AND THERE IS A CHEAP TEST WORTH RUNNING**: the "one dispatched zgemm" arm needs
   `QCHEM_BLAZE_BLAS=ON`, which was correctly refused on 2026-08-15 because the system BLAS was then
   NETLIB.  ⚠ It is not any more — `libblas`/`liblapack` now resolve to `openblas-pthread` — so the 34
   GFlop/s path is available again, **single-threaded, composing with the pin and with the levels above**.
   That is a rebuild and one run.
3. **The \f$V_H\f$ field build at 0.89×** — that bucket is 6.6% of the threaded wall and it is mostly
   `SymmetrizeGMap`, the IBZ star-average over the point group (48 ops), which is a serial walk over a
   \f$\{G\}\f$ map.  It got LOUDER, not quieter, when §5f's lever A removed the gathers around it.

### 7d. The Amdahl residual, and the NaF anomaly that is still open

Solve \f$S+P/12=t_{12}\f$ against \f$S+P=t_1\f$:

| row | Amdahl \f$S\f$ | as % of serial | what 7c attributes it to |
|---|---|---|---|
| MnO ALL DEFAULTS | **104 s** | 26% | ⚠ **re-attributed 2026-09-06**: the ~59 s "unbucketed" is the **Hamiltonian ctor's 46 s exclusive half** (×2, once per anneal stage) + ~21 s of non-scaling buckets + imbalance |
| MnO `BECKE_XC=0` | 32 s | 22% | same shape, no mesh build |
| NaF SR2 Γ | 9 s | 40% | ⚠ **matches its SERIAL setup buckets (9.8 s) to ~7%** |

⚠ **NaF REMAINS THE ANOMALY.**  Its inferred serial fraction still equals its setup, yet the very same
setup buckets thread at 6–8× on MnO.  The likeliest reading is now SIZE, not structure: NaF is a 2-atom
cell whose mesh build is 7 s and whose Φ tables are 0.35 s — too little work per thread to amortise the
regions.  ▶ Cheap to settle: the same bucket table as 7c, taken on NaF.

⚠ **CPU inflation is 1.41–1.58×** on our side even with the barrier spin gone — thread management, load
imbalance, memory bandwidth.  CP2K's is 1.8× threaded and **8.3× under MPI**, so this is not a gap against
them; it is a cost of the last speedup increment.

---

## 8. WHAT IS STILL MISSING

★ **AS OF 2026-09-04, in the user's priority order** (`doc/OpenWork.md` carries the same list as the
forward queue):
0. ✅ **THREADED ROWS — DONE 2026-09-06, both sides (§7).**  CP2K at 12 threads was the missing piece; it
   turns out its parallelism does not help on these cells (0.82–1.09× on OMP, 1.44× on 12 MPI ranks), so
   the remaining threading work is entirely OUR OWN: ~59 s of unbucketed serial SCF work, a dense XC GEMM
   on blaze's serial kernels, and the IBZ star-average.
1. **CLOSE BIN 1 AND BIN 2 SINGLE-THREADED** — per-iteration CPU and pre-SCF setup, against CP2K.  ✅ §5a
   is FILLED as of 2026-09-05 (nine rows, both codes serial, setup split out): bin 1 is ahead on seven of
   nine and comes down to the `CP2K_COMPAT=1` row at 2.05×; bin 2 is the Becke mesh build.  The whole-run
   table is NOT informative (mixed thread states, 3× different iteration counts) and is kept only for
   \f$E_{tot}\f$ and RSS.
2. ✅ **THE BUSY-WAIT BARRIER** — fixed, and everything re-run at 12 threads (§7).  Cause was identified 2026-09-04:
   nothing in the tree sets `KMP_BLOCKTIME` or `OMP_WAIT_POLICY`, and LLVM's libomp spins **200 ms** after
   every parallel region — with per-shell-pair regions that is mostly spin.  It matches the inflation
   already recorded: **663 s threaded CPU against 500 s serial for the same work**.
3. ✅ **The MnO FM row** (footnote ⁵) — re-taken 2026-09-05: 411 s CPU / 474 MB, 0.83× CP2K per SCF step.
4. ✅ **NaF's CP2K iteration counts** — read off its own logs (SR2 16 steps, full-SR 27); both rows now have
   a per-iteration column in §5a.
5. **An attribution A/B** — the 09-04 deltas are "current vs banked" across everything since 08-28, with
   no parent-commit A/B on the VA recipe.

⏸ **PARKED, deliberately** (user, 2026-09-04: *"defocusing"*): the
\f$\lVert V_{xc}-V_{xc}^{fit}\rVert\f$ fit-quality study.  It stays worth doing one day — every
Becke-vs-uniform cost number here is taken at UNKNOWN-EQUAL accuracy, so §6's "Becke is a negative
acceleration" is a COST statement and not a verdict — but it is not what bins 1 and 2 need.




- **★ REPEAT THE WHOLE TABLE AT 12 THREADS** (user, 2026-08-19).  This cut is the SERIAL baseline — CP2K
  genuinely serial, qchem serial-except-BLAS — which is the right reference for "how much work does each
  code do", and the wrong one for "how fast is it on this box".  Both sides scale from the environment:
  CP2K's `psmp` takes `OMP_NUM_THREADS`, qchem takes `GPW_OMP_THREADS` for the pair loops plus whatever
  blaze does with the BLAS.  Run every row at `OMP_NUM_THREADS=12` / `GPW_OMP_THREADS=12` and keep BOTH
  cuts: the serial one is the algorithmic comparison, the threaded one is the user-facing time, and the
  ratio between them is qchem's parallel efficiency — which nothing currently measures.
  **Deliberately deferred until the qchem times come down** (user): re-measuring an 20-minute row now, when
  Steps 2–3 are about to change it by an order of magnitude, buys a number with a short shelf life.
- **The multi-k MnO row** (`MNO_KMESH=2`, ALL-TRIM so the whole k-sum runs real) — cost unmeasured, and CP2K
  needs a matching `&KPOINTS` deck.  At 20 min for the Γ row, budget for it before starting.
- **A k-point CP2K deck for NaF**, so the 2×2×2 row becomes a comparison instead of a qchem-only datapoint.
- **The Si shifted-MP row**, which needs the fractional-k regression fixed first (footnote ¹).
- **An A/B on NaF's XC route.**  Today's direct runs agree with CP2K to +0.88 mHa (SR2) and +1.35 mHa (SR),
  whereas July's numbers were 0.19 / 0.10 mHa — but those came from the GRID-CONTINUATION test
  (`DISABLED_NaFGridContinuation`, coarse→fine seeding, uniform-XC era) and not from this one, so it is not
  a like-for-like drift.  One A/B would say whether ~1 mHa is NaF's honest agreement at today's defaults or
  the price of a recipe change; until then do not quote NaF as "0.1 mHa class".


---

## 9. STANDING OBSERVATIONS


- **Si Γ agrees to 1.1e-5 Ha** at 1.5× CP2K's CPU time and 1.8× its RAM.  On the small cells the gap is
  modest; the diffuse and magnetic cells are where it opens up.
- **NaF Γ agrees to 0.9 mHa (SR2) / 1.3 mHa (full SR)** at 13.2× / 2.1× the CPU time and 3.3× / 16.6× the
  RAM.  The full span costs us 3090 MB where CP2K needs 186 MB — the same RAM signature as MnO, and the
  reason the diffuse cells are where this table bites.  (The SR2 row's 13.2× against the full span's 2.1× is
  not a typo and is not yet explained: the SMALLER basis is the worse ratio.  Both are Γ, same cell, same
  code path — worth one look, because whatever it is scales the wrong way.)
- **MnO is where both columns hurt — but Steps 2 and 3 took a large bite out of both (2026-08-19/20).**  The
  AFM row went 6.0× → **4.2×** CPU, 23× → **6.2×** RAM (4947 → 1350 MB against CP2K's 217 MB) and 3.2× →
  **1.38×** wall on an identical 118-function span; the FM row still carries the pre-Step-2 numbers.  **The
  RAM half was the streams** (Step 2), so Step 4 was answered by Step 2 rather than by a campaign of its
  own, exactly as predicted.  CP2K caches nothing and its kernels are just fast — what is left on our side
  is the Φ-shaped ρ GEMM and the Becke partition, in that order for wall and the reverse for CPU.
- **MnO carries a ~100 mHa configuration-BLIND offset and a −37.15 mHa configuration-SELECTIVE one** on that
  identical span (`OpenWork` Step 5) — and qchem sits BELOW a variational reference there, which convicts an
  operator or convention rather than a basis.  The run now prints its own term breakdown
  (`Ekin/Een/Eee/Exc/Enn/E_alphaZ`), which is where Step 5's term-by-term comparison starts.
- **The FREE production MnO run still folds NOTHING** — `[fold]` prints `NONE` at all three sites on a cell
  whose magnetic group has 12–24 ops (`OpenWork` Step 2).  Any MnO runtime row taken from a free run is
  measuring an unfolded run, and must say so.  Arming T3 does not change this: it is a `MNO_IMPOSE`
  decision, i.e. physics, and the 1.5× wall / 3.7× RAM measured above is what a free run is paying for the
  freedom.
---

## 10. STANDING MEASUREMENTS / OPEN PERF LEVERS (moved from `OpenWork.md` §4b, 2026-10-01)

Verbatim archive of the old tracker: `doc/Records/OpenWork_History5.md`.  Each row: what is open · next concrete action · record.

**Where it stands (numbers from §5a and the rule-3f table):** per SCF iteration serial qchem is AHEAD of CP2K on seven of nine §5a rows; the like-for-like parity row (`CP2K_COMPAT=1`, fixed-point stage) is **1.13×**, 1.72× for the capped two-stage probe, and its whole residual is ONE gather (lever B, behind N4).  Rule-3f whole-run, same measure/threshold/loop: MnO AFM-II **5m42s vs 6m14s wall** (33 vs 44 iterations); Si 2×2×2 8.6 s vs 5.6 s.  Peak RAM on the parity routes 113–132 MB vs CP2K's 217 MB (whole-run MnO 262 vs 217 MB).  Threaded (12 cores) 4.14 s/step vs CP2K 5.90 (they run 0.82–1.09× their own serial).  History: the MnO per-iteration cost started at **67× CP2K's** (573 s vs 8.5 s, `BenchmarkHistory.md` §cache; the earlier "100×" was an instrument artefact).  No QE timing exists.

**Parked levers:** B (one gather per spin) behind N4; C (GDM trial densities) behind OT (OpenWork §2).

| row | what is open | next concrete action · record |
|---|---|---|
| **Linear D-mixing DIVERGES on free Si with no Fock accelerator** (found 2026-09-20 by `scripts/retake5a`) | `GPW_Si.Γ_CP2K` under `CP2K_COMPAT=1 GPW_ACC=null` (plain linear D-mixing, `ProductionGates` α=0.30): descends to −7.11505 by iteration 19, then diverges (E −7.0886 at 60, [F,D] 2e-3 → 0.2).  The adaptive relax raised α 0.30 → 0.45 at iteration 2 and never re-damped.  CP2K's `DIRECT_P_MIXING` at α=0.4 converges the same cell to 1e-7 in 12 steps; our Kerker G0=1 alone takes 54, Kerker + Pulay(8) 13 (E −7.115067447, the anchor).  DIIS on the Fock side had been masking it | two questions: (a) why does the adaptive controller (V1.18) not re-damp on a rising energy — a rising-E, rising-[F,D] trajectory is exactly its trigger; (b) is direct P mixing at fixed α=0.4 stable on our loop (a `GPW_ALPHA` knob + adaptive off)?  Cheap: seconds per run.  Not a physics question; until answered the no-Fock-side route on our side is Kerker + Pulay · `doc/Benchmark.md` rule 3f, `scripts/retake5a` |
| **Re-take Benchmark §5a on CP2K's convergence measure** (rule 3f, 2026-09-20) | every `q iters` / `c steps` / `total ×` cell in §5a compares a 1e-3–1e-5 mixer residual with CP2K's `EPS_SCF` on max\|ΔP\| — two different questions; the s/it columns survive, the counts and totals do not | re-run each §5a row with `Δρmeasure=MaxΔD`, `MinΔρ=EPS_SCF`, `NMaxIter=MAX_SCF` off its deck, the deck-shaped loop (Null accelerator, Pulay 8, no MOM) where the deck mixes Broyden/Pulay, and `CP2K_COMPAT=1`; one sitting, ~1 h serial (Si rows seconds, NaF minutes, MnO ×2 the bulk).  The two MnO rows are DONE (§5 ⁸); the Si Γ example is in rule 3f · `doc/Benchmark.md` rule 3f |
| **Size the Becke grid** (was item 1 / bin 2) | MnO's setup is 184 s = 47% of the default run (CP2K 8.1 s), 136.6 s of it TWO Becke mesh builds — but the RECIPE (nR=40, degree 29) is **over-generous on every system**: 3.5× on Si/Al, ~25× on NaF/Mn (`BeckeLadder`, scored by \f$\max|\Delta V_{xc}(i,j)|\f$ — score the MATRIX, the energy cancels error the operator does not) | ▶ **a POLICY CALL for the user**: at an absolute \f$\max|\Delta V_{xc}|\le\f$ 1e-4 the answer is nR=40, GL-17…21 (2–3× cheaper than production).  ⚠ no default flips on ladder evidence alone — Al is non-monotonic and a frozen ladder understates the self-consistent shift on a metal ⇒ a converged A/B on Al first.  Two items survive whatever the calibration says: (a) why TWO builds for one cell (each anneal stage builds one); (b) the build threads at 8.2× on MnO — NaF's serial fraction is still its setup (likeliest SIZE).  Then the per-element radial scaling, the coarse-end routing calibration, Becke+IBZ (the real-space star-average is untested on this route) · History4 "1. BIN 2" |
| **BM(3) — the stock Lebedev rule is already site-invariant** | Lebedev-29 (302 dirs) is EXACTLY invariant under Si's \f$T_d\f$ site group (0 unmatched of 7248) because Lebedev rules are octahedral orbits; W2b nonetheless builds a site-adapted rule at 886 dirs/atom — a **2.9× sitting unclaimed** | test the STOCK rule for site invariance first and reuse it when it passes (the test must stay EXACT); keep W2b as the fallback; measure on MnO before spending anything — it vanishes as site symmetry drops · History4 row BM |
| **The ρ̃-mixed sampling bucket** (the actual per-iteration XC lever) | `FourierMixCD.C:65` samples \f$\rho(r)=\sum_G\tilde\rho(G)e^{iGr}\f$ by DIRECT SUMMATION at every mesh point, every Kerker/Pulay iteration: 35.0 s / 6 iterations serial against the DM GEMM's 1.70 s.  ⚠ It is NOT the low-rank D bucket (that GEMM is nearly bypassed on ρ̃-mixed recipes — the 7–8× rank win is real but bites only on DM-backed routes) and NOT Φ-sparsity (Φ is 48% dense on MnO; per-atom batching = 1.05×) | attribute the cost INSIDE that bucket before optimising: an FFT to a coarse uniform grid + interpolation to the mesh (O(N log N + npts)) is legitimate only for the cusp-free CORRECTION, not the full ρ; an adaptive G-ball keyed to \f$|\tilde\delta|\f$ is cheap late and exact at convergence.  N4 decides what is sampled · History4 "Vxc MUST BE FED THE DM ρ(r)" |
| **The k-scaling gap** (cross-k gather memo) | Si per-iteration cost rises 6.9× from Γ to 8 k-points; CP2K's rises 1.03× — we start 6.7× ahead and spend it on k.  Cause: the per-offset reductions \f$B_{ij}(n)\f$ are k-INDEPENDENT but `IntegrateMemo` is bypassed whenever a density screen is passed (32 of 51 gathers on Si 2×2×2 are the SAME FIELD).  ⛔ two attempts refuted (History4): withholding the screen perturbs the SEED trajectory; flooring vanishing weights breaks the fold's orbit invariance | the trajectory-exact route: memoize \f$B\f$ over the UNION of active sets and let each block contract only its own; raise `kMaxIntegrateMemos` (4) with it.  Pin 19 already removed the gather's D-screen, which is the other half.  CP2K also folds 8 → 4 by TIME REVERSAL on the shifted mesh — a separate 2× · History4 "THE k-SCALING GAP" |
| **Step 2's remainder — the per-iteration G-space folds** | the {G}-star fold is wired at two STATIC sites; the per-iteration consumers (ρ̃, the Poisson multiply, the V_xc gathers, G_ERI3 columns, seed structure factors) are UNFOLDED: 12–24× on MnO's magnetic group, 48× cubic, unclaimed; **T3.4b** multi-k per-block arming of the pair-stream fold (Γ-only is armed by default) | extend the ball fold to the per-iteration sites (the FFT itself does not fold trivially); T3.4b = union-of-reps stream caches or the star-summed joint scatter.  Space-group collocation reduction route (b) IMPOSES the symmetry — say so · History4 "Step 2", `doc/Records/SymmetryUpgradePlan.md` |
| **\f$B_{ij}(R)\f$ k-independent 1E memo** | "keep k out of the key"; payoff only on multi-k | time with row KP (OpenWork §2) · `GPWPlan1.md` item 5 |
| **Becke partition, the loop itself** | `-march=native` alone measured 1.13× with the loop still SCALAR; `BeckeImage` is array-of-structs with a data-dependent `P>0` exit that saves only 2.3% of work | vectorise first (chunk the exit, SoA); the `norm()` table (~2 MB gather vs a 15-cycle hardware sqrt) is DOUBTFUL and points the opposite way — re-measure after vectorising, do not assume.  Both are small beside V1.22 · History4 "Becke partition, what is LEFT" |
| **`Eval`/`EvalGradient` still run their own image loop** | the per-point callers (KB quadrature) were not migrated to the `LatticeSum1E` point-set seam because the seam re-derives its offset list per call; so the code's *"ONE remaining explicit image list"* is narrowed, not retired | cache the offsets on the seam side, then migrate; it retires the last place GPW enumerates images itself · History4 "THE Φ BLOCH POINT-SUM SEAM" |
| **Φ-sparsity for LARGE cells** | ⛔ refuted on MnO (48% dense) — the item's own caveat "the win grows with cell size" was the whole story | re-measure with `GPW_PHI_SPARSITY=1` on a battery-scale supercell BEFORE building anything |
| **The singles-vs-pairs census** | the factored density turns pairs into singles (8778 → 118 on MnO) but the comparable count is (i,R) singles vs (i,j,R) pairs at the same ε, weighted by box volume — pairs SCREEN far harder (Gaussian product theorem), which is why CP2K collocates pairs | build the singles-side census over the orbital reach before any singles-route work; the crossover is system-dependent and may favour BOTH routes · History4 "THE FACTORED FORM CHANGES THE OBJECT" |
| **Parity deviation #8** | `CP2K_COMPAT=1` is our best KNOWN parity, not proven; the deviation list grew 4 → 7 every time anyone looked | finding the next one is part of any parity claim; the like-for-like row is a `GPW_MNO_NMAX`-capped probe — step one is an uncapped or equal-cap comparison · `doc/Benchmark.md` §2, §5a |
| **Residual 1.33–1.45× CPU inflation at 12 threads** | load imbalance / memory bandwidth; ~4.1 s left in Phase 1, deliberately not taken | revisit on a 2×2×2 LiMn₂O₄ supercell where the buckets reshuffle — optimising the tail of a development cell is fitting to the wrong system |
~~~~
