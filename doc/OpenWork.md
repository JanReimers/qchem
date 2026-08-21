# Open Work — the live tracker (v2, cut 2026-08-19)

**READ THIS AT SESSION START.**  One place for everything in flight, in the order we intend to do it.
The previous tracker — closed threads A–E, the 2026-06-30 and 2026-08-15 orderings, and the
runtime-campaign rounds 1–4 as written at the time — is retired to `doc/OpenWork_History1.md`.  Only
still-open work lives here.

---

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
**★★ MEASURED 2026-08-19 — THE MOMENT GAP WAS THE INSTRUMENT.**  With `XC_GridEngine` now reporting the
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

**0c — MAKE THE ORDER LEGIBLE.  THE CONSOLE CANNOT TELL YOU *WHEN* (user, 2026-08-20).  OPEN.**

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
- **Becke partition, what is LEFT in that loop** (none of it ε-related; all found during the census):
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
- **★★ THE ρ GEMM — LOW-RANK D.  The successor to the refuted Φ-sparsity item, and the only remaining plan
  for the largest per-iteration bucket.**  `IrrepCD_Core::DM_RhoAtPoints` forms
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
  **This is now the best-measured un-taken lever in the code: ~7–8× on the LARGEST per-iteration bucket.**
  ⚠ Instrument note: the PSD test must be on a RELATIVE floor.  The first cut tested \f$\lambda_{\min}<0\f$
  and screamed "NOT PSD" on every run at \f$\lambda_{\min}\approx-1.8\times10^{-15}\f$ against
  \f$\lambda_{\max}=13.2\f$ — one ulp.  An eigensolver ALWAYS returns O(eps·λmax) negatives for a
  numerically-PSD matrix; a real negative eigenvalue from an extrapolated D would sit orders of magnitude
  above that floor.
  *Not a clever trick:* evaluating ρ in ORBITAL space rather than DM space is what most codes do; the only
  observation here is that this one took the DM route, which is a reasonable choice that happens to cost n/r.
- **★ NEW DEFECT, found by the instrument above: THE IMPOSED XC MESH LOSES ITS SITE BLOCKS — and the
  Step 0a site-moment instrument is SILENTLY DEAD on every imposed run.**  `MakePeriodicBeckeMesh` records
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

### Step 4 — RAM  ·  read it off Step 1's table, then decide

**✅ LARGELY ANSWERED BY STEP 2, as predicted — do not open it as a track.**  Arming the T3 fold took the
MnO AFM row's peak RSS **4947 → 1349 MB** (23× CP2K → 6.2×) on top of round 3's 5.78 → 3.70 GB, because the
RAM was the streams and folding cuts their *demand* by the orbit factor.  What remains is the Φ tables
(Step 3's screening item, which is a RAM lever as much as a time one) and the fact that a FREE run folds
nothing and therefore still pays the full 4947 MB.  Re-read the number off the table after Step 3; reopen
this only if it still binds then.

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
- Δρ/N convergence gate (`doc/SCFStrategyPlan.md`); GDM fallback-diagonalize breadcrumb (run 59's silent
  +302 mHa hop); the per-channel ortho duplication above; the fingerprint's overconfident
  "raise NMaxIter" advice.

---

## Parked threads (real work, deliberately not in the plan above)

- **A. Spherical SALC — S3b, the libcint-spherical extractor**  ·  `doc/SphericalSALCPlan.md`.  The ONLY
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
