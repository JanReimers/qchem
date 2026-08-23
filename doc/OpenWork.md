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
  **What happens now instead.**  `XC_GridEngine::RhoPol` takes the DM route only when the density it is
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

---

## NEW IDEA — fast evaluation of DM-ρ(r) by FACTORING D

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

### ✅ SESSION RESULTS 2026-08-21/22 — THE Vxc REPAIR IS BUILT, MEASURED, AND OPT-IN

**⛔ FIRST, TWO CLAIMS BELOW ARE NOW WRONG AND ARE CORRECTED HERE.**
1. *"at the fixed point ρ[D]=ρ_mix, so the converged answer is unchanged"* (route (b)) — **NO.**  ρ̃ is a
   band-limited FIT PROJECTION, so ρ̃_mix(r) ≠ ρ[D](r) pointwise even at convergence; the two have DIFFERENT
   fixed points.  MEASURED: NaF **139 μHa**, MnO individual terms **~100–150 mHa** (Ekin +110, Een −148) while
   the TOTAL moves 8 μHa — the variational signature.  Route (b) is therefore MORE ACCURATE, not neutral.
   ⚠ Consequence for **Step 5**: its term-by-term breakdown vs CP2K moves by ~100 mHa with this flag, so the
   flag must be PINNED before that comparison means anything.  (It does NOT explain the ~100 mHa offset —
   the total barely moves.)
2. The cost framing predating CALL COUNTS.  `report::Timed` now counts entries and prints `[xN, s/call]`.

**THE MEASUREMENTS** (MnO AFM-II imposed, 41 it, 97160 Becke points, n=122, r=17):

| bucket | calls | s/call | total | share |
|---|---|---|---|---|
| XC-mesh ρ sampling (matrix-free ρ̃) | 42 | **5.19** | **218.0 s** | **51%** |
| collocate density (pair scatter) | 330 | 0.105 | 34.7 s | 8% |
| XC-mesh quadrature H_xc | 168 | 0.175 | 29.4 s | 7% |
| integrate-back (pair gather) | 166 | 0.116 | 19.2 s | 5% |
| XC-mesh ρ sampling (DM GEMM, low-rank) | 41 | **0.042** | 1.7 s | 0.4% |
| (same, full-rank) | 41 | 0.296 | 12.1 s | — |

**⇒ the ρ̃ sampling was THE largest cost in the run, and the DM route is ~124× cheaper per sampling.**
Low-rank on the GEMM: **7.05×**, against the predicted n/r = 122/17 = 7.2.

**★ AND IT IS AN ACCURACY FIX, NOT ONLY A SPEED ONE.**  `GPW_RHO_NEGATIVE=1` census:

| route | points with ρ<0 (of 97160) | negative mass | min ρ |
|---|---|---|---|
| matrix-free ρ̃ | **8.5–18%, every iteration, persistent to convergence** | −0.003…−0.007 e | **−0.147** |
| DM / DM-source | **0**, over 200+ samplings and 2 cells | 0 | 0 |

A truncated Fourier series RINGS NEGATIVE at the cusps — and `SlaterExchange::GetVxc` guards `if (ro>0)` and
returns 0, so up to a fifth of the atom-centred quadrature contributed NOTHING to v_xc/E_xc, silently.
ρ=‖L†Φ‖² cannot do that.

**BUILT (all opt-in; default path byte-identical, suite 756/756):**
- `FactoredRho<Leaf>` — derived leaf, memoized per density serial, Tr(D) re-check against a stale memo.
  `IrrepCD_Factory(..., RhoRoute)`; `EigenTrim` throws (Cholesky's pivots ARE the spectrum).
  ⚠ The memo is the point as much as the GEMM: `LowRankFactor` used to run an O(n³) `pstrf` on EVERY call.
- `tDM_Sourced_CD` — the mixer retains the D its field was built from (SHARED ownership; aliasing
  `shared_ptr` for the polarized split) + **`EffectiveAlpha()`**.  `GPW_XC_DM_SOURCE=1`.
- Instruments: `Timed` call counts, `GPW_RHO_NEGATIVE`, `GPW_KERKER_SPECTRUM`, α_eff in the ρ_mix trace
  column (`Ker 0.45>0.33`), `GPW_XC_DM_BOOST`.

### ★★ α_eff — READ THE DAMPING OFF THE MIX INSTEAD OF SETTING IT

**The defect first:** route (b) as first wired fed XC ρ[D_out] UNDAMPED while Hartree kept Kerker's α·f(G) —
a HALF-DAMPED map.  NaF DIIS went 34 → 100+ iterations.  That is NOT evidence that ρ̃ is a better XC input
(user's challenge, and it was right); it is a wiring error.  Restoring damping fixes it.

\f[ \alpha_{\rm eff}=\frac{\lVert\tilde\rho_{mix}-\tilde\rho_{in}\rVert_2}{\lVert\tilde\rho_{out}-\tilde\rho_{in}\rVert_2}
   = \alpha\sqrt{\frac{\sum_G f(G)^2|\delta\tilde\rho|^2}{\sum_G|\delta\tilde\rho|^2}} \f]

— α times the residual-power-weighted RMS of the filter.  **MEASURED, not modelled**: two accumulators in the
loop `KerkerMix` already runs, so it needs no knowledge of WHICH preconditioner ran (exact for Kerker, = α
for a linear mix, may EXCEED 1 for Pulay — extrapolation overshoots on purpose).  Deposited on the field, so
the XC channel matches the damping without a SECOND independently-settable mixing policy.

| NaF stage 1 (DIIS) | baseline | undamped | flat α=0.25 | **α_eff** |
|---|---|---|---|---|
| iterations | 34 | **100+, FAILED** | 55 | **43** |

α_eff BEATS the hand-set flat value.  Converged energies are damping-INDEPENDENT (2e-8 across α) while
sitting 139 μHa from baseline — path vs fixed point, cleanly separated.

**THE BOOST (`GPW_XC_DM_BOOST`, default 1) IS JUSTIFIED AND IS NOT A DUPLICATE OF α.**

| iterations | baseline | undamped | α_eff ×1 | **α_eff ×2** |
|---|---|---|---|---|
| MnO | 41 | 46 | 53 | **39** (208.7 s vs baseline 425 s = **2.04×**) |
| NaF | 34 | 100+ | 43 | 47 |

No universal multiplier — MnO wants 2, NaF wants 1.  **Raising α instead does NOT work**: MnO at α=0.9
OSCILLATES with AND without the DM route (E=−45.2, 80 it), i.e. HARTREE has no headroom at 0.9 while XC ran
happily at 0.71.  α moves both channels; the boost moves only XC, which is the whole point — Kerker's
low-G damping is medicine for the 4π/G² divergence that XC's finite kernel never had.

### ★★ THE RESIDUAL SPECTRUM (`GPW_KERKER_SPECTRUM=1`) — AND AN MnO CONVERGENCE FINDING

Binned in x=G/G₀ (the filter's own argument, so cells compare):

| | 0.5–1 (f≈0.30) | 1–2 (f≈0.64) | 2–3 (f≈0.84) | 3–5 (f≈0.93) |
|---|---|---|---|---|
| **NaF**, iter 1 | **0.0%** | 21.6% | 43.5% | 31.9% |
| **NaF**, final | **0.0%** | 47.1% | 47.4% | 3.6% |
| **MnO**, iter 1 | 3.8% | 60.0% | 24.5% | 8.0% |
| **MnO**, final | **75.4%** | 17.5% | 5.2% | 1.5% |

**NaF's residual has NO power below x=1, ever** — Kerker barely filters it, so α_eff≈0.8α is a faithful
summary and there is nothing for a boost to recover.  **MnO's residual MIGRATES DOWN**: by the end 75% sits
where f=0.30, dragging α_eff 0.35 → 0.197.  That is exactly why ×2 helps MnO and not NaF — one mechanism,
both cells, visible in the data.
**★ INDEPENDENT OF THIS TRACK: late MnO convergence is rate-limited by the |G|≈0.65 CHARGE mode damped to an
effective 0.45×0.30 = 0.135.**  `GPW_SCF_UT.C` predicted that arithmetic; this measures which mode dominates,
which turns the `MNO_KERKER_G0` sweep from "a PREDICTION to be swept" into a measurement.

### ✅ GDM AT kT=0 IS UNAFFECTED (user's canary)

Inside the geodesic the Fock comes from `itsCD` and every line-search energy from `GetTotalEnergy(cdt)`, both
D-backed — so route (b)'s branch (which needs a density with NO DM face) cannot fire.  MEASURED on the NaF
stage-2 GDM: 9 it (α_eff) / 8 it (undamped) / 25 it (baseline), **zero** convention-shift events, no
line-search pathology.  It also REMOVES a discontinuity: mixed steps used to feed XC ringing ρ̃ while geodesic
steps fed exact ρ[D]; now both feed ρ[D].

### ★★★ MESH ⊥ ρ-ROUTE — AND THE MISSING CELL WAS NEVER MISSING (user, 2026-08-22)

ρ for XC has two INDEPENDENT axes: WHICH POINTS (uniform / Becke) and HOW ρ is produced (PAIR collocation /
SINGLES Φ-table, optionally factored).  `GpwOptions::vxcFit` and `xcMesh.cellKind` ALREADY express exactly
that (§6a fit/grid separation), and `VxcFit::Delta + UnitCellKind::Uniform` is supported, documented and
prints `[XC quadrature] DELTA fit on the uniform cell mesh`.  It had simply never been RUN — `Auto` gives the
historical pairing (PW on uniform, Delta on Becke).  PW+Becke is impossible (no G raster); polarized forces
Delta.  **Si Γ, Delta+Uniform, `SI_VXCFIT_DELTA=1`:**

| | lowrank=1 | lowrank=0 |
|---|---|---|
| U-4x | −7.115059004, **11 it CONVERGED** | −7.115059004, **11 it CONVERGED** |
| U-sel | −7.114980442, 60 it (degenerate manifold, benign) | −7.114980442, 60 it |

**BIT-IDENTICAL** — the clean exactness test, third material, uniform grid.  U-4x also beats the `Auto`
baseline (11 vs 14 it).  Si rank: **r=4–12 of n=24** (r=4 = its 4 occupied pairs), so low rank now holds on
ionic-magnetic AND covalent.  Low-rank speedup on ρ 1.26× at n=24 and H_xc UNCHANGED — the fourth
confirmation that **the rank helps the density side and never the operator side**.

⚠ **A FAILED ATTEMPT WORTH RECORDING (and its lesson).**  Before finding the above I added a `GPW_XC_SINGLES`
flag sourcing ρ from a Φ GEMM inside `PWFittedVxc` while leaving H_xc on the screened pair gather.  Si went
14 → 60 iterations and the energy moved 35 μHa.  **Vxc and Exc must reach D by ONE discrete path** (the
`GPW.RawXCConsistencyFD` gate is what enforces it); ρ from an unscreened GEMM paired with a screened adjoint
is two functionals.  `XC_GridEngine` owns ρ AND its quadrature over one Φ table, making the mismatch
unrepresentable; `PWFittedVxc` splits those owners, which is what permitted the bug.  Reverted.

### ★★★ DESIGN ITEM — SEPARATION OF CONCERNS IN THE XC TERMS (user, 2026-08-22)

**The ask:** `fit basis {Delta,PW,aux} x quad grid {Uniform,Becke} x op assembly {pair,singles} x
functional {LDA,GGA,hybrid}` should be freely combinable — *"the software should not hard code certain
combos and not others"* — **and every combination must be variational, so GDM just works.**  What follows is
what a day of measuring says the axes actually are; the short version is that one of the four is not an axis,
one coupling is mathematical and must be made structural, and one of our two XC terms is misnamed.

**1. ⛔ THE `PWFittedVxc` NAME IS A LIE ON THE PRODUCTION PATH.**  Since 0.5(f2) the RAW branch is the
default and it performs **no plane-wave fit at all** — the code's own comment says *"No ball fit anywhere"*.
It evaluates \f$v_{xc}\f$ pointwise on the fit raster and assembles H through `applyRawAdjoint`.  So it is a
**δ-basis quadrature**, exactly like `DeltaFittedVxc`; `VxcFit::PlaneWave` names a route it does not take.
(Discovered the hard way: a `GPW_FIELD_SPECTRUM` census of the fit coefficients printed nothing, because the
coefficients are never computed.)  The genuinely-fitted BALL branch survives only as a dormant fallback.

**2. AXIS C (pair vs singles) IS AN ALGORITHM, NOT A SEMANTICS.**  On the SAME points both compute
\f[ H_{ij}=\sum_g w_g\,v(r_g)\,\chi_i(r_g)\chi_j(r_g) \f]
— the pair route materialises \f$\chi_i\chi_j\f$ into multigrid boxes, the singles route contracts through
Φ.  Same formula, different evaluation order; they differ ONLY through the pair route's ε-screening and
multigrid truncation.  So it is a COST/TRUNCATION strategy, not a user-facing choice.
**⚠ But it is not a free choice either, and neither route always wins** (user's question): the pair cost
scales with the SCREENED pair count and wins where screening bites (large cells, diffuse bases); the singles
GEMM wins where n is small against npts or screening is weak — MnO's 4-atom cell measured a Φ-sparsity
ceiling of only ~2×, i.e. screening barely helps there.  So it is a **crossover**, and the right shape is the
one `UnitCellKind::Auto` already uses: an automatic selector with an explicit override, **latched for the
run** (switching mid-SCF changes the truncated operator, i.e. the functional — the same hazard the RAW/BALL
latch exists for).

**3. ★ THE ONE MATHEMATICAL COUPLING: the ρ-route and the H-route MUST BE ADJOINT.**  \f$H=\partial
E/\partial D\f$ holds only if the integrate-back is the EXACT adjoint of the collocation *on the same
truncated operator* — what `GPW.RawXCConsistencyFD` gates.  So these are **not two axes; they are one object
with two faces.**  `XC_GridEngine` already is that (`Rho()` and `Matrix()` over one Φ).  `PWFittedVxc` splits
them across two owners, which is exactly why a mismatch was expressible there — MEASURED 2026-08-22:
unscreened singles ρ + screened pair H sent Si from 14 to 60 iterations and moved E by 35 μHa.

**4. A FITTED REPRESENTATION IS VARIATIONAL IFF E IS DEFINED THROUGH THE FIT.**  \f$H_{ij}=\sum_k c_k\langle
c_k|\chi_i\chi_j\rangle\f$ is fine — with c the STATIONARY point of the least-squares fit in the metric that
also defines E (the robust/Dunlap form), the \f$dc/dD\f$ terms vanish and H is exactly \f$\partial
E/\partial D\f$.  This is standard in molecular DFT.  Today's BALL route is non-variational because E comes
from grid quadrature while H comes from the fit — TWO DEFINITIONS OF ONE ENERGY, not a defect of fitting.

**⇒ THE REAL AXES ARE THREE, AND THEY COMMUTE:**

| axis | what it is |
|---|---|
| **Quadrature** | points + weights + the adjoint-paired ρ/H routes; pair-vs-singles is INTERNAL strategy |
| **Vxc representation** | δ (identity) · PW · Gaussian aux — fitted ones use the ROBUST form |
| **Functional** | LDA · GGA · hybrid (already abstract: `ExFunctional`) |

**⇒ THE ABSORPTION.**  `PWFittedVxc`'s production behaviour IS `DeltaFittedVxc` over a raster-backed mesh
with the pair strategy.  Its points are already expressible (wrap `G_Quadrature::GridPoints()` with weight
Ω/N).  So:
- `XC_GridEngine` gains the strategy — same two faces, internally Φ-backed (singles) or 3C-tensor-backed
  (pair).  ⚠ One structural wrinkle: the pair strategy needs `Overlap3C(*fitBasis)`'s
  `applyRaw`/`applyRawAdjoint`, so a pair-strategy engine takes the FIT BASIS where a singles engine takes
  only the mesh.  Different ctors, same faces.
- The quadrature XC term (today's `DeltaFittedVxc`, itself misnamed — its "fit basis" is the identity)
  becomes THE term.  `PWFittedVxc` has nothing left on the production path.
- What remains of it is the BALL branch — and THAT is what deserves the "Fitted" name, rebuilt in the robust
  form.  `itsRhoIsRaw`'s latch, the BALL fallback and the real-TRIM throw re-home there.
- `VxcFit` then selects REPRESENTATION (identity vs fitted), not secretly the grid AND the assembly.

**⇒ AND "GDM JUST WORKS" BECOMES STRUCTURAL**, resting on exactly two things: ρ and H from one adjoint-paired
object (a mismatch becomes unrepresentable), and any fit defined in the robust form.  The FD gate then
becomes a test PARAMETERISED OVER ALL COMBINATIONS, so a new one cannot ship non-variational.

**THE CELLS THAT ARE UNREACHABLE TODAY — and note two of them need no Becke mesh at all.**  "Uniform" is TWO
point sets in this code: the fit raster (corners \f$i/N\f$, N from the fit G-ball) and
`UnitCell::CreateIntegrationMesh` (midpoints \f$(i+\frac12)/n\f$).  Each is welded to ONE assembly:

| points | pair | singles |
|---|---|---|
| fit raster (corners) | ✅ `PWFittedVxc` | ❌ |
| UnitCell mesh (midpoints) | ❌ | ✅ Delta+Uniform |
| Becke | ❌ (no box structure to gather into — the genuinely hard one) | ✅ Delta+Becke |

**MEASURED INPUTS ALREADY IN HAND.**  `VxcFit::Delta + UnitCellKind::Uniform` is supported, documented (§6a
fit/grid separation) and had simply never been RUN — Si there is BIT-IDENTICAL under factoring and converges
11 vs 14 iterations.  `GPW_MESH_ORTHO=M` measures \f$T(\Delta G)=\sum_g w_g e^{i\Delta G\cdot r_g}\f$, the
whole non-orthogonality of a PW basis under a mesh: uniform ~1e-15 (exact), Becke 1e-3 rising to 2.6e-2 by
\f$|\Delta G|>6\f$, and L=29→35 improves the high-\f$|\Delta G|\f$ end 2.7× and the low end not at all —
the angular rule running out.  So a PW fit on Becke is WELL CONDITIONED where it matters.  ⚠ And the fit
basis is built at `relCutoff * {G}_rho` with `GridCutoffFactor()` = 1.0 for every functional, so v_xc uses
the FULL density ball (24k G on MnO) — never swept below 1.  At 1/3 the radius that is ~900 vectors and the
metric solve becomes trivial; the honest way to pick it is E_xc convergence, not a rule of thumb.

**⇒ THE MECHANISM: `CreateVxcFitBasisSet` MAKES ALL THE DECISIONS, AND KEEPS ITS NAME (user, 2026-08-22).**
The basis already has TWO factories with the SAME inputs returning different halves of one answer --
`CreateVxcFitBasisSet(st,mp)` and `CreateXCQuadrature(st,mp)` (the latter already a BUNDLE:
`{mesh, fold, sigmas, flipFixed}`).  `Ham_PW_DFT::BuildTerms` then re-decides on `becke`/`polarized`/`delta`
and picks a TERM CLASS -- which is where the combinations get hard-coded.  Fix: the quadrature is ABSORBED
into the fit basis and `BuildTerms` stops switching, because **under the weight-vector view a fit basis is
not defined without points** (\f$W_c[g]=w_g c^*(r_g)\f$), so the mesh is CONSTITUTIVE of the basis, not a
sibling of it.  `PlaneWaveFit_IBS` already works this way (it IS a `G_Quadrature`).
- ⚠ NOT an "XCEngine" (user): *"Engine just means something that does some work, so it adds no info"* --
  and `XC_GridEngine` is already an instance of that non-name.  The factory keeps its name and its job.
- **The δ case stops being a null-object smell.**  The 2026-08-01 objection was to a ZERO-FUNCTION
  pseudo-basis; a **delta basis** is a different object -- n_pts genuine δ functions, exact on that point
  set, metric DIAGONAL (\f$\langle\delta_g|\delta_{g'}\rangle=w_g\delta_{gg'}\f$), fit = identity.  Under
  the old framing there was nothing to represent; under "a fit basis is a family of weight vectors" it is the
  most natural basis there is, and the doc's *quadrature is the δ-basis special case of fitting* becomes a
  TYPE instead of an analogy.
- **The assembly strategy needs NO extra input**: both routes are already functions of (orbital basis, fit
  basis) -- pair is `orb.Overlap3C(*fitBasis).applyRawAdjoint`, singles is Φ at the fit basis's points.
  Capability decides what is POSSIBLE, policy (beside `ResolveXCMesh`) decides which.

**⇒ THE ACTUAL CHANGE SET IS SMALL AND LOCAL -- `FIT_SF_ABS<T>` NEEDS NO WIDENING.**  ⛔ (I first claimed a
five-basis blast radius; wrong.)  That face is thin (`isOrtho` + `SymmetrizeRaster`) and the quadrature is
already reached by CAPABILITY CROSS-CAST, not through it (`OrthoScalarFitter::FitGrid()` casts to
`G_Quadrature`).  So the molecular bases (Atom, SymmetryAdapted, the generic default) are untouched.  What
actually blocks the missing cells is that **`G_Quadrature` BUNDLES two different things**:

| member | universal? |
|---|---|
| `GridPoints()`, `Integral(f)` | ✅ any point set / any weights |
| `RhoOnGrid`, `ForwardFFT`, `GridCoeff`, `FieldCoeffs` | ❌ **RASTER ONLY** |

A δ-basis-over-Becke can provide the first row and cannot provide the second, and today they are one face,
so it cannot answer at all.  **And the narrow face already exists: points+weights IS `qcMesh::Mesh`** (user's
standing guidance to use it as much as possible; NB its weights need not be quadrature weights -- they can
carry PW phases, which is what makes the weight-vector view a type and not a metaphor).  So:
1. split `G_Quadrature` → a `qcMesh::Mesh` accessor (universal) + the raster transform face (`PW_Grid_Evaluator`);
2. `OrthoScalarFitter::FitGrid()` splits the same way (points/Integral off the mesh, transforms off the raster);
3. a δ fit-basis type carrying a `qcMesh::Mesh`;
4. `CreateVxcFitBasisSet` picks basis + mesh; `CreateXCQuadrature`'s bundle folds into the returned basis.
⚠ One sizing trap: a δ basis has n_pts "functions", so anything reasoning about a fit basis by COUNTING
functions (grid sizing, memory reports, the `relCutoff` arithmetic) must not choke on 97160.

#### ✅ DONE 2026-08-22 — all four steps + THE ABSORPTION.  756/756; the two-route Si gate BIT-FOR-BIT

Steps 1–4 as written, then the absorption the section calls for (user: "also do the absorption").  What
landed, in the order it was built and tested:

| | |
|---|---|
| **1** | `G_Quadrature` → **`BasisSet::Quadrature`** (`Mesh()`, universal — new leaf module `qchem.BasisSet.Quadrature`) + **`G_RasterTransform`** (`RhoOnGrid`/`ForwardFFT`/`GridCoeff`/`FieldCoeffs`, honestly raster-only).  `PeriodicGridEvaluator` now STORES its raster as a `qcMesh::Mesh` (weights Ω/N) and `GridPoints()` reads it, so there is one copy, not two.  `qcMesh::Integrate(mesh,f)` added beside `SiteIntegrals`. |
| **2** | `GriddedScalarFitter::Grid()` → `Mesh()` + `Raster()`; `OrthoScalarFitter` cross-casts each half separately. |
| **3** | **`BasisSet::FIT_SF_Delta<T>`** (face) + **`DeltaFit_IBS`** (concrete, Bloch Γ) — the δ basis carries the `XCQuadrature` bundle (mesh + fold + σ tags), `Mesh()` from it, `isOrtho`, and `op(r)`/`Gradient` **throw** (user ruling: the value is a distribution and nothing wants it; `MakeOverlap` never appears because it lives on `FIT_SF_NonOrtho`, not on the face δ implements). |
| **4** | `VxcFit` moved DOWN into `qcMesh::MeshParams` (alias kept in `qchem::Hamiltonian`) so `CreateVxcFitBasisSet` can read it: ONE factory call returns δ-over-mesh or PW-over-raster.  `BuildTerms` resolves only `Auto` — the one decision that needs to know the run is polarized — and stamps it. |
| **absorption** | **`XC_Quadrature`** = the abstract two-faced object (ρ, and its exact adjoint H).  Two strategies: **`XC_SinglesQuadrature`** (ex `XC_GridEngine`, Φ table) and **`XC_PairQuadrature`** (ex `PWFittedVxc`'s guts — raw collocation + `applyRawAdjoint`, the BALL fallback, the R2.16 route latch).  `MakeXCQuadrature(fitBasis)` picks by CAPABILITY.  `PWFittedVxc` is GONE; `DeltaFittedVxc*` → **`Vxc_Quadrature` / `Vxc_QuadraturePol` / `Vcorr_QuadraturePol`**, one term set for every combination. |

**MEASURED, not asserted.**  `GPW_SCF.DeltaFitUniformGridMatchesPWFit_SiGamma` drives BOTH strategies
end-to-end; stashed/rebuilt/reran to compare against the pre-change binary:

| route | before | after | iters |
|---|---|---|---|
| Si PW-fit (now the PAIR quadrature) | −7.115067665 | **−7.115067665** | 11 → 11 |
| Si Delta uniform (now SINGLES) | −7.115059008 | **−7.115059008** | 11 → 11 |

Full sweep 756/756 (169 s), including the FD variational gates and the real/complex bitwise gate on the
pair route.  ⚠ One honest caveat: E_xc on the absorbed pair route is now accumulated as
`Σ (w·ε)·ρ` where it was `Σ w·(ε·ρ)` (the shared term body), and step 1 replaced the raster's
`(Σf)·Ω/N` with `Σ w_g f_g` — algebraically identical, last-digit different in principle; nothing
measurable moved at the 1e-9 the gate prints.

**FREE SIDE-EFFECT of the absorption:** the raw collocation now runs ONCE per iteration, not once per
TERM — the exchange and correlation terms share one `XC_Quadrature` the way the Becke pair always did.

**WHAT DID NOT CHANGE, deliberately.**  `CreateXCQuadrature` survives as the basis-side mesh BUILDER that
the whole-set `CreateVxcFitBasisSet` calls for the δ case (and that several tests drive directly); it is
simply no longer on the Hamiltonian's path — no caller re-decides anything.  The robust/Dunlap rebuild of
the genuinely-fitted BALL branch (item 4 above) is still open, as is PW-fit-on-Becke (I3) and the
pair-vs-singles auto-selector: capability picks the strategy today, so an override knob is the remaining piece.

#### REVIEW FIXES, same day (user, on reading the above)

**⛔ `VxcFit` DOES NOT BELONG IN qcMesh.**  It rode there as a `MeshParams` field so the factory could read
it — but `MeshParams` describes POINTS and `VxcFit` describes FUNCTIONS, which are *precisely the two axes
this whole section exists to separate*, so folding one into the other's parameter block re-welded them at
the type level while the code was busy unwelding them.  Moved to `qchem.BasisSet.Fit_IBS`, beside the
fit-basis types it selects between, and passed as its own ARGUMENT: `CreateVxcFitBasisSet(cl, mp, fit)`.
qcMesh is pure geometry again.

That argument sits on the WHOLE-SET `tBasisSet` factory only — the level the Hamiltonian actually asks at.
A Bloch BLOCK builds only its own lineage's fitted basis (the per-block δ branches are reverted), and
`DeltaFit_IBS` moved from qcLattice_BS up to qcBasisSet, since a δ basis over a mesh is not periodic-specific.
Molecular callers pass `Auto` — "give me your own representation" — and never interpret the other two;
`tBasisSet<double>` throws on an explicit `Delta`, because there is no real δ basis to build.

**⛔ `GriddedScalarFitter` HAD NO CONSUMERS — DELETED.**  Its whole purpose was "the XC term borrows the
grid THROUGH the fitter, one owner".  After the absorption the quadrature object holds the fit basis and
reads `BasisSet::Quadrature` directly, so `Raster()` had **zero** callers and `Mesh()` had one, redundant
with the caller's own.  A face nobody consumes is not an abstraction: `Factory(cFIT_SF_ABS)` now returns
plain `FunctionFitter_Scalar<dcmplx>`, and where the ortho fitter samples is its own business — exactly as
the AO fitter's Becke mesh is.

**⛔ `FIT_SF_Delta` DID NOT EARN ITS EXISTENCE — it does now.**  As first built its entire content was
`XCQuad()`: one getter, one consumer (`MakeXCQuadrature`), zero behaviour — a type whose only job was to
carry a struct across a cast, i.e. the 2026-08-01 null-object objection in a new costume (a basis with no
FUNCTIONS then, a basis with no BEHAVIOUR now).  **A fit basis is for DOING the fit and delivering E and
the H matrix, not for holding its functions — in the δ case its grid — and giving them away** (user).  So
the getter is gone and the face answers operations, every one of them over a field SAMPLED AT ITS OWN
POINTS, in its own order, so no caller learns where those points are or what kind of mesh they came from:

| | |
|---|---|
| `Integrate(f)` | \f$\int f=\sum_g w_g f_g\f$ — the E_xc quadrature |
| `SiteIntegrals(f)` | \f$\int w_A f\f$ per site block; EMPTY when the mesh has no atomic partition |
| `Sample(field)` | a field that evaluates itself anywhere, sampled at my points (matrix-free densities) |
| `Symmetrize(f)` | the orbit-mean projector; no-op on a free run, so the caller never asks whether symmetry was imposed |
| `SymmetrizeSpin(ρ,m)` | the Shubnikov (ρ even, m odd, flip-fixed zeros) projection — was a three-way branch in the consumer |

`XC_SinglesQuadrature` is now built ON that basis and keeps only POLICY (which source ρ comes from, the
per-serial caches, the spin channels, the DM-source damping).  The bundle's invariants and its `[fold]`
announcement moved into the basis ctor, where the object is made — providers self-report.

**⛔ AND "QUADRATURE" IS AN IMPLEMENTATION DETAIL, SO IT LEFT `BuildTerms`** (user).  The builder now says
*what the model is* and asks for XC terms: `MakeVxcTerms(exch, corr, XFitBasis, polarized)` returns the
pair, and the strategy choice, the sharing of one quadrature across the pair (ρ and the Φ tables built once
per ITERATION, not once per term), and the adjoint pairing are all inside.  The TERMS lost the mesh too:
they integrate through `XC_Quadrature::Integrate` and report through `NumPoints()`, so no term touches
weights or points.  `MakeXCQuadrature` picks by CAPABILITY — has the basis got `G_RasterTransform`? — never
by asking what the basis IS.

**WHAT IS LEFT, and where it lives.**  Exactly one coordinate escape survives: `Quadrature::Mesh()`, needed
because the ρ GEMM is initiated by the DENSITY (`cDM_CD::DM_RhoAtPoints`, which owns D and lives in a
library above qcBasisSet), so points and the Φ cache must still reach it.  It is confined to the two
`XC_*Quadrature` strategies — no term, no Hamiltonian, no fit-basis client sees it — and closing it means
flipping `DM_RhoAtPoints` to take the δ basis and ask it per block, which brings Φ and BOTH contractions
into the basis and touches qcChargeDensity's IrrepCD/composite tree.  That is the only piece of this design
item still outstanding; everything else about δ-vs-fitted and pair-vs-singles is settled.

⇒ **AND IT IS NOT ITS OWN ITEM: it closes as an instance of `doc/CleanupCandidates.md` R1.0** (the "FAKE
RADIAL `op(r)`" cure — *consumers stop touching `op(r)`; code that wants an INTEGRAL asks the basis for
it*).  Φ_gi = ⟨δ_g|χ_i⟩/w_g is a cross-basis overlap, so the same discipline that fixes the atomic fake
radial turns `DM_RhoAtPoints(points, tables)` into `DM_RhoAtPoints(basis)`.  R1.0 now carries the scope
(step (1) only — step (2), honest Y_lm `op(r)`, is explicitly out, user 2026-08-22), the enabling move
(`VectorFunction<T>` off `IrrepBasisSet<T>`, onto `Orbital_1E_IBS<T>` — four pointwise sites, three of
them orbital), and the naming trap (a δ basis's `Sample` is the FIT coefficient vector, NOT
`MakeOverlap`, since ⟨δ_g|f⟩ = w_g f(r_g)).

Re-verified after all three: 756/756, and the two-route Si gate still prints −7.115067665 / −7.115059008 at 11/11 iterations.

---

## ★★★ ONE FIT-BASIS INTERFACE (user, 2026-08-22) — ✅ DONE 2026-08-23 (the spec follows, annotated)

#### ✅ WHAT LANDED.  758/758 (count UP by one gate, not down); Si two-route gate BIT-UNMOVED at 11/11.

`FIT_SF_ABS<T>` lost `NumPoints` / `Sample(f)` / `Integrate(values)` and gained the two integrals every
fit basis can answer **per FUNCTION** — `Overlap(f)` = ⟨f_a|f⟩ (moved UP from `FIT_SF_NonOrtho`) and
`Integrals()` = ⟨f_a|1⟩ (on the Gaussian side, literally `Charge()`).  The count is `GetNumFunctions()`,
which the removed `NumPoints` was shadowing.  Each fitter applies its own metric to that one projection:
`BasisSet::OrthogonalFit(b,f)` = ⟨f_a|f⟩/⟨f_a|f_a⟩ is the δ fitter's whole `DoFit` **and** the default
`DM_RhoAtPoints`, which is where `OverlapDiagonal()` finally got a production consumer.  Full record,
including both risks measured rather than assumed: `doc/CleanupCandidates.md` R1.0 increment 3.

Three things the spec below did not anticipate, all recorded there in full:
- **⚠ the plane-wave fitter must NOT fit through `Overlap`, and not for bitwise reasons** — the XC
  assembly looks ṽ up at orbital index differences that run to twice the orbital ball, outside the fit
  ball, where the honest projection is ZERO and the raster's `GridCoeff` aliases.  Routing it through the
  ball would delete most of v_xc.  So `G_RasterTransform` grew `RasterSize()`/`Sample(f)`/`Integral(v)` —
  point vocabulary belongs on the face that is honestly about voxels.  (Not a name clash: the PW fit basis
  carries no self-overlap ⟨i|j⟩ at all — that lives on the orbital `EPW_Orbital1E_IBS` tier, which an
  auxiliary fit basis does not ride.)
- the ⚠ BITWISE pin was real (δ's fit is now `fl(w·f)/w`, not `f`) but ~1e-13 Ha against a 1e-10 gate.
- **`Values` → the 3-centre overlap (row 4) ⛔ first written off, then LANDED the same day.**  "It needs an
  `Overlap3C(δ)` overload on every orbital lineage" follows only if the ORBITAL basis receives; with the
  FIT basis receiving (user: δ does \f$H_{xc}\f$ itself, needing only `op(r)`) no orbital lineage is
  touched.  `Values` and `Quadrature(orb,v)` are gone, replaced by `Overlap3C(orb)`; `Integrals()` became
  `Charge()` per the house convention (⟨i|j⟩/⟨i|f|j⟩/⟨i|f⟩ = Overlap, ⟨i|1⟩ = Charge, 2e = Repulsion).
  The fitter holds the coefficients and contracts — the same two calls the Gaussian fitter makes.
  Then TWO more left the δ face (user): `Overlap3C` went up to `FIT_SF_ABS<T>` with the default
  `return orb.Overlap3C(*this)` — Gaussian and PW inherit their existing tensor untouched, δ overrides —
  and `SymmetrizeSpin` went up with the default `{Symmetrize(rho); Symmetrize(m);}`, which is bit-identical
  to the branch δ already ran without σ tags, leaving δ only the Shubnikov case.  **`FIT_SF_Delta` is down
  to three declarations**: `isOrtho()`, the mixed real-TRIM `Overlap3C` overload, and `SiteIntegrals` —
  and then `SiteIntegrals` went too, **by INJECTION**: `CreateVxcFitBasisSet` — which BUILDS the mesh —
  now also hands the same `shared_ptr<const qcMesh::Mesh>` to the XC strategy, which takes the moments
  where ρ_σ is already cached.  No getter; a creator gave its creation to two collaborators.
  **`FIT_SF_Delta` is down to ONE declaration**: the mixed-run real-TRIM `Overlap3C` overload.
  ⛔ **And the probe that verified the injection found a DEFECT — tracked as R1.0d:** an IMPOSED-symmetry
  Becke mesh silently loses its site blocks (the orbit-consistency filter re-emits into a fresh
  `MeshBuilder` without `BeginSite`), so per-site moments vanish on exactly the runs that want them — free
  1-atom run reports `NSites()==1` and prints moments; both imposed 2-atom runs report `0` and are silent.
  The MnO AFM campaign is imposed by construction.
  ⚠ `ForwardFactoredT<dcmplx>` is exercised by NO enabled test (measured) — pre-existing, not a
  regression, but live code on the default route.  Tracked as **R1.0c** in doc/CleanupCandidates.md with a
  cheap unit gate (`applyRawFactored(L)` == `applyRaw(LL†)`, both scalars).

#### ★ NEXT INCREMENT, SPECCED NOT BUILT (user, 2026-08-23) — SEPARATE THE METRIC AXIS INTO FACES

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

### The spec as written (user, 2026-08-22)

**THE GOAL, in one sentence.**  A δ fit basis is a family of FUNCTIONS — n_pts genuine δ functions with a
diagonal metric — so it must present EXACTLY the interface a Gaussian auxiliary fit basis presents, and
every "point" word must come off the fit-basis faces.

**WHY (the finding that opened it).**  The 2026-08-22 work put `NumPoints()` / `Sample(field)` /
`Integrate(values)` on `FIT_SF_ABS<T>`.  Those are the IMPLEMENTATION leaking through: they describe a
quadrature, not a family of functions.  They looked right only because for a δ basis n_functions ==
n_points, so the wrong accessor gives the right number.  The user's ruling: *"a DeltaFitBasis is also a
family of functions ... as such it should have exactly the same interface as the Gaussian fit basis set."*

### THE INTERFACE DELTA

| remove (point vocabulary) | replace with (function vocabulary) | δ's realization |
|---|---|---|
| `FIT_SF_ABS<T>::NumPoints()` | `IrrepBasisSet<T>::GetNumFunctions()` (already there) | n_pts |
| `FIT_SF_ABS<T>::Sample(f)` → per-POINT values | **`Overlap(const ScalarFunction<double>&)` → `vec_t<T>`, per FUNCTION** — moved UP from `FIT_SF_NonOrtho`, where it already exists, onto `FIT_SF_ABS<T>` | \f$\langle\delta_g\|f\rangle = w_g f(r_g)\f$ |
| `FIT_SF_ABS<T>::Integrate(values)` | the coefficients dotted with the per-function integrals \f$\langle f_a\|1\rangle\f$ (`Charge()` already has that shape on the CD side; add its SF sibling) | \f$\int\delta_g = w_g\f$, so \f$c\cdot w = \sum_g w_g f_g\f$ |
| `FIT_SF_Delta::Values(orb)` (Φ tables) | the 3-CENTRE overlap \f$\langle\delta_g\|\chi_i\chi_j\rangle\f$ every fit basis provides; its adjoint contraction is today's `Quadrature(orb,v)` | \f$w_g\chi_i(r_g)\chi_j(r_g)\f$ |

**`Overlap(field)` unifies across all three representations**, which is the load-bearing claim of the whole
increment — verify it early:
- **Gaussian**: projects on its own Becke mesh (the existing `Fit_IBS::Overlap(Sf)` body, unchanged).
- **δ**: \f$w_g f(r_g)\f$.
- **plane wave**: the FORWARD FFT — literally what `OrthoNormalScalarFitter::DoFit` does by hand today
  (sample, `ForwardFFT`, `FieldCoeffs`).  Moving it onto the basis is what makes the fitters uniform.

Each FITTER then applies its own metric to that one projection: \f$S^{-1}\f$ solve (Gaussian), divide by
`OverlapDiagonal()` (δ), or nothing (orthonormal PW).  **This is where `OverlapDiagonal()` finally gets a
consumer** — it has none today beyond its gate.

### THE XC TERM RE-PLUMB (what makes the deletions possible)

The molecular term is already the template (`Imp/FittedVxc.C`): compose a FIELD \f$v_{xc}\circ\rho\f$, hand
it to `DoFit`, contract with `Overlap(orbitalBasis)`.  Port the periodic terms to the same two calls:
- **H**: `ContractAdjoint(Overlap3C(fitBasis), c)` — for δ that IS \f$\Phi^\dagger\mathrm{diag}(wc)\Phi\f$.
- **E_xc**: \f$\sum_a e_a\langle\rho|f_a\rangle\f$ — the ε-fit coefficients against the density's projection
  ONTO the fit basis.  For δ, \f$\langle\rho|\delta_g\rangle=w_g\rho(r_g)\f$, so it stays \f$O(n_{pts})\f$
  and reuses the ρ samples already in hand.  (Do NOT route it through an E-MATRIX + `DM_Contract`: that is
  the molecular shape and it costs an extra \f$O(n_{pts}n^2)\f$ GEMM per iteration here.)

**THE ONE IRREDUCIBLE POINT-NESS, and it does not belong in a signature.**  \f$v_{xc}\f$ is pointwise
nonlinear, so something must evaluate ρ, apply the functional, and project.  That stays INSIDE: the term
composes the field, the BASIS samples it inside `Overlap(field)`, the field's BULK evaluation asks the
density, and the density gets its Φ tables by asking the basis (`DM_RhoAtPoints(q)`, already landed).  The
user's D-GEMM ruling survives intact and no coordinate appears in any interface.

### PINS AND TRAPS

- ⚠ **BITWISE.**  δ's fit is \f$c_g=\langle\delta_g|f\rangle/w_g = w_g f_g/w_g\f$.  Today `Sample` returns
  \f$f_g\f$ DIRECTLY because the \f$w_g\f$ cancels exactly; going through multiply-then-divide can move the
  last bit of every coefficient, and the pinned gates print 10 digits.  Decide deliberately: either keep the
  cancellation explicit inside the δ fitter (recommended — same interface, exact arithmetic) or accept the
  ulp move and re-pin.  **Measure first**: `GPW_SCF.DeltaFitUniformGridMatchesPWFit_SiGamma` must still
  print −7.115067665 (pair) and −7.115059008 (singles) at 11/11 iterations.
- **MIXED SCALARS.**  A real TRIM block wants `hmat_t<double>` while its complex siblings want
  `hmat_t<dcmplx>` (doc/RealComplexPlan.md 3c-3).  That is what `Fitting::FitContraction<U>` exists for —
  do not collapse it back into a single-scalar face.
- **NO GETTERS.**  Nothing may hand out points, weights, or a mesh.  The current tree has zero such
  escapes; keep it that way (the whole 2026-08-22 arc was closing them one at a time).
- **NAMING.**  `Sample` ≠ `Overlap`: the first is per-POINT values, the second per-FUNCTION integrals.  They
  coincide in length only for δ.  Do not blanket-rename onto the `MakeOverlap` family without re-deriving
  each one (R1.0 records this trap).
- **`SiteIntegrals` does not belong on a fit basis at all** — it is an atomic-moment OBSERVABLE.  Move it to
  whatever owns the partition; it is the last non-fit-basis question on `FIT_SF_Delta`.

### EXPECTED END STATE

`FIT_SF_Delta` mostly dissolves: δ becomes a fit basis whose metric is diagonal, distinguished by
`isOrtho()` + `OverlapDiagonal()` rather than by having its own face.  What legitimately remains δ-specific
is `SymmetrizeSpin` (the magnetic (ρ,m) pair projection — δ is the only representation a polarized run can
use, a plane-wave fit having no per-channel collocation).

### VERIFICATION

`ninja allTests && scripts/memsafe ctest -j8` — **757 tests, all passing**, and the count must not DROP
(a vanished test is a silent regression; see CLAUDE.md on `_NOT_BUILT`).  Beyond the suite: the Si
two-route gate above (both numbers, both iteration counts), `M_DFT.*` for the molecular lineage, and
`GPW.OverlapDiagonalPerRepresentation` which pins the three metric diagonals.

### SUGGESTED ORDER

1. `Overlap(field)` onto `FIT_SF_ABS<T>`; implement for PW (forward FFT) and δ; Gaussian inherits its
   existing body.  Suite green here — nothing consumes it yet.
2. Fitters take their projection from it and apply their own metric.  `OverlapDiagonal` gets its consumer.
3. Periodic XC terms onto the two-call shape (H via `Overlap3C`, E via coefficients·integrals).
4. Delete `Sample` / `Integrate` / `NumPoints` / `Values`, and whatever of `FIT_SF_Delta` is left empty.
5. Then the ruled MOLECULAR conformance (ρ through the D-GEMM; sample ONCE and derive both \f$v_{xc}\f$ and
   \f$\epsilon_{xc}\f$ from it) — worth MEASURING first, since molecular meshes are small and the win may
   be uniformity rather than time.

Background and the full 2026-08-22 record: `doc/CleanupCandidates.md` R1.0.


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
the XC mesh already runs":** Φ-cache + per-iteration GEMM is the existing `XC_GridEngine` pattern, so the
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
   =(\Phi^\dagger W\Phi)_{ij}\f$ is already a Φ-shaped GEMM in `XC_GridEngine` — so the XC side is ALREADY
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
