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
instrument is a table we redo).

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
**Outcome:** either the moment gap collapses toward the ~1.5× the Fourier measure already suggests, or it
survives as a sharp question — for a fraction of one overnight run.  Historical `m_stag` numbers in the
plan docs stay readable only as collapse/survival evidence; they are not moments and must not be quoted as
such.

**0b — make the FOLDS visible.**  The user's standing complaint, still true: `GPW_STREAM_FOLD` is opt-in
(`src/BasisSet/Lattice_3D/Imp/BasisSet.C:202`) and a grep for a fold factor on cout/cerr returns
nothing.  Print a fold-factor line in every grid/stream report — orbit order actually used, representative
vs total pair counts, the {G}-star reduction where it is armed.  `grids.stream` already carries
`foldOps`/`repPairs`; this is about the console and about the per-iteration sites that report nothing at
all.  Needed by Step 1 (a runtime table must say whether folding was on) and by Step 2 (which is measured
in fold factors).

### Step 1 — THE HEAD-TO-HEAD TABLE, built as a standing benchmark  ·  plan: new, see below

Si, NaF, MnO — **same basis in both codes**, ≥1 multi-k row — reporting **wall time, peak RAM and energy**
for qchem and CP2K side by side.

**This is the instrument, not a report.**  Everything optimised in rounds 3–4 was measured against a
single hand-run disabled test with `GPW_MNO_NMAX=3`.  Steps 2–4 have no finish line without this.

- **Like-for-like IS available and does NOT wait on Step 6** (user, 2026-08-19): CP2K can be forced down
  to the spherical spans qchem holds at full rank.  The decks and basis entries already exist —
  `IntegrationTests/CP2K/mno_{afm2,fm}_gpw_v{a,b}.inp` and the `VALENCE-LOWQ-V{A,B}` entries in
  `IntegrationTests/CP2K/VALENCE-LOWQ-BASIS`, variants **VA (N=118)** and **VB (N=128)** — and qchem has
  run VA at full rank ("kept 118 of 118", runs 61/62).  Si and NaF need the same treatment: one shared
  span per material, held full-rank by both codes, or the row is not a comparison.
- **The MnO accuracy row is already banked** (`doc/SphericalLatticePlan.md`, the VA verdict table); what is
  missing everywhere is the **runtime and RAM columns**, the Si/NaF rows, and the multi-k row.
- Multi-k cost is currently unmeasured — `MNO_KMESH=2` is ALL-TRIM, so its whole k-sum runs real.
- Report fold state, stream coverage (`pts64/pts32/dropped`, `meanRun`), thread count and BLAS mode
  alongside each number, or the table will not be reproducible.

### Step 2 — ARM THE SYMMETRY FOLDS  ·  plan: `doc/SymmetryUpgradePlan.md`

**The biggest single runtime multiplier left, by an order of magnitude** — this is Step 3's work done
properly, not a separate track.  Inventory (2026-08-15, re-verified 2026-08-19):

- the {G}-star fold is wired at exactly TWO static sites (the local-PP sweeps, imposed runs only);
- the **per-iteration G-space consumers are UNFOLDED** — ρ̃, the Poisson multiply, the V_xc gathers, the
  G_ERI3 columns, the seed structure factors: **12–24× on the MnO magnetic group, 48× cubic**, unclaimed;
- the REAL-space T3 pair-stream orbit fold (**71× measured on the cache**) is BUILT but still opt-in
  (`GPW_STREAM_FOLD=1`; the open-shell D-asymmetry slosh retraction is why).

Work: default-arm T3 where safe + auto-arm at multi-k (T3.4b), extend the ball fold to the per-iteration
sites, and report every fold factor (Step 0b).  Caveats: the FFT itself does not fold trivially; the
dominant per-iteration cost is real-space, so T3 is the big one.  For scale: rounds 3–4 bought 1.4× by
hand-tuning kernels.

### Step 3 — RUNTIME, CONTINUED  ·  plan: `doc/GPWPlan1.md`

Measure against Step 1, after Step 2 (folding changes what is hot).

- **Φ-table SCREENING — the O(N)-XC increment, and THE next structural lever.**  Φ is stored dense
  (npts × n) but the true object is SPARSE: batch the mesh (per atom / per radial shell) with a per-batch
  significant-function list and every Φ-shaped cost (build, ρ GEMM, H_xc GEMM) becomes O(npts·n_sig²).
  **The win grows with cell size** — MnO's 4 atoms understate it.
- **Becke mesh build** (~12–32 s): `BeckeCutoff` alone was 11.5% of the round-3 profile.
- **NOT worth re-attacking without new information:** the per-iteration scatter/gather.  Round 4 measured
  it at 61% irreducible per-(pair,point) emit; it needs fewer TERMS (looser `GPW_DENSITY_EPS`, the T3
  fold, a smaller span), not a faster kernel.

### Step 4 — RAM  ·  read it off Step 1's table, then decide

Not a campaign yet.  Stream RAM already fell 5.78 → 3.70 GB (round 3's run-length encoding) and Step 2's
folding cuts stream *demand* by the orbit factor.  CP2K caches nothing — its kernels are just fast — so
our answer is "fold, then re-tier the budgets", which is Steps 2–3 with a different readout.  Open it as
its own track only if Step 1 still says it binds.

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

**Next concrete move (cheap, does NOT wait on Steps 1–4):** the term-by-term breakdown —
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
