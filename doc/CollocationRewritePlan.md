# The collocation rewrite — from a per-point walk to a separable contraction

**Goal: close the remaining ~22× against CP2K on the on-the-fly path**, and take the 4218 MB pair-stream
cache down with the same change.  The measured record that motivates it is in `doc/OpenWork.md`
("THE ON-THE-FLY BOX WALK") and `doc/Benchmark.md`.

> ✅ **COMPLETE 2026-08-27 (steps 0–8).**  Both halves of the goal landed: the cache is DELETED (3915 → 155 MB
> on the MnO probe, 1323 → 463 MB on the imposed benchmark row) and `GPW_CONTRACT_CUBE` is the DEFAULT
> collocation kernel, with a (shell pair, offset) TASK LIST in the cache's place.  MnO's box walk went
> 1139.8 → 78.9 s (**14.5×** over the campaign); against CP2K the standing is now 2.6× CPU / 2.1× RAM on
> MnO, and Si Γ and NaF full-SR BEAT CP2K on CPU outright.  792/792 on both kernel settings.
> **This file is now a RECORD, with one forward item: §3c-bis, the \f$(L_a,L_b)\f$ batching.**

## 0. Vocabulary and the measurement rule (both learned the hard way, 2026-08-26/27)

- **"BOX WALK", never "kernel"** (user).  *Kernel* is overloaded — M&D integral kernels, XC kernels, GPU
  kernels.  The box walk is `NR_Evaluator::ForShellPairBox` plus the two lambdas that consume it,
  `scatterShell` (in `CollocateDensity`) and `integrateShell` (in `IntegratePotential`): the code that,
  for one shell pair at one cross-cell offset, walks the grid points of the product's compact box.
- **"DENSITY-MATRIX SUB-BLOCK", never "D-block"** (user).  In a codebase whose hardest case is Mn *d
  shells*, "D-block" reads as angular momentum.  \f$D_{ij}\f$ over \f$i\in\f$ shell I, \f$j\in\f$ shell J.
- **ORTHO vs NON-ORTHO METRIC** — and the metric is CELL geometry, not a crystal-system label
  (user, 2026-08-27, correcting an earlier claim here that it was "unrelated to the crystal system").
  The grid metric IS the cell's metric tensor, rescaled per axis:
  \f[ h_{ab}=s_a\!\cdot\!s_b=\frac{G_{ab}}{N_aN_b},\qquad G_{ab}=\mathbf a_a\!\cdot\!\mathbf a_b,
      \qquad s_a=A\hat e_a/N_a \f]
  and \f$G\f$ is exactly what fixes the cell's lengths and angles.  So it has EVERYTHING to do with
  crystallography.  What comes apart is this: **the crystal SYSTEM classifies the LATTICE, while the metric
  belongs to the CELL chosen to represent it.**

  | cell as rastered | \f$|a|\f$ | \f$\alpha=\beta=\gamma\f$ | \f$G\f$ off-diag | volume |
  |---|---|---|---|---|
  | Si / NaF `FCCUnitCell(10.26)` — PRIMITIVE | 7.2549 | **60.000°** | 26.32 | 270.0 |
  | ...the same lattice's CONVENTIONAL cell | 10.2600 | 90.000° | 0.00 | 1080.0 |
  | MnO AFM-II \f$(a,\tfrac a2,\tfrac a2;\ldots)\f$, \f$a=8.40\f$ | 10.2879 | **33.557°** | 88.20 | 296.4 |
  | ...MnO's own FCC primitive, for the ratio | 5.9397 | 60.000° | 17.64 | 148.2 |

  ⇒ **ALL THREE PRODUCTION CELLS ARE RHOMBOHEDRAL** — Si and NaF at 60° (the FCC primitive of a CUBIC
  lattice), MnO at 33.557° (and its volume 296.35 is exactly 2× the FCC primitive's 148.18, i.e. the AFM-II
  magnetic doubling).  "Non-ortho metric" here means *rhombohedral*, specifically.  The only ortho-metric
  grids in the suite are the atom-in-box tests (Na, O₂, Mn₂, Na₂).

  ★ **AND THE ORTHO FAST PATH IS PURCHASABLE — at 2×, not 4×** (corrected 2026-08-27; the first cut of
  this paragraph quoted the CONVENTIONAL cubic cell).  The smallest ORTHOGONAL-metric cell of the FCC
  lattice is **body-centred tetragonal at 2× the primitive volume**: \f$a(0,0,1)\f$,
  \f$\tfrac a2(1,-1,0)\f$, \f$\tfrac a2(-1,-1,0)\f$, lengths \f$a,\tfrac a{\sqrt2},\tfrac a{\sqrt2}\f$
  (enumerated, not recalled).  So the trade is 2× the grid points to buy a path measured ~1.5–2× faster per
  point — roughly BREAK-EVEN, not the clear loss first claimed here.  The conclusion survives, but only
  just.  ⇒ We stay in the general case, and the non-ortho kernel is the target rather than a fallback.
  (CP2K calls the two paths `ortho` and `general`.)

  ⛔ **AND FOR MnO IT IS NOT A CHOICE AT ALL.**  Its cell is rhombohedral because the MAGNETIC ORDER is —
  AFM-II stacks along [111], reducing the lattice to trigonal — so no low-index orthogonal cell exists.
  That settles it independently of the Si arithmetic.

### 2b. ★ THE PRIMITIVE-vs-CONVENTIONAL DICHOTOMY — the symmetry half is a NON-issue (user, 2026-08-27)

> *"Primitive: less atoms, fake low symmetry. FCC: 4x the atoms, proper full symmetry.  There must be some
> (group theory) technique that allows you to work with less atoms and maintain the full symmetry…"*

**The premise does not hold, and that is the good news: the primitive cell HAS the full symmetry.**  All
**48/48** cubic point-group operations are INTEGER matrices in the primitive basis (\f$P^{-1}WP\f$ integral
for every \f$W\in O_h\f$ — enumerated).  Only the cell's SHAPE is rhombohedral; the LATTICE is cubic and
every operation maps it to itself.  The code already relies on this: the Si Γ run prints
`[IBZ] point group |ops|=48` and `[stream fold] 48/48 ops folded`, on the FCC PRIMITIVE cell.

**And the sought technique IS the primitive cell.**  FCC's four "extra" centring translations are not
non-lattice — they ARE lattice translations, and look extra only because the conventional description
carries four lattice points.  The primitive cell is precisely the quotient by them.  So *fewer atoms +
full symmetry* is already achieved; what it costs is the diagonal metric, which is geometry, not group
theory.  (SALC is a different axis: it reduces the BASIS, block-diagonalising H and S by irrep — see
doc/MolecularSymmetryPlan.md — not the grid.)

⚠ **ONE IDEA THAT MIGHT ESCAPE THE TRADE — parked, not scheduled.**  Decouple the GRID from the CELL: use
cubic voxels commensurate with the CONVENTIONAL cell (diagonal metric) but store only a primitive cell's
worth, wrapping through the centring translations.  Metric diagonal, storage primitive.  The catch is the
Poisson solve — the FFT domain stops being a simple box, and an FCC-periodic function has nonzero Fourier
components only on the BCC reciprocal sublattice (3/4 of the conventional box's G-vectors are identically
zero).  Real, non-trivial, and orthogonal to this plan.
- **EVERY speed number is the sum of exactly two timing-ledger buckets** — `scf: integrate-back (pair
  gather)` + `scf: collocate density (pair scatter)`.  They are EXCLUSIVE and DISJOINT, so no setup
  subtraction ever appears.  `GPW_REPORT=1` prints them.

## 1. ★★★ THE TESTING-LEVEL PIN (user, 2026-08-27) — this plan's most important process constraint

> *"this project probably has too much (expensive) testing at the integration level, and not enough at the
> (cheaper, better isolated) unit level."*

**The box walk today has ZERO unit tests.**  Grep `src/*/tests/` for `ForShellPairBox`, `ForPairBox`,
`CollocateDensity`, `IntegratePotential`: nothing.  Every scrap of coverage for the hottest,
most load-bearing loop in the project arrives through a full SCF in `IntegrationTests/`.

**What that cost, concretely.**  The four box-walk edits of 2026-08-26 were validated only end-to-end, so
each edit-measure cycle was a 3-minute NaF SCF or a 9-minute MnO probe, and the evidence for
"bit-identical" was *"the printed energy did not move"* — an end-to-end proxy for a POINTWISE property.
A unit-level oracle would have been faster AND strictly stronger.

⇒ **STEPS 0–3 BELOW ARE UNIT TESTS**, in `src/BasisSet/Molecule/tests/`, added to the EXISTING
`UTMolecule_BS` exe (a new source file in its `CMakeLists.txt` list — deliberately NOT a new exe, so the
`allTests` DEPENDS trap in CLAUDE.md does not apply).  Integration tests enter only at step 8, as
ACCEPTANCE, not as the development loop.

## 2. The algebra — and we already own it

The M&D collapse CP2K rebuilds per task (`cab_to_cxyz`) exists here and is CACHED:
`src/BasisSet/Molecule/Evaluators/PG_Cart_MnD/GaussianRF.C`, `struct Ω : public Cacheable2`, interned in
the process-global `Cache2` on the primitive pair.  It carries \f$\alpha_p=a+b\f$ (CP2K's `zetp`), the
product centre \f$P\f$ (`rp`), the prefactor `Eij` (`prefactor`), and `H2`.

\f[ \chi_a(r)\,\chi_b(r-R) \;=\; E_{ij}\sum_{NLM} E^x_{N n_a n_b} E^y_{L l_a l_b} E^z_{M m_a m_b}\;
    \Lambda_{NLM}(r-P;\alpha_p) \f]

and `Hermite2` stores those coefficients as `d`,`e`,`f` indexed \f$(N,n_a,n_b)\f$, \f$(L,l_a,l_b)\f$,
\f$(M,m_a,m_b)\f$ — **already per-direction separable**, which is exactly the property that makes
\f$\Lambda_{NLM}(r-P)\f$ factor into three 1-D functions of one Cartesian coordinate each.

⇒ **Our starting point is BETTER than CP2K's, not worse.**  They re-derive the collapse per task; we have
it cached as geometry.

### 2a. Do NOT copy CP2K's `pab` fold — \f$D\f$ is the only per-iteration object here

\f$\Omega\f$, \f$\alpha_p\f$, \f$P\f$, `Eij`, `H2`, the boxes and the chords are ALL GEOMETRY.  CP2K folds
the density matrix into `cab` before building the cube because it rebuilds the cube every task anyway.  We
do not, so the split is:

| object | depends on | rebuilt |
|---|---|---|
| separable 1-D/2-D tables | \f$\Omega\f$ + the grid | **never** (geometry; valid across iterations AND k-blocks) |
| the coefficient tensor \f$c_{NLM}\f$ | \f$\Omega\f$ + \f$D\f$ | per iteration — \f$O(n_I n_J n_{NLM})\f$ per pair, **no grid points** |
| the grid contraction | both | per iteration |

Three things this settles:
1. **The no-cut discipline and the D-aware screen need not move for \f$D\f$'s sake at all** — \f$D\f$ never
   enters the geometry object.
2. The weight is \f$\mathrm{fold}\cdot\mathrm{Re}[D_{ij}\overline{e^{ikR_n}}]\f$: a REAL scalar per (pair,
   **offset**), because the Bloch phase rides on the image.  Any weighting is per-offset.
3. **The integrate direction has no \f$D\f$ whatsoever** — it PRODUCES \f$h_{ij}\f$; `screenD` there is a
   screening magnitude, not a weight.

⚠ \f$D\f$ is over CONTRACTED, normalized components (`ns[i]`); \f$\Omega\f$ is over PRIMITIVE pairs.  Every
current basis is uncontracted so they coincide TODAY — `gi[p]` exists for a reason; do not assume it.

## 3. The steps

### ✅ Step 0 — THE GATE: **PASSED 2026-08-27, 6.4× on the MnO-representative cube**
`src/BasisSet/Molecule/tests/M_PG_BoxWalk.C`, in `UTMolecule_BS`.  Cubic 16 a.u. cell, 96³ grid, two
centres 4.0 a.u. apart (the MnO Mn–O contact), \f$\varepsilon=10^{-10}\f$.  The contraction pays its OWN
setup every repetition — the collapse AND the tables — so this is the honest per-(pair, offset) cost, not
a table-reuse flatter.

| case | walk | contract | ratio |
|---|---|---|---|
| **Mn-like d × O-like p (the MnO contact)** | 20.7 ms | 3.25 ms | **6.36×** |
| diffuse d × diffuse p (the big box) | 211.2 ms | 25.0 ms | **8.43×** |
| d × d | 4.92 ms | 0.76 ms | **6.48×** |
| s × s (the cheap end) | 12.5 ms | 2.97 ms | 4.23× |

⇒ **Above the ≥5× bar on every case with \f$L>0\f$.**  The s×s row is expected and is not a
counter-example: a 1-FMA innermost loop has no structural win to show, and s×s pairs are not where the
time goes.  ★ The DIFFUSE case gains most (8.4×), which is the right way round — diffuse pairs own the
biggest boxes and therefore most of the cost.
⛔ **AND THIS ROW IS NOT THE ANSWER: it is an ORTHO-METRIC grid, which NO production cell has.**  Quote
the non-ortho gate below instead.

<details><summary>the original statement of the gate</summary>
Time a single representative (shell pair, offset) from MnO (an Mn *d* × O *p* pair at a real level) two
ways with no SCF anywhere: today's `ForShellPairBox`, and a hand-written separable contraction over the
same cube.

⛔ **IF IT IS NOT ≥5× AT THE CUBE LEVEL, STOP.**  The walk is ~90% of the uncached run, so anything less
cannot survive Amdahl.

★ **WHY THIS STEP EXISTS.**  The exp-recurrence experiment (doc/OpenWork.md) removed 20% of the profile
and bought 3%, after a full gated implementation.  A cube-level number would have said so in an hour.
**A hot symbol's share is an upper bound on the win, not an estimate of it** — and the only cheap way to
tell the difference is to build the replacement kernel in isolation FIRST.
</details>

### ✅ Step 0b — THE GATE THAT COUNTS: **non-ortho metric, 3.5–7×** (2026-08-27)
`SkewSeparableContractionGate`, on the REAL metrics: Si `FCCUnitCell(10.26)` at 24³ and MnO
\f$(8.4,4.2,4.2;\ldots)\f$ at 32³, centres at the true contact distances (Si–Si 4.44 a.u., Mn–O 4.2 a.u.).
Mathieu's three 2-D exponential tables + the polynomial re-expanded in grid-index powers; innermost loop
\f$l_p{+}1\f$ FMAs and three table lookups.  Ratios over three runs (absolute times swing ~2× on this
shared box; the ratios hold to ~20%):

| non-ortho case | ratio |
|---|---|
| Si FCC primitive, d × p | 2.6–3.6× |
| MnO rhombohedral, d × p | 3.6–4.5× |
| **MnO rhombohedral, DIFFUSE d × p** | **6.3–7.2×** |
| MnO rhombohedral, d × d | 4.2–4.6× |

⇒ **Lower than the ortho row's 6.4–8.4×, as expected**: with an ortho metric the exponential tables are
\f$O(n)\f$ (1-D); with a non-ortho one they are \f$O(n^2)\f$ (three 2-D tables), so the transcendental
saving drops from \f$\sim n^2\f$ to \f$\sim n\f$.
★ **The DIFFUSE pair gains most (6–7×), and it is the one that matters** — its cube costs 114–128 ms
against 16–34 ms for the mid pairs, so it dominates the cost-weighted answer.
⇒ **VERDICT: the gate is cleared**, but on a smaller margin than the ortho number advertised, and the
projected 5–15× in §4 should be read as 4–7× until measured in situ.

★ **KNOWN HEADROOM, and it resurrects the parked branch AT THE RIGHT LEVEL.**  The prototype builds its
2-D tables with \f$3n^2\f$ DIRECT `exp` calls; CP2K builds the same tables with the multiplicative
recurrence (`grid_cpu_collint.h:752`, seeded symmetrically outward from the cube centre).  That is exactly
the identity rejected in doc/OpenWork.md — rejected there because it sat inside an \f$O(n^3)\f$ point
loop and was z-anisotropic.  On an \f$O(n^2)\f$ TABLE it is neither.  So these ratios are a FLOOR.

⚠ **THREE PROTOTYPE DEFECTS FOUND AND FIXED WHILE GETTING THIS NUMBER — all of them made the algorithm
look worse than it is, and the first two would have been quoted as fact:**
1. `IndexPoly::MulLinear` zeroed and copied its full \f$9^3=729\f$ slots ~375× per cube; degree-bounding
   every loop was worth ~2×.
2. The innermost loop called a helper doing THREE modulo operations per point — the identical defect that
   cost 13.9% in the production walk, reintroduced in the thing meant to measure it.
3. Two MnO cases had the centres 7.27 a.u. apart, where the prefactor screen kills the whole box: the walk
   timed **0.00 ms** and the ratio read 0.00×.  A DEGENERATE case, not a fast one.  The gate now asserts
   the cube is non-trivial before timing it.
⇒ **A prototype's own overhead is the easiest way to measure the wrong thing.**  First reading was 2.6×
across the board; after the fixes, 3.5–7×.

### ✅ Step 1 — THE ORACLE: DONE 2026-08-27 (two tests, same file)
`SeparableCollapseMatchesTheWalkPointwise` sweeps \f$L_a,L_b\in\{0,1,2\}\f$ × three \f$\alpha_I\f$ × two
\f$\alpha_J\f$ (36 shell pairs) and checks the NAIVE collapse against the walk at every point the walk
visits — so a failure convicts the ALGEBRA, not the contraction.  `ContractionMatchesTheWalkOverTheCube`
then checks the fast three-step contraction over the whole cube.

⚠ **AND THE ORACLE'S FIRST CUT WAS WRONG IN A WAY WORTH RECORDING.**  It used a RELATIVE per-point error
and reported deviations of \f$10^{12}\f$ on a correct collapse.  A shell pair summed over components with
mixed-sign weights passes through ZERO at points inside the cube (a d shell's \f$x^2\f$-like combinations
do it routinely), so the relative criterion divides by ~0 and converts a \f$10^{-16}\f$ error into a
catastrophic-looking one.  **The collocated value's meaning is absolute** — the walk screens it at
absolute \f$\varepsilon\f$ — so the tolerance is absolute, scaled to the cube's own peak.  Same family of
mistake as the drift audit that could not see its own worst case: *the criterion has to match what the
quantity means.*

<details><summary>the original statement of the step</summary>
`ForShellPairBox` STAYS in the tree as the reference implementation.  A new gtest asserts the new cube
equals the old walk POINTWISE across a case matrix — orthorhombic AND triclinic cells, \f$L=0\ldots3\f$,
near and far offsets, sharp × diffuse pairs, boxes that wrap the cell.  Written and green against the old
path first, so it is a real instrument when the kernel lands rather than a retrofit to whatever the kernel
happens to do.
⚠ **The oracle must be pointwise**, not an integrated total: an integral hides sign-cancelling errors, and
this is exactly the "integrated observables" rule's dual — for a CORRECTNESS oracle on a field, the field
is the observable.
</details>

### ✅ Step 2 — RESOLVED 2026-08-27: **binomial monomials about P, with \f$\alpha_p,P,E_{ij}\f$ from the cached \f$\Omega\f$**

⛔ **THE CRITERION THIS PLAN PROPOSED DOES NOT DISCRIMINATE.**  It said to choose on "which keeps
integrate-back an exact adjoint".  Adjointness follows from using the SAME coefficients in both directions
with the contraction reversed — available in ANY polynomial basis.  `ContractedCubeReproducesThePairIntegrals`
measures it: against an arbitrary pseudo-random field \f$V\f$, every component pair's
\f$\sum_g V\chi_i\chi_j\f$ from the contracted cube matches the walk's to **2.1e-19 on a scale of 3.2e-5**
(relative 7e-15, i.e. machine precision), and to 9.5e-13 including the boundary points the contraction adds
— under the 1e-10 collocation \f$\varepsilon\f$.  ⇒ Adjointness is available either way; decide on other
grounds.

**THE GROUNDS THAT DO DECIDE:**
1. **On the grid, Hermite buys nothing.**  \f$\Lambda_{NLM}\f$ is *polynomial × the same Gaussian*, so it is
   a change of POLYNOMIAL BASIS, not of structure.  And the non-ortho path re-expands the polynomial in
   GRID-INDEX powers regardless — so starting from Hermite would add a Hermite→monomial conversion that
   the binomial route never needs.
2. **\f$\Omega\f$ is reused for the part that matters anyway** — \f$\alpha_p\f$, \f$P\f$, \f$E_{ij}\f$ and
   the `Cache2` identity are the structural, cached quantities; the per-axis binomial is ~10 lines of
   geometry-only arithmetic.
3. The binomial route is what CP2K runs, so it is proven at scale.

⚠ The \f$\Omega\f$/`H2` convention was CONFIRMED from the code, not assumed:
`Overlap2C = (\pi/\alpha_p)^{3/2} E_{ij} H2(0,p_a,p_b)` (`Imp/GaussianRF.C:335`) pins
\f$\Lambda_{000}=e^{-p|r-P|^2}\f$ and \f$E_{ij}=e^{-ab|AB|^2/p}\f$ — the standard M&D convention, and the
one the prototype uses.

<details><summary>the original statement of the step</summary>
Hermite (reuse `H2`'s `d`/`e`/`f`) or binomial re-expansion to monomials about \f$P\f$ (CP2K's route).
**The criterion is NOT arithmetic count.**  It is which one makes integrate-back come out as the EXACT
ADJOINT with the same coefficients: \f$\int\rho V = \mathrm{Tr}(Dh)\f$ to machine precision is load-bearing
(it is what makes the GPW energy variational), and it is far easier to preserve than to repair.
</details>

### ✅ Step 3 — RESOLVED 2026-08-27: **the screen DISSOLVES.  It was not the blocker.**

**THE OBSERVATION THAT SETTLES IT, and it is in the walk's own code:** the BOX is already sized by the
UNION tolerance (`ForShellPairBox` is handed `epsUnion`).  So the per-component `epsHere[k]` test does NOT
shrink the walk — it only declines to ACCUMULATE a value it has already computed.  A contracted cube has
nothing to decline: the cube is formed for the whole shell pair regardless.
⇒ **Dropping the per-point component screen costs NO work and REMOVES a truncation** — the contracted form
is strictly MORE accurate, not less.

Measured (`PerComponentScreenDropsSubEpsOnly`, weights spread over four decades so the per-component
tolerances differ widely — the case the screen exists for): the screen drops **48322 of 60984 (component
pair, point) terms — 79% of them** — for a worst pointwise loss of **4.4e-10** against a cube peak of
5.1e-5.  That is 4.4× \f$\varepsilon\f$, inside the \f$n_In_J\varepsilon\f$ bound a magnitude screen
guarantees.  ⇒ It drops a great deal and costs almost nothing either way.

**WHAT SURVIVES, and it is where the real saving always was:** the per-(pair, offset) PREFACTOR pre-filter
(`PairPrefactorExp` ≥ −ln ε kills the whole term).  It touches no grid points, it is what shrank the
D-aware boxes, and it carries over to the contracted form unchanged.

<details><summary>the original statement of the step</summary>
Today's D-aware tolerance is PER COMPONENT PAIR (`epsHere[k]`, applied to \f$|val|\f$ at each point), and
pre-filtering on it is what made the shell hoist pay at all — measured 1.03–1.08× without it against
2.13× with.  A contracted cube has no per-component identity: one cube serves the whole shell pair.
- **Try the cheap option first:** keep the per-component pre-filter deciding only WHICH PAIRS ENTER the
  contraction.  It already runs once per offset with no grid points involved, so it may carry over
  unchanged.
- **Only if that fails:** derive a bound on \f$|c_{NLM}|\f$ and screen the coefficient tensor.
</details>

### ✅ Step 3c — THE TASK LIST: **yes, build it — but NOT for the reason it looks like** (user, 2026-08-27)

> *"Since we are planning to re-evaluate at every iteration rather than cache, would it make sense to build
> a list of i,j pairs (or i,j,R triples) that pass the screen?"*

**This is exactly CP2K's TASK LIST** (`cp2k/src/grid/grid_task_list.h`): a task is
`(level, iatom, jatom, iset, jset, ipgf, jpgf, border_mask, block_num, radius, rab)` — the (shell pair,
primitive pair, image) triple with the reach precomputed and the level resolved.  Same object, same reason.

⛔ **BUT THE OBVIOUS JUSTIFICATION IS FALSE, and measuring it first is what makes this worth doing.**
`OffsetEnumerationIsNotTheCost`, on a 4-atom MnO-shaped cell, 36 shells, 32³ grid:

| | |
|---|---|
| enumerate + screen all 4482 (shell pair, offset) tasks | **2.26 ms** |
| actually walk them (102,450,921 grid points) | **2294 ms** |
| ⇒ enumeration's share of the pass | **0.10 %** |

So hoisting the enumeration saves **nothing**.  `ForImageOffsets` rebuilds a `CellsInSphere` vector per
shell pair per call and it does not matter.  The same goes for the per-task geometry (reach, \f$P\f$,
\f$\alpha_p\f$, \f$E_{ij}\f$, the fractional bounding box): ~100 flops against ~100 M grid points.

★ **WHAT THE TASK LIST IS ACTUALLY FOR — and it is the whole point of the rewrite:**

| | measured on the fixture above |
|---|---|
| the equivalent VALUE cache (what we store today) | **781.6 MB** |
| a TASK LIST over the same set (48 B/task) | **0.205 MB** |
| ratio | **3810×** |

and on the REAL NaF counts (406 pairs, 5514 offsets, 57,444,703 cached points): 459 MB of values against
265 KB of tasks, **~1700×**.

⇒ **The task list is what makes "re-evaluate every iteration" a DESIGN rather than a regression.**  Today
"uncached" means re-deriving everything and paying 4218 MB to avoid it.  With a task list the GEOMETRY is
derived once into a few MB, and only the ARITHMETIC repeats — which the rewrite makes 4–7× cheaper.  That
is precisely CP2K's posture, and it is the first configuration that beats them on BOTH axes at once.

**The other three reasons, in honest order (none of them enumeration):**
1. **BATCHING BY \f$(L_a,L_b)\f$.**  Sorting the list groups tasks with the same total degree, so the
   contraction's innermost `for a<=lp` runs with \f$l_p\f$ known at compile time.  CP2K dispatches exactly
   this way (`grid_cpu_collocate.c`, `always_inline` low-\f$l_p\f$ variants).  ⚠ Unmeasured here.
2. **A FLAT, SORTABLE PARALLEL LOOP.**  Today there are two nested `schedule(dynamic)` loops (pairs, then
   shell pairs); one list sorted longest-first balances better and is trivially chunkable.
3. **THE D-AWARE FILTER BECOMES AN \f$O(1)\f$ PER-TASK PREDICATE**, cleanly outside the walk — and the
   static list is a strict SUPERSET, since \f$\varepsilon_{ij}=\varepsilon/|c_{ij}|\ge\varepsilon\f$ means
   a per-iteration \f$D\f$ can only ever remove tasks, never add one.  So the list is built once at
   \f$\varepsilon\f$ and filtered per iteration.

⚠ **Do NOT put \f$D\f$ in the task list** (§2a): the list is geometry, valid across iterations AND
k-blocks; the density is a per-iteration filter and a per-iteration coefficient contraction.

### ▶ 3c-bis — DESIGN NOTES FOR REASON 1, the \f$(L_a,L_b)\f$ batching (specced 2026-08-27, NOT built)

Reason 1 above is now the next lever on this path (step 7's own closing line), and the questions it raises
were worked through with the user the day step 7 landed.  Recorded here so the next session inherits the
answers instead of re-deriving them.

**★ THE TEMPLATE PARAMETER IS \f$l_p\f$ ON THE GRID SIDE AND \f$(L_a,L_b)\f$ ON THE SETUP SIDE — because
`MakePairPoly` is a COLLAPSE, and the two halves sit on opposite sides of it.**

- **`ContractCube` / `GatherCube` / `ToGridPoly` see only \f$l_p\f$.**  By the time they run, all
  \f$n_I\times n_J\f$ component pairs have been summed — with their density weights — into ONE coefficient
  tensor `q.c[sx][sy][sz]`; there is no loop bound left in them that distinguishes d×s from p×p.  Nor does
  \f$(L_a,L_b)\f$ tighten the tensor: a Cartesian shell of total \f$L\f$ contains components with any single
  exponent from 0 to \f$L\f$, so \f$s_x\f$ reaches \f$L_a{+}L_b\f$ on its own and the tensor fills the full
  simplex \f$s_x{+}s_y{+}s_z\le l_p\f$.  ⇒ `template<int LP>`, 9 instantiations (0..`kMaxPoly`−1).
  Templating these on \f$(L_a,L_b)\f$ would emit 25 copies of 9 distinct functions — strictly worse, and
  worse for the batching argument, which wants ONE resident specialization at a time.
- **`MakePairPoly` / `MomentsToPairs` see \f$(L_a,L_b)\f$ separately.**  \f$n_I\f$, \f$n_J\f$ become
  compile-time \f$(L{+}1)(L{+}2)/2\f$, the per-component monomial exponents become constants, and
  `Binom1D`'s degree with them.  ★ The concrete win there is FOOTPRINT, not flops: `wd[kMaxShell*kMaxShell]`
  is **16 KB of stack per task** to hold the 36 doubles a d×d pair needs (same for `bAll`), `GridPoly` is
  \f$9^3\f$ = 729 doubles where the live simplex at \f$l_p{=}4\f$ is 35, and `MulLinear` allocates a SECOND
  one per call.  All of that collapses to tens of bytes.
- **They compose — there is no choice to make.**  A batch dispatch keyed on \f$(L_a,L_b)\f$ knows both, so
  one `switch` yields `MakePairPolyN<LA,LB>` and `ContractCubeN<LA+LB>` from the same arm.

**⛔ AND IT SPLITS INTO A BIT-IDENTICAL HALF AND AN ANCHOR-MOVING HALF.  This decides the SCHEDULING.**

| increment | bit-identical? | when it can land |
|---|---|---|
| **1. unroll with `LP` fixed, KEEPING the `E1` table** | **YES** — same terms, same summation order, just no loop | any time, gated on the existing anchors |
| **2. Horner, i.e. DROP `E1`** | **NO** | only inside an anchor re-bank |

`E1[a]` is built by repeated multiplication (`e *= e1`) and summed low-to-high; Horner re-associates the
whole polynomial, so every GPW anchor moves again.  Horner is the more attractive half — it deletes an
entire \f$(l_p{+}1)\times n_0\f$ table, its build loop and its memory traffic from the innermost path — but
it belongs in the SAME re-bank as A2–A6 on `doc/OpenWork.md`'s anchor-moving roster, or it makes exactly
the masking problem that roster exists to prevent.  ⇒ **Two increments, not one.**

**⛔ §5 OPEN QUESTION 2 CONSTRAINS THE BATCH KEY.**  `MakePairPoly` declines contracted radials today and
falls through to the walk.  If that question is ever answered CP2K's way, a task becomes
**(shell pair, primitive pair, offset)** and §5 question 3 moves `PairLevel` with it from a shell-pair to a
primitive-pair property.  ⇒ Key the batch on \f$(L_a,L_b)\f$, but do NOT hard-wire "one cube per shell pair"
into the dispatch interface; let the TASK identity stay opaque.

**⛔ THE GENERIC ARM MUST SURVIVE.**  The kernel declines \f$l_p\ge\f$ `kMaxPoly` and contracted shells, and
since 2026-08-27 `GPW_CONTRACT_CUBE=0` (the reference walk) is the arm nothing exercises by default.  A
templated dispatch keeps the generic path as its default arm; it does not replace it.

**✅ §4b IS ADJACENT BUT DISJOINT.**  The parked \f$\Omega\f$ fold is to be revisited "AFTER step 7, when the
setup buckets are a larger share of what is left", which is the same argument as the \f$(L_a,L_b)\f$ SETUP
half — but they do not overlap: §4b measured that collocation and integrate-back never touch \f$\Omega\f$ at
all.  Two independent levers on the same rising share.

**▶ FIRST ACTION, and it is a MEASUREMENT (§1).**  Time `MakePairPoly` + `ToGridPoly` against `ContractCube`
over the existing MnO-shaped fixture in `M_PG_BoxWalk.C` — the same shape as `OffsetEnumerationIsNotTheCost`
— and let it decide whether the setup half is worth templating at all.  ⚠ The envelope says a few percent,
and this campaign's envelopes have been wrong twice in the same direction: enumeration was "obviously worth
hoisting" and measured **0.10%**, the \f$\Omega\f$ fold was elegant and measured **0.6% of an MnO run**.
Measure before templating anything on \f$(L_a,L_b)\f$; `LP` on the grid side needs no such permission.

### Step 4 — The ORTHO-METRIC kernel — ⚠ DEMOTED: a stepping stone that ships NOTHING
Three 1-D tables, three-step contraction.  ⛔ **No production cell has an ortho metric** — only the
atom-in-box tests do.  So this buys simpler bring-up of the contraction machinery and nothing else; do it
only if step 5 proves hard to land in one go, and never quote its numbers as the result.

### ✅ Steps 5 + 6 — BUILT AND MEASURED 2026-08-27.  **4.98× on the MnO acceptance case.**
Gated on `GPW_CONTRACT_CUBE`, **default OFF** (see the verdict below).  `ForShellPairBox` remains the
reference implementation; `BoxGeom`/`MakeBoxGeom` is shared by both so they cannot drift.

| box walk (two ledger buckets), streams disabled | walk | contract | |
|---|---|---|---|
| Si Γ | 1.474 s | 0.372 s | **3.96×** |
| NaF Γ | 140.3 s | 38.9 s | **3.61×** |
| **MnO AFM-II Γ (acceptance)** | **392.7 s** | **78.9 s** | **4.98×** |

MnO wall 429 → 113 s, CPU 470 → 154 s, **peak RSS unchanged** (166.2 → 168.5 MB), trajectory identical
(iters, lastΔρ, Eamp all unchanged), `Etot` −59.69580383 → −59.69580301 (8.2e-7 Ha).

⇒ **SESSION TOTAL on the MnO box walk: 1139.8 → 78.9 s = 14.5×**, i.e. 569.9 → 39.4 s per iteration.
Against CP2K's 8.5 s/iteration the standing goes **67× → ~5×**.

**ADJOINTNESS IS THE ACCEPTANCE TEST, and it passes** (`CollocateIntegrateAreExactAdjoints`).
\f$|\int\rho V-\mathrm{Tr}(Dh)|\f$ relative: walk/walk 5.5e-15, **contract/contract 4.4e-15**, mixed
1.4e-14.  `GatherCube` is the literal transpose of `ContractCube` — same tables, same chord, same
coefficients — so this holds by construction, not by tuning.

⛔ **THE ONE REAL BUG, and it was mine, not the method's.**  The chord bracket used
`floor(root1)-1 .. ceil(root2)+1`.  The correct bracket is `ceil(root1) .. floor(root2)` — a point is kept
iff the quadratic is ≤0 there.  The extra points all carry the SAME SIGN for a given pair, so they did not
cancel in an integral:

| | over-wide bracket | exact bracket |
|---|---|---|
| pointwise relative error | 2.4e-11 | **2.4e-14** |
| INTEGRATED relative error | ~5e-10 | **~7e-15** |

★ This retracts an earlier claim in this file that the deficit was *"intrinsic to the non-ortho route"*
(index-unit polynomials times cancelling 2-D tables).  It was not intrinsic; it was a sloppy bracket.
⇒ **MEASURE THE INTEGRATED DIFFERENCE, NOT ONLY THE POINTWISE ONE.**  Systematic pointwise errors do not
cancel, and the SCF sees the integral.  The pointwise number looked acceptable at every stage while the
integrated one was 30× worse.

**WHAT STILL DIFFERS between the two paths is the per-component screen** (step 3), which the contracted
cube has no analogue for and needs none.  Measured integrated cost of that screen: 1.5e-7 absolute, 5.8e-5
relative on a pair whose weights span four decades.  ⇒ **THE CONTRACTION IS THE MORE ACCURATE OF THE TWO**;
Si's residual 1.8e-7 Ha is the WALK's truncation, not the kernel's.

### ⚠ THE DEFAULT: still OFF, on one piece of evidence
`ctest` is **783/783 with the flag off** and **782/783 with it on**.  The single failure is
`GPW_SCF.ImposedOrderLostIsAPostconditionFailure_Na2Box`, whose own comment records that its recipe is
marginal — *"alpha=0.5 and a generous cap are LOAD-BEARING, not decoration"*, converging in **66 of a
100-iteration cap**.  The kernel's perturbation pushes it past the cap; its ORDER behaviour is still
correct (the moment dies at step 65, which is what the gate is about), but `DidConverge()` is asserted
first.
⇒ Same fragility class as the exp recurrence's Na-doublet flip: a marginal SCF, a small perturbation, a
different trajectory.  **Not evidence that the kernel is wrong — but not nothing either.**  Defaulting on
should wait for either a robust Na2 fixture or an understanding of the trajectory sensitivity.

### Step 5 (original statement) — The NON-ORTHO kernel — **THE ONLY ONE WITH PRODUCTION USERS**.
"Mathieu's trick" (`cp2k/src/grid/cpu/grid_cpu_collint.h:532,752`): the quadratic form factors EXACTLY into
three 2-D tables, \f$e^{-\zeta_p Q(i,j,k)}=T_{ij}T_{jk}T_{ki}\f$ with
\f$T_{ij}=e^{-\zeta_p(d_i^2h_{ii}+2d_id_jh_{ij})}\f$, plus re-expansion of the polynomial in grid-index
powers.  Per grid point the Gaussian is then three lookups and two multiplies.
⚠ **Seed the table recurrences SYMMETRICALLY OUTWARD FROM THE CUBE CENTRE**, as CP2K does.  We already know
why from the 2026-08-26 underflow analysis: a recurrence seeded in the tail can start at (or below) the
underflow floor and stay zero through the vertex, where the value is significant.

### Step 6 — The integrate-back mirror.
Same tables, reversed contraction (grid → coefficients → block).  The adjointness gate from step 2 is the
acceptance test.

### ✅ Step 7 — DONE 2026-08-27: **the cache is deleted, the task list replaces it, and the kernel is the default.**

**FIRST, THE RE-MEASUREMENT THE PLAN DEMANDED — both sides, same binary, same probe, one session.**
`MNO_SKIP_FM=1 GPW_MNO_NMAX=2 GPW_REPORT=1 GPW_CONTRACT_CUBE=1 OMP_NUM_THREADS=1 GPW_OMP_THREADS=1`,
AFM arm, symmetry FREE, no fold (`[fold] … NONE`), 99% CPU on both rows — i.e. serial, unfolded, and the
two runs differ ONLY in the stream budget (`GPW_STREAM_BUDGET_PTS=0 GPW_STREAM_BUDGET_PTS_F32=0`).

| MnO AFM-II Γ probe | cache ON | cache OFF | what the 3.9 GB bought |
|---|---|---|---|
| `scf: collocate density` (×18) | 15.60 s | 48.11 s | 3.08× |
| `scf: integrate-back` (×10) | 8.75 s | 22.80 s | 2.61× |
| **the two box-walk buckets** | **24.35 s** | **70.92 s** | **2.91×** |
| `setup: collocation stream build` | 34.23 s | 14.16 s ¹ | — |
| wall / CPU | 1m35.9s / 95.3 s | 2m02.0s / 121.4 s | **1.27×** |
| **peak RSS** | **3915 MB** | **155 MB** | **0.040×** |
| `Efinal` | −59.69580308 | −59.69580301 | 7e-8 Ha apart |

¹ the counting pass the build ran even at budget 0 — pure waste that deleting the cache removes outright.
Net it out and the honest post-deletion figures are ~1m48s wall / ~107 s CPU, so **the 3.9 GB bought about
1.1× on the run**, not the ~18× it bought before the contraction kernel landed.  ⇒ The trade the cache
existed for is gone: **25× the RAM for ~1.1× the CPU.**  Deleted.

★ **AND THE 2.91× IS NOT WHAT IT LOOKS LIKE.**  It is a replay-vs-recompute ratio on the buckets *alone*,
and those buckets are now ~58% of the run — which is exactly why the whole-run number is 1.1×.  Quoting the
bucket ratio as "the cache was worth 3×" is the mistake this table exists to prevent.

**WHAT LANDED**

1. **The value cache is gone** — `PairOffsetStream`, `PairStreams`, `StreamCache`, `EnsureStreams`,
   `StreamBudgetHeadroom`, `ReleaseStreams` (and its `LatticeSum1E` face + three forwarders), the two
   `GPW_STREAM_BUDGET_PTS*` knobs, the fp64/fp32 tiering and the run-length encoding.  ~370 lines.
   With them go the whole class of defects they carried: budget starvation across grid stages
   (doc/GPWPlan.md 0.5(b)), the self-heal rebuild, the fp32 overflow tier's 6e-8 relative replay noise, and
   `GPW_Evaluator`'s dtor having to hand a global budget back.
2. **The TASK LIST replaces it** (§3c): `BoxTask{n, Roff, pf}` per (shell pair, offset), derived ONCE per
   basis instance and never invalidated — the cell and the centres are that instance's own data, and the
   list is fold-independent (the T3 fold is a per-iteration filter, not a build-time one).  On MnO it is
   35 k tasks ≈ 1.3 MB against 3.9 GB.  Reported as `grids.boxTasks`; built in `setup: box task list`
   (measured **0.6 ms** on Si Γ).
3. **`GPW_CONTRACT_CUBE` DEFAULTS ON** — A1, and it is not separable from A7: with the walk as the default
   the cache-less path is ~18× slower, so the two had to land together.  `=0` opts back onto
   `ForShellPairBox`, which stays as the reference implementation the unit oracle checks against.
4. **The collocation and integrate-back paths are now ONE route each**, not two.  The per-component-pair
   replay lambdas (`scatter`, `integratePair`) are deleted; every term goes through the shell-blocked
   `scatterShell` / `integrateShell`, which is where the contraction kernel already lived.  The threaded
   loops lost a whole pass each and now run LONGEST-FIRST off the task list (§3c reason 2); the serial
   loops keep the natural order, because reordering changes the cross-pair grid reduction and serial is
   the bit-anchor path.
5. **Unit coverage first, per §1**: `M_PG_BoxWalk.TaskListIsExactlyTheEnumerationItReplaces` asserts the
   hoisted list is bitwise the same offsets, in the same order, with the same `pf`, as the un-hoisted
   enumeration — and that a per-iteration weight can only ever remove tasks from it.

**★★★ AND IT FOUND TWO LATENT DEFECTS THE CACHE HAD BEEN HIDING.**  Both were live in code that only
OVER-BUDGET pairs used to reach; routing every pair through it fired them at once.  This is the strongest
argument for the deletion that was not in the plan when it was written: **a cache that makes one of two
routes rare does not remove the second route, it removes the test coverage of it.**

1. ⛔ **THE SHIFTED-MP DEFECT, ALIVE IN A SECOND COPY.**  `integrateShell`'s D-aware screen still read
   \f$\mathrm{Re}[D_{ij}\overline{e^{ikR_n}}]\f$ where the per-pair path had long since been corrected to the
   magnitude \f$|D_{ij}|\f$ (doc/Benchmark.md footnote 1).  A real part is not a magnitude: at a
   quarter-integer \f$k\f$, \f$e^{2\pi ikn}=i^n\f$ is purely imaginary for every ODD offset, so against the
   real, positive `screenD` the test discarded every odd-offset term and \f$H\f$ came out exactly real.
   Symptom the moment the cache went: `GPW_SCF.SR_2x2x2ShiftedMP_vs_CP2K` converged to −3.735 instead of
   −7.868 — **4.1 Ha**, NON-AUFBAU.  Fixed; the row now reproduces its banked −7.868473428 to all 10 s.f.
2. ⛔ **THE WALK'S PER-COMPONENT \f$|v|<\varepsilon_{ij}\f$ SCREEN BROKE THE ADJOINT.**  With the cache, both
   directions replayed ONE frozen stream, so their truncations were identical by construction.  Without it
   the collocation walked a D-aware box while the gather did not, and `GPW.AnalyticIntegrateBackAdjoint`
   plus both XC finite-difference gates sat ~1.2e-8 relative off — scaling linearly with `GPW_DENSITY_EPS`
   (green again at 1e-12).  The screen is DELETED from both walk fallbacks: it never shrank the walk (the
   box is already sized by \f$\varepsilon_{union}\f$, so it only declined to accumulate a value it had
   already computed — step 3), and removing it makes the walk and the contracted cube truncate IDENTICALLY.
   ⇒ The reference implementation is now a truer reference, and the adjoint holds on either route.

**⚠ WHAT IS STILL UNMEASURED**: §3c reason 1, batching by \f$(L_a,L_b)\f$ so the innermost `for a<=lp` runs
with \f$l_p\f$ known at compile time (CP2K's `always_inline` low-\f$l_p\f$ variants).  The task list is the
prerequisite and it now exists; the dispatch does not.  **That is the next lever on this path** — and
`doc/Benchmark.md`'s re-taken MnO row says why it is needed: on a long imposed run the two box-walk buckets
are 477 s of 976 s CPU, so that row LOST 1.47× when the cache went even as its RAM fell 2.9×.

### Step 7 (original statement) — Reconcile the three consumers of the current geometry.
- **The stream cache** — this is where the RAM win lands.  It costs **4218 MB** because it stores per-point
  VALUES, \f$O(n^3)\f$ per (pair, offset); separable tables are \f$O(n)\f$ / \f$O(n^2)\f$.  It may become
  unnecessary outright.  ★ First route that can deliver CP2K's memory profile AND its CPU — which is what
  `doc/Benchmark.md` said was needed and what caching never gave.
- **The T3 orbit fold** — keyed per (pair, offset), reads the orbit-projected \f$D\f$; needs re-derivation
  against a per-shell-pair cube.
- **The static local-PP sweep** — same walk, different tolerance rule (absolute \f$\kappa\f$, explicit
  `pairLevels`).  It is ~2% of a run but it shares the code, so it must not be forgotten.

### ✅ Step 8 — ACCEPTANCE, PASSED 2026-08-27.

**`ctest -j8`: 792/792 on BOTH kernel settings** — default (contracted) and `GPW_CONTRACT_CUBE=0` (the
reference walk).  ⚠ The polarity of the rot warning has REVERSED: the WALK is now the arm nothing exercises
by default, so `GPW_CONTRACT_CUBE=0 ctest -j8` is the run to keep making at breakpoints.  Test count is
unchanged (−1 `GPW.StreamCacheReleaseUnstarvesLaterGrid`, whose entire subject was deleted;
+1 `M_PG_BoxWalk.TaskListIsExactlyTheEnumerationItReplaces`).

**Four gates moved, and each moved for a reason that is recorded rather than tuned away:**

| gate | what happened | resolution |
|---|---|---|
| `GPW.StreamCacheReleaseUnstarvesLaterGrid` | its subject (budget, residency, starvation, self-heal) no longer exists | DELETED, with a tombstone comment pointing at the task-list unit test |
| `GPW_SCF.SR_2x2x2ShiftedMP_vs_CP2K` | 4.1 Ha off, NON-AUFBAU | a real BUG, fixed — the `Re[]`-vs-magnitude screen above |
| `GPW.StreamFoldReducedMatchesFull_SiDiamond_HalfK` | reduced-vs-full 1.18e-10 against a 1e-12 tolerance | the tolerance was asserting the CACHE's bit-identity.  At Γ the two arms still agree to machine precision (same union tolerance); at general k the offset phase varies across an orbit so they agree at the SCREENING tier.  Tiers are now k-dependent and the unweighted h gate keeps the tight one as the control |
| `GPW_SCF.ImposedOrderLostIsAPostconditionFailure_Na2Box` | did not converge inside its cap | NO stable mixing step exists (α 0.65 ✔ / 0.7 ✔ / 0.75 ✘ / 0.8 ✔) — the fixture is chaotic in Δρ while the physics it gates is not.  The gate now SWEEPS four measured-good steps.  ⇒ The durable fix is the Δρ/N convergence gate (A4): the run reaches a state its own fingerprint calls *"DENSITY-DEGENERATE … benign"* while `DidConverge()` calls it a failure |

**`doc/Benchmark.md` is re-taken** — Si Γ, Si 2×2×2 (both meshes), NaF SR2 Γ, NaF SR2 2×2×2 and NaF full-SR
Γ, plus the MnO AFM-II VA row, each through `scripts/bench` with the run banner beside it.  Headline: **two
rows now beat CP2K on CPU outright** (Si Γ 0.48×, NaF full-SR 0.43×), **every** qchem row is under CP2K's
RAM for the first time, NaF full-SR went 3090 MB / 219 s → **58 MB / 43 s**, and the MnO row went
1323 → 463 MB while LOSING 1.47× CPU.  Read that last one: the cache was still earning its keep on the
longest row, and the case for deleting it rests on RAM plus the two defects, not on a free lunch.

## 4. Expected payoff, stated as a range with its uncertainty

Per-point work goes from *2 exps + 6 power tables + \f$(n_I{+}n_J)\f$ products + a loop over live component
pairs* to *\f$l_p{+}1\f$ fused multiply-adds against a table*.  That SUGGESTS 5–15× on the box walk and
4–8× overall.

⛔ **This is a projection from a cost law, and the last projection made from a profile share in this
campaign was off by 7×.**  It is not a commitment, and step 0 exists precisely to replace it with a
measurement before anyone spends a week.

## 4b. ⏸ PARKED — folding the \f$\Omega_{ij}\f$ list by symmetry (user, 2026-08-27)

> *"Is it possible to symmetry fold the \f$\Omega_{ij}\f$ list?"*

**YES, and the structure is unusually favourable — but MEASUREMENT says it is not where the time is.**

**THE ALGEBRA.**  With \f$p=a+b\f$, \f$P-A=-\tfrac bp\overline{AB}\f$, \f$P-B=\tfrac ap\overline{AB}\f$ and
\f$E_{ij}=e^{-\frac{ab}p|\overline{AB}|^2}\f$ — so \f$E_{ij}\f$ and every `H2` coefficient depend on the
centres ONLY through the difference vector \f$\overline{AB}=A-B\f$.  Only \f$P\f$ carries absolute
position.  ⇒ **\f$\Omega\f$'s content is TRANSLATION-INVARIANT**, a function of
\f$(a,b,L_A,L_B,\overline{AB})\f$, and the expensive part (the E-coefficient recursion) is shared by every
pair with the same separation vector.

**AND UNDER A CUBIC POINT OP THE FOLD IS AN INDEX REMAP, NOT A RECOMPUTATION.**  The 48 ops of \a m-3m are
SIGNED AXIS PERMUTATIONS in Cartesian coordinates, and `Hermite2` already stores `d`,`e`,`f` PER AXIS
indexed \f$(N,n_a,n_b)\f$ — so \f$\Omega(W\overline{AB})\f$ is \f$\Omega(\overline{AB})\f$ with the axis
tables permuted and odd powers sign-flipped.  (A general op needs the Cartesian representation matrices,
which is the machinery the T3 fold already carries.)

⛔ **WHY IT IS PARKED.**  Measured 2026-08-27 with the contraction kernel on, the \f$\Omega\f$-consuming
buckets (analytic 1E lattice sums + local-PP SHORT + vet basis) are **6% of a Si run and 0.6% of an MnO
run**.  Collocation and integrate-back do not touch \f$\Omega\f$ at all — the walk reads `radials`/`pols`
directly and `MakePairPoly` derives \f$\alpha_p,P,E_{ij}\f$ itself.  ⇒ Folding it PERFECTLY would buy less
than the exp recurrence did, and that bought 3%.

★ **WHERE IT WOULD PAY, and the second one is the interesting case:**
1. **The molecular / 4-centre path**, where \f$\Omega\f$ IS the central object (`Repulsion4C` takes two by
   reference) — see `doc/ERI4Rework.md`, which is already about banking pairs.
2. **The lattice path PRECISELY BECAUSE IT CANNOT CACHE.**  `GeometryCacheBudget(0)` deliberately evicts
   every \f$\Omega\f$ immediately on lattice paths (the 2026-07-22 OOM: ~2M Hermite3 tables, 11 GB), so it
   rebuilds per (pair, offset) BY DESIGN.  A symmetry fold gives **reuse WITHOUT storage** — the same
   "task list, not value cache" insight as §3c, applied to \f$\Omega\f$.  That is the version of this idea
   worth revisiting, and the moment to revisit it is AFTER step 7, when the setup buckets are a larger
   share of what is left.

## 5. Open questions to answer before step 4

1. **The screen** (step 3) — the only one that can force a redesign.
2. **Contracted bases.**  Every current basis is uncontracted, so \f$\Omega\f$'s primitive pair IS the
   shell pair.  A contracted basis needs a sum over primitive pairs per shell pair, each with its own
   \f$\alpha_p\f$, \f$P\f$ and box.  Decide whether the cube is per PRIMITIVE pair (CP2K's choice) or per
   shell pair with an inner primitive loop.
3. **Level assignment.**  `PairLevel` is a shell-pair property today; per-primitive-pair cubes would make
   it a primitive-pair property, which is what the REL_CUTOFF rule actually wants (the doc there already
   says so: *"For an uncontracted basis this IS the per-primitive-product assignment"*).
