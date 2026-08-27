# The collocation rewrite — from a per-point walk to a separable contraction

**Goal: close the remaining ~22× against CP2K on the on-the-fly path**, and take the 4218 MB pair-stream
cache down with the same change.  Forward plan; the measured record that motivates it is in
`doc/OpenWork.md` ("THE ON-THE-FLY BOX WALK") and `doc/Benchmark.md`.

## 0. Vocabulary and the measurement rule (both learned the hard way, 2026-08-26/27)

- **"BOX WALK", never "kernel"** (user).  *Kernel* is overloaded — M&D integral kernels, XC kernels, GPU
  kernels.  The box walk is `NR_Evaluator::ForShellPairBox` plus the two lambdas that consume it,
  `scatterShell` (in `CollocateDensity`) and `integrateShell` (in `IntegratePotential`): the code that,
  for one shell pair at one cross-cell offset, walks the grid points of the product's compact box.
- **"DENSITY-MATRIX SUB-BLOCK", never "D-block"** (user).  In a codebase whose hardest case is Mn *d
  shells*, "D-block" reads as angular momentum.  \f$D_{ij}\f$ over \f$i\in\f$ shell I, \f$j\in\f$ shell J.
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
⚠ **ORTHORHOMBIC only, i.e. the OPTIMISTIC case.**  This gate could only ever KILL the plan; passing it
does not prove the triclinic path pays.  Step 5 needs its own number.

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

### Step 2 — Choose the expansion basis.  **Unit level.**
Hermite (reuse `H2`'s `d`/`e`/`f`) or binomial re-expansion to monomials about \f$P\f$ (CP2K's route).
**The criterion is NOT arithmetic count.**  It is which one makes integrate-back come out as the EXACT
ADJOINT with the same coefficients: \f$\int\rho V = \mathrm{Tr}(Dh)\f$ to machine precision is load-bearing
(it is what makes the GPW energy variational), and it is far easier to preserve than to repair.

### Step 3 — Settle the SCREEN.  This is the real blocker, not the density matrix.  **Unit level.**
Today's D-aware tolerance is PER COMPONENT PAIR (`epsHere[k]`, applied to \f$|val|\f$ at each point), and
pre-filtering on it is what made the shell hoist pay at all — measured 1.03–1.08× without it against
2.13× with.  A contracted cube has no per-component identity: one cube serves the whole shell pair.
- **Try the cheap option first:** keep the per-component pre-filter deciding only WHICH PAIRS ENTER the
  contraction.  It already runs once per offset with no grid points involved, so it may carry over
  unchanged.
- **Only if that fails:** derive a bound on \f$|c_{NLM}|\f$ and screen the coefficient tensor.

### Step 4 — The ORTHORHOMBIC kernel, gated.
Three 1-D tables, three-step contraction.  Simpler, and it gets the contraction machinery right without
cross terms.  Validated against the step-1 oracle, then the cubic-box integration cases (Na, O2, Mn2, Na2).

### Step 5 — The TRICLINIC kernel — what Si/FCC, NaF and MnO actually run.
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

### Step 7 — Reconcile the three consumers of the current geometry.
- **The stream cache** — this is where the RAM win lands.  It costs **4218 MB** because it stores per-point
  VALUES, \f$O(n^3)\f$ per (pair, offset); separable tables are \f$O(n)\f$ / \f$O(n^2)\f$.  It may become
  unnecessary outright.  ★ First route that can deliver CP2K's memory profile AND its CPU — which is what
  `doc/Benchmark.md` said was needed and what caching never gave.
- **The T3 orbit fold** — keyed per (pair, offset), reads the orbit-projected \f$D\f$; needs re-derivation
  against a per-shell-pair cube.
- **The static local-PP sweep** — same walk, different tolerance rule (absolute \f$\kappa\f$, explicit
  `pairLevels`).  It is ~2% of a run but it shares the code, so it must not be forgotten.

### Step 8 — Acceptance (and ONLY here does the integration level appear).
Si → NaF → MnO uncached probe (`MNO_SKIP_FM=1 GPW_MNO_NMAX=2`) → full `ctest -j8` → a `doc/Benchmark.md`
row with the run banner copied beside it.

## 4. Expected payoff, stated as a range with its uncertainty

Per-point work goes from *2 exps + 6 power tables + \f$(n_I{+}n_J)\f$ products + a loop over live component
pairs* to *\f$l_p{+}1\f$ fused multiply-adds against a table*.  That SUGGESTS 5–15× on the box walk and
4–8× overall.

⛔ **This is a projection from a cost law, and the last projection made from a profile share in this
campaign was off by 7×.**  It is not a commitment, and step 0 exists precisely to replace it with a
measurement before anyone spends a week.

## 5. Open questions to answer before step 4

1. **The screen** (step 3) — the only one that can force a redesign.
2. **Contracted bases.**  Every current basis is uncontracted, so \f$\Omega\f$'s primitive pair IS the
   shell pair.  A contracted basis needs a sum over primitive pairs per shell pair, each with its own
   \f$\alpha_p\f$, \f$P\f$ and box.  Decide whether the cube is per PRIMITIVE pair (CP2K's choice) or per
   shell pair with an inner primitive loop.
3. **Level assignment.**  `PairLevel` is a shell-pair property today; per-primitive-pair cubes would make
   it a primitive-pair property, which is what the REL_CUTOFF rule actually wants (the doc there already
   says so: *"For an uncontracted basis this IS the per-primitive-product assignment"*).
