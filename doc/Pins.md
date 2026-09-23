# Durable pins — the invariants a session must not violate

**Cut 2026-09-08 out of `doc/OldPlans/GPWPlan.md`'s "Durable pins / invariants" section**, because that file is a
RECORD of the 2026-07 campaign and the pins are not: they govern work all over the tree, and burying
project-wide invariants inside a finished campaign's plan is how they get missed.  Two user rulings that
lived only in session memory are folded in and marked as such.

These are **rulings, not preferences.**  Each is here because violating it produced a wrong number, a wrong
interface, or a retracted verdict at least once.  If you think one is wrong, say so and get it changed —
do not work around it quietly.  **Three files, three questions** (ruled 2026-09-16): `CLAUDE.md` answers
*how do I work here* (conventions, build/test/box discipline, tool paths, the doc system); this file answers
*what must the code obey* (physics, numerics AND design invariants — one paragraph each, earned by a wrong
number, a wrong interface or a retracted verdict, pointing at the RECORD that earned it); a RECORD in
`doc/` answers *why is it this way* (the evidence and the rejected alternatives).  A pin is the distillate;
the argument stays in the record.

---

## 1. THERE IS NO CUT IN r SPACE — Gibbs ringing is like a wrecking ball

> **"THERE IS NO CUT in r space, Gibbs ringing is like a wrecking ball."** (user, 2026-07-16; wording
> sharpened 2026-09-08 — the *reason* is now part of the pin, because "no cut" alone reads as fastidiousness
> and it is not.)

A real-space lattice sum is an **ε-CONVERGED SERIES for a FIXED operator**.  Magnitude screening is the ONLY
truncation mechanism.  **A radius must never appear as a parameter, member, or concept in any interface** —
not user-facing, not internal.

- A truncation radius yields a **DIFFERENT operator**, not "the operator to ε".  Sharp truncation in r is
  multiplication by a step function, whose transform rings — and that ringing lands on everything downstream.
  Measured: the Rcut=2a NaF metric lost **2.25 e** per mid-slosh loading.
- A radius must **never be a conditioning crutch**.  That job belongs to the basis, or to rank reduction.
- ⚠ **The G direction is different IN KIND.**  The Ecut ball is a PROJECTION onto a finite auxiliary
  subspace: variational (adjoint-exact), exponentially controlled, systematically improvable.  That is a
  legitimate resolution dial, not a cut.
- **End state: ONE knob per direction** — ε in r (convergence tolerance), Ecut in G (projection resolution).

## 2. Everything is a fit — and a fit is (integration GRID) × (FIT BASIS)

> (user, 2026-09-06.)  **A quadrature IS a fit.**  Never say a route has "no fit".

The grid and the basis are **ORTHOGONAL AXES**.  A grid does not imply a basis.  Which pairings are worth
using is high-level **POLICY** and must never be hard-coded — "Becke = delta basis" is itself a conflation.
The uniform route is {uniform raster} × {plane-wave \f$\{G\}\f$}; the Becke route is {Becke mesh} × {delta
basis}.  Both are fits of the same \f$v_{xc}\f$, which is what makes one scoreboard legitimate.

## 3. Fit quality is grid-convergence of ρ, NEVER \f$\Delta E_{total}\f$

The fit is non-variational, so a total-energy difference does not bound its error and can flatter it.
Score a fit by the convergence of ρ (or of the operator that is actually diagonalised —
\f$\max|\Delta V_{xc}(i,j)|\f$), against a fine reference in the same family.

## 4. Report an INTEGRATED observable, never a point sample of a field

> (user, 2026-08-19.)

MnO's \f$m_{stag}\f$ sampled at one point 0.7 bohr along +x is a spin **DENSITY**, not a moment; it was never
derived, and its direction-dependence IS the d-occupation confounder.  Valid as a **collapse detector only**.
⚠ A frozen-density point probe also understates the self-consistent shift on a metal — grid error feeds back
through the density, the Fermi level and the occupations.

## 5. Spin-polarized is the native formulation

Unpolarized is the \f$\zeta=0\f$ collapse, not the base case.  New terms are written spin-native
(`FittedVxcPol` / `FittedVcorrPol`), and this governs everything still to come — GGA, +U, non-collinear.

## 6. PW fitting must look IDENTICAL to molecular fitting at interface level

> (user.)  No "Fourier" in an abstract face.  Any fit/aux basis comes from the orbital basis via
> `Create{CD,Vxc}FitBasisSet(...)` — **the factory is the seam even when it is trivial**, and inlining the
> trivial G-space fit was REJECTED for exactly this reason.  **Never assume `orbital == fit`.**

## 7. Ill-conditioning is a BASIS problem

Not a solver bug and not a code bug.  "LASolver" symptoms are basis conditioning (SIPP diffuse → SIPP_SR).
Use well-conditioned bases for SCF; the accuracy levels are Low/Medium/High, and the old N3/N5 pools were
test-only and no longer exist.

## 8. Two self-consistent lattice schemes — do NOT mix them

**(A)** complete-Bloch analytic single-sum matrices (what GPW has; correct as \f$R_{cut}\to\infty\f$).
**(B)** truncated-Bloch collocation Gram matrices (always PSD).
Scheme-B overlap with scheme-A analytic kinetic gave \f$E_{kin}=-300\f$.  Stay in scheme A.

## 9. Pseudopotential smoothness is what makes GPW work

All-electron cores are too sharp; validate with a well-conditioned GTH valence basis, never all-electron.
GAPW is out of scope.  Relatedly, the **pseudo-wall is an asymptote**: change a basis interface only for a
NEW INTEGRAL TYPE.

## 10. Regression style for periodic energies

Periodic/GPW energies are **"did-E-move" anchors** — pin the converged value, no `Converged()` guard.  Where
a real-space-on-lattice quantity must equal its finite counterpart, assert **bit-consistency** (`L_PP`-style)
rather than an absolute oracle.

**A moved anchor is RE-JUDGED against an INDEPENDENT route, never merely refreshed** (added 2026-09-16 from
`doc/Records/TestSuitePlan.md` §2): KP-0 re-pinned \f$-7.45137\to-7.45294\f$ only after the band-folding-equivalent Γ
supercell agreed.  And the two kinds of failure are not alike: an ENERGY anchor can go stale and must be
judged; a failing CHARGE, count or weight sum is physics and cannot (user, 2026-09-09).

## 11. An explicit phase beats an automatic one — the `UseChargeDensity` lesson

> **USER, 2026-09-08:** *"this code used to have exactly `tDynamic_HT::UseChargeDensity(cd)`, and I was too
> clever by half trying to make it all 'automatic' which ended being a source of bugs and finally blocking
> irrep omp … lesson learned!"*

The term stack once had an EXPLICIT *"here is the density, prepare yourself"* call.  It was replaced by
automatic, self-correcting, per-object density-serial guards — each one locally correct, and collectively:

- a class of bugs (a memo believing it was fresh; the `itsRho`/`itsXCMix` aliasing that made
  \f$\alpha=0.25\f$ and \f$\alpha=1.0\f$ produce **bit-identical** runs; the DM-source staleness guard
  that had to be added back as a LIVE check because its `assert` was compiled out in Release), and
- an architectural block: lazy-fill-on-first-touch turned every k-independent memo into a
  write-on-first-touch, which is what stopped the per-block loop being threadable at all.

`tHamiltonian::RefreshForDensity` (2026-09-08) is `UseChargeDensity` returning by another name.  ⇒ **When
work must happen once per density / per iteration / per geometry, give it a PHASE and a name.**  An
automatic guard hides the phase structure in call order, and call order is not a thing anyone can see.

⚠ **The same suspicion now falls on the \f$H_{ij}\f$ cache** (user, same message): `tDynamic_HT_Imp` stores
its result in `mutable CacheMap itsCache` keyed by `Irrep`, filled during the block loop, purely so the
ENERGY pass (`GetEMatrix` → `IrrepCD::DM_Contract`) does not recompute what the Fock pass just built.  Same
shape, same smell — and it is the one remaining write inside the loop.  ⛔ **`DB_Cache` is NOT the answer, and the
reason is LIFETIME rather than the key**: it is a process-wide store that **never evicts**, built for
cross-run sharing (its own header: *"allow data sharing between separate runs"*), while the
\f$H_{ij}\f$ memo turns over **every SCF iteration** and is never reusable across runs.  Twenty iterations
would leave twenty generations of every block in it.  *Ask what a cache EVICTS before asking what it keys
on.*  ▶ The fix in the spirit of this pin is an explicit **per-iteration scope** that owns the matrices and
dies with the iteration — which also retires the last write-shaped obstacle in the block loop, because the
slots are CREATED in the phase and the loop only fills nodes that already exist.  Filed as `doc/CleanupCandidates.md` R1.0h.

## 12. No grad-student knobs

Policy enums, not numeric dials.  A number a user has to tune is a design failure looking for somewhere to
live.

## 13. A symmetry op acts on GRID INDICES as \f$DUD^{-1}\f$, never as \f$U\f$

A grid point is \f$k=(i+s)/N\f$ **componentwise**, so an op \f$U\f$ induces on the index lattice the
CONJUGATED map \f$M = D U D^{-1}\f$ with \f$D=\mathrm{diag}(N)\f$, i.e. \f$M_{ab}=N_a U_{ab}/N_b\f$ and
\f$i' = M(i+s)-s\f$.  \f$M=U\f$ **only when the mesh is isotropic**, which is why applying \f$U\f$ to the
indices is right on every \f$n\times n\times n\f$ mesh and silently wrong on \f$2\times1\times1\f$.

⛔ **AND A MOD-\f$N\f$ WRAP IS NOT A LICENCE TO PROCEED.**  Reducing a stray image back into \f$[0,N)\f$
always yields *a* grid point, so a map that is not a mesh symmetry looks like one.  The test is
\f$M\f$ INTEGRAL (then \f$\det M=\det U=\pm1\f$ makes it unimodular over \f$\mathbb{Z}\f$, so the action is
a PERMUTATION) plus \f$(M-I)s\f$ integral for a shifted mesh — both properties of the OP, checked once,
not of the point.  An op failing either is not a symmetry of that mesh and must be dropped whole.

**What it cost:** the IBZ stars overlapped on Si \f$2\times1\times1\f$, \f$\Sigma w=1.5\f$, and the SCF
carried 12 electrons in an 8-electron cell (KP-0, 2026-09-09; record in `doc/Records/OpenWork_History3.md`).
▶ Corollary, still open: the group that symmetrizes \f$\rho\f$ must then be intersected with the mesh
symmetries, or the density is projected into a symmetry the sampling does not have
(`doc/CleanupCandidates.md` R1.0r).

## 14. An IrrepBasisSet carries ONE irrep of G; libraries follow the FAMILY, modules carry the GROUP

`BasisSet = ⊕_irreps IrrepBasisSet`, the irrep a label of G, the symmetry group of the Hamiltonian.
Carriers are not invariants (\f$Y_{lm}\f$, \f$e^{i\mathbf{k}\cdot\mathbf{r}}\f$, a SALC each *transform as*
an irrep); k is an irrep label of the translation group exactly as \f$l\f$ labels O(3).  Three orthogonal
axes: **G** (block labels, `qcSymmetry`), **family** (the analytic seed — Gaussian, Slater, BSpline, PW,
delta; = the integral ENGINE), **construction** (subduce \f$G_{big}\downarrow G\f$ or induce site \f$\uparrow G\f$;
derived, never free).  Placement: **a LIBRARY is an engine** (its mass is its integrals, which factorise
over the family, so `UnitCell` inside the Gaussian engine is legitimate and a cut on the G axis is ruled
out); **a MODULE carries the G** (`…Gaussian.Point.*` vs `…Gaussian.Lattice.*`), enforced by the ctest grep
*no `.Point.` module imports `qchem.UnitCell` or `qchem.Symmetry.Lattice_3D.*`* (`scripts/audit-basisset-gtags`).
Spin is a factor of G (\f$G_{spatial}\times SU(2)\f$) until a double group dissolves it — which is why
Pol/UnPol is an imposed SUBGROUP (V1.37), not a type.  **What it cost:** the tree had been cut on the wrong
axis (`qchem.UnitCell` imported inside `Molecule/`); V1.33 re-cut it in eleven commits.  Record:
`doc/Records/BasisSetTaxonomyPlan.md` §1; Doxygen `\ref basisset_taxonomy`.

## 15. Smearing needs kT ABOVE the frontier splitting, and GDM as built is fixed-occupation

Two measured facts from the Fermi-smearing build (`doc/Records/GPWPlan1.md`, 2026-07-26): **(a)** kT must EXCEED the
frontier splitting or the occupations slosh-rotate instead of converging (NaF: 1e-2 converges, 1e-3 does
not); **(b)** the GDM direct minimiser DIVERGES under smearing because its geodesic direction is the
fixed-occupation \f$[F,D]\f$, not the free-energy gradient (which carries an occupation-response term).
⇒ smeared runs use the fixed-point stage (DIIS/Kerker/Pulay); GDM tail-polishes only at kT=0 until it has a
smearing-aware direction.  ⚠ This is a limit of OUR parameterisation, not of the method — CP2K's OT has the
same axis scaffolded; **never conflate hold-the-block with don't-smear** (`doc/Records/SCFStrategyPlan.md`).

## 16. A basis SPAN can reverse a magnetic ordering — match spans before comparing to an oracle

MnO with Cartesian d (its \f$r^2e^{-\alpha r^2}\f$ s-contaminants): FM below AFM by 40 mHa.  The same
cell through the spherical-d view: AFM below FM by 45.5 mHa — **no code bug, the span**.  Mechanism: the
extra s-like functions let the density rearrange s-character away from the l=0 KB projectors (an l=0
repulsion DODGE worth 37% of the weak basin's reward), a freedom the oracle's spherical basis structurally
lacks.  And the earlier "8 mHa agreement" with CP2K was contaminants-vs-diffuse COMPENSATION.  ⇒ Before
any energy is compared to an oracle, the two spans are matched exponent-for-exponent
(`valence_lowq_sph` v2 = the CP2K transcription), and a Cartesian-d basis is never used for a d-metal
ordering question.  Record: `doc/Records/SphericalLatticePlan.md` I0–I2.

## 17. A class reports CONTEMPORANEOUSLY with its own activity — console order == execution order

`CurrentReport` is a global sink; each class emits at the moment it does the thing, and never tells
another class to emit.  An `Emit*()` method on an abstract face, or a "reporter" that PULLS state out of
objects after the fact, is the defect (user, 2026-09-11; V1.5 deleted the `Emit*()` faces).  Corollary for
trace columns: a printed number is either physics or a gate the run CONSUMES — printing \f$\alpha_{eff}\f$
implied it was used, and it was deleted for that reason (user, 2026-09-13).  Record: `doc/Records/RunReportPlan.md`
(the design), `doc/CleanupCandidates.md` V1.5.

## 18. XC is fed the MIXER'S density; a separately-damped XC feed destroys Kerker's mode selectivity

Four collapsed MnO states (−45.5, −46.3, −56.4, −38.5 against the converged −61.403) came from ONE cause,
measured as a monotone dose-response in \f$E_{ee}\f$ (13.5 → 29.0 → 35.1): handing \f$V_{xc}\f$ a density
damped by a FLAT \f$\alpha_{eff}\f$ un-damps the low-G CHARGE mode 2.4× while barely touching the AFM mode
— **Kerker's low-G charge-slosh damping is what holds the AFM basin**, and the moment death is a
consequence of the charge runaway, not a spin effect.  So the XC feed is never a second, independently
damped copy of the density.  The ρ≥0 GOAL survives (the DM route gives 0 negative points against 15%
for the band-limited ρ̃): the form that keeps it is the **cusp deficit**,
\f$\rho_{XC}=\rho_{mix}+(\rho[D]_{exact}-\rho[D]_{BL})\f$ — Hartree's own mixed array plus the sharp
content only the DM can supply, with no \f$\alpha_{eff}\f$ to choose (N4, still to be measured).
Two facts that ride with it: **\f$\tilde\rho_{mix}\f$ and \f$\rho[D]\f$ have DIFFERENT fixed points**
(ρ̃ is a band-limited fit projection; NaF 139 μHa, MnO terms ~100 mHa at 8 μHa total), so "at convergence
they agree" is false; and **\f$E_{ee}\f$ is a validated charge-slosh detector** (T3).  Record:
`doc/Records/OpenWork_History4.md` "ITEM 1 MEASURED" + "N4".

## 19. Never density-screen the GATHER — \f$h_{ij}\f$ is diagonalised, not traced

Dropping a term because \f$D_{ij}=0\f$ is sound for the ENERGY (\f$\mathrm{Tr}(Dh)\f$ is blind to it) and
WRONG for the Fock matrix, which is diagonalised to make the next density: zeroing \f$h_{ij}\f$ wherever
the density vanishes is a SELF-FULFILLING truncation — a pair with no density can never acquire any, and
with a DIAGONAL SAD seed that is every off-diagonal element.  CP2K does not D-screen at all (checked in
`task_list_methods.F`: one global `eps_rho_rspace`, geometry only).  The fix that "worked" (floor a
vanishing weight) BROKE the stream fold's orbit invariance (0.14 against 2.4e-8) — so the gather's D-screen
was REMOVED (2026-09-04) and the D-aware tolerance survives only on the COLLOCATION, where the weight
really is the scatter weight.  ⚠ `cij==0` had meant BOTH "structurally absent" and "zero density"; only
the first may be excluded.  Record: History4 "THE k-SCALING GAP" + "ATTEMPT 2".

## 20. Ask what a matrix MEANS before symmetrizing it

The first cut of the T3 stream fold orbit-averaged the integrate-back's `screenD` — a matrix of
\f$|D_{ij}|\f$ MAGNITUDES, not a density matrix.  Signed averaging cancelled mixed-σ orbits to ~0, the
D-aware screen dropped live terms, and the imposed O₂ triplet collapsed by 2.3 Ha.  A screen is reduced by
the orbit **MAX** (`FoldScreenMax`), a density by the orbit PROJECTION (`FoldProjectedD` — reading the
representative's own \f$D_{ij}\f$ SAMPLES the orbit, and sampling equals projecting only if D is already
symmetric).  Gate: the dimer-in-a-box cell — single-atom Si cells miss this bug entirely.  Record:
History4 "Step 2 — ARM THE SYMMETRY FOLDS".

## 21. Pivoted Cholesky where D is PSD — and a Cholesky that FAILS is the canary, never a silent fallback

Standing preference (user, 2026-08-20): factor a density matrix by pivoted Cholesky (greedy on the
diagonal, truncation bounded by the trailing diagonal, no rotational noise), not a trimmed
eigendecomposition.  Its LIMIT, narrowed the same day: the objection to trimmed eigen comes from ORBITAL
work where the factor is INVERTED (\f$S^{-1/2}\f$, \f$1/\lambda\f$ amplification); where the factor is only
MULTIPLIED (\f$\rho_g=\|L^\dagger\Phi_g\|^2\f$) the error is bounded and eigen is admissible.  **D is PSD in
this tree only because no mixer EXTRAPOLATES D** (`LinearMixer` is convex, α∈[0,1]; Pulay/Broyden act on
\f$\tilde\rho\f$) — a property of today's mixer set, not a theorem.  It dies with a density-space Pulay,
α>1, or MP/cold smearing (negative occupations).  ⇒ a failing pivoted Cholesky is exactly the signal that D
left the cone: make it LOUD and route to the eigen split.  ⚠ PSD tests need a RELATIVE floor — an
eigensolver always returns O(ε·λmax) negatives.  (USPP/PAW augmentation charges can drive ρ<0 with a PSD
D; not this tree, which is norm-conserving.)  Record: History4 "THE ρ GEMM — LOW-RANK D".

## 22. Basis AUTO-trimming is a VET-stage, SYMMETRY-EQUIVARIANT, WHOLE-ORBIT decision on S — never a per-function filter at ortho time

Three user rulings (2026-08-14/15/26): **(a) not display-only** — the trim happens BEFORE anything
downstream is built (grid ladder, collocation task lists, KB projections all fall out of the surviving
function list; filtering at ortho time does the dropped functions' work for nothing); **(b) the rank
decision is a property of S, i.e. of the BASIS, made ONCE** — not re-derived per spin channel;
**(c) drop whole ORBITS under the (magnetic) space group, never individual AOs** — greedy per-function
pivoting resolves symmetry-tied pivots by numerical noise (runs 58–60 dropped O₁'s p(0.18) but O₂'s
s(0.15)), and a partial orbit is a symmetry-BROKEN basis that costs both site equivalence and the run's
ability to converge at all (*"would sometimes remove only 3/4"*).  Report the decision as a BASIS (species / shell /
exponent), not bare indices.  Open work: the vet-stage trim itself (`doc/OpenWork.md`).  Record: History4
"Continuous — CLEANUP".

**Addendum 2026-09-23 (user, on the word AUTO).**  The vet-stage trim is an **automatic** trim and stays
one — what is refuted is auto-trimming in the other two places it was ever tried: the SCRIPT that rewrote
the committed `.bsd` (retired; VA/VB exist because a run could not say which span it used) and the
ORTHO-TIME per-function drop.  ⇒ the ortho-time path gets exactly TWO behaviours: **(1) shut up and work,
or (2) make noise with enough information to fix the BASIS** — never a silent third option that quietly
edits the span.  ⚠ Today it does the third: `LASolverCholeskyPivoted` drops and prints `dropped AO index 47`,
a bare index, from a linear-algebra layer that does not know what a shell or an exponent is.
★ **The cure is an EXCEPTION, not a decorated return type** (user: *"this may be one situation where
exceptions are a good design.  Probably cleaner than decorating the LASolver return types with fallible
flags and other index info"*) — `throw` is already this tree's marker (CLAUDE.md), the throw carries the
indices and pivots, and the layer that OWNS the basis catches it and names species/shell/exponent.
**Two things make the vet-stage trim harder than it looks, and the unit tests must cover both** (user):
(a) under CARTESIAN d/f a shell is not a pure \f$l\f$ — the \f$l-2\f$ contaminants (s inside d, p inside f)
mean "drop a whole orbit" has to reckon with functions that carry two characters at once, which is the same
defect that produced the SR span's two-exponent s window; (b) the rank decision depends critically on
LATTICE SPACING, so a trim validated on one cell says nothing about a denser one — the tests need a
spacing axis, not a single geometry.

## 23. +U is ORBITAL-RESOLVED — U is a vector over (site, shell, site-group irrep); the manifold is an INPUT, never Mn-d by assumption

Ruled 2026-09-16 (user, on Macke et al. JCTC 2024 and ACBN0).  In the eigenbasis of the site occupation
matrix \f$E_U=\sum_i\tfrac{U_i}{2}\lambda_i(1-\lambda_i)\f$; the shell-averaged Dudarev form is the special case
\f$U_i=U\f$ — the same shape as pin 5 (unpolarized is the ζ=0 collapse), applied to +U.  The (t2g, e_g)
split is the site-point-group irrep decomposition, so the labels are fixed by symmetry; eigenvalue
tracking is only for a site symmetry lower than the split.  **Why it is a ruling and not a preference:**
shell-averaging suppresses intrashell screening (perturbing t2g and e_g together zeroes the channel that
screens them: FeS₂ U 7.37 → 3.29/2.16 resolved), and the WRONG manifold is worse than the wrong U — the
correction that opened β-MnO₂'s gap was on **O-p_z**, not Mn-d, and correcting FeS₂'s hybridised e_g at all
broke its structure.  User: *"I have seen other examples where O played an unexpected role in TMOs."*  ⇒ the
term takes a LIST of (site, shell, irrep, U); no code path may assume the Hubbard atom is the transition
metal.  Projector = Löwdin OAO on the site block.  U values are never hand-tuned in production (pin 12):
ACBN0-style from our own on-site ERIs, checked against QE `hp.x`.  Record: `doc/OpenWork.md` §1 step 5.
**Addendum 2026-09-21 (increment 2, earned by a wrong table):** the labels have TWO groups and neither is
the cell's.  The SITE group is the declared decoration's Shubnikov stabiliser (σ=None) — the order splits
t2g → a1g + e_g and the labels must see it.  The PARENT group that names a site level "e_g < t2g" is the
point group of the site's **coordination environment** (`Lattice_3D::SiteEnvironmentRotations`), NOT the
(super)cell's grey stabiliser: on the rhombohedral AFM-II MnO cell the latter is D_3d (12 ops) with or
without decoration and names nothing — measured, after the tree had asserted O_h for a week.  And the
occupation matrix is NEVER symmetrised: symmetry NAMES the eigenvectors of the density's own n (isotypic
projectors, `purity` printed), it does not edit them — a free run's broken symmetry must keep its own
occupations, and the functional must stay dE/dD.  Inside a degenerate cluster the eigenbasis is rotated
to the projectors (n is unchanged); inside a NEARLY degenerate one (four λ≈0.999 on a full majority
shell) the names are ill-conditioned by nature — that is Macke's tracking problem, and the printed
`parentage` says so rather than hiding it.

---

**Where these came from.**  1, 3, 5, 7, 8, 9, 10, 12 were `doc/OldPlans/GPWPlan.md`'s pins section (2026-07).  11 is the user's `UseChargeDensity` post-mortem (2026-09-08).
13 is the KP-0 multi-k defect (2026-09-09).  14–17 were harvested 2026-09-16 when their plan files went RECORD
(`BasisSetTaxonomyPlan`, `GPWPlan1`, `SphericalLatticePlan`, `RunReportPlan`); pin 10's anchor rule came from
`TestSuitePlan` the same day.  18–22 were harvested from the v2 `OpenWork.md` when it was rebuilt as v3
(2026-09-16, `doc/Records/OpenWork_History4.md`) — the ⛔ findings that were durable rather than in the weeds.  23 is the DFT+U ruling of the same day.
2, 4, 6 are user rulings recorded in session memory (`feedback_everything_is_a_fit`,
`feedback_integrated_observables`, `feedback_pw_fitting_uniform_interface`) and had no home in the repo
until now.
