# Durable pins — the invariants a session must not violate

**Cut 2026-09-08 out of `doc/GPWPlan.md`'s "Durable pins / invariants" section**, because that file is a
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
`doc/TestSuitePlan.md` §2): KP-0 re-pinned \f$-7.45137\to-7.45294\f$ only after the band-folding-equivalent Γ
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
carried 12 electrons in an 8-electron cell (KP-0, 2026-09-09; record in `doc/OpenWork_History3.md`).
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
`doc/BasisSetTaxonomyPlan.md` §1; Doxygen `\ref basisset_taxonomy`.

## 15. Smearing needs kT ABOVE the frontier splitting, and GDM as built is fixed-occupation

Two measured facts from the Fermi-smearing build (`doc/GPWPlan1.md`, 2026-07-26): **(a)** kT must EXCEED the
frontier splitting or the occupations slosh-rotate instead of converging (NaF: 1e-2 converges, 1e-3 does
not); **(b)** the GDM direct minimiser DIVERGES under smearing because its geodesic direction is the
fixed-occupation \f$[F,D]\f$, not the free-energy gradient (which carries an occupation-response term).
⇒ smeared runs use the fixed-point stage (DIIS/Kerker/Pulay); GDM tail-polishes only at kT=0 until it has a
smearing-aware direction.  ⚠ This is a limit of OUR parameterisation, not of the method — CP2K's OT has the
same axis scaffolded; **never conflate hold-the-block with don't-smear** (`doc/SCFStrategyPlan.md`).

## 16. A basis SPAN can reverse a magnetic ordering — match spans before comparing to an oracle

MnO with Cartesian d (its \f$r^2e^{-\alpha r^2}\f$ s-contaminants): FM below AFM by 40 mHa.  The same
cell through the spherical-d view: AFM below FM by 45.5 mHa — **no code bug, the span**.  Mechanism: the
extra s-like functions let the density rearrange s-character away from the l=0 KB projectors (an l=0
repulsion DODGE worth 37% of the weak basin's reward), a freedom the oracle's spherical basis structurally
lacks.  And the earlier "8 mHa agreement" with CP2K was contaminants-vs-diffuse COMPENSATION.  ⇒ Before
any energy is compared to an oracle, the two spans are matched exponent-for-exponent
(`valence_lowq_sph` v2 = the CP2K transcription), and a Cartesian-d basis is never used for a d-metal
ordering question.  Record: `doc/SphericalLatticePlan.md` I0–I2.

## 17. A class reports CONTEMPORANEOUSLY with its own activity — console order == execution order

`CurrentReport` is a global sink; each class emits at the moment it does the thing, and never tells
another class to emit.  An `Emit*()` method on an abstract face, or a "reporter" that PULLS state out of
objects after the fact, is the defect (user, 2026-09-11; V1.5 deleted the `Emit*()` faces).  Corollary for
trace columns: a printed number is either physics or a gate the run CONSUMES — printing \f$\alpha_{eff}\f$
implied it was used, and it was deleted for that reason (user, 2026-09-13).  Record: `doc/RunReportPlan.md`
(the design), `doc/CleanupCandidates.md` V1.5.

---

**Where these came from.**  1, 3, 5, 7, 8, 9, 10, 12 were `doc/GPWPlan.md`'s pins section (2026-07).  11 is the user's `UseChargeDensity` post-mortem (2026-09-08).
13 is the KP-0 multi-k defect (2026-09-09).  14–17 were harvested 2026-09-16 when their plan files went RECORD
(`BasisSetTaxonomyPlan`, `GPWPlan1`, `SphericalLatticePlan`, `RunReportPlan`); pin 10's anchor rule came from
`TestSuitePlan` the same day.
2, 4, 6 are user rulings recorded in session memory (`feedback_everything_is_a_fit`,
`feedback_integrated_observables`, `feedback_pw_fitting_uniform_interface`) and had no home in the repo
until now.
