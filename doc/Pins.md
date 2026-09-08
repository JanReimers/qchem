# Durable pins — the invariants a session must not violate

**Cut 2026-09-08 out of `doc/GPWPlan.md`'s "Durable pins / invariants" section**, because that file is a
RECORD of the 2026-07 campaign and the pins are not: they govern work all over the tree, and burying
project-wide invariants inside a finished campaign's plan is how they get missed.  Two user rulings that
lived only in session memory are folded in and marked as such.

These are **rulings, not preferences.**  Each is here because violating it produced a wrong number, a wrong
interface, or a retracted verdict at least once.  If you think one is wrong, say so and get it changed —
do not work around it quietly.  Coding conventions (naming, includes, ownership, style) live in
`CLAUDE.md`; this file is physics and numerics.

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

## 11. No grad-student knobs

Policy enums, not numeric dials.  A number a user has to tune is a design failure looking for somewhere to
live.

---

**Where these came from.**  1, 3, 5, 7, 8, 9, 10, 11 were `doc/GPWPlan.md`'s pins section (2026-07).
2, 4, 6 are user rulings recorded in session memory (`feedback_everything_is_a_fit`,
`feedback_integrated_observables`, `feedback_pw_fitting_uniform_interface`) and had no home in the repo
until now.
