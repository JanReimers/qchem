# Durable pins — invariants worth knowing before you touch the code

Each pin is here because ignoring it produced a wrong number, a wrong interface, or a retracted verdict at
least once.  They are the project's current best understanding, not scripture: everything is open to debate,
and if one looks wrong, say so and change it rather than working around it.  Pin numbers are cited from
code, docs and memory, so they are stable; a pin that was folded into another stays as a one-line pointer.

**Three files, three questions.**  `CLAUDE.md`: how do I work here (conventions, build/test, tool paths, the
doc system).  This file: what must the code obey (physics, numerics, design invariants; one short section
each, pointing at the record that earned it).  A RECORD in `doc/`: why is it this way (evidence, rejected
alternatives).

---

## 1. No cut in r space

> "THERE IS NO CUT in r space, Gibbs ringing is like a wrecking ball." (user, 2026-07-16)

A real-space lattice sum is an ε-converged series for a FIXED operator.  Magnitude screening is the only
truncation, and a radius must not appear as a parameter, member or concept in any interface.

- A truncation radius gives a DIFFERENT operator, not "the operator to ε": a sharp cut is a step function
  in r, its transform rings, and the ringing lands on everything downstream.  Measured: an Rcut=2a NaF
  metric lost 2.25 e per mid-slosh loading.
- A radius is never a conditioning crutch; that job belongs to the basis or to rank reduction.
- The G direction does have a cut (Ecut), but a different kind: a projection onto a finite auxiliary
  subspace, which is variational, exponentially controlled and systematically improvable.  What is banned
  is the real-space cut and its Gibbs ringing.
- End state: one knob per direction — ε in r, Ecut in G.

## 2. Densities, Hartree and Vxc are always fits (grid × fit basis)

For software consistency, every DFT route to the charge density, the Hartree term and Vxc goes through the
fitting framework, even when the fit is trivially exact.  Examples: plane-wave orbitals with a plane-wave ρ
(Parseval makes the fit exact), a fit basis with unit metric, or a trivial projection ⟨ab|fit⟩ — all are
still "a fit".  A quadrature grid is a fit whose basis is (finite-element) delta functions, so a Becke mesh
is {Becke mesh} × {delta basis} and the uniform route is {uniform raster} × {plane waves}.  Both fit the
same \f$v_{xc}\f$, which is what makes one scoreboard legitimate.

Interface consequences (formerly pin 6): PW fitting must look identical to molecular fitting at interface
level.  No "Fourier" in an abstract face.  Any fit/aux basis comes from the orbital basis via
`Create{CD,Vxc}FitBasisSet(...)`; the factory is the seam even when it is trivial, so inlining the trivial
G-space fit was rejected.  Do not assume `orbital == fit`.

## 3. Fit quality is grid-convergence of ρ, never \f$\Delta E_{total}\f$

The fit is non-variational, so a total-energy difference does not bound its error and can flatter it.
Score a fit by the convergence of ρ (or of the operator actually diagonalised,
\f$\max|\Delta V_{xc}(i,j)|\f$) against a fine reference in the same family.

## 4. Report an integrated observable, not a point sample of a field

(user, 2026-08-19.)  MnO's \f$m_{stag}\f$ sampled at one point 0.7 bohr along +x is a spin density, not a
moment; it was never derived, and its direction-dependence is the d-occupation confounder.  Use it only as a
collapse detector.  A frozen-density point probe also understates the self-consistent shift on a metal,
because grid error feeds back through the density, the Fermi level and the occupations.

## 5. Spin-polarized is the native formulation

Unpolarized is the \f$\zeta=0\f$ collapse, not the base case.  New terms are written spin-native
(`FittedVxcPol` / `FittedVcorrPol`); this governs GGA, +U and non-collinear work to come.

## 6. (folded into 2)

## 7. Ill-conditioning is a basis problem — read `[basis trim]`

Ill-conditioning is not a solver bug.  "LASolver" symptoms have been basis conditioning every time.
`valgen` readily produces valence bases with diffuse functions; put those in a unit cell at ordinary bond
lengths and S(k) becomes ill-conditioned.  **Always read the `[basis trim]` output and take it seriously.**
Use well-conditioned bases for SCF (accuracy levels Low/Medium/High; the old N3/N5 pools were test-only and
are gone).  This is where the project has lost the most time.  Two related rules:

- **Auto-trimming is a vet-stage, symmetry-equivariant, whole-orbit decision on S** (user, 2026-08-14/15/26,
  formerly pin 22).  (a) It happens before anything downstream is built (grid ladder, collocation task
  lists, KB projections follow the surviving functions).  (b) The rank decision is a property of S, made
  once, not per spin channel.  (c) Drop whole orbits under the (magnetic) space group, never single AOs:
  per-function pivoting breaks symmetry-tied pivots by numerical noise (runs 58–60 dropped O₁'s p(0.18) but
  O₂'s s(0.15)), and a partial orbit is a symmetry-broken basis.  Report the decision as species / shell /
  exponent, not bare indices.  Two things make it harder than it looks, and unit tests should cover both:
  under Cartesian d/f a shell is not pure l (s inside d, p inside f), and the rank depends on lattice
  spacing, so a trim validated on one cell says nothing about a denser one.
- **The ortho-time path has two behaviours: work, or fail with enough information to fix the BASIS** — never
  a silent third that edits the span.  Today `LASolverCholeskyPivoted` drops and prints a bare index
  ("dropped AO index 47") from a layer that knows nothing about shells.  The proposed cure is an exception
  carrying indices and pivots, caught by the layer that owns the basis (`throw` is this tree's marker).
  Open work: the vet-stage trim itself (`doc/OpenWork.md`).  Record: `OpenWork_History4.md` "Continuous — CLEANUP".
- **Pivoted Cholesky where D is PSD; a failing Cholesky is the canary** (user, 2026-08-20, formerly pin 21).
  Factor density matrices by pivoted Cholesky rather than a trimmed eigendecomposition, but only where the
  factor is not inverted (where it is only multiplied, eigen is fine).  D is PSD here only because no mixer
  extrapolates D (`LinearMixer` is convex, Pulay/Broyden act on \f$\tilde\rho\f$) — a property of today's
  mixers, not a theorem; density-space Pulay, α>1 or MP/cold smearing would break it.  So a failing
  Cholesky means D left the cone: make it loud and route to the eigen split.  PSD tests need a relative
  floor, since an eigensolver always returns O(ε·λmax) negatives.

## 8. One scheme for periodic lattice matrices

GPW builds its lattice matrices from complete-Bloch analytic single sums (scheme A, correct as
\f$R_{cut}\to\infty\f$).  The alternative, truncated-Bloch collocation Gram matrices (scheme B, always PSD),
is a different operator.  Mixing them (scheme-B overlap with scheme-A kinetic) gave \f$E_{kin}=-300\f$.  Take
all matrices of one calculation from one scheme.  (Kept as a one-line historical guard; drop it if scheme B
is gone from the tree.)

## 9. Change a basis interface only for a new integral type

Correctness first, then efficiency, then end-user convenience, then developer convenience, then readability.
(Pseudopotential-smoothness and GAPW remarks that lived here were history; GAPW is out of scope.)

## 10. Two kinds of energy pin — say which, and cite the reference

An energy pin in a test is one of two things, and the comment should say which:
1. **Absolute**: compared to an independent oracle.  Give the reference (literature, or an in-house
   CP2K/QE/ABINIT run with its deck).
2. **Relative ("did E move")**: compared to the previous commit's converged value, with no `Converged()` guard.

A relative pin that moves is re-judged against an independent route and never merely refreshed (KP-0 was
re-pinned −7.45137 → −7.45294 only after the band-folding-equivalent Γ supercell agreed;
`doc/Records/TestSuitePlan.md` §2).  Energy anchors can go stale; a failing charge, count or weight sum is
physics and cannot (user, 2026-09-09).

## 11. Caching: static things in the DB cache, per-iteration things in an explicit phase

Caching trades RAM for runtime.  Cache only quantities that are static across SCF iterations (integrals
keyed by geometry and basis), and use the existing framework: `src/BasisSet/Internal/{DB_Cache,Cache2,Cache3,Cache4}.C`.

The one known exception is the \f$H_{ij}\f$ matrices, which the Fock pass builds and the energy pass reuses
(`tDynamic_HT_Imp::itsCache`).  `DB_Cache` is not the home for them: it never evicts and is built for
sharing across runs, while \f$H_{ij}\f$ turns over every iteration, so twenty iterations would leave twenty
generations in it.  Ask what a cache evicts before asking what it keys on.  The intended fix is a
per-iteration scope that owns those matrices and dies with the iteration (`CleanupHistory3.md` R1.0h).
I know of no other exception; if you add one, say why at the declaration.

Underlying lesson (user, 2026-09-08, the `UseChargeDensity` post-mortem): replacing an explicit "here is the
density, prepare yourself" call with automatic per-object freshness guards produced bugs (a memo believing
it was fresh; `itsRho`/`itsXCMix` aliasing that made α=0.25 and α=1.0 bit-identical) and blocked threading
the per-block loop (write-on-first-touch).  When work must happen once per density, iteration or geometry,
give it a named phase rather than hiding it in call order.

## 12. The code computes physics numbers; the user supplies scope (provisional)

A number the physics determines should never be a user dial.  The origin: an early DFT+U asked the user for
U, which would stall a new student for weeks.  The answer is that the code learns to compute U (and J, V)
self-consistently; this is hard, but the know-how accumulates in the code.  What the user still supplies is
SCOPE: (1) which atoms/shells get a correction, (2) which energy window defines the correlated subspace.
That is still a real decision, but a much smaller one.  See pin 23.

Boundary between compiled code and the (planned) GUI agent, as a working proposal rather than a ruling:
- Code owns anything deterministic and reproducible that a test can pin: computing U, mixing parameters,
  convergence decisions.
- The agent/user owns choices that depend on scientific intent.
- For scope choices the code should still *surface the evidence* as structured output (candidate windows
  from a DOS minimum or a projected-character gap, candidate sites by projected d/f weight, each with its
  reason), so the agent and the GUI both consume one API and a recommendation is reproducible.  The code
  recommends; the user confirms.

## 13. A symmetry op acts on grid indices as \f$DUD^{-1}\f$, not as \f$U\f$

A grid point is \f$k=(i+s)/N\f$ componentwise, so an op \f$U\f$ acts on the index lattice by
\f$M=DUD^{-1}\f$ with \f$D=\mathrm{diag}(N)\f$ (\f$i' = M(i+s)-s\f$).  \f$M=U\f$ only on an isotropic mesh,
which is why applying \f$U\f$ to indices works on every \f$n\times n\times n\f$ mesh and silently fails on
\f$2\times1\times1\f$.  A mod-N wrap does not validate a map: reducing a stray image into \f$[0,N)\f$ always
yields some grid point.  Test the OP once: \f$M\f$ must be integral (then it is unimodular, so a
permutation), and \f$(M-I)s\f$ integral for a shifted mesh.  An op failing either is not a symmetry of
that mesh and must be dropped whole.  Cost: IBZ stars overlapped on Si \f$2\times1\times1\f$,
\f$\Sigma w=1.5\f$, 12 electrons in an 8-electron cell (KP-0; `OpenWork_History3.md`).  Still open: the
group that symmetrizes ρ must be intersected with the mesh symmetries (`OOD-SOLID-Cleanup.md` R1.0r).

## 14. Libraries follow the basis FAMILY, modules carry the symmetry GROUP

`BasisSet = ⊕_irreps IrrepBasisSet`, each carrying one irrep label of G, the Hamiltonian's symmetry group.
Three orthogonal axes: **G** (block labels, `qcSymmetry`), **family** (Gaussian, Slater, BSpline, PW, delta;
this is the integral engine), **construction** (subduce or induce; derived).  A library is an engine (its
integrals factorise over the family, so `UnitCell` inside the Gaussian engine is legitimate); a module
carries the G (`…Gaussian.Point.*` vs `…Gaussian.Lattice.*`), enforced by `scripts/audit-basisset-gtags`
(no `.Point.` module imports `qchem.UnitCell` or `qchem.Symmetry.Lattice_3D.*`).  Spin is a factor of G until
a double group dissolves it, so Pol/UnPol is an imposed subgroup (V1.37), not a type.  Cost: the tree was
cut on the wrong axis and V1.33 re-cut it in eleven commits.  Record: `BasisSetTaxonomyPlan.md` §1.

## 15. Smearing needs kT above the frontier splitting; GDM as built is fixed-occupation

(a) kT must exceed the frontier splitting or the occupations slosh instead of converging (NaF: 1e-2
converges, 1e-3 does not).  (b) The GDM minimiser diverges under smearing, because its direction is the
fixed-occupation \f$[F,D]\f$ and not the free-energy gradient.  So smeared runs use the fixed-point stage
(DIIS/Kerker/Pulay) and GDM polishes only at kT=0 until it has a smearing-aware direction.  This is a limit of
our parameterisation, not of the method (CP2K's OT has the same axis scaffolded): do not conflate
hold-the-block with don't-smear.  Records: `GPWPlan1.md` (2026-07-26), `SCFStrategyPlan.md`.

## 16. Match basis spans before comparing to an oracle

MnO with Cartesian d (s-contaminants): FM 40 mHa below AFM.  Same cell with spherical d: AFM 45.5 mHa below
FM.  No code bug, only the span.  The extra s-like functions let the density dodge the l=0 KB projectors,
which the oracle's spherical basis cannot do, and an earlier "8 mHa agreement" with CP2K was compensation
between contaminants and diffuse functions.  So: match both spans exponent-for-exponent before comparing
energies (`valence_lowq_sph` v2 is the CP2K transcription), and don't use a Cartesian-d basis for a d-metal
ordering question.  Record: `SphericalLatticePlan.md` I0–I2.

## 17. A class reports contemporaneously with its own activity

`CurrentReport` is a global sink; each class emits when it does the thing, so console order is execution
order.  An `Emit*()` method on an abstract face, or a reporter that pulls state out of objects afterward,
is the defect (user, 2026-09-11; V1.5 deleted the `Emit*()` faces).  Corollary: a printed trace number is
either physics or a gate the run consumes; printing \f$\alpha_{eff}\f$ implied it was used, so it was
deleted (user, 2026-09-13).  Records: `RunReportPlan.md`, `CleanupHistory3.md` V1.5.

## 18. A separately damped XC feed destroys Kerker's mode selectivity

Four collapsed MnO states (−45.5, −46.3, −56.4, −38.5 vs the converged −61.403) had one cause, a monotone
dose-response in \f$E_{ee}\f$ (13.5 → 29.0 → 35.1): giving \f$V_{xc}\f$ a density damped by a flat
\f$\alpha_{eff}\f$ un-damps the low-G charge mode 2.4× and barely touches the AFM mode.  Kerker's low-G
charge-slosh damping is what holds the AFM basin; the moment collapse follows from the charge runaway and
is not a spin effect.  So the XC feed is never a second, independently damped copy of the density.  The
ρ≥0 goal survives (the DM route gives 0 negative points against 15% for the fitted ρ̃) via the cusp
deficit \f$\rho_{XC}=\tilde\rho_{mix}+(\rho[D]-\tilde\rho_{fit})\f$, with no \f$\alpha_{eff}\f$ to choose (N4,
still to be measured; pin 28 makes it the RULE on a field mixer).  Also: \f$\tilde\rho_{mix}\f$ and \f$\rho[D]\f$ have different fixed points (ρ̃ is a
fit projection; NaF 139 μHa), so "at convergence they agree" is false, and \f$E_{ee}\f$ is a
validated charge-slosh detector (T3).  Record: `OpenWork_History4.md` "ITEM 1 MEASURED" + "N4".

## 19. Never density-screen the gather — \f$h_{ij}\f$ is diagonalised, not traced

Dropping a term because \f$D_{ij}=0\f$ is sound for the energy (\f$\mathrm{Tr}(Dh)\f$ is blind to it) and
wrong for the Fock matrix, which is diagonalised to make the next density.  Zeroing \f$h_{ij}\f$ where the
density vanishes is self-fulfilling: a pair with no density can never acquire any, and with a diagonal SAD
seed that is every off-diagonal element.  CP2K does not D-screen at all (`task_list_methods.F`: one global
`eps_rho_rspace`, geometry only).  Flooring a vanishing weight broke the stream fold's orbit invariance
(0.14 vs 2.4e-8), so the gather's D-screen was removed (2026-09-04); the D-aware tolerance survives only on
the collocation, where the weight really is the scatter weight.  `cij==0` had meant both "structurally
absent" and "zero density"; only the first may be excluded.  Record: History4 "THE k-SCALING GAP" + "ATTEMPT 2".

## 20. Ask what a matrix means before symmetrizing it

The first T3 stream fold orbit-averaged `screenD`, a matrix of \f$|D_{ij}|\f$ magnitudes, not a density
matrix.  Signed averaging cancelled mixed-σ orbits to ~0, the D-aware screen dropped live terms, and the
imposed O₂ triplet collapsed by 2.3 Ha.  A screen is reduced by the orbit MAX (`FoldScreenMax`), a density
by the orbit projection (`FoldProjectedD`); reading the representative's own \f$D_{ij}\f$ samples the orbit,
which equals projecting only if D is already symmetric.  Gate: the dimer-in-a-box cell (single-atom Si
cells miss this bug).  Record: History4 "Step 2 — ARM THE SYMMETRY FOLDS".

## 21, 22. (folded into 7)

## 23. DFT+U lessons (summary; evidence in `doc/Records/HubbardUHistory.md`)

Plan: `doc/HubbardUPlan.md`.  The earlier long texts of pins 23–25 are archived in the history record.

- **+U is orbital-resolved.**  U is a vector over (site, shell, site-group irrep); the term takes a list of
  (site, shell, irrep, U) and no code path may assume the Hubbard atom is the transition metal.  Shell-averaging
  suppresses intrashell screening (FeS₂ U 7.37 → 3.29/2.16 resolved), and the wrong manifold is worse than the
  wrong U: β-MnO₂'s gap opened with a correction on O-p_z, not Mn-d.  Projector = Löwdin OAO on the site block.
- **The manifold and the energy window are inputs, not derivations.**  With entangled bands (NiO Ni-3d/O-2p)
  there is no window-independent "correlated subspace": cRPA and LRT disagree 16× on Sr₂FeO₄ from the window
  alone (Carta et al., arXiv:2505.03698).  The code should recommend a window where a natural one exists
  (pin 12), not hide the choice.  LRT and cRPA belong behind the same estimator face as ACBN0.
- **Labels: two groups, neither the cell's.**  The site group is the declared decoration's Shubnikov
  stabiliser; the parent group that names "e_g < t2g" is the point group of the site's coordination
  environment (`Lattice_3D::SiteEnvironmentRotations`), not the cell's (D_3d on rhombohedral AFM-II MnO,
  which names nothing).  The occupation matrix is never symmetrised: symmetry names its eigenvectors, it
  does not edit them.  Nearly degenerate clusters have ill-conditioned names (Macke's tracking problem).
- **A linear-response U is conditioned on the state it linearises about.**  NiO's hp.x 5.267 eV is
  U_LR(U_in=3 eV) (5.434 at U_in=0); quote U_in beside every response value.  When χ₀ ≈ χ (MnO d⁵, ZnO d¹⁰)
  same-site U is a difference of near-equal small numbers and not an oracle.
- **Two comparisons, never merged.**  Ours ÷ published-ACBN0 (1.8–2.8×) measures projector completeness
  (their minimal PAO keeps ~60% of the norm), not screening.  Only ours ÷ an independent matched-PP oracle
  (hp.x) tests screening; label each oracle row matched-PP or different-PP.  Never build a gate on an
  ACBN0-derived target (retracted 2026-09-23; screened-ACBN0 refuted 2026-09-25).
- **An iteration cap is not a verdict.**  Three wrong conclusions here were a capped run read as "does not
  converge".

## 24, 25. (folded into 23)

## 26. CP2K DFT+U reports `trq` scaled by `fspin`

`dft_plus_u.F` multiplies the reported manifold trace by `fspin` (0.5 for RKS, 1.0 for UKS), so
`E_DFT+U = alpha*trq` read at face value gives χ 2× too small on a closed-shell oracle run; double it (exact
when U_MINUS_J=0).  Earned by the CK-alpha Si cross-check (χ −14.5164 vs ours −14.5117 only after the fix).
Records: `HubbardUPlan.md`, `OpenWork_History5.md` §2.

## 27. A Kerker-type filter `g²/(g²+G0²)` must define its G=0 value

With a leaf's `G0=0` (the undamped m channel of (ρ,m) mixing) it is 0/0 at G=0 and the NaN rides into the
rebuilt channels and v_xc.  Guard: `g2>0 ? … : 1.0` (full mixing, not frozen); gate
`GPW_Si.Γ_Imp_Pol_Kerker_eqUnpol`.  Record: `OpenWork_History5.md` ("(ρ,m) Kerker on an exact SINGLET").

## 28. XC is evaluated from ρ[D], never from a fitted ρ

> "Calculating XC from any 'fitted' ρ has caused massive problems in the past." (user, 2026-10-22)

\f$V_{xc}\f$ and \f$E_{xc}\f$ see the density of the density matrix, \f$\rho[D](r)=\sum_{ij}D_{ij}\chi_i(r)\chi_j^*(r)\f$
evaluated pointwise on the XC quadrature -- never the FIT's representation of it, whatever the fit basis (\f$\tilde\rho_{fit}\f$,
or a mixer's \f$\tilde\rho_{mix}\f$).  A fitted ρ misses whatever its basis cannot hold -- on a plane-wave fit, the core cusp,
with Gibbs ringing (pin 1): 15% of MnO's XC points came out negative (min −0.154) against none on the ρ[D] route -- and it has
a DIFFERENT fixed point: NaF 139 μHa, the O₂ box on O q6 **+7.72 mHa** (−30.70246 vs −30.71018; 2026-10-22).  At kT = 0 the
factored D (low rank, `policy.dmLowRank`) makes ρ[D] cheap -- 130× cheaper per sampling than the matrix-free fitted ρ on MnO.
(Pin 2 is a different statement: it routes every density through the fit INTERFACE; this pin says which ρ XC must see.)

**Two exceptions (user, 2026-10-22), both transient:**
1. **The SAD seed** has no D -- one iteration.
2. **A G-dependent (Kerker-type) preconditioner** exists only in the fit's {G} representation, so its step cannot be taken on D.
   During such a step XC sees \f$\rho_{XC}=\tilde\rho_{mix}+(\rho[D]-\tilde\rho_{fit})\f$ (the CUSP DEFICIT, N4): the damped fitted
   part plus what the fit cannot hold, taken from the current D.  At the fixed point \f$\tilde\rho_{mix}=\tilde\rho_{fit}\f$, so it IS
   ρ[D] (O₂ box: −30.71018173, the ρ[D] value to all digits).  Bare \f$\tilde\rho_{mix}\f$ is NOT this exception -- that is the
   +7.72 mHa.  Nor is the raw output ρ[D] alone: it strips the filter's damping from the XC channel (MnO AFM-II collapsed to
   −45.5 Ha vs −61.403; pin 18).  The part outside the fit needs no damping: no mixer acts outside the fit space.

- A mixer with NO G-dependent filter mixes D itself, history and all -- \f$\rho[\sum_i c_iD_i]=\sum_i c_i\rho[D_i]\f$ pointwise -- so XC
  sees \f$\rho[D_{mix}]\f$ with no correction.  Linear D-mixing and direct minimisation (GDM) comply (O₂ box: linear −30.71018167/173).
  A Pulay history without Kerker (G0 = 0) is in this class too, but today's `PulayMixer` still mixes the fitted ρ̃ (OpenWork §3).
- The cusp deficit is wired on Kerker and Kerker-preconditioned Pulay, unpolarized and in both spin coordinates, (ρ↑,ρ↓) and
  (ρ,m), and is the DEFAULT (`SCFParams::XCCuspDeficit` = true, 2026-10-22); false is a control, never a recipe.  The raster (pair)
  XC route meets it with its raw shadow (ρ[D] collocated pointwise, mixed on the raster) and refuses the flag on its BALL (fitted)
  fallback.  Free MnO AFM-II holds its order with it (moment 4.453); the flip moved no anchor's verdict (`OpenWork_History5.md`).

Records: `OpenWork_History4.md` "ITEM 1 MEASURED" (MnO negatives, the collapse), `OpenWork.md` §3 "Kerker/Pulay as the
DEFAULT mixer" and N4 (the O₂ and MnO arms).
