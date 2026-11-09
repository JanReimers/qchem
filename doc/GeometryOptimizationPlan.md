# Geometry optimization: molecules first, on analytic forces (DESIGN NOTE)

Lane B, 2026-11-07, rewritten the same day after the user's rulings (§12).  Status: a design for the user to read before any code.

It covers:
- the OpenWork rows "A cheap host for geometry-optimization experience" (the user's sequence: molecules, then rutile TiO₂, LiCoO₂,
  Li_xTiS₂) and level 1 of "Geometry optimization: imposed relaxation, then DISCOVER…";
- the forces seam of `doc/Records/CalcNetworkReview.md` §3.3 (the `PotentialEnergySurface` face, §3.1).

Evidence it cites:
- the seed test (`doc/Records/SeedTestReport.md`);
- a PySCF probe on H₂O, `~/Code/qchem6-runs/H2O/fdforces/`.

The first version proposed a finite-difference force product.  The user ruled it out (§12.3), so this version designs the ANALYTIC
forces instead.  FD survives only as a per-term unit-test helper at fixed D (§6).

## 0. What exists, and the shape

**Already there:**
- `qchem.PotentialEnergySurface`: `ConvergeAt(cell, R, want) → Outcome<PESPoint, PESFailure>`, forces optional in `PESPoint`.
- `MoleculeSession` and `SolidSession` implement it.  `MoleculeSession` carries each point's D block by block onto the next point.
- A forces request returns `NoGradient` today.

```
   GeometryOptimizer ── q, E, dE/dq ──► CoordinateMap ──► PotentialEnergySurface ──► MoleculeSession │ SolidSession
   (BFGS, trust radius;                 (molecule: Cartesian;   (energy, forces)          │ at the converged D:
    model-surface unit tests)            solid: cell + Wyckoff)                           ▼
                                                         F = − Σ_terms ∂E_t/∂R |_D  − Tr(W ∂S/∂R)    (§4)
```

## 1. The optimizer, behind a coordinate map (ruled Q1)

**The optimizer sees only a coordinate vector q, E(q) and g = dE/dq.**  It never sees R or a cell.

A `CoordinateMap` turns q into a geometry (cell, R) and turns a geometry's forces into g = −Jᵀ F, with J = ∂R/∂q:
- **Molecule:** Cartesian.  q = R, J = 1, with the 6 (5 if linear) external modes projected out of g and of the step.
- **Solid, level 1:** the cell plus the free Wyckoff parameters of the imposed group.  J is the group's own map from parameters to
  positions, so symmetrized forces come for free.  The cell part needs stress (the variable-cell seam, `CellChange` today).
- **Later:** redundant internals could be a second molecular map, and nothing above the map changes.

**The method.**
- Quasi-Newton BFGS, started from a model Hessian (Lindh's force-field form, which needs only Z and R).
- A rational-function or trust-region step, starting at 0.3 bohr, shrunk on an energy rise and that step rejected.  No line search.
- A `PESFailure` at a point stops the optimization with an Outcome naming the point and why.

**Convergence** is max|g|, RMS g, max step and ΔE.  The force threshold is DERIVED from ε through a named calibration (ruled Q2,
§8 has the calibration).

**Tests:** analytic model surfaces only (a Morse or LJ cluster through a test PES, a quadratic with a known Hessian, a map with a
nontrivial J).  They check convergence, step counts, invariance and refusals.  **No real optimization runs before analytic forces
exist** (ruled §12.3).

## 2. Grids fixed across one optimization, picked once from the two-step run

The rule (the user's fixed-settings and one-grid-set-per-series rulings, applied to one optimization):
1. **Step 1:** seed (SAD), ε = 1e-4, converge at the start geometry.
2. **Audit at that D on the H measure:** max|ΔH_ij| against the next rung up, the pilot's rule.  A fixed fine reference costs
   more than it buys: the seed test's addendum replays the rule, 6 of 6 meet ε.  The ε here is the one the force target needs (§8).
3. **Fix the RULE for the whole optimization:** (nR, L, angular kind, radial map, Becke order) per species, the Hartree fit basis
   and every raster.  The mesh POINTS move with the atoms; its parameters never change.  The basis recipe and the PP set are frozen,
   as `MoleculeSession` already freezes them.
4. **Never resize mid-optimization.**  A resize makes E(R) discontinuous and the BFGS update meaningless.
5. **Report health per point:** δq of the mesh, free at every point.  At the converged geometry, audit again and report whether
   the pick differs.  Acting on a different pick is the NEXT run's business (pin 31 catch 2).

Evidence: on all five seed-test systems the 1e-5 pick made at the 1e-4 D is the converged D's pick.  That test did not cover a
geometry change between the steps; step 5's final audit checks it, for every optimization.

## 3. The force: what it is, and what it needs

At a converged density, with the energy E[D, R] and the orthonormality constraint, for each atom A:

    F_A = − Σ_t ∂E_t/∂R_A |_D  +  Tr(W ∂S/∂R_A)

- W = Σ_σ Σ_i f_iσ ε_iσ c_iσ c_iσᵀ is the energy-weighted density, spin-summed.  The formula is spin-native, with D = D↑ + D↓
  wherever a term reads the total.
- Each ∂E_t/∂R_A |_D is the term's OWN derivative at frozen D: its operator's centres and its basis functions' centres move, and
  D's AO coefficients do not.  That is exactly what a per-term finite difference at fixed D measures, which is why it is each
  term's unit oracle (§6).
- The formula holds only where the SCF stopped at a stationary point of E in D.  So **every term must give H = ∂E/∂D first** (pin 32):
  a term in an inconsistent shape has no analytic force, only `NoGradient`.

## 4. Every force term, its oracle, and its estimate

Routes as the molecular factory builds them (`src/Hamiltonian/Internal/Imp/Hamiltonians.C`).  **Default stack** = the first route of
each row.  Estimates are lane-days, with the shared layer of §5 paid once (2.5 days).

| term (class) | the gradient | oracle: unit (FD at fixed D, no SCF) + integrals | whole-molecule oracle | pin 32 | estimate |
|---|---|---|---|---|---|
| **Pulay** (basis, through S) | +Tr(W ∂S/∂R): ⟨∂χ\|χ⟩ | FD of Tr(W S(R)) at fixed W; libcint `int1e_ipovlp` | PySCF | — | 0.5 |
| **Kinetic** | −Tr(D ∂T/∂R) | FD of Tr(D T(R)); `int1e_ipkin` | PySCF | static, consistent | 0.25 |
| **Ven** (all-electron nuclear) | basis-centre part + operator part (the nucleus moves) | FD; `int1e_ipnuc` (centres) + `int1e_iprinv` (operator, origin on the nucleus) | PySCF | static | 0.5 |
| **IonIon** | Coulomb between Zion cores, closed form | FD of the sum; translational invariance | PySCF | none | 0.1 |
| **PP local** (GTH, closed form when every block `CanAssembleField`) | three-centre: two basis centres + the PP centre (the PP centre's part = − the sum of the other two, by translation invariance) | FD at fixed D; no libcint GTH, so our engine's l±1 identity (∂/∂A_x of a Cartesian Gaussian = 2α G(l+1) − l G(l−1)) on the existing closed form | PySCF GTH (`pseudo=gth-pade`; worked on H₂O, §8) | static | 1.5 |
| **PP local, quadrature route** (a block that cannot assemble the field) | needs the PP mesh's weight derivatives | FD at fixed D | — | static | `NoGradient` first; with the XC weight derivative, 0.5 |
| **PP non-local** (KB projectors) | Σ D_ij h_pq ∂(⟨χ_i\|p⟩⟨p\|χ_j⟩): two-centre ⟨χ\|p⟩ derivatives | FD at fixed D; l±1 on the closed-form ⟨χ\|p⟩ | PySCF GTH | static | 1.0 |
| **Hartree, `Vee_DFT_Dunlap`** (default) | 2 cᵀ ∂C/∂R · D − cᵀ ∂G/∂R c, plus the charge constraint's λ ∂N_fit/∂R (fit functions move with their atoms) | FD at fixed D, with c re-solved at R±h (stationary in c, so dc/dR drops out analytically: the FD checks that); `int3c2e_ip1` (orbital centres), `int3c2e_ip2` (fit centre), `int2c2e_ip1` | PySCF DF-J with the SAME auxiliary basis (Coulomb metric, no charge constraint: a near oracle, its gap measured once) | **consistent** (Cᵀ G⁻¹ C) | 1.5 |
| Hartree, `Vee_DFT` δ (radial Poisson / direct sum) | qᵀ K ∂q/∂R + ½ qᵀ ∂K/∂R q: the points and weights move, and K (radial Poisson) moves with the sites | FD at fixed D | PySCF (exact J) at a converged grid | consistent (Oᵀ K O), first order | `NoGradient` first; 2-3 |
| Hartree, multipole | as δ, plus the channel expansion | FD at fixed D | — | **implied second derivative not yet checked** | `NoGradient`; check first |
| **XC, `Vxc_DFT` δ on the Becke mesh** (default) | (a) Σ_g w_g v_xc ∂ρ(r_g)/∂R: basis-function derivatives at the points, and the points ride with their atom; (b) the GRID RESPONSE Σ_g ∂w_g/∂R e_xc(ρ_g): the Becke partition's derivative (Johnson–Gill–Pople 1993), through its iterated s(μ) (order 2 = three iterations) | FD at fixed D on the moving mesh, (a) and (b) together; (b) alone as a qcMesh unit test: ∂w_g/∂R by FD of the weights.  Σ_A F_A = 0 is the invariance check (b) is needed for | PySCF with `grid_response=True` where its grid can match ours (original Becke, no radii adjustment), else at a converged grid (the H₂O probe agreed to 3.5e-7 Ha/bohr) | **consistent** (pin 28: H = Oᵀ f′(OD)) | (a) 1.0, (b) 2.0 |
| XC, A1 auxiliary fit (`FittedVxc`) | through the fit | FD at fixed D | — | **suspect** (an overlap-metric fit of v_xc): check before any gradient | `NoGradient` |
| XC, multipole | through the channels | — | — | not checked | `NoGradient` |
| HF `Vee_HF`, `Vxc_HF` (J, K) | 4-centre derivatives, contracted with D and D | FD at fixed D; `int2e_ip1` | PySCF RHF/UHF gradients | exact (no fit) | 1.5 (after DFT) |
| +U (`Hubbard_U`) | projector-occupation derivatives | FD at fixed D | QE / ABINIT for solids | to check | later |

**Default stack total:**
- 2.5 shared (§5);
- 0.5 Pulay, 0.25 kinetic, 0.1 ion-ion, 1.5 PP local, 1.0 PP non-local, 1.5 Hartree Dunlap, 3.0 XC;
- 1.0 gates (§6) and 1.0 the optimizer (§1).

That is **about 12 lane-days to an LDA, GTH, DZVP-MOLOPT molecule optimized on analytic forces.**  The all-electron Ven adds 0.5;
HF adds 1.5.

## 5. Where the code lives, and how it meets lane S's sprint

**The shared layer (2.5 days), touching no tier-1 hub:**
- **A derivative-integral evaluator** built from `AoShellSource` (the shells the basis already exposes).  It packs them for libcint
  as `PG_LibCint` already does, and returns ∂M/∂R_A as 3 matrices per atom for ovlp, kin, nuc/rinv, 3c2e and 2c2e.  Our in-house
  MnD engine gains the l±1 identity for the GTH integrals libcint lacks, and libcint checks it on every integral the two share.
  Building it beside the basis, not on the `Orbital_*` hubs, keeps their faces unchanged.
- **W from the converged wave function**, spin-summed.
- **The force assembly:** `ConvergeAt(forces)` on `MoleculeSession` sums the terms' gradients, adds the Pulay term, reports per term
  (one `report::Row` per term: max|F|, Σ_A F_A), and returns `NoGradient` naming the first term without one.

**The per-term gradient is a CAPABILITY face** (CalcNetworkReview §3.3).
- `GeometryGradient { std::vector<rvec3_t> dE_dR(…) const; }` lives in its own module and is reached by the sanctioned
  abstract→abstract cross-cast.  A term that does not implement it is a `NoGradient` with its name, never a silent zero.
- **How the assembly finds the terms is a tier-1 `Hamiltonian` hub question**: it needs one query "each term, or its gradient face".
  That needs a design discussion with the user (CLAUDE.md, the protected hubs).  §11 asks it.

**The sprint mapping (lane S, `doc/HamiltonianSprintPlan.md`):**
- **Step 3b (landed):** `Vee_DFT_Dunlap` is its own class.  The Dunlap gradient is Dunlap's own formula on Dunlap's integrals, which
  is exactly what V1.1b asked of the two Hartree versions.
- **Step 4, V1.41 option C:** one current state per term, filled at the refresh.  A gradient READS that state at the converged D
  (Dunlap's c and λ, the XC term's v_xc and e_xc at the points), so no term recomputes a fit to give its force.  **The term-side
  gradient code lands after step 4.**  Before it, it would be written against the `*_RealBlock` caches step 4 deletes.
- **Step 5, V-ONE-OF-N:** the coefficient record.  Dunlap's c is an `rvec_t` today and its gradient contracts it with ∂C/∂R, so
  it is written against the record when step 5 has landed, or moves with it.  The δ Hartree route's q on the points is the record's
  other realization, which is one more reason that route waits.
- **Step 5, V-NONLINEAR-MESH:** the term loses mesh ownership.  So the Becke weight derivative belongs to the MESH, not the term: it
  lives in qcMesh as a pure function (∂w_g/∂R_B for an atom-centred mesh), unit-tested against FD of the weights, whoever owns the
  mesh after step 5.
- **Before step 4 lands**, these can be built with no Hamiltonian edit:
  - the derivative-integral evaluator (libcint vs MnD gates);
  - W;
  - the Becke weight derivative in qcMesh;
  - the FD-at-fixed-D test helper;
  - the optimizer and its coordinate maps.

## 6. Tests and oracles

1. **Per-term FD at fixed D** (the helper that replaces the FD product):
   - Freeze D's AO coefficients (and W for Pulay).
   - Evaluate ONE term's energy at R ± h along each Cartesian, h = 2e-3 bohr.
   - Compare with the term's analytic gradient.
   - No SCF, so a unit test in seconds.  Truncation is (h²/6)E‴: about 1e-6 Ha/bohr on H₂O's stretch at that h (§8).
2. **Invariance:** Σ_A F_A = 0 and zero torque, per term where it holds and in total.
   - The XC term alone breaks Σ F = 0 without the grid response (b).  That is the test that (b) is right.
3. **libcint against our engine:** every derivative integral the two share, element by element.
4. **Whole molecule against PySCF:** H₂O, NH₃, CH₄, then benzene, at displaced geometries.
   - Same Cartesian basis, GTH set, functional; DF-J with the same auxiliary basis for the Dunlap route.
   - Then the optimized geometries against PySCF's own optimizer, and our step counts.
5. **The optimizer** on model surfaces only (§1).

## 7. What the molecular facade lacks today

| gap | today | needed |
|---|---|---|
| ε | `CalcOptions::accuracy` is consumed by no term ("consumed by: none yet"); ε is wired on the solid facade only (V-FACADE-ACCURACY) | the molecular XC mesh sized from ε on H: the audit below, its pick fed back as the setup's `mesh` |
| a molecular grid audit | none in our code (the seed test's TeF₆ ladder was PySCF; `SolidCalculation::XCMeshAudit` is solid-only) | an audit at fixed D against the next rung up (the `qchem.Mesh.BeckeClimb` rule, lane B 2026-11-07) |
| two steps across setups | `MoleculeSession` carries D between points of ONE setup | carry D from the loose setup into the tight one (same basis, so D carries) |
| the deck | refuses molecules | molecules in the deck, and an `optimize` post-action beside `displace` |
| basis | the transcribed DZVP-MOLOPT `.bsd` is H C N O only (gated against CP2K 2026-11-02) | enough for the molecular phase; `scripts/cp2k2bsd` for Li, Ti, Co, S, F with the solids |
| forces | `NoGradient` | §4 |

## 8. ε for forces: the H₂O measurement (evidence, kept)

PySCF probe, H₂O LDA, GTH, DZVP-MOLOPT, our molecular rule (MHL m3 α2, Becke order 2), atom-centred meshes moving with the atoms;
force on one H.  The grid force error is F on the rung minus F on (250, 59), both by central differences at the same h:

| rung (nR, L) | points | \|ΔE\| at R₀ (μHa) | max\|ΔH_ij\| | grid force error, stretch / bend (Ha/bohr) | force / H |
|---|---|---|---|---|---|
| (40, 17) | 13,164 | 1.2 | 2.0e-4 | 2.1e-3 / 6e-7 | 10 |
| (60, 23) | 34,864 | 2.0 | 6.8e-6 | 9.5e-5 / 3.0e-5 | 14 |
| **(80, 29) default** | 72,282 | 0.09 | 3.6e-7 | 4.9e-6 / 2.7e-6 | 14 |
| (100, 35) | 129,722 | 0.13 | 6.2e-8 | 4.5e-7 / 4e-8 | 7 |
| (150, 47) | 345,292 | 0.004 | 2.3e-9 | ~1e-8 | — |

What it shows:
1. **On a mesh that moves with the atoms, the grid's force error is a smooth BIAS:** it is the same at h = 1e-3 to 1e-2 to 1%.  An
   analytic force on that mesh carries the same bias, since it differentiates the same E(R).
2. **The force error follows max|ΔH_ij| (7–14×) and not ΔE.**  (40, 17) has the smaller ΔE of the first two rungs and the 20× larger
   force error.  That is pin 31's measure clause, measured on a force.
3. **The calibration (ruled Q2):** f_tol ≈ η·ε_H, with η = 15 measured on ONE molecule and two coordinates.
   - It is a named calibration with this record at its declaration, printed beside its result (CLAUDE.md, rigorous vs empirical;
     the XC mesh has no rigorous bound).
   - With a margin of 3: f_tol = 45·ε_H, so 4.5e-4 Ha/bohr at ε = 1e-5, the conventional default threshold.
   - A tight 1.5e-5 Ha/bohr needs ε_H ≈ 1e-6.
   - η is re-measured on the analytic forces of NH₃, CH₄ and benzene before it is trusted.
4. **SCF tolerance:** the analytic force is first order in the SCF residual (Hellmann–Feynman holds only at the stationary point).
   The gate is on `MinFD`, and its value comes from the same budget: ‖[F, D]‖ well below f_tol.

## 9. Toward the solids (named so the molecular design does not block them)

- **TiO₂, level 1 (imposed):** a, c and u through the solid coordinate map.  The cell part needs stress (`CellChange` today).
- **GPW terms:** each has its own gradient.
  - the collocation of ρ on the uniform raster, and the Hartree potential's gather;
  - V_loc long;
  - the Ewald ion-ion;
  - the periodic Becke partition (an ε-converged image series: its derivative is a series too).
- **Uniform rasters bring the egg-box:** an error that oscillates with the raster spacing and does not move with the atoms.  §8's
  "bias, not noise" fails there; measure it on TiO₂ before a force is trusted on a solid.

## 10. Steps and estimate

| step | content | needs | estimate |
|---|---|---|---|
| 1 | this note (user review); the hub question of §11 | — | done |
| 2 | the optimizer + coordinate maps on model surfaces | nothing | 1 |
| 3 | the derivative-integral evaluator (libcint + MnD l±1), W, the FD-at-fixed-D helper, the Becke weight derivative in qcMesh | nothing (no Hamiltonian edit) | 3.5 (2.5 shared + 1.0 weight derivative tests) |
| 4 | the term gradients of the default stack: Pulay, kinetic, ion-ion, PP local, PP non-local, Dunlap, XC (a)+(b) | lane S's step 4; §11's hub ruling | 6 |
| 5 | gates against PySCF; η re-measured; first optimizations (H₂O, NH₃, CH₄, benzene) | 4 | 1.5 |
| 6 | the molecular audit + two-step pick + cross-setup carry; molecules in the deck + `optimize` | V-FACADE-ACCURACY (coordinate with lane S) | 3 |

Total about 15 lane-days.  Steps 2-3 can start now.

## 11. Questions for the user

1. **The hub query.**  The force assembly needs each term's gradient face.  Should it come from one new query on the tier-1
   `Hamiltonian` hub (each term, or a `GeometryGradient` cross-cast per term)?  Or should it stay off the hub, with the session
   holding the terms' gradient faces when the factory builds them (the factory-concentrates-concrete-wiring rule)?  Lane B leans to
   the factory route: no hub change.
2. **Routes without a gradient.**  Is `NoGradient` the right first answer for every non-default route (δ and multipole Hartree, A1
   and multipole XC, the PP quadrature), each with its own row?
3. **The Dunlap oracle.**  PySCF DF-J has no charge constraint.  Is a near oracle with its gap measured once acceptable, beside the
   exact per-term FD?

## 12. Rulings (user 2026-11-07)

Verbatim; its section numbers are the first version's (§1's FD layer is now §6's test helper, §5's measurement is now §8).


1. **Cartesian**, behind a coordinate map, so the optimizer sees only q and dE/dq; molecules map Cartesian, solids map cell + Wyckoff parameters; redundant internals could be a second molecular map later.
2. **The force threshold is derived from ε** through a NAMED calibration (the 15× is one system): *"If we learn that we need a fudge factor that is fine."*
3. **No finite-difference force product.**  FD gradients spoil the superlinear convergence of BFGS / Newton near the minimum (Numerical Recipes); libcint and PySCF are the oracles; FD would delay the real thing: *"at some point soon we need to take on the real version."*  So the next step is the ANALYTIC-forces design note (molecules first, every term with its oracle: libcint derivative integrals, PySCF whole-molecule gradients, or a per-term finite difference at fixed D as a unit test).  The optimizer is built against analytic model surfaces in unit tests; §1's `FiniteDifferenceForces` survives only as a test helper.  §5's H₂O grid-force measurement stays as evidence for ε.

**§11 answered (user 2026-11-07):**
1. **Forces are a whole-Hamiltonian answer, like the energy**: one `Forces(...)` query on the `Hamiltonian` hub, which sums its own terms' gradients.  No production client needs per-term forces (the optimizer consumes the total); the per-term breakdown is written by each term into the run report, and the per-term oracles are tests (allowed to reach internals).  A hub change, so its exact surface comes to the user before code.
2. **`NoGradient` for every non-default route** (δ / multipole Hartree, A1 / multipole XC, the PP quadrature): refused with an Outcome naming the term and the route BEFORE an optimization starts, each route with its own row.
3. **The PySCF DF-J near oracle for Dunlap** is accepted beside the exact per-term finite difference at fixed D.

**Automatic differentiation spike (user 2026-11-07, "a very sound use of two days"):** before any term gradient is hand-derived, 1-2 days with Enzyme (an LLVM plugin, built from source against our LLVM 21) on the Becke-mesh δ XC energy at fixed D, H₂O: against the per-term finite difference, Σ_A F_A = 0, cost vs one energy, tape memory, friction with modules / Blaze / OpenMP; then dE/dD by AD against the term's H (a pin 32 checker).  Never differentiate through the SCF iterations: forces are dE/dR at the converged D plus the Pulay term.  The result decides the architecture: libcint for integral derivatives, AD for grid and fit terms (and Hessian-vector products for level 2), or hand derivation.

**AD SPIKE VERDICT (lane B, 2026-11-08; evidence, numbers and friction: `doc/Records/EnzymeSpike.md`).  AD works wherever Enzyme sees the body: a leaf differentiated in its own defining TU, or inline / template code it imports.**
- **Exact.** The δ XC kernel (`qchem.AD.XCKernel`, a pure E_xc(R, D) on plain arrays) equals the library's term to 8e-13.
  - Its AD dE/dR equals the per-term FD at fixed D to the FD's own error (≤ 2.3e-9), and Σ_A F_A ≈ 1e-16.
  - Its AD dE/dD equals the term's H to 2e-12 (pin 32 holds).
  - Its Becke-weight gradient equals the hand-derived `BeckePartitionGradient` to 1e-13.
  - All are deterministic, and threaded: OpenMP over independently differentiated serial blocks.
- **Cost** (untuned, against one energy): dE/dR 17× serial and 21× on 4 threads; dE/dD 6.5× / 7.7×.  The tape is < 1 MB, because it is per block.
  A hand-coded XC gradient is typically 3–4× an energy, so AD is about 5× slower than hand-coded here.  For a force, asked once per geometry,
  that is affordable.
- **The danger is silent.**
  - Enzyme's loop caching returns a WRONG reverse-mode gradient with no diagnostic, at any `-O`, for a loop-local array filled through
    an inner loop with a data-dependent trip count.  Minimised: `doc/Records/enzyme_loopcache_repro.cpp`.  The kernel met it as "wrong at
    `-O3`".
  - At `-O2 -march=native` (256-bit vectors) the spike's blocked kernel segfaulted.
  - So AD code compiles at `-O2 -mprefer-vector-width=128`.  No flag set makes it safe, so **every AD gradient is gated by a
    finite-difference test**, and the Enzyme commit is pinned (`1a87ccc`).
  - Draft upstream issues: record §6 (not filed).
- **ARCHITECTURE (user, 2026-11-08): each Imp TU differentiates ITS OWN leaf functions in place, and the derivatives COMPOSE by the
  chain rule across the OOD faces.**
  - In place means the leaf function and its `__enzyme_autodiff` registration sit side by side, in the TU that defines the leaf.
  - For the XC force:
    - the basis Imp delivers ∂χ/∂R;
    - the Becke Imp delivers ∂w/∂R;
    - the XC term's Imp delivers ∂E/∂ρ_g = w_g v_g and ∂E/∂w_g = ρ_g ε_g;
    - the face-level code multiplies them, F_A = Σ_g ρ_g ε_g ∂w_g/∂R_A + Σ_g w_g v_g ∂ρ_g/∂R_A.  That is composition, not a hand derivation of
      a formula.
  - The spike's `qchem.AD.XCKernel`, a separate library restating the formulas, was the DEMONSTRATION only.  Under this architecture
    `qcAD` keeps the typed wrapper (`ReverseInto` on Enzyme's marker entry, `NoGradient`), and the restated kernel retires once the leaves
    land.
  - What the model COSTS (measured 2026-11-08, scratch probes; record §4 items 1 and 9):
    - **Visibility.** Imported `inline` and template code IS visible to Enzyme, instantiated into the TU's IR.  A non-inline function
      defined in another unit is NOT ("No reverse pass found"), and neither is an in-class member of an exported class: in a module that
      is not implicitly inline.  Registered in its OWN defining unit, either one works, and the importer needs no plugin.  So the rule
      is: register where the body is, or make the leaf `inline` in its interface.
    - **Leaves that exist today.**
      - The Becke weight: `qcMesh::BeckePartitionWeight` is a free function.  Its `BeckeCutoff` is defined in another unit (`Mesh.C`), so
        it moves `inline` into `qchem.Mesh`'s interface, or into the same unit.
      - The LDA functionals: `VWN_Correlation`'s and `SlaterExchange`'s energy functions sit in their own units.
      - libxc is external C, with its own analytic derivatives; not ours to maintain.
    - **Leaves that need extracting.** The basis has none.  The PG_Cart point sweep (`IrrepBasisSet::operator()`) calls `GaussianRF` /
      `PrimGaussian` objects defined in another unit (`Evaluators/PG_Cart_MnD`), through pointers.  It needs a plain shell evaluator
      (centre, exponents, coefficients, Cartesian powers, norms → component values) that both the sweep and its AD call: one
      implementation, and its registration beside it.  The same holds for PG_Spherical.  The libcint delivery has its own derivative
      evaluators (`GTOval_ip_*`).
    - **The composition.** It is new code on §13's term faces, the face-level products above.  Each composition is gated against the
      per-term FD, because a leaf's AD cannot see a mistake in how the leaves are combined.
    - **Flags.** A TU that registers derivatives compiles at `-O2 -mprefer-vector-width=128`, a per-file CMake property: it is
      .pcm-compatible both ways, and `-march` is unchanged.  Its ENERGY speed showed no measurable cost:
      - the XC kernel energy: 0.0290 s at the tree's `-O3 -march=native`, 0.0288 s at the AD flags;
      - a strict-order φᵀDφ contraction: 7.84 s against 7.91 s.  An FP reduction is not autovectorised without fast-math, and Blaze's own
        SIMD is chosen by `-march`, which the vector-width preference does not touch.
      - So per-file flags suffice.  A separate AD-only translation of the same source is the fallback if a hot TU ever measures slower:
        the leaf `inline` in its interface, the energy TU at the tree's flags, an AD TU at the AD flags.
- **DERIVATIVES ARE ENZYME'S (user, 2026-11-08: "either 1) hand code all derivatives (and maintain them!!!) xor 2) Use Enzyme ... not some
  blend").**  A build without Enzyme has energies and no forces; the XC force refuses with `NoGradient`.  That is the intent.
  - **The principle: ONE derivative method for every derivative, Enzyme** (user, 2026-11-08: "Clean (consistent) code dictates that all
    derivatives are done the same way.  BUT then Enzyme becomes non-optional.  So I say that temporarily we leave gradients as they are,
    simply because we still need them even if Enzyme+Forces are turned off").
  - **New force derivatives are Enzyme only.**
  - **The ENERGY path's existing hand derivatives stay, FOR NOW**, because a build without Enzyme still needs them: the basis `Gradient(r)`
    faces (`PrimGaussian`, `GaussianRF`, `IrrepBasisSet`; ∇_r χ) and our functionals' `GetVxc` / `GetFxc` (`SlaterExchange`,
    `VWN_Correlation`; libxc's are libxc's).  They move to Enzyme when Enzyme becomes non-optional: an OpenWork row holds that deferred
    decision.
  - **`BeckePartitionGradient` (2026-11-07) is retired from the production path** (user, 2026-11-08).  It was a force derivative, which no
    energy path uses.  The production `qcMesh::BeckePartitionGradient` is now Enzyme's AD of `BeckePartitionWeight`, registered in place
    in `BeckePartition.C`.  The hand derivation stays ONLY as the test oracle in `src/Mesh/tests/BeckePartitionTests.C`: it gates the AD
    gradient at 1e-13, where FD resolves ~1e-9, and AD's failures are silent.

## 13. The `Forces` surface on the hub (DESIGN, for the user before code; lane B 2026-11-07)

Ruling §11.1 puts forces on the `Hamiltonian`, beside the energy.  The hub already has the shape twice:
- `tHamiltonianAt<T>` (V1.41) is the Hamiltonian AT one density, and the energy is asked of it (`GetTotalEnergy()`).
- `MakeResponseKernel` folds the dynamic terms' `tResponse_HT` capability faces "so the caller never sees the term list", and
  names every term that lacks one.

So forces follow both precedents.

**The hub additions (tier-1, the exact surface):**

```cpp
template <class T> class tHamiltonian            // qchem.Hamiltonian, existing
{
    ...
    //! \brief The terms that cannot give a geometry gradient, each as "<term>: <route>" (e.g. "Vee_DFT: δ fit, radial Poisson"):
    //! EMPTY = Forces is available.  A static property of the term set, so a session asks it BEFORE the first SCF of an optimization
    //! and refuses with PESFailure::NoGradient listing them (ruling §11.2).  Folded over the terms like IsVirialValid.
    virtual std::vector<std::string> MissingGradients() const=0;
};

template <class T> class tHamiltonianAt          // existing: the Hamiltonian at one density
{
    ...
    //! \brief −dE/dR_A for every atom, Ha/bohr, the structure's atom order: Σ over the terms of their own derivative at this density,
    //! plus the orthonormality (Pulay) term +Tr(W ∂S/∂R_A).  \a W is the energy-weighted density Σ f ε c cᵀ, in the same blocks and
    //! spin channels as this value's density.  Each term writes its contribution into the run report as it is summed (a
    //! report::Row "forces": term, max|F_A|, |Σ_A F_A|).  THROWS when MissingGradients() is not empty: a caller that did not ask
    //! first is a defect, not an outcome (the Outcome was the session's, before the optimization).
    std::vector<rvec3_t> Forces(const tDM_CD<T>* W) const;
};
```

**Behind the hub, no tier-1 change:**
- **The term face, a new module** `qchem.Hamiltonian.GeometryGradient`:
  - `template <class T> class tGeometryGradient { virtual void AddForces(std::vector<rvec3_t>& F, const tbs_t<T>* wholeBasis, const tDM_CD<T>* D) const=0; };`
  - Implemented by `Kinetic`, `Ven`, `IonIon`, `Ven_PP_Local` and `Ven_PP_NonLocal` (closed form), `Vee_DFT_Dunlap`, and `Vxc_DFT`
    on the δ basis.
  - `tHamiltonianImp` reaches each by the sanctioned abstract→abstract cross-cast, exactly as `MakeResponseKernel` reaches
    `tResponse_HT`.
  - A term without the face is a `MissingGradients()` entry, and a route-dependent term (Hartree, XC) names its route.
- **The basis face, a new module** `qchem.BasisSet.GeometryDerivatives`, beside `AoShellSource`:
  - the capability a Gaussian block gives: dS, dK, dV per atom;
  - Σ D c ∂(ij|k)/∂R on an orbital block, and Σ c c ∂(k|l)/∂R on a fit block.
  - PG_Cart implements it through the libcint evaluator on its own `PGData` (the integrals committed today, gated against MnD).
  - Terms and the Pulay sum reach it by cross-cast from the orbital block they already hold.
- **W is the facade's.**  `MoleculeCalculation` holds every level (ε, f, orbital), so it builds W as a density-matrix value in the
  same blocks as D.  `WaveFunction` (tier-1) gains nothing.

**The first increment's limits, each refused by name, never silently zero:**
- symmetry-blocked runs (SALC blocks): a displacement breaks the blocks, so the contraction needs the AO-level D.  It is a
  `MissingGradients()` entry until the AO assembly lands; the optimizer's symmetrized map does not need SALC blocks.
- HF exchange; +U.
- every non-default route (§11.2).

**Sprint order (§5, unchanged):** the term faces read each term's current state at the converged density, so they land after lane
S's step 4 (V1.41 option C).  `tHamiltonianAt::Forces` is the natural home once stage 2 moves that state into the value.

**Built and gated already on lane-b (not landed), with no hub change:**
- the optimizer + Cartesian map;
- the Becke weight gradient;
- dS / dK / dV;
- the contracted Dunlap 3c / 2c gradients.

**BUILT (lane B 2026-11-08, 9f4a79423):**
- The hub surface: `tHamiltonian::MissingGradients(wholeBasis, st)` / `MakeForces`, and `tHamiltonianAt::Forces(W, st)`.  Both take the
  Structure, and MissingGradients the basis, for the per-face `Can…` queries.
- The term face `tGeometryGradient_HT` (qchem.Hamiltonian.Internal.GeometryGradient), with Kinetic, Ven and IonIon (closed form; periodic
  refuses).
- The density walk `tDM_CD::DM_Pullback(tStatic_PC, st)`.
- The Pulay sum, Tr(W ∂S), in `Imp/HamiltonianForces.C`.
- Gated on E1 (bare nuclei): the ion-ion force; the fixed-D per-term FD; and the FD of the RE-CONVERGED energy, which checks the Pulay term.
- Still missing: the Dunlap Hartree, the δ XC (after lane S's 3b/3c: V-XC-ONE-TERM), the PP terms, HF exchange, and the PW
  structure-factor nuclear force.  `MissingGradients()` names each.

### 13.1 A value-with-derivatives carrier on the derivative-aware faces (lane B 2026-11-08)

**APPROVED (user, 2026-11-08: "Looks good").**  The approved shape:
- The sparse-in-atom first-order jet is the currency on the derivative-aware faces.
- Each face offers both directions: forward jets, and a pullback of a given adjoint.
- **The name is `PositionJet<T>`** (user: "I like PositionJet<T> now that I know it's a diff geom term").  Its declaration says what it
  is: the first-order jet in the atom positions, in Ehresmann's sense, the same object as a Rall number or a dual number.
- Second order is deferred.
- Its matrix and mesh-array realizations form one V-ONE-OF-N record (CLAUDE.md), if more than one exists.
- `qchem.BasisSet.GeometryDerivatives` sits behind a tier-1 hub (`BasisSet`), so its exported surface goes to the controller before that
  part is committed.

Code starts after the Becke leaf lands.  The proposal as written follows.


The user's input (2026-11-08): *"I used to work with Rall algebra where every function evaluation returns {f, df/dx1, df/dx2, ...}.  With
all that info in one blob, chaining is trivial and it makes no difference whether it is a cross library or inter planetary function
call."*  So the currency on every derivative-aware face is a typed carrier, {value, ∂value/∂R}: forward mode at the interface.  Enzyme fills
it INSIDE each Imp, and chaining across a face is ordinary arithmetic on carriers, so library and vtable boundaries stop mattering.  This is
first-order forward AD by operator overloading at the face level (Rall numbers, dual numbers, jets), with Enzyme doing the leaves.

**The type.** A first-order jet in the atom positions, SPARSE in atoms: a value plus a short list of (atom index, ∂/∂R_atom) entries.
- Its arithmetic writes the chain rule ONCE: sums and products merge the atom lists; `exp`, `sqrt`, a functional's ε(ρ) and so on carry f′.
- A matrix-valued jet ({S, ∂S/∂R_A}) holds a libcint integral's derivative.  Its non-zeros are the pairs on atom A, so it is sparse too.
- The name is open.  Name it for what it IS: `PositionJet<T>`, or `WithGradient<T>`.  One family, with the matrix and mesh-array
  realizations as the V-ONE-OF-N record.

**(1) Cost.**
- A DENSE carrier (3N components on every value) is plain forward mode, ∝ 3N per value: ~10× the energy on H₂O (3N = 9) and ~60× on 20
  atoms.  That is no.
- Our derivatives are SPARSE in atoms:
  - a basis function χ_i(r_g) depends on its own centre, and on the point's owner (the point rides with it): at most 2 atoms;
  - ρ_g = Σ D_ij χ_i χ_j depends on the k_g atoms whose functions reach g (`PGData::Reaches`, the magnitude screen pin 1 allows);
  - the Becke weight w_g depends on every atom: no exact sparsity, only magnitude.
- So the carrier is sparse in atoms, and filling one costs:
  - **χ:** one 3-vector, ∇χ, serves both atoms.  MEASURED: Enzyme vector-forward (width 3) on a d shell (6 components, 3 primitives)
    gives value + 3 derivatives in 83.7 ns against 36.2 ns for the value alone, **2.31×**, exact against FD, with no tape.
  - **ρ:** the value costs O(nbf_g²), dominated by Dχ.  The jet is ∂ρ/∂R = 2 Σ_i (Dχ)_i ∂χ_i, an extra O(3 nbf_g), so ≈ 1.1–1.3× the value.
    The width k_g never multiplies the cost, because each χ_i carries one atom.
  - **w_g:** a scalar of all n atoms, so REVERSE mode fills its jet, at a few × the O(n²) weight (§2 below).  For 20 atoms that is 400
    cutoff evaluations per point, against ~10⁴ for ρ at nbf_g ~ 100: a small share.
- ESTIMATED from these pieces, not measured end to end: an XC force at about 2.5–3.5× one XC energy, on H₂O and on 20 atoms alike.  The
  spike's reverse-mode kernel MEASURED 17–21× on H₂O.  The jet form also keeps no tape: each point's jets are transient.

**(2) Where reverse mode stays right, and both directions on one face.**
- **A scalar of many inputs:**
  - the Becke weight w_g(R_1..R_n);
  - the final E, which is a sum over points of a scalar times jets, so it needs no adjoint pass: it is linear at the top.
- **A large linear operator between ρ and E:** the δ-ρ Hartree, E_H = ½ ρᵀKρ.
  - Pushing ρ's jets through K forward costs one solve per derivative direction.
  - The adjoint is ONE solve, and it already exists: v = Kρ.  Likewise the Dunlap fit, whose adjoint is its coefficients.
- So a face offers both directions:
  - FORWARD: return jets of its outputs;
  - PULLBACK: given ∂E/∂(its outputs), return ∂E/∂R as one jet.
  - The term picks per output.  XC: forward χ jets into ρ jets, w_g's jet by reverse, then a pullback of E through the scalar functional.
    Hartree on the mesh: the potential v_g as ρ_g's adjoint, then ρ_g's jets contracted with it.
  - The carrier holds the forward derivative.  The adjoint is a plain value the term already has (v_xc, v_H, fit coefficients).
- **This also sidesteps Enzyme's silent defect.** Forward mode is EXACT on the loop-cache repro (`doc/Records/enzyme_loopcache_repro.cpp`)
  on both 1a87ccc and 871de19f, at `-O2` and `-O3`; the defect is in reverse mode's caching.  Reverse is left only for small scalar
  leaves (the Becke weight), each FD-gated.

**(3) Enzyme's forward mode as the in-Imp filler.**
- `__enzyme_fwddiff(f, enzyme_width, k, enzyme_dupv, stride, x, dx, ..., enzyme_dupv, stride, out, dout)` gives k tangents in one call, a
  vector-forward width of 3 per atom.  That fills χ's jet in the basis Imp.
- The Becke Imp fills w_g's jet by `__enzyme_autodiff`: one scalar, n inputs.
- The functionals keep their hand `GetVxc` FOR NOW (the ruling above), so ε(ρ)'s jet is ε′ = (v − ε)/ρ.  When Enzyme becomes non-optional,
  forward mode replaces it in place.

**(4) What it replaces in §13.**
- The hub surface is unchanged: `MissingGradients()`, `tHamiltonianAt::Forces(W)`.
- `tGeometryGradient::AddForces` stays as each term's face.  Behind it, the hand composition the AD verdict described (the face-level
  products F_A = Σ ρε ∂w/∂R_A + Σ w v ∂ρ/∂R_A) is no longer written per term.  A term evaluates its OWN energy expression on jets, or calls
  a pullback, and the chain rule lives once, in the jet type.
- The basis face `qchem.BasisSet.GeometryDerivatives` returns jets: χ at points (forward, Enzyme) and matrix jets for the integrals
  (libcint's derivative integrals, as committed).  The Pulay term is Tr(W S) with S a matrix jet.
- Gates:
  - the jet type's own unit tests (its arithmetic against FD);
  - each term's per-term FD gate at fixed D, unchanged;
  - the hand Becke oracle.

**Open for the user before code:**
1. The carrier, sparse in atoms, as the face currency.
2. The two directions per face: forward jets, and a pullback of a given adjoint.
3. The name.
4. Second order (a Hessian-vector product for level 2) is a second-order jet or a forward-over-pullback.  Not proposed now.

### 13.2 The derivatives live ON the integral faces (user ruling 2026-11-08; supersedes §13's separate basis face)

The user, on a separate geometry-derivative face and on parallel `_dR` twins: *"why not just add the derivative version right into
Integrals_Overlap<T>, they (mostly) go together"* and *"yes a large number of these classes are all going to be augmented with differential
versions.  I accept that in advance."*  That pre-approval is in CLAUDE.md's hub bullet (main 75b87fd28): a derivative method beside a hub
face's existing method needs no design discussion, and the HubInterfaceAudit snapshot commit cites the ruling.  So:
- **Each integral face gains its derivative beside its method,** as a pullback: given the adjoint A = ∂E/∂M, it returns
  Σ_ij Re(A_ij* ∂M_ij/∂R_atom) per atom, in the structure's site order, as `Outcome<rvec3vec_t, AD::NoGradient>`.
  - `Integrals_Overlap<T>::MakeOverlapPullback(W, st)` beside `MakeOverlap`: the Pulay term.
  - `Integrals_Kinetic<T>::MakeKineticPullback(D, st)`.
  - `Integrals_Nuclear<T>::MakeNuclearPullback(D, st)`: the basis centres AND the nuclei moving.
  - `Evaluatable_IBS<T>::JetsOn(pts, owner, st, PointJets&)` beside `operator()(r)`.  It is FORWARD, on the reviewed slab shape, and its
    gradient source is still with the user.
  - The fit faces gain theirs when Dunlap needs them: `Orbital_DFT_IBS`'s 3-centre repulsion, and the fit metric (k|l).  ∫f_k does not
    move, so the charges need none.
- **Every family answers each method:**
  - Gaussian: libcint's derivative integrals (gated in `M_LibCint.C`).
  - Plane waves: a real implementation.  The functions do not move, so Overlap, Kinetic and the jets are EXACT ZEROS.  The operators that
    move with the atoms (Nuclear, the PP) are the structure-factor Hellmann–Feynman force.  The user notes that the zero is where the
    one-of-N idea helps later, not now.
  - A family with no derivative yet: a one-line `NoGradient` refusal.  The methods are pure virtual, so the compiler makes every family
    decide (compile-time over run-time).
- **No factory cross-casts.** A term asks the face it already holds.
- **The capability query is per FACE, named for what its client does with the answer** (controller, 2026-11-08: one block-wide bool is
  too coarse).  The faces' queries are `CanPullbackOverlap()`, `CanPullbackKinetic()`, `CanPullbackNuclear()` and `CanGiveJets()`, each
  false exactly when its method would refuse.  `MissingGradients()` is each term asking its own faces; it never calls a pullback to find
  out, and needs no shared capability base.
- **Built in increment 1:**
  - PG_Cart: libcint pullbacks, through the Gaussian point mixin, chosen at compile time by its evaluator type.
  - Plane waves: exact zeros for Overlap and Kinetic.
  - Everything else refuses in one line: the other Gaussian evaluators, Bloch Gaussians, PW nuclear, LAPW, radial, Dirac, SALC and the fit
    blocks.
  - `JetsOn`: real on PG_Cart; every other family refuses.
- **`JetsOn`'s gradient source: the existing hand `Gradient(r)`** (user ruling 2026-11-08, "Use Gradient(r)").  ∂φ/∂R_centre = −∇_r φ,
  so there is no second (Enzyme) derivative of φ, the jets work without Enzyme, and a Gaussian Cartesian block never refuses.  When Enzyme
  becomes non-optional (the OpenWork row), `Gradient(r)` moves to Enzyme and `JetsOn` follows for free.  The plain shell-evaluator
  extraction (the Enzyme χ leaf of the AD verdict) is NOT needed for this increment.
- **Every derivative method has its own unit gate, per family that implements it** (the user, 2026-11-08: "we will need to make sure
  these are all correct").  Each gate is a central FD of the SAME quantity at displaced geometry, with the step and the truncation-error
  estimate written at the test, plus translation invariance where it must hold: `src/BasisSet/Gaussian/tests/M_GeometryDerivatives.C`
  and `src/BasisSet/Lattice/tests/PlaneWaveUT.C`.
  - **The first gate on the hand `Gradient(r)` found it WRONG.** `PrimGaussian::Gradient` had its sign flipped: the s functions read
    exactly −FD, so every PG_Cart and PG_Spherical ∇φ was wrong.  The fix was one line.  Its consumers then were ∇ρ / ∇φ outputs that no
    LDA energy reads, which is how it stayed hidden.  This is the first entry in the "Enzyme non-optional" row's tally.
- §13's separate `qchem.BasisSet.GeometryDerivatives` module, and the `_dR` twins, are WITHDRAWN.  `PositionJet<T>` (committed, qcMath)
  and `PointJets` (now a value type beside `Evaluatable_IBS`) stand.

## 14. ε on the molecular facade: step 6's design (DESIGN, for agreement before code; lane B 2026-11-07)

Pin 31 rule 1 ("ε is structure neutral") on `MoleculeCalculation`.  Today `CalcOptions::accuracy` is consumed by no term (the run
says "none yet").  This section says what ε sizes on a molecule, how the two-step run of §2 works there, and what the audit is.
No Hamiltonian hub change: the audit reads `MakeTermEnergies` and `MakeTermMatrices`, which the solid audit already uses.

### 14.1 What ε sizes, and with what measure

The measure is `Consumes()`.  A run built on its own answers `Energy`, so it is judged on |ΔE| ≤ ε·N_atoms.  A `MoleculeSession`
point (a surface, forces, an optimization) answers `MatrixElements`, so it is judged on max|ΔH_ij| ≤ ε, in the normalised basis.  This
is the solid's rule.  `MoleculeCalculation` gains the override; today it inherits `Energy` from the base.

| what | today | under ε | class after |
|---|---|---|---|
| **Becke XC mesh (nR, L)** | stated default (MHL m3 α2, 80, Lebedev 29, order 2); F₂ GTH misses ε by 2.7× on E | climbed by `qchem.Mesh.BeckeClimb` at a fixed D, on the run's measure | ε-AUDITED + raised |
| **δ-ρ Hartree rung** (`CoulombFit::Delta`, non-default) | sits on the XC mesh (or `coulomb.nRadial`) at `ResolvedLMax` 16 | rides the SAME climb: its distance is the max over every mesh-borne term, so the mesh meets ε for all of them | ε-AUDITED (the lMax cut stays FIXED, flagged) |
| **Coulomb fit set** (Dunlap, the default) | the orbital basis's own fit set, unaudited | audited at fixed D against EXACT J (ERI4, four-centre: no fit at all), on the run's measure | ε-AUDITED, WARN on a miss (see Q2) |
| screens (`screenEps` 1e-10, PGData, Hermite, ERI4 floors) | FIXED | not in this item; they stay rows of `doc/Records/EpsilonInventory.md` §3 | FIXED |
| SCF gate | stated per deck | geometry runs take `MinFD` from f_tol (§8 item 4); energy runs stay the user's open call (inventory §4 ⚑) | — |

The audits compare term by term, so a term that does not change between the two meshes reads 0, with no special case.  That covers
the analytic terms: kinetic, Ven, the GTH local and non-local terms, Dunlap J, and ion-ion.

**Where the climb starts: today's default.**  `CalcOptions::mesh` becomes unstated by default (the solid's `xcMesh Auto`).  An
unstated mesh enters the ladder at today's (80, 29) and only ever goes UP.  So an anchor moves only where today's mesh misses ε,
toward the converged value.  Each move is re-judged against the (150, 47) reference (pin 10), and those are the anchors pin 31 calls
wrong.  Entering lower, which saves points, is a separate lever with its own measurement: the radial-rule row (m = 2, the outer-shell
drop) and its pending ruling.

**A STATED mesh** is audited and reported against ε.  If it is OVER, the report says by how much and it is never raised (pin 31
rule 3; the solid's stated-Becke rule).

### 14.2 The two-step run (§2 applied)

1. **Converge at the entry rung.**  An energy run uses the run's own ε.  A geometry point uses ε = 1e-4 for this step only (§2 step 1).
2. **Climb at that D.**  Call `ClimbBeckeLadder(entry, target, cal, false, distance, 0)`.  `distance(a, b)` builds the mesh-borne terms
   on meshes a and b over the SAME basis and the SAME D (a molecule's basis does not depend on its mesh, so nothing is restored),
   refreshes them, and returns |ΔE| or max|ΔH| per `Consumes()`.  It caches by rung, because each step's reference is the next
   step's rung.
3. **If the climb raised:** rebuild the point on the settled mesh and re-converge ONCE from the carried D (`MakeCarriedDensity`, this
   branch).  Then confirm with one more audit (`maxRaises 1`).  This is the solid pilot's shape, `kPilotRounds` included.
4. **Freeze.**  `SizedSetup()` returns the setup with the settled mesh stated in it and marked resolved.  A `MoleculeSession` steps
   every later point on it, so the mesh parameters never change mid-optimization (§2 rule 4).  At the final geometry the run audits
   once more and REPORTS whether the pick differs.  Acting on that is the next run's business (§2 rule 5).
5. **Cost:** report it, never hide it (`PilotCost()`: audit rounds and re-convergence iterations).  On water at (80, 29) one climb
   step is two XC builds at about 1.3-1.7× the run's points, below one SCF iteration each.  No cost policy fires.

**The Coulomb fit audit** runs once, after step 1 (it has no ladder to climb).  It measures |ΔE_J| and max|ΔH_J| between Dunlap J
and exact J at the same D.  Exact J is O(N⁴) through ERI4.  Above a size that setup can compute from the basis count, the audit is
SKIPPED as a WARN-AND-RUN policy row that gives its predicted cost and the reason.

**Every decision is a `policy` row** (CLAUDE.md, every policy reports itself), in the scope where it fires:
- the mesh's origin (stated / climbed from / settled at) and each audit's delivered value against its target;
- the fit audit's verdict, or why it was skipped;
- the `[accuracy] consumed by:` banner, which changes from "none yet" to the terms that consume ε.

### 14.3 The audit's code shape

- **`MoleculeCalculation::MeshAudit(const MeshParams& reference)`** and **`CoulombFitAudit()`**: both are fixed-D, term-by-term
  comparisons.  They return the solid's `GridAuditResult` (per term ΔE and max|ΔH| with the element named, plus the totals).  They
  report under a `gridAudit` scope opened by their caller.
- **The result type and the normalised-basis fold (`FoldΔH`)** move out of `SolidCalculation` into a facade-neutral
  `qchem.GridAudit` module in qcCalculation.  The molecular side uses it first.  `SolidCalculation` adopts it AFTER lane A's R1.0r
  lands, as one small follow-up commit, so this item does not touch lane A's file.
- **Per species: measure before building (Q1).**  `MakeMolecularMesh` takes ONE `MeshParams` for every atom, as the solid's mesh
  does, and the climb climbs one.  Making (nR, L) per species is a qcStructure interface change, plus a climb with one axis pair per
  species.  The banked TeF₆ frozen-D data (`~/Code/qchem6-runs/RadialRule/`; the numpy mesh reproduces ours to 1e-14) can answer
  what it buys without any library code: climb Te and F separately and compare points against the uniform climb.  The radial-rule
  report already found that α per element buys nothing; (nR, L) per element is the untested half.

**Q1 MEASURED (lane B 2026-11-08, the ruled replay): per-species (nR, L) is NOT built.**  TeF₆ at its frozen D, XC only, today's rule
(MHL m3 α2, order 2), the C++ ladder (radials 40-140, degrees 29-65 with the Lebedev menu's own rounding), ε = 1e-5, against the
(400, 77 GL) reference.  Each species' points form their own block, so every (Te rung, F rung) pair is a sum; a pick must pass with
every rung above it passing too.  Cheapest mesh meeting BOTH |ΔE| ≤ ε·N and max|ΔH| ≤ ε:

| system | uniform | per species (Te, F) | saving |
|---|---|---|---|
| PP | 168,688 at (80, 29) | the same | 0% |
| AE | 301,212 at (100, 35) | 282,474 at Te (80, 29), F (100, 35) | 6.2% |

On E alone the PP pair saves 21.5% (F (60, 29)), but H decides the mesh, as on LiF.  6% is below the rung step, so the climb stays one
(nR, L) for every atom, and `MakeMolecularMesh` keeps its one `MeshParams`.  Data: `~/Code/qchem6-runs/RadialRule/species_replay.{py,log}`.

### 14.4 Files, for sequencing with lane G

| file | change | overlap |
|---|---|---|
| `src/Calculation/MoleculeSetup.C` | `mesh` unstated by default (an Auto marker, as `SolidCalcOptions::xcMesh`); `MoleculeSetup` gains the resolved mesh + `meshResolved` | **lane G** (policy rows) |
| `src/Calculation/Imp/MoleculeSetup.C` | resolve the entry rung; the stated/unstated policy row | **lane G** |
| `src/Calculation/MoleculeCalculation.C`, `Imp/MoleculeCalculation.C` | `Consumes()`, `MeshAudit`, `CoulombFitAudit`, the pilot, `SizedSetup`, `PilotCost`, the banner | **lane G** if the "molecular bridge" includes `Outcome_` (survey row 25) |
| `src/Calculation/PotentialEnergySurface.C` / `MoleculeSession` (in `MoleculeCalculation.C`) | answers `MatrixElements`; pilots at the first point, steps on `SizedSetup` | as above |
| new `src/Calculation/GridAudit.C` (+ CMake) | the facade-neutral audit result and fold | none |
| `src/Calculation/tests/` (UTCalculation) + `IntegrationTests/` | water/F₂: the climb raises F₂ and leaves water; a stated mesh is reported OVER and not raised; the fit audit against ERI4 | none |
| **not touched** | the Hamiltonian hub, `SolidCalculation` (R1.0r), `Deck.C` (molecules in the deck stay step 6's second half), `pybind/`, `XCPolicy.C` (no new calibration: the entry rung is today's default) | — |

**Order and estimate:** GridAudit lift + `MeshAudit` + `CoulombFitAudit` (1), then the pilot and `SizedSetup` + `Consumes` (1), then
the session (0.5), then the per-species replay (0.5, measurement only).  About 3 lane-days.  The first piece touches only the
MoleculeCalculation files and the new module, so lane G's `MoleculeSetup.C` edits can land first.

### 14.5 Questions

1. **Per species:** run the TeF₆ replay first, and build per-species (nR, L) only if it saves more than the rung step (≈ 20% of the
   points)?
2. **The Coulomb fit set, MODEL or NUMERIC** (inventory §3 ⚑)?  If NUMERIC, a miss must raise it, which needs a ladder of fit sets
   (we have one set per orbital basis).  If MODEL, a miss WARNs (pin 31: model choices are exempt, but the code warns where it sees
   a poor one).  Lane B proposes MODEL with the audit as the WARN, until a fit-set ladder exists.
3. **Energy runs pay one audit by default** (two XC builds at fixed D).  The solid today audits its Becke mesh only on geometry runs
   (inventory §2, known miss 1).  Lane B proposes that both facades audit every run, and that the solid's energy-run gap become its
   own row.

**§13 and §14 APPROVED (user 2026-11-07, "All approved"):** the §13 `Forces` surface on the `Hamiltonian` hub (`tHamiltonian::MissingGradients()`, `tHamiltonianAt::Forces(W)`; term faces, the basis `GeometryDerivatives` face and W behind the hub); the §14 molecular-ε design with the controller's answers: (1) replay the banked TeF₆ fixed-D data before building per-species (nR, L); (2) the Coulomb fit set is NUMERICAL (pin 31 lists fit-basis incompleteness as non-model), so a miss WARNs with a policy row until a fit-set ladder exists, which gets its own row; (3) the sizing climb runs by default on energy runs on both facades, the separate post-SCF confirmation audit stays opt-in.  The Enzyme build is approved too.

