# Linear response (DFPT / CPHF) — scoping and interface design for A7

**Born 2026-09-27**, the dedicated A7 design session asked for in `doc/HubbardUPlan.md` §0 (row A7 and "A7
scoping insights").  **Status: DESIGN RULED 2026-09-27** (user: *"This looks very good"*; D1–D6 agreed, §6),
and **REVISED the same day** for the user's second round of requirements (§4b: PW/Sternheimer extensibility, SO
coupling, forces; §4c: self-consistent U; §6: the Li-supercell question, U/V/J).  **Stage R0 is DONE and
validated (§5b).**  It stays at the top of `doc/` while A7 is executed, and moves to `doc/Records/` when A7 lands.

## ▶ START HERE (next session, written 2026-09-28)

**Where A7 stands.**  R0 — χ₀(q) by sum over states — is built (new library `src/Response/`, qcResponse),
unit-gated (`UTResponse`: a brute-force ring, insulator + Fermi metal), and validated on NiO against hp.x:
χ₀(q)/χ₀(R=0) on Ni1 3d agrees to 0.1–1 % (§5b).  ctest 932/932 at `2306f636`.

**NEXT = stage R1** (user, 2026-09-28: *"continue with the DFPT implementation"* — so R1 goes BEFORE CK-1; CK-1
checkpoint/restart is still owed before the Oct 6–20 unattended window, `OpenWork.md` §2).  R1 = molecular
CPHF/CPKS static polarisability, the first stage with a KERNEL; its oracle is PySCF (`~/Code/pyscf-env`) on H₂O at
HF and LDA.  In order:
1. **The abstract faces C1/C2** (§3): R0 used a concrete MO-basis value type (`Response::BlockPairs`) on
   purpose; R1 is where `TransitionDensity` (qcChargeDensity) and `TransitionFock` (qcHamiltonian Types) become
   abstract faces, because H1's signature must name them.  Decide their homes FIRST — that is the interface
   review the user asked for up front.
2. **H1 `ResponseKernel`** via `tHamiltonian::MakeResponseKernel()`, and H2 `tResponse_HT` on the HF J/K terms
   (they are already linear in D).  The **finite-difference kernel lives in `src/Response/tests/`** (D6, ruled:
   test tree only, friend hooks via `src/forward.H` if ever needed).
3. **M1**: `LinearOperator<T>` + a matrix-free Krylov solver in qcMath/qcLASolver (`Apply` takes a tolerance).
4. **A `DipoleProbe`** (P1's first non-Hubbard concrete) and the solver loop `(1 − R₀K)δD = R₀V`.
5. ⚠ **Generalise the block PAIRING before any molecular run**: R0's `Reference` is Bloch-only — its partner
   map casts to `BlochQN` (`IsShiftOf`) and its D5 check reads the k-star.  A molecule's blocks are point-group
   irreps, and a dipole couples DIFFERENT irreps (A₁→B₁ in C₂v): that is S1's selection rule
   ("Partners(ket irrep, perturbation irrep)"), of which the k-mesh `MeshShift` is one instance.  Start R1 on
   a C₁ (no symmetry) H₂O if the product table is not ready.

**Recipes and traps banked this round** (details §5b, §7): every multi-k run needs `GPW_OMP_THREADS` (the serial
default idled 15 of 16 cores); a MAGNETIC imposition keeps the full k-mesh, so `<P>_IMPOSE=1` satisfies D5; NiO
needs the pin-22 vet trim `NIO_VET=1 NIO_ORTHO_TOL=1e-3` (raw shells Ni s 0.06 + d 0.18) — without it a KB GHOST
state (−36.5 Ha at one k) or a gapless state appears.  The open question "why does a near-null direction host an
attractive KB ghost at all" is parked in `OpenWork.md` §4a (user, 2026-09-28: DFPT first).

★ **The constraints that set the design** (user, 2026-09-26, `doc/HubbardUPlan.md` after the A7 insights):
*"I am mostly concerned about changes to abstract interfaces … If we are adding behaviour to the Hamiltonian
interfaces they should be general perturbation theory (PT) interfaces, not DFT specific.  So if possible the
same PT interfaces would work for MP2 corrections of an HF Hamiltonian."*  And (2026-09-27): *"the design
should also handle perturbations like SO coupling and forces for geometry optimization"*, and a PW code must
be able to join *"by simply extending (new derived implementation with a Sternheimer solver) an existing
interface"*.  §3 is the list of interface changes, §4 checks it against HF/CPHF and MP2, and §4b checks it
against PW, SOC and forces.

Read before this file: Timrov, Marzari, Cococcioni PRB **98**, 085127 (2018) §III–IV
(`~/Code/DFPT3-timrov2018.pdf`); QE `HP/src/hp_main.f90`, `hp_solve_linear_system.f90`,
`LR_Modules/response_kernels.f90`; CP2K `qs_linres_methods.F` (`apply_op` = `apply_op_1` + `apply_op_2`).

---

## 1. The physics in one screen, rewritten for a GAUSSIAN basis

Perturb Hubbard channel \f$J\f$ with \f$\alpha_J\hat P_J\f$ and measure the occupation response
\f$\chi_{IJ}=dn_I/d\alpha_J\f$.  The bare (Kohn–Sham) response is \f$\chi_0\f$, the self-consistent one is
\f$\chi\f$, and \f$U_I=(\chi_0^{-1}-\chi^{-1})_{II}\f$ (Timrov eq 16).  While the system responds, the +U
potential is held at its ground-state value (Timrov eq 20, "keep \f$\hat V_{\rm Hub}\f$ fixed").  A
monochromatic perturbation \f$\alpha\hat P_Je^{i\mathbf q\cdot\mathbf R}\f$ couples Bloch block \f$k\f$ to
block \f$k+q\f$.  The \f$N_q\f$ mesh of such perturbations is exactly equivalent to an \f$N_q\f$ supercell
(`HubbardUPlan.md` §2).

Written as block-pair matrices, with \f$\delta F\f$ the first-order Fock matrix and \f$\delta D\f$ the
first-order density matrix, the equations are:
\f[
\delta D_{k+q,k}=\mathcal R_0\big[\delta F\big]
  =\sum_{nm}C_{m,k+q}\,\frac{f_{nk}-f_{m,k+q}}{\varepsilon_{nk}-\varepsilon_{m,k+q}}\,
   \big(C^\dagger_{m,k+q}\,\delta F_{k+q,k}\,C_{nk}\big)\,C^\dagger_{nk},
\qquad
\delta F=V_{\rm pert}+\mathcal K[\delta D].
\f]
- \f$\mathcal R_0\f$ is the **independent-particle response**.  It depends only on the reference state
  (orbitals, energies and occupations).
- \f$\mathcal K\f$ is the **response kernel**, \f$\partial F/\partial D\f$ at \f$D_0\f$ (Hartree plus
  \f$f_{xc}\f$ for KS, \f$J-K\f$ for HF, zero for any term that is held fixed).  It is a property of the
  Hamiltonian.
- \f$\chi_0\f$ is \f$\mathcal R_0\f$ with no kernel.  \f$\chi\f$ solves the linear fixed point
  \f$(1-\mathcal R_0\mathcal K)\,\delta D=\mathcal R_0V_{\rm pert}\f$.

★★ **The biggest difference from hp.x: we have every virtual orbital.**  A plane-wave code solves a
Sternheimer equation with a \f$\hat P_c=1-\hat O\f$ projector because it cannot afford the empty states
(Timrov eq 27, `sternheimer_kernel`).  A Gaussian basis diagonalises the whole space at every k, so
\f$\mathcal R_0\f$ is an exact dense **sum over states**.  It needs no inner CG, no \f$\alpha\hat O\f$ shift and
no NSCF run at k+q (`hp_run_nscf`), provided the reference already holds the k+q block.  Two things follow:
- \f$\chi_0(\mathbf q)\f$ is cheap at every q.  It is a validation target on its own, with no kernel
  involved (stage R0, §5).
- The cost of the method is the **kernel application**: roughly one Fock build of δρ per response
  iteration.  Insight 6's timings (55 min for SrVO₃, 1 h 48 min for LiCoO₂) are mostly Sternheimer and NSCF
  work, and we do not do that work.

---

## 2. The decomposition: four objects, one per reason to change (SRP)

This refines insight 1.  hp.x's four stages become four objects, each owned by the library that already
knows the thing it needs:

```
   REFERENCE  (qcResponse)            PROBE  (Hubbard: qcHamiltonian face)       KERNEL  (qcHamiltonian)
   orbitals C, ε, f per (k,σ)         channel J -> V_pert (Adjoint)              ∂F/∂D at D0, per term;
   + the run's OccupancyRule          δD -> n_I       (Forward)                  a frozen term answers 0
   + the k -> k+q block map           (the SAME LowdinProjector as E_U)
          │  R0 : δF -> δD                   │                                        │  K : δD -> δF
          └──────────────┬───────────────────┴──────────────────┬─────────────────────┘
                         ▼                                      ▼
              LinearResponseSolver  (qcResponse):  (1 - R0 K) δD = R0 V_pert    (Krylov, qcMath)
                         │   χ0_IJ(q) = n_I(R0 V_J)      χ_IJ(q) = n_I(δD_J)
                         ▼
              ResponseMatrices: Σ_q -> χ_IJ(R), invert, U = diag(χ0⁻¹ − χ⁻¹), V_IJ = off-diagonal
                         ▼
              LR_HubbardU : HubbardUEstimator  (a strategy beside ACBN0)
```

| object | its ONE reason to change | built from (already exists) |
|---|---|---|
| **Reference** \f$\mathcal R_0\f$ | how an unscreened system responds (occupancy statistics, degeneracies, the k+q pairing) | `tWaveFunction` const face, `TOrbitals::GetCoeff`, `OccupancyRule` (insight 2), `BlochQN` (N, ik) |
| **Probe** (perturb / measure) | what we poke and what we read (Hubbard channel, dipole, irrep slot) | `LowdinProjector` Forward/Adjoint pair (insight 4: this is the SAME object that builds E_U) |
| **Kernel** \f$\mathcal K\f$ | how the Hamiltonian's potential depends on the density | the term lists, folded like `IsVirialValid` / `SiteMoments` |
| **Solver** | the numerics of a linear fixed point | a matrix-free Krylov solver over a `LinearOperator` face (§3, M1) |

★ **The perturbation and the measurement are the Adjoint and Forward halves of one pair.**  Our Hubbard term
was built on the MatrixForward/MatrixAdjoint pair (user, 2026-09-11): \f$n=T^\dagger DT\f$ in one direction
and \f$V=TWT^\dagger\f$ in the other.  Perturbing channel J means Adjoint(\f$e_J\f$), and measuring
\f$n_I\f$ means Forward(\f$\delta D\f$).  So \f$\chi_{IJ}=\mathrm{Forward}_I\,\mathcal R\,\mathrm{Adjoint}_J\f$, and
it cannot come out non-symmetric because of a projector mismatch.  **Insight 4, projector consistency, is
enforced by construction and does not depend on discipline.**

---

## 3. The abstract-interface changes, every one of them

Each row gives the owner, the SOLID argument and the stage that first needs it.  **Anything not listed
here does not change.**

### Hamiltonian (qcHamiltonian) — the rows the user asked to see first
| # | change | why this shape |
|---|---|---|
| **H1** | **NEW capability `ResponseKernel<T>`**, handed out by `tHamiltonian::MakeResponseKernel()` (same idiom as `MakeHubbardUEstimator`: capability out, term list hidden).  ONE question: `TransitionFock InducedFock(const TransitionDensity&) const`, i.e. \f$\delta F=(\partial F/\partial D)\,\delta D\f$. | **General PT, not DFT.**  It is the Hamiltonian's linearisation, the object that CPHF, CPKS, TDHF/TDDFT (Casida A±B), MP2's Z-vector and DFPT all use.  It is named for what the client consumes, "what Fock change does this density change induce", not after DFPT or U (CLAUDE.md naming rule). |
| **H2** | **NEW term capability `tResponse_HT<T>`**, a cross-cast face on dynamic terms like the RealBlock faces: `RefreshForResponse(const TransitionDensity&)` + `ResponseMatrix(bra, ket)` over one ALLOWED BLOCK PAIR (S1: (k+q,σ′) ← (k,σ) in general, same-spin k+q for U).  The kernel folds it over the terms.  **A term reads δρ only through the TransitionDensity's SAMPLING faces** (δρ_σσ′ on a mesh, δρ(G+q)) and never through its matrices; that is what lets a PW TransitionDensity arrive without any term changing (§4b).  **Static terms are never asked** (they are D-independent by definition).  A dynamic term without the face **fails loudly** (the RealBlock idiom), and is never silently zero. | This mirrors the existing `RefreshForDensity`/`GetMatrix` pair (pin 11, explicit phase).  The δρ-derived state (δV_H, \f$f_{xc}\delta\rho\f$ on the mesh) is k-independent, exactly like its ground-state twin. |
| H3 | `ExFunctional` gains the **second derivative** \f$f_{xc}^{\sigma\sigma'}(\rho_\uparrow,\rho_\downarrow)\f$ (libxc provides `fxc`; VWN/Slater get it analytically) | Needed only by the analytic XC kernel.  Spin-native from the first line (pin 5). |
| H4 | **Hubbard's frozen mode goes onto `HubbardProjection`** (today `FreezeOccupations` exists only on the concrete `Hubbard_U`) | Freezing is how Timrov eq 20 is expressed, and **the kernel of a frozen term is zero by its own state**.  So H1 needs no "exclude +U" flag, and no DFT-specific parameter reaches the generic face. |
| H5 | **`HubbardUEstimator` ISP split**: the strategy face keeps `Evaluate`/`Write`.  `Apply` (writing U into the term) moves to a separate `HubbardUTarget` face, which is really `HubbardProjection::SetU`.  `HubbardEstimate` keeps (site, l, U_eff) plus an optional inter-site row, and the ACBN0-only fields (`UbarBare`, `Nup`, …) go into that estimator's own `Write` (pin 17). | Today's face is ACBN0-shaped.  Its input is a feed of orbitals, and ACBN0 is a functional of orbitals, while LR needs the Hamiltonian.  Each strategy's INPUT becomes its own construction business; the facade chooses the strategy (DIP).  Estimating a U and applying a U are two responsibilities. |

### Below the Hamiltonian
| # | owner | change | why |
|---|---|---|---|
| C1 | qcChargeDensity | **NEW ABSTRACT face `TransitionDensity<T>`**: the perturbation's symmetry label (P1) plus the SAMPLING faces the terms need: δρ_σσ′ on the XC mesh, and at \f$\mathbf G+\mathbf q\f$ for Hartree.  **The representation is a concrete's business** (revised 2026-09-27 for §4b): `AO_TransitionDensity` holds block-pair matrices \f$\delta D_{k+q,k}\f$ (Gaussian, dense PW); a future `Orbital_TransitionDensity` holds the first-order orbitals \f$\{\psi_v,\delta\psi_v\}\f$ (a large-cutoff PW code, where an \f$N_{\rm pw}^2\f$ matrix cannot exist).  **It is NOT a `tChargeDensity`.** | **LSP**: it is complex, has zero trace, is not normalisable and is not mixable along an SCF trajectory.  If it derived from `tChargeDensity`, `GetMatrix(…, δD)` would compile and return nonsense (v_xc of a transition density).  As a separate type that mistake is a build error (compile-time over runtime).  At q=0 it may *contain* a `tDM_CD` as its implementation (composition, reusing the collocation), but it does not *become* one. |
| C2 | qcHamiltonian (Types) | **NEW ABSTRACT face `TransitionFock<T>`**, a distinct type from C1.  Its faces are the two ways a potential is consumed: `Matrix(bra, ket)` (Gaussian R₀) and `Apply(bra, ket, ψ)` (Sternheimer, \f$\delta\hat V|\psi\rangle\f$) | Pin 20, ask what a matrix MEANS: δD is contravariant and δF is covariant, and their contraction is an energy.  One container type for both would let the solver add a density to a potential. |
| E1 | qcElectronConfiguration | **`OccupancyRule::ResponseWeight(ε_n,f_n,ε_m,f_m)`**, plus the q=0 Fermi-level shift δμ (keeps Tr δD = 0).  **Contract revised 2026-09-27 (review comment 2):** the Integer rule's weight is well defined only if every COUPLED occupied→empty pair has \f$\Delta_{nm}=\varepsilon_m^{\rm empty}-\varepsilon_n^{\rm occ}\f$ resolved above the reference's own eigenvalue noise \f$\delta\varepsilon\f$ (the SCF's last `GetEigenValueChange`).  The Reference computes the **response gap** \f$\Delta_{\min}=\min\Delta_{nm}\f$ over the pairs the perturbation actually couples, and returns an **Outcome**: FAIL if \f$\Delta_{\min}\le0\f$ (INVERTED: the reference is not a gapped aufbau state across that pair) or \f$\Delta_{\min}\lesssim\delta\varepsilon\f$ (UNRESOLVED), naming (k, k+q, σ, n, m).  It **always reports** \f$\Delta_{\min}\f$ and the implied bound \f$\delta\chi_0/\chi_0\lesssim\delta\varepsilon/\Delta_{\min}\f$ beside χ₀, and never silently. | Insight 2 made concrete.  The rule that decided \f$f(\varepsilon)\f$ is the only thing that knows \f$f'(\varepsilon)\f$; the response reads the SAME composed rule the ground state used, so there is no second metal-detection path (pin 15).  ★ **Why exact degeneracy was the wrong test, and why the danger is worse than "a near-touch somewhere": `Crystal_EC`'s INSULATOR mode fills a FIXED count per k-block, with NO aufbau across blocks.**  Within one block, occupied ≤ empty by construction (unless MOM holds it).  Across a (k, k+q) pair nothing orders them, so an occupied band at k+q may sit ABOVE an empty one at k, and the denominator crosses ZERO between mesh points.  So the q≠0 pairs are where it bites and the q=0 pairs are safe, which is exactly why the check runs over COUPLED pairs.  ⚠ **What it must NOT do: flag a small POSITIVE gap as corruption.**  A small-gap insulator's large χ₀ is correct physics.  The only thing that makes it wrong is a gap not resolved above δε, so the threshold is a MEASURED quantity and not a tunable (pin 12).  `FermiOccupancy` needs no gate: its weight is bounded by \f$|f'|\le1/4kT\f$, and it is computed in the cancellation-safe form near \f$\Delta\to0\f$. |
| S1 | qcSymmetry | **Block-pair selection: `Partners(ket irrep, perturbation irrep) → bra irreps`**, i.e. the selection rule bra ∈ Γ_pert ⊗ ket.  **The k-mesh case, `Shifted(k-block, q) → k′-block`, is one instance** (revised 2026-09-27: q IS the perturbation's irrep of the translation group; a spin flip is ΔM_s = ±1; a molecular dipole in C₂v couples A₁ to B₁).  Abelian first, like the molecular symmetry plan.  **Commensurability contract (review comment 3): an off-mesh q is UNREPRESENTABLE, not checked.**  q is never a `rvec3_t`.  It is a `MeshShift` (integer \f$\Delta ik\f$ on the k-mesh's own grid N) that only the k-mesh can construct, so `Shifted` cannot be handed an off-mesh q, and for any MeshShift k+q IS a mesh point (shifted MP meshes included, since q is a difference of mesh points).  A **q-mesh that does not divide the k-mesh** (\f$N_k\bmod N_q\ne0\f$) is a user CONFIGURATION error: the facade returns it as an Outcome at setup, before any SCF is paid for.  A partner block that is not STORED (an IBZ-reduced reference, against D5) is a broken precondition: `Reference` construction throws. | Insight 3, generalised so that SOC, spin-flip and molecular symmetry do not each invent their own pairing.  The contract prefers a type that cannot express the error over a check that catches it (build failure over runtime failure, user 2026-08). |  **In our Bloch gauge it is a pure index map**: the phase is \f$e^{i\mathbf k\cdot\mathbf R_n}\f$ over integer lattice offsets (`GPW_Evaluator.C`, `LatticeSum1E.C`), so \f$\phi_{k+G}\equiv\phi_k\f$ and no G-phase bookkeeping exists (QE needs its `ikqs` map because its PW basis is \f$e^{i(\mathbf k+\mathbf G)\cdot\mathbf r}\f$).  Requires a **full, unreduced mesh** at first (§6 D5). |
| **P1** | qcResponse (face), concretes in their owners | **NEW `Perturbation<T>`** (added 2026-09-27 for SOC and forces, §4b): its symmetry label (the irrep for S1); the one-body operator derivative \f$h^{(1)}\f$ as a `TransitionFock`; and **`MovesBasis()`**, true when the basis depends on the parameter (nuclear displacement; GIAO field).  Only then does the response need the metric \f$S^{(1)}\f$ and the derivative integrals of every term.  The Probe (§2) is a Perturbation paired with its measurement | The Probe only knew "channel J", so a displacement or an SOC operator had no place to sit.  `MovesBasis()` is named for what the client does with the answer (add the \f$S^{(1)}\f$ terms), not after the cause ("is geometric"). |
| M1 | qcMath / qcLASolver | **`LinearOperator<T>` face + a matrix-free Krylov solver** (GMRES; MINRES/CG where the problem is Hermitian).  It works on an ABSTRACT vector (D3).  `Apply` takes a tolerance, which an exact operator ignores and an inexact one (Sternheimer, §4b) honours via FGMRES | The response equation is LINEAR, and Krylov is the optimal method for a linear problem (Pulay applied to a linear map is GMRES in disguise).  **The density mixers are NOT reused**: they are `GField`-typed, SCF-trajectory-shaped and live behind `.Internal.`.  This keeps the numerics seam free of any electronic-structure vocabulary. |

### Listed so nothing is a surprise later, but OUT of A7's scope
| # | owner | change | when |
|---|---|---|---|
| H6 | qcHamiltonian | **term gradient face `tGradient_HT<T>`**: \f$\partial E_{\rm term}/\partial\lambda\f$ at fixed D for a basis-moving λ (nuclear forces), with the \f$-{\rm Tr}[W S^{(1)}]\f$ Pulay term from the energy-weighted density W | **the Forces row** (`OpenWork.md` §2: *"design note first (which terms need a gradient face)"*).  It is a FIRST-order object and never touches the response solver (§4b) |
| H7 | qcHamiltonian | the kernel's **transverse (spin-flip) block**, \f$f_{xc}^{\perp}\f$ | only with a non-collinear ground state (§4b SOC) |

### Deliberately unchanged
- **`tWaveFunction` / `tSCFWaveFunction`.**  The const face is already documented as the face for *"future
  property/post-HF code"*, and the response reads only that.  The occupation policy comes from the facade,
  which owns it (V1.11 moved it off the WF on purpose).
- **`tChargeDensity` and its mixable/DM faces.**  Nothing is added (see C1).
- **`tStatic_HT`, `tDynamic_HF_HT` signatures.**  An HF term's `GetMatrix` is already linear in D, so its
  `tResponse_HT` is its own J/K build applied to δD.
- **`HubbardProjection`'s existing methods.**  The probe consumes `LowdinProjector` exactly as E_U does,
  generalised from a block to a block PAIR: \f$V_{k+q,k}=T_{k+q}\,W\,T_k^\dagger\f$ and
  \f$n=\sum_kT_k^\dagger\,\delta D_{k,k+q}\,T_{k+q}\f$.

### The new library
**qcResponse** sits above qcWaveFunction and qcHamiltonian and below qcCalculation.  It holds `Reference`,
`ResponseProbe` (+ `HubbardChannelProbe`, `DipoleProbe`), `LinearResponseSolver`, `ResponseMatrices`,
`LR_HubbardU`, and later `MP2`.  Following the IrrepCD preference, **`Reference` answers operations
(`ApplyR0`, `ToMO`) and exposes no `GetC()`/`GetEigenvalues()`**.  `ApplyR0` delegates to an
`IndependentResponse` strategy: `SumOverStates` now, and `Sternheimer` when a large-cutoff PW code arrives
(§4b).  The FD kernel is NOT in this library: it lives in `src/Response/tests/` (D6).

---

## 4. Does it serve HF and MP2?  (the user's test)

| client | Reference / R0 | Kernel K (H1) | Probe | Solver |
|---|---|---|---|---|
| **DFPT Hubbard U** (A7) | yes | Hxc; +U frozen ⇒ 0 | Hubbard channel (Adjoint/Forward) | yes |
| **CPHF static polarisability** (HF) | yes (integer occupancy) | \f$J[\delta D]-K[\delta D]\f$, already linear | dipole | yes |
| **CPKS polarisability** (LDA) | yes | \f$J+f_{xc}\f$ | dipole | yes |
| **MP2 energy** | yes: \f$C,\varepsilon,f\f$ per block (ToMO); pair denominators are a *sibling* of R0 | **not used** | — | — |
| **MP2 relaxed density / gradient** (Z-vector) | yes | **the same \f$J-K\f$ kernel** | RHS = MP2 Lagrangian | **the same solver** |
| TDHF / TDDFT (Casida), later | yes (the A−B diagonal) | the same K | — | eigen- instead of linear solve |

★ **Verdict: H1 is the only thing the Hamiltonian learns, and it is theory-neutral.**  The MP2 energy never
touches the Hamiltonian at all.  It needs the reference and a two-body integral source, and the two-body
source belongs to the basis (`BareCoulombSource`/`ERI4`), not to the Hamiltonian.  MP2 *properties* reach
H1 through the Z-vector, which is the same CPHF solve as a polarisability.  **R1's molecular CPHF gate is
what turns this table from a claim into evidence.**

⚠ **Where the generality stops, on purpose.**  A geometric perturbation (phonons, forces, MP2 gradients)
moves the Gaussians, and that adds an \f$S^{(1)}\f$ metric term to the RHS and to \f$\mathcal R_0\f$ (the
Pulay/orbital-response overlap term of CPHF).  U, fields and the MP2 energy have a **fixed** basis and fixed
projectors (Timrov: *"the localized orbitals are a fixed basis set"*).  The RHS type (`TransitionFock`)
leaves room for an \f$S^{(1)}\f$ companion, but none is built.

---

## 4b. Does it extend to PW, SO coupling and forces?  (user, 2026-09-27)

### PW with a Sternheimer solver: yes, by adding concretes, after three revisions made today
The user's test: a PW code joins *"by simply extending (new derived implementation with a Sternheimer solver)
an existing interface"*.  That holds if these three things are true:
1. **\f$\mathcal R_0\f$ is a STRATEGY behind the Reference.**  Face `IndependentResponse`, with two concretes:
   - `SumOverStates`, for Gaussian and for **our own dense PW**, which already diagonalises fully, so it
     needs nothing new.
   - `Sternheimer`, for a large-cutoff PW code.  It solves
     \f$(\hat H-\varepsilon_v+\alpha\hat O)|\delta\psi_v\rangle=-\hat P_c\,\delta\hat V|\psi_v\rangle\f$ over the
     occupied states only.  For metals it uses the smeared form of Baroni et al. RMP 73, 515 (2001) §II.C,
     and its occupation weights come from E1's rule again.
2. **C1/C2 are abstract, and their representation belongs to the concrete** (the §3 revision).  A term reads
   δρ only through sampling faces, and a potential is consumed as `Matrix` or as `Apply`, so no term knows
   which Reference produced δρ.
3. **The solver works on an abstract vector** (M1; D3 corrected).

H1, H2, P1, S1 and E1 do not change.  QE's own seams map one-to-one: `sternheimer_kernel` → the Sternheimer
concrete; `dv_of_drho` → H1; `mix_potential` → M1; `hp_dnsq` → the Probe's measurement.
⚠ **One consequence to design in:** a Sternheimer \f$\mathcal R_0\f$ is an INEXACT inner solve.  QE tightens it
as the outer loop converges (`thresh = min(0.1·√dr2, 1e-2)`), and an outer Krylov method over an inexact
operator needs a flexible variant (FGMRES).  So `LinearOperator::Apply` should take a tolerance from day
one.  An exact operator ignores it.

### SO coupling: two different things, and only one of them is perturbation theory
- **SOC IN the ground state** (self-consistent, non-collinear spinors; magnetic anisotropy, heavy elements)
  is a Hamiltonian TERM plus a two-component basis and SpinGroup.  That is ground-state work, and it is
  anticipated already: `PreservesReal()` names SOC as the term that answers no.  It adds no PT interface.
  Response on top of it needs spin-mixed block pairs (S1) and the transverse kernel (H7).
- **SOC AS a perturbation** (second-order SOC energies, anisotropy by PT, g-shifts) is a P1 concrete with
  \f$h^{(1)}=\xi(r)\,\mathbf L\cdot\mathbf S\f$.  Its \f$L_\pm S_\mp\f$ parts couple (k,↑) to (k,↓), which is
  exactly S1's selection rule with \f$\Delta M_s=\pm1\f$; \f$L_zS_z\f$ is same-spin.  The second-order energy is
  \f$\langle\psi_0|h^{(1)}|\psi^{(1)}\rangle\f$ from \f$\mathcal R_0\f$ (plus the transverse kernel if the response
  is self-consistent).  **No new interface beyond P1, S1 and H7**, and these are the reasons P1 exists and
  S1 was generalised today.
- ⚠ In our PP world SOC comes from the PSEUDOPOTENTIAL.  HGH's relativistic tables carry spin-orbit
  coefficients (the \f$k_{ij}\f$).  **Checked 2026-09-27: `src/Pseudopotential/Data/gth_potentials.json` does NOT
  carry them** (each channel has `h`, `l` and `r` only).  So SOC in any form starts with a PP-data increment.

### Forces: they do NOT go through the response solver, and that is the physics, not a gap
Forces on a variational SCF energy (HF, KS, DFT+U at fixed U) are FIRST-order.  By Hellmann–Feynman and
Wigner's 2n+1 theorem they need only the ground state:
\f[
F_A=-{\rm Tr}\big[D\,h^{(A)}\big]-(\text{each term's explicit }\partial E/\partial R_A)+{\rm Tr}\big[W\,S^{(A)}\big].
\f]
So forces use **P1** (a displacement, with `MovesBasis()` true) and a per-term **gradient face (H6)**.  They
use neither H1 nor the solver.  H6 is the Forces row's design note (`OpenWork.md` §2), and its hardest term
is ours: **the +U projector moves with its atom, and an ortho-atomic one depends on S of every atom** (QE
needed a dedicated paper for ortho-atomic +U forces, Timrov et al. PRB 102, 235159 (2020)).  The response
machinery enters geometry work in exactly three places, and each is a reuse:
1. **Hessian / phonons** = the derivative of the forces = the solver with a `MovesBasis()` perturbation,
   whose \f$S^{(1)}\f$ is the metric term the RHS type left room for.
2. **Gradients of a NON-variational energy** (MP2, or a U that depends on the density) = the Z-vector: the
   same solver, with a Lagrangian RHS.
3. **Relax with a self-consistent U**: the outer loop of §4c alternates relax and U.  U is constant during
   each relaxation, so the forces are the plain first-order ones.

## 4c. How U (and V, J) become self-consistent, and why it does not touch DIIS, GDM or mixing

**Never per SCF iteration.  Always as an OUTER loop around a converged SCF.**  Timrov et al. PRB 103,
045141 (2021) do exactly this: converge at \f$U_{\rm in}\f$, run LR about that state with \f$V_{\rm Hub}\f$ frozen
(H4), take \f$U_{\rm out}\f$, and repeat until \f$|U_{\rm out}-U_{\rm in}|\f$ is below tolerance, optionally
relaxing the structure inside each cycle.  That is the loop `SolidCalculation::ConvergeHubbardU` already
runs for ACBN0, and `LR_HubbardU` (R4) plugs into the same loop.
- **Inside each SCF, U is a constant.**  The problem is an ordinary variational one, and DIIS, GDM/OT,
  Kerker/Pulay mixing and MOM see nothing unusual.
- **Between outer steps the SCF restarts as a fresh STAGE.**  This was measured the hard way on
  2026-09-21: a continued Iterate let eight near-zero Pulay residuals extrapolate straight back onto the old
  density and report "converged" in one iteration.  The LR solve warm-starts from the previous outer step's
  δD, which is cheap because the state moved only a little.

**Why per-iteration updating is wrong, and not merely expensive.**  ACBN0 as published is a "pseudo-hybrid"
functional: U[ρ] is recomputed from the current density.  Then the Fock matrix built with the current U is no
longer the gradient of the energy unless it also carries \f$\partial E/\partial U\cdot\partial U/\partial D\f$.
- Fixed-point methods (density mixing, Pulay, commutator DIIS) still converge, but to the fixed point of an
  INCONSISTENT functional.
- Energy-based direct minimisers (GDM, OT) break, because their line search needs a gradient that is
  consistent with the energy.
- For LR it is also meaningless: \f$\chi\f$ is the linearisation about a CONVERGED state, and a χ taken from a
  mid-SCF density linearises about nothing.

A density-dependent U also puts \f$\partial U/\partial R\f$ into the forces (§4b), which is item 2 of that
list, the Z-vector, for a quantity nobody wants to differentiate.  ⇒ **The outer loop is the design, for
ACBN0 and LR alike.**

---

## 5. Stages: each with its own oracle, cheapest first (insight 6)

★ **Design for general q; q=0 is the efficiency special case.**  This matches the project's Pol/UnPol stance.
`TransitionDensity` carries q from day one, and every face takes a block PAIR.  q=0 is the case where the
pair is (k, k).

| stage | delivers | interface rows | oracle (a wrong number is a bug in NEW code) |
|---|---|---|---|
| **R0** ✅ machinery (§5b) | \f$\chi_0(\mathbf q)\f$ by sum over states over the DECLARED Hubbard channels only (D2 scope), primitive cell, full mesh, **no kernel**, reported WITH its response gap \f$\Delta_{\min}\f$ and bound (E1) | E1, S1, C2 (as the probe RHS), qcResponse skeleton (`Reference`, `HubbardChannelProbe`) | hp.x's printed χ₀ on the A6 matched-PP decks: SrVO₃ χ₀(V,V) = −1.7822 (metal: exercises Fermi ResponseWeight + δμ); NiO χ₀ = −0.113 at U_in = 3 eV (insulator: Integer).  Same PP, projector, k-mesh and q-mesh (`IntegrationTests/QE/README.md`).  **No Hamiltonian change at all.** |
| **R1** | CPHF/CPKS, molecular, finite field replaced by response | H1, H2 (Hartree/J, K), C1 at q=0, M1, `DipoleProbe`; **H1 first backed by a finite-difference kernel** \f$[F(D_0+h\delta D)-F(D_0-h\delta D)]/2h\f$ built from the PUBLIC `GetMatrix`, which needs no term code | PySCF (`~/Code/pyscf-env`) static polarisability of H₂O at HF and LDA, and our own finite-field SCF.  **The FD kernel then stays permanently as the unit-test oracle for every analytic `tResponse_HT`**, the same pattern as `Hamiltonian/tests/GPW_XC_FD.C`. |
| **R2** | periodic q=0 self-consistent χ: analytic Hartree + LDA \f$f_{xc}\f$ through GPW | H2 (periodic Hartree/XC), H3, H4 | (a) FD kernel vs analytic on a solid; (b) in a **supercell**, R2 *is* LR-cDFT, checked against a finite-difference cDFT run (perturb with a static \f$\alpha\hat P_J\f$, re-converge; Timrov §III).  **C1's 32-atom MnO measurement says whether (b) is affordable.** |
| **R3** | q ≠ 0 kernel: δρ collocated from (k+q, k) pairs, Hartree at \f$\mathbf G+\mathbf q\f$ | C1 at q≠0 (the collocation pair loop takes a per-image complex weight) | hp.x U: SrVO₃ 6.2502 eV (q 2×2×2), NiO 5.267 eV **at U_in = 3 eV** (frozen +U, H4) |
| **R4** | `LR_HubbardU` behind the split estimator face; `ConvergeHubbardU` drives it; inter-site V_IJ and orbital-resolved channels **reported** | H5 | self-consistent U vs Timrov 2021; the ACBN0 loop still runs unchanged.  ★ **The U_0 / U_SC table (user, 2026-09-27: *"the chemist in me wants a feel for how these behave"*) comes from ONE run**: start the loop at U_in = 0 — outer step 1 IS U_0, the last is U_SC.  The loop prints one greppable `[U table]` line per manifold: U_0, U_SC, outer steps, and the GAP and SITE MOMENT at both ends — because U_0 linearises the U=0 state, and a large U_0→U_SC change usually means U changed the state's character (NiO loses AFM-II at U=0, trap 3) |

**R3 is the cost centre and the only stage whose size is not yet known.**  The ground-state collocation
already applies a per-image Bloch phase.  A transition density needs \f$e^{-i(\mathbf k+\mathbf q)\cdot\mathbf R'}e^{i\mathbf k\cdot\mathbf R}\f$
on the image pair.  Read the collocation pair loop and size R3 before committing to it (ruling D1).

---

### 5b. R0 execution record (2026-09-27) — the machinery VALIDATED; the ground state is the open item
**Code:** `d7c95c92` (+ `bd4bba7a`: an unmeasured eigenvalue noise is NaN and says so).  Unit gates: the
brute-force ring (insulator + Fermi metal, every real-space element), the response weights, MeshShift; ctest
927/927.  Logs: `~/Code/qchem6-runs/a7_r0_nio/`.

**The NiO gate** (`gpwprobe nio`, AFM-II, U_in = 3 eV on Ni 3d, `orthofull`, k 2×2×2 full mesh, q 2×2×2):
```
GPW_OMP_THREADS=12 GPW_SPHERICAL=1 NIO_KMESH=2 NIO_IMPOSE=1 NIO_ORTHO_TOL=1e-3 NIO_U=3 NIO_U_RADIAL=orthofull
NIO_SKIP_FM=1 NIO_ANNEAL=5e-3,0 NIO_ACC=Null NIO_MOM=0 NIO_PULAY=8 NIO_PULAY_START=5 NIO_MEASURE=maxdd
NIO_EPS=1e-6 NIO_CHI0=2   gpwprobe nio
```
(`NIO_IMPOSE=1` keeps the FULL mesh for a MAGNETIC imposition — `DetectPointOps`: the Shubnikov k-fold is
Γ-only for now — so it satisfies D5.  Converges in 19 + 10 iterations.)

★ **THE SHAPE MATCHES hp.x TO 1–2 %** — χ₀(q)/χ₀(R=0) on Ni1 3d, independent of the overall magnitude:

| q-star (weight) | ours | hp.x `NiOgO` |
|---|---|---|
| Γ (1) | 0.918 | 0.935 |
| 3-star (3) | 1.048 | 1.041 |
| 3-star (3) | 0.982 | 0.982 |
| (½,½,½)-type (1) | 0.989 | 0.996 |

Star members agree to 1e-5 (the symmetry is respected), χ is Hermitian and the real-space block is real.  This
validates exactly what R0 adds: the k+q pairing, the phase convention and the Fourier sum.

**The magnitude is the GROUND STATE's, not the response's:** on-site χ₀(Ni1) −0.1476 eV⁻¹ vs hp.x −0.1130,
ratio **1.31**; our gap 2.18 eV vs hp.x 2.86 eV, ratio **1.31**.  An independent-particle response scales as
1/gap, so the residual is the gap, not R0.

⛔ **THE GAP IS A BASIS-CONDITIONING QUESTION, filed as `OpenWork.md` §4a "NiO VA: the diffuse Ni s".**  At the
default `orthoTol=1e-4` every run (free Ladder+MOM, free deck-shaped, imposed deck-shaped) converged GAPLESS —
one empty level at Γ below occupied levels at other k — and **E1's gate refused it** (`INVERTED coupled pair`)
instead of printing a χ₀ for a metal: the review round's comment 2, caught on its first real material.  At
`1e-3` the dropped AOs are index 0/47 = each Ni's FIRST function, the α = 0.06 diffuse s; the state becomes a
2.18 eV insulator **0.36 Ha HIGHER** in energy (−106.3729 vs −106.7348).  A function "reproducible by the kept
set" lowering a converged energy by 0.36 Ha smells like the GPW near-dependence dive (the MnO 136-span saga);
GPW collocation is not strictly variational, so that is a measurement to make, not a verdict.  **Follow-up the same day:** the uniform
(QE-like) and a fine Becke quadrature reproduce the gapless state (±24 µHa), and QE's own band list puts the
CBM at Γ with a 3.5 eV direct gap where ours is ~0.13 eV — a state BELOW the PW reference, which the Ritz bound
forbids for an incomplete basis: ill-conditioning × integral error, so the 1e-3 insulator is likely the
physical state (full reasoning and the running eps=1e-14 test: `OpenWork.md` §4a row).  ★★ **RESOLVED 2026-09-28 — R0 VALIDATED ON A CLEAN
STATE.**  The pathology was a KB ghost in a near-null direction; the pin-22 VET-STAGE trim at orthoTol=1e-3
(`NIO_VET=1 NIO_ORTHO_TOL=1e-3`: raw shells Ni s {0.06} + Ni d {0.18} from both sites; min eig S 9.8e-7 → 4.7e-3;
116/116 in every block) gives a clean insulating NiO: E −106.2021, gap 1.30 eV, 15+6 iterations.  χ₀(q)/χ₀(R=0)
on Ni1 3d vs hp.x: Γ **0.936/0.935**, 1.038/1.041, 0.989/0.982, 0.986/0.996; on-site −0.1617 vs −0.1130 eV⁻¹ — the
residual is the ground-state gap (1.30 vs 2.86 eV), a basis/physics comparison.  Log:
`~/Code/qchem6-runs/a7_r0_nio/nio_k222_imposed_VET1e-3_chi0.log`.  **The NiO gate recipe is §5b's line plus
`NIO_VET=1 NIO_ORTHO_TOL=1e-3`.**

### 5a. Timeline, with the infrastructure it leans on (2026-09-27, user: fold in KP and checkpointing)
1. **R0** ✅ — code landed (`d7c95c92`); VALIDATED on NiO 2026-09-28 (§5b).  Lesson already banked: **`GPW_OMP_THREADS` is
   part of every multi-k recipe** (serial default ⇒ ~1.2 of 16 cores; `OpenWork.md` row KP, measured).
2. **CK-1 checkpoint/restart** (`OpenWork.md` §2 row "SCF checkpoint/restart") — ⚠ REORDERED 2026-09-28: AFTER R1
   (user: DFPT first), but still before the Oct 6–20 window: every A6/A7 material's converged state saved once and reused; CK-2 then lets χ₀ run on a
   stored state with no SCF.
3. **R1** (molecular CPHF, cheap, needs neither of the above).
4. **The cross-k gather memo, then KP** (row KP's measured order) — when R2/R3's kernel applications, which
   are gather-shaped per block PAIR, become the wall.  R3's q-points are also embarrassingly parallel at the
   PROCESS level (as hp.x's `start_q/last_q`), which needs no code.
5. **R2 → R3 → R4** on the stored materials.

## 6. Rulings (all six ruled 2026-09-27)

| | fork | recommendation |
|---|---|---|
| **D1** | **q ≠ 0 (R3) vs supercell-at-q=0 (R2 in a supercell).**  Both give the same χ_IJ(R). | ✅ **RULED: agreed.**  Build R0 at general q (cheap, no kernel), take R1–R2 at q=0, decide R3 after C1 and after sizing the pair loop.  The user's prior (earlier analysis): multi-q beats q=0-in-supercells for our design, and the next paragraph agrees.  **The Li-configuration question** (user: *"Li_0.125Mn2O4 in a 2x2x2 super cell: in the q=0 method do we then need a super-super-cell 4x4x4?  Or will the 2x2x2 then serve dual purpose?"*): the perturbation's EFFECTIVE supercell is always **(the cell the ground state runs in) × (the q-mesh)**.  So the configuration cell serves dual purpose only as far as it is already big enough that the perturbed Mn's images have decayed.  If it is not big enough, q=0 needs the super-super-cell (8× the atoms, ~512× the dense-diagonalisation cost), while multi-q stays in the 2×2×2 and adds a q-mesh on its small BZ.  Only the q route scales, which is D1's answer again.  ★ **But the programme may never need an LR run in a configuration cell at all**: `HubbardUPlan.md` §3 computes U at the three end members in their own small cells and ASSIGNS it per Mn site by local oxidation state, so the CE training supercells need no new U.  An LR run on a Li_x cell is the TRANSFERABILITY CHECK (does a charge-ordered LiMn₂O₄ reproduce the end members' per-site U?), and there the 2×2×2-plus-q-mesh answer applies. |
| **D2** | **Ambition: same-site U only, or the full matrix (inter-site V, orbital-resolved channels)?** (insight 5) | ✅ **RULED: agreed.**  **U, V and J are three different quantities, not two naming conventions** (user asked).  **U** = on-site Coulomb, the DIAGONAL of the CHARGE response \f$(\chi_0^{-1}-\chi^{-1})_{II}\f$.  **V** = inter-site Coulomb (I≠J), the OFF-DIAGONAL of the same matrix.  It is an interaction, not a hopping: DFT+U+V's energy \f$-\tfrac12\sum_{I\ne J}V_{IJ}{\rm Tr}[n^{IJ}n^{JI}]\f$ acts on the inter-site occupation \f$n^{IJ}\f$, which measures hybridisation, so V FAVOURS hopping-like mixing (hence the association).  **J** = on-site Hund's exchange, from the SPIN response: perturb with \f$\alpha(\hat P_I^\uparrow-\hat P_I^\downarrow)\f$, measure the site magnetisation, and take J from the same inverse-response difference in the magnetisation channel (Linscott, Cole, Payne, O'Regan, PRB 98, 235157 (2018); read the sign and normalisation from the paper before coding, not from here).  In this design **J is one more Probe** (a spin-antisymmetric channel), not an interface change.  It gives A5's missing J oracle a native route, and it is the only LR J we could compare against ACBN0's \f$\bar J\f$.  Dudarev's \f$U_{\rm eff}=U-J\f$ is where the "U, J" pairing comes from.  ⚠ Notation: V also names potentials here (\f$V_{\rm pert}\f$, \f$V_{\rm Hub}\f$); write \f$V_{IJ}\f$ with its indices always.  **Compute the full \f$\chi_0^{-1}-\chi^{-1}\f$ from day one**: the channel list is the Probe's input (pin 23 generalised one level) and the matrix is inverted anyway.  **Consume only diagonal U** in the +U term.  A +V term is a new *term*, not an interface change, and gets its own OpenWork row.  Orbital-resolved (irrep-slot) channels come for free from the probe, but flag them as numerically delicate (small-χ inversion, the MnO lesson).  ★ **SCOPE, confirmed 2026-09-27 (review question): "full" means full over the DECLARED channel set**, i.e. the run's input manifold list × the q-mesh's cell images.  It does not mean every possible site, irrep slot or spin channel in the crystal.  R0's first PR: the Hubbard manifolds the deck lists (SrVO₃: V 3d; NiO: both Ni 3d sites, as `hp.x` has them), the χ₀ matrix over those channels per q, Fourier to R.  Nothing else.  Irrep-slot channels and J's spin channels are later Probe concretes, each its own increment.  **UNCERTAINTY IS PROPAGATED, NOT ASSUMED (review comment 4).**  Every U, V and J carries an error bar from first-order perturbation of the inverse, \f$\delta(\chi^{-1})=-\chi^{-1}\,\delta\chi\,\chi^{-1}\f$, with δχ taken from the Krylov residual and δχ₀ from E1's \f$\delta\varepsilon/\Delta_{\min}\f$ bound.  It is amplified by \f$1/\chi^2\f$: large exactly where χ is small (the MnO cancellation, χ₀ ≈ χ).  **The gate is on the CONSUMED quantity**: Dudarev's \f$U_{\rm eff}=U-J\f$, with \f$\delta U_{\rm eff}=\sqrt{\delta U^2+\delta J^2}\f$, amplified by \f$U/(U-J)\f$ (about 1.25 for a typical TMO, and unbounded as U → J), never "U converged" and "J converged" separately.  **And `ConvergeHubbardU`'s `tolU_eV` must exceed \f$\delta U_{\rm eff}\f$, or the outer loop chases noise and cannot converge legitimately**: the iteration-cap trap one level up. |
| **D3** | Unknown in the solver: δD (block-pair matrices) or δV (grid potential, as QE mixes `dvscf`)? | ✅ **RULED 2026-09-27: keep everything general** (user).  ⚠ **CORRECTED the same day: the proposal's "δD" was a Gaussian-ism.**  A large-cutoff PW code cannot form an \f$N_{\rm pw}^2\f$ δD at all.  The solver works on an ABSTRACT vector (M1) whose representation the Reference chooses (C1's concretes).  Gaussian: δD block pairs.  PW-Sternheimer: δψ, or δV on the grid. |
| **D4** | Strategy selection and the H5 split | ✅ **RULED: agreed.**  The facade picks `{ACBN0, LinearResponse}`.  `SolidCalculation::EstimateHubbardU`'s orbital-feed loop becomes ACBN0's adapter (it moves into qcResponse, beside `Reference`, which walks the same blocks). |
| **D5** | IBZ-reduced reference meshes | ✅ **RULED: agreed.**  **Full mesh only** for the first version: the k+q block must be *stored*.  Unfolding by rotating orbitals is a later optimisation, and is where C2's shifted-MP fold defect would matter. |
| D6 | Keep the FD kernel after the analytic ones exist? | ✅ **RULED 2026-09-27: yes, as a test oracle, and it LIVES IN THE TEST TREE** (user: *"ideally it gets moved into the unit test or IntegrationTest area.  If it needs hooks into the production code we should make those hooks private with friend access for the testing harness (see src/forward.H)"*).  It starts there on day one, not later: it is built from the PUBLIC `tHamiltonian::GetMatrix` plus a density factory, so it should need no hooks.  Home: `src/Response/tests/` (unit over integration).  R1's production kernel for HF is the analytic one, which is trivial because J and K are linear in D.  At q≠0 an FD kernel cannot exist (\f$D_0\pm h\,\delta D_{q\ne0}\f$ is not lattice-periodic). |

---

## 7. Traps already paid for, carried forward
- **A response is a function of the state it linearises about** (`HubbardUPlan.md` §1 trap 2).  NiO's
  hp.x 5.267 eV is \f$U_{\rm LR}(U_{\rm in}=3\,{\rm eV})\f$ with V_Hub frozen.  R3's NiO gate must run at the
  same U_in with H4 frozen.
- **d⁵ MnO is not an oracle** (χ₀ −0.045 → χ −0.043, the inverse difference cancels).  Keep it out of R0–R3.
- **The iteration cap has caused 3 wrong conclusions.**  The Krylov solver reports its residual, and a
  non-converged χ is an Outcome, never a number (CLAUDE.md: a call that can legitimately fail returns an
  Outcome).
- **A number with no error bar invites being over-read** (review round, 2026-09-27).  χ₀ reports its response
  gap \f$\Delta_{\min}\f$ (E1), χ reports its residual, and U/J/U_eff report the propagated uncertainty (D2).  A
  gate that compares two values without their error bars is a gate on noise.
- **Matched everything, or the comparison is meaningless.**  Use the same UPF, projector flavour, k-mesh and
  q-mesh as the hp.x deck.  If R0 disagrees, suspect the projector first (insight 4) and the ResponseWeight
  second.
