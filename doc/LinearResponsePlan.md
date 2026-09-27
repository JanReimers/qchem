# Linear response (DFPT / CPHF) — scoping and interface design for A7

**Born 2026-09-27**, the dedicated A7 design session asked for in `doc/HubbardUPlan.md` §0 (row A7 and "A7
scoping insights").  **Status: a PROPOSAL.**  §6 lists the rulings it needs from the user.  No code has been
written against it.  It stays at the top of `doc/` while A7 is executed, and moves to `doc/Records/` when A7
lands.

★ **The constraint that sets the design** (user, 2026-09-26, `doc/HubbardUPlan.md` after the A7 insights):
*"I am mostly concerned about changes to abstract interfaces … If we are adding behaviour to the Hamiltonian
interfaces they should be general perturbation theory (PT) interfaces, not DFT specific.  So if possible the
same PT interfaces would work for MP2 corrections of an HF Hamiltonian."*  Section 3 is the list of interface
changes.  Section 4 checks that list against HF/CPHF and MP2.

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
| **H2** | **NEW term capability `tResponse_HT<T>`**, a cross-cast face on dynamic terms like the RealBlock faces: `RefreshForResponse(const TransitionDensity&)` + `ResponseMatrix(bra k+q, ket k, Spin)`.  The kernel folds it over the terms.  **Static terms are never asked** (they are D-independent by definition).  A dynamic term without the face **fails loudly** (the RealBlock idiom), and is never silently zero. | This mirrors the existing `RefreshForDensity`/`GetMatrix` pair (pin 11, explicit phase).  The δρ-derived state (δV_H, \f$f_{xc}\delta\rho\f$ on the mesh) is k-independent, exactly like its ground-state twin. |
| H3 | `ExFunctional` gains the **second derivative** \f$f_{xc}^{\sigma\sigma'}(\rho_\uparrow,\rho_\downarrow)\f$ (libxc provides `fxc`; VWN/Slater get it analytically) | Needed only by the analytic XC kernel.  Spin-native from the first line (pin 5). |
| H4 | **Hubbard's frozen mode goes onto `HubbardProjection`** (today `FreezeOccupations` exists only on the concrete `Hubbard_U`) | Freezing is how Timrov eq 20 is expressed, and **the kernel of a frozen term is zero by its own state**.  So H1 needs no "exclude +U" flag, and no DFT-specific parameter reaches the generic face. |
| H5 | **`HubbardUEstimator` ISP split**: the strategy face keeps `Evaluate`/`Write`.  `Apply` (writing U into the term) moves to a separate `HubbardUTarget` face, which is really `HubbardProjection::SetU`.  `HubbardEstimate` keeps (site, l, U_eff) plus an optional inter-site row, and the ACBN0-only fields (`UbarBare`, `Nup`, …) go into that estimator's own `Write` (pin 17). | Today's face is ACBN0-shaped.  Its input is a feed of orbitals, and ACBN0 is a functional of orbitals, while LR needs the Hamiltonian.  Each strategy's INPUT becomes its own construction business; the facade chooses the strategy (DIP).  Estimating a U and applying a U are two responsibilities. |

### Below the Hamiltonian
| # | owner | change | why |
|---|---|---|---|
| C1 | qcChargeDensity | **NEW `TransitionDensity<T>`**: per (k, σ) a `mat_t<T>` for block pair (k+q, k), plus q.  It supplies δρ where the terms need it: on the XC mesh, and at \f$\mathbf G+\mathbf q\f$ for Hartree.  **It is NOT a `tChargeDensity`.** | **LSP**: it is complex, has zero trace, is not normalisable and is not mixable along an SCF trajectory.  If it derived from `tChargeDensity`, `GetMatrix(…, δD)` would compile and return nonsense (v_xc of a transition density).  As a separate type that mistake is a build error (compile-time over runtime).  At q=0 it may *contain* a `tDM_CD` as its implementation (composition, reusing the collocation), but it does not *become* one. |
| C2 | qcHamiltonian (Types) | **NEW `TransitionFock<T>`**: the same container shape, a distinct type | Pin 20, ask what a matrix MEANS: δD is contravariant and δF is covariant, and their contraction is an energy.  One container type for both would let the solver add a density to a potential. |
| E1 | qcElectronConfiguration | **`OccupancyRule::ResponseWeight(ε_n,f_n,ε_m,f_m)`**, plus the q=0 Fermi-level shift δμ (keeps Tr δD = 0) | Insight 2 made concrete.  The rule that decided \f$f(\varepsilon)\f$ is the only thing that knows \f$f'(\varepsilon)\f$.  `IntegerOccupancy` returns \f$(f_n-f_m)/(\varepsilon_n-\varepsilon_m)\f$ and **throws** when a partially-occupied pair is degenerate (that is hp.x's "should NOT be treated as a metal" condition, and here it becomes a diagnosis).  `FermiOccupancy` takes the \f$-f(1-f)/kT\f$ limit.  The response reads the SAME composed rule the ground state used, so there is no second metal-detection path (pin 15). |
| S1 | qcSymmetry (Lattice_3D) | **k-mesh query `Shifted(k-block, q) → k'-block`** | Insight 3.  **In our Bloch gauge it is a pure index map**: the phase is \f$e^{i\mathbf k\cdot\mathbf R_n}\f$ over integer lattice offsets (`GPW_Evaluator.C`, `LatticeSum1E.C`), so \f$\phi_{k+G}\equiv\phi_k\f$ and no G-phase bookkeeping exists (QE needs its `ikqs` map because its PW basis is \f$e^{i(\mathbf k+\mathbf G)\cdot\mathbf r}\f$).  Requires a **full, unreduced mesh** at first (§6 D5). |
| M1 | qcMath / qcLASolver | **`LinearOperator<T>` face + a matrix-free Krylov solver** (GMRES; MINRES/CG where the problem is Hermitian) | The response equation is LINEAR, and Krylov is the optimal method for a linear problem (Pulay applied to a linear map is GMRES in disguise).  **The density mixers are NOT reused**: they are `GField`-typed, SCF-trajectory-shaped and live behind `.Internal.`.  This keeps the numerics seam free of any electronic-structure vocabulary. |

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
(`ApplyR0`, `ToMO`) and exposes no `GetC()`/`GetEigenvalues()`**.

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

## 5. Stages: each with its own oracle, cheapest first (insight 6)

★ **Design for general q; q=0 is the efficiency special case.**  This matches the project's Pol/UnPol stance.
`TransitionDensity` carries q from day one, and every face takes a block PAIR.  q=0 is the case where the
pair is (k, k).

| stage | delivers | interface rows | oracle (a wrong number is a bug in NEW code) |
|---|---|---|---|
| **R0** | \f$\chi_0(\mathbf q)\f$ for Hubbard channels by sum over states, primitive cell, full mesh, **no kernel** | E1, S1, C2 (as the probe RHS), qcResponse skeleton (`Reference`, `HubbardChannelProbe`) | hp.x's printed χ₀ on the A6 matched-PP decks: SrVO₃ χ₀(V,V) = −1.7822 (metal: exercises Fermi ResponseWeight + δμ); NiO χ₀ = −0.113 at U_in = 3 eV (insulator: Integer).  Same PP, projector, k-mesh and q-mesh (`IntegrationTests/QE/README.md`).  **No Hamiltonian change at all.** |
| **R1** | CPHF/CPKS, molecular, finite field replaced by response | H1, H2 (Hartree/J, K), C1 at q=0, M1, `DipoleProbe`; **H1 first backed by a finite-difference kernel** \f$[F(D_0+h\delta D)-F(D_0-h\delta D)]/2h\f$ built from the PUBLIC `GetMatrix`, which needs no term code | PySCF (`~/Code/pyscf-env`) static polarisability of H₂O at HF and LDA, and our own finite-field SCF.  **The FD kernel then stays permanently as the unit-test oracle for every analytic `tResponse_HT`**, the same pattern as `Hamiltonian/tests/GPW_XC_FD.C`. |
| **R2** | periodic q=0 self-consistent χ: analytic Hartree + LDA \f$f_{xc}\f$ through GPW | H2 (periodic Hartree/XC), H3, H4 | (a) FD kernel vs analytic on a solid; (b) in a **supercell**, R2 *is* LR-cDFT, checked against a finite-difference cDFT run (perturb with a static \f$\alpha\hat P_J\f$, re-converge; Timrov §III).  **C1's 32-atom MnO measurement says whether (b) is affordable.** |
| **R3** | q ≠ 0 kernel: δρ collocated from (k+q, k) pairs, Hartree at \f$\mathbf G+\mathbf q\f$ | C1 at q≠0 (the collocation pair loop takes a per-image complex weight) | hp.x U: SrVO₃ 6.2502 eV (q 2×2×2), NiO 5.267 eV **at U_in = 3 eV** (frozen +U, H4) |
| **R4** | `LR_HubbardU` behind the split estimator face; `ConvergeHubbardU` drives it; inter-site V_IJ and orbital-resolved channels **reported** | H5 | self-consistent U vs Timrov 2021; the ACBN0 loop still runs unchanged |

**R3 is the cost centre and the only stage whose size is not yet known.**  The ground-state collocation
already applies a per-image Bloch phase.  A transition density needs \f$e^{-i(\mathbf k+\mathbf q)\cdot\mathbf R'}e^{i\mathbf k\cdot\mathbf R}\f$
on the image pair.  Read the collocation pair loop and size R3 before committing to it (ruling D1).

---

## 6. Rulings needed before any code

| | fork | recommendation |
|---|---|---|
| **D1** | **q ≠ 0 (R3) vs supercell-at-q=0 (R2 in a supercell).**  Both give the same χ_IJ(R). | Build R0 at general q regardless: it is cheap and has no kernel.  Take R1–R2 at q=0.  **Decide R3 after C1** (supercell wall/RSS) and after sizing the pair loop.  The interfaces are q-general either way, so neither answer forces an interface change. |
| **D2** | **Ambition: same-site U only, or the full matrix (inter-site V, orbital-resolved channels)?** (insight 5) | **Compute the full \f$\chi_0^{-1}-\chi^{-1}\f$ from day one**: the channel list is the Probe's input (pin 23 generalised one level) and the matrix is inverted anyway.  **Consume only diagonal U** in the +U term.  A +V term is a new *term*, not an interface change, and gets its own OpenWork row.  Orbital-resolved (irrep-slot) channels come for free from the probe, but flag them as numerically delicate (small-χ inversion, the MnO lesson). |
| **D3** | Unknown in the solver: δD (block-pair matrices) or δV (grid potential, as QE mixes `dvscf`)? | **δD.**  It is basis-generic: molecules and solids use the same solver, and it is the representation the Probe measures.  δV-on-grid is a GPW-only optimisation, and can come later behind the same `LinearOperator`. |
| **D4** | Strategy selection and the H5 split | The facade picks `{ACBN0, LinearResponse}`.  `SolidCalculation::EstimateHubbardU`'s orbital-feed loop becomes ACBN0's adapter (it moves into qcResponse, beside `Reference`, which walks the same blocks). |
| **D5** | IBZ-reduced reference meshes | **Full mesh only** for the first version: the k+q block must be *stored*.  Unfolding by rotating orbitals is a later optimisation, and is where C2's shifted-MP fold defect would matter. |
| D6 | Keep the FD kernel after the analytic ones exist? | **Yes, as a permanent test oracle**, never as a production path (at q≠0 it cannot exist: \f$D_0\pm h\,\delta D_{q\ne0}\f$ is not lattice-periodic). |

---

## 7. Traps already paid for, carried forward
- **A response is a function of the state it linearises about** (`HubbardUPlan.md` §1 trap 2).  NiO's
  hp.x 5.267 eV is \f$U_{\rm LR}(U_{\rm in}=3\,{\rm eV})\f$ with V_Hub frozen.  R3's NiO gate must run at the
  same U_in with H4 frozen.
- **d⁵ MnO is not an oracle** (χ₀ −0.045 → χ −0.043, the inverse difference cancels).  Keep it out of R0–R3.
- **The iteration cap has caused 3 wrong conclusions.**  The Krylov solver reports its residual, and a
  non-converged χ is an Outcome, never a number (CLAUDE.md: a call that can legitimately fail returns an
  Outcome).
- **Matched everything, or the comparison is meaningless.**  Use the same UPF, projector flavour, k-mesh and
  q-mesh as the hp.x deck.  If R0 disagrees, suspect the projector first (insight 4) and the ResponseWeight
  second.
