# OpenWork_History5 — verbatim archive of OpenWork v3 post-§1-shrink as of 2026-10-01 (RECORD, no open work)

> ⚠ **rundeck / gpwprobe timings before 2026-10-31 had THREADED BLAS under them** (fixed c95c1fa4): neither main pinned OpenBLAS, so it ran its default 16-thread pool, busy-waiting beside our OpenMP regions, while the banner printed `BLAS=1`.  `gpwprobe` dates from 2026-09-15, `rundeck` from 2026-10-03.  **Re-measure before quoting any such wall/CPU time as a serial comparison.**  Numbers from ctest / ITMain / the UT exes and from scfrun are UNAFFECTED: those mains always pinned.  The pin now runs when qcCommon loads, for every executable, and the banner's `BLAS=N` reports the count OpenBLAS actually runs.

# Open Work — the live tracker (v3, rebuilt 2026-09-16)

**READ THIS AT SESSION START.**  Four sections, in the order a session should read them: **§1 NEXT** (the one
queued action), **§2 MAJOR FEATURES** (capabilities the code does not have yet, against the battery north-star),
**§3 NON-OOD CLEANUP** (tooling, build, hygiene — the SOLID/OOD debt is `doc/CleanupCandidates.md` v2, not here),
**§4 REMAINING TODO** (measurements, performance and accuracy items that are neither).  Then the parked and
descoped lists.  **If it is not in this file or `CleanupCandidates.md`, it does not need doing.**

**Where the history went.**  The v2 tracker (cut 2026-08-19, 1815 lines) is `doc/Records/OpenWork_History4.md`,
verbatim and in full — every measurement, retraction and refuted attempt is still there, cited below by section
title.  Earlier arcs: `OpenWork_History1.md` (threads A–E, runtime rounds 1–4), `_History2.md` (the Vxc repair +
fit-basis interface arc), `_History3.md` (Stage A 2026-09-08: the instruments, the head-to-head table, KP-0, TE).
The durable ⛔ findings left for `doc/Pins.md` (pins 18–22) rather than staying here as prose.

**The three rules this file lives by.**  (1) A ✅ goes to the current history file the day it is written, with a
one-line stub here.  (2) A row states its NEXT CONCRETE ACTION and the ONE record to read; the argument lives
in the record.  (3) A ⛔ that is durable becomes a pin; one that is in the weeds stays in the history — neither
stays here.

---

## 1. NEXT — the queued programme (agreed 2026-09-09)

Steps 1–4 are CLOSED: (1) KP-0 multi-k (2026-09-09); (2) the cleanup sweep + sprint S (2026-09-15); (3) TE test
suite as a product of axes (2026-09-15); (4) the plan polish (2026-09-16).  Records:
`doc/Records/OpenWork_History3.md`, `_History4.md` §"THE QUEUED PROGRAMME".  **A fresh session starts here:**

### ▶ 5. DFT+U — programme step 5 (`doc/Records/ParallelAndOraclePlan.md` Phase 3)

**Status 2026-10-01.**  Orbital-resolved +U (pin 23) is BUILT and CP2K-gated (increments 1–3); our ACBN0
estimator is refuted as a source of U values (pin 25).  The route is **our own linear response (A7,
`doc/LinearResponsePlan.md`)**; the oracle table (matched-PP hp.x: NiO 5.267, FeS2 7.535, SrVO3 6.25, ...) is
banked.  Open: A7 R3 step 4 (facade q loop + Si supercell gate), C1 32-atom MnO supercell sizing, Gate 2 spinel
sizing, C2 shifted-MP fold defect, Li q3 + basis-vs-PP check.

**Read `doc/HubbardUPlan.md` ▶ START HERE** for the exact next actions, gates and gotchas; the full history is
`doc/Records/HubbardUHistory.md` (including this section as it stood before 2026-10-01).

## 2. MAJOR FEATURES — capabilities the code does not have yet

Ordered by the battery roadmap (`doc/BatteryMaterialsRoadmap.md`: Tier 1 = GPW + GGA + spin + DFT+U + forces +
CE/MC on NC PPs; Tier 2 = USPP/PAW + k-point throughput).  `state` is what exists in the tree today.

| feature | why | state 2026-09-16 | next concrete action · record |
|---|---|---|---|
| **DFT+U** | essential for localised 3d — plain LDA/GGA gets voltages badly wrong | NOT STARTED; oracle validated (CP2K); the term faces (`MatrixForward/Adjoint`) exist | **§1 above** |
| **SCF checkpoint/restart (wavefunction + ρ to disk)** | crash-resume (the Oct 6–20 unattended window), a LIBRARY of converged states per material (user 2026-09-27: *"restart any material from a good state"*), the self-consistent-U outer loop across PROCESSES; AIMD's warm start later | **CK-1 ✅ DONE 2026-09-28** (design + execution record → `doc/Records/OpenWork_History4.md` §"CK-1"): `SolidCalculation::SaveState(path)`; `SolidCalcOptions::saveStateTo` (auto-save after EVERY `Converge` — each anneal stage, each ACBN0 outer step — converged or not; written aside and renamed into place); `SolidCalculation::Restart(path, lat, mol, opts, params)` → `Outcome<unique_ptr<SolidCalculation>, RestartRefusal>`, EXACT RESUME vs WARM START with every difference named.  Format `qchem-solid-state 1`: serial HDF5 through `qchem.HDF5` (src/Common), layout + fingerprint classes in `src/Calculation/SolidState.C`.  Gates: UTCalculation `SolidState.*` (7; Si 2×2×2: exact resume ΔE ≤ 2e-15 in 3 iterations Unpol-real and Pol-complex; grid-continuation warm start's FIRST iterate 2e-10 Ha from the converged answer; real↔complex TRIM blocks both ways; refusals on k-mesh / spin group / missing file), UTCommon `HDF5.*` (3) | **CK-2**: a WaveFunction read from disk (`tWaveFunction` without `tSCFWaveFunction` — anticipated by `SCFWaveFunction.C`'s ISP comment) so χ₀/ACBN0/gaps run on a stored state with NO SCF; every CK-1 file already carries C/ε/f for ALL orbitals.  **CK-3**: warm start onto a DIFFERENT k-mesh via the real-space \f$D(\mathbf R)=\sum_kw_kD_ke^{i\mathbf k\cdot\mathbf R}\f$.  CK-1 residuals, each its own small decision: (a) an exact resume still pays ≥ 3 Fock builds, and a NEAR-converged start CRAWLS — ρ_mix halves to 0.5/0.4 and DIIS resets on 1e-12 noise (the Ecut=30 warm start: 9 iterations from 2e-10 Ha) — an accelerator heuristic, not the restart; CK-2 removes the need for the exact-resume SCF; (b) the +U occupation ECHO (n recomputed from the restored D, checked against a saved copy) was NOT built: n is a function of D and the Hamiltonian recomputes it — build it the first time a state disagrees, or with frozen-occupation mode (where n IS independent state and must be saved); (c) no `Restart` overload for an annealed schedule and no mid-stage checkpoint (every N iterations) — add when a run needs one; (d) `h5py` is absent from the system Python, so the Python-side read is untested (the PySCF venv can take it).  States live OUTSIDE git: `~/Code/qchem6-runs/states/<material>/`.  QE `.save` samples: `IntegrationTests/QE/checkpoints/` |
| **DFT+U+J** (Hund's unlike-spin term) | the fractional-SPIN error — magnetic coupling in open-shell TMOs; +U alone leaves it (Macke outlook) | NOT STARTED; ACBN0 delivers J beside U from the same on-site ERIs (§1 step 3) | same term seam as +U, one more scalar per manifold; land after the +U anchor |
| **DFT+U+V** (intersite Hubbard V) | hybridised / charge-transfer insulators where an on-site term cannot restore the bond | NOT STARTED; Macke §5: orbital-resolved U already does most of what +V was added for, so it is a follow-on not a prerequisite | two-centre occupation numbers on the same projector; decide after the O-p manifold question on MnO is answered |
| **Linear response on REAL TRIM blocks** (found 2026-09-28, A7 R2) | every k of a Γ-only or Γ-centred 2×2×2 mesh is TRIM, so with the default LDA stack those blocks run REAL — i.e. almost every run we do (NiO/MnO/SrVO₃ k222, all the A6 materials).  The response path (R1/R2: `tResponse_HT<T>::GetMatrix(const tobs_t<T>*…)`, `OrbitalFrame<T>`, `TransitionBlock<T>`, `AO_TransitionFock<T>`) serves blocks of the RUN's scalar only, so a complex run with real blocks is REFUSED — loudly (`MakeOrbitalFrame` throws), never a wrong number.  Workaround today: the ground state with `forceComplex` — same physics (the 3c-3 gate pins ON == OFF), but complex arithmetic and memory on every block.  ⚠ R0's χ₀ is NOT affected (its Reference is scalar-agnostic): only the self-consistent χ (and so U₀) needs the workaround | NOT STARTED.  Blocks the U₀-vs-hp.x series only by COST, not correctness | two routes — (a) **real-block siblings now**: a `ResponseRealBlock` face per the `Dynamic_HT_RealBlock` idiom (Hartree/XC/+U each get a `GetMatrixR` over their existing scalar-generic `MakeMatrixT` body, which the δ-state already feeds), `FrameBlock` as a per-block scalar variant (as `cd_child_t`), `TransitionBlock` likewise (the composite already takes `tDM_CD<double>` children), `AO_TransitionFock` storing real blocks upcast; ~a day, but it ADDS to the 31 `*R` methods V1.35 exists to delete; (b) **ride V1.35** (`CleanupCandidates.md`: the two-parameter `tDynamic_HT<TBlock,TRun>`), after which the response face is one template with no sibling.  Recommendation: (a) only if the U₀ series is cost-bound on `forceComplex`; measure one material both ways first · `doc/LinearResponsePlan.md` §5d |
| **A7 R3 deferred: efficiency + clean-ups** (logged 2026-09-29 under the user's numbers → efficiency → clean-up order) | R3 lands the RIGHT NUMBERS first; these wait: (a) q-STAR reduction (one q per star, needs the site map for χ_IJ; ≤2× on NiO); (b) the Becke XC gather, the largest bucket per kernel application (~5 of ~11 s on NiO); (c) block-GMRES / warm start across channels and q; (c′) the transition gather at 3.2× the ground state (two `MomentsToPairs` passes; scatter 2.2×); (d) H5's full split (ACBN0's `Apply` onto `HubbardUTarget`, R4). (e) the raster (uniform-XC) route collocates each pair TWICE per kernel application -- once for the Hartree ball, once for the raw XC feed -- where the ground state's CollocMemo replays the second (R3 step 3); (f) the Becke transition forward/adjoint are two SERIAL, unscreened GEMMs per pair (the ground-state ones are threaded).  Also deferred, with their own rows: the real-pair route ("Linear response on REAL TRIM blocks"), CK-2 ("SCF checkpoint/restart").  (Retiring R2's periodic forwarding is DONE: R3 step 3.) | NOT STARTED — by design, until R3's NiO and supercell numbers are in | pick up after R3 step 5; measure first (the §3d sizing table is the baseline) · `doc/LinearResponsePlan.md` §3d "Deferred ledger" |
| **LRT and cRPA as first-class Hubbard-U estimation strategies** | today ACBN0 (on-site ERIs + a screening functional) is the ONLY in-tree route; hp.x and ABINIT's `ucrpa` are external oracles we run and compare against by hand.  User (2026-09-25, reading `doc/HubbardUPlan.md` §7's Carta et al. arXiv:2505.03698): no problem running every material in these papers and supporting BOTH methods natively | NOT STARTED; the abstract seam already exists (`HubbardProjection`/`HubbardUEstimator`, increment 3 — `ACBN0` is one concrete strategy behind it, doc/Pins.md pin 23's 2026-09-25 addendum) | LRT needs a perturb-and-respond capability (dn/dα on the manifold, doc/HubbardUPlan.md §5's DFPT-adjacent open question) we do not have yet; cRPA needs χ0 in a product basis (screening's own §2 route (d)) — both large, so this is a DIP placeholder not a near-term build.  Land ONE new concrete strategy behind the existing face when either becomes buildable, never a special case |
| **A THIRD independent U oracle: CP2K finite-α LRT (CK-alpha)** | A7's DFPT and QE's `hp.x` are both validated against each other, but a bug shared by two Sternheimer-family implementations wouldn't be caught by either — a THIRD oracle built on the ORIGINAL (pre-DFPT) Cococcioni–de Gironcoli recipe (converge at α=±1e-3 Ha, χ=Δn/Δα from the definition, no solver/kernel at all) is a genuinely independent check, and it's cheap: CP2K already had an unused `alpha` keyword, just wired to the wrong backend | ✅ **PATCH LANDED 2026-09-29**, `~/Code/cp2k` branch `ck-alpha` (commit `b5430fa`): CP2K's `&DFT_PLUS_U ALPHA` keyword existed but was `#if defined(__SIRIUS)`-gated — parsed only for the SIRIUS/pwdft backend, dead for the native Quickstep GPW path everyone (including our own CP2K comparisons) actually uses.  Un-gated the one keyword-read (`qs_kind_types.F`), threaded `alpha` through both GPW projector flavours (`dft_plus_u.F` `lowdin`/`mulliken`: `v_block(isgf,isgf) += alpha`, energy `+= alpha*trq`, and the `IF (u_minus_j==0) CYCLE` guard now also checks alpha so a PURE alpha-probe run at U_in=0 — the case this recipe needs — doesn't skip the atom).  `mulliken_charges` (a separate reporting/forces routine) deliberately left untouched, out of scope for a single-point probe.  Built as `cp2k_ckalpha.ssmp` (renamed via `OUTPUT_NAME`, never touches the pristine `/usr/bin/cp2k.psmp` oracle `doc/Records/CP2KBuild.md` validates every row against; fully static, no `.so` of its own).  **Smoke-tested, not yet numerically validated**: a native-PADE (no libxc) single-Ni deck with `U_MINUS_J=0, ALPHA=0.5 eV` gives a nonzero, iteration-varying DFT+U energy (0.0794→0.0849 Ha across SCF steps) — proof the mechanism fires every iteration, not a one-off; before the patch this would have been exactly 0.0 | ✅ **NUMERICALLY VALIDATED 2026-09-30**: matched Si deck (`IntegrationTests/CP2K/ckalpha/si_ckalpha_a{0,p,m}.inp`, FCC a=10.26 bohr, SIPP-SR/GTH-PADE-q4, Γ only, native `PADE` XC, `PLUS_U_METHOD LOWDIN`, U_MINUS_J=0 + `ALPHA [hartree] ±1.0E-3` on the Si1 p-shell (L=1) only — byte-for-byte the same system as `src/Response/tests/Kernel.C ResponsePolarizability.GPW_Si_Chi_eqFiniteDifferenceCDFT`).  **A UNIT TRAP found and resolved on the way**: CP2K's `dft_plus_u.F` scales its reported `trq` by `fspin` (0.5 for a restricted/RKS run, 1.0 for UKS — `src/dft_plus_u.F:317-319`), so `E_DFT+U = alpha*trq_reported` reads as HALF the true manifold occupation for our closed-shell Si run; read naively, χ comes out 2× too small.  Doubling (`n_true = 2·E_DFT+U/alpha`, exact since U_MINUS_J=0 leaves no other term) gives **n_true(+α) = 2.985484, n_true(−α) = 3.014516** — matching our own gtest's printed `n(+a) 2.98549  n(-a) 3.01451` to 5 significant figures — and **χ_CP2K = −14.5164 Ha⁻¹ vs our χ_FD = χ_LR = −14.5117 Ha⁻¹, a 0.033% relative difference**.  Two independent codes, independent DFT+U implementations (Dudarev/LOWDIN in both, but unrelated code), and — critically — one via a DFPT solver cross-checked against its own finite difference, the other via nothing but the finite difference itself, agree to 3 parts in 10⁴.  **This is the third-oracle validation the row exists for.**  Logs: `IntegrationTests/CP2K/ckalpha/si_ckalpha_a{0,p,m}.out` (not committed — regenerate from the decks; `cp2k_ckalpha.ssmp` and its build tree live outside this repo) | **NEXT**: (1) land a committed gtest or `gpwprobe` arm that automates this cross-check (today it was three manual `cp2k_ckalpha.ssmp` runs + hand arithmetic for the fspin correction) so it can be repeated on NiO/MnO without re-deriving the unit trap; (2) run the SAME recipe at U_in ≠ 0 with V_Hub FROZEN in the ±α runs (LinearResponsePlan.md's matching rules: freeze V_Hub, `doc/LinearResponsePlan.md` L674-692) — CP2K's own +U re-occupies every SCF step, so this needs either a CP2K-side freeze (not built) or accepting the unfrozen-vs-frozen gap measured on our own code (Si 2 eV: −17.65 unfrozen vs −14.42 frozen) as the expected offset; (3) separately, **a real but UNRELATED bug found on the way** (2026-09-29, still open): `-DCP2K_USE_LIBXC=ON` against qchem6's own patched libxc fork (`~/opt/libxc-qchem-patches`) crashes CP2K's `XC_FUNCTIONAL` section builder at startup — `"trying to add a subsection with a name (LDA_C_1D_CSS) that was already used"` (`input_section_types.F:590`).  Every committed CP2K oracle deck in this repo uses libxc-style split functional names (`LDA_X`+`LDA_C_VWN`), so THIS build can't reproduce them verbatim yet (native `PADE` works fine, no libxc needed, and is enough for the alpha probe) — worth its own investigation, separate from CK-alpha |
| **PBE / GGA** (then PBEsol, BLYP) | THE materials workhorse; #1 for the north-star | NOT STARTED; `Hamiltonian::Model` can list `PBE` with a "not wired" throw; collocation emits ρ only, not {ρ, ∇ρ} | build the ∇ρ collocation + the `∇·` term in the potential; prerequisite: retire the `GetEpsXc()=0.75*GetVxc()` base default (exact for Dirac exchange only — silent-wrong the day a GGA forgets to override; `CleanupCandidates.md` I.1 residual).  Spin-native from day one (pin 5) · `doc/OldPlans/FacadeDFTPlan.md` |
| **Hybrid functionals** (PBE0, B3LYP, HSE06) | molecular gold standard / solid band gaps | NOT STARTED; HF exchange exists for atoms/molecules (`Vxc`/`VxcPol` on `tDynamic_HF_HT`), **no HF for solids** (`IrrepCD<dcmplx>::AccumulateDirect` asserts out) | needs periodic exact exchange first — a real track, not an enum; design the canonical-pair scatter so periodic HF inherits it (`doc/OldPlans/ERI4Rework.md` §9) |
| **GW** (G₀W₀ quasiparticle corrections; RPA/ACFDT total energies ride the same χ₀ and W) | PARKED — a SPECTRAL capability (gaps, band alignment, photoemission), which the voltage curve does not ask for; RPA correlation is the total-energy cousin that WOULD matter for energetics.  Asked 2026-09-16 | NOT STARTED; what it needs: the unoccupied manifold (FREE in a Gaussian basis — PW codes pay for hundreds of empty bands), χ₀ in a product basis (= our density-fit basis, the RI basis GW codes use), Σ_x = periodic EXACT EXCHANGE (**missing**, the same prerequisite as hybrids), ε(q,ω) frequency integration + the q→0 head/wings (the genuinely GW-specific part).  O(N⁴) conventional, N³ low-scaling — 1–2 orders above +U on the same cell | stays behind +U, GGA, forces; moves up only if RPA total energies are wanted for energetics.  Oracle exists: CP2K periodic G₀W₀ in GPW (Wilhelm et al., low-scaling).  Dependency chain: periodic exact exchange → hybrids → GW/RPA |
| **LibXC-polarized** | the functional zoo in two spin channels | `Libxc_LDA` is UNPOLARIZED-ONLY by construction (never passes two channels; `Factory` throws for `SpinGroup::Polarized` + `XC::LibXC`) | pass libxc's `XC_POLARIZED` contract through the spin-native `ExFunctional` face; gate against `VWN5PolarizedMatchesLibxc` |
| **Forces** (Hellmann–Feynman + Pulay, incl. dV_PP/dR) | Tier 1; relaxations and phonons — and the PIVOT of every temperature-dependent property: τ(T) for band transport, polaron barriers, Li migration barriers, AIMD, lattice expansion all hang off it | NOT STARTED; roadmap hooks: keep `dV_PP/dR` computable, the collocation adjoint must support position derivatives | design note first (which terms need a gradient face) · `doc/BatteryMaterialsRoadmap.md` stage 3 |
| **Boltzmann transport (CRTA) + Seebeck** — the one transport capability that is nearly free for us | σ_αβ(T,μ) = e²∫dε(−∂f/∂ε)Σ_αβ(ε), Σ = Σ_nk v⊗v τ δ(ε−ε_nk).  In CRTA σ/τ and the τ-INDEPENDENT Seebeck S(T) come out; T enters through the Fermi window and μ(T) | NOT STARTED, but the ingredients exist: ε_nk at ANY k from the real-space lattice sums H(R), S(R) that `LatticeSum1E` enumerates (no Wannier interpolation — the step PW codes pay for); v_nk = i Σ_R R e^{ik·R}⟨ψ\|H(R)−εS(R)\|ψ⟩ is one more contraction; μ(T) = `TakeElectronsFermi` | store H(R)/S(R) (the B_ij(R) memo, row KP), a dense-k eigensolver loop (10⁴–10⁵ k), the velocity contraction, the Fermi-window integral.  Oracle: BoltzTraP2 fed our ε_nk.  ⚠ A real τ(T) = electron–phonon (phonons ⇒ forces; e-ph matrix elements = EPW-class) — behind forces.  Asked 2026-09-16 |
| **Polaron-hopping conductivity** — what LiMn₂O₄ / LiFePO₄ / MnO actually do | small polarons: the carrier localises on one Mn (Mn³⁺/Mn⁴⁺, the Jahn–Teller e_g electron) with its lattice distortion and hops thermally, σ(T) = (σ₀/T) exp(−E_a/kT); band transport is the wrong model | NOT STARTED.  Needs, in order: **+U with occupation CONTROL** to put the carrier on a chosen site (the frozen/constrained-occupation mode already in step 5's day-one interface; orbital resolution matters — the polaron IS the e_g electron) → **forces + relaxation** (the polaron is electron + distortion) → the hopping barrier E_a by NEB between site A and site B, or Marcus/Holstein from the two constrained diabatic states → attempt frequency (a phonon) | after step 5 and forces; the constrained-occupation half is the part we get first.  User 2026-09-16: *"I was not aware of 2"* — worth a note in the +U design that occupation control has TWO clients (LR-cDFT and polarons) |
| **Ionic (Li⁺) conductivity** | σ_ion(T) = n q² D(T)/kT (Nernst–Einstein), D(T) from migration barriers E_m; the occupancy/vacancy side is the lattice-gas / CE / MC stage | NOT STARTED; behind forces (NEB) and the lattice gas.  ★ **USER PIN 2026-09-16: E_m is a SET, not a number** — *"you have to try various routes (that are not identical under symmetry folding).  Each route has its own E_m (and multiplicity) … if two or more E_m are close you need to account for all of them."*  ⇒ enumerate the symmetry-INEQUIVALENT migration paths as ORBITS of (site → neighbour-site) pairs under the space group (our `Fold`/`SymOp` machinery is exactly the enumerator), NEB each representative once, carry each path's multiplicity z_i, and combine as D = (1/6)Σ_i z_i ν_i a_i² exp(−E_{m,i}/kT) when paths are independent — or **kinetic Monte Carlo on the lattice** when they couple (percolation, concentration dependence), which is the same lattice-gas machinery.  Never report one barrier | NEB is a forces client; the path-orbit enumeration is a symmetry client we could build now; AIMD is the brute-force cross-check (all paths at once, but needs high T and thousands of force evaluations) |
| **Molecular dynamics (BOMD / XL-BOMD)** — and AIMD as its consumer | Li diffusion by brute force, finite-T structure, the AIMD cross-check on the NEB path-orbit picture | NOT STARTED.  ★ DESIGN NOTE (user + Claude, 2026-09-16): Car–Parrinello is a PLANE-WAVE-native trick (fictitious electron mass; cheap H·ψ by FFT; 10⁵ coefficients to avoid diagonalising) — in a Gaussian basis H(k) is small and DENSE (118 → ~1000 on a 32-atom cell) and building it dwarfs diagonalising it, which is why CP2K's own GPW engine runs Born–Oppenheimer MD with density/orbital EXTRAPOLATION (ASPC) or extended-Lagrangian BOMD (Niklasson), not CP.  ⇒ what we would add, WITHOUT touching the abstract faces: (1) an MD driver above the SCF (Verlet + thermostat); (2) **a dedicated `tSCFIterator` type** whose SEED is the propagated/extrapolated density and whose loop-face runs a few iterations to a loose tolerance — the four role seams (`SCFStrategyPlan.md` §2) stay as they are; (3) **a warm-started `tLASolver` concrete** (Davidson/LOBPCG seeded with the previous step's eigenvectors, behind the existing `Solve`) — Lanczos per se buys little on dense n≈10³, WARM STARTING is what pays; OT/GDM warm-started from the previous orbitals is the direct-min equivalent (row OT).  ⚠ GPW-specific: the EGG-BOX effect — atoms moving across the fixed uniform fit grid ripple the forces; the atom-centred Becke mesh moves with the atoms and is immune | behind forces; then the iterator type is small.  Oracle: CP2K MD (same engine class) |
| **USPP / PAW** (augmentation facet) | Tier 2 — VASP-class affordability | NOT STARTED, deliberately (the roadmap's reframe: correctness physics first) | do not bake `S=I` or `ρ=ΣDχχ` into new PP or density code; pin 21 records why augmentation breaks "ρ<0 ⇒ D not PSD" |
| **OT — orbital transformation minimiser** | the minimiser CP2K ships; the ONLY way to time a minimiser against theirs; gates lever C | NOT STARTED; the role seams anticipate another direct minimiser (`SCFStrategyPlan.md` §7); GDM as built is fixed-occupation (pin 15) | build OT beside GDM under the accelerator/loop-driver seam, WITH the smearing-aware search direction (OT+smearing is a coupled orbital+occupation minimisation, not forbidden); then time OT-vs-OT and re-open lever C · `doc/Records/SCFStrategyPlan.md` §7, `doc/Records/OTNotes.md` |
| **k-point parallelism** (row KP) | the one parallel axis every other code has (CP2K `PARALLEL_GROUP_SIZE`, confirmed in source; QE `-nk`, VASP `KPAR` believed, not verified) | ✅ the pre-warm exists (`tHamiltonian::RefreshForDensity`, 2026-09-08 — every k-independent memo is warmed before the block loop); ~~one write remains INSIDE the loop (`itsCache` keyed by `Irrep`)~~ **STALE, checked 2026-09-27: R1.0h's `HT_SlotOwner::PrepareSlots` pre-creates those slots before the loop** (`HamiltonianTerm.C`); `itsByL/itsByLSeen` is the opt-in `GPW_NL_PER_L` diagnostic only | shared prologue → read-only parallel k loop → density reduction.  ★ **MEASURED 2026-09-27 — the embarrassment has arrived (A7 R0, NiO AFM-II +U, k 2×2×2 = 16 blocks, `gpwprobe nio`):** `perf` on the SERIAL run: block-loop Fock build (gather) ~44 % of cycles, collocation `Contract` (OUTSIDE the block loop) ~38 %, ~1.2 of 16 cores busy.  **Lever 0 was free: `GPW_OMP_THREADS=12`** (default 1) → ~8.5 cores busy, setup 258 s → 79 s, ~40 s/iteration.  Order proposed: (1) threads ON for every multi-k run (recipe), (2) the cross-k gather memo (§4b "k-scaling gap": 16 blocks over only 2 distinct fields ⇒ the gather is ~8× redundant — work REMOVAL), (3) KP for what remains of the block loop — its ceiling on this box is now the ~2× of idle cores, and it must be designed against the inner pair-loop OMP (pin: nested parallel sections).  Not started · History4 row KP |
| **Space-group irreps as block labels** (T⋊P — the next G row of pin 14) | Tier-A BZ reduction beyond the point-group fold; the IBZ weights are per-k already | the fold under the MESH-symmetry subgroup exists (KP-0); irreps of the little group do not | lands in `qcSymmetry.Lattice_3D` + the `Gaussian.Lattice` container, touching no library boundary · `doc/Records/SymmetryUpgradePlan.md` §9 |
| **SSB descent** — the methodical symmetry route | today MnO's AFM-II is ASSUMED, never derived; the free run is first-class but blind | design in `SymmetryUpgradePlan.md` §3b; exists: `SymmetryDefects`, the ops chokepoint; does NOT exist: `Impose::Subgroup`, subgroup closure, crystal irreps, **CD persistence (nothing at all)** | step 3 cannot be a single free iteration (the symmetric solution is a stationary point; SSB is second-order) — it must measure GROWTH or CURVATURE · `SymmetryUpgradePlan.md` §3b |
| **The second magnetic material** — and the derived order parameter | only the CELL is material-specific (user, 2026-08-25); `m_stag` hardcodes MnO's two Mn sites and a 0.7-bohr offset | `siteSpins` decoration + integrated site moments both live on the run now, so \f$m_{order}=\sum_A\sigma_A\mu_A/\sum_A|\sigma_A|\f$ is derivable for ANY collinear ordering | derive it when the SECOND material arrives (one cell cannot tell material-specific values from MnO's accidents) · History4 "ONLY THE CELL IS MATERIAL-SPECIFIC" |
| **Fermi smearing, finished** | metals and the anneal | Fermi–Dirac per block + `GlobalMu`/`ShFermi` reservoirs exist; kT is a hand knob | (a) principled kT tied to the gap / DOS / a target entropy (pin 15 is the constraint); (b) Gaussian / MP / cold flavours + the ½(E+A) T→0 extrapolation reported beside E, −TS, A — decide per battery need (FD is the true finite-T physics; MP/cold are numerically nicer for a T→0 answer; ⚠ MP/cold give NEGATIVE occupations, pin 21's canary) · `doc/Records/GPWPlan1.md` "Future considerations" |
| **Molecular spin-resolved SAD seed** | the Hund-split tables exist; only the PW `SeedCD` reads them | molecular `NumericCD` still hands a polarized run ρ/2 per channel | channel-aware `NumericCD` assembly + an O₂-triplet gate (Hund-split seed vs ρ/2: same basin, fewer iterations).  Feature wish, not a defect |
| **Diffuse-basis ACTUATOR** / the vet-stage symmetric trim | the detector (`PivotedCholeskyDrops`, `basis.removed` report) landed; nothing ACTS | no `Prune(indices)` exists; ortho-time pivot filtering is the only trim | decide auto-prune (`IrrepBasisSet::Prune` + `BasisSet::Prune`, RKB prunes paired large+small) vs report-only-forever (the 80/20: the user reruns either way).  **Pin 22 governs the shape**: vet-stage, a property of S made once, whole ORBITS never single AOs, reported as a basis not indices; ortho-time filtering becomes the fallback.  Free controlled experiment: `VALENCE_LOWQ_VA` under Cartesian d is rank-deficient by 10 and auto-drops to exactly SR's 122 — a null control; runs 58–60 (O₁ p(0.18) vs O₂ s(0.15)) are the real evidence.  ★ **DECIDED 2026-09-23 (user, pin 22 addendum): AUTO-PRUNE, not report-only** — the vet-stage trim is an automatic trim and stays one; what is refuted is auto-trimming in the SCRIPT (retired) and at ORTHO TIME.  Ortho-time gets two behaviours only — shut up and work, or make noise enough to fix the BASIS — and the cure for its bare-index report is an **exception** carrying indices/pivots, caught where the basis is owned, NOT a decorated `LASolver` return type (user: cleaner than fallible flags).  ⚠ **Needs a slew of unit tests, and two axes they must cover**: (a) under CARTESIAN d/f a shell is not a pure *l* — the *l*−2 contaminants (s in d, p in f) mean a "whole orbit" carries two characters at once; (b) the rank decision depends critically on LATTICE SPACING, so a trim validated on one cell says nothing about a denser one · History4 "Continuous — CLEANUP" |
| **libcint lattice engine** (`PG_LibCint` realises `LatticeSum1E`) | the faster 1E/KB/3C engine for GPW | = the ISP-split `LatticeSum1E` that is V1.38's prerequisite | `CleanupCandidates.md` V1.38 (stashed) |
| **Spherical SALC S3b** — the libcint-spherical extractor | the one empty cell of the molecular grid | S1–S5 done; the in-house spherical SALC ships without it | must match libcint's real-harmonic ORDER and NORMALISATION (a foreign convention); libcint-spherical presents AS a `PGData` with spherical components (a trap).  Genuinely separable · `doc/OldPlans/SphericalSALCPlan.md` |
| **valgen `--auto`** — rule-based valence-window refinement | the next d-metal basis | plan, not built | build when a second d-metal basis is needed · `GPWPlan1.md` "valgen --auto" |
| **The run report for the GUI** | the binding side's consumer | `qchem.Reporting` complete; the wishlist waits on a consumer | `meta` section (cheap, one `EmitSection` in `Converge`); field-metadata registry (code key / terse label / hover text, render-side only); detail-level FILTER; `Renderer` DIP split; `RollingFileSink`; `basis.removed` named by L/α/atom; `schemaVersion`; structure/symmetry/irreps/Hamiltonian sections; HDF5 sidecar; a **"literature units" display toggle** (idea, user 2026-09-24, prompted by ABINIT reporting +U's U/J in eV but PW cutoffs in Ry — every code picks its own convention per QUANTITY, not one global unit system) — per-field unit choice in the metadata registry above, not a second internal unit system · `doc/Records/RunReportPlan.md` |
| **Lattice gas / cluster expansion / MC** | Tier 1's last stage — the actual voltage curve | specced, deliberately not built | `doc/LatticeGasPlan.md` |
| **A second oracle** (QE first; Elk for all-electron; GPAW/SIESTA only as GPW-like timing peers) | a question ONE oracle cannot answer | `~/Code/q-e`, `~/Code/abinit`, VASP all built (CLAUDE.md); nothing run against them yet | ⚠ TRIGGERED ONLY by such a question — today that is §4 row "MnO accuracy" (−99.7 mHa vs CP2K, operator not named).  Also cheap and their regime: CP2K's own 32-atom MnO supercell, to close the `Benchmark.md` §7b caveat (is their 1.09× OMP a never-parallelised route or too few tasks) · `doc/Records/ParallelAndOraclePlan.md` Phase 4 |

---

## 3. NON-OOD CLEANUP — tooling, build, hygiene

(The SOLID/OOD debt list — faces, casts, term architecture — is `doc/CleanupCandidates.md`.  This is the rest.)

| item | what | next concrete action |
|---|---|---|
| **Module toolchain** | `import std;` + a modular Blaze fork — banish the preprocessor; Clang-only modules | `doc/ModuleToolchainPlan.md` (not started); `doc/Records/cmakenotes.md` has the build notes |
| **Env-knob graduation** (GPWPlan1 item 1; N5's remainder) | **~45 `GPW_*` `getenv` knobs remain in `src/`** (`grep -rn 'getenv("' src`); `raster`/`cutoffFactor` are TYPED options the `CP2K_COMPAT` policy does not reach | grade each: verification instrument / ops valve (documented env, keep) vs a setting a user might tune (a typed `SolidCalcOptions` field, or the `CP2K_COMPAT` deviation list).  Pin 12 is the criterion.  ⚠ `GPW_BECKE_L/NR/ALPHA` are consulted ONLY for arguments passed `<0` — they silently do nothing against a caller-supplied degree |
| **`GPW_OMP_THREADS` naming and semantics are confused** (found 2026-09-30, CK-alpha/MnO Mn-3d session: a user-observed "500% CPU" on a run whose own console banner said `GPW_OMP_THREADS=1 (BLAS pinned to 1)`) | Thread control is not a GPW/basis-set setting — the name says otherwise, and it is not even the WHOLE story: `qchem::WorkerThreads()` (`src/Common/Parallel.C`) gates the GPW pair loops, the PW Hartree/XC/PP terms and the response GEMM (all default to 1, honoured), but `src/Structure/Imp/UnitCell.C`'s Becke-mesh construction is a SEPARATE, DELIBERATELY-always-parallel region (bit-identical at any thread count, so the authors let it default to the OpenMP runtime's own team size — this box's 16 — and only CAPS it if `GPW_OMP_THREADS` is set) — which is exactly what fired during the MnO run's setup.  User (2026-09-30): we want ONE flag, `QCHEM_OPENMP_THREADS`, and `=1` must mean 1 EVERYWHERE, including inside BLAS/LAPACK (today's `PinBlasToOneThread` calls `openblas_set_num_threads`, unverified against the reference `liblapack.so.3` this binary also links — a routine that happens to be single-threaded by construction, not by any explicit pin, is a latent trap the next linked LAPACK could break).  **Open: what should `QCHEM_OPENMP_THREADS=0` mean** — decide before renaming, since `GPW_OMP_THREADS`'s unset-vs-explicit-0 distinction does not obviously generalise | rename `GPW_OMP_THREADS` → `QCHEM_OPENMP_THREADS` (one flag, no material/module-scoped variants); audit every `#pragma omp parallel` site (`grep -rn "pragma omp parallel" src`) against it, starting with `UnitCell.C`'s Becke-mesh region; verify the BLAS/LAPACK pin actually holds at `=1` (a quick `ps -L`/thread-count check on a run that does real Cholesky/eigensolve work, not just GEMMs); rule on `=0`'s meaning (candidates: "no cap, OpenMP's own default" [today's unset-GPW_OMP_THREADS behaviour for the Becke region] vs "explicitly serial, same as =1" vs "reject the value") and write it down here once decided |
| **External run-data organization** (found/fixed 2026-09-30: user noticed `IntegrationTests/QE/tmp/` showing dirty in every `git status`, which raised the bigger question) | Four codes (CP2K, QE, ABINIT, our own `qchem6`) each generate decks/logs/checkpoints, and the storage had grown inconsistent: `~/Code/cp2k-runs/` was one flat pile of ~70 files for 6+ materials (plus a 460 MB duplicate `.kp` restart pair); `~/Code/abinit-runs/` mixed loose files with ad hoc subfolders; there was no `~/Code/qe-runs/` sibling at all, so QE output split three ways — committed decks (correct), untracked scratch left inside the repo (the `git status` trigger), and one run (`hp_NiOgO`) dumped under `qchem6-runs` though the compute was QE's; `qchem6-runs/` itself mixed campaign-code folders (`a7_r2`, `a7_r3`), ad hoc names (`gate1`, which turned out to hold TWO materials), and material names as siblings, with only `states/<material>/` following a sound pattern | ✅ **DONE 2026-09-30 (this session): the convention is now `~/Code/{cp2k,qe,abinit,qchem6}-runs/<Material>/<descriptive-name>.{ext}`**, material spelled consistently with `materials.json`/`doc/HubbardUPlan.md`'s table, elemental/PP-validation runs not tied to one material under `_atoms/`, checkpoints under `<Material>/states/`.  Applied: `.gitignore` now covers `IntegrationTests/QE/tmp/` + `*.UPF` (closes the git-dirty trigger); `cp2k-runs`, `abinit-runs`, `qchem6-runs` all reorganized into this shape (pure `mv`, nothing deleted); `qe-runs/` created and `hp_NiOgO` moved into it as `qe-runs/NiO/`; `gate1` split into `LiMn2O4/` and `MnO2/`.  **Still open**: (1) the 460 MB duplicate `mno_afm2_gpw_sr_222g-RESTART.kp`/`.kp.bak-1` pair in `cp2k-runs/MnO/` — a deletion, deliberately NOT done without asking; (2) `qchem6-runs/batch/`'s job-queue infrastructure (`queue/done/lock/logs/work/`) should stay pure orchestration scratch, but today's jobs left their real output sitting in `batch/work/<job>/` rather than "graduating" to `<app>-runs/<material>/` on completion — either do this by hand each time or teach the job scripts (or the driver) to copy the deliverable files to their proper home when a job lands in `done/`; (3) no enforcement — a future session can still drop files in the wrong place, this is a convention not a guard.  **Also raised, unresolved**: the user recalled a web archive for computational-materials-science run data and wondered about using it instead of trusting a local SSD with everything — candidates found (not yet acted on, needs the user's own account, not something to sign up for on their behalf): [NOMAD](https://nomad-lab.eu) (code-agnostic, parses raw input/output from many codes including CP2K/QE/ABINIT, zip-upload) and [Materials Cloud](https://www.materialscloud.org/) (AiiDA/QE-centric, same EPFL/Marzari group behind `hp.x`; Carta et al.'s own KCuF₃/Sr₂FeO₄ inputs already reached us via a `materials_cloud_submission` folder, so there is already a live precedent for using it).  If adopted, it would sit ALONGSIDE this local convention (an archive for finished/citable results, not a replacement for day-to-day scratch), not instead of it |
| **`MinΔρ` compares three different quantities** (the A4 Δρ/N gate) | `LinearMixer::Mix` normalises by charge, `ReDampMix` does not, `KerkerMixer::MixField` takes an ∞-norm over G-coefficients — one threshold, three scales, so `MinΔρ` means something different per recipe and per system size; and Δρ measures the MIXER'S STEP, not the distance to the answer (two mixers floored 8e-7 Ha apart at Δρ values 5× apart) | one intensive gate, Δρ/N, in one place (`doc/Records/SCFStrategyPlan.md` "Δρ convergence gate should be intensive"); anchor-moving ⇒ sprint S (§4).  Related: `DidConverge()` disagreeing with the run's own "DENSITY-DEGENERATE (benign)" fingerprint, and the fingerprint's overconfident "raise NMaxIter" advice |
| **The SCF banner names the wrong mixer** | `Calculation.C:122` prints `Pul` when Pulay depth > 0 and never mentions Kerker, though `MakeGSpaceMixer` composes them (`[Pulay] ENABLED … G0=1` has the truth); a benchmark reader would conclude the preconditioner was swapped out.  With `KerkerG0=0` the fallback to linear D-mixing prints no mixer line at all | one line: print the composition (`Ker+Pul(8)`), and print the fallback |
| **`_Long` test budget + ctest label** | TestSuitePlan ruling 5 proposed 60 s CPU (de-facto ceiling 74 s, `GPW.XCPotentialConsistencyFD`); no `ctest -L long` label exists; the 5-minute MnO AFM-II gate is `DISABLED_` for that reason | set the budget, add the label, promote `GPW_MnO.DISABLED_Γ_Shub_Pol_Smear_Anchor` under it |
| **`GPW_CONTRACT_CUBE=0` can rot** | the reference box walk is an investigation opt-out nothing exercises by default (792/792 both ways on 2026-08-27) | run `GPW_CONTRACT_CUBE=0 ctest -j8` beside the plain sweep at any breakpoint that touches the collocation path |
| **`M_PG_BoxWalk.WhereTheContractionSpendsItsTime` flakes under `-j8`** | it asserts on wall-clock | assert on a ratio of counts, or make it an instrument (`gpwprobe`) |
| **`GPW_NaF.DISABLED_Γ_GridContinuation`** | needs a facade grid-continuation face | small facade addition, then re-enable |
| **The PW basis family has no facade** | `SolidCalculation` is built over a Gaussian basis; the 2 PW SCF drivers live in `IntegrationTests/PW/Harness.C` | grow a basis-family axis on the facade, then retire the harness |
| **`MakeWater()` etc. inline in the molecular tests** | `qchem.Materials` has the 4 molecules; the molecular tests do not read them yet | mechanical, no `src/` consumer, not urgent |
| **Source comments cite the old flat `doc/` paths** | `doc/GPWPlan.md` ×35, `doc/SymmetryUpgradePlan.md` ×59, `doc/RealComplexPlan.md` ×44 … now under `Records/` or `OldPlans/` | user: fine as is (findable by name).  If ever done: one sed, its own commit, when nothing is building |
| **Low-rank D: qualify the code's own claim** | the doc string says *"the rank is the same from tol 1e-6 to 1e-12, so the cut is not a tuning decision"* — a kT=0 statement; at kT=5e-3 the pivot spectrum has a thermal tail at ~7e-6, a 4–5 decade gap, not 12 | qualify, do not delete; make the pivoted-Cholesky failure LOUD and route to the eigen split (pin 21) |
| **Build flags** | ✅ settled: Release `-O2` is a DELIBERATE choice (CMakeLists.txt), `-march=native` the Release default since 2026-09-04 (it licenses FMA contraction — bit-identical anchors are re-checked, not assumed); `QCHEM_ARCH_EXPERIMENT` is the opt-in A/B | nothing |

---

## 4. REMAINING TODO — measurements, performance, accuracy

Each row: what is open · the next concrete action · the record.  Priority is roughly top-down inside each group.

### 4a. Accuracy

| row | what is open | next concrete action · record |
|---|---|---|
| **Shifted-MP fold lowers Si by 1.02 mHa** (found 2026-09-20 by the rule-3f re-take) | `GPW_Si.k222s_Imp_CP2K`: imposed −7.868473429 under EVERY criterion and both XC grids; FREE −7.867452508 = CP2K's −7.867436530 to 16 µHa.  The k=±¼ shifted 2×2×2 MP mesh is the suite's only fractional-k (non-TRIM) sampling, and the point-group fold of it — star weights, the ρ star-average, or the TRIM/complex block split — moves the energy by 1 mHa DOWNWARD.  Charge is 8 in both, so it is not KP-0's weight sum | first PIN which half: (a) `QCHEM_IMPOSE_SYMMETRY=1` with the ρ star-average disarmed vs armed (R1.0r's question — the shifted mesh is invariant under the cubic group, but `DetectPointOps` star-averages ρ under the full 48 whatever the mesh); (b) the folded k-set and weights printed against the unfolded 8 (pin 13's integrality test passes — does the STAR SIZE?); (c) the same A/B on the Γ-centred 2×2×2 (all TRIM), which agrees imposed/free.  A 1 mHa fold error on the fractional-k mesh is a sampling defect that every future k-mesh run inherits; fix before KP-1 · `doc/Benchmark.md` §5a re-take, `feedback_complex_type_vs_value` |
| **NiO VA: the diffuse Ni s decides insulator vs metal** (found 2026-09-27 by A7 R0) | NiO AFM-II U=3 eV, k 2×2×2, `orthofull`: at `orthoTol=1e-4` (pivoted Cholesky keeps pivots to 2e-4) every recipe — free Ladder+MOM, free deck-shaped, Shubnikov-imposed — converges GAPLESS (one empty Γ level below occupied levels elsewhere; E = −106.7348); at `1e-3`, which drops each Ni's FIRST AO (the α = 0.06 s of the VA block `s 7 × {0.06..28}`), it is a 2.18 eV insulator at E = −106.3729, **0.36 Ha HIGHER**.  hp.x on the same deck: insulating, 2.86 eV.  ⚠ HubbardUPlan's NiO ACBN0 rows were taken at 1e-4 (Γ-only, where the Γ ghost has no other k to invert against) ★ **2026-09-27 follow-up: the XC QUADRATURE IS EXONERATED** — at 1e-4, imposed: Becke default −106.7348334, Becke nR 80 / L 35 −106.7334395, uniform 96 Ha −106.7334158, all gapless with the same Γ level (ε ≈ 0.2099); fine Becke = uniform to 24 µHa at under half the points (default Becke sizing 1.4 mHa off — a separate data point for "Size the Becke grid").  ★ **The Γ level IS QE's CBM, 3.4 eV too low**: hp.x's `NiOgO.scf.2` has its conduction minimum at Γ (band 17, 13.33 eV; Γ direct gap 3.5 eV) — the dispersive Ni-4s/O-3s band — ours has a Γ direct gap of ~0.13 eV.  A state BELOW the plane-wave reference cannot be basis incompleteness (Ritz upper bound at fixed H): it is ILL-CONDITIONING × INTEGRAL ERROR — the state lives in S's near-null direction (pivot 2e-4) where ε ≈ c†Hc/c†Sc amplifies tiny H/S errors (pin 7); the 0.36 Ha drop is the same ghost (GPW is not strictly variational).  ⇒ the 1e-3 insulator (2.18 eV vs QE 2.86) is likely the PHYSICAL state ★★ **2026-09-27, the VET-STAGE TRIM landed (`c0b1f4fd`, pin 22) and exposed a SUBSET DIVE.**  The vet removes Ni s {0.06} from both sites (ortho then keeps 126/126 in every block — a consistent basis at last), and **"VA minus Ni s 0.06" DIVES: E = −134.27 (27.5 Ha BELOW the full VA), Een −189.9 vs −159.0**, 23 Ha low already at iteration 1 from the same IonicSAD seed.  A SUBSET must be variationally HIGHER (Ritz).  Excluded by A/B (3 iterations each, identical to 2 mHa): the vet's in-process trial builds / any cache (a STATED trim built once, `NIO_TRIM=28:0:0.06`, reproduces −125.517135479 to all digits); the lattice-sum reach (`GPW_SCREEN_EPS=1e-14`, 339 cells); the V_loc-long G-ball (`GPW_LONG_SWEEP=1`).  eps 1e-14 on the FULL basis moves E ~20 mHa and leaves the Γ ghost.  ⇒ the SAME CLASS as "The 136-function span" (MnO 132-of-136 dove 20.7 Ha): **some term's matrix elements are not SUBSET-INVARIANT** — ⟨φ_i|V|φ_j⟩ between two KEPT functions changes when a third leaves.  The orthoTol 1e-3 "insulator" never met it only because its per-block drops were inconsistent (index 1 = Ni s 0.167 in one block, 0 in another).  Logs `~/Code/qchem6-runs/a7_r0_nio/` ★★★ **RESOLVED 2026-09-28 — it was a PSEUDOPOTENTIAL GHOST STATE in a near-null direction, and the vet trim at orthoTol=1e-3 cures it.**  The −134.27 "dive" console showed an OCCUPIED level at **−36.55 Ha at ONE k, (0,½,½)** (EenNL/Eloc blown), NOT a subset-invariance failure: the trimmed basis at 1e-4 still had min eig S = 3.2e-4 at the same zone-boundary k-family (the worst k is the ZONE BOUNDARY, (½,½,0): 9.8e-7 before any trim — not Γ).  **`NIO_VET=1 NIO_ORTHO_TOL=1e-3` trims RAW SHELLS Ni s {0.06} then Ni d {0.18}** (the Mn d 0.18 of MnO's `va`), min eig S 9.8e-7 → **4.7e-3**, ortho 116/116 in EVERY block: E = **−106.2021** (Ritz-consistent: above the full basis), Eloc/EenNL −85.14/−75.25 (full VA −85.16/−74.49), lowest level −0.42 (the full basis's lone −1.29 Γ level was the same artefact), **INSULATOR 1.30 eV** (QE 2.86), CBM at Γ as in QE, 15+6 iterations.  MnO's history, for the record (the .bsd spans): its trims were O p 0.18 + O s 0.15 (+ Mn d 0.18 in `va`); **NiO's Ni block was NEVER vetted** (`va` Ni = `sph` Ni).  Remaining: the 1.30 vs 2.86 eV gap is now a basis/physics comparison, not a pathology | (0) ~~a SUBSET-INVARIANCE unit gate~~ (superseded: the console named a ghost, not a subset defect); OPEN: why the near-null direction hosts a −36 Ha KB ghost at all — a trimmed basis should be merely worse, not attractive (a unit gate: the KB nonlocal matrix's lowest generalized eigenvalue vs λ_min(S) on the 1e-4-trimmed NiO) — per term (S, T, local PP, KB, Hartree/XC at a FIXED density), the shared-function block of full-basis vs subset matrices must be IDENTICAL; names the guilty term with no SCF, and is the permanent precondition for the vet trim being safe (it is committed but must not be used on NiO VA until this is fixed); (1) DONE: 1e-4 with `GPW_SCREEN_EPS=GPW_DENSITY_EPS=1e-14` (was 1e-10) — if the Γ level rises and E climbs toward −106.37, precision × conditioning is confirmed and the cure is pin 22's vet-stage trim of the VA Ni s, NOT a looser orthoTol; if nothing moves, audit the KB nonlocal PP matrix elements on the diffuse s; (2) name the Γ level's character (s-weight on the α = 0.06 functions); (2) is the 0.36 Ha REAL — re-run the 1e-4 state at a finer density grid (`NIO_ECUT` up) and under `CP2K_COMPAT`; a collocation artefact moves, physics does not; (3) a CP2K single point with the same basis if (2) is ambiguous; (4) if it is the near-dependence dive, the cure is pin 22's vet-stage trim of the VA Ni s, and the Ni block's valgen record gets the note · `doc/LinearResponsePlan.md` §5b, logs `~/Code/qchem6-runs/a7_r0_nio/` |
| ~~**Linear response on an IMPOSED run symmetrizes δρ**~~ ✅ **CLOSED 2026-09-29 by A7 R3 step 3** | the R2 composite route it describes is retired; the pair route has no fold; imposed Γ raster χ == free χ | → `doc/Records/OpenWork_History4.md` §"Linear response on an IMPOSED run" |
| **GDM declines every UNPOLARIZED (folded-doublet) run** (found 2026-09-29; user: *"there must be a way to just make it work for unpol systems"*) | Both engagement checks in `SCFAcceleratorGDM.C` (`UseFD`) assume occupations 0/1: idempotency tests Tr D′ = Tr D′² and block occupancy tests Tr(D′P) = N_occ, but a folded doublet has D′ = 2P (log: `Tr(D')=8 vs Tr(D'^2)=16`).  So GDM, and a `Ladder`'s GDM rung, silently fall back to diagonalising steps on every unpolarized run -- a restart rule "use GDM or DIIS->GDM" would read as satisfied while it is not. | test D′/g (g = the level capacity: 2 folded, 1 per channel -- the Reference's own `g`); then CHECK whether the step itself (gradient ∝ g[F′,P], trust, energy model) needs g too. Gate: a closed-shell unpolarized run converged by GDM == its polarized singlet twin (the Pol == UnPol shape); include a molecular closed shell · `SCFAcceleratorGDM.C` |
| **Print S(k) conditioning on EVERY run** (logged 2026-09-29; user rule proposal: "always check the [basis trim] output for conditioning problems") | `[basis trim]` is printed only by the vet stage (`gpwprobe <P>_VET=1`), so unit tests, the gtest harness and most facade runs never show S(k) -- the Si response test (min eig 0.027, cond 98) had to be measured by hand. | the facade prints one line per k-block at setup: min eig S(k), max, condition number (one eigen-decomposition per block), beside the `[ortho]` pivot line; then the CLAUDE.md rule is checkable on every run · `SolidCalculation` setup |
| **MnO accuracy — name the operator** (was Step 5 / item 4) | the VA (N=118) exact-span table, both codes at full rank, zero span/symmetry/ensemble excuses: qchem AFM −61.40298 / FM −61.44158, CP2K −61.30333 / −61.30478 ⇒ a **−99.7 mHa configuration-BLIND offset** (below a variational reference on an identical span ⇒ an operator or convention, not a basis; suspects: G=0 / alignment conventions, XC quadrature, V_loc bookkeeping) and a **−37 mHa configuration-SELECTIVE, span-independent FM-favouring bias** | ⚠ first PIN the XC ρ source (`GPW_XC_DM_SOURCE`) — individual terms move ~100 mHa with it while the total moves 8 μHa.  Then the CHEAP move: term-by-term against CP2K's energy blocks (Ekin/Eee/Exc/E_loc/E_NL, `Een` split V_loc/V_nl) on the **Mn ATOM** first (seconds; oracle banked at −14.2440 restricted / −14.658 polarized), then the crystal; a configuration-blind offset should show on one atom.  Then `MNO_KMESH=2` (k-convergence moves the ordering more than the physical 6J₁+12J₂ ≈ 4 mHa).  One control: the AFM and FM arms may not be converged under equally safe conditions if Kerker's charge damping is load-bearing (pin 18).  ⚠ re-read every "weak-moment basin" conclusion against the INTEGRATED moment (pin 4) · `doc/Records/SphericalLatticePlan.md`, History4 "Step 5" |
| **N4 — the cusp-deficit XC feed** (pin 18's open half) | \f$\rho_{XC}=\rho_{mix}+(\rho[D]_{exact}-\rho[D]_{BL})\f$ is a DESIGN ARGUMENT, not a result; `SCFParams::XCCuspDeficit` exists, OFF (sprint item A6).  It also carries a bin-1 prize: lever B (one gather per spin) is blocked ONLY by XC needing the raw ρ≥0 feed — routing \f$V_H\f$ through the raw adjoint moves the Hartree block 6e-5 relative (ball vs raw truncate differently), so B is refuted UNTIL N4 lands | two measurements before it is believed: (a) is \f$\rho_{mix}+\Delta_{cusp}\f$ pointwise non-negative (`GPW_RHO_NEGATIVE` census answers it); (b) is \f$\Delta_{cusp}\f$ as iteration-static as the core-electron argument claims.  Gating instrument: the radial spectrum of \f$(\alpha f-1)\tilde\delta\f$ vs |G| per iteration — computable entirely inside the mixer.  Anchor-moving ⇒ sprint S · History4 "N4" |
| **N3 — charge and spin mix in different channels** ✅ **PROMOTED 2026-09-20** | Kerker per spin channel damps the spin channel too, and the spin channel has NO 4π/G² to justify it (cf. VASP `AMIX_MAG`); `QCHEM_MIX_RHO_M=1` selects the (ρ,m) basis with Kerker on ρ and plain linear on m | MEASURED on MnO AFM-II with the deck's loop shape (`gpwprobe mno`, Null accelerator, no MOM, Pulay 8, max\|ΔD\|<1e-6): **(ρ,m) 18 / 21 iterations vs (up,dn) 22 / 24** (U=0 / U=4 eV), energies identical to 1e-10 Ha, m̃(q) 3.13/3.16 e.  Pin 18's "will LOOK like a regression" did not materialise — that warning was written against the two-history recipe.  ⇒ `QCHEM_MIX_RHO_M` qchem default flipped to ON (still a declared deviation, off under `CP2K_COMPAT`).  The G=0 carve-out that keeps the TOTAL moment mixable stays (`FourierMixCD.C:36`) · History4 "N3" |
| **Sprint S — the anchor-moving batch** | A1/A5/A7 landed (A5 moved nothing: a converged run lands on the same number from either seed).  Left: **A2** V1.22 (Becke partition once per orbit representative — also up to ~23× on the largest CPU bucket: 98816 points, 4290 orbits; imposed runs only), **A3** §K (fit-{G} densification, deferred by ruling), **A4** the Δρ/N gate (§3), **A6** `XCCuspDeficit` (N4) | do them in ONE re-bank window so each delta is attributable; V1.22 and §K are `CleanupCandidates.md` rows · History4 "THE ANCHOR-MOVING SPRINT" |
| **The 136-function span** (was Step 6; time-boxed) | why CP2K holds the full diffuse span and we must strip it.  ⛔ the screen-discipline hypothesis is REFUTED (run 64: 132-of-136 at eps 1e-12 still dove to −82.19, 20.7 Ha below CP2K on a SUBSET of its span) | surviving candidates: (1) SVD/eigen-consistent F/S filtering; (2) the symmetry-INEQUIVARIANT AO drop — the vet-stage trim (§2) is the near-prerequisite; (3) project the near-null directions out of F as well as S.  Time-boxed research; do not let it grow into a track.  Re-run with `GPW_MNO_VERBOSE=1` · History4 "Step 6" |
| **Na₂'s polarized singlet will not converge from an AFM seed at α=0.3** | Δρ parks at ~1e-2 and oscillates while E is flat to 1e-9, on DIIS/GDM/Kerker/Pulay/smearing alike; α=0.5 converges in 66; UNPOLARIZED converges in 30 at α=0.3 | "the two-channel run of a ζ=0 system is much harder to converge than its unpolarized twin" is worth explaining or fixing before a magnetic campaign pays for it · History4 "TWO SMALLER LOOSE ENDS" |
| **The ‖V_xc − V_xc_fit‖ fit-quality study** | every Becke-vs-uniform cost number is taken at unknown-equal accuracy, so "Becke is a negative acceleration" is a COST statement, not a verdict | ⏸ parked (user: *"defocusing"*); it collides with the Becke grid-size item below, so whoever does that will be standing next to it |

### 4b. Performance — the standing measurements

The head-to-head instrument is `doc/Benchmark.md` (COPY the run command from §5a).  Where it stands: **per SCF
iteration ahead of CP2K on seven of nine rows**; the like-for-like parity row is **1.13× on the comparable
fixed-point stage, 1.72× for the two-stage probe**, and its whole residual is ONE gather (lever B, behind N4);
lever C (GDM's trial densities) is unjudgeable until OT exists; **peak RAM: solved, we win** (113–132 MB on the
parity routes vs 217); **iteration counts: RE-JUDGED 2026-09-20 under Benchmark rule 3f — on the same measure, threshold and loop shape we take 22 / 24 (imposed) and 33 / 37 (free, parity) against CP2K's 44 / 104**; the old "DOCUMENT, do not chase" verdict compared different convergence measures (History4 "BIN 4").  Threaded (12 cores) we beat CP2K on this class of system
(4.14 s/step vs 5.90; they are 0.82–1.09× their own serial) — no parallel opportunity CP2K exploits is missing
on 4-atom cells; a 100-water box would invert it.

| row | what is open | next concrete action · record |
|---|---|---|
| ✅ **ctest ran NO Γ-named SCF test from 2026-09-15 to 2026-09-20** (found 2026-09-20 when a gate that failed by hand read "Passed 0.02 s" in ctest) | CMake 4.2 discovers through `--gtest_output=json`; googletest escapes non-ASCII BYTES as `\u00XX`; `Γ` (CE 93) came back as `Î` + U+0093 and the filter matched nothing — gtest exit 0, ctest "Passed".  48 tests, five days of green sweeps that never ran them — including the N3-promotion sweep, which is how the (ρ,m) singlet NaN (next row) got through | FIXED the same day: `DISCOVERY_EXTRA_ARGS --gtest_output=` (text-listing discovery; CLAUDE.md Tests).  The Γ spelling STAYS (user ruling 2026-09-15) — the tooling was wrong, not the name.  Habit: an SCF test cannot pass in 0.0x s · CLAUDE.md |
| ✅ **(ρ,m) Kerker on an exact SINGLET produced a NaN Fock matrix** (found 2026-09-20 by the +U polarized Kerker gate — the first Kerker singlet the promoted default met) | `RasterKerker`'s filter `g²/(g²+G0²)` with the m leaf's `G0=0` (linear, undamped) is 0/0 at G=0; the NaN rode the raw-raster shadow into the rebuilt channels and v_xc.  `KerkerStep` (G-space) had the guard; its raster twin did not | FIXED: the same `g2>0 ? … : 1.0` guard; gate `GPW_Si.Γ_Imp_Pol_Kerker_eqUnpol` (a no-U polarized Kerker singlet == the unpolarized run to 1e-6 under the default basis) · `FieldMixer.C` |
| **Linear D-mixing DIVERGES on free Si with no Fock accelerator** (found 2026-09-20 by `scripts/retake5a`) | `GPW_Si.Γ_CP2K` under `CP2K_COMPAT=1 GPW_ACC=null` (plain linear D-mixing, `ProductionGates` α=0.30): descends to −7.11505 by iteration 19, then diverges (E −7.0886 at 60, [F,D] 2e-3 → 0.2).  The adaptive relax raised α 0.30 → 0.45 at iteration 2 and never re-damped.  CP2K's `DIRECT_P_MIXING` at α=0.4 converges the same cell to 1e-7 in 12 steps; our Kerker G0=1 alone takes 54, Kerker + Pulay(8) 13 (E −7.115067447, the anchor).  DIIS on the Fock side had been masking it | two questions: (a) why does the adaptive controller (V1.18) not re-damp on a rising energy — a rising-E, rising-[F,D] trajectory is exactly its trigger; (b) is direct P mixing at fixed α=0.4 stable on our loop (a `GPW_ALPHA` knob + adaptive off)?  Cheap: seconds per run.  Not a physics question; until answered the no-Fock-side route on our side is Kerker + Pulay · `doc/Benchmark.md` rule 3f, `scripts/retake5a` |
| **Re-take Benchmark §5a on CP2K's convergence measure** (rule 3f, 2026-09-20) | every `q iters` / `c steps` / `total ×` cell in §5a compares a 1e-3–1e-5 mixer residual with CP2K's `EPS_SCF` on max\|ΔP\| — two different questions; the s/it columns survive, the counts and totals do not | re-run each §5a row with `Δρmeasure=MaxΔD`, `MinΔρ=EPS_SCF`, `NMaxIter=MAX_SCF` off its deck, the deck-shaped loop (Null accelerator, Pulay 8, no MOM) where the deck mixes Broyden/Pulay, and `CP2K_COMPAT=1`; one sitting, ~1 h serial (Si rows seconds, NaF minutes, MnO ×2 the bulk).  The two MnO rows are DONE (§5 ⁸); the Si Γ example is in rule 3f · `doc/Benchmark.md` rule 3f |
| **Size the Becke grid** (was item 1 / bin 2) | MnO's setup is 184 s = 47% of the default run (CP2K 8.1 s), 136.6 s of it TWO Becke mesh builds — but the RECIPE (nR=40, degree 29) is **over-generous on every system**: 3.5× on Si/Al, ~25× on NaF/Mn (`BeckeLadder`, scored by \f$\max|\Delta V_{xc}(i,j)|\f$ — score the MATRIX, the energy cancels error the operator does not) | ▶ **a POLICY CALL for the user**: at an absolute \f$\max|\Delta V_{xc}|\le\f$ 1e-4 the answer is nR=40, GL-17…21 (2–3× cheaper than production).  ⚠ no default flips on ladder evidence alone — Al is non-monotonic and a frozen ladder understates the self-consistent shift on a metal ⇒ a converged A/B on Al first.  Two items survive whatever the calibration says: (a) why TWO builds for one cell (each anneal stage builds one); (b) the build threads at 8.2× on MnO — NaF's serial fraction is still its setup (likeliest SIZE).  Then the per-element radial scaling, the coarse-end routing calibration, Becke+IBZ (the real-space star-average is untested on this route) · History4 "1. BIN 2" |
| **BM(3) — the stock Lebedev rule is already site-invariant** | Lebedev-29 (302 dirs) is EXACTLY invariant under Si's \f$T_d\f$ site group (0 unmatched of 7248) because Lebedev rules are octahedral orbits; W2b nonetheless builds a site-adapted rule at 886 dirs/atom — a **2.9× sitting unclaimed** | test the STOCK rule for site invariance first and reuse it when it passes (the test must stay EXACT); keep W2b as the fallback; measure on MnO before spending anything — it vanishes as site symmetry drops · History4 row BM |
| **The ρ̃-mixed sampling bucket** (the actual per-iteration XC lever) | `FourierMixCD.C:65` samples \f$\rho(r)=\sum_G\tilde\rho(G)e^{iGr}\f$ by DIRECT SUMMATION at every mesh point, every Kerker/Pulay iteration: 35.0 s / 6 iterations serial against the DM GEMM's 1.70 s.  ⚠ It is NOT the low-rank D bucket (that GEMM is nearly bypassed on ρ̃-mixed recipes — the 7–8× rank win is real but bites only on DM-backed routes) and NOT Φ-sparsity (Φ is 48% dense on MnO; per-atom batching = 1.05×) | attribute the cost INSIDE that bucket before optimising: an FFT to a coarse uniform grid + interpolation to the mesh (O(N log N + npts)) is legitimate only for the cusp-free CORRECTION, not the full ρ; an adaptive G-ball keyed to \f$|\tilde\delta|\f$ is cheap late and exact at convergence.  N4 decides what is sampled · History4 "Vxc MUST BE FED THE DM ρ(r)" |
| **The k-scaling gap** (cross-k gather memo) | Si per-iteration cost rises 6.9× from Γ to 8 k-points; CP2K's rises 1.03× — we start 6.7× ahead and spend it on k.  Cause: the per-offset reductions \f$B_{ij}(n)\f$ are k-INDEPENDENT but `IntegrateMemo` is bypassed whenever a density screen is passed (32 of 51 gathers on Si 2×2×2 are the SAME FIELD).  ⛔ two attempts refuted (History4): withholding the screen perturbs the SEED trajectory; flooring vanishing weights breaks the fold's orbit invariance | the trajectory-exact route: memoize \f$B\f$ over the UNION of active sets and let each block contract only its own; raise `kMaxIntegrateMemos` (4) with it.  Pin 19 already removed the gather's D-screen, which is the other half.  CP2K also folds 8 → 4 by TIME REVERSAL on the shifted mesh — a separate 2× · History4 "THE k-SCALING GAP" |
| **Step 2's remainder — the per-iteration G-space folds** | the {G}-star fold is wired at two STATIC sites; the per-iteration consumers (ρ̃, the Poisson multiply, the V_xc gathers, G_ERI3 columns, seed structure factors) are UNFOLDED: 12–24× on MnO's magnetic group, 48× cubic, unclaimed; **T3.4b** multi-k per-block arming of the pair-stream fold (Γ-only is armed by default) | extend the ball fold to the per-iteration sites (the FFT itself does not fold trivially); T3.4b = union-of-reps stream caches or the star-summed joint scatter.  Space-group collocation reduction route (b) IMPOSES the symmetry — say so · History4 "Step 2", `doc/Records/SymmetryUpgradePlan.md` |
| **\f$B_{ij}(R)\f$ k-independent 1E memo** | "keep k out of the key"; payoff only on multi-k | time with row KP (§2) · `GPWPlan1.md` item 5 |
| **Becke partition, the loop itself** | `-march=native` alone measured 1.13× with the loop still SCALAR; `BeckeImage` is array-of-structs with a data-dependent `P>0` exit that saves only 2.3% of work | vectorise first (chunk the exit, SoA); the `norm()` table (~2 MB gather vs a 15-cycle hardware sqrt) is DOUBTFUL and points the opposite way — re-measure after vectorising, do not assume.  Both are small beside V1.22 · History4 "Becke partition, what is LEFT" |
| **`Eval`/`EvalGradient` still run their own image loop** | the per-point callers (KB quadrature) were not migrated to the `LatticeSum1E` point-set seam because the seam re-derives its offset list per call; so the code's *"ONE remaining explicit image list"* is narrowed, not retired | cache the offsets on the seam side, then migrate; it retires the last place GPW enumerates images itself · History4 "THE Φ BLOCH POINT-SUM SEAM" |
| **Φ-sparsity for LARGE cells** | ⛔ refuted on MnO (48% dense) — the item's own caveat "the win grows with cell size" was the whole story | re-measure with `GPW_PHI_SPARSITY=1` on a battery-scale supercell BEFORE building anything |
| **The singles-vs-pairs census** | the factored density turns pairs into singles (8778 → 118 on MnO) but the comparable count is (i,R) singles vs (i,j,R) pairs at the same ε, weighted by box volume — pairs SCREEN far harder (Gaussian product theorem), which is why CP2K collocates pairs | build the singles-side census over the orbital reach before any singles-route work; the crossover is system-dependent and may favour BOTH routes · History4 "THE FACTORED FORM CHANGES THE OBJECT" |
| **Parity deviation #8** | `CP2K_COMPAT=1` is our best KNOWN parity, not proven; the deviation list grew 4 → 7 every time anyone looked | finding the next one is part of any parity claim; the like-for-like row is a `GPW_MNO_NMAX`-capped probe — step one is an uncapped or equal-cap comparison · `doc/Benchmark.md` §2, §5a |
| **Residual 1.33–1.45× CPU inflation at 12 threads** | load imbalance / memory bandwidth; ~4.1 s left in Phase 1, deliberately not taken | revisit on a 2×2×2 LiMn₂O₄ supercell where the buckets reshuffle — optimising the tail of a development cell is fitting to the wrong system |

---

## Parked — real work, deliberately not queued

- **Lever B** (one gather per spin) — behind N4.  **Lever C** (GDM's trial densities) — behind OT.
- **The GPW (orbital family × fit family) pairing is FROZEN in the engine** (`GPW_Evaluator` owns its PW
  density-fit grid; pin 2 says pairings are POLICY at the factory) — becomes a defect the day a second fit
  family wants the same orbital engine.
- **The basis-side nullable vendor** (retire `GetRhoOnGrid`'s empty-vector signalling) — needs a ruling
  (`CleanupCandidates.md` R1.0n/R1.0o).
- **A `MixingPolicy` derived from (ordering, cell, functional)** — "the user says AFM MnO, not `G0=1.0,
  Pulay=0`"; expert overrides stay behind a named policy so a bad combination becomes a deliberate act
  (N1/T4's shape; the detectors T1–T3 that make it safe are built).

## Deferred & descoped — recorded so they are not re-litigated

- Fold `QchemTester` + the pybind bridge onto the facade — test-harness/binding cleanup, not lib surface.
- Container utils to `src/` (`sample_scalar`/`sample_gradient`, `Structure::BoundingBox`) — binding convenience.
- `SCFParams` ASCII rename — DROPPED; solved by C++20 designated initializers (`34ccf302`).
- `MolecularSym_EC` → `FixedIrrepOcc_EC` rename — belongs with the queued symmetry-naming cleanup.
- The geometry hoist of the box walk (§5 of ScreeningPlan) — REFUTED on measurement (13.2% ceiling vs a 14.5%
  price); the exp recurrence — REJECTED (anisotropic, flips a degenerate basin); column-major Φ fill — a
  net loss; chasing CP2K's OMP — settled by `Benchmark.md` §7.  All in History4; do not re-propose.

## ✅ GDM declined every UNPOLARIZED run (closed 2026-10-01)
Cause: `UseFD`'s idempotency test (Tr D′ = Tr D′²) and block-occupancy test (Tr D′P = nocc) assumed g=1, and `Create`'s `occ` is the reservoir ELECTRON count (`GetN`), used as an orbital count. Fix: g = round(Tr D′²/Tr D′) read off D′ itself (no interface change); idempotent iff g·Tr D′ = Tr D′²; nocc = N/g refreshed each UseFD; block check against g·nocc. The gradient and the diagonal quadratic model both omit the common factor g, so step length is unchanged. Gate: `SCFAcceleratorGDM.AFoldedDoubletEngagesAndStepsLikeItsPolarizedTwin` (orbitals identical to 1e-10) and `M_Calculation.WaterPureGDM` (11 geodesic steps, DIIS anchor). Full ctest 983/983.

## ✅ Shifted-MP fold lowered Si by 1.02 mHa (closed 2026-10-01)
`GPW_Si.k222s_Imp_CP2K`: imposed −7.868473 vs FREE −7.867453. Cause: a shifted MP mesh on fcc is not Oh-invariant; `FoldGrid` skips ops that do not map the mesh onto itself (8 → 2 irreducible via the subgroup), but `DetectPointOps` handed ρ's G-space/raster star-average all 48 ops, so ρ was symmetrized over an Oh-orbit k-set bigger than the one the bands were summed over. Fix: filter `recipDensity`/`directDensity` through `Fold::MapsMeshOntoItself` (BasisSet.C). Result −7.867453736 = FREE to µHa; test tolerance tightened 3e-3 → 1e-4. Diagnostic (a)-(c) were not needed. Γ-centred and TRIM meshes keep all 48 ops (full ctest 983/983).

## ✅ valgen could not mint a fully-stripped ion (closed 2026-10-21)
`--electrons 0` was the "unset" sentinel (`r.electrons > 0 ? r.electrons : Zion`), so Sr²⁺ at q2 silently came back
NEUTRAL (charge 2.0003, found 2026-10-05 minting SrVO₃'s seeds).  Fix: `ValenceBasisRecipe::electrons` defaults to
−1 (= neutral, Zion); 0 is a fully-stripped ion and is REFUSED (`std::invalid_argument` from Generate{ValenceBasis,
SeedDensity,FloorEnergy}; valgen exits 1 with the reason).  Decision on "what is an empty-valence seed entry":
NONE is needed -- `SeedCD`'s species rule seeds `target<=0` as nothing (scale 0) without consulting the library, and
`MagneticDecoration` calls it non-magnetic; there is nothing to validate a basis in either (validate on the
neutral).  The library's Sr entry (Nelec 2, q 2) is the genuine neutral, untouched.  Gate:
`UTCalculation ValenceBasisGen.ElectronsUnsetIsNeutralAndZeroIsRefused`.

## ✅ l=4 (g) shells in a lattice basis (closed 2026-10-21)
TiO₂ + Ti g 1.2 failed three ways (vet "drops that map to no nameable shell"; CholeskyPivoted dropping ALL 18 Ti g;
`Hubbard_U: shells of one l with different component counts`).  ONE cause: `Math::SphericalShell(l)` tabulated
l=0..3 and `assert`ed above -- compiled out in Release, so l=4 returned `{}`.  `BuildCartToSphere` still counted
2l+1 = 9 columns per g shell but filled none (9 all-zero view functions per shell → the Cholesky drops), and the
view's `GetAoShells` built a 0-component `SphericalShellRep` (shifted offsets → unnameable drops, mismatched
counts).  Fix: `Math::GenerateSolidHarmonics(l)` (Rodrigues/Legendre expansion, any l, m=−l..l sin/cos) as the
`default:` of `SphericalShell`; the l≤3 tables are kept byte-identical and the generator reproduces them up to a
positive scale per m.  Gates (UTGaussian_BS): `M_Spherical.{g,h}_orthogonal`,
`generated_harmonics_match_tables_and_are_harmonic` (l≤6), `M_SphericalView.g_matches_native_family` (view ≡ native
spherical family, eigen/pencil/span-projector probes), `g_shell_keeps_nine_unit_components`, and the user's
"first check MnD vs libcint for Cartesian g": `M_LibCint.matrix_1E_matches_scalar_g` (S/T/V, LMax=4) PASSES -- the
MnD g integrals were never the problem.  End-to-end: the CP2K/QE-gated TiO₂ Γ ortho-atomic deck on
`valence_lowq_tm.bsd` + Ti g 1.2 (scratch) now vets clean (min eig S 8.3e-4, nothing trimmed), builds all twelve
+U manifolds, and converges in 35 iterations to −71.135687 Ha (−0.63 mHa vs no g).  The χ₀/U effect of g on b6 is a
REMAINDER, kept on the OpenWork row.

## ✅ LibXC-polarized (closed 2026-10-21)
`Libxc_LDA` was scalar-only and `Factory(XCFunctional)` threw for Polarized + LibXC.  Now the wrapper inits each id
twice (XC_UNPOLARIZED for the scalar face and exact ζ=0 -- byte-identical for every unpolarized run -- and
XC_POLARIZED for the two-channel face); a negative channel is clamped to 0 on its own (libxc's rule).  ⚠ The
per-channel SPLIT of ε matters, not just Σρ_σε^σ: the energy is Σ_σ Tr(D_σ ε^σ_fit), so with libxc's one-ε split for
Dirac exchange polarized N came out 0.146 mHa off the in-house DiracVWN although every pointwise value matched to
1e-9.  An exchange id (`info->kind==XC_EXCHANGE`) therefore uses the exact spin-scaling relation (the base class's
channel-separable face, SlaterExchange's split); correlation ids use libxc's ε for both channels (VWN5's split).
Result: polarized N LibXC SVWN5 vs DiracVWN 4.23 µHa, smaller than the same pair UNPOLARIZED (7.47 µHa, the
implementation floor) and unchanged under a 1000× tighter SCF gate.  Gates: `UTHamiltonian
LDA_XC.LibxcWrapperPolarizedFace` (raw libxc, VWN5 and Slater faces incl. ρ 1e-8 and ζ→1; the only disagreement is
v↓ at EXACTLY ζ=1, ~1e-5, zero-weight), `ITMain A_DFT_U_Consistency.PolarizedLibxcSVWN5MatchesInHouseDiracVWN`.

## ✅ Ladder restart starts on the DIIS+Kerker rung (closed 2026-10-21)
USER DECIDED (2026-10-21): `SolidCalculation::Restart` enters the Ladder at its GDM rung.  Mechanism:
`SolidAcceleratorOptions::ladderStartOnGDM` → `SCFAcceleratorLadder(..., start)` (announced on the console);
`Restart` sets it for ITS stage only and restores the caller's value afterwards, so later stages on the same
calculation start as asked.  A smeared restart cannot engage GDM (D′ not idempotent): the Ladder's existing RETREAT
rule drops to DIIS at the first projection step, loudly.  Gates: `GPW_NaF.Γ_Imp_eqExactResume` and
`Γ_Imp_eqColdStart` now run the PRODUCTION Ladder (their GDM override removed): 2 iterations each (were 9 and 12 on
the DIIS-first Ladder), exact resume dE −4e-10 Ha, first iterates within the 1e-7 / 1e-6 pins;
`UTCalculation SolidState.LadderRestart_Smeared_RetreatsToDIIS` (Si k222, kT 0.02: enters GDM, retreats, lands on
the cold smeared Ladder energy to 1e-6).  `rundeck`'s restart path goes through the same `Restart`.

## ✅ Response: REPORT orbital-resolved channels + inter-site V_IJ (closed 2026-10-21)
V_IJ had already landed (`a18ef3dd`: `QMeshResponse::Interaction`, `SolidCalculation::InterSiteHubbardV`, the `[V]` lines); the
row was stale on that half.  The resolved channels: the response channels can now be the +U term's own U SLOTS (site irrep x
grey parent, `HubbardChannels::SlotProjectors` -- so a resolved U names a slot `Uirrep` can carry) or the occupation
EIGEN-ORBITALS inside each slot (spin- and k-summed n, a multi-dimensional irrep's partners grouped: QE orbital_resolved /
Macke JCTC 20, 4824 (2024)).  A channel is amplitudes Xᵀℓ on an orthonormal basis X of its subspace, so nothing below the
probe changed.  Deck: `hubbardLinearResponse` `"resolve": "manifold"|"slot"|"eigen"`; `perturb` still names manifolds (all
their channels perturbed).  Printed: a `[response] channel` table (label, n, composition on the named projector functions);
U full (hp.x) / 1×1 (1/χ₀_cc − 1/χ_cc: one orbital alone in the q-mesh supercell, what a ±α FD on one orbital measures) /
orbital-diagonal (χ₀, χ off-diagonals between one manifold's channels zeroed, inter-site kept -- Macke's pointwise
inversion); the record carries `U1x1_eV`, `UorbitalDiagonal_eV`, `manifold`, V with labels.  Gate (UTResponse
`ResponsePolarizability.GPW_Si_dManifold_SlotChannelsSpanTheManifold`, Si k211 + one d shell): Td d = E {z2, x2-y2} + T2
{xy, yz, xz}; Σ_cc' χ_cc' over the slots = the manifold's χ (χ₀ and χ); eigen = slot when each irrep occurs once; 1×1 =
orbital-diagonal for one site.  Measured there: U(E) 61 eV at n = 0.003 vs U(T2) 7.4 and the shell-averaged 7.1 eV -- a
near-empty channel's resolved U is not a physical number.  TiO₂ Γ ortho-atomic smoke (scratch, `resolve: eigen`, Ti0 3d):
D₂h d = 2A_g + B1g + B2g + B3g resolves into 5 channels (the two A_g mix xy/z² on the global axes; n 0.25 / 0.63), O p into
z vs (x,y); full U 3.0–7.0 eV, 1×1 −0.6–2.2 eV (a 6-atom Γ cell: no isolation; the oracle needs a supercell/q-mesh).

## ✅ A carried density is re-filled at kT = 0 (closed 2026-10-21, anchor-moving batch item 1)
Cause (`SCFIterator::Initialize`): every seed's iteration-0 Fock-diagonalize-and-fill ran on the occupation policy's DEFAULT
(integer aufbau, kT = 0, no MOM) because the run's own policy is only assembled in `Iterate`; for a CARRIED converged density
(a restart, an anneal / ±α / U-loop stage) iteration 1 started from that kT = 0 re-fill.  Fix: the caller says what a seed IS
(`SCFIterator::SeedRole`, explicit-seed ctor): `Guess` (every strategy seed; unchanged) or `WorkingDensity` (set only by
`SolidCalculation::BuildStage` and `Restart`).  A carried seed is held to the first `Iterate`, and when the run's policy is
not the default (kT > 0 or MOM) its iteration-0 fill is REDONE under that policy -- a restored fixed point refills to itself.
At kT = 0 without MOM the path is bit-identical.  Two dead ends, recorded so they are not retried: (1) "a matrix-backed seed
is a carried one" is FALSE (the Uniform and plane-wave strategy seeds are matrix-backed, in other block layouts -- 11 tests
threw on mismatched composites); (2) making the carried density the working density OUTRIGHT loses a perturbed stage's full
first step: the ±α SCFs limit-cycled once the linear mixer grew relax 0.2 → 1.0.  Also: the ±α / restore stages of
`HubbardFiniteDifferenceChi` take a restart recipe (LINEAR mixing's relax capped at 0.3 -- a Kerker/Pulay α is left alone: the
k211 FD gate's ±α runs converge at 0.45 and not at 0.3; at kT = 0 a Ladder enters at its GDM rung), and
`Restart` enters the GDM rung at kT = 0 ONLY -- a smeared stage restored from an integer state hands GDM an idempotent D′,
which engages and never applies the smearing (Si k222 kT 0.02: 6 mHa off; gate `SolidState.LadderRestart_Smeared_StartsOnDIIS`).
Gate: `SolidState.SmearedExactResume_StartsAtTheSavedState` (Si k222, kT 0.02): without the fix the first iterate is 5e-6 Ha
off the saved state and the run needs > 3 iterations; with it, it starts on the state.  Sweep vs the batch baseline (HEAD
`69a58211`): 1051/1051, ONE test moved -- `GPW_Al.Γ_Imp_Anneal_Anchor` 48 → 28 iterations, energy unchanged to 1e-9; the FD
gate's χ is bit-identical (−14.41588683).
**CORRECTION 2026-10-22 -- the "orbital-diagonal" U is REMOVED** (bug report from the β-MnO₂ oracle session: TiO₂ Γ, b6 basis,
`resolve: eigen`, 22 channels -- every orbital-diagonal U −22 to −846 eV, while 1×1 sat at 1.3–4.3 eV).  Not mis-assembly: numpy on
the logged χ₀/χ reproduces every printed number.  The CONSTRUCTION was ill-posed (it was my reading of Macke's "pointwise", which
zeroes ALL off-diagonals = our 1×1): zeroing a site's inter-orbital off-diagonals while keeping the inter-site ones does not keep a
negative-definite response definite -- the screened χ's same-site off-diagonals are large, and the modified χ had eigenvalues +0.44,
+0.55, so its inverse's diagonal went large and wrong-signed, worst on χ's near-null channels (exactly where the full inversion is
largest).  Output, record field `UorbitalDiagonal_eV` and `QMeshResponse::Interaction/U(group)` are gone; full and 1×1 stay (1×1 on
a q-mesh supercell is the isolated single-orbital number QE's orbital-resolved ±α FD measures).  The Si d gate now asserts the
resolved χ₀ and χ are negative definite.  If inter-site screening in a resolved U is wanted: a PRINCIPAL submatrix (one orbital
type on every site), which keeps definiteness.

## ✅ N4 — the cusp-deficit XC feed, made the DEFAULT (closed 2026-10-22; pin 28)

The row as it stood before closing:

| **N4 — cusp-deficit XC feed** (pin 18; REQUIRED by pin 28, user 2026-10-22) | ρ_XC = ρ̃_mix + (ρ[D] − ρ̃_fit) (ρ̃_fit = the fit's representation of ρ[D], any fit basis) is how a field mixer meets pin 28.  WIRED 2026-10-22 on Kerker AND Pulay, (ρ↑,ρ↓) and (ρ,m) (`PulayMixer`: correction against the LATEST out, not the extrapolated one; (ρ,m): per spin, formed by `PolarizedDensityMixer`; UT `KerkerMix.CuspDeficitIsPerSpinOnThePulayArmInEitherChannelBasis`); the PAIR XC sampler still refuses it.  `SCFParams::XCCuspDeficit` stays OFF.  MEASURED (Pulay + Null): O₂ box −30.70246 → **−30.71018173** (= ρ[D]; 24 → 25 it, both channel bases); Si k222 unchanged (12 it, 10 digits); free MnO AFM-II (`decks/MnO_AFM2_free.json`) 29 → 28 it, −61.41608505 → −61.41613923, order SURVIVED (4.455 → 4.453).  (a) ρ_XC ≥ 0? NOT mid-iteration (O₂ ≤ 13.5 % of points, min −4.5e-4; MnO ≤ 25 %, min −0.038) but → 0 at convergence (MnO 0 points, O₂ 1e-7 e), against the fitted ρ̃'s 20 % / −4.6 e (O₂) and 8.6 % at MnO's convergence.  (b) iteration-static: not measured | NEXT: make `XCCuspDeficit` the DEFAULT (pin 28) -- anchor-moving: sweep, re-judge every moved Kerker/Pulay anchor against its oracle (pin 10; the `_CP2K` Kerker cells compare to CP2K's grid-collocated ρ, so read them first), then the campaign decks.  Then lever B (needs this) · H4 "N4" |

**Closed 2026-10-22 (user: "make XCCuspDeficit the default (pin 28)").**  `SCFParams::XCCuspDeficit` defaults to TRUE; false is the
pre-pin-28 control (`GPW_Si.Γ_Becke_Imp_Kerker_FitXC_CP2K`, the former plain cell; the former `..._CuspDeficit_CP2K` is now the
default-named `Γ_Becke_Imp_Kerker_CP2K` -- grammar token `CuspDeficit` retired, `FitXC` added).  The pair (raster) XC route
needed no correction map: a Kerker/Pulay-mixed density answers RAW with its raster SHADOW -- ρ[D] collocated pointwise and mixed
on the raster -- which is ρ[D] at every quadrature point at the fixed point; that route now accepts the flag and refuses it only on
the BALL (fitted-ρ) fallback.  Hence the sweep moved only singles-route (Becke) Kerker runs.

Sweep vs the 2026-10-22 baseline (`ctest -j4` + `-L long`, `scripts/anchordiff --tol 1e-10`): 1061/1061 + 4/4 pass; 7 movers,
every one re-judged against its oracle, no verdict changed:

| test | before → after (Ha) | oracle | gap to oracle |
|---|---|---|---|
| GPW_Si.Γ_Becke_Imp_Kerker_CP2K | −7.1149381244 → −7.1149381646 | CP2K −7.11505788 | +0.12 mHa (unchanged) |
| GPW_NaF.Γ_Becke_Imp_eqUni | −24.4312534648 → −24.4312534677 | CP2K SR2 −24.4312134 | −0.04 mHa |
| GPW_NaF.Γ_Imp_eqColdStart / Γ_Imp_Anchor | ≤ 1.8e-9 | same | +0.81 mHa |
| GPW_NaF.Γ_Imp_eqExactResume | 0 (39 → 36 it) | — | — |
| GPW_MnO.Γ_Shub_Pol_Smear_Anchor_Long | −61.4145468 → −61.4145709 | CP2K AFM-II −61.470570 | 56.02 → 56.00 mHa; moment 4.4506 → 4.4490 |
| GPW_MnO.Γ_U_Shub_Pol_Smear_CP2K_Long | dE 0.6001933 → 0.6001958; E_U 0.5873195 → 0.5869510 | CP2K dE 0.6173543, E_U 0.6095092 | dE −17.16 mHa (unchanged); E_U 22.19 → 22.56 mHa (gate 0.05) |

The renamed Kerker pair reproduces the old cells to all digits.  Measured before the flip (Pulay + Null): O₂ box −30.70246 →
−30.71018173 (= ρ[D]); free MnO AFM-II −61.41608505 → −61.41613923, order SURVIVED; ρ_XC not pointwise ≥ 0 mid-iteration but → 0
at convergence.

## ✅ MnO accuracy — name the operator (closed 2026-10-23)

The row as it stood when closed (the full hunt, rung by rung):

| **MnO accuracy — name the operator** | VA exact-span: −99.7 mHa configuration-BLIND offset + −37 mHa FM-favouring selective bias vs CP2K.  **RE-MEASURED 2026-10-23 under pin 28** (decks + logs `~/Code/qchem6-runs/MnO/va_pin28_*`, VA spherical, the +U test's CP2K-shaped loop): (1) the Mn ATOM in a 16-bohr box (VA, periodic, `IntegrationTests/CP2K/mn_box16_gpw_va.inp`) agrees with CP2K GPW to **−0.71 mHa** (E_xc +0.10 mHa) -- the blind offset is NOT a one-Mn property; (2) the two codes' MnO arms were NOT the same STATES: with the moment free (common μ) our FM settles at moment 8.24 (CP2K 9.998) and our AFM carries −TS −15.4 mHa (CP2K −1.0) -- the old "selective" number compares different states; (3) FM with the moment FIXED at 10 (multiplicity 11, no shared μ) matches CP2K's state (moment 10.0008, −TS −7.66 vs −7.06 mHa) and E_xc to +14.2 mHa, and the offset stays **−143.9 mHa** (A −61.44864 vs −61.30478) -- so ~−158 mHa sits in kinetic + PP + Hartree; (4) Becke vs uniform XC agree to 0.9 mHa on FM (the Becke tail drop is not it); our moment-8.24 FM is 13.4 mHa BELOW the high-spin one in our code, a state CP2K never reached.  Also seen: Kerker/Pulay + Null did NOT converge the AFM arm on the UNIFORM XC route in 200 iterations **BOXES 2026-10-23** (16-bohr, periodic, same basis/PP, `IntegrationTests/CP2K/{o,mno}_box16_gpw_v*.inp`, materials `O_box16` / `MnO_box16`): the O ATOM (triplet, moment fixed) agrees to **+0.14 mHa** at 144 Ha and −0.99 mHa at 72 Ha (E_xc ≤ 1.1 mHa) once our density cutoff matches CP2K's -- the first +26 mHa was our AUTO 36 Ha (row below); the Mn–O DIMER at the crystal bond length (moment fixed at 5, kT = 0): CP2K OT −30.34790; our only converging recipe (Kerker/Pulay fixed 0.3 + Null) was at −30.3491 when stopped at iteration 107, i.e. within ~1 mHa, while GDM stalls at −30.34225 and linear/DIIS 2-cycles (rows below).  ⇒ the −144 mHa needs the CRYSTAL: Mn, O and one Mn–O pair are each clean.  (CP2K `DETAILED_ENERGY` is an energy-decomposition option, not a term print)   **SUBLATTICES 2026-10-23** (the crystal's cell with one species removed, FM moment fixed, 72 Ha both, materials `MnO_AFM2_Mnsub` / `MnO_AFM2_Osub`, decks `IntegrationTests/CP2K/{mn,o}sub_gpw_v*.inp`): Mn₂ **+0.13 mHa** (E_xc +0.38); O₂ ≈ −2.6…−3.4 mHa (CP2K did not converge in 200 Broyden steps, wandering ±0.4 mHa; ours 23 iterations).  ⇒ every one-species / two-atom rung is within ~3 mHa: the offset appears ONLY with Mn AND O in the periodic cell.  Lead: the crystal's collocation charge leak ρ_lost/N 5.1e-5 (1.3 m e) is ~15–20× the boxes'.  EXONERATED 2026-10-23 (user's suspect): the D-aware collocation screen (`policy.dAwareScreen=false` = CP2K's geometry-only screen) moves the crystal by −0.27 µHa (FM m10) / −0.53 µHa (AFM), leak unchanged.  **GRID LADDER 2026-10-23** (FM moment 10 FIXED in both codes, `~/Code/cp2k-runs/MnO/mno_fm_gpw_va_m10_c{144,288}.out`): offset **−143.86 mHa at 72 Ha, −144.13 mHa at 144 Ha** (CP2K −61.30477589 / −61.30438787, ours −61.44863782 / −61.44851959; E_xc +14.4 / +13.9 mHa) -- grid-CONVERGED in both codes, and our leak 6.22e-5 is identical at both cutoffs (not a resolution effect).  E_alphaZ is consistent across the rungs: α_Mn/Ω and α_O/Ω from the sublattices reproduce the crystal's 3.17254 Ha to 2e-5 (a cross-species slip would be Ha-sized).  ⇒ a converged operator/convention difference that needs BOTH species in the periodic cell   **LATTICE SCAN 2026-10-23** (FM moment 10, 72 Ha, materials `MnO_AFM2_a{9p5,11,13}`, decks `IntegrationTests/CP2K/mno_fm_m10_a*.inp`): offset −143.9 (a 8.40, Mn–O 4.20) · ≈ −40.5 (a 11, ours −61.0389 not fully converged, CP2K −60.99837) · **−30.35 mHa** (a 13, ours −60.82359711 vs −60.79325075, same state: −TS −27.37/−27.25 mHa, Mn 5.15/5.13; E_xc +1.4 mHa); a 9.5 CP2K did not converge (±0.05 Ha).  Neither c + k/Ω nor c + B/d fits (both predict ≈ −61…−73 at 8.40); c + A·e^{−κd} through the three: **κ ≈ 1.7 / bohr, asymptote ≈ −28 mHa** ⇒ TWO components: a SHORT-RANGE one (≈ −116 mHa at the crystal, overlap range) that a single Mn–O pair does NOT show (dimer clean) -- so it needs TWO O on one Mn: a mixed-species THREE-centre term, ⟨χ_O₁|V_Mn|χ_O₂⟩, i.e. the local-PP short-range lattice sum (the separable KB is built from the dimer-tested two-centre overlaps; cf. CleanCode D13, single-species HGH_LocalPotential ignoring its Z) -- and a LONG-RANGE tail ≈ −28 mHa (CP2K Mulliken ±0.27 e charge transfer persists at a 13)   **FOUND + FIXED 2026-10-23 (user: "unit tests directly on the short range part of V would be most efficient")**: Mn q7 has NO short-range local PP, so the live 3-centre short term is O's well between functions on other atoms.  New gates `GPW.DiffuseDVshortMixedSpeciesOracle` (raster oracle from CP2K's GTH parameters; Mn-Mn-on-O / Mn-O / O-O classes all ≤ 1e-5 after the fix) and `GPW.VshortIsTranslationInvariant` showed the two translation-equivalent Mn of MnO AFM-II with V_short blocks **2.0e-3 Ha apart**: `BuildOpTerms` (PG_Cart_MnD) kept each operator image by |R_off| ≤ reach, a sphere about the OPERATOR's home position, while the lattice sum fixes χ_i in the home cell -- so each χ_i lost the O-well images on the far side of its neighbours, differently per site.  Fixed: keep an image iff it is within reach of SOME home basis centre.  Effect on VA (exact span, 72 Ha, both arms now land on CP2K's state by themselves): **AFM −61.35040091 / FM(m10) −61.35183415 → offsets −47.08 / −47.06 mHa (were −108.0 / −143.9); Δ(AFM−FM) +1.43 vs CP2K +1.45 mHa -- the −37 mHa configuration-SELECTIVE bias is GONE**; AFM −TS −0.976 vs −0.977 mHa, E_xc within 1 mHa.  Sweep: SR anchor −61.41455 → −61.35052 (re-judged, MnO.C header), +U CP2K test passes (moved +78 mHa on the summed arms), Si multi-atom −4.5e-5, NaF ~1.6e-7; NaF coarse→fine first iterate 4e-8 → 2.3e-6 (not convergence slop; WHY is open).  Remaining: a configuration-BLIND −47 mHa   **AFTER THE FIX, 2026-10-23**: lattice scan (FM m10, 72 Ha) −47.1 (a 8.40) · ≈ −36.4 (a 11, ours oscillating ±0.1 mHa about −61.0347) · **−37.43 mHa** (a 13, ours −60.83068209) -- the remainder does NOT decay: ≈ −37 mHa FLAT from Mn–O 5.5 to 6.5 bohr (neither overlap- nor Madelung-like) plus ≈ −10 mHa at the crystal spacing.  Ruled out since: a non-magnetic O (moment fixed 0, symmetry imposed: −0.95 mHa, −TS identical −19.0954 mHa -- the expanded cells force O ~non-magnetic); every ONE-ELECTRON lattice sum is equivariant under the Mn1↔Mn2/O1↔O2 translation (`GPW.LatticeSumsAreTranslationEquivariant`: S, T, V_short, KB ≤ 5e-14; V_long 8.7e-7, grid class) | NEXT: the flat ≈ −37 mHa needs BOTH species in a periodic cell at any spacing while every isolated piece in every state is clean -- suspects: the mixed-species electrostatic conventions (point-charge Ewald + erf V_long + E_alphaZ against CP2K's Gaussian core charges, self and overlap energies) and the many-electron terms on the collocated density; a = 16 / 20 points to confirm it is truly constant (CP2K may need OT there); CHARGED ATOMS RULED OUT 2026-10-23 (CP2K `CHARGE` is integer, so ±1 = ~3.7× the crystal's 0.27 e; 16-bohr box, moment fixed, `IntegrationTests/CP2K/{mnp,om}_box16_gpw_v*.inp`): O⁻ −0.91 mHa (neutral O −0.99), Mn⁺ −4.90 mHa of which −4.40 is the charged-cell E_alphaZ convention (row below; cancels in any neutral cell) -- corrected −0.50, as neutral Mn's −0.71.  The offset does not follow the charge; re-run the FM-below-AFM gate (`DISABLED_Γ_Shub_Pol_Smear_Ordering_Long`) and every campaign anchor built on MnO/NiO/oxide energies (they all carried this defect); why NaF's coarse/fine fixed points separated.  Earlier plan: pin `QCHEM_XC_DM_SOURCE` first (DONE: pin 28); then term-by-term vs CP2K blocks on the Mn ATOM (oracle −14.2440 / −14.658), then crystal, then `MNO_KMESH=2`; read weak-moment conclusions against the INTEGRATED moment (pin 4) · `Records/SphericalLatticePlan.md`, H4 "Step 5" |

**Closed by the Ewald fix (2026-10-23).**  An independent Python Ewald (η-independent to 8 digits, rocksalt Madelung 1.747565 exact)
gave −57.2917665 for the MnO AFM-II cell's Zion 7,7,6,6 point ions against our Enn −57.3390044: −47.24 mHa, scaling exactly as 1/a
(Δ·a = −0.3968 at a 8.40 and 13).  `EwaldEnergy` enumerated real-space images |R| ≤ rcut instead of the PAIR distance |τa−τb+R| ≤ rcut;
the O at (¾,¾,¾) sits 21.8 bohr from the Mn at the origin against rcut 21.9 (sublattices, with 14.5 bohr home separations, lost only
0.01 mHa).  Fixed (enumerate to rcut + max home-pair separation, keep r ≤ rcut); gate `StructureTests.EwaldFarApartHomeAtoms`.
VA after both fixes: AFM −61.30316303 / FM(m10) −61.30459627 vs CP2K −61.30332518 / −61.30477589 → +0.16 / +0.18 mHa.  +U test: E(0)
+0.16, E(U) +0.27, dE +0.11, E_U +0.009 mHa.  SR anchor −61.41455 → −61.35052 (V_short) → −61.30328 (Ewald), 0.05 mHa from CP2K VA.
Charged-cell aside (its own row): E_alphaZ uses N_e, CP2K's equivalent ΣZ.

## ✅ A Γ supercell breaks its OWN translation symmetry (closed 2026-10-23)

| **A Γ supercell breaks its OWN translation symmetry (~2e-4)** | free Si 2×1×1 at Γ, uniform XC, converged to Δρ 1e-7: χ₀ of sites 1 and 3 (translation images) differ by 1.7e-4 relative, U by 6.6e-4 eV, while the q-mesh route is translation-symmetric by construction.  It is the FLOOR of the supercell-equivalence gate (5e-4) | find which ground-state ingredient is not translation-covariant (per-site occupations first; then `dAwareScreen off`, the raster's divisions vs 2× the primitive's) · `Records/LinearResponsePlan.md` §3e step 4 |

**Closed:** the V_short operator-image defect (pin 29).  Supercell gate after the fix: χ₀R rel 9.25e-6, χR 1.16e-5; U per site 2.079615 /
2.079605 / 2.079615 / 2.079605 eV (q-mesh 2.079559 / 2.079520).  Tolerances tightened to 5e-5 relative and 3e-4 eV.
Left unconverged and stopped (not needed once the crystal agreed): the linear O–Mn–O and Mn–O–Mn triatomics in a 16-bohr box (materials
`OMnO_box16` / `MnOMn_box16`, CP2K OT decks `IntegrationTests/CP2K/{omno,mnomn}_box16_ot.inp`: CP2K −46.13211978 / −45.03771777).  Ours,
with the V_short fix but before the Ewald fix, Pulay 0.2 + Null at kT = 0: O–Mn–O −46.1333 (Δρ 5e-3) after 99 iterations, Mn–O–Mn
−44.979 (Δρ 5e-3, clearly another state) -- box SCFs with a near-degenerate open d frontier converge badly here (OpenWork §3 GDM row).

## ✅ +U projector independent of the orbital basis (OpenWork §2; closed 2026-10-25) -- the row as it stood

| **+U projector independent of the orbital basis** (user 2026-10-24) | increments 1–3 of `doc/HubbardProjectorPlan.md` LANDED 2026-10-24 (opt-in deck `"projector": "pseudoAtomic"`; the O 2p projector = QE's `PP_CHI` to 1e-8).  **Increment 4 LANDED 2026-10-25** (`b4cfc789`: ACBN0 over the pseudo-atomic set via its own exact one-centre integrals; TiO₂ b6 vs hp.x U(Ti) +4.2 % → +2.3 %, occupations = QE's; NiO cell +11 meV).  NEXT: increment 5 (name `siteShells` on the CP2K-gated manifolds first, then flip the default).  The gate FAILED and refuted its premise: the β-MnO₂ Γ gap does not move with the projector (VA 0.128 → 0.128 eV, VA + O p 0.18 0.237 → 0.234; QE 0.325) — the gap is the ORBITAL basis.  What the projector sets is E_U: VA + O p 0.18 pseudo-atomic 0.02214 Ha vs QE 0.02216 (site-shell 0.01767) | **USER 2026-10-24: keep the flexibility; look elsewhere for the gap.**  Increment 4 is done (above); increment 5 (make it the default, retire `radial`/`atomicRadial`, re-judge +U anchors) is in progress.  The β-MnO₂ gap discrepancy (0.231 vs 0.325 eV on VA + O p 0.18) is NOT the projector, the basis (flat under every addition), convergence or pin 28: at U = 0 every band is QE's to 0.04 eV, and under U only the U'd O-normal band (falls 0.09 eV less) and the Mn d LUMO (falls 0.13 eV more) respond differently — the self-consistent response to U (HubbardUHistory 2026-10-24 "gap residual") **RESOLVED 2026-10-25 by the plane-wave route: on QE's PP and a PW basis the U response is QE's to 1–2 % (Mn spin response +0.0385 vs +0.0391, Mn d band −0.769 vs −0.773 eV), so the Gaussian residuals (+0.023, −0.852) are the GAUSSIAN ORBITAL BASIS** — next: a Mn-d-response-aimed Gaussian basis judged against the PW U − U0 table (HubbardUHistory 2026-10-25 "THE PLANE-WAVE VERDICT") · record `Records/HubbardUHistory.md` 2026-10-24 (pseudo-atomic projector) |

Closed by increment 5 (`095f0bec` + `94ed536f`, both tiers 0 anchors moved vs the pre-flip sweep).  Plan text: `Records/HubbardUHistory.md` "doc/HubbardProjectorPlan.md -- retired 2026-10-25".

## ✅ The magnetic decoration is DERIVED FROM THE SEED (OpenWork §3; +U half closed 2026-10-26) -- the row as it stood

| **The magnetic decoration is DERIVED FROM THE SEED** (2026-10-23) | With `siteSpins` empty, `SolidCalculation` labels each atom magnetic iff the atomic library holds a spin pair at the SEED's per-species target (`MagneticDecoration`, `src/Calculation/Imp/SolidCalculation.C` ~348 and the Hubbard site-group path ~447).  So the seed strategy picks the imposed ordering: MnO AFM-II with `SeedStrategy::SAD` labels neutral O (2p⁴) UP on both sites → decoration Mn↑ Mn↓ O↑ O↑, a FERRIMAGNET → 12 Shubnikov ops instead of 24 → the run lands 0.78 mHa lower with Mn +4.4554/−4.4418, O −0.0068 (IonicSAD: O²⁻ non-magnetic, 24 ops, Mn ±4.4506, O exactly 0).  User ruling: the ordering is chemistry, stated like the Hubbard subspace, never inferred; a seed change must never change the imposed symmetry | make `siteSpins` MANDATORY for a polarized run that imposes or carries a Hubbard manifold (throw in `SolidCalculation`, same check in `deck::Vet` so `scripts/vet-decks` catches decks); retire the auto path.  **36 decks** (polarized + Hubbard, 13 also imposed: NiO, MnO, MnO₂, LiMn₂O₄, FeS₂, TiO₂) state none: judge each (a non-IonicSAD seed may have run on a mislabelled decoration — its recorded number may move), don't copy the derived one.  The SR MnO anchor already states `{+1,-1,0,0}` · `src/Calculation/Imp/SolidCalculation.C` |

**Closed for +U (2026-10-26).**  `RequireStatedMagneticOrdering` (exported beside `SolidCalcOptions`) throws for a polarized run with a
Hubbard manifold and no `siteSpins` (`greyImposition` exempt: it is the stated "no decoration" control), and for a `siteSpins` of the
wrong length or with an entry outside {−1,0,+1}.  The facade calls it first in `BuildBasis` (so `Restart` too); `deck::Resolve` calls it
so `scripts/vet-decks` refuses the deck and `Run` refuses before claiming a revision.  The +U site-group path no longer derives from the
seed (`BuildRun`).  UT `Deck.APolarizedPlusUDeckMustStateItsMagneticOrdering`.

**The 36 decks, judged.**  All 36 ran `seed: IonicSAD`.  What IonicSAD derived, measured per structure on the current code (and what neutral
SAD would have derived -- the defect's other arm):

| structure | decks | IonicSAD derived | SAD would derive | stated now |
|---|---|---|---|---|
| NiO_AFM2 | 22 (13 imposed) | Ni +1 −1, O 0 0 | O +1 +1 (ferrimagnet) | `[1,-1,0,0]` = derived |
| MnO_AFM2 | 4 | Mn +1 −1, O 0 0 | O +1 +1 | `[1,-1,0,0]` = derived |
| MnO2_lambda_spinel | 4 | Mn +1 ×4, O 0 | all +1 | `[1,1,1,1,0×8]` = derived |
| LiMn2O4_spinel | 4 | **all 0** (the per-species Mn target has no library spin pair) | all +1 | `[0,0,1,1,1,1,0×8]` (materials.json's deliberate FM gate-1 order) |
| FeS2_pyrite | 1 | all 0 | all 0 | all 0 (low-spin d⁶) |
| TiO2_rutile | 1 | all 0 | O +1 | all 0 (d⁰) |

No recorded number moves: the stated decoration equals the derived one everywhere except LiMn₂O₄, and those four decks are unimposed, where
the decoration only selects the +U site groups -- the σ=None stabiliser of a uniform (FM) decoration is the grey one, as for all zeros.  The
LiMn₂O₄ result is evidence for §3 "IonicSAD targets are per SPECIES" / "The polarized IonicSAD seed loses its spin".  Code-side callers
updated: `GPW_MnO.Γ_U_Shub_Pol_Smear_CP2K_Long` (`{+1,-1,0,0}`), three `GPW_Si` pol-singlet +U tests and `UTResponse` Kernel (`{0,0}`).

**REMAINDER (kept in the live row):** an imposed polarized run WITHOUT +U still derives its decoration from the seed.

## ✅ IonicSAD targets are per SPECIES, not per SITE (OpenWork §3; closed 2026-10-26) -- the row as it stood

| **IonicSAD targets are per SPECIES, not per SITE** (2026-10-22) | `IonicSADTargets` returns `std::map<Z,int>` and the loop overwrites it atom by atom, so every atom of an element gets the LAST one's target; `MagneticDecoration` and `SeedCD`/`PolarizedSeedCD` consume that map.  Harmless once equivalent atoms get equal charges (the tie-break above), WRONG for inequivalent atoms of one element: a charge-ordered cell (§3's LiMn₂O₄, Mn³⁺/Mn⁴⁺; any mixed-valence cathode) collapses to one Mn charge state.  MEASURED 2026-10-26: on LiMn₂O₄ the per-species Mn target has NO library spin pair, so IonicSAD's derived decoration was ALL ZEROS (Mn non-magnetic; λ-MnO₂ gets Mn +1) -- the spinel's seed spin is suspect too (row "The polarized IonicSAD seed loses its spin") | per-SITE targets (`std::vector<int>` parallel to the structure's atoms) through `IonicSADTargets` → `MagneticDecoration` → the seed ctors; a UT on a two-valence cell.  Before the spinel campaign (HubbardUPlan S5/S6) · `src/ChargeDensity/Imp/{Seed,SeedCD}.C` |

**Closed 2026-10-26.**  `IonicSADTargets` returns `std::vector<int>`, one target per atom (atom-index order), and the targets sum to
the cell's charge.  `MagneticDecoration(st, functional, ionicNvalBySite)` reads each site's own target.  `SeedCD` builds one RADIAL
CLASS per (Z, target, flipped-in-this-channel) and partitions the G assembly into groups: group k holds each species' k-th class,
ranked unflipped first and then by target.  A single-valence cell therefore makes exactly the old calls: one full-structure sweep, or
the unflipped/flipped pair of a channel seed.  UTs `IonicSADSeed.*` (UTChargeDensity) on an Mn₃O₄-like box (Mn³⁺ Mn³⁺ Mn²⁺): the targets
conserve charge site by site; the decoration reads each site's target; and the seed's ρ̃(G), total and polarized, equals the SUM of two
seeds that each hold one charge state (1e-12 relative).  The seed charge matches the site targets only to 0.016 e (45.016 for 45: the
library radials, the diffuse O²⁻ entry); the per-species map was 2 e heavy.

**What moved.**  Full `ctest -LE long` sweep against the 0f353ff5 baseline: `scripts/anchordiff --tol 1e-10` reports **0 of 1087 tests
moved** (energy and iteration count).  An audit of every material in `materials.json` finds only ONE whose per-site targets differ from
the old per-species ones: **LiMn₂O₄_spinel**, Mn 3,3,4,4 (2 Mn⁴⁺ + 2 Mn³⁺) where every Mn used to get 4.  Its 4 decks
(`campaigns/LiMn2O4/*_U{0,2,4}`, `probes/LiMn2O4_spinel_gate1`) now run on a charge-ordered seed, so their banked numbers are stale.
(Al_fcc and CsI_cscl have no library entry and were skipped by the audit.)

**REMAINDER (kept in the live row):** the electronegativity tie-break puts the two Mn charge states on symmetry-equivalent 16d sites by
atom index, so the charge order is CHOSEN rather than derived; state it the way `siteSpins` is stated if the spinel campaign needs a
specific one.

## ✅ Pulay without Kerker mixes D + the default mixer (OpenWork §3; closed 2026-10-26) -- the rows as they stood

| **Pulay without Kerker mixes the FITTED ρ̃; it should mix D** (pin 28, 2026-10-22; corrected 2026-10-24) | `PulayMixer` keeps its history as fitted ρ̃ even at `kerkerG0 = 0`.  ⚠ XC does NOT see the bare fitted ρ there any more: since N4 became the default (2026-10-22) it sees the cusp-deficit ρ, which IS ρ[D] at the fixed point -- measured 2026-10-24 on Si Γ (δ route): pure Pulay and Kerker-Pulay converge to the same −7.1150587901 Ha (on the raster route pure Pulay already mixes the collocated ρ[D], `shadow`).  What is still wrong is the FORM: the cusp deficit is pin 28's excused route only beside a G-dependent filter, and with none the step needs no {G} at all; the trace marks it `cusp*` (the `XC ρ` column).  A trajectory / convergence question, not a converged-number one: no campaign deck uses G0 = 0 (all NiO / SrVO₃ / CrossCode decks run Kerker-Pulay, G0 1.0).  ρ[D] is linear in D (ρ[Σcᵢ Dᵢ] = Σcᵢ ρ[Dᵢ]), so with no filter the whole Pulay history can live on density matrices and XC sees ρ[D_mix] exactly, with no cusp deficit | §3 ORDER step 2, WITH the Kerker/Pulay default pair (row below; one re-bank window).  ✅ its density face LANDED 2026-10-26 (`6833216b`): `tMixableDensity::Axpy` is the one binary mutation, `MixIn` = `ReScale(1-c)` + `Axpy(c)` with exact no-op shortcuts (c = 0 / factor 1); energies within 1e-10 both tiers, three iteration counts moved by 3 (GPW_SiBox ×2, GPW_Al: the old mix was FMA-fused).  NEXT: a D-space Pulay mixer (history of D_in/D_out, residual inner product on the density -- the run's quadrature, as D-MINDRHO's gate) chosen by the factory when G0 = 0; gate: O₂ box at G0 = 0 equals the linear / ρ[D] energy (−30.71018) · `src/ChargeDensity/Internal/PulayMixer.C` |

| **Kerker/Pulay as the DEFAULT mixer** | USER 2026-10-04: *"in general Kerker/Pulay should be the default mixer."*  Today `SCFParams` defaults to linear D-mixing (relax 1.0, no G-space preconditioner) + DIIS; measured the same day, free Si k211 at those defaults ran into its 200-iteration cap while Kerker/Pulay (`pulayDepth 8, pulayStart 5, kerkerG0 1.0, startingRelaxRo 0.45`) converged in 24 (also LinearResponsePlan §3e finding 5: restarts).  But the SAME day free Si at Γ went the other way: linear+DIIS 40 iterations, Kerker/Pulay 68 (pulayStart 3, relax 0.3), 72 (relax 0.2), never at relax 0.45 -- so the default needs its relax/start chosen, not just the switch flipped | PLAN (2026-10-22): an `SCFParams` mixing mode {Auto, Linear, Kerker, Pulay}, default Auto -- the PERIODIC iterator resolves Auto to Kerker/Pulay, molecules to linear (`StartingRelaxRo`'s 1.0 is LINEAR's, so Auto carries its own recipe relax); choose relax/start from a sweep of a few recipes over the GPW anchors with `scripts/anchordiff`.  ⚠ since D-MINDRHO the Kerker gate is HONEST (it used to stop 50–120× early), so Kerker iteration counts measured before 2026-10-22 are not comparable · `src/SCFParams`, `Records/LinearResponsePlan.md` §3e step 4 **MEASURED 2026-10-22 (R0 sweep: `Auto` → Kerker/Pulay 0.45/G0 1/depth 8/start 5, the accelerator UNCHANGED): a regression.**  19 tests fail, test time 1182 → 3025 s: restarts (`SolidState.*`) 8–10 → 80+ iterations, the Response FD SCFs (Si Γ, DIIS) DIVERGE (200 cap, −6.465 vs −7.115 Ha), O₂ box 19 → 60, PW k222 throws (its test builds the iterator WITHOUT the Structure, so no cell / reciprocal lattice -- not a PW limitation; PW Γ ran Pulay): Auto must ask `CanMixInGSpace`.  CAUSE: Pulay on ρ AND DIIS on F extrapolate at once.  Probes (rundeck, same E in every converged arm): Si k222 free linear+DIIS 18 / Pulay+DIIS cap / **Pulay+Null 12**; Si Γ free tight linear+DIIS cap / **Pulay+Null 14**; Al k222 kT 0.01 linear+DIIS 23 / Pulay+DIIS 23 / **Pulay+Null 11** (and the banked MnO free recipe is Pulay+Null).  ⇒ the default is a PAIR (mixer + accelerator), not a mixer flip: Kerker/Pulay wants NO Fock extrapolator on a seeded run; restarts keep the Ladder at its GDM rung (mixer bypassed).  USER 2026-10-22: build the pair.  **P1 sweep (2026-10-22: `SolidCalcOptions::accelerator` = `Auto` -- a seeded stage beside Kerker/Pulay gets `Null`, a carried kT = 0 stage the Ladder on its GDM rung; `Auto` mixing stays LINEAR beside any Fock extrapolator, `SCFAccelerator::ExtrapolatesFock()`): 2 failures (R0: 19), restarts fine, Response FD 2–4× faster (101→28, 111→42, 59→14 iterations), free Si k222 18→12 -- but NOT yet a default**, four blockers: (1) ACCURACY -- XC on the band-limited ρ̃_mix (pin 18): O₂ box +7.72 mHa (−30.70246 vs −30.71018 for linear, either accelerator); XC from ρ[D] (−30.71018173, 28 it) and Kerker + N4 cusp deficit (−30.71018173, did not converge in 120) both remove it exactly, so the cure is N4 on the PULAY arm (today Kerker-only; the (ρ,m) channel basis refuses it too); Si ≤ 6e-6, Mn₂ box 4.7e-5.  (2) ROBUSTNESS -- Si Γ +U 2 eV ground state DIVERGES under Pulay + Null (`ResponsePolarizability.GPW_Si_U2_FrozenChi_eqFiniteDifferenceLRT`, 200 cap, −7.0614); linear + DIIS converged it.  (3) COST -- the ρ̃ path per iteration on a big vacuum box: Mn box 1 → 22 s/iteration (test 10 → 410 s), O₂/Si/Na boxes 2–4×.  (4) the U loop on Ladder@GDM carried stages does not converge U in 12 outer steps (`GPW_Si.Γ_U_Kerker_Pulay_Cplx_LinearResponseUIsAFixedPoint`; DIIS stages did).  **LANDED NEUTRAL 2026-10-22** (`openwork-sweep`): `SCFParams::Mixer` {Auto, Linear, Kerker, Pulay} + `ResolveMixing`, `SCFAccelerators::Type::Auto` (resolves to DIIS today, `SolidCalculation::MakeAccelerator` is the seam), `SCFAccelerator::ExtrapolatesFock()`, `ChargeDensity::CanMixInGSpace`, deck keys `scf.mixer` / `"Auto"` accelerator; every default resolves to linear + DIIS as before.  Blocker (1) is GONE (N4 is the default, 2026-10-22).  NEXT (§3 ORDER step 2, after `seed-tiebreak`, before D-ENV-RERUN): switch `MakeAccelerator` to the P1 pair (+ a carried-stage linear step capped at 0.3), re-sweep against the 2026-10-22 baseline, and work blockers (2)–(4). |

**Closed 2026-10-26.**  (1) `tMixableDensity::Axpy` is the one binary mutation; `MixIn` = `ReScale(1-c)` + `Axpy(c)` with exact
no-op shortcuts (`6833216b`: energies within 1e-10 both tiers; three iteration counts moved by 3 -- the old one-expression mix was
FMA-fused).  (2) `MakeCopy` + `Dot` (block Frobenius) and `PulayDMixer`: history on D_in / D_out copies, B bilinear in (in,out),
the next D formed in place; `Mixing::Pulay` at `KerkerG0` = 0 routes to it.  O₂ box: −30.7101817 (the ρ[D] oracle −30.71018173),
`XC ρ = DM` on every iteration.  (3) PIN 28 IS A CRASH (user: "simply crash the run if cusp* or fitted! appears"): `RequirePin28`
in the SCF loop; the old fitted-XC control is now `GPW_Si.Γ_Becke_Imp_Kerker_FitXC_Pin28Throws`.  (4) THE DEFAULT (user chose D-Pulay
over Kerker/Pulay): `Auto` = Pulay on D beside Null on seeded and carried stages; Ladder@GDM / DIIS only on a RESTORED stage.
The Kerker pair's sweep (2026-10-26, same baseline): 2 failures (Si Γ +U diverged; the U loop false-converged on Ladder@GDM --
OpenWork §3 "GDM false-converges"), vacuum boxes 20 → 708 s.  The D-Pulay sweep: both pass; short tier 554 → 471 s before the PW
fix below, 536 s after (the PW XC grid is 8× the points); long-tier k211 FD 537 → 85 iterations.  Anchors: the Si Γ cluster
−5.6e-6 is a CORRECTION (converged tight, linear+DIIS and Pulay both −7.1150678444); k222 gate slop either side of −7.7784728;
the rest capped / expect-throw tests.  `GPW_Si.Γ_U_ACBN0_OuterLoopFeelsTheNewU`: its "tolU 1e-6 cannot converge in 3" premise was
loose-SCF slop (the +U potential commutes with D on Si's half-filled p) -- tolU 0 now forces maxOuter.  (5) Found by the crash:
the plane-wave XC ball (§3 "Every plane-wave number predates the 4× XC ball"): `G_OrbitalPairBall` (stated by `PlaneWave_IBS`), the
`fit=` source (an EXACT fit, user's definition: ρ always goes through the fit), the three PW unit tests rewritten to the density grid.


### ✅ 2026-10-26 β-MnO₂: finish the U-free α·P test on PW (closed same day; verdict: `Records/HubbardUHistory.md` 2026-10-26 "second half: PW")

| **β-MnO₂: finish the U-free α·P test on PW** (2026-10-26) | the Gaussian half is in (HubbardUHistory 2026-10-26: Δm / Δn per α nearly basis-insensitive, control vs p5 f3), and the Gaussian U spin response climbs monotonically with basis completeness (+0.0229 → +0.0303 vs PW +0.0385, O d the first mover).  The PW α run (`~/Code/qchem6-runs/MnO2/alpha/alpha_pw.json`, 80 Ry, α 0.5 eV on site-0 Mn d) SLOSHED (59 iterations at Δρ 1–6e-2) and was stopped | **ON HOLD for the PW XC-density fix** (user 2026-10-26: PW XC sees ρ truncated to the wavefunction bandwidth, not the full ρ[D] — so the PW U − U0 table, +0.0385, is re-JUDGED after the fix (pin 10), not refreshed).  NO seed could be banked: a PW state cannot be saved yet (row "PW states cannot be saved" below; the U0h rerun `alpha/pw_U0h_save` reproduced −92.6795194 Ha in 31 iterations and saved nothing).  Ready deck `alpha/alpha_pw2.json` (α 0.5 eV, startingRelaxRo 0.2; drop its `restartFrom` until PW save lands).  After the fix: rerun the U0h reference, then α (and the U arm) on the fixed binary, rerun with startingRelaxRo ~0.2 and/or α 0.1–0.2 eV, seeded from the PW U0h state (`state.save` → `restartFrom`); compare Δm(Mn0), Δ\|m\|(Mn1), Δn per unit α against the GPW control / p5 f3 rows — matching ⇒ the U deficit is the basis's RESPONSE, not the U functional.  Also: O d + dense Mn d together on p5 f3 (the remaining 0.008) · `Records/HubbardUHistory.md` 2026-10-26 |

## ✅ Same-exponent shells of different l merged silently (OpenWork §3; orbital half closed 2026-10-26) -- the row as it stood

| **Two shells of DIFFERENT l sharing one exponent on one atom are MERGED -- silently** (bit 2026-10-24) | the known, flagged-not-fixed inverted condition in `PG_Cart/Imp/IrrepBasisSet.C` ~76 (and the duplicate check in `PG_Spherical`/`PG_LibCint`): `GaussianRF::operator==` ignores L, so a same-exponent shell is a "duplicate", the new l is never added, and the higher-Lmax radial REPLACES the old one.  β-MnO₂ scan: VA + O d 1.2 (VA already has O p 1.2) loaded 172 functions instead of 184 (Cartesian +12 not +24) and Etot came out +0.75 Ha ABOVE the smaller basis; O d 2.0 loads 184 correctly.  No error, no warning: any hand-built or scanned .bsd can hit it (CleanCode R1.0b is blocked on the same merge).  The code comment: the naive flip fixes ORBITAL bases but moves the DENSITY-FIT bases (Si₂ size-consistency 0.0026 → 0.024 Ha) | NOW: make it LOUD for orbital bases (throw on a same-exponent different-l pair on one atom, naming both shells) while the fit-basis path is understood; then the real fix (keep each l its own shell; the fit path's non-additivity) with a reader test that counts functions · `src/BasisSet/Gaussian/PG_Cart/Imp/IrrepBasisSet.C` |

**Closed for orbital bases 2026-10-26.**  `ReadShells(Reader&, const Structure&, ShellMerge)` replaces the three verbatim copies of the
reading loop (the libcint copy also erased-then-inserted through an invalidated iterator).  `KeepEachL` (every `Orbital_IBS`): a shell is a
duplicate only when its radial AND l-set match; a same-exponent shell of another l is its own block.  `LegacyMerge` (the default, so every fit
basis): the historical loop verbatim.  UTs `Gaussian94Reader.ReadShellsKeepsEachL…` / `AnOrbitalBasisLoadsEveryShellOfASharedExponent`
(H s + p at one exponent: 4 functions on MnD Cartesian, MnD spherical and libcint -- was 1).  Found live: NiO's step-3 scan v2 (Ni p 0.30 /
0.90 on the TM file, whose Ni d window starts at 0.30) loaded 146 functions for 156; v3 188 for 198.  Sweep: 0 of 1092 moved.  Shipped-file
audit: every shared-exponent pair is in a fit file (A1/A2 coul+exch, cd, x) except `dzvp_test.bsd` Be s/p.


## ✅ xcMesh recipe fields dropped under cellKind Auto (OpenWork §3; closed 2026-10-26) -- the row as it stood

| **`xcMesh` knobs are DROPPED SILENTLY when `cellKind` is Auto** (2026-10-24) | `qcMesh::ResolveXCMesh(mp)` returns `BeckeXCParams()` (nR 40, L 29) whenever `cellKind==Auto` — the default — discarding any `nRadial` / `angularDegree` / radial kind the deck set; a β-MnO₂ mesh-convergence run came back bit-identical before this was noticed | honour the given fields when resolving Auto to Becke (Auto picks the KIND, not the recipe), or throw when a recipe field is set under Auto; a deck test that a set `nRadial` reaches the `[Becke grid]` line · `src/Mesh/XCPolicy.C` ~423 |

**Closed 2026-10-26** by refusal, not by honouring: the recipe fields' defaults (`MeshParams{}`) are not the Becke recipe's, so "honour what
was set" cannot tell a stated field from a default.  `RequireNoRecipeUnderAuto` (`src/Mesh/XCPolicy.C`) throws `std::invalid_argument`
naming every recipe field that differs from `MeshParams{}` under Auto, at the top of both `ResolveXCMesh` overloads; deck `Resolve` calls
it so `rundeck --check` refuses the deck (the run-time throw alone passed `--check`).  Rundeck's echoed records carry every field AT its
default, so they still parse.  UTs `XCPolicy.AGridFieldSetUnderAutoThrowsNamingIt`, `Deck.AnXCMeshRecipeFieldUnderAutoIsRefusedAtResolve`.

## ✅ The cusp-deficit guard refused an exact-fit basis (OpenWork §3; closed 2026-10-26) -- the row as it stood

| **The cusp-deficit guard refuses an EXACT-fit basis** (2026-10-25) | `SCFParams::XCCuspDeficit` (default on) makes `PairDensitySampler` throw at iteration 2 of a plane-wave run ("asks for XC from ρ[D] … this mixed density has no raster shadow") — but on plane waves the BALL ρ IS ρ[D] (the {G} ball holds every G_i − G_j), so the deficit is identically zero.  Workaround: deck `scf.xcCuspDeficit false` | the fit basis answers what the sampler consumes — "is my fitted ρ exactly ρ[D]?" (true for PW, false for GPW) — and the guard asks it; gate: the PW-1 smoke deck converges with the default `SCFParams` · `src/ChargeDensity/Internal/Imp/DensitySampler_Pair.C` |

**Closed 2026-10-26.**  `RequireRawForCusp(cd, isRaw, ballHoldsPairs)` (`DensitySampler_Pair.C`) skips the throw when the fit basis
answers `G_OrbitalPairBall::HoldsEveryOrbitalPair()` -- plane waves since the 4× XC ball (same day) -- because the BALL ρ there IS ρ[D]
and the deficit is identically zero; the XC column reads `fit=`.  Gate: the PW-1 smoke deck (`~/Code/qchem6-runs/MnO2/pw/MnO2_beta_A1AFM.r007.json`)
with `scf.xcCuspDeficit true` (the default) converges.

## Δρ gate "below the quadrature floor" (OpenWork §3, 2026-10-23) -- REFUTED 2026-10-26; the row as it stood

| **A Δρ gate below the quadrature floor can never pass, silently** (2026-10-23) | β-MnO₂ (VA, Γ, GDM stage at kT 0): the energy fell monotonically to −93.0237884 Ha while Δρ (`minDeltaRhoPerElectron`) sat at 2e-5…7e-5 — its minimum, 2.34e-5, is EXACTLY the `ρ_lost/N` column, the run quadrature's own error integrating the charge, constant all run.  The L1 gate is measured on that quadrature, so a `1e-6` gate on this cell only ever waits for `NMaxIter`.  Same family: a deck stage that does not set the [F,D] gates inherits `MinFD 1e-5` / `MinΔFD 1e-7`, which a non-variational GPW GDM stage (‖[F,D]‖ ~5e-4) cannot meet either | confirm the floor (Δρ of a converged state re-sampled: ≥ ρ_lost/N?), then make it LOUD: the iterator (or `deck::Vet` once the mesh is known) warns when a gate sits below ρ_lost/N, and the trace marks a Δρ that has hit it.  Gate such stages on ΔE (the `MinΔE` comment already says so for GPW) · `src/SCFIterator/Imp/SCFIterator.C` ~467 |

**The floor half is refuted (2026-10-26), from the code and the sweep log, no new run.**  The gate (`ResidualPerElectron` →
`Vxc_Quadrature::RealSpaceL1Distance`) is ∫|ρ_a−ρ_b| sampled on ONE quadrature: a constant charge-integration error (ρ_lost) is common
to both densities and cancels, so it is no floor.  Measured on `~/Code/qchem6-runs/anchors/r11_short.log` (every gridded trace row with
ρ_lost/N ≠ 0): GDM 28 of 49 rows below ρ_lost/N, min ratio 2e-4; PulD 126 of 984 and Null 131 of 1235, min ratio 1e-7; DIIS 0 of 30.  The
β-MnO₂ stall (Δρ 2e-5…7e-5, minimum 2.34e-5 = ρ_lost/N) is the GDM crawl the direct-min driver measures as successive-iterate distance.
The REMAINDER (the inherited [F,D] gates a GPW GDM stage cannot meet) stays in the tracker under its own name.


## Inherited [F,D] gates on a GDM stage (OpenWork §3, 2026-10-23) -- REFUTED 2026-10-27; the row as it stood

| **A deck stage that does not set the [F,D] gates inherits `MinFD 1e-5` / `MinΔFD 1e-7`** (2026-10-23; the row's Δρ-floor half REFUTED 2026-10-26) | a non-variational GPW GDM stage (‖[F,D]‖ ~5e-4) cannot meet them, so it only waits for `NMaxIter`.  The Δρ "quadrature floor" was not one: the gate is ∫\|ρ_a−ρ_b\| on ONE quadrature, so a constant charge-integration error cancels, and in the 2026-10-26 sweep log Δρ falls through ρ_lost/N on GDM rows (down to 2e-4 × ρ_lost/N) and on Pulay/Null rows (1e-7 ×); β-MnO₂'s stall at 2.34e-5 was GDM crawling (rows "GDM crawls…", "GDM false-converges…") | warn (iterator, at stage start) when a direct-min stage carries an [F,D] gate it cannot meet, or default a GDM stage's [F,D] gates off and gate it on ΔE · `src/SCFIterator/SCFParams.C` ~31/46 · history `Records/OpenWork_History5.md` |

**Refuted 2026-10-27 (both halves; the Δρ half on 2026-10-26).**  There is no GPW [F,D] floor under GDM: in the `fecad6be` sweep log
(`~/Code/qchem6-runs/anchors/r10r23_short.log`) every NaF GPW GDM stage is still falling ~3×/iteration when its ΔE/Δρ gates stop it
(`Γ_Imp_Anchor` cold fine: 1.39e-2 → 1.88e-3 → 7.76e-4; coarse: 3.6e-2 → 7.9e-3 → 2.5e-3 → 5.3e-4), and the Si GPW restarts
(`SolidState.ExactResume_*`, `WarmStart_OtherGrid_*`) reach 1e-7 … 1e-12.  ⚠ CORRECTED 2026-10-27: those restart numbers say how tightly the SAVED state was converged (the tests gate Δρ 1e-7 / ΔE 1e-10 before saving), not what GDM can reach -- and both are GTH pseudopotential runs, as is NaF.  MEASURED the same day instead (`~/Code/qchem6-runs/NaF/naf_fdprobe_*`, NaF SR2 Γ, GPW + GTH): the Kerker+Pulay state converged to max|dD| 1e-6 / ΔE 1e-11 has ‖[F,D]‖ 2.8e-7; GDM from it, gated on minFD alone, reaches 5e-11 in 58 iterations (wandering 3e-8…9e-8 on the way); GDM from a loose state (ΔE/E 1e-4) reaches 1.4e-6 in 10 iterations -- and is then kicked back by the stale-field fallback (OpenWork row "GDM's mixed FALLBACK step…").  No PP [F,D] floor above 1e-10 on this cell.  `minFD 1e-5` is a few iterations further at that rate
(extrapolated, not run).  The β-MnO₂ stage that ran to `NMaxIter` (`basis_scan/m0_free_twostage_prefix.log` stage 2, open gap 2.3e-2)
decays LINEARLY, ~7 %/iteration -- the crawl of row "GDM crawls…", where the evidence now sits.  And the [F,D] gate is the RIGHT one for a
direct minimiser: row "GDM false-converges…" wants GDM gated on its gradient, not on its own step.  Nothing to change in the defaults.


## ✅ The polarized IonicSAD seed lost its spin (OpenWork §3, 2026-10-23) -- RESOLVED by re-measure 2026-10-27; the row as it stood

| **The polarized IonicSAD seed loses its spin** (2026-10-23) | MnO AFM-II SR Γ, 24-op group, seed-tiebreak's TRUE ionic seed (Mn²⁺ d⁵ ×2 + O²⁻ ×2, 26 e): iterate 1 is −42.20 Ha with Mn moment **0.30** (neutral SAD: −59.85 Ha, 3.93); the moment peaks at 2.86 (iteration 14), dies by iteration 30, and the run ends non-magnetic at −60.97043 Ha, unconverged at iteration 82.  Two Mn²⁺ Hund pairs should give ~5 after one diagonalization.  Two more data points (2026-10-23): the VA +U gate's U=0 arm (`Γ_U_Shub_Pol_Smear_CP2K_Long`, same seed) also starts at 0.30, keeps the order but takes 38 iterations (18 on the old seed; energy unmoved, 3e-10); and the OLD per-species seed (neutral Mn + O²⁻, 30 e for 26) started −52.75 Ha / 1.80, neutral SAD −59.32 / 3.68 — the true ionic seed is the WORST start of the three | check the seed's OWN m(r) before the first Fock build (is the Mn²⁺ pair at target 5 found and its spin channels kept?), then whether O²⁻ (unbound in LDA, the diffuse charge-state entry) wrecks the first iterate.  BEFORE β-MnO₂ (HubbardUPlan S2): the seed-tiebreak fix exists to make IonicSAD seed its AFM · `src/ChargeDensity/Imp/{Seed,SeedCD}.C` |

**Re-measured 2026-10-27 on `3941b9e3`** (decks + logs `~/Code/qchem6-runs/MnO/mno_afm2_shub_seed_{IonicSAD,SAD}*`: `decks/MnO_AFM2_free.json`
with `siteSpins [1,-1,0,0]`, `imposeSymmetry true` → the 24-op Shubnikov group, SR, Γ, Kerker + Pulay 8/5, kT 5e-3):

| seed | seed moment | iterate 1 | converged |
|---|---|---|---|
| IonicSAD (Mn²⁺ ×2 + O²⁻ ×2) | 4.794 | −48.96 Ha, m 1.06 | 15 iterations, −61.3032807795, m 4.459 |
| neutral SAD | 4.781 | −59.54 Ha, m 3.85 | 14 iterations, −61.3032807795, m 4.459 |

The loss is gone.  The 2026-10-23 measurement predates the V_short image fix (pin 29, whose bias favoured FM), the Ewald fix, the per-SITE
IonicSAD targets (`aca82fd4`) and the stated ordering (`0f353ff5`); which of them removed it was not bisected.  The ionic seed's first
iterate is still the worse start (1.06 vs the seed's own 4.79 -- one diagonalization in the Mn²⁺/O²⁻ potential), but the order is back by
iteration 3 and costs one iteration.  Not gated by a new test: it would be a ~10-minute `_Long` twin of `Γ_Shub_Pol_Smear_Anchor_Long`.


## ✅ GDM's mixed fallback used the field mixer's stale density (OpenWork §3; found and fixed 2026-10-27) -- the row as it stood

| **GDM's mixed FALLBACK step builds its Fock from the field mixer's STALE density** (2026-10-27) | when GDM's line search is exhausted (`RejectStep` below `TrustMin` -- at a small gradient the energy decrease drops under the E noise), `DirectMinStep`'s `mixedStep` builds the Fock from `itsMixer->FockDensity(*itsCD)` -- and the Kerker/Pulay mixers IGNORE the argument and return their own running field `itsIn`, which nothing updates while GDM moves D.  So the "stable degradation" jumps back to the density of the mixer's LAST mix (at a restart: the restored state).  Measured, NaF SR2 Γ (GPW + GTH), GDM resumed from a loose state (`~/Code/qchem6-runs/NaF/naf_fdprobe_{L,D}*`): ‖[F,D]‖ 1.4e-1 → 1.4e-6 in 10 iterations, then iteration 14 (XC column `shadow` = the mixer's field) jumps +2.6 mHa to −24.42781, back at the loose start (−24.42979), [F,D] 2e-1; again near 29 and 44 -- three `direct-min MIXED fallback` calls, no console message.  From a TIGHT state (probe C) the same fallbacks are small (the stale field is near the answer) and [F,D] still reaches 5e-11 | re-seat a field mixer on the CURRENT density before the fallback uses it (a `tDensityMixer` face: present the working D as the running field, drop the history; a D-mixer needs nothing), and say so on the console; then decide whether an exhausted line search at a gradient below the E-noise floor should END the stage (converged at the noise floor) rather than step at all · `src/SCFIterator/Imp/SCFIterator.C` `DirectMinStep`, `src/ChargeDensity/Internal/{Kerker,Pulay}Mixer.C` |

**Fixed 2026-10-27.**  `tDensityMixer::Reseat(const cd_t& current)` (no-op default for a D-mixer, whose `FockDensity` IS the working D;
returns whether anything was re-seated) and `tFieldMixer::ReseatField(const GField&)`: `KerkerMixer` takes the current field and
re-presents it; `PulayMixer` also drops its history and restarts priming (residuals measured about densities the run has left would
extrapolate toward a mix of stale states); `PolarizedDensityMixer` re-seats each leaf in its channel basis (one `LeafFieldsOf` helper now
shared with `Mix`) and rebuilds the (ρ,m) channel pair.  `DirectMinStep`'s `mixedStep(why)` re-seats before building the Fock and prints
`[gdm] iteration N: <reason> -> one MIXED fixed-point step from the current D (field mixer re-seated on it, history dropped)` -- always,
not only under `gdm_trace`.  UT `JointPulay.ReseatFieldMakesAFreshMixerAtTheCurrentDensity` (a re-seated Pulay / Kerker steps exactly as
a fresh mixer seeded there).  GATE, `~/Code/qchem6-runs/NaF/naf_fdprobe_D*` (NaF SR2 Γ, GPW + GTH, GDM resumed from the loose state L,
gated on minFD 1e-7 alone): before, ‖[F,D]‖ 1.4e-6 at iteration 10, then kicks at 14 / ~29 / ~44 (+2.6 mHa, back to the loose start),
NOT converged at 60; after, two announced fallbacks (iterations 14, 19) with E flat to 1e-14 relative, CONVERGED at 23 (6.1e-8),
Etot −24.4303949715 = stage A's.  Short sweep: 1099 pass, 0 of 1097 moved (no existing gate falls back mid-GDM).


## ✅ After a stage that ended on a GDM step, the live orbitals and ε were not a valid reference (OpenWork §3; found and fixed 2026-10-27) -- the row as it stood

| **⛔ After a stage that ENDS ON A GDM STEP, the live orbitals and ε are not a valid reference: χ₀ is garbage** (2026-10-27) | NiO k222 (v3 basis, U_in 3), restarted from the converged overnight state (`states/NiO_AFM2.r010.h5`, EXACT RESUME, 1 iteration): χ₀(site 0, q=Γ) = **−0.00703** /Ha with the restart's default Ladder@GDM (−0.00714 on a second run: not even reproducible), **−2.89940** when the stage's last step is a diagonalization (`accelerator Null`) -- the overnight non-restarted run gave −2.89938.  Energy, moments, response gap identical in all three.  Mechanism, in part: GDM's `OrbitalsAt` returns the geodesically ROTATED orbitals with the PREVIOUS Fock's pseudo-canonical energies attached (`itsSeo`/`itsSev`), so the wave function's (C, ε) are not eigenpairs of F[D]; the sum-over-states R0 (and any other consumer of orbitals: gaps, ACBN0's occupied set, band printouts) assumes they are.  UNEXPLAINED: why a near-zero final step costs a factor 400.  The SAVED FILE is not affected (D is rebuilt from the written C at save; a restart reads only D) | (1) explain the factor 400 (compare C†FC with diag ε per block after the GDM step); (2) end every direct-min stage with a PSEUDO-CANONICALIZATION (diagonalize F within the occupied and within the virtual space: D unchanged to the bit, ε become F's) -- or a full diagonalization, the user's call; (3) make the response REFUSE a reference whose blocks are not eigenpairs of F to a stated tolerance (loud, like pin 28's throw); gate: the restart probe's χ₀ == the diagonalizing arm's · decks/logs `~/Code/qchem6-runs/NiO/kaxis/p0_{diag,gdm,t6}*` · `src/SCFAccelerator/Internal/Imp/SCFAcceleratorGDM.C` ~302 |

**The factor 400, explained.**  Audit (`ResponseKernel.GPW_Si_k211_GDMFinishedStage_OrbitalsAreEigenpairs`, per block: C†SC−I, the
occupied-virtual block of C†FC, diag(C†FC)−ε) on Si k211: a diagonalizing stage 1e-10 everywhere; the GDM-finished stage had exact
OCCUPIED ε, a 1e-8 gradient -- and VIRTUAL ε off by up to **2.7 Ha**.  `OrbitalsAt` moved the occupied block along the geodesic, then
COMPLETED the virtuals by an SVD of the projected old virtuals: their singular values are all ~1, so the SVD returned an arbitrary
rotation and reordering of the virtual space, and `e[no+a]=itsSev[a]` labelled each virtual with another one's energy.  Harmless to
GDM (each step re-derives its energies from F); fatal to every consumer of (C, ε): R0's denominators ε_occ−ε_vir came out far too
large, so χ₀ came out too small -- 5.4× on Si (−1.348 vs −7.271), 400× on NiO (−0.00703 vs −2.8994), and not even reproducible.

**Fixed.**  (1) `SCFAcceleratorGDM.C` `OrbitalsAt`: the virtuals ride the same unitary exp(A) as the occupied block,
X_v(θ) = X_v + U(cos θ − 1)U†X_v − Y V sin θ U†X_v, so column a stays pseudo-canonical virtual a and keeps its energy to O(θ):
Si virtual ε error 2.7 Ha → 8e-9, χ₀ −7.2705593 (diagonalizing arm −7.2705592); GDM 31 → 26 iterations (its virtual basis no
longer jumps).  (2) `Response::EigenpairResidual` (max|C†FC − diag ε| over the complex blocks) and `SolidCalculation::
HubbardLinearResponse` REFUSES (Outcome, Why::NotConverged) above **1e-4 Ha** -- χ₀'s relative error ~ residual/gap, so <~1e-3 on a
0.1 Ha gap while a mislabel is refused by >1000×.  Measured residuals: Si 1e-10 (diag) / 8e-9 (GDM); every response test in the
sweep ≤ 3.9e-8; NiO k222 GDM restart **1.03e-5** (its [F,D] noise 2.4e-6) -- an initial 1e-5 bound refused that healthy state.
The end-of-stage PSEUDO-CANONICALIZATION the row proposed was NOT built: with the root cause gone the residual is O(final step),
and the refusal measures it on every response.  Gates: the Si eigenpair test; NiO restart probe
(`~/Code/qchem6-runs/NiO/kaxis/p0_gdm_fixed2.log`) χ₀ −2.8993982 / χ −1.8254405 == the diagonalizing arm (`p0_diag.log`) to the
printed digits.  The SAVED FILE was never affected (D is rebuilt from the written C at save; a restart reads only D).
Sweep: 1102 pass; moved: `GPW_NaF.Γ_Imp_eqColdStart` −8.9e-9 Ha (38 → 37 iterations) and `Γ_Imp_eqExactResume` −2e-10 Ha --
the GDM rung's trajectory, within convergence noise, judged expected.


## A CARRIED density's first iteration fed XC the FITTED ρ (OpenWork §3, 2026-10-28) -- FIXED 2026-10-29; the row as it stood
| **A CARRIED density's first iteration feeds XC the FITTED ρ, so an "EXACT RESUME" does not start at its own fixed point** (2026-10-28) | TiO₂ k222 Kerker+Pulay (Null accel, kT 5 mHa) restart from a state converged at Δρ 1e-6: iteration 1 opens at Δρ 6.2e-4 and E 1.9e-5 Ha off, XC ρ column `fitted`, then `cusp` from iteration 2; the fresh Pulay history at α 1 amplifies it (Δρ to 2e-2) and the resume takes 15 iterations.  The saved D is fine: `SCFIterator`'s `SeedRole::WorkingDensity` path keeps the D-backed density (`itsCarried`), but the mixer's first field has no raster shadow, so `PairSource` returns `Fitted` -- legal under pin 28 only for a seed WITHOUT a D.  Hits every carried density: restarts, anneal stages, ±α, the U loop | seed the mixer FROM the carried D (`tDensityMixer::Reseat`: the D becomes the running field, so iteration 1 samples ρ[D] / cusp); gate: an exact resume's iteration 1 at Δρ ≤ the saved one · `src/SCFIterator/Imp/SCFIterator.C` l.236–246 / 370 |

**Verdict.**  Not a restart-only defect: EVERY field-mixer run's first Fock (and every GDM fallback re-seat) fed XC the fitted ρ
on the singles (Becke) route, because that route reaches ρ[D] only through the field's N4 cusp correction (ρ_mix − ρ_out,
`ExactSourceOf`), which the mixers formed only inside `Mix` -- a field taken straight from a D (the seed field, `ReseatField`, the
(ρ,m) channels' first rebuild) carried none.  A restart showed it worst: the saved fixed point re-opened off itself.  FIX: such a
field is presented WITH its correction, identically zero (`PresentFromD`, `FieldMixer.C`; the polarized (ρ,m) mixer forms it against
the working channels), and the working D is seated as the first field's DM source (`SCFIterator`, after `CreateMixer`).  The
proposed `Reseat` route was not needed.  Gates: UTCalculation `SolidState.ExactResume_{Pulay,KerkerPulay}_Null_StartsAtTheSavedState`
(Si k222 polarized, Becke, Null accelerator): Kerker-Pulay resume `fitted` / Δρ 8.6e-5 / 6 iterations BEFORE, `cusp` / Δρ 7.3e-10 /
2 AFTER; TiO₂ Γ deck iteration 1 `fitted` -> `cusp`.  Sweep 1114/1114; moved 5, all Kerker runs whose iteration 1 changed:
`GPW_NaF.Γ_Becke_Imp_eqUni` +2.9e-8 (12 -> 10 iterations), `GPW_NaF.Γ_Imp_Anchor` +3.7e-9, `GPW_Si.Γ_U_Imp_Pol_Kerker_eqUnpol`
+8e-10, `GPW_NaF.Γ_Imp_eqExactResume` +3e-10, `GPW_Si.Γ_U_Kerker_Pulay_Cplx_LinearResponseUIsAFixedPoint` (67 -> 66) -- the SCF
path, within convergence; banked `~/Code/qchem6-runs/anchors/firstfield_short.log`.

## ✅ An IMPOSED Shubnikov group did not hold the Mn0 ↔ Mn1 spin flip (OpenWork §3, 2026-10-28) -- FIXED 2026-10-29; the row as it stood

| **An IMPOSED Shubnikov group does not hold the Mn0 ↔ Mn1 spin flip: β-MnO₂ sites drift apart under imposition** (2026-10-28) | `~/Code/qchem6-runs/MnO2/uscan/control_U1.log` (A1AFM, 16 ops = 8 spatial + 8 spin-flip imposed, U on O p Γ2 = 1 eV, Pulay+Kerker 0.45, kT 1e-3): the [+U] Mn d occupations are site-equal to 1e-4 for the first iterations, then the two Mn sites' charges differ by up to 0.1 e and moments by 0.43 (+3.65 / −3.40 at iteration 300); U = 0.5, 2 eV and the p5f3_Od basis likewise, and converged U = 0 passes through 1e-4.  An imposed group should hold that to roundoff on every iterate.  Suspect: the G-space mixed density (the m channel, `ρ_mix` feeding XC and Hartree) or the [+U] occupations from D_out are not symmetrised under the spin-flip ops, so a near-soft ferrimagnetic mode grows | find the first iterate with a nonzero site asymmetry and which object (ρ̃ mixed, D_in, D_out, W_U) breaks it first; then rerun the uscan small-U leg · `Records/HubbardUHistory.md` 2026-10-28 |

**Diagnosis and fix (2026-10-29, `3b2a163a`).**  A 60-iteration trace of `uscan/control_U1` on a local build printing the [+U]
occupations to 12 digits (`~/Code/qchem6-runs/MnO2/uscan/leak/trace_U1.log`): at iteration 1 the FLIP-partner residual
(n(Mn0,↑) − n(Mn1,↓), n(O2,↑) − n(O3,↓)) is already 4e-9 while the σ=None partners (O2 ↔ O5, same spin) sit at 1e-11; it then grows
1e-8 (it 10) → 1e-6 (27) → 6e-2 (56), and the σ=None residual follows late (1e-11 → 1e-5) through the unprojected +U occupations.
SEED: the multigrid's coarse levels 15³ and 5³ (odd N) cannot hold τ = (½,½,½), which every Mn0 ↔ Mn1 flip op carries, while all
eight σ=None ops have τ = 0.  AMPLIFIER: by design (`BasisSet/Lattice/Imp/BasisSet.C`, the 2026-08-11 MnO run-37 note) only the XC
engine's (ρ,m) star-average imposed the flip ops; every per-channel consumer (Hartree, the mixer, the +U occupations) got the σ=None
subgroup -- "a PARTIAL projector is legal".  It is not, when the unprojected direction is SCF-unstable: intermediate U on O p_z
softens A1AFM toward the ferrimagnetic state (11.7 mHa below it at U = 0), so the seed grew.  FIX: `tBasisSet::ProjectSpinPair`
-- the GPW basis builds M_g on the Γ block's AOs (BuildOperationRep's shell permutation × angular rep, image centre matched modulo
the lattice) and `tCompositeWF` projects each block's (D↑, D↓) pair, D_s ← |G|⁻¹ Σ_g M_g D_{σ_g(s)} M_gᵀ, before building the
density, so ρ, m, the mixer history and the +U occupations are all imposed at once.  Re-traced: every flip and σ=None residual is
EXACTLY 0 over 60 iterations.  Gates: `UTLattice_BS ShubnikovPair.*` (two Mn bcc AFM, spherical d: the projected pair obeys
ρ↑(r) = ρ↓(r + a(½,½,½)) and a σ=None 4-fold pointwise; idempotent; charge conserved; a free run untouched).  Multi-k magnetic runs
keep the partial projection.  The leak was NOT why U = 0.5–2 eV fail to converge: with it fixed, Kerker+Pulay still cycles at
U = 1 eV (reaches −93.3663 at Δρ 4.5e-4, then a smooth ~10-iteration runaway, 0.1 Ha UP) -- a mixer problem, `HubbardUHistory.md`
2026-10-29.


## An EXACT RESUME of a field-mixer state walked away (β-MnO₂ A1AFM r015, U = 0; reported 2026-10-29) -- FIXED 2026-10-29
Probe `~/Code/qchem6-runs/MnO2/uscan/leak/resume_U0.json` (Kerker/Pulay, Null accel, kT 1e-3, imposed A1AFM): iteration 1 A −93.3798733
(saved −93.3837082), Δρ 9.8e-3, then −93.130 after 3 iterations.  Ruled out: the restore (E[D_saved] = −93.3809516 == the saved run's
E(internal); the 2.76 mHa difference to the saved A is its −TS) and the Shubnikov pair projection (the saved pair is symmetric; projecting
it changes nothing).  TWO causes, both in the restart's first steps: (1) a carried density was never the working density -- the
iteration-0 fill diagonalized F[D_saved] and THAT seeded the mixer, so the first mixed iterate sat two undamped Fock steps from the saved
state (a 6.4e-6 residual -> 1e-2 on this soft magnetic system); (2) the fresh mixer PRIMED with plain Kerker (pulayStart 5) AT the fixed
point and walked 0.29 Ha away before Pulay brought it back (23 iterations).  FIX: `SeedRole::ExactResume` (a restore with no WARM
difference; `SolidCalculation::Restart` passes it): the carried D IS the first working density; and an exact resume engages Pulay from
iteration 1 (`tSCFIterator::Iterate`, announced on the console).  r015: iteration 1 A −93.3837059 (2.3 μHa off), converged in 4
iterations to 1.3e-8.  NOT met, and cannot be: "iteration 1 at Δρ <= the saved state's" -- the saved D is the last OUTPUT density, one
residual from the fixed point, so its first Fock measures J·δ (2.4e-4 here, ~40x the saved 6.4e-6); the energy holds.  The Pulay
history is still not saved (a WARM restart's old residuals extrapolate onto the old fixed point; an exact resume rebuilds it in 2
steps).  Gates: UTCalculation `SolidState.ExactResume_{Pulay,KerkerPulay}_Null_StartsAtTheSavedState` (+ the immediate-Pulay line).
Sweep 1119/1119; moved vs `firstfield_short`: the restart tests (fewer iterations) and four Shubnikov box tests that move identically
on the 06260200 build without this change (the merge's, e.g. `GPW_Na2Box.Γ_Becke_Shub_Pol_OrderLostThrows` 105 -> 2835 iterations).
Banked `~/Code/qchem6-runs/anchors/exactresume_short.log`.  The saved state's own convergence belongs in the restart line: OpenWork row
"A restart's fingerprint is incomplete", item (4).


## [F,D] has no single home (OpenWork §3, 2026-10-28) -- DONE 2026-10-29; the row as it stood
| **[F,D] has no single home: the Null accelerator computes none, and DIIS and GDM each compute their own and COMBINE them differently** (2026-10-28) | (1) `SCFAcceleratorNull::GetError()` returns 0, so `FD < MinFD` always passes and Δρ alone decides convergence on every Kerker/Pulay run (the NiO / SrVO₃ / TiO₂ campaign decks); the trace's threshold cell still prints the default `(1e-05)` when `minFD` is unset, reading as a live criterion.  (2) Two copies of [F′,D′]: `tSCFIrrepAcceleratorDIIS::UseFD` and `tSCFIrrepAcceleratorGDM::UseFD`, same formula -- but the managers combine the irreps differently: DIIS √Σ‖E_k‖², GDM max_k ‖E_k‖.  So `minFD` is a different gate depending on the active rung, the Ladder's column jumps at the DIIS→GDM hand-off, and on N k-blocks DIIS's number is up to √N larger.  One place is EASY (user wanted it, 2026-10-28): every irrep accelerator, Null included, is already handed (F, D′) through `UseFD` | one free function `OrbitalGradient(F′, D′)` in `qchem.SCFAccelerator` used by all three `UseFD`s, and ONE irrep combination on the manager base: **max_k ‖E_k‖** (ruled by the user 2026-10-28: a convergence metric is judged on its WORST irrep, for any metric -- not for k-mesh independence, which a 1/N_k would give √Σ as well; DIIS's √Σ goes); Null keeps its irreps and reports it; then the `----` / `(n/a)` special cases in `SCFIterator.C` (l.779, 917) go away.  SCOPE (user 2026-10-28): one CONSISTENT treatment for every convergence metric, not just [F,D] -- the per-irrep ones ([F,D]; Δρ, today a whole-cell ∫|ρ_out−ρ_in|, measured per block from ΔD) are GATED on the worst irrep; the whole-cell scalars (ΔE/E, ρ_lost/N) have no per-irrep piece and get one normalization/threshold convention instead.  Display: the max; an RMS beside it only if it earns its column (each trace column must earn its width -- an error-entropy column was considered and declined as unreadable for most users).  Expect GDM-run iteration counts to move (the gate changes) -- anchordiff the sweep.  Then the RESPONSE-REFERENCE gate (today `kResponseReferenceMinΔρ` 1e-9 in `Deck.C`, refused for littleGroupK, a `--check` note otherwise; user 2026-10-28) becomes the worst irrep's ‖[F,D]‖ -- the orbital residual the symmetry measurement and R0 actually need, with no quadrature in it · `src/SCFAccelerator/Internal/Imp/SCFAccelerator{DIIS,GDM,Null}.C` |

**Verdict.**  One gradient (`OrbitalGradient` / `GradientNorm` in `qchem.SCFAccelerator`, used by the DIIS, GDM and Null `UseFD`s) and ONE
irrep combination, the WORST irrep (`WorstIrrep`; DIIS's √Σ gone, GDM already max).  The Null accelerator measures and reports its [F,D]
(its blocks register with the manager, as DIIS's do); no Null special case is left in the iterator, and its runs gate on `MinFD` like
any other (user: the real fix, no carve-out).  WHY IT MATTERED: `LinearMixer::Adapt` sets its step from [F,D] -- on a Null run it read 0, so
alpha grew x1.5 at iteration 1 (clamped to 1.0) and never moved again; the "coin-flip" `GPW_Na2Box.Γ_Becke_Shub_Pol_OrderLostThrows` was
an undamped symmetric 2-cycle that roundoff asymmetry used to break and the exact Shubnikov pair projection (3b2a163a) no longer did
(2835 iterations; now 210, single alpha, and every alpha of its old sweep converges in 208-211).  KNOCK-ONS: the atom DIIS engage gate
rescaled Z²·0.1/32 -> Z²·0.08/32 (the [F,D] unit changed; Kr Libxc went chaotic unrescaled); Δ[F,D] RETIRED as a default gate (user:
it is a stagnation measure, 1e-7 silently demanded [F,D] ~ 1e-7 and overrode MinFD); default NMaxIter 20 -> 100.  And the ruling "an SCF
test fails if it hits NMaxIter" (user), enforced at the main level (`qchem.SCFCap` + `BoundedSCF`; CLAUDE.md): 54 offenders found, 46
fixed or honestly bounded, 8 open (the row above).


## Eight SCF tests failed the cap guard (OpenWork §3, 2026-10-29) -- FIXED 2026-10-29; the row as it stood
| **Eight SCF tests fail the cap guard: they never converged** (2026-10-29; the guard `qchem.SCFCap` found 54, the rest fixed or honestly bounded) | `ValenceBasisGen.{Fluorine_q7, FluorineSeedDensityAnionIsDiffuse, AssembleValenceLowqFile}` (ONE F⁻ pseudo-atom SCF: [F,D] stalls at 9e-3, Δρ 3e-4); `OSi_PP_U.MultiSpeciesRouting` ([F,D] 2.1, diverging); `M_DFT.OxygenTripletBelowSinglet` (Δρ ~0.13, open-shell π*²); `P/A_PG_DFT.Energy/Z5` (B, Δρ 0.17, open 2p¹); `ResponsePolarizability.GPW_Si_U2_FrozenChi_eqFiniteDifferenceLRT` (±α SCFs diverge: 140 -> 264 iterations, 3 mHa); `f_U/A_HF_dfPin.SelfPin/BSpline_Z92` (converges 1.4e-8 from a 1e-8 pin its own gates cannot resolve).  All "passed" before by measuring an unconverged iterate | work each to a converged fixed point (the next session step); the U pin: tighten its gates until its energy is a fixed point, then re-judge the pin (pin 10) · this row, `CLAUDE.md` (the cap rule) |

**Verdict, one cause each -- none was the SCF code being wrong, each was a test (or recipe) asking for something its run could not deliver:**
- **U self-pin** (`A_HF_dfPin`): a 1e-8 pin under gates ([F,D] < Z·1e-5) that only resolve ~1e-8.  Recipe -> [F,D] < Z·1e-9, Δρ/N < Z·1e-8:
  every one of the nine cases sits on its fixed point and **every pin stands as written** (worst 4.7e-9, Eu; U 1.6e-9).  No re-pin.
- **OSi +U**: a closed-shell SINGLET for two open-shell atoms (O 2p⁴ ³P, Si 3p² ³P) at R = 20 bohr -> the quintet; 13 iterations.
- **O₂ singlet / B atom**: an integer determinant over an exactly DEGENERATE partly filled level (π*ˣ/π*ʸ; B 2p¹) flips every
  iteration -> smearing kT 1e-4 (equal fractions at any kT for a degenerate level; NIST's spherical atom).  Needed the
  converge-once `Calculation(st, opts, SCFParams, acc)` constructor (CleanCode D-CALC-DOUBLE-SCF, done for these callers): the old
  constructor ran its own SCF on DEFAULT params before the test's recipe could apply.
- **F⁻ pseudo-atom** (3 ValenceBasisGen tests): the pseudo-atom Hamiltonian is FITTED, so [F,D] floors (3e-5) while E is flat to 1e-9;
  the generator now gates its validating SCFs on the ENERGY (ΔE/E 1e-10), and a recipe can name its accelerator -- F⁻ (closed shell,
  2p unbound in the window, ε +0.07) runs GDM, where DIIS jittered 5e-6 Ha forever.
- **Si frozen-χ FD**: the ±α restarts run linear D-mixing on the Null accelerator; with its [F,D] now real, `LinearMixer::Adapt`'s
  ×0.8 / ×1.5-with-a-clamp rule returned to EXACTLY α = 1.0 and the run sat on a period-3 limit cycle.  Two fixes: (1) the user's --
  the shrink/grow ratio is IRRATIONAL (grow ×φ, and a step that would cross the ceiling closes only 1/φ² of the gap, so α never
  repeats); (2) `SCFParams::MaxRelaxRo`, the CEILING of the adaptive step, which the facade's restart recipe caps at 0.3 beside the
  start (CLAUDE.md's restart rule -- capping only the start let it grow back to 1.0).  Deck key `maxRelaxRo`.  Measured and rejected:
  a ceiling at the last retracted α (43 atom/molecule failures); √2 growth with a 1/φ approach (two atoms too slow).
Sweep (cap guard on) 1119/1119; no energy moved by more than 1e-8 except the FD test returning to its pre-change value (264 -> 160
iterations); iteration counts moved on PW_Si/NaF/CsI anchors, GPW_Al anneal and the Na₂ box (the new adaptive step).  Banked
`~/Code/qchem6-runs/anchors/capguard_short.log`.

## C0 tracker sweep (2026-10-29): closed rows moved verbatim from `doc/OpenWork.md`

Each row as it stood on main `5143eda7`.  Open REMAINDERs became new OpenWork rows or are owned elsewhere (see the stub at the end of OpenWork §3).

| row | what | record |
|---|---|---|
| ✅ Response: REPORT orbital-resolved channels + inter-site V_IJ | DONE 2026-10-21: `hubbardLinearResponse` `"resolve": "slot"|"eigen"` (channels = the +U term's U slots / occupation eigen-orbitals inside them, with compositions), U printed full / 1×1; V_IJ was already printed (`a18ef3dd`).  REMAINDER: the β-MnO₂ matched-QE oracle (HubbardUPlan S2); compositions are on GLOBAL axes (no local-octahedron frame); a resolved U is reported, not fed to `hubbardLoop` (`Uirrep`) | `Records/OpenWork_History5.md`, `doc/HubbardUPlan.md` S2/S6 |
| ✅ +U projector independent of the orbital basis | DONE 2026-10-25 (HubbardProjectorPlan increments 1–5; plan retired verbatim to `Records/HubbardUHistory.md`): the pseudo-atomic projector is the DEFAULT, `siteShells` names CP2K's LOWDIN, `radial`/`atomicRadial` retired; TiO₂ vs hp.x U(Ti) +4.2 % → +2.3 %, occupations = QE's.  REMAINDER: CleanCode D-RK-ONCE, D-CONTRACTED-SHELL; the +U campaign re-runs are §3 ORDER step (3) | `Records/OpenWork_History5.md` |
| **PW-3** ✅ 2026-10-25 → `Records/HubbardUHistory.md` "THE PLANE-WAVE VERDICT" | remainder: E_U +2 % vs QE, plausibly our TRUE set overlap O vs QE's represented B†B (opt-in B†B would settle it) | · `src/BasisSet/PlaneWave/Imp/PlaneWave_IBS.C` |
| ✅ LibXC-polarized | DONE 2026-10-21: XC_POLARIZED face (exchange ids by spin scaling); polarized N LibXC ≡ DiracVWN to 4 µHa | see `Records/OpenWork_History5.md` |
| ✅ GDM on UNPOLARIZED runs | DONE 2026-10-01: capacity g read off D′ (Tr D′²/Tr D′), nocc=N/g; gradient/model step need no g (ratio invariant). Still open: DIIS→GDM ladder restart on unpolarized *solids* (GPW) not re-measured | see `Records/OpenWork_History5.md` |
| ✅ MnO accuracy — name the operator | DONE 2026-10-23: the −99.7/−37 mHa VA offsets were TWO image-list defects of one class (pin 29): the short-range local-PP operator images (`BuildOpTerms`, −61/−97 mHa AFM/FM + the whole −37 mHa selective bias) and the Ewald real-space images (`EwaldEnergy`, −47.24 mHa).  VA now **+0.16 / +0.18 mHa** from CP2K (AFM / FM), Δ(AFM−FM) +1.43 vs +1.45, +U dE / E_U +0.11 / +0.009 mHa.  REMAINDER: the expanded a = 13 cell is −6.9 mHa off (CP2K capped at 200 steps there; an artificial ionic-metallic state); every MnO/NiO/oxide CAMPAIGN density, χ₀ and U carried the V_short defect -- re-run them (the Ewald fix shifts total energies only) | `Records/OpenWork_History5.md` |
| ✅ N4 — cusp-deficit XC feed | DONE 2026-10-22: wired on Kerker AND Kerker-preconditioned Pulay, both spin coordinates, and the DEFAULT (`SCFParams::XCCuspDeficit` = true, pin 28); the pair (raster) XC route serves it with its raw shadow and refuses it on the BALL fallback.  Moved anchors re-judged against their oracles, none changed verdict (largest: MnO +U E_U −0.37 mHa, 22.6 mHa from CP2K).  REMAINDER: (b) is the deficit iteration-static -- not measured; lever B (one gather per spin) is now unblocked by it | `Records/OpenWork_History5.md` |
| ✅ Ladder restart starts on the DIIS+Kerker rung | DONE 2026-10-21: `Restart` enters the Ladder at its GDM rung (`ladderStartOnGDM`); NaF exact resume 9→2, coarse→fine 12→2; smeared restarts retreat to DIIS.  REMAINDER: a restart whose density moved more (geometry, U) is unmeasured | `Records/OpenWork_History5.md` |
| ✅ **IonicSAD targets are per SPECIES, not per SITE** — DONE 2026-10-26 | `IonicSADTargets` returns one target per ATOM; `MagneticDecoration` and `SeedCD`/`PolarizedSeedCD` read them site by site (a radial class per (Z, target, flip); the G assembly sweeps one sub-cell per class rank, so a single-valence cell keeps today's calls).  Sweep: **0 of 1087 tests moved** (anchordiff 1e-10, iterations included).  Audit of every library material: ONLY **LiMn₂O₄_spinel** changes (Mn 3,3,4,4 = 2 Mn⁴⁺ + 2 Mn³⁺; was 4 for all, 2 e heavy) — its 4 decks (+U 0/2/4, gate1 probe) re-run on a charge-ordered seed.  REMAINDER: the tie-break puts the two charge states on SYMMETRY-EQUIVALENT 16d sites by atom index — a chosen charge order, not a derived one; state it like `siteSpins` if the spinel campaign needs a specific one | `Records/OpenWork_History5.md` "✅ IonicSAD targets are per SPECIES" · `src/ChargeDensity/Imp/{Seed,SeedCD}.C` |
| ✅ **The magnetic decoration is DERIVED FROM THE SEED** — +U half DONE 2026-10-26 | a polarized run with a Hubbard manifold now THROWS without `siteSpins` (`RequireStatedMagneticOrdering`: facade + `deck::Resolve`, so `scripts/vet-decks` and a revision-less refusal); the seed-derived path for the +U site groups is retired.  The 36 decks state theirs; **no recorded number moves** (all 36 ran IonicSAD, whose derivation was measured per structure; every difference from it is a uniform decoration, whose σ=None site stabiliser is the grey one, on unimposed decks).  REMAINDER: an IMPOSED polarized run WITHOUT +U still derives its decoration from the seed (the MnO SAD ferrimagnet of the original row is that path) — extend the throw to it, judging the integration tests that rely on it | `Records/OpenWork_History5.md` "✅ The magnetic decoration is DERIVED FROM THE SEED" · `src/Calculation/Imp/SolidCalculation.C` |
| ✅ **The polarized IonicSAD seed lost its spin** — RESOLVED by re-measure 2026-10-27 | on today's code the true ionic seed (Mn²⁺ ×2 + O²⁻ ×2, Shubnikov 24-op, SR Γ) keeps AFM-II: seed 4.79, iterate 1 moment 1.06, 15 iterations to −61.3032807795 (neutral SAD: 14, same energy to 1e-10).  History: `Records/OpenWork_History5.md` |
| ✅ **`xcMesh` knobs were DROPPED SILENTLY when `cellKind` is Auto** — DONE 2026-10-26 | Auto picks the KIND and its calibrated recipe, so a recipe field set beside it now THROWS naming the field (`RequireNoRecipeUnderAuto`, both `ResolveXCMesh` overloads; deck `Resolve` runs it, so `--check` refuses the deck).  History: `Records/OpenWork_History5.md` |
| ✅ **Two shells of DIFFERENT l sharing one exponent on one atom are MERGED -- silently** — orbital half DONE 2026-10-26 | `ReadShells` (`Point/Reader.C`), the ONE copy of the reading loop the three readers (PG_Cart, PG_Spherical, PG_LibCint) each carried: `ShellMerge::KeepEachL` for orbital bases (a duplicate = same radial AND same l-set), `LegacyMerge` (bit-identical) for the density-fit bases.  Bit NiO's step-3 basis scan the same day (Ni p 0.30 beside the Ni d window at 0.30: +2 functions, not +12).  0 of 1092 tests moved.  REMAINDER: the fit bases' merge (the inverted condition; Si₂ size-consistency 0.0026 → 0.024 Ha if flipped) -- every shared-exponent pair in the shipped files is in a FIT file (A1/A2/cd/x); `dzvp_test.bsd` Be s/p is the one orbital one | `Records/OpenWork_History5.md` · `src/BasisSet/Gaussian/Point/Reader.C` |
| ✅ A carried density is re-filled at kT = 0 | DONE 2026-10-21 (anchor batch item 1): a carried seed (`SeedRole::WorkingDensity`: Restart, BuildStage) is refilled under the RUN's policy when it is not the default (kT > 0 or MOM); FD stages take a restart recipe (linear relax ≤ 0.3, Ladder at GDM at kT = 0) | `Records/OpenWork_History5.md` |
| ✅ l=4 (g) shells in a lattice basis | DONE 2026-10-21: `SphericalShell(4)` was EMPTY in Release; general solid-harmonic generator, MnD g ≡ libcint; TiO₂ + Ti g now vets, builds +U, converges.  REMAINDER: the χ₀/U effect of Ti g on the b6 basis (b-series bases were scratch; rebuild b6 + Ti g 1.2) to finish TiO₂'s l-convergence | `Records/OpenWork_History5.md`, `Records/HubbardUHistory.md` §"TiO₂ ortho-atomic at Γ" |
| ✅ **After a stage that ended on a GDM step, the live orbitals and ε were not a valid reference** — FIXED 2026-10-27 | GDM's `OrbitalsAt` completed the virtuals by an SVD that scrambled them, then labelled each with another virtual's energy (2.7 Ha off on Si k211): χ₀ 5× (Si) / 400× (NiO) too small after any GDM-finished stage, restarts included.  The virtuals now ride the same geodesic; the response measures the reference's eigenpair residual and REFUSES above 1e-4 Ha.  NiO GDM restart: χ₀ −0.00703 → −2.8993982 (= the diagonalizing arm).  REMAINDER: re-judge any χ/U/gap/ACBN0 number produced straight after a GDM-finished stage before today (the NiO overnight run did not end on GDM).  History: `Records/OpenWork_History5.md` |
| ✅ **GDM's mixed fallback built its Fock from the field mixer's STALE density** — FIXED 2026-10-27 | `tDensityMixer::Reseat` (Kerker, Pulay, the polarized pair: the current D becomes the running field, history dropped); `DirectMinStep`'s fallback re-seats first and says so on the console.  NaF GDM from a loose state: 60 iterations unconverged with three +mHa kicks → 23 iterations to ‖[F,D]‖ 6e-8.  REMAINDER: whether a line search exhausted below the E-noise floor should END the stage instead of stepping (design, user's call) · history `Records/OpenWork_History5.md` |
| ✅ **Pulay without Kerker mixes D** — DONE 2026-10-26 | `PulayDMixer`: the history on density MATRICES (`tMixableDensity::Axpy` + `MakeCopy` + `Dot`), chosen by `Mixing::Pulay` at `KerkerG0` = 0 for any run; V_xc sees ρ[D_mix] -- O₂ box −30.7101817 with `XC ρ = DM` every iteration (`GPW_O2Box.Γ_Imp_M3_Pulay_Anchor`).  Pin 28 is now a THROW (`RequirePin28`), gated by `GPW_Si.Γ_Becke_Imp_Kerker_FitXC_Pin28Throws` | `Records/OpenWork_History5.md` |
| ✅ **The DEFAULT mixer: Pulay on D** (was "Kerker/Pulay as the DEFAULT") — DONE 2026-10-26 | user 2026-10-26: Pulay on D, not Kerker/Pulay (the Kerker pair left Si Γ +U diverging, the U loop unconverged and the vacuum boxes 20–35× slower per iteration).  `Auto`: seeded or CARRIED stage → Null + Pulay on D (α 0.45, depth 8, start 5); RESTORED stage (Restart) → Ladder@GDM at kT = 0 / DIIS smeared, linear capped at 0.3.  Sweep vs 6833216b: short tier 554 → 536 s, every blocker of 2026-10-22 gone, anchors judged (Si Γ cluster −5.6e-6 = a CORRECTION to the tight fixed point −7.1150678444; k222 gate slop either side of −7.7784728).  REMAINDER: molecules' `Auto` stays linear (D-Pulay is available to them by request); Kerker is a stated choice | `Records/OpenWork_History5.md` |
| ✅ A Γ supercell breaks its OWN translation symmetry (~2e-4) | DONE 2026-10-23: it was the V_short operator-image defect (pin 29) -- the Si 2×1×1 image sites now agree to the printed digits (U sites 1/3 2.079605 / 2.079605 eV, were 6.6e-4 apart), χ₀ q-mesh vs supercell 1.6e-4 → 9.3e-6; gate tightened 5e-4 → 5e-5 (U 1e-3 → 3e-4 eV) | `Records/OpenWork_History5.md` |
| ✅ valgen cannot mint a fully-stripped ion | DONE 2026-10-21: unset is −1, 0 is refused loudly (no SeedCD entry needed) | see `Records/OpenWork_History5.md` |
| ✅ **The cusp-deficit guard refused an EXACT-fit basis** — DONE 2026-10-26 | the guard asks the fit basis `G_OrbitalPairBall::HoldsEveryOrbitalPair()` (plane waves at the ≥ 4 Ecut ball); the PW-1 smoke deck converges with the default `SCFParams`.  History: `Records/OpenWork_History5.md` |
| ✅ An IMPOSED Shubnikov group did not hold the Mn0 ↔ Mn1 spin flip -- FIXED 2026-10-29 (`3b2a163a`, history: `doc/Records/OpenWork_History5.md`, end) | -- |
| ✅ [F,D] has no single home -- DONE 2026-10-29 (history: `doc/Records/OpenWork_History5.md`, end) | -- |
| ✅ Eight SCF tests failed the cap guard -- all FIXED 2026-10-29 (history: `doc/Records/OpenWork_History5.md`, end) | -- |
| ✅ A carried density's first iteration fed XC the FITTED ρ -- FIXED 2026-10-29 (history: `doc/Records/OpenWork_History5.md`, end) | -- |

Also moved (relocated, not closed): "Second oracle" is a RULE (CLAUDE.md's oracle table), not a feature row; "Seed library keyed by `q`" is folded into `doc/CleanCode.md` D-STRUCTDATA step 2 -- its scoping text is kept here:

| row | what | record |
|---|---|---|
| **Second oracle** (QE → ABINIT/VASP all built) | TRIGGERED ONLY by a question one oracle cannot answer (today: MnO −99.7 mHa below) | `Records/ParallelAndOraclePlan.md` Phase 4; also CP2K's 32-atom MnO supercell for the §7b OMP caveat |
| **Seed library must be keyed by the PP variant `q`** | every `atomic_*_densities.json` entry already carries `q`; the lookup does not take it, so a q3 run cannot ASK for the q3 seed -- today an ambiguous (Z,functional,Nelec) only THROWS (neutral since 9f926e37, explicit `Nval` since 5ddc6990).  A throw stops a wrong pick; it does not let the caller pick the right one.  Scoped 2026-10-02, NOT started (user: make it an OpenWork item): the seed code sees only a `Structure`, which has no `q`, so it is a plumbing change, ~9 sites: `AtomicDensity` Get/Has × {Density, SpinPair} gain a trailing `int q=-1` (filter entries to that `q`; unique neutral is then `Nelec==q`); `IonicSADTargets`, `MagneticDecoration`, `MakeSeedDensity<T>` (+2 explicit instantiations), `SeedCD`/`PolarizedSeedCD` ctors gain a `qByZ` map (default empty = today's behaviour); the facade (`SolidCalculation.C`, 3 call sites) builds it from `SolidCalcOptions::species`; `SCFIterator.C:151` passes none.  ⚠ a species whose `q` has NO library entry then stops silently using another variant's seed (HasAtomicDensity false → the existing neutral-amplitude fallback or a throw): expect the sweep to flag any test whose `species` q disagrees with the library's | WHO supplies `q` is the same question as `doc/CleanCode.md` D-STRUCTDATA step 2 (`pseudopotentials.json`) and `OOD-SOLID-Cleanup.md` V-CALCNET -- build it so the `qByZ` map is what step 2's PP set hands over, then do the plumbing once; fold D-BASIS-PP in · `src/ChargeDensity/Imp/{AtomicDensity,Seed,SeedCD}.C` |


## ✅ `u_bands` after a GDM stage: the virtuals carried ~0 slot weight (OpenWork §3, 2026-10-24) -- CLOSED 2026-10-29 by re-measure (C1 B0)
| **`u_bands` after a GDM stage: virtual bands carry ~0 weight on every +U slot** (2026-10-24) | `m0_imp_Op018_PA.log` (stage 2 GDM, kT 0): the virtual levels' ENERGIES match the smeared run's (LUMO→+1.92→+4.14 eV) but every virtual's slot weights sum to ≲0.01 (the stage-1 run: LUMO Mn d 0.65 / O p 0.24 / O s 0.08).  Occupied bands are fine.  So after GDM the virtual orbitals' coefficients do not belong to their eigenvalues (stale, or not in the AO frame) -- a diagnostic reading them is wrong, and so would be anything else that reads GDM's virtuals (a response, a band plot) | **FIRST (C0 2026-10-29): very likely the GDM-virtuals defect FIXED in `7fb46357` -- re-run `m0_imp_Op018_PA` stage 2 and close the row if the virtuals now carry slot weight.**  Then: check what GDM leaves in the virtual block's coefficients vs its eigenvalues (`tSCFIrrepAcceleratorGDM`), then make the virtuals canonical at the end of a GDM stage · `src/Calculation/Imp/SolidCalculation.C` `ReportBandSlotWeights` |

**Verdict (2026-10-29, lane A C1 B0): fixed by `7fb46357`** (GDM's virtuals ride the same geodesic as the occupied block, History5
"After a stage that ended on a GDM step").  Re-measured on lane-a (`c7d2a939` binary): `m0_imp_Op018_PA` re-run cold (both stages,
`~/Code/qchem6-runs/MnO2/basis_scan/m0_imp_Op018_PA_b0recheck.{json,log}`, −93.33526355, order SURVIVED, Mn ±3.412), then its GDM stage
alone restarted from that state with `QCHEM_DIAGNOSTICS=u_bands` (`m0_imp_Op018_PA_b0gdm.{json,log}`, 7 GDM iterations, ends on a GDM
step): the LUMO is Mn d Γ0 0.651 / O p Γ0 0.243 / O s 0.084, LUMO+1 Mn d 0.53 / O p 0.46, the next ones 0.86 / 0.88 / 0.75 on Mn d --
the stage-1 (smeared, diagonalizing) characters the row quoted, energies unchanged (LUMO→+1.92→+4.14 eV).  Seen on the way, for the
panel: that GDM restart stage is the linear crawl ([F,D] 1.40e-3 → 6.29e-4 in 6 steps, ~14 %/step) stopped by its ΔE/E gate (OpenWork
"GDM crawls…", the open-gap witness).

## ✅ GDM false-converges on a carried stage whose Hamiltonian CHANGED (OpenWork §3, 2026-10-26) -- CLOSED 2026-10-29 by re-measure (C1 B1)

VERDICT: the C1 B1 panel case #6, `GPW_Si.Γ_U_Ladder_LinearResponseUIsAFixedPoint_Panel` (the row's own setup: the Si +U
linear-response loop under an EXPLICIT `Ladder`), no longer false-converges: the loop runs **outer = 3, CONVERGED, U_SC
1.2399 eV** -- the DIIS route's value -- with the carried stages at 11 and 10 iterations (the cold first stage 79).  The row's
symptom was 12 outer steps NOT converged with U never moving the state.  Consistent with 1e0b2069 (2026-10-29: the default gate
carries [F,D] = MinFD, which is the row's own proposed fix, "gate a GDM stage on its GRADIENT"); NOT bisected.  The case stays on
the panel, GREEN, budget 125 (100 measured); log `~/Code/qchem6-runs/anchors/panel-07af0909-lane-a.log`.

The row as it stood:

| **GDM false-converges on a carried stage whose Hamiltonian CHANGED** (2026-10-26) | D-MINDRHO's Δρ = ∫|ρ_out−ρ_in| measures a DIRECT MINIMISER's own step, not its distance to the fixed point: a U-loop outer stage on Ladder@GDM (`GPW_Si.Γ_U_Kerker_Pulay_Cplx_LinearResponseUIsAFixedPoint` under the Auto pair) starts at [F,D] 4.6e-7 (Si p λ = ½: the new U's potential is small), GDM's preconditioned steps decay ~0.7×/iteration and the Δρ < 1e-7 gate passes in 2–10 iterations with Ekin/Een/Eee/Exc unchanged to 1e-10 -- U never moves the state, the loop runs 12 outer steps NOT converged (DIIS stages: 3 outer, U_SC 1.2399).  MITIGATED 2026-10-26: the `Auto` policy puts GDM only on a RESTORED stage (same Hamiltonian); an EXPLICIT Ladder on a carried stage still has it | gate a GDM stage on its GRADIENT ([F,D]) or on the aufbau density of F[ρ_in] rather than on its step; meanwhile the Auto pair puts GDM only on a true RESTART (same Hamiltonian) · `src/SCFIterator/Imp/SCFIterator.C`, `src/SCFAccelerator/` |

## Closed 2026-10-30 (lane A, B2 R2.5): the imposed no-+U decoration

Verdict: `RequireStatedMagneticOrdering` refuses an imposed polarized run without `siteSpins`; the seed derivation in `SolidCalculation::BuildBasis` is deleted; the 10 short-tier tests that derived now state the same values (`doc/SCFRobustnessPlan.md` "B2 as built").  The row, verbatim:

| **An IMPOSED polarized run WITHOUT +U still derives its magnetic decoration from the seed** (remainder of the closed "The magnetic decoration is DERIVED FROM THE SEED", 2026-10-26) | the +U path throws without `siteSpins`; the no-+U imposed path still derives it (the MnO SAD ferrimagnet of the original row) | extend `RequireStatedMagneticOrdering` to it, judging the integration tests that rely on it · `src/Calculation/Imp/SolidCalculation.C` |

## Closed 2026-10-30 (lane A, B2): the serial long tier timeout

Verdict: user ruling 2026-10-30 -- the routine long tier is `ctest -L long -LE panel` (default timeout); the panel runs on demand and its `_Panel_Long` passes carry TIMEOUT 14400 (root `CMakeLists.txt`).  The row, verbatim:

| **The serial long tier exceeds the ctest timeout** (lane A, 2026-10-30) | `GPW_MnO2.Γ_U_Shub_Pol_Smear_Kerker_Converges_Panel_Long` TIMED OUT at 7200 s serial (4268 s at 4 threads; its U2eV twin passed serial in 6410 s) -- so the serial long baseline the campaign now compares against cannot be green | raise that test's TIMEOUT (or the panel tier's) to ~3 h, or let the long tier run threaded with its own baseline · log `~/Code/qchem6-runs/anchors/lane-a-b2r21_long.out` |

## "A restart's fingerprint is incomplete" (OpenWork §3, 2026-10-29) -- gap (1) DONE 2026-10-30 (C4 step 5); the row as it stood

The functional is recorded from the run's `PseudopotentialSet` and sits in WARM; `pp.set` is REFUSE; atom positions and shell
centres are WARM (a moved atom prints ONE line with its max displacement); a state from before step 5 is migrated on comparison
(`Migrated`, `src/Calculation/Imp/SolidState.C`).  UTs `SolidState.AFunctionalChangeAndAMovedAtomAreWarmAndSayHowFar`,
`AStateFromBeforeStep5StillLoads`.  Gaps (2)–(4) stay open.

| **A restart's fingerprint is incomplete, and its WARM line names keys, not consequences** (user 2026-10-29: grids, fit bases, U, the functional all move the converged ground state -- a restart across them MUST iterate, and should say so) | Today (`qchem.SolidState`, `MakeStateFingerprint`): WARM records +U (site/l/U/α/projector/U per irrep), densityEcut, cutoffFactor, the XC mesh kind+size, multiplicity, ortho kind/tol, kT, block realness; each difference prints `differs: <key>: saved <a>, now <b>`.  GAPS: (1) **the functional is hard-coded `"LDA"`** (`SolidCalculation::Fingerprint`, l.1110) AND sits in REFUSE -- a PBE run records itself as LDA, and an LDA→PBE restart (a normal start) would be refused; it belongs in WARM, recorded truthfully.  (2) NOT recorded: the XC fit basis (`VxcFit`), the raster (BallOnly/AliasFree), the GPW tolerances (screenEps, mgridEcuts, relCutoff), the Becke eps.  (3) the WARM line gives the values but not WHAT they mean: user wants specifics, e.g. "XC fit basis changed (Delta -> PWFit) between the saved state and this run: the ground state moves, expect more iterations to reach the new one".  (4) (user 2026-10-29) the saved state's OWN convergence: the summary keeps the iteration count but not the final Δρ / [F,D] or which mixer converged it -- a state saved at Δρ 6e-6 by a Kerker/Pulay mixer on a soft system (β-MnO₂ A1AFM r015) re-opens one undamped Fock away at Δρ ~40× that, so even an EXACT resume should say "saved at Δρ 6.4e-6 by a mixed run: the first Fock amplifies that residual; expect iterations" | add the missing keys to WARM (functional moved REFUSE→WARM, recorded from the run); give each WARM key a one-line meaning + consequence beside its saved/now values (a table keyed on the fingerprint name, so an unknown key still prints raw); gate: a restart across each new key prints its line and converges · `src/Calculation/Imp/SolidState.C` (MakeStateFingerprint, Differences), `src/Calculation/Imp/SolidCalculation.C` (Fingerprint, the WARM print l.~1167) |

## ✅ "PW states cannot be saved (or restarted)" (OpenWork §3, 2026-10-26) -- FIXED 2026-10-30 (lane B); the row as it stood

`MakeStateFingerprint` records `basis.family` (REFUSE; a state without it migrates as "gaussian" -- no plane-wave state could be
written before) and, for a plane-wave basis, `basis.planeWaveEcut` in place of the AO shells: with the cell and each block's k
and size, already REFUSE, the cutoff fixes the G-set.  The density travels as before, per-block D (`WriteSolidState` /
`RestoreDensity`) -- the row's ρ(G) carrier was NOT built: D is nG² per spin block (334 MB at nG = 4571, β-MnO₂ 80 Ry), so a
large PW state is a large file; ρ(G) stays the lever if that bites.  `BuildPoint` computes the fingerprint whenever
`saveStateTo` is set, so an unsaveable run throws before its SCF (the deck's `state.save` sets it).  A restored `Auto` stage with
`numBands > 0` starts on DIIS (GDM needs every virtual).  Gates (Si Γ, 4 Ha): exact resume in 2 iterations, full and partial
eigensolve; a PW state refused by a Gaussian run on `basis.family`.

| **PW states cannot be saved (or restarted) -- and the deck says so only AFTER the SCF** (2026-10-26) | `state.save` on a `basis.family "planeWave"` deck: `SolidCalculation::SaveState` → `MakeStateFingerprint` throws "the Bloch block does not report its AO shells (AoShellSource)"; the CK-1 net catches it and prints `WARNING: state NOT saved` after the SCF.  β-MnO₂ 80 Ry U0h: 52 min, converged, nothing banked.  So no PW run can seed another (a site-asymmetric α·P from the α = 0 state, a cutoff ladder, the U arm), and `restartFrom` on PW cannot work either | a PW fingerprint (cell, Ecut, the G-set / k, spin) beside the AO one, and save/restore of the PW density (ρ(G) is the natural, size-agnostic carrier — D_ij is nG², 0.3 GB+ per spin at 40 Ha, cf. PW-2); meanwhile make it LOUD UP FRONT: `deck::Vet` / the facade refuse `state.save` on the PW family before the SCF · `src/Calculation/Imp/SolidCalculation.C` (SaveState, Fingerprint), `qchem.SolidState` |

| **A threaded NaF exact resume STALLS GDM at iteration 1, and its mixed fallback feeds XC the FITTED ρ** (lane A, 2026-10-30, found by the reproducibility fix) | `GPW_NaF.Γ_Imp_eqExactResume_Panel` at `QCHEM_OPENMP_THREADS=4`, now deterministic: 5 iterations vs budget 3 every run.  Iteration 1's GDM line search STALLS on the restored (4-thread-converged) state → `[handoff] GDM -> one MIXED fixed-point step (line search STALLED)` → that step's XC ρ column reads `fitted` although the banner says `rho_mix + cusp deficit` → first iterate 1.3e-4 Ha off the saved E (serial: no stall, first iterate at the saved E, ≤ 3 iterations).  Two questions: why does the line search find no descent from a converged state only when threaded (threaded collocation vs serial-order integrate-back no longer exact adjoints at ULP level?), and why does the RE-SEATED presentation route XC to the fit (`PresentFromD` with cusp should carry the zero N4 correction, pin 28; the 2026-10-28 carried-density fix) | print E(t) of iteration 1's line search at 1 vs 4 threads (`QCHEM_DIAGNOSTICS=gdm_trace`), and the XC source of the re-seated `FockDensity` · `src/SCFIterator/Imp/SCFIterator.C` `DirectMinStep`, `src/ChargeDensity/Internal/FieldMixer.C` `PresentFromD` |

**✅ VERDICT 2026-10-30 (lane A, 7a016d0b).** `gdm_trace` at 1 vs 4 threads: the cold runs are identical, and at the restored state ([F,D] 7e-15) BOTH thread counts find E(t) == E(0) TO THE BIT at all 12 backtracks -- there is no descent to find at a fixed point.  The retry then landed on a rounding coin-flip (-3.6e-15 serial = "descent"; 0 at 4 threads = STALLED), so threading only chose the side.  Two defects: (1) the STALL -- strict descent is ill-posed at the fixed point; a t=1 trial whose rise is <= 64 eps |E0| is now taken (`DirectMinStep`); (2) the FITTED ρ -- `HandOff`'s re-seat built a new presentation with no DM provenance, and the fallback's Fock ran before the next Mix seated it, so XC sampled the fit (pin 28); `HandOff` now seats the current D after a re-seat.  (2) was the walk-off on its own: with only (2) the stalled fallback reads `cusp`, stays at the answer and converges in 2.  Both fixes: 2 iterations at 1 and 4 threads.  No sweep test re-seats (grep `re-seated` in the main-2a3852c9 logs: none), so (2) has only this measurement.

**AMENDED 2026-10-30 (8650f2d2): fix (1) REVERTED.**  The rounding-level accept fired long before the fixed point: once |g| ~ 1e-7 the true descent g² ~ 1e-14 is already under 64 eps|E| ~ 1e-13, so GDM took every full step unbacktracked and limit-cycled ([F,D] 2e-7 ↔ 4e-7; the landing sweep caught `ResponseKernel.GPW_Si_k211_GDMFinishedStage` 41 → 220 and `GPW_Si.Γ_U_Ladder_LinearResponseUIsAFixedPoint_Panel` 179 → 400 = cap).  Fix (2), the re-seat seat, stands alone and is sufficient: the fixed-point stall still happens at 4 threads, but the fallback now stays on the answer.  `eqExactResume_Panel` passes at 4t; the four cases are 0 moved vs main-2a3852c9 serially.  Lesson: a "no descent is possible" test must look at the GRADIENT (or predicted descent), not at |ΔE| against the energy's rounding.

## ✅ "valgen cannot add a seed entry diff-friendly" (OpenWork §3, 2026-10-22) -- FIXED 2026-10-30 (lane B, d10d38cc); the row as it stood

`valgen --seed-into <library>` inserts or replaces ONE entry by (Z, q, Nelec, functional) as TEXT (`SpliceSeedEntry`,
`qchem.ValenceBasisGen`): a match is replaced in place, otherwise the entry is appended after the last one, in the file's tab
indentation; every other byte is kept (UT `ValenceBasisGen.SeedIntoChangesTheEntryAndNoOtherByte`, on the shipped 12.7k-line
library).  `GenerateSeedDensity` writes `q`.  Refused: an entry with no `q`, a library holding the identity twice, and replacing
a SPIN-RESOLVED entry with an unpolarized one -- the CLI smoke on Na would have silently dropped the shipped Na entry's
`rho_up`/`rho_dn` (the IonicSAD seed reads them); regenerate with `--spin`.

| **valgen cannot add a seed entry diff-friendly** (2026-10-22, the Mn⁴⁺ seed) | `valgen --seed-out` writes a standalone entry WITHOUT `q` (the PP variant every library entry carries, and the key the row above needs), and nothing inserts it into `atomic_valence_densities.json`: re-serializing the library re-formats all 12.7k lines (a 26k-line diff for one entry), so the Mn⁴⁺ entry went in by a hand-written text splice | `valgen --seed-into <library>`: insert or replace ONE entry by (Z, q, Nelec, functional) as TEXT, keeping the file's tab indentation and number format (diff = the entry only), always writing `q`; a UT that the round trip leaves every other byte unchanged · `CLIapps/valgen.C`, `src/ChargeDensity/Data/` |

## ✅ "A restart's fingerprint is incomplete" gaps (2)-(4) -- DONE 2026-10-30 (lane B, ce6b1f94); the row as it stood

(2) WARM now records the RESOLVED v_xc fit (`xc.fit`), the raster (`grids.raster`), `grids.ladderFactor`, the Becke eps
(`grids.beckeEps`, Becke mesh only) and every GPW tolerance (`gpw.<name>`; Gaussian family only); a state written before them
compares as agreeing (`Migrated`, keyed on the absent `xc.fit`).  (3) every WARM line ends in its meaning and consequence, from a
table keyed on the fingerprint name (prefix match; an unknown key prints raw), e.g. `grids.densityEcut: saved [20], now [30] --
the density grid changed: the converged ground state moves -- expect iterations to reach the new one`.  (4) the summary records
the final Δρ (`densityChange`, format 1.1) and the restart prints it, saying when a density-MIXED run (no [F,D]) saved it.  Found
on the way: `Same()` compared entries below 1 to an ABSOLUTE 1e-10, so screenEps 1e-10 vs 1e-12 read as the same run (now 1e-10
relative, 1e-14 floor).  Gate: `SolidState.ARestartAcrossTheLadderPrintsItsLineAndConverges` (the raster arm costs ~230 s on the
UT cell, so the ladder carries the gate), `TheNumericsKeysAreWarmAndAnOlderStateIsStillExact`.  REMAINDER (new row): the run
POLICY is not recorded.

| **A restart's fingerprint is incomplete, and its WARM line names keys, not consequences** (user 2026-10-29: grids, fit bases, U, the functional all move the converged ground state -- a restart across them MUST iterate, and should say so) | Today (`qchem.SolidState`, `MakeStateFingerprint`): WARM records +U (site/l/U/α/projector/U per irrep), densityEcut, cutoffFactor, the XC mesh kind+size, multiplicity, ortho kind/tol, kT, block realness; each difference prints `differs: <key>: saved <a>, now <b>`.  GAPS: (1) DONE 2026-10-30 (C4 step 5: functional WARM and truthful, `pp.set` REFUSE, positions WARM with max displacement; `Records/OpenWork_History5.md`).  (2) NOT recorded: the XC fit basis (`VxcFit`), the raster (BallOnly/AliasFree), the GPW tolerances (screenEps, mgridEcuts, relCutoff), the Becke eps.  (3) the WARM line gives the values but not WHAT they mean: user wants specifics, e.g. "XC fit basis changed (Delta -> PWFit) between the saved state and this run: the ground state moves, expect more iterations to reach the new one".  (4) (user 2026-10-29) the saved state's OWN convergence: the summary keeps the iteration count but not the final Δρ / [F,D] or which mixer converged it -- a state saved at Δρ 6e-6 by a Kerker/Pulay mixer on a soft system (β-MnO₂ A1AFM r015) re-opens one undamped Fock away at Δρ ~40× that, so even an EXACT resume should say "saved at Δρ 6.4e-6 by a mixed run: the first Fock amplifies that residual; expect iterations" | add the missing keys of (2) to WARM; give each WARM key a one-line meaning + consequence beside its saved/now values (a table keyed on the fingerprint name, so an unknown key still prints raw); gate: a restart across each new key prints its line and converges · `src/Calculation/Imp/SolidState.C` (MakeStateFingerprint, Differences), `src/Calculation/Imp/SolidCalculation.C` (Fingerprint, the WARM print l.~1167) |

## ✅ "A free run's XC mesh follows the group detected at EACH geometry, so a surface jumps" (F, opened and closed 2026-10-30, lane B)

FOUND by the first run of PES gate (c), `GPW_NaF.Γ_PathIsContinuous_Long` (F walked along x, δ = −0.095 … +0.095, h = 0.01 bohr,
Γ, VALENCE_LOWQ_SR2, Becke(20,2,24), free).  Both carries (AO, Löwdin) agreed to ~1e-11, so the SCF was not the cause.  The path:
the 16 points with |δ| ≥ 0.025 sat on a smooth sextic to 2e-9, the two at ±0.005 1.1e-6 Ha ABOVE it; E(δ) − E(−δ) grew exactly
linearly (an odd term of 9.4e-5 Ha/bohr at a site where rocksalt's mirror through F forbids one).  Cold starts (no carry):
δ = 0 sat 2.5e-5 Ha BELOW the curve (1.15e-4 at k222); the default Becke mesh (Lebedev 29) shrank both effects ~20× (step ~1.5e-6,
asymmetry 6e-8…1e-7); the UNIFORM XC mesh was exactly mirror-symmetric (1e-13) and smooth; the AliasFree raster changed nothing
(asymmetry 1.00e-6 vs BallOnly 9.85e-7, the same δ = 0 step).  CAUSE: KAxisPlan ruling 1 (a free run gets symmetric NUMERICS)
built the free run's XC mesh invariant under the space group DETECTED at its geometry -- 48 ops / 16752 points on the site, 8 ops /
17536 (52696 on the default mesh) off it -- so the quadrature switched as F moved, and off the site the 8-op group's site-adapted
angular set is not mirror-symmetric, so ±δ were quadratured on sets that are not mirror images.  None of the C2 rows (uniform-raster
equivariance, BallOnly aliasing, Kerker equivariance, the low-rank D factor, R1.0r) covers it: they are about one geometry.
FIX: `SolidCalcOptions::xcMeshSymmetric` (default true = ruling 1, unchanged for a standalone run; false = the plain atom-centred
mesh, `GPWParams::freeMeshInvariant`), resolved into the setup, so a session and `SolidCalculation(setup)` agree bit for bit;
`SolidSession` THROWS on a free Gaussian setup with it true (no silent override); the fingerprint's `grids.xcMesh` gains " plain"
(states written before keep their string).  The session ops are not knowable in advance for an arbitrary path, so the plain mesh is
the only frozen choice no path can invalidate -- and the cheaper one: 7.6k points against 16.7–17.5k on NaF.  AFTER: mirror 5e-10
over the whole path; the k-th differences fall geometrically (Δ³ 4.5e-7, Δ⁴ 1.0e-7, Δ⁵ 2.0e-8), the signature of a smooth curve.
The gate's jump measure moved from Δ³ (which conflates a jump with the curve's own E'''h³) to Δ⁵ < 1e-7 (a step s shows as up to
10 s); on the OLD code it reads 5.9e-6 (mirror 1.7e-5), so the gate catches this defect ~60× over.  Gates: the path test;
`SolidSession.RefusalsAreOutcomesAndAWrongCountThrows` (the throw).  REMAINDER (own row): the F site is a MAXIMUM at this level.

## PES / run-setup remainders (1)–(3) -- DONE 2026-10-30 (lane B); the items as they stood

(1) gate (c) -- E along a 20-point path has no jump above the SCF tolerance -- and the NaF Γ instance: `GPW_NaF.Γ_PathIsContinuous_Long`
(3.8 min serial, 71 MB), the item above.  (2) the Q3 `BasisConditioning` refusal UT: `SolidSession.AStepOutOfTheTrimWindowIsRefusedAsBasisConditioning`
(9f1889f9; Si Γ at orthoTol 1e-3: min eig S 0.027 at the reference, 4.5e-4 at a quarter-bond squeeze -> refused; a 0.05 bohr step
converges).  (3) the LÖWDIN carry of C, built and measured (3c940a3c, `SolidSession::Carry::Lowdin`): Si Γ 4 carried steps 41 vs 42
iterations (AO), energies equal to 1e-11; NaF Γ 19 carried points 110 vs 96 iterations (AO), energies equal to ~1e-11 -- no lever,
so the default stays AO.

## ✅ "The SCF observer has two homes" (OpenWork §3, 2026-10-30) -- DONE 2026-10-30 (lane B); the row as it stood

`SolidCalcOptions::onIteration` is deleted.  The per-iteration observer is the last argument of both `SolidCalculation`
constructors, of `SolidCalculation::Restart` and of `SolidSession` (before `Carry`), exactly as the molecular
`Calculation(setup, scf, observer)` / `MoleculeSession` take it: an observer is not a resolved decision, so it has no place in
the setup.  Attached before iteration 1 as before (stage 0 of an annealed recipe, a restart's first iterate, every session
point); `OnIteration` still REPLACES it for a later `Converge`.  Callers migrated: 3 `SolidState` UTs, 4 MnO, 1 O₂ box, 2 NaF
restart gates; the deck never serialised it.  No number moves (the observer only reads).

| **The SCF observer has two homes** (controller 2026-10-30, from lane G's molecular hook) | `Calculation` / `MoleculeSession` take the per-iteration observer as a CONSTRUCTOR argument (lane G: an observer is not a resolved decision); the solid side still carries `SolidCalcOptions::onIteration` inside the resolved options | move it to a `SolidCalculation` / `SolidSession` constructor argument, migrate every caller, delete the field (lane B) · `src/Calculation/` |

## ✅ "A restart's fingerprint does not record the run POLICY" — DONE 2026-10-30 (lane B, 21fe6e27)

`Deviation` gains `movesState`: the routes that change the Hamiltonian or its quadrature (dmLowRank, streamFold, xcFromDM,
beckeXC, dAwareScreen, hubbardEigen) are recorded as `policy.<route>` WARM keys from the run's OWN `RunPolicy(o.policy)`;
mixRhoM (mixing channels, identical energies) and imposeSymmetry (a permission; the imposition itself is resolved) are not,
so they cannot turn an exact resume into a warm start.  A flipped route prints "policy.xcFromDM: saved 'off', now 'on' -- a
declared CP2K-deviation route changed the Hamiltonian"; a state saved before the keys existed resumes exactly (`Migrated`).
UT `SolidState.ThePolicyRoutesAreWarmAndAnOlderStateIsStillExact`.  Not covered: the RESOLVED imposition (`setup->imposed`)
is not itself a fingerprint key.

| **A restart's fingerprint does not record the run POLICY** (remainder of "A restart's fingerprint is incomplete", gaps 1–4 DONE 2026-10-30) | the declared CP2K-deviation routes (`SolidCalcOptions::policy`: xcFromDM, beckeXC, dmLowRank, ...) change the Hamiltonian, so they move the converged ground state like a grid does, but `MakeStateFingerprint` records none of them: a restart across a policy change is called an EXACT resume | record the RESOLVED route table WARM (one key per route, `policy.<route>`, so a banner reformat is not a difference) with its consequence line; a state without them compares as agreeing · `src/Calculation/Imp/SolidState.C`, `Records/OpenWork_History5.md` |

## ✅ "Spinel charge order is CHOSEN by atom index" — DONE 2026-10-30 (lane B, 4badf3eb)

`SolidCalcOptions::siteCharges` (deck `solid.siteCharges`), one integer per site like `siteSpins`: empty = derived (the
electronegativity tie-break by atom index, unchanged); Resolve refuses a wrong count, a non-neutral list, or a seed other
than IonicSAD.  It reaches `SeedKey::siteCharges`, so `IonicSADTargets` feeds the seed AND the magnetic decoration from it.
`decks/probes/LiMn2O4_spinel_gate1.json` states the order the tie-break derived ([1,1,4,4,3,3,-2×8]: Mn⁴⁺ on Mn0/Mn1), so no
number moves; λ-MnO₂ has no tie.  UTs `IonicSADSeed.AStatedChargeOrderReplacesTheTieBreak`,
`Deck.AStatedChargeOrderIsOnePerSiteNeutralAndIonicSAD`.

| **Spinel charge order is CHOSEN by atom index** (remainder of the closed "IonicSAD targets are per SPECIES", 2026-10-26) | the tie-break puts LiMn₂O₄'s Mn³⁺/Mn⁴⁺ on symmetry-equivalent 16d sites by atom index -- a chosen charge order, not a derived one | state it in the deck like `siteSpins` when the spinel campaign needs a specific order · `src/ChargeDensity/Imp/Seed.C` |

## ✅ "The molecular facade's default XC mesh (30 radial / Lebedev 5) is 0.5 Ha/atom off" — DONE 2026-10-31 (lane B, 8a8e39b8)

`CalcOptions::mesh` is now MHL 80 / Lebedev 29 (302 directions; m 3, α 2, Becke 2 unchanged).  Scan (single thread),
reference MHL 150 / Lebedev 35, error per atom at 80/29: water LDA 1 μHa, N₂ Xα 3, F₂ GTH q7 27, Si atom 0.1, Si₂ at
R 20 bohr 53 (the slowest case: the far-separated pair).  The 30/5 error was ANGULAR -- water LDA 30/29 is 0.7 mHa from
converged, 80/5 still 65 mHa -- and NOT variational.  Cost: water LDA SCF 0.14 s → 2.8 s (20x); 100/35 would buy ~1 μHa/atom
for 1.75x more.  60/17 measured as good on these four (≤ 18 μHa/atom, 1/3 the cost) but non-monotone at the 50 μHa level
(F₂ 80/23 +74 μHa), so the default is the controller's 80/29.

Every moved anchor, re-judged against PySCF (`scripts/pyscf_molecule_oracle.py`: exact J and XC, the same Cartesian dzvp, 200/590 grid; PySCF HF reproduces our
HF water and O₂ anchors to 1e-8, so basis and geometry are identical):

| anchor | 30/5 (old) | 80/29 (new) | 150/35 | PySCF exact | new − PySCF |
|---|---|---|---|---|---|
| M_Calculation.WaterLDA (+ M_DFT.WaterPolarizedLDA) | −75.932462 | −75.872652 | −75.872648 | −75.877299 | +4.65 mHa |
| M_DFT.Water (+ WaterPolarizedSAD), Xα 0.74 | −76.149301 | −76.100897 | −76.100894 | −76.104768 | +3.87 |
| M_DFT.WaterSpherical, Xα 0.74 | −76.148556 | −76.100121 | −76.100117 | −76.104104 (sph) | +3.98 |
| M_DFT.N2, Xα 0.75197 | −109.217602 (pin −109.21679 at 2e-3 rel) | −109.225960 | −109.225954 | −109.227199 | +1.24 |
| M_DFT.OxygenTripletLDA | −149.256288 | −149.263716 | −149.263710 | −149.299178 (UKS) | **+35.5** ⚠ |
| Si2_PP_U E1 (Si pseudo-atom) | −3.759439 | −3.761718 | −3.761718 | -- | additivity |
| Si2_PP_U E2 (Si₂, R 20) | −7.516293 | −7.523336 | −7.523443 | -- | E2−2E1 2.6 → 0.10 mHa |
| F₂ GTH q7 valence_lowq_vb (lane G, no test) | −49.236805 | −48.192616 | −48.192561 | −48.200545 (CP2K −48.20058) | +8.0 |

The closed-shell rows sit 1-5 mHa ABOVE exact: PySCF's A1_coul J-fit alone moves 2 mHa DOWN, so the fitted Vxc (`A1_exch`)
carries it -- and the F₂ "PP route" row reads the same offset (re-attributed in OpenWork §3).  The O₂ triplet's +35.5 mHa
(+28.3 exchange-only) is a separate defect of the POLARIZED path: new §3 row.  The 30/5 values sat 10-55 mHa BELOW exact.
Test-side consequences: M_DFT.N2/Water drop the 2e-3 RELATIVE tolerance (0.2 Ha) that hid the mesh; M_Calculation.WaterLDA
gains a 6 mHa PySCF oracle line; M_Sym's same-mesh invariance runs keep 30/5 explicitly (107 + 175 s on the default) and drop
the retired Δ[F,D] gate (never met at the converged mesh's 1e-8 [F,D] wobble); OSi_PP_U gets kT 1e-4 (on a converged mesh the
O 2p triple is exactly degenerate and integer aufbau flipped it to the cap -- 30/5's 12 directions had split it); the
UTResponse LDA polarizability gate is 22 s on the new default (+0.6 / +1.4 / +0.7 percent vs PySCF CPKS).  A_PG_DFT (atoms vs
NIST) unchanged in pass/fail.

| **The molecular facade's default XC mesh (30 radial / Lebedev 5) is 0.5 Ha/atom off** (lane G S4, 2026-10-30) | `CalcOptions::mesh` (`src/Calculation/MoleculeSetup.C`: MHL, m 3, α 2.0, Becke 2; doc calls it "the proven molecular values" -- proven on HF only?).  F₂ LDA / GTH q7 / valence_lowq_vb: default −49.236805 vs converged −48.19259 (± 3e-5 from 80/29 up) = 1.044 Ha LOW; all-electron LDA water dzvp: −75.932462 vs −75.872652 = 60 mHa low.  Too LOW = the quadrature error is not variational; HF unaffected | recalibrate the default (80/29 converges both to < 2e-5/atom, measure the cost), gate it with an LDA molecule anchor at a converged mesh; re-judge any molecular DFT anchor banked on 30/5 (pin 10) · lane B · `src/Calculation/MoleculeSetup.C` |

## ✅ "Molecular PP route is 4 mHa/atom off CP2K on F₂" and "Molecular POLARIZED DFT is 35 mHa above PySCF on the O₂ triplet" — VERDICT 2026-10-31 (lane B, the δ-fit gate)

Neither route is defective.  Both numbers are the molecular A1 Gaussian FIT, in two parts: (1) the A1 Vxc fit, (2) the fit
bases' missing D shells (`LegacyMerge`, now OpenWork §3 "The molecular FIT bases drop their D shells").  The instrument was the new
δ Vxc fit basis on the molecular path (`DeltaFit_IBS<double>`, `CalcOptions::vxcFit = Delta`): exact XC on the mesh, so with the
SAME A1_coul Coulomb fit our energy must equal PySCF's density-fit route.  Measured (80/29; δ mesh-converged to < 1 μHa by 300/35;
energies in Ha, PySCF from `scripts/pyscf_molecule_oracle.py`):

| case | A1 (today) | δ (today) | A1, D kept | δ, D kept | PySCF A1_coul J-fit | PySCF exact |
|---|---|---|---|---|---|---|
| H₂O LDA dzvp | −75.872652 | −75.883203 | −75.880837 | −75.879455 | −75.879447 | −75.877299 |
| H₂O Xα 0.74 | −76.100897 | −76.110816 | −76.108315 | −76.106988 | −76.106980 | −76.104768 |
| N₂ Xα 0.75197 | −109.225960 | −109.248504 | −109.221696 | −109.228113 | −109.228107 | −109.227199 |
| F₂ GTH q7 vb | −48.192616 | −48.264376 | −48.197111 | −48.201937 | −48.201880 | −48.200545 (CP2K −48.20058) |
| O₂ LDA triplet | −149.263716 | −149.305761 | −149.301882 | −149.300684 | −149.300675 | −149.299178 |
| O₂ Xα 0.7 triplet | −148.902318 | −148.937147 | −148.933472 | −148.932130 | −148.932120 | −148.930624 |

("D kept" = the fit bases read with `KeepEachL`, an experiment, NOT landed: it breaks the Si pseudo-atom -- see the new row.)
With the D shells kept the δ route is PySCF's to 6-10 μHa (57 μHa on F₂): the molecular PP operators and the polarized machinery
are right.  The O₂ triplet's 35 mHa excess was the A1 per-spin Vxc fit on a d-less fit basis; with D kept, A1 is within 1-6 mHa
of PySCF on every case (the A1 Vxc-fit error proper).  The N₂ per-term check (δ, today): E_xc −13.39074 vs PySCF −13.39001, the
Hartree energy 26 mHa below -- the Coulomb fit, not XC.

COST, δ vs A1 (single thread, SCF wall / iterations; "today" | "D kept"):  H₂O LDA 2.67/13 vs 2.84/13 | 2.36/12 vs 2.80/12;  H₂O Xα
8.30/39 vs 3.32/16 | 2.71/14 vs 3.87/17;  N₂ 1.09/6 vs 1.11/6 | 1.10/6 vs 1.46/6;  F₂ 1.97/8 vs 1.75/8 | 1.97/8 vs 2.13/8;  O₂ LDA
10.0/22 vs 7.9/16 | 13.5/28 vs 14.7/24;  O₂ Xα 10.2/25 vs 19.8/41 | 15.6/29 vs 15.7/25.  Per iteration δ costs the same as A1 or
less on these molecules at 80/29.  The molecular DEFAULT stays A1 (Auto) in this chunk -- the user's call (pin 30 keeps A1 either way).

| **Molecular PP route is 4 mHa/atom off CP2K on F₂** (lane G S4, 2026-10-30; RE-ATTRIBUTED lane B 2026-10-31: probably NOT the PP) | PySCF with the SAME GTH-PADE-q7 + valence_lowq_vb, EXACT J and XC (no fit), gives −48.200545 (lda,vwn): CP2K's −48.20058 to 37 μHa, so the PP data and basis agree across two codes.  Ours −48.19256 (converged mesh) is 8.0 mHa ABOVE -- the same sign and size as all-electron water with NO PP (+4.65 mHa above PySCF exact, +6.8 above PySCF's A1_coul J-fit route): the DGauss A1 fitted-Vxc route (`A1_exch`), which the crystal route (GPW, no A1 fit) does not use.  Closed-shell sizes: water LDA +4.65, water Xα +3.9, N₂ Xα +1.2, F₂ +8.0 mHa | split the molecular E_xc: evaluate E_xc of the converged density ON THE GRID (no Vxc fit) beside the fitted value -- if that closes the 8 mHa the row becomes 'the fitted-Vxc route is 1-4 mHa/atom high', else gate E_loc/E_NL on the isolated F atom as first planned · oracle `scripts/pyscf_molecule_oracle.py` · `~/Code/cp2k-runs/F2/` |

| **Molecular POLARIZED DFT is 35 mHa above PySCF on the O₂ triplet** (lane B, 2026-10-31) | dzvp O₂ (R 2.282 bohr), multiplicity 3, converged mesh: LDA −149.263716 vs PySCF UKS lda,vwn exact −149.299178 (+35.5 mHa; +37.0 vs its A1_coul J-fit); exchange-only Xα 0.7 −148.902318 vs −148.930624 (+28.3) -- so exchange/fit, not VWN.  Our UHF triplet = PySCF to 1e-8 (state, basis, geometry fine); closed-shell water carries only +4 mHa, so ~30 mHa is the POLARIZED path (`FittedVxcPol`, the per-spin fit) | same split as the F₂ row: grid E_xc of the converged spin densities vs the fitted value; then the per-spin Vxc fit against a closed-shell ζ=0 run of the same molecule · oracle `scripts/pyscf_molecule_oracle.py`; `IntegrationTests/M_DFT.C` OxygenTripletLDA carries the ⚠ |

## ✅ "The molecular FIT bases drop their D shells" (D-FIT-MERGE) + the δ Vxc DEFAULT — DONE 2026-10-31 (lane B, f8e5848d)

**User ruling 2026-10-31: "δ default, A1 kept."**  The molecular and atom Vxc default is the δ basis (exact XC quadrature on
the run's own mesh); the fitted auxiliary basis stays selectable as `VxcFit::Auxiliary` (A1 Gaussians on a molecule, the radial
set on an atom) -- pin 30 keeps it as a strategy that holds the framework honest.  The numbers behind the ruling are the δ-vs-A1
table above.  `PseudoAtomOrbital` (the +U projector every +U solid consumes) STATES `Auxiliary`, so no +U anchor moves; the SAD /
IonicSAD seeds use the stored atomic densities and never run an atom SCF; ValenceBasisGen follows the default (an offline accuracy
tool: its UTs pass unchanged).

**D-FIT-MERGE.**  `ReadShells` has one rule for every Gaussian basis: a shared-exponent shell of another l is its own block;
`ShellMerge::LegacyMerge` is deleted.  N₂'s A1_coul now loads 68 functions (it loaded 32), UT
`Gaussian94Reader.TheFitBasesCarryEveryShellOfTheirFile`.  The Si fit metrics are well-conditioned with the D shells (normalised
min eigenvalue 1.2e-4 Coulomb, 2.7e-3 overlap): no conditioning rule was needed.  The "Si₂ size-consistency 0.0026 → 0.024 Ha" that
kept the merge was an ill-posed SCF: Si 3p² at integer occupation has no fixed point -- PySCF agrees with EXACT integrals -- and only
the d-less (spherical) Hartree hid it.  Degenerate open shells (the Si pseudo-atom, the O-Si quintet, the smeared O₂ singlet of
panel case 11) run at **kT 1e-2**, stated at each test: kT 1e-4 / 1e-3 slosh under DIIS, GDM and Ladder alike once the Hartree is
anisotropic (the rule GPW_SiBox already stated for the box).  Si₂ additivity at 20 bohr has no converging recipe (nor in PySCF):
an OpenWork §3 row; its routing check reads Enn on a bounded run.  Panel case 11's iteration budget 33 → 55 (46 measured at
kT 1e-2) -- lane A's panel, flagged.

**THE ANCHORS, judged once (pin 10).**  Independent route: PySCF with the SAME A1_coul Coulomb fit and exact XC
(`scripts/pyscf_molecule_oracle.py`); for the atoms, CP2K's ATOM code.

| gate | before | now | independent route | Δ |
|---|---|---|---|---|
| M_Calculation.WaterLDA (+ M_DFT.WaterPolarizedLDA) | −75.8726517 | −75.8794551 | −75.8794474 | 8 μHa |
| M_DFT.Water (+ WaterPolarizedSAD), Xα 0.74 | −76.1008972 | −76.1069880 | −76.1069799 | 8 μHa |
| M_DFT.WaterSpherical | −76.1001208 | −76.1075981 | −76.1075818 (Cartesian aux) | 16 μHa |
| M_DFT.N2, Xα 0.75197 | −109.2259604 | −109.2281129 | −109.2281074 | 6 μHa |
| M_DFT.OxygenTripletLDA | −149.2637158 | −149.3006843 | −149.3006746 | 10 μHa |
| M_DFT.F2_PP (was the δ gate −48.264376) | -- | −48.2019367 | −48.2018797 (CP2K −48.20058) | 57 μHa |
| M_DFT.WaterLDA_AuxiliaryFit (A1) | -- | −75.8808369 | −75.8794474 | −1.39 mHa |
| M_DFT.OxygenTripletLDA_AuxiliaryFit (A1) | -- | −149.3018824 | −149.3006746 | −1.21 mHa |
| M_DFT.F2_PP_AuxiliaryFit (A1) | -- | −48.1971113 | −48.2018797 | +4.77 mHa |
| Si_PP_U.Medium (Slater atom) | −3.7470781 | −3.7470398 | CP2K ATOM −3.747045 | 5 μHa (was 33) |
| Si_PP_U.Polarized | −3.7733032 | −3.7732782 | -- | |
| O_PP_U.SlaterHigh | −15.7478640 | −15.7478193 | CP2K ATOM −15.748119 | 300 μHa (was 255) |
| Si2_PP_U E1 (Si atom) | −3.7617185 (unsmeared artefact) | A −3.7818599, E −3.743669 (kT 1e-2) | −3.7436684 (smeared) | 1 μHa |
| GPW_O2Box finite twin (δ) | −30.0116800 | −29.9994810 | −29.9994536 (CP2K box −29.998761) | 27 μHa |
| GPW_MnBox finite twin (δ) | −14.6382750 | −14.6365281 | GPW box −14.6380455 / converged −14.63650 | 1.5 / 0.03 mHa |
| GPW_SiBox finite twins | −3.7617 (unsmeared) | internal −3.743669 (kT 1e-2) | -- | |

Also: the LDA FD polarizability (UTResponse) is 0.6 / 0.4 / 0.1 percent from PySCF CPKS (was 0.6 / 1.4 / 0.7); its GMRES tol is 1e-8
(the h = 1e-4 FD noise floor sits just above 1e-9 with the larger fit basis).  Measured in a scratch tree: ITMain molecular / atom /
PP / box set 156/156, UTCalculation 96/96, UTGaussian_BS 85/85, UTResponse LDA/HF 2/2; the solid anchors are the landing sweep's.

| **The molecular FIT bases drop their D shells: the Coulomb fit is 4-62 mHa LOW** (lane B 2026-10-31, the δ-fit gate; promoted from CleanCode D-FIT-MERGE) | `ReadShells(..., ShellMerge::LegacyMerge)` (`src/BasisSet/Gaussian/Point/Reader.C`) never adds a new l to a shared-exponent radial, and every A1 file has `SP` + `D` shells on the SAME exponents: N₂'s A1_coul loads 32 functions, not 68.  MEASURED with the δ Vxc route (exact XC on the mesh) against PySCF with the SAME A1_coul Coulomb fit: LegacyMerge -3.8 (H₂O) / -20 (N₂) / -62 (F₂ q7) mHa; KeepEachL 6-10 μHa (H₂O, N₂, O₂ triplet) and 57 μHa (F₂) -- quadrature level.  BUT KeepEachL breaks the Si pseudo-atom (sipp, A1 and δ alike: 100 iterations unconverged, E -3.7436 vs -3.7617) and Si₂ (E -2.4 / -5.6): the old 'size-consistency' symptom is a Si-specific CD-fit pathology, still unexplained | diagnose the Si fit (conditioning of the A1_coul Si block's Coulomb metric with its d shells; ConstrainedFF) before flipping the fit merge; then re-anchor every molecular DFT test against `scripts/pyscf_molecule_oracle.py` · CleanCode D-FIT-MERGE |

## ✅ "The SCF iteration row has two homes" — DONE 2026-10-31 (lane B)

ONE `IterationRow{iteration, cells, eb}` per iteration, built from the ONE column list (`tSCFIterator::RowColumns`) the
console formats: the invariant frame (Etotal, [F,D], ρ_mix, accel, svMin, cfg, gap, m_site) became `ColumnData` cells beside
the variable ones (their text writers are the old fixed writers' exact output), and `ColumnData` gained the cell's number and
plot scale.  The observer receives `IterationCell{title, tolerance, value (NaN for a "----"), text (padding stripped), scale}`
for every column, in display order; `IterationRow::ValueOf(title)` reads one by name.  Built when the run prints OR observes.
`SCFProgress` is retired with no shim, and its absolute |ΔE| with it (the row is exactly the console's, controller ruling): the
GPW harness computes its own.  Byte-identical console proven by goldens committed BEFORE the refactor
(`SCFTrace.MolecularConsoleIsByteIdentical`, `SCFTrace.SolidConsoleIsByteIdentical`: the trace block plus the tables printed after
it, which inherit the stream state the rows leave) and `SCFTrace.TheObserverReceivesTheConsoleRow` (titles in header order, each
cell text in its console line).  The console HEADER lines are still written by their own writers (unchanged, byte-identical) --
the row's cells carry the same titles.  ⚠ `pybind/qchem_bridge.cpp` still names `SCFProgress` (2 sites): lane G migrates the
bridge (flagged).  The ρ_lost/N "dead column" row is untouched: that column is one list entry, so if it goes, it goes from both.

| **The SCF iteration row has two homes** (lane G, 2026-10-31; for the GUI's live per-column plots, user 2026-10-31) | the console row (`IterationTrace` / `DisplayColumns`, molecular and solid layouts) and the observer's `SCFProgress{iteration, energy, dE, commutator, drho, order, eb}` are separate: the observer lacks Δ[F,D], ρ_lost/N, the XC-ρ flag, virial, ρ_mix, accel, svMin, cfg, gap/HOMO/LUMO and every column's title/tolerance; its `dE` is absolute while the solid console's is ΔE/E; the trace is built only under `Verbose` | ONE `IterationRow{iteration, cells[{title, tol?, value, text, scale}], eb}` built once per iteration when Verbose OR an observer is set; `DisplayColumns` formats it, the observer receives it; `SCFProgress` retired (≈23 consumers incl. pybind, which lane G migrates); the absolute \|ΔE\| is DROPPED (the row is exactly the console's).  Lane B, after `deck::Build`; design in lane G's note `lane-g-iteration-row.md` · `src/SCFIterator/Imp/SCFIterator.C` |

### anchordiff double-counts a solid stage's iterations (X; closed 2026-10-31, lane B)

`scripts/anchordiff` summed every `iters=` a test printed.  A `SolidCalculation` stage prints it TWICE: the result line
(`[label] CONVERGED iters=N ...`) and the per-stage summary (`[label stage k/n] kT=... iters=N A=E-TS=...`).  A test
that gained a stage summary between two sweeps therefore read as a doubled iteration count -- the 8ab42909 panel's only
"move" (a MnO2 resume, 5 -> 10) was this, not the SCF.  Fix: lines matching ` stage \d+/\d+\] ` are not counted (the result
line always prints).  Both logs are parsed by the same rule, so old baselines stay comparable.

## ✅ "CP2K oracle: is `ADDED_MOS` converged on the smeared TMO decks?" (lane G 2026-10-30) -- CLOSED 2026-10-31 (lane H); the row as it stood

| **CP2K oracle: is `ADDED_MOS` converged on the smeared TMO decks?** (lane G 2026-10-30) | the GUI's CP2K translator first sized `ADDED_MOS` at 16 vs the banked MnO AFM-II +U deck's 20 and came out 15 μHa high: too few added MOs truncate the Fermi-Dirac occupations.  Nobody has shown 20 is itself converged | re-run the banked smeared decks (MnO, NiO, β-MnO₂ +U) at 2× `ADDED_MOS`; a move > 1 μHa re-banks per pin 10 · `doc/Records/CP2Kresults.md` |

**Verdict.**  MnO AFM-II GPW VA (kT as banked), `~/Code/cp2k-runs/MnO_AFM2/`:

| deck | EPS_SCF | ADDED_MOS | steps | E (Ha) | E_U (Ha) |
|---|---|---|---|---|---|
| U = 0 | 1e-6 | 20 (banked) | 44 | −61.30332518 | — |
| U = 0 | 1e-6 | 40 | 44 | −61.303325178 | — |
| +U 4 eV | 1e-6 | 16 (GUI, r001) | — | −60.6859558 | — |
| +U 4 eV | 1e-6 | 20 (banked) | 104 | −60.68597088 | 0.60950923 |
| +U 4 eV | 1e-6 | 40 | 103 | −60.68595142 | 0.60950893 |
| +U 4 eV | 1e-8 | 20 | 136 | **−60.68588816** | 0.60950778 |
| +U 4 eV | 1e-8 | 40 | 117 | −60.68588676 | 0.60950785 |

U = 0 does not move at 2× ADDED_MOS.  At EPS_SCF 1e-6 the +U SCF (Broyden) stops on a ±15 μHa step-to-step oscillation, so the 16/20/40 spread does not follow ADDED_MOS; tightened to 1e-8, 20 and 40 agree to 1.4 μHa and both sit **83 μHa above the banked value**.  Re-judged per pin 10 against our own route (`GPW_MnO.Γ_U_Shub_Pol_Smear_CP2K_Long`, dE = 0.6174645): 110 μHa from the banked CP2K shift, 27 μHa from the tight one -- the independent route supports the move, so the +U anchor is re-banked at the 1e-8 ADDED_MOS 20 value (same deck as banked, EPS_SCF only) in `CP2Kresults.md` and `IntegrationTests/GPW/MnO.C` (gate tolerance 0.05 Ha, unchanged).  NiO has no banked CP2K energy; β-MnO₂'s oracle is QE (no CP2K deck) -- nothing else to re-run.  Lesson for the translator: a +U CP2K oracle quoted at μHa needs EPS_SCF ≤ 1e-8.


### Closed 2026-10-31 (lane A): GPW verdict lines repeat every SCF refresh

| **GPW verdict lines repeat every SCF refresh** (lane A, 2026-10-31) | Al fcc k222 smeared (`VALENCE_LOWQ_SR`): the `[grid audit] UNDER: V_loc long` verdict prints ~10× (once per iteration, as a console line + its `·` progress twin); the verdict is a property of the raster, fixed at setup | emit it once at setup (pin: emit contemporaneously with the activity, and the activity is the sizing) · `src/BasisSet/Gaussian/Lattice/Imp/GPW_Evaluator.C` |

**✅ FIXED f3f6ccdc:** `EmitGridAudit` speaks once per run report (the first k-block evaluator; the raster is k-independent): Al fcc k222 8 identical `UNDER: V_loc long` lines → 1.  The same commit has `gridAudit` state the reference raster's verdicts (its fence spans the reference's lazy first evaluation).


### Closed 2026-10-31 (lane A): `A_PG_DFT` asserts the SIGNED relative error

| **`A_PG_DFT` asserts the SIGNED relative error** (lane A, 2026-10-31) | `IntegrationTests/A_DFT.C` `EXPECT_LT(RelativeError(E, ref), tol)`, but `RelativeError` = (E_ref − E)/E_ref is signed and its doc (`src/Common/PeriodicTable.C:96`) says callers assert on `fabs`.  `P/A_PG_DFT.Energy/Z5` (B, Xalpha, dzvp, kT 1e-4) prints **−1.23 %** (−24.644 Ha, 0.30 Ha below its reference) on lane-a and −1.28 % on lane-b, and PASSES: any energy below the reference passes | wrap in `std::fabs`, then judge each `A_PG_DFT` case against a basis-matched reference (dzvp is not the NIST limit: the 1.2 % is either the reference's functional/α or a real defect) · `IntegrationTests/A_DFT.C` |

**✅ FIXED (the commit after 0bd60fa2): a pin-10 RE-JUDGEMENT, not a refresh.**  The defect was wider than the row said: EVERY Slater-Xalpha atom gate (`A_DFT_atom` SG/SL U+P and `A_PG_DFT` U+P, ~30 assertions) compared against NIST's atom table, which is **LSDA** (Dirac + VWN) -- a different functional -- and sat 0.005–12 % BELOW it (H −12.2 %, Li −1.9 %, B −1.3 %), passing ONLY because `EXPECT_LT` took the signed error.  The LSDA cases themselves match NIST to ≤ 600 ppm.  Independent route: PySCF UKS/RKS with E_x scaled by α/(2/3) (`1.5*alpha*SLATER`) at the same per-Z Schwarz α, uncontracted ANO-RCC, non-relativistic (`scripts/pyscf_xalpha_atoms.py`, run data `~/Code/qchem6-runs/Atoms/pyscf_xalpha_unc-ano-rcc.txt`); the reference's own basis error is ~1e-4 Ha at B (vs unc-aug-cc-pVQZ −24.661717 and vs our radial Slater High).  Two STATES at a degenerate open shell (`doc/SCFRobustnessPlan.md` §C1-DEG): the radial code's spherical M_L ensemble (PySCF frac_occ: B −24.653941; ours −24.653999) and the molecular path's pure integer component (PySCF B −24.661837; dzvp −24.656287, +198 ppm above, variational).  After: SG ≤ 121 ppm, SL ≤ 603 ppm, dzvp +135…+381 ppm one-sided; tolerances tightened to ~2.5× the worst.  Sweep 0 moved (assertions only).

## ✅ Benchmark §10 "More threads than asked on a vacuum box" — RESOLVED 2026-10-31 (lane H); the row as it stood

| row | measurement | next action |
|---|---|---|
| **More threads than asked on a vacuum box** (seen 2026-10-23) | `rundeck` on `Mn_box16` (VA, linear + DIIS) with `QCHEM_OPENMP_THREADS=4` showed ~14 threads at 70–80 % CPU in `top -H`; D-THREADS' audit says every `omp parallel` region takes `WorkerThreads()`, so the suspects are a library pool (OpenBLAS / FFT) spinning or doing work with `QCHEM_BLAS_THREADS` unset | sample it (`gdb -p … -ex "thread apply all bt 6"`) during the SCF loop and name the pool; if it is real work, decide whether it belongs under the one knob · `src/Common/Parallel.C` |

**Verdict (2026-10-31, c95c1fa4).**  The pool was **OpenBLAS**.  A perf profile of a β-MnO₂ PW `rundeck` run (6 OpenMP threads asked) showed ~12.5 cores busy inside zheev (`zhemv`/`zgemm`), with ~5% of samples in `exec_blas_async_wait`.  `FixBlasThreads` was called only by scfrun and the gtest mains; `rundeck` (2026-10-03) and `gpwprobe` (2026-09-15) never called it, so OpenBLAS sized its own pool (16 threads) while the banner printed the REQUESTED `BLAS=1`.  Fix: `FixBlasThreads` + `StopOmpThreadsBusyWaiting` now run from a load-time initialiser in `src/Common/Imp/Parallel.C`, the ONE home, and every main's copy was deleted.  The banner prints `openblas_get_num_threads()` with a WARNING if it differs from the request.  Checked: every exe links qcCommon, and rundeck reports `BLAS=1` by default and `BLAS=4` under `QCHEM_BLAS_THREADS=4`.

## ✅ PW-5: Γ Shubnikov imposition on the plane-wave block — DONE 2026-10-31 (lane H); the row as it stood

| row | why | next action |
|---|---|---|
| **PW-5: Γ Shubnikov imposition on the plane-wave block** (lane H, 2026-10-31) | PW refuses `imposed` (`SolidSetup.C`: "the plane-wave family has no symmetry imposition yet"); at U = 1 eV β-MnO₂ A1AFM is a SADDLE free (GPW free from the converged imposed state leaves it in 20 iterations, OpenWork §3 PW row), so the PW side of the GPW-vs-PW spin-response comparison at intermediate U cannot be measured without it | the PW analogue of `tBasisSet::ProjectSpinPair` (pair average of (D↑, D↓) over the magnetic Γ group, σ-flip ops swap the channels) acting on G-vectors; gate = the `ShubnikovPairUT` real-space checks on a PW block, then β-MnO₂ U = 1 eV PW imposed converges A1AFM · `Records/HubbardUHistory.md` (2026-10-31 DECISION) |

**Verdict.**  `PW_BasisSet::ProjectSpinPair` landed (08637da0, on main via 3d49c0c6).  It applies the phased G-permutation e^{-2πi m'·τ}, m' = W⁻ᵀm, per Shubnikov op and swaps the channels under σ-flip.  Its gates: `ShubnikovPairUT` (two PW tests: the flip and the rotation in real space; a projector that conserves charge and leaves a free run alone).  β-MnO₂ PW 40 Ha, U = 1 eV, imposed (r021): CONVERGED A1AFM in 87 iterations, E −92.6627400205, Mn ±3.4869.  The physics result is in `HubbardUHistory.md` (2026-10-31, "The PW U = 1 eV point").  Remainder: its ~30 s serial cost per iteration is the OpenWork §4 bullet.


## ✅ "⛔ GPW `BallOnly` raster + the C = 2 density floor ALIAS tight pairs" -- CLOSED 2026-11-01 by C2 step 3 (lane A, landed on main 14be0052); the row as it stood

Its ruling (2026-10-29: the raster N is SIZED FROM ONE ε, `BallOnly`/`AliasFree` retire as names) is delivered: GPW rasters are demands sized by the pilot audit on what the client consumes (pin 31, `doc/GridEpsilonPlan.md` step 3); 11 anchors re-banked against independent routes.  The REMAINDER -- the auto floor ignores a local-PP field sharper than the basis -- stays in `doc/OpenWork.md` §3.

| **⛔ GPW `BallOnly` raster + the C = 2 density floor ALIAS tight pairs** (2026-10-25; was "local-PP matrix elements wrong for a tight primitive") | β-MnO₂ + one Mn d at α = 76.7 collapsed to −6917 Ha "converged" (Eloc −8764 Ha).  Unit tier (UTLattice_BS `GPW.TightPrimitiveVlongOracle` / `…CollocationOracle`, closed forms at any α): every SIZING rule follows α (density floor, V_loc-long harmonic ball, completion rung at 2·E_fine) — the defect is the real-space SAMPLING: a pair is routed to E_cut ≈ 2p, and on a `BallOnly` raster (Nyquist = ball edge) it aliases at e^{−G_N²/p} = e^{−4}: a tight pair's density 4–16 % off at low G, its V_loc-long diagonal ±1.5–2.7 Ha (non-variational ⇒ the collapse).  `AliasFree` at the same C = 2: ≤ 4e-5 / 3e-4 Ha.  The raster doc's premise ("CP2K's bet, controlled by Ecut") holds only with CP2K's REL_CUTOFF-scale cutoffs, not with C = 2.  Second, latent: the auto floor ignores a local-PP FIELD sharper than the basis (α = 2 alone: 0.03 Ha even alias-free).  **Grid audit landed (d3a487e1): `grids.auditDensity` / `grids.auditVlocLong` + an `[grid audit] UNDER` console line** — the β-MnO₂ baseline reads UNDER on both (alias 1.8e-2 / 7.1e-3).  **It also breaks EQUIVARIANCE** (2026-10-28, TiO₂ k221 converged to Δρ 1e-9, KAxisPlan §4a): the orbital-image residual of the P4₂/mnm ops is 9.96e-7 on BallOnly vs 9.0e-9 AliasFree, and Etot moves 0.54 mHa -- one more argument for the AliasFree default.  **On an ANCHOR material at the mHa scale** (C2 M1, 2026-10-30, `doc/GridEpsilonPlan.md` §5): `decks/bench/NaF_SR2_Gamma_parity` E_AliasFree − E_BallOnly = **−941 μHa** (BallOnly audit UNDER, density alias 1.5e-2 at the p = 80 Na pair; MnO free −23 μHa, alias 9.4e-3; Si Γ/k222/k222s, audit clean, +9/+14/+15 μHa); AliasFree costs 2.2–6.2× wall, so BallOnly stays default and the ε-sized N is the cure -- M2a (2026-10-31): a FIXED-DENSITY quadrature error, ≈ −864 μHa of it the uniform XC raster (Becke XC: BO−AF −77 μHa), so NaF's §5a parity gap vs CP2K is mostly BallOnly's (+0.82 vs AliasFree −0.12 mHa) | **RULED 2026-10-29 (user, `doc/QueueCampaignPlan.md` §7): the raster N is SIZED FROM ONE accuracy target ε (alias ≤ ε per multigrid level), so `BallOnly`/`AliasFree` retire as names; consistent error across terms; lane A chunk C2.**  Earlier framing: **DECISION (user):** default raster `AliasFree` (8× raster points, every GPW anchor moves) vs `BallOnly` with a raster-aware floor (E_cut ≳ 5.5 p_max); the field floor max(2α_max, the PP field's need); a term-specific raster if the audit shows terms need different grids.  Then re-bank the anchors (`scripts/anchordiff`), and make an SCF whose energy falls below any physical bound FAIL · `src/BasisSet/Gaussian/Lattice/Imp/GPW_Evaluator.C` (EnsureLevels, MakeLocalPPLong), `SolidCalcOptions::raster` |

## ✅ "A bare `xcMesh.cellKind Becke` is the struct default" -- CLOSED 2026-11-01 (lane A, with step 3b); the row as it stood

Fix: `qcMesh::XCMeshDefault()` = `BeckeXCParams()` with `cellKind Auto` is `SolidCalcOptions::xcMesh`'s default (one source of the Becke numbers; the deck overlays on it), so stating only the kind runs (40, 2.0, Lebedev 29).  The `MeshParams` struct defaults are UNCHANGED (molecular, projector and fit meshes default-construct one).  `RequireNoRecipeUnderAuto` treats a field equal to either default as not stated (records written before the fix carry the struct values under Auto).  Moved (anchordiff vs main 14be0052, |dE| > 1e-10): no integration anchor; the Si unit tests on a bare stated Becke rose +74 to +79 mHa (the 12-direction rule over-bound): `ResponseKernel` q13 Hermitian / thread-identity, `ResponseLittleGroup.GPW_Si_k222_q211_Becke_ReducedEqualsFull`, `OrbitalImage` free/imposed Becke carry-over, the two `SolidState` field-mixer exact-resume panel cases (reference SCF 7→12 / 9→10 iterations; the resume claims hold) -- all identities, re-judged by their own claims.  The Fock thread-identity gate and `SessionSizesItsRastersOnMatrixElements` now STATE the old (30, 1.0, 5) mesh (the first: a new §3 row; the second: its subject is the rasters, 107 s → 17 s).

| **A bare `xcMesh.cellKind Becke` is the struct default (nR 30, Lebedev degree 5), not the calibrated `BeckeXCParams()` (40, 29)** (2026-11-01, lane A) | deck-reachable: `xcMesh: {cellKind: "Becke"}` with no recipe fields keeps `MeshParams{}` (nRadial 30, angularDegree 5 = the 12-direction rule), and so may the in-tree callers that set only the kind (`grep -rn "cellKind=qcMesh::UnitCellKind::Becke;"`: `SolidSession.C`, `SolidState.C`, `Kernel.C`, `Boxes.C` ×2, the Deck and Structure tests -- check each for recipe fields set after).  Measured: Si Γ max|ΔH_xc| 1.8e-2 per element against (40, 29) -- a silent substitution of a far coarser grid than the one `Auto` would build.  Not fixed in step 3b (it moves those anchors) | resolve a stated Becke whose recipe fields are all at the struct defaults to `BeckeXCParams()` in `ResolveXCMesh` (or make the struct's Becke defaults the calibrated ones), re-judge the moved anchors; the deck writer should then record the resolved recipe · `src/Mesh/XCPolicy.C`, `src/Mesh/Mesh.C` (MeshParams) |

## ✅ β-MnO₂: the U spin deficit (O p → Mn d cross response) — RESOLVED 2026-11-01 (lane H); the row as it stood

| row | why | next action |
|---|---|---|
| **β-MnO₂: the U spin deficit is NOT in the Mn d → Mn d response (charge, spin AND orbital match PW) — probe the O p → Mn d CROSS response** (2026-10-27) | Mn d α (94–105 %), α↑ 0.05 eV (107 %), per-slot α_k 0.05 eV (code `ef056880`, deck `alphaSlot_eV`: every slot, spin, site within the 1e-4 print) — all PW's.  But the deficit's U arm is U = 4.62 eV on the O p Γ2 slot (plane normal), U = 0 on Mn d: the Mn d moment deficit (60–79 %) is the CROSS response O p_z → Mn d moment, which no Mn d probe touches | 2026-10-28 O p_z probe (α −0.1 eV, both spins): the CHARGE cross response is PW's (3–5 %) and converges to it with the basis; the Mn \|m\| response is OPPOSITE in sign to U's, 1.7–1.8× PW's and does NOT move with the basis — so the basis-sensitive part is U's SPIN-dependent part −U(n_σ − n̄).  The leak was fixed (`3b2a163a`), and the 2026-10-29 plain-Kerker ladder (0, 0.5, 1, 2, 4.62 eV) shows the Mn moment is NON-MONOTONIC in U: it DIPS (−0.0076 / −0.0061 at 1 eV, = 3.5 × the α probe, i.e. the spin-averaged channel) and then rises as O p_z fills and its spin polarization reverses (+0.023 / +0.030 at 4.62 eV).  The basis difference grows from 0.0015 to 0.0074.  **PW at U = 1 eV (2026-10-31, PW-5 imposed): Δm_d −0.0035** (predicted −0.005; GPW −0.0075 / −0.0061).  So PW dips too, and the deficit builds up STEADILY across U: 0.0040 + 0.0114 (control, 0→1→4.62 eV), 0.0026 + 0.0055 (p5f3).  It is not the dip carried to large U.  The O p_z moment and the Mn d charge track PW at every U; only the Mn d MOMENT per unit O p_z change is short.  Next: a GPW Mn-basis ladder on Δm_d/ΔO p_z m at U = 1 eV (imposed, about 80 iterations a point); the control U = 2 eV limit-cycle recipe is now optional · `Records/HubbardUHistory.md` 2026-10-25/27/28/29 |

**Verdict.**  The U = 1 eV Mn-basis ladder (main 0fab43d6, AliasFree pinned; `~/Code/qchem6-runs/MnO2/uscan/ladder/`) measured the cross ratio Δm_d / ΔO p_z m: control 0.69, + O d 0.55, + Mn p5 f3 0.62, + both 0.52, + Mn p5 f3 + O d4 f1 0.45, PW 0.30.  It falls monotonically as O polarization is added; every other channel (O p_z moment, Mn d charge, band responses) is already on PW's.  The PW/GPW deficit is orbital-basis incompleteness in the Mn–O bond region, dominated by O polarization (C2v O site: dipole + quadrupole crystal field ⇒ O 2p mixes into d and f); no code path is implicated (user goal 2026-11-01: understand the discrepancy and rule out code bugs).  Full table: `HubbardUHistory.md` 2026-11-01.  Side finding → new §3 row "Dense Mn d midpoints make β-MnO₂ GPW DIVERGE".

## ✅ "rundeck cannot run a displacement path/pair as one SolidSession" -- CLOSED 2026-11-01 (lane B, 7af24457); the row as it stood

| row | why | next action |
|---|---|---|
| **rundeck cannot run a displacement path/pair as one SolidSession, so a user's two-point curvature sizes each point as a standalone Energy and can be silently 40 % wrong** (lane B, 2026-11-01; F: blocks forces/geometry from a deck) | pin 31 sizes a SolidSession on max|ΔH| and freezes the pilot's sizing across points (`SizedSetup`, `src/Calculation/SolidCalculation.C:651`), but only `src/Calculation/tests/SolidSession.C` drives it -- no deck key or postSCF action does.  MnO q7 VA ΔE(0.1 bohr): +1.12e-4 on the auto 40³ raster vs +1.565e-4 at densityEcut 160 (egg-box −3.9e-5 vs +3.7e-6) | a deck form -- `structure.path [displacements]` or a postSCF `displace` action -- that runs the points through ONE SolidSession (H-graded rasters + Becke sizing, frozen across points); gate: the MnO q7 VA pair from a deck matches the hand-sized 160+ Ha value.  The front end of the geometry optimizer will need the same path · `src/Calculation/Imp/Deck.C`, `src/Calculation/SolidCalculation.C` |

**Verdict.**  The postSCF action `displace` (`{"points":[{"atoms":[i,...], "delta":[x,y,z]}, ...]}`, `atoms` empty = every atom, the rigid
twin) re-resolves the deck's setup UNSIZED and runs the reference and each point through ONE `SolidSession`: its points consume
`MatrixElements`, so the pilot starts at the rigorous bound and the rasters + basis recipe are frozen across the points.  Pre-flight
refuses a free GPW run on the symmetric XC mesh, an atom past the structure, a zero delta, no points.  GATE (MnO AFM-II q7 VA, Γ,
the +m Mn +0.1 bohr along x, Becke (60, 2, 48) plain, the AUTO raster, `~/Code/qchem6-runs/MnO/msite/MnO_msite_q7va_fm_displace.*`):
the main run sized as an Energy on 60³; the session sized at the bound, x (19.35, 21.58, 19.35), 125³ / 135³; ΔE **+1.5522e-4**
(15 iterations, carried) vs the hand-sized ladder +1.5655e-4 / +1.5703e-4 / +1.5726e-4 (densityEcut 160 / 240 / 320) and CP2K
+1.5577e-4 (1200 Ry) -- the standalone auto-raster Energy pair gave +1.12e-4.  The rigid twin: **−4.4e-11** (the hand-sized egg-box
at 160 Ha was +3.7e-6).  Cost: SCF 1869 s at 4 threads, peak RSS 2.7 GB (a heavy-slot run).  Si_diamond plumbing check: one atom
+0.05 bohr +4.52e-4, rigid +1.2e-7.  Not reused: the main run's density (the session seeds afresh); a deck that only wants the
points still pays one Energy SCF first -- a `structure.path` form could skip it if that ever matters.


## ✅ Dense Mn d midpoints make β-MnO₂ GPW DIVERGE from iteration 1 (closed 2026-11-01, lane A)

**The row, verbatim:**

| **Dense Mn d midpoints make β-MnO₂ GPW DIVERGE from iteration 1 on main 0fab43d6 (regression)** (lane H, 2026-11-01) | `basis_scan/pol/p5f3_dD.bsd` (Mn p5 f3 + Mn d 0.56 … 24.66, vet off, CholeskyPivoted 1e-7 keeps 296/296, min pivot 4e-6) converged on 2026-10-26 (first iterate −93.10 Ha, AliasFree).  On main 0fab43d6 the FIRST iterate is +4.19e4 Ha at x = 8 (`accuracy` 0, rasterX/vlocRasterX/xcRasterX 8 -- the stated AliasFree) and +845 Ha at the default ε; 300 iterations end at +4.7e3 Ha (`uscan/ladder/p5f3_dD_U0.log`, `probe/dD_default.log`).  The same runs on control / Mn_p5f3 / p5f3_Od / p5f3_Od4f (no dense d) reproduce the old numbers to 1e-7 Ha.  So a Hamiltonian term breaks on a tight d set (the 24.66 Mn d is just above `routingBeta` 24); the ladder table carries the extra 144 Ha rung in both.  Suspect C2 step 3 (ε-sized rasters) or a commit between 2026-10-26 and 0fab43d6; an SCF whose first energy is +4e4 Ha should FAIL loudly (cf. the BallOnly row's "below any physical bound" ask) | bisect with `p5f3_dD_U0.json` at NMaxIter 1 (first Etot −93.10 vs +4e4 is the gate; ~1 min a step) between the 2026-10-26 run's commit and 0fab43d6; then the term-by-term energy split of iteration 1 · `src/BasisSet/Gaussian/Lattice/Imp/GPW_Evaluator.C`, `~/Code/qchem6-runs/MnO2/uscan/ladder/` |

**Verdict: the RASTER, amplified by the overlap's conditioning -- not a Hamiltonian term, not a commit.**  The basis is the same 296 functions
in both runs (not the reader); λ_min(S) = 6.75e-7.  First iterate vs the raster alias exponent x (90³ unless noted; `~/Code/qchem6-runs/MnO2/dDregress/`):
x 8 (64³) +41888 Ha (Eee +38770, Ekin +3440); x 11 −80.75302; x 13 / 16 −80.75311.  The 2026-10-26 run sat on AliasFree, which is x ≈ 16 on
this cell, NOT 8: AliasFree sized G² = 4·p_h·ln(1e6) with p_h = pβ/(p+β), ≈ 13.8p under HartreeXC routing (β = p/3), rounded up by the FFT
menu -- so lane H's "x 8 = the stated AliasFree" deck was half its alias exponent (doc/GridEpsilonPlan.md row 1, corrected).  At the default ε the
first SCF starts on the ball's own rasters (x ≈ 4: Eelec +923 Ha) and the pilot (C2 step 3) only audits after convergence, which never came.

**Fix (lane A):** (1) `SCFIterator` REFUSES a first iterate past its layout's bound (`FirstEnergyRefusal`): a FINITE system's energy above the structure constants, E − Σ{Constant} (rest mass and −TS included), above 0 -- the electrons-at-infinity theorem (the TOTAL can be positive, e.g. compressed H₂); a PERIODIC cell, which has no such limit (Ewald: E_el +1.59 Ha on healthy Si Γ; even E_tot is +0.047 Ha on healthy Na q1), only on a stated SCALE, E_tot above +1 Ha per electron (here +1102 / +22).  The first commit be857a17 used `GetElectronicEnergy` (Kinetic + Potential) for both: 93 short tests refused (every GPW cell, and DHF, whose rest mass that sum leaves out); fixed in the next code commit -- `SCFFailure::Why::UnphysicalEnergy`, both facades, the message names the raster/conditioning cause.  (2) `SolidCalculation::Converge`: a refused GPW run on UNRESOLVED rasters restarts ONCE from the
seed on rasters with e^(−x) ≤ λ_min(S), x = max(ln(1/λ_min), the ball's) -- 14.21 here; a STATED raster is never overridden.  Measured: the
stated-x 8 deck refuses at iteration 1; the default deck refuses, restarts at x 14.21, converges in 24 (Etot −93.4419391947) and the pilot SETTLES
with 0 re-convergences (every raster term meets 6e-5 Ha).  Gates: `GPW_MnO2.Γ_U_Shub_Pol_Smear_Kerker_Pulay_RefusesAnAliasedFirstIterate_Long`,
`…_ConvergesThroughARefusedFirstIterate_Panel_Long` (the basis vendored as `decks/panel/MnO2_p5f3_dD.bsd`).  Lane H's pilot ask (45³ vs 90³): on
this cell the pilot settled on 90³ at x 14.21 after the restart; the unrefused case is the remainder row.  pybind `why_name` lacks the new Why
(flagged to lane G).

## ✅ "The linker does not enforce library dependencies: no `-Wl,-z,defs` on our shared libraries" -- CLOSED 2026-11-01 (lane B); the row as it stood

| row | why | next action |
|---|---|---|
| **The linker does not enforce library dependencies: no `-Wl,-z,defs` on our shared libraries** (user, 2026-11-01; S) | CLAUDE.md says "the linker enforces no cyclic dependencies between libraries"; today only CMake does, and only for DECLARED `target_link_libraries` cycles among shared libraries.  A library that calls into another without declaring it links fine (undefined symbols are allowed in a .so by default) and is caught only, if at all, when an executable links.  C++20 module imports catch most of it, but not all (a non-module symbol, a template instantiated across libraries) | next free lane: add `-Wl,-z,defs` to every project shared library (not submodules), fix each undeclared dependency it reports by DECLARING it (never by allowing undefined), full allTests build + short/long sweep (link flags only: 0 anchor moves expected); then make the CLAUDE.md sentence true as written · root `CMakeLists.txt` |

**Verdict.**  The root `CMakeLists.txt` (end of file) walks the project's own trees -- the root, `src/`, `IntegrationTests/`,
`CLIapps/` -- for `SHARED_LIBRARY` targets and adds `LINKER:-z,defs` to each, so a new library is covered with no list to keep;
submodules and `pybind/` are untouched.  NOTHING was undeclared: all 26 libraries (qchem, the qc* family, qcPW_TestFields and the two
IntegrationTests harnesses) linked clean on the first try, and again on a full clean rebuild after merging main 57912b4a (lane A's
Fold / xcMeshAudit / `GPW/LiF.C` code included); `ldd -r` reports no undefined symbol in any; a throwaway .so calling an undefined
function fails under the same flag (the flag bites).  Sweep: short 1201/1201, long 5/5, 0 anchor moves vs main-1d19aad6 (a link flag
cannot move numbers).  The CLAUDE.md sentence "the linker enforces no cyclic dependencies between libraries" is now true of
undeclared dependencies too; a sanitizer build (its runtime resolves only in the executable) would have to drop the flag.

## ✅ "`multiplicity: 1` on an ODD electron count silently runs a DOUBLET" -- CLOSED 2026-11-01 (lane B); the row as it stood

| row | why | next action |
|---|---|---|
| **`multiplicity: 1` on an ODD electron count silently runs a DOUBLET** (lane H, 2026-11-01) | `src/Calculation/Imp/SolidSetup.C` DECISION 2: `twoS = multiplicity>1 ? multiplicity-1 : nel%2`, so 1 on odd Nelec gives nUp=nel/2+1, nDown=nel/2 and the parity check ("catches a singlet asked of an odd electron count") can never fire.  The `SolidSetup.C` doc says 1 = the EXPLICIT two-channel singlet (nUp=nDn).  Bit bcc Li (1 e/cell, kT 1e-3, per-spin Fermi levels): the run went fully spin-up (moment 1.1, ↓ band empty 0.13 Ha below the ↑ Fermi level, MOM-guard NON-AUFBAU) and printed Etot −0.2196 Ha against the unpolarized −0.2727 -- 53 mHa, reported CONVERGED (`~/Code/qchem6-runs/Li/mu/badspin/`) | decide what 1 means on odd Nelec: REFUSE with the parity message (CLAUDE.md: a user-reachable choice is an Outcome/throw, never a substitution), or (metal) require `spinsShareFermi` so the counts float; add the case to `src/Calculation/tests/Deck.C` beside the siteSpins rows |

**Verdict.**  REFUSED (CLAUDE.md: a user-reachable impossible request is a refusal, never a substitution).  Decision 2 is now one
exported function, `SpinChannelCounts(opts, nel)` (`src/Calculation/SolidSetup.C`): multiplicity 1 is the explicit two-channel singlet,
nUp = nDown, and on an odd count it throws naming the count and the three ways out -- multiplicity 2 (the doublet), 0 (minimal spin),
or `spinsShareFermi` true (one Fermi level sets the channels, so the counts are only the starting split and keep nel%2 there).  The
parity / 2S > N check now reaches multiplicity 1 too.  `deck::Resolve` calls it on the RESOLVED Nelec, so `rundeck --check` refuses
the deck before any SCF.  UT `Deck.MultiplicityOneOnAnOddCountIsRefusedUnlessTheSpinsShareAFermiLevel` (NaF: 8 electrons pass, 7
refused naming "odd Nelec 7", 7 with spinsShareFermi / multiplicity 2 / minimal spin pass, a doublet of 8 and 2S = 9 hit the parity
message).  No shipped deck asked for it (`scripts/vet-decks`: 52 of 56 OK, the 4 failures are the panel-deck basis-path row).



## ✅ THREADED multi-k run HANGS FOREVER: a data race on `DeltaScalarFitter`'s lazy `mutable` caches (closed 2026-11-02, lane A, 0e99a8a4)

**The row, verbatim:**

| **THREADED multi-k run HANGS FOREVER: a data race on `DeltaScalarFitter`'s lazy `mutable` caches** (F; found 2026-11-02) | lane H's `Li_bcc_q1_k8_u` (bcc Li, k8, Becke XC, QCHEM_OPENMP_THREADS=8) printed the header and sat at iteration 1 for 10 h at 200% CPU.  gdb: two threads spinning in `std::map::find` → `_M_lower_bound` on the SAME `IntCache` (`src/Fitting/Internal/DeltaFunctionFitter.C:181` `Integrator`: unlocked find-then-emplace), called from `tCompositeWF::AssembleFock`'s parallel block loop via `tVxc_Quadrature::MakeMatrix` → `DensitySampler_Singles` → `Overlap`/`Contract`.  A concurrent `emplace` corrupted the red-black tree, so lookups never terminate.  AssembleFock's contract ('the serial cold pass warms every cache; later passes do density-derived work only') does not hold for this cache: the Becke-quadrature integrators are first touched in a parallel pass.  Serial sweeps can never see it (every sweep is serial) | fix: guard `Integrator`/`Handle` (and every lazily-filled `mutable` map reachable from AssembleFock) with a mutex, or warm them in the cold pass, and ASSERT the contract (a cache miss in a parallel pass throws); then a 2-thread multi-k Becke UT that would hang/race (run under TSan once); audit the other `mutable SymMap` caches in `src/Fitting` · gdb stacks in the controller session (2026-11-02) |

**Verdict.**  The row's diagnosis held, with one correction to its cure.  The XC sampler builds ONE `DeltaScalarFitter` per field value v
(`tSinglesDensitySampler::FitFor`) and every Bloch block of the parallel Fock pass contracts against it, so the fitter's `Integrator`/`Handle`
maps are FIRST filled inside that pass: no serial cold assembly can warm an object born after it, which also makes the row's "make a
parallel-pass cache miss throw" FALSE here -- that miss is legitimate.  The guard is a lock: `DeltaScalarFitter::Integrator` (Handle is reached
only from it), and the same per-block find/emplace pattern in `ConstrainedFF::Integrator` (the Gaussian fitter) and `DeltaFit_IBS::Table` /
`Tensor` (two locks, tensor then table; Table never takes the tensor's).  Audit of the rest of AssembleFock's reach: the terms' `itsCache` /
`itsRealCache` are pre-created per block by the refresh (each block writes its own node, by design); the pair sampler's `itsAdj`, Hubbard's
projectors and the PP `itsByL` already lock.  AssembleFock's comment now says the warm-up is an optimisation, never the guard.  Gate:
`ResponseKernel.GPW_Si_k311_FockIsBitIdenticalAcrossThreads_Becke_Pol` (threaded pass FIRST, Becke × Pol = two fitters, the hanging
configuration).  TSan (`build/TSan`, the three threaded Fock gates): **0 warnings; with the fitter lock removed, 29** (the same pass and the
`DeltaFit_IBS` closures it reaches).  **THE SILENT HALF (lane H, 2026-11-02).**  Before the lock a threaded run did not only hang: bcc Li q1 k444 (imposed Oh, Becke XC, kT 1 mHa)
at 8 threads CONVERGED to Etot −0.272725 vs −0.252040 serial (−20.7 mHa, almost all E_xc: −0.20183 vs −0.17985; first iterate ρ_lost/N 1.0e-1
vs −2.7e-3), with no warning -- and the Si k311 Fock bit-identity gates passed through it.  On the locked build (rundeck, same deck): 1 / 2 / 8 / 8
threads −0.2520403163 / −0.2520403162 / −0.2520403180 / −0.2520403180, 10 iterations each: the same thread count bit-identical, different
counts within 1.8e-9 (a different collocation partition, row 160's thread-count remainder).  Gated to convergence by
`GPW_Li.k444_Becke_Imp_Smear_ThreadedRunIsReproducibleAndMatchesSerial_Long` (8 = 8 exactly, 8 vs 1 within 1e-8).

## ✅ "A molecule's mesh never reports itself: `MakeMolecularMesh` emits nothing to the run report" -- CLOSED 2026-11-02 (lane B); the row as it stood
| **A molecule's mesh never reports itself: `MakeMolecularMesh` emits nothing to the run report** (S; user 2026-11-02) | lane G found 'a molecule's real mesh size only appears once the run has started/finished'.  Cause: the PERIODIC Becke builder files its own `grids.becke` as it builds (`src/Structure/Imp/UnitCell.C:384`, the pattern we want: each section adds its own info to the CurrentReport as it runs, so it is there BEFORE the SCF), but the molecular/atom builder (`src/Structure/Imp/Molecule.C:39` `MakeMolecularMesh`, reached via `Structure::CreateIntegrationMesh` from Orbital_DFT_IBS/BasisSet/PP_Local/PP_NonLocal) files nothing.  User ruling: no 'one class tells another to report' (EmitReport-style) fix -- the builder reports itself | fix, per the user's design (2026-11-02): THE CALLER SETS THE CONTEXT, THE CALLEE REPORTS ONE LEVEL IN.  Each caller of `CreateIntegrationMesh` (the XC fit quadrature, PP_Local, PP_NonLocal, ...) opens a cursor scope naming ITSELF just before the call; the mesh builder writes its fields with cursor-relative `report::Set` (natom, radial, nRadial, angular, L, beckeOrder, points) and never chooses a key or an ID.  This is the cursor design already in `src/Common/Reporting.C:172`, but it has only a TOP-LEVEL `Section` and an array `Row`: add the missing nested-object scope (RAII, descends into `cursor[key]`).  WHY (user, 2026-11-02): with cursor-relative reporting the console report IS the run sequence -- unexpected nesting shows a peculiar run sequence at once (here: the Hamiltonian built inside the `basis` scope, Calculation.C:157, which root-anchored EmitAt calls had been hiding).  So: move the Hamiltonian build out of the basis scope; root-anchored `EmitAt`/`EmitSection` are the pattern being retired.  Then convert the periodic builder too: `UnitCell.C:384`'s root-anchored `EmitAt("grids","becke")` picks its own key, the same collision waiting to happen; UT: a molecule run's report has `grids.becke` before iteration 1 (an observer reads it at the first row) · lane G's run-spec notes |
**Verdict (lane B, 2026-11-02).**  Built to the user's design (caller sets the context, the callee reports one level in; scopes
purely RELATIVE, no root anchor).  `report::Scope(key)`, the nested-object RAII the cursor lacked: it descends into
`cursor[key]` and on close renders only the non-object members it added or changed, so an outer scope whose inner scopes
rendered stays quiet (UT `Reporting.AScopeDescendsOneLevelAndRendersOnlyWhatItChanged`).  The finding that decided the
shape: `Calculation::Converge` built the Hamiltonian INSIDE `rpt::Section basis`, and its comment said the root-anchored
`EmitAt` calls existed to escape that.  Now the molecule and atom facades build the basis in a silent `Scope("basis")`, the
Hamiltonian in `Scope("grids")`, and the iterator's conditioning rows in the rendering `Section("basis")` (the basis still
renders once, complete); the solid facade wraps its Hamiltonian build in `Scope("grids")`, and the grid audit's reference
Hamiltonian in `Scope("gridAudit" | "eggBox")` so it cannot overwrite the run's.  Every `CreateIntegrationMesh` caller opens a
scope naming itself (`xcQuadrature` in the three quadrature factories and the route's `kind`, `ppLocal`, `ppNonLocal` (also the
GPW KB mesh), `auxiliaryFit`, `residualMesh`, `dipoleProbe`); `MakeMolecularMesh`, `MakePeriodicBeckeMesh` and the uniform
cell mesh write `report::Set` fields and choose no key (`grids.becke` and `grids.unitCellUniform` are gone).  Paths measured on
H₂O LDA / PP at iteration 1: `grids.xcQuadrature`, `grids.auxiliaryFit`; PP: `basis.ppLocal`, `basis.ppNonLocal` (lazy, the
row's remainder).  UTs `MolecularMesh.TheMeshFilesItsFieldsOneLevelIntoTheCallersScope`,
`CalculationFacade.AMoleculesXCMeshIsInTheRunReportByIterationOne`.

## ✅ "CONTRACTED shells make the GPW gather crawl: CH₄ on MOLOPT DZVP spends > 35 min in one V_H gather" -- FIXED 2026-11-02 (lane B, 320affe8); the row as it stood
| **CONTRACTED shells make the GPW gather crawl: CH₄ on MOLOPT DZVP spends > 35 min in one V_H gather** (lane B, 2026-11-01; S, blocks the H/C/N gate) | CH₄ in a 15 bohr box, Γ, DZVP-MOLOPT-GTH (`scripts/cp2k2bsd`, 33 functions, 7 primitives per contraction down to α 0.0217): setup 247 s (V_loc-long audit on 96³), then the first `Vee_Hartree` gather (`ForShellPairBox` inside `IntegratePotential`) ran > 35 min at 4 threads without finishing (killed; gdb stack in the lane-B session).  Our own bases are single-primitive, so the box walk has never seen a 7×7 primitive pair with a 0.02 tail.  `~/Code/qchem6-runs/CH4/` | measure first (one gather, `perf`): per-primitive-pair boxes sized by the most diffuse primitive?  CP2K contracts on the grid per primitive pair with its own radius; compare the task count, then decide (per-primitive boxes, or a contracted pair's box from its own decay) · `src/BasisSet/Gaussian/Evaluators/PG_Cart_MnD/` (`ForShellPairBox`) |
**Verdict (lane B, 2026-11-02).**  Measured by reading, not perf: the separable-contraction kernel DECLINED every contracted
shell (`MakePairPoly` returns live = false when a radial has more than one primitive), so the ground-state collocate and gather
fell back to the reference walk `ForShellPairBox`.  It sizes ONE box from the shell pair's most diffuse primitives (MOLOPT
shares 7 exponents down to 0.0217 across every shell of an element, so every shell pair's box is ~24 bohr in reach and wraps
the 15 bohr cell's images: 16 209 tasks) and evaluates all 7 × 7 primitives at every point.  The transition collocation
(`ForTransitionBoxes`, R3) already summed primitive pairs through the cube kernel, each with its own box.
`ContractPrimitivePairs` / `GatherPrimitivePairs` give the ground state the same route (the shell's level kept, so the
numbers are the walk's to the screening tier).  Gate: `M_TransitionCollocation.WalkAndContractionRoutesAgree` now also runs H₂O
dzvp, walk vs cubes 4e-13 (ρ) / 1.3e-13 (h).  CH₄ DZVP-MOLOPT then converged (pilot 12 iterations, production 11).  The
remaining cost (~80 s/iteration, ~300 s V_loc audit) is the open P row: primitive pairs on their own grid levels.

## ✅ "A molecule has no build-only step: the `Calculator` constructor runs the whole SCF" -- CLOSED 2026-11-02 (lane B); the row as it stood
| **A molecule has no build-only step: the `Calculator` constructor runs the whole SCF** (S; 2026-11-02) | crystals have `deck::Build` (its own report run '<name> build', so grids/basis show BEFORE the SCF); a molecule does not, so a client cannot see the setup without running it.  Consequence found by lane G's report audit: `pybind/qchem_bridge.cpp` `qcb_basis_info` REBUILDS the basis outside any run (facade default options, no SALC blocking, no trim) to get S eigenvalues / nFunctions -- a second, possibly disagreeing, source for what the report's `basis.perIrrep` / `basis.removed` already record | user (2026-11-02): "the GUI needs a trigger to read all setup panels and start the actual work" -- i.e. TWO calls, Build then Run, for molecules as for crystals.  Add a molecular Build (resolve + build basis + Hamiltonian under its own report run, no SCF), the facade's `Converge`/Run then reuses it; the GUI reads every setup panel (basis, grids, PP) from the build run, then starts the SCF; GUI shape (user, 2026-11-02): a BUILD button beside Run -- Build does all the pre-SCF steps and fills the setup panels; Run alone is also fine (it builds if nothing is built, reuses a Build whose spec is unchanged, rebuilds after an edit); lane G then deletes `qcb_basis_info` and reads `basis.perIrrep` from the build run · lane G's audit (2026-11-02) |
**Verdict (lane B, 2026-11-02).**  `Calculation::Prepare(setup, onIteration)` and `AtomCalculation::Prepare(Z, charge, opts, acc,
onIteration)`, the twins of `SolidCalculation::Prepare`.  Both facades' `Converge` was split into `BuildGraph` (the accelerator,
the basis in a silent `Scope("basis")`, the Hamiltonian in `Scope("grids")`, the SCF iterator in the rendering
`Section("basis")`, whose wave-function build files basis.perIrrep / basis.removed) and the iterate.  `Prepare` runs
`BuildGraph` inside its own report run "<structure ID> build", closed before it returns; the FIRST `Converge` reuses the graph
(its own run then holds only the SCF), a later one rebuilds as before; the constructors are unchanged (they build inside the
converge run, as they always did).  Done together with the PP mesh's run sequence (same defect): `Ham_PP` builds ONE PP
quadrature in its own build phase and injects it into `PP_Local` and `PP_NonLocal` (it was built lazily per term per block
inside the SCF iterator's construction, so it filed under `basis.*`), now `grids.ppQuadrature`.  UTs
`CalculationFacade.PrepareFilesTheWholeSetupInItsOwnBuildRunBeforeAnySCF` (H₂O LDA + PP: the filed build run holds
basis.perIrrep and grids.xcQuadrature / ppQuadrature / auxiliaryFit, no scf header, the observer saw no row; Prepare + Converge
equals the constructor bit for bit), `CalculationFacade.AnAtomsPrepareFilesItsSetupBeforeAnySCF`,
`CalculationFacade.APseudopotentialRunsPPQuadratureIsInGridsNotBasis`.  Found on the way: the seed build files `scf.xcRhoRoutes`
(root-anchored, `DensitySampler_Singles.C:209`) into the build run -- listed on the mesh-report remainder row.

## ✗ RETRACTED "ninja does not recompile the importers of a changed module interface (stale objects)" -- 2026-11-02 (lane B, the reporter); the row as it stood
| **ninja does not recompile the importers of a changed module interface (stale objects)** (S, a FALSE-GREEN risk; lane B 2026-11-02) | lane B changed `Calculation`'s exported class layout; the first sweep showed 23 failures that were STALE OBJECTS: A_PP.C / L_PP.C (which would not even compile against the new interface) and the facade tests (segfault on the old layout) were never rebuilt -- the C++20 module dyndep missed them.  Touching the importers fixed it.  The danger is the converse: a stale object can also PASS, so a sweep can be green on code that no longer compiles | reproduce on a toy (change an exported class's layout, check `ninja -n` lists every importer); find whether it is CMake's dyndep scanning (the `.ddi`/`.modmap` of importers) or a clang-scan-deps cache; until fixed, a lane that changes an exported LAYOUT does a clean rebuild of affected targets before its landing sweep · CLAUDE.md 'ninja rebuilds were flaky' |
**Retraction (lane B, 2026-11-02): ninja was right; the BUILD FAILED and the sweep ran anyway.**  The build log of that sweep
(lane B's background job, never read past its tail) opens with `FAILED: IntegrationTests/CMakeFiles/ITMain.dir/A_PP.C.o ...
no matching constructor for initialization of 'Hamiltonian::PP_NonLocal'` (and `L_PP.C`): the importers WERE recompiled -- and
did not compile.  ninja stopped there, so ITMain and the other test exes were never RELINKED, and the sweep, chained after the
build with `;` instead of `&&`, ran the OLD executables against the NEW shared libraries: the "undefined symbol
...PP_NonLocalC1...MeshParams" and the segfaults on `Calculation`'s new layout.  Fixing the two test files and rebuilding is what
cleared it (the touch was incidental).  Checked both ways: a toy (CMake 4.2.3, ninja 1.13.2, clang 21.1.6; a library exporting a
class, an executable importing it; change the layout and the signature, build the library alone then the exe, or both at once)
recompiles and relinks the importer every time; on the real tree, touching `src/Hamiltonian/Internal/Terms.C` makes `ninja -n`
list every one of its importers (A_PP, L_PP, the PW and GPW files).  **The real false-green risk is a sweep over a failed
build**, and it is now refused where every sweep passes: `scripts/snaprun` runs a ctest command only when `ninja -n allTests`
has no work to do (exit 3 otherwise, naming the pending work; tested on a built and a deliberately stale tree).


## ✅ Accelerator knobs unreachable from the facades (closed 2026-11-02, lane A, 536c15f3 + the Ladder-default commit)

**The row, verbatim:**

| **Accelerator knobs unreachable from the facades: the Ladder's DIIS→GDM switch (molecules) and DIIS/ladder keys (crystals)** (S; lane G / user GUI pass, 2026-11-02) | molecules: `Calculation` hands `SCFAccelerators::Factory` only {NProj, EMax, EMin, SVTol}, so the Ladder's ethresh/stall/floor/switchat cannot be set; crystals: the deck has no keys for `SolidAcceleratorOptions` (diis.Nproj, SVTol, EngageBelowFD, ladder*), so the GUI refuses them | one option struct per facade carries the ladder fields (`AcceleratorOptions`, MoleculeSetup.C:97; `SolidAcceleratorOptions` + deck keys, round-trip UT), the facade forwards them; lane G wires the GUI · lane G's note |

**Verdict.**  Molecules: `AcceleratorOptions` gains `trust`, `ladderEThresh`, `ladderStall`, `ladderFloor`, `ladderSwitchAt`, `ladderSignal`;
`ToAcceleratorJson` is the one translation `Calculation` and `AtomCalculation` hand the factory (they passed only NProj/EMax/EMin/SVTol);
`scfrun` maps its flags onto it.  The json door's tail hand-off defaulted OFF (`switchat` 0) and a converging DIIS never stalls, so a molecular
Ladder WAS plain DIIS (lane G: water and ethene LDA/dzvp, GDM never engaged); the default is now |ΔE/E| < 1e-6 (`EnergyChange`, CLAUDE.md's
rule), gated by `SCFTrace.AMolecularLadderHandsOffToGDMBelowTheEnergyTail`.  Crystals: deck key `acceleratorOptions` = `SolidAcceleratorOptions`
(`diis` {nProj, engageBelowFD, svTol, svTolRel}, `gdm` {engageBelowFD, trust, trustBackoff, trustMin, pcFloor, aufbauTol, aufbauMaxRepicks},
`ladder` {gdm {…}, eThresh, stall, floor, switchAt, signal}), forwarded to `SolidSetup::Resolve`, written to every record; `ladderStartOnGDM`
stays the facade's per-stage decision.  UTs: `Deck.AcceleratorOptionsAreDeckKeysAndRoundTrip`, `MoleculeSetup.EveryAcceleratorKnobReachesTheFactoryDoor`.

## ✅ The two AUTO-ball rows -- CLOSED 2026-11-02 (lane A, the density-ball pilot + V_loc long sized from eps); the rows as they stood

| **The AUTO density cutoff under-resolves O q6 in an O-only cell** (2026-10-23) | `densityEcut` auto = 2·α_max of the ORBITAL basis: 36 Ha for an O box (VB O, α_max 18), while the O GTH-PADE-q6 local part is sharper (the log's own `alpha_pp=8.15`).  O atom in a 16-bohr box: 36 Ha is **+27 mHa** above the converged value; 72 Ha −0.99 mHa, 144 Ha +0.14 mHa vs CP2K.  Hidden in every Mn-containing cell (Mn sets α_max = 36 → 72 Ha).  **WORSE ON A SOFT BASIS (lane B 2026-10-31, the landing sweep):** O₂ in a 16-bohr box on sipp O (α_max 2) -- auto = 4 Ha, the grid audit prints UNDER (alias 3.9e-3 vs 1e-6) and the run carries on -- lands **0.711 Ha LOW**, non-variational (−30.71018); 60 Ha −29.998793, 120 / 240 Ha −29.998767 = CP2K −29.998761 to 6 μHa (`~/Code/cp2k-runs/O2/O2_sipp_box16_*`, same basis + PP + box; PySCF −29.998770).  `GPW_O2Box.Γ_Imp_M3_eqFinite` 'matched' it only through the old 30/5 molecular mesh (−30.736); it now states densityEcut 120 with a CP2K oracle line (a54bf84e).  ⚠ panel #4 `GPW_O2Box.Γ_Imp_M3_Pulay_Anchor_Panel` still anchors −30.71018 at auto (forwarded to lane A).  The fix is lane A's C2 (the raster sized from ε, honouring the audit's UNDER) | §3 ORDER step 0 FIRST: read the auto `densityEcut` each campaign cell printed (TiO₂, FeS₂, SrVO₃, NiO, MnO₂, spinels) -- if any is O-limited, this row lands BEFORE the re-runs, else it is independent.  Fix: fold the PP's own sharpness into the auto rule (max of the orbital and PP-projector/local widths), gate: O box at auto within 1 mHa of the 144 Ha value; then re-check every O-only / light-element anchor (O₂ box: `GPW_O2Box.Γ_Imp_M3_eqFinite`) · `~/Code/qchem6-runs/O_box16/` |

| **The auto density floor ignores a local-PP FIELD sharper than the basis** (2026-10-25; the REMAINDER of "⛔ GPW `BallOnly` raster + the C = 2 density floor ALIAS tight pairs", closed 2026-11-01 by C2 step 3 -- `Records/OpenWork_History5.md`) | the raster half is cured: C2 step 3 sizes every raster from ε (pin 31, the pilot audit), `BallOnly`/`AliasFree` retired.  Still open: the auto `densityEcut` floor follows the BASIS (max 2α), not a local-PP field sharper than it -- α = 2 alone is 0.03 Ha off even alias-free; O2Box after the re-bank is still 0.106 Ha above CP2K (−30.1051713 vs −29.998761, box16, 400 Ry), attributed to this | measure the O2Box gap against densityEcut at a fixed sized raster; if it closes, add the PP field's need to the floor (max(2α_max, the PP field's G-tail at ε)) and re-bank · `src/BasisSet/Gaussian/Lattice/Imp/GPW_Evaluator.C`, `IntegrationTests/GPW/Boxes.C` |

**Verdict.**  The density ball (the G-ball the density is collocated into) was never judged: the pilot's `x + 8` reference raised the
rasters at the run's ball, so its truncation sat on both sides and cancelled.  `SolidCalculation::CalibrateRasters` now runs a nested
ladder like the Becke mesh's: the run against twice the ball at the converged D (`CompareAtFixedD(..., ballScale 2)`), the ball doubling
(`cutoffFactor`) while the Hartree term is over pin 31's target, a stated `densityEcut` honoured and reported OVER.  Every other term the
audit reads is judged too (the honest gate, controller ruling): one over target that a bigger ball cannot cure blocks SETTLED.  That
exposed the V_loc long routing floor -- the per-pair harmonic rule's fixed `vlocEps` 1e-5 re-placed over the multigrid levels as the ball
doubled, alternating ±6e-5 Ha on MnO (target 4e-5) -- now sized from the target, `GridAccuracy::VlocEps()` = eps/10 (measured margin:
the delivered term error is 1.0-1.7x vlocEps per element, 1.4-6x per cell energy, at 1e-5 and 1e-7 alike).  Evidence: Si Γ geometry,
V_loc long row 1.4e-5 at vlocEps 1e-5, 3e-8..5e-7 at 1e-7; lane B's stated-ball MnO series 250 → 320 Ha drifted 169 μHa at 1e-5 and
0.75 μHa at 1e-7 (`~/Code/qchem6-runs/MnO/lanea_vloc7_e{250,320}`), so MnO's free deck has the reference −61.303112 (±1 μHa) and its
auto cell lands −61.3031131.  O₂ box16: 106 mHa → 0.11 mHa from CP2K.  Two pre-existing pilot defects fixed beside it: the convergence
tail estimate d/(1−d/d_prev) blew up when a round re-measured an unchanged rung (now only between two different rungs), and
`PilotCost::settled` read true at the round cap.  Cost (MnO free, 4 threads): energy 463 s / 1.0 GB vs 183 s / 796 MB at the old ball,
almost all the bigger ball (the per-SCF cost scales N_G); geometry 1h22m / 3.3 GB (§4 row: the XC-mesh sampler's per-point G-sum is about
half the growth).  Remainders: the O tail row, the Becke-mesh-on-MnO-geometry row, the MnO VA (+39 μHa) and NaF (+0.82 mHa) converged CP2K-gap rows (§3).

**The third row, closed 2026-11-03 at lane A's landing (main's §3 row as it stood, verbatim):**

| **The density G-ball and V_loc-long were never judged against ε: the AUTO `densityEcut` truncates every sharp pair by a fixed 1.8 %, and the fixed `vlocEps` 1e-5 scatters ±80 μHa** (F/V; lane B 2026-11-02) -- **fix in flight, lane A** | `densityEcut` auto = 2·α_max, so the audited ball tail of the tightest pair is e^{−4} = 0.0183 in EVERY cell (MnO prints `ballTail 0.0183156` at 72 Ha exactly as H₂O does at 24 Ha).  The pilot sized the rasters to ε (pin 31) but never the ball: `auditDensity` ran with `judgeBall=false` under the comment "Hartree's 1/G² kernel suppresses the ball tail … by design" (`GPW_Evaluator.C:1160-1163`), refuted here (XC sees the truncated ρ with no 1/G²).  **Molecules** (DZVP-MOLOPT-GTH, 15 bohr Γ, Becke (60, GL 48), vs CP2K 1200 Ry, auto 24 Ha → 100 Ha): H₂O **+330 → +2.5 μHa**, NH₃ +66 → +0.5, CH₄ +11 → −0.4 (CP2K converged: 1600 Ry +1.1 μHa; XC mesh (100, GL 64) 0.1 μHa; a 20 bohr box keeps +336 μHa, so not electrostatics; decks `decks/probes/{H2O,NH3,CH4}_dzvp_molopt.json`, 3628ddcd).  **MnO AFM-II free** (`decks/MnO_AFM2_free.json`, VALENCE_LOWQ_SR; 15 iterations, ±4.4587 every arm): auto 72 / 100 / 160 / 250 / 320 / 400 Ha = −61.30328078 / −61.30319112 / −61.30311021 / −61.30314294 / −61.30331206 / −61.30317248 -- from 250 Ha the ball tail is < 1e-6, yet E scattered ±80 μHa with SCF ΔE 1e-13 and Exc within 16 μHa (Ekin −340 μHa 250 → 320): the V_loc-long routing tolerance (fixed `vlocEps` 1e-5).  Lane A reran 250 / 320 at `vlocEps` 1e-7: −61.3031121 / −61.3031128 (0.75 μHa) ⇒ **reference −61.303112 ± 1 μHa; auto is 169 μHa low** (ball + routing together).  The msite row's q7 VA drift (160 / 240 / 320 Ha) is the same scatter.  Distinct from the O-q6 row above (PP sharpness; its fix would not cover MOLOPT H₂O, where the orbital pairs set α_max) | lane A: every term gated against ε (the ball judged, the pilot sizing `densityEcut`), `vlocEps` derived as ε/10 instead of fixed.  Gates: the three probe decks with `densityEcut` removed within 3 μHa of CP2K; `decks/MnO_AFM2_free.json` at default within ε of −61.303112; re-judge every anchor that moves (pin 10).  Close to history at lane A's landing · logs `~/Code/qchem6-runs/{H2O,NH3,CH4}/*_dzvp_molopt*.log`, `~/Code/qchem6-runs/MnO/ball_MnO_AFM2_free_*.log`, `lanea_vloc7_e*`, CP2K `~/Code/cp2k-runs/{H2O,NH3,CH4}/` |

⚠ **CORRECTION to that row's wording:** "refuted here (XC sees the truncated ρ with no 1/G²)" is WRONG on the mechanism.  XC reads ρ[D]:
under the cusp-deficit default that holds at the fixed point, and the linear D-mixing and GDM routes read it directly.  The ball enters
through the G-space quadratures, Hartree and V_loc long.  The `x + 8` audit missed it because its reference kept the run's own ball, so
the truncation sat on both sides and cancelled.  The false comment (`GPW_Evaluator.C:1160-1163`, and the DensityGridAudit doc comment) is
gone; CleanCode D-FALSE-CONTRACT (1) records it.

**Gates.**  `decks/MnO_AFM2_free.json` at default: −61.30311314, 1.1 μHa from the reference −61.303112 (ε 4e-5): MET.  The probe decks with
`densityEcut` removed (lane A runs `~/Code/qchem6-runs/<m>/lanea_AB_<m>_dzvp_molopt.{json,out}`, ε = 1e-5 default) -- the gate table below.  `vlocEps` is now `GridAccuracy::RigorousVlocEps()` (ε, the per-pair RELATIVE bound of the harmonic routing rule, no S prefactor)
divided by `vlocRoutingEta` (the measured calibration: max 6 per cell, 1.7 per element, × ≥ 1.5 margin, rounded to 10; deck
`solid.accuracyCalibration.vlocRouting`), combined in `VlocEps()` only; a stated `tolerances.vlocEps` wins; the run report prints
rigorous, calibration, result and source (`Deck.AVlocRoutingCalibrationReachesTheEvaluatorAndAStatedVlocEpsWins`).

| probe (DZVP-MOLOPT-GTH, 15 bohr Γ, Becke (60, GL 48)) | auto ball (old) | stated 100 Ha | **lane A auto** | ball settled | ball's own delivered (cf vs 2·cf) | wall / RSS (4 threads) |
|---|---|---|---|---|---|---|
| H₂O vs CP2K −17.1824875 | +330 μHa | +2.5 | **−13.3 μHa** (−17.18250083) | cf 4 (48 Ha) | 1.47e-5 | 1h06 / 0.75 GB |
| NH₃ vs CP2K −11.7016477 | +66 | +0.5 | **−0.75** (−11.70164845) | cf 4 | 1.25e-6 | 2h00 / 0.69 GB |
| CH₄ vs CP2K −8.0351296 | +11 | −0.4 | **−1.2** (−8.035130806) | cf 4 | 8.3e-7 | 2h45 / 0.79 GB |

NH₃ and CH₄ meet the 3 μHa gate; **H₂O does NOT**: it lands 13.3 μHa below CP2K, inside its own target (ε·N = 3e-5) and the size the
ball audit itself predicts (1.47e-5 at cf 4 vs 8; the stated-100 Ha run sits 15.8 μHa above it).  At the default ε the pilot is RIGHT to
stop there: the 3 μHa gate asks for ε ≈ 1e-6 on H₂O.  (Wall: the three ran serially beside a sweep; per-SCF cost at cf 4 is the §4
sampler row.)


## ✅ "No GTH valence basis covers H, C or N: organic molecules cannot take the pseudopotential path" -- CLOSED 2026-11-03 (lane B); the row as it stood

**The row, verbatim:**

| **No GTH valence basis covers H, C or N: organic molecules cannot take the pseudopotential path** (user via lane G, 2026-11-01: benzene + PP is refused) | the sipp* / valence_lowq* sets in BasisSetData cover O, Si, F, Na, Al, Mn, Li, Ni, V, Sr, Ti, Fe, Co, Cr, Cu, K, S only, so SOLID² refuses benzene with PP ticked (correctly: a refusal, not a substitution) | lane B: import MOLOPT DZVP/TZV2P-GTH for H, C, N (and check O against the existing set) into BasisSetData (a CP2K → .bsd converter, the reverse of bsd2cp2k); gate: benzene or CH₄ on the PP path matches CP2K with the same basis + GTH PP.  GUI side (lane G): the refusal names the uncovered elements.  **2026-11-01 (lane B):** `scripts/cp2k2bsd` (88a89bca) and H q1 / C q4 / N q5 in `gth-lda-valence` + their seeds (f21ad3c6) are built; **2026-11-02, GATE PASSED for H and C:** CH₄ in a 15 bohr box, Γ, DZVP-MOLOPT-GTH + GTH-PADE q1/q4, spherical -- CP2K −8.0351296 (1200 Ry; 800 Ry −8.0351292), ours −8.0351407 on a converged Becke XC mesh ((60, GL 48) −8.0351412, (60, GL 64) −8.0351412, (100, GL 48) −8.0351407): 11 μHa absolute, 2 μHa/atom, 8 μHa of it in Exc.  The default Becke (40, site-adapted 29) is 52 μHa short on CH₄ (−8.0350893).  Committed: `BasisSetData/dzvp_molopt_gth.bsd` (H C N O, the gate in its preamble), `decks/probes/CH4_dzvp_molopt.json`.  **N and O GATED 2026-11-02 at a stated densityEcut 100 Ha** (3628ddcd): CH₄ −0.4, NH₃ +0.5, H₂O +2.5 μHa vs CP2K; the old CH₄ 11 μHa and the auto-cutoff gaps were the density ball (§3 row 'The density G-ball and V_loc-long were never judged').  **TZV2P-MOLOPT GATED 2026-11-02** (f4060639, `BasisSetData/tzv2p_molopt_gth.bsd`, same recipe): CH₄ −0.4, NH₃ +0.6, H₂O +2.6 μHa vs CP2K.  **BENZENE GATED 2026-11-03** (af3b9999, 18 bohr box, same recipe): −37.6867871 vs CP2K −37.6867742, −13.0 μHa = 1.1 μHa/atom (Exc −9.4 of it) · `BasisSetData/`, `scripts/bsd2cp2k`, `scripts/cp2k2bsd` |

**Verdict.**  DZVP-MOLOPT-GTH and TZV2P-MOLOPT-GTH for H C N O are in `BasisSetData/` (`scripts/cp2k2bsd`), with H q1 / C q4 / N q5 in
`gth-lda-valence`.  Same-basis gates vs CP2K 1200 Ry (15 bohr Γ, GTH-PADE, spherical, Becke (60, GL 48), stated densityEcut 100 Ha):
DZVP CH₄ −0.4, NH₃ +0.5, H₂O +2.5 μHa; TZV2P CH₄ −0.4, NH₃ +0.6, H₂O +2.6 μHa; benzene (DZVP, 18 bohr) −13.0 μHa = 1.1 μHa/atom.
Probe decks `decks/probes/{CH4,NH3,H2O,C6H6}_dzvp_molopt.json`, `{CH4,NH3,H2O}_tzv2p_molopt.json`.  The auto densityEcut gap
is the open §3 density-ball row (fix in flight, lane A).  Remainder (lane G): the SOLID² refusal names the uncovered elements.

## ✅ "The cfg column flags occupation jumps on CH₄ almost every iteration" -- FIXED 2026-11-03 (lane B, b0b78a74); the row as it stood

**The row, verbatim:**

| **The cfg column flags occupation jumps on CH₄ almost every iteration** (user via lane G, 2026-11-03) | CH₄ run in C2v: t2 → A1 + B1 + B2, exactly degenerate; the `*` (occupation jump) marks most iterations, the user expects none.  Lane G's reading: reordering among degenerate cross-irrep levels flips the per-irrep occupation.  ⚠ VERIFY FIRST: CH₄ is closed-shell a1² t2⁶, so all three t2 partners are FULLY occupied and a swap among them cannot change per-irrep occupations -- if `*` still fires, the detector is comparing something else (level ORDER, or a degenerate VIRTUAL pair).  CONFIRMED 2026-11-03 (lane G, from the user's console orbital dump): the run was neutral, UNPOLARIZED HF, so a spin-polarized run is ruled out | reproduce on the CH₄ probe and print what the detector compares; then either (a) the detector ignores swaps within the degeneracy tolerance (pin to `MergeTol`, see the ε inventory §4), or (b) the aufbau orders degenerate cross-irrep levels stably.  A real occupation change must still be flagged · `src/WaveFunction/` (occupation / aufbau), `src/SCFIterator/Imp/SCFIterator.C` (the cfg cell) |

**Verdict.**  Reproduced with `scfrun --mol CH4 --model HF --sym` (C2v; `--mol CH4` and `--sym` added for it): `*` on most rows.  What
the detector compared: an ENERGY-ORDERED label string, e.g. `(1a₁)²(2a₁)²(3a₁)²(1b₂)²(1b₁)²` → `(1a₁)²(2a₁)²(1b₁)²(3a₁)²(1b₂)²` --
the three full t2 partners (3a₁, 1b₁, 1b₂) reordered on noise, no occupation moved.  Neither reading in the row: it was LEVEL ORDER
among occupied levels.  The O₂ triplet golden's two `*` (iterations 2-3) were the same artefact, an ↑ level crossing a ↓ level.  Fix
(neither (a) nor (b)): the cfg compares the aufbau's own decision, per irrep block (wave-function block order) each orbital's rounded
occupation in block order (`OccupationsByBlock`, `SCFIterator.C`); the merged level table is not used (a MergeTol host is whichever
level came first).  Tests: `SCFTrace.ADegenerateClosedShellAcrossIrrepsNeverFlagsCfg` (CH₄ C2v, no row flags),
`SCFTrace.AnElectronMovingBetweenIrrepsFlagsCfg` (NO doublet: the π* electron moves b₁ ↔ b₂ on iterations 2-3 and still flags).
Notes: no natural run tried flags now (N₂, CO, O₂ triplet / singlet, C₂, B₂ with symmetry; atoms B C O Al Sc V Cr Fe Cu polarized),
and a kT = 0 metal does not hop either: the solid fills every k-block with the same count (Na q1 fcc k333, 1 electron per block).


## ✅ "GDM never engages on a symmetry-blocked molecule whose electrons span several irreps" -- FIXED 2026-11-04 (lane B, da50a531); the row as it stood

**The row, verbatim:**

| **GDM never engages on a symmetry-blocked molecule whose electrons span several irreps** (user, 2026-11-03, CH₄ C2v HF via the GUI) | the ladder hands DIIS → GDM at |ΔE/E| 1.9e-9 and GDM declines in EVERY block: "the occupied density is NOT this minimizer's leading 5 orbitals (Tr(D'P_block)=6)" (a₁), "=2" (b₁, b₂), then retreats to DIIS for the rest of the run.  Cause (read, `SCFAccelerator/Internal/Imp/SCFAcceleratorGDM.C:69,99`): each irrep's GDM gets `occ` = the RESERVOIR's electron count (10, the whole molecule: the aufbau floats electrons across irreps), so `itsNocc = itsNel/g` = 5 in every block, while the a₁ block actually holds 3 doubly-occupied orbitals and b₁/b₂ one each.  The check is right; its expected count is wrong.  ⚠ The message also names the WRONG cause ("a MOM or Fermi-smeared occupation") -- neither is present (CleanCode D-FALSE-CONTRACT candidate).  The hops the user suspected were the cfg-column artefact (fixed 6aa96d65), unrelated | take the block's occupied count from the IDEMPOTENT D' itself (`itsNocc = round(Tr D' / g)` when `idem`), not from the reservoir; keep the leading-block check (it still catches MOM / smearing); correct the message to say which count it expected and from where.  UT: CH₄ C2v HF engages GDM at the hand-off and converges; a MOM run still declines.  Sweep + panel (SCF path) · `src/SCFAccelerator/Internal/Imp/SCFAcceleratorGDM.C` |

**Verdict.**  Three defects, the second and third hidden behind the first:
1. **The block's occupied count** (the row's cause): `itsNocc` came from the Create `occ` = the RESERVOIR count (10 for CH₄).  Now
   `UseFD` reads it off the idempotent D' (`itsNocc = round(Tr D'/g)`, beside the capacity g it already read); the reservoir count
   is no longer used.  The leading-block check stays and still catches a MOM (non-aufbau) pick; its message now says what it
   compared ("Tr(D')/g = N occupied, but Tr(D'P_leading) = x of g·N") and names the non-aufbau pick as the cause, not "MOM or smearing".
2. **An empty (or full) block vetoed the step**: CH₄ C2v has an empty a₂ block, and the composite steps only when EVERY block steps.
   A block with nocc = 0 or nocc = n now steps PASSIVELY (diagonalise F', no rotation), and a composite RejectStep counts only the active blocks.
3. **The held line-search fill used the reservoir count** (latent; reached for the first time once 1 and 2 let GDM engage): the
   trial fill `HeldOccupationPolicy` put PrescribedBlockN (10) in every block, giving E +55.8 Ha and "Matrix sizes do not match".
   The orbitals cannot supply the count because `UpdateOrbitals` rebuilds them unoccupied (reading them gave a dead-gate [F,D] = 0).
   `CompositeWF` now records each block's electron count after every unheld fill (`itsLastBlockN`), and a held fill reuses it.
   A held fill with no record throws.

**Result.**  CH₄ C2v HF (`scfrun --mol CH4 --model HF --sym`): the Ladder hands DIIS → GDM at iteration 8-9, stays on GDM to
convergence with the energy non-increasing, E = −40.194456980543 (DIIS-only −40.194456980544), and no DECLINING line.
Tests: `SCFTrace.AMolecularLadderEngagesGDMOnASymmetryBlockedMolecule` (Ladder vs DIIS on CH₄ C2v: all GDM after the hand-off,
energy monotone, |ΔE| < 1e-8); `SCFAccelerator.TheBlockCountIsReadOffTheDensityNotTheReservoir`,
`.AnEmptyBlockIsPassiveAndDoesNotVetoTheStep`, `.AMomOccupationStillDeclines` (D' on orbitals 0 and 2 still declines).
Sweep: short 1263/1263, long 10/10, anchordiff 0 moved.  Panel+periodic 24/25: the one red is TiO₂ k222 at 62 vs a budget of 49, the
open pilot-budget row (66 on main-4b4a6641).  Its accelerator is Null, so this change is not involved.

**The CLAUDE.md claim "GDM now engages on unpolarized (occ=2) runs too" (2026-10-01) was verified on single-block runs only**:
`WaterPureGDM` (no symmetry, energy only, the hand-off "checked by hand"), the folded-doublet UT (one `Irrep()` block), and B UHF.
Before da50a531, no symmetry-blocked molecule with electrons in more than one irrep had ever run GDM.


## The ρ̃ sampler, step 1 (lane A, landed ac66be2d, 2026-11-04) -- the design note, verbatim

### Lane A — ρ̃ sampler design note (OpenWork §4 step 1), 2026-11-04

lane-a = main 3400d501 merged (clean, no conflicts).  Design only; no code yet.

#### What is actually being sampled (the call path, measured)

Both of the MnO-geometry buckets are the SAME kernel, `FourierMixCD::operator()(const rvec3vec_t&)`, a direct
per-point G-sum (folded half-space, per-axis phase tables, one complex MAC per folded G per point):

| bucket (`lanea_ballcost_geometry.out`) | calls | s/call | total | what is sampled |
|---|---|---|---|---|
| `XC-mesh rho sampling (matrix-free density)` | 122 | 11.6 | 1415 s | the WHOLE mixed ρ̃, both spins (`DensitySampler_Singles.C` matrix-free branch: a mixed density with no exact source) |
| `XC-mesh rho sampling (all iterations)` self | 281 | 2.45 | 689 s | ρ[D] by the Φ GEMM **plus** the N4 cusp-deficit correction ρ̃_mix − ρ̃[D], which `FieldMixer::Present` builds as a `FourierMixCD` over the same ball and `FieldProjection().Project(*corr)` G-sums at every Becke point |

N_r = 31 664 Becke points (after the tail drop), and N_G is the density ball (72 → 160 Ha takes the cost from 1.82 to 6.51 s/call, ∝ N_G).
So cost = N_r · N_G/2 complex MACs per channel, which is O(N_G·N_r).

Side question to settle while I'm in there: **why 122 calls take the matrix-free branch at all**, since the cusp-deficit route should give every mixed
density an exact source.  Candidates: Pulay-extrapolated densities, or the pilot's audit/ladder SCFs.  If it's the Pulay path, its exact source
is missing and that's a pin-28 question (XC reads ρ[D]).  I'll report it as a finding and won't fold it into the sampler change.

#### Route: type-2 NUFFT.  Reusing the levels' collocated ρ is REJECTED

**(b) The levels' collocated ρ: rejected.**
- The multigrid levels collocate ρ[D] (the pair scatter) on RASTER nodes.  What the XC mesh needs is the MIXED field, or its correction,
  at BECKE points, which are off-raster.
- Getting a raster value onto an off-node point is an interpolation.
  - The only exact interpolant of a band-limited raster is the G-sum itself.
  - A local polynomial or Lagrange interpolant has an error that is not band-limited and has no ε-bound: it rings at the cusps, which is pin 1's
    wrecking ball.
- The levels also don't hold ρ_mix: the mixer's field lives in G.  So (b) is either inexact or collapses into (a).

**(a) Type-2 NUFFT: chosen.**
ρ(r) = Σ_dm c(dm) e^{2πi dm·s}, with s = fractional coordinates of r (G = B dm, so G·r = 2π dm·s).  That's a standard 3-D type-2 NUFFT on the
integer lattice dm, for ANY cell (triclinic included), with no Cartesian assumption.  Steps:
1. **Deconvolve:** c'(dm) = c(dm)/φ̂(dm).  O(N_G).
2. **FFT** onto an oversampled grid n_a = σ·(2·m_a+1), σ = 2.  This uses pocketfft through `qchem.FFT` (`src/Math/FFT.C`), with real output because ρ is real.
3. **Interpolate:** each point gathers the w³ grid values nearest it, weighted by the SEPARABLE kernel φ(s_1)φ(s_2)φ(s_3).  The per-point weights are
   3·w numbers and are cached per mesh: the Becke points are fixed for the run, and so is the ball.

Kernel: the "exponential of semicircle" φ(z) = exp(β(√(1−z²) − 1)), β = 2.30·w (FINUFFT's kernel).  The error is ≈ 10^{1−w} relative to ‖c‖₁.
The Gaussian kernel is the fallback: a closed-form bound, but about 2× the width for the same digits.

Pin 1 check:
- Nothing truncates a lattice sum.  The output is the SAME band-limited series as the G-sum, to a stated ε.
- The kernel's compact support is a property of the kernel by construction, not a cut of a field.  Its error is its Fourier alias tail, and that
  is exactly what the bound measures.
- It adds no Gibbs content beyond what the series already has.  The cusp-deficit route already keeps the cusps on ρ[D] (the Φ GEMM) and samples
  only the smooth correction.

#### Dispatch: the point set answers for itself (added 2026-11-04, user requirement)

**A fact first, because it changes what "untouched" can mean.**
- Today the uniform cell mesh has NO inverse-FFT path for a ρ̃-mixed density.
- `DeltaFit_IBS::Overlap(f)` (`src/BasisSet/DeltaFit_IBS.C:144`) calls `f(mesh.Points())` on BOTH meshes, so the uniform XC mesh
  (`UnitCell::CreateIntegrationMesh`, cellKind Uniform: n³ fractional MIDPOINTS) went through the same direct G-sum as Becke.
- The paths that ARE inverse FFTs are the PW fit (`PlaneWaveFit_IBS::Overlap`, G-space) and the raster XC route (`GetRhoOnGrid`).
  Neither calls `FourierMixCD::operator()(points)`, so neither is touched by this change.
- So the uniform mesh GAINS the exact inverse FFT here.  It is not "kept": it did not have one.

**Who holds what.**
The sampler's client is the δ fit basis (`DeltaFit_IBS`, via `FieldProjector::Project` → `OrthogonalFit` → `Overlap(f)`).  It holds a
`FitQuadrature` (= `qcMesh::FoldedMesh`), which holds a `shared_ptr<const qcMesh::Mesh>`.  It does NOT know the mesh kind, and it must not.
The field it is handed is a plain `ScalarFunction<double>`.  So the dispatch needs two halves: the POINT SET makes the sampler, and the
FIELD chooses whether to use one.

**The design (two virtuals, no kind switch anywhere):**
1. **qcMath, a new data-free face `PointSampling`:**
   - `virtual const rvec3vec_t& Points() const = 0;`
   - `virtual std::unique_ptr<FourierSampler> MakeFourierSampler(const FourierBox&) const;`
   The DEFAULT returns the NUFFT sampler on `Points()`.  `FourierBox` = the reciprocal basis B (G = B·dm) and the dm index box.
   `FourierSampler::Sample(terms) → rvec_t` is values on MY points, in MY order.
2. **`qcMesh::Mesh` implements `PointSampling`.**
   - The uniform raster's override cannot be a subclass: `CreateIntegrationMesh` returns `Mesh` BY VALUE, and copies are everywhere,
     so a subclass slices.
   - So the override travels as a SAMPLING STRATEGY OBJECT the BUILDER installs.  `Mesh::MakeFourierSampler` forwards to it, and a
     null strategy means the default (NUFFT).
   - `UnitCell::CreateIntegrationMesh`, the uniform branch, installs `RasterFourierSampling{n, half-voxel offset}`.  That strategy
     returns an inverse-FFT sampler.
   - The inverse FFT is EXACT on the midpoint raster for ANY n: e^{2πi k·(j+½)/n} is n-periodic in k, so each k is accumulated into
     bin (k mod n) with its half-voxel phase.  It is the same series, no approximation, and no σ or w.
   - Becke and every other builder install nothing, so they take the NUFFT default.
3. **The field's half, `ScalarFunction<double>`:** `virtual rvec_t Sample(const PointSampling& ps) const`.  Its default is
   `(*this)(ps.Points())`.  `FourierMixCD` overrides it with `ps.MakeFourierSampler({B, box})->Sample(itsRho)*itsScale`.  So nothing
   ever asks a mesh what kind it is, or asks a field whether it is a Fourier series.
4. **The client:** `DeltaFit_IBS::Overlap(f)` calls `f.Sample(*mesh)` instead of `f(mesh->Points())`.  That is a one-line change, and
   every non-Fourier field takes the default, bit-identical to today.

**Dependency direction:**
- qcMath holds `PointSampling`, `FourierSampler`, the NUFFT and the FFT.
- qcMesh implements the face (qcMesh → qcMath, as today).
- Structure (`UnitCell`) installs the raster strategy (it already builds the mesh).
- BasisSet (`DeltaFit_IBS`) calls the face.
- ChargeDensity (`FourierMixCD`) overrides the field's half.
- Nothing points upward, and `FourierMixCD` names only qcMath faces.

**Design point for the user:** `Mesh.C` documents itself as "a plain value type with no virtuals ... No polymorphism".  This gives it
one abstract face (`PointSampling`) and one held strategy.  The strategy keeps it a value type (copies keep their sampling, nothing
slices), but it does gain a vtable.  The alternative, to put the face on `FoldedMesh` (what `DeltaFit_IBS` holds), only moves the
vtable, because the raster layout still has to be recorded where the points are built.

**Gates (unit tests):**
- NUFFT vs the direct G-sum on a Becke mesh: the `KerkerMix.BatchedSamplerEqualsTheDirectSeries` pattern, moved onto a real
  `MakePeriodicBeckeMesh` point set.
- The uniform raster sampler vs the direct G-sum on `CreateIntegrationMesh`'s uniform mesh, to round-off.  It must also be the FFT
  sampler: a counting stub sampler proves the strategy, not the NUFFT, answered.
- **"The uniform path is untouched (bitwise)" holds for every NON-Fourier field** (default `Sample` = `f(Points())`, bitwise).  For a
  ρ̃-mixed field the uniform mesh changes from the direct sum to an exact FFT, so equality there is to round-off, not bitwise.  A bitwise
  uniform gate for ρ̃ is only possible if the uniform mesh keeps the direct sum.  Say if that is what the user wants.

**Where step 1 stands (honest):**
- I had coded step 1 before this requirement arrived.  The NUFFT now sits INSIDE `FourierMixCD::operator()(points)`, so it serves
  both meshes, and the sweep of that is running now.
- The NUFFT module (`src/Math/NUFFT.C`) and its unit tests carry over unchanged.  Its width rule is verified to deliver its tolerance
  at w = 4..14, after a sizing fix: n is sized from max|k|, not from the box's extent.
- The dispatch is a restructure of where the call is made, so I'll rework it into the design above and resweep before committing.
  The current sweep stays as evidence that the NUFFT delivers 0 moves.

#### Accuracy against ε (pin 31, CLAUDE.md "rigorous ≠ calibration")

- **Rigorous**, in its own function with its derivation in the comment: ‖δρ‖_∞ ≤ E_ES(w,σ) · ‖c‖₁, where E_ES is Barnett's 2021 aliasing bound for
  the ES kernel.  Then the consumed measure:
  - single-point: |ΔE_xc| ≤ ‖δρ‖_∞ · Σ_i w_i |v_xc(r_i)|, which must be ≤ ε per atom;
  - forces/geometry: max|ΔH_ij| ≤ ‖δρ‖_∞ · max_ij Σ_i w_i |f_xc(r_i)| |φ_iφ_j|, which must be ≤ ε per element.
  Both prefactors are integrals the quadrature already answers, so no new face is needed.  From them, invert for the smallest w.
- **Calibrated**, a named η = delivered/bound measured on the first call of a run by a direct G-sum on a point subsample.  It's an audit, and the run
  report prints the rigorous w, η, the w used and the delivered error.  It's exposed as an expert knob, beside `vlocRoutingEta` in `GridAccuracy`
  (`GPWTolerances.C`), and combined with the rigorous value in ONE visible place.
- A STATED width or tolerance is honoured and reported against ε.

#### Landing in two commits, so the move batch is attributable

1. **NUFFT at near-machine accuracy** (w sized for a relative 1e-13).  Expect 0 anchor moves at 1e-10 against the fc106a15 baselines.  Gates:
   - a unit test of the NUFFT against the direct sum (random coefficients and points, triclinic cell, error ≤ bound) in `src/Math/tests`;
   - a unit test of `FourierMixCD` batched vs pointwise in `src/ChargeDensity/tests`;
   - the MnO gate bit-close (≤1e-10) to the G-sum.
   The cost even at this width is about 31 664·13³ ≈ 7e7 MACs.
2. **w sized from ε** (Rigorous + η).  This is the real ε-sizing, and its moves at ~ε are re-judged per pin 10 as a batch, separately from the speedup.

#### Expected cost (MnO, 160 Ha ball, per channel per call)

| | today | NUFFT, w = 13 (step 1) | NUFFT, w ≈ 7 (ε=1e-5, step 2) |
|---|---|---|---|
| per-point work | N_G/2 ≈ 15k MAC | w³ = 2197 | 343 |
| FFT, about 74³ oversampled | — | ~10 ms | ~10 ms |
| estimate | ~3 s | ~0.15 s | ~0.03 s |

That's 20–100× on these two buckets.  The 1415 + 689 s on the geometry run would drop to tens of seconds, which is about half of the per-SCF growth
the density ball cost.  Bit-identical across thread counts: the interpolation is a per-point gather, and the FFT is serial or deterministic.  The
gradient comes for free when a GGA needs it (a type-2 NUFFT with iG weights), so `FourierMixCD::Gradient` stops being a throw then, not now.

#### Files

- NEW `src/Math/NUFFT.C` (`qchem.NUFFT`, beside `qchem.FFT`): `Type2Plan(index box, σ, w)` and `Sample(points as fractional coordinates, coefficients)`.
  Weights are cached per point set.  Also the ES kernel, its φ̂, and the Rigorous bound function.
- NEW `src/Math/tests/NUFFT.C`, added to an existing UT exe (`UTMath`), with no new exe.
- `src/ChargeDensity/Imp/FourierMixCD.C`: the batched `operator()` goes through the plan.  The direct sum stays as the reference (the audit, and
  the pointwise op).
- ~~`FieldMixer.C` hands a plan in~~: not needed.  The plan is cheap to build per call (deconvolution is O(N_G), and the per-point weights
  are 3·w kernel evaluations), so nothing is cached and nothing is threaded through the mixer.
- Dispatch (above): `src/Math` (the `PointSampling` / `FourierSampler` faces), `src/Mesh/Mesh.C` (implements them, holds the strategy),
  `src/Structure/Imp/UnitCell.C` (the uniform branch installs the raster strategy), `src/BasisSet/DeltaFit_IBS.C` (`f.Sample(*mesh)`),
  `src/Math/ScalarFunction.C` (the default `Sample`).
- `src/BasisSet/Gaussian/Lattice/GPWTolerances.C`: the NUFFT width as Rigorous + Calibrated (step 2).
- `src/ChargeDensity/Internal/Imp/DensitySampler_Singles.C`: report keys (w, bound, η, delivered), written cursor-relative.
- Keep-out files: none touched (no pybind/, Pins, or QueueCampaignPlan).

Unless you or the user redirect, I start step 1 on your go.

#### Side finding: the 122 matrix-free calls are NOT a pin-28 defect

They come from the convergence residual, `LoopDriver.C:90`: `ResidualDensities` presents the in/out pair WITHOUT an exact source by
design, and measures Δρ on the Becke quadrature every iteration.  That is 60 iterations × 2 plus the seeds, which makes 122.  XC never sees
these densities, so there is no OpenWork row.  (A 2× lever is to sample in − out once, since the residual is linear, but the NUFFT makes
it moot.)

## ✅ The run record files the fit bases (2026-11-04, lane B, 6a0b913f)

**The row, verbatim:**

| **The run record does not file the FIT BASES** (user via lane G, 2026-11-04, GUI Basis tab) | a molecular DFT run (CH₄ LDA, DZVP) files `basis/{name, nFunctions, perIrrep}` for the orbital basis and `grids/auxiliaryFit` for the fit's mesh, but neither the Hartree ρ fit basis (A1_coul, `BasisSet/Gaussian/PG_Cart/Imp/IrrepBasisSet.C:228`) nor the Vxc fit (δ on the XC mesh by default, or A1_exch) appears anywhere in the record | file them under `basis` as `fits: [{role, name, nFunctions}]` (role "Hartree ρ fit" / "Vxc fit"), written CURSOR-RELATIVE by the code that creates each fit, inside the caller's basis Section (report = run sequence); every facade (atom / molecule / solid: the GPW plane-wave fits too).  The GUI's Basis tab (qchem-GUI 775dfe2) already reads `basis.fits` · `src/BasisSet/` (CreateCDFitBasisSet / CreateVxcFitBasisSet), the facades' basis sections |

**Verdict.**  The row's placement (`basis.fits`) was overruled twice the same day.  First: the Hamiltonian builds the fits, because only it knows whether a run needs one, so the sequence stays and the basis context does not get them.  Then the user found `grids.hamiltonian.fits` odd and ruled the tree: a top-level `hamiltonian`, each mesh filed in construction order just before the term that needs it, and each term filing its own fit basis.

- The three facades' `Scope("grids")` around the Hamiltonian build is `Scope("hamiltonian")`.  A molecule or atom now files no `grids` at all; a solid's `grids` keeps what is built before the Hamiltonian (`ladder`, `densityEcut`, `cutoffFactor`, `vlocEps`, the audit rows, `latticeSums`) and the still root-anchored `boxTasks`.
- `BuildTermFit` (`Hamiltonians.C`) opens `Hartree` / `Vxc_DFT` around the fit build and files `hamiltonian.<term>.fitbasis: [{role, name, nFunctions}]`.  A fit's own mesh or grid lands in the term (`hamiltonian.Hartree.auxiliaryFit`, `.densityFit`, `hamiltonian.Vxc_DFT.vxcFitGrid`), which also ended two auxiliary fits sharing one key under Vxc = Auxiliary.
- `UnitCell::EmitUniformGridReport` is cursor-relative, so `densityFit` / `vxcFitGrid` move with their context.  That exposed two fit bases that are not the Hamiltonian's: the plane-wave seed's and the G-space mixer's, now `seedDensity.densityFit` and `densityMixer.vxcFitGrid`.  On a Becke run the mixer's was the only `grids.vxcFitGrid`, which read as the Hamiltonian's.  The HF SAD bootstrap's sibling Hamiltonian files under `seedHamiltonian`, a path no check has exercised at run time.
- `Name()` on a Gaussian auxiliary fit is its file (`A1_coul`, `A1_exch`); `BasisSetID` is untouched.

Tests: `FacadeMatrix` (every cell: each `hamiltonian` member is a mesh with points or a term with complete `fitbasis` rows), `CalculationFacade.EachAuxiliaryFitFilesInItsOwnTerm`.  Sweep: 1279/0 short, 11/0 long, 0 anchor moves against the branch's previous sweep (`GPW_NaF.Γ_GeometryRunSizesItsBeckeMesh_Long` flickers 33/34 iterations at dE = 0 between main sweeps, independent of this change).  `qviz/gridinfo.py` moves in lockstep on its own branch (lane G).

## ✅ HF needs no mesh (2026-11-04, lane B, 73e92020)

**The row, verbatim:**

| **HF needs no mesh: its Δρ gate should come from D** (user 2026-11-04: "HF should have no mesh!"; found by lane G in the GUI's Grids tab) | D-MINDRHO (2026-10-22) made Δρ the L1 distance per electron on a real-space quadrature; a Hamiltonian with none (HF, 1-electron) falls back to a default-MeshParams mesh built only for the gate (`SCFIterator.C:288-320`; CH₄ 30 radial × L5 = 1,692 points, no ε link) | USER RULING: when the Hamiltonian has no real-space quadrature, the gate is the S-METRIC D norm, √Tr(ΔD S ΔD S) per electron, and NO mesh is built (`residualMesh` goes, with the D-MINDRHO fallback).  Basis-invariant: it equals the Frobenius ΔD norm in an orthonormal basis, the invariant form of the RMS(ΔD) other codes print.  The Δρ columns need NOT be comparable between HF and DFT, so label the column for what it is, and re-judge the HF gate's default value (L2-type, so its numbers differ from L1); HF anchors may stop an iteration earlier or later, re-judged as stopping-tolerance moves.  DFT unchanged · `src/SCFIterator/Imp/SCFIterator.C` (ResidualPerElectron), `src/SCFIterator/SCFParams.C` |

**Verdict.**  Built to a second user ruling the same day (V-INTERFACE-TIERS: `DeltaRho` + `DeltaRhoLabel` replace `RealSpaceL1Distance`).  The first design kept a capability query on the hub; the user rejected it as Ask, not Tell.  `tHamiltonian::DeltaRho(in, out)` answers the per-electron residual and `DeltaRhoLabel()` titles it; the NaN signal and `SCFIterator::ResidualPerElectron` are gone, and the iterator neither branches nor sees a mesh.

- A Hamiltonian holding a real-space quadrature measures the L1 distance on it, titled `Δρ`, unchanged.  One without (HF, the 1-electron and Dirac test models, Vxc=Auxiliary DFT) measures √Tr(ΔD S ΔD S)/N, titled `‖ΔD‖`, and builds no mesh.  The Vxc=Auxiliary case includes the pseudo-atom SCFs behind the +U projectors, so a solid's log can show both titles.
- `tDensityMatrixDistance` (`qchem.ChargeDensity`) is the D-side face, `tRealSpaceResidual` (`Hamiltonian/Internal`) the quadrature side.  Blocks' squares add.
- **The default stays 1e-4, no new knob.**  Near the gate the D norm reads 0.4-2.3× the L1 on HF H₂O, CH₄, N₂, CO, O₂ triplet, Be, B, Ne, Fe; at the first L1 < 1e-4 it reads 3e-6..1.2e-4.  Record at `SCFParams`.

**Measured against main-e54c99e7** (short 1283/0, long 11/0).  Counted from SCF trace rows, because `anchordiff` printed 0 (it could not see these until `042e06d3`).  Stopping-tolerance moves only: 9 tests take one more iteration, none fewer; over 481 runs compared on their last row the largest energy move is 1e-12 Ha, one print digit.

| test (instances) | iterations | energy, Ha |
|---|---|---|
| `A_HF_U` He (BSpline Low, Medium; BSpliner Medium) | 12 → 13 | -2.8616774411 .. -2.8616799918 |
| `A_HF_U` Be (same three) | 13 → 14 | -14.5729965418 .. -14.5730231252 |
| `A_HF_P` Li (BSpline, BSpliner Medium) | 19 → 20 | -7.4327508213, -7.4327509009 |
| `MoleculeSession.ACarriedStep…` water HF, 3 runs | 94 → 95 in all | -76.0229320021 |

No solid and no quadrature DFT run moved.  No panel case is HF.  The `_Panel` cases in the short tier did not move; the `_Panel_Long` cases were not run.  Their solid SCFs keep the quadrature L1, and their only changed gate is the +U pseudo-atom's (Mn and O), which did not move in the long tier (`GPW_MnO2…RefusesAnAliasedFirstIterate_Long`).

**The trap it found.**  `dynamic_cast` from a virtual base to a public virtual face returns null, silently, on libstdc++ when the face is declared BEFORE a non-public non-virtual base sharing that virtual base (reproduced standalone; libc++ is correct).  It made every DFT run print `‖ΔD‖`, and only the solid golden caught it.  The face now sits after the base, with the reason at `PWTerms.C`.  A text scan of the 318 classes in `src` found no other instance; it is a heuristic, not a compiler check, validated against the pre-fix shape.

## ✅ `anchordiff` counts molecular and atom iterations (2026-11-04, lane B, 042e06d3)

It summed `iters=`, which only solids print, so 73e92020's 9 iteration moves read as 0 moved.  It now also compares the SCF trace rows, keeping both counts because they differ: the rows add seed and pseudo-atom SCFs, and the `PW_*` prototypes print no trace.  The summary line gives how many tests carry an energy and an iteration count (283 of 1277 short tests, was 117).

**Did an earlier "0 moved" hide anything?  No.**  Rerun over the 2026-11-04 baselines, tol 1e-10: the six steps of the main chain from main-6fa69b0d to main-e54c99e7, plus nine lane sweeps against their base, 15 short and 13 long comparisons: no test moved on rows alone.  All seven mains give 7141 trace rows over 281 tests, so the count is live on every log and the iteration counts did not change.  What the old tool did flag is unchanged: `GPW_NaF.Γ_GeometryRunSizesItsBeckeMesh_Long` flickering 33/34, and `GPW_NaF.Γ_Becke_Imp_eqUni` (1.8e-9 Ha), `GPW_NaF.Γ_Imp_eqExactResume_Panel` and `GPW_Si.Γ_Becke_Imp_Kerker_CP2K` (1-2e-10 Ha).

## ✅ The ρ̃ sampler, step 2: the NUFFT width sized from ε (2026-11-04, lane A, a183196b)

**What landed.**  The width is the smallest w with η·`Type2::Bound`(w)·‖c‖₁ ≤ τ, combined in `NUFFTSampler::Load`.  τ is set once from the seed, as the minimum over the consumers the run has: XC energy (target/∫|v_xc|), XC matrix (Cauchy-Schwarz on |f_xc|) and the SCF residual (κ_R·gate·N/2∫1).  A run with no target (eps 0) keeps the fixed step-1 width.

**Rulings (user and controller, 2026-11-04).**  Set everything before the iterations and report problems during them: τ never changes inside an SCF.  The `samp/τ` column is the per-Load audit against the direct series, flagged `!` over 1.  The ONE tightening happens at the next Hamiltonian build, tighten-only and once per run, and is reported.  A stage that keeps its Hamiltonian, or a restart into a new object, is report-only.  κ_R = 0.1 is a margin, not a measurement.

**The calibration.**  The first cut, η = 1/30, came from a synthetic series (delivered/bound 1/80..1/47).  The full sweep then read `samp/τ` peaks of 4.4 (Si Γ Becke Kerker), 2.9 (NaF), 2.2 (MnO) and 1.4 (TiO₂ k222), so routine runs went over target and tightened, and Si ended over.  The column saw-tooths: it falls with ‖c‖₁ at fixed w and jumps when w drops a step, and at that jump the sizing is exactly met, so delivered/bound = η × peak = 0.147 at most.  Times a 1.5 margin, rounded up: η = 1/4.  At 1/4 the sweep peaks at 0.81 with no over-target run, so the tightening is exercised only by `SolidState.AnOverTargetSamplerIsFlaggedAndTightenedOnceAtTheNextSCF`.

**Anchors, pin 10** (merged tree against main-0f69df88; logs `lane-a-step2d_*`).  Short, 2 moves, both Becke + ρ̃: `GPW_NaF.Γ_Becke_Imp_eqUni` −1.8e-9 Ha and `GPW_NaF.Γ_Imp_eqExactResume_Panel` +2e-10 Ha, each about 4.4e-6 Ha by bound (`samp/τ` · τ · P_E).  Long: `GPW_NaF.Γ_GeometryRunSizesItsBeckeMesh_Long` 33 → 35 iterations, ΔE 0.  Uniform, PW and molecular paths: 0 moves.  The panel has no same-commit baseline; against the older 4-thread lane-a-ball5 log (pre step 1) 3 tests moved, and the largest, β-MnO₂ Kerker+Pulay −5.3e-7 Ha, is an eps-0 run stopped at its 40-iteration bound with `samp/τ` at `----`, so the width did not size it.

## The Becke ladder's ~3e-5 is not a floor: it is L=47's angular error, and the ladder ends too low (2026-11-05, lane A)

**Question** (controller): the pilot's Becke climb reads about 3e-5 and never settles on MnO and LiF.  Is that a floor, and from what?  It drives the 23 mesh rebuilds of the MnO geometry run.

**Answer.**  Not a floor, on MnO.  At a fixed converged D, the error is all in the angular degree L, and it falls to about 1.3e-6 by L = 65.  The ladder's top rung is (120, GL 53), whose own error is about 1e-5, so the pilot's rung test cannot pass with margin.  LiF has no split yet: the numbers here are MnO only.

**Tolerance and radial split.**  MnO geometry cost deck, `accuracy` 0, one converged D at (100, GL 47), fixed-D `xcMeshAudit`, max|ΔH_ij| of the XC term:

| reference | beckeEps 1e-6 | beckeEps 1e-10 |
|---|---|---|
| (160, 47) radial only | 5.11e-7 | 5.08e-7 |
| (100, 65) angular only | 2.391e-5 | 2.391e-5 |
| (160, 65) both | 2.391e-5 | 2.391e-5 |
| (120, 53) | 2.457e-5 | 2.457e-5 |
| (140, 59) | 3.052e-5 | 3.052e-5 |

So `beckeEps` is ruled out for this quantity (unchanged to 4 digits), and nR = 100 is converged to 5e-7.  That also caps anything non-angular, sampler noise included, at about 5e-7.

**The angular series.**  nR fixed at 100; `Tr(DΔH)` summed over both spin channels; `sNormDH` as defined at `GridAuditTerm`.  The first L = 65 row measures L = 53's error, not L = 65's:

| base | ref | max\|ΔH\| | Tr(DΔH) | sNormDH | max element |
|---|---|---|---|---|---|
| 47 | 53 | 2.46e-5 | -4.7e-6 | 1.05e-4 | Mn d α = 7.92, diagonal, \|D_ij\| 2.1e-5 |
| 47 | 59 | 3.05e-5 | -3.2e-6 | 1.11e-4 | same, on its AFM partner |
| 47 | 65 | 2.39e-5 | -8.2e-7 | 1.05e-4 | same |
| 47 | 77 | 2.33e-5 | -3.2e-7 | 1.12e-4 | same |
| 65 | 53 | 9.98e-6 | -3.9e-6 | 9.4e-5 | Mn d α = 3.72, diagonal, \|D_ij\| 4.0e-4 |
| 65 | 77 | 1.29e-6 | +5.0e-7 | 1.33e-5 | same |
| 65 | 89 | 1.63e-6 | +4.4e-7 | 1.26e-5 | same |
| 65 | 101 | 1.29e-6 | +4.1e-7 | 1.13e-5 | same, on its AFM partner |

Reading:
- Error against the finest references: about 2.3e-5 at L = 47, 1e-5 at L = 53, 1.3e-6 at L = 65; roughly ×0.4 per ΔL = 6, from those three points.  L = 59 is at least 7e-6, inferred from one element's differences with assumed signs.
- Any floor is at or below about 1.3e-6, under the 1e-5 target.
- At base 47 all four audits pick one element up to the AFM mirror, so 2.46 / 3.05 / 2.39 is not a max over different elements.  It is real scatter of about 0.7e-5 on top of a much larger L = 47 error.  An early reading of "slow, irregular convergence" was wrong.
- The D-weighted trace is 5 to 70× below max|ΔH| at L = 47.  The geometry gate stays max|ΔH| ≤ ε per element; `trDdH` and `sNormDH` are report-only (pins 19/20).

**The partition matters.**  beckeOrder 4 at (100, 47) against L = 65 / 77 reads 3.71e-5 (up from 2.39e-5 at order 3), `sNormDH` 3.7e-4 (up from 1.05e-4), and the largest element moves to a diffuse Mn d (α = 0.38, majority spin, |D_ij| 0.13) that the density does see.  ΔE against L = 65 grows from 6e-7 to 6.9e-6 Ha.

**Cause: not established.**  At order 3 every max element is a diagonal Mn d in that Mn's minority spin channel (8 of 8 audits), which suggested v_xc^σ ~ ρ_σ^(1/3) being sharp in angle where ρ_σ is small.  The order 4 result shows a second, partition-driven source, so that is at most part of the story.  A first guess, tight functions on neighbouring atoms, was wrong: both functions sit on one atom.

Runs: `~/Code/qchem6-runs/MnO/lanea_ladderfloor_*`.  Measures: d0e8bcb0.

## The Becke mesh's error, axis by axis: rule alignment, radial floor, partition order (2026-11-05, lane A)

Continues the section above.  MnO AFM-II and LiF, ε = 1e-5 per element, one converged D per table, runs `~/Code/qchem6-runs/{MnO,LiF}/lanea_*`.  Plain meshes (`xcMeshSymmetric:false`), no fold.  Lebedev above degree 35 ran on a measurement patch (`MnO/lanea_lebedev_measure.patch`), not on main.

**Aligned Lebedev carries an extra error that rotation removes.**  MnO's mesh was Gauss-Legendre above 35 because `LebedevMenu` stops at 434 directions, not for a symmetry reason.

| nR 100, order 3 | points | angRot 0 | 0.5 | 1.1 |
|---|---|---|---|---|
| Lebedev 35 | 115k | 4.0e-4 | 8.6e-5 | |
| Lebedev 41 | 157k | | 1.7e-5 | |
| Lebedev 47 | 204k | 6.3e-5 | 8.1e-6 | 1.3e-5 |
| Lebedev 53 | 259k | | 6.0e-6 | 5.8e-6 |
| Lebedev 65 | 386k | 4.7e-6 | 2.4e-6 | 1.9e-6 |
| GL 47 | 303k | 2.4e-5 | 1.9e-5 | |
| LiF Lebedev 47 | 100k | 1.0e-4 | 7.0e-6 | |
| LiF Lebedev 53 | 127k | | 4.6e-6 | |

Two sources: an alignment alias (largest element a diffuse s) and a tight-d angular error (Mn d, α = 7.92 or 3.72, |D_ij| small) that remains and sets the rung.  A 0.05 bohr displaced Mn at Lebedev 53 reads 5.96e-6, unchanged.

**A cubic-harmonic measure ranks the energy error, not max|ΔH|.**  M(θ), the rotated rule's quadrature of unit A1g harmonics above its degree, weighted e^{-0.1(l-L)}, needs no SCF (`MnO/cubicM.py`).  Lebedev 47:

| θ | √M | \|ΔE\| | sNormDH | max\|ΔH\| |
|---|---|---|---|---|
| 0 | 0.097 | 3.6e-5 | 7.9e-4 | 6.3e-5 |
| 1.1 | 0.0245 | 8.0e-6 | 1.6e-4 | 1.3e-5 |
| 0.5 | 0.0167 | 1.4e-6 | 7.8e-5 | 8.1e-6 |
| 0.32 | 0.0113 | 1.2e-6 | 5.2e-5 | 9.5e-6 |
| 0.12 | 0.0097 | 2.7e-7 | 6.3e-5 | 1.1e-5 |

θ = 0.32 and 0.12 were predicted before they were run.  The trace with D is close to site-invariant and a single tight-d element is not, which fits.  A free rotation axis reached the same best M as the fixed axis.

**Radial floor at L = 65** (|nR 100 − reference|):

| ref nR | MnO max | MnO sNormDH | LiF max | LiF sNormDH |
|---|---|---|---|---|
| 40 | 1.1e-4 | 4.3e-3 | 3.9e-4 | 2.9e-3 |
| 60 | 2.8e-5 | 3.8e-4 | 1.1e-5 | 1.8e-4 |
| 80 | 1.9e-6 | 3.9e-5 | 9.3e-6 | 1.9e-5 |
| 120 | 2.3e-7 | 3.6e-6 | 2.2e-7 | 1.8e-6 |

`sNormDH` falls about ×10 per step on both; LiF's max jumps with the worst pair.  In L the two measures are about equally smooth, and `sNormDH` is the more angle-sensitive.

**Partition order.**  `BeckeCutoff` iterates k+1 times, so the code's order 2 is Becke's published three-iteration scheme.  Order 2 at (100, 47) reads 2.42e-5 / 2.36e-5, no gain over order 3, and equivalent sites differ by 6e-7 relative at 40/29.  Order 4's excess is angular, not a floor: GL 65 reads 2.2e-6 on MnO, but rotated Lebedev reads 1.2e-5 (53) and 2.4e-6 (65) on MnO, 2.9e-5 and 1.35e-5 on LiF.

Mesh build, one build each, counts independent of load:

| mesh | order | images/pt | P-set/pt | `BeckeCutoff` calls | build |
|---|---|---|---|---|---|
| 100/65 | 3 | 1020 | 28.9 | 2.97e10 | 146, 148 s |
| 100/65 | 4 | 309 | 4.8 | 1.31e9 | 8.3, 6.0 s |
| 40/29 | 2 | 17270 | 615 | 1.13e11 | 458 s |
| 40/29 | 3 | 1051 | 30.6 | 1.67e9 | 6.9 s |
| 40/29 | 4 | 307 | 4.9 | 7.15e7 | 0.36 s |

The cause is the shell loop's `newPB<eps` stop and `inPSet` test (`src/Structure/Imp/UnitCell.C:243-249`).  No weight is lost: at order 4, `beckeEps` 1e-8 leaves Sum(w) unchanged to 8 digits, and Sum(w) reaches V/4 = 74.088 at 120/89 (74.08803; order 3 74.08800), so its 100/65 deficit (74.08758) is quadrature.

**Symmetric runs.**  `angRot` never reaches `MakeInvariantAngularMesh`.  Imposed MnO at degree 47: 2342 directions, 622416 points, 26554 after a 23.4× fold, 3.8e-6 against its own degree 53.  The rule exists at 53 and fails at 59 for this 12-op site group.

### The ladder's gates (2026-11-05, lane A)

Both gates are full geometry runs from Auto at ε = 1e-5, 4 threads, one displaced atom: `~/Code/qchem6-runs/{MnO,LiF}/lanea_ladder_geometry.{json,out}`.

| | MnO | LiF |
|---|---|---|
| settles on | 100 / Lebedev 53 rotated | the same |
| pilot rounds, re-convergence iterations | 3, 6 | 1, 0 |
| angular, L 53 against L 65 | 5.3e-6 | 3.4e-6 |
| radial, nR 100 against 120 | 3.1e-7 | 1.3e-7 |
| Becke builds in the run | 22 | 9 |
| wall | 4354 s | 1381 s |

The old walk on MnO took 3327 s and 23 builds and never settled: it ended NOT SETTLED on 100 / GL 47.  **The ladder's estimate was wrong: it said about 15 builds and about 2550 s.**  The run is 31% dearer in wall and 34% in CPU than the old walk.  The estimate counted one build per rung, but a pilot round builds the run's mesh, both references and the ball's reference, and a re-convergence builds the mesh again.  Most of the old walk's builds were also far cheaper than the entry rung, at 32k to 185k points against 259k, so mesh build per call went from 31 s to 62 s.  The rounds are now set by the density ball (factor 2, 4, 8), not by XC.  OpenWork carries the row on the repeated builds.

The displaced point moves by 1.8e-6 Ha between the two meshes (E −61.3031147 against −61.3031165), and dE by 4.7e-7.

**Lebedev 59 was measured and left out.**  MnO 5.4e-6 / 4.1e-6 at the two angles beats L 53 only slightly, and LiF at 8.7e-6 is worse than L 53's 4.6e-6, so it buys nothing over 53.  The ladder uses 53, 65, 77 and 89 (GL above the Lebedev menu), and the audit reference for L 53 is L 65.

The first LiF attempt stopped in pre-flight, because its deck stated rasters beside an accuracy target.  With those removed it ran as above.

### The pilot's two references (2026-11-05, lane A)

The row "rebuilds an identical Becke mesh 8 times" named the repeated build.  The loop was the defect, so no cache was added.  A round built four references: raster, angular rung, radial rung, ball.  The raster and ball references built the run's own mesh and compared it with itself, and the two rung references held rasters and ball fixed.  On a Becke mesh only the delta fit is admitted and it samples rho from D, so XC depends on the mesh alone (every Exc row of those two audits read exactly 0 in all 3 MnO rounds), and every other term depends on the rasters and ball alone.  So the raster reference now carries the angular rung and the ball reference the radial one.  Every audit still runs at the final D.  An admitted PlaneWave fit on Becke would break the separation; the comment at the site says so and Si's gate checks the raster half.  The ball half has no test of its own: it rests on the Exc rows above and on the gates' identical output.

| gate | mesh lines in the run | wall | pilot output |
|---|---|---|---|
| MnO geometry | 22 -> 16 (stage 15 -> 9) | 4354 -> 2665 s | identical, run id aside |
| LiF geometry | 9 -> 7 | not comparable (ran beside the panel) | identical |

**The 1689 s is mostly load, not this change.**  The two MnO runs saw different boxes: the identical work before the first geometry audit took 1284 s then, 879 s now, and collocate 8.96 against 6.67 s a call.  Scaling by that 1.35-1.45 puts the change's own saving at about 450-620 s, under the 900 s predicted beforehand.  No A/B on one load was run.  Sweep against main-9ed28468: short 1299/0, long 11/0, 0 moved; `SolidSession.GeometryRunPicksAndSizesTheBeckeMesh_Long` 2113 -> 1573 s, `GPW_NaF.Γ_GeometryRunSizesItsBeckeMesh_Long` 1978 -> 1495 s, same caveat.  Peak RSS unchanged at 3.1 GB, though each reference now holds a raised raster and a finer mesh together.

**Molecular PP terms in closed form ✅ (2026-11-05, lane S, `lane-s-pp-analytic`).**  `PP_Local` and `PP_NonLocal` now ask the block's `Orbital_PP_IBS` (Gaussian short part, erf long part as a 3-centre Coulomb integral against the Gaussian core charge, Gaussian-Y_lm projectors), in Cartesian and spherical, a SALC block forwarding to its raw basis; a Gaussian molecular PP run builds no PP mesh.  The mesh stays for atomic radial blocks.  TeF₆ GTH q6/q7 against the converged -153.4084029984: 80_29 from -6.11e-4 to -2.21e-5, which is the audited XC error (-2.19e-5), and 150_47_gl from -5.7e-6 to -4.0e-7; wall time and RSS unchanged.  Pin 10: `M_DFT.F2_PP` moved +4.2e-5, from 5.7e-5 to 1.5e-5 off PySCF, and `F2_PP_AuxiliaryFit` the same.  `Si2_PP_U.LargeSeparation` -9.6e-2 is a two-iteration bounded iterate of a pair that never converges: its Si atom agrees to 1.4e-11 and the pair's first iterate to 1.1e-5, and iteration 2 amplifies that.  `OSi_PP_U` runs 26 iterations, not 16, and lands 1.4e-5 lower, from a first iterate 2.3e-4 apart.  It did NOT close the spherical NaN row: that run still dies with no PP mesh.

Original row: | **Molecular PP terms are quadrature, solids' are analytic: a structure gap** (user 2026-11-05) | GPW assembles the local PP from the analytic form factor and the KB projections analytically (`GPW_Evaluator`, `ProjectOntoGaussianYlm`); the molecule's `PP_Local` / `PP_NonLocal` integrate both on a Becke mesh (`MolecularPseudopotentialPlan.md` §3, "Path A", which claimed the mesh "nails it").  TeF₆ on GTH-q6: the default mesh is −6.1e-4 Ha off its best value, worse than all-electron; lane B's fixed-density audit will say how much is these terms | if the audit confirms it: take the molecular terms to the same analytic forms (Gaussian short part, erf long part, Gaussian-Y_lm projectors), which removes a grid from ε's inventory.  They are static, so this is accuracy, not run time · `src/Hamiltonian/Internal/Imp/Hamiltonians.C`, `GPW_Evaluator.C` |


## 40/29 on an energy run, and the box's 9× (lane A, 2026-11-06)

**Which runs reach the aligned rule.**  `ResolveImposition` (`Imposition.C`) decides.  An imposed run, and a free run with `xcMeshSymmetric` true (the default), take the site-adapted mesh, never the aligned rule.  A free run with `xcMeshSymmetric` false takes the plain mesh, aligned unless `angRot` is stated: every `SolidSession`, every deck `displace` path (both require it), `IntegrationTests/GPW/LiF.C:27`, `NaF.C:300,345`, `src/Calculation/tests/SolidSession.C`, and the 7 `decks/probes/*_molopt.json`.  The ladder rotates only above degree 35, so its rungs 29 and 35 are aligned.  Si-like systems under Auto, and the Becke veto, run the uniform raster and have no angular rule.  Molecules: plain, aligned, order 2, `mhl_m` 3.

**40/29 at each run's own converged D**, against a GL reference ((100,53) on MnO and NaF, (150,53) on TeF₆), ε = 1e-5:

| run | points | max\|ΔH\| | ΔE | ΔE / ε·N | wall |
|---|---|---|---|---|---|
| MnO plain | 31,664 | 1.062e-3 | −4.857e-4 | 12 | 8:57 |
| MnO rotated 0.5 | 31,840 | 1.118e-4 | +1.618e-4 | 4.0 | 8:31 |
| MnO site-adapted | 97,256 | 2.847e-4 | +3.964e-5 | 1.0 | 23:41 |
| NaF plain | 15,928 | 1.120e-4 | −6.829e-5 | 3.4 | 3:21 |
| NaF rotated 0.5 | 16,016 | 5.981e-5 | −4.482e-5 | 2.2 | 3:14 |
| NaF site-adapted | 41,440 | 6.401e-5 | −5.079e-5 | 2.5 | 5:42 |
| TeF₆ plain | 63,620 | 1.435e-4 | −1.333e-4 | 1.9 | 15:02 |
| TeF₆ rotated 0.5 | 64,152 | 5.404e-5 | −3.910e-5 | 0.56 | 13:35 |
| TeF₆ site-adapted | 202,184 | 7.023e-5 | −5.934e-5 | 0.85 | 34:09 |

ε·N is 4e-5, 2e-5, 7e-5.  Walls ran beside other jobs, so read them as rough.  Rotation is free and beats the plain rule on every row.  Site-adapted costs 2-3× the wall and does not beat rotation on max\|ΔH\| anywhere; its one better ΔE, on MnO, is a single sign-changing number.  Decks: `~/Code/qchem6-runs/{MnO,NaF,TeF6}/energy40/`.

**The 9× is the Becke order.**  The molecule ran order 2 and `mhl_m` 3 (`MoleculeSetup.C`), the box order 3 and `mhl_m` 2 (`BeckeXCParams`).  One change at a time on the box at (80,29), which is the row's question:

| box at (80,29) | points | vs (200,35) max\|ΔH\|, ΔE | vs (150,53 GL) max\|ΔH\|, ΔE |
|---|---|---|---|
| as run: order 3, m 2 | 127,980 | 6.704e-5, −6.744e-5 | 8.754e-5, −8.791e-5 |
| `beckeEps` 1e-10 | 129,432 | 6.704e-5, −6.731e-5 | 8.754e-5, −8.788e-5 |
| Cartesian | 127,980 | 6.658e-5, −6.550e-5 | 8.694e-5, −8.538e-5 |
| `mhl_m` 3 | 115,210 | 6.146e-5, −6.357e-5 | 8.196e-5, −8.417e-5 |
| rotated 0.5 | 128,908 | 1.534e-5, +5.270e-6 | 1.249e-5, +1.033e-6 |
| **order 2** | 134,504 | **7.647e-6**, −1.471e-6 | 9.201e-6, −3.497e-6 |
| order 2, m 3 (the molecule's mesh) | 121,038 | 7.436e-6, −2.052e-5 | 8.991e-6, −2.249e-5 |
| order 2, m 3, rotated | 120,812 | 4.741e-6, −1.945e-5 | 4.763e-6, −1.923e-5 |

The molecule read 7.36e-6, ΔE −2.19e-5, as an XC-only fixed-density audit (the PP matrices never entered, so the analytic PP landing does not change it).  The box on the molecule's mesh reads 7.4e-6, ΔE −2.05e-5: no box effect is left.  The dropped tails were never the cause: at `beckeEps` 1e-10 the box still drops 39,688, because those points carry no weight.  The error needs both the aligned rule and order 3, since rotation or order 2 alone removes most of it; why they combine is not measured.

**Order 2 is not a general cure.**  NaF plain 40/29: 7.63e-5, −5.46e-5 (order 3 1.12e-4, −6.83e-5); rotated 5.79e-5, −3.80e-5.  MnO plain: 8.11e-4 but ΔE −7.13e-4, worse than order 3's −4.86e-4 (first audit only; stopped, as the second audit was not needed).  And the periodic partition pays for it: TeF₆ (80,29) 1,264 s to build against about 58 s, NaF 40/29 1:32:39 of wall against 3:21, because the softer cutoff widens the image series the screen must sum.  `beckeOrder` stays 3.

**The audit's reference carried the same error.**  Every mesh above that is converged gives one number: run minus ΔE against (150,53 GL) is −153.403944 on order 3 plain, order 3 rotated, and order 2.  The aligned (200,35) sits 2.05e-5 off it, the rotated one 4e-6.  `xcRef` copies the run's `angRot`, so a plain run's Lebedev reference is aligned.  This confirms the audit's sign too, on all three materials.

**Box against molecule, closed.**  Box at converged XC −153.4039436; add back V_long −5.2e-5 and Hartree −2.8e-6 from the pilot's audit: −153.4038888.  The exact-J spherical molecule at 200_53_gl reads −153.4038870.  They agree to 1.8e-6, where `TeF6EpsilonReport.md` §6 had 2.5e-5 from the aligned reference.  That molecule number predates the analytic PP and was not rerun.

Original rows, verbatim:

| **Does the energy default carry the aligned-Lebedev error?** (lane A's ladder note, 2026-11-05) | aligned Lebedev is 5× worse than rotated at L=35 on MnO, ΔE 2.4e-4 against ε·N = 4e-5, and 40/29 and 60/35 are aligned Lebedev.  But a free run with the default `xcMeshSymmetric` true uses the site-adapted rule, and rotation only reaches a plain mesh, so the affected runs are not yet known | first list which runs reach the aligned rule, and measure 40/29 itself.  Only then decide on rotating it, perhaps with an M-chosen angle (M ranks \|ΔE\|); that moves anchors, so its own batch under pin 10 · `BeckeXCParams`, `AngularFactory.C` |

| **A box's XC error at a nominal Becke mesh is 9× the molecule's at the same mesh, and the raster pilot does not size that mesh** (lane B, 2026-11-05, from TeF₆ in a 20 bohr cube, GTH-LDA, DZVP-MOLOPT-SR-GTH, ε 1e-5) | The pilot settles every raster term inside ε·N and leaves `XC L0` at 0, so `xcMesh` stays at its default (40,29): against (200,35) that is ΔE −1.13e-4 (1.6× ε·N) and max\|ΔH\| 1.23e-4 (12× ε), a miss only a stated `xcMeshAudit` shows.  At a stated (80,29) the box reads max\|ΔH\| 6.7e-5, worst on Te d; the molecule at (80,29) reads 7.4e-6 against 400_77_gl.  The two are not the same mesh or basis: the box drops 24% of its points as tails (`beckeEps` 1e-6, 127,980 against 168,688) and is spherical, the molecule Cartesian.  Cause unmeasured | rerun the box at (80,29) with `beckeEps` 0 (or far below 1e-6) and Cartesian, one change at a time; then decide whether the pilot should size the XC mesh.  Record: `~/Code/qchem6-runs/TeF6/TeF6_report.md` §6 |

**Multipole's ρ̃ leg and its `lMax` guard ✅ (2026-11-06, lane S, `28eae5f9b`).**  Cause: ρ̃ synthesized onto the mesh lost charge, 30× the mesh's error, unmoved by L and falling 10× with nR.  ρ̃_A is now integrated on its own complete shells with unpartitioned weights (`SiteMultipoles::OwnShellOverlap`), where Y_lm is exactly orthonormal, so ρ̃'s charge is the mesh's and the diagnostic is right by construction.  Water 80/29: +530 → +32 μHa, δ's level; TeF₆ PP 80/29 +540 μHa (δ +542).  The bound is rigorous: 2·lMax ≤ the angular degree, `MultipoleLMaxBound`, refused by `AdmitCoulomb`.  The two loosened oracle tolerances are tightened to 5e-5.  The rows' claim "flat in l" was premature: see the open row on multipole Vxc.

Original rows:

| **Multipole's coarse-rung error is located but not explained, and its charge diagnostic cannot see it** (user question 2026-11-06) | same potential and SCF as δ × Poisson: on TeF₆ PP 80/29 Ekin, Een and Exc agree to every printed digit, and only Eee differs, 347.15030 against 347.14392 (6.4 mHa, before the radial-map fix `b0f1ea7e`), so it is all the energy's ρ̃ leg.  Flat in l (water 8-20 to 0.03 mHa), so not l truncation.  `HartreeChargeError` reads the δ mesh charge, 3.347e-5 on both routes, so a multipole run reports ratio 6.9 while its real error is about 13× that | **Multipole stays a real option, and its Eee must be right (user 2026-11-06: δ being better is not known to be universal).**  Lane S, after its TeF₆ table: measure ρ̃ − ρ at the mesh points (its integral, where it sits, by shell and by site) at fixed D, then nR alone and L alone; fix the cause; report ρ̃'s own charge on the multipole route · `src/BasisSet/MultipoleFit_IBS.C`, `src/Mesh/SiteMultipoles.C` |

| **Multipole `lMax` has no guard against the angular degree** (lane S 2026-11-06) | on `CoulombFit::Multipole`, `lMax` too high for the mesh aliases: water at 80/29 is flat to 0.03 mHa for l = 8-20, then 5 mHa out at l = 24, silently.  The route's own ρ̃ leg needs the shells to resolve it, so a user can reach this | derive the bound from the rule's degree, keep any measured margin a separate named calibration, refuse above it with an Outcome in `AdmitCoulomb` · `src/Hamiltonian/Imp/Factory.C` |

**The aligned rule leaves the default ✅ (2026-11-06, lane A).**  `BeckeXCParams` now rotates a Lebedev rule by 0.5, every ladder rung rotates, and `XCAuditReference` never gives an aligned reference (a run stating `angRot` 0 keeps its rule and is still audited against a rotated one).  It reaches only the plain mesh; site-adapted and imposed runs are bit-identical.  (40,29) stays the energy ballpark with no energy-stage climb, under the "H before E" ruling.

What rotation gives up, measured on the rule against its own C4z image at the converged D (temporary probe, not committed; the aligned control reads 6.7e-16): NaF 40/29 max|ΔH| 9.6e-6, ΔE 3e-12; NaF 60/35 5.5e-6, 3e-10; TeF₆ box 40/29 1.65e-5, 2e-9.  In a NaF session with F at ±0.1 bohr on each axis, the rotated rule splits the axes by 2.0e-7 Ha on a curvature of 2.7e-5, and the ± mirror holds to 2e-9 (aligned 9e-10).  So exact cubic equality is lost at about ε on H at this mesh, under the rule's own error (5.4e-5 TeF₆, 6e-5 NaF).

Anchors, by anchordiff against main-0b0626aa at tol 1e-10, final sweep: 0 of 1317 short tests and 1 of 14 long ones moved (a first sweep moved one short test, below).  `GPW_LiF.Γ_Becke_XCMeshAuditReadsTheCurvatureError_Long` re-judged under pin 10 against GL (80,48): curvature minus the audited mesh error is 1.22e-5 aligned, 1.245e-5 rotated, and Etot minus ΔE agrees to 1e-6; max|ΔH| falls 2.1e-3 to 2.0e-4 and the mesh error on the curvature 1.97e-5 to 6.8e-6.  The gate now pins that difference too.  `SolidSession.SessionSizesItsRastersOnMatrixElements` states a deliberately coarse Lebedev 5, which inherited the rotation (Etot +8e-2 at 12 directions): it now states `angRot` 0 and matches the baseline.  A deck stating a low Lebedev degree inherits the rotation the same way; below 29 the rotated rule is unmeasured.

The health audit is dropped.  One rung up, rotated, on the TeF₆ box is about 60 s, 48 s of it the mesh build: a third of a 3:00 stated-raster run, about 4% of the 23:41 pilot run.  That misses the 5% bar on the run it would cost most, and a stated `xcMeshAudit` already does it.  It read well: (60,35) gave 5.37e-5, −3.66e-5 where (150,53 GL) gave 5.40e-5, −3.91e-5.

Original row: | **The energy default still carries the aligned-Lebedev error on the plain route, and `xcMeshAudit`'s Lebedev references carry it too** (lane A, 2026-11-06) | 40/29 measured on MnO, NaF and the TeF₆ box.  On the plain route rotation costs nothing and cuts max\|ΔH\| 2-9× and \|ΔE\| 1.5-3.4×; site-adapted and imposed runs do not reach the aligned rule.  An audit reference at degree 29-35 copies the run's alignment, so (200,35) was 2e-5 off on TeF₆.  No route meets max\|ΔH\| ≤ ε at 40/29 (5-106×), which under "H before E" is not a defect of the energy default | proposal with the controller: rotate the plain rungs below degree 36, rotate audit references, report one health audit, no energy-stage climb · `XCPolicy.C`, `Deck.C` `xcRef`; record `doc/Records/OpenWork_History5.md`, "40/29 on an energy run" |

## CP2K MnO at matched accuracy: no cutoff meets the H target (lane A, 2026-11-06)

MnO AFM-II, the 118-function VA span, free, kT 5e-3, both codes serial and both stopping at max|ΔP| 1e-8 (rule 3f).  CP2K `mpirun -np 1`, `~/Code/cp2k-runs/MnO/lanea_time_c{576,2304}`, REL_CUTOFF 60, Broyden.  qchem `~/Code/qchem6-runs/MnO/match/`, Kerker/Pulay, `QCHEM_OPENMP_THREADS=1`, 99% CPU.

**CP2K has no matched cutoff for H.**  At one fixed density, max|ΔH_ij| against 3200 Ry reads 1.08e-3 at 1152 Ry, 2.66e-4 at 1600, 4.34e-4 at 1800, 3.71e-4 at 2000, 2.19e-4 at 2304: no decay past 1600, a floor 20 to 40 times ε.  REL_CUTOFF changes nothing.  The worst element is the diagonal of index 89/106, O's last s by my reading, not verified.  Off the diagonal the maximum is 1 to 2e-4.  This is the same floor our uniform XC raster showed on Si, and CP2K's XC is on a uniform raster, but the cause in CP2K is unmeasured.  So the rows below are the energy match and the best H CP2K reaches.

| run | max\|ΔH\| | Etot (Ha) | steps | wall | setup | s/step | peak RSS | load |
|---|---|---|---|---|---|---|---|---|
| CP2K 576 Ry, the energy match | 1.4e-2 ᵃ | −61.30305970 | 116 | 1587 s | 9 s | 13.6 | 0.64 GB | 5–7 |
| CP2K 2304 Ry, the best H | 2.2e-4 | −61.30305343 | 119 | 2749 s | 12 s | 23.0 | 3.9 GB | 2–5 |
| qchem, grids stated | 5.3e-6 | −61.30305473 | 18 | 1523 s | 246 s | 70.9 | 2.8 GB | 1 |
| qchem, the pilot sizes the grids from nothing | 5.3e-6 | −61.30305473 | 18+10+4 ᵇ | 7508 s | | | 3.5 GB | 7 → 1 |

ᵃ Against 2304 Ry, r240, at fixed density.  ᵇ The H pilot's own SCFs; its 7508 s includes the 777 s energy run the deck does first, and 3 audit rounds.  The rest of that run is one displaced point, 1328 s, 17 iterations.

CP2K's energy ladder 144/288/576 Ry gives −61.30332518, −61.30309458, −61.30305970, limit about −61.303054: 576 Ry is 5.7 μHa off it, and 2304 Ry is 1.3 μHa from us.  The qchem H figure is the audit's XC angular term against one rung up; Hartree is 1.1e-7, the density ball 3.4e-7.  The two H columns are measured differently: CP2K's against its own higher cutoff, ours against the pilot's reference.

**Reading.**  With grids known, a qchem H-sized run costs what CP2K's energy match costs, and 0.55 of its best-H run, at an H error 40 times below that floor.  The win is entirely the step count.  Per step CP2K is 3 to 5 times cheaper, and 116 Broyden steps against 18 Pulay steps is a mixer difference I did not try to tune away, so a better CP2K mixer would move this.  A single run that pays for its own sizing costs 2.7 times CP2K's 2304 Ry wall; the pilot pays once for a series that shares a grid set.  qchem's own energy default is not in the table because it does not match: 40/29 here is −61.3028931, 1.6e-4 off, 4 times ε·N, as already measured above (777 s).

Walls are not quiet-box numbers for CP2K.  qchem's 61 s/iteration held from load 7 to 1 in the pilot, and an earlier CP2K 576 run took 2015 s at unknown load against 1587 s here, so read CP2K's walls to about 25%.

**The stated run first came out 5.8 μHa high.**  `accuracy` 0 drops the derived `vlocEps` 1e-6 back to the no-target routing, −61.30304896.  Stating `vlocEps` 1e-6 gives −61.3030547298, the pilot's value to every printed digit, cold from IonicSAD.  So the settled grids are a complete description only with `vlocEps`.

Original row: | **CP2K MnO timing at MATCHED accuracy** (user, 2026-11-03; step 3 of the sequence above) | the §5a parity rows compare CP2K at 288 Ry, whose own cutoff ladder still moves E by 35 μHa (144/288/576 Ry), against our ε-sized grids; raise CP2K's CUTOFF/REL_CUTOFF until its max|ΔH_ij| (or the energy, for an energy run) meets the same 1e-5, and time that, serial, `mpirun -np 1` · `doc/Benchmark.md` §5a, `~/Code/cp2k-runs/MnO/`

## ✅ The pilot prints what restates its grids (lane A, 2026-11-07, 13d4aff4a)

The raster pilot now prints, and files under `rasterPilot.restate`, the deck keys that reproduce the run it settled: `accuracy` 0, `cutoffFactor`, the three raster demands, `vlocEps`, and the Becke mesh's keys, each double written round-trip exact.  Short sweep green, 0 of 1319 moved by anchordiff against main-bc15bfd5 at 1e-10.

Original row: | **A deck restating a pilot's settled grids also needs its `vlocEps`** (lane A, 2026-11-06) | `accuracy` 0 with the settled rasters and mesh, on MnO VA, falls back to the no-target V_long routing and lands 5.8 μHa off the pilot's energy; stating `vlocEps` 1e-6 reproduces it to every printed digit.  The pilot's SETTLED line names the rasters and the mesh but not `vlocEps`: print it there, and check what lane H's Li grid set restates · `doc/Records/OpenWork_History5.md`, "CP2K MnO at matched accuracy"

## ✅ Grid-pilot beta findings (b)-(d) (lane A, 2026-11-07, b1f150ac8)

(b) the sizing line prints ε in %g (`sizing the rasters to 2e-05 Ha`), and GatherSharpness reads α_max off the ladder's own resolution (`ExponentAtCutoff`), so the line and the ladder agree.  (c) a postSCF `gridAudit` / `eggBox` / `xcMeshAudit` sets ok = (measured ≤ the run's target): max|ΔH_ij| ≤ ε on a geometry run, |ΔE| ≤ ε·N on an energy run; the record carries `target` and `measure`; accuracy 0 is not judged and says so.  Recorded ok flips (true → false), every test and deck carrying an audit checked (no deck in `decks/` carries one): only `GPW_LiF.Γ_Becke_XCMeshAuditReadsTheCurvatureError_Long`, both geometries, |ΔE| 7.39e-5 / 6.72e-5 vs 2e-5 Ha; it asserts the audit's numbers, not ok, and stays green.  `GPW_Si.Γ_GridAuditSeparatesTheRasterTerms`, `Γ_EggBoxLeavesTheAnalyticTermsInvariant`, `Γ_GridAuditReusesThePilotsLastRound` and `Deck.AK222RunFilesEachRasterVerdictOnceAndTheAuditAsksTheBasisForIt` meet their targets.  (d) a postSCF `gridAudit` at the pilot's margin (8) reuses the pilot's last round when that round audited the same D (a D generation counter) and was not a Becke angular rung; the reference reads `(the pilot's last round, not re-measured)` on the console and in the report.  Gate `GPW_Si.Γ_GridAuditReusesThePilotsLastRound` (reuse at margin 8, re-measured at 12, the restated deck re-measures cold and agrees to the SCF gate).  Sweep: short 1319/1319, long 13/13, anchordiff vs main-2bdb1256 0 moved in both tiers.

Original row: | **Grid-pilot beta-test findings on Li (lane H, 2026-11-03)** | (a) FIXED 2026-11-05: the record carries the console's verdict; (b) solids print `sizing the rasters to 0.000 Ha` (fixed-format ε: print in %g / scientific); (c) `xcMeshAudit` reports ok=true at 44× ε (bcc q3 Becke 40/29 → 60/35 moves 0.44 mHa; the radial warning predicted it), and the q1 solids auto-chose a 512-point Uniform XC mesh the pilot calls converged while Becke 60/35 differs by 26-35 μHa -- the audit's ok criterion is not ε; (d) COST: bcc q3 aug SCF 590 s vs the pilot ~17,000 s (round-1 audit 11,150 s), and the postSCF gridAudit repeats the pilot's last audit (2,460 s); (e) the q3 atom L28 was OOM-killed (global) in a postSCF audit at 4.3 GB RSS | lane A, with the sampler work (much of (d) is the ρ sampler, OpenWork §4 sequence): (b) is a small fix; (c) joins the ε row's known miss 1 (energy-run Becke audit) and the Becke-ladder-does-not-settle row -- ok must mean ≤ ε; (d) skip the postSCF gridAudit when the pilot's last round audited the same grids at the same D; (e) an RSS budget for the audit's reference basis (memsafe-style refusal, not a kill) · `src/Calculation/Imp/SolidCalculation.C` (raster pilot, postSCF audits) |

## ✅ The short-range local PP, rebuilt per build and per k (lane A, 2026-11-07, 3bdde5b30 + 7611f7610 + 852e8af8a)

Two levers, as the controller sequenced them.  (1) The analytic short-PP block was cached under the basis's full `BasisSetID`, which
carries the rasters, ladder and tolerances, so every rebuilt basis missed; it now keys on what the analytic sum reads (molecular basis,
k, cell, home-only, screenEps; `GPW_Evaluator::LocalPPShortIDFragment`; the grid fallback keeps the full key).  Si k222 with the pilot:
`local-PP SHORT [x8]` (one per k block) against `local-PP LONG [x24]`.  (2) The MnD evaluator builds each 1E operator's R-resolved
series {(n, K_ij(n))} once per (operator, cell, sites[, local operator]) and phases it per k; one fix for overlap, kinetic, nuclear and
the local Gaussian, with the `LatticeSum1E` face unchanged.  The phase sum repeats the per-k arithmetic term for term, so the routes
agree bit for bit (gate `LatticeSeries.TheRResolvedRouteEqualsThePerKRoute_{Si,MnO}`, all four operators at five k).  Measured in one
binary, interleaved, 1 thread pinned, 8 k: Si local Gaussian 3.06 / 0.376 s, S+T+V 0.374 / 0.052 s; MnO valence 20.6 / 2.58 s,
1.50 / 0.21 s: the cost becomes one build.  Cosmetic: the reused postSCF audit names "pilot reference (x + 8) = x(...)" (852e8af8a).
Sweep lane-a-shortpp vs main-95927d22: short 1343/1344 before 852e8af8a (the one failure was that name, fixed and re-run: the
audit/pilot/raster subset 31/31), long 13/13; anchordiff 0 moves in both.  Open: re-measure the q3 aug pilot (the prediction: the
2,610 s bucket falls to about one build per basis, ~6 s x builds).

Row as it stood:

· **The short-range local-PP block is rebuilt on every basis build, and per k** (lane A, 2026-11-07, measured on the bcc Li q3 aug pilot, main 09eb3caa, 4 threads, `~/Code/qchem6-runs/Li/mu/pilot2/Li_bcc_q3_aug_pilot.log`).  `setup: local-PP SHORT (analytic lattice sum)` is the run's largest bucket: 2,610 s over 432 calls, 6 s each, = 6 builds (run setup, the re-convergence rebuild, 4 audit references) × 72 IBZ blocks.  `MakeLocalPPShort` is a pure function of (molecular basis, structure, PP, k phase) and never reads the rasters, and `CompareAtFixedD` hands each reference the run's own basis and PP (only the egg-box remakes it).  Two fixes, in order (controller 2026-11-07): (1) build once and hand the blocks to every rebuild and reference (allowed under the lazy rule: fixed inputs), saving ~5/6 = 2,175 s (36% of that run), expect 0 moves; (2) the real-space lattice sum is k-independent, so sum it once R-resolved and phase it to each k (`LatticeSum1E::MakeLocalGaussian`; check whether S, T and the other 1E lattice sums share the per-k pattern, which would make it one fix for every 1E term), gated by a unit test that the R-resolved route equals the per-k one block by block to roundoff on a small cell with several k; the summation order changes, so name and re-judge every test that moves (pin 10) · `src/BasisSet/Gaussian/Lattice/Imp/GPW_Evaluator.C`, `src/Calculation/Imp/SolidCalculation.C`

## ✅ Grid-pilot beta-test finding (e) on Li: the audit's reference RSS refusal (lane A, 2026-11-08, d54bec926)

The q3 atom L28 run was OOM-killed (globally) in a postSCF audit at 4.3 GB RSS.  Measured (lane A 2026-11-07, a temporary
VmRSS/VmHWM probe around each reference build, malloc_trim at entry; tables in the shared notes `laneA_e_rss_calibration.md`): an
audit reference's cost is per k BLOCK and grid-sized.  On Li q3 aug k12 (72 blocks) the Becke audit alone added 2.0 GB to a
1.3 GB run.  Fix (controller ruling 2026-11-07): `qchem.AuditFootprint` sizes every reference before it is built. It keeps three
named pieces apart: the byte count (Phi tables on the XC mesh plus one raster field per block); the eta/kappa calibration (1.09,
156 B per raster point, from 20 real audits on Li k444/k666/k12 and NaF k222, with the record at its declaration); and the band
(0.78 - 1.50; 0.78 rather than the ruled 0.79 because the record's minimum is 0.788, accepted by the controller).  REFUSE ON RESOURCES when the low
end exceeds the memory available (the nearest cgroup limit minus its usage, else MemAvailable; the row names its source), WARN AND
RUN when only the high end does.  The pilot is checked before the first SCF (no iteration spent).  A mid-pilot refusal stops it NOT SETTLED with its
record kept.  The postSCF audits return `Outcome<GridAuditResult, ResourceRefusal>`.  Policy row `AuditReferenceMemory`.  Expert
deck key `solid.auditFootprint`.  Not calibrated: a polarized run, a PlaneWave Vxc fit, a non-imposed Becke mesh.  The
mid-pilot refusal paths (round > 1, a climb rung) have no deterministic test; the two doors that do (preflight, postSCF) are
`AuditFootprint.*` in UTCalculation.

## ✅ D-DEGEN-FILL: a partly filled degenerate level at kT = 0 (lane A, 2026-11-08, 675bf14b8)

Row opened by lane S (2172f88d7, on its 3b branch): at kT = 0 a merged degenerate level that nothing stated was filled orbital by
orbital in the eigensolver's order, so D depended on roundoff (`PW_NaF.Γ_Anchor` 4 e⁻ into a 3-orbital level as 2+2+0; on lane
S's 3b tree the run never converged).  User rulings 2026-11-08: keep the OccupationPolicy framework -- each option its own
`OccupancyRule`, chosen by a value in `OccupationConfig`, reachable from the deck, the default a facade policy reporting itself as
one `choose` row; HF defaults to the integer fill ("does one of B give the lower symmetry, lower energy? If so then B": Boron
UHF/dzvp shared −24.40769 vs integer −24.525932, the same-spin self-interaction a fractional HF occupation pays); option C (an HF
warning row) dropped ("C will confuse users"); a STATED configuration (MOM, a held fill, Atom_EC per m) wins exactly.

Landed: `IntegerOccupancy` / `SharedDegenerateOccupancy` (`OccupationConfig::degenerate`, width = `SCFParams::MergeTol`); the one
kernel `qchem::Orbitals::ShareDegenerateLevels` for the per-block fill and the ranked molecular reservoir (a level split over
irreps is one level); `SCFParams::DegenerateLevels` = deck key `degenerateLevelFill` (Auto | OrbitalByOrbital | Shared), Auto
decided in `qchem.DegenerateFillReport` (DFT and solids share, HF integer).  Found on the way: GDM's `ComputeStep` took one held
integer step on a non-idempotent D' (now refused, like `Ready`); a declined direct minimiser's adaptive linear fallback is capped at
`DeclinedDirectMinRelax` 0.25; the iterator's constructor fills a fresh seed under the default rule, so a Shared run refills it at
the first Iterate (re-running Init instead pushed the seed Fock into DIIS twice: NaF "converged" at iteration 1).  Two fill-order
traps the landing sweep caught: the reservoir hands orbitals over irrep by irrep (water C2v HF filled 5A1 over 1B1, +3.24 Ha), so
it states its aufbau as a −energy rank; and a rankless fill must keep the STORED order, which for GDM is its leading occupied block
(an energy re-sort moved the MnO dimer off GDM's manifold).

Gates: `PW_NaF.Γ_Anchor` IonicSAD 16 / Uniform 28 iterations, one state (1.2e-9 Ha apart, test tightened to 1e-7);
`GPW_Al.Γ_Imp_Converges` −1.9212152506 in 9 (main −1.9212127589 in 14, the state lane S's 3b reaches); Boron UHF Ladder / PureGDM
on the −24.5259316 anchor; `M_Calculation.BoronLDASharedLevelGDMDeclines` (the GDM interplay).  Sweep (vs lane-a-e): short 1369/1369,
moves Si2_PP_U.LargeSeparation −0.605 (a BoundedSCF 2-iteration iterate; its claim reads Enn), GPW_Al −2.49e-6, GPW_NaF Becke
+1.9e-9 and ExactResume −2e-10 (noise), iterations only on PW_NaF 108→88, GPW_SiBox eqUnfolded 122→75 / eqFinite 74→52; long 13/13,
iterations only NaF Becke-mesh 31→29, LiF XCMeshAudit 42→54; panel at 4 threads vs main-c153c0fb: only the MnO dimer GDM case
moved, restored bit-identical by its deck stating OrbitalByOrbital; TiO₂ k222 59 vs 49 as on main (now on QCHEM_PANEL_RED).
Remainders → OpenWork §3: HF + stated Shared starts from an integer sibling seed; a stated GDM on a solid meets the shared default.

## ✅ A STATED GDM meets the shared degenerate fill: guide the user (lane A, 2026-11-09, caba39d09)

Row opened from D-DEGEN-FILL: the solid default shares a partly filled degenerate level at kT = 0, GDM cannot rotate a
non-idempotent D', so a deck stating `accelerator: GDM` never ran GDM.  User ruling 2026-11-08: "Many users will have no idea
that GDM only works for integers. So we need to guide them."  Landed in `qchem.DegenerateFillReport`, both facades: a STATED GDM
with the fill unstated takes `OrbitalByOrbital`, its `choose` row teaching why and how to get the symmetric state; GDM and
`Shared` both stated file a `warn` row `DegenerateLevelFillUnderGDM` naming the fallback that runs instead, at its relax
(`DescribeDeclinedFallback`).  Unit tests `DegenerateFillPolicy.*` (six, molecule and solid).

Pulay-D fallback cap: NOT capped, measured.  Shared + GDM (declining every step), stated Pulay-D depth 8, alpha 0.45 / 0.3 / 0.25 /
0.2: Si triplet box (`Si_box16`, SIPP_SR) 12 / 11 / 11 / 11 iterations; MnO dimer (`MnO_box16`) 67 / cap 200 / cap 200 / 57.  The
step does not decide a Pulay history's convergence, and a 0.25 ceiling would fail the MnO case.  The linear-fallback cap
(`DeclinedDirectMinRelax` 0.25, CLAUDE.md) stays as it is.  Logs: `~/Code/qchem6-runs/{Si_box16,MnO_box16}/degen_gdm_shared_statedpulayd_a*.log`.

Seen on the way: the MnO dimer under Shared at 0.45 / 0.2 converges to −30.3485166 / −30.3485171 Ha, below CP2K OT −30.34790; the
two capped runs end lower still (−30.34914, not converged): a different state, relevant to the "GDM crawls" row.  The solid banner
prints the unresolved alpha (0.5) while the Auto Pulay-D runs at 0.45 → new OpenWork §3 row.  Sweep vs main-8473b2e0: short
1418/0, long 14/0, 0 moves (panel not run: no iterator change, and its one GDM deck states its fill).

## ✅ The raster star-average's ops live on the fit BASIS (lane S, 2026-11-09, d75fed558 + the seed fix)

The row as it stood in `doc/OpenWork.md` §3:

| **The raster star-average's ops live on the fit BASIS** (lane S, 2026-11-08; restated after the user's ruling that an ID names the basis, never its role) | On GPW the Hartree CD fit basis and the XC fit basis are the SAME basis (same {G} ball, so the same `BasisSetID`, and every cached table is rightly shared), built as two `PlaneWaveFit_IBS` instances that differ in ONE member: `itsDirectOps`, the crystal direct ops, given only to the XC instance (`tGPW_IBS::CreateVxcFitBasisSet` passes `SymmetryOps()`, `CreateCDFitBasisSet` none).  The ops are read only by `G_RasterTransform::Symmetrize`, called from the k-composite's `GetRhoOnGrid` (`CompositeCD.C`) on the raster ρ_DM.  So two instances of one basis behave differently, which an ID-keyed cache cannot see; harmless today only because no cached 3-centre build reads the ops (`MakeOverlap3C` slices to the `PW_Grid_Evaluator`) | Split what the basis KNOWS from what the run CHOOSES (user 2026-11-08): the basis keeps the TRANSFORM RULE, "how do my functions transform under g", as a capability that takes the op as an argument; only the stored ops LIST (which group is in use) leaves the basis, owned by the crystal / k-composite and passed in.  The seam already exists: the tier-1 hub `qchem.Symmetry.Representation` (`Image(g,·)` / `AdjointImage(g,·)`), realized by `DenseRepresentation::MakeRepMatrix(g)` for Gaussian sets (shell rotations × centre permutation) and `MonomialRepresentation::MakeMonomial(g)` for plane waves (`PlaneWave_IBS`: G → W·G with the τ phase, a phased permutation).  The PW FIT basis today is `IdentityRepresentation` (E only, "no matrix-image client of a fit set").  So: give `PlaneWaveFit_IBS` the plane-wave `MakeMonomial(g)` (the same rule `PlaneWave_IBS` uses, shared, not copied); the k-composite averages over ITS ops through that one rule, for the raster ρ (`Symmetrize`, a voxel permutation + glide, is the same action in real space) and for ρ̃ (D2's `PairProjector::StarAverage`); character-table row assignment keeps asking the same face.  Then the two GPW instances are identical and the second can go.  ✅ APPROVED by the user 2026-11-08 ("Yes, that looks very well thought out"); lane S schedules it in or right after V-ONE-OF-N 3c, as its own commit.  **Refined (user 2026-11-08): this is the FIRST PIECE of `doc/OOD-SOLID-Cleanup.md` V-SYM-SEQUESTER-OPS** (op lists stop travelling).  So the k-composite holds the crystal's GROUP (a qcSymmetry object), not a vector of ops, and asks it to star-average: `group.StarAverage(field, rep)`, passing the field and the basis's `Representation` rule; the group iterates its own ops.  D2's `PairProjector::StarAverage(Projections&, const std::vector<ReciprocalOp>&)` changes the same way (it takes the group, or the group calls the projector's rule).  Dependency direction as V-SYM-SEQUESTER-OPS says: qcMesh/qcBasisSet → qcSymmetry, never the converse.  Record: this row |

**Done (user ruling B, 2026-11-09: a group plus a field face; "for the purpose of interface sanity we can think of molecules as k=0 systems").**  `qchem.Symmetry.StarGroup` holds its own structure-neutral `SymOp`s and runs `StarAverage(StarField&)`; the field supplies its image under one op (`StartStar` / `AddImage(g)` / `FinishStar(1/|G|)`).  ρ̃ is a `GMapStarField` (deferred ops, the old m-outer loop); the raster is `BasisSet::G_RasterStarField::MakeStarField` (a new non-hub face, realized by `PlaneWaveFit_IBS`).  One G-index rule, `Symmetry::Lattice_3D::ImageOfGIndex`, serves the ρ̃ image and `PlaneWave_IBS::MakeMonomial`.  The k-composite holds `SymmetryImposition::DensityGroup()`; `recipDensity` is gone.  HUB TRIM (user 2026-11-09, "Delete it"): `G_RasterTransform::Symmetrize` left the G_FieldEvaluator hub; the snapshot was regenerated.  The approved route through `Representation` was replaced by this face, because the hub acts on matrices, and a star-average acts on a one-index field that can reach outside any function list.

**The one move, traced:** the uniform seed's composite (`Seed.C`) held no ops, yet its raster was averaged (the basis carried the ops) and its ρ̃ was not.  With the ops off the basis it lost both; it now holds the run's group, so both are averaged.  Sweep vs main-daaaa917: short 1430/1430, long 15/15 with 0 moves; 13 GPW imposed Γ tests move by at most 7.4e-7 Ha, by ±1-2 iterations (one iteration's ΔE at the gates' convergence; the baseline spreads 1.1e-6 between eqUnfolded's two arms).  Every other piece is bit-identical by construction (UT `PlaneWaveFields.TheRasterAndGMapImagesOfAnOpAgree`).  Remainders: OpenWork §3 "Star-average remainders".


## ✅ E is evaluated on the idempotent D_out (lane S, 2026-11-09, d21377785 + bab40220d + a4afb9989)

The row as it stood in `doc/OpenWork.md` §3:

| **PW ball-route XC is non-variational: ±1e-8 Ha noise floor** (lane S, 2026-11-08) ⛔ | On a cold, all-full plane-wave run the ball-route XC energy is not variational in the density, so the total energy wanders ±1e-8 Ha between runs that differ only by roundoff (measured on the V-ONE-OF-N 3a binary too: pre-existing, not 3b).  Two 1e-8 gates sit on that floor and are NOT loosened (user ruling 2026-11-08): `SolidState.ExactResume_PlaneWave` and `Deck.ASavedStateRestartsAndTheRevisionNamesItsParent`, both listed in `QCHEM_KNOWN_RED` (root CMakeLists.txt, `ctest -L knownRed`).  When either goes red (passes its gate) it comes off the list and this row closes | owner V-XC-ONE-TERM: make the ball-route E_xc the functional the V_xc is the derivative of (one sampling for both), then re-run both gates · record: this row (lane S 3b sweep classification, group B) |

**Re-diagnosed: the XC was never the cause.**  `PW_Terms.TheXCMatrixIsItsEnergysDerivative` (715c81b9f) checks, on a complex PW
block at fixed D, that FD of E equals Re Tr(HΔ) to the FD error: Dirac 7e-10, LDA 1.3e-10, Hartree 3.6e-11.  The ball route IS
variational.  The floor came from WHERE E was evaluated: under D-mixing (`LinearMixer`, `PulayDMixer`) `Mix` overwrote the working
D with the mix, and `tSCFIterator` then called `TotalEnergy(itsCD)` on it.  A mix of D_out's is not idempotent, so E there is FIRST
order in the mixing error: ±1e-8..1e-7 Ha at Δρ ~1e-8.  `Deck.ASavedState...` (a GPW run, not PW) ended at -7.1150678405; its own
restart at [F,D] 7e-15 read -7.1150678654, 2.5e-8 lower.  A field mixer (Kerker/Pulay ρ̃) always kept its own running field, so its
E was already E[D_out].

**Fix (d21377785).**
- `OwnedDIn<T>` (`LinearMixer.C`): a D-mixer keeps the mixed D_in as its own untracked copy and hands it to the next Fock through
  `FockDensity` only for the D_out it mixed (a `Version()` match); otherwise the working density.  At relax == 1 nothing is stored.
- STRUCTURAL GUARD: `tDensityMixer::Mix(const cd_t& working, ...)` -- no mixer can write the SCF's density.  `LinearMixer::Adapt`
  re-mixes its own D_in on a retraction (algebraically the old quarter step), and the iterator no longer recomputes E.
- ASSERT (live in RelChecked): `tSCFIterator::EnergyDensityIsTheOutput()` -- the density E is evaluated on equals this iteration's
  C f C† output (‖ΔD‖ ≤ 1e-10 N).  With in-place mixing restored it fires on `Deck.ASavedState...` (exit 134); silent on the fix.
  It compares densities through the existing `GetChangeFrom`, so the protected `ChargeDensity` hub is unchanged.
- E CONSUMERS all read E_out now, because itsCD is never mixed: the ΔE gates and `SetEnergy` (the Ladder's ΔE/E hand-off), `GetEnergy()`
  (the facades' `Energy()`, `ReportResult`, the saved-state summary), and GDM's start.

**Gates.**  Both knownRed cases PASS and left the list: `Deck.ASavedState...` 9 iterations to -7.1150678654 = its restart (main 15);
`SolidState.ExactResume_PlaneWave` 16 iterations (main 74, a ±1e-7 stall).  `SCFTrace` goldens re-captured (bab40220d): only E and its
derived columns move, [F,D] / Δρ byte-identical; E[D_out] is second order (Si Γ iteration 4: -7.1150545 vs -7.1140720 on the mix).

**Anchor moves** (anchordiff vs `main-08d9536c`, tol 1e-10; controller accepted): short 68 moved, all |ΔE| ≤ 1.7e-8 Ha, all on
D-mixing runs (atoms and molecules at relax < 1 D-mix too); molecular moves ≤ 2e-10 at unchanged iteration counts.  Iteration counts
mostly DROP (ExactResume_PlaneWave 74→16, PartialEigensolve 31→15, ten SolidState gates 15→12, LittleGroup / OrbitalImage 28→12,
LR 163→138).  Long: 6 moved, ≤ 5e-10.  Two cases it EXPOSED (no numerical defect; rows in §3, a4afb9989):
`GPW_Si.Γ_U_Ladder_LinearResponseUIsAFixedPoint_Panel` (the Ladder's ΔE/E hand-off: panel-red) and `A/A_SL_HF_ion.A/60` (a noise-floor
[F,D] gate on cond(S) = 1e15: knownRed, trim first).

**Cost, and the follow-up that removed it (controller ruling 2026-11-09, option B).**  As first landed, a D-mixing run refreshed
the Hamiltonian at BOTH D_out (the energy) and the mixed D_in (the next Fock): one extra COLLOCATION of ρ[D_out] per spin per
iteration and nothing else (serial, the two builds side by side): MnO AFM-II Γ PulD 6.1 → 8.1 s/iteration (+33 %, collocations
28 → 48 over 12 iterations), Si k222 PulD and every Kerker run unchanged (a field mixer's Fock is its own presented field, so it
always collocated D_out alone; the §5a rows and the banked solid recipe are Kerker).  An energy-only evaluation at D_out would not
help (E_H and E_xc need ρ[D_out], which IS the cost); carrying a combined field through the mixer's Axpy needed a projector
identity on a hub face, and recovering the mixer's coefficients inside the collocation memo by least squares was rejected (a
tolerance, a guess at a relation the mixer knew, and the memo's bit-identical replay contract broken).

The ruling: when a D-mixer's NEXT Fock is not D_out, the trace reports the first-order energy
E_T = E[D_in] + Re Tr(F[D_in](D_out − D_in)) (`tSCFIterator::TaylorEnergy`; F[D_in] reaches the density's own static contraction
through a `tStatic_CC` adapter on the Hamiltonian at D_in, so spin, k and real TRIM blocks are the density's own walk), and the
loop refreshes ONCE at its end for the exact E_KS[D_out] (final E, saved state, the RelChecked assert, post-SCF consumers), printing
|E_T − E_KS|.  Where the next Fock IS D_out (a field mixer, a relax-1 passthrough, direct-min) the exact energy is free and stays.
**E_T steers, E_KS judges** (controller, 2026-11-09).  STEER (read E_T on a D-mixing step): the MinΔE gate, the Ladder's |ΔE/E|
hand-off (`SetEnergy`), the virial gate, the trace and the observer, and through the rows the run's factor-scale TREND checks (charge
slosh: Eee against its run minimum; order collapse: the site moment below 1 % of its peak, sustained), which a mixed D_in following
the D_out's carries too.  JUDGE (exact E_KS[D_out], refreshed where it judges): the iteration-1 refusal (`FirstEnergyRefusal`: one
extra refresh per run; E_T is a few Ha off there by construction, O₂ HF 3.4 Ha), the final E, the saved state's summary, the RelChecked
assert.  Neither: the GDM line search (its own E(t=0)) and the linear mixer's `Adapt` (only [F,D]).  E_T is formed in the driver
while both densities are live and ACTIVE (before `installNew` supersedes D_in: the Hamiltonian refuses a superseded density, which
the RelChecked tree caught), and only when the mixer will not pass D_out straight through (`tDensityMixer::PassesOutputThrough`:
a relax-1 linear mixer, every molecular default, keeps the free exact energy and pays nothing).
- GATE `UTSCFIterator TaylorEnergy.*` (diagnostic `energy_taylor` prints E_KS[D_out] beside E_T, changing no number): Si k222 PulD
  slope 2.09 over 8 iterations, gap/‖D_out−D_in‖² 0.032–0.058 while the gap falls 3.8e-3 → 4e-10 Ha; O₂ HF triplet (polarized,
  exact exchange) slope 1.86, ratio 0.05–0.22 over 3.4 Ha → 5e-11 Ha; final |E_T − E_KS| 2.7e-15 / 5.7e-14 Ha.  At Si Γ the step
  direction rotates between iterations (ratio 0.0006–0.056): second order as a BOUND, not a clean line -- why the gate is k222.
- TIMING (MnO PulD, serial): collocations 26 per 12 iterations (28 before), XC samplings 26 (27), 6.5 s/iteration against 6.1.
  The remainder is one XC-matrix quadrature per spin: the adapter re-assembles F[D_in] blocks the wave function already holds
  (OpenWork §4 row).  `E_T` is evaluated BEFORE the residual: the residual samples D_out on the XC mesh, whose sampler caches one
  density, and E[D_in] after it re-sampled D_in (27 → 38).

## ✅ Fixed-size tables indexed by a runtime L: a THROW at the entry, not an assert (lane H, 2026-11-09, 98e91b4bb)

The row, verbatim:

| **Fixed-size tables indexed by a runtime L are guarded by ASSERT only** (F, contract; lane H 2026-11-04, sibling survey of the Hermite3 overrun) | `Boys::AuxillaryFJ::theLookUp[121][23]` (assert `Fj.size() < thejMax+6` is looser than the table: jmax ≤ 16 is the real bound); `Wigner3j` (constexpr LMax 4); `RelWigner3j` (KMax).  All dead in Release (NDEBUG), so an over-L caller reads past the table silently, exactly as Hermite3 did.  `PG_Cart_MnD::kMaxPoly` is SAFE (declines at run time); Hermite1 now throws at construction | Make each a THROW at construction/entry naming L and the bound (Hermite1's pattern), or size it from L (Hermite3's); one UT per table that asks for L = bound+1 · HubbardUHistory.md 2026-11-04 |

Landed (CLAUDE.md "Where to check, and with what": the bound checked ONCE at the entry, thrown, never assert-only):
- **Boys `AuxillaryFJ`** takes jmax at construction and throws `std::out_of_range` for jmax outside [0, 16].  The real bound is the
  table's: the Taylor step reads columns jmax..jmax+6 of `theLookUp[121][23]`, so jmax <= 16 (a `static_assert` pins the width).
  The old `assert(Fj.size()<thejMax+6)` was dead in Release AND looser than the table (it let jmax reach 20, column 26 of 23).
  `RNLM` constructs it first, before indexing its static `theRjs[thejMax-Max]`.  Reachable by Ltotal = La+Lb+Lc+Ld > 16 (four h shells).
- **`Wigner3j`** (constexpr table, LMax 4): `Wigner3j::CheckL(l)` throws, called once from the `NR_Angular` and `RKB_Angular`
  constructors (one per irrep basis), so the per-lookup asserts in the angular loops stay asserts.  Found on the way: `RelWigner3j`'s
  KMax = 5 admits κ = +5, whose l = 5 indexes the integer 3j table past LMax on the RKB path -- the RKB constructor's `CheckL(l(κ))`
  closes it.  A `static_assert` ties LMax to the `Fact` table the Racah core indexes (n = 4 LMax + 2 at worst).
- **`RelWigner3j`** has no table of its own; its κ/k bounds are the `Fact` table's reach (`static_assert`).  `RelWigner3j::Check`
  throws at each lookup (not in a hot loop: only the unit tests call it today).
UTs: `M_BoysTable.RefusesJBeyondTheTable` (and `RNLM(17, ...)` throws), `M_BoysTable.ValuesAtTheBound` (F_j(0) = 1/(2j+1) and
F_j(5) against the defining integral, j <= 16), `RelWigner3jTests.TablesRefuseBeyondTheirBound`.

### ✅ HF with a STATED shared degenerate fill starts from an integer seed (lane A, 2026-11-09)
FIXED in `tSCFIterator::Initialize`: an HF run seeded through its LDA sibling (a matrix-free SAD seed) now sets
`itsFreshSeedFill`, so a stated Shared rule refills the sibling's integer-filled orbitals at the first Iterate, as it
does for any fresh seed.  Gate `M_Calculation.BoronUHFStatedSharedTakesTheSharedState` (Boron UHF/dzvp, SAD seed
STATED): −24.4076937683 (7 iterations) with the fix, −24.5259316138 (8, the integer state) without it, so it
discriminates.  The seed must be stated: HF's default seed is the core guess (MoleculeSetup's AutoSeed row), which
never takes the sibling path, so the first gate drafted on the default seed passed with the fix disabled.
Seen on the way: O₂ singlet HF + Shared from SAD lands on the broken −149.5438874 either way (its π_g splits beyond
MergeTol at that fixed point); from the core guess it reaches the shared −149.2962220 → new §3 row "A stated Shared fill
does not guarantee the symmetric state".

### ✅ A solid run's banner prints the UNRESOLVED mixer alpha (lane A, 2026-11-09)
The diagnosis was WRONG: the banner already received the resolved params.  The cause was a cout precision leaked by an
earlier `setprecision(1)`, which printed the Auto Pulay-D step 0.45 as `alpha=0.5`.  FIXED: `EmitSCFBanner` states its
own format (defaultfloat, precision 6) and restores cout's.  The remaining leak sites → new §3 row "cout precision
leaks between prints".
Sweep for both (vs main-cf1c8366 short/panel, main-08d9536c long): short 1443/0, 0 moved (+1 new gate); long 15/0,
0 moved; panel|periodic at 4 threads 25/25, 0 moved.

### ✅ A molecule's wavefunction cannot be saved or restarted (user 2026-11-09; lane A, V-FACADE-STATE)
DONE as ruled (doc/OOD-SOLID-Cleanup.md V-FACADE-STATE, now in CleanupHistory3).  The shared half is the new module
`qchem.SCFState` (src/Calculation/SCFState.C): the framing, the base fingerprint keys, the comparison, the writer and reader, the
Löwdin carry (generalised from the solid's) and the restore checks.  The abstract `Calculation` gained a public non-virtual
`SaveState(path)` and `JudgeState(path)` (the refusal or the warm list, without restarting: the binding's ask), a public
`RestartedFrom()` (exact / warm lines / displacement / Löwdin-carried / the electrons the restored D holds on the run's basis), the
protected template method `RecordState(StateRecord&, kT)` (the base writes the structure, electron and block keys, each facade
overrides, calls it and adds its own), and `JudgeRestart` (kind check, comparison, carry, the EXACT RESUME / WARM START lines).
`MoleculeCalculation::Restart` / `PrepareRestart` and `AtomCalculation::Restart` / `PrepareRestart` mirror the solid's: an exact
resume seeds `SeedRole::ExactResume` and a Ladder starts on its GDM rung at kT = 0 (a `startOnGDM` key on the json accelerator
door); a molecular setup with an unstated, unsized XC mesh adopts the saved state's mesh when the eps agrees (a `policy` row).
FORMAT: `qchem-state 2.0`, root attr `kind`, per block an `irrep` label and the overlap `S`; the reader still takes the solid's
`qchem-solid-state 1.x` (as kind solid; no S, so a moved-atom restart from one carries verbatim, as it always did), so the
vendored decks/states/ files resume unchanged.
Gates: FacadeMatrix SaveRestartRoundTrips green for Atom (O, LDA) and Molecule (H2O PP LDA, the fresh setup adopting the mesh);
M_Calculation `BoronUHFExactResumeStartsAtTheSavedAnswer` and `WaterLDASymmetryExactResumeStartsAtTheSavedAnswer` (first iterate
within 1e-9 of the saved E, <= 2 iterations); `WaterRestartAtMovedAtomsCarriesTheDensityThroughTheLowdinFrame` (one H moved
0.064 bohr: Tr(D S') = 10 to 1e-12, then the fresh run's energy to 1e-8; DISCRIMINATES: the verbatim carry holds 9.99506);
`ASavedStateIsJudgedByClass` (another basis REFUSED by name, HF -> LDA a named WARM start, a molecule state to an atom refused by
KIND); UTCalculation `SolidState.AVendoredFormat1StateReadsThroughTheSharedReader`.  Found on the way: a by-value
`RestartedFrom()` dangled in a range-for over `->warm` (segfault): it returns a reference.
Follow-ups (OpenWork §2/§3): the binding + GUI (lane G) and the molecular deck keys (V-FACADE-OPTS); the geometry session's
verbatim carry; a SALC basis identified by name, irreps and sizes only.

**✅ The moved-atom carry, MEASURED: the verbatim AO carry stays, for file restarts and geometry sessions alike (lane A,
2026-11-09; was the OpenWork row "A geometry session carries D verbatim: Tr(D S') != N").**  The row's premise was wrong.  The
V-FACADE-STATE landing (035a9900f) carried a moved-atom restart through the Löwdin frame to hold Tr(D S') = N, and the same was
asked of `MoleculeSession`.  Measured on H2O dzvp, one step from the converged reference, first-iterate E1 - E_final (Ha):

| carry | Tr(D S') | LDA O 0.02 | LDA H 0.064 | LDA H 0.137 | HF O 0.02 | HF H 0.064 | HF H 0.137 |
|---|---|---|---|---|---|---|---|
| verbatim AO | 10.006 / 9.995 / 9.995 | 1.3e-5 | 4.6e-6 | 2.0e-5 | 5.8e-7 | 1.7e-6 | 1.6e-5 |
| Löwdin frame, A = S'^-1/2 S^1/2 | 10 | 1.4e-4 | 9.0e-6 | 1.4e-4 | 3.2e-6 | 4.3e-6 | 4.5e-5 |
| natural orbitals re-orthonormalised in S' (= McWeeny of the verbatim D) | 10 | 1.2e-4 | 2.4e-5 | 1.2e-4 | 2.1e-6 | 1.5e-6 | 1.3e-5 |
| projection S'^-1 <new\|old> C, re-orthonormalised (controller's request) | 10 | 1.5e-4 | 1.9e-4 | 7.8e-4 | 6.3e-6 | 5.5e-6 | 3.3e-5 |

Every N-exact carry starts LDA about 10x further from the answer; HF is a near-tie with the natural orbitals.  The fixed-N
Gaussian Coulomb fit does NOT explain it: with the Hartree term on the mesh (`CoulombFit::Delta`) verbatim still wins by the same
factor (1.4e-5 vs 1.2e-4 to 1.6e-4).  The projection, which keeps the density where it WAS in space, is worst: the density
follows its nuclei, and so do the AO coefficients.  Iteration counts were recorded too but are accelerator noise (a DIIS period-2
cycle under LDA; 11 to 85 for the same start quality).  The solid side had found the same earlier (`SolidSession::Carry`, AO
default, 2026-10-30: Si 41 vs 42, NaF 110 vs 96 -- no lever).  Landed: one shared `CarryToMovedAtoms` (qchem.SCFState: D kept as
AO coefficients, the new S attached; the table at its declaration) used by `Calculation::JudgeRestart` and by the session's
step (`MoleculeCalculation::CarryOntoThisBasis`, which also files the step in `RestartedFrom()`); `RestartSummary::lowdinCarried`
is gone, `electrons` stays as a diagnostic; `LowdinCarry` remains as the solid session's `Carry::Lowdin` option.  Gates on the
FIRST ITERATE (within 1e-5 Ha of the converged energy; measured 1.7e-6): UTCalculation
`MoleculeSession.EveryStepsFirstIterateStartsNearItsAnswer` (two steps) and ITMain `M_Calculation.WaterRestartAtMovedAtomsStartsNearItsAnswer`;
with the carry switched off the first iterate is 4.1 Ha out and the session gate fails.  The earlier V-FACADE-STATE entry above
(its Löwdin carry and Tr = N gate) is superseded by this one.

## ✅ GeometryOptimizer rejected a roundoff-level energy rise after the gradient had converged (lane H, 2026-11-09, 07f5e89e5 + 7bbeb3647)

The row as it stood:
| **`GeometryOptimizer` rejects a roundoff-level energy rise after the gradient has converged** (lane H, 2026-11-09, the `WyckoffMap` rutile model) | at max\|g\| 2e-9 the next step's ΔE is below E's roundoff (E −4.53 Ha); a 1-ulp "rise" rejects the step and cuts the trust radius ×¼ ten times until `TrustCollapsed`, on a point whose gradient (8e-16) had already converged.  A real SCF point's energy noise is its convergence tolerance, far above 1 ulp, so any maxStep tighter than the step that noise allows fails the same way | accept a step whose rise is below an energy noise floor (the surface's own ε, or `energyChange`) when the new point meets the force criteria, or test convergence at the trial point before the rise check; a UT: the quadratic bowl at maxStep 1e-7 converges.  For V-OPTIMIZE-DIP (`doc/OOD-SOLID-Cleanup.md`), which owns the algorithm half · `src/Calculation/GeometryOptimizer.C` (Run) |

**Closed by V-OPTIMIZE-DIP §8 steps 1-2** (`doc/Records/VOptimizeDIP_Design.md` §4).  The algorithm moved to qcOptimize's `BFGS_RFO`
with three noise rules:
1. a rise within δ(here) + δ(there) is accepted;
2. the trust ratio is judged only when the predicted change exceeds the noise;
3. every evaluated point is judged by the client's rule.  A REJECTED trial whose gradient meets the rule is the new
   `Why::GradientDisagrees`, never a radius collapse.

δ is the cost's `ValueNoise`, or 16 ε_mach |value| when the cost states none.  The geometry client's `PESCost` states
`OptimizerOptions::energyNoise`, which `FromEps` sets to eps.

Gates:
- `WyckoffMap.RutileHasOneFreeParameterAndTheOptimizerFindsItsMinimum` runs at maxStep 1e-6, without the 1e-4 workaround;
- UTOptimize `Minimize.ARiseWithinTheValueNoiseIsNotARise`: a variational value noise of 1e-6.  Declared, it converges in 23
  steps; the undeclared twin collapses at step 30.


## ✅ cout precision leaks between prints (lane A, 2026-11-10, 64db52231)

The row as filed:

| **cout precision leaks between prints** (lane A, 2026-11-09, from the solid banner) | several prints set `setprecision(1)` (some with `scientific`) and never restore it, so a later print inherits it: the solid banner read `alpha=0.5` for 0.45 (fixed: it now states its own format), and `[lattice sums] ... reach=3e+01` still prints at the leaked precision.  Sites: `SCFIterator.C` :116 :942 :1055, `Response/Imp/LittleGroup.C` :31, `FieldMixer.C`, `SCFAcceleratorDIIS.C`, `SphericalQNs.C` | every print that changes cout's format restores it (a small RAII `FormatGuard` in qchem.Reporting), or prints through an `ostringstream`; the console renderer of the run report already does · the sites above |

**Closed.**  `qchem.Reporting` gained `report::FormatGuard`, an RAII object that saves a stream's flags, precision and fill
and puts them back when its scope ends.  One is opened at the top of every scope that sends a format manipulator to cout,
or to an `ostream&` it was handed, about 35 sites in 27 files.  Prints that saved and restored precision by hand
(`prec0`) also got one, because `prec0` missed the flags and the early returns.

Seen in the sweep:
- `[lattice sums] ... alpha_max=4e+01 ... reach=3e+01` now prints `alpha_max=36 ... reach=30.3485`.
- Two lines had been blessed WITH the leak, and now state their own format:
  - `[energy] final E_KS[D_out]` read `-1.496295470276e+02`, because the trace row's `scientific` leaked onto it.  It now
    prints defaultfloat at 13 significant figures.  The SCFTrace goldens were re-blessed for that line only.
  - `GPW_Mn2Box.Γ_Becke_Shub_Pol_Smear_KeepsOrder`'s own result line read 10 digits only while an earlier print leaked
    `setprecision(10)`.  anchordiff saw a 1.8e-5 "move" that was 6-digit rounding.

A print that relies on cout's ambient format is the same defect from the other side: a value line states its format.

Not done:
- A runtime detector, for example a gtest listener that fails a test which leaves cout's format changed.  It would
  catch a new leak only on paths the tests reach.
- Math's own `operator<<` writers already save and restore by hand, and are below qcCommon, so they keep their own code.

## ✅ The Ladder's |ΔE/E| hand-off fires on a lucky tie, starves on a cycle (lane A, 2026-11-10, 61c7c7a94)

The row as filed:

| **The Ladder's \|ΔE/E\| hand-off fires on a lucky tie, starves on a cycle** (lane S, 2026-11-09, exposed by E[D_out]) | `SCFAcceleratorLadder` hands DIIS → GDM on \|ΔE/E\| < 1e-6 alone.  `GPW_Si.Γ_U_Ladder_LinearResponseUIsAFixedPoint_Panel` (Si +U, linear D-mixing α 1): its [F,D] / Δρ trace is identical before and after E[D_out] through iteration 26.  On main it handed off at 27 because two MIXED-D energies tied to 2e-7 while [F,D] was still 1.3e-2 (a lucky tie), then GDM finished at 78.  On E[D_out] the DIIS limit cycle (period ~20: iterations 7-12 ≈ 27-32) shows its real ±1e-4 swings, the hand-off waits to 55, and GDM's ~50 linear iterations meet the 100 cap.  **On the panel-red list** (`QCHEM_PANEL_RED`, a4afb9989).  ⚠ The 1e-6 hand-off is the USER's own CLAUDE.md rule ("switches to GDM when ΔE/E<1e-6"): changing its criterion goes to the user with lane A's measurement | owner SCF robustness (lane A, after V-FACADE-STATE): gate the hand-off on [F,D] (or require both) and measure on the panel; take the case off `QCHEM_PANEL_RED` when it goes red · record: this row; traces `~/Code/qchem6-runs/anchors/lane-s-dout_short.log` vs `main-08d9536c_short.log` |

**Closed.**  The tail criterion was not the defect.  The STALL trigger, which exists for "DIIS ran out of steam", was gated on
`DIIS::Exhausted()` (B singular or not extrapolating for 3 steps), and a DIIS that extrapolates fine while going nowhere never reaches
it.  User ruling 2026-11-10 ("yes that sounds like a good exit policy"):
- STALL fires when the best [F,D] has not improved by 0.1 % for 5 DIIS steps, with |ΔE/E| > 1e-8 still required, and without `Exhausted()`.
- The best-[F,D] tracker starts at the first REAL error.  It had been seeded from the err = 0 read before any orbitals, so
  `itsNoImprove` counted from iteration 1.  `Exhausted()` masked this; un-gating without the fix (arm S) handed off at iteration 5
  everywhere and broke the SCFTrace golden.
- The CLAUDE.md tail rule (|ΔE/E| < 1e-6) is unchanged.
- `GPW_Si.Γ_U_Ladder_LinearResponseUIsAFixedPoint_Panel` left `QCHEM_PANEL_RED`.
- New UTs: `SCFAcceleratorLadder.ALimitCycleHandsOffAfterStallStepsWithoutABetterError`, `AStalledErrorOnASettledEnergyIsAPlateauNotAStall`,
  `TheErrorBeforeAnyOrbitalsIsNotABest`.

The measurement (main c70781405 plus a probe; logs `~/Code/qchem6-runs/anchors/ladder-handoff/arm_*.log`).  Summed SCF iterations
per test; `!` = the 100 cap.  Arms: E = today (|ΔE/E| < 1e-6); F x = [F,D] < x alone; B x = both; T = E plus the un-gated stall (landed).

| case | E | F 1e-2 | F 1e-3 | B 1e-2 | B 1e-3 | T |
|---|---|---|---|---|---|---|
| Si Γ U Ladder LR (panel-red) | 136 (h66, cap) | 79 (h14) | 100! | 100! | 100! | **90** (h8) |
| NaF Imp eqColdStart | 40 | 45 | 57 | 45 | 57 | 39 |
| NaF Imp eqExactResume | 38 | 44 | 56 | 44 | 56 | 37 |
| NaF Imp Anchor | 21 | 21 | 25 | 21 | 25 | 20 |
| NaF Becke Imp eqUni | 12 | 12 | 19 | 12 | 19 | 12 |
| Boron UHF Ladder GDM | 11 | 11 | 7 | 11 | 7 | 11 |
| SCFTrace blocked molecule | 22 | 23 | 23 | 22 | 22 | 22 |

The other Ladder cases (LiBox, SiBox ×2, SolidSession, FacadeMatrix, the SCFTrace energy tail) were identical in every arm.  On Si the
DIIS rung is in a period-~20 limit cycle with [F,D] between 1.3e-3 and 1e-1, so no threshold on that signal is a reliable trigger:
|ΔE/E| fires on a lucky tie, [F,D] < 1e-2 on a lucky trough, and [F,D] < 1e-3 or "both" never fire.  [F,D] alone also costs NaF 5-17
iterations, because GDM finishes well from [F,D] 1e-2..4e-2.

Landing sweep `lane-a-ladder1`: short 1498/1498, long 15/15, panel|periodic (4 threads) 25/25; RelChecked builds.  Every WILL_FAIL case
failed for its own known reason (the iteration cap on Na2Box, MnO2 ×2, MnOBox and A_SL_HF_ion 60; Si one-spin m ≠ 0; MnO ordering;
TiO2 53 > 49, the same on main).  anchordiff against main c70781405 moves only the four Ladder cases above, by at most 1.5e-9 Ha.

Caveat seen in arm T, not the hand-off's fault, logged beside "An SCF energy can jump with Δρ blind to it": NaF Anchor's cold SCF stopped
8.9e-7 Ha above arm E's stop, on a GDM step with |ΔE/E| 1.9e-10 while [F,D] was 3e-3.  The restart lands both arms on the same energy to 1.5e-9.

## ✅ Three SCF gate rows, landed 2026-11-02 in d0567558a, closed 2026-11-10 (lane A)

The rows as filed:

| **SCF convergence gates have no real "off": 1e30 is a sentinel** (user via lane G, 2026-11-01: "the real fix is to actually have a proper off flag for each rather than some heuristic guess") | `SCFParams` MinFD / MinΔFD / MinΔE / MinVirial / MinΔρPerElectron use 1e30 for off; `IterationRow` carries it as each column's tolerance, so the GUI draws a gate line at 1e30, and `SCFIterator.C:295` / `RequireLiveGates` / the SCFTrace golden headers test `< 1e30` | lane A: each gate becomes `std::optional<double>` (nullopt = off); the column factories (FDCol etc.) pass it straight to `ColumnData::tolerance`, which already prints no threshold and sends none across the bridge; replace every `< 1e30` test; deck/JSON "off" spelled explicitly; re-bless the SCFTrace headers.  GUI needs no change (no GUI-side filtering -- ruled out as a heuristic).  pybind/ may need a flag to lane G if SCFParams is bound · `src/SCFIterator/SCFParams.C`, `src/SCFIterator/Imp/SCFIterator.C` |
| **The molecular SCF iteration row still prints Δ[F,D]; it should carry ΔE/E** (user via lane G, 2026-11-01: "delta [F,D] should be replaced by delta E/E for iteration log and SCF graphs") | Δ[F,D] was retired as a default gate (1e0b2069), but the molecular console layout still prints it; the periodic layout already prints ΔE/E.  The GUI plots whatever cells `IterationRow` carries, so it needs no change | lane A, with the 1e30 off-flag row: replace the Δ[F,D] column with ΔE/E in the molecular layout, re-bless the SCFTrace goldens · `src/SCFIterator/Imp/SCFIterator.C`, `src/SCFIterator/tests/SCFTrace.C` |
| **Linear D-mixing has no floor on α; a stopped SCF must propagate** (user via lane G, 2026-11-01, fighting a benzene SCF: "the max/min mixing factor (alpha)") | `SCFParams` has StartingRelaxRo / MaxRelaxRo (both in the GUI, lane-g 62071c53) but `LinearMixer` multiplies α by 0.8 on every [F,D]-worsening step with no lower bound.  Separately, the GUI Stop button throws a bridge-local exception from the row sink through `Calculation`/`SolidCalculation::Converge`, which today lets it through unchanged | lane A, with the SCF gate rows: `SCFParams::MinRelaxRo` (default 0 = no floor) through `LinearMixerFactory`, clamped where `itsRelax *= 0.8`; deck key `minRelaxRo`; the GUI adds an "α min" row.  Add a UT that an exception thrown from the row sink propagates out of Converge untouched (Converge must never swallow std::exception).  Open question for later: should a user cancel be an Outcome (a legitimate stop) rather than an exception? · `src/ChargeDensity/Internal/LinearMixer.C`, `src/SCFIterator/SCFParams.C` |

**Found done** while picking lane A's next work (2026-11-10): d0567558a ("SCFParams gates become std::optional; ΔE/E on the
molecular console; MinRelaxRo") did all three.
- `Gate = std::optional<double>` with `Meets`, and decks spell an off gate "off".
- ΔE/E is every layout's column.
- `MinRelaxRo` floors `LinearMixer`'s retraction.
- A UT checks that an observer exception leaves both facades' Converge untouched.

The rows were never moved.  One leftover was fixed on closing: `FacadeMatrix.C`'s SolidProbe still assigned `1e30` to three Gates,
which with `std::optional` made them LIVE gates at 1e30 rather than off; they are now `std::nullopt`.
The third row's open question (should a user cancel be an Outcome?) stays in the live tracker as its remainder.

## ✅ SALC state identity: a molecular SALC state identifies its basis by name, irreps and sizes only (lane A, 2026-11-10, 76bf2f986)

The row as filed:

| **A molecular state identifies a SALC basis by name, irreps and sizes only** (lane A, 2026-11-09, from V-FACADE-STATE) | an AO-blocked molecular basis records its shells (REFUSE), but a symmetry-adapted one records the basis name, the block irrep labels and sizes: a `.bsd` edited under the same name with the same block sizes would restart without a refusal | record the UNSYMMETRIZED basis's shells beside the SALC blocks (the SALC builder holds them) · `MoleculeCalculation::RecordState` |

**Closed.**  `BuildBasis` hands back the raw basis that `PG::SymmetryAdapt` adapts (the shared_ptr the adapted basis already
keeps alive).  `MoleculeCalculation` holds it as `itsAoBasis`, and `RecordState` records its shells exactly as an AO run's:
REFUSE `basis.shellL` / `shellComponents` / `shellPrimitives` / `exponents` / `coefficients`, WARM `shellCenters`.  No new
face on any basis class: a SALC block cannot honestly answer `AoShellSource` (its functions are combinations), so the facade,
which built both, holds the pair.

Coordination with lane G's `basis.fileChecksum` (basisFile runs): the shells come from the basis the run BUILT, so a stated
file records that file's shells and the two keys agree.  An edited file is ONE refusal listing both keys (the checksum also
refuses a comment-only edit the shells cannot see: it is the stricter of the two, never the opposite verdict).
A SALC state saved between 2026-11-09 and this commit refuses on the new keys; none is vendored.

Gate: `M_Calculation.ASymmetryBlockedStateIdentifiesItsBasisByItsShells` (H2O dzvp: the SALC and AO states carry equal shell
keys; before the fix the SALC state carried none).  Sweep `lane-a-salcid1`: short 1513/1513, long 15/15, 0 anchor moves.
