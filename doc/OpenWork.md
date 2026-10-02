# Open Work — the live tracker (v4, slimmed 2026-10-01)

**READ THIS AT SESSION START.**  §1 NEXT (DFT+U pointer) · §2 MAJOR FEATURES · §3 ACCURACY / OPEN DEFECTS ·
§4 PERFORMANCE (summary; levers live in `doc/Benchmark.md` §10).  Debt lives elsewhere: SOLID/OOD →
`doc/OOD-SOLID-Cleanup.md`; hygiene/tooling/build → `doc/CleanCode.md`.  **If it is in none of these, it does
not need doing.**  Rules: a row = what · NEXT concrete action · the ONE record; a ✅ goes to a history file the
same day; a durable ⛔ becomes a pin.  History: v3 verbatim = `doc/Records/OpenWork_History5.md`; v2 =
`OpenWork_History4.md` (cited below as "H4" by section title); `_History1-3` older.

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

Battery roadmap order (`doc/BatteryMaterialsRoadmap.md`: Tier 1 = GPW + GGA + spin + DFT+U + forces + CE/MC on NC
PPs; Tier 2 = USPP/PAW + k-point throughput).  Full prose of every row: `OpenWork_History5.md` §2.

| feature | state | next action · record |
|---|---|---|
| **DFT+U (+J, +V)** | +U built + CP2K-gated; J/V not started (same term seam, one scalar / two-centre occupations) | §1.  +J after the +U anchor; +V only after the MnO O-p manifold question |
| **Linear response (A7) + U estimators** | R0–R3 steps 1–3 done; LRT/cRPA as native `HubbardUEstimator` strategies not built | `doc/LinearResponsePlan.md` ▶ START HERE; deferred efficiency ledger = its §3d |
| **LR on REAL TRIM blocks** | response path serves the RUN's scalar only; a complex run with real blocks is REFUSED (loud); workaround `forceComplex` (R0 χ₀ unaffected) | measure one material both ways; build real-block siblings only if cost-bound, else ride V1.35 · `LinearResponsePlan.md` §5d |
| **CP2K finite-α LRT (3rd U oracle)** | patch + Si validation DONE (χ −14.5164 vs ours −14.5117, 0.033%; `~/Code/cp2k` branch `ck-alpha`, `cp2k_ckalpha.ssmp`); ⚠ CP2K RKS `trq` is scaled by fspin=0.5 (pin 26) | (1) automate the cross-check as a gtest/gpwprobe arm; (2) U_in≠0 needs a CP2K-side V_Hub freeze (not built) · `IntegrationTests/CP2K/ckalpha/` |
| **SCF checkpoint/restart** | CK-1 ✅ (SaveState / saveStateTo / Restart, HDF5 `qchem-solid-state 1`; states in `~/Code/qchem6-runs/states/<material>/`) | **CK-2** WaveFunction read from disk (χ₀/ACBN0/gaps with NO SCF); **CK-3** warm start onto a different k-mesh · H4 "CK-1" |
| **PBE / GGA** (then PBEsol, BLYP) | NOT STARTED; `Model::PBE` throws; collocation emits ρ only | ∇ρ collocation + `∇·` term; retire `GetEpsXc()=0.75*GetVxc()` base default first (OOD §I.1); spin-native (pin 5) · `doc/OldPlans/FacadeDFTPlan.md` |
| **LibXC-polarized** | `Libxc_LDA` unpolarized-only; Factory throws for Polarized+LibXC | pass `XC_POLARIZED` through the spin-native `ExFunctional` face; gate vs `VWN5PolarizedMatchesLibxc` |
| **Hybrids / periodic exact exchange** | HF exists for atoms/molecules; none for solids (`IrrepCD<dcmplx>::AccumulateDirect` asserts) | periodic EXX first (a track, not an enum); design the canonical-pair scatter so it inherits · `doc/OldPlans/ERI4Rework.md` §9 |
| **GW / RPA** | PARKED (spectral, not voltage) | behind +U, GGA, forces; needs EXX + χ₀ product basis + ε(q,ω); oracle CP2K G₀W₀ |
| **Forces** (HF + Pulay, dV_PP/dR) | NOT STARTED; the pivot for relaxations, phonons, NEB, MD, τ(T) | design note first (which terms need a gradient face) · `BatteryMaterialsRoadmap.md` stage 3 |
| **Transport / polarons / Li⁺ conduction / MD** | NOT STARTED; all behind forces (+ constrained-occupation +U for polarons; lattice gas for Li⁺) | CRTA+Seebeck is nearly free (H(R), S(R) from `LatticeSum1E`; oracle BoltzTraP2).  ★ USER PIN: Li⁺ E_m is a SET of symmetry-inequivalent path barriers + multiplicities (path-orbit enumeration via `Fold`/`SymOp`), never one number.  MD = BOMD/XL-BOMD above the SCF with a warm-started `tSCFIterator`/`tLASolver` — NOT Car–Parrinello; mind the egg-box effect on the uniform grid · H4 rows "Boltzmann", "Polaron", "Ionic", "Molecular dynamics" |
| **USPP / PAW** | Tier 2, deliberately not started | do not bake `S=I` or `ρ=ΣDχχ` into new PP/density code (pin 21) |
| **OT minimiser** | NOT STARTED; gates perf lever C | build beside GDM under the accelerator seam, smearing-aware · `Records/SCFStrategyPlan.md` §7, `Records/OTNotes.md` |
| **k-point parallelism (KP)** | pre-warm exists; cross-k gather memo = perf lever (Benchmark §10); embarrassment measured: NiO 16 blocks ~1.2 of 16 cores | threads ON for every multi-k run (`QCHEM_OPENMP_THREADS=12`), then gather memo, then KP block loop (mind nested-OMP pin) · H4 row KP |
| **Space-group irreps as block labels** | KP-0 fold under the mesh subgroup exists; little-group irreps do not | `qcSymmetry.Lattice_3D` + `Gaussian.Lattice` · `Records/SymmetryUpgradePlan.md` §9 |
| **SSB descent / second magnetic material** | MnO AFM-II is ASSUMED; `m_stag` hardcodes MnO | `Impose::Subgroup`, CD persistence (CK-1 now exists), growth/curvature test; derive order parameter when a SECOND material arrives · `SymmetryUpgradePlan.md` §3b |
| **Fermi smearing, finished** | FD + reservoirs exist; kT a hand knob | principled kT; MP/cold flavours (⚠ negative occupations, pin 21) · `Records/GPWPlan1.md` |
| **Molecular spin-resolved SAD seed** | tables exist; only PW `SeedCD` reads them | channel-aware `NumericCD` + O₂-triplet gate |
| **Diffuse-basis ACTUATOR** | detector landed (`PivotedCholeskyDrops`); nothing acts.  DECIDED 2026-09-23: AUTO-PRUNE at the vet stage (pin 22) | `Prune(indices)` (whole orbits); ortho-time failure = an exception carrying indices; ⚠ unit tests must cover Cartesian d/f l−2 contaminants AND lattice-spacing-dependent rank · H4 "Continuous — CLEANUP" |
| **Spherical SALC S3b** | libcint-spherical extractor, the one empty grid cell | match libcint's real-harmonic order/normalisation · `OldPlans/SphericalSALCPlan.md` |
| **valgen `--auto`** | not built | when a second d-metal basis is needed · `GPWPlan1.md` |
| **Run report for the GUI** | `qchem.Reporting` complete; wishlist waits on a consumer (meta section, field-metadata registry incl. per-field "literature units" toggle, detail filter, HDF5 sidecar) | `Records/RunReportPlan.md` |
| **Lattice gas / CE / MC** | specced, deliberately not built | `doc/LatticeGasPlan.md` |
| **`MixingPolicy` from (ordering, cell, functional)** ⏸ | "user says AFM MnO, not `G0=1.0, Pulay=0`"; detectors T1–T3 built | N1/T4 shape; expert overrides behind a named policy · H5 Parked |
| **Second oracle** (QE → ABINIT/VASP all built) | TRIGGERED ONLY by a question one oracle cannot answer (today: MnO −99.7 mHa below) | `Records/ParallelAndOraclePlan.md` Phase 4; also CP2K's 32-atom MnO supercell for the §7b OMP caveat |

---

## 3. ACCURACY / OPEN DEFECTS

| row | what is open | next action · record |
|---|---|---|
| ✅ Shifted-MP fold (Si −1.02 mHa) | DONE 2026-10-01: ρ was star-averaged over all 48 ops while the k-fold used only the mesh subgroup; both now use `MapsMeshOntoItself`. Imposed = FREE = −7.867454 (CP2K −7.867437) | `Records/OpenWork_History5.md` |
| **NiO VA: pseudopotential ghost in a near-null S direction** | RESOLVED in practice 2026-09-28: `NIO_VET=1 NIO_ORTHO_TOL=1e-3` trims Ni s{0.06}, Ni d{0.18}; min eig S 9.8e-7→4.7e-3; insulator 1.30 eV (QE 2.86), E −106.2021.  OPEN: why a trimmed basis hosts a −36 Ha KB ghost at all | unit gate: KB nonlocal lowest generalized eigenvalue vs λ_min(S) on the 1e-4-trimmed NiO; then the 1.30 vs 2.86 eV basis/physics comparison · `LinearResponsePlan.md` §5b, H5 row |
| ✅ GDM on UNPOLARIZED runs | DONE 2026-10-01: capacity g read off D′ (Tr D′²/Tr D′), nocc=N/g; gradient/model step need no g (ratio invariant). Still open: DIIS→GDM ladder restart on unpolarized *solids* (GPW) not re-measured | see `Records/OpenWork_History5.md` |
| **Linear D-mixing diverges on free Si, no Fock accelerator** | `CP2K_COMPAT=1 GPW_ACC=null`: adaptive relax raises α 0.30→0.45 and never re-damps | why V1.18 does not re-damp on rising E; is fixed-α 0.4 direct-P stable? · Benchmark §10 |
| **MnO accuracy — name the operator** | VA exact-span: −99.7 mHa configuration-BLIND offset + −37 mHa FM-favouring selective bias vs CP2K | pin `GPW_XC_DM_SOURCE` first; then term-by-term vs CP2K blocks on the Mn ATOM (oracle −14.2440 / −14.658), then crystal, then `MNO_KMESH=2`; read weak-moment conclusions against the INTEGRATED moment (pin 4) · `Records/SphericalLatticePlan.md`, H4 "Step 5" |
| **N4 — cusp-deficit XC feed** (pin 18) | ρ_XC = ρ_mix + (ρ[D]_exact − ρ[D]_BL) is a design argument; `SCFParams::XCCuspDeficit` OFF; blocks lever B | (a) is ρ_mix+Δ pointwise ≥0 (`GPW_RHO_NEGATIVE` census); (b) is Δ iteration-static?  Anchor-moving ⇒ re-bank window · H4 "N4" |
| **The 136-function span** (time-boxed) | why CP2K holds the full diffuse span; screen hypothesis REFUTED | candidates: SVD-consistent F/S filtering; symmetry-inequivariant drop (vet trim is the prerequisite); project near-null out of F · H4 "Step 6" |
| **Na₂ polarized singlet won't converge from AFM seed at α=0.3** | Δρ ~1e-2 oscillation, E flat; unpol converges in 30 | explain/fix before a magnetic campaign pays · H4 "TWO SMALLER LOOSE ENDS" |
| **‖V_xc − V_xc_fit‖ study** ⏸ | Becke-vs-uniform costs are at unknown-equal accuracy | parked (user "defocusing") |
| **Ladder restart starts on the DIIS+Kerker rung** | `GPW_NaF.Γ_Imp_eqExactResume`: from a CONVERGED NaF state the Ladder's rung 0 (DIIS, Kerker α=.25) steps the density before the first energy is read -- first iterate 1.3e-4 Ha off, 4 iterations until |dE/E|<1e-6 hands off to GDM, 9 in all; standalone GDM from the same file starts at 4e-10 and converges in 2.  SAME on the coarse(Ecut 40)->fine restart: GDM 2 iterations (first iterate 4e-8), Ladder 12 (first 1.5e-4; 5 DIIS iterations wander at 1e-5 before GDM takes over); one material, one grid change -- a restart whose density moved more (geometry, U) is unmeasured  The saved state is fine (diagnosed 2026-10-02).  User: restarts should use the Ladder/GDM; "immediately" was too emphatic, but a restart need not spend 4 iterations on DIIS | decide whether `SolidCalculation::Restart` should enter the Ladder at its GDM rung (or loosen `ladderSwitchAt` for a restart); the test pins GDM today · `IntegrationTests/GPW/NaF.C`, `src/SCFAccelerator/Internal/Imp/SCFAcceleratorLadder.C` |
| **Seed library must be keyed by the PP variant `q`** | every `atomic_*_densities.json` entry already carries `q`; the lookup does not take it, so a q3 run cannot ASK for the q3 seed -- today an ambiguous (Z,functional,Nelec) only THROWS (neutral since 9f926e37, explicit `Nval` since 5ddc6990).  A throw stops a wrong pick; it does not let the caller pick the right one.  Scoped 2026-10-02, NOT started (user: make it an OpenWork item): the seed code sees only a `Structure`, which has no `q`, so it is a plumbing change, ~9 sites: `AtomicDensity` Get/Has × {Density, SpinPair} gain a trailing `int q=-1` (filter entries to that `q`; unique neutral is then `Nelec==q`); `IonicSADTargets`, `MagneticDecoration`, `MakeSeedDensity<T>` (+2 explicit instantiations), `SeedCD`/`PolarizedSeedCD` ctors gain a `qByZ` map (default empty = today's behaviour); the facade (`SolidCalculation.C`, 3 call sites) builds it from `SolidCalcOptions::species`; `SCFIterator.C:151` passes none.  ⚠ a species whose `q` has NO library entry then stops silently using another variant's seed (HasAtomicDensity false → the existing neutral-amplitude fallback or a throw): expect the sweep to flag any test whose `species` q disagrees with the library's | WHO supplies `q` is the same question as `doc/CleanCode.md` D-STRUCTDATA step 2 (`pseudopotentials.json`) and `OOD-SOLID-Cleanup.md` V-CALCNET -- build it so the `qByZ` map is what step 2's PP set hands over, then do the plumbing once; fold D-BASIS-PP in · `src/ChargeDensity/Imp/{AtomicDensity,Seed,SeedCD}.C` |
| **Anchor-moving batch (re-bank once)** | left: V1.22 (OOD), §K (CleanCode), `MinΔρ` gate (CleanCode D-MINDRHO), N4 | do in ONE window so each delta is attributable · H4 "THE ANCHOR-MOVING SPRINT" |

---

## 4. PERFORMANCE — summary (detail and levers: `doc/Benchmark.md` §10, instrument §5a)

qchem started at **67× CP2K's per-iteration cost on MnO** (573 s vs 8.5 s).  Now, serial, one thread each: ahead
per SCF iteration on 7 of 9 §5a rows; like-for-like parity row **1.13×** (≈ within 13% of CP2K on CPU; whole-run
MnO AFM-II 5m42s vs 6m14s wall at the same convergence measure); **peak RAM wins** (113–132 MB vs 217 MB on
the parity routes; 262 vs 217 whole-run); threaded 4.14 vs 5.90 s/step.  QE: no timing documented.
Open levers (one-liners; Benchmark §10): lever B (one gather/spin, behind N4) · lever C (behind OT) · cross-k
gather memo (Si 6.9× Γ→8k vs CP2K 1.03×) · Becke grid sizing (policy call) · ρ̃-sampling bucket · per-iteration
G-space folds · re-take §5a on rule 3f.  Rule: COPY the run command from §5a, never reconstruct it.
