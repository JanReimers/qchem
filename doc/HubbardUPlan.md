# Self-consistent orbital-resolved (U, J) — the live plan

Born 2026-09-23, **condensed 2026-10-01** to: what is left, and the lessons.  Everything else (the full
argument, every run, every retraction) is the verbatim archive `doc/Records/HubbardUHistory.md` — cited below as
**H §N** (= section N of the old plan, which the archive keeps) or **H·OW** (the OpenWork §1 block).  North star:
`doc/BatteryMaterialsRoadmap.md` (Li\f$_x\f$Mn\f$_2\f$O\f$_4\f$ voltage curve).  Ruling underneath: `doc/Pins.md`
pins 12, 23, 24, 25.  A7's own plan: `doc/LinearResponsePlan.md`.

★ **Why this exists (once).**  Orbital-resolved +U with U determined from the code's OWN on-site ERIs, Gaussian
basis, transition-metal oxides, 14 GB desktop.  A supercell-free U is the only kind that scales to a
**composition sweep**, which is what a voltage curve is.

---

## 1. Status (2026-10-01)

**Built** (OpenWork step 5 increments 1–3; H·OW): `Hubbard_U` term (scalar-generic, spin-native, real-TRIM,
`FreezeOccupations`); manifold = an INPUT list `HubbardManifold{site,l,U,...}` (pin 23); `ManifoldSymmetry`
per-site-irrep U vector; three projectors (Löwdin column / contracted atomic / ortho-atomic with U=0
spectators); `BasisSet::BareCoulombSource` + `ERI4Block`; `Hamiltonian::ACBN0` estimator;
`SolidCalculation::ConvergeHubbardU` outer loop; `gpwprobe mno|nio|gate1`.  CP2K shell-averaged +U oracle gate
PASSED (MnO U=4: ΔE +0.6002 vs +0.6174 Ha; 24/22 iters vs CP2K 104/44).  ctest was 908/908 at increment 3.

**Routes to a U** (H §2): (a) ACBN0 as-is — ⛔ REFUTED (overshoots 2.6× on NiO d, drives MnO U away from the
literature).  (b) ACBN0 with a screened bulk kernel — ⛔ REFUTED 2026-09-25 (needed ω's 50.6 % apart on Ni 3d vs
O 2p against matched-PP hp.x; localisation sign backwards).  (c) hp.x values as input — legitimate under pin 12,
one QE run per material per composition.  **(d) our own linear response (DFPT/finite-difference) — THE ROUTE,
ruled 2026-09-25 as A7.**  What survives of (a)/(b): the instruments (bare on-site ERIs, projected occupations,
manifold/outer-loop machinery).  No ACBN0 number is quoted as physics.

**Gates** (H §4): Gate 1 magnetic robustness — PRELIMINARY PASS (order survived 6/6 arms, none converged at
NMAX=80).  Gate 2 run sizing — NOT DONE.  Gate 3 oracle set — DONE (became A6).  Gate 4 supercell-vs-q-mesh +
k-parallelism — NOT DONE.

**Headline oracle numbers** (matched-PP hp.x LRT, U in eV, U_in = 0 unless noted; k 4³ / q 2³ except NiO/MnO k 2³;
`ortho-atomic` projector unless noted; full decks/logs in `IntegrationTests/QE/README.md` §A6):

| material | manifold | U | note |
|---|---|---|---|
| NiO | Ni 3d | **5.267** (U_in=3) / 5.434 (U_in=0) | the one SOLID gate point; different U_in ⇒ different number, never average |
| NiO | O 2p | 8.514 | replaced an old cRPA *bound* ≳4 |
| MnO | Mn 3d | 0.198 (atomic) / 0.986 (ortho) | NOT USABLE: χ/χ₀ = 0.98 / 0.955 (d⁵ Pauli-saturated) — proven not a projector artifact |
| MnO | O 2p | ~~26.56~~ → 11.150 (ortho) | atomic projector was the artifact here; redone 2026-09-30 |
| SrVO₃ / KCuF₃ / Sr₂FeO₄ | V / Cu / Fe 3d | 6.250 / 8.163 / 8.012 | KCuF₃, Sr₂FeO₄ geometry from Carta's own input files (`~/Code/reprints/materials_cloud_submission/`); their U is MLWF-based, not comparable |
| LiMO₂ M=V,Cr,Fe,Co,Ni | M 3d | 5.953 / 5.811 / 7.592 / 7.307 / 9.173 | idealised R-3m template; LiNiO₂ U_in≈0 `.save` archived `IntegrationTests/QE/checkpoints/linio2_U0.save/` |
| TiO₂ rutile | Ti 3d | 4.637 | d⁰ gap, no 2-step |
| FeS₂ pyrite | Fe 3d | **7.5351** | identical to 4 s.f. on all 4 Fe sites; ~5 h run; original console log lost to a driver-bug rerun (`batch/recovered/` holds the Hubbard_parameters.dat) |
| ZnO | Zn 3d | 35.34 | NOT AN ORACLE: closed d¹⁰ shell, χ₀≈χ≈0 |
| Si (CP2K CK-alpha) | 3p | χ only | χ_CP2K −14.5164 vs ours −14.5117 Ha⁻¹ (0.033 %) — validates our method |

**Our own LRT on MnO Mn-3d** (free AFM-II, Γ, U_in=0; logs `~/Code/qchem6-runs/MnO/`): manifold-definition
trend is monotone — ours `every` χ/χ₀ 0.693, U 2.62 eV; ours `ortho` 0.837, 1.23 eV; hp.x ortho-atomic 0.955,
0.99; hp.x atomic 0.981, 0.20.  CP2K CK-alpha MnO χ = −1.894 Ha⁻¹ (broad manifold, same order as ours).  The
remaining 0.837 vs 0.955 gap is OPEN (basis/projector detail, not chased).  Interpretation: the closed-shell
cancellation is REAL physics of the same-site channel (FD and DFPT compute the same dn/dα, measured 3e-4 on Si),
not a solver artifact.

## 2. ▶ START HERE — what to do next

Every action appears here once; if it is not here it is a finding, not a task.

**Track A — the U value (decides whether the method is real)**
- **A7 (the active build) — `doc/LinearResponsePlan.md` ▶ START HERE.**  R0 (χ₀(q) sum-over-states, NiO vs
  hp.x), R1 (molecular CPHF, HF α == PySCF 1e-6), R2 (periodic q=0 kernel, χ_LR == dn/dα to 2.7e-6) and CK-1
  (SaveState/Restart, HDF5) DONE; **R3 steps 1–3 DONE 2026-09-29** (kernel at any q); **NEXT = R3 step 4**:
  `HubbardLinearResponse(qmesh)` loops `Reference::QMesh`, prints χ₀(q), χ(q) with Krylov residuals, forms
  χ(R) and U = (χ₀⁻¹−χ⁻¹)_II over the q-mesh supercell; gate = primitive Si with k,q 2×1×1 == the 2×1×1 supercell
  at Γ (MEASURE its tolerance first, LRP §3d gate table).  Oracles: hp.x NiO 5.267 (U_in=3, frozen +U), SrVO₃
  6.2502.  NiO states saved under `~/Code/qchem6-runs/` (check `ls`; use the FREE one for FD cross-checks); NiO
  free frozen χ(q=0) = −2.0700 Ha⁻¹.  Design constraints (user): interface changes are general perturbation
  theory (must also serve MP2/Z-vector), not DFT-specific; plan interfaces before coding.  Intended scope
  decisions still to make: same-site U only (hp.x-equivalent) vs inter-site V from day one (pin 23 generalises
  one more level); wire the response weight through the existing `OccupationPolicy` (Integer vs Fermi(kT)),
  never a second metal-detector; response density projected through the SAME projector object as the ground
  state; k+q→k′+G map as a first-class mesh query (C2 adjacent).
- **A5** (optional, any time): ABINIT `lruj` (`~/Code/abinit`, `mpirun --force-mpirun`) for J — hp.x gives none.
- **A6** is complete (table above).  Remaining nothing queued; β-MnO₂ (Macke's second material) never scoped.

**Track B — can we RUN the spinels (needs no oracle)**
- **B3 Gate 2**: one converged SCF per composition (λ-MnO₂ 12 atoms, LiMn₂O₄ 14, Li₂Mn₂O₄ 16 — the last NOT
  yet in `materials.json`, JT-tetragonal, deferred) at Γ, wall + peak RSS logged.  Estimate to beat: 3–8× MnO's
  ~6 min.  Measure, do not plan on the estimate.  The gate-1 arms did not converge at NMAX=80: re-run with
  more iterations or the MnO-style mixer recipe before quoting any energy (probe: `gpwprobe gate1 <material>
  [U_eV]`, `GPW_SPHERICAL=1`; logs in `~/Code/qchem6-runs/MnO2/`, `LiMn2O4/`).
- **B4**: Li q3 basis (`valence_semicore.bsd` + a `BasisSetData` enum value) **landing together with the
  basis-vs-PP check (§4 item 2)**, then the q1-vs-q3 discriminator on the Li intercalation energy
  E[LiMn₂O₄]−E[λ-MnO₂]−E[Li].  q1 is committed and working; B1–B3 do not wait.  Li q3 recipe:
  `valgen --q 3 --shell 0:8:0.05:60 --nmax 60 --floor`, E = −4.23584 Ha (validated, uncommitted).

**Track C — infrastructure**
- **C1**: one 32-atom MnO 2×2×2 supercell SCF at Γ, wall + peak RSS — says whether the supercell route exists
  for us at all (14 GB; untested assumption).  Cheap.
- **C2**: the 1.02 mHa shifted-MP fold defect (`doc/OpenWork.md` §4) — fix before KP-1; multi-k U inherits it.
- **C3**: k-parallelism (`doc/OpenWork.md` §2 row KP) — ACBN0/LRT need k-mesh SCFs (~10 % MnO U sensitivity Γ
  vs 2×2×2); this plan is the "suitably embarrassed" trigger.

**Long unattended runs (user away Oct 6–20)**: only what B3/C1 justify — the three compositions × the outer
loop, k-mesh arms, a calibration LR run if C1 says the supercell fits.  ⛔ `scripts/memsafe -p`, never bare.

## 3. The target and its open risks

Li\f$_x\f$Mn\f$_2\f$O\f$_4\f$: compute (U, J) at λ-MnO₂ (Mn⁴⁺ d³), Li₂Mn₂O₄ (Mn³⁺ d⁴), no supercells; assign
**U per Mn SITE by local oxidation state** (pin 23 already takes a per-site list), never interpolate on x.
LiMn₂O₄ is a **transferability CHECK** (per-site U's from a charge-ordered run must reproduce the end members),
not a third node.  **U(O 2p) at all three compositions is a first-class target** (ligand-hole/oxygen redox
moves the voltage; the `orthofull` arm gives it free from the same SCF).  Collinear magnetism is enough
(user 2026-09-23).  Risks in biting order:
1. Magnetic robustness — pyrochlore Mn sublattice is frustrated: many near-degenerate collinear states, and
   different Li configurations may land in different ones, polluting CE energy differences.  Goal =
   reproducible, not true-ground-state.
2. Mn³⁺ d⁴ is Jahn–Teller (e_g¹); Macke's FeS₂ warning: correcting hybridised e_g wrecked the lattice
   parameter.  Do U at FIXED geometry first.
3. LiMn₂O₄ may not charge-order at the DFT level (then the transferability check is weak — say so).
4. O 2p is not a spectator by assumption: carry it at U=0 and measure.

## 4. Open questions / unfinished (these rot silently — clear or re-park at session start)

- Does **J** transfer like U?  Only U_eff = Ū−J̄ is quotable (our bare J̄≈7.5 eV carries eq-13 self-terms).
- Is per-site-oxidation-state U stable when two Mn are crystallographically equivalent but electronically
  not (charge ordering)?  Check imposed-symmetry machinery does not average them.
- **NOTHING CHECKS THE BASIS AGAINST THE PSEUDOPOTENTIAL**: `SolidCalcOptions::species {"Li",3}` + a `.bsd`
  keyed by element only ⇒ a q3 run can silently get the q1 `LI` block.  Fix = machine-readable per-element q
  provenance in the `.bsd` header, factory throws on mismatch; same family as OOD-SOLID row D-SEED1 (seed
  library keyed (Z, functional), first match wins).  Lands WITH B4.
- GGA before any VALUE comparison with the PBE literature (ACBN0 7.63/3.0, Macke, Carta are PBE).
- Spherical atom resolves a degenerate shell by picking orbitals; matters for the +U atomic radial on a broken
  shell (measure first; captured norm prints: NiO 0.9999, MnO VA 0.991).
- ABINIT `ucrpa` Ni-3d crashes with a Fortran integer overflow at optdriver=4 (`~/Code/abinit-runs/ni3d_noU/`);
  off the critical path (Ni-3d has hp.x).  ABINIT is paused (A2) and PAW-only (different-PP oracle).
- Follow-ons the same term seam takes (user): DFT+U+V, +U+J, and resolved versions — all from LR-cDFT through
  DFPT (`~/Code/timrov2018-DFPT.pdf`, `~/Code/DFPT1-gonze1989.pdf`, `DFPT2-baroni2001.pdf`); any empirical
  U-tuning must be planned in a much wider parameter-fitting context.  Natively supporting both LRT and cRPA
  behind `HubbardUEstimator` is an OpenWork §2 feature row.  Literature candidates list: H §7.

## 5. Gotchas / lessons (every trap that cost a wrong conclusion; evidence in H)

**Oracle discipline**
- **An hp.x number is conditioned on the state it linearises about.**  NiO's 5.267 eV is U_LR(U_in=3 eV)
  (decks inherit `U Ni-3d 3.0` from QE's benchmark); MnO decks carry none.  Read the `HUBBARD` block; quote
  U_in with the value.  (pin 24; H §1 trap 2)
- **Two comparisons, never merged.**  ours÷published-ACBN0 (1.8–2.8) measures PROJECTOR COMPLETENESS (their
  PAO-3G keeps ~60 % of the norm), not screening.  Only ours÷an INDEPENDENT oracle (hp.x, cRPA — never ACBN0)
  tests screening.  A target must never be ACBN0-derived (gate 3 was circular once; retracted 2026-09-23).
  ACBN0's N̄² is not a screening model — it vanishes as the basis completes.  (pin 25; H §1, §4)
- **Label every oracle row `matched-PP` or `different-PP`.**  ABINIT is PAW-only ⇒ second opinion; hp.x with
  `gth2upf` (`CLIapps/gth2upf.C`) UPFs is matched (our occupations agree with hp.x to 1–2 %, so the NiO
  factor-2.6 gap is the FUNCTIONAL, not manifold/projector/PP).
- **A bound is not a point.**  The "6 % one-factor agreement" lived on an O-2p *lower bound*; the real value
  (8.51) broke it.  A withdrawn objection that was conditioned on a wrong number gets re-examined when the
  number changes.
- **Same-site LRT is ill-posed for Pauli-saturated/closed shells** (MnO d⁵, ZnO d¹⁰): χ₀≈χ, so χ₀⁻¹−χ⁻¹ is a
  difference of near-equal small numbers.  Do not quote MnO-d or ZnO-d hp.x values as oracles (pin 24).
- **Matched manifold before comparing**: our `every` (CP2K every-shell) vs `ortho` (single Löwdin 3d) vs hp.x
  atomic/ortho-atomic move MnO U 2.62→0.20 eV monotonically; `Hubbard_projectors='atomic'` is fragile (MnO O-2p
  26.56 → 11.15 under ortho-atomic).  'file' (Wannier) is a different object (Carta's KCuF₃).
- **cRPA is unreliable in entangled bands** (Carta et al. arXiv:2505.03698: Sr₂FeO₄ Fe-3d 0.42 cRPA vs
  6.9–7.3 LRT).  Our ABINIT NiO O-2p cRPA 1.2 eV is suspected of that pathology; trust LRT.  The window/manifold
  is an INPUT (pin 23 addendum).  cRPA, ACBN0 and LRT are three different objects called "U".
- **ABINIT ucrpa**: run optdriver=4 SERIAL (`-np 1`; MPI segfaults); `plowan_*` vars are COMMON across datasets,
  unsuffixed; `getwfk` off-by-one; under `nsppol=2` the "Average U and J" summary is WRONG — trust the per-block
  "Hubbard cRPA interaction" line only.  OpenMP made it slower.  ~100 min/NiO.
- **Bare F⁰ comparisons**: PySCF's plain (2l+1)⁻²Σ(mm|m'm') = 24.29 eV and our eq-10 Ū_bare = 27.3 eV are
  different quantities by construction (pair-count denominator excludes self-pairs; `scripts/a1_bare_f0_check.py`
  predicts it to 0.1–0.5 %).  Our ERIs agree with PySCF.
- **CP2K CK-alpha with `MAX_SCF 1` does NOT give χ₀** (energy block reports the INPUT occupation ⇒ α·n(α=0),
  identical trq for ±α).  Dead end; `queue/03_mno_ckalpha_chi0.sh` kept for the record.
- **An oracle number within tolerance is not a mechanism check** (the PolarizedMixCD defect passed the CP2K gate
  with n=0 on alternate refreshes; found by the `[+U]` trace).  CP2K's +U form is DIAGONAL populations, not
  Dudarev eigen-decomposition: `RunPolicy::HubbardEigen`, knob `QCHEM_U_EIGEN` (CP2K_COMPAT member).

**SCF / running**
- **THE ITERATION CAP HAS CAUSED THREE WRONG CONCLUSIONS** (IonicSAD NiO "not converged" at 80 needed 114; the
  Ni²⁺ "cannot be seeded" claim was the cap; Mn³⁺/Mn⁴⁺ ditto).  Probe default now 200.  Gate-1's six arms also
  hit NMAX=80.  Never conclude from a capped run.
- **LDA NiO at U=0 LOSES AFM-II — physics, not seeding** (SAD 79 iters −109.2693031; SAD+MOM bit-identical;
  IonicSAD 48 iters −109.2693702; N↑/N↓ 4.289/4.289).  So ACBN0's outer loop on NiO is a non-result (sites decouple
  to U_eff −1.05 eV, Hartree 2.04× floor); MnO's monotone loop was a property of d⁵.  NiO's U is one-shot only.
  hp.x NiO starts at U_in=3 for this reason.  **Magnetic robustness under a changing U is a gate.**
- **The seed does not change the converged number** (SAD vs IonicSAD NiO: 31 µHa, U_eff within 0.7 %) — but
  record the seed with the table anyway (the VA/VB span lesson, pin 16).
- Free MnO: the banked Kerker+Pulay+Null recipe converges (CLAUDE.md); default Ladder/GDM/MOM does not.  Match
  CP2K's measure (max|ΔD|), ONE history, no MOM (Benchmark rule 3f).
- **Löwdin coefficients + AO-basis integrals** gave Ū=182 eV: ACBN0's P̄ is the AO-basis density on the manifold's
  functions; Löwdin enters only the charges.  Si at Γ cannot test the renormalisation (s-only/p-only orbitals).
  eq 10c carries no N̄ in the denominator — that asymmetry IS the screening.
- A Cartesian d (6 comps, s-contaminant) is refused by `Hubbard_U` — use the spherical view (`GPW_SPHERICAL=1`);
  the SR Mn span under the spherical view is the wrong basis (VA is the default for spherical MnO).
- The spherical view's `AoShell::norm` "all ones" was FALSE for raw real harmonics (1:3:4) and made a d shell
  under O_h come out as three irreps; and the AFM-II supercell's grey stabiliser is D_3d (12), not O_h — parentage
  is a property of the coordination polyhedron (pin 23 addendum).
- **MnO CP2K LOWDIN manifold = every shell of that l** (40×40 Löwdin block per Mn per spin); E_U=0.61 Ha is a
  MECHANISM number, not physics.  A physically meaningful +U names ONE d manifold.

**QE / hp.x recipes** (details `IntegrationTests/QE/README.md` §A6; `hp.x` = `~/Code/q-e/bin/hp.x`, always `mpirun`,
**never nq=1** — it returns non-results)
- hp.x needs a **q-mesh = the supercell size of the method**; a primitive-cell perturbation replicates onto all
  equivalent sites and returns Σ_J χ_IJ.  Mesh q is mathematically a supercell (supercell vs q-mesh differ in
  cost: q-mesh is N_q small problems, RAM-friendly, keeps symmetry — Gate 4/C1).
- Insulators (even nonmagnetic LiCoO₂) need the 2-step smeared→`occupations='fixed'` recipe (the GAP, not
  magnetism); metals one smeared step; TiO₂/ZnO straight fixed.  Hubbard atoms must be listed FIRST in
  `ATOMIC_POSITIONS`.  GTH V q5 needs 300 Ry; Cu q11 450 Ry; GTH 3d in NiO 280 Ry (100 Ry was 0.19 Ha off).
- **Symmetrise a real relaxed cell before debugging anything else** (hp.x `D_S (l=2) not orthogonal` on a cell
  matching to 8 figures; Sr₂FeO₄).  `gth2upf` limits: one shell per l (Sr q10 semicore silently integrates to 2 e⁻
  — use Sr q2); Cu 3d¹⁰4s¹ aufbau limit cycle fixed by a kT-anneal fallback.  Atomic-valence mapping: light-valence
  PPs (Mn q7, Ni q10, Sr q2).
- Batch driver: `~/Code/qchem6-runs/batch/run_queue.sh` once raced a duplicate rerun over a finished job; back up
  `*.Hubbard_parameters.dat` before anything can touch it.
- Run cost (serial, 300 Ry, 2×2×2 q): SrVO₃ 55 min; LiCoO₂ 1h48; FeS₂ ~5 h; Sr₂FeO₄ 11 h.
- Every deck: QE symmetry ON (`nosym` false) throughout the table; k and q grids are not uniform across rows
  (NiO/MnO k 2³; later rows k 4³).

**Seeds / basis**
- IonicSAD available for Mn³⁺ (⟨r⟩ 1.078), Mn⁴⁺ (1.021), Ni²⁺ (1.008); a partially filled MINORITY shell under a
  filled majority is a long descent, a partial majority over an empty minority is not.
- PySCF: `~/Code/pyscf-env` (venv — Python 3.14 + PEP 668 forbid system pip); `scripts/gate3_screening_test.py`,
  `gate3_omega_sensitivity.py` consume PySCF's range-separated ERIs.
