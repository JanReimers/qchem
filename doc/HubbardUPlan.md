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

## 1. Status (2026-10-01; §2 updated 2026-10-21)

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

**Headline oracle numbers** (U in eV, U_in = 0 unless noted; full decks/logs in `IntegrationTests/QE/README.md` §A6).  Three
METHODS, never merged: **DFPT** = analytic linear response (hp.x's Sternheimer over a q-mesh; ours `hubbardLinearResponse`, a Krylov
solve of the same equations), **FD** = finite difference (±α SCFs; CP2K CK-alpha, ours `hubbardFiniteDifference`).  Columns:

| material | manifold | QE hp.x U — DFPT, matched-PP, ortho-atomic unless noted, k 4³ / q 2³ (NiO/MnO k 2³) | CP2K U — FD, Γ only, LOWDIN-all · ours — DFPT, same Γ/basis/LOWDIN | note |
| --- | --- | --- | --- | --- |
| NiO | Ni 3d | **5.267** (U_in=3) / 5.434 (U_in=0) | not converged ⛔ (U_in=0 ⇒ nonmagnetic; U_in=3 AFM never settles, 7 recipes) | the one SOLID gate point; different U_in ⇒ different number, never average |
| NiO | O 2p | 8.514 | — | replaced an old cRPA *bound* ≳4 |
| MnO | Mn 3d | 0.198 (atomic) / 0.986 (ortho) | **4.94** (χ/χ₀ 0.656) | NOT USABLE: χ/χ₀ = 0.98 / 0.955 (d⁵ Pauli-saturated) — proven not a projector artifact |
| MnO | O 2p | ~~26.56~~ → 11.150 (ortho) | **7.77** (χ/χ₀ 0.439) | atomic projector was the artifact here; redone 2026-09-30 |
| SrVO₃ / KCuF₃ / Sr₂FeO₄ | V / Cu / Fe 3d | 6.250 / 8.163 / 8.012 | SrVO₃ **4.43** (χ/χ₀ 0.041; Γ-only metal) · ours 4.41–4.65 (Γ fixed points, parked) / — / — | KCuF₃, Sr₂FeO₄ geometry from Carta's own input files (`~/Code/reprints/materials_cloud_submission/`); their U is MLWF-based, not comparable |
| LiMO₂ M=V,Cr,Fe,Co,Ni | M 3d | 5.953 / 5.811 / 7.592 / 7.307 / 9.173 | — | idealised R-3m template; LiNiO₂ U_in≈0 `.save` archived `IntegrationTests/QE/checkpoints/linio2_U0.save/` |
| TiO₂ rutile | Ti 3d | 4.637 | **4.86** (χ/χ₀ 0.286) · ours **4.870** (0.286) | d⁰ gap, no 2-step.  **At Γ, ortho-atomic, both DFPT: hp.x 5.019 (χ/χ₀ 0.330) vs ours 5.742 (0.248) — χ₀ +31 %, χ −1.6 %: OPEN corner case** (§2) |
| FeS₂ pyrite | Fe 3d | **7.5351** | **6.01** (χ/χ₀ 0.119) · ours **6.007** (0.120) | identical to 4 s.f. on all 4 Fe sites; ~5 h run; original console log lost to a driver-bug rerun (`batch/recovered/` holds the Hubbard_parameters.dat) |
| ZnO | Zn 3d | 35.34 | — | NOT AN ORACLE: closed d¹⁰ shell, χ₀≈χ≈0 |
| Si (CP2K CK-alpha) | 3p | χ only | χ only (not a U) | χ_CP2K −14.5164 vs ours −14.5117 Ha⁻¹ (0.033 %) — validates our method |

**CP2K CK-alpha on the other oracle materials** (2026-10-05; Γ only, PADE, `VALENCE-LOWQ-TM` basis, α on ONE site, LOWDIN, U_in=0;
decks/jobs `~/Code/qchem6-runs/batch/queue/07–09`, shared `lib_ckalpha_generic.sh`; runs in `~/Code/cp2k-runs/<Material>/`; each whole
pipeline takes 20–40 s): TiO₂ Ti 3d 4.86 eV vs hp.x 4.637 (+5 %); FeS₂ Fe 3d 6.01 vs 7.535 (−20 %); SrVO₃ V 3d 4.43 vs 6.250 (−29 %,
RKS because d¹ is odd-electron, trq ×2).  χ₀ is α-linear (re-measured at α/4: −13.99 / −33.44 / −144.4 Ha⁻¹, identical).  Differences from
hp.x that are NOT yet separated: Γ-only vs k 4³/q 2³ (SrVO₃ is a metal — its χ₀ is a few-level DOS at Γ, the least trustworthy of the
three), LOWDIN-all-functions vs ortho-atomic manifold, unfrozen V_Hub does not apply at U_in=0.
**✅ CROSS-CODE GATE PASSED 2026-10-21** (ours on the SAME cell/basis file/LOWDIN manifold/α, `decks/campaigns/CrossCode/`; record
H §"A7 cross-code gate"): TiO₂ χ₀/χ/U = −13.917/−3.987/4.870 vs CP2K −13.993/−3.997/4.864; FeS₂ −33.270/−3.987/6.007 vs
−33.457/−3.990/6.007; d traces 2.981 vs 2.987 and 7.692 vs 7.692.  ⇒ the TiO₂ +5 % / FeS₂ −20 % gaps to hp.x are SETUP (Γ vs mesh,
LOWDIN vs ortho-atomic), not method.  SrVO₃ at Γ is PARKED like CP2K/NiO (a degenerate t2g metal: our unpolarized run lands on
several fixed points within 2e-3 Ha, U 4.41–4.65 eV, CP2K's 4.43 inside the band; polarized it is a different, magnetic state).

**Our own LRT on MnO Mn-3d** (free AFM-II, Γ, U_in=0; logs `~/Code/qchem6-runs/MnO/`): manifold-definition
trend is monotone — ours `every` χ/χ₀ 0.693, U 2.62 eV; ours `ortho` 0.837, 1.23 eV; hp.x ortho-atomic 0.955,
0.99; hp.x atomic 0.981, 0.20.  CP2K CK-alpha (LOWDIN, spans every Mn function — the broadest manifold; Γ, PADE, α on Mn1 only, `VALENCE-LOWQ-VA`): χ = −1.894, χ₀ = −2.885 Ha⁻¹ ⇒ χ/χ₀ **0.656**, U **4.94 eV** — continues the monotone trend below our `every` (0.693); compare the RATIOS, U is a small difference of two large inverses.  O 2p (α on O1, L=1): χ −1.966, χ₀ −4.477, χ/χ₀ 0.439, U **7.77 eV** (hp.x ortho 11.15).  The
remaining 0.837 vs 0.955 gap is OPEN (basis/projector detail, not chased).  Interpretation: the closed-shell
cancellation is REAL physics of the same-site channel (FD and DFPT compute the same dn/dα, measured 3e-4 on Si),
not a solver artifact.

## 2. ▶ START HERE — what to do next

Every action appears here once; if it is not here it is a finding, not a task.

**Track A — the U value (decides whether the method is real)**
- **A7 (our own linear-response U) — `doc/LinearResponsePlan.md` ▶ START HERE.**  R0–R4, CK-1 and H5 are DONE
  (2026-10-05): `hubbardLoop` with `"estimator":"linearResponse"` runs Timrov's outer loop over a q-mesh and prints
  `[U table]`.  NiO U_SC 5.727/5.744 eV (2×2×2 q-mesh) vs hp.x 5.267; SrVO₃ U(V 3d) 6.818 vs the FD-matched hp.x
  6.177 (+10 %) with χ/χ₀ matched to 1.2 % — the excess is χ₀'s magnitude, as on NiO.  Still open, none blocking:
  inter-site V_IJ and orbital-resolved channels REPORTED (χ's real-space off-diagonals already hold them); the NiO
  ground-state gap (1.30 vs 2.86 eV at U_in=3) that would move U toward hp.x; one q-mesh outer step is ~2.2 h if
  U_SC becomes routine.  Design constraint unchanged: interface changes are general perturbation theory, not DFT-specific.
- ✅ **Cross-code gate vs CP2K CK-alpha — PASSED 2026-10-21** (§1; record H §"A7 cross-code gate"): TiO₂ and FeS₂ agree
  to 6 meV / <1 meV; SrVO₃ at Γ parked.  The goal stays CORNER CASES where our U or χ disagrees with CP2K/QE, without
  getting bogged down in convergence heuristics (user 2026-10-21; why NiO-CP2K and SrVO₃-Γ are parked).
- **CORNER CASE FOUND 2026-10-21 — TiO₂ ortho-atomic at Γ vs hp.x at Γ** (both DFPT, same state, same PP and same
  pseudo-atom radials: `gth2upf` writes OUR pseudo-atom orbitals as PP_CHI; record H §"TiO₂ ortho-atomic at Γ"): d
  occupation 2.405 vs 2.208 (+9 %), **χ₀ −14.404 vs −10.969 Ha⁻¹ (+31 %)**, χ −3.566 vs −3.625 (−1.6 %), U 5.742 vs
  5.019 eV.  Not a projector SCALE (χ/χ₀ 0.248 vs 0.330).  Lead: the ortho-atomic SPECTATORS are truncated in
  `valence_lowq_tm.bsd` (Ti 4s captures 0.9675 of its norm, O 2p 0.9762; Ti 3d 0.9952), so Ti 3d is under-orthogonalised
  against O 2p.  **Discriminated the same day:** + O p 0.18 (vet keeps it in rutile; Ti s 0.05 is vet-TRIMMED) takes O 2p capture to
  0.9994, the occupation to 2.230 (+1 %), χ₀ to +7 %, U to 5.561 eV — the truncated O 2p spectator was most of it.  But one O p
  shell is worth 0.215 Ha: `valence_lowq_tm` is 288 mHa above QE's plane-wave E on this cell, so the CP2K gate validated the
  method, never the basis.  **RESOLVED the same day — NO method disagreement** (H §"TiO₂ ortho-atomic at Γ", CORRECTION):
  with f on Ti and O (`b6`, 6.9 mHa above QE's E) the occupation matches to 0.1 %, χ to 0.3 %, χ₀ +3.1 %, U 5.125 vs 5.026 (+2 %),
  and U tracks the energy gap toward ≈ 5.01 eV.  The cause was the truncated O 2p spectator plus missing angular momentum (a
  PW basis has no l limit).  Next: the k/q mesh (hp.x 4.637) on a basis with f, and S1's O 2p (per-irrep) check on the same basis.
- **A5** (optional, any time): ABINIT `lruj` (`~/Code/abinit`, `mpirun --force-mpirun`) for J — hp.x gives none.
- **A6** is complete (table above).  β-MnO₂ (Macke's second material) is now item S2 below.

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
- ✅ **C2** the 1.02 mHa shifted-MP fold defect — DONE 2026-10-01 (`doc/OpenWork.md` §3 stub; `Records/OpenWork_History5.md`).
- **C3**: k-parallelism (`doc/OpenWork.md` §2 row KP) — ACBN0/LRT need k-mesh SCFs (~10 % MnO U sensitivity Γ
  vs 2×2×2); this plan is the "suitably embarrassed" trigger.

**Step-3 readiness — the ORDER before (U, J) at λ-MnO₂ / Li₂Mn₂O₄ (§3), agreed 2026-10-21.**  Not ready yet: several items below
change what a Mn U MEANS, so doing §3 first risks numbers that are thrown away.  Items already listed elsewhere are named, not repeated.
- **S1 — finish the QE method checks**: the TiO₂ corner case above (basis scan running 2026-10-21), then **U(O 2p) on TiO₂ at Γ vs hp.x
  at Γ** (same deck, perturb an O 2p manifold; the only O data point today is MnO, CP2K 7.77 vs hp.x ortho 11.15, unexplained).
- **S2 — a Mn⁴⁺ d³ proxy in a small cell: β-MnO₂** (rutile, 6 atoms; needs a `materials.json` entry, a Mn hp.x deck, collinear AFM).
  Tests pin 24 on λ-MnO₂'s exact d count: t2g³ is a fully spin-polarised half shell, the MnO-d⁵ saturation shape.  Measure χ/χ₀ for
  the total-d AND the resolved t2g/e_g channels before believing either.  **Literature: Macke et al., JCTC 20, 4824 (2024)
  §4.3** (`~/Code/reprints/ct3c01403.pdf`; PBEsol, SSSP PPs, OAO projectors, A1-AFM in a 48-atom 2×2×2 supercell, geometry RE-RELAXED
  each outer step): orbital-resolved U(Mn t2g) 1.59 eV (gap 0), + U(O pz) 4.62 eV → gap 1.01 eV and V₀ within 1 %; shell-averaged
  U(Mn 3d) 6.34 (gap 0.02), U+V 6.76 / V(Mn–O) 0.99, 1.10 (their ref 59); expt a 4.40, c 2.88 Å.  The resolved U is well posed there
  because they invert POINTWISE (χ, χ₀ off-diagonals between orbitals zeroed first — intrashell screening kept); a full inversion
  is the pin-24 singular case, so our resolved channels must use the same inversion to compare.  The ligand that matters is
  O **pz** (normal to the Mn₃O plane, π to t2g) — for us a separate irrep of the O site's C₂ᵥ (the planar-environment fix
  `305450f4` is what makes that site group computable).  NOT a matched oracle (PBEsol, different PP, relaxed geometry): build one
  as A6 did — QE 7.6 (has `orbital_resolved` + orbital-resolved U/α in the HUBBARD card) with our `gth2upf` UPFs, LDA, EXPERIMENTAL
  geometry, A1-AFM: hp.x DFPT for shell-averaged Mn 3d, pw.x ±α (FD) for Mn t2g and O pz with 1×1 inversion; then ours.
- **S3 — fix `doc/OpenWork.md` §3 "A carried density is re-filled at kT = 0"** before any U_SC campaign: every `hubbardLoop` outer
  step is a `BuildStage` from the carried density, so on a smeared / near-degenerate magnetic state each step kicks the state
  (SrVO₃ landed on different fixed points, 2026-10-21).
- **S4 — GGA/PBE (`doc/OpenWork.md` §2 row)**, then re-validate ONE U against hp.x in PBE.  U is functional-dependent (dn/dα about an
  LDA state is an LDA number) and the voltage literature is PBE+U, so the §3 campaign runs in PBE or is redone.
- **S5 — B3 Gate 2 sizing** (Track B) decides whether **C3 k-parallelism** comes before §3: Γ-only is too crude (MnO U ~10 % Γ vs
  2×2×2) and one NiO q-mesh outer step is ~2.2 h.
- **S6 — then §3**: U per Mn SITE (pin 23) + U(O 2p) at the end members, FIXED experimental geometry (§3 risk 2), U_eff = Ū − J̄
  only (no LR J: it needs a spin perturbation we do not have; second order for voltages).  Resolved t2g/e_g U is REPORTED, applied
  only if S2 shows it well-posed.  **Print V_IJ** from the real-space χ off-diagonals (`doc/LinearResponsePlan.md` START HERE open
  item); implement the +V Hamiltonian term only if U-only voltages need it (Mn–O V is what fixed voltages in the extended-Hubbard
  literature).
- **S7 — forces (`doc/OpenWork.md` §2 row) in parallel from S4 on**: needed by the CE/voltage stage (every Li configuration
  relaxed), NOT by the U step.

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
- **TEST-NAME GRAMMAR vs the +U tests** (moved from the clean-code pass, 2026-10-02; for a Hubbard-focused session to
  decide): `scripts/testgrid --check` (grammar: CLAUDE.md "Tests — naming", `doc/Records/TestSuitePlan.md`) rejects
  three `GPW_Si` tests in `IntegrationTests/GPW/Si.C` -- `Γ_U_ACBN0_Imp_Pol_eqUnpol`, `Γ_U_Atomic3p_Imp_Pol_eqUnpol`,
  `Γ_U_ACBN0_OuterLoopFeelsTheNewU` -- because `ACBN0` / `Atomic3p` (the U FUNCTIONAL and the manifold PROJECTOR
  variant, named after the model token `U`) are in no axis vocabulary.  Either the grammar grows a "U flavour" axis
  (after `<Model>`), or the names drop the flavour (and it moves into the claim or the body).  NOT decided: it is not
  clear which is right, and it depends on how many U flavours the product space will really carry (ACBN0, screened
  ACBN0, LRT-U, cRPA-U, +V, +J ...).  Decide that first; then fix `scripts/testgrid` AXES + CLAUDE.md, or rename.
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
  identical trq for ±α).  **FIX (2026-10-05): `MAX_SCF 2` + `DIRECT_P_MIXING ALPHA 1.0`, restart from the α=0 state — iteration 2's input density IS the single diagonalization's output; use the LAST `DFT+U energy` line** (jobs `~/Code/qchem6-runs/batch/queue/03_mno_chi0_Mn3d.sh`, `04_mno_O2p.sh`, shared `lib_mno_ckalpha.sh`).  **NiO in CP2K does NOT converge**: U_in=0 reaches a NONMAGNETIC fixed point (Pulay α 0.1, 88 steps, Ni moments 0.000 — agrees with 'LDA NiO loses AFM-II at U=0'); U_in=3 gives an AFM state (Mulliken ±0.6–1.2) but 4 recipes (Pulay 0.1/0.05, kT 0.005/0.01, Ni d 0.30/0.38) swing 0.05–0.1 Ha for 100 iterations; cause (CP2K SCF vs basis noise vs the +U state) NOT separated.  Parked; trial decks `~/Code/cp2k-runs/NiO/try_[A-G].inp`.
- **Bases for the other oracle materials** (Ti V Cr Fe Co Cu Ni Li K Sr S): `BasisSetData/valence_lowq_tm.bsd` + CP2K `IntegrationTests/CP2K/VALENCE-LOWQ-TM-BASIS`, conditioned on λ_min(S(k)) ≥ 4e-4 in idealised host crystals (`scripts/basis_conditioning/`); NOT registered in `Factory.C`; committed Ni (s .06, d .18) is 1e-6 on NiO; method + valgen improvement = `doc/OpenWork.md` §3 row 'valgen is optimistic about diffuse functions'.
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
