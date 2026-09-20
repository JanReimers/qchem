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

Steps 1–4 are closed (KP-0 2026-09-09; the CleanupCandidates sweep + sprint S 2026-09-15; TE 2026-09-15; the
plan polish 2026-09-16 — records in `doc/Records/OpenWork_History3.md` and `_History4.md` §"THE QUEUED
PROGRAMME").  **A fresh session starts here:**

### ▶ 5. DFT+U — programme step 5, `doc/Records/ParallelAndOraclePlan.md` Phase 3

**RULED 2026-09-16 (user, on the two papers in `~/Code`: Macke/Timrov/Marzari/Colombi Ciacchi, JCTC 2024,
"Orbital-Resolved DFT+U" `ct3c01403.pdf`; Agapito/Curtarolo/Buongiorno Nardelli, ACBN0, `1406.3259v3.pdf`):
+U is built ORBITAL-RESOLVED as the formulation — U is a vector over (site, shell, site-group irrep) — and
the shell-averaged Dudarev form is the special case where every U_i is equal.**  Pin 23.  The deciding
observation (user: *"pretty much decides it for me … I have seen other examples where O played an unexpected
role in TMOs"*): in Macke et al. the correction that opened β-MnO₂'s gap was on **O-p_z, not Mn-d at all**
(their §4.3, Table 5), and correcting FeS₂'s hybridised e_g was what wrecked its lattice parameter.  ⇒ the
Hubbard MANIFOLD is an input, never an assumption: the term takes a list of (site, shell, irrep, U) and Mn-d
is one entry, not the design.

★ **Write it SCALAR-GENERIC** — MnO at Γ is a real-TRIM run, so +U must serve the mixed corner (real block, complex density) like every periodic term; make it ONE `template<class TBlock>` body (as `Kinetic<T>`) with the real-block mixin as one line, so `CleanupCandidates.md` V1.35 collapses it for free rather than converting it (ruled 2026-09-19).  ★ **Write it against `MatrixForward<T>` and `MatrixAdjoint<T>` from the start** (user, 2026-09-11;
`MatrixIntegrator` itself was DELETED `7a41cca6`): +U is a `Dynamic_HT`, its occupation-matrix forward and its
potential adjoint are exactly that pair, and orbital resolution changes only the scalar per eigenvalue in
between — in the eigenbasis of the site occupation matrix \f$E_U=\sum_i \tfrac{U_i}{2}\lambda_i(1-\lambda_i)\f$
(Macke eq 6), the potential is built diagonal and rotated back, so mixing and forces need no adaptation.
**Our site groups do the resolution for free**: the (t2g, e_g) split IS the site-point-group irrep
decomposition (the same machinery as the site-adapted meshes), so the labels are fixed by symmetry and the
paper's eigenvalue-TRACKING algorithm is needed only when the site symmetry is lower than the split wanted.
Projectors: **Löwdin OAO** = \f$S^{-1/2}\f$ on the site block (`LASolver` forms it; Macke shows OAO beats NAO
consistently; truncation spheres are a plane-wave artefact we do not have).

> **▶ EXECUTED — increment 1 LANDED 2026-09-20 (`6aa8170d`; spec written 2026-09-19 after the tree
> reconnaissance, built as specified).  What is in the tree:**
> - **`Hubbard_U`** (`src/Hamiltonian/Internal/Hubbard.C`, module `qchem.Hamiltonian.Internal.Hubbard`): a
>   `cDynamic_HT` with the real-TRIM corner, spin-native, ONE `template<class TBlock> MakeMatrixT` body
>   (`MakeMatrix`/`MakeMatrixR` are the two one-line instantiations, as `Vee_Hartree`) — the scalar-generic shape,
>   so V1.35 collapses it for free.  `IsVirialValid()=false`.
> - **Input = a manifold list** (pin 23): `Hamiltonian::HubbardManifold{site,l,U}` on `Hamiltonian::Factory` /
>   `Ham_PW_DFT`; at the facade `SolidCalcOptions::hubbard` + `HubbardU(site,l,U_eV)` (eV outside, Hartree in).
>   `l` is never assumed 2, `site` never assumed the TM.  Increment 1 = ONE U per manifold (Macke eq 6 with all
>   \f$U_i\f$ equal); the per-site-irrep vector is the same code with a vector.
> - **The manifold's functions**: the block answers the new abstract face **`BasisSet::AoShellSource`** ("I am
>   built from atom-centred shells": `Gaussian::Point::Orbital_1E_IBS`, `tGPW_IBS` through `MolecularBlock()`, and
>   the spherical view with ITS OWN column layout) and `ShellRep::L()` names the shell's l.  ★ **CP2K's LOWDIN
>   manifold, read off `src/dft_plus_u.F`, is EVERY shell of angular momentum l on the atom** (all `nsb`
>   contractions: 8 d shells ⇒ a 40×40 Löwdin block per Mn per spin on the VA span), and `Columns` reproduces
>   that for parity — so CP2K's E_U=0.61 Ha at U=4 eV is a MECHANISM number (\f$\sum_i q_{ii}(1-q_{ii})\approx8.3\f$
>   over the 40 POPULATIONS of a diffuse 8-zeta span — see the form finding below), not physics; a physically meaningful +U names ONE d manifold, which is exactly the
>   orbital-resolved/ACBN0 direction.  ⚠ A Cartesian d (6 components, the s-contaminant among them) is REFUSED
>   with the fix in the message: run the spherical view (`GPW_SPHERICAL=1`).
> - **Projector = Löwdin, born on the pair**: per block `LowdinProjector<TBlock> : MatrixForward<TBlock>,
>   MatrixAdjoint<TBlock>`, \f$T=S^{1/2}[:,M]\f$ (`blazem::eigen` on the block overlap, geometry-fixed);
>   FORWARD \f$n=T^\dagger DT\f$ (D carries \f$w_kf_\nu\f$, the k-sum is the composite's block sum), ADJOINT
>   \f$V=TWT^\dagger\f$; `NumCoefficients` = \f$\sum_M m_M^2\f$.  The density gets the forward through its own
>   `ProjectOnto` (the term realises `Fitting::ScalarProjector`, as `DeltaScalarFitter` does), the term keeps the
>   adjoint.  **`PrepareSlots` is the geometry phase**: every block's pair is built there, because the composite
>   sizes its sum by `NumCoefficients()` BEFORE any block runs (the first bug the Si gate found).
> - **Occupations from \f$D_{out}\f$**: `EnsureOccupations` takes the `cDM_CD` itself or the DM-backed source a
>   mixed density retains (`cDM_Sourced_CD::DMSource`), per channel through `ChannelOf`; `SpinGroup::None` =
>   \f$n_{tot}/2\f$ in both channels (pin 5).  A matrix-free SEED before any block: \f$E_U=0\f$, \f$V=0\f$, version
>   NOT stamped (second bug the gate found: an empty W reached the adjoint).  DECLARED on the facade's run banner
>   (`+U: LOWDIN, shell-averaged, occupations from D_out*`) — CP2K mixes P; same fixed point, different trajectory.
> - Energy/potential in the occupation eigenbasis (`Analyse`): \f$E_U=\sum\tfrac U2\lambda(1-\lambda)\f$,
>   \f$W=\sum U(\tfrac12-\lambda)vv^T\f$, ONCE per density in `RefreshForDensity` (the eager phase), served per
>   block from the cache.  **`FreezeOccupations(bool)`** from day one (LR-cDFT perturbation runs; polaron
>   occupation control).  `Occupations(M,σ)` / `HubbardEnergy()` for the probes.  `QCHEM_U_TRACE=1` prints the
>   per-refresh `[+U]` line (N per manifold, max λ, E_U) to stdout (pin 17; gtest attaches no console).
> - **Gates**: `UTHamiltonian` `LowdinProjector.*` ×4 (adjointness, idempotent D ⇒ λ∈{0,1}, trace = charge,
>   NumCoefficients); `GPW_Si.Γ_U_ZeroUIsExactlyAbsent` (U=0 ≡ no term to 1e-10) and
>   `Γ_U_IsAPositiveSelfConsistentFunctional` (E_U = 1.5U exactly: three half-filled p eigenstates per site per
>   spin — the functional is honest on a case with a closed-form answer).  `scripts/testgrid` grew the **`model`**
>   axis (`U` = DFT+U, default LDA) after `fit`.  ctest 891/891.
> - **The oracle gate `GPW_MnO.DISABLED_Γ_U_Shub_Pol_Smear_CP2K`** (VA span + spherical view, the Shub anchor's
>   recipe, two arms U=0 / U=4 eV on both Mn d, ~6 min each): compares the U-INDUCED SHIFT ΔE and E_U against
>   CP2K's +0.61735 / 0.60951 Ha (the 100 mHa absolute offset between the codes is U-independent to first order).
>   ★★ **FIRST RUN (Dudarev form): E_U = 0.080 Ha against CP2K's 0.610 — and the cause is CP2K's FORM, not
>   ours**: `dft_plus_u.F` keeps ONLY THE DIAGONAL of the Löwdin block (`IF (isgf == jsgf)`), i.e. the 40
>   populations, no eigen-decomposition — not the rotationally-invariant Dudarev functional.  On this span the
>   populations are all fractional (Σq(1−q)≈8.3) and the eigenvalues near-integer (N_maj=4.80, max λ=0.997).
>   **USER RULING 2026-09-20: the diagonal-population form is a `CP2K_COMPAT` member** — `RunPolicy::HubbardEigen`,
>   knob `QCHEM_U_EIGEN` (default on = Dudarev; off under the umbrella), `doc/Benchmark.md` §2 row 8, the first
>   PHYSICS deviation on that list; read once by `Hubbard_U`'s constructor; the gate sets it through the N5
>   hatch.  Both numbers are banked in `doc/Records/CP2Kresults.md`.  **VERDICT on CP2K's form (2026-09-20): PASSED** — ΔE = +0.6002 Ha (CP2K +0.6174), E_U = 0.5873 (CP2K 0.6095), ΔE−E_U 0.0129 (CP2K 0.0078), **24/22 iterations with/without U on the deck's own measure and loop shape (CP2K 104/44)**; tolerance 50 mHa on a 0.6 Ha effect with the codes 100 mHa apart absolutely.  Step 5.1 is BANKED.
>   ⚠ **The first pass PASSED WITH A DEFECT** (same day): `PolarizedMixCD` — the polarized ρ̃-mixed density the Fock build sees — carries no D and answers no source face; the term asked the TOTAL, got nothing, zeroed n on every Fock build and applied V = U/2·P.  The `[+U]` trace showed it (n=0 on alternate refreshes) and the oracle agreement did not.  Fixed: the DM is resolved PER CHANNEL (view → `cDM_CD` → `DMSource()`), a source-less channel keeps its n; pinned by `GPW_Si.Γ_U_Imp_Pol_Kerker_eqUnpol` (1 s, fails on the old term).  **Lesson for the record: an oracle number within tolerance is not a mechanism check.**
>   **Run-time at parity (`CP2K_COMPAT=1`, user's ask, 2026-09-20)**: setup ~10 s vs CP2K ~5 s; **9.6 s/iteration vs CP2K 9.05** (1.06×; the +U refresh itself costs us ~15 ms, CP2K's full-matrix S½PS½ ~0.75 s/iter); neither of our free arms CONVERGED (120-iteration cap, a period-2 "ρ rotates" cycle) where CP2K's deck converges the same free cell in 44/104 with `BROYDEN_MIXING ALPHA 0.2 BETA 1.5 NBUFFER 8` = Kerker-preconditioned Broyden on ρ̃.  Our recipe keeps its history on the Fock side (DIIS) and `PulayDepth=0`.  ⇒ the like-for-like CONVERGENCE gap is on the free run, and Kerker+density-history is the suspect (user); next probe = the parity arms with `PulayDepth=8` (our `PulayMixer` IS Kerker-preconditioned density Pulay), U=0 arm first — this is step 5's item 2 (N3) territory.
>   ★★ **RESOLVED the same day — three things were fighting, and none of them was the physics** (user: *"we usually end up iterating toward 'convergence' in parameter settings … anytime CP2K finishes many minutes before ITMain my antenna goes up"*; targets (1) very close energies (2) similar convergence rates (3) similar runtime/RAM).  **(a) The MEASURE**: CP2K's `EPS_SCF` is max|ΔP_ij| on the AO density matrix (1e-6); our `MinΔρ` on Kerker was the G-space residual (1e-5 ≈ 1e-4 relative, ~100× looser) — on the anchor recipe max|ΔD| sat at **3e-3** for 200 iterations while the residual read "converged" at 43.  `SCFParams::Measure::MaxΔD` puts CP2K's measure on the loop (`GetMaxChangeFrom` on the mixable face); Benchmark rule 3f.  **(b) TWO HISTORIES**: Fock-side DIIS (the Ladder) and density-side Pulay extrapolate each other's output — with both on, max|ΔD| never left 1e-2.  **(c) MOM** hands a smeared degenerate d frontier unequal occupations that swap every iteration: D rotates at constant ρ.  The deck has ONE history (Broyden on ρ̃, `BETA 1.5 NBUFFER 8`) and no MOM.  **With the deck's loop shape and measure** (`MNO_ACC=Null MNO_MOM=0 MNO_PULAY=8 MNO_MEASURE=maxdd MNO_EPS=1e-6`, probe knobs added): imposed Becke recipe **22 / 24 iterations** (U=0 / U=4 eV) to max|ΔD| < 1e-6, energies unchanged to 4e-7 Ha, vs CP2K **44 / 104**; at `CP2K_COMPAT=1` the FREE run now CONVERGES (33 iterations, AFM-II held, m̃(q) 3.13 e): setup 1.9 s vs 8.1, **10.3 s/iter vs 8.29** (3 gathers + 2 collocations against 2 + 2 — §5f lever B, the whole residual), **wall 5:42 vs 6:14**, RSS 262 vs 217 MB.  The oracle gate now carries this recipe.  ⇒ N3 (item 2) was then MEASURED the same day and PROMOTED (§4 row): (ρ,m) 18 / 21 vs (up,dn) 22 / 24 iterations, energies identical.

In order:
1. ✅ **DONE 2026-09-20 (the block above). The ORACLE ROW FIRST — shell-averaged, because that is all CP2K has.**  `&DFT_PLUS_U` per `&KIND` with
   `U_MINUS_J` and `PLUS_U_METHOD MULLIKEN | LOWDIN` (verified in the installed 2025.2 input reference).
   MnO AFM-II, the deck we trust (`IntegrationTests/CP2K/mno_afm2_gpw_va.inp`), one `U_MINUS_J` on the Mn
   kind, run as `doc/Benchmark.md` §5a says.  Ours: the per-irrep vector with all U_i equal, LOWDIN, declared on
   the `RunPolicy` deviation line.  This banks shell-averaged +U AND validates the orbital-resolved plumbing
   in one anchor.  ⚠ Mulliken and Löwdin give different occupation matrices for the same density — match the
   flavour before comparing a number.
2. ✅ **DONE 2026-09-20 (measured and promoted, §4 row N3). N3 lands WITH +U, not after it** — charge and spin need separate preconditioning (§4 row N3): +U on an
   antiferromagnet is spin-channel-sensitive and today's mixer takes charge medicine in the magnetisation
   channel; on top of that mixer a +U bug and a mixing bug are indistinguishable.  ⚠ N3 will LOOK like a
   regression on MnO (pin 18) — the N1 detectors are what make that judgeable.
3. **THEN THE U VALUES — ACBN0-style self-consistent (U, J) from OUR OWN on-site ERIs**, the "no knob" route
   (pin 12: a per-irrep vector of hand-set U's is a grad-student knob squared).  ACBN0 evaluates the full
   Anisimov on-site HF energy (bare ERIs on the Hubbard centre, occupations renormalised by the site's
   projected charge) so it is orbital-resolved BY CONSTRUCTION and delivers Hund's **J** as well.  ★ Where a
   Gaussian-basis code has an unfair advantage: Agapito et al. had to project plane waves onto a fitted
   "PAO-3G" minimal basis to get ERIs at all; **we already own the one-centre d-shell ERIs** (the atomic
   `Cache4`/Rk machinery is exactly the four-index (mm′|m″m‴) on one centre).  ⚠ Their renormalisation is a
   MULLIKEN charge, which is basis-sensitive with diffuse functions — our whole 136-span story — so take the
   Löwdin renormalisation and MEASURE the sensitivity to the valence basis explicitly before trusting a U.
4. **The U ORACLE IS QUANTUM ESPRESSO** (`~/Code/q-e`, built; `mpirun` always): `hp.x` implements the
   LR-cDFT / DFPT determination of orbital-resolved U (Macke's method lives in `pw.x`/`hp.x`).  Compare
   ACBN0's (U, J) against `hp.x` on the same cell; a disagreement is a FINDING to record, not a bug — the
   literature already shows cRPA, LR-cDFT and ACBN0 do not agree with each other.  This is exactly the PAR
   row's trigger: a question one oracle cannot answer.
5. **Do NOT build LR-cDFT ourselves** unless 3 and 4 disagree inexplicably.  **Why it is expensive TODAY and
   why that changes** (user asked for this to be explicit): the METHOD is cheap in principle — the
   perturbation is \f$\alpha\,\hat P_{manifold}\f$, i.e. the SAME projector the +U term already owns, times a
   scalar; and the responses \f$\chi_0\f$ (first non-self-consistent iteration) and \f$\chi\f$ (converged) are
   ordinary SCF runs.  What makes it expensive is the RUN COUNT × RUN COST: (a) the perturbation must not see
   its own periodic images, so the classic recipe is a **2×2×2 supercell** (96 atoms for FeS₂/β-MnO₂ in
   Macke; 32 for MnO), (b) several α values per manifold, per site type, plus a self-consistency loop over U
   (3–4 rounds), i.e. tens of converged supercell SCFs, and (c) inverting the response matrices with the
   off-diagonal intrashell elements zeroed (Macke eq 18 — the step that keeps intrashell screening).  On our
   tree TODAY a 4-atom MnO cell is ~7 min converged, the 32-atom supercell is UNTESTED (§4 row "Size the Becke
   grid": setup is 47% of the run and scales with the mesh; PAR 2.2 / Phase 2.1 answered supercell
   CONVERGENCE, not cost), and there is no k-point parallelism (§2 row KP) — so the honest estimate is
   days of wall per material, on a code path nobody has profiled at that size.  TOMORROW: the projector comes
   free with step 5.1, the supercell cost is what items "Size the Becke grid" and KP are already attacking
   (2–3× and the k axis), and once a 32-atom cell runs in tens of minutes LR-cDFT is a SCRIPT over the +U
   term, not a capability.  What stays expensive to BUILD is **DFPT** — the monochromatic-perturbation
   linear-response solver that lets QE do it in the primitive cell without supercells; that is a genuine
   solver increment (a Sternheimer/response machinery we have no seam for) and is what `hp.x` gives us for
   free as an oracle.  ⇒ the order above: use QE for the values, keep the supercell-finite-difference route
   as the fallback we could script once the run cost is down, never build DFPT for this.

★ **Does the OOD/SOLID structure get in the way of LR-cDFT?  Checked against the tree 2026-09-16 — NO, it is
what makes "a script over the +U term" true** (user's concern, worth answering once): the perturbation
\f$\alpha\hat P_{manifold}\f$ is DENSITY-INDEPENDENT ⇒ a `tStatic_HT`, `tHamiltonian::Add(tStatic_HT*)`, nothing
else changes; \f$\chi_0\f$/\f$\chi\f$ are the +U term's own `MatrixForward` read at iteration 1 and at
convergence through `SolidCalcOptions::onIteration` (live FROM CONSTRUCTION, so iteration 1 is seen); the
zeroed-off-diagonal inversion is a `CLIapps/` probe.  **The one thing to build INTO the +U term from day one:**
a FROZEN-OCCUPATION mode — Macke fixes the Hubbard potential at its unperturbed self-consistent value during
the perturbation runs so the measured curvature is DFT-only; the potential is then built from a STORED n,
not the current one.  It is the term's own state; put it in the interface now, not later.

**Contingency (user, 2026-09-16):** if QE turns out much faster than us on 2×2×2 supercells, that is
ANOTHER optimisation campaign to run-time parity with QE — while `doc/Benchmark.md` §5a keeps its CP2K rows
as the regression anchor so one parity is never traded for the other.  ⚠ QE is PW/PAW with a per-k cost
model, ours is per-pair: the like-for-like number is "wall to a converged U on the same 32-atom cell", not a
per-routine comparison; `CP2K_COMPAT`'s deviation-list discipline gets a QE sibling.

⚠ Two caveats standing: every number in Macke et al. is PBE/PBEsol, so a like-for-like comparison of
orbital-resolved VALUES waits on GGA (§2) — LDA+U on MnO still tests the mechanism and the CP2K anchor; and
MnO is charge-transfer-leaning, so an O-p entry in the manifold list is a live question for it too, not
only for β-MnO₂.  Follow-ons that the same term seam takes: Hund's **+J** (unlike-spin term; what the
fractional-SPIN error needs — Macke's outlook, and the magnetic-coupling question) and intersite **+V**
(DFT+U+V for hybridised/charge-transfer cases; Macke §5 argues orbital-resolved U already does most of what
+V was added for).  Both are §2 rows.

⏸ **Parked decision that does not block it:** the basis-side nullable vendor that would retire `GetRhoOnGrid`'s
empty-vector-means-no-route signalling (`CleanupCandidates.md` R1.0n/R1.0o).

---

## 2. MAJOR FEATURES — capabilities the code does not have yet

Ordered by the battery roadmap (`doc/BatteryMaterialsRoadmap.md`: Tier 1 = GPW + GGA + spin + DFT+U + forces +
CE/MC on NC PPs; Tier 2 = USPP/PAW + k-point throughput).  `state` is what exists in the tree today.

| feature | why | state 2026-09-16 | next concrete action · record |
|---|---|---|---|
| **DFT+U** | essential for localised 3d — plain LDA/GGA gets voltages badly wrong | NOT STARTED; oracle validated (CP2K); the term faces (`MatrixForward/Adjoint`) exist | **§1 above** |
| **DFT+U+J** (Hund's unlike-spin term) | the fractional-SPIN error — magnetic coupling in open-shell TMOs; +U alone leaves it (Macke outlook) | NOT STARTED; ACBN0 delivers J beside U from the same on-site ERIs (§1 step 3) | same term seam as +U, one more scalar per manifold; land after the +U anchor |
| **DFT+U+V** (intersite Hubbard V) | hybridised / charge-transfer insulators where an on-site term cannot restore the bond | NOT STARTED; Macke §5: orbital-resolved U already does most of what +V was added for, so it is a follow-on not a prerequisite | two-centre occupation numbers on the same projector; decide after the O-p manifold question on MnO is answered |
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
| **k-point parallelism** (row KP) | the one parallel axis every other code has (CP2K `PARALLEL_GROUP_SIZE`, confirmed in source; QE `-nk`, VASP `KPAR` believed, not verified) | ✅ the pre-warm exists (`tHamiltonian::RefreshForDensity`, 2026-09-08 — every k-independent memo is warmed before the block loop); ⛔ one write remains INSIDE the loop (`tDynamic_HT_Imp::GetMatrix`'s `mutable CacheMap itsCache` keyed by `Irrep`, k-dependent) + `itsByL/itsByLSeen` | per-block storage (or no cache under a parallel loop) for that memo, then shared prologue → read-only parallel k loop → density reduction.  ⏸ BY AGREEMENT not now: every row we run is Γ (width 2), so the payoff cannot be measured until multi-k rows exist — *"it should stay on the list until we are suitably embarrassed"* · History4 row KP |
| **Space-group irreps as block labels** (T⋊P — the next G row of pin 14) | Tier-A BZ reduction beyond the point-group fold; the IBZ weights are per-k already | the fold under the MESH-symmetry subgroup exists (KP-0); irreps of the little group do not | lands in `qcSymmetry.Lattice_3D` + the `Gaussian.Lattice` container, touching no library boundary · `doc/Records/SymmetryUpgradePlan.md` §9 |
| **SSB descent** — the methodical symmetry route | today MnO's AFM-II is ASSUMED, never derived; the free run is first-class but blind | design in `SymmetryUpgradePlan.md` §3b; exists: `SymmetryDefects`, the ops chokepoint; does NOT exist: `Impose::Subgroup`, subgroup closure, crystal irreps, **CD persistence (nothing at all)** | step 3 cannot be a single free iteration (the symmetric solution is a stationary point; SSB is second-order) — it must measure GROWTH or CURVATURE · `SymmetryUpgradePlan.md` §3b |
| **The second magnetic material** — and the derived order parameter | only the CELL is material-specific (user, 2026-08-25); `m_stag` hardcodes MnO's two Mn sites and a 0.7-bohr offset | `siteSpins` decoration + integrated site moments both live on the run now, so \f$m_{order}=\sum_A\sigma_A\mu_A/\sum_A|\sigma_A|\f$ is derivable for ANY collinear ordering | derive it when the SECOND material arrives (one cell cannot tell material-specific values from MnO's accidents) · History4 "ONLY THE CELL IS MATERIAL-SPECIFIC" |
| **Fermi smearing, finished** | metals and the anneal | Fermi–Dirac per block + `GlobalMu`/`ShFermi` reservoirs exist; kT is a hand knob | (a) principled kT tied to the gap / DOS / a target entropy (pin 15 is the constraint); (b) Gaussian / MP / cold flavours + the ½(E+A) T→0 extrapolation reported beside E, −TS, A — decide per battery need (FD is the true finite-T physics; MP/cold are numerically nicer for a T→0 answer; ⚠ MP/cold give NEGATIVE occupations, pin 21's canary) · `doc/Records/GPWPlan1.md` "Future considerations" |
| **Molecular spin-resolved SAD seed** | the Hund-split tables exist; only the PW `SeedCD` reads them | molecular `NumericCD` still hands a polarized run ρ/2 per channel | channel-aware `NumericCD` assembly + an O₂-triplet gate (Hund-split seed vs ρ/2: same basin, fewer iterations).  Feature wish, not a defect |
| **Diffuse-basis ACTUATOR** / the vet-stage symmetric trim | the detector (`PivotedCholeskyDrops`, `basis.removed` report) landed; nothing ACTS | no `Prune(indices)` exists; ortho-time pivot filtering is the only trim | decide auto-prune (`IrrepBasisSet::Prune` + `BasisSet::Prune`, RKB prunes paired large+small) vs report-only-forever (the 80/20: the user reruns either way).  **Pin 22 governs the shape**: vet-stage, a property of S made once, whole ORBITS never single AOs, reported as a basis not indices; ortho-time filtering becomes the fallback.  Free controlled experiment: `VALENCE_LOWQ_VA` under Cartesian d is rank-deficient by 10 and auto-drops to exactly SR's 122 — a null control; runs 58–60 (O₁ p(0.18) vs O₂ s(0.15)) are the real evidence · History4 "Continuous — CLEANUP" |
| **libcint lattice engine** (`PG_LibCint` realises `LatticeSum1E`) | the faster 1E/KB/3C engine for GPW | = the ISP-split `LatticeSum1E` that is V1.38's prerequisite | `CleanupCandidates.md` V1.38 (stashed) |
| **Spherical SALC S3b** — the libcint-spherical extractor | the one empty cell of the molecular grid | S1–S5 done; the in-house spherical SALC ships without it | must match libcint's real-harmonic ORDER and NORMALISATION (a foreign convention); libcint-spherical presents AS a `PGData` with spherical components (a trap).  Genuinely separable · `doc/OldPlans/SphericalSALCPlan.md` |
| **valgen `--auto`** — rule-based valence-window refinement | the next d-metal basis | plan, not built | build when a second d-metal basis is needed · `GPWPlan1.md` "valgen --auto" |
| **The run report for the GUI** | the binding side's consumer | `qchem.Reporting` complete; the wishlist waits on a consumer | `meta` section (cheap, one `EmitSection` in `Converge`); field-metadata registry (code key / terse label / hover text, render-side only); detail-level FILTER; `Renderer` DIP split; `RollingFileSink`; `basis.removed` named by L/α/atom; `schemaVersion`; structure/symmetry/irreps/Hamiltonian sections; HDF5 sidecar · `doc/Records/RunReportPlan.md` |
| **Lattice gas / cluster expansion / MC** | Tier 1's last stage — the actual voltage curve | specced, deliberately not built | `doc/LatticeGasPlan.md` |
| **A second oracle** (QE first; Elk for all-electron; GPAW/SIESTA only as GPW-like timing peers) | a question ONE oracle cannot answer | `~/Code/q-e`, `~/Code/abinit`, VASP all built (CLAUDE.md); nothing run against them yet | ⚠ TRIGGERED ONLY by such a question — today that is §4 row "MnO accuracy" (−99.7 mHa vs CP2K, operator not named).  Also cheap and their regime: CP2K's own 32-atom MnO supercell, to close the `Benchmark.md` §7b caveat (is their 1.09× OMP a never-parallelised route or too few tasks) · `doc/Records/ParallelAndOraclePlan.md` Phase 4 |

---

## 3. NON-OOD CLEANUP — tooling, build, hygiene

(The SOLID/OOD debt list — faces, casts, term architecture — is `doc/CleanupCandidates.md`.  This is the rest.)

| item | what | next concrete action |
|---|---|---|
| **Module toolchain** | `import std;` + a modular Blaze fork — banish the preprocessor; Clang-only modules | `doc/ModuleToolchainPlan.md` (not started); `doc/Records/cmakenotes.md` has the build notes |
| **Env-knob graduation** (GPWPlan1 item 1; N5's remainder) | **~45 `GPW_*` `getenv` knobs remain in `src/`** (`grep -rn 'getenv("' src`); `raster`/`cutoffFactor` are TYPED options the `CP2K_COMPAT` policy does not reach | grade each: verification instrument / ops valve (documented env, keep) vs a setting a user might tune (a typed `SolidCalcOptions` field, or the `CP2K_COMPAT` deviation list).  Pin 12 is the criterion.  ⚠ `GPW_BECKE_L/NR/ALPHA` are consulted ONLY for arguments passed `<0` — they silently do nothing against a caller-supplied degree |
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
parity routes vs 217); **iteration counts: DOCUMENT, do not chase** (different mixers, different orbital steps,
different convergence measures — History4 "BIN 4").  Threaded (12 cores) we beat CP2K on this class of system
(4.14 s/step vs 5.90; they are 0.82–1.09× their own serial) — no parallel opportunity CP2K exploits is missing
on 4-atom cells; a 100-water box would invert it.

| row | what is open | next concrete action · record |
|---|---|---|
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
