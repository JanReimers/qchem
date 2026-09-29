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

> **✅ INCREMENT 2 LANDED 2026-09-21 (code `37c6501c` + the gates/fixes commit of this date) — the
> per-site-irrep U vector, the orbital resolution pin 23 is about.**  What is in the tree:
> - **`Hamiltonian::ManifoldSymmetry`** (exported from `qchem.Hamiltonian.Internal.Hubbard`, unit-testable
>   without an SCF): the site group on one manifold as \f$D(g)=\bigoplus_{\rm shells}(N_a/N_b)\mathrm{Rep}(b,a)\f$,
>   its irreps CLUSTERED off a generic symmetrised matrix and named by character vector, isotypic projectors
>   \f$P_k\propto\sum_g\chi_k(g)D(g)\f$ (scale read off \f$\mathrm{Tr}A^2/\mathrm{Tr}A\f$, so a complex-type pair is
>   safe), and **the U SLOTS from group theory alone**: (site irrep k, parent irrep p) with
>   \f$\dim=\mathrm{Tr}(P^{site}_kP^{grey}_p)\f$, sorted (dim, parent, irrep).  `Label(n)` eigen-decomposes the
>   density's OWN n (never symmetrised), rotates the eigenbasis INSIDE each degenerate cluster to the
>   projectors (the seed n = 0 is one cluster; n is unchanged), names each eigenvector by the slot carrying
>   most of it, and reports `purity` (min site weight) and per-slot `parentage` (min parent weight).
>   `DudarevInEigenbasis(λ, v, slot, Uirrep, U, W)` is the one-U-per-slot functional `Analyse` now calls.
> - **`Lattice_3D::SiteEnvironmentRotations(atom, shells=3)`** — the point group of the site's coordination
>   (all images within the first three neighbour distances; ops enumerated as \f$R=BA^{-1}\f$ over
>   same-species first-shell triples, kept when orthogonal and a symmetry of the whole environment).  The
>   facade fills `HubbardManifold::greyOps` from it; `siteOps` stays `SiteRotations(site, decoration)`.
> - The facade banner declares `shell-averaged` | `ORBITAL-RESOLVED*` with each manifold's `Uirrep` and
>   the two group orders; the TERM prints the slot table once at `PrepareSlots` (stdout + console, pin 17)
>   and `by slot: … purity … parentage …` on every `[+U]` refresh line.  `gpwprobe mno` knob
>   `MNO_U_IRREP=a,b,c` (eV per slot; a wrong count throws with the table in the message).
> - **Gates**: `UTHamiltonian ManifoldSymmetry.*` ×4 — a d shell under O_h (48 signed permutations) splits
>   {2, 3} with the textbook characters (C_4: 0/−1, C_3: −1/0, i: 2/3, C_2: 2/−1), projectors idempotent /
>   orthogonal / complete; an O_h-symmetric n with (U_eg, U_t2g) reproduces \f$E_U\f$ and \f$W\f$ by hand and
>   equal `Uirrep` is bit-identical to the scalar U; a d shell on a D_3d site inside O_h has slots {1, 2, 2}
>   = a1g<t2g, e_g<e_g, e_g<t2g (8 shells ⇒ {8,16,16}); the seed n = 0 is named with purity 1.
>   `UTStructure SiteGroups.*` ×2 — MnO AFM-II: `SiteRotations` 12 with AND without decoration,
>   `SiteEnvironmentRotations` 48 (Mn and O), site ⊂ environment; diamond Si: 24 both ways.  `GPW_Si.Γ_U_*`
>   unchanged (T_d p: one slot, dim 9 = 3 shells × 3).  ctest 901/901.
> - **MnO AFM-II (`gpwprobe mno`, VA span, spherical, U = 4 eV)**: `site group of 12 ops (grey 48); U slots:
>   [0] dim 7 (a1g) < t2g  [1] dim 14 (e_g) < e_g  [2] dim 14 (e_g) < t2g` — 7 d shells × {1, 2, 2};
>   `purity 1.000` on every refresh; a1g parentage 100 %; the two e_g slots 84–85 % on the MINORITY channel
>   (the trigonal field really mixes the two e_g copies) and 50–60 % on the MAJORITY, where four λ ≈ 0.999 are
>   degenerate to 1e-3 and the eigenvectors are ill-conditioned by nature.  Minority occupation sits in
>   e_g<e_g (0.22–0.30 e) not e_g<t2g (0.01–0.04 e): the σ-covalent back-donation from O-p, as it should.
>   **`MNO_U_IRREP=4,4,4` is bit-identical with the shell-averaged run** (Etot −61.07915339 both, E_U per
>   refresh identical to 8 digits); `4,5,3` moves E_U (0.1407 → 0.1458 on the first refresh).
> - ★★ **Three defects found and fixed on the way** — (a) the slot table was being built by labelling a
>   RANDOM matrix's eigenvectors, so a rank-1 a1g slot could be missed and a runtime eigenvector would
>   append a slot past `Uirrep`'s end: now group theory; (b) **the spherical lattice view's `AoShell::norm`
>   said "all ones", which is FALSE for the raw real harmonics (angular norms² 4π/15 : 4π/5 : 16π/15, a
>   1 : 3 : 4 ratio) — the site-group rep was non-orthogonal and a d shell under O_h came out as THREE
>   irreps**; `BuildCartToSphere` now returns each view column's normalisation and the table reports it, so
>   `BuildOperationRep`'s \f$N_a/N_b\f$ convention holds on the view too (the I4 lattice-SALC consumer would
>   have hit the same); (c) **the AFM-II SUPERCELL's grey stabiliser is D_3d (12), not O_h** — the tree's
>   `SiteRotations` comment asserted O_h for a week; parentage is a property of the coordination polyhedron
>   (pin 23 addendum).  The spec's "{3, 2}" is {2, 3} in the code (sorted by dimension); the printed table
>   is the order `Uirrep` follows.
> - **Not built, stated**: eigenvalue TRACKING across near-degenerate clusters (Macke's algorithm) — with
>   equal U per slot inside a degenerate shell nothing depends on it; with unequal U it is the physics that is
>   ill-posed, not the code.  The spec is in `doc/Records/OpenWork_History4.md` §"THE QUEUED PROGRAMME ·
>   increment 2 spec".  **NEXT = increment 3: ACBN0 (U, J) from our own on-site ERIs** (item 3 below).


> **▶ INCREMENT 3 — REMAINDER (a) DONE, (b) DECIDED: OUR ACBN0 ON NiO AGAINST THE hp.x ORACLE, 2026-09-22.**
> NiO is the clean test because hp.x is healthy there (5.27 eV) where MnO's d⁵ made it unreliable.  Built
> first: a valgen-validated **Ni q10 valence basis** (s 7 × {0.06..28}, d 8 × {0.18..44}; pseudo-atom
> converged, 0.56 mHa above the pool floor — the SAME block in `va` and `sph`, because the VA/VB trims are
> MnO-vs-CP2K rank statements and NiO's oracle is not CP2K), **`materials.json` `NiO_AFM2`** at a = 7.88 bohr
> (QE's benchmark cell — an estimator is compared with an oracle on the ORACLE's geometry), and
> **`gpwprobe nio`**: `RunMnO`/`MnO()` became `RunTMO`/`TMO` over a `TmoSpec`, so MnO and NiO are one arm with
> two specs and MnO's `MNO_*` knobs are untouched.  **MnO re-run as the refactor's regression check:
> U_eff 10.774 / 10.692 eV against the banked 10.77, E = −61.4112412 against the banked imposed anchor.**
> - ★★ **THE MEASUREMENT (NiO AFM-II, LDA, Γ, imposed, deck-shaped recipe, \f$U_{\rm in}=3\f$ eV):**
>
>   | route | U(Ni 3d) | U(O 2p) | \f$N_{3d}\f$ ↑/↓ |
>   |---|---|---|---|
>   | ours, ACBN0, ortho-atomic (Ni only) | **14.51 / 14.76 eV** | — | 4.920 / 3.445 |
>   | ours, ACBN0, full ortho-atomic set | **13.80 / 13.99 eV** | 7.27 eV | 4.894 / 3.204 |
>   | ACBN0 paper (PBE, Mulliken, PAO-3G, self-consistent) | 7.63 | 3.0 | — |
>   | **hp.x linear response** (our GTH UPF, ortho-atomic, k/q 2×2×2) | **5.27** | — | 4.980 / 3.368 |
>
>   Spectators (full set): Ni 4s 0.09 eV (nearly empty → the \f$d^0\f$ limit), O 2s 23–28 eV (full shell →
>   bare-like) — the same pattern MnO showed.
> - **SEED, and it does not matter** (re-measured 2026-09-23): the rows above are `IonicSAD`, the probe's
>   current default after the Ni²⁺ retraction; the same runs on `SAD` gave 14.60/14.67 and 13.89/13.91, i.e.
>   within 0.7 % (O 2p within 0.02 %), with total energies 31 µHa apart.  The number is a property of the
>   converged state, not of the guess.  ⚠ IonicSAD needed 114 iterations where SAD needed 63 and so read
>   "NOT converged" at the old 80 default — the probe default is now 200.
> - ★★ **THE MANIFOLD AGREES WITH THE ORACLE; ONLY THE FUNCTIONAL DISAGREES.**  Our ortho-atomic 3d
>   occupations are within **1.2 % / 2.1 %** of hp.x's on the same cell with the same UPF and the same
>   projector (4.921/3.440 vs 4.980/3.368) — as MnO's were (4.93/0.45 vs 4.988/0.481).  So the factor 2.6 to
>   hp.x is not the projector, not the manifold, not the pseudopotential: it is the U functional.
> - ★★ **SLICE D's VERDICT NOW HOLDS ON A SECOND MATERIAL AND A SECOND MANIFOLD.**  ours ÷ published ACBN0:
>   MnO d **2.31**, MnO O-2p **2.75**, NiO d **1.82**, NiO O-2p **2.42** — all in 1.8–2.8, none near 1.  ⚠ That
>   spread is too loose to call "a common factor 2.5", so it is quoted as a range; what it does show is a
>   PROJECTOR-COMPLETENESS effect (which is manifold-independent in sign and rough size) rather than a physics
>   disagreement (which would vary by manifold).  On NiO, where hp.x is healthy, ours is 2.6× linear response
>   and the paper's own ACBN0 is 1.45× it.
> - ⛔ **THE ORACLE ROW WAS CONDITIONED AND THE CONDITION WAS NOT WRITTEN DOWN.**  hp.x's 5.267 eV is the
>   response around a ground state **already at \f$U_{\rm in}=3\f$ eV**: `IntegrationTests/QE/NiOg.scf.*.in`
>   carry `HUBBARD {ortho-atomic} / U Ni-3d 3.0`, inherited from QE's `hp_insulator_us_magn` benchmark.  It is
>   neither a U=0 response nor a self-consistent U.  The MnO decks carry no `HUBBARD` U, so those rows ARE
>   \f$U_{\rm in}=0\f$ — the two hp.x rows were never on the same footing.  `IntegrationTests/QE/README.md`
>   now carries the qualifier and a `U_in` column.  **A response is a function of the state it linearises
>   about; quote the state with the number.**
> - ⛔ **LDA NiO AT U=0 LOSES THE AFM-II ORDER, AND THE OUTER LOOP CANNOT BE RUN ON NiO.**  At U=0 the SCF
>   CONVERGES (79 iterations, Δρ → 0) to a NON-MAGNETIC state: the integrated site moment goes 2.107 e →
>   1.6e-6 e, dead from iteration 45, with the collapse written on its face as \f$N\uparrow/N\downarrow\f$ =
>   4.2811/4.2811.  Its U_eff = 15.9 eV is paramagnetic NiO and is NOT banked; the order guard refused it.
>   **MOM does not rescue it** — MOM on gives a BIT-IDENTICAL energy (−109.2693031), i.e. it faithfully holds
>   the non-magnetic pattern it is given.  This is the textbook LDA-NiO failure and is exactly why the QE
>   benchmark starts at U = 3.  And the ACBN0 OUTER LOOP on NiO is a NON-RESULT: after the first U update
>   every SCF failed, the two sites decoupled (site 1 reached U_eff = **−1.05 eV**), the order died and the
>   Hartree term ran away (Eee 19.6 → 39.9 Ha, 2.04× its own floor).  Trajectory 14.6 → 17.7 → 4.6 → 16.6 →
>   17.4 → 19.1 — thrashing, not converging.  ⇒ **MnO's monotone 8-step loop was a property of d⁵, not of the
>   loop**: a fragile antiferromagnet does not survive a changing U on today's recipe.  One-shot is all NiO
>   has, and that is a prerequisite for any U-functional work on it.
> - ★★ **(b) THE U-FUNCTIONAL DECISION: SCREENED ACBN0, WITH ORACLE-U AS AN INTERIM BRIDGE.**
>   **(a) ACBN0 as is — REFUTED.**  It overshoots on both materials and both manifolds, and its outer loop
>   drives U the WRONG WAY (MnO 10.8 → 16.8 monotonically, away from the literature).  No ACBN0 U is banked.
>   **(c) oracle-computed U per material as INPUT — legitimate but not the destination.**  Pin 12 permits it
>   (a computed input is not a hand-set knob) and it is what a production +U run should use TODAY, declared on
>   the run banner.  It does not scale to the north-star: a cathode voltage curve needs a U per COMPOSITION,
>   i.e. an hp.x supercell run at every Li concentration.  **(b) ACBN0 with a screened interaction — the
>   direction.**  The diagnosis is specific: the bare \f$F^0\approx27\f$ eV for Ni 3d is RIGHT (it is a bare
>   integral on the atomic Slater scale), and \f$\bar N^2\approx0.64\f$ is measuring BASIS COMPLETENESS, not
>   dielectric screening — so the missing physics is \f$\varepsilon^{-1}\f$, which for a TMO is ~1/4 and is
>   the right size to close a factor 2.6.  We already own the machinery (`BasisSet::BareCoulombSource` +
>   `ERI4Block`), and a screened kernel is a NEW INTEGRAL TYPE, which the pseudo-wall pin explicitly allows.
>   ⚠ **The screening length must come from the density** (Thomas-Fermi on the valence ρ, or an RPA
>   \f$\varepsilon\f$), never be set by hand — a hand-set λ is pin 12 all over again, one layer down.
> - **NEXT:** (1) the screened-kernel slice (the testable claim: ONE dielectric factor moves all four
>   measured manifolds, and if it does not, the effect is not screening); (2) the NiO magnetic state has to
>   survive a U change before any loop or trajectory on NiO means anything — today it does not; (3) GGA before
>   any value comparison with the PBE literature (both the ACBN0 paper's 7.63/3.0 and the 4–7 eV range are PBE).

> **▶ INCREMENT 3 — (3) THE hp.x ORACLE: NUMBERS IN, 2026-09-21 (`4f992092`; recipe, decks, outputs in
> `IntegrationTests/QE/`).**  The hp.x build reproduces QE's own NiO benchmark (7.0521 vs 7.0514 eV).  The
> matched route works: NiO with OUR GTH UPFs (gth2upf: our PP, PP_CHI = our pseudo-atom 3d) gives
> **U(Ni 3d) = 5.27 eV** (LDA, ortho-atomic, 280 Ry — the cutoff GTH 3d needs; 100 Ry was 0.19 Ha off) beside
> QE's 7.05 (PBEsol, US-PP).  **MnO AFM-II, same route: U(Mn 3d) = 0.96 eV** (0.20 with the plain atomic
> projector + O 2p in the inversion), stable against cutoff and projector: χ₀ = −0.045 → χ = −0.043 — the
> d⁵ high-spin shell's bare response is 2.5× weaker than NiO's and essentially unscreened (χ/χ₀ 0.96 vs 0.62),
> the closed-shell-per-spin weakness of linear-response U.  pw.x agrees with us on the manifold itself
> (atomic 3d 4.988↑/0.481↓ vs our 4.93↑/0.45↓; ortho-atomic 4.980↑/0.253↓).
> - ★★ **THE ORACLE VERDICT ON ACBN0 (MnO, LDA):** ACBN0 one-shot 10.8 eV, self-consistent 16.8 eV; linear
>   response 0.96 eV; literature screened U 4–7 eV.  The two routes bracket the physical range from opposite
>   sides by a factor ~5 each: ACBN0's \f$\bar N^2\f$ renormalisation is not a screening model (slice D), and
>   hp.x's linear response is unreliable for a filled-majority/empty-minority d⁵ shell in LDA.  **Neither is
>   banked as MnO's U; the +U anchors keep U = 4 eV (CP2K parity).**  What IS banked: the projector, the
>   occupations, the bare on-site integrals, and the machinery (estimator, outer loop, UPF writer) — the
>   instruments a better U-functional will be judged with.
> - **NEXT (increment 3 remainder):** ✅ (a) and (b) are DONE — the block above (2026-09-22).  (c) GGA before
>   any value comparison with the PBE literature: still open.
> - ▶ **THE EXECUTION PLAN IS NOW `doc/HubbardUPlan.md`** (born 2026-09-23): the background above in brief,
>   the four routes with the SUPERCELL column that decides between them, the Li_xMn2O4 target (U per Mn SITE,
>   not interpolated on x), and the gate order — magnetic robustness, run sizing, and the one-dielectric-factor
>   test that can REFUTE the screened route before a line of kernel is written.
> - Cost record: pw.x MnO 280 Ry ≈ 4 min serial; hp.x Mn-only 2×2×2 at 280 Ry 49 min on 4 ranks; the every-
>   thing-in-one 100 Ry atomic run was a 42-min non-result.  ⛔ mpirun always; nq = 1 never.

> **▶ INCREMENT 3 — (3) THE hp.x ORACLE, IN PROGRESS 2026-09-21 (`8f5aefea`, `0b1476d9`).**
> `CLIapps/gth2upf` writes a QE UPF for one of OUR GTH pseudopotentials with PP_CHI = OUR pseudo-atom's
> orbitals, so pw.x/hp.x run the same PP and the same `atomic` projector as `HubbardU_Atomic`.  VALIDATED:
> the isolated Mn q7 atom in pw.x (spherical 3d⁵4s² fixed) gives E −14.2433 Ha vs CP2K ATOM −14.2414 / ours
> −14.2442, eigenvalues shifted by one common +0.17 eV (the box's G=0 alignment), 3d–4s splitting to 4 meV.
> **MnO AFM-II in pw.x (LDA sla+vwn, 100 Ry, k 2×2×2):** E −61.518 Ha, moments ±4.51 μB, gap 1.32 eV, atomic
> 3d occupations **4.988↑/0.481↓ vs our atomic projector 4.93↑/0.45↓** — the two codes agree on what the
> manifold holds.  hp.x needs the 2-step magnetic-insulator recipe (smearing → fixed occupations,
> tot_magnetization 0) and a q-MESH: **nq=1 returned U(Mn 3d) −0.41 eV / U(O 2p) 24.4 eV — a non-result**
> (the perturbation repeats with the 4-atom cell); the 2×2×2 q-mesh run is the oracle number (decks and
> recipe in `IntegrationTests/QE/`).
> - ⚠ Our own MnO one-shot ACBN0 and the loop are still Γ-only; the k222 sensitivity was ~10 % (slice C).

> **▶ INCREMENT 3 — SLICE D LANDED 2026-09-21 (`4a00d797`): THE OUTER LOOP AND THE ORTHO-ATOMIC PROJECTOR —
> and the verdict on ACBN0 as a screening model.**  `SolidCalculation::ConvergeHubbardU(params, {maxOuter,
> tolU_eV})`: estimate → `HubbardUEstimator::Apply` (the term's `SetU`: U, occupation version AND the matrix
> caches reset — the \f$V_k\f$ blocks are cached by density serial and the density has not changed) → a NEW
> STAGE through `BuildStage` (a continued `Iterate` let the converged Pulay history extrapolate the first
> post-U step back onto the old density and report "converged" in one iteration) → re-converge → repeat; the
> trajectory is the result.  `HubbardManifold::orthoAtomic` / `HubbardU_OrthoAtomic`: the flagged
> contracted manifolds Löwdin-orthogonalised among themselves; a SPECTATOR is a manifold with U=0, so O 2s/2p
> + Mn 4s at U=0 reproduce QE's ortho-atomic set (and get their own estimates).  `gpwprobe mno MNO_ACBN0=n`,
> `MNO_U_RADIAL=every|atomic|ortho|orthofull`.  Gates `LowdinProjector.OrthoAtomicManifoldsAreOrthonormal
> AsASet`, `GPW_Si.Γ_U_ACBN0_OuterLoopFeelsTheNewU`.  ctest 908/908.
> - ⚠ **BASIS CORRECTION (`6284e8ff`, same day, user):** every MnO number in the slice C/D blocks BEFORE this
>   line was taken on the **SR span under `GPW_SPHERICAL`** — the SR Mn block has only two s exponents because
>   its s span lives in the CARTESIAN d contaminants, gone under the spherical view: the wrong basis for a
>   spherical arm (the probe now defaults a spherical MnO arm to VA and refuses sr).  And the "atomic 3d" had
>   been the pseudo-atom RUN IN THE SITE'S OWN SHELLS, which (VA trimmed of the 0.18 d shell; SR with no 4s)
>   is a too-compact 3d — eps(3d) −0.18 (VA) / −1.04 Ha (SR) against CP2K's ATOM code −0.257.  In a 16+16
>   pool our pseudo-atom reproduces CP2K ATOM to 3/0.3/0.1 mHa (E/eps 4s/eps 3d; gate
>   `ValenceBasisGen.MnQ7PseudoAtomInALargePool`), so the radial is now the POOL orbital PROJECTED onto the
>   site's shells (captured norm printed: MnO VA 3d 0.991).  **Corrected numbers (VA, LDA, Γ, U=0):** atomic
>   3d \f$\bar U\f$ 15.0 / \f$\bar J\f$ 4.3 / \f$U_{\rm eff}\f$ **10.77 eV** (bare 24.3/6.2; \f$N_{d\uparrow}\f$ 4.93 →
>   4.02 renormalised), every-shell 8.95 eV; **self-consistent atomic 3d 10.77 → 15.37 → 16.57 → 16.75 →
>   16.785 eV in 7 outer steps** (the two Mn differ by ~1 %, an asymmetry to watch).  The verdict below is
>   unchanged in substance; the SR-span figures are superseded.
> - ★★ **THE NUMBERS (MnO AFM-II, LDA, Γ, deck-shaped recipe, imposed) — SR SPAN, SUPERSEDED, kept for the record.**  One-shot at U=0, \f$U_{\rm eff}\f$(Mn 3d):
>   10.87 (atomic) / 10.87 (ortho, Mn–Mn) / 10.75 eV (full ortho set); spectators: Mn 4s 0.27 eV (nearly
>   empty → the \f$d^0\f$ limit), O 2s 32 eV (full shell → bare-like), **O 2p 7.36 eV (paper 2.68)**.
>   **Self-consistent, atomic 3d, from U=0: 10.9 → 17.2 → 18.9 → 19.2 → 19.26 eV in 8 monotone outer steps**
>   (SCFs of 13/12/10/10 iterations after the first).  +U localises d → \f$\bar N\f$ rises → U rises.
> - ★★ **VERDICT: ACBN0's renormalisation is a projector-completeness effect, not a screening model.**  With
>   a compact atomic projector on a complete Gaussian basis, \f$\bar N\f$ is 0.75–0.85 per d state and eq 12
>   cannot bring a 29 eV bare average below ~11 eV; the paper's Mn d (4.67) and O 2p (2.68) sit a COMMON
>   factor ~2.5 below ours, i.e. their PAO-3G-projected plane-wave states carry ~60 % of their norm — the
>   "screening" is what the minimal projection basis drops.  The literature's screened U for MnO (cRPA,
>   LR-cDFT: 4–7 eV) is a response quantity; **hp.x (linear response) is the value oracle, as ruled** — the
>   ACBN0 machinery stays as the on-site-ERI instrument (bare \f$F^0\f$, \f$J\f$, the projected occupations)
>   and as the loop scaffold for whatever U-functional replaces the renormalisation.  No ACBN0 number is
>   quoted as physics.
> - **NEXT = (3) hp.x on MnO with a MATCHED pseudopotential and projector**: write the UPF ourselves
>   (`gth2upf`: the GTH q7/q6 local + separable parts we already carry, PP_CHI = OUR pseudo-atom 4s/3d and
>   2s/2p, PP_RHOATOM from the same atom) so QE's `atomic`/`ortho-atomic` projector IS our χ; validate the
>   UPF on the isolated pseudo-atom's eigenvalues vs `AtomCalculation`; then `pw.x` AFM-II LDA + `hp.x`
>   (nq 2×2×2) → U(Mn 3d), U(O 2p).  `hp.x` is BUILT (`~/Code/q-e/bin/hp.x`, 2026-09-21); QE's own
>   `test-suite/hp_insulator_us_magn/NiO.*` is the deck template (the same AFM-II cell).

> **▶ INCREMENT 3 — SLICE C LANDED 2026-09-21 (`10a4c83b`): THE CONTRACTED MANIFOLD, and the first
> screened U.**  `HubbardManifold::radial` (one coefficient per \f$l\f$-shell on the site; empty = CP2K's
> every-shell convention) + `atomicRadial`: the facade runs the GTH pseudo-atom (LDA, unpolarized) in EXACTLY
> the site's shells (valgen's recipe) and takes its lowest occupied \f$l\f$ orbital as the contraction, with a
> normalisation check across the two codes' radial conventions (`[+U radial]` line; `AoShell` now carries its
> radial so a shell can be recognised).  `LowdinProjector`: a contracted manifold's \f$\chi_m=\phi[:,c]\tilde V\f$
> are S-orthonormalised within the manifold and projected ATOMICALLY, \f$T=S[:,c]\tilde V\f$
> (\f$T^\dagger c=\langle\chi|\psi\rangle\f$, QE's `atomic`); a column manifold stays Löwdin, bit-identical.
> `HubbardU_Atomic(site,l,U)`; `gpwprobe mno MNO_U_RADIAL=atomic`; gates `LowdinProjector.AContractedManifold
> ProjectsAtomically`, `GPW_Si.Γ_U_Atomic3p_Imp_Pol_eqUnpol`.  ctest 906/906.
> - ★ **Two more things the paper had to teach on the way.**  (a) Contracting the LÖWDIN-orthogonalised d AOs
>   with the raw-frame radial (\f$S^{1/2}[:,c]\tilde V\f$) is a different function on a strongly overlapping
>   7-exponent span — a 3d charge of 0.45 on MnO; the frame-independent object is \f$\langle\chi|\psi\rangle\f$.
>   (b) **eq 10c carries no \f$\bar N\f$**: the pair-count denominators use the UNRENORMALISED populations while
>   the numerator carries \f$\bar P\f$ twice — that asymmetry IS the screening (\f$\bar U\propto\bar N^2\f$); my
>   first version weighted both and cancelled it, which is why slices A+B read "near-bare".
> - ★★ **THE MEASUREMENT (MnO AFM-II, LDA, Γ, U=0, deck-shaped recipe, imposed):**
>
>   | manifold (projector) | \f$\bar U\f$ | \f$\bar J\f$ | \f$U_{\rm eff}\f$ | bare \f$\bar U/\bar J\f$ | \f$N_{d\uparrow}\f$ (renorm.) |
>   |---|---|---|---|---|---|
>   | pseudo-atom 3d (atomic) | 15.3 | 4.4 | **10.9 eV** | 28.8 / 7.5 | 4.83 (3.62) |
>   | every d shell, 35 fn (Löwdin) | 11.5 | 3.3 | 8.1 eV | 19.5 / 5.3 | 5.02 (3.93) |
>
>   Screening ≈ \f$\bar N^2\approx0.55\f$ on both.  The paper's Mn value is 4.67 eV (PBE, Mulliken, PAO-3G, dense
>   k, self-consistent) — a factor ~2 below ours, which is what `hp.x` (item 4) is for.  ⚠ \f$\bar J\approx7.5\f$ eV
>   bare is NOT Hund's J (~1 eV): eq 13's numerator includes the \f$m_1=m_2=m_3=m_4\f$ self-terms, which is why
>   the paper quotes only \f$U_{\rm eff}=\bar U-\bar J\f$ where they largely cancel.  On Si the nearly unbound
>   pseudo-atom 3p (ε = −0.019 Ha, 95 % on α = 0.16) OVER-COUNTS — renormalised charge 2.24 > bare 1.83, the
>   two sites' atomic χ overlap — the known weakness of non-orthogonalised atomic projectors and why hp.x
>   prefers **ortho-atomic**; Mn 3d is compact and unaffected.
> - **k-mesh sensitivity MEASURED (same day):** 2×2×2 (8 k, no symmetry reduction, `|ops|=0`): atomic 3d
>   \f$U_{\rm eff}\f$ = 10.06 / 9.77 eV on the two Mn (\f$\bar U\f$ 13.9/13.5, \f$\bar J\f$ 3.85/3.75; bare 28.9/7.45
>   unchanged, as it must be; \f$N_{d\uparrow}\f$ 4.86 → 3.50 renormalised) against 10.9 at Γ — ~10 %, not the
>   factor 2 to the paper; the residual is functional/projector/self-consistency, i.e. the hp.x question.
> - **REMAINDERS / NEXT:** (1) ✅ the k-mesh sensitivity (above); (2) the **outer loop**
>   (\f$U^{(n)}\to U^{(n+1)}\f$, from 0, to \f$10^{-4}\f$ eV) as a facade driver rather than by hand; (3) the
>   **ortho-atomic** projector (Löwdin among the atomic functions of all sites) — a third projector kind; (4)
>   **hp.x on MnO** with a matched projector (item 4): the value oracle.  Until (4), no ACBN0 U is quoted as
>   physics; the +U anchors keep U = 4 eV on the every-shell manifold (CP2K parity).

> **▶ INCREMENT 3 — SLICES A + B LANDED 2026-09-21 (`b580203b`, `090ba17d`, + the AO-basis fix); THE FIRST
> MnO NUMBER IS IN, AND IT DECIDES THE NEXT SLICE.**  What is in the tree: `BasisSet::BareCoulombSource`
> (+`ERI4Block::Transform`) realised in the Gaussian ERI4 mixin over `FourC`, forwarded by `tGPW_IBS`,
> transformed by the spherical view — gate `UTGaussian_BS BareCoulomb.*` (one d shell reproduces
> \f$F^0\f$, \f$F^0+4F^2/49+36F^4/441\f$ and the REAL-basis pair exchange \f$(5/98)(F^2+F^4)\f$ against an
> independent radial quadrature — ⚠ the quoted \f$(F^2+F^4)/14\f$ is Anisimov's complex-basis average, 7/5 of
> it); `Hamiltonian::ACBN0` on the public `HubbardUEstimator` face via `tHamiltonian::MakeHubbardUEstimator`
> (the `SiteMoments` pattern: the term and integrals stay `.Internal.`; the composite finds the term through the
> abstract `HubbardProjection` face); `TOrbital::GetCoeff` + `TOrbitals::GetBasisSet`;
> `SolidCalculation::EstimateHubbardU()`; `gpwprobe mno MNO_ACBN0=1`; gate `GPW_Si.Γ_U_ACBN0_Imp_Pol_eqUnpol`.
> ctest 904/904.
> - ★ **A defect the Si gate could not see and MnO showed at once:** Löwdin-basis coefficients paired with
>   AO-basis integrals gave \f$\bar U=182\f$ eV.  The paper's \f$\bar P\f$ (eq 9) is the AO-basis density
>   matrix on the manifold's functions; Löwdin enters ONLY in the charges (per-orbital \f$\bar N_i\f$ over the
>   same-(species,l) set; per-function populations in the pair-count denominators).  And Si at Γ is NO test of
>   the renormalisation: Γ₁ is s-only and Γ₂₅′ p-only by symmetry, so every occupied orbital's p charge is 0 or 1.
> - ★★ **THE MEASUREMENT (MnO AFM-II, VA span, 7-shell d manifold, deck-shaped recipe, imposed):**
>   at U=0: \f$\bar U=20.0\f$, \f$\bar J=5.5\f$, \f$U_{\rm eff}=14.5\f$ eV — bare (unrenormalised) 19.5 / 5.3;
>   at U=4 eV: 20.4 / 5.3 / **15.1** (the outer loop drifts UP as +U localises d).  Renormalised d charge
>   3.93↑/0.28↓ per Mn.  **The renormalisation barely bites, and that is structural**: on a COMPLETE 7-shell
>   span every d-like KS state is ~93 % inside the manifold (\f$\bar N_i\approx1\f$), and for a localised
>   \f$d^5\f$ shell eq 12 collapses to \f$\approx\tfrac{25}{20}F^0\f$ — the bare shell average.  ACBN0's
>   screening IS \f$\bar N_i<1\f$: a manifold the KS states do not fully live in (the paper's minimal PAO-3G,
>   4.67 eV for Mn; PBE, Mulliken).  So gate (4) of the spec — the CONTRACTED single-3d manifold — is not a
>   sensitivity check but the decisive measurement, and it is the same object increment 1 named as the
>   physically meaningful +U manifold.
> - **NEXT SLICE (C): a CONTRACTED manifold.**  `HubbardManifold` gains an optional radial contraction (one
>   coefficient per shell of the site's \f$l\f$ shells — the atom's own 3d from `AtomCalculation`/the SAD
>   machinery, or a caller's vector); the projector generalises from a column selector to
>   \f$T=S^{1/2}V\f$ with \f$V^\dagger SV=I\f$ (the χ's S-orthonormal), `ManifoldSymmetry` sees 5 functions,
>   `BareCoulomb` transforms the shell integrals through \f$V\f$ on all four indices (`ERI4Block::Transform`
>   already exists).  Then: ACBN0 on the contracted manifold vs the 7-shell one (the sensitivity the user
>   asked for, measured), the outer loop, and `hp.x` as the value oracle (item 4).  A U from the 7-shell
>   manifold is NOT to be quoted as physics.

> **▶ INCREMENT 3 — ACBN0 (Ū, J̄) FROM OUR OWN ON-SITE ERIs — SPEC 2026-09-21 (written after reading
> `~/Code/1406.3259v3.pdf` eqs 1–13 and the tree; user: LAPACK is fair game for any library).**
> - **The formula (paper eqs 8–13, spin-unrestricted, ONE manifold M of functions {m}; write \f$m_1..m_4\f$):**
>   \f$\bar U=\dfrac{\sum_{m_1m_2m_3m_4}\sum_{\sigma\sigma'}\bar P^\sigma_{m_1m_2}\bar P^{\sigma'}_{m_3m_4}(m_1m_2|m_3m_4)}{\sum_{m\ne m'}N^\alpha_mN^\alpha_{m'}+\sum_{mm'}N^\alpha_mN^\beta_{m'}+\sum_{mm'}N^\beta_mN^\alpha_{m'}+\sum_{m\ne m'}N^\beta_mN^\beta_{m'}}\f$,
>   \f$\bar J=\dfrac{\sum_{m_1m_2m_3m_4}\sum_\sigma\bar P^\sigma_{m_1m_2}\bar P^\sigma_{m_3m_4}(m_1m_4|m_3m_2)}{\sum_{m\ne m'}N^\alpha_mN^\alpha_{m'}+\sum_{m\ne m'}N^\beta_mN^\beta_{m'}}\f$,
>   \f$U_{\rm eff}=\bar U-\bar J\f$ (Dudarev), with \f$(m_1m_2|m_3m_4)=\int\phi_{m_1}\phi_{m_2}\,r_{12}^{-1}\,\phi_{m_3}\phi_{m_4}\f$ the
>   BARE integrals over the manifold's functions in the CENTRAL CELL ONLY (paper: \f$g=l=m=0\f$ — no lattice sum,
>   no screening; the renormalisation IS the screening), \f$N^\sigma_m=\bar P^\sigma_{mm}\f$.
> - **The renormalised density matrix, LÖWDIN not Mulliken (user, item 3 above):** per k and spin, each occupied
>   orbital \f$i\f$ has Löwdin coefficients \f$\ell_i=T^\dagger c_i\f$ in the manifold (\f$T=S^{1/2}[:,M]\f$, the
>   projector `Hubbard_U` already owns), and its RENORMALISED occupation \f$\bar N_i=\sum_{M'\sim M}\|T_{M'}^\dagger c_i\|^2\f$
>   — its Löwdin charge in EVERY manifold of the same (species, l) in the cell (the paper's \f$\{\bar m\}\f$: both Mn
>   on AFM-II MnO); then \f$\bar P^\sigma_M=\sum_kw_k\sum_if_{ki}\bar N_{ki}\,\ell_{ki}\ell_{ki}^\dagger\f$.  \f$\bar N_i=1\f$
>   for every orbital gives back \f$n\f$, the +U occupation matrix — the first gate.  This needs the ORBITALS,
>   not D: `TOrbital<T>::GetCoeff()` (the AO-basis coefficients beside the existing `GetCoeffPrime`).
> - **The integrals — a NEW INTEGRAL TYPE, so a new basis face (the pseudo-wall pin allows exactly this):**
>   `BasisSet::BareCoulombSource::BareCoulomb(cols)` → `ERI4Block` (the \f$m^4\f$ tensor over a chosen function
>   subset).  Realised ONCE in the `Gaussian::Orbital_ERI4_IBS<E>` mixin over the evaluator's `FourC` (PG_Cart and
>   PG_Spherical get it together; libcint's matrix-delivery engine throws "not this increment"); FORWARDED by
>   `tGPW_IBS` to its molecular block (a Bloch sum of an AO is the AO in the central cell); TRANSFORMED by the
>   spherical view through its own \f$T\f$ (cart→sphere on all four indices).  MnO VA: 42 Cartesian d functions
>   on one Mn → \f$42^4/8\f$ M&D integrals, seconds, geometry-fixed.
> - **Where ACBN0 lives:** qcWaveFunction sits ABOVE qcHamiltonian, so the estimator cannot take a `WaveFunction`.
>   `Hamiltonian::ACBN0` takes plain orbital data — per (block, spin): weight \f$w_k\f$, coefficient columns,
>   occupations — and asks `Hubbard_U` for `LowdinCoefficients(block, C)` (what the client CONSUMES: the Löwdin
>   coefficients of given orbitals in each manifold).  The FACADE extracts that from `GetOrbitals(irrep)` after
>   Converge and runs the OUTER LOOP the paper uses: SCF at \f$U^{(n)}\f$ → \f$(\bar U,\bar J)^{(n+1)}\f$ → repeat
>   until \f$|\Delta U|<10^{-4}\f$ eV, from \f$U^{(0)}=0\f$ ("the true variational solution" is their future work
>   too).  Per-slot \f$\bar U_k\f$ (the same sums restricted to a slot's functions) is REPORTED as a diagnostic;
>   wiring it into `Uirrep` is one flag once the shell-averaged value is trusted.
> - **Gates:** (1) `UTGaussian_BS`/`UTHamiltonian`: ONE normalised d Gaussian shell — \f$\bar U_{\rm bare}=F^0\f$ and
>   \f$\bar J_{\rm bare}=(F^2+F^4)/14\f$ against Slater integrals from an independent 2-D radial quadrature
>   (the textbook shell averages: the whole chain FourC → c2s → normalisation in one number), rotational
>   invariance, \f$(ii|jj)\ge0\f$, the 8-fold ERI symmetry; (2) \f$\bar N_i\equiv1\f$ reproduces `Hubbard_U`'s
>   \f$n\f$ bit-for-bit; (3) MnO AFM-II: \f$U^{(n)}\f$ converges, printed per outer iteration with
>   \f$(\bar U,\bar J,N_m)\f$; the paper's PBE/Mulliken/PAO-3G value is 4.67 eV (Mn) — ours is LDA/Löwdin/7-shell,
>   so agreement is NOT the claim; (4) **basis sensitivity measured explicitly** (user): \f$U\f$ on the VA span vs a
>   contracted single-3d manifold, before any U is trusted.  Oracle for the VALUE = item 4, `hp.x`.
> - ⚠ **A physically meaningful manifold is ONE radial d function** (increment 1's finding: CP2K's all-shells
>   manifold is a mechanism number).  ACBN0 over 35 functions gives an average over diffuse shells too; the
>   contracted-3d manifold (the atom's own 3d, `AtomCalculation`) is the comparison arm of gate (4), and the
>   likely production manifold.  Decide from the measurement, not now.


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
| **SCF checkpoint/restart (wavefunction + ρ to disk)** | crash-resume (the Oct 6–20 unattended window), a LIBRARY of converged states per material (user 2026-09-27: *"restart any material from a good state"*), the self-consistent-U outer loop across PROCESSES; AIMD's warm start later | **CK-1 ✅ DONE 2026-09-28** (design + execution record → `doc/Records/OpenWork_History4.md` §"CK-1"): `SolidCalculation::SaveState(path)`; `SolidCalcOptions::saveStateTo` (auto-save after EVERY `Converge` — each anneal stage, each ACBN0 outer step — converged or not; written aside and renamed into place); `SolidCalculation::Restart(path, lat, mol, opts, params)` → `Outcome<unique_ptr<SolidCalculation>, RestartRefusal>`, EXACT RESUME vs WARM START with every difference named.  Format `qchem-solid-state 1`: serial HDF5 through `qchem.HDF5` (src/Common), layout + fingerprint classes in `src/Calculation/SolidState.C`.  Gates: UTCalculation `SolidState.*` (7; Si 2×2×2: exact resume ΔE ≤ 2e-15 in 3 iterations Unpol-real and Pol-complex; grid-continuation warm start's FIRST iterate 2e-10 Ha from the converged answer; real↔complex TRIM blocks both ways; refusals on k-mesh / spin group / missing file), UTCommon `HDF5.*` (3) | **CK-2**: a WaveFunction read from disk (`tWaveFunction` without `tSCFWaveFunction` — anticipated by `SCFWaveFunction.C`'s ISP comment) so χ₀/ACBN0/gaps run on a stored state with NO SCF; every CK-1 file already carries C/ε/f for ALL orbitals.  **CK-3**: warm start onto a DIFFERENT k-mesh via the real-space \f$D(\mathbf R)=\sum_kw_kD_ke^{i\mathbf k\cdot\mathbf R}\f$.  CK-1 residuals, each its own small decision: (a) an exact resume still pays ≥ 3 Fock builds, and a NEAR-converged start CRAWLS — ρ_mix halves to 0.5/0.4 and DIIS resets on 1e-12 noise (the Ecut=30 warm start: 9 iterations from 2e-10 Ha) — an accelerator heuristic, not the restart; CK-2 removes the need for the exact-resume SCF; (b) the +U occupation ECHO (n recomputed from the restored D, checked against a saved copy) was NOT built: n is a function of D and the Hamiltonian recomputes it — build it the first time a state disagrees, or with frozen-occupation mode (where n IS independent state and must be saved); (c) no `Restart` overload for an annealed schedule and no mid-stage checkpoint (every N iterations) — add when a run needs one; (d) `h5py` is absent from the system Python, so the Python-side read is untested (the PySCF venv can take it).  States live OUTSIDE git: `~/Code/qchem6-runs/states/<material>/`.  QE `.save` samples: `IntegrationTests/QE/checkpoints/` |
| **DFT+U+J** (Hund's unlike-spin term) | the fractional-SPIN error — magnetic coupling in open-shell TMOs; +U alone leaves it (Macke outlook) | NOT STARTED; ACBN0 delivers J beside U from the same on-site ERIs (§1 step 3) | same term seam as +U, one more scalar per manifold; land after the +U anchor |
| **DFT+U+V** (intersite Hubbard V) | hybridised / charge-transfer insulators where an on-site term cannot restore the bond | NOT STARTED; Macke §5: orbital-resolved U already does most of what +V was added for, so it is a follow-on not a prerequisite | two-centre occupation numbers on the same projector; decide after the O-p manifold question on MnO is answered |
| **Linear response on REAL TRIM blocks** (found 2026-09-28, A7 R2) | every k of a Γ-only or Γ-centred 2×2×2 mesh is TRIM, so with the default LDA stack those blocks run REAL — i.e. almost every run we do (NiO/MnO/SrVO₃ k222, all the A6 materials).  The response path (R1/R2: `tResponse_HT<T>::GetMatrix(const tobs_t<T>*…)`, `OrbitalFrame<T>`, `TransitionBlock<T>`, `AO_TransitionFock<T>`) serves blocks of the RUN's scalar only, so a complex run with real blocks is REFUSED — loudly (`MakeOrbitalFrame` throws), never a wrong number.  Workaround today: the ground state with `forceComplex` — same physics (the 3c-3 gate pins ON == OFF), but complex arithmetic and memory on every block.  ⚠ R0's χ₀ is NOT affected (its Reference is scalar-agnostic): only the self-consistent χ (and so U₀) needs the workaround | NOT STARTED.  Blocks the U₀-vs-hp.x series only by COST, not correctness | two routes — (a) **real-block siblings now**: a `ResponseRealBlock` face per the `Dynamic_HT_RealBlock` idiom (Hartree/XC/+U each get a `GetMatrixR` over their existing scalar-generic `MakeMatrixT` body, which the δ-state already feeds), `FrameBlock` as a per-block scalar variant (as `cd_child_t`), `TransitionBlock` likewise (the composite already takes `tDM_CD<double>` children), `AO_TransitionFock` storing real blocks upcast; ~a day, but it ADDS to the 31 `*R` methods V1.35 exists to delete; (b) **ride V1.35** (`CleanupCandidates.md`: the two-parameter `tDynamic_HT<TBlock,TRun>`), after which the response face is one template with no sibling.  Recommendation: (a) only if the U₀ series is cost-bound on `forceComplex`; measure one material both ways first · `doc/LinearResponsePlan.md` §5d |
| **A7 R3 deferred: efficiency + clean-ups** (logged 2026-09-29 under the user's numbers → efficiency → clean-up order) | R3 lands the RIGHT NUMBERS first; these wait: (a) q-STAR reduction (one q per star, needs the site map for χ_IJ; ≤2× on NiO); (b) the Becke XC gather, the largest bucket per kernel application (~5 of ~11 s on NiO); (c) block-GMRES / warm start across channels and q; (d) H5's full split (ACBN0's `Apply` onto `HubbardUTarget`, R4). Also deferred, with their own rows: the real-pair route ("Linear response on REAL TRIM blocks"), CK-2 ("SCF checkpoint/restart"), retiring R2's periodic forwarding (§4a "Linear response on an IMPOSED run symmetrizes δρ") | NOT STARTED — by design, until R3's NiO and supercell numbers are in | pick up after R3 step 5; measure first (the §3d sizing table is the baseline) · `doc/LinearResponsePlan.md` §3d "Deferred ledger" |
| **LRT and cRPA as first-class Hubbard-U estimation strategies** | today ACBN0 (on-site ERIs + a screening functional) is the ONLY in-tree route; hp.x and ABINIT's `ucrpa` are external oracles we run and compare against by hand.  User (2026-09-25, reading `doc/HubbardUPlan.md` §7's Carta et al. arXiv:2505.03698): no problem running every material in these papers and supporting BOTH methods natively | NOT STARTED; the abstract seam already exists (`HubbardProjection`/`HubbardUEstimator`, increment 3 — `ACBN0` is one concrete strategy behind it, doc/Pins.md pin 23's 2026-09-25 addendum) | LRT needs a perturb-and-respond capability (dn/dα on the manifold, doc/HubbardUPlan.md §5's DFPT-adjacent open question) we do not have yet; cRPA needs χ0 in a product basis (screening's own §2 route (d)) — both large, so this is a DIP placeholder not a near-term build.  Land ONE new concrete strategy behind the existing face when either becomes buildable, never a special case |
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
| **Linear response on an IMPOSED run symmetrizes δρ** (found 2026-09-29 by reading, A7 R3 design note; NOT run) | A Hubbard probe (one site) or a q ≠ 0 wave BREAKS the imposed group, but two density paths still symmetrize what they are handed: (1) the T3 stream fold, armed on imposed Γ-only runs on the SHARED molecular block — `CollocateDensity`'s `FoldProjectedD` projects δD onto the group, and `IntegratePotential` fills partner elements assuming a symmetric δV; (2) the raster star-average in `Composite_Fourier::GetRhoOnGrid` (the fit basis's `G_RasterTransform::Symmetrize`), which the pair sampler's `Sample(δ, σ)` goes through on the uniform-XC route.  The Becke `FoldedMesh` path is correctly skipped.  Every R2 gate ran FREE, so neither path has been exercised; NiO's R3 recipe (k 2×2×2, Becke) misses both. | R3 step 3 retires the R2 periodic forwarding and routes q = 0 through the new never-folding transition entry points; until then, `HubbardLinearResponse` refuses (Outcome) an imposed run on a raster XC route or with a Γ-folded collocation.  A gate to prove it first: the free-vs-imposed χ on Si Γ (they must agree) · `doc/LinearResponsePlan.md` §3d finding 5 |
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
