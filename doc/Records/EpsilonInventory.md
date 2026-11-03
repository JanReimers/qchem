# ε inventory — every numerical approximation a run makes, against the one accuracy target

RECORD, step 0 of the OpenWork §3 row "ε must govern EVERY numerical approximation" (pin 31, extended 2026-11-03).
Taken by lane B on 2026-11-03 from `main` f400f63d. For the pilot and `vlocEps`, the source is **lane A's uncommitted
tree** (`~/Code/qchem-lane-a`, the density-ball ladder and the derived `vlocEps`), which was about to land. Re-check every
line number against the tree before acting on it.

**Classes** (pin 31):
- **MODEL**: a choice of model, exempt. Covers the functional, orbital basis, PP, smearing kT, structure, and the +U projector flavour.
- **ε-DERIVED**: sized from ε by a formula. "R/C separated" means the rigorous bound and the named calibration live apart
  (CLAUDE.md Design, `GridAccuracy::Rigorous` / `Calibrated`).
- **ε-AUDITED**: the run measures it against ε (which pilot check is named).
- **STATED**: the user states it. It is honoured and reported against ε.
- **FIXED**: a constant that ignores ε. This is a **defect** under pin 31, however small it measures today.
- **n/a**: not an approximation of the converged answer: an algorithm control (the converged answer does not depend on it),
  a roundoff identity, or a reporting threshold. Listed so the survey can be checked as complete.
- ⚑ **USER'S CALL**: the class itself is a ruling the user has not yet made.

ε defaults to 1e-5: Ha/atom on an energy run, Ha per matrix element on a geometry run (`GridAccuracy::eps`,
`GPWTolerances.C:134` lane A). `accuracy 0` means no target: the rasters are as stated.

## 0. Scope by facade: ε is one concept for atom, molecule and solid (pin 31, 7ecf191f)

**ε exists today on the SOLID facade ONLY**: `SolidCalcOptions::accuracy`, a `GridAccuracy` (deck `solid.accuracy`;
OOD-SOLID-Cleanup V-FACADE-ACCURACY). `MoleculeCalculation` and `AtomCalculation` carry no ε. Every tolerance they use is
FIXED by construction, so no row in the M or A column can be ε-derived or ε-audited until V-FACADE-ACCURACY lands.

⚠ **`AtomCalcOptions::accuracy` is NOT ε.** It is a `BasisSetAccuracy` {Low, Medium, High} exponent-pool PRESET
(`AtomCalculation.C:50`): it chooses the orbital basis, which is a MODEL choice and exempt. It shares the name only, and
V-FACADE-ACCURACY renames it.

A = AtomCalculation, M = MoleculeCalculation, S = SolidCalculation (GPW unless marked PW). **ε?** = the facades where the
tolerance consumes ε today, by derivation or audit.

| tolerance (row below) | A | M | S | ε? |
|---|:-:|:-:|:-:|---|
| raster demands `rasterX` / `vlocRasterX` / `xcRasterX` (§1) | | | ● | S |
| `vlocEps` (§1) | | | ● | S (lane A) |
| density G-ball `densityEcut` / `cutoffFactor` (§1) | | | ● | S (lane A ladder) |
| routing `fieldSharp` / `relFieldSharp` / `ladderFactor` / κ (§1) | | | ● | — |
| `screenEps` | | ● (`ToleranceBinding`) | ● | — |
| `densityEps` | | | ● | — |
| `beckeEps` (periodic partition images) | | | ● | — |
| periodic Becke mesh (§2) | | | ● | S on geometry runs only |
| molecular XC mesh (80, L29) (§2) | | ● | | — |
| atom XC mesh (50, 1 dir) (§2) | ● | | | — |
| PW XC ball 4× (§2) | | | ● PW | — |
| Vxc fit basis: Delta default, A1 opt-in (§2) | ● | ● | ● | n/a at Delta |
| Coulomb fit set, FittedVee (§3) | ● (radial fit set) | ● | | — |
| PGData, Hermite and ERI4 floors (§3) | (radial bases: own integrals, not surveyed) | ● | ● (via PG_Cart_MnD) | — |
| SCF gates `MinΔρPerElectron` / `MinFD` (§4) ⚑ | ● | ● | ● | — |
| `MergeTol` (§4) | ● | ● | ● | — |
| `orthoTol` / vet trim (§4) | | | ● | MODEL |
| LASolver fallback floor (§4) | ● | ● | ● | — |
| k-mesh ⚑, Ewald (§5) | | | ● | — |
| +U `tolU_eV`, LR Krylov `tol` (§5) | | | ● | — |
| `BasisSetAccuracy` pool preset | ● | | | MODEL (not ε) |

The A and M columns are where the FIXED list (§7) grows once ε reaches those facades. The molecular XC mesh (rank 2) and
the Coulomb fit set (rank 4) have no ε to judge against today, not just no audit.

## 1. The GPW solid path

| name · where set | default | approximates | class | evidence of size today |
|---|---|---|---|---|
| `GridAccuracy::eps` · `GPWTolerances.C:134`, deck `solid.accuracy` | 1e-5 | — (the target itself) | — | — |
| `rasterX` (density ladder: Hartree, the uniform XC's input) · `GPWTolerances.C:63` | 0 → pilot | alias of the collocated ρ̃ on the raster | **ε-DERIVED**: rigorous `kDensityTail` (S 0.35, n 3) + self-calibrated η per run (`MeasuredDemand`); R/C separated ✓ | NaF BallOnly −941 μHa was the XC raster (pin 31) |
| `vlocRasterX` · `:64` | 0 → pilot | alias of V_loc-long's integrate-back | **ε-DERIVED**: `kVlocLongTail` (S 2.35, n 3) + pilot η; R/C separated ✓ | MnO V_long 0.26 Ha in H vs 5 μHa in E (pin 31) |
| `xcRasterX` (uniform XC on L0) · `:65` | 0 → pilot | v_xc[ρ] sampled on the raster | **ε-AUDITED**: pilot `MeasuredPowerDemand` (no bound; labelled MEASURED) | Si max\|ΔH_xc\| 1.0e-2 → 4.8e-4 at x 21 → 54 |
| `vlocEps` (V_loc-long per-pair routing) · `:29`; derived in `SolidSetup.C:236-241` (lane A) | 0 → ε / η_routing = ε/10 | each (pair, offset)'s spectral tail on its multigrid level | **ε-DERIVED** (lane A): `RigorousVlocEps()` = ε per pair, calibration `vlocRoutingEta` 10 (measured 1.0–6); R/C separated ✓; stated > 0 wins | fixed 1e-5 scattered MnO ±80 μHa across densityEcut 250/320/400; at 1e-7 it was 0.75 μHa (lane B + lane A, 2026-11-02) |
| `kVlocEpsNoTarget` · `:20` (lane A) | 1e-5 | the same, on an `accuracy 0` run | **STATED** (no target is a stated grid); reported "(no-target fallback)" | as above |
| density G-ball: `densityEcut` / `cutoffFactor` · `SolidSetup.C:100-101` | auto = 2·α_max (ball 2× that) | the density's own Fourier ball (the tail of the tightest pair, e^{−4} = 1.8 % at auto) | **ε-AUDITED** (lane A): nested ladder, run vs twice the ball at fixed D, judged on Hartree; auto doubles, stated is honoured and reported OVER | H₂O DZVP-MOLOPT +330 μHa at auto vs 2.5 at 100 Ha; NH₃ 66; CH₄ 11; MnO auto 169 μHa below the 1e-7 reference |
| `localPPRelCutoff` κ · `:39` | 30 Ha (tail e^{−15} = 3e-7) | the local-PP sweeps' absolute pair→level rule | **FIXED** | "60 is the self-convergence check" (comment); no recorded number |
| `fieldSharp` · `:51` | 2/3·α_max | the KS-field core exponent: which pairs route to which raster | **FIXED** (the raster audit holds it on BOTH sides) | none; the `vlocEps` precedent shows routing error hides from raster audits |
| `relFieldSharp` · `:42` | 1/3·α_max | the HartreeOnly routing floor | **FIXED** (HartreeOnly route only) | β = 0 diverges (+904 Ha), comment |
| `ladderFactor` · `SolidSetup.C:102` | 4 | multigrid gradation (which level a pair lands on) | **FIXED** (the audits keep it on both sides) | none |
| `relCutoff` · `:54` | 0 (the evaluator's own rule) | CP2K REL_CUTOFF parity | **STATED** (>0 only) | — |
| `mgridEcuts` · `:46` | empty (auto ladder) | an explicit CP2K CUTOFF/3^i ladder | **STATED** | — |
| `screenEps` · `:49` | 1e-10 | magnitude screen of the analytic 1E lattice sums (S, ⟨p²⟩, V_loc) and the collocation boxes' reach | **FIXED** | "harmless at ε 1e-5 by every measurement so far" (OpenWork row); unmeasured at 1e-7 |
| `densityEps` · `:57` | 1e-10 | the collocation task-list floor (the lowest any screener answers) | **FIXED** | as `screenEps` |
| pilot reference margin · `SolidCalculation.C:448` | x_ref = x_run + 8 | the audit's own reference raster | **n/a** (the audit's instrument). Its adequacy bounds every ε-AUDITED row above | none |
| pilot round cap · `:449` | 4 | — | **n/a** (cost cap; at the cap the run reports NOT SETTLED, lane A) | — |
| `beckeEps` · `Mesh.C:174` | 1e-6 | the periodic Becke partition's image (competitor) series | **FIXED** | weights move ~1e-6 relative; 1e-5 breaks the site-share gate (comment) |

## 2. XC quadrature (solid, molecule, atom)

| name · where set | default | approximates | class | evidence |
|---|---|---|---|---|
| periodic Becke mesh, **geometry** run · `BeckeXCParams()` `XCPolicy.C:273-287` | (40, MHL 2, Lebedev 29) | E_xc / V_xc quadrature | **ε-AUDITED**: the pilot's nested rung ladder on max\|ΔH_xc\| (`BeckeLadderAbove`; stated mesh honoured, reported OVER) | NaF/LiF (20, 2, 24) 4–10× wrong on a 0.1 bohr curvature |
| periodic Becke mesh, **energy** run · same | same | same | **FIXED** (known miss 1: the rung audit runs only `if (onH && …Becke)`, `SolidCalculation.C:1459` lane A) | CH₄ DZVP-MOLOPT: default 52 μHa short against a 5-atom target of 50 μHa; the MOLOPT probes pass only because they state (60, GL 48) |
| Uniform-vs-Becke selector margin · `XCPolicy.C:187` | `kUniformMargin` 2 | picks the XC route by cost | **n/a** (picks a route; each route's accuracy is its own row) | — |
| plane-wave path XC ball · `PlaneWave` evaluator (pin 28 floor) | 4 × the wavefunction ball | the XC fit on {G} | **FIXED**. The PW family has no pilot at all (`SolidSetup.C`: the raster block is `!s->planeWave`) | at 1× it was 60 μHa on β-MnO₂ (OpenWork PW row) |
| `planeWaveEcut` · `SolidSetup.C:95` | stated (required) | the orbital basis | **MODEL** (the orbital basis) | — |
| molecular XC mesh · `MoleculeSetup.C:63-66` | (MHL 80, Lebedev 29) | E_xc quadrature on a molecule | **FIXED** | measured against 150/35, per atom: water 1 μHa, N₂ 3, **F₂ GTH q7 27 μHa** (2.7× ε 1e-5) |
| atom XC mesh · `AtomCalculation.C:84` | (MHL 50, 1 direction) | radial quadrature (atoms spherical) | **FIXED** | "the proven atom values"; no number at hand |
| Vxc fit basis (solid, molecule, atom) · `VxcFit::Auto` | Delta (v_xc sampled on the quadrature's own points) | — | **n/a** at the default (no fit beyond the grid). `Auxiliary` (A1) is opt-in, and then **FIXED** (fit basis unaudited) | Si pseudo-atom δ vs A1: 5 vs 33 μHa from CP2K ATOM |

## 3. The molecular path's own screens and fits

| name · where set | default | approximates | class | evidence |
|---|---|---|---|---|
| Coulomb fit (FittedVee) auxiliary basis · `Hamiltonian/Factory.C` (the lineage's fitted auxiliary basis) | the orbital basis's own fit set | J by density fitting | **FIXED** (no audit). ⚑ whether a fit set is MODEL or NUMERIC is open | none recorded on the molecular path |
| `screenEps` (molecular basis) · `ToleranceBinding`, `PG_Cart_MnD/Evaluator.C:110-125` | 1e-10 | pair magnitude screen | **FIXED** | as the GPW row |
| PGData pointwise floor · `PGData.C:41`; `DeltaFit_IBS.C:58` | 1e-10 | a primitive's support radius | **FIXED** | none |
| Hermite recursion floor · `Hermite2.C:71` | 1e-13 | Hermite coefficients dropped | **FIXED** (near roundoff) | none |
| ERI4 screen · `ERI4.C:57` | 5e-16 | four-centre integrals dropped | **FIXED** (at roundoff) | none |
| `symmetryTol` · `MoleculeSetup.C:76` | 1e-4 | point-group detection | **MODEL** (the structure) | — |

## 4. SCF, ortho, basis conditioning

| name · where set | default | approximates | class | evidence |
|---|---|---|---|---|
| `MinΔρPerElectron` (+ `Δρmeasure`) · `SCFParams.C:31` | 1e-4, ∫\|Δρ\|/N | an unconverged D, i.e. an energy error | ⚑ **USER'S CALL** (today: STATED per deck, never judged against ε) | the old mixer-residual measure was 50–120× loose for Kerker/Pulay (comment `SCFParams.C:33-36`) |
| `MinFD` · `SCFParams.C:62` | 1e-5 | orbital gradient [F,D] | ⚑ **USER'S CALL** (as above) | E error ~ ‖[F,D]‖² near convergence; unmeasured against ε |
| `MinΔE`, `MinΔFD`, `MinVirial` · `SCFParams.C:47-52` | off (1e30 sentinel, OpenWork row) | — | ⚑ **USER'S CALL** | — |
| `NMaxIter` · `SCFParams.C:27` | 100 | — | **n/a** (a safety cap; a cap hit FAILS the test, `qchem.SCFCap`) | — |
| `MergeTol` · `SCFParams.C:71` | 1e-4 Ha | which levels count as degenerate (occupation sharing) | **FIXED**, or MODEL if it is read as the occupation rule. Flagged | none |
| Ladder / DIIS / GDM controls (`ladderEThresh`, `ladderSwitchAt`, `PCFloor`, `eMax`, `svTol`, …) · `SCFAccelerator/Factory.C:35-37`, `SCFAccelerator.C:170-216`, `MoleculeSetup.C:101-115` | various | the path to the fixed point, not the fixed point | **n/a** (algorithm controls) | — |
| `orthoTol` · `SolidSetup.C:189` | 0 (no drop) | which near-dependent functions the ortho drops, i.e. the span | **MODEL** (the span), set by a number and reported (`[basis trim]` when `basis.vet`) | pin 22 |
| vet-stage trim · `VetStageTrim`, `BasisSet/Lattice/BasisSet.C:175` | off unless `basis.vet` | the same, per element and l | **MODEL**, reported as `[basis trim]` | — |
| LASolver singular floor · `LASolverLapack.C:50` | 1e-8 absolute | the no-clean-gap ortho fallback | **FIXED** (fallback path only) | none |

## 5. Lattice sums, k-mesh, +U

| name · where set | default | approximates | class | evidence |
|---|---|---|---|---|
| **k-mesh** · `Deck.C:146` | {1,1,1} (Γ) | the BZ integral | ⚑ **USER'S CALL**: a BZ quadrature is a numerical approximation of the crystal model. Today STATED, never judged against ε | campaign decks state their own mesh |
| `kShift` | deck | MP offset | follows the k-mesh row | — |
| smearing `kT` | deck | — | **MODEL** (pin 31 lists it) | — |
| Ewald truncation · `Ewald.C:53` | p = 34 (e^{−34} ≈ 2e-15), η = √π/Ω^{1/3} | real- and reciprocal-space series cutoffs | **FIXED**, but trivially derivable (p = ln(1/ε')) | 2e-15: negligible |
| +U projector (pseudo-atomic orbital) · `PseudoAtomOrbital.C` | its pool basis | the projector function | **MODEL** (the projector flavour). The PAO's own solve is 3 / 0.3 / 0.1 mHa vs CP2K ATOM (E, ε_4s, ε_3d; comment) | — |
| +U self-consistent loop `tolU_eV` · `Deck.C:103` | 1e-4 eV | \|ΔU\| per manifold | **FIXED** (a derived U feeds E) | none |
| LR χ Krylov `tol` · `Deck.C:108`, `Krylov.C:58` | 1e-8 (LR), 1e-10 (generic) | relative residual of χ₀/χ | **FIXED** | none |

## 6. Not approximations (checked, excluded)

These are roundoff and identity tolerances, or reports:
- `kZeroTol` 1e-12 (`KPlusG.C:19`)
- `kMeshMatchTol` 1e-8 and direction `dtol` 1e-9 (`SymmetrizeMesh.C:55, 413`)
- the fractional-site `tol` 1e-6 (`BasisSet/Lattice/Imp/BasisSet.C:394`)
- `LittleGroup` `tol` 1e-6
- the occupation and hole report thresholds 1e-6 (`SCFIterator.C:71, 104`)
- `kOrderFloor` 0.05 (the order detector)
- `GaussLegendre` eps 3e-14 (node roots)

## 7. The FIXED rows, ranked by likely impact at ε = 1e-7 Ha/atom

1. **Periodic Becke mesh on an ENERGY run.** It is unaudited. CH₄ measured 52 μHa at the default, about 100× the 5-atom target at 1e-7. Fix: run the rung audit on energy runs too (OpenWork known miss 1).
2. **The molecular XC mesh (80, Lebedev 29).** F₂ is 27 μHa per atom, water 1 μHa per atom: 10–270× over at 1e-7. It needs the same nested-rung audit as the periodic Becke mesh.
3. **The multigrid ROUTING constants: `fieldSharp`, `relFieldSharp`, `ladderFactor`, and κ `localPPRelCutoff`.** Every raster audit holds them on both sides, and `vlocEps` showed that a routing error hides from a raster audit (±80 μHa on MnO). κ = 30 already leaves a 3e-7 relative tail. Measure each against a doubled value at fixed D, as the ball ladder does.
4. **The molecular Coulomb fit basis.** It is unaudited, and there is no number on the molecular path. ⚑ The user must rule whether a fit set is MODEL or NUMERIC.
5. **The plane-wave path's XC ball (4×), and the PW family generally.** It has no pilot and no target, and it moved 60 μHa once.
6. **`beckeEps` 1e-6.** The partition weights move at 1e-6 relative, against E_xc of order 10 Ha per cell. The cancellation in Σw = 1 is unmeasured.
7. **`screenEps` and `densityEps` 1e-10.** Harmless at 1e-5 by every measurement so far. They are the candidates to scale with ε (ε·10⁻³ or a pilot term).
8. **The +U loop `tolU_eV` and the LR Krylov `tol`.** They feed E only through U. They should be derived from ε through dE/dU (the E-vs-U slope).
9. **`MergeTol`.** It decides degeneracy, and so occupations. It is a class question as much as a size question.
10. **Negligible at 1e-7, but FIXED by the letter of pin 31:** Ewald p = 34, the ERI4 screen 5e-16, the Hermite floor 1e-13, PGData 1e-10, and the LASolver fallback floor.

⚑ **Open for the user (classes not decided here):**
- The **k-mesh**: should a BZ quadrature be judged against ε, e.g. a k-ladder audit like the ball's, or stay a stated model input?
- The **SCF gates** (`MinΔρPerElectron`, `MinFD`): should the gate derive from ε, i.e. stop when the estimated E error of an unconverged D is below ε·N?
- **Whether an auxiliary fit set is a model choice.**
