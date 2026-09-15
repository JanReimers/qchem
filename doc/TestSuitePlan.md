# The Test Suite as a PRODUCT SPACE — axes, naming, file breakdown, and what falls out (item **TE**)

*Drafted 2026-09-15 as THE QUEUED PROGRAMME's step 3 (`doc/OpenWork.md`).  Status: PLAN — rulings 1–6 GIVEN
2026-09-15 with amendments, ruling 7 (materials = data at the Calculation level, lattices in qcStructure) the same
day (§11).  **Phase 1 DONE 2026-09-15 (suite wall 207 s → 71 s); phase 2 (the harness) is next.***

## 0. The ask, and what "done" looks like

**USER, 2026-09-08:** *"for the test review I am also concerned about organization. For example for SCF can we
define a number of axes {Basis Set: PW, GPW, LAPW...} × {material: Si, NaF, MnO, Na, Al,...} × {real space
Grids: Uniform, Becke} × {Gamma, multi k} × {Symmetry imposed: yes, no} × {kT: 0, anneal} etc... And then
decide file breakdown {By basis set?} and layout chosen permutations in a consistent order with a consistent
test naming convention."*

⇒ **An SCF test is a POINT in a product space.  Name it as one, file it as one, and the coverage grid — which
cells are tested, which are holes, which are tested twice — becomes something a reader (or a script) can SEE.**
The layout is the deliverable; the deletions, re-enablements and cost wins are consequences.

**Done means:**
1. `doc/` states the axis vocabulary (§1), the claim kinds (§2) and the naming grammar (§3) — and a test name
   that does not parse under the grammar is a naming defect, checkable by a 40-line script (§5).
2. Every solid-state SCF test lives at `IntegrationTests/<Basis>/<Material>.C` under suite `<Basis>_<Material>`,
   ordered inside the file by the axis order (§4).
3. ONE harness: every solid test drives `qchem::SolidCalculation`, as every molecular test already drives
   `qchem::Calculation` (§6).  `RunGPW` / `RunGpw` / `RunGpwAnnealed` are gone.
4. `PolarizedRunKeepsItsSpin` is a mixer unit test in `src/ChargeDensity/tests` and the 217–251 s
   integration run is deleted (§7).  Suite wall time roughly halves.
5. Every `DISABLED_` test has a recorded verdict — promoted, re-enabled, or deleted (§8).
6. The unit-level tests parked in `IntegrationTests/` are re-homed under `src/<lib>/tests` (§9).
7. `doc/TestFacadeMigrationPlan.md` retires to `doc/OldPlans/` — it describes a scaffold that no longer exists.

## 1. The axes — the vocabulary

The user's list is open-ended by design; this is what the suite ACTUALLY varies today (census of every
`GPW_SCF`, `PlaneWaveDFT` and `RealComplexTerms` body, 2026-09-15), with a TOKEN per value.  Tokens are
C++ identifiers, and UTF-8 identifiers are fine — `Γ` is spelled `Γ` (user 2026-09-15; verified through gtest's
macros, `--gtest_list_tests` and `--gtest_filter` with the project's clang 21, and `SymQNTests.Ωκ_SequenceIndex`
already goes through ctest and TestMate today).  Each axis has a DEFAULT that is **elided from the name** — the
default is the shipped `SolidCalcOptions` default, so a name with no tokens at all IS the facade's own recipe.
**Exception: k is ALWAYS named** — it is the first token, so every name starts from a known place and reads
left-to-right; `GPW_Si.Γ_CP2K`, never `GPW_Si.CP2K`.

| # | axis | tokens | default (elided) | notes |
|---|---|---|---|---|
| 1 | **Basis family** | `GPW`, `PW`, `LAPW` | — always named (suite prefix) | the ENGINE, per `doc/BasisSetTaxonomyPlan.md` §1 |
| 2 | **Material** | `Si`, `Al`, `Na`, `NaF`, `MnO`, `CsI`, `Jellium`, `Cosine` | — always named (suite suffix) | a cell + a formula; the CELL is the only material-specific thing (`doc/OpenWork.md` "ONLY THE CELL") |
| 2b | **…in a box** | `SiBox`, `NaBox`, `O2Box`, `MnBox`, `Mn2Box`, `Na2Box` | | `CellImages::HomeCellOnly` is not a knob, it is a different MATERIAL: an atom or molecule in a periodic box whose oracle is the MOLECULAR facade |
| 3 | **k sampling** | `Γ`, `k211`, `k222`, `k222s` (shifted MP), `k311` | `Γ` — but ALWAYS named (see above) | lower-case `k` (user 2026-09-15); `s` = the CP2K shifted-MP convention; `k311` is the mixed real/complex mesh (has a non-TRIM k — `feedback_complex_type_vs_value`) |
| 4 | **XC grid** | `Uni`, `Becke` | `Auto` (the facade's cost selector, V1.26) | the INTEGRATION grid |
| 5 | **v_xc fit basis** | `PWFit`, `DeltaFit` | `Auto` | ORTHOGONAL to 4 — `doc/Pins.md` "everything is a fit" |
| 6 | **Symmetry** | `Imp` (space group), `Shub` (Shubnikov, polarized), `Grey` (negative control) | free | `Imp` on an AFM density is the V1.28 hazard; `Shub` is the only imposition a magnetic run may ask for |
| 7 | **Spin** | `Pol` (explicit two-channel singlet), `M2`/`M3`/`M6` (multiplicity), `ShFermi` (spins share μ) | unpolarized (ζ=0 collapse) | `Pol` ≠ unpolarized: it is the cross-check that the polarized machinery collapses to the unpolarized anchor |
| 8 | **Occupation** | `Smear`, `Anneal` (staged kT), `GlobalMu` (one μ over the k-mesh) | kT=0 aufbau, per-block filling | `Anneal` is a SCHEDULE, not a temperature |
| 9 | **Convergence machinery** | `GDM`, `Ladder`, `Kerker`, `Pulay`, `MOM`, `SeedMOM` | DIIS, linear D-mixing, no MOM | the `doc/SCFStrategyPlan.md` role seams |
| 10 | **Ansatz** | `Cplx` | real TRIM blocks | `forceComplex` — the downgrade direction (`doc/RealComplexPlan.md` §1) |
| 11 | **Seed** | `UniSeed`, `SpinSeed` | `IonicSAD` (V2.2) | `SpinSeed` = the spin-SAD channel seed |
| 12 | **Route** (implementation, not physics) | `Unfolded`, `Singles`/`Pairs`, `DMSource` | the shipped route | appears ONLY inside a twin claim (`eqUnfolded`): a route is never a test's point, only its control arm |

Axes 1–8 are the PHYSICS axes the user named (plus the fit-basis split the pins demand).  Axes 9–12 are
MACHINERY axes; they exist in the suite today and the grid has to show them, but a name carries them only
when the test is ABOUT them.

**Molecular and atomic tests already have a proto-convention** — `A_HF_U.Energy`, `M_DFT.WaterPolarizedLDA`,
`M_Sym.water_HF_spherical_polarized`: prefix `A_`/`M_`/`L_` = system class, then model `HF`/`DFT`/`DHF`, then
`U`/`P` = spin.  Their axes are `{class: A, M} × {radial basis: SL, SG, BS} × {model: HF, DFT, DHF} × {spin: U, P}
× {angular: Cart, Sph} × {engine: MnD, LibCint} × {symmetry}`.  ▶ **They are NOT re-cut in this plan** — they
are cheap (`A_*`+`M_*` = 4.6% of test CPU, `CTestCostData.txt` 2026-09-15), their names already read as a
grid, and the harness collapse the solids need (§6) was done for them in 2026-08
(`IntegrationTests/CMakeLists.txt`: "the QchemTester scaffold is RETIRED").  Two exceptions ride along:
`M_Sym`'s snake_case is regularised, and **S3b** gets its cell — the one EMPTY cell of the molecular grid:
`{symmetry: yes} × {engine: LibCint} × {angular: Spherical}`.  `M_HF_U.WaterLibCintSpherical` (no symmetry) and
`M_Calculation.WaterSymmetryLibCint` (Cartesian) both exist; the SALC extractor for libcint's real-harmonic
ORDERING is what is guarded out (`doc/OldPlans/SphericalSALCPlan.md` S3b — the TE row's item (c)).

## 2. The three CLAIM kinds — what a test asserts, not what it runs

Reading the 35 enabled `GPW_SCF` bodies, every test asserts exactly one of three things, and the KIND is the
first thing a reader wants to know:

| kind | token | the assertion | how the anchor is JUDGED (KP-0's rule) |
|---|---|---|---|
| **Oracle anchor** | `CP2K` | E(point) == a number from a CP2K deck in `IntegrationTests/CP2K/` | a CP2K number is re-run tight-eps + converged density before it is quoted (`CLAUDE.md`); a move is a PHYSICS question |
| **Did-E-move anchor** | `Anchor` | E(point) == OUR banked number | a move is re-judged against an INDEPENDENT route (KP-0 re-pinned -7.45137 → -7.45294 via the supercell ladder), never merely refreshed |
| **Twin** (one-axis invariance) | `eq<Token>` | E(point) == E(point with ONE axis moved to `<Token>`) | no anchor at all; the tolerance is the gate resolution.  `eqFinite` = the MOLECULAR facade is the twin |
| **Property** | a verb | a non-energy claim: `Converges`, `Stalls`, `Schema`, `RunsReal`, `Deterministic`, `KeepsOrder`, `OrderLostThrows`, `MomentRelaxes`, `SeedMirror` | |

★ **The twin naming is what makes coverage MACHINE-READABLE.**  `GPW_Si.k222_Becke_Imp_eqFree` says, without
opening the file: basis GPW, Si, 2×2×2, Becke grid, symmetry imposed, and the claim is that it matches the
same point with symmetry free.  A script can list every edge of the product graph the suite pins.  Today
that test is `BeckeXC_IBZ_SiDiamond` and the claim has to be read out of the `EXPECT_NEAR` message.

⚠ Two twins today assert against the SHARED CONSTANT rather than the twin RUN
(`PolarizedSingletMatchesUnpolarizedSiGamma` and `PolarizedSeedSingletMatchesUnpolarizedSiGamma` both
pin −7.11506 instead of running the unpolarized arm).  That is a legitimate cost trade (the unpolarized arm IS
`Γ_CP2K`, run anyway) and the name still says `eqUnpol`; the body comment must say "asserted via the shared
anchor `Γ_CP2K`".

## 3. The naming grammar

```
TEST(<Basis>_<Material>,  [<k>_][<Grid>_][<Fit>_][<Sym>_][<Spin>_][<Occ>_][<Machinery>_][<Ansatz>_][<Seed>_]<Claim>)
      └── suite ──┘        └───────────── the POINT: axis tokens in AXIS ORDER, defaults elided ────────────┘ └ §2 ┘
```

Rules:
1. **Axis order is fixed** (the table order in §1).  `Becke_k222` is a violation; `k222_Becke` is the spelling.
2. **Defaults are elided — except k.**  `GPW_Si.Γ_CP2K` is Si at Γ on the facade's own recipe.  Because the defaults are
   the FACADE's defaults, "no tokens" always means "what a user gets".  (A test that pins a default explicitly
   for emphasis — `kT0` — does not get a token; the emphasis goes in the body comment.)
3. **One claim per test.**  A body that asserts an anchor AND a property (`NaPseudoAtomInBoxDoublet`:
   `eqFinite` + a did-E-move pin) names the PRIMARY claim and lists the secondary in its first comment line.
   Two primary claims = two tests.
4. **A twin names the MOVED axis by its token**: `eqUni`, `eqFree`, `eqUnpol`, `eqCplx`, `eqPWFit`,
   `eqAufbau`, `eqUnfolded`, `eqFinite`.  The point named on the left is the ARM UNDER TEST; the control arm
   is the default (or the named token).
5. **`DISABLED_` is not a token.**  After §8 there are no disabled tests; an instrument is a `CLIapps/` binary.
6. **`_UT`/`UT` suffixes go.**  The directory says integration; the suite says the point.
7. Unit tests under `src/<lib>/tests` keep their class-named suites (`KerkerMix`, `BeckeMesh`) — the product
   grammar is for SCF tests, whose subject IS a point in the product.  (§9 moves the misfiled ones.)

### 3a. Every enabled solid SCF test, today → proposed

`GPW_SCF` (35 enabled), in the file order they will be RE-FILED in (suite, then axis order):

| today | proposed | kind | note |
|---|---|---|---|
| `SiliconGammaConverges` + `SolidCalculationMatchesTheSiAnchor` | `GPW_Si.Γ_CP2K` | oracle | ★ **the same test twice** once the harness is the facade (§6) — merge; the `ResolvedXCMesh` asserts ride along |
| `GridsReportSchema` | `GPW_Si.Γ_Schema` | property | |
| `RealTRIMBlocksRunRealInReport` | `GPW_Si.Γ_RunsReal` | property | the `Cplx` control arm is inside the body |
| `CrossRunFirstRunAnomalyProbe` | `GPW_Si.Γ_Deterministic` | property | three identical facade runs, bitwise |
| `RealTRIMBlocksMatchComplex_SiGamma` | `GPW_Si.Γ_eqCplx` | twin | |
| `SmearingInertOnGap` | `GPW_Si.Γ_Smear_eqAufbau` | twin | −TS = 0 on a gapped cell |
| `PolarizedSingletMatchesUnpolarizedSiGamma` | `GPW_Si.Γ_Imp_Pol_eqUnpol` | twin | via the shared anchor (§2 ⚠) |
| `PolarizedSeedSingletMatchesUnpolarizedSiGamma` | `GPW_Si.Γ_Imp_Pol_SpinSeed_eqUnpol` | twin | idem |
| `SharedFermiLevelLetsTheMomentRelax` | `GPW_Si.Γ_Imp_M3_ShFermi_Smear_MomentRelaxes` | property | held-vs-shared μ arms inside |
| `StreamFoldImposedGamma_SiDiamond` | `GPW_Si.Γ_Imp_eqUnfolded` | twin | route twin: reduced streams vs full |
| `BeckeXCMatchesUniformXC_SiGamma` | `GPW_Si.Γ_Becke_eqUni` | twin | asserts E_xc, ρ_lost and the v_xc FIELD, not E |
| `DeltaFitUniformGridMatchesPWFit_SiGamma` | `GPW_Si.Γ_Uni_DeltaFit_eqPWFit` | twin | the fit-basis axis alone |
| `PolarizedSingletMatchesUnpolarized_PWFitRaster` | `GPW_Si.Γ_Uni_PWFit_Imp_Pol_eqUnpol` | twin | 1e-6 — the polarized PW v_xc route (V2.3) |
| `SiliconMultiKPlumbing` | `GPW_Si.k211_Anchor` | did-E-move | KP-0's re-judged −7.45294 |
| `SR_2x2x2ShiftedMP_vs_CP2K` | `GPW_Si.k222s_CP2K` | oracle | |
| `SiDiamondIBZ_NonSymmorphic` | `GPW_Si.k222_Imp_CP2K` | oracle | −7.77846 is the Γ-centred CP2K deck |
| `BeckeXC_IBZ_SiDiamond` | `GPW_Si.k222_Becke_Imp_eqFree` | twin | |
| `RealTRIMBlocksMatchComplex_SiMixedMesh` | `GPW_Si.k311_eqCplx` | twin | the only non-TRIM k in the suite |
| `RealTRIMBlocksWithMOMMatchComplex_SiMixedMesh` | `GPW_Si.k311_Uni_MOM_eqCplx` | twin | R2.21 |
| `SiPseudoAtomInBoxMatchesFinite` | `GPW_SiBox.Γ_Uni_eqFinite` | twin | vs the molecular facade |
| `SmearingConvergesDegenerateShell` | `GPW_SiBox.Γ_Imp_Smear_eqFinite` | twin | the internal energy vs finite |
| `StreamFoldOpenShellMatchesUnfolded_SiAtomInBox` | `GPW_SiBox.Γ_Uni_Imp_Anneal_eqUnfolded` | twin | |
| `AlFCCDegenerateShellAufbauStalls` | `GPW_Al.Γ_Stalls` | property | `EXPECT_FALSE(converged)` — the negative control for `Anneal` |
| `AlFCCAnnealedMetal` | `GPW_Al.Γ_Anneal_Anchor` | did-E-move | |
| `AlFCCMetalGlobalMu` | `GPW_Al.k222_Smear_GlobalMu_Anchor` | did-E-move | |
| `AlFCCMetalIBZExact` | `GPW_Al.k222_Imp_Smear_GlobalMu_eqFree` | twin | 1e-4: IBZ == full mesh |
| `NaFCCMetalGlobalMu` | `GPW_Na.k222_Imp_Smear_GlobalMu_Anchor` | did-E-move | |
| `NaPseudoAtomInBoxDoublet` | `GPW_NaBox.Γ_Imp_M2_eqFinite` | twin | secondary: did-E-move −0.141933 (V2.2's basin gate) |
| `O2TripletInBoxMatchesFinite` | `GPW_O2Box.Γ_Imp_M3_eqFinite` | twin | |
| `MnAtomInBoxDChannel` | `GPW_MnBox.Γ_M6_Smear_eqFinite` | twin | secondary: did-E-move −14.6380; `GPW_MN_SPHERICAL` arm stays an env A/B |
| `PolarizedRunKeepsItsSpin` | **→ `KerkerMix.PolarizedDensityMixesPerChannel` (unit, §7); integration test DELETED** | | |
| `ImposedShubnikovHoldsAFMThroughSCF_Mn2Box` | `GPW_Mn2Box.Γ_Becke_Shub_Pol_Smear_KeepsOrder` | property | THE integration-scale collapse detector |
| `ImposedOrderLostIsAPostconditionFailure_Na2Box` | `GPW_Na2Box.Γ_Becke_Shub_Pol_OrderLostThrows` | property | N1/T3 |
| `MnOSeedSublatticesAreEqualAndOpposite` | `GPW_MnO.Γ_Pol_SeedMirror` | property | seed-level, no SCF |
| `MnOSeedVxcMirrorOnBeckeMesh` | `GPW_MnO.Γ_Becke_Pol_SeedVxcMirror` | property | seed-level, no SCF |
| `MnOImposedShubnikovKeepsTheSeedStaggering` | `GPW_MnO.Γ_Shub_Pol_SeedDecoration` | property | seed-level, no SCF |

`PlaneWaveDFT` — the SCF-level ones (the integral-level ones move out, §9):

| today | proposed | kind |
|---|---|---|
| `ScfJelliumUniform` | `PW_Jellium.Γ_Converges` | property |
| `ScfWeakCosineSelfConsistent` | `PW_Cosine.Γ_Converges` | property |
| `ScfSiliconDiamondConverges` | `PW_Si.Γ_Anchor` | did-E-move |
| `ScfSiliconBZSampled` | `PW_Si.k222_Anchor` | did-E-move |
| `FrameworkSiliconGammaMatchesPrototype` | `PW_Si.Γ_eqPrototype` | twin (route: the SCFIterator framework vs the standalone loop) |
| `FrameworkSiliconGammaThroughSCFIterator` | `PW_Si.Γ_Converges` | property |
| `FrameworkSilicon2x2x2ThroughSCFIterator` | `PW_Si.k222_Converges` | property |
| `FrameworkNaFThroughSCFIterator` | `PW_NaF.Γ_Anchor` | did-E-move |
| `FrameworkCsIThroughSCFIterator` | `PW_CsI.Γ_Anchor` | did-E-move |
| `PolarizedSeedAFMStaggering` | `PW_MnO.Γ_Pol_SeedStaggered` (check the cell) | property |
| `ItemK_Explore_ScfDensity`, `ItemK_RelCutoffDensifiesAndConvergesVxc` | campaign names — re-read and re-cut, or delete if the verdict is banked in `doc/OldPlans/PlaneWavePlan*.md` | |

## 4. The file breakdown

**Recommendation: DIRECTORY per basis family, FILE per material, ONE exe.**

```
IntegrationTests/
  gtestmain.C
  GPW/  Harness.C   Si.C  Al.C  Na.C  NaF.C  MnO.C  Boxes.C      # Boxes = SiBox NaBox O2Box MnBox Mn2Box Na2Box
  PW/   Harness.C   Si.C  NaF.C  CsI.C  Model.C                   # Model = Jellium, Cosine
  Molecule/  M_HF_U.C  M_DFT.C  M_Sym.C  M_Calculation.C  M_Umbrella.C
  Atom/      A_HF_U.C  A_HF_P.C  A_HF_dfPin.C  A_DHF.C  A_DFT.C  A_DFT_U.C  A_PP.C  L_PP.C
  CP2K/      (the oracle decks, unchanged)
```

Why this and not the alternatives:
- **One file per basis** (`GPW.C`) — the user's first suggestion — leaves a 4000-line file; the material is
  the second axis in the suite name anyway, so it is the natural file unit.  `ctest -R GPW_` still selects a
  family, `-R _Si\.` a material across families.
- **One exe per basis** (`ITGPW`, `ITPW`) — ctest already load-balances per TEST, so a split buys nothing at
  `-j8`, and every new exe must be added to `allTests` DEPENDS or it silently never runs (`CLAUDE.md`).  ONE
  `ITMain` stays; the subdirectories are source lists.  `.vscode` TestMate needs no change.
- **Suite = `<Basis>_<Material>`**, not `<Basis>` with the material as the first name token: gtest's suite is
  the unit `--gtest_filter` and TestMate group by, and "all of Si" is the question asked far more often than
  "all k222".  The grid script (§5) answers the other cuts.

**Inside a file, tests are laid out in AXIS ORDER** — k first (Γ block, then k211, k222, k222s, k311), grid
inside k, and so on — exactly the order in §3a.  A file opens with its coverage block: the list of its points
and claims, one per line, in that order.  That block is the human-readable grid; §5's script is the check that
it matches the code.

**`Harness.C`** per family is what is LEFT of `GPW_SCF_UT.C`'s first 970 lines after §6: the basis factories
(`MakeBasisSR`, `MakeBasisLowQ`), the production-shaped `SCFParams` (`ProductionGates()`, `TightGates()`) and
the twin helper (`ExpectRealComplexTwins`).  Nothing that runs an SCF — and NO cells: the cells come from
`qchem.Materials` (§6b), because the probe binaries and eventually the GUI need the same ones.

## 5. The coverage grid — and the holes it exposes TODAY

A script `scripts/testgrid` (Python, ~40 lines: `ctest -N` or `ITMain --gtest_list_tests` → parse every
`<Basis>_<Material>.<tokens>_<Claim>` → a table, one row per test, one column per axis, defaults filled in).
It has two jobs: (1) print the grid; (2) exit non-zero on a name that does not parse — the naming convention
becomes a CI check the day the re-file lands.  `ctest`'s `CTestCostData.txt` gives it a cost column for free.

Laying §3a out that way ALREADY shows what the tracker could not (each is a finding, not an action — the
actions are the rulings in §11):

| hole / duplicate | evidence | what it means |
|---|---|---|
| **`GPW × NaF` = ZERO enabled tests** | all 7 NaF tests are `DISABLED_` | the suite's best-validated CP2K comparison (0.10–0.19 mHa, `doc/GPWPlan.md`) has NO standing gate.  The cost was the reason (the full-SR OOM campaign); the SR2 Γ run is the candidate to re-enable under a cost budget |
| **`GPW × MnO` at SCF level = ZERO** | 3 seed-level tests only; the converged AFM-II (run 38) is a hand run | the north-star material has no converged gate.  Long-tagged candidate once cost is judged |
| **machinery axis beyond DIIS = ZERO enabled** | `GDM`, `Ladder`, `Kerker`, `Pulay`, `SeedMOM` occur ONLY in disabled bodies (plus `Kerker` in the test §7 removes) | after §7, Kerker has no INTEGRATION coverage at all.  A cheap `GPW_Si.Γ_Kerker_eqDIIS` fills it — and if it runs the SINGLES route, it also fills the DM-source XC hole (`OpenWork.md` step 3 ⚠, V1.18e) |
| **`Becke × k≠Γ`** = one cell (`k222_Becke_Imp`) | | fine for now; note it |
| **`Anneal`** on Al and SiBox only; `LAPW` has NO SCF test | | the LAPW row of the basis axis is empty by construction (no LAPW SCF yet) — the grid says so explicitly instead of nobody noticing |
| **duplicate**: `SiliconGammaConverges` ≡ `SolidCalculationMatchesTheSiAnchor` after §6 | same cell, same gates, same anchor | merge |
| **near-duplicate**: `PolarizedSingletMatchesUnpolarizedSiGamma` vs `..._PWFitRaster` | differ only in `Uni_PWFit` vs `Auto` | both stay — the fit-basis axis is the point of the second — but the names now SAY that |

## 6. The harness collapse — `RunGPW` / `RunGpw` / `RunGpwAnnealed` → `qchem::SolidCalculation`

Today `GPW_SCF_UT.C` runs SCF through FIVE doors: `RunGPW` (positional, 11 call sites), `RunGpw(GpwOptions)`
(37), `RunGpwAnnealed` (3), the facade (14) — and **`RunMnO`**, a lambda inside `DISABLED_MnO_AFM2_RhombohedralGamma`
(line 4341) that already drives `SolidCalculation` but is a driver all the same: it assembles the cell
(with three env-knob discriminators, `MNO_SWAP_SUBLATTICE` / `MNO_SWAP_ORDER` / `MNO_SHIFT`, plus `MNO_KMESH`),
the recipe and a `GpwReport` bracket, and returns its own result struct (`MnOArm`).  It collapses like the
others: the cell and its decoration move to the materials module (§6b), the recipe becomes the test body's
one `SolidCalcOptions` block, the discriminators become flags on the promoted probe (§8), and `MnOArm` +
`Instrumentation()` go with them — a campaign's instruments belong in the probe binary, not the gate.  `SolidCalcOptions` already carries every `GpwOptions` knob
(the 2026-08-25 harvest: `spinsShareFermi`, `greyImposition`, `momFromSeed`, `siteSpins`, `onIteration`;
`realTRIMBlocks` ≡ `!forceComplex`).  The molecular half did exactly this collapse in 2026-08 and the
scaffold evaporated — that is the precedent, and `doc/TestFacadeMigrationPlan.md` (which planned it) retires.

**Audit before the move** — what a test reaches through the driver that the facade may not expose:
- `GpwResult.E` (`EnergyBreakdown` by term: `E["MinusTS"]`, `E["Kinetic"]`) — facade has `LastIterateTerms()`;
  confirm `Result()->Terms()` or equivalent on the CONVERGED state.
- `GpwReport` / `FpRow` telemetry — the fingerprint observer; `onIteration` covers it.
- `XCProbe` (the v_xc FIELD on the grid, `BeckeXCMatchesUniformXC_SiGamma`) and `RouteProbe` — read the
  converged density through the neutral `ScalarFunction` face the facade exposes (`Result()->Density()`);
  the v_xc field itself may need one accessor.
- `ReportSymmetryFound` — the symmetry banner; `Calculation` has the report cursor.
- Annealing (`RunGpwAnnealed`'s stage schedule `{kT}` × `{accelerator}`) — `SolidCalcOptions` comments say
  "an annealed recipe … stage 0", so a schedule exists; confirm the `{"DIIS","GDM"}` per-stage accelerator.
- `MakeGpwAccelerator("Ladder")` — `SCFAccelerators::Type::Ladder` exists on the facade enum.

Each gap is ADDITIVE on the facade (the molecular rule: fill the facade, never keep the scaffold).  Anchors
stay byte-for-byte; a moved anchor during this step is a defect in the move, not a re-pin.

### 6b. ONE place per material — the lattice/structure set-up (user 2026-09-15)

**Census, 2026-09-15:** the Si diamond cell (`FCCUnitCell(10.26)` + two atoms) is rebuilt by hand 32 times in
`GPW_SCF_UT.C`, 9 in `GPW_UT.C`, 3 in `PlaneWaveDFTUT.C`, 3 in `src/BasisSet/Gaussian/tests/M_PG_BoxWalk.C`, 4
in `src/Structure/tests/SymmetrizeMeshUT.C`, 1 in `src/SCFIterator/tests/SCFTrace.C`; the MnO AFM-II cell
exists in four spellings (`GPW_SCF_UT.C` ×2 incl. `RunMnO`'s env-parameterised one, `GPW_UT.C`,
`PlaneWaveDFTUT.C`, `StructureTests.C`).  Every one of them carries its own lattice constant, its own atom
order and its own species/valence list — which is exactly how "same material, different `a`" or "same
material, Mn added second" drifts into an anchor silently.

**RULING (user, 2026-09-15): the GENERIC and the CONCRETE live at different heights of the DAG.**
- A *lattice type* (FCC, BCC, hexagonal, … — all 14 Bravais lattices) is structure and belongs in
  `qcStructure`, beside `UnitCell` where `FCCUnitCell` already is.  Today FCC is the only one.
- A *material* (Si diamond at a = 10.26 with `{{"Si",4}}`; MnO AFM-II with its decoration) is a USE of a
  structure plus what the pseudopotential says about it — it belongs at the `SolidCalculation` level
  (`qcCalculation`, the top of the DAG), NOT in `qcStructure`.  The GUI will want the same list — "pre-defined
  materials (and molecules) for users to try out" — so it is a LIBRARY concern, not a test helper, and it is
  DATA, not C++: a JSON file, on the pattern the tree already uses twice (`src/Pseudopotential/Data/
  gth_potentials.json` + `PSEUDO_DATA_PATH`; `src/ChargeDensity/Data/atomic_valence_densities.json`).

So the one place is TWO things, each an `OpenWork.md` row of its own (added 2026-09-15: **BL** and **MD**):

**(i) `qcStructure` — the 14 Bravais lattices (row BL).**  `FCCUnitCell(a)` generalises to the family:
`BravaisCell(Bravais::FCC, a)`, `BravaisCell(Bravais::Hexagonal, a, c)`, `BravaisCell(Bravais::Rhombohedral,
a, α)`, … one constructor per lattice system's free parameters, returning a `UnitCell`; the rhombohedral
AFM-II MnO cell that `RunMnO` writes out as a `Matrix3D` by hand is `Rhombohedral` from `Cubic-F` doubled along
[111] — a named construction, not nine numbers.  This is what a test or the data file names when it says
"rocksalt": `FCC + a two-atom basis`.  Scope: the lattice VECTORS; symmetry detection already exists
(`src/Symmetry/Lattice_3D/SpaceGroup.C`) and stays where it is.

**(ii) `qcCalculation` — `src/Calculation/Data/materials.json` (+ `molecules.json`) and a `qchem.Materials`
module that reads them (row MD).**  One entry per material: the Bravais type + parameters, the atom basis
(species, fractional positions, per-site spin decoration), the PP species/valence list — and NOTHING about
the run (no k-mesh, no grid, no kT: those are the test's own axes, `SolidCalcOptions`' job).

```json
{ "Si_diamond": { "lattice": {"type": "FCC", "a": 10.26},
                  "atoms":   [{"Z": 14, "frac": [0,0,0]}, {"Z": 14, "frac": [0.25,0.25,0.25]}],
                  "species": [["Si", 4]] },
  "MnO_AFM2":   { "lattice": {"type": "Rhombohedral", "a": 8.40, "from": "FCC doubled along [111]"},
                  "atoms":   [{"Z": 25, "frac": [0,0,0], "spin": +1}, {"Z": 25, "frac": [0.5,0.5,0.5], "spin": -1},
                              {"Z": 8,  "frac": [0.25,0.25,0.25]},    {"Z": 8,  "frac": [0.75,0.75,0.75]}],
                  "species": [["Mn", 7], ["O", 6]] } }
```
```
namespace qchem::Materials
{
    struct Material { std::shared_ptr<UnitCell> cell; std::vector<std::pair<std::string,int>> species; int Nelec; };
    Material Get(const std::string& name, double aOverride = 0.0);   // "Si_diamond", "MnO_AFM2", ...
    std::vector<std::string> Names();                                // the GUI's pick-list
    Material AtomInBox (int Z, int valence, double a);               // SiBox/NaBox/MnBox: not data, a recipe
    Material DimerInBox(int Z, int valence, double a, double d, const std::vector<int>& siteSpins = {});
}
```
- `Nelec` is DERIVED (Σ valence × count), never stored — a stored count that disagrees with the species list
  is the kind of defect a data file invites.
- `RunMnO`'s geometry discriminators: `siteSpins` is the data file's per-site `spin`; a rigid `shift` and the
  atom-ORDER swap are probe flags applied to the built cell (translate / permute), not material data.
- The k-mesh is NOT in the entry: `Lattice_3D(*m.cell, k)` at the call site; it is axis 3.
- The molecular `MakeWater()` etc. go the same way (`molecules.json`, same module) — the GUI wants both lists.
- **Consumers below `qcCalculation` in the DAG** (`src/Structure/tests`, `src/BasisSet/*/tests`,
  `src/SCFIterator/tests`) cannot import `qchem.Materials`.  They get (i) — `BravaisCell(FCC, 10.26)` plus two
  `AddAtom` lines — which is what they build by hand today minus the lattice-vector arithmetic.  The
  duplication that MATTERS (lattice constant, atom order, species list, decoration) is gone from every place
  that runs an SCF; a unit test of a mesh or a box walk is entitled to its own two-atom cell.

Phase 2 builds (i) for the lattice types the suite actually uses (FCC, simple cubic for the boxes,
rhombohedral for MnO — the other 11 land with row BL on their own schedule) and (ii) for the seven materials
in §1, then repoints every SCF-running call site; `grep -c 'FCCUnitCell cell(' IntegrationTests` → 0.

## 7. `PolarizedRunKeepsItsSpin` → a mixer unit test (EARLY: it is 27% of the suite by itself)

**The claim, read from the body:** "asking for Kerker on a polarized density must not change the physics" —
i.e. `KerkerMixerFactory` handed a polarized seed must compose ONE ρ̃ mixer PER CHANNEL
(`PolarizedDensityMixer`), never a single-map mixer over the spin-blind total.  That is a FACTORY property
(`src/ChargeDensity/Imp/DensityMixer.C::ComposePeriodic`), and the negative control is already a valve there
(`QCHEM_SPINBLIND_KERKER`).  It is not an SCF property; the 12-iteration Mn sextet on a Becke mesh was the
only instrument available on 2026-08-07.

**The unit tests** (in `src/ChargeDensity/tests/KerkerMix.C`, `UTChargeDensity` — an existing exe, no
`allTests` edit):
1. `KerkerMix.PolarizedSeedComposesPerChannel` — a hand-built polarized `FourierDensity` pair with m ≠ 0
   (the `JointPulay.C` `ΔG_Map` helpers build such fields with no SCF); the factory returns a
   `PolarizedDensityMixer` whose channel views are live, and one `Mix` step preserves the INTEGRATED moment
   (the G=0 component of ρ↑−ρ↓) exactly as the per-channel filters predict.
2. `KerkerMix.SpinBlindValveCollapsesTheMoment` — under the valve, the same step returns a leaf and the
   moment is gone.  Documents the negative control so the valve never rots silently.
   (Needs a periodic basis for `CreateVxcFitBasisSet`: the cheapest GPW Si cell at Γ, built, never iterated.)

**DONE 2026-09-15** — three tests, not two: the factory composition, the per-channel STEP (driven through a
20-line `TwoChannelWorking` test double holding two `FourierMixCD` channels, so the composed `Mix` runs with
no SCF), and the valve.  The Mn box is a plane-wave basis at Ecut=4 — no Gaussian basis, no PP, no
Hamiltonian; the seed is the library's Mn Hund pair.  41 ms for the three.

**Then the integration test is DELETED** — not shortened.  What it also claimed ("the S=5/2 moment is there
every iteration", "the energy is the polarized one") is covered by `GPW_MnBox.Γ_M6_Smear_eqFinite` (same cell,
linear mixing) and by `GPW_Mn2Box.Γ_Becke_Shub_Pol_Smear_KeepsOrder` (order sustained THROUGH an SCF, 17 s).
The `(b)` remark in the TE row about the 0.7-bohr point probe is already moot: the probe was replaced by the
integrated site moment on 2026-09-14 (`SolidCalcOptions` note).

## 8. The `DISABLED_` class — every one gets a verdict

23 tree-wide (17 `GPW_SCF`, 4 `GPW`, 1 `Reporting`, 1 `SCFTrace`).  The zero-assert ones went on 2026-09-09;
these all assert.  The rule for the verdict (TE row (d), KP-0): a disabled test is one of —
**(P) an INSTRUMENT** (a ladder/sweep a human runs with env knobs and reads) → promote to a `CLIapps/` probe
binary, out of gtest; **(R) a gate that was parked for COST** → re-enable under a cost budget, `_Long`-tagged
if needed; **(D) a campaign whose verdict is BANKED in `doc/`** → delete; **(F) a real claim that is
currently FAILING or unfinished** → an open tracker row, not a disabled test.

| test | verdict | why |
|---|---|---|
| `GPW_SCF.DISABLED_SiSupercellLadder` | **P** → `CLIapps/gpwprobe --ladder` | `SI_LADDER=n,n,n` env-driven; it is the independent route KP-0 judged an anchor against — keep it RUNNABLE, not in ctest |
| `DISABLED_SR_2x2x2GammaCentred_vs_CP2K` | **R** → `GPW_Si.k222_CP2K` | an oracle anchor (−7.77846, deck `si_fcc_gpw_222_gamma.inp`), parked at "~4 min" in 2026-07; `SiDiamondIBZ_NonSymmorphic` pins the same number IMPOSED in 31 s, so the free arm is the missing `eqFree` twin — time it post-box-walk before ruling |
| `DISABLED_SingleKSweepProbe` | **P** | a k sweep |
| `DISABLED_TermTranslationInvariance` | **R** → `GPW_Si.Γ_TranslationInvariant` (property) | a real invariance (the Rcut>0 KB fix's guard, 1e-3 vs a 1.7 Ha bug); one-electron terms only, so cheap — find out why it was parked |
| `DISABLED_NaFixedDensityTermProbe` | **D** | a term-by-term probe of the Na doublet; V2.2 closed the doublet (`doc/CleanupHistory.md`) |
| `DISABLED_NaFImposedGDMSmearProbe` | **P** or **D** | the NaF imposed+GDM+smear recipe probe; if `GPW_NaF.Γ_CP2K` re-enables (below), this is a variant of it |
| `DISABLED_NaFRocksaltGamma` | **R** → `GPW_NaF.Γ_CP2K` (SR2 basis, `NAF_KMESH=1`) | THE NaF oracle (0.2 mHa vs CP2K).  Parked for cost during the full-SR OOM campaign; the SR2 Γ arm is what `doc/Benchmark.md` times.  Budget it (`memsafe`), and it becomes the NaF row's anchor |
| `DISABLED_NaFGridContinuation` | **P** | the coarse→fine continuation recipe (`GC_*` env knobs) — an instrument |
| `DISABLED_NaFFullBasisRankReduction`, `DISABLED_NaFFullBasisEigenTol` | **D** | `doc/GPWPlan.md` §1: rank reduction DEMOTED to automation, λ~1e-6 runs clean via seeded aufbau — banked |
| `DISABLED_BeckeRecipeLadder_{SiGamma,NaF,MnSextet,AlFCC}` | **P** → `CLIapps/gpwprobe --becke-ladder` | the `BeckeLadder` harness `doc/Benchmark.md` names as the grid-sizing instrument (user 2026-09-06: size the Becke grid) — it is needed, and it is not a test |
| `DISABLED_RotatedLebedevXCProbe_SiGamma` | **D** | R2.15 flipped Lebedev degree-gated ≥29 and landed (`project_concurrent_cleanup_branch`); the rotation probe's verdict is banked |
| `DISABLED_BeckeXCMatchesUniformXC_NaFSR2` | **R** → `GPW_NaF.Γ_Becke_eqUni` | the NaF twin of the Si Becke gate; same cost budget as `GPW_NaF.Γ_CP2K` |
| `DISABLED_MnO_AFM2_RhombohedralGamma` | **F** → tracker | its header still describes the d-channel KB as "WRONG by ~12.6×" — that was EXONERATED on 2026-09-14 (sprint S).  Re-run; if run 38's converged AFM-II reproduces, it becomes `GPW_MnO.Γ_Shub_Pol_Smear_CP2K` (−61.470570), `_Long` if the cost says so; the header is rewritten either way |
| `GPW.DISABLED_IllConditionedChargeProbe` | **D** | an N3/N5-style ill-conditioned pool probe (`feedback_scf_accuracy_levels`) |
| `GPW.DISABLED_DiffuseD{PairVlong,VlongSharpField,KB}Oracle` | **D** or **R** | the diffuse-d oracles from the 2026-08-13 "V_long defect" that was RETRACTED 2026-08-14 (the HGH single-species gotcha); if they pass now they are cheap `GPW.*` basis-level gates, else delete |
| `Reporting.DISABLED_VisualDump` | keep | a renderer PREVIEW, ruled 2026-09-09 |
| `SCFTrace.DISABLED_SolidNonPP_ShouldBringTheVirialColumnBack` | **F** → tracker | a claim about a feature that does not exist yet; a disabled test is the wrong place to remember it |

`CLIapps/gpwprobe` is ONE binary with sub-commands, sharing `IntegrationTests/GPW/Harness.C` (so a probe and a
test build the same cell the same way).  The env-knob style (`SI_LADDER=`, `GC_KERKER_G0=`) becomes flags.

## 9. Unit tests parked in `IntegrationTests/` — re-home them (`feedback_unit_over_integration_tests`)

| file | what it is | home |
|---|---|---|
| `GPW_UT.C` (31 enabled: collocation, lattice sums, stream fold, general-k matrices, KB vs mesh, Shubnikov GMap) | basis-set and term-level unit tests — NO SCF | `src/BasisSet/Gaussian/tests/` (`UTGaussian_BS`) for the basis ones; the two Shubnikov tests to `src/Symmetry/…/tests` |
| `PlaneWaveDFTUT.C`, the non-SCF 18 (`BasisIntegralScalar`, `HartreeSingleCosineMatchesPoisson`, `DensityGridTransformRoundtrip`, …) | PW basis + Hamiltonian-term unit tests | `src/BasisSet/Lattice/tests/PlaneWaveUT.C` exists — extend it |
| `RealComplexTermsUT.C` (10) | Hamiltonian-term real/complex twins | `src/Hamiltonian/tests/` (`UTHamiltonian`) |
| `EigenSolverUT.C` | LA | `src/LASolver/tests/` (`UTLASolver`) |
| `ValenceBasisGen_UT.C` | the generator | its own library's `tests/` |

Existing exes, so no `allTests` edit — but **check ctest's total N goes up by the number moved and down in
`ITMain` by the same number** (the `UTSCFAccelerator` lesson).

## 10. Running order and acceptance

Each phase is one commit series, green at every step (`ctest -j8`, 860 today), and the anchors do not move
in ANY phase — this is a refactor of WHERE and WHAT-NAMED, not of physics.

| phase | does | acceptance |
|---|---|---|
| **0 — paper** | this file; the rulings in §11 | user thumbs-up |
| **1 — the cost win** ✅ 2026-09-15 | §7: THREE `KerkerMix` unit tests (`PolarizedSeedComposesPerChannel`, `PolarizedStepMovesEachChannelAtAlpha`, `SpinBlindValveCollapsesTheChannels` — 41 ms, teeth checked: the first two FAIL under the valve); `PolarizedRunKeepsItsSpin` deleted | **`ctest -j8`: 207 s → 71 s wall** (862 enabled tests, was 860); N: 885 total incl. 23 disabled.  One unrelated timing flake seen under load (`M_PG_BoxWalk.WhereTheContractionSpendsItsTime`, passes alone) |
| **2 — the harness** | §6: facade gaps filled; `BravaisCell` for the three lattice types in use + `materials.json`/`qchem.Materials` for the seven materials (§6b) and every SCF-running call site repointed; every `GPW_SCF` body on `SolidCalculation`; `RunGPW`/`RunGpw`/`RunGpwAnnealed`/`RunMnO` deleted; `TestFacadeMigrationPlan.md` → `OldPlans/` | identical anchors; `GPW_SCF_UT.C` shrinks by ~900 lines; `grep -c 'FCCUnitCell cell(' IntegrationTests` → 0 |
| **3 — the re-file** | §3a renames + §4 directories + per-file coverage blocks + `scripts/testgrid` | `testgrid` parses 100% of `GPW_*`/`PW_*`; ctest N unchanged; the first `Γ_` name is confirmed in `ctest -N` and the TestMate tree before the rest are renamed |
| **4 — the disabled class** | §8 verdicts: `CLIapps/gpwprobe`, re-enables under budget, deletes, two tracker rows | zero `DISABLED_` outside `Reporting.VisualDump`; N goes up by the R's |
| **5 — re-home the unit tests** | §9 | ctest N conserved across exes |
| **6 — fill the first holes** | `GPW_Si.Γ_Kerker_eqDIIS` (Kerker + singles route ⇒ DM-source XC coverage); S3b: `M_Calculation.WaterSymmetryLibCintSpherical` (needs the extractor — real work, not a test) | the two zero-coverage rows in §5 turn green; the molecular grid has no empty cell |

Phase 1 goes first because every later phase's dev loop runs the suite.  Phase 2 before 3 because renaming a
test and rewriting its body in the same diff is unreviewable.  Phase 6 is the first NEW physics coverage and
is last on purpose — the user's ruling: organisation first, cost second, holes are a consequence.

**Out of scope, on purpose:** the molecular/atomic re-cut (§1, cheap and already a grid); the Becke grid
SIZING (`doc/Benchmark.md`, a physics question the promoted ladder serves); DFT+U tests (step 5 of the
programme — they get cells in this grid, `GPW_MnO.Γ_Shub_Pol_Smear_U_CP2K`, when +U exists).

## 11. Rulings — GIVEN 2026-09-15 (user: "I agree with all six rulings"), with amendments

From the same review: **`Γ` is spelled `Γ`**, **k tokens are lower-case** (`k222`) and **k is always named** —
applied above; `RunMnO` is the fifth driver in §6; and the cell set-up is split generic/concrete, §6b (ruling 7).

1. **The token spellings** (§1).  In particular: `Uni`/`Becke`; `Imp`/`Shub`/`Grey`; `Pol` vs `M1` for the
   explicit two-channel singlet; `Smear`/`Anneal`/`GlobalMu`; `k222s` for shifted MP.
2. **Defaults elided = the facade's defaults** (§3 rule 2) — this ties the naming to `SolidCalcOptions`, which
   is the point, but it also means a facade default change RENAMES tests (the IonicSAD flip would have turned
   every `IonicSAD` token into `UniSeed` on the old ones).  Accept that coupling?
3. **Directory per basis, file per material, one exe** (§4) vs one file per basis.
4. **`CLIapps/gpwprobe`** as the home of the promoted instruments (§8), sharing the test harness — vs deleting
   the ladders outright and re-deriving when the Becke sizing question is picked up.
5. **The two cost-gated re-enables** — `GPW_NaF.Γ_CP2K` and `GPW_MnO.Γ_Shub_Pol_Smear_CP2K`: what is the per-test
   budget above which a test is `_Long` (excluded from the default `ctest -j8`, run by `ctest -L long`)?  Today
   the de-facto ceiling is 74 s (`GPW.XCPotentialConsistencyFD`).  Proposal: 60 s CPU.
6. **Phase order** (§10), specifically 2-before-3.
7. **GIVEN 2026-09-15 — and CORRECTED**: the first draft put `qchem.Materials` in `src/Structure/`; the user
   ruled that concrete materials do NOT belong in `qcStructure` — they are a `SolidCalculation`-level concern,
   eventually JSON data the GUI also reads; only the GENERIC (`FCCUnitCell` → all 14 Bravais lattices) belongs
   in `qcStructure`.  §6b now says so; the two halves are `OpenWork.md` rows **BL** and **MD**.
