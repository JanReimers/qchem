# Step 2 — the remaining CleanupCandidates rows, and why each is not obvious

> ⛔ **RETIRED 2026-09-14 — sprint S (V2.2 + V2.5) landed, 860/860, zero re-pins.**  Groups A, C and D are
> closed; E was never work; the two group-B rows that remain (**V1.34**, **R1.0b**) are blocked on things
> outside the sweep and live as rows in `doc/CleanupCandidates.md`.  The programme's step 2 hands off to
> step 3 (**TE**) per `doc/OpenWork.md`.  Kept as a RECORD of how the rows were grouped and what each
> grouping turned out to be worth (group D's finding; group C's "anchor-moving" rows moved nothing).

Cut 2026-09-09 after batches 1–3 of the sweep; done-markers added 2026-09-10 (`851/851`); re-cut 2026-09-14
(`857/857`): groups A and D are CLOSED, C is what is left; **sprint S closed C the same day (`860/860`).**
This is a READING AID for `doc/CleanupCandidates.md`, not a second tracker — every row here is live in
that file and the verdicts belong there.  It exists because "~24 open items" is not a plan, and because
the rows differ less in size than in WHAT IS BLOCKING THEM, which the alphabetical ordering hides.

## ▶ STATE, 2026-09-14 — WHAT CLOSED SINCE THE CUT

| row | verdict | commit |
|---|---|---|
| **R1.0j** | ✅ renamed `XC_Quadrature` → `DensitySampler`; forward AND adjoint migrated to `MatrixForward`/`MatrixAdjoint` on BOTH routes; `Symmetrize` forwarders deleted.  ⏳ the three leftover tenants go with R1.0h | `d180277b` `a7be9b35` `c7272a22` |
| **R1.0e** | ✅ library home SETTLED and executed — `qcChargeDensity`, interface+factory public, concretes `.Internal.`, tests through the factory | `bd11b903` |
| **R1.0h** | ✅ CLOSED 2026-09-14 — slot pre-creation (`56f3db12` `0b753d1e`) + `ChargeBreakdown{lost, siteMoments}` as the site-moment OWNER; the owning SCOPE object DECLINED (no payload once the tenants resolved: `Matrix` stays, `Integrate` stays, `SiteMoments` left).  Found in passing: `charge.lost` was never set on polarized GPW runs | `018aa3e2` |
| **V1.37** | ✅ CLOSED 2026-09-14, all three steps: `SpinGroup`, ONE `tCompositeWF` + ONE `tComposite_CD` with channel VIEWS (steps 1–2); ONE `Vxc`/`FittedVxc`/`Vxc_Quadrature` built FOR the group, spin-native `ExFunctional`, eleven `Ham_*_U/_P` → six, `SymMap` (step 3).  User rulings on the way: clean code over 1e-16 anchors; `Pol` gone, not aliased | `b0310692` `35d811ea` |
| **V2.1** | ✅ CLOSED BY V1.37 step 3 — the "one PolarizedVxc" the user asked for on 2026-08-08 is the one term per operator; ⚠ the row's premise that the scalar ρ cache would DIE was wrong-way: on an SU(2) run the folded doublet's ONE raster is exactly the efficiency the fold buys, so `Rho` stays and the term halves it | `35d811ea` |
| **(the point probe)** | ✅ `SetOrderParameter` / `orderProbe` DELETED; `m_site` = the integrated site moment is the intrinsic order column (user, on `DISABLED_MnO_AFM2_RhombohedralGamma`) | `0563217b` |
| **V1.35a** | ✅ the second slot hook DISSOLVED via `HT_SlotOwner<TRun>` — a diamond done correctly | `1ac9ae50` |
| **V1.36** | ✅ `FittedVcorrPol` memoizes its \f$v_c\f$ fits (two fitters, one per spin) | `3b88f260` |
| **V1.20c** | ✅ `Projector3` promoted out of `Internal` (a written-rule violation, and my `GMap` fix had been half a fix) | `ee50a5a1` |
| **V1.20d** | ✅ the two remaining sites (`DBCacheClient` promoted alone; Slater's re-export → plain import); closing the first EXPOSED a live cross-family reach from `qcCalculation` into the cache mechanism, now `EmitIntegralsCacheReport()`; the audit is a ctest test (`InternalReexportAudit`).  ⚠ Had never been a row — only a paragraph saying "filed as V1.20d" | `0e07933e` |
| **V4.1 / V4.2** | ✅ CHECKED — neither trigger has fired | — |
| **V1.33** | ✅ the taxonomy re-cut, all five plan steps: `qcRadial_BS` / `qcPlaneWave_BS` / `qcGaussian_BS{Point,Lattice}` / `qcLattice_BS`, ctest `BasisSetGTagAudit`, `PG_Cart` untagged as the G=1 seed.  ⚠ pybind breaks (flagged) | `c2cb79a3`..`d5ddb1a5` |
| **V1.17** | ✅ `GetSpinDensity` off the base and onto `tSpinResolvedWF<T>`, a data-free cross-cast face on the `tSpinResolved_CD` model; the raw `new` went with it | `95640bca` |
| **V1.32** | ✅ `IrrepCD<T>` → non-template `FiniteIrrepCD`; three one-live-branch conditionals collapsed, `IrrepHF_PairBase` deleted.  `PeriodicIrrepCD<T>` asked and DECLINED | `9a073f39` |
| **V1.14** | ✅ both `Emit*()` faces deleted under the USER'S REPORTING RULING (each class reports at its OWN trigger; nobody tells another when); basis built INSIDE the run so shells self-announce.  ⛔ the row's "pulling reporter + toggles on SCFParams" fix was the wrong direction; the `bool&` toggles are the design | `fe78682a` `d5d42fb4` |
| **V1.12** | ✅ `EnergyBreakdown` → keyed role-tagged contributions + unsummed diagnostics + `ChargeBreakdown` seed; band form reachable (throws where TrDV is not claimed).  User rulings: −TS is an energy; second map OK; roles not prefixes | `e1ac8527` |
| **V1.2** | ✅ `Orbital_PP_IBS` + the `SpeciesField` vocabulary in qcBasisSet, `Math::Gaussian` in qcMath; the PP models implement the faces; `qcLattice_BS → qcPseudopotential` REMOVED, `qcPseudopotential → qcBasisSet` added.  Bit-identical | `fd7f8099` |
| **V2.2 / V2.5** | ✅ **SPRINT S, 2026-09-14** — GPW seeds `IonicSAD` by default (`Uniform` opt-in; Al states it; a missing species throws with the way out named) and the KB mesh fallback floors its cutoff at `C α_max + α_β` via a new `SpeciesProjectorSet_R::SharpnessR`.  ★ The two "anchor-moving" rows moved NO anchor: 860/860 with zero re-pins (converged energies agree at printed precision; Si Γ runs 17 → 8 iterations).  ★ V2.5's sweep EXONERATED the analytic l=2 KB and re-enabled the five-week-disabled d-channel gate | (this session) |
| **V1.18** | ✅ WIDENED into the density-mixer reorganisation (user code review, 12 points) and executed in SEVEN increments: module split, pure faces, α_eff deleted, adaptive step → one cross-cast method (re-fetch dropped), DM source → provenance seated by the driver, `ΔG_Map` operators, three named factories with param structs, mixers own their field.  ⏳ (g) the PolarizedRunKeepsItsSpin unit test is TE work | `e60087bd`…`ead8bfcb` |

★ **AND THREE ROWS WERE CREATED BY THIS WORK**, all live in `CleanupCandidates.md`: **R1.0r** (ρ is
star-averaged under a bigger group than the k-mesh has), **V1.35** (the axis fusion — needs a PLAN, not a
session), ~~**V1.20d**~~ (closed 2026-09-13 — it had never actually been filed as a row).

⇒ **GROUP A IS CLOSED (2026-09-14): R1.0j, R1.0h and V1.37 all landed; V1.12 and V1.33 before them; V1.38 is a
STASH, not a question.  Group C is untouched (and must stay whole — see its note).  GROUP D IS CLOSED (2026-09-13): V1.17
and V1.32 closed 2026-09-10 exactly the way their rows said; V1.14 closed 2026-09-11 the OPPOSITE way — the
user's reporting ruling reversed the row's proposed fix; V1.18 closed 2026-09-13 as something much BIGGER
than its row — the user's code review of the file turned it into a seven-increment reorganisation.**  ★ That is the
finding worth carrying: group D was labelled *"genuinely open interface questions — just judgement"*, but
the judgement had ALREADY BEEN MADE in each row and then left un-executed.  The two closed here needed no
new decision — one pointed at an idiom already in the tree, the other carried its own pre-ruling on how far
NOT to go.  **Before treating a group-D row as open, check whether its verdict is already written in it.**

▶ **The rows that ARE obvious are deliberately not listed.**  If a row is a one-liner, do it; it does not
need a page.

## ▶ WHAT IS LEFT, AND WHAT TO DO NEXT (2026-09-14)

Groups A and D are closed and E was never work, so this file's remaining content is **B (blocked) and C
(anchor-moving)** — five rows, and none of them is a design question any more:

| row | group | state | what unblocks it |
|---|---|---|---|
| **V2.3** | B | ✅ CLOSED 2026-09-14 — the gate was run: no throw, 6e-9 Ha off the unpolarized answer in the same 17 iterations (`PolarizedSingletMatchesUnpolarized_PWFitRaster`) | — |
| **V2.2** | C | ✅ CLOSED 2026-09-14 (sprint S) — default flipped; nothing to re-pin | — |
| **V2.5** | C | ✅ CLOSED 2026-09-14 (sprint S) — the floor was missing altogether, not just its \f$\alpha_{pp}\f$ term; the d-channel gate is back | — |
| **V1.34** | B | half-realised `FitContraction<U,TFit>`; `bad_cast` in Release on a real TRIM block via the ball route | **N4**'s verdict on whether the ball-fit route survives |
| **R1.0b** | B | SP/"L" shells in the Gaussian94 reader | the `PG_Cart::IrrepBasisSet` same-exponent merge bug |

**In order:**
1. ~~**V2.3 first**~~ ✅ done — it deleted the row (the fix had been in the tree since 2026-08-28).
2. ~~**The sprint S = V2.2 + V2.5 together**~~ ✅ done 2026-09-14 — and cheaper than even the re-cut
   feared: the one full `ctest` pass re-pinned NOTHING.  A converged run lands on the same number from
   either seed, and the KB floor binds on no production run.  What "anchor-moving" had really been
   describing was the COST OF FINDING OUT, which is one sweep.  Records in `doc/CleanupHistory.md`.
3. **V1.34 and R1.0b stay blocked** on things outside this sweep (N4; the reader bug).  Neither is worth
   forcing: V1.34's honest fix depends on a route decision, R1.0b's payoff is the S3b spherical lineage.

⇒ **2 is done, so this file is retired** (`doc/OldPlans/`), and the programme's step 2 hands off to step 3
(**TE**, the test-suite axes — `PolarizedRunKeepsItsSpin` early) and then DFT+U, per `doc/OpenWork.md`.

---

## A. Blocked on a design question nobody has answered yet

### R1.0j — ✅ DONE (2026-09-09/10, closed with R1.0h 2026-09-14) — "XC quadrature" is misnamed and does too much
Renamed `DensitySampler` (`c7272a22`), relocated to `qcChargeDensity` (`bd11b903`), both assembly halves
through `MatrixForward`/`MatrixAdjoint`.  The three tenants resolved WITHOUT a scope (see R1.0h): `SiteMoments`
LEFT (the sampler keeps the quadrature op `SiteIntegrals(f)`; the term computes the observable; `ChargeBreakdown`
carries it); `Matrix` STAYS (the adjoint pairing `LatchRoute` guards); `Integrate` STAYS — the singles strategy can
hold NO mesh, so `qcMesh::Integrate` would add a branch.  ★ What R1.0j(4) called "the density↔operator seam" is
what the class now is.

### R1.0h — ✅ DONE 2026-09-14 (`018aa3e2`) — the \f$H_{ij}\f$ cache, and the scope that was not needed
Part 1 (`56f3db12` `0b753d1e`): the block loop performs no map INSERTION; five cache holders, pure hooks.
Part 2: the row's "OWNING SCOPE" was tested against the tree and had NO PAYLOAD — its memory case was ruled
out on 2026-09-09, and its tenants resolved by themselves (above).  ⇒ Declined, with the record saying so;
if it ever returns it returns for the k-parallel axis, merged with a term rewrite, never alone.  ⚠ Found in
passing: `Vcorr_QuadraturePol::GetEnergy` never set `charge.lost`, so ρ_lost/N read 0 on every polarized
GPW run since 2026-09-04 — fixed.

### V1.12 — ✅ DONE 2026-09-13 (`e1ac8527`) — `EnergyBreakdown`'s 13 public data members
Two of them are not energies: `GridChargeLost` is a GPW health DIAGNOSTIC (its own comment says so) and
`MinusTS` is WF-side entropy.  `E_alphaZ` is lattice-only, sitting in a structure-neutral struct.

**Was not obvious because** "keyed contributions + a small fixed set of roles for the totals" touches every
term, every `operator+=`, and every Display — it happened ONCE, the user's review supplied the roles and the
two-number entry (E, Tr(D·V)) that makes the band form reachable, and +U now adds one entry.  ⚠ Two of
my premises were corrected on the way: −TS IS an energy (dimensions), and a second non-summed map is a
data-structure choice the reporting ruling does not forbid.

### V1.33 ✅ DONE 2026-09-13 — see the STATE table; record in `CleanupHistory.md`, plan `doc/BasisSetTaxonomyPlan.md` (RECORD)

### V1.38 — Point spec in the core + one thin IBS class per (G, engine) (filed 2026-09-14, STASHED)
The `BasisSetTaxonomyPlan.md` §5 sequel, sized 2–4 sessions.  **Not obvious because** it should come AFTER an
ISP split of `LatticeSum1E` (the monster face it would otherwise mixin-forward), and nothing on the battery
path needs it yet.  Triggers and full row in `CleanupCandidates.md`.

### V1.37 — ✅ DONE 2026-09-14 (`b0310692` steps 1–2, `35d811ea` step 3) — Pol/UnPol are imposed subgroups, not types
Executed in one day.  ★ **The addendum's "13 `IsPolarized()` sites" undercounted the wrong thing twice**: steps
1–2 had eleven MORE casts to the abstract polarized FACE than the five concrete ones it listed; step 3's thirteen
were nine DECLARATIONS and four bool→enum conversions, and the real job was the TYPE SPLIT behind them (five Pol
term classes, five Pol Hamiltonian classes, a Spin-tagged exchange functional).  Two user rulings changed the
landing: **clean code over 1e-16 anchors** (the spin-grouped sums that made the first landing bit-identical were
removed; totals unchanged at printed precision) and **`Pol` gone, not aliased**.  Dirac Hamiltonians answer
`Polarized` always (spin inside the double group — the taxonomy's last row).  Records in `CleanupHistory.md`.

---

## B. Blocked on something else landing first

### R1.0b — shared-radial (SP / "L") shells in the Gaussian94 reader
**Blocked on** the flagged `PG_Cart::IrrepBasisSet` bug: the reader MERGES same-exponent shells across
\f$l\f$, which is why valgen carries the standing "keep exponents disjoint across l" rule.  Multi-l shells
must be REPRESENTABLE before shared radials can be emitted.

⚠ **Do not assume the obvious payoff.**  Sharing exponents does not by itself remove the d-contaminant
redundancy: the contaminant is \f$r^2e^{-\alpha r^2}\f$ (n=2) while an s at the same \f$\alpha\f$ is
\f$e^{-\alpha r^2}\f$ (n=0) — independent functions.  The MEASURED cause was the NUMBER of s functions
whose span mimics the contaminant.  Shared radials buy compactness and CP2K-comparable structure;
conditioning is an experiment (`GPW_SCF.MnAtomInBoxDChannel` + the vet's λmin/cond readout).

★ **The bigger prize is noted in the row and is a different job:** wiring the SPHERICAL lineage into the
lattice sums (the parked S3b work) removes the contaminant entirely — spherical d has none — and makes the
CP2K comparisons apples-to-apples (its log for our own shell list: 55 Cartesian vs 47 spherical functions).

### V2.3 — ✅ CLOSED 2026-09-14 — polarized plane-wave Vxc does not throw, and now a gate says so
The row's mechanism (`PW_XC::itsRhoGrid` keyed on the aliased `Version()`) had been fixed since the raster
route grew its per-spin pair cache on 2026-08-28; nothing had run it.  `PolarizedSingletMatchesUnpolarized_PWFitRaster`:
6e-9 Ha off the unpolarized PW-fit answer, same iteration count.  Record in `CleanupHistory.md`.

### V1.34 — the fitter's contraction face is templated but half-realised
`Fitting::FitContraction<U,TFit>` exists for two scalars and is implemented for ONE, so a real TRIM block
on the ball-fit XC route gets `std::bad_cast` — in Release, mid-SCF.  It carried a hand-written throw whose
stated reason was *"nothing real-block reaches it"*, an observation about REACHABILITY rather than about
the interface, and it stopped being true the moment a polarized run could take that route.

**Not obvious because** the fix is the DESIGN one (*give capabilities only to types that have them*), not
an added overload — and which way it resolves depends on whether the ball-fit route survives **N4**.
⚠ **CHECKED 2026-09-10: the adjoint migration did NOT dissolve this.**  I suspected it might; it did not.
`FitContraction::Overlap` is still the adjoint path for the pair route's BALL FALLBACK and for the three
molecular terms, so the half-realised face is still reachable.  Row stands as written.

---

## C. Anchor-moving — cheap to code, expensive to land

Each is a few lines, but every one re-seeds or re-sizes something pinned energies depend on.  ⇒ They want
ONE measured re-bank, not three separate ones.  That is what item **S** (the anchor-moving sprint) exists
to schedule; do not land these piecemeal.

- ~~**V2.2 — GPW's seed defaults to `Uniform`.**~~  ✅ CLOSED 2026-09-14: defaults to `IonicSAD`, `Uniform`
  is the explicit opt-in (Al says so — no library entry), and the re-bank moved nothing.
- ~~**V2.5 — `PPMeshParams()` sizes its uniform mesh with no \f$\alpha_{pp}\f$ term.**~~  ✅ CLOSED 2026-09-14.
  The integrand named here had already left (local PP in G-space; KB analytic); the one consumer is the
  mesh ORACLE of the KB gates, and what it lacked was ANY floor — it inherited an explicit 20 Ha density
  cutoff on an α_max=36 basis, which is the whole of the "l=2 disagreement" that had a gate disabled.
- ~~**V2.1 — collapse `Delta_XC`×2 into the polarized pair at \f$\zeta=0\f$.**~~  ✅ CLOSED BY V1.37 step 3
  (2026-09-14), and better than the row asked: not "one PolarizedVxc forwarding to two", but ONE term per
  operator that asks the density for its channels — no forwarding, no pair, and NO 2× cost on closed shells
  (an SU(2) run keeps the folded doublet's single raster and hands the functional \f$\rho/2\f$ twice; the
  spin face at exact \f$\zeta=0\f$ is the scalar path bit for bit).  Moved anchors: only the fit-of-a-sum
  where `Ham_DFTcorr` had two `FittedVxc`s (roundoff).  Not an anchor-moving item after all.

---

## D. Genuinely open interface questions — no blocker, just judgement

- **V1.17 — ✅ DONE 2026-09-10 (`95640bca`), and it was the smallest row in this file, not the biggest.**
  `tSpinResolvedWF<T>` now carries `GetSpinDensity`, `tUnPolarizedWF` cannot be asked at all, and the
  owning return finishes V1.25.  ★ **The row's own sentence "the correct idiom exists one library over"
  WAS the implementation plan** — `tSpinResolved_CD` is the same shape solved the same way, and its
  comment even states the rule being violated.  Nothing had to be designed, only noticed.
  ⚠ Two lessons banked in the history entry: the measured scope was **2 implementors, 1 client** (the
  row read far heavier than it was — measure before deferring), and the client's unconditional assign
  had been correct only BY ACCIDENT of the null return, which the cross-cast forced into the open.
  ⚠ It had also been parked since 2026-08-17 for a real-TRIM session that has since finished — **a park
  note outlives its reason silently.**
- **V1.14 — ✅ DONE 2026-09-11 (`fe78682a` `d5d42fb4`), and the row's fix was BACKWARDS.**  User ruling:
  `CurrentReport` is a global so nothing is threaded; each class reports CONTEMPORANEOUSLY with its own
  activity (console order == execution order); a class telling another to emit is the defect.  So: no
  pulling reporter, no toggles on `SCFParams` — the `Emit*()` faces were deleted and each provider announces
  at its own trigger (`FillOrbitals` for usage; shell CONSTRUCTION for exponents, which meant building the
  basis INSIDE the run bracket in both facades).  `EmitGridReport` was already gone.  The `bool&` toggles
  are the design (Reporting.C names them beside the sink) and stay.
- **V1.18 — ✅ DONE 2026-09-13, seven increments (`e60087bd` … `ead8bfcb`).**  The row was the ALGEBRA half
  of what the user's code review of `DensityMixer.C` found; the review supplied the other half (layout, faces,
  factories) and became the spec.  The straddle is gone because the MIXERS OWN THEIR FIELD and `FourierMixCD`
  is a presentation built whole; the ISP half is `Kerker/PulayMixerFactory` taking a GENERIC seed.  ⚠ Two
  review premises were wrong and were resolved with the user (α_eff WAS used → deleted anyway, not physics;
  ReDamp is NOT a line-search failure → one adaptive method on a cross-cast face).  Full record in
  `doc/CleanupHistory.md`; (g) the `PolarizedRunKeepsItsSpin` mixer unit test is left for TE.
- **V1.32 — ✅ DONE 2026-09-10 (`9a073f39`).**  Small and self-contained exactly as advertised, and it
  compiled first try.  The parameter was holding up THREE conditionals with one live branch each; the
  `IrrepHF_PairBase` alias died with it.  ⚠ **One thing that looks like a fourth dead branch is not:**
  `IrrepCD_Factory`'s own `if constexpr` STAYS — a guard that PREVENTS an instantiation is not a branch
  that SERVES one, even spelled identically.  And `PeriodicIrrepCD<T>` was asked and DECLINED per the
  row's pre-ruling: the two leaves sit adjacent and look like a symmetry begging to be completed, but
  one leaf's scalar is an accident of history and the other's is physics.
- **V1.2 — ✅ DONE 2026-09-13 (`fd7f8099`).**  The probe held; the PP faces were ALREADY neutral in structure,
  so the job was vocabulary + rename + delete `Integrals_Pseudo`.  Two rulings taken in a short design
  discussion (faces live WITH the service in qcBasisSet, qcStructure untouched; the range split is an
  argument), one nit (`Math::Gaussian`, not `RadialGaussianTerm`).  The edge did not invert so much as
  reverse: qcLattice_BS no longer links qcPseudopotential, which now links qcBasisSet.  ⚠ Lesson: a row
  parked as "never started" for five weeks was a two-hour job once the probe was re-checked.

---

## E. Not tasks (so they stop reading as work)

- **V4.1 / V4.2 are WATCH TRIGGERS, not items.**  Split `CollocMemo` when a THIRD consumer appears;
  promote `SolveSPD`/NNLS to qcMath when a SECOND consumer of small dense LS/NNLS appears.
  ✅ **Checked 2026-09-09, still true 2026-09-10: neither trigger has fired** — NNLS still has exactly one
  consumer (`SymmetrizeMesh.C`), and `CollocMemo` still has two.
- **V1.24(i)/(iii) and V1.30 are assigned to the MnO dev**, not to this sweep.

---

## ▶ R1.0j RE-MEASURED — the `MatrixIntegrator` question, and what became of it

*(User, 2026-09-09: "XC quadrature has recently been worked on … I think a lot of its responsibilities have
been moved elsewhere, in particular `qcMesh::MatrixIntegrator<T>`.  Who are the consumers?  And what are its
current responsibilities?")*

⛔ **THE 2026-09-09 ANSWER WAS "NOTHING HAS MOVED YET" — and that is no longer true.  Kept because the
measurement is what drove everything after it.**  At the time, `grep MatrixIntegrator src/Hamiltonian/`
returned zero hits: the face was BUILT (R1.0l) and a second realization proved it (R1.0m), but the rewiring
was blocked and the engine carried every responsibility it had.

### ✅ WHAT ACTUALLY MOVED, 2026-09-09/10

| responsibility | then | now |
|---|---|---|
| FORWARD \f$D\to\rho(r_g)\f$ | the engine's own `Rho`/`RhoPol` | `qcMesh::MatrixForward`, via `ScalarProjector::Forward` |
| ADJOINT \f$v\to\langle i|v|j\rangle\f$ | the engine reached into `applyRawAdjoint` on the concrete tensor | `qcMesh::MatrixAdjoint`, **both routes** (`d180277b`, `a7be9b35`) |
| the integral rule | two differently-ordered `Integrate`s were possible | ONE: the integrator takes the RULE at construction |
| symmetry projection | two members, then two forwarders | gone — the client calls its own `FoldedMesh` (`d180277b`) |
| the name | `XC_Quadrature` (zero functional references, 2 of ~13 members quadrature) | `DensitySampler` (`c7272a22`) |
| the library | `qcHamiltonian` | `qcChargeDensity`, interface+factory public (`bd11b903`) |

⛔ **AND `MatrixIntegrator` ITSELF IS DELETED (`7a41cca6`)** — user ruling: nobody needs both halves, so a
named pair face only added confusion.  The census agreed: nothing held it.  ★ **The pairing guarantee turned
out to be a property of CONSTRUCTION, not of a type** — one object is built once and its two halves go to the
two clients that each need one.

⛔ **THAT ALSO REFUTED R1.0m's BLOCKER, WHICH WAS MINE.**  R1.0m said the rewiring could not be done because
"`MatrixIntegrator` assumes ONE owner holds both directions" while this tree owns them on opposite sides of a
boundary, and proposed two designs to bridge it.  Neither was needed: nothing ever needed to hold both, and
the "boundary" was just the boundary between two CLIENTS — which is what handing out two halves is for.

### ✅ WHAT WAS LEFT IN THE ENGINE, AND WHERE IT WENT (closed 2026-09-14 with R1.0h)

Of the nine responsibilities the 2026-09-09 census found, one was always legitimate (the integral rule), six
had left by 2026-09-10, and the last three resolved on 2026-09-14 — **without the owning scope the row said they
were waiting for**:

- **`Matrix`** stays: the forward/adjoint pairing that `LatchRoute` guards (RAW collocated \f$\rho_{DM}\f$ vs
  the BALL round trip minimise different functionals; a per-DENSITY decision) is the class's reason to exist.
  ⚠ `LatchRoute` is NOT redundant — an earlier guess that `BlockAdjoint` had made it so was wrong.
- **`Integrate`/`NumPoints`** stay: the singles strategy can hold NO mesh (`MakeDensitySampler(fb)` with a
  default `FitQuadrature`, used by tests), so delegating to `qcMesh::Integrate` adds a branch, not removes one.
- **`SiteMoments`** left.  The sampler keeps `SiteIntegrals(f)` — a quadrature question, the atom-partitioned
  sibling of `Integrate`; the spin-native XC term computes \f$\mu_A=\int w_A(\rho_\uparrow-\rho_\downarrow)\f$
  in its ENERGY pass and writes it into `ChargeBreakdown::siteMoments`; the SCF trace emits it once per
  iteration and the facade's detectors read it off `SCFProgress`.  The "fire exactly once per new density"
  coupling turned out to be the ENERGY PASS, which already fires exactly there.

★★★ **THE SHARING SURVIVED, and is now structural rather than a discipline.**  The XC term is ONE object per
Hamiltonian (`MakeVxcTerm` over a composite functional) with ONE sampler, so the collocation is shared between
its Fock pass and its energy pass by construction — there is no second term left to share it WITH.  ⛔ Never
"eliminate" the sampler by pushing \f$\rho\f$ back into the terms: that is how the 4.8 s/iteration on NaF
comes back.  `itsSrcVersion` (the LIVE staleness check on the DM-source route) is untouched.
