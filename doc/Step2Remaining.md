# Step 2 — the remaining CleanupCandidates rows, and why each is not obvious

Cut 2026-09-09 after batches 1–3 of the sweep; **done-markers added 2026-09-10 (`851/851`)**.
This is a READING AID for `doc/CleanupCandidates.md`, not a second tracker — every row here is live in
that file and the verdicts belong there.  It exists because "~24 open items" is not a plan, and because
the rows differ less in size than in WHAT IS BLOCKING THEM, which the alphabetical ordering hides.

## ▶ STATE, 2026-09-10 — WHAT CLOSED SINCE THE CUT

| row | verdict | commit |
|---|---|---|
| **R1.0j** | ✅ renamed `XC_Quadrature` → `DensitySampler`; forward AND adjoint migrated to `MatrixForward`/`MatrixAdjoint` on BOTH routes; `Symmetrize` forwarders deleted.  ⏳ the three leftover tenants go with R1.0h | `d180277b` `a7be9b35` `c7272a22` |
| **R1.0e** | ✅ library home SETTLED and executed — `qcChargeDensity`, interface+factory public, concretes `.Internal.`, tests through the factory | `bd11b903` |
| **R1.0h** | ⚗️ **HALF DONE** — slot pre-creation landed, the block loop performs no map insertion.  ⏳ the owning scope remains (now unblocked) | `56f3db12` `0b753d1e` |
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
| **V1.18** | ✅ WIDENED into the density-mixer reorganisation (user code review, 12 points) and executed in SEVEN increments: module split, pure faces, α_eff deleted, adaptive step → one cross-cast method (re-fetch dropped), DM source → provenance seated by the driver, `ΔG_Map` operators, three named factories with param structs, mixers own their field.  ⏳ (g) the PolarizedRunKeepsItsSpin unit test is TE work | `e60087bd`…`ead8bfcb` |

★ **AND THREE ROWS WERE CREATED BY THIS WORK**, all live in `CleanupCandidates.md`: **R1.0r** (ρ is
star-averaged under a bigger group than the k-mesh has), **V1.35** (the axis fusion — needs a PLAN, not a
session), ~~**V1.20d**~~ (closed 2026-09-13 — it had never actually been filed as a row).

⇒ **Group C is untouched (and must stay whole — see its note).  GROUP D IS CLOSED (2026-09-13): V1.17
and V1.32 closed 2026-09-10 exactly the way their rows said; V1.14 closed 2026-09-11 the OPPOSITE way — the
user's reporting ruling reversed the row's proposed fix; V1.18 closed 2026-09-13 as something much BIGGER
than its row — the user's code review of the file turned it into a seven-increment reorganisation.**  ★ That is the
finding worth carrying: group D was labelled *"genuinely open interface questions — just judgement"*, but
the judgement had ALREADY BEEN MADE in each row and then left un-executed.  The two closed here needed no
new decision — one pointed at an idiom already in the tree, the other carried its own pre-ruling on how far
NOT to go.  **Before treating a group-D row as open, check whether its verdict is already written in it.**

▶ **The rows that ARE obvious are deliberately not listed.**  If a row is a one-liner, do it; it does not
need a page.

---

## A. Blocked on a design question nobody has answered yet

### R1.0j — ✅ LARGELY DONE 2026-09-09/10 — "XC quadrature" is misnamed and does too much
Measured: the engine touches a functional **zero times** — no `ExFunctional`, no `GetVxc` anywhere in the
interface or either implementation unit; the functional lives in the TERMS (`Vxc_Quadrature` holds it and
maps it over the points).  And of `XC_SinglesQuadrature`'s members only **two** are quadrature
(`Integrate`, `NumPoints`).  ⇒ It is named for its CLIENT, not its responsibility.

✅ **RENAMED `DensitySampler` / `SinglesDensitySampler` / `PairDensitySampler` (`c7272a22`) and RELOCATED to
`qcChargeDensity` (`bd11b903`).**  Renamed IN PLACE first, ahead of the relocation this row wanted to bundle
it with, on the user's ruling — the misleading name cost reading time every day and a later `git mv` was
cheap.  Both halves of the assembly now go through `MatrixForward`/`MatrixAdjoint`.
⏳ **WHAT REMAINS:** three tenants (`Matrix`, `Integrate`/`NumPoints`, `SiteMoments`) that cannot leave
separately — see R1.0h below, and the re-measured appendix at the end of this file.

### R1.0h — ⚗️ HALF DONE 2026-09-09 — the \f$H_{ij}\f$ cache
`tDynamic_HT_Imp::GetMatrix` memoizes per-`Irrep` INSIDE the block loop, and it is the **last remaining
write** in that loop after the eager-refresh phase landed.  It exists because the ENERGY pass re-asks for
the same block (`GetEMatrix` → `IrrepCD::DM_Contract`).

⛔ **`DB_Cache` was ruled out on LIFETIME, not keying** — `DB_Cache` never evicts (its own header: "allow
data sharing between separate runs"), while this memo turns over every SCF iteration, so twenty iterations
would leave twenty generations of every block.  *Ask what a cache EVICTS before asking what it keys on.*

✅ **PART 1 LANDED (`56f3db12`): the block loop performs no map INSERTION in the ordinary path.**
`RefreshForDensity` now takes the basis as well as the density and has two duties — pre-create this
iteration's per-irrep slots, then pre-warm the k-independent memos.  `operator[]` mutates the tree only when
the key is ABSENT, so pre-creating the nodes is the whole fix; the bodies became FILL-IF-EMPTY with a 0×0
matrix as the sentinel.  ⚠ There were **five** cache holders, not one.
✅ **AND THE HOOKS ARE PURE (`0b753d1e`)** — which made the compiler name three terms that were silently
skipping a phase they needed (`FittedVee`/`FittedVxc` refitting inside the loop; `FittedVxcPol` never
forwarding to its children).
⏳ **WHAT REMAINS: the OWNING SCOPE**, and it is where `DensitySampler`'s three tenants go.  ⚠ Measured
before choosing: its bounded lifetime reclaims ~6 MB on MnO against a ~500 MB run, so the memory argument is
weak — the value is a home for the tenants.  **Now unblocked**: the library move is done.

### V1.12 — ✅ DONE 2026-09-13 (`e1ac8527`) — `EnergyBreakdown`'s 13 public data members
Two of them are not energies: `GridChargeLost` is a GPW health DIAGNOSTIC (its own comment says so) and
`MinusTS` is WF-side entropy.  `E_alphaZ` is lattice-only, sitting in a structure-neutral struct.

**Was not obvious because** "keyed contributions + a small fixed set of roles for the totals" touches every
term, every `operator+=`, and every Display — it happened ONCE, the user's review supplied the roles and the
two-number entry (E, Tr(D·V)) that makes the band form reachable, and +U now adds one entry.  ⚠ Two of
my premises were corrected on the way: −TS IS an energy (dimensions), and a second non-summed map is a
data-structure choice the reporting ruling does not forbid.

### V1.33 ✅ DONE 2026-09-13 — see the STATE table; record in `CleanupHistory.md`, plan `doc/BasisSetTaxonomyPlan.md` (RECORD)

### V1.37 — Pol/UnPol are imposed subgroups, not types (filed 2026-09-13)
Spin is a factor of G (SU(2) imposed = UnPol = `Spin::None` doublet; U(1)_z = Pol = Up/Down; nothing = spinors).
ONE composite over full `Irrep`s for WF and CD, `GetChannel(Spin)` a VIEW; `tPolarized_CD`'s two-level tree goes.
Forward-incompatible otherwise with the double-group rows.  A campaign (53 files); **V1.33 landed 2026-09-13, so this is unblocked**.  Full row in
`CleanupCandidates.md`.

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

### V2.3 — polarized plane-wave Vxc throws
`PW_XC`'s `itsRhoGrid` is keyed on `cd->Version()` alone, which a POLARIZED density ALIASES across
channels (a polarized density's `Version()` forwards to its Up child).  Needs per-spin rho-grid caches —
the same trap the engine's `RhoPol` pair-cache already fixes.
Mechanically clear; downstream of the XC engine's ownership question (R1.0j / R1.0e).

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

- **V2.2 — GPW's seed defaults to `Uniform`.**  The Na-doublet campaign showed a STABLE WRONG BASIN for
  electron-sparse systems: a lone-electron doublet converged 72 mHa high with every health metric green.
  The molecular facade already defaults DFT to SAD.  Candidate: default GPW to `IonicSAD`, `Uniform`
  becomes explicit opt-in.
- **V2.5 — `PPMeshParams()` sizes its uniform mesh with no \f$\alpha_{pp}\f$ term.**
  `mp.eCut = C\,\alpha_{max}\f$, but the integrand is \f$\langle\chi_i|V_{short}|\chi_j\rangle\f$ with
  exponent \f$2\alpha_{max}+\alpha_{pp}\f$.  The floor is simply missing.  One consumer today (the
  KB-projector grid fallback), so exposure is bounded.
- **V2.1 — collapse `Delta_XC`×2 into the polarized pair at \f$\zeta=0\f$.**  User's restatement is the
  target: *"there should be only one PolarizedVxc that simply stores two abstract Vxc pointers and does
  the obvious function forwarding."*  Costs ~2× XC pointwise work for closed shells; decide with a perf
  measure.  Falls out for free if it lands: the engine's second (scalar) rho cache dies.

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

### ⏳ WHAT IS LEFT IN THE ENGINE, AND WHERE IT GOES

Of the nine responsibilities the 2026-09-09 census found, one was always legitimate (the integral rule) and
six have left.  Three remain, and they leave **together, with R1.0h's owning scope** — not separately:

- **`Matrix` + `Integrate`** cannot leave on their own without splitting the forward/adjoint pairing that
  `LatchRoute` guards.  ⚠ **`LatchRoute` is NOT redundant** — an earlier guess of mine that the per-block
  `BlockAdjoint` had made it so was WRONG.  It guards the FORWARD's route (RAW collocated \f$\rho_{DM}\f$ vs
  BALL round trip, which minimise different functionals) against changing mid-SCF: a per-DENSITY decision.
- **`SiteMoments`** needs an observable owner AND the "fire exactly once per new density" coupling that only
  the sampler knows (`EmitSiteMoments` fires inside `RhoPol`'s serial-advance branch).  `PartitionedMoments`
  is already a null-guard plus `qcMesh::SiteIntegrals` — nothing to move there.

★★★ **AND THE SHARING IS NOT NEGOTIABLE.**  The class exists so the exchange and correlation terms share ONE
collocation — without it the pair re-evaluated the Bloch image sums pointwise four times per iteration, 4.8
s/iteration on NaF, essentially the whole Becke premium (user, 2026-09-09: *"very important"*).  ⚠ Since the
2026-09-04 one-gather change `MakeVxcTerms` builds ONE term, so the surviving sharing is between that term's
FOCK pass and its ENERGY pass — **the same shape and cause as R1.0h's \f$H_{ij}\f$ cache**, which is why the
two are one job.  ⛔ Never "eliminate" the sampler by pushing \f$\rho\f$ back into the terms.
⚠ And do not drop `itsSrcVersion` on the way out: a deliberately LIVE staleness check, not an assert.
