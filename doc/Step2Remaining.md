# Step 2 — the remaining CleanupCandidates rows, and why each is not obvious

Cut 2026-09-09, after batches 1–3 of the sweep (nine rows closed, three filed, `850/850`).
This is a READING AID for `doc/CleanupCandidates.md`, not a second tracker — every row here is live in
that file and the verdicts belong there.  It exists because "~24 open items" is not a plan, and because
the rows differ less in size than in WHAT IS BLOCKING THEM, which the alphabetical ordering hides.

▶ **The rows that ARE obvious are deliberately not listed.**  If a row is a one-liner, do it; it does not
need a page.

---

## A. Blocked on a design question nobody has answered yet

### R1.0j — "XC quadrature" is misnamed and does too much
Measured: the engine touches a functional **zero times** — no `ExFunctional`, no `GetVxc` anywhere in the
interface or either implementation unit; the functional lives in the TERMS (`Vxc_Quadrature` holds it and
maps it over the points).  And of `XC_SinglesQuadrature`'s members only **two** are quadrature
(`Integrate`, `NumPoints`).  ⇒ It is named for its CLIENT, not its responsibility.

**Not obvious because** the fix is a rename PLUS a relocation, and the destination is a parked decision
(R1.0e: `qcChargeDensity` vs a new leaf library).  Renaming before that is churn.
▶ See the dedicated section at the end of this file — the responsibilities and consumers were
re-measured 2026-09-09 and the `MatrixIntegrator` question has a sharper answer than this row records.

### R1.0h — the \f$H_{ij}\f$ cache
`tDynamic_HT_Imp::GetMatrix` memoizes per-`Irrep` INSIDE the block loop, and it is the **last remaining
write** in that loop after the eager-refresh phase landed.  It exists because the ENERGY pass re-asks for
the same block (`GetEMatrix` → `IrrepCD::DM_Contract`).

⛔ **`DB_Cache` was ruled out on LIFETIME, not keying** — `DB_Cache` never evicts (its own header: "allow
data sharing between separate runs"), while this memo turns over every SCF iteration, so twenty iterations
would leave twenty generations of every block.  *Ask what a cache EVICTS before asking what it keys on.*

**Not obvious because** the right shape is an explicit PER-ITERATION SCOPE that owns the matrices and dies
with the iteration (`doc/Pins.md` pin 11) — a new lifetime concept, not a relocation.  Also on the
critical path for **KP** (k-point parallelism): it is the one write that stops the block loop being
read-only.

### V1.12 — `EnergyBreakdown`'s 13 public data members
Two of them are not energies: `GridChargeLost` is a GPW health DIAGNOSTIC (its own comment says so) and
`MinusTS` is WF-side entropy.  `E_alphaZ` is lattice-only, sitting in a structure-neutral struct.

**Not obvious because** "keyed contributions + a small fixed set of roles for the totals" touches every
term, every `operator+=`, and every Display — and it wants to happen ONCE, ideally carrying +U's new terms
in rather than landing just before them.

### V1.33 — the `BasisSet` taxonomy is on the wrong axis
`src/BasisSet/{Atom, Molecule, Lattice_3D}` classifies by PHYSICAL SYSTEM; the contents are classified by
BASIS KIND (`Radial` / `Polarized{Cartesian|Spherical}` / GPW / PW).
**Evidence it has already failed:** `qchem.UnitCell` is imported at **five** sites inside
`BasisSet/Molecule/`.  A molecule has no unit cell.

**Not obvious because** it is a directory-and-target reorganisation of the largest library in the tree,
and it interacts with V1.20's family ruling (`.Internal.` marks the FAMILY boundary).

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

- **V1.17 — `tWaveFunction::GetSpinDensity()` returns null as the unpolarized answer.**  A capability half
  the hierarchy lacks, declared on the BASE, with every client null-checking a raw pointer (and the raw
  `new` behind it).  The correct idiom exists one library over: `tSpinResolved_CD` as a cross-cast face.
  ★ This one directly contradicts the spin-native-is-primary bias — the polarized WF is the PRIMARY type,
  not a special case bolted on through a nullable getter — so it is the most overdue of the set.
- **V1.14 — report-emission creep on neutral faces.**  `EmitBasisUsage`, `EmitRadialReport`,
  `EmitGridReport` (PURE — it forces every implementor), plus function-local-static
  `bool& ReportBandGap()` / `ReportGridCharge()` process-globals that LEAK STATE BETWEEN TESTS (the
  `SCFIterator` comment admits it).  Fix: a reporter/visitor that PULLS; toggles on `SCFParams`.
- **V1.18 — `FourierMixCD` tell-don't-ask + `MakeDensityMixer` ISP.**  `RhoTilde()` hands out the raw
  `ΔG_Map` and `PulayMixer` runs the whole DIIS algebra OUTSIDE the density; `SetRawRho` + an external
  `RasterKerker` is a get/compute/set straddle.  Separately, `MakeDensityMixer` takes `const tDM_CD*` but
  uses only `GetTotalCharge` + a Fourier cast — excluding the matrix-free seeds BY TYPE, not by intent.
- **V1.32 — de-template the finite `IrrepCD<T>` → `FiniteIrrepCD`.**  After lineage-as-class the finite
  leaf has exactly one instantiation and the factory's `if constexpr` makes finite-complex
  UNREPRESENTABLE, so the parameter is vestigial and the name no longer says the load-bearing thing:
  *Finite* is the identity, not the scalar.  ⚠ The row already rules that `PeriodicIrrepCD<T>` must NOT
  follow — its T is load-bearing (real TRIM vs general k).  Small and self-contained.
- **V1.2 — `Orbital_PP_IBS`, the structure-neutral PP-integral face.**  Would invert the
  `qcBasisSet(qcLattice_BS) → qcPseudopotential` dependency edge.  User APPROVED attempting it 2026-08-05;
  the feasibility probe is banked in the file's appendix.  Not blocked — just never started.

---

## E. Not tasks (so they stop reading as work)

- **V4.1 / V4.2 are WATCH TRIGGERS, not items.**  Split `CollocMemo` when a THIRD consumer appears;
  promote `SolveSPD`/NNLS to qcMath when a SECOND consumer of small dense LS/NNLS appears.
  ✅ **Checked 2026-09-09: neither trigger has fired** — NNLS still has exactly one consumer
  (`SymmetrizeMesh.C`).
- **V1.24(i)/(iii) and V1.30 are assigned to the MnO dev**, not to this sweep.

---

## ▶ R1.0j RE-MEASURED, 2026-09-09 — the `MatrixIntegrator` question

*(User: "XC quadrature has recently been worked on … I think a lot of its responsibilities have been moved
elsewhere, in particular `qcMesh::MatrixIntegrator<T>`.  Who are the consumers?  And what are its current
responsibilities?")*

⛔ **NOTHING HAS MOVED YET.  `grep -rn "MatrixIntegrator\|MatrixForward\|MatrixAdjoint" src/Hamiltonian/`
returns ZERO hits.**  The face was BUILT (R1.0l) and a second realization proved it
(`ScreenedMatrixIntegrator`, R1.0m), but the rewiring was **blocked** and the engine still carries every
responsibility it had.  See R1.0m: `MatrixIntegrator` assumes ONE owner holds both directions, and in this
tree the FORWARD is driven by the ChargeDensity (which contracts its own private \f$D\f$ and aggregates
over blocks) while the ADJOINT is driven by the BASIS, per block.  The pairing `XC_Quadrature` enforces is
enforced ACROSS an ownership boundary — which is why it needed a bespoke class, and why `LatchRoute`
exists at all.

### Consumers — there are only three, and all three are TERMS

`MakeXCQuadrature` is called from exactly one production site (`PWTerms.C:364`, inside `MakeVxcTerms`) and
from three test sites.  The `shared_ptr<const XC_Quadrature>` it returns is held by:

| consumer | file | what it uses |
|---|---|---|
| `Vxc_Quadrature` | `Imp/PWTerms_XC.C` | `Rho`, `Integrate`, `Matrix`, `NumPoints` |
| `Vxc_QuadraturePol` | `Imp/PWTerms_XC.C` | `RhoPol`, `Integrate`, `Matrix`, `NumPoints` |
| `Vcorr_QuadraturePol` | `Imp/PWTerms_XC.C` | `RhoPol`, `Integrate`, `Matrix`, `NumPoints` |

★ **All three are XC/correlation terms, and ONE engine is SHARED by a pair** — that sharing is the whole
performance reason the class exists (without it the pair re-evaluated the Bloch image sums pointwise four
times per iteration: measured 4.8 s/iteration on NaF).
▶ **So the client list does not justify the name either**: a Hartree term or a \f$+U\f$ projector would
want the same object, and neither is XC.

### Current responsibilities, by category

Interface: 8 virtuals.  Implementations: `XCQuadrature.C` 429 lines + `_Singles.C` 477 + `_Pair.C` 319.

| # | responsibility | members | is it quadrature? |
|---|---|---|---|
| 1 | **The integral rule** | `Integrate`, `NumPoints`, `FunctionIntegrals` | ✅ YES — this is the whole of it |
| 2 | **FORWARD: \f$D\to\rho(r_g)\f$** | `Rho`, `RhoPol`, `Projector`, `SampleOne`, `Refresh`, `RefreshPol` | ⛔ this is `MatrixForward` |
| 3 | **ADJOINT: \f$v\to\langle i|v|j\rangle\f$** | `Matrix`×2, `MatrixT` | ⛔ this is `MatrixAdjoint` |
| 4 | **Per-density CACHING + staleness** | `itsRhoVersion`, `itsPolVersion`, `itsSrcVersion`, the cross-invalidation rule | ⛔ policy, not quadrature |
| 5 | **Symmetry projection** | `Symmetrize`, `SymmetrizeSpin` | ⛔ belongs to `FoldedMesh` (and now IS there — these are forwards) |
| 6 | **An OBSERVABLE and its reporting** | `SiteMoments`, `PartitionedMoments`, `EmitSiteMoments` | ⛔ a physics observable that happens to be free here |
| 7 | **Cache-warming policy** | `WarmForDensity` | ⛔ the eager-refresh phase (KP) |
| 8 | **DM-source damping** | `itsXCMix`, `itsXCMixUp/Dn`, `GPW_XC_DM_MIX` | ⛔ an SCF knob living in a quadrature |
| 9 | **Route latching** (pair only) | `LatchRoute`, `itsRouteLatched`, `itsLatchedRaw` | ⛔ exists BECAUSE of the ownership split |

⇒ **Categories 2 and 3 together ARE `MatrixIntegrator`.**  That is the sharpest statement of R1.0j
available: the row said "the engine is named for its client"; the re-measurement says **the engine is a
`MatrixIntegrator` wearing six other hats**, and five of those hats (4, 6, 7, 8, 9) are POLICY that has no
business in a quadrature at all.

### What that makes the real increment

Not a rename, and not "wire it onto a `MatrixIntegrator`".  In order:
1. **Resolve the ownership split** (R1.0m gives two designs: a block-aggregating decorator that owns the
   loop the density currently owns, or promoting `LatchRoute` from a runtime check to a type).  Everything
   else is downstream — category 9 EXISTS only because of it.
2. **Then** categories 2+3 collapse onto `MatrixForward`/`MatrixAdjoint` and the class loses ~half its
   surface.
3. **Then** the rename is obvious and cheap, because what is left really is a quadrature plus a small,
   nameable set of policies — and 6 (site moments) can move to whoever owns observables, 8 (DM damping) to
   the SCF knobs where it belongs.
