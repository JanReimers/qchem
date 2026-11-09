# Hamiltonian cleanup sprint

PLAN, born 2026-11-06 (user: geometry optimization is Hamiltonian centred, so clear the H rows of `doc/OOD-SOLID-Cleanup.md` first).  Under
execution after the Saturday 2026-11-07 reboot.  One lane in sequence, lane S: these rows share hot files.  Retires to `doc/Records/` when
step 5 lands.  The rows stay the record of what is open; this file holds only the order, the gates and the rulings not yet on a row.

## Order

| step | rows | why here | gate |
|---|---|---|---|
| 1 | `§I.1`, `D13`, `BM(1)` | independent, small | sweep; a move is re-judged against the oracle |
| 2 | `V1.34` | ruled shape (3) not yet realised, and it blocks the rest | sweep, 0 moves |
| 3 | `R1.0q`: `FittedVxc`, then `Vee_Hartree` and the quadrature XC terms | needs step 2 | sweep, 0 moves · ✅ closed FOR NOW 2026-11-07 (`Vee_Hartree` stays on `Contract`/`ContractAdjoint`); the general return-type question is open with the user, design nothing for it (`doc/Records/CleanupHistory3.md`, R1.0q) |
| 3b | `V1.1b`: `Vee_DFT`, Dunlap and normal | see the ruling below; touches the term step 3 passes through | design note first; sweep, 0 moves on the default route |
| 4 | `V1.35` with `V1.41` stages 2-3, planned together | the `*_RealBlock` caches are in both; both change the tier-1 `Hamiltonian` hub | design notes first; sweep + panel at 4 threads; one timed MnO or TeF₆ run before and after · ✅ 2026-11-07, 0 moves (`doc/Records/CleanupHistory3.md`, V1.41 + V1.35) |
| 5 | `V-ONE-OF-N` (the coefficient record; `Vee_Hartree` joins the pair, reopening `R1.0q`), `V-XC-ONE-TERM`, then `V-NONLINEAR-MESH` | the term loses mesh ownership in the second | `VxcProductSpace`; moves re-judged against PySCF; one timed run |

**Step 4's inner order (lane S's notes `lane-s-v135-blast-radius.md`, `lane-s-v141-stages-2-3.md`):** the two `V1.41` rulings; stage 2 for `Vee_Hartree`, `Vxc_Quadrature`, `Hubbard_U`; `V1.35` as deletions (once `MakeMatrix(bs, s, DensityID)` is in, no block signature names the run scalar, so `TRun` may not be needed); then the rest of stage 2, and stage 3.  `FittedVee` and `FittedVxc` need coefficients as values, a tier-1 `Fitting.FunctionFitter` change, and do not block `V1.35`.

**`V1.41` Decision 1, RULED (user 2026-11-06): option C.**  One current state per term, filled by the refresh; the block query reads it; the value passes the `DensityID` only as a check against the stamp, live until confidence allows its removal (no client sees the ID, so it can go with no client edit).  No keyed map, no per-term value objects.  This rewords stage 3: no IMPLICIT "density I saw last", not "no state".  Two values alive at once is given up: nothing uses it (`lane-s-v141-stages-2-3.md`, addendum).  Decision 2, DECIDED (controller, delegated by the user 2026-11-06): a source-less density gets n = 0 unless frozen, frozen n being a parameter.  If the proof below hits a real path, the fallback is handing the previous occupations over explicitly, and that comes back to the controller.  Proved by making `Hubbard.C:821`'s non-zero case throw across the `_U_` GPW tests and the +U panel.

**`V-ONE-OF-N` joins the sprint (user 2026-11-07; placed at step 5, but the timing is the controller's and lane S's call).**  One coefficient record (`rvec_t` or `ΔG_Map`, one populated) is what `V-XC-ONE-TERM` needs for `FittedVee`/`FittedVxc` coefficients as values, and `MatrixForward` returns it, so `Vee_Hartree` joins the pair and `R1.0q` closes on its original target.  Lane S writes the design note during step 4's long sweeps; it touches tier-1 `Fitting.FunctionFitter` and maybe `qchem.BasisSet` (`Projector3`), so the user sees it first.  The note answers: the record's name and members; clients reach it only through dispatching functions (`IrrepCD` reads `dense` today); a mismatch throws (`Contract`'s assert); whether `Overlap`'s dense / diagonal / identity rides now or waits for sparse PW.  Pin 32 applies to every fitted term it touches.

**Step 4 hub surface APPROVED (user 2026-11-07)** as lane S proposed it: `DensityID` moves to `qchem.ChargeDensity`; `tDynamic_CC<T>` gains `GetMatrix(bs, s, DensityID)`, and `GetEMatrix` takes `DensityID`; the `*_RealBlock` classes and the 31 `*R` methods go; in stage 3, `tDynamic_HF_HT<T>::GetMatrix(bs, s, DensityID, wholeBasis)`.  **Added (user 2026-11-07, option (a)):** `tDynamic_HF_HT<T>::RefreshForDensity(wholeBasis, cd)`, so the HF J/K is built in the refresh and not lazily on the first block; the wholeBasis latch moves there.  HF stage 2 runs after `Hubbard_U`.  **Corrected (user 2026-11-07):** `RealBlockEnergy_CD::DM_ContractE` (its own face in `qchem.ChargeDensity`, not on `Dynamic_CC_RealBlock`) becomes `DM_ContractE(const tDynamic_CC<double>*, DensityID) const`; the class stays, `tChargeDensity::DM_Contract` unchanged.  **Stage 3 (user 2026-11-07):** `PrepareSlots(bs, cd)` becomes create-AND-fill (every block for the current density, built at the Fock refresh; the block query a pure read), RENAMED for what it does; `InvalidateCache` goes; the energy pass's `GetEMatrix` of the Hartree fit terms computes from the energy refresh's fit, no memo (controller).

**Step 5 also adds `VxcFit::Multipole` as a column of `VxcProductSpace` (user 2026-11-07):** `AdmitVxc` admits it and `GridHartree.C` tests it, but the 27-cell gate never sees it.

**`V-ONE-OF-N` design APPROVED (user 2026-11-07)** on lane S's note (`doc/Records/VOneOfNDesign.md`; its §9 is the sub-step order): non-template `FitCoefficients` (`rvec_t` or `ΔG_Map`, one populated, no public member, a mismatch throws); `MatrixForward::Forward` returns it and `Adjoint` takes it (`Vee_Hartree` joins the pair); `tChargeDensity<T>` gains `Coefficients(fit)` / `RepulsionCoefficients(fit)`, and the `FourierDensity` face and `FourierDensityBase` go (D2); `GetRhoOnGrid`'s clients take `Coefficients(fit)` (R1.0n/o); `Projector3` → `PairProjection<T>`; `Overlap`'s identity / diagonal waits for sparse PW.  `FittedVxc` is NOT pin-32 consistent (H from the v fit, E from the ε fit): step 5's `VxcProductSpace` work, not settled by the record.

Not in the sprint: `V1.39` (its own trigger is GGA), `V-BLOCKED-MATRIX` (it may fall out of step 4; decide there).

**While in each file (user 2026-11-06, `OOD-SOLID-Cleanup.md` header):** concrete classes not seen by clients are low priority, but look for the same sequence of code appearing in two structure-specific classes, and log or fold it.

## How a step is validated

No row needs a material campaign.  Those came from numerical questions; a cleanup asks whether the numbers stayed the same.

- **0 moves expected:** short + long, `scripts/anchordiff --tol 1e-10`.  About 1.3 hours.
- **Inside the Fock loop** (`V1.41`, the `R1.0r` magnetic branch): also the panel at 4 threads.  Hours, unattended.
- **Moves on purpose** (`V-XC-ONE-TERM`): re-judged against PySCF on small molecules (pin 10).
- **Not seen by a sweep:** wall time, so one timed run where caches move; threads (baselines are serial); opt-in routes, which need an
  energy gate, not only a pointwise one (the radial Poisson bug, `b0f1ea7e`).

## Design notes the user sees first

`Hamiltonian` is a tier-1 hub, so these gate the sprint more than validation does.

1. **`V1.35`: the blast radius (user 2026-11-06).**  `tDynamic_HT<TBlock,TRun>` looks clean by itself.  The real question: do the two
   parameters propagate through existing code, or do many classes become templates that were not?  That may be fine, but see it first.
   The note counts it: which classes change, which become templates, which stay concrete because the parameters stop at a boundary, how
   many module interface units are touched, and what that does to build time.  `*RealBlock` is named in about 20 files today.
2. **`V1.41` stages 2-3:** the `MakeMatrix(bs, s, DensityID)` term interface (ruled), and what moves per term.
3. **`V1.1b`:** below.

## Ruling: `V1.1b` (user 2026-11-06, to be copied to its row)

The Dunlap fit is not named explicitly today, and it should be.  Target: two versions of the Hartree term (`FittedVee`, or `Vee_DFT` under
the new naming), one Dunlap and one normal.  Each creates the correct fitter implementations and calls the correct 3-centre integrals.
Today one `FittedVee` serves both and tells them apart by `itsGridCharge`, and the energy formula 2·EeeFit − EeeFitFit is shared without saying
what it needs: that is the defect.  It needs a real ρ̃ with a solved potential, so Dunlap and multipole qualify and δ does not (user 2026-11-06: 'only
Dunlap' was too narrow a framing).  The two versions are the shapes in `doc/Records/HartreeFitConsistency.md`: Dunlap Cᵀ G⁻¹ C, second order;
normal Oᵀ K O, first order, with `K` a `CoulombOperator`.

**Open for the note: can they be structure neutral?**  Wanted, if achievable.  The likely obstacle is G = 0: a periodic Coulomb metric
diverges for a charged fit function, so a periodic Dunlap fit needs its neutrality stated, where a normal fit in G-space drops G = 0 in `K`.
Say what each version does on a molecule, an atom and a solid, and where neutrality truly fails, name it.

**Accepted (user 2026-11-06, on lane S's note `lane-s-v11b-two-hartree-terms.md`):** the normal term asks K for its energy (δ gives ½q·v, multipole its two-leg form, 0 moves), since 2·EeeFit − EeeFitFit needs a real ρ̃ with a solved potential, not the Coulomb metric as such.  Dunlap is refused on solids, a second stated exception to "every combination" beside multipole × periodic; on the plane-wave fit basis the two are one formula.

## Term names (user 2026-11-06: `Vee_Hartree` is redundant, `Vxc_Quadrature` likewise; RULED, the table below)

`<Operator>_<Theory>[_<Variant>]`: the theory is HF or DFT, and a variant appears only where two realizations differ in FORMULA, never for structure or grid.

| today | ruled |
|---|---|
| `Vee` | `Vee_HF` |
| `Vxc` | `Vxc_HF` (ruled) |
| `FittedVee` on a Gaussian fit | `Vee_DFT_Dunlap` |
| `FittedVee` on a grid fit, `Vee_Hartree` | `Vee_DFT` |
| `FittedVxc`, `Vxc_Quadrature`, `tVxc_Quadrature` | `Vxc_DFT` (ruled) |
| `PP_Local`, `PP_NonLocal`, `Ven_PP_Short`, `Ven_PP_Long`, `Ven_PP_NonLocal` | `Ven_PP_Local`, `Ven_PP_NonLocal`, keeping `_Short` / `_Long` only while the split is real |

**When the XC and Hartree renames can land (lane S 2026-11-06).**  A table row with several sources is several live classes in one namespace, so the name waits for them to merge.  `tVxc_Quadrature<T>` takes `Vxc_DFT` when V1.35 deletes the `Vxc_Quadrature` subclass, with `FittedVxc` following in step 5.  `Vee_Hartree` takes `Vee_DFT` when it merges with the molecular term after step 4.

Unchanged: `Ven`, `Kinetic`, `DiracKinetic`, `RestMass`, `IonIon`, `Hubbard_U`.  Renames ride with the step that touches the class, not as a separate sweep.
