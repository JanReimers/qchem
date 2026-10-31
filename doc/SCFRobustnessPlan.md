# SCF robustness — the design note for chunk C1 (DRAFT 2026-10-29, lane A; awaiting the user's ruling)

Born at the top of `doc/` because C1 executes it (`doc/QueueCampaignPlan.md` §3, §4 C1).  Retires to `doc/Records/` when C1
merges.  **No code changed in this step.**  §1 is the evidence, §2 the policies as checkable rules, §3 the regime matrix,
§4 the panel, §5 the build order, §6 the questions only the user can answer.

Sources read: `Records/OpenWork_History5.md` (all closed rows since 2026-09), `Records/HubbardUHistory.md` (2026-09-20 … 2026-10-29),
the open §3 rows of `OpenWork.md`, `Pins.md` 15/18/23/27/28/29, `Records/SCFStrategyPlan.md` §3–§4, §10 and its appendices, OOD D9 /
V-MIXER-GBALL / V1.27, CleanCode V1.24 / D-DIIS-FDMIN / D-CALC-DOUBLE-SCF, and the code they name (`SCFParams.C` `ResolveMixing`,
`SolidCalculation::MakeAccelerator`, `SCFIterator.C` `RequirePin28`, `PolarizedDensityMixer.C`, `SCFCap.C`).

---

## 1. Incident inventory — 53 incidents, classified by ROOT CAUSE

"Late" = how long the defect cost before it was named (iterations of a run, or calendar time from the first symptom).  "How" = how
it was found.  "Caught by" = the policy of §2 that would have stopped it at iteration 1 (or at construction), if any; **—** = a
physics fact or an open question that no policy catches, only the panel can track.

### A. Two components keep histories of different objects (P1)

| # | incident (date found) | late | how diagnosed | caught by | record |
|---|---|---|---|---|---|
| A1 | MnO free: Fock DIIS (Ladder) + density Pulay extrapolate each other's output; max\|ΔD\| never left 1e-2 (2026-09-20) | weeks of "MnO is hard"; 120-iteration caps | user's "CP2K finishes minutes earlier" antenna → compare loop SHAPE with CP2K's deck | P1 (refuse the pair) | HubbardUHistory 2026-09-20 ★★ |
| A2 | Pulay on ρ beside DIIS on F: free Si k222 stalls, Si Γ diverges (2026-10-22) | one probe session | A/B probes for the Auto mixer default | P1 | `SolidCalculation.C` `accelerator` doc |
| A3 | Polarized mixer ran ONE Pulay per spin channel: c↑≠c↓, synthesised moments, MnO run-27 ejections (2026-08-10) | days | E/m_stag trajectory, then reading the class | P1 (one extrapolator over all channels) | `PolarizedDensityMixer.C` header |
| A4 | U outer loop: a continued `Iterate` let the converged Pulay history extrapolate the first post-U step back onto the old density → "converged in 1 iteration" (2026-09-21) | caught on the first loop | suspicious 1-iteration convergence | P1 (a stage hand-off drops history) | HubbardUHistory incr. 3 slice D |
| A5 | GDM's mixed fallback built its Fock from the field mixer's STALE running field (+2.6 mHa kicks, three silent fallbacks per run) (2026-10-27) | months latent; one probe | NaF fdprobe trace, XC column `shadow` | P1 (hand-off re-seats) + P2 (density lineage) | History5 "GDM's mixed fallback" |
| A6 | Exact resume of a field-mixer state: a fresh Pulay primed with plain Kerker AT the fixed point walked 0.29 Ha away (2026-10-29) | 23 iterations per resume | probe from a saved state, ruling out restore and projection first | P1 (hand-off from "restored" to "mixer" is explicit) | History5 "An EXACT RESUME … walked away" |
| A7 | MOM on a smeared degenerate d frontier: unequal occupations swapped every iteration, D rotates at constant ρ (2026-09-20) | weeks (same as A1) | same CP2K loop-shape comparison | P4a (MOM off in the magnetic row) | HubbardUHistory 2026-09-20 |

### B. XC fed the wrong ρ (pins 18, 28)

| # | incident | late | how | caught by | record |
|---|---|---|---|---|---|
| B1 | XC fed a separately damped ρ (α_eff): four collapsed MnO states, −45.5 … −38.5 vs −61.403 (2026-08) | weeks | dose-response in E_ee | pin 18 / P2 (already a throw) | pin 18 |
| B2 | Fitted ρ into XC: O₂ box +7.72 mHa; 15 % negative XC points on MnO (2026-10-22) | months latent | the `XC ρ` trace column | pin 28 `RequirePin28` (throw since 2026-10-26) | pin 28 |
| B3 | Pulay without Kerker mixed the fitted ρ̃ (2026-10-22/24) | 2 days | pin 28 audit | `RequirePin28` | History5 "Pulay without Kerker" |
| B4 | PW XC ball at 1× Ecut: every PW number missing the G_max…2G_max shell (+60 μHa) (2026-10-26) | all PW history | **`RequirePin28`'s crash** — the throw found it, not a person | `RequirePin28` (worked as P2 intends) | pin 28 |
| B5 | A carried density's (and EVERY field-mixer run's) first Fock fed XC the fitted ρ (2026-10-28) | all field-mixer history | restart opening off its own fixed point | P2 (pin 28 at iteration 1, including the seed-from-D case) | History5 "A CARRIED density's first iteration" |
| B6 | +U on the polarized ρ̃-mixed density: the term found no D, applied n = 0 on alternate refreshes; the oracle PASSED (2026-09-20) | one day | the `[+U]` trace, not the oracle | P2 (density lineage: every consumer can name its D) | HubbardUHistory 2026-09-20 "PASSED WITH A DEFECT" |

### C. The grid / raster (P3; most of the fix is C2)

| # | incident | late | how | caught by | record |
|---|---|---|---|---|---|
| C1 | `BallOnly` + C = 2 aliases tight pairs: β-MnO₂ + Mn d α 76.7 "converged" at −6917 Ha (2026-10-25) | hours | unit closed forms at any α | P3 (raster's alias estimate vs need) + a physical-bound check | OpenWork §3 "⛔ GPW BallOnly" |
| C2 | Same raster breaks equivariance 1e-6 vs 9e-9: TiO₂ littleGroupK, Etot 0.54 mHa (2026-10-28) | "the AliasFree day" | KAxisPlan P4 orbital-image residual | P3 (littleGroupK declares 1e-6) | same row |
| C3 | Uniform-raster XC not equivariant even IMPOSED (2.8e-4 / 7.4e-3 Ha) (2026-10-27, OPEN) | unknown; energies hide it | KAxisPlan P4 | P2 (group residual of F on an imposed run) | OpenWork §3 |
| C4 | FREE run's Becke mesh not group-invariant (5e-4 → 2.6e-9) (2026-10-27) | unknown | KAxisPlan P4 | P2 | OpenWork §3 (cost remainder open) |
| C5 | AUTO density cutoff misses O q6 (+27 mHa O box) (2026-10-23) | hidden in every Mn cell | O box vs CP2K | P3 (PP field sharpness as a declared need) | OpenWork §3 |

### D. Restarts and carried stages (P4a rows 4–5)

| # | incident | late | how | caught by | record |
|---|---|---|---|---|---|
| D1 | Iteration-0 re-fill at kT = 0 for a carried density: FeS₂ first iterate 0.93 Ha off, SrVO₃ 14 Ha (2026-10-21) | every carried stage since the iterator existed | first-iterate Δρ printed on FD stages | P2 (a carried stage's iteration 1 within the saved residual) | History5 "re-filled at kT = 0" |
| D2 | Ladder restart entered on DIIS+Kerker: 9/12 iterations where GDM takes 2 (2026-10-02 → ruled 10-21) | 3 weeks | NaF restart A/B | P4 (Restart recipe) | CLAUDE.md, History5 |
| D3 | A smeared restart from an integer state handed GDM an idempotent D′: engaged, never applied smearing, 6 mHa off (2026-10-21) | one gate | Si k222 kT 0.02 gate | P4a (GDM is kT = 0 only, pin 15) | History5 |
| D4 | relax 1.0 diverged a restart (Si Γ 2026-09-29); the cap of the START let the adaptive α grow back to 1.0 (2026-10-29) | a month apart, same defect twice | FD test period-3 cycle under the cap guard | P4 (recipe caps the CEILING) | History5 "Eight SCF tests" |
| D5 | Exact resume: the carried D was never the working density (two undamped Fock steps) (2026-10-29) | 23 iterations per resume | β-MnO₂ r015 | P2 (iteration 1 of an exact resume ≈ saved E) | History5 |
| D6 | Restart fingerprint incomplete; functional hard-coded "LDA" (gaps 1–4 DONE 2026-10-30; the run-POLICY remainder is an OpenWork row) | — | user review | P3 (producer reports what it delivered) | OpenWork §3 |

### E. Imposed symmetry that leaked (P2)

| # | incident | late | how | caught by | record |
|---|---|---|---|---|---|
| E1 | Shubnikov flip ops projected only in XC; odd coarse multigrid levels seeded a flip residual 4e-9 at iteration 1 that grew to 6e-2 (2026-10-28/29) | **found at iteration 300**; 6 runs × 1–1.7 h; ~2 days | a 12-digit `[+U]` trace on a local build | **P2: group residual of every density consumer at iteration 1 (4e-9 vs 1e-11)** | History5 "IMPOSED Shubnikov group" |
| E2 | V_short image list centred on the operator: translation-equivalent Mn got blocks 2e-3 Ha apart; free β-MnO₂ broke its own symmetry; MnO FM 38 mHa below AFM (2026-10-23) | months | site-resolved oracle comparison | P2-like unit guard (pin 29: translation-equivalent sites get identical blocks) | pin 29 |
| E3 | Γ supercell broke its own translations ~2e-4 (2026-10-23) | — | χ₀ of translation images | pin 29 guard | History5 |
| E4 | Imposing the CHEMICAL group on an AFM density demagnetises it (V1.28; standing warning) | — | design review | P2 (stated decoration vs the imposed group) | `SolidCalculation.C` `imposeSymmetry` |

### F. GDM (P1, P2, P4a)

| # | incident | late | how | caught by | record |
|---|---|---|---|---|---|
| F1 | GDM silently declined EVERY unpolarized run (g = 2 not handled); "use GDM" read as satisfied (2026-09-29) | since GDM existed | reading `UseFD`'s engagement tests | P2 (an asked accelerator that declines says so at iteration 1) | History5 |
| F2 | Stale virtuals after a GDM-finished stage: virtual ε off 2.7 Ha, χ₀ 400× too small, not reproducible (2026-10-27) | months latent | NiO χ₀ restart A/B, then an eigenpair audit | P3 (the response declares an eigenpair residual — now `EigenpairResidual` refuses > 1e-4) | History5 |
| F3 | GDM false-converges on a carried stage: Δρ measures its own decaying step; U never moves the state (2026-10-26, OPEN, mitigated) | 12 outer steps | Si U loop: energies unchanged to 1e-10 | P4 (a direct-min stage gates on [F,D], never Δρ) | OpenWork §3 |
| F4 | GDM crawls: linear tail −5…−10 %/step at a degenerate frontier AND at an open gap (2026-10-23/27, OPEN) | 100–200 iterations per run | trace reading; aufbau re-pick tried and refuted | panel (budget 15) | OpenWork §3 |
| F5 | Aufbau re-pick (built 2026-10-24): different state, still crawling; kept OFF | one day | dimer + m0 gates | — (refuted premise: no integer aufbau solution) | HubbardUHistory 2026-10-24 |
| F6 | Ladder exhaustion: DIIS ground 30+ iterations after GDM was vetoed (2026-08-13) | 30 × 30 min | run 49's tail | P4 (stop on exhaustion, done) | SCFStrategyPlan appendix |
| F7 | Line search always starts at t = 1: 3–4 Fock builds per GDM step (2026-08) | standing | NaF trace | panel (cost) | SCFStrategyPlan "GDM" note |

### G. The linear / adaptive mixer's control law

| # | incident | late | how | caught by | record |
|---|---|---|---|---|---|
| G1 | Linear D-mixing diverges on free Si without a Fock accelerator: α 0.30 → 0.45, never re-damps (2026-09-20, OPEN) | DIIS masked it for months | `scripts/retake5a` | panel | OpenWork §3 |
| G2 | Mn–O dimer linear + DIIS 2-cycles −29.4 ↔ −18.3 Ha, α flapping 0.8 ↔ 1.0 (2026-10-23) | one run | trace | panel | OpenWork §3 (same row) |
| G3 | Null accelerator reported [F,D] = 0, so `LinearMixer::Adapt` grew α to 1.0 at iteration 1; Na2Box an undamped 2-cycle, 2835 iterations (2026-10-29) | months latent | the [F,D]-one-home work | P2 (a gate quantity that is identically 0 is a defect) | History5 "[F,D] has no single home" |
| G4 | ×0.8 / ×1.5 with a clamp returned EXACTLY to α = 1.0: period-3 limit cycle (2026-10-29) | — | cap guard | panel (fixed: irrational ratio) | History5 "Eight SCF tests" |

### H. Gates that measured the wrong thing (P4)

| # | incident | late | how | caught by | record |
|---|---|---|---|---|---|
| H1 | Our G-space residual ~100× looser than CP2K's max\|ΔP\|: "converged" at 43 while max\|ΔD\| sat at 3e-3 for 200 (2026-09-20) | weeks | CP2K loop comparison | P4 (the stop gate is part of the recipe) | Benchmark rule 3f |
| H2 | DIIS combined irreps √Σ, GDM max: `minFD` meant a different thing per rung (2026-10-28) | — | reading | done (worst irrep, one home) | History5 |
| H3 | **An iteration cap read as a verdict — three wrong physics conclusions** (pin 23); 54 tests "passed" on an unconverged iterate (2026-10-29) | months | the `qchem.SCFCap` guard | done (cap fails the test) | pin 23, History5 |
| H4 | Δ[F,D] default 1e-7 silently demanded [F,D] ~ 1e-7, overriding MinFD (2026-10-29) | — | the one-home work | done (retired) | History5 |
| H5 | Δρ "quadrature floor" and "inherited [F,D] gates": both REFUTED — the stalls were the GDM crawl (2026-10-26/27) | 4 days of a wrong hypothesis | sweep-log statistics | — (the panel would have shown the crawl is the cause) | History5 |
| H6 | `MinΔρ` compares three different quantities across mixers (H5 row, open in CleanCode lineage) | — | reading | P4 (one gate definition per recipe) | History5 §3 |

### J. Seeds and the landscape itself (physics: only the panel tracks these)

| # | incident | late | how | caught by | record |
|---|---|---|---|---|---|
| J1 | IonicSAD targets per SPECIES (last atom wins) + no Mn⁴⁺ entry: β-MnO₂ seeded non-magnetic, 200 iterations (2026-10-22/26) | days | reading the seed | P2 (seed moment vs stated `siteSpins` at iteration 0) | HubbardUHistory 2026-10-22 |
| J2 | Decoration derived from the seed, not stated (2026-10-26; non-+U remainder OPEN) | — | review | P2 (`RequireStatedMagneticOrdering`) | OpenWork §3 |
| J3 | LDA NiO at U = 0 loses AFM-II (physics; 3 runs, 2 seeds, 1 fixed point) (2026-09-21) | — | order guard refused it | `OrderLost` (exists) | pin 23 |
| J4 | β-MnO₂ U = 0: A1AFM is a LOCAL minimum, ferri 11.7 mHa lower (QE agrees) — it is the soft direction E1 and K1 fall along (2026-10-27) | — | QE twin | — | OpenWork §3 |
| J5 | Degenerate frontier at kT = 0 has NO integer aufbau solution (dimer 3e-5…1e-3 Ha pair; m0 2.3 mHa) (2026-10-24) | — | the re-pick experiment | P4a (that regime does not run a kT = 0 stage) | HubbardUHistory 2026-10-24 |
| K1 | Kerker+Pulay overshoots β-MnO₂'s U = 1 eV basin; plain Kerker holds it; control U = 2 limit-cycles under plain Kerker (2026-10-29, OPEN) | 6 × 1–1.7 h runs, then 2 days of ladder | trace + A/B | panel (RED cases) | OpenWork §3 |
| K2 | Na₂ polarized singlet from an AFM seed at α 0.3: Δρ parks at 1e-2 on every accelerator (2026-08, OPEN) | parked 2 months | — | panel (RED) | OpenWork §3 |
| K3 | Polarized Γ Si jumps to a spin-polarized state under α = 1e-3 Ha (2026-10-26, OPEN) | — | UTResponse fixture | panel (RED) | OpenWork §3 |

### What the inventory says

1. **47 of 53 are ASSUMPTION failures** (classes A–H plus the seed defects J1–J2: one component assumed something of another), as
   §3 of the campaign plan guessed.  Six (J3–J5, K1–K3) are the landscape or still unexplained.  None of the 53 was a missing algorithm.
2. **The two most expensive were caught by nothing until a human traced digits**: E1 (≈ 2 days, found at iteration 300, a residual
   of 4e-9 at iteration 1) and A1 (weeks, found by comparing with CP2K's loop shape).  Both are exactly P2 / P1 in their first form.
3. **The one cheap catch was a throw**: B4 was found by `RequirePin28` crashing, with no person looking.  That is the strongest
   evidence for P2's "fail, don't mark" stance (the user's "simply crash the run", pin 28).
4. **Hand-offs are where state goes stale** (A4, A5, A6, B5, D1, D5, F2): 7 incidents, all at a boundary between two objects
   (stage → stage, rung → rung, accelerator → mixer, saved → live).  Each was fixed locally; no single rule names them.
5. **Gates are a recurring blind spot** (H1–H6, F3): a stop criterion that does not measure the fixed point turns a defect into a
   "converged" number.  The cap rule (H3) closed the loudest version; F3 is still open.
6. **The open rows are mostly class G/F/K** — mixer control laws and GDM — and are exactly the ones with no panel to measure a fix
   on.  Every fix so far has been verified on ONE case, which is how "fix one, break another" happens (e.g. the restart relax cap
   moved iteration counts on PW_Si / NaF / CsI / Al / Na₂, 2026-10-29).

---

## 2. The policies as checkable rules

Each rule: the statement, how to check it (a test or a throw), where it lives, what it replaces.  ⚠ = the evidence does not
support the policy as worded in `QueueCampaignPlan.md` §3; the change is stated.

### P1 — one history per run, and every hand-off is explicit

- **R1.1 At most ONE extrapolator is live at any iteration.**  Extrapolators: DIIS on F, Pulay (on D, or Kerker-preconditioned on
  ρ̃), a future Broyden, GDM's own CG/line-search memory.  Damping (linear α) and preconditioning (Kerker) are NOT histories and
  compose freely.  *Check:* construction THROWS on a Fock extrapolator (DIIS, or a Ladder with a DIIS rung) beside a density
  history (`PulayDepth > 0`), naming both and the evidence (A1, A2).  *Where:* `qchem::ResolveMixing` (`SCFParams.C`) — today its
  explicit `Mixing::Pulay` branch and the "hand-set PulayDepth is kept" branch both let the pair through beside DIIS; only `Auto`
  avoids it.  *Replaces:* the "belt-and-suspenders facade WARN" of `SCFStrategyPlan.md` §6, which was never built.
- **R1.2 A multi-channel density has ONE joint extrapolator** (one B, one c over all channels).  Already true
  (`MixJointly`, A3); the rule is the regression guard: a UT that a polarized Pulay's coefficients are channel-independent.
- **R1.3 Every hand-off goes through ONE function and prints ONE line.**  `HandOff(from, to, why)`: drops the outgoing object's
  history, re-seats the incoming one on the CURRENT working D (`tDensityMixer::Reseat` already exists for one case), and prints
  `[handoff] iteration N: <from> -> <to> (<why>)`.  The seven hand-offs in the tree today: Ladder rung change, Ladder retreat, GDM
  mixed fallback, Pulay priming → history, new stage (`BuildStage`), restart → first iteration, exact-resume Pulay-from-1.
  *Check:* grep — no `Reset()`/`Reseat()` call outside `HandOff`; a UT per hand-off that the incoming object steps exactly as a fresh
  one seeded at the current D (the pattern of `JointPulay.ReseatFieldMakesAFreshMixerAtTheCurrentDensity`).
  *Replaces:* the ad-hoc prints in `SCFIterator.C` (`[restart] exact resume…`, `[gdm] … MIXED …`) and the Ladder's own.
- **R1.4 D9: residual → per-channel preconditioner → the one extrapolator.**  `PulayMixer` today OWNS Kerker; split it into a
  `KerkerPreconditioner` decorator over a preconditioner-free extrapolator (SCFStrategyPlan §3 "a decorator, not a hand-off").  Then
  Broyden is a drop-in sibling and "plain Kerker" is (Kerker, linear) with no special class.
- **R1.5 V-MIXER-GBALL** (MOVED TO C2, controller ruling 2026-10-30: C2 re-sizes the density/pair ball from the one accuracy
  target and re-judges every moved anchor once -- the mixer's G-space is chosen there, with the FourierDensity face it needs):
  the mixer's G-space comes from the density's own representation through a face named for what the mixer
  consumes, not from `CreateVxcFitBasisSet`.  (Structural, rides R1.4.)
- ⚠ **Not "unrepresentable at the facade" yet.**  Making the illegal pair a type error means merging `SolidCalcOptions::accelerator`
  and `SCFParams`' mixer fields into one recipe value; that is P4's named-recipe work (R4.1).  Until then R1.1's throw is the guard.
  This is a sequencing note, not a rejection.
- **Extend P1 to the occupation seam?**  A7 (MOM on a smeared degenerate frontier) is a history in disguise — MOM's reference is
  the previous iterate's orbitals.  Proposed: MOM is refused beside kT > 0 (it has no meaning when occupations are fractional anyway).
  **User's call (§6 Q5).**

### P2 — invariants measured every iteration, failing at iteration 1

The set is SMALL and each costs ≤ one matrix product per block per iteration.  Every one THROWS (pin 28's rule: a mark is read by
nobody), naming the object that broke it and the iteration.

| rule | invariant | tolerance | would have caught | where |
|---|---|---|---|---|
| R2.1 | **Imposed group**: on an imposed run, the density EVERY consumer reads (Hartree's ρ̃, the mixer's input, XC's ρ, the +U occupations, the (D↑,D↓) pair) is invariant under the DECLARED group, flip ops included | 1e-10 relative | E1 at iteration 1 (4e-9 vs 1e-11), C3, E4 | one `RequireImposedGroup` in `SCFIterator`, asking each consumer through a face named for it ("the density you will read") |
| R2.2 | **pin 28** (exists: `RequirePin28`), extended to the first iteration of a carried/re-seated density | as today | B5 | `SCFIterator.C` |
| R2.3 | **Density lineage**: the density the Fock is built from is the loop's working density (a serial check), and every density consumer can name its D | exact | A5, B6 | `SCFIterator` + `tDensityMixer::FockDensity` |
| R2.4 | **Charge**: Tr(D S) = N per spin reservoir | 1e-10 · N | — (cheap, standing) | `tCompositeWF` |
| R2.5 | **Stated decoration**: an imposed or +U polarized run has `siteSpins` stated and the seed's site moments match their signs at iteration 0 | sign | J1, J2 | `RequireStatedMagneticOrdering` (extend to non-+U, the open row) |
| R2.6 | **No silent decline**: an accelerator the caller ASKED for that does not engage says so once, at the first iteration it declines | — | F1 | `tSCFAccelerator::Engageable` callers |
| R2.7 | **No dead gate**: a gate quantity that is identically zero (or not computed) is an error, not a pass | — | G3 | the iterator's column/gate code |

⚠ **Two items of the plan's P2 list are NOT per-iteration invariants and move:**
- **"Aufbau holes" is not an invariant** — at a degenerate frontier at kT = 0 there is no integer aufbau solution (J5, measured
  twice), and GDM is fixed-occupation by design (pin 15).  A hole is REPORTED at stage end ("ended non-aufbau: pair at 2.3 mHa")
  and is a P4a matter (that regime should not run a kT = 0 stage).  Recommend: drop from P2.
- **"The reference eigenpair residual" is a hand-off contract, not a per-iteration one** — mid-GDM the (C, ε) are legitimately not
  eigenpairs of F[D].  It is P3's first instance (`Response::EigenpairResidual` already refuses > 1e-4).  Recommend: move to P3.

### P3 — precision contracts

- **R3.1 A consumer DECLARES what it needs** as a small value `PrecisionNeed {worstIrrepFD, ΔρPerElectron, eigenpairResidual,
  equivariance, aliasEstimate}` (unset = does not care).  Consumers today: littleGroupK (equivariance 1e-6 ⇒ reference Δρ/N 1e-9,
  `a2e664f1`), the linear response (eigenpair residual 1e-4, `7fb46357`), forces (to be declared by the forces note).
- **R3.2 A producer REPORTS what it delivered**: the SCF stage (its final worst-irrep [F,D], Δρ/N, which recipe, aufbau or not), the
  raster (`[grid audit]`'s alias estimate, d3a487e1), the restart (the saved state's own stop level — OpenWork "fingerprint" item 4).
- **R3.3 Pre-flight compares need with promise** (deck `Vet` before the SCF, where the stage's gates are known) and either
  ESCALATES the producer (the stage tightens its own gate; later the raster refines itself) or REFUSES with both numbers.  After the
  run, the consumer re-checks the DELIVERED value (the existing refusals stay as the post-check).
- **Scope split with C2**: C1 builds the value types, the pre-flight and the SCF-side producer (generalizing `kResponseReferenceMinΔρ`
  and `EigenpairResidual` into it).  The raster/grid producer side is C2's — and C2's ruling (raster N sized from ONE accuracy
  target ε) makes most of the grid escalation automatic.  **User to confirm (§6 Q6).**

### P4 — named recipes, gated by a panel

- **R4.1 A recipe is ONE value**: (accelerator, mixer/extrapolator, preconditioner, step α and its ceiling, occupation policy, stop
  gate, cap).  Names: `Auto` (chooses a row of §3 and PRINTS which), `Insulator`, `MagneticInsulator`, `Metal`, `Restart`,
  `CarriedStage`.  Deck key `scf.recipe`; an expert override is a STATED field beside the name, printed on the banner as a deviation.
  *Replaces:* the separate `accelerator` + `Mixer`/`PulayDepth`/`KerkerG0`/`StartingRelaxRo` knobs as the primary interface, the
  parked §2 row "`MixingPolicy`", and `SolidCalculation::MakeAccelerator`'s inline Auto branch (it becomes the `Auto` row chooser).
- **R4.2 The stop gate belongs to the recipe and must measure the fixed point**: a direct-min stage gates on worst-irrep [F,D]
  (never Δρ alone — Δρ is its own step, F3); a fixed-point stage gates on Δρ/N AND [F,D].  Construction refuses a direct-min
  stage whose only gate is Δρ.
- **R4.3 Every recipe change is judged on the whole panel (§4) before it lands**: iteration counts AND fixed points
  (`scripts/anchordiff` on the panel's log against the last panel baseline).

### P4a — the regime matrix is the `Auto` spec (§3)

### P5 — the campaign stop rule

- **R5.1** A campaign session that has spent **~2 h of wall** on one run's convergence, OR has tried **three recipe variants** of
  one deck without a converged fixed point, STOPS: it files a §3 row (deck path, log path, record id, what was tried) and hands the
  deck to lane A, which adds it to the panel's RED list.  Evidence: K1 (β-MnO₂ U = 0.5–2: 6 × 1–1.7 h runs before the leak was
  traced), H5 (4 days on a wrong floor hypothesis), and pin 23's "three wrong conclusions from a cap".
- **R5.2** A campaign deck names a recipe (R4.1); a deck that hand-tunes knobs carries a `_doc` line saying which panel case
  justifies the deviation.  The `rundeck` banner prints it.
- **R5.3** The controller enforces it at intake: a new CONVERGENCE row without a panel case is not accepted (`QueueCampaignPlan` §6).

---

## 3. The regime matrix (P4a) — the `Auto` spec

Each cell cites its evidence or says **UNMEASURED**.  "Pulay-D" = Pulay on D, no Kerker (`theAutoMixing` 0.45 / 8 / 5).
"K+P" = Kerker-preconditioned Pulay on ρ̃ with the cusp deficit (pin 28 exception 2).  All rows: one extrapolator (R1.1), no MOM
unless stated.

| regime | accelerator | mixer / step | preconditioner | occupation | stop gate | evidence |
|---|---|---|---|---|---|---|
| **1. insulator, kT = 0, closed or simple open shell** (NaF, Si, TiO₂, O₂ box) | Null | Pulay-D | none | aufbau, integer | worst-irrep [F,D] + Δρ/N | Auto default 2026-10-26 ("beat linear + DIIS on every probe"); Si Γ K+P 13 vs linear divergent (G1); O₂ box `Γ_Imp_M3_Pulay_Anchor` |
| 1′. same, polish | GDM, only from Δρ/N ≲ 1e-5 | — | GDM's own (PCFloor) | fixed | [F,D] only (R4.2) | NaF: GDM from a tight state 5e-11 in 58 (History5 "Inherited gates"); **from a seed: UNMEASURED as a recipe** |
| **2. magnetic insulator, near-degenerate d frontier** (MnO AFM-II, NiO, β-MnO₂ A1AFM, Mn–O dimer) | Null | K+P 8/5, (ρ,m) channels | Kerker on ρ only; m undamped (pin 27 guard) | Fermi kT 1e-3…5e-3, ONE μ; **no MOM** | Δρ/N + [F,D]; free energy | MnO free recipe (CLAUDE.md, `decks/MnO_AFM2_free.json`); N3 (ρ,m) 18/21 vs 22/24; MOM fights (A7); pin 15a kT ≥ frontier splitting |
| 2′. same, kT = 0 polish | **do NOT** (refuse in Auto) | — | — | — | — | J5: no integer aufbau solution; F4 crawl; F5 refuted |
| 2″. same, U on a soft direction (β-MnO₂ U(O p_z) 0.5–2 eV) | Null | **plain Kerker α 0.1** where K+P overshoots | Kerker | kT 1e-3 | Δρ/N | K1: holds U = 1; **control U = 2 limit-cycles: NO KNOWN RECIPE** |
| **3. metal / smeared, kT > 0** (Al, SrVO₃) | Null (GDM forbidden, pin 15b) | K+P | **Kerker** (long-wavelength sloshing ∝ L²) | Fermi, one μ (`perKFermiLevels` false) | free-energy ΔF + Δρ/N | `GPW_Al.Γ_Imp_Anneal_Anchor` (28 iter); SrVO₃ K+P decks.  ⚠ **Auto today gives a metal Pulay-D with NO Kerker: UNMEASURED** on a large metallic cell |
| **4. restart, same Hamiltonian** (exact resume, no WARM difference) | kT = 0: Ladder@GDM; kT > 0: Null | kT = 0: linear ≤ 0.3; field mixer: Pulay from iteration 1 | as the saved run | the saved run's | as the saved run | NaF 9 → 2 (D2); r015 23 → 4 (A6/D5); Si kT 0.02 `SmearedExactResume` |
| **5. carried stage, CHANGED Hamiltonian** (U loop, ±α, kT anneal, WARM restart incl. coarse → fine grid, geometry step) | Null; GDM allowed ONLY with R4.2's [F,D] gate | Pulay-D (fresh history, R1.3) or K+P as the parent stage; linear ceiling 0.3 | as the parent stage | the stage's own (refilled under its policy, D1) | Δρ/N + [F,D] | Si U loop K+P 66 iterations; F3 (GDM false-converges, Δρ gate); ±α Kerker-Pulay 0.45 (k211 FD); coarse → fine GDM 2 iterations (NaF, CLAUDE.md).  **Geometry step: UNMEASURED** (the forces note) |

**Inconsistencies the matrix exposes in today's `Auto`:**
- A WARM restart (coarse → fine, a changed U) is treated as row 4 (Ladder@GDM at kT = 0) though its Hamiltonian changed — row 5.
  It works today on NaF because the change is small and GDM's [F,D] is then near zero; F3 shows what happens when the gate is Δρ.
  Fix = R4.2 (gate GDM on [F,D]), after which rows 4 and 5 may share GDM safely.
- Row 3 vs row 1: `Auto` does not know a metal from an insulator (it sees kT, not a gap).  Either the deck says `Metal`, or `Auto`
  reads the stage-1 gap.  **User's call (§6 Q4).**
- Row 2′: `Auto` today will put a kT = 0 restart of a magnetic oxide on GDM (row 4).  Correct only when the saved state has an open
  gap ≳ 10 mHa; at a 2 mHa pair it crawls (F4).

A cell earned by a wrong number becomes a pin (the plan's rule).  Candidates once C1 measures them: 2′ (no kT = 0 stage at a
degenerate frontier), R4.2 (a direct minimiser gates on its gradient).

---

## 4. The panel (P4)

### How it runs

- A test whose name ends in **`_Panel`** gets the ctest label `panel` (a third discovery pass in `qchem_discover_tests`, exactly
  like `_Long`, with `FAIL_REGULAR_EXPRESSION "no tests were run"`).  Run: `flock ~/Code/qchem6-runs/.heavy.lock scripts/memsafe ctest -L panel`
  from `build/Release`, keep `LastTest.log` as `~/Code/qchem6-runs/anchors/panel-<sha>.log`, and `scripts/anchordiff` it against
  the previous panel log (per test: final Etot and summed SCF iterations — the instrument already exists).
- **Two tiers**: the FAST cases (≲ 10 min total at 4 threads; run on every C1 commit that touches the SCF path) and the FULL
  tier (run before a C1 merge in the heavy slot, target ≲ 90 min at 8 threads).  How the full tier is labelled is §6 Q1.
- **The budget is an ASSERT, not a reading**: a small RAII `qchem::IterationBudget budget(N, "<row or recipe>")` on the existing
  `qchem.SCFCap` hook counts the SCF iterations in its scope and `ADD_FAILURE`s above N.  Every panel case also asserts its FIXED
  POINT (Etot to its anchor tolerance, and the order parameter / moment where the case is magnetic) — an iteration count without
  the state is how F5 "improved" the dimer onto a different state.
- **GREEN cases** pass today; budget = baseline iterations × 1.25, rounded up, set by C1's first panel run.  **RED cases** are the
  open rows; each names its row and its TARGET budget, and is marked `WILL_FAIL` in CMake so the panel is green while it fails and
  goes red the day it is fixed — forcing the bookkeeping (remove `WILL_FAIL`, close the row).  (§6 Q2.)
- Panel decks live IN THE TREE (`decks/panel/`), self-contained: every case below that today reads a basis or a state from
  `~/Code/qchem6-runs` gets that `.bsd` vendored (or registered) and starts cold or from a committed small state.

### The cases (provisional budgets: measured counts where a record has them; C1 step B1 re-measures all)

| # | case | regime | source today | status | budget | tier | row |
|---|---|---|---|---|---|---|---|
| 1 | NaF exact resume | 4 | `GPW_NaF.Γ_Imp_eqExactResume` | GREEN (2) | 3 + first iterate pin | fast | — |
| 2 | NaF coarse → fine | 5 (WARM) | `GPW_NaF.Γ_Imp_eqColdStart` (+ `SolidState.WarmStart_OtherGrid_StartsAtTheAnswer`) | GREEN (37 total, 2 on fine) | 47 / 3 | fast | — |
| 3 | field-mixer exact resume | 4 | `SolidState.ExactResume_{Pulay,KerkerPulay}_Null_StartsAtTheSavedState` | GREEN (2) | 3 | fast | — |
| 4 | O₂ box, Pulay-D | 1 | `GPW_O2Box.Γ_Imp_M3_Pulay_Anchor` | GREEN | ×1.25 | fast | — |
| 5 | Si +U loop | 5 | `GPW_Si.Γ_U_Kerker_Pulay_Cplx_LinearResponseUIsAFixedPoint` | GREEN (66) | 83 | fast | — |
| 6 | Si +U loop under an explicit Ladder@GDM | 5 | new twin of #5 | RED (F3: false-converges) | #5's | fast | "GDM false-converges" |
| 7 | free Si, linear D-mixing, Null | 1 | `GPW_Si.Γ_CP2K` with linear + Null | RED (G1: diverges) | 15 (CP2K direct-P 12) | fast | "Linear D-mixing" |
| 8 | Na₂ box polarized singlet, AFM seed, α 0.3 | 1/2 | new (H4 "TWO SMALLER LOOSE ENDS" setup) | RED (K2) | 38 (unpolarized twin 30) | fast | "Na₂ AFM seed" |
| 9 | Na₂ box Shubnikov order-lost | 2 | `GPW_Na2Box.Γ_Becke_Shub_Pol_OrderLostThrows` | GREEN (210 — itself a symptom) | 263 | fast | — |
| 10 | polarized Γ Si under α = 1e-3 Ha | 1 | UTResponse `ConvergedSi(Polarized)` fixture | RED (K3: jumps) | m = 0 state kept | fast | "Polarized Γ Si jumps" |
| 11 | B atom / O₂ singlet (degenerate partly filled level, kT 1e-4) | 2 (molecular) | `P/A_PG_DFT.Energy/Z5`, `M_DFT.OxygenTripletBelowSinglet` | GREEN | ×1.25 | fast | — |
| 12 | Al fcc anneal | 3 | `GPW_Al.Γ_Imp_Anneal_Anchor` | GREEN (28) | 35 | fast | — |
| 13 | MnO AFM-II IMPOSED | 2 | `GPW_MnO.Γ_Shub_Pol_Smear_Anchor_Long` | GREEN (14–15) | 19 | full | — |
| 14 | MnO AFM-II FREE | 2 | `decks/MnO_AFM2_free.json` | GREEN (33 at CP2K_COMPAT, Benchmark CP2K-LOOP) | 42 + AFM-II held | full | — |
| 15 | Mn–O dimer, kT 0, mult. 6, GDM | 2′ | `~/Code/qchem6-runs/MnO_box16/mno_box16_va.json` | RED (F4: crawls, 44–200) | 15 at the SAME state (Mn 5.35) | full | "GDM crawls" |
| 16 | β-MnO₂ A1AFM imposed, U(O p_z) = 1 eV, K+P | 2″ | `~/Code/qchem6-runs/MnO2/uscan/control_U1*.json` (K+P arm) | RED (K1: overshoot) | 40 (U = 0 / 4.62 converge in 16–24) | full | "Kerker+Pulay overshoots" |
| 17 | same, plain Kerker α 0.1 | 2″ | `uscan/control_U1a.json` | GREEN (converged; ≈ 60) | ×1.25 | full | — |
| 18 | β-MnO₂ control U = 2 eV | 2″ | `uscan/control_U2a.json` | RED (limit cycle, 300) | 80 (p5f3_Od converged in 78) | full | same row |
| 19 | β-MnO₂ r015 exact resume | 4 | `uscan/leak/resume_U0.json` | GREEN (4) | 5 | full | — |
| 20 | TiO₂ k222 little-group reference to Δρ/N 1e-9 | 1 | `decks/campaigns/TiO2/tio2_k222_q222_littleGroupGate.json` (SCF stage only) | GREEN (count unrecorded) | ×1.25 | full | — |
| 21 | NiO AFM-II U_in 3 eV, cold | 2 | `~/Code/qchem6-runs/NiO/step3_qmesh.json` ground stage, at Γ or k111 to fit the tier | UNMEASURED cost | — | full | — |

Notes.  #9's 210 iterations on a 2-atom box is GREEN only because the gate's subject is the order loss; it is in the panel because
it is a cheap witness of the linear step (G3).  NiO (#21) is in on the plan's list but no cold NiO run fits the tier as banked
(the k222 states are overnight runs); B1 measures a Γ / k111 arm and drops it to a `_Long`-only manual case if it does not fit.
**The plan's candidate list had no METAL**; #12 is added because row 3 of the matrix is otherwise unmeasured, and a larger metallic
cell (SrVO₃ k444) is the case that would test "Pulay-D without Kerker" on long-wavelength sloshing — UNMEASURED cost, full tier if it fits.

### B1 as built — the first baseline (2026-10-29, lane-a at 07af0909 + B1; log `~/Code/qchem6-runs/anchors/panel-07af0909-lane-a.log`)

GREEN budget = ⌈1.25 × baseline⌉.  A budget counts every SCF stage inside the case; the facade's atomic solves (the seed's, and the
+U projector's pseudo-atom) are SETUP and run under `qchem::SetupSCF`, uncounted — the first baseline counted them, +122 on each β-MnO₂ case.

| # | test (`_Panel` = fast tier, `_Panel_Long` = long) | baseline | budget | status |
|---|---|---|---|---|
| 1 | `GPW_NaF.Γ_Imp_eqExactResume_Panel` | 2 | 3 | GREEN |
| 2 | `GPW_NaF.Γ_Imp_eqColdStart_Panel` (the fine stage only) | 3 (4 on the re-check) | 5 | GREEN |
| 3 | `SolidState.ExactResume_{Pulay,KerkerPulay}_Null_StartsAtTheSavedState_Panel` | 2, 2 | 3 | GREEN |
| 4 | `GPW_O2Box.Γ_Imp_M3_Pulay_Anchor_Panel` | 25 | 32 | GREEN |
| 5 | `GPW_Si.Γ_U_Kerker_Pulay_Cplx_LinearResponseUIsAFixedPoint_Panel` | 36 (15+3+2+16) | 45 | GREEN |
| 6 | `GPW_Si.Γ_U_Ladder_LinearResponseUIsAFixedPoint_Panel` | 100 (cold 79 + 11 + 10) | 125 | **GREEN — RED claim not reproduced**: outer = 3, U_SC 1.2399 (= the DIIS route); the row closes |
| 7 | `GPW_Si.Γ_Null_Linear_Converges_Panel` (now under `policy.cp2kCompat`, the row's condition) | 17 | 22 | **GREEN — RED claim not reproduced**: the re-damp holds (α 0.89 ↔ 0.50); the row's Mn–O dimer witness stays open |
| 8 | `GPW_Na2Box.Γ_Pol_Null_Linear_SpinSeed_Converges_Panel` | 38 = cap, unconverged (E 3e-6 above the twin) | 38 | RED |
| 9 | `GPW_Na2Box.Γ_Becke_Shub_Pol_OrderLostThrows_Panel` | 232 | 290 | GREEN |
| 10 | `ResponsePolarizability.GPW_Si_Γ_Pol_OneSpinAlpha_KeepsZeroMoment_Panel` (UTResponse) | Si-p moment 0.191 e (jumped) | 1 run | RED (built 2026-10-30, B2: readout `SolidCalculation::MakeHubbardMoments`) |
| 11 | `M_DFT.OxygenTripletBelowSinglet_Panel` (the B atom arm not built) | 26 | 33 | GREEN |
| 12 | `GPW_Al.Γ_Imp_Anneal_Anchor_Panel` | 15 (13+1+1) | 19 | GREEN |
| 13 | `GPW_MnO.Γ_Shub_Pol_Smear_Anchor_Panel_Long` | 14 | 18 | GREEN |
| 14 | `GPW_MnO.Γ_Pol_Smear_Kerker_Pulay_Converges_Panel_Long` | 15 | 19 | GREEN (E −61.3032807797) |
| 15 | `GPW_MnOBox.Γ_M6_GDM_Converges_Panel_Long` | 15 = cap, Δρ/N 6.7e-2 | 15 (provisional: ruling 7, the CP2K OT state comparison first) | RED |
| 16 | `GPW_MnO2.Γ_U_Shub_Pol_Smear_Kerker_Pulay_Converges_Panel_Long` | 40 = cap, Δρ/N 5.2e-2 | 40 | RED |
| 17 | `GPW_MnO2.Γ_U_Shub_Pol_Smear_Kerker_Converges_Panel_Long` | 143 (E −93.36627496) | 179 (deck cap 180) | GREEN — but 143, not the record's ≈ 60: the record's run was WARM |
| 18 | `GPW_MnO2.Γ_U_Shub_Pol_Smear_Kerker_ConvergesAtU2eV_Panel_Long` | 80 = cap | 80 | RED |
| 19 | `GPW_MnO2.Γ_U_Shub_Pol_Smear_Kerker_Pulay_Saved_eqExactResume_Panel_Long` (built 2026-10-30, lane A: r015 vendored in `decks/states/`) | 5 (4 threads; iteration 1 at −93.3837059, 2e-6 off the saved E; iteration 2 steps 3e-4 away and returns) | 5 | GREEN, at its budget |
| 20 | `GPW_TiO2.k222_Pol_Smear_Kerker_Pulay_Converges_Panel_Long` | 39 | 49 | GREEN |
| 21 | — NiO cold | — | — | NOT BUILT (no arm fits the tier as banked) |

Wall: the whole panel ran 1 h 48 min at `-j2` (fast tier ≈ 2 min; #17 alone 82 min, #18 44 min).

**The metal measurement (ruling 4)**: Al fcc k444 (one atom), kT 10 mHa, Null, Pulay 8/5, α 0.45, Uniform seed, imposed —
Pulay on D **11** iterations, Kerker+Pulay **14**, same fixed point (−2.0774499094 / −2.0774499091).  A one-atom cell has no
long-wavelength sloshing for Kerker to damp, so this says only that Pulay-D is not WORSE on a small metal; the SrVO₃ k444 arm
(the larger cell) is unmeasured.  Decks: `~/Code/qchem6-runs/Al/al_k444_metal_{pulayD,kerkerpulay}.{json,log}`.

---

## 5. Build order inside C1 (panel first)

Each step: one commit (or a short series), sweep-green, panel re-run, `anchordiff` vs the previous panel log.

| step | what | rows closed / advanced (§9 tags) | gate |
|---|---|---|---|
| **B0** | Re-measure the check-first row: `m0_imp_Op018_PA` stage 2 virtual slot weights after `7fb46357` | "`u_bands` after a GDM stage" (F) — close if fixed | virtuals carry slot weight |
| **B1** | **The panel**: `_Panel` discovery pass + label; `qchem::IterationBudget`; vendor the out-of-tree bases/decks into `decks/panel/`; write cases 1–21; first baseline run sets the GREEN budgets; RED cases `WILL_FAIL` | instrument for everything below | panel green with the RED list as expected |
| **B2** | **P2 invariants** R2.1–R2.7 (each a throw; R2.1 first — it is E1) | "an IMPOSED polarized run without +U derives its decoration" (S) via R2.5 | each throw has a UT that breaks the invariant on purpose (the `Pin28Throws` pattern); panel unmoved |
| **B3** | **P1**: R1.1 throw in `ResolveMixing`; R1.3 `HandOff`; R1.4 Kerker as a decorator (D9); ~~R1.5~~ (→ C2, 2026-10-30); delete `DIISParams::FDMin` (D-DIIS-FDMIN); rename `FDMax` → `EngageBelowFD` (V1.24(i)) | D9, V-MIXER-GBALL, D-DIIS-FDMIN, V1.24(i) (S) | bit-identical on the sweep except where a test used an illegal pair (each judged); panel unmoved |
| **B4** | **R4.2** — a direct-min stage gates on worst-irrep [F,D] | "GDM false-converges" (F) → panel #6 GREEN | #6 within #5's budget at #5's fixed point |
| **B5** | **P4 / P4a**: the recipe value (R4.1) + `Auto` as the §3 chooser, printed; deck key `scf.recipe`; refuse row 2′ | §2 "`MixingPolicy`" (F); the CLAUDE.md recipe paragraphs become pointers to §3 | panel: every GREEN case ≤ budget under `Auto` |
| **B6** | Work the RED cases, one at a time, each judged on the WHOLE panel: (a) GDM crawl — `PCFloor` 1e-2/1e-3 and a line search that extends t (#15); (b) "GDM: should an exhausted line search END the stage?" (after Q3); (c) linear D-mixing re-damp on rising E (#7, G1/G2); (d) K+P overshoot — residual-norm safeguard / Pulay restart on rising E (#16, #18); (e) Na₂ AFM seed (#8); (f) polarized Si jump (#10) | "GDM crawls", "GDM exhausted line search", "Linear D-mixing" remainder, "Kerker+Pulay overshoots", "Na₂ AFM seed", "Polarized Γ Si jumps" (F) | each RED → GREEN with no GREEN case over budget |
| **B7** | **P3**: `PrecisionNeed` / delivered values, deck pre-flight; fold `kResponseReferenceMinΔρ` and `EigenpairResidual` into it | (prepares C2's raster producer) | the two existing refusals unchanged in behaviour, now pre-flighted |
| **B8** | Mechanical: D-CALC-DOUBLE-SCF migration; V1.27 `TraceColumns` from what the run IS, with the `ρ_lost/N` dead column (rides V1.27) | D-CALC-DOUBLE-SCF, V1.27, "`ρ_lost/N` dead column" (S) | sweep wall time down (no double SCF), panel unmoved |

Order rationale: B1 before everything (every later step is measured on it); B2 before B3 because the invariants catch a P1
refactor's mistakes at iteration 1; B4 before B5 because it removes the only reason `Auto` must tell row 4 from row 5; B6 after
B5 so each fix lands as a change to ONE recipe cell; B7 late because its main producer is C2's.

### B2 as built (started 2026-10-29 on 28fe4a31; R2.1 measured + reported, R2.4 + R2.5 + R2.7 landed, #10 built)

**R2.1 is measured where a projector already runs, not by a second projection.**  Every consumer that imposes the group
already applies a projector P to its input; \|f − Pf\| / \|ρ\| on that input is nearly free (one raster copy) and is exactly
the number that found E1 (the 4e-9 flip residual at iteration 1).  The D pair is NOT a measuring point: `tCompositeWF`
projects (D↑, D↓) by construction (3b2a163a), and a raw D_out at a partly filled degenerate shell is legitimately not
invariant.  Wiring: `DensitySampler::GroupResidual()` (singles: measured before `Symmetrize` / `SymmetrizeSpin`) →
`Vxc_Quadrature` records it with the Fock build's pin-28 source → `ChargeBreakdown::xcGroupResidual` → the iterator.
**Two stages**: (1) MEASURE: `QCHEM_DIAGNOSTICS=group` prints it every iteration; run it on panel #13 (MnO imposed) and on
β-MnO₂ A1AFM imposed, Kerker and Pulay-D; (2) THROW at the tolerance those runs set.  The plan's 1e-10 is a guess, and if
the odd coarse multigrid levels still leave a τ-break at ~1e-9 on every iteration (3b2a163a removed the AMPLIFIER, not
necessarily the SEED), a 1e-10 throw would fail every imposed magnetic Kerker run.
**Finding while wiring it**: the raster (pair) XC route does not star-average ρ at all.  An imposed run on that route
(`cp2kCompat` vetoes Becke) feeds XC an UNPROJECTED density, and R2.1 can only answer "not measured" there.  Whether that is
a hole (imposition silently partial on the CP2K-parity route) is the first question the measurement stage answers.
**R2.1 MEASURED 2026-10-30** (B2 build on 83f04c23, `QCHEM_DIAGNOSTICS=group`, every iteration):

| case | route / mixer | residual max\|f−Pf\|/max\|ρ\| |
|---|---|---|
| `GPW_NaBox.Γ_Imp_M2_eqFinite` | Becke singles, Pulay-D | 2e-15 … 7.6e-15 |
| `GPW_O2Box.Γ_Imp_M3_Pulay_Anchor_Panel` (#4) | Becke singles, Pulay-D | 1.6e-15 … 1.9e-15, **one 7.6e-11 at iteration 12** |
| `GPW_MnO.Γ_Shub_Pol_Smear_Anchor_Panel_Long` (#13) | Becke singles, Kerker+Pulay field | 4e-15 (it. 1), **3.15e-10 (it. 2)**, then decaying with the step to 1.5e-12 |
| `GPW_MnO2.Γ_U_Shub_Pol_Smear_Kerker_Converges_Panel_Long` (#9) | Becke singles, Kerker field | 2.5e-12 … 9.9e-11 throughout (143 SCF iterations; the 122 pseudo-atom solves report "not measured") |
| `GPW_Si.Γ_Imp_Pol_eqUnpol` | raster | not measured (the route does not project) |

Readings: (a) the D-space mixer is exact (roundoff 1e-15); its one spike is the low-rank D factor — at O₂ iteration 12 the
Pulay-D mix grows a 5e-9…1e-13 eigenvalue tail with DEGENERATE pivot pairs (`dm_rank`: rank 7 → 20 of 24), and the
pivoted-Cholesky cut through such a pair splits an orbit; 7.6e-11 is the size of the dropped tail.  (b) The Kerker field
mixer is NOT exactly equivariant: it leaves 1e-12…3e-10 every iteration, scaling with the step — the τ-break the 3b2a163a
fix left behind as a SEED (the AMPLIFIER is gone: nothing grows).  (c) The E1 leak read 4e-9.  **Stage 2 TRIED AND WITHDRAWN (2026-10-30)**: a throw at 1e-9 failed 8 CORRECT short-tier
runs (anchors unmoved) -- `GPW_NaF.Γ_Imp_*` read 0.10…0.51 from iteration 2 on, `GPW_Si.k222_Becke_Imp_*` 0.21 at iteration 1.
Split on NaF: the cusp-deficit correction field is clean (1e-14), ρ[D] itself is 0.51 off.  Not a defect: those runs EVALUATE
REDUCED -- the T3 stream fold scatters orbit-representative pair triples with their orbit multiplicity, "and the caller's dense
group-average ... reproduces the full collocation" (`PG_Cart_MnD/Evaluator.C`), and an IBZ k-sum is a star-representative sum
the same way -- so the pre-projector field is incomplete BY CONSTRUCTION and P is part of the evaluation, not a guard.  The
four measured cases above happen to evaluate full.  **So the pre-projector residual is a valid measure only on a full
evaluation**; R2.1 stays a `QCHEM_DIAGNOSTICS=group` report (`ReportImposedGroup`), and the throw waits for a measure that
holds on both routes: the D-side residual \|D − P D\| / \|D\| where the reduced route already forms P D (`FoldProjectedD`) --
a symmetric D is exactly the precondition the reduced scatter states.  OpenWork row "R2.1 as a throw".  (a), (b) and the
unmeasured raster route are OpenWork rows too.
**R2.4 (charge per reservoir)**: `RequireReservoirCharge` at the end of `tCompositeWF::FillOrbitals` -- each occupation
reservoir's filled electrons (k-weighted Σ occ) equal its prescribed N to 1e-10·max(N,1) (NaN throws); UT
`ReservoirCharge.R2_4_ALostElectronThrows`.  **R2.7 (no dead gate)**: `RequireLiveGates` at the iteration that declares
convergence -- a live [F,D] gate reading exactly 0.0 was never computed (it arrives from the accelerator's `GetError()`;
evidence G3) and throws; UT `SCFTrace.R2_7_ADeadGateThrows`.  **Not Δρ/N**: the first build also checked it and failed 13
correct one-electron atom tests (`A_*_E1`/`DE1`) -- a density-independent Hamiltonian reproduces its seed bit for bit, so an
exact 0 there is a verdict.  Value-based detection is the stopgap; the durable form is a gate that was not computed SAYING so
(an optional), which is the "[F,D] has no single home" remainder.
**R2.5 (stated decoration), 2026-10-30**: `RequireStatedMagneticOrdering` now refuses an IMPOSED polarized run without
`siteSpins` (+U or not; `greyImposition` exempt), and the seed-derived decoration in `SolidCalculation::BuildBasis` is GONE.  The
short tier was first run with an instrument in the derive branch: 10 tests derived, every one a sensible ordering (Si {0,0} ×6,
O₂ {+1,+1} ×2, Na {+1}, Na₂ / Mn₂ {+1,−1}); each now STATES exactly what it derived (anchors cannot move), and every shipped
deck already vets.  The plan's second clause (the seed's site moments match the stated signs at iteration 0) is NOT built.
**#10 (the polarized Si jump)**: the readout is `SolidCalculation::MakeHubbardMoments()` -- each +U manifold's Tr n↑ − Tr n↓
through the term's own projectors, hoisted from `ConvergeHubbardU`'s [U table] snapshot (one walk, both callers), because the
UTResponse Si fixture runs a UNIFORM XC mesh (no Becke basins, so `RunDiagnostics::OrderFinal` reads 0 by absence -- the first
draft of the case caught that).  RED as recorded: 0.191 e (n↑ 1.59 / n↓ 1.41), listed in `QCHEM_PANEL_RED`.

---

## 6. Questions that need the user's ruling

1. **Panel label mechanics.**  `_Panel` suffix → label `panel`, and the regular sweep becomes `ctest -j8 -LE 'long|panel'` (a
   CLAUDE.md change), with the full tier as `_Panel_Long`?  Or should fast panel cases stay in the regular sweep too (they are
   mostly existing tests) and only the full tier be excluded?
2. **Known-failing panel cases: `WILL_FAIL` (run every time, turn red when fixed) or `DISABLED_` (the current rule for a failing
   claim with an open row)?**  Recommendation: `WILL_FAIL` — a fix is then announced by the panel itself.
3. **GDM line search exhausted below the energy-noise floor** (OpenWork row): end the stage?  Recommendation: end it with a
   DISTINCT outcome (`NoiseFloor`, not `Converged`), which fails a test unless the gate was met or a `BoundedSCF` says why.
4. **Metals in `Auto`**: should `Auto` detect a metal (stage gap / kT) and add Kerker, or must a deck say `Metal`?  Unmeasured
   either way; B1 measures Pulay-D vs K+P on Al and, if affordable, SrVO₃ first.
5. **MOM under P1**: refuse MOM beside kT > 0?
6. **P3 scope**: C1 builds the contract types + the SCF producer; the grid producer is C2's (with its ε-sized raster).  Confirm.
7. **The Mn–O dimer target state**: CP2K OT (−30.34790) minimizes at fixed N without aufbau; is the panel's fixed point OUR crawl
   state (Mn 5.35) or CP2K's?  The panel needs one before #15 can turn green.
8. **Reject or keep**: I recommend dropping "aufbau holes" from P2 (it is a report, J5) and moving "reference eigenpair residual"
   from P2 to P3.  Everything else in P1–P5 is supported by the evidence, with the sequencing note on P1 (type-level
   unrepresentability waits for R4.1) and the added "three variants" trigger on P5.

### B3 as built (2026-10-30, lane A, on the lane-b pre-merge 3f22991c; NOT YET SWEPT -- sweeps with the twins after lane B lands)

- **R1.1** (`00cd6f98`): `ResolveMixing` throws (`RequireOneExtrapolator`) on Mixing::Pulay, or Auto with a hand-set
  `PulayDepth`, beside an accelerator that `ExtrapolatesFock()`.  Plain Kerker / linear α still compose with DIIS.  UT
  `MixingMode.R1_1_ADensityHistoryBesideAFockHistoryThrows`.  Static scan: no test or deck uses the pair.
- **V1.24(i) + D-DIIS-FDMIN** (`6240130e`): `FDMax` → `EngageBelowFD` (DIIS and GDM), `FDMin` deleted (never read); every
  positional `DIISParams`/`GDMParams` initialiser now designated.
- **R1.3** (`7a452e8d`, UTs `95d0b23e`): `qchem.HandOff` (qcCommon) -- `HandOffRequest` + the one line `[handoff] iteration N:
  <from> -> <to> (<why>)`.  Deciders that cannot act REPORT through their own face: `SCFAccelerator::TakeHandOff` (the ladder's
  advance / retreat; its own `***` prints are gone), `tDensityMixer::TakeHandOff` (Pulay priming → history, field and D-space;
  the polarized pair forwards its up leaf).  `tSCFIterator::HandOff` prints them and is the ONLY `Reseat` call (the GDM mixed
  fallback); the exact-resume Pulay-from-1 line goes through it; `SolidCalculation::BuildStage(why, …)` prints every stage
  boundary (anneal stages, the +U loop, the FD-LRT arms + restore).  Of the seven hand-offs, "restart → first iteration" keeps
  its `[restart] EXACT RESUME / WARM START` banner (a report of the restore, not a change of extrapolator).  Console-only:
  `SolidState`'s exact-resume test now reads the new line (`38496064`).
- **R1.4** (`02de08b2`): `KerkerPreconditioner` (FieldMixer) holds relax, G0, the step, the raster step and the G-space;
  `KerkerMixer` = (filter, no history), `PulayMixer` = (filter, Pulay history).  Bit-identical Etot sequences vs the B2 sweep
  on Si Kerker pol, Si +U Kerker+Pulay, O2Box Pulay-on-D, the SolidState Kerker+Pulay exact resume.
- **R1.5 MOVED TO C2** (controller ruling 2026-10-30, so B3 stays bit-identical).  The preconditioner OWNS the mixer's G-space
  (one constructor argument); the proposed owner (the density's own Coulomb/CD fit) is a different type (`FIT_CD_ABS`) from the
  face the mixer reads ρ̃ through (`FIT_SF_ABS`, `FourierDensity::GetFourierDensity`), and a different ball moves every
  Kerker/Pulay anchor -- C2 re-sizes that ball from the one accuracy target and re-judges the moved anchors ONCE, designing the
  face with the ball it serves.

### Rulings (user via the controller, 2026-10-29: "Everything looks good")

1. RULED 2026-10-29: fast panel cases (`_Panel`) STAY in the regular sweep; the full tier is `_Panel_Long`, which is `long`-labelled and so excluded from `-LE long`.  Both carry the label `panel` (`ctest -L panel`).
2. RULED 2026-10-29: a known-failing panel case is `WILL_FAIL` with its open tracker row, not `DISABLED_` (CLAUDE.md's DISABLED_ rule says so).
3. RULED 2026-10-29: an exhausted GDM line search below the noise floor ends the stage with a distinct `NoiseFloor` outcome, which fails a test unless the gate was met or a `BoundedSCF` says why (B6(b)).
4. RULED 2026-10-29 (pending a measurement): B1 measures Pulay-D vs K+P on Al (SrVO₃ if affordable) and brings it back; leaning: `Auto` detects a metal.
5. RULED 2026-10-29: MOM beside kT > 0 is refused.
6. RULED 2026-10-29 ("the split sounds very reasonable"): C1 builds the contract mechanism (need vs delivered, the pre-flight comparison) and the SCF producer, absorbing the response-reference Δρ refusal and the eigenpair-residual refusal; C2 wires the grid producer with its ε-sized raster.  B7 is unblocked.
7. RULED 2026-10-29 (provisional): compare our crawl state's moments and occupations with CP2K OT's first; provisionally CP2K's lower-energy state is #15's target; report the comparison before #15's budget is fixed.
8. RULED 2026-10-29: "aufbau holes" leave P2 (an end-of-stage report); "reference eigenpair residual" moves to P3.

### C1-DEG: degenerate open shells after D-FIT-MERGE (lane A, 2026-10-31; first measurements, the lane-b part still to run)

**Question (controller, 2026-10-31).**  With the fit bases' D shells back (lane B, f8e5848d), the molecular Hartree is anisotropic.  The degenerate open-shell PP atoms (Si singlet, O–Si quintet, the smeared O₂ singlet, panel #11) slosh at kT 1e-4 / 1e-3 under DIIS, GDM and Ladder, and lane B lands them at a stated kT 1e-2.  This section asks two things.  Is the slosh an occupation FEEDBACK loop (occupations SWAP; gain ~1/kT; a stiff fixed point exists), or a symmetric SADDLE (occupations settle EQUAL while [F,D] stalls; cure = a symmetry-breaking kick)?  And what is the finite-kT target?

**Step 1, measured on lane-a (pre-merge: d-less, SPHERICAL fit).**  B atom, dzvp, polarized LDA, multiplicity 2, 40 capped reruns read through `Calculation::GetLevels()` (probe kept out of the tree: scratch `boron_probe.C`):
- kT = 0: the three p↑ levels are degenerate to 1e-5 Ha at EVERY iteration, even with the electron in one p (occupation {1,0,0} by energy order).  So the d-less fitted Hartree + Vxc see only the spherical part of ρ.  The integer fill picks an arbitrary rotation inside the exactly degenerate triple each iteration, and E wanders over 8 mHa (−24.3393 … −24.3470) for 40 iterations without converging.  This is NEITHER a swap NOR a saddle: it is a zero mode of the fill (no restoring force selects the occupied p).
- kT = 1e-4, 1e-3, 1e-2: occupations 1/3, 1/3, 1/3 from iteration 1 (exact degeneracy ⇒ equal fractions at any kT); converges in 22 iterations, monotone.  The pre-merge "kT 1e-4 cures it" was this: a spherical problem made spherical.

**Step 1, lane-b binary (anisotropic Hartree), as-is test only.**  `P/A_PG_DFT.Energy/Z5` (B, Xalpha, kT 1e-4) CONVERGES in 24 iterations (lane-a: 43).  The level table shows the p↑ electron WHOLLY in one p (occupation 1/1).  The anisotropic Hartree splits the triple by ≫ kT, so smearing at 1e-4 is integer occupation in practice: a spontaneously broken symmetry that holds for B.  So B with ONE p electron does not reproduce the slosh.  The sloshing cases all put TWO electrons (or a hole pair) in a degenerate set: Si 3p², O₂ π*², O–Si.  **The per-iteration occupation trace on those, on lane-b code, is NEXT**, once D-FIT-MERGE is on main and merged here (lane-b's tree is lane B's; a separate build costs ~1 h of the build lock).  Probe: rerun the scratch probe for Si (sipp, gth-lda-valence, singlet) and the O₂ singlet at kT 1e-4 / 1e-3, recording the three (two) frontier occupations AND the splitting each iteration.  Read-out: occupations alternating between two orbitals with a split that flips sign = feedback (gain ~ 1/kT, since df/dε = −f(1−f)/kT); equal occupations with a split ≈ 0 while [F,D] stalls = saddle.

**Step 2 plan (MOM at kT = 0, stated occupations).**  Ruling 5 (MOM beside kT > 0 is refused) leaves MOM at kT 0 on a STATED state.  Each test names the state it pins:
- B {p_x↑};
- Si triplet p_x↑p_y↑;
- O–Si quintet;
- not "Si 3p²": closed p_z² and open p_x p_y differ in energy.

Judge each against PySCF at the SAME fixed occupation (`scripts/pyscf_molecule_oracle.py` has the fixed-spin Si triplet).  Use IMOM (Barca–Gilbert–Gill 2018): the reference is the INITIAL occupied set, so the run cannot drift back to aufbau, which plain MOM's moving reference allows.  Today's `MOMStartIter` delayed capture is a delayed IMOM, so the capture iteration must be the seed for a stated state (the seed already IS the stated occupation).

**Step 3: finite-kT candidates (reprints in `~/Code/reprints/MOM/`), measured on B first, then the step-1 cases:**
- **Ensemble direct minimisation** (Marzari–Vanderbilt–Payne PRL 79, 1337 (1997); Freysoldt–Boeck–Neugebauer PRB 79, 241103(R) (2009)).  Minimise the free energy over orbitals AND an occupation pseudo-Hamiltonian together.  The occupations move with the orbitals instead of being refilled from ε each step, so the 1/kT gain of the refill never enters the loop.  This removes the feedback loop by construction and is the principled kT > 0 target.  It is also the natural extension of our GDM: GDM's held Grassmann block becomes the orbital half, and a pseudo-Hamiltonian gradient adds the occupation half (SCFStrategyPlan: our parameterisation, not the method, is what bars GDM at kT > 0).
- **Occupation annealing / pseudo-fractional occupation** (Rabuck–Scuseria JCP 110, 695 (1999)).  Start at a high kT and lower it on a schedule.  This is our `schedule` annealed driver already (stages at falling kT).  It is a recipe, not a new method: cheap to try first, and it says whether a kT 1e-4 fixed point exists below the kT 1e-2 one.
- **DO-MOM vs smeared MOM** (Levi–Ivanov–Jónsson JCTC 16, 6968 (2020)).  Direct optimisation (exponential transformation) with MOM held on the occupied set converges the degenerate/excited cases where SCF+MOM oscillates.  That is GDM+MOM at kT 0, which is step 2's machinery with GDM as the engine.  Smeared MOM (`MOMSmearPenalty`) stays refused by ruling 5.

**Order.**  (a) the step-1 trace on lane-b code decides feedback vs saddle; (b) step 2 (IMOM, stated states, PySCF-judged), which is needed whatever (a) says, because a stated integer state is the honest gate; (c) annealing as a recipe on the kT 1e-2 cases; (d) ensemble DM only if (a) says feedback and (c) cannot reach kT ≤ 1e-3.

**C1-DEG correction (user via the controller, 2026-10-31): the trace must fill BY ENERGY across irreps.**  `Atom_EC` / `PseudoAtom_EC` (the `AtomCalculation` path) assign occupations BY IRREP: B's p electron sits in the m = −1 sub-irrep at every iteration.  That is a maximally stretched determinant, and it removes the degenerate-shell problem by construction, so it is no reproduction.  The trace runs through `Molecule_EC` (one ranked reservoir per spin spanning the point-group irreps: `spansSpatial`, `ranksIntegerFill`).  CHECKED: every kT 1e-2 case in f8e5848d is a molecular-path `Calculation` (`MoleculeSetup::Resolve`): `Si2_PP_U` si1, the O–Si quintet, the `GPW_SiBox` finite twins, the O₂ singlet.  The `AtomCalculation` cases there are unchanged and unsmeared.  The B result above (`A_PG_DFT`, `Calculation`) and the lane-a probe (`Calculation`; levels label "No Symmetry", i.e. C1, one irrep) are molecular-path too.

**Step 2 has two routes to a STATED state, chosen by whether the group separates the partners:**
- **Irrep-assigned occupation** (Atom_EC-style), where the imposed point group puts the degenerate partners in DIFFERENT irreps (an atom in D2h: p_x b3u, p_y b2u, p_z b1u; O₂ D2h: π*_x b2g, π*_y b3g).  The state is stated as counts per irrep: deterministic, no reference to drift, nothing to tune.  B {b3u↑}, Si triplet {b3u↑ b2u↑}, Si closed {b1u↑↓}, and the O₂ singlet in one π* irrep (the stretched determinant, ≠ the smeared even singlet: say which) are each a stated state.  This is the PREFERRED route; it needs the facade to take per-irrep counts on a Molecule_EC (today only the multiplicity is stated) and the group imposed.
- **MOM/IMOM**, only where symmetry does not separate them: a C1 run, a low-symmetry molecule, or a stated state INSIDE one irrep (two partners of one degenerate irrep, e.g. an atom's p set in a group where p is one 3-d irrep, if we ever impose O_h/I_h).
- **Measure both on one case** (B in D2h by irrep vs B in C1 by IMOM): same energy to the convergence floor, and the iteration counts compared.

**What a real-irrep stated state IS (controller, 2026-10-31).**  Route (a) states a DETERMINANT in real orbitals.  That is not always the spectroscopic term:
- **B {b3u↑}** is a pure ²P component: any one component of ²P is one determinant, and a superposition over m that differs in ONE spin-orbital collapses to one determinant.  NIST's spherical B is the M_L ENSEMBLE (1/3 each, the smeared fill), not a pure state, so judge {b3u↑} against PySCF at the same occupation, not against NIST.
- **Si closed {b1u↑↓}** and **the O₂ singlet in one π\* irrep** are NOT pure singlet terms in real orbitals.  π\*_x² is an equal ¹Δg/¹Σg⁺ mixture and lies ABOVE ¹Δg, and Si p_z² stands to ¹D the same way.  Both remain valid stated states: anchor each against PySCF at the SAME occupation, and **never label them ¹Δg / ¹D** in a test name or comment.
- **A pure ¹Δg (Λ = 2: (π\*₊₁)²) or Si ¹D M_L = 2 (p₊₁α p₊₁β)** IS one integer determinant, but only in angular-momentum-resolved (complex, m) irreps, the kind Atom_EC's m sub-irreps are.  The molecular D2h path cannot state it.  This is a LIMIT of route (a), recorded here, not a defect queued.
- **Genuinely multi-determinant states** (an M_S = 0 open-shell singlet, …) are out of scope for KS stated states.

**C1-DEG-DETECT: the degenerate-frontier detector (user request 2026-10-31; OpenWork row "No single-determinant solution", tagged X).**  User: the software should *"(eventually) have a robust way of detecting situations where they are asking for the real single det. solution where none exists and advise the user what to do (switch to complex basis set or specify occupations by irrep, or other?)"*.  The case is a graduate student who asks for an O₂ singlet or Si 3p² in real orbitals and gets a slosh, or a silent term mixture.  **The detector IS the step-1 trace instrument: build it with the trace; the advice and the report wiring come later.**
- **DETECT** (per iteration, from data the SCF already holds: the level list with ε, occupation, irrep, spin):
  - Group the levels of one spin into DEGENERATE SETS: the same irrep row set, |Δε| ≤ τ_deg.  τ_deg is relative to the set's width at the seed, never a fixed number: the anisotropic Hartree splits a set by the very fill under test.
  - (i) **Aufbau violation inside a set:** an occupied member above an empty partner, or an unequal integer fill of a set that straddles E_F.
  - (ii) **Swap counter:** the filled member of a set changes between iterations.  Count it from the occupied-subspace projection onto the set, not from orbital indices, which reorder freely inside a degenerate set; on lane-a's spherical fit this counter reads the zero mode as "fill rotates every iteration".
  - (iii) **Ensemble at convergence:** a smeared run ends with fractional occupations within kT of E_F on a degenerate set.  The answer is then an ENSEMBLE, and its energy is the ensemble's.
- **REPORT** (later): through the SCF's `Outcome`, as the failure reason or as a note on a converged run, never a throw.  It names the set (irreps, ε, fill) and what was asked (kT, symmetry, multiplicity), and fires BEFORE `NMaxIter` where (i)/(ii) persist for a set number of iterations.
- **ADVISE** (later), offering only the options the system has:
  - (a) state occupations by irrep, where the group separates the partners;
  - (b) complex m / Λ-resolved orbitals for a pure M_L / Λ component (¹Δg = (π\*₊₁)²); only the atom path has them today, so say so;
  - (c) accept the ensemble (smearing, equal fill) and know its energy is the ensemble's;
  - (d) IMOM in C1 on a stated state;
  - (e) a broken-symmetry seed.
  - Never label a real-orbital mixture as a pure term.
- **CONSUMERS:** the console, the run report json, and later the GUI's in-app agent (advice → a deck edit).
