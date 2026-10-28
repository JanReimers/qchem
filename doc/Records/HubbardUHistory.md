# HubbardUHistory — verbatim archive as of 2026-10-01 (RECORD, holds no open work)

Live plan: doc/HubbardUPlan.md. Below: the full HubbardUPlan.md as of 2026-10-01, then OpenWork §1.

---

# Self-consistent orbital-resolved (U, J) from our own on-site ERIs — the plan

**Born 2026-09-23**, at the top of `doc/` because it is being executed (`doc/OpenWork.md` §1 step 5,
increment 3 remainders).  Retires to `doc/Records/` the day the U-functional lands.  The north star above
it is `doc/BatteryMaterialsRoadmap.md`; the ruling under it is `doc/Pins.md` pin 23.

★ **What is new here, and worth saying once.**  Orbital-resolved DFT+U, determined self-consistently from
the code's *own* on-site ERIs, in a **Gaussian** basis, on transition-metal oxides, on a 14 GB desktop.
Each piece exists somewhere; the combination does not.  ACBN0 was built on plane waves and had to invent a
PAO-3G projection to get four-index on-site integrals at all — we own them natively (`Cache4`/Rk).  CP2K
has the Gaussians and hand-sets U.  QE has the U determination and is plane-wave.  Nobody has the
orbital-resolved, self-contained, Gaussian combination, and the reason to care is not novelty: it is that
a supercell-free U is the only kind that scales to a **composition sweep**, which is what a voltage curve
is.

---

## 0. START HERE — what to do next

The plan below is the ARGUMENT; this is the QUEUE.  **Two tracks, and they do not block each other** — A
decides whether the U functional is real, B decides whether the material runs at all.  Every action in the
plan appears here exactly once; if it is not in this list it is a finding, not a task.

**★ Track A's screening question is ANSWERED, and it answered NO.**  A1–A3 are done: route (b) [ACBN0 with
a screened bulk kernel] is refuted on real, matched-PP oracle data for both Ni-3d and O-2p (§4).  A4 does
not start.  **The routing decision is MADE (user, 2026-09-25): (d), our own finite-difference/DFPT linear
response — queued as A7, in-house and material-agnostic once built.**  It is queued BEHIND A6: broadening
the matched-PP hp.x oracle set across oxides/sulfides/fluorides first, both because it needs no new
capability (cheap, real progress now) and because it is real convergence-recipe practice this plan will
need repeatedly regardless of which route wins.  (A2's ABINIT cRPA route stays paused, superseded by A2b's
matched-PP hp.x route, which is what A3 was actually run against.)

### Track A — is the SCREENED route real?  (needs no spinels)
| | action | state |
|---|---|---|
| A1 | ~~Resolve the 11 % bare-\f$F^0\f$ discrepancy~~ | ✅ **RESOLVED 2026-09-23** (§6.1): NOT a bug — the two numbers are different quantities by definition, and comparing them was the error, not either computation |
| A2 | ~~ABINIT `ucrpa` on NiO (O 2p first), then MnO~~ | ⚠ **PAUSED 2026-09-25, redirected to A2b** — §7's Carta et al. read shows cRPA is unreliable in exactly NiO's hybridised-band regime (16× errors documented on their own materials), so the ≈1.2 eV dp-dp number is suspected of being that pathology, not a value to chase further via model-convention sweeps.  Ni-3d crash (§6.3) and MnO run both deprioritised with it |
| A2b | ~~Extend hp.x (matched-PP LRT) to O-2p~~ | ✅ **DONE 2026-09-25**: U(Ni 3d)=5.4343, U(O 2p)=8.5139 eV, both matched-PP (`IntegrationTests/QE/README.md`, `NiOgO.*`) |
| A3 | ~~Re-run gate3 scripts against A2b's real target~~ | ✅ **DONE 2026-09-25 — ROUTE (b) REFUTED.**  ω's 50.6 % apart (was 6 % on a wrong number), localization sign backwards again |
| A4 | ~~the screened kernel on `BareCoulombSource`, with \f$\varepsilon\f$ COMPUTED~~ | ⛔ **DO NOT START — A3 refuted, not held.**  §2's routing table now points at (c) or (d); which one is a plan-shape decision, not this row's to make |
| A5 | ABINIT `lruj` for **J** (hp.x gives none) — answers §5's first open question | optional, any time |
| **A6** | **Broaden the matched-PP hp.x oracle set** — get GOOD at converging oxides/sulfides/fluorides across the papers' own benchmark materials (SrVO₃, KCuF₃, Sr₂FeO₄, the LiMO₂ series M=V–Ni, TiO₂, ZnO, FeS₂ — §7; MnO's own O-2p redone with `ortho-atomic` instead of the fragile `atomic` projector already used once) | ⛔ **IN PROGRESS 2026-09-25.**  Element coverage checked: Ti/V/Fe/Co/Sr/K/F/Zn all convert cleanly on the first attempt (`gth2upf`); **Cu (3d10-4s1, needed for KCuF₃) hit a genuine aufbau limit cycle, fixed** with a kT-anneal-then-cold-MOM fallback (hot smear + damped mixing, then MOM-held re-converge) — lands on the correct 4s¹3d¹⁰ ground state, E=−47.926315 Ha.  Pure fallback, zero effect on elements that already converge (Mn/Ni/O/Li regression-checked identical).  **SrVO₃ (material queue item 1) DONE 2026-09-25: U(V 3d) = 6.2502 eV, matched-PP LRT** (`IntegrationTests/QE/README.md` §A6) — found and routed around a real `gth2upf` limitation on the way: Sr q10 (semicore 4s4p5s) silently integrates to 2 electrons not 10 (the pseudo-atom EC caps one shell per `l`, can't hold two s-shells at once); used Sr q2 instead, the same light-valence convention already used for Mn/Ni.  **KCuF₃/Sr₂FeO₄ (queue items 2–3) DONE 2026-09-28/29: U(Cu 3d)=8.1629 eV, U(Fe 3d)=8.0116 eV** — Dr. Carta replied with his actual input files, removing all sourcing uncertainty; see the row's own detail below and `IntegrationTests/QE/README.md` §A6.  **LiCoO₂ (from queue item 4, the LiMO₂ series) DONE 2026-09-26: U(Co 3d) = 7.3070 eV**, matched-PP LRT — picked ahead of V/Cr to shake out the recipe on a well-characterized member first.  Structure sourced and bond-length-verified from Pinsard-Gaudart et al. 2011 (a real citation, not memory).  Found a genuine `hp.x` requirement along the way: Hubbard atom(s) must be listed FIRST in `ATOMIC_POSITIONS` or it refuses to run.  Also found LiCoO₂ needs the SAME 2-step (smeared → fixed-occupation) recipe as MnO/NiO despite being NONmagnetic — the 2-step is about the GAP, not about magnetism, correcting this plan's earlier "magnetic insulator" framing.  **LiVO₂ DONE 2026-09-26: U(V 3d) = 5.9526 eV** — geometry from a real citation (user-supplied, Mat. Res. Bull. 27, 555 (1992)) after a self-relaxed (QE `vc-relax`, our own PPs) stand-in was tried and then DISCARDED the moment the citation existed.  Idealized (untrimerized) LiVO₂ comes out METALLIC under nonmagnetic LDA — SrVO₃'s recipe (no 2-step), not LiCoO₂'s.  Full detail: `IntegrationTests/QE/README.md` §A6.  **LiCrO₂ DONE 2026-09-26: U(Cr 3d) = 5.8111 eV** — the flagged aufbau risk did NOT materialize (`gth2upf` converged Cr q6 to 3d⁵4s¹ cleanly, no anneal fallback needed); structure from Garg et al., *Crystals* 9(1), 2 (2019), bond-length verified (Cr–O 2.002 Å, Li–O 2.114 Å vs the paper's 2.003/2.113 Å); also METALLIC under nonmagnetic LDA, same pattern as V.  **LiFeO₂ DONE 2026-09-27: U(Fe 3d) = 7.5915 eV** — real LiFeO₂ isn't R-3m (correctly flagged from memory before sourcing), idealized into the row's template anyway, valid for validating our own methodology; geometry from Materials Project mp-19419, bond-length-verified; also metallic, not a repeat of MnO's d⁵ trap (that was spin-polarized-response-specific).  **LiNiO₂ DONE 2026-09-27: U(Ni 3d) = 9.1730 eV — LiMO₂ ROW COMPLETE (V, Cr, Fe, Co, Ni)** — real LiNiO₂ is
also JT-active (confirmed via Materials Project's own DFT+U-relaxed entry coming back genuinely monoclinic,
not idealized), same idealize-into-the-row decision as Fe; geometry from Seo et al., *J. Electrochem. Soc.*
165 (2018) A2554, bond-length verified; metallic, largest U of the row.  **First use of the new checkpoint
practice** (`doc/OpenWork.md` SCF-restart row): archived the converged U_in≈0 `.save` to
`IntegrationTests/QE/checkpoints/linio2_U0.save/` before `hp.x` touched it, so a later self-consistent-U
rerun can warm-start instead of reconverging from scratch — the five earlier row members do not have this.
Row summary table + all detail: `IntegrationTests/QE/README.md` §A6.  **TiO₂ (rutile) DONE 2026-09-27: U(Ti 3d) = 4.6368 eV** — a genuinely different chemistry (Ti⁴⁺ d⁰, no partially-filled shell at all); geometry from QE's own `PP/examples/example08` (authored by Iurii Timrov, one of the `hp.x` method's own authors), cross-checked against CP2K's independently-sourced Wyckoff-cited rutile cell (<0.01% agreement); no 2-step needed (`occupations='fixed'`, a genuine d⁰ gap, no occupation ambiguity).  **ZnO DONE 2026-09-27 — NOT A USABLE ORACLE: U(Zn 3d) = 35.3358 eV**, an order of magnitude outlier.  Zn²⁺ d¹⁰ is a genuinely CLOSED shell; χ₀(Zn,Zn)=−0.00347 → χ(Zn,Zn)=−0.00309, both tiny and nearly equal — the SAME closed-shell pathology already named for MnO's d⁵ (χ₀⁻¹−χ⁻¹ nearly cancels), reached by a different mechanism (genuine shell closure, not spin-cancellation).  Flag ZnO/Zn-3d alongside MnO/Mn-3d as ill-posed for same-site LRT U; do not quote 35 eV as a value.  Ground state itself is fine, checkpoint archived.  Next: FeS₂ (user has its citation: pyrite, cubic Pa-3̄, a=5.418 Å, from `ct3c01403.pdf`).  **KCuF₃ DONE 2026-09-28: U(Cu 3d) = 8.1629 eV** — Dr. Carta replied with his ACTUAL input files (`~/Code/reprints/materials_cloud_submission/`); used his exact geometry (a=4.066704097 Å cubic Pm-3̄m, the SAME template as SrVO₃) but OUR OWN recipe (his Hubbard-parameter determination is MLWF-based, still not directly comparable to our ortho-atomic `hp.x` route even with matched geometry).  His deck also corrected the SI's stated dual convention (his ecutrho/ecutwfc=8, not the SI's "four times").  Cu q11 needed a real cutoff scan (450 Ry, not the session's usual 300) and exercised the aufbau fallback again, bit-for-bit reproducing the logged E_atom.  Metallic, no 2-step.  **Sr₂FeO₄ DONE 2026-09-28/29: U(Fe 3d) = 8.0116 eV** — K₂NiF₄-type body-centered tetragonal, geometry given in exact Cartesian form in his deck (bond-length verified).  `hp.x` CRASHED on his raw relaxed cell (`D_S (l=2) ... not orthogonal` — the three cell vectors match to only ~8 figures, a real relaxation tolerance, and QE's symmetry-finder needs exact invariance for the l=2 Wigner rotation build); fixed by symmetrizing the cell (averaging the vectors' common components, same atomic fractions, identical total energy to the last printed digit after the fix — confirms it was pure noise removal).  **General lesson for any future real-relaxed-cell input**: symmetrize before debugging anything else.  Longest run this session (11h4m, the larger lower-symmetry 7-atom cell).  Close to LiFeO₂'s 7.59 eV despite different oxidation state (Fe⁴⁺ d⁴ vs Fe³⁺ d⁵) and host — sane cross-check.  BOTH of Dr. Carta's materials now done; KCuF₃/Sr₂FeO₄ queue items CLOSED.  Full detail: `IntegrationTests/QE/README.md` §A6.  Remaining A6 queue: FeS₂ (pyrite, cubic Pa-3̄, a=5.418 Å, citation in hand) |
| **A7** | **Build route (d): our own finite-difference/DFPT linear response**, material-agnostic (§5's open question; CP2K `qs_linres_*` as the Gaussian-basis reference) | ⛔ **QUEUED, AFTER A6** (user, 2026-09-25).  Does NOT reopen route (b) — a single uniform screening length is refuted regardless of how ε is computed; DFPT gives the general, manifold-resolved case (b) was a crude shortcut for, and SUPERSEDES it (once built, U comes from the response directly, no bare-tensor-times-kernel construction).  DIP-based estimator strategy, `doc/OpenWork.md` §2 row.  **DESIGN RULED 2026-09-27, STAGE R0 DONE + VALIDATED 2026-09-28, STAGE R1 (molecular CPHF, the first KERNEL; HF α == PySCF to 1e-6) DONE 2026-09-28: `doc/LinearResponsePlan.md` (read its ▶ START HERE)** (interface changes H1–H5/C1–C2/E1/S1/P1/M1, stages R0–R4, D1–D6 all ruled; extends to PW-Sternheimer, SOC, forces, §4b; U self-consistency = outer loop only, §4c).  **User intends to plan/scope/SOLID-design A7 in a DEDICATED session** (2026-09-26) — read "A7 scoping insights" just below the material queue FIRST: six concrete things this session's `hp.x` runs surfaced (architecture decomposition, occupation-policy reuse, k/q commensurability, projector consistency, a same-site-vs-intersite scoping fork, real cost data) |

### Track B — can we RUN the material?  (needs no oracle)
| | action | state |
|---|---|---|
| B1 | ~~Three spinel structures into `materials.json`~~ | ✅ **2/3 DONE 2026-09-25**: λ-MnO₂ (12 atoms) and LiMn₂O₄ (14) landed, ase-generated + bond-length-verified (§4 prerequisites).  Li₂Mn₂O₄ (16, JT-tetragonal) deliberately deferred — gate 1 doesn't need it |
| B2 | ~~Gate 1 — magnetic robustness~~ | ✅ **PRELIMINARY VERDICT 2026-09-25: order SURVIVED in all 6 arms** (λ-MnO₂ and LiMn₂O₄ × U=0,2,4 eV) — see §4 |
| **B3** | **Gate 2 — run sizing**: one converged SCF per composition at Γ, wall + peak RSS logged | after B2 |
| B4 | `valence_semicore.bsd` + a `BasisSetData` enum value (§6.3) **landing together with the basis-vs-PP check** (§6.2), then the q1-vs-q3 discriminator on the Li intercalation energy | ⚠ Li **q1 is already committed and working** — B1–B3 do NOT wait on this |

### Track C — infrastructure both tracks eventually need
| | action | state |
|---|---|---|
| C1 | **Gate 4's measurement**: one 32-atom MnO 2×2×2 supercell SCF at Γ, wall + peak RSS — says whether the supercell route exists for us at all | not started; cheap; settles an assumption the tracker has carried untested |
| C2 | The **1.02 mHa shifted-MP fold defect** (`doc/OpenWork.md` §4) — *"fix before KP-1"*; any multi-k U inherits it | prerequisite for C3 |
| C3 | **k-parallelism** (`doc/OpenWork.md` §2 row KP) — parked "until we are suitably embarrassed"; this plan is the embarrassment | after C2 |

**Oct 6–20 (user away) = the long unattended runs**, and only what B2/B3/C1 have justified: the three
compositions × the outer loop, the k-mesh arms, and the calibration LR run if C1 says the supercell fits.
⛔ `scripts/memsafe -p`, never bare.

### A6's material queue (session handoff, 2026-09-25)

★ **This is A6, not A2b** — A2b was NiO O-2p specifically and is done/closed.  A6 is the broader
"get good at converging matched-PP hp.x oracles across chemistries" item, and this is its ordered
material list for a fresh session to pick up.  Element coverage for every element below is already
checked (`gth2upf`, this session) — Ti/V/Fe/Co/Sr/K/F/Zn convert cleanly; Cu needed and got a kT-anneal
fallback fix (aufbau limit cycle, 3d10-4s1 vs 3d9-4s2 near-degeneracy — see A6's row above and the
`gth2upf.C` commit).  **No element in this list is blocked.**

★ **Reference table (started 2026-09-30, user: "there should be lots of reference numbers... we should
probably start a table").**  Every `hp.x` matched-PP LRT value found so far, pulled out of the prose below
into one place.  All are U (or U(eff)=U−J) in eV, matched pseudopotential/projector (our GTH, converted with
`gth2upf`) unless flagged otherwise; **U_in is the ground state the response is linearised about** — 0 eV
unless noted (`doc/Pins.md`: "a response is a function of the state it linearises about").  This is U_0
only — no U_SC (self-consistent) column yet, since that needs A7 stage R4's outer loop (not built); add it
alongside U_0 the day a material has both, per NOTES' own suggested shape (state at both ends: gap + moment
at U=0 and at U_SC too, not just the two U numbers).

**Checked directly against the decks (2026-09-30), not assumed:**
- **Symmetry is uniform, not a mix**: no deck sets `nosym`/`noinv`, and every checkpoint XML sampled
  (`kcuf3`, `linio2`, `sr2feo4`, `tio2`, `zno`) confirms `<nosym>false</nosym> <noinv>false</noinv>` — QE's
  default full space-group detection/reduction is ON in every run in this table (a different axis from our
  own code's FREE-vs-Shubnikov-imposed choice, but constant across the whole table either way).
- **The 2-step recipe is smeared → fixed, not literally kT>0 → kT=0**: step 1 sets `occupations='smearing'`
  with `degauss` > 0 (Ry) — Gaussian for MnO/NiO/SrVO₃, Marzari-Vanderbilt `'cold'` for the LiMO₂ row +
  KCuF₃ + Sr₂FeO₄ — a NUMERICAL broadening, not a physical electronic temperature, though it plays the same
  role (fractional occupations so SCF can find the right band ordering before a gap opens).  Step 2 reruns
  with `occupations='fixed'` and no `degauss` at all (true integer occupations, `nbnd` carried over from
  step 1).  Materials marked "no 2-step" below have only ONE smeared deck (real metals, nothing to land
  on); TiO₂/ZnO go straight to `occupations='fixed'` with no smearing step (clean d⁰/closed-shell gap,
  nothing to help).
- **`k`-grid (the SCF's own BZ sampling) is NOT uniform, but the `q`-grid (hp.x's linear-response
  perturbation mesh — the thing the README's own warning calls "the supercell size of the method") IS**:
  every `*.hp.in` sets `nq1=nq2=nq3=2` with no exception.  The underlying `pw.x` `K_POINTS` differs: MnO/NiO
  (A2b, earlier) ran 2×2×2; every A6 material from SrVO₃ onward ran 4×4×4 (finer, and arguably the more
  consequential choice for the actual metals in the LiMO₂/KCuF₃/Sr₂FeO₄ row, where k-sampling matters more
  than for an insulator).
- **Hubbard projector: `ortho-atomic` everywhere except MnO, which was `atomic`** — this is the input-deck
  confirmation of the MnO O-2p row's "fragile projector" flag, and ✅ **REDONE 2026-09-30** (user: "does it
  make sense to re-run MnO with ortho-atomic... at least then table is consistent") — see the row and the
  physics discussion right after the table.

| material | manifold | U_LRT (eV) | U_in | k-grid / q-grid | projector | ground state | citation (geometry) | status |
|---|---|---|---|---|---|---|---|---|
| NiO | Ni 3d | **5.2670** | 3 eV | 2³/2³ | ortho-atomic | AFM-II insulator, 2-step | `NiOgO.*` (A2b) | the one SOLID gate-3 oracle point (§4) |
| NiO | Ni 3d | 5.4343 | 0 | 2³/2³ | ortho-atomic | AFM-II insulator, 2-step | `NiOgO.*` (A2b) | different Hubbard-channel-set convention than the row above (Q10 effect) — NOT the same number, do not average |
| NiO | O 2p | 8.5139 | 0 | 2³/2³ | ortho-atomic | AFM-II insulator, 2-step | `NiOgO.*` (A2b) | matched-PP, replaces an earlier literature *bound* (cRPA ≳4 eV) |
| MnO | Mn 3d | 0.198 (atomic) / **0.9856** (ortho-atomic) | 0 | 2³/2³ | atomic → ortho-atomic | AFM-II insulator, 2-step | `mno.hp.in` / `mno_oa.hp.in` | **STILL EFFECTIVELY NOT USABLE, and now PROVEN not a projector artifact**: χ/χ₀ = 0.955 under ortho-atomic (only 4.5% screening) vs 0.981 under atomic (1.9%) — the number moved 5× but the pathology (weakest screening of any manifold in this table by a wide margin) did not go away under a completely different projector.  See physics discussion below the table |
| MnO | O 2p | ~~26.56~~ → **11.1503** | 0 | 2³/2³ | atomic → **ortho-atomic** | AFM-II insulator, 2-step | `mno_oa.hp.in` (this session) | ✅ **REDONE 2026-09-30, fragile-projector flag CONFIRMED REAL**: the number changed by more than 2× and landed in-family with the rest of the table (cf. NiO O-2p 8.51 eV) — unlike Mn-3d, O-2p's earlier number really was mostly a projector artifact, not physics.  Logs: `mno_oa.{scf,scf2,hp}.out` (not committed, regenerate from the decks) |
| SrVO₃ | V 3d | **6.2502** | 0 | 4³/2³ | ortho-atomic | metal, no 2-step | ABINIT `tucalc_crpa_1.abi` cell | χ/χ₀≈0.081 (far more screened than NiO/MnO — a real metal) |
| KCuF₃ | Cu 3d | **8.1629** | 0 | 4³/2³ | ortho-atomic | metal, no 2-step | Carta et al. (author's own input files) | our ortho-atomic projector vs their MLWF — not the same convention even with matched geometry |
| Sr₂FeO₄ | Fe 3d | **8.0116** | 0 | 4³/2³ | ortho-atomic | metal, no 2-step | Carta et al. (author's own input files) | close to LiFeO₂'s 7.59 despite Fe⁴⁺ d⁴ vs Fe³⁺ d⁵ — a sane cross-check, not identical chemistry |
| LiCoO₂ | Co 3d | **7.3070** | 0 | 4³/2³ | ortho-atomic | insulator (1.56 eV gap), 2-step | Pinsard-Gaudart et al. 2011 | 2-step needed despite being NONmagnetic — it's about the gap, not magnetism |
| LiVO₂ | V 3d | **5.9526** | 0 | 4³/2³ | ortho-atomic | metal (idealized, untrimerized), no 2-step | Mat. Res. Bull. 27, 555 (1992) | |
| LiCrO₂ | Cr 3d | **5.8111** | 0 | 4³/2³ | ortho-atomic | metal, no 2-step | Garg et al., *Crystals* 9(1), 2 (2019) | |
| LiFeO₂ | Fe 3d | **7.5915** | 0 | 4³/2³ | ortho-atomic | metal, no 2-step | Materials Project mp-19419 | real LiFeO₂ isn't R-3m; idealized into the row's template (user decision) |
| LiNiO₂ | Ni 3d | **9.1730** | 0 | 4³/2³ | ortho-atomic | metal, no 2-step | Seo et al., *JES* 165 (2018) A2554 | largest U of the LiMO₂ row; U_in≈0 checkpoint archived (`checkpoints/linio2_U0.save/`) for a future U_SC warm start |
| TiO₂ (rutile) | Ti 3d | **4.6368** | 0 | 4³/2³ | ortho-atomic | insulator (d⁰ gap), no 2-step | QE `PP/examples/example08` (Timrov) | cross-checked vs CP2K's independently-sourced cell to <0.01% |
| ZnO | Zn 3d | 35.3358 (outlier) | 0 | 4³/2³ | ortho-atomic | insulator, no 2-step | — | **NOT A USABLE ORACLE**: Zn²⁺ d¹⁰ closed shell, χ₀≈χ≈0 by a different mechanism than MnO's (shell closure, not spin-cancellation) — do not quote 35 eV as a value |
| Si (CK-alpha, not hp.x) | Si 3p | *n/a — a χ cross-check, not a U* | 0 | Γ/Γ (CP2K, not QE) | n/a (LOWDIN) | insulator | `si_ckalpha_a{0,p,m}.inp` (this session) | χ_CP2K −14.5164 Ha⁻¹ vs our χ_FD/χ_LR −14.5117 Ha⁻¹ (0.033%) — validates the METHOD (our DFPT/FD), not a material U; `doc/OpenWork.md` "CK-alpha" |
| FeS₂ (pyrite) | Fe 3d | **7.5351** | 0 | 4³/2³ | ortho-atomic | nonmagnetic (LS Fe²⁺ d⁶), no 2-step | Macke et al., `ct3c01403.pdf` (their own cited experimental structure) | identical to 4 s.f. across all 4 symmetry-equivalent Fe sites; ~5h run, the longest `hp.x` this session — see A6 queue note for a driver-bug recovery caveat (number is solid, original console log is not preserved) |

✅ **FeS₂ DONE (2026-09-30, batch queue, `queue/02_fes2_hpx.sh`)**: structure pulled directly from
`ct3c01403.pdf` (Macke et al., the orbital-resolved-U paper itself — cubic Pa-3̄, a=5.418 Å, the paper's own
cited EXPERIMENTAL structure, not their PBE-relaxed one, since qchem6 has no GGA/relaxation yet to reproduce
that); S 8c coordinates hand-derived from Pa-3̄ symmetry at x=0.386, cross-checked against CP2K's own
`c_1_FeS2.inp` regression deck to <0.003 (that deck's coordinates are intentionally perturbed for a
symmetry-finder test).  **S converts cleanly via `gth2upf`** (E_atom −10.069183 Ha) — the "S untested"
open item from A6's scope is resolved.  SCF (nonmagnetic, `occupations='fixed'`, no 2-step) **converged in
39 iterations**, E = −330.04024036 Ry — confirms the diamagnetic low-spin Fe²⁺ d⁶ (t₂g⁶eₘg⁰) prediction, no
smearing needed.  `hp.x` (shell-averaged Fe-3d, `ortho-atomic`, q 2×2×2) took ~5 hours (the longest `hp.x`
run this session — lower symmetry, 12-atom cell, harder GMRES convergence than MnO's 4-atom cell) and gave
**U(Fe 3d) = 7.5351 eV**, IDENTICAL to 4 significant figures across all 4 symmetry-equivalent Fe sites — a
strong internal consistency check (Pa-3̄ really does treat them all the same, and the response converged
cleanly, not to 4 different noise floors).
⚠ **Recovery note**: a bug in the batch driver's halt/resume logic (fixed same session,
`batch/run_queue.sh`) raced this job into an unnecessary ~full rerun right as the original finished — the
ORIGINAL run's own `fes2.Hubbard_parameters.dat` survived (backed up before anything could touch it,
`batch/recovered/`) and is the number quoted above, but its `.scf.out`/`.hp.out` console logs were
overwritten by the aborted duplicate before recovery, so the full iteration-by-iteration log for THIS
specific run is not preserved (the deck + the Hubbard_parameters.dat table are; that is enough to reproduce
it, just not to re-read the original convergence trace).  The redundant rerun was killed once recognized, no
data lost beyond that log.

Remaining A6 queue after FeS₂: none named (the paper's own material set — pyrite + β-MnO₂ — is now both
started; β-MnO₂ was never in our A6 queue, it's Macke et al.'s SECOND material, not previously scoped here).

The independent-oracle *ratio* table (ours vs hp.x/cRPA,
the thing gate 3's screening test actually consumes) is separate and already exists at §4's "manifold / ours /
independent oracle / ratio" table below — this table is the raw hp.x values feeding it, not a replacement.

★ **How `Hubbard_projectors` maps onto "pick atoms" vs "pick bands in an energy window" (user question,
2026-09-30).**  QE's `Hubbard_projectors` (`~/Code/q-e/PW/src/ldaU.f90:134`) takes exactly three values:
`'atomic'`, `'ortho-atomic'`, `'file'`.  **Both `atomic` and `ortho-atomic` are the SAME choice — pick
atoms** — they differ only in whether the chosen pseudo-atomic orbitals (the PP's `PP_CHI` radial functions,
the same object our own `gth2upf` writes from the GTH database) are used as-is (`atomic`: non-orthogonal,
can overlap between neighbouring atoms or between shells on one atom — the "fragile" choice, MnO's O-2p
row) or Löwdin-symmetrically-orthogonalized within the manifold first (`ortho-atomic`: \f$O^{-1/2}\f$
applied per `force_hub.f90`/`stres_hub.f90` — an honest orthonormal projector, everyone else's default,
and the convention our own code's LOWDIN manifold already matches, per the `[+U]` console banner).
**"Pick bands within an energy window" is the OTHER option, `Hubbard_projectors='file'`**: QE does not build
that projector itself — `'file'` loads externally-constructed ones, and the standard source is a
disentangled, energy-window-selected, maximally-localized Wannier function from Wannier90.  That is exactly
what KCuF₃'s row above already flags: **Carta et al.'s Hubbard projector is an MLWF, not an atomic
orbital** — the row you're asking about is a live example of option (2), sitting right next to our own
option-(1) `ortho-atomic` number for the same material, which is why the two U(Cu 3d) values are not
expected to agree even with the geometry matched exactly.

★ **Is the χ₀≈χ "closed-shell problem" a finite-difference/LRT artifact DFPT would fix, or real physics?
(user question, 2026-09-30, prompted by the MnO/ZnO rows above.)**  REAL PHYSICS, not a solver artifact —
and the ortho-atomic rerun above is direct evidence either way (§ below).
- **Not an artifact**: finite-difference LRT and DFPT (Sternheimer) compute the IDENTICAL quantity, χ =
  dn/dα of the SCF ground state — DFPT is just an analytic α→0 solve of the same derivative a finite step
  estimates.  This is not asserted, it is MEASURED: this session's CK-alpha result (above, Si) has our own
  DFPT (`chi_LR`) and an independent finite-difference (`chi_FD`, and now CP2K's own separate finite-
  difference implementation) agreeing to 3e-4 relative.  A different solver cannot return a different answer
  for a well-defined derivative it is computing exactly.
- **Why the response is genuinely small for a high-spin d⁵ (or closed d¹⁰) shell**: Dudarev's α shifts the
  on-site potential UNIFORMLY across the whole manifold, both spins.  For Mn-3d in MnO, the MAJORITY channel
  is already at occupation 1 in every orbital — Pauli-saturated, nowhere to put more charge, dn/dα≈0 by
  construction.  The MINORITY channel sits deep in the exchange gap (Δ_ex ~ several eV) — a 1e-3 Ha (27 meV)
  shift is nowhere near enough to pull a state across it, so its KS-polarizability energy denominator
  (~1/Δ_ex) suppresses its contribution too.  With BOTH channels inert, there is very little to redistribute
  at any level, bare or screened — χ₀ and χ are each small AND close together not because screening happens
  to cancel them, but because there is barely a response for screening to act on.  U = χ₀⁻¹−χ⁻¹ then divides
  by the difference of two near-equal small numbers: ill-conditioned by construction.  ZnO's d¹⁰ closed
  shell reaches the same place by a different route (no partial occupation anywhere in the manifold at all).
- **The reassuring part**: this is roughly where DFT+U's own physical motivation is weakest — U corrects
  delocalization/self-interaction error in fractionally- or near-degenerately-occupied states, and a rigidly
  filled-majority/empty-minority shell has comparatively little of THAT error to correct via THIS specific
  same-site reoccupation channel.  It does not mean MnO needs no U (LDA famously gets its gap wrong without
  one) — it means the error U fixes shows up more through Mn–O hybridization/charge-transfer than through
  "how far does Mn-d's own occupation move when you poke Mn-d directly," which is exactly why cRPA (a
  differently-probed, methodologically independent quantity — see below) is the queued fix for this row, not
  a reason to distrust DFT+U itself.
- ✅ **DIRECT EVIDENCE, same session**: the ortho-atomic rerun above is effectively a same-code,
  different-projector test of this claim.  χ/χ₀ went from 0.981 (atomic) to 0.955 (ortho-atomic) for Mn-3d —
  moved, but stayed the weakest screening of any manifold in the whole table by a wide margin (everywhere
  else is 0.08–0.6) — while O-2p's number changed by more than 2×.  **A projector swap that leaves one
  manifold's pathology essentially intact while correcting another's is exactly what "real physics in one
  case, projector artifact in the other" looks like.**
- ✅ **RUN 2026-09-30 — RESULT DOES NOT CANCEL, and this is a genuinely open question, not a closed one.**
  Free AFM-II MnO, Γ, complex, U_in=0, site0 Mn-3d perturbed (`gpwprobe mno`, the banked Kerker+Pulay+Null
  recipe, `GPW_OMP_THREADS=14`) — first attempt (default Ladder/GDM/MOM accelerator) did NOT converge in 200
  iterations; this recipe converged cleanly in 18–29 iterations for all four SCFs (ground state + two ±α +
  the restore).  Self-consistent LR: χ₀ = −4.6021304, χ = −3.1905158 (**χ/χ₀ = 0.693, 31% screening**), U =
  2.61606 eV.  Finite-difference cross-check on the SAME manifold: χ_FD = −3.2259878 — **1.1% from the LR
  value** (a real internal cross-check, looser than Si's 2.7e-6 but on a much harder, magnetic, multi-SCF
  case).  **This is NOT the near-total cancellation hp.x shows** (0.955–0.981 under either projector) — it
  sits in the same "healthy screening" range as every other manifold in the table (0.08–0.6).
  ⚠ **NOT YET AN APPLES-TO-APPLES COMPARISON**: this run used `MNO_U_RADIAL`'s DEFAULT, `"every"` — CP2K's
  every-shell manifold (every d-type primitive shell on the site, un-orthogonalized against each other),
  which is a BROADER object than hp.x's single contracted, Löwdin-orthogonalized 3d orbital (`ortho-atomic`)
  or the even more minimal `atomic` one.  A broader manifold spanning multiple radial shells may simply have
  more room to redistribute charge among ITSELF under the α shift than one tightly-defined atomic orbital
  does — which would make this a manifold-DEFINITION effect, not evidence against the Pauli-
  saturation/exchange-gap argument above.  Log: `/home/janr/Code/qchem6-runs/MnO/mno_afm2_free_U0_radEvery_chi_20260930.{log,cmd}`.
  ✅ **MATCHED-MANIFOLD RERUN, same session: `MNO_U_RADIAL=ortho` (the two TM 3d sets Löwdin-orthogonalised
  against each other — our closest analogue to hp.x's `ortho-atomic`).**  Converged the same way (18–29
  iters).  χ₀ = −4.3056558, χ = −3.603541 (**χ/χ₀ = 0.837, 16.3% screening**), U = 1.23137 eV.  FD
  cross-check: χ_FD = −3.631901 — 0.78% from LR (tighter than the `every` run).  **This is a real,
  reproducible, MONOTONIC trend across four points, not noise:**

  | manifold definition | χ/χ₀ | screening | U (eV) |
  |---|---|---|---|
  | ours, `every` (broad, multi-shell) | 0.693 | 31% | 2.62 |
  | ours, `ortho` (single Löwdin orbital, matched to hp.x) | 0.837 | 16% | 1.23 |
  | hp.x `ortho-atomic` | 0.955 | 4.5% | 0.99 |
  | hp.x `atomic` (non-orthogonalized) | 0.981 | 1.9% | 0.20 |

  **Narrowing the manifold toward a single atomic-like orbital moves the ratio monotonically toward hp.x's
  near-total cancellation.**  This REFINES, rather than refutes, the Pauli-saturation/exchange-gap argument
  above: it isn't that the d⁵ configuration has literally zero accessible response at any resolution — it's
  that a manifold built from ONE tightly-localized orbital per site has the least internal freedom to
  redistribute charge among itself under the α shift, so it saturates hardest; a manifold spanning more
  radial freedom (still l=2, still on the same site) has more room and screens more like a normal shell.  Both
  measurements are "real" — they are honest dn/dα of two DIFFERENTLY-DEFINED manifolds, not one right answer
  and one wrong one.  **The remaining gap (0.837 vs 0.955–0.981) is still open** — plausibly basis/projector
  detail (our `ortho` still isn't bit-identical to QE's `ortho-atomic`: different pseudopotential radial
  functions, different contraction), not yet chased further.  State checkpoints (CK-1 HDF5, warm-startable via
  `SolidCalculation::Restart`) + exact command lines + logs for both runs: `/home/janr/Code/qchem6-runs/MnO/`.

  ✅ **A FIFTH point, from a THIRD code: CP2K CK-alpha on MnO Mn-3d (2026-09-30, batch queue).**  Same recipe
  as the Si validation (`doc/OpenWork.md` "CK-alpha") — `cp2k_ckalpha.ssmp`, native PADE, U_MINUS_J=0,
  `ALPHA [hartree] ±1.0E-3` on the Mn1 `&KIND` only, `PLUS_U_METHOD LOWDIN`, same AFM-II UKS broken-symmetry
  cell/basis (`VALENCE-LOWQ-VA`) as the banked `mno_afm2_gpw_va.inp` oracle.  All three SCFs converged (94,
  96, 104 steps).  **χ = −1.8938 Ha⁻¹** (UKS: `fspin=1.0`, so `trq = E_DFT+U/α` needs no ×2 correction, unlike
  the RKS Si case — `src/dft_plus_u.F:317-319`).
  ⚠ **This is χ ONLY, not a χ/χ₀ ratio** — the two-point ±α finite-difference recipe gives the SCREENED
  (self-consistent) response directly; unlike `hp.x`'s DFPT, it has no companion χ₀ (bare/frozen-potential)
  output built in.
  ⛔ **TRIED 2026-09-30, DEAD END: `SCF_GUESS RESTART` + `MAX_SCF 1` from the converged α=0 state does NOT
  give χ₀.**  Both ±α runs returned IDENTICAL `trq` to all printed digits (5.5181333169) regardless of sign
  — the telltale sign of an artifact, not a response.  Diagnosis: CP2K's step-1 "Total energy" summary block
  reports the DFT+U energy using the occupation matrix of the INPUT (restarted, unperturbed) density, not
  the density that comes OUT of that step's diagonalization — so `E_DFT+U = α·n(α=0)`, trivially linear in α
  with a FIXED n, giving χ₀ = 0 by construction.  It measures nothing; it echoes the old occupation back.
  (This also explains why the first, crashed attempt's number looked suspiciously close to the converged χ
  value — `α·n(α=0)` is numerically close to `α·n(α, relaxed)` because n≈5.5 dominates over the small
  response Δn, not because χ₀≈χ.)  A genuine χ₀ would need the OUTPUT occupation of step 1 specifically
  (e.g. `&PRINT &PLUS_U`'s per-iteration occupation table, if CP2K populates it before the step-1 summary —
  unverified) rather than the final energy block — a real investigation, not a quick fix.  Not chased
  further; open.  Decks: `/home/janr/Code/qchem6-runs/batch/lane2/queue/03_mno_ckalpha_chi0.sh` (kept for the
  record, not for its number).
  **Still informative without the ratio**: CP2K's `PLUS_U_METHOD LOWDIN` on `VALENCE-LOWQ-VA` spans ALL of
  Mn's d-type shells (the basis note: "Mn keeps 7s+8d"), i.e. the SAME broad, multi-shell manifold philosophy
  as our own `every` convention, not hp.x's single contracted orbital.  χ = −1.894 sits in the SAME ORDER OF
  MAGNITUDE as our own χ for that manifold breadth (`every` −3.191, `ortho` −3.604) — nothing like hp.x's
  two-orders-of-magnitude-smaller numbers (−0.043, −0.071).  A independent THIRD code, with its own basis and
  projector implementation, again lands in the "healthy response" regime for a broad manifold rather than
  hp.x's near-cancellation — consistent with, not a refutation of, the manifold-breadth story above.  Logs:
  `/home/janr/Code/qchem6-runs/batch/logs/01_mno_ckalpha.log`, decks regenerated by
  `/home/janr/Code/qchem6-runs/batch/queue/01_mno_ckalpha.sh`.

★ **Is cRPA the same idea as ACBN0's renormalized occupancies? (user question, 2026-09-30.)**  No — three
genuinely different objects, all called "U":
- **ACBN0**: a single-SCF, mean-field ANSATZ.  Computes the BARE on-site four-index Coulomb integrals (real
  ERIs over the localized projector basis — the same kind of object our own `BareCoulombSource`/`ERI4Block`
  builds), then rescales them by an occupation-matrix-derived factor \f$\bar N^2\f$ meant to stand in for
  screening.  No polarizability, no response calculation, no frequency dependence is ever computed — it is
  read off ONE converged ground state's own density matrix.  This is exactly why the plan already states
  "ACBN0's \f$\bar N^2\f$ renormalisation is not a screening model" (§4) and why it overshoots hp.x by 2–3×
  on every manifold that HAS a usable hp.x number.
- **cRPA**: a genuine ab initio linear-response calculation.  Computes the full RPA polarizability of the
  crystal, then CONSTRAINS it — removes polarization channels that are transitions WITHIN the target
  correlated subspace (e.g. Mn-3d→Mn-3d), since those are what the Hubbard model itself is meant to handle,
  not something DFT should already be screening away.  Everything else (O-2p, other bands, interstitial
  states) still screens the bare Coulomb matrix element into \f$W=v/(1-P_rv)\f$.  **This is why cRPA
  sidesteps MnO's cancellation problem**: it never asks "how much does the Mn-3d shell's OWN occupation move
  when you poke Mn-3d" (the question that goes ill-conditioned when that shell is nearly inert) — it asks
  "how much does everything ELSE in the crystal screen the bare Mn-3d matrix element," which stays
  well-posed regardless of how inert Mn-3d itself is.
- **hp.x's LRT/DFPT (our A7)** is a third distinct thing again: an effective-model U backed out of an actual
  observable, dn/dα, of the interacting KS system — not a rescaled bare integral, not a screened matrix
  element.
All three get called "U," but they answer different physical questions, so none is obligated to agree
numerically with the others — consistent with ACBN0's already-logged 2–3× overshoot, and the reason every
oracle row in this doc is labelled by method/PP-match rather than averaged.

1. ~~**SrVO₃ FIRST**~~ ✅ **DONE 2026-09-25**: U(V 3d) = **6.2502 eV**, matched-PP LRT, 300 Ry (GTH V-q5
   is hard, converged to 0.94 mRy at 300), k 4×4×4 / q 2×2×2, χ₀(V,V)=−1.7822 → χ(V,V)=−0.1436
   (χ/χ₀≈0.081 — far more strongly screened than MnO's 0.96 or NiO's 0.62, as expected: SrVO₃ is a real
   metal).  Geometry from ABINIT's own `tests/tutoparal/Input/tucalc_crpa_1.abi` (acell 7.2605 bohr cubic).
   No AFM/2-step recipe needed, as predicted — one plain smeared `scf` step sufficed.  Found (and routed
   around, not fixed) a real `gth2upf` gap: Sr q10 semicore silently integrates to 2 e⁻ not 10 (the
   pseudo-atom EC allows only one shell per `l`); used Sr q2 instead.  Full detail + literature comparison
   still open: `IntegrationTests/QE/README.md` §A6.
2. ✅ **KCuF₃ DONE 2026-09-28: U(Cu 3d) = 8.1629 eV** — needed Cu specifically (now fixed).  **Corrected 2026-09-26 (user)**:
   Carta et al. actually use the HIGH-SYMMETRY CUBIC perovskite for KCuF₃, not the Jahn-Teller-distorted
   one ("for the purpose of this work, we consider KCuF₃ in the high symmetry cubic perovskite structure...
   for ease of computation" — main text) — so this is the SAME Pm-3m template as SrVO₃, no Wyckoff-table
   sourcing problem at all.  The one real gap: the paper says the cell (including cell parameters) is
   FULLY RELAXED but never states the numerical lattice constant, in the main text or the SI. Email sent to
   the corresponding author (Dr. Carta) asking for it, plus their functional/pseudopotential/k-mesh — most
   of which turned out to already be in the SI (PBE, PseudoDojo norm-conserving, 84 Ry, cold smearing
   0.01 Ry, spin-unpolarized; only the cell constant is genuinely missing), so the email was trimmed to ask
   only for that.  ⚠ **Their Hubbard projector is a Wannier function (MLWF via Wannier90), not an atomic
   projector** — even with their exact cell, our `hp.x` ortho-atomic number won't be on the same convention
   as theirs (same caveat as MnO/NiO's "quote the condition with the value").  Elements: K, Cu, F (all
   checked).  Parked, not blocking A6 — see item 4 for what ran instead while this waits.
3. ✅ **Sr₂FeO₄ DONE 2026-09-29: U(Fe 3d) = 8.0116 eV** — the other Carta et al. benchmark, and their own "d-only" entangled-band
   cautionary case (Fe-3d U_cRPA=0.42 vs U_LRT=6.94–7.29 eV) — directly relevant to re-testing our own NiO
   finding on a second material.  Tetragonal K₂NiF₄-type layered perovskite (I4/mmm), which — unlike cubic
   KCuF₃'s single lattice constant — has a free INTERNAL atomic coordinate (apical-anion Wyckoff 4e
   z-parameter) that symmetry alone doesn't fix.  Same email as KCuF₃ now also asks for *a*, *c*, **and**
   that z.  Elements: Sr, Fe, O (all checked).  Parked pending reply, same as item 2.
4. ✅ **The `LiMO₂` series (M=V–Ni)** from the cRPA-comparison paper (§7) — same rocksalt-derived
   framework across the row, so once one is running the rest are cheap geometry swaps.  Elements: Li, V,
   Cr, Mn, Fe, Co, Ni, O (Cr untested — likely needs the SAME Cu-style anneal fix, d5-4s1 near-degenerate
   with d4-4s2, same mechanism, not yet confirmed).  **LiCoO₂ DONE 2026-09-26: U(Co 3d) = 7.3070 eV**
   (`IntegrationTests/QE/README.md` §A6) — R-3m structure from Pinsard-Gaudart et al. 2011 (real
   single-crystal XRD, Co–O/Li–O bond lengths verified to <0.001 Å against literature after the
   hexagonal→rhombohedral-primitive conversion).  Needed the 2-step (smeared→fixed-occupation) recipe
   despite being nonmagnetic — LiCoO₂'s low-spin Co³⁺ d⁶ is a real band gap (LDA 1.56 eV), and the 2-step
   recipe turns out to be about the GAP, not about magnetism specifically.  Also surfaced a real `hp.x`
   rule: the Hubbard atom must be listed FIRST in `ATOMIC_POSITIONS`.  **LiVO₂ DONE 2026-09-26:
   U(V 3d) = 5.9526 eV** — R-3m structure from a real citation (Mat. Res. Bull. 27, 555 (1992),
   a=2.8388 Å, c=14.828 Å, z(O)=0.25749), bond-length-verified (V–O 1.988 Å, Li–O 2.121 Å).  A
   self-consistently `vc-relax`ed stand-in geometry (our own GTH-LDA via QE's own relaxer — `qchem` itself
   has no relax capability) was built and then discarded the moment the citation came in — never keep a
   defensible guess once a real number exists.  Comes out METALLIC under nonmagnetic LDA in this idealized
   (untrimerized) cell — no 2-step recipe needed, unlike Co.  **LiCrO₂ DONE 2026-09-26: U(Cr 3d) = 5.8111 eV**
   — the flagged aufbau risk did NOT materialize (`gth2upf` converged Cr q6 to 3d⁵4s¹ cleanly on
   the first try, no anneal fallback needed).  Structure from Garg et al., *Crystals* 9(1), 2 (2019),
   bond-length verified (Cr–O 2.002 Å, Li–O 2.114 Å vs the paper's 2.003/2.113 Å).  Also METALLIC
   under nonmagnetic LDA (Cr³⁺ half-filled t₂g³, same reasoning as V), no 2-step needed.  **LiFeO₂ DONE 2026-09-27: U(Fe 3d) = 7.5915 eV** — real LiFeO₂ is NOT R-3m (user correctly recalled this from memory before sourcing began: high-spin Fe³⁺ d⁵ favours cubic-disordered or orthorhombic phases), idealized into the row's SAME untrimerized R-3m template anyway (user decision: valid for validating our own U methodology, and matches the cRPA-comparison paper's own apparent convention — its k-mesh table groups LiFeO₂ with LiCrO₂/LiCoO₂ under the same isotropic 13×13×13 mesh, unlike LiMnO₂'s distorted-cell-shaped 6×13×6).  Geometry from Materials Project mp-19419 (GGA+U-relaxed, ICSD-backed), already primitive-rhombohedral, bond-length-verified.  Also METALLIC under nonmagnetic LDA — NOT a repeat of MnO's d⁵ cancellation trap, which was specific to MnO's spin-polarized AFM response, absent here in a spin-restricted (nspin=1) treatment.  **LiNiO₂ DONE 2026-09-27: U(Ni 3d) = 9.1730 eV — ROW COMPLETE.**  Real LiNiO₂ is also JT-active (Materials Project's own DFT+U-relaxed entry comes back genuinely monoclinic, confirming it), same idealize-anyway decision; geometry from Seo et al., *J. Electrochem. Soc.* 165 (2018) A2554, bond-length verified; metallic, largest U of the row (9.17 eV).  First use of the new checkpoint-archiving practice (`doc/OpenWork.md`): U_in≈0 `.save` kept at `IntegrationTests/QE/checkpoints/linio2_U0.save/` for a future self-consistent-U warm start.
5. **TiO₂, ZnO, FeS₂** — the remaining ACBN0-paper benchmark set (rutile TiO₂ and wurtzite ZnO already
   in ACBN0's own four-material study alongside MnO/NiO; FeS₂ is Macke's e_g-hybridisation warning case,
   §3 risk 2).  Elements: Ti, Zn, S (S untested — not yet checked in `gth_potentials.json`), Fe (checked).

### A7 scoping insights (session handoff, 2026-09-26)

★ **Not a plan for A7 — six things this session's `hp.x` runs (SrVO₃, LiCoO₂, plus MnO/NiO earlier)
learned about what building our own DFPT engine actually involves, for whoever scopes/designs it next.**
None of this is guesswork: it comes from reading `hp.x`'s own behaviour and output, and Carta et al.'s
supplemental derivation (`~/Code/supplementary.pdf` §I, and §III "Computational details").

1. **The architecture decomposes cleanly, and `hp.x` shows the seam.**  Every run this session printed
   the SAME four-stage structure: an outer loop over (Hubbard atom, q-point-in-star); an inner
   Sternheimer/CG solve per perturbation (the `chi: iter# ... residue` lines); a response-occupation-matrix
   accumulation (`hp_dnsq`); and a final χ₀/χ assembly + inversion (`U = χ₀⁻¹ − χ⁻¹`, supplementary Eq. 36).
   That is four separably-testable objects, not one monolith: a perturbation driver, a Sternheimer/CG
   numerical kernel, a response-density accumulator, and a χ-matrix assembler.  Design each as its own seam
   before writing the outer loop.
2. **The metal/insulator numerics are NOT a QE quirk to route around — they are the real content of the
   response weight, and this project already has the seam that handles it.**  Three materials across this
   plan (LiCoO₂ this session, MnO/NiO earlier) needed a smeared-then-fixed-occupation 2-step before `hp.x`
   would run at all, failing outright otherwise ("DOS at Fermi level too small... should NOT be treated as
   a metal"); SrVO₃, an actual metal, needed no such step.  The reason is the `(f_n−f_m)/(ε_n−ε_m)`
   response-weight sum: near-degenerate `ε_n≈ε_m` terms are only well-behaved under a Fermi(kT) occupancy
   (which broadens the denominator) or under an exact integer gap — never under smearing applied to a
   system that doesn't actually have a smearable Fermi surface.  `src/ElectronConfigurations/OccupationPolicy.C`
   ALREADY has exactly this axis — `occupancy {Integer, Fermi(kT)}`, composed once at assembly, not
   re-asked per fill (R2.21).  A7's response-weight construction should be built as a NEW reader of that
   SAME existing `OccupationState`/`OccupationPolicy` pair, not a second, parallel metal-detection path —
   the ground-state SCF already decided Integer vs Fermi(kT) once; the DFPT response should just ask it,
   the same way the ground state's own energy sum does.
3. **The k/q-mesh commensurability constraint is real and will need a k+q→k′+G lookup.**  Carta et al.'s
   own SI says it plainly: "our implementation is limited to q point grids that are commensurate with the
   k point grid... the k+q point is mapped onto another k′ point within the original grid, modulo a
   reciprocal lattice vector."  Whatever holds our Bloch/k-mesh state needs that map as a first-class
   query, not an afterthought discovered mid-implementation — C2's shifted-MP fold defect (`doc/OpenWork.md`
   §4) is adjacent territory and worth reading before designing this.
4. **Projector consistency is non-negotiable, and this is the one lesson the whole A1–A6 arc has hammered
   on repeatedly** (matched-PP oracle work, `IntegrationTests/QE/README.md` throughout): the response
   density MUST be projected through the SAME projector-flavour object (Löwdin / atomic / ortho-atomic)
   already used for the ground-state occupation matrix — never a second, parallel projection implementation
   for the response side.  If A7 ever disagrees with A6's own `hp.x` numbers on a shared material, the
   FIRST suspect must be "did the response side use a different projector," not "is DFPT wrong."
5. **A real scoping-ambition decision, not an implementation detail — decide it up front.**  `hp.x`'s
   "coarse-grained" LRT (supplementary §I.E) only ever produces a SAME-SITE U per Hubbard atom (it
   averages away everything else, Eq. 35).  The GENERALIZED LRT it's coarse-grained FROM (supplementary
   §I.B, restricted to a target subspace via Eq. 18) is naturally inter-site and multi-orbital from the
   start — and our own on-site `BareCoulombSource`/`ERI4Block` machinery is already RICHER than `hp.x`'s
   (native 4-index integrals, not a Wannier-projected proxy standing in for them).  Pin 23's "the manifold
   is an input, never a derivation" already generalizes one level for site/shell/irrep; A7 is the natural
   place to ask whether it should generalize ONE MORE level, to an inter-site V from day one, rather than
   building "hp.x-equivalent" first and bolting V on later.
6. **Real cost data to scope against, not a guess.**  This session's `hp.x` wall times (serial, 300 Ry,
   4–5-atom primitive cells, 2×2×2 q-mesh): SrVO₃ (metal) 55 min; LiCoO₂ (insulator, 2-step) 1h48m; MnO/NiO
   (insulator, 2-step, from earlier sessions) similar.  A native DFPT engine's FIRST validation target
   should be the cheapest one with an already-trusted external number to check against — SrVO₃ or NiO, not
   a from-scratch material — so a wrong answer is a bug in the new code, not a confound from an unfamiliar
   material.
6. **MnO's own O-2p, redone** with `ortho-atomic` instead of the fragile, non-orthogonalised `atomic`
   projector already used once (`IntegrationTests/QE/README.md`'s `mno.hp.in`, U(O 2p)=26.56 eV, flagged
   fragile at the time) — cheapest item on this list, no new geometry or elements, just a deck edit
   mirroring what A2b already did for NiO.

Each material needs the FULL recipe A2b/NiO went through: `gth2upf` per element, a sourced geometry, a
`pw.x` ground-state recipe hunt (ecut/k-mesh/smearing — MnO and NiO each took real iteration to get
right, expect the same here), then `hp.x`.  **A digression to fix something the recipe hunt exposes in
our own code is explicitly in scope, not a distraction** (user, 2026-09-25) — the Cu fix above is the
first example.

User comment: I have two papers on general Sternhaimer DFPT: ~/Code/{DFPT1-gonze1989.pdf,DFPT2-baroni2001.pdf}. If you already know how to do this no need to read them.  There is also a modern paper on precicely our application: ~/Code/DFPT3-timrov2018.pdf which is probably worth reading.  In order to keep our code framework clean I am mostly concered about changes to abstract interfaces.  So we should look at that before coding.  If we are adding behaviour to the Hamiltonian interfaces they should be general perturbation theory (PT) interfaces, not DFT specific.  SO if possible the same PT interfaces would work for MP2 corrections of an HF Hamiltonian.

---

## 1. Background — where we stand (2026-09-22)

**Built and banked** (increments 1–3, `doc/OpenWork.md` §1 step 5): the `Hubbard_U` term (scalar-generic,
spin-native, real-TRIM corner); the manifold as an **input** list of (site, shell, irrep, U) with no code
path assuming the Hubbard atom is the transition metal (pin 23); `ManifoldSymmetry` slots from the site's
coordination point group; three projector flavours (Löwdin column, contracted **atomic**, **ortho-atomic**
with U=0 spectators) matching CP2K's and QE's conventions; `BasisSet::BareCoulombSource` + `ERI4Block` for
the on-site bare integrals; the `ACBN0` estimator; and `SolidCalculation::ConvergeHubbardU`, the outer loop.

**Measured.**  On MnO and NiO, with the same UPF and projector as the oracle:

| | our ACBN0 | seed | ACBN0 paper | hp.x (linear response) | our occupations vs hp.x |
|---|---|---|---|---|---|
| MnO (d⁵) U(Mn 3d) | 10.77 eV | IonicSAD | 4.67 | 0.96 ⚠ | 4.93/0.45 vs 4.988/0.481 |
| NiO (d⁸) U(Ni 3d), ortho | **14.51–14.76** | IonicSAD | 7.63 | **5.27** | 4.920/3.445 vs 4.980/3.368 |
| NiO (d⁸) U(Ni 3d), full ortho set | **13.80–13.99** | IonicSAD | 7.63 | **5.27** | 4.894/3.204 |
| NiO U(O 2p) | 7.27 | IonicSAD | 3.0 | — | |

★ **The seed column exists because the seed CHANGED under these numbers** (`NiOSpec` went SAD → IonicSAD on
2026-09-23 when the "Ni²⁺ cannot be seeded" claim was retracted), and a table that cannot say which
configuration produced it is the failure mode the VA/VB span files exist to prevent.  **Re-measured on both
seeds, the answer does not move**: ortho 14.60/14.67 (SAD, 63 iterations, E = −109.2322219) vs 14.51/14.76
(IonicSAD, 114 iterations, E = −109.2322528) — 31 µHa apart, U_eff within 0.7 %, occupations to four
decimals; the full ortho set agrees to 0.7 % on Ni 3d and to **0.02 %** on O 2p.  ⇒ the banked numbers are a
property of the converged state, not of the guess that found it.  ⚠ The IonicSAD ortho run needed **114**
iterations where SAD needed 63, and reported "NOT converged" at the old 80-iteration default — **the third
time an iteration cap produced a wrong conclusion this week**; the probe default is now 200.

★★ **The finding that matters: the manifold is right and the functional is wrong.**  Our projected
occupations agree with hp.x to 1–2 % on the same cell, same pseudopotential, same projector — so the
factor 2.6 in the *value* is the U functional alone, not the projector, manifold or PP.

⛔ **TWO DIFFERENT COMPARISONS LIVE IN THIS TABLE AND THEY ANSWER DIFFERENT QUESTIONS.  Never merge them.**
- **ours ÷ published ACBN0** (1.8–2.8 across MnO-d, MnO-O2p, NiO-d, NiO-O2p) measures **PROJECTOR
  COMPLETENESS** — the paper's PAO-3G-projected plane-wave states carry ~60 % of their norm and ours carry
  ~100 %.  It is *not* evidence about screening, and a kernel fitted to close it would be fitting out a
  basis artefact.
- **ours ÷ an INDEPENDENT oracle** (hp.x, cRPA — never ACBN0) is the only thing that can test the
  **SCREENING** hypothesis.  Gate 3 holds that comparison, and an earlier version of gate 3 mixed the two,
  which is the circularity that had to be retracted on 2026-09-23.
**ACBN0's \f$\bar N^2\f$ renormalisation is not a screening model** — it vanishes as the basis approaches
completeness, which is exactly the regime a Gaussian code lives in.

**⚠ Three traps already paid for** (do not re-pay them):
1. **hp.x's MnO number is not usable** — d⁵ high-spin is filled-majority/empty-minority, \f$\chi_0=-0.045\to\chi=-0.043\f$
   (96 % of the bare response survives), so \f$\chi_0^{-1}-\chi^{-1}\f$ nearly cancels.  NiO (d⁸) is the healthy
   oracle: \f$\chi_0=-0.113\to\chi=-0.071\f$.
2. **An hp.x number is conditioned on its starting U.**  The NiO 5.267 eV is \f$U_{\rm LR}(U_{\rm in}=3\,{\rm eV})\f$
   (the decks inherit `U Ni-3d 3.0` from QE's benchmark); the MnO decks carry none.  Read the `HUBBARD` block
   before comparing.  A response is a function of the state it linearises about.
3. **LDA NiO at U=0 loses AFM-II — and it is PHYSICS, not seeding** (strengthened 2026-09-23).  Three runs
   now reach the same non-magnetic fixed point by different paths: SAD seed 79 iterations to
   E = −109.2693031, SAD+MOM *bit-identical* (so MOM holds the collapse rather than preventing it), and
   **IonicSAD 48 iterations to −109.2693702** — 67 µHa apart, same state, collapse signature
   \f$N\uparrow/N\downarrow\f$ = 4.2889/4.2889 in both.  ⇒ **two different seeds, different trajectories, one
   fixed point**: LDA NiO at Γ on this cell genuinely has the non-magnetic solution as its SCF fixed point,
   so gate 1 is asking a physics question and no amount of seeding will answer it.  The
   **ACBN0 outer loop is a non-result on NiO** — sites
   decouple to \f$U_{\rm eff}=-1.05\f$ eV, Hartree runs to 2.04× its floor.  MnO's clean monotone 8-step loop
   was a property of d⁵, not of the loop.  ⇒ **magnetic robustness under a changing U is a gate, not a detail.**

---

## 2. The decision, and why the APPLICATION decides it

Four routes to a U, not three — the fourth is `doc/OpenWork.md` step 5 item 5's parked fallback:

| | route | needs a supercell? | status |
|---|---|---|---|
| (a) | ACBN0 as is | **no** | ⛔ REFUTED — overshoots on both materials and both manifolds; its loop drives U *away* from the literature (MnO 10.8 → 16.8) |
| (b) | **ACBN0 with a screened interaction** | **no** | ★ THE ROUTE |
| (c) | hp.x values as input | **yes** (q-mesh ≡ supercell) | a bridge: legitimate under pin 12 (computed input ≠ hand-set knob), but one QE run per material *and per composition* |
| (d) | our own finite-difference linear response | **yes** | parked on cost; scriptable over the projector we already own |

★★ **The supercell column is the whole argument, and it is a TRANSLATIONAL-SYMMETRY statement** (user
asked, 2026-09-23).  (c) and (d) measure a *response* to \f$\alpha\hat P_I\f$ — a potential shift on ONE
site's manifold — and want the matrix \f$\chi_{IJ}=dn_I/d\alpha_J\f$ including \f$J\neq I\f$, because U
comes from a difference of INVERSES.  Under PBC the perturbation replicates coherently onto every equivalent
site, so a primitive cell returns \f$\sum_J\chi_{IJ}\f$ (the q=0 component) — one number where a matrix was
needed.  An isolated perturbation is genuinely not translationally symmetric, and there are exactly two
cures.  **Supercell**: break it explicitly, to a period long enough that the response to the nearest image
has decayed (the classic 2×2×2).  **DFPT**: never break it — perturb monochromatically with
\f$\alpha\hat Pe^{i\mathbf q\cdot\mathbf R}\f$, which IS Bloch-periodic and lives in the primitive cell,
then Fourier-sum \f$\chi_{IJ}(\mathbf R)=N_q^{-1}\sum_{\mathbf q}\chi_{IJ}(\mathbf q)e^{-i\mathbf q\cdot\mathbf R}\f$.
An \f$N_q\f$ mesh is MATHEMATICALLY EQUIVALENT to an \f$N_q\f$ supercell — the same dichotomy as frozen-phonon
vs DFPT phonons, which is what DFPT was invented for.  ⇒ our own `nq=1` non-result was q=Γ only, i.e. all
sites perturbed IN PHASE (`IntegrationTests/QE/README.md`: *"the q-mesh is the supercell size of the
linear-response method"*).  **Either cure costs a supercell-equivalent per composition.**

★ **(a) and (b) PERTURB NOTHING.**  They evaluate on-site integrals on the ground-state density in whatever
cell you already ran — no perturbation, no images, nothing to break, so the question never arises.  NiO's U
came out of a 4-atom cell in ~10 minutes.  A voltage curve needs U at many compositions; only the ACBN0
family delivers that without a supercell each time.  (a) is refuted ⇒ **the application selects (b)**, and
it selects it independently of the physics argument.

**What (b) actually requires, honestly.**  The bare \f$F^0\approx27\f$ eV is right (it is a bare integral on
the atomic Slater scale) and \f$\bar N^2\approx0.64\f$ measures basis completeness, so the missing physics is
the dielectric screening \f$\varepsilon^{-1}\f$.  ⛔ **\f$\varepsilon\f$ must be COMPUTED, never entered** — a
hand-set \f$\lambda\f$ or \f$\varepsilon\f$ is pin 12 one layer down, and the "~1/4 for a TMO" that motivated
this route was a *sizing* argument (does the gap have the right magnitude to be screening?), not a proposed
value.  Three ways to get it, in increasing order of cost and correctness:
- **Thomas–Fermi** from the valence density: no empirical input, but it is a *metallic* screening model and
  these are insulators — it will over-screen.  Cheap; expect it to be wrong, and measure *how* wrong.
- **\f$\varepsilon_\infty\f$ from our own response** to a field.
- **cRPA**: the principled route, the one the literature uses for screened U, and a solver increment.

★ **CP2K is the reference implementation for the machinery, not for the value.**  `dft_plus_u.F` only
*applies* U (`u_ramping` ramps toward a hand-set `u_minus_j_target` as an SCF aid — not a determination),
so there is nothing to learn there about values.  But `qs_linres_*.F`, `response_solver.F` and especially
`qs_linres_polar_utils.F` ("Polarizability calculation by dfpt … Berry phase operator … periodic Raman")
are **periodic DFPT response in a Gaussian basis**.  That is exactly the "solver increment we have no seam
for", implemented in our own basis type rather than QE's plane waves — worth reading before we conclude
anything about what a computed \f$\varepsilon\f$ costs us.

⚠ **The honest size of (b).**  ACBN0's renormalisation *is* the paper's screening model; replacing it with
a real screened kernel is not tuning ACBN0, it is substituting a different mechanism and keeping the
bookkeeping.  If the screened kernel ends up needing cRPA, ask explicitly whether we are building cRPA with
extra steps — and if so, whether (d) is cheaper after all.  **Step 3 below is designed to answer that
before we commit.**

---

## 3. The target: a Li\f$_x\f$Mn\f$_2\f$O\f$_4\f$ voltage curve without a U per configuration

The user's scheme: compute (U, J) at the three ordered compositions — **λ-MnO₂** (Mn⁴⁺, d³),
**LiMn₂O₄** (nominally Mn³·⁵⁺), **Li₂Mn₂O₄** (Mn³⁺, d⁴) — with no supercells, and use those for every
Li\f$_x\f$Mn\f$_2\f$O\f$_4\f$ configuration the cluster expansion needs.  Feasible exactly because (b) is
supercell-free (§2).

★ **One refinement, and it falls out of what is already built.**  Do not interpolate U on the global
composition \f$x\f$ — **assign U per Mn SITE by its local oxidation state**.  U is a property of the local
electronic configuration, not of a composition average, and pin 23 already makes the term take a per-site
list, so this needs no new capability:
- λ-MnO₂ → \f$U({\rm Mn}^{4+}), J({\rm Mn}^{4+})\f$
- Li₂Mn₂O₄ → \f$U({\rm Mn}^{3+}), J({\rm Mn}^{3+})\f$
- LiMn₂O₄ → **a transferability CHECK, not a third interpolation node**: if it charge-orders into Mn³⁺ +
  Mn⁴⁺, the per-site U's it produces must reproduce the two end members' values.  That is a falsifiable
  prediction, which an interpolation on \f$x\f$ is not.

★★ **AND U(O 2p) AT ALL THREE COMPOSITIONS, AS A FIRST-CLASS TARGET — not a spectator** (user, 2026-09-23).
Three reasons it is on the critical path for the VOLTAGE specifically, not just for the band structure:
(1) pin 23 exists because β-MnO₂'s decisive correction was on **O-p_z, not Mn-d**; (2) we measured
U(O 2p) = **7.27 eV** on NiO — against the paper's 3.0, and nowhere near negligible; (3) delithiation is
formally Mn³⁺→Mn⁴⁺ but a real fraction of the hole lands on **O 2p (ligand hole / oxygen redox)**, so
U(O 2p) moves the computed voltage directly.  The O sublattice also stops being equivalent once Li is
partially removed.  ★ **This costs nothing extra**: the `orthofull` arm already gives every listed manifold
its own estimate from ONE SCF — the NiO run returned Ni 3d, Ni 4s, O 2s and O 2p together.  List O 2p (and
O 2s) at U=0 and read all of them off each of the three runs.

Then each Mn in a CE training supercell takes its U from its own local environment, and **the CE training
runs need no new U calculations at all**.

**Risks, in the order they are likely to bite:**
1. ⛔ **Magnetic robustness — and it is NOT an argument for non-collinear.**  Ruled 2026-09-23 (user):
   **collinear is enough for a room-temperature voltage curve**; non-collinear order is not on this path.
   The risk that remains is narrower and real: the spinel Mn sublattice is the **pyrochlore lattice —
   geometrically frustrated**, so there are many near-degenerate COLLINEAR states, and an SCF can slide
   between them.  Two ways that bites: (a) the imposed order collapses under a U change, as it did on NiO
   (trap 3) — and NiO was an *unfrustrated* rocksalt AFM, so this is strictly harder; (b) worse for our
   purpose, different Li configurations land in DIFFERENT magnetic states, which pollutes exactly the
   energy DIFFERENCES the cluster expansion is fitted to.  ⇒ what matters is that the magnetic state be
   **consistent and reproducible across runs**, not that it be the true ground state.  **Gate 1 below
   exists for this and nothing else.**
2. **Mn³⁺ d⁴ high-spin is Jahn–Teller active** (e_g¹) — it is the whole story of LiMn₂O₄'s structural
   transition, and Li₂Mn₂O₄ is tetragonally distorted because of it.  Orbital-resolved U on e_g is both the
   best test of orbital resolution and the most dangerous place to apply it: Macke's FeS₂ warning is that
   correcting a hybridised e_g wrecked the lattice parameter.  Do U at **fixed geometry** first; U ↔ JT
   distortion is a coupled problem and must not be entered accidentally.
3. **LiMn₂O₄ charge ordering is delicate at the DFT level.**  If it does not charge-order, the
   transferability check is weaker — say so rather than reading the average as a third node.
4. **O 2p is not a spectator by assumption.**  We measured 7.27 eV on NiO (paper 3.0), and pin 23 exists
   because β-MnO₂'s decisive correction was on O-p_z.  Carry O 2p in the manifold list at U=0 and *measure* it.

---

## 4. The work, in order

Gates 1–3 are cheap enough for now-to-Oct-5; the long unattended runs are sized for the **Oct 6–20** window
(user away).  ⛔ Unattended runs go through **`scripts/memsafe -p`** (cgroup + OOM shield), never bare.

**Prerequisites (no physics).**  ⚠ **Read the STATUS, not the prose** — the Li entry below documents a
test that is DESIGNED but NOT RUN, and it was mistaken for finished work once (2026-09-23).  The queue in
§0 is authoritative: Li q1 is **done and committed**, the q1-vs-q3 discriminator is **B4 and not blocking**,
the spinel structures are **B1 and are the real prerequisite**.
- **A Li valence basis** — ✅ **q1 DONE AND COMMITTED**; ⛔ **the q1-vs-q3 DISCRIMINATOR IS NOT RUN**
  (queue item B4, and it blocks nothing else).  q1 vs q3 is a REAL TEST, not a formality (user, 2026-09-23).
  There is no
  `LI` block in any `valence_lowq_*.bsd`.  GTH LDA offers Li **q1** (2s¹ only; 1s frozen into the core) and
  **q3** (1s²2s¹ explicit).  The tension is specific to a cathode: Li is nearly fully ionised, so q1's frozen
  core is being asked to describe an ion whose valence electron has LEFT — exactly the regime where a frozen
  core is least justified, because the 1s sees a different potential once 2s is gone.  q3 cannot have that
  problem but costs functions on every Li site in a 14-atom cell.  ⇒ mint BOTH with
  `valgen --nmax 60 --floor`, and settle it on a number the voltage cares about: **the Li intercalation
  energy** (E[LiMn₂O₄] − E[λ-MnO₂] − E[Li]) computed both ways.
  ✅ **BOTH MINTED AND VALIDATED 2026-09-23.**  q1: `--q 1 --shell 0:5:0.03:2`, converged, E = −0.189367 Ha,
  0.11 mHa below the pool floor — committed to `valence_lowq_{va,sph}.bsd` as the VALENCE variant, matching
  those files' own convention (Mn is q7 not q15, Na is q1 not q9).  q3: `--q 3 --shell 0:8:0.05:60`,
  converged, E = −4.23584 Ha, gap 0.107 mHa — validated but NOT committed, because a `.bsd` block is keyed
  by ELEMENT and the two cannot coexist in one file.  ⛔ **Two pieces of plumbing the discriminator needs
  before it can run**, both found while minting: (i) a basis-file variant that can carry q3's `LI` block
  (a `BasisSetData` enum value + its two map entries — cheap, but wire it with a real consumer, not
  speculatively); (ii) `CleanupCandidates.md` row **D-SEED1** — the atomic seed library is keyed by
  (Z, functional) and returns the FIRST match, so "neutral Li" is ambiguous between q1's 1 electron and
  q3's 3, decided silently by file order.  The A/B is not trustworthy until that throws or is keyed on q.  If q1 and q3 agree there, q1 is free
  throughput for the CE training set; if they do not, q3 is mandatory and we have learned why.  Seed density
  too (Li⁺ is a stripped cation, so `HasAtomicSpinPair` correctly calls it non-magnetic).
- **The Mn seed ions are already checked and they are fine** (2026-09-23): Mn³⁺ (d⁴) and Mn⁴⁺ (d³) both
  generate cleanly at `--nmax 60` — converged, moments 4.000 / 3.000, ⟨r⟩ 1.078 / 1.021 bohr, shrinking with
  charge as a cation should.  ⇒ **IonicSAD is available for all three spinel compositions**, which matters
  for gate 1: the ionic seed is the basin chooser, and a neutral-superposition seed would be a much worse
  guess for Mn⁴⁺ (four electrons away in a 7-electron PP) than it was for Ni²⁺ (two).
  ⚠ **A retraction that belongs here** (it produced a work item that no longer exists): "Ni²⁺ cannot be
  seeded, high-spin d⁸ is minority-d³ in a five-fold shell and the atom occupies whole irreps" was WRONG.
  The evidence — a non-aufbau run with ⟨r⟩ = 2.91 bohr for a cation — predated the iteration-cap fix made in
  the same session; at `--nmax 120` Ni²⁺ converges cleanly (charge 8.001, moment 2.000, ⟨r⟩ 1.008, between
  Mn³⁺'s 1.078 and Ni³⁺'s 0.938).  NiO now uses IonicSAD like MnO.  What is TRUE is milder: a partially
  filled MINORITY shell under a filled majority is a long descent (two orbitals 2e-6 Ha apart straddle the
  boundary), while a partially filled MAJORITY over an empty minority is not — Mn³⁺, Mn⁴⁺, Ni³⁺ (d⁷) and
  Co³⁺ (d⁶) all converge by 60.  **Both halves of that were the cap, twice.**
- ⛔ **The three spinel structures in `materials.json` — NOT STARTED, and this is queue item B1: the real
  blocker for every run in track B.**  Primitive cells: λ-MnO₂ 12 atoms (4 Mn, 8 O),
  LiMn₂O₄ 14 (2 Li, 4 Mn, 8 O), Li₂Mn₂O₄ 16 — plus whatever magnetic decoration gate 1 settles on.
  ⚠ Lattice constants are anchors: take them from a named source and say which.

**Gate 1 — does the magnetic state survive a U change?**  (the NiO lesson; blocks everything downstream.
⚠ **Seeding is NOT the lever** — trap 3 now has three runs and two different seeds converging to the same
non-magnetic fixed point, so a better seed will not buy the order.  The levers are the mixer's
magnetisation channel (§4 row N3), kT, and U itself.)
Run λ-MnO₂ and LiMn₂O₄ at fixed U = 0, 2, 4 eV and watch the integrated site moment.  If the order dies as
it did on NiO, no loop on this material means anything and the fix (mixer preconditioning in the
magnetisation channel, §4 row N3; or a different ordering) comes first.  **Cheap: three short SCFs each.**

★ **RUN 2026-09-25 — PRELIMINARY VERDICT: ORDER SURVIVED IN ALL SIX ARMS.**  `gpwprobe gate1 <material>
[U_eV]` (new command, generic over any `materials.json` entry: finds the Mn sites, puts the same U on
every one, reports `SolidCalculation`'s own `RunDiagnostics` — already material-agnostic, nothing new to
build).  Ferromagnetic seed (materials.json's decoration), `GPW_SPHERICAL=1` (the Cartesian-d Hubbard_U
guard fires otherwise, same as MnO/NiO), `NMAX=80`, Γ-only:

| material | U (eV) | seed→peak→final (e) | verdict | Eee end/floor |
|---|---|---|---|---|
| λ-MnO₂ | 0 | 4.9 → 5.005 → 4.31 | SURVIVED | 1.30 |
| λ-MnO₂ | 2 | 4.9 → 4.9 → 4.262 | SURVIVED | 1.64 |
| λ-MnO₂ | 4 | 4.9 → 4.9 → 4.459 | SURVIVED | 1.71 |
| LiMn₂O₄ | 0 | 4.832 → 5.159 → 3.374 | SURVIVED | 2.58 |
| LiMn₂O₄ | 2 | 4.832 → 5.126 → 3.646 | SURVIVED | 2.64 |
| LiMn₂O₄ | 4 | 4.832 → 5.426 → 4.739 | SURVIVED | 2.64 |

Every final value sits within ~30 % of its own seed/peak — nothing resembling NiO's collapse to <2 % of
peak.  ⇒ **no magnetic-order redesign is needed before gate 2/3.**  ⚠ **Caveat, and it is real: NONE of
the six converged at NMAX=80** — this is a trajectory reading (the instrument `SolidCalculation` runs
regardless of convergence), not a converged-energy verdict.  λ-MnO₂'s Hartree sloshing (1.3–1.7×) is
mild; LiMn₂O₄'s (2.58–2.64×) is not, consistent with §3 risk 1's own flag that the mixed-valence,
geometrically frustrated (pyrochlore) Mn sublattice would be the harder of the two to settle — LiMn₂O₄'s
run here also uses a single Mn3+-everywhere approximation for the formally 3.5+ average, not the true
charge-ordered state, which is a plausible contributor to the extra sloshing.  **Not yet production
numbers for gate 2/3** — needs either more iterations or a mixer/schedule tuned the way MnO/NiO's own
recipe was (doc/OpenWork.md N3 row), before an energy or a U from these cells is quoted.  Evidence:
`~/Code/qchem6-runs/gate1/*.log`.

**Gate 2 — size the runs.**  One converged SCF per composition at Γ, timed and RSS-logged.  Estimate to
beat: ~316 basis functions for LiMn₂O₄ against MnO's 118, so expect 3–8× MnO's ~6 min ⇒ 20–50 min per SCF,
⇒ an 8-step outer loop is an overnight run per composition.  **Measure it; do not plan on the estimate.**

**Gate 3 — BUILD THE ORACLE SET, then test the screening hypothesis.**  ⚠ **This gate replaces an earlier
version that was CIRCULAR** (caught 2026-09-23): it proposed fitting one \f$\varepsilon^{-1}\f$ to move all
four measured manifolds onto "their targets", but two of those targets were the ACBN0 paper's OWN values
(MnO O-2p 2.68, NiO O-2p 3.0).  Fitting to those fits out PAO-3G incompleteness — precisely the thing we
established is NOT screening.  **A target must be independent of ACBN0.**  Against independent oracles:

| manifold | ours | independent oracle | ratio |
|---|---|---|---|
| NiO Ni-3d | 13.89 | `hp.x` **5.27** (linear response) | **2.64** ← the only SOLID point |
| MnO Mn-3d | 10.77 | literature 4–7 eV (cRPA / LR-cDFT) | 1.5–2.7 (a range) |
| NiO O-2p | 7.27 | cRPA \f$\bar U_{2p}\gtrsim4\f$ eV, "independent of the TM" (ACBN0 paper's own citation [97]) | ≲1.8 (a bound) |
| MnO O-2p | 7.36 | same | ≲1.8 (a bound) |

⛔ **Two things fall out, both unwelcome.**  (1) We have **ONE solid oracle point**, not four — a range and
two bounds otherwise — so "does a single factor fit?" is UNDERDETERMINED and cannot refute anything yet.
(2) What structure there is looks like **TWO CLUSTERS, not one factor**: d ≈ 2.6, p ≈ 1.8.

⛔ **And a physics argument runs the WRONG WAY.**  A bulk \f$\varepsilon\f$ screens everything equally.  Make
it q-resolved — \f$U\sim\sum_q|\rho_\phi(q)|^2\varepsilon^{-1}(q)v(q)\f$ — and a MORE localized orbital has
broader \f$\rho_\phi(q)\f$, samples larger q, where \f$\varepsilon^{-1}(q)\to1\f$: so TM 3d should need
**less** screening correction than O 2p.  We observe it needing **more** (2.6 vs 1.8).  That is backwards for
a dielectric picture, and it is not explained by orbital differences — our occupations match the oracle's to
1–2 %.  ⇒ **this is the specific observation that could refute (b)**, and it is available for the price of
assembling the oracle set.

★★ **AND THE POINTS ARE ALREADY BUILT: ABINIT DOES cRPA *AND* LINEAR-RESPONSE J** (user asked 2026-09-23;
verified in `~/Code/abinit/src`).  Two capabilities QE does not give us, each hitting one of gate 3's two
weaknesses:
- **`ucrpa`** (+ `ucrpa_bands`, `ucrpa_window`) — **constrained RPA**, which is *methodologically
  independent of linear response*.  This is the one that matters: the O-2p "oracle" above is a literature
  **bound** (cRPA \f$\gtrsim4\f$ eV) and cRPA is exactly the method that produced it, so `ucrpa` turns a bound
  into a VALUE on our own materials — and gives an independent Mn-3d number on MnO, where hp.x is broken by
  the d⁵ shell (trap 1).  Two of the four soft rows become solid.
- **`lruj`** (`src/98_main/lruj.F90`, "Linear Response U **and J**", citing Cococcioni & de Gironcoli PRB 71
  035105) driven by `macro_uj` / `pawuj_det` ("Determine U (or J) parameter", with/without a compensating
  charge bath).  **hp.x gives U only**, which is precisely why §5's first open question — *does J transfer
  the way U does?* — has had no oracle.  ABINIT answers it.
- ⚠ **BUT IT IS A DIFFERENT-PP ORACLE, AND THE TABLE MUST SAY SO.**  ABINIT's Hubbard is **PAW-only** (every
  Hubbard file lives under `src/65_paw`; the keyword is literally `usepawu`), and our GTH is
  norm-conserving.  The discipline that made the QE comparison worth anything was `gth2upf` — the oracle ran
  OUR pseudopotential and OUR projector, which is why the occupations matched to 1–2 % and the disagreement
  could be pinned on the functional.  We cannot do that here without a PAW dataset.  ⇒ ABINIT numbers are a
  SECOND OPINION, not a matched comparison; cRPA's methodological independence is worth more here than PP
  matching, but **label every oracle row `matched` or `different-PP`** or the next reader will average them.
- Cost note: ABINIT is an MPI build like the rest — `mpirun` always, and it additionally needs
  `--force-mpirun` (`CLAUDE.md`).  cRPA needs bands/windows chosen, which is a real input-convergence
  question of its own; treat the first run as a recipe hunt, not a number.

★ **PySCF is the OTHER install, and it is for the KERNEL, not for the values** (so: second, not first).
Free, `pyscf.pbc` gives periodic Gaussians with k-points and density fitting, and — the part we would
actually use — a four-index ERI engine with **range-separated (`omega`) integrals**, i.e. the exact object
route (b)'s screened kernel needs checking against.  It buys nothing for gate 3.
⚠ **Not `pip install` on this box as it stands** (checked 2026-09-23): Python is **3.14.4** and
`/usr/lib/python3.14/EXTERNALLY-MANAGED` is present, so PEP 668 makes apt-managed Python refuse installs
into system site-packages.  `python3-pip` provides BOTH `pip3` and `pip` (they are the same thing here —
there is no python2), but the name is not the issue; the venv is:
```
sudo apt install python3-venv python3-pip
python3 -m venv ~/Code/pyscf-env && ~/Code/pyscf-env/bin/pip install pyscf
```
Inside a venv `pip` and `pip3` are identical, so the question stops mattering.  ⚠ **Python 3.14 is new
enough that a PySCF manylinux wheel may not exist yet** — if pip falls back to building from source it will
want cmake + libcint, which is where the project's standing "source builds are the DEFAULT" policy takes
over anyway (and is arguably what we want, since the point of PySCF here is to read its integral engine).
`python3-numpy` / `python3-scipy` are in apt if a source build needs them outside the venv.

★★ **FIRST RESULT, 2026-09-23 — THE ONE-FACTOR TEST PASSES ITS FIRST REAL TRIAL, AND IT CORRECTS THE
ARGUMENT ABOVE.**  Run on paper with the newly installed PySCF (2.14.0, `~/Code/pyscf-env`), using its
range-separated (`omega`) four-index integrals on OUR OWN contractions read off the `[+U radial]` banner —
`scripts/gate3_screening_test.py`:

| manifold | bare \f$F^0\f$ | independent target | needed ratio | needed \f$\omega\f$ |
|---|---|---|---|---|
| Ni 3d | 24.29 eV | 5.27 (hp.x, **matched PP**) | 0.217 | **1.039 a.u.** |
| O 2p | 20.97 eV | ≳4 (cRPA **bound**) | 0.191 | **0.975 a.u.** |

**The two manifolds want the same screening length to 6 %.**  ⇒ ⛔ **MY "TWO CLUSTERS, AND THE LOCALIZATION
ARGUMENT RUNS THE WRONG WAY" OBJECTION WAS APPLIED TO THE WRONG QUANTITY** and is withdrawn.  It compared
ratios of our **\f$U_{\rm eff}\f$**, which already carries ACBN0's \f$\bar N^2\f$ renormalisation — i.e. the
projector-completeness effect — so it was measuring completeness and screening mixed together.  Against the
**BARE \f$F^0\f$**, which is what a screened kernel actually modifies, the picture is one factor.  And the
localization physics comes out RIGHT, not backwards: the same \f$\omega\f$ screens the compact Ni 3d less
(ratio 0.217) than the diffuse O 2p (0.191), which is the required direction.
⚠ **Three things keep this from being a verdict.**  (1) \f$1/\omega\approx1\f$ bohr is a SHORT length —
inside the 3d orbital itself and well inside the 3.94 bohr Ni–O bond — so whatever this is, calling it
"screening by the medium" needs an argument; an `erfc` cut-off at 1 bohr reshapes the on-site
self-interaction rather than dressing it.  (2) TWO points, one of which is a BOUND, is not a fit.  (3) ✅ **RESOLVED 2026-09-23 (A1) — NOT a
cross-check failure.**  PySCF's plain \f$(2l+1)^{-2}\sum_{mm'}(mm|m'm')\f$ = 24.29 eV and our banner's
eq-10 \f$\bar U_{\rm bare}\f$ = 27.3 eV are two DIFFERENT quantities by construction — eq 10's numerator is
an unrestricted double sum over the density matrix while its pair-count denominator (eq 10c) excludes
same-orbital/same-spin self-pairs, so \f$\bar U_{\rm bare}\f$ is inflated above the plain average by a
factor set by the spin-resolved occupation, not by any integral disagreement.  Derived analytically
(diagonal \f$P\f$, uniform-within-spin \f$N_m\f$) and checked on a live `gpwprobe nio` run (`NIO_ACBN0=1
GPW_SPHERICAL=1 NIO_U_RADIAL=ortho`, U=0): predicted \f$\bar U_{\rm bare}=26.99\f$ eV from the run's own
bare \f$N_\uparrow/N_\downarrow\f$ (4.9669/4.9667 and 2.9597/2.8819 on the two sites) against the actual
banner 27.018 / 27.131 eV — **0.1–0.5 % apart**, the residual being the real (non-uniform,
crystal-field-split) per-orbital occupation the approximation ignores.  `scripts/a1_bare_f0_check.py`.
⇒ our ERI4Block integrals and PySCF's agree; gate 3's ω values (built on PySCF's plain F0,
the correct object for a kernel that screens the bare integral TENSOR) are unaffected.  ⚠ This same
`gpwprobe nio` run independently reproduced **trap 3** (U=0 LDA NiO's AFM-II order decaying to 1.7% of its
peak, "DENSITY-DEGENERATE... benign") — a diagnostic run only, no U quoted from it.
⛔ **AND THE 6 % IS FRAGILE — IT LIVES OR DIES ON THE O 2p TARGET, WHICH IS A LOWER BOUND** (measured
2026-09-23, same script).  The needed \f$\omega\f$ for O 2p against a range of assumed targets, compared with
Ni 3d's 1.039:

| O 2p target | needed \f$\omega\f$ | disagreement with Ni 3d |
|---|---|---|
| 3.0 eV | 1.185 | 14 % |
| **4.0 eV** (the cRPA bound) | **0.975** | **6 %** |
| 5.0 eV | 0.828 | 20 % |
| 6.0 eV | 0.715 | 31 % |
| 7.0 eV | 0.624 | 40 % |

The literature figure is \f$\gtrsim4\f$ eV — a **lower bound** — so the true value may well be 5 or 6, and the
one-factor agreement degrades fast above 4.  ⇒ **the ABINIT `ucrpa` run is DECISIVE, not decorative**: it is
the difference between this hypothesis standing and falling, and it is the single highest-value action in
the whole plan.  Do not build a kernel before it.

★★ **FIRST ABINIT `ucrpa` RESULT, 2026-09-24 (NiO O-2p, "dp-dp" model) — a real number, with a real bug
isolated along the way.**  Decks in `~/Code/abinit-runs/` (a scratch working directory outside the repo,
mirroring `~/Code/cp2k-runs/` — the `.abi` recipes themselves are worth banking, the multi-GB WFK/DEN
restart files are not).  Cell + PAW: ABINIT's own validated `tests/tutorial/Input/tlruj_1.abi` (AFM-II,
\f$a=7.88\f$ bohr, matching `NiOSpec`) with `Psdj_paw_pbe_std/{Ni,O}.xml` (**PBE, not LDA** — no LDA Ni
PAW is bundled locally, so this run is a different-PP **and** different-XC second opinion, labelled as
such).  Fatbands (`pawfatbnd`) on a first GS+NSCF pass showed NiO's Ni-3d/O-2p valence complex (bands
11–26) is **thoroughly hybridised** — unlike SrVO3's cleaner d/p separation, there is no clean "pure O-2p"
band range, so the "dp-dp" model (Amadon2014: Wannier + screening exclusion both spanning the whole
hybridised complex) is close to forced, not chosen.
- ⛔ **Three real bugs hit and fixed on the way, worth recording so a future session does not re-pay
  them.**  (1) `abinit`/`--dry-run` run bare instead of under `mpirun` hangs (CLAUDE.md's warning,
  re-confirmed the hard way).  (2) The `plowan_bandi/bandf/natom/iatom/nbl/lcalc/projcalc` variables are
  **COMMON across datasets, unsuffixed** — constructed in the Wannier dataset, read back in the screening
  and effective-interaction datasets; suffixing them `2` (dataset-2-only, the natural-looking choice)
  makes dataset 3 read zeros ("Lower and upper values of the selected bands 0 0") with no other complaint.
  Caught by diffing against ABINIT's own validated `tests/tutoparal/Input/tucalc_crpa_2.abi`, not by
  guessing. (3) `getwfk3 -2` (relative) is an off-by-one for a 4-dataset recipe — dataset 2 (the
  well-diagonalised NSCF + Wannier construction) is `-1` away from dataset 3, not `-2`.
- ⛔ **The optdriver=4 (effective-interaction) step CRASHES under MPI parallelism** (`-np 4`: segfault in
  `m_prep_calc_ucrpa.F90` after k-point 1/64, no diagnostic) but runs to completion under `-np 1`
  (serial) — matching, independently, `tucalc_crpa_2.abi`'s own `TEST_INFO` comment: *"results with 24
  procs are non-reproducible at present! ... There must be a bug."*  **Always run this dataset serial.**
  OpenMP threading (`OMP_NUM_THREADS>1`, `-np 1`) was tried as a speed lever and made things SLOWER for
  this cell size (iteration 9 in 13.5 min vs. plain serial's iteration 21 in 5.75 min) — overhead
  dominates at this problem size; plain serial is the fastest working recipe found. A serial NiO run
  (4 datasets) costs **~100 minutes wall**, dataset 4 alone taking the majority of it.
- ⛔⛔ **A SECOND, genuine ABINIT bug, isolated (not just worked around): the "Average U and J" summary
  table is WRONG under `nsppol=2`.**  Every per-block calculation ("Hubbard cRPA interaction for w=1,
  U=...") printed a sensible, ecuteps-responsive number (bare 4.1456 eV, cRPA 1.2359 eV at
  `ecuteps=4`; cRPA 1.2141 eV at `ecuteps=7` — a small, correctly-signed convergence trend, screening
  reducing from bare as physics demands).  But the "Average U and J as a function of frequency" table
  printed immediately after — which in ABINIT's own SrVO3 tutorial output is documented to just ECHO that
  same per-block number — instead printed a DIFFERENT, LARGER value (4.9437 eV at `ecuteps=4`, 4.8565 eV
  at `ecuteps=7`), **identical across every atom and every one of the four Up-Up/Up-Down/Down-Up/Down-Down
  spin combinations**, in the SAME run.  Removing the `usepawu`/`dmatpuopt` "for printing" block (copied
  from the tutorial) changed nothing — that hypothesis is REFUTED.  **Decisive isolation**: reran ABINIT's
  own `tucalc_crpa_2.abi` (SrVO3, non-magnetic, `nsppol=1`) UNMODIFIED and reproduced its published
  reference numbers **bit-for-bit** (bare 15.3789 eV, cRPA 2.7546 eV, J 0.5997 eV, "Average U and J" table
  MATCHING the per-block number exactly, as documented) — confirming our build and methodology are sound.
  SrVO3's tutorial never exercises the 4-way spin-combination average at all (`nsppol=1` prints only ONE
  "Up-Up" block, nothing to combine) — exactly the code path NiO's AFM-II cell forces open.  ⇒ **the
  `nsppol=2` cross-spin averaging in ABINIT's ucrpa "Average U and J" print is broken**; this is the same
  code region the tutorial's own test-suite comment already flags as unreliable, now shown broken in a new
  (magnetic) way.  **Trust the per-block "Hubbard cRPA interaction" number, never the "Average U and J"
  summary, on any `nsppol=2` ucrpa run.**
- ⛔ **A THIRD bug: Ni-3d (l=2) crashes with a Fortran integer overflow** at the same point (start of the
  optdriver=4 k-point loop) that the O-2p (l=1) run gets past cleanly — same deck otherwise (bands 11-26,
  same cell).  Not yet root-caused; flagged in §6.3, not on the critical path (Ni-3d already has hp.x's
  matched-PP oracle at 5.27 eV; ABINIT's Ni-3d row was always going to be a second opinion, not decisive).
- ★ **Best current number: \f$U_{\rm cRPA}({\rm O\ 2p},\,\omega=0,\,dp\text{-}dp\text{ model})\approx1.21\f$–\f$1.24\f$
  eV** (ecuteps 7→4 Ha; \f$J\approx0.33\f$ eV), from the per-block print, trusted per the isolation above.
  ⚠ **This is well BELOW the \f$\gtrsim4\f$ eV literature bound gate 3's 6 % agreement leans on** — but the
  SrVO3 tutorial's OWN convergence table (§5 there) shows the SAME material/orbital's U swinging from 1.6
  to 12.0 eV **purely from the choice of screening-exclusion model** (\f$t_{2g}\f$-\f$t_{2g}\f$ vs \f$dp\f$-\f$dp\f$ vs
  \f$d\f$-\f$dp\f$(a) vs \f$d\f$-\f$dp\f$(b)), so a 1.2 eV "dp-dp" number is not necessarily in tension with a
  literature bound quoted under a DIFFERENT, unstated model convention — **the two are not yet
  comparable**, and making them comparable (matching whatever model convention the literature's
  \f$\gtrsim4\f$ eV figure actually used) is the next real step, not a parameter sweep on this one.
  ⛔ **Not yet quotable against gate 3** until that model-matching is done — recorded here as progress, not
  as gate 3's answer.

★★ **A3 RUN 2026-09-25 — WITH A REAL O-2p ORACLE, THE ONE-FACTOR TEST FAILS, AND ROUTE (b) IS REFUTED.**
A2b (`IntegrationTests/QE/README.md`) gave a matched-PP hp.x value, U(O 2p) = 8.5139 eV, replacing the old
4 eV **bound** the 6 % agreement above was built on.  Re-ran `scripts/gate3_screening_test.py` +
`gate3_omega_sensitivity.py` against it:

| manifold | bare \f$F^0\f$ | independent target | needed ratio | needed \f$\omega\f$ |
|---|---|---|---|---|
| Ni 3d | 24.29 eV | 5.27 eV (hp.x, matched PP) | 0.217 | **1.039 a.u.** |
| O 2p | 20.97 eV | **8.51 eV** (hp.x, matched PP — no longer a bound) | 0.406 | **0.513 a.u.** |

**The two needed screening lengths are now 50.6 % apart** (the sensitivity scan's own worst case at the old
bound's upper end, 7.27 eV, already read 42 %; the real value is past even that).  ⛔ **AND THE
LOCALIZATION SIGN HAS FLIPPED BACK TO WRONG.**  With the real target, O 2p's needed ratio (0.406) is now
LARGER than Ni 3d's (0.217) — the diffuse orbital needs LESS reduction from bare than the compact one.
That is backwards for a dielectric picture (a compact orbital samples large \f$q\f$, where
\f$\varepsilon^{-1}(q)\to1\f$, so it should need LESS correction, not more) — the exact objection raised
2026-09-23 and WITHDRAWN that same day when the 4 eV bound made the two manifolds agree to 6 %.  **The
withdrawal was itself conditioned on a wrong number.**  Both halves of the plan's own stated refutation
criterion (§2: *"if d and p need systematically different factors AND the localization sign stays wrong,
the residual is not bulk screening and (b) is refuted"*) are now met, on real oracle data for BOTH
manifolds, not a bound on one of them.  **⇒ ROUTE (b) [ACBN0 with a screened bulk kernel] IS REFUTED.**
⚠ **What is NOT refuted**: the underlying instrument (bare on-site ERIs, the manifold/outer-loop machinery,
§7's pin 23 addendum on LRT/cRPA as pluggable strategies) — only the SPECIFIC hypothesis that a single
basis-independent bulk screening length, applied uniformly to the bare integral tensor, reproduces both
oracles.  Per §2's own routing table, refuting (a) [ACBN0 as-is] and now (b) leaves **(c)** [hp.x/QE values
as input, per-composition — legitimate under pin 12, but a supercell-equivalent per composition] and **(d)**
[our own finite-difference/DFPT linear response] as the surviving routes, with §2's own observation that
(b) and (d) share much of their cost if a computed \f$\varepsilon\f$ needs DFPT anyway — **this is a plan
shape question for the user, not a call to make alone.**

★ **A cross-check that is independent of the kernel FORM.**  A flat \f$1/\varepsilon\f$ (no length scale at
all) would need \f$\varepsilon = 4.61\f$ for Ni 3d, against NiO's experimental \f$\varepsilon_\infty\approx5.7\f$.
The needed value sitting slightly BELOW \f$\varepsilon_\infty\f$ is the physically right ordering — an on-site
U samples large q where \f$\varepsilon^{-1}(q)\to1\f$, so it must be screened LESS than the macroscopic limit.
That is mildly supportive of a dielectric picture *whatever* kernel shape turns out to be right.

⇒ this RAISES the value of the ABINIT runs below: turning the O 2p bound into a cRPA VALUE, and adding
Mn 3d / Mn O-2p rows, is what turns a suggestive two-point coincidence into a test.

**So the gate is: get more independent points BEFORE building a kernel.**  ★ The reframing that makes this
affordable: **(c)/(d) are unaffordable PER COMPOSITION but perfectly affordable ONCE.**  Use them for what
they are good at — a calibration set on a handful of materials where linear response is healthy (`hp.x`
run properly, i.e. with its \f$U_{\rm in}\f$ declared) — then test the screened kernel against it, then run
(b) in production across the composition sweep.  No amount of cleverness about (b) fixes a one-point oracle
set.  **Refutation criterion:** if d and p need systematically different factors AND the localization sign
stays wrong, the residual is not bulk screening and (b) is refuted; (d) becomes the route and this plan
changes shape.

**Gate 4 — SUPERCELL vs q-MESH IN OUR CODE, and the k/irrep parallel axis** (user, 2026-09-23:
*"they may not be equivalent on our code"*).  §2's equivalence is MATHEMATICAL.  Computationally the two
diverge, and on this box they diverge in the same direction for three independent reasons:
- **Asymptotics favour the q-mesh by \f$N_q^2\f$.**  A supercell is one \f$N_qN\f$ problem — \f$O((N_qN)^3)\f$ in
  the dense-algebra part.  A q-mesh is \f$N_q\f$ problems of size \f$N\f$ — \f$O(N_qN^3)\f$.  At \f$N_q=8\f$ that
  is 64× in the cubic term.
- ⛔ **RAM is our binding constraint, and the supercell route is RAM-hostile.**  Grids scale with cell
  volume, matrices as \f$(N_qN)^2\f$; the tracker's own estimate says **the 32-atom MnO supercell is
  UNTESTED**.  A q-mesh runs one small cell at a time and fits trivially.  On 14 GB this may be the
  difference between slow and *impossible* — which is a fact to MEASURE, not to reason about.
- ⛔ **A supercell perturbation breaks the symmetry we lean on.**  Our imposed runs fold the Becke mesh
  11.55× on MnO AFM-II (12 ops); an isolated \f$\alpha\hat P_I\f$ destroys that by construction.  The q-mesh
  route keeps each primitive-cell run's symmetry.
- ★ **And the q-mesh is embarrassingly parallel over q** — which is the k/irrep parallel axis we already
  own a row for.

★★ **k-PARALLELISM IS ON THE CRITICAL PATH FOR EVERY ROUTE, (b) INCLUDED.**  Not just for (c)/(d): the
ACBN0 estimator consumes converged orbitals, we measured **~10 % k-sensitivity** in MnO's U between Γ and
2×2×2, and a trustworthy production U therefore needs a k-mesh SCF per composition.  §2 row **KP** is parked
"*BY AGREEMENT not now … it should stay on the list until we are suitably embarrassed*" — the parking
condition was "no multi-k rows exist to measure the payoff on".  **This plan is the embarrassment.**  What
KP needs is already scoped there: the pre-warm exists, and one k-dependent write remains inside the loop
(`tDynamic_HT_Imp::GetMatrix`'s `mutable CacheMap itsCache` keyed by `Irrep`, plus `itsByL/itsByLSeen`).
⚠ **And there is a PREREQUISITE DEFECT**: the shifted-MP fold lowers Si by 1.02 mHa (§4, found 2026-09-20)
— "a 1 mHa fold error on the fractional-k mesh is a sampling defect that every future k-mesh run inherits;
fix before KP-1".  Multi-k U values inherit it too.

**The measurement this gate asks for** (cheap, and it settles an assumption the tracker has carried
untested): run **one 32-atom MnO 2×2×2 supercell SCF at Γ**, logging wall and peak RSS.  That single number
says whether route (d)-by-supercell is available to us at all, or only in principle.

**Then, and only then: the screened kernel.**  A new integral type on `BareCoulombSource` (the pseudo-wall
pin allows exactly this), with \f$\varepsilon\f$ computed — Thomas–Fermi first as the cheap bound, reading
`qs_linres_polar_utils.F` before committing to anything more.

**The Oct 6–20 long runs** (only what the gates have justified): the three compositions × the outer loop,
plus a NiO re-run as the control once the magnetic-robustness fix exists, plus the k-mesh arms gate 4 says
we need (Γ-only was ~10 % on MnO's U) — and, if gate 4's 32-atom measurement says the supercell route fits,
the calibration LR run that turns gate 3's one solid oracle point into several.  ⛔ `scripts/memsafe -p`,
never bare: these are unattended.

---

## 5. Open questions (write the answer here when it is earned)

- Does **J** transfer the way U does?  ACBN0 gives J for free; hp.x does not give us a J to check it
  against.  The atomic limit (\f$J\approx1\f$ eV for 3d) is the only oracle we have — and our bare
  \f$\bar J\approx7.5\f$ eV is NOT Hund's J (it carries eq 13's self-terms), which is why only
  \f$U_{\rm eff}=\bar U-\bar J\f$ is quotable.
- Is the per-site-oxidation-state assignment stable when two Mn sites are crystallographically equivalent
  but electronically inequivalent (charge ordering)?  That is a symmetry-breaking question and the
  imposed-symmetry machinery has an opinion — check it does not average the two.
- If a computed \f$\varepsilon\f$ needs DFPT, we will have BUILT the machinery that makes a DFPT-based (d)
  nearly free — hp.x's own method with \f$\alpha\hat P\f$ in place of the electric field, CP2K's
  `qs_linres_*` as the Gaussian reference.  ⇒ **(b) and (d) share most of their cost**, which makes the
  choice between them a hedge rather than a bet.  `doc/OpenWork.md` step 5 item 5's "never build DFPT for
  this" was costed against QE's plane-wave implementation and assumed DFPT would be built ONLY for U; if (b)
  forces it anyway, reopen that ruling rather than inherit it.
- The spherical atom resolves a partially-filled degenerate shell by picking orbitals, not by occupying the
  shell uniformly — energetically converged but symmetry-broken.  Harmless for a seed (spherically averaged
  anyway); worth a thought for the **+U atomic radial**, which takes "the lowest occupied l orbital" and on a
  broken shell that need not be the spherical average.  Measure before caring: the captured-norm line already
  prints (NiO 0.9999, MnO VA 0.991).
- GGA before any value comparison with the PBE literature (the paper's 7.63/3.0 and the 4–7 eV range are
  both PBE).  Still open, still gating the *value* comparisons, not the *method* work.

---

## 6. UNFINISHED — started and not completed (distinct from §5, which is things we do not KNOW)

⚠ **These will rot silently if nobody looks.**  Each says what exists, what is missing, and where the
evidence is.  A fresh session should clear or re-park them before starting new work.

1. ✅ **RESOLVED 2026-09-23 — the 11 % bare-\f$F^0\f$ "discrepancy" was a units mismatch, not a defect.**
   PySCF's plain \f$(2l+1)^{-2}\sum_{mm'}(mm|m'm')\f$ (24.29 eV) and our banner's eq-10 \f$\bar U_{\rm bare}\f$
   (27.3 eV) are different quantities BY DEFINITION: eq 10's pair-count denominator (10c) excludes
   same-orbital/same-spin self-pairs while its numerator does not, so \f$\bar U_{\rm bare}\f$ is inflated
   above the plain average by a factor set by the spin-resolved occupation alone.  Verified two ways: (a)
   analytically, for a diagonal density matrix with uniform-within-spin \f$N_m\f$; (b) against a live
   `gpwprobe nio` run's own bare \f$N_\uparrow/N_\downarrow\f$, predicting \f$\bar U_{\rm bare}\f$ to
   0.1–0.5 %.  Full detail in gate 3's caveat (3) above.  **No fix needed anywhere** — our ERI4Block
   integrals and PySCF's agree; gate 3's ω values were already built on the correct (plain-F0) object and
   are unaffected.
2. ⛔ **NOTHING CHECKS THE BASIS AGAINST THE PSEUDOPOTENTIAL** (added 2026-09-23 — this one had been said
   aloud and never written down, which is how it nearly got lost).  A run declares its PP variant through
   `SolidCalcOptions::species` (`{"Li",3}`) while the basis comes from a `.bsd` file whose blocks are keyed
   by ELEMENT only.  So a q3 run can be handed the q1 `LI` block — a basis built and validated for a
   one-electron valence, describing three — **silently**, with no diagnostic anywhere.  It is latent today
   only because every element in `valence_lowq_*` happens to be the variant its runs use.  The shape of the
   fix: the `.bsd` files already carry rich headers, so a machine-readable per-element provenance line
   (invisible to a Gaussian94 reader) lets the factory assert each block's q against the run's declared
   valence, and THROW on a mismatch.  Same family as `CleanupCandidates.md` **D-SEED1** (the seed library
   keyed on (Z, functional)), and it becomes live the moment `valence_semicore.bsd` exists.
3. **Li q3 is validated but uncommitted.**  `--q 3 --shell 0:8:0.05:60`, converged, E = −4.23584 Ha, gap
   0.107 mHa.  It cannot share a file with q1, so it needs `valence_semicore.bsd` + a `BasisSetData` enum
   value and its two map entries.  Left unwired deliberately (no consumer yet, and a dead file is worse than
   a recorded command) — the command is in the `.bsd` header and in §4's prerequisites.  ⇒ item 2 above
   should land WITH it, not after.
4. **ABINIT Ni-3d `ucrpa` crashes with a Fortran integer overflow.**  Same deck as the working O-2p run
   (`~/Code/abinit-runs/ni3d_noU/`), same cell, same bands (11-26) for the Wannier construction and
   screening exclusion — only the correlated orbital swapped (l=2/projector=5 on the Ni sites instead of
   l=1/projector=3 on the O sites).  Crashes at the same point O-2p gets past cleanly: the start of
   optdriver=4's k-point loop ("begin of Ucrpa calc for a nb of kpoints of: 64").  Not root-caused.  Not on
   the critical path — Ni-3d already has hp.x's matched-PP oracle (5.27 eV, §4 gate 3 table); this run was
   only ever going to be a second opinion, unlike O-2p where it turns a bound into the only value we have.
   Evidence: `~/Code/abinit-runs/ni3d_noU/run.log`.

---

## 7. Literature hunt — candidates for a physics-based screening replacement (UNVETTED, 2026-09-25)

★ **Not findings — a queue for the user to sift** (web search only, no paper read in full).  The ACBN0
paper itself is arXiv Oct 2014 / PRX Jan 2015 — eleven years of follow-on work exists and none of it has
been checked against what we need: a REPLACEMENT for eq 12/13's Mulliken-type \f$\bar N\f$ renormalisation
(§ this session's finding: it is a basis-completeness artefact, not a screening model) with something
whose screening is a real, basis-independent dielectric quantity.  User is downloading and annotating
inline; this list is the starting set, not the final one.

**Most promising — direct hits on our exact problem:**
- [Comparative analysis of methods for calculating Hubbard parameters using cRPA](https://arxiv.org/abs/2503.11142) (Phys. Rev. B, May 2025) — systematically compares cRPA projection/Wannierisation schemes specifically for **entangled bands**, exactly NiO's Ni-3d/O-2p hybridisation problem (§4 gate 3's model-convention gap); benchmarks on LiMO₂ (M=V–Ni) and SrMO₃ (M=Mn,Fe,Co) — SAME element family as our own MnO/NiO.
- User comments: This paper defines and compares three methods of defining the polarizability function: 1.
Band method, 2.Disentanglement method, 3. Weighted method.  As such one would hope that in the conclusions section they would recommend one of these as being superior ... no such luck!  There are no conclusions of any sort on the "Conclusions" section, just a summary of what they did, and suggestiosn for future work.
  However the paper contains a substantial amount of calculated results for U/J accross SrMO3 and LiMO2 series, for us to compare with.  Methods 2&3 are based on Wannier functions.  They use https://github.com/wannier-developers/wannier90.  We can download an install if useful.
  Ultimately everything boils down choosing the bands (energy window) from which to compute the Wannier functions.  d only, d + Ox-p, d-frontier only etc.  As with Mulliken pop analysis we are over interpreting a single determinent approximate wave function and its defining basis set. Strictly speaking there is no such thing as "d-band", but there is a "mostly-d-band" but only within the context of our product WF and its associated basis set.

- [Bridging constrained random-phase approximation and linear response theory for computing Hubbard parameters](https://arxiv.org/abs/2505.03698) (2025) — connects cRPA (our ABINIT route) and linear-response (our hp.x route) methodologically; could directly bear on why our two oracle types disagree in scale.
- User comments: Yes this paper does exactly what the title says.  They are able to get agreement between cRPA and LRT "using well-defined Wannier projectors not only allows for a systematic comparison between LRT and cRPA (and potentially other methods to calculate U ), but also offers greater transferability across  different implementations."
- ★★ **READ IN FULL 2026-09-25 — and it redirects A2, not just explains it.**  Carta, Timrov, Beck & Ederer
  formally bridge LRT and cRPA (their Eq. 5) for an ISOLATED set of bands: once you account for (1) cRPA
  typically dropping the xc-kernel response that LRT naturally includes, and (2) cRPA's coarse-graining to
  a purely monopolar response discarding excitation channels INSIDE the interacting subspace that LRT
  keeps, the two agree to a few percent (KCuF₃ Cu-3d: 10.01 vs 9.97 eV).  **But for an ENTANGLED
  interacting/screening split — their "d-only" case — cRPA becomes ambiguous and collapses to an
  unphysically small U while LRT "remains largely unaffected": Sr₂FeO₄ Fe-3d gives U_cRPA = 0.42 eV vs.
  U_LRT = 6.94–7.29 eV, a 16× gap, SAME orbital, SAME material, from the window/method choice alone.**
  ⇒ **NiO's Ni-3d/O-2p complex (bands 11–26, no clean separation, per this doc's own fatbands finding) is
  exactly their "entangled" case.**  Our ABINIT O-2p cRPA number (≈1.2 eV, dp-dp model) is therefore
  SUSPECTED of being this SAME cRPA-in-a-hybridised-subspace pathology, not new physics about O-2p
  screening — and no `ucrpa_bands` convention search fixes an intrinsically ill-posed calculation.  Their
  own conclusion points the other way: trust LRT in the entangled regime.  **⇒ REDIRECT: the next O-2p
  action is extending hp.x (matched-PP LRT, our existing trusted oracle) to O-2p, not further ABINIT model
  sweeps** — added to §0 as A2b.  Mechanism, briefly: when D and R overlap in energy, screening channels
  that physically belong inside the correlated subspace get misattributed to the screening subspace,
  driving cRPA's U artificially low — the SAME shape of error as A1's Mulliken lesson (a quantity that
  looks like it measures physics but is actually measuring how the method's own bookkeeping handles a
  subspace split/basis choice).  Their fix for comparing methods at all was not dissolving the window
  choice but forcing BOTH methods onto ONE EXPLICIT, SHARED projector (Wannier) — pin 23's 2026-09-25
  addendum (doc/Pins.md) generalises this: the manifold/window is an INPUT, never a derivation, the same
  ruling pin 23 already makes for site/shell/irrep, one level up.  User: happy to run every material these
  papers use and to support BOTH LRT and cRPA natively (DIP, behind the existing `HubbardUEstimator` face)
  — filed as a `doc/OpenWork.md` §2 feature row, not a near-term build (each needs real new capability:
  LRT a perturb-and-respond mechanism, cRPA a χ0 in a product basis).

**On ACBN0's basis dependence specifically (this session's finding, independently):**
- [Orbital-Resolved DFT+U for Molecules and Solids](https://pubs.acs.org/doi/10.1021/acs.jctc.3c01403) (JCTC, 2023/2024, arXiv:2312.13580) — explicitly compares Mulliken vs Löwdin-orthogonalised projectors for the renormalised occupation, reports Löwdin improves self-consistency stability (we already made the same Löwdin-not-Mulliken choice, `Hamiltonian.C`'s design note item 3 — worth checking whether they also diagnose the basis-completeness failure mode).
- User comments: We processed this paper before in a different session durind DFT+U planning.  That is where the idea of diagonalizing the orbital occupation matrix comes from. They reference ~/Code/timrov2018-DFPT.pdf for the DFPT method.  Some important quotes: 1) "Recalling that the main motivation of Hubbard U corrections lies in the mitigation of local SIE (self interaction errors) through recovery of PWL (piecewise linearity) of the total energy, the Hubbard manifold should contain those and only those states
that substantially contribute to the former. Oftentimes, self-interaction occurs in partially occupied d and f shells due to their high electron count and localization; hence, these are the
traditional targets of Hubbard U corrections. Nevertheless, self-interaction can also manifest itself in s and p shells, ..." 2) "After all, the correction of all magnetic quantum orbitals within a given shell using the same scalar U parameter is inherently a simplistic approximation."

- [Pseudo-hybrid density functional ACBN0 for Hubbard U correction in a numeric atom-centered orbital basis](https://arxiv.org/abs/2609.12198) (2026, very recent) — NAO basis, i.e. the same "how much does the projection basis distort U" question in a different localised-basis code.
- User comments: I think paper uses ACBN0 as is, without acknowledging or addressing the shortcomings we have identified for that method.  They use PBE and SCAN xc functions which we don't have working yet.  But if we ever want band gaps, magnetic moments and U values to compare with they do present some good tables of numbers for Cr2 O3 , Cu2 O, CuO, MnO, NiO, and CoO.  

**DFT+U+V (intersite) follow-ons from the ACBN0 lineage:**
- [Efficient First-Principles Approach with a Pseudohybrid Density Functional for Extended Hubbard Interactions](https://arxiv.org/abs/1911.05967) (2019) — the ACBN0→ACBN0+V extension (intersite Hubbard V), likely by overlapping authors.
- [DFT+U+V is equivalent to DFT+U with density-dependent hybridized projectors](https://arxiv.org/abs/2607.18071) (2026) — theoretical reformulation, very recent.

- User comment:  We need to consider 1) DFT+U+V, 2) DFT+U+J, 3) 1&2 with Resolved U (and J?, and V?) but all from the standpoint of LR-cDFT calculated throught DFPT (~/Code/timrov2018-DFPT.pdf).

**Screened-kernel FORM (relevant to route (b)'s actual kernel construction, not just the U value):**
- Analytical treatment of the Yukawa screened Coulomb interaction in a plane-wave basis (2025; ScienceDirect/ADS) — closed-form matrix elements for a Yukawa/Thomas-Fermi-screened kernel in a PW basis; the Gaussian-basis analogue is what `BareCoulombSource` would need for route (b)'s A4.
- Calculation of Effective Coulomb Interaction for Pr³⁺, U⁴⁺, UPt₃ (cond-mat/9501111, classic) — early Yukawa/Thomas-Fermi-screened Slater-integral fit; background, not current.

**Background / reviews:**
- Hubbard-corrected DFT energy functionals: the LDA+U description of correlated systems (Himmetoglu, Marzari, Cococcioni; Int. J. Quantum Chem. 2014, arXiv:1309.3355) — the standard broad review, pre-dates ACBN0-specific critique but frames the double-counting/screening landscape it sits in.
- DFT+U within the framework of linear combination of numerical atomic orbitals (arXiv:2202.05409) — same LCAO-basis-dependence territory as our finding, different code family.

**A different philosophy, for contrast (probably not what we want, but worth knowing it exists):**
- Machine learning the Hubbard U parameter in DFT+U using Bayesian optimization (npj Comput. Mater., 2020) — fits U empirically against a target property rather than deriving it from screening physics; the opposite direction from route (b).
- user comment: If we do decide to support emprical methods (tune parameters {a,b,c...} in order optimize agreement measured properties {A,B,C,...}) I would to plan it in a much wider context than just tuning U.

---

---

# OpenWork §1 (DFT+U step 5) as of 2026-10-01

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



---

## Pins.md pins 23-25 verbatim as of 2026-10-01 (archived when Pins.md was condensed)

## 23. +U is ORBITAL-RESOLVED — U is a vector over (site, shell, site-group irrep); the manifold is an INPUT, never Mn-d by assumption

Ruled 2026-09-16 (user, on Macke et al. JCTC 2024 and ACBN0).  In the eigenbasis of the site occupation
matrix \f$E_U=\sum_i\tfrac{U_i}{2}\lambda_i(1-\lambda_i)\f$; the shell-averaged Dudarev form is the special case
\f$U_i=U\f$ — the same shape as pin 5 (unpolarized is the ζ=0 collapse), applied to +U.  The (t2g, e_g)
split is the site-point-group irrep decomposition, so the labels are fixed by symmetry; eigenvalue
tracking is only for a site symmetry lower than the split.  **Why it is a ruling and not a preference:**
shell-averaging suppresses intrashell screening (perturbing t2g and e_g together zeroes the channel that
screens them: FeS₂ U 7.37 → 3.29/2.16 resolved), and the WRONG manifold is worse than the wrong U — the
correction that opened β-MnO₂'s gap was on **O-p_z**, not Mn-d, and correcting FeS₂'s hybridised e_g at all
broke its structure.  User: *"I have seen other examples where O played an unexpected role in TMOs."*  ⇒ the
term takes a LIST of (site, shell, irrep, U); no code path may assume the Hubbard atom is the transition
metal.  Projector = Löwdin OAO on the site block.  U values are never hand-tuned in production (pin 12):
ACBN0-style from our own on-site ERIs, checked against QE `hp.x`.  Record: `doc/OpenWork.md` §1 step 5.
**Addendum 2026-09-21 (increment 2, earned by a wrong table):** the labels have TWO groups and neither is
the cell's.  The SITE group is the declared decoration's Shubnikov stabiliser (σ=None) — the order splits
t2g → a1g + e_g and the labels must see it.  The PARENT group that names a site level "e_g < t2g" is the
point group of the site's **coordination environment** (`Lattice_3D::SiteEnvironmentRotations`), NOT the
(super)cell's grey stabiliser: on the rhombohedral AFM-II MnO cell the latter is D_3d (12 ops) with or
without decoration and names nothing — measured, after the tree had asserted O_h for a week.  And the
occupation matrix is NEVER symmetrised: symmetry NAMES the eigenvectors of the density's own n (isotypic
projectors, `purity` printed), it does not edit them — a free run's broken symmetry must keep its own
occupations, and the functional must stay dE/dD.  Inside a degenerate cluster the eigenbasis is rotated
to the projectors (n is unchanged); inside a NEARLY degenerate one (four λ≈0.999 on a full majority
shell) the names are ill-conditioned by nature — that is Macke's tracking problem, and the printed
`parentage` says so rather than hiding it.

**Addendum 2026-09-25 (user, reading Carta, Timrov, Beck & Ederer, arXiv:2505.03698 — the manifold
question is not a defect to engineer away):** the same "manifold is an INPUT" ruling extends from WHICH
site/shell/irrep to WHICH ENERGY WINDOW, and for the identical reason.  Carta et al. formally bridge cRPA
and linear-response U for an ISOLATED set of bands (their Eq. 5: the two agree once the xc-response and
the coarse-graining's dropped intra-subspace channels are both accounted for) — but for an ENTANGLED
interacting/screening split, cRPA becomes ambiguous and can collapse to an unphysically small U while LRT
"remains largely unaffected" (their Sr₂FeO₄ Fe-3d: 0.42 eV cRPA vs. 6.94–7.29 eV LRT, SAME orbital, SAME
material, a 16× gap from the WINDOW CHOICE alone).  NiO's Ni-3d/O-2p complex (bands 11–26, no clean
separation — `doc/HubbardUPlan.md` A2) is exactly their "entangled" case, so our own O-2p cRPA number is
suspected of being this same pathology, not new screening physics — a finding earned by the paper, not by
our own retraction discipline this time, but the same shape as A1's Mulliken lesson: **there is no
basis-independent "the d-band" any more than there is a window-independent "the correlated subspace."**
A grad student forced to draw an energy-window cut on a DOS plot is not doing something wrong — it is the
only thing that CAN be done — and Carta et al.'s own fix was not to dissolve the choice but to make it
ONE EXPLICIT, SHARED PROJECTOR (Wannier) so two methods could even be compared.  ⇒ two consequences, not
one: (i) the manifold/window stays an input, never a recommended default masquerading as a derivation;
(ii) **the code should RECOMMEND a window when a natural one exists** (a DOS minimum, a projected-character
gap) rather than only accept one — surfacing the choice is not the same as making it, and a tool that can
show *why* a cut is natural is strictly better than one that is silent.  Architecturally: LRT and cRPA
belong in the tree as ADDITIONAL CONCRETE STRATEGIES behind the same abstract estimator face `ACBN0`
already sits behind (`HubbardProjection`/`HubbardUEstimator`, increment 3) — DIP, not a special case for
each — tracked as a `doc/OpenWork.md` §2 feature row, not decided here.

## 24. A linear-response number is conditioned on the state it linearises about — and same-site LRT is ill-posed for Pauli-saturated shells

An hp.x (or our DFPT) U is dn/dα at a SPECIFIC reference state: NiO's 5.267 eV is U_LR(U_in = 3 eV), not a
property of the material, and it differs at U_in = 0 (5.434).  Quote U_in beside every response value and read
the `HUBBARD` block of a deck rather than assuming zero.  And when χ₀ ≈ χ (MnO d⁵ half-filled, ZnO d¹⁰ closed
shell) the same-site U = (χ₀⁻¹ − χ⁻¹)_II is a difference of near-equal small numbers: it is not an oracle, and no
projector choice rescues it (`doc/Records/HubbardUHistory.md`, A6 table).  Also: an iteration CAP is not a
verdict — three wrong conclusions in this project were a capped run read as "does not converge".

## 25. Two comparisons, never merged: projector completeness vs screening — and no target may be derived from the method under test

Ours ÷ published-ACBN0 (factors 1.8–2.8) measures PROJECTOR COMPLETENESS (their minimal PAO keeps ~60 % of the
norm), not screening.  Only ours ÷ an INDEPENDENT oracle (matched-pseudopotential hp.x, never ACBN0 or an
ACBN0-derived target) tests the screening model; label every oracle row matched-PP or different-PP.  ACBN0's N̄²
renormalisation vanishes as the basis completes, so it is not a screening model, and an oracle gate built on an
ACBN0 target is circular (retracted 2026-09-23; screened-ACBN0 refuted 2026-09-25).

## A7 cross-code gate vs CP2K CK-alpha (2026-10-21)

**Verdict: PASSED on TiO₂ and FeS₂; SrVO₃ at Γ PARKED.**  Ours = `rundeck decks/campaigns/CrossCode/<m>.json` (`7f411887`; SrVO₃ deck
fixed the same day, below): Γ only, the CP2K runs' cell and basis file (`BasisSetData/valence_lowq_tm.bsd` via `basis.file`), LOWDIN
all-shells manifold on site 0 (35 functions: 7 uncontracted d primitives), U_in = 0, kT 0.005, α = 1e-3 Ha.  CP2K = the 2026-10-05
CK-alpha pipelines (`~/Code/qchem6-runs/batch/lib_ckalpha_generic.sh`, logs `batch/logs/07–09_*.log`).  Records
`~/Code/qchem6-runs/{TiO2/TiO2_rutile.r003, FeS2/FeS2_pyrite.r002, SrVO3/SrVO3.r007–r010}.json`.  χ in Ha⁻¹.

| | d trace (α=0) | χ₀ | χ (LR) | χ (FD) | χ/χ₀ | U eV |
|---|---|---|---|---|---|---|
| TiO₂ ours | 2.981 | −13.917 | −3.98688 | −3.98567 | 0.2865 | **4.870** |
| TiO₂ CP2K | 2.987 | −13.993 | — | −3.99651 | 0.2856 | **4.864** |
| FeS₂ ours | 7.692 | −33.270 | −3.98700 | −3.93534 | 0.1198 | **6.007** |
| FeS₂ CP2K | 7.692 | −33.457 | — | −3.98992 | 0.1193 | **6.007** |
| SrVO₃ ours, unpol r009 / r010 | 4.096 / — | −143.90 / −134.82 | −5.915 / −5.603 | −5.771 / −5.671 | 0.0411 / 0.0416 | 4.412 / 4.655 |
| SrVO₃ CP2K (RKS) | 4.098 | −144.28 | — | −5.89684 | 0.0409 | **4.426** |

Agreement on the insulators: χ₀ 0.5–0.6 %, χ 0.07–0.3 %, U to 6 meV (TiO₂) and <1 meV (FeS₂); same-state checks pass (d traces,
FeS₂ ±α SCFs ΔE = ±7.69e-3 Ha = α·n_d).  ⇒ The TiO₂ +5 % and FeS₂ −20 % gaps CP2K-vs-hp.x are SETUP (Γ vs k 4³/q 2³, LOWDIN-all vs
ortho-atomic), not method.  Cost: TiO₂ 21.6 min wall / 715 MB (6 threads, LR + FD); FeS₂ 38.9 min / 2.6 GB (6 threads; SCF 28 iterations, LR 15 kernel applications).  OPEN, not chased: FeS₂'s OWN FD χ is 1.3 % below its LR (TiO₂: 3e-4) although the ±α SCFs stay on the ground-state branch and the restore returns to it (ΔE 5e-8 Ha); CP2K's FD agrees with our LR to 0.07 %, so the gap is in our FD path, not the response (candidates: the iteration-0 refill below, the FD SCF tolerance at Δn = χα ≈ 4e-3).
**RESOLVED 2026-10-22 (anchor-moving batch):** rerun of `decks/campaigns/CrossCode/fes2_Fe3d_G.json` after the carried-density refill
(±α stages refilled under the run's kT) and D-MINDRHO (the honest convergence gate; the old Kerker/Pulay gate stopped 50–120× early):
χ_FD −3.98893 vs χ_LR −3.98723 (+0.04 %; was −1.3 %), CP2K's FD −3.98992 (now 0.025 % apart); ground state and LR unchanged
(Etot −164.9113610, 28 iterations, U 6.007 eV); ±α SCFs 33 iterations each.  Which of the two changes closed it is not separated (an
item-1-only rerun would).  Run data `~/Code/qchem6-runs/FeS2/fes2_Fe3d_G_newgate_20261022.{log,json}`.

**SrVO₃ — why parked (user 2026-10-21: look for corner cases where U/χ disagree, do not get bogged down in convergence heuristics).**
(1) The committed deck had `multiplicity 1`; on 25 e that is a 13↑/12↓ DOUBLET (`twoS = Nelec % 2`), which at Γ polarizes (V +2.60,
O −0.53 each, net +1.0; 0.16 Ha BELOW the unpolarized state; U 1.68 eV, χ/χ₀ 0.93) and is multi-stable under the FD α (+α landed 0.1 Ha
up, V moment 0.95).  `spinsShareFermi` (moment free) also polarizes (m_V 1.0–1.6, unsettled at 56 iterations).  CP2K ran RKS, so the
deck is now `multiplicity 0` (ζ=0: one Spin::None channel, 25 e, smeared).  (2) Unpolarized, Γ SrVO₃ puts 1 e in a 6-fold
t2g×spin level at E_F; four converged fixed points within 2e-3 Ha give χ −5.60…−5.91, U 4.41–4.65 — CP2K's 4.43 is inside the band,
and no sharper statement is honest.  (3) A restart lands on a DIFFERENT fixed point than the one it left (r010: the "restored"
α = 0 state is 1.8e-3 Ha below the run's own ground state), so its FD χ cannot test linearity.

**Found on the way (OpenWork §3 row "A carried density is re-filled at kT = 0"):** every SCF stage built from a carried, converged
density does an iteration-0 Fock-diagonalize-and-fill on the occupation policy's DEFAULT state (integer aufbau, kT = 0, no MOM;
`src/SCFIterator/Imp/SCFIterator.C` `Initialize`), and iteration 1 starts from THAT, not the carried density.  First FD iterate:
TiO₂ (gap 1.8 eV) Δρ 4e-4; FeS₂ (gap 0.73 eV, kT 0.136 eV) Δρ 3.4e-3, E 0.93 Ha off, swinging to −116 Ha before returning; SrVO₃
(metal) Δρ 2.7e-2, E 14 Ha off.  Pulay at relax 1.0 (the deck's ground-state params, passed straight to the FD re-converge) then
amplifies it.  The LR numbers never pass through a restart, so the gate stands on them.

## The cross-code gate re-run after the pin-29 image-list fixes (2026-10-24)

OpenWork §3 ORDER step 1 (the CHEAP tier).  Same decks (`decks/campaigns/CrossCode/{tio2,fes2}_*_G.json`, TiO₂ b6 6×6
`~/Code/qchem6-runs/TiO2/basis_scan/b6_TiO_postpin29.json` = `b6_TiO.json` with the retired `minDeltaRho` → `minDeltaRhoPerElectron 1e-6`),
binary = `openwork-sweep` at `30a41775` (image fixes + anchor batch + N4 default; NOT seed-tiebreak).  Logs `*_postpin29_20261024.log`
in `~/Code/qchem6-runs/{TiO2,TiO2/basis_scan,FeS2}/`; records `TiO2_rutile.r006`, `basis_scan/TiO2_rutile.r002`, `FeS2_pyrite.r003`, the
control `FeS2_pyrite.r004`.

| | E_tot Ha | χ₀ | χ LR | χ FD | U eV |
|---|---|---|---|---|---|
| TiO₂ LOWDIN, banked 10-21 | −71.1350613 | −13.917 | −3.98688 | −3.98567 | 4.870 |
| TiO₂ LOWDIN, 10-24 | −71.1350964 | −13.9727 | −3.99553 | — | **4.863** |
| TiO₂ CP2K | | −13.993 | | −3.99651 | 4.864 |
| TiO₂ b6 6×6, banked 10-21 | −71.41653 | | | | Ti 3.051 / O 6.735 (Ti 1×1 5.125) |
| TiO₂ b6 6×6, 10-24 | −71.4165332 | | | | Ti 3.0504 / O 6.7345 (Ti 1×1 5.1265) |
| FeS₂, 10-22 (anchor tree; XC from ρ_mix) | −164.9113610 | −33.2807 | −3.98723 | −3.98893 | 6.0070 |
| FeS₂, 10-24 CONTROL (`xcCuspDeficit=false`) | −164.9113609 | −33.2793 | −3.98720 | −3.99322 | 6.0070 |
| FeS₂, 10-24 (cusp deficit, the default) | −164.9114008 | **−32.7800** | −3.98059 | −3.98469 | **6.0059** |
| FeS₂ CP2K | | −33.457 | | −3.98992 | 6.007 |

**The V_short fix did NOT move these cells.**  FeS₂ is the clean witness: the control (today's code, the 10-22 XC route) reproduces the
10-22 run to 1e-7 Ha, χ₀ 4e-5 relative, U to 2e-5 eV.  Rutile: E within 0.04 mHa, U within 2 meV (an upper bound -- the 10-21 runs predate
the anchor batch too).  Pin 29's own prediction: the defect grows with the home-atom separation, so compact cells see little of it --
MnO AFM-II's −61 mHa is the other end.  ⚠ So a small move here does NOT clear NiO: rocksalt AFM-II is MnO's own cell type.

**What did move FeS₂ is the N4 cusp deficit (pin 28): χ₀ −1.5 % (0.53 % → 2.0 % from CP2K), χ −0.17 %, U −1 meV** -- the banked gate was
measured with XC on the bare ρ̃_mix (the log's "XC rho source: rho_mix", a pin-28 violation), and moving XC to ρ[D] moved χ₀ AWAY from
CP2K while U, a ratio, held.  Rutile's χ₀ moved TOWARD CP2K (0.5 % → 0.15 %), so it is not one sign.  Open (OpenWork §3 re-run row):
CP2K's GPW evaluates XC on its REGULAR-GRID density, not on an atomic quadrature, so a 1–2 % χ₀ difference may be the XC-density
convention; the discriminator is FeS₂ on our raster (pair) route, `--set solid.policy.beckeXC=false` -- collocated ρ[D] on a uniform grid, CP2K's way.

**RESOLVED 2026-10-25 -- it IS the XC quadrature** (`FeS2_pyrite.r005`, log `fes2_Fe3d_G_raster_20261025.log`; trace `shadow` / `Pul+K` on all
135 rows, pin-28 clean): raster E −164.9077261 Ha (3.7 mHa above Becke), **χ₀ −33.3333 (0.37 % from CP2K; Becke + cusp 2.0 %)**, χ LR
−3.98767 (0.06 % from CP2K's FD), U 6.0075 eV; own FD χ −3.97862 (0.23 % from its LR; Becke 0.10 %).  Like-for-like with CP2K (ρ[D] on
a uniform grid at the density cutoff), χ₀ agrees; the Becke atomic grid resolves the core density a uniform grid cannot, and χ₀ -- not
U -- feels that at the 2 % level.  The fitted-ρ̃ control sat near CP2K for the same reason (band-limited).  U is route-independent:
6.0059 / 6.0070 / 6.0075 eV on Becke+cusp / Becke-fitted / raster vs CP2K 6.007.  ⇒ a χ₀ comparison against CP2K is a RASTER-route
comparison; a Becke χ₀ is compared with CP2K only through U.

**NiO at the cell level, same day** (`rundeck decks/campaigns/NiO/nio_q0_step3_newroute.json` -- the exact banked deck: free, k222, q = 0,
U_in 3 eV frozen, ortho-atomic, Kerker-Pulay; restart from the saved FREE state, which re-converged in 11 iterations; binary `f3cce7ef`;
record `~/Code/qchem6-runs/NiO/NiO_AFM2.r008`, log `nio_q0_step3_newroute_postpin29_20261024.log`):

| NiO AFM-II, q = 0 | banked 2026-10-04 | 2026-10-24 |
|---|---|---|
| E_tot / E_nn / E_elec Ha | −106.2023644 / −95.7965 / −15.5775 | −106.1384955 / −95.7245 / −15.5856 |
| χ₀ Ni0 / Ni1 | −4.1181 / −4.0644 | −4.02564 / −4.02549 |
| χ Ni0 / Ni1 | −2.07002 / −2.05166 | −2.04332 / −2.04329 |
| U (2×2) Ni0 / Ni1 eV | 6.4502 / 6.4793 | **6.4677 / 6.4677** (1×1 6.5577) |

**Unlike rutile and pyrite, NiO's Hamiltonian moved** (E_elec −8.1 mHa; E_nn +72.0 mHa is the Ewald constant) -- rocksalt AFM-II, MnO's
cell type -- **and the defect's signature is in the response: the two Ni, equivalent by the AFM-II symmetry, had χ₀ 1.3 % and U 29 meV
apart; now 4e-5 and 2e-5 eV** (pin 29: equivalent sites got different V_short blocks).  The MEAN U moved +3 meV (6.4648 → 6.4677), the
mean χ₀ −1.6 % (not separated between V_short and the N4 cusp deficit, which moved FeS₂'s χ₀ by the same amount; a
`--set scf.xcCuspDeficit=false` control would split them).  So the banked NiO U numbers (cell 6.09-era loop, q-mesh 5.73) are not
re-run now: their mean moves by meV, and step (3) re-runs them anyway after the projector decoupling (2b) -- which changes exactly these
ortho-atomic projectors.  SrVO₃ (one V per cell, no equivalent-site pair to split) is skipped by the same rule.

## TiO₂ ortho-atomic at Γ vs hp.x at Γ (2026-10-21) — an OPEN corner case

Ours: `rundeck decks/campaigns/CrossCode/tio2_Ti3d_G_ortho.json` = the CP2K-validated Γ deck (same state: Etot −71.1350613 Ha, as
the LOWDIN run) with QE's ortho-atomic set (every atom's PSWFC: Ti 4S+3D, O 2S+2P as U = 0 spectators, Ti0 3d perturbed), multiplicity 0;
record `~/Code/qchem6-runs/TiO2/TiO2_rutile.r005.json`, 1 min 44 s.  QE: the archived `IntegrationTests/QE/tio2.{scf,hp}.in` with k 1³,
q 1³ (`~/Code/qe-runs/TiO2/tio2_G.*`, `mpirun -np 4`, 1 min 13 s); its UPFs' PP_CHI are OUR pseudo-atom orbitals (`CLIapps/gth2upf.C`),
so both codes project onto the same radials.  QE prints the Ti1 column only; by symmetry the 2×2 inversion gives hp.x's 5.0191 eV, the
1×1 scalar 5.026 (the matrix-vs-scalar difference is 7 meV, so our single-site χ is compared directly).

| Ti 3d, Γ, ortho-atomic | ours | QE hp.x | ours/QE |
|---|---|---|---|
| d occupation (both spins) | 2.405 | 2.208 | +8.9 % |
| χ₀(1,1) Ha⁻¹ | −14.404 | −10.969 (−0.40310 eV⁻¹) | +31 % |
| χ(1,1) Ha⁻¹ | −3.566 | −3.625 (−0.13321 eV⁻¹) | −1.6 % |
| χ₀(2,1), χ(2,1) eV⁻¹ | — | +0.02430, +0.00221 | |
| χ/χ₀ | 0.248 | 0.330 | |
| U eV (scalar) | 5.742 | 5.026 (2×2 5.019) | +14 % |

Almost all of the 0.72 eV is the χ₀⁻¹ term; not a projector SCALE (that would leave χ/χ₀ alone).  Lead: the ortho-atomic spectators are
truncated by `valence_lowq_tm.bsd` (captured norm Ti 4s 0.9675, O 2p 0.9762 on 4 p shells; Ti 3d 0.9952, O 2s 0.9971), so Ti 3d is
under-orthogonalised against O 2p — consistent with a larger occupation and a larger bare response.  The discriminator (same deck on
`VALENCE_LOWQ_VA`) could not run: that file has no Ti.  CP2K's LOWDIN-all has no spectators, which is why the CP2K gate could not see it.
Found on the way: `SiteEnvironmentRotations` threw on rutile's coplanar O environment (fixed `305450f4`), and `Gaussian94Reader` reports a
missing element as a garbage-byte angular-momentum label (✅ fixed the same day, `a833d5f1`: `FindAtom` throws naming element and file).

**Discriminator, same day — the O 2p spectator is most of it.**  `valence_lowq_va.bsd` cannot help: its O block IS the TM O block, and a
VA-style Ti would only make Ti s more diffuse (TM's header: in the 3d series the diffuse s hurts conditioning most).  Scratch basis =
TM + Ti s 0.05 + O p 0.18 (`basis.file`, `vet` on; scratch, not committed).  In rutile the vet TRIMS Ti s 0.05 (min eig S 6.2e-6 → 3.3e-5
after) — the free atom would have kept it — and KEEPS O p 0.18: O 2p captured norm 0.9762 → 0.9994 (Ti 4s stays 0.9675), min pivot
3.4e-3 → 9.0e-4.

| Ti 3d, Γ, ortho-atomic | TM | TM + O p 0.18 | QE hp.x |
|---|---|---|---|
| E_tot Ha | −71.1351 | −71.3499 | −71.4235 (−142.84694 Ry) |
| d occupation | 2.405 | 2.230 (+1.0 %) | 2.208 |
| χ₀(1,1) Ha⁻¹ | −14.404 | −11.785 (+7.4 %) | −10.969 |
| χ(1,1) Ha⁻¹ | −3.566 | −3.458 (−4.6 %) | −3.625 |
| U eV (scalar) | 5.742 | 5.561 (+10.6 %) | 5.026 |

The occupation lands on QE and the χ₀ excess drops 31 % → 7 %: an under-orthogonalised Ti 3d (truncated O 2p spectator) was most of the
corner case.  The larger lesson is the energy: ONE O p shell is worth 0.215 Ha, and TM sits 288 mHa above the plane-wave answer (74 mHa
with the tail).  CP2K uses the same basis, so the CP2K gate validated the METHOD on an incomplete basis, never the basis.  The remaining
+0.54 eV is not yet split between the Ti 4s spectator (vet-limited in rutile) and general basis incompleteness.

**Basis convergence scan, same day (user: "converge a TiO₂ basis toward the PW energy before chasing the last +0.54 eV").**  Every
variant = `valence_lowq_tm.bsd` + added uncontracted shells, `vet` on, judged in rutile (scratch decks/bases, not committed; runner
`mkvariant.py` / `runvariant.sh` in the session scratchpad).  ⚠ the density cutoff is AUTO (2·α_max), so variants ran on 72, 90 and
220 Ha grids; U is stable across all three.  QE: E = −71.42347 Ha (300 Ry; both codes Dirac + VWN5).

| variant (cumulative) | E Ha | gap mHa | d occ | χ₀ / χ Ha⁻¹ | U eV | vet |
|---|---|---|---|---|---|---|
| TM | −71.13506 | 288 | 2.405 | −14.404 / −3.566 | 5.742 | — |
| b1 = + O p 0.18 | −71.34993 | 73.5 | 2.230 | −11.785 / −3.458 | 5.561 | Ti s 0.05 TRIMMED |
| + O d 1.0 (on b1) | −71.34845 | 75.0 | 2.230 | −11.764 / −3.529 | 5.398 | trims TM's OWN Ti s 0.10 |
| + Ti p 0.3 (on b1) | −71.36460 | 58.9 | 2.233 | −11.829 / −3.520 | 5.431 | |
| + O s 0.15 (on b1) | −71.36717 | 56.3 | 2.226 | −11.632 / −3.477 | 5.486 | pivot 4.3e-4 |
| b2 = b1 + O s 0.15 + Ti p 0.3 | −71.37453 | 48.9 | 2.238 | −11.825 / −3.532 | 5.403 | |
| b3 = b2 + Ti p 0.9 | −71.38258 | 40.9 | 2.234 | −11.739 / −3.539 | 5.371 | |
| b3 + tight Ti (s 60, d 90) | −71.38294 | 40.5 | 2.235 | −11.766 / −3.543 | 5.368 | |
| b3 + tight O (p 20, s 45) | −71.39622 | 27.3 | 2.234 | −11.781 / −3.546 | 5.365 | |
| b4 = b3 + tight O + O d 2.0 | −71.39783 | 25.6 | 2.234 | −11.789 / −3.566 | 5.323 | pivot 4.1e-4 |
| b4 + O p 50, s 110 | −71.39775 | 25.7 | 2.233 | −11.762 / −3.562 | 5.325 | (grid 220 Ha) |
| b4 + Ti s 0.07 | = b4 | | | | 5.323 | Ti s 0.07 TRIMMED |
| b4, Ti 4s DROPPED from the ortho-atomic set | = b4 | | 2.242 | −11.799 / −3.568 | 5.320 | |

**Verdict: U is converged with respect to the basis, the grid and the Ti 4s spectator, at 5.32 ± 0.03 eV, against hp.x 5.026 (scalar;
5.019 2×2) — +6 %, carried by χ₀ +7.5 % (χ −1.6 %, occupation +1.2 %).**  The O 2p spectator WAS the 31 % → 7 % step; nothing since
moves χ₀ toward QE.  The vet is the right judge: it trimmed every diffuse Ti s a free atom would keep (0.05, 0.07; and TM's own 0.10
once O d 1.0 joined).  The energy gap (25.6 mHa) is basis in untried directions (O f, more O d, intermediate Ti p/d), not grid.
Not systematic with the other materials: SrVO₃'s χ₀ was SMALLER than hp.x's (−1.317 vs −1.426 eV⁻¹) — a "χ₀ excess" is not a
common sign.  Next suspects (not chased): what each code's bare response at Γ contains (hp.x iteration-1 Sternheimer vs our
Krylov over the basis' virtuals), and the ortho-atomic projector construction itself (QE's O^{-1/2} on its radial mesh vs ours).

**⚠ CORRECTION, same day — the "+6 % at a converged basis" verdict above was PREMATURE: the basis was converged in s/p/d, not in
ANGULAR MOMENTUM** (user: "PW has no l limit — do we need l=3 on Ti?").  Continuing the scan with f shells (all vet-clean):

| variant (cumulative on b4) | E Ha | gap mHa | d occ | χ₀ / χ Ha⁻¹ | U eV |
|---|---|---|---|---|---|
| b4 | −71.39783 | 25.6 | 2.234 | −11.789 / −3.566 | 5.323 |
| b5 = + Ti f 0.6 | −71.41028 | 13.2 | 2.220 | −11.517 / −3.586 | 5.226 |
| b5 + Ti f 1.8 | −71.41061 | 12.9 | 2.219 | −11.512 / −3.588 | 5.221 |
| b5 + O f 1.0 | −71.41403 | 9.4 | 2.214 | −11.357 / −3.584 | 5.196 |
| b5 + Ti f 0.25 | −71.41371 | 9.8 | 2.215 | −11.439 / −3.616 | 5.146 |
| **b6 = b5 + O f 1.0 + Ti f 0.25** | **−71.41653** | **6.9** | **2.211** | **−11.312 / −3.613** | **5.125** |
| b6 + O f 2.5 | −71.41664 | 6.8 | 2.211 | −11.315 / −3.614 | 5.125 |
| QE hp.x | −71.42347 | 0 | 2.208 | −10.969 / −3.625 | 5.026 (2×2 5.019) |

**Verdict: NO method disagreement.**  At b6 the occupation matches QE to 0.1 % and χ to 0.3 %; χ₀ is +3.1 % and U +2.0 %, and U
tracks the energy gap (~0.016 eV per mHa over the last steps) — extrapolating to zero gap gives ≈ 5.01 eV, on QE's 5.02–5.03.  The
corner case was (1) the truncated O 2p ortho-atomic spectator (+31 % → +7.5 %) and (2) missing f on Ti AND O (+7.5 % → +3.1 %).
My own prediction (a Hylleraas bound: a bigger virtual space can only RAISE |χ₀|) was wrong in sign — f moved |χ₀| DOWN through the
occupied states and the projector, the escape clause of that argument.  Also measured: the XC mesh is not it (b4 on Becke nR40/L29,
Uniform, Becke nR60/L35: χ₀ −11.789 / −11.742 / −11.743, U 5.323 / 5.326 / 5.331).  Ti g (l=4) could not be tested: three separate
failures (OpenWork §3 row "l=4 shells in a lattice basis").  LESSON: a Gaussian-vs-PW χ₀ comparison needs the basis converged in l
as well as in exponent range, and a CP2K match on the same basis cannot show it.  Cost: an f-bearing TiO₂ Γ run is ~10–20 min (5 threads).

**U(O 2p) on TiO₂ at Γ vs hp.x — PASSED, same day** (S1).  QE: `~/Code/qe-runs/TiO2/tio2_GO.*` (HUBBARD ortho-atomic, U Ti-3d AND O-2p
1e-8, k 1³, q 1³; hp.x inverts the full 6×6 Hubbard-site matrix).  Ours: the b6 basis, all six manifolds perturbed (deck indices 0, 1, 5,
7, 9, 11), ~15 min.  χ in Ha⁻¹ (QE eV⁻¹ × 27.2114).

| | ours | QE | Δ |
|---|---|---|---|
| occupation Ti 3d / O 2p | 2.211 / 4.937 | 2.208 / 4.931 | +0.1 % / +0.1 % |
| χ₀ (Ti1,Ti1) / (O3,O3) | −11.313 / −7.194 | −10.969 / −7.093 | +3.1 % / +1.4 % |
| χ₀ (Ti1,O4) / (Ti1,O3) / (Ti1,Ti2) / (O3,O6) | 4.134 / 0.975 / 0.588 / 2.250 | 4.057 / 0.837 / 0.661 / 2.215 | |
| χ (Ti1,Ti1) / (O3,O3) | −3.6135 / −2.3762 | −3.6249 / −2.3675 | −0.3 % / +0.4 % |
| χ (Ti1,O4) / (Ti1,O3) / (O3,O6) / (O3,O4) | 1.1075 / 0.4464 / 0.4906 / 0.1297 | 1.1068 / 0.4324 / 0.4905 / 0.1344 | |
| **U(O 2p)** eV | **6.735** | **6.758** | **−0.3 %** |
| **U(Ti 3d)** eV, 6×6 | **3.051** | **2.927** | +4.2 % |

The screened χ matrix matches element by element; the residual is χ₀'s Ti block (+3 %, the f-basis tail above), which the 6×6
inversion amplifies in U(Ti) (+2 % scalar → +4 %).  **And U(Ti 3d) is NOT a property of Ti alone: adding O 2p to the Hubbard set takes
hp.x's U(Ti) 5.019 → 2.927 eV (ours 5.125 → 3.051)** — the Ti–O off-diagonals enter the inversion.  Per-irrep O pz (Macke) is not
reported by our response yet (`doc/LinearResponsePlan.md` START HERE: orbital-resolved channels REPORTED).

**The k/q mesh vs hp.x — PASSED, same day (S1's last item).**  hp.x re-run at k 2³ / q 2³ (`~/Code/qe-runs/TiO2/tio2_k2.*`, 70 min on 4
ranks, 6 q): **U = 4.6373 eV — identical to its k 4³ / q 2³ 4.6368 (0.5 meV)**, so TiO₂'s k-mesh is converged at 2³ and the Γ → 4.637
change is the q-mesh (supercell size).  Ours: b6 basis, k 2³ / q 2³ FREE (D5: an imposed rutile group would reduce the mesh), kT = 0 (QE
`occupations='fixed'`; the pre-flight refuses a smeared q ≠ 0 response with per-k Fermi levels), Ti0 3d perturbed; record
`TiO2_rutile.r027` (scratch dir), 1 h 50 min on 12 threads, peak RSS 3.2 GB, 13–15 kernel applications per q.

| k 2³ / q 2³, home cell R = 0 | ours | hp.x | Δ |
|---|---|---|---|
| E_tot Ha | −71.60463 | −71.61246 (−143.22493 Ry) | 7.8 mHa (the b6 basis gap, as at Γ) |
| d occupation | 2.297 | 2.283 | +0.6 % |
| χ₀ Ha⁻¹ | −11.418 | −11.360 (−0.41746 eV⁻¹) | +0.5 % |
| χ Ha⁻¹ | −3.809 | −3.867 (−0.14211 eV⁻¹) | −1.5 % |
| **U eV** | **4.757** | **4.637** | **+2.6 %** |
| Γ → q 2³ shift eV | −0.368 | −0.389 scalar / −0.382 2×2 | 0.02 eV |

The q ≠ 0 machinery reproduces hp.x's finite-size shift to 0.02 eV; the residual is the Γ basis offset carried over (+0.10 → +0.12 eV).
One definition difference left: ours perturbs the Ti0 sublattice only (8 of hp.x's 16 sites; 7 meV at Γ).  **TiO₂ is closed as a
cross-code point: CP2K-matched on CP2K's basis, QE-matched on an l-converged one, at Γ, with O 2p, and on the q-mesh.**

## HubbardUPlan §2 as it stood on 2026-10-22, before the handoff rewrite (verbatim)

The Track-A corner-case bullet and S1/S2, archived when §2 was condensed to a START HERE for a new session.

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
  PW basis has no l limit).  ✅ **O 2p at Γ PASSED** (6×6 Hubbard-site matrix, b6): U(O 2p) 6.735 vs hp.x 6.758 eV (−0.3 %), U(Ti)
  3.051 vs 2.927 (+4 %), screened χ element by element (H, same section).  ✅ **k/q mesh PASSED**: hp.x k2/q2 = 4.6373 ≡ its k4/q2
  4.6368 (k converged at 2³); ours (b6, k2/q2, kT 0) **4.757 vs 4.637 (+2.6 %)**, χ₀ +0.5 %, χ −1.5 %, the Γ→q2 shift −0.368 vs −0.389 eV.
  **TiO₂ is CLOSED as a cross-code point.**  Left in S1: per-irrep O pz.  **First `resolve: eigen` run (TiO₂ Γ, b6, 2026-10-21)**: the 1×1 (Macke
  pointwise) column is well behaved — Ti channels 1.77–2.06 eV, O p 3.96–4.31 eV (the [1̄10] plane normal is the MIDDLE channel,
  n 1.629, 4.04 eV — not the most occupied, as it is in β-MnO₂) — and full inversion is the
  ill-posed 17–36 eV Macke warns of.  The ORBITAL-DIAGONAL variant (−22 to −846 eV) was ill-posed by DEFINITION (zeroing a site's
  inter-orbital off-diagonals leaves χ indefinite) and is REMOVED (`db09a3f4`, merged `1d33184c`, ctest 1050/1050).  ⚠ **QE
  orbital-resolved FD (pw.x ±α_m restarts, 6-atom cell, Γ; `~/Code/qe-runs/TiO2/orfd/`) is NOT yet a reliable oracle**: channel
  definitions agree (QE Ti eigenvectors xz/yz · 0.23 z²+0.77 xy · x²−y² · … vs ours 0.27/0.73), but on the pure x²−y² channel QE's
  first iteration shows NO response (identical to 8 digits for ±α) and χ ends with the wrong sign — QE's eigenstate tracking appears
  to start from a different basis at a restart, which would also shrink |χ₀| on k1/k2 (QE −0.418/−1.266 vs ours −0.455/−1.587 Ha⁻¹;
  U₁ₓ₁ QE 2.29/+0.27 vs ours 2.03/−0.20 eV).  The −0.20 eV channel is anti-screened in ours (|χ| 1.606 > |χ₀| 1.587); QE's k2
  is barely screened (1.251 vs 1.266) — both a difference of near-equal inverses.  Next: find how QE seeds the tracked eigenvectors
  on a restart (`PW/src` tracking, ref 127 of Macke) or run the ±α from scratch, before trusting any QE resolved number.

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
  **Set up 2026-10-21** (`024bd682`; `IntegrationTests/QE/README.md` §S2): `materials.json` `MnO2_beta_A1AFM` (Bolzan 1993 neutron cell,
  a 4.4041 / c 2.8765 Å, u 0.3046; corner Mn +m, body-centre −m), QE step 1 at Γ and k 4×4×6.  ⛔ **Plain LDA has NO gapped reference**:
  Γ is a 10 meV two-level artifact at E_F (the parked-SrVO₃ shape), k446 a semimetal (Mn ±3.0 μB, 4 of 1680 occupations fractional).
  So hp.x's fixed 2-step is out.  **RULED 2026-10-21 (user): route (b)** — wait for the OpenWork §2 row "Response: REPORT
  orbital-resolved channels" to land (being built in a parallel session), then open the gap with U(O pz) at U_in, as Macke does
  (PBEsol: U(O pz) 4.62 + U(Mn t2g) 1.64 eV → 1.01 eV gap), and run the matched oracle from that gapped reference.  Rejected:
  (a) a matched smeared oracle — a semimetal's χ₀ is set by the smearing width.  ✅ **Unblocked 2026-10-21** (resolved channels merged,
  `16528a14`; ctest 1050/1050).  **QE: U(O pz) 4.62 eV ALONE opens a 0.842 eV gap** at k 4×4×6 (+ U(Mn t2g) 1.64: 0.856 eV) —
  `IntegrationTests/QE/README.md` §S2 has the eigen-indices.  O pz is a PURE irrep (plane normal of O's C₂ᵥ), so our `Uirrep_eV` can
  set the same U_in; Mn t2g cannot be (one t2g state shares D₂h's a_g with an e_g one) and is not needed for the gap.  Next: the
  matched oracle from that reference — QE pw.x ±α_m FD (orbital-resolved; hp.x is shell-averaged only) vs ours `resolve: eigen`, on a
  k-mesh both can afford, after the Mn/O basis gets the TiO₂ treatment (O p tail + f, vetted in rutile).  Our side
  resolves the structure (`--check`: 38 e); its Mn basis is VA's MnO-conditioned block — the TiO₂ lesson (O p tail + f) applies first.

**Corrections and findings of 2026-10-22, after the text above:**
- ⚠ The "[1̄10] plane normal is the MIDDLE channel, n 1.629" correction above was WRONG — QE's eigenvector table was read by ROWS
  where QE prints COLUMNS.  Signed compositions (`69778563`) settle it: on every rutile O site slot Γ0 = the σ orbital along the C₂
  axis (O at (u,u,0): (x+y)/√2, n 1.629, 1×1 U 4.04 eV), Γ1 = z, Γ2 = the plane normal ((x−y)/√2, n 1.749, 1×1 U 3.96 eV) —
  matching QE read by columns (σ 1.626, normal 1.738).  The normal is the MOST occupied O p channel in TiO₂ AND β-MnO₂.
- ⚠ Consequently the β-MnO₂ m0 runs of 2026-10-22 put U(O pz) 4.62 eV on slot 0 = the σ BOND, not the normal; the deck
  (`~/Code/qchem6-runs/MnO2/basis_scan/m0.json`) now targets slot 2.
- β-MnO₂ seeding: the library had no Mn⁴⁺ (q7, 3 e) entry AND `IonicFormalCharges` drained the first of two equal donors (Mn+7 /
  Mn+1), so IonicSAD seeded one Mn as nothing and the other non-magnetic — the AFM order never formed (m0: non-magnetic, oscillating,
  200 iterations).  Tie-break + Mn⁴⁺ entry make the AFM form (Mn +3.04/−3.74, SURVIVED) but that run (U on σ) did not converge.
  Both changes are PARKED on branch `seed-tiebreak` (`2c4a7773`): they make `GPW_MnO.Γ_Shub_Pol_Smear_Anchor_Long` (MnO AFM-II SR Γ
  imposed, smeared) LOSE its order (moment dead at iteration 30; passes at HEAD) — HEAD's pass rests on the buggy seed.
  ⚠ CORRECTION (same day): that seed is NOT "+4/0" — `IonicSADTargets` keeps one target per SPECIES (last atom wins), so HEAD seeds
  MnO as NEUTRAL Mn d⁵4s² ×2; the fix gives the Mn²⁺ d⁵ the anchor's comment claims.  (The parked commit's message says "+4/0".)
- QE β-MnO₂ at Γ with U(O pz) 4.62 eV only: gap 0.325 eV, E = −186.85322133 Ry (`~/Code/qe-runs/MnO2/mno2b_G_OpzU.scf.out`).
- QE orbital-resolved FD on TiO₂ PARKED by the user (restart eigenstate tracking: x²−y² shows no first-iteration response).

## 2026-10-23 — the S2 merge sequence, executed

The spec as written (verbatim from `doc/HubbardUPlan.md` §2 S2):

  **What blocks the first converged run — and THE MERGE SEQUENCE (spec'd 2026-10-22; user: batch FIRST, incremental merges preferred):**
  ⛔ IonicSAD cannot seed β-MnO₂'s AFM on main (no Mn⁴⁺ entry + the formal-charge tie bug + per-species targets; OpenWork §3 rows
  "IonicFormalCharges drains…" and "IonicSAD targets are per SPECIES").  The fix is PARKED on `seed-tiebreak` (`2c4a7773`) because it
  gives MnO the Mn²⁺ d⁵ seed its SR Γ anchor claims, and that arm then loses AFM.  The seed fix MOVES ANCHORS (every IonicSAD cell with
  repeated cations), so it lands ON TOP of the anchor batch's baseline and MnO is re-judged ONCE on the combined code:
  1. **Merge `anchor-batch` into qchem6 main** (`git fetch ~/Code/qchem1 anchor-batch`; head `7591ec69` on 2026-10-22: refill-on-restart
     `32b9d359`, D-MINDRHO `b60c9722`, D-FERMI-DEFAULT `45b4487d`, Becke/screening `99b3d1b6` `0aece9ca` `81dbb79e`, temp names
     `7ebc4420`).  Only its DONE items: the Kerker/Pulay-default window is not needed (incremental merges are fine).  Expect doc
     conflicts (OpenWork, HubbardUPlan, records: keep both).  D-FERMI-DEFAULT RETIRED the deck key `globalFermi`: run
     `scripts/vet-decks` and migrate any deck still using it (e.g. `decks/campaigns/SrVO3/srvo3_k444_chi0.json`).  Build `allTests`;
     `scripts/memsafe ctest -j8 -LE long` AND `ctest -j3 -L long` green (the batch's STATUS says nothing needed a re-judge).
  2. **Rebase `seed-tiebreak` onto that main.**
  3. **`IntegrationTests/GPW/MnO.C` `Γ_Shub_Pol_Smear_Anchor_Long`: `o.seed = SeedStrategy::SAD`** (neutral polarized SAD: Hund pairs +
     flips, and the imposition decoration then uses neutral targets — what the arm really ran on), and correct its comment: the recipe's
     "IonicSAD seed (Mn²⁺ d⁵ + O²⁻)" was in fact neutral Mn d⁵4s²; the true Mn²⁺ seed loses AFM on this SR Γ arm (moment dead at
     iteration 30, −60.97044 Ha) — add that as a remainder on OpenWork §3 "MnO ordering at Γ/SR: FM lies BELOW AFM-II".
  4. **Run that test ALONE**: expect E −61.41455 ± 2e-3 and the integrated Mn moment > 4.0 — bit-identical to main if SAD really is
     main's effective seed.  If not: STOP and take it to the user (the fallback is a re-judge, not a recipe tune).
  5. **Full sweep, both tiers** (the count goes up by the new UT `IonicFormalCharges.EquivalentDonorsShareTheTransfer`), note any
     anchor the new seed moved in other IonicSAD tests (the VA MnO +U gate passed with it), then merge `seed-tiebreak` into main.

**How it ran:**
1. `anchor-batch` (7591ec69) merged into main as `cacea828`; OpenWork conflict resolved as in qchem1's 6f7c6578.  `scripts/vet-decks`
   49/49 OK (no deck still used `globalFermi`: the SrVO3 hit was a comment string).  Fast 1052/1052, long 4/4.
2. `seed-tiebreak` rebased (no conflicts).
3–4. ⚠ Step 4 FAILED its premise first: with `SeedStrategy::SAD` the anchor landed at −61.41533 (0.78 mHa lower, Mn +4.4554/−4.4418,
   O −0.0068) on a 12-op group.  Cause: the decoration is DERIVED from the seed's per-species targets, and neutral O (2p⁴) is labelled
   UP on both sites → a ferrimagnet.  Main's IonicSAD had O²⁻ (non-magnetic) → 24 ops.  The CORRECTION of 2026-10-22 was half right:
   main's IonicSAD seed was neutral Mn ×2 **+ O²⁻ ×2** (30 e for a 26 e cell), not neutral SAD.  User ruling: the ordering is
   STATED, never inferred (OpenWork §3 row).  The anchor got `siteSpins {+1,-1,0,0}` on main (`52735109`, bit-identical: −61.4145467920,
   17 iterations); on the branch, SAD + stated spins → −61.4145467920 in 15 iterations, Mn ±4.4506, O 0 — the same fixed point from a
   different first iterate (pin 10).  The TRUE ionic seed (Mn²⁺ ×2 + O²⁻, 24 ops) loses AFM: iterate-1 moment 0.30, −60.97043 Ha,
   unconverged at 82 iterations (OpenWork §3 row).
5. Branch sweep: fast 1053/1053 (+1 UT), long 4/4; `anchordiff` at 1e-10 Ha: no energy moved; iterations moved only in the SR anchor
   (17 → 15) and the VA +U gate's U=0 arm (18 → 38, the fixed IonicSAD seed's 0.30 start).  Merged as `beffb23b`.

## 2026-10-23 — β-MnO₂ m0 (VA, Γ, U(O pz) 4.62 eV on slot Γ2), after the V_short + Ewald fixes (main 768a4f70)

Recipe found on the way (decks `~/Code/qchem6-runs/MnO2/basis_scan/m0{,_imp}.json`, logs beside them): `siteSpins` STATED; stage 1 Kerker+Pulay (8/5, α 0.45) at kT 1e-3 to Δρ 3e-4; stage 2 standalone GDM at kT 0 gated on |ΔE/E| < 1e-9 with Δρ 1e-4 and the [F,D] gates off.  Rejected: Kerker+Pulay alone (its history goes unstable near Δρ 1e-4: `m0_pulay_diverged.log`), relax 1.0 (`m0_relax1_stopped.log`), Ladder at kT 0 from the seed (aufbau 2-cycle: `m0_ladder_kT0_2cycle.log`).  Before the fixes the FREE run broke the symmetry (Mn +2.87/−3.79, O non-equivalent: `m0_free_twostage_prefix.log`) — a symptom of the V_short defect.

| arm | stages | Etot (Ha) | Mn | O moments | O p slot Γ2 occ. up/dn | frontier at the end |
|---|---|---|---|---|---|---|
| imposed (16 ops) | 20 + 9 | −93.01214291 | ±3.3581 | ±0.1008 | 0.954 / 0.972 | hole on every other GDM row; aufbau (gap 1.9 eV) only at the LAST step, a configuration change |
| free | 28 + 15 | −93.01215953 | +3.3586/−3.3592 | ±0.10 (4e-3) | 0.954 / 0.972 | NON-AUFBAU: hole, gap −2.35 mHa |

QE (same U, Γ): Γ2 0.961 / 0.972, gap 0.325 eV, E −93.42661 Ha (ours 0.414 Ha high: the basis scan).  Stage 1 ends with −TS = −1.6 mHa at kT 1e-3: the smeared frontier is fractionally occupied, and the kT = 0 GDM stage keeps the occupied block it was handed — OpenWork §3 "GDM crawls at a near-degenerate frontier" (same mechanism as the Mn–O dimer).  No gap from a kT = 0 GDM stage is trustworthy until that is fixed.

## 2026-10-24 — GDM aufbau re-pick: built, measured worse, kept OFF

User asked to let GDM re-pick its occupied block by aufbau when the frontier crosses (OpenWork §3 "GDM crawls…").  Built in `tSCFIrrepAcceleratorGDM::ComputeStep`: after the pseudo-canonical split, if the block's lowest virtual is > `AufbauTol` below its highest occupied level, re-choose the lowest-energy `no` orbitals, restart CG; cap `AufbauMaxRepicks` = 8 per block, then a loud give-up.  Unit tests `ACrossedFrontierIsRepickedByAufbau`, `TheRepickCapLeavesTheBlockAlone`.

| gate | without | with re-pick |
|---|---|---|
| Mn–O dimer (VA, kT 0, mult. 6, GDM from seed) | −30.342242 at iter 44, crawling, Mn 5.35 | −30.344618 at 200, Mn **4.54** (other state), still crawling; 8 + 8 re-picks then give-up; ends non-aufbau |
| β-MnO₂ m0 imposed (stage 2) | 9 iterations, −93.01214291, aufbau only on the last step | 18 iterations, −93.0121424, ends non-aufbau (0.21464 below 0.21694) |

Late swaps are a near-degenerate pair that crosses back after every swap (no integer aufbau solution at kT 0); early swaps on m0 moved up to 5 orbitals across 0.4 Ha, from stage 1's rounded fill.  User: keep the code, DEFAULT OFF (`AufbauTol = 0`).

## 2026-10-24 — β-MnO₂ Γ gap: the Hamiltonians AGREE at U = 0; the difference is in the U term

| Γ, smeared kT 1e-3 Ha | QE | ours (VA, imposed) |
|---|---|---|
| U = 0: frontier pair (bands 19/20) | 8.3259 / 8.3362 eV (Δ 0.010 eV), occ 0.61/0.39 per spin, −TS −1.09 mHa | 0.2460 / 0.2464 Ha (Δ 0.011 eV), occ 0.6/0.4 — `m0_imp_U0.log` (`QCHEM_DIAGNOSTICS=frontier`) |
| U = 0: next band | +1.75 eV | +1.78 eV |
| U = 0: Etot | −93.46773 | −93.07477 (0.393 Ha: basis) |
| U(O pz) 4.62: which orbital | eigenvector of ns with the LARGEST occupation, (−p_x+p_y)/√2 = plane normal | slot Γ2 (the same orbital) |
| U(O pz): O p occupations up / dn | 0.791 0.831 **0.961** / 0.640 0.808 **0.972** | 0.774 0.816 **0.954** / 0.637 0.793 **0.972** |
| U(O pz): O moment | 0.163 | 0.100 |
| U(O pz): frontier | GAPPED 0.3255 eV, integer, −TS = 0 | pair 2.3 mHa apart, fractional (−TS −1.6 mHa in stage 1) |
| U(O pz): E_U | 0.02216 Ha (= Dudarev Σ U/2 λ(1−λ) on the U'd eigenvalue, checked by hand: 0.02197) | NOT PRINTED in the energy breakdown |

So: same LDA Hamiltonian, same U orbital, near-identical occupations, yet QE's U splits the pair and ours does not.  Suspects left: the U POTENTIAL's projector (QE: ortho-atomic over all UPF atomic wavefunctions; ours: our atomic radial, Löwdin over the deck's set) and its energy form; next: print our E_U and compare the frontier pair's O pz projection in both codes (QE `projwfc.x` on ortho-atomic).

**Band characters (QE `projwfc.x`, ortho-atomic; run on copies of the saves with `<dftU>` stripped — projwfc 7.6 cannot read the orbital-resolved block) and our levels (`QCHEM_DIAGNOSTICS=frontier=16`), eV from the Mn d band (QE band 20):**

| QE band (character) | 13 | 14 | 15 (O z) | 16 (O normal) | 17 (in-plane) | 18 (in-plane) | 19 HOMO (O z) | 20 Mn d | 21 | 22 |
|---|---|---|---|---|---|---|---|---|---|---|
| QE U(O pz) | −3.10 | −2.94 | −2.39 | −1.06 | −0.54 | −0.45 | −0.33 | 0 | +1.89 | +4.11 |
| ours U(O pz) | −2.88 | −2.56 | −1.89 | −0.67 | −0.27 | −0.13 | −0.13 | 0 | +1.95 | +3.98 |

At U = 0 the same bands agree to ≤ 0.04 eV (QE: band 19 = pure O normal π, band 20 = Mn d, 0.010 eV apart).  With U, QE's normal band drops ~1.05 eV and its O z band rises 0.76 eV relative to Mn d; ours does the same thing ~0.3 eV less: EVERY occupied O p band sits 0.2–0.5 eV higher relative to Mn d.  E_U 0.02424 vs 0.02216 Ha, same Dudarev eigenbasis form, same orbital.  Suspect: the projector radial (QE: the UPF PSWFC; ours `atomicRadial`) and its Löwdin tails; QE's O moment 0.163 vs ours 0.100.  Also: the earlier "fractional frontier" was mostly FD smearing (kT 1e-3 Ha) of our 0.128 eV gap.

**The projector tail (same day).**  Our O 2p projector = the GTH pseudo-atom's 2p (the SAME function the UPF's PP_CHI holds — `gth2upf` wrote it) PROJECTED onto the O site's p shells: VA's stop at 0.465, so it captures 0.976 of the norm; ⟨r⟩ 1.186 vs 1.265 bohr, norm beyond 2.5 bohr 2.1 % vs 6.6 %, normalised overlap 0.988.  Adding ONE diffuse O p shell (0.18; decks `m0_imp_{s1,U0}_Op018.json`, `vaOp018.bsd`; vet: nothing trimmed, min eig S 2.7e-4):

| Γ, U(O pz) 4.62 | VA | VA + O p 0.18 | QE |
|---|---|---|---|
| projector captured | 0.976 | 0.9994 | — |
| gap | 0.128 eV (fractional at kT 1e-3) | 0.237 eV (integer) | 0.325 eV |
| O normal band − Mn d | −0.67 | −0.88 | −1.06 |
| E_U (Ha) | 0.02424 | 0.01767 | 0.02216 |
| Etot U / U=0 (Ha) | −93.01214 / −93.07477 | −93.34058 / −93.38443 | −93.42661 / −93.46773 |

The remainder (0.09 eV of gap, 0.1–0.2 eV on the O p bands) is what a SITE-ONLY contracted projector cannot hold at any exponent: the tail's overlap with the NEIGHBOURS' functions.  → OpenWork §2 row "+U projector independent of the orbital basis".  And O p 0.18 is the basis scan's first big win: 0.33 of the 0.41 Ha to QE, both arms.


## 2026-10-24 — The pseudo-atomic +U projector (HubbardProjectorPlan increments 1–3): it matches QE's, and the gap does not move

Built (commits on main, this day): an atom-centred function set and the basis's overlap with it (`qchem.BasisSet.AtomCentredSet`; GPW through the KB projection seam, the set's own Bloch overlap through a new engine pair kernel `LatticeSum1E::MakeOverlap(phase,A,g,h)`), `LowdinProjector::Atomic(B,O,…)` (ortho-atomic in the TRUE metric O), the pseudo-atom builder `qchem.PseudoAtomOrbital` (moved out of the facade), and the deck field `"projector": "pseudoAtomic"`.  The O 2p projector is now the pseudo-atom orbital in its own 16-exponent pool: against the O UPF's `PP_CHI.2` (`~/Code/qe-runs/SrVO3/O.pz-gth-q6.UPF`) the normalised overlap is **1.00000000** and the norm beyond 2.5 bohr is 6.59 % in both.

The increment-3 gate (decks `~/Code/qchem6-runs/MnO2/basis_scan/m0_imp_s1{,_Op018}_PA.json`, stage 1, kT 1e-3, `QCHEM_DIAGNOSTICS=frontier=16`):

| Γ, U(O pz) 4.62 eV | VA site-shell | VA **pseudo-atomic** | VA+O p 0.18 site-shell | VA+O p 0.18 **pseudo-atomic** | QE |
|---|---|---|---|---|---|
| O 2p captured by the block (B†S⁻¹B/O) | 0.976 | 0.994 | 0.9994 | 1.000 | — |
| gap (HOMO→LUMO, eV) | 0.128 | **0.128** | 0.237 | **0.234** | 0.325 |
| E_U (Ha) | 0.02545 | 0.03014 | 0.01767 | **0.02214** | 0.02216 |
| O p slot 2 (normal) up / dn | 0.949 / 0.972 | 0.940 / 0.966 | 0.968 / 0.978 | 0.960 / 0.972 | 0.961 / 0.972 |
| Etot (Ha) | −93.01334 | −93.00879 | −93.34058 | −93.33536 | −93.42661 |

**Verdict: the gate as written FAILS, and its premise is refuted.**  The projector is now QE's own function and is basis-independent (captured ≈ 1 in both bases), yet on plain VA the gap is unchanged to 0.001 eV and on VA + O p 0.18 to 0.003 eV.  So the 0.128 → 0.237 eV move that one diffuse O p shell made was the ORBITAL basis (the KS states, Etot −0.33 Ha), not the projector's tail; the 2026-10-24 "remainder is the tail's overlap with neighbours' functions" reading above is wrong.  What the projector DOES set is E_U: with the better basis the pseudo-atomic E_U is QE's to 0.1 % (site-shell route: 20 % low), and the O p occupations agree with QE's to ≤ 0.001 on the U'd slot.  The remaining gap and energy difference to QE (0.09 eV, 0.09 Ha) is the orbital basis — the S2 basis scan's job.  The basis-independence check of the occupations (VA vs VA + O p 0.18 to ~1e-3) cannot be read off these runs: the DENSITY changes with the basis (0.02 on slot 2 in both routes alike).

## 2026-10-24 — β-MnO₂ diffuse/polarization basis scan (pseudo-atomic projector; decks + logs `~/Code/qchem6-runs/MnO2/basis_scan/scan/`)

Base VA + O p 0.18 (gap 0.234 eV, Etot −93.33536); each row adds ONE shell; gap from the frontier window (stage 1, kT 1e-3 Ha, imposed):

| added | gap (eV) | ΔE (mHa) | vet |
|---|---|---|---|
| Mn p 0.5 | 0.231 | −11.4 | nothing trimmed |
| Mn d 0.18 | 0.229 | −0.45 | trimmed Mn s 0.1 instead |
| O p 0.08 | 0.234 | 0 | trimmed O p 0.08 itself |
| O s 0.15 | 0.242 | −20.9 | trimmed Mn s 0.1 |
| O d 0.8 | 0.226 | −14.2 | nothing |
| O d 1.25 | 0.220 | −8.9 | nothing |
| O p 0.08, vet OFF, orthoTol 1e-7 (diagnostic) | **0.234** | −25.4 | 176 kept, min pivot 6.4e-6 |
| (O d 1.2: INVALID — merged with O p 1.2 by the same-exponent bug, 172 functions, +753 mHa) | | | |

The gap is FLAT (0.220–0.242) under every addition, including the diffuse O p kept past the conditioning wall.  So beyond O p 0.18 the residual 0.09 eV is NOT the Gaussian orbital basis.  Band by band (eV from Mn d; QE bands 13–22): ours −3.04 −2.95 −2.30 **−0.84** −0.41 −0.32 −0.23 | +1.92 +4.14 vs QE −3.10 −2.94 −2.39 **−1.06** −0.54 −0.45 −0.33 | +1.89 +4.11: within ~0.1 eV except the U'd O-normal band (+0.22: pushed down 80 % as far as QE's) and the in-plane bands that carry normal character (+0.13).  With projector, occupations and E_U matching QE, the residual is HOW STRONGLY U ACTS PER BAND — next: our bands' weights on the U'd projector vs QE's projwfc (band 16: 0.987).

**Per-band +U slot weights vs QE `projwfc.x` (same day; `QCHEM_DIAGNOSTICS=u_bands=10`, VA + O p 0.18, pseudo-atomic projector; log `m0_imp_s1_Op018_PA_ubands.log`).**  Slots: O p Γ2 = plane normal (U'd), Γ1 = z, Γ0 = in-plane σ.

| QE band | character | QE | ours |
|---|---|---|---|
| 15 | O z | 0.991 | Γ1 0.993 |
| 16 | O normal (U'd) | 0.987 | Γ2 0.988 |
| 17 | in-plane / normal | 0.706 / 0.290 | 0.689 / 0.308 |
| 18 | in-plane / normal | 0.738 / 0.257 | 0.724 / 0.273 |
| 19 HOMO | O z | 0.994 | Γ1 0.996 |
| 20 LUMO | Mn d / in-plane / O s | 0.649 / 0.245 / 0.082 | 0.652 / 0.242 / 0.083 |

Characters agree to ~0.01 on every band, so the first-order U shift per band is QE's (the "U acts more weakly per band" hypothesis is WRONG).  Energies from the deep pure-O-z band 15 (untouched by U): O normal +1.46 vs QE +1.33, HOMO O z +2.07 vs +2.06, LUMO Mn d +2.30 vs +2.39.  The gap is short because the Mn d LUMO sits 0.09 eV LOW (and the U'd band 0.13 eV high): same characters, different energies — the self-consistent response to U, not its action.  Open.


## 2026-10-24 — The gap residual is the SELF-CONSISTENT RESPONSE of two bands to U; convergence and pin 28 cleared

Decks `~/Code/qchem6-runs/MnO2/basis_scan/m0_imp_{,U0_}Op018_PA.json` (stage 1 + banked GDM stage 2, kT 0) and `m0_imp_s1t_{,U0_}Op018_PA.json` (stage 1 only, kT 1e-3 like QE, gated at Δρ 3e-5 just above the quadrature floor); VA + O p 0.18, pseudo-atomic projector, `QCHEM_DIAGNOSTICS=u_bands=10`.

1. **Convergence: cleared.**  The GDM-converged U state (kT 0, |dE/E| < 1e-9) moves the bands ≤ 0.03 eV from the loose stage-1 run (O normal +1.456 vs +1.459, HOMO +2.069 vs +2.067, LUMO +2.274 vs +2.301, eV from the deep O z band).
2. **Pin 28 (XC from ρ[D]): cleared.**  GDM mixes D itself, so its fixed point has XC on ρ[D_mix] with no fitted ρ̃ and no cusp-deficit correction; the residual is still there.
3. **Trap: a U = 0 arm cannot take the kT-0 GDM stage.**  U = 0's frontier is a near-degenerate pair (QE occupations 0.61 / 0.39); integer occupations put it 20 mHa ABOVE the smeared state (−93.36180 vs −93.38443).  The U / U=0 difference is therefore read on the tight smeared pair.
4. **Per band, matched by character, eV from the deep pure O z band (U-free):**

| band | ours U=0 | QE U=0 | ours U | QE U | ΔU ours | ΔU QE |
|---|---|---|---|---|---|---|
| O normal (U'd slot) | +3.141 | +3.107 | +1.459 | +1.332 | **−1.682** | **−1.775** |
| O z (HOMO under U) | +2.031 | +2.033 | +2.068 | +2.068 | +0.037 | +0.035 |
| in-plane pair | +2.620 / 2.668 | +2.596 / 2.635 | +1.894 / 1.979 | +1.856 / 1.944 | −0.726 / −0.689 | −0.740 / −0.691 |
| Mn d (LUMO) | +3.151 | +3.117 | +2.299 | +2.393 | **−0.852** | **−0.724** |

At U = 0 every band is QE's to ≤ 0.04 eV on this basis (a near-uniform +0.03).  Under U two bands respond differently: the U'd O-normal band falls 0.09 eV LESS than QE's and the Mn d LUMO (no U on it) falls 0.13 eV MORE — exactly the gap deficit (0.034 − 0.128 = −0.094 eV).  The U'd slot's occupations are QE's (0.960 / 0.972 vs 0.961 / 0.972), the characters agree to 0.01 (previous entry), the ortho-atomic set is QE's (Mn 4s+3d, O 2s+2p per site, from `mno2b_G_OpzU.projwfc.out`) and E_U agrees to 0.1 %, so the FIRST-order shift U(½−λ) is QE's.  What differs is the second-order, self-consistent response of the rest of the Hamiltonian.  Hint (partition-dependent, not conclusive): the O moment falls 0.0276 under U in ours (Becke cells, 0.1467 → 0.1191) vs 0.0144 in QE (spheres R = 0.176 alat, 0.1499 → 0.1355); Mn rises +0.0034 vs +0.0084.

Side defect found: after a GDM stage `u_bands` prints ~0 slot weight on every VIRTUAL band (energies right) → OpenWork §3 row.

**The occupation response to U, on the shared ortho-atomic projectors (same day).**  QE: `mno2b_U0_tok` and `mno2b_G_OpzU_Mn` (`~/Code/qe-runs/MnO2/`; token U = 1e-8 on Mn 3d and O 2p so ns is printed — orbital-resolved form, QE refuses to mix it with shell-averaged U; energies unchanged to 1e-8 Ry).  Ours: the tight smeared pair `m0_imp_s1t_{U0_,}Op018_PA` (`[+U]` lines).

| | U=0 QE / ours | U QE / ours | Δ QE / ours |
|---|---|---|---|
| Mn d total | 5.0411 / 5.0708 | 4.9187 / 4.9423 | −0.122 / −0.129 |
| Mn d moment | 3.514 / 3.518 | 3.551 / 3.541 | +0.037 / +0.023 |
| O p total | 4.924 / 4.911 | 5.002 / 4.994 | +0.078 / +0.083 |
| O p moment | 0.180 / 0.183 | 0.163 / 0.161 | −0.017 / −0.022 |

The CHARGE response (Mn d → O p) agrees to 5–7 %; the Mn SPIN response is ~60 % of QE's (the O-moment hint, now on shared projectors).  Order of magnitude: neither difference alone moves the Mn d LUMO 0.13 eV (spin ~0.01 eV via an Mn exchange splitting ~0.9 eV/μB; charge only through an unscreened on-site Coulomb) — so not an occupation effect; the clearly different ingredient is the Mn SPIN response → the Mn XC (exchange) potential's response is the lead suspect.  Open.

**The three XC checks (user's list, same day; decks + logs `~/Code/qchem6-runs/MnO2/basis_scan/xc/`, U / U = 0 pairs on the tight smeared decks, bands named by `u_bands` character, eV from the deep O z band):**

| U arm | gap | O normal | O z top | Mn d | Mn d moment |
|---|---|---|---|---|---|
| baseline (Becke nR 40, L 29, Delta fit, ρ_mix + cusp deficit) | 0.231 | +1.459 | +2.068 | +2.299 | +3.5407 |
| 1. mesh nR 80, L 35 (425 k pts, ~3×) | 0.232 | +1.459 | +2.067 | +2.299 | +3.5406 |
| 2. V_xc fit PlaneWave | 0.235 | +1.457 | +2.067 | +2.302 | +3.5409 |
| 3. `policy.xcFromDM` ("XC rho source: rho[D] WHOLESALE") | 0.231 | +1.459 | +2.068 | +2.299 | +3.5407 |

The U = 0 arms likewise (Mn d +3.151, moment +3.5176–3.5178).  ALL THREE RULED OUT: the mesh, the fit and the ρ route move nothing at the 1e-3 eV / 1e-4 e level (check 3 is identical to 1e-8 — at the fixed point the cusp-deficit route IS ρ[D]).  The Mn spin response (+0.023 vs QE +0.037) and the Mn d band's extra 0.13 eV drop remain unexplained; the decisive instrument is the plane-wave route (OpenWork §2 PW-1…4): same PP as QE, so any difference would be a Hamiltonian difference by construction.


## 2026-10-25 — QE cutoff ladder: both residuals survive at a cutoff a DENSE plane-wave solve can reach

Question: can the plane-wave instrument (OpenWork §2 PW-1…4) run DENSE, at a cutoff matched with QE, instead of waiting on a sparse/matrix-free route (PW-2)?  Only if QE's own response to U is already converged at that low cutoff.  Decks + outputs `~/Code/qe-runs/MnO2/ladder/` (`mno2b_lad<Ry>_{U0,U}`; ecutrho = 4 × ecutwfc; otherwise the 300 Ry `mno2b_U0_tok` / `mno2b_G_OpzU_Mn` decks verbatim, token Mn 3d U so ns is printed); extractor reproduces the banked 300 Ry numbers exactly.  Moments from the ortho-atomic ns (up − down), bands from `projwfc.x` (spin 1, eV from the deep pure O z band; O normal = the odd one of the three O pxy bands, Mn d = the Mn d + O s frontier band).

| ecutwfc (Ry) | N_PW (Γ) | Δ Mn d moment | Δ Mn d total | Δ O p total | Δ O normal (eV) | Δ Mn d band (eV) | gap U (eV) |
|---|---|---|---|---|---|---|---|
| 60 | ~2.9 k | +0.058 | −0.120 | +0.075 | −1.767 | −0.629 | 0.601 |
| 80 | ~4.5 k | +0.039 | −0.131 | +0.083 | −1.736 | −0.773 | 0.271 |
| 100 | ~6.4 k | +0.037 | −0.127 | +0.081 | −1.751 | −0.752 | 0.289 |
| 150 | ~11.7 k | +0.037 | −0.124 | +0.079 | −1.768 | −0.732 | 0.318 |
| 300 (banked) | 33 k | +0.037 | −0.122 | +0.078 | −1.774 | −0.723 | 0.325 |
| **ours** (VA + O p 0.18) | — | **+0.023** | −0.129 | +0.083 | **−1.682** | **−0.852** | 0.231 |

From 80 Ry up both residuals are QE's at the 300 Ry level within a fraction of their size: the Mn spin response is +0.037–0.039 (ours +0.023), the Mn d band falls 0.72–0.77 eV (ours 0.85), the U'd O-normal band 1.74–1.77 (ours 1.68).  60 Ry is NOT converged (even the U = 0 bands are 0.3 eV off).  ⇒ **a dense PW run at 40–50 Ha (4.5–6.4 k PWs, real TRIM block) matched to QE at 80–100 Ry discriminates basis vs Hamiltonian; PW-2 is off this diagnosis's critical path.**

**Trap: QE's orbital-resolved U has several self-consistent fixed points.**  From the atomic start the 100 and 150 Ry U arms converged to OTHER states — spin-down O ns off-diagonal sign flipped (+0.026 vs −0.082), the U'd orbital's occupation 0.83–0.88 instead of 0.97, E_U doubled (0.091 vs 0.044 Ry), at 150 Ry even the O moment reversed — 23 / 9 mHa ABOVE the state reached by seeding from the 300 Ry U run (`startingpot='file'` with its charge-density + `occup.txt` + `Hubbard_m_order.txt`).  The 80 Ry atomic start did not converge in 100 iterations.  All U rows above are SEEDED (9–13 iterations); the atomic-start outputs are kept as `*.atomicseed.scf.out` / `*.capped100.scf.out`.  A cross-code U comparison must name the basin (E_U and the U'd slot's ns), not just Etot.

**First plane-wave run through the facade (same day; PW-1 first increment, 7b89611c; deck + logs `~/Code/qchem6-runs/MnO2/pw/pw_smoke_E20_U0.*`).**  β-MnO₂ A1-AFM, U = 0, polarized, free, kT 1e-3, Kerker+Pulay, `basis.family planeWave`, 20 Ha: 1591 PWs (QE at 40 Ry: 1591), 13 iterations, **A = −87.8140996 Ha vs QE `mno2b_lad40_U0` −175.62822603 Ry = −87.8141130 Ha (Δ 13 μHa)**.  So at a matched basis the shared Hamiltonian stack (GTH local + KB incl. Mn d channel, Hartree, polarized LDA on the FFT raster, Ewald) is QE's.  Two breaks on the way, both OpenWork rows: the cusp-deficit guard throws on an exact-fit basis (run needed `xcCuspDeficit false`), and the dense cost (34 s of eigensolve per iteration, 1.84 GB peak at 1591 PWs ⇒ ~15 GB at the 80 Ry match) stops the 80 Ry arm until the memory is attributed.  Still needed for the U residuals: PW-3 (+U overlap face on `PlaneWave_IBS`).

**The matched-cutoff U = 0 arm (same day; footprint commit 8e2e02de; deck + log `~/Code/qchem6-runs/MnO2/pw/pw_E40_U0.*`).**  β-MnO₂ A1-AFM, U = 0, polarized, free, kT 1e-3, Kerker+Pulay, plane waves at 40 Ha = QE ecutwfc 80 Ry: **4571 PWs, identical to QE's**, `numBands` 30 (QE nbnd), 31 iterations, 8.09 GB peak, 78 min serial.  **Internal energy E = −92.6767990 Ha vs QE `mno2b_lad80_U0` E = F + TS = −92.6767752 (Δ −24 μHa)**; the free energies differ by 1.6 mHa, all of it −TS (ours −2.72 mHa Fermi-Dirac, QE −1.10 mHa `smearing='gaussian'` at the same width — the near-degenerate U = 0 frontier is fractionally filled here, unlike at 40 Ry).  Spin-up bands 1–20 from band 1: ≤ 7.5 meV from QE's, largest at the 19/20 frontier pair.  ⇒ at the cutoff where QE's U response is converged (the ladder above), the shared Hamiltonian on plane waves is QE's.  The U arm needs PW-3 (the +U overlap face on `PlaneWave_IBS`); a like-for-like −TS needs Gaussian smearing on our side, or QE at Fermi-Dirac.

## 2026-10-25 — THE PLANE-WAVE VERDICT: the U residuals are the GAUSSIAN ORBITAL BASIS

PW-3 landed (6c60c623 + 74dd6732: the +U pseudo-atomic overlap face on `PlaneWave_IBS`).  Decks + logs `~/Code/qchem6-runs/MnO2/pw/pw_E40_{U,U0h}.*`: β-MnO₂ A1-AFM at Γ, plane waves at 40 Ha = QE ecutwfc 80 Ry (4571 PWs, identical to QE's), the 12 ortho-atomic pseudo-atomic manifolds of `m0_imp_s1t_Op018_PA` (U = 4.62 eV on O p slot Γ2 = the plane normal, U = 0 elsewhere; the U=0 twin carries the same set at U = 0), free, kT 1e-3, Kerker+Pulay, `numBands` 30.  U: 16 iterations; U=0: 31 (A identical to the no-+U run, −92.6795194); both ~9.05 GB peak.  Projectors captured by the 80 Ry block (B†S⁻¹B/O): Mn d 0.9997, Mn s 1.0000, O s 0.9995, O p 0.9991.  QE = the cutoff ladder's seeded 80 Ry pair.  Bands from the deep pure O z band (`u_bands` characters).

| U − U0 at 80 Ry | QE | **ours, plane waves** | ours, Gaussian (VA + O p 0.18, 2026-10-24) |
|---|---|---|---|
| Mn d moment | +0.0391 | **+0.0385** | +0.023 |
| Mn d total | −0.1314 | **−0.1307** | −0.129 |
| O p total | +0.0829 | **+0.0825** | +0.083 |
| O p moment | −0.0187 | **−0.0189** | −0.022 |
| O-normal band (eV) | −1.736 | **−1.725** | −1.682 |
| Mn d band (eV) | −0.773 | **−0.769** | −0.852 |

Levels too: U=0 Mn d moment 3.4902 / total 5.0568 vs QE 3.4905 / 5.0569; U arm O normal +1.472, Mn d +2.446 eV vs QE +1.465 / +2.436, gap 0.281 vs 0.271 eV; the U'd O p slot occupations (0.788, 0.824, 0.957) up / (0.804, 0.639, 0.969) dn vs QE's ns eigenvalues (0.788, 0.826, 0.958) / (0.805, 0.641, 0.969).  Energies: U = 0 internal −24 μHa from QE; under U +0.59 mHa, of which E_U 24.19 vs 23.69 mHa (+2%) is the bulk — and a 1e-3 occupation shift on the U'd slot IS 2% of (n − n²): the projector convention (our TRUE set overlap O vs QE's represented B†B, captured 0.9991) is the natural source; not chased.

**Verdict.**  With QE's pseudopotential and a plane-wave basis, the same Hamiltonian stack (local + KB, Hartree, polarized LDA, Ewald, ortho-atomic +U) responds to U as QE does, to 1–2 %.  The two residuals of the Gaussian runs — the Mn spin response at 60 % of QE's and the Mn d band falling 0.08–0.13 eV too far — are therefore the GAUSSIAN ORBITAL BASIS.  This REVISES the 2026-10-24 basis-scan reading ("beyond O p 0.18 the residual is not the Gaussian orbital basis"): that scan added one diffuse/polarization shell at a time and watched the GAP, which stayed flat; the deficit is in the Mn d radial freedom to respond to SPIN and to U, which a single added shell did not supply.  Next: a Gaussian basis aimed at the Mn d response (uncontracted / tight d, a second d radial), judged on the U − U0 table above against the PW column — the PW run is now the reference at this PP.

## 2026-10-25 — the basis bisection: neither O p (either side) nor diffuse Mn d moves the residual; sharp Mn d is BLOCKED by a local-PP defect

Mn vs O × sharp vs diffuse, each variant = VA + O p 0.18 plus shells on ONE element on ONE side of its ladder, as a U / U0 pair on the banked imposed recipe (`m0_imp_s1t_{,U0_}Op018_PA`), judged on U − U0 against the PW 80 Ry column.  Decks, bases, logs: `~/Code/qchem6-runs/MnO2/basis_scan/bisect/`; extractor = the `[+U]`-line / `u_bands` parse that reproduces the baseline row exactly.

| U − U0 | Mn d moment | Mn d total | O p total | O p moment | O normal (eV) | Mn d band (eV) |
|---|---|---|---|---|---|---|
| PW 80 Ry (reference) | +0.0385 | −0.1307 | +0.0825 | −0.0189 | −1.725 | −0.769 |
| baseline VA + O p 0.18 | +0.0229 | −0.1285 | +0.0833 | −0.0217 | −1.682 | −0.852 |
| + Mn d 0.56, 1.19, 2.55, 0.18 (no trim, pivoted 1e-7) | +0.0245 | −0.1307 | +0.0846 | −0.0216 | −1.679 | −0.852 |
| + O p 1.93, 4.98, 20.6 | +0.0230 | −0.1284 | +0.0832 | −0.0218 | −1.686 | −0.847 |
| + O p 0.29, 0.75, 0.11 (vet trimmed O p 0.11, 0.18) | +0.0235 | −0.1283 | +0.0831 | −0.0219 | −1.679 | −0.852 |

**O is ruled out in both directions; the diffuse / mid-low Mn d side is ruled out.**  Two runs were not interpretable and taught two defects (OpenWork §3, 2026-10-25): (1) the vet trim at orthoTol 1e-4 deleted the VALENCE core of the densified Mn d ladder (sharp arm: d 0.38, 0.82, 1.74, 3.72; diffuse arm: s 0.1, d 0.18, 0.38, 0.56 — the diffuse U energy went UP 94 mHa after ADDING shells), hence the no-trim rerun above; (2) **a Mn d primitive at α = 76.7 makes the SCF COLLAPSE** — −6917 Ha "converged" after the trim, −15 000 → −21 000 Ha with the full valence d kept — with **Eloc = −8764 Ha** (EenNL +4.8): the GPW LOCAL-PP matrix elements are wrong for a tight pair, which no variational basis can do with correct integrals.  So the sharp Mn d side is UNTESTED, and the same defect is now a SUSPECT for the residual itself: it acts on exactly the tight Mn d pairs, and a milder error at the ladder's α = 36 (pair 72) could bias the Mn d level and its response while leaving the U = 0 bands within 0.04 eV.  Next: isolate the local-PP defect in a unit test (one Mn site, one d primitive, α = 36 vs 77, GPW Ven_PP vs the analytic molecular GTH local PP in a box) — it unblocks the sharp arm and tests the suspicion on the baseline in one step.

**The grid is NOT the residual (same day; decks + logs `~/Code/qchem6-runs/MnO2/basis_scan/ecut/base_E300_*`).**  The tight-pair aliasing (OpenWork §3; `Records/GridSharpTerms.md`) put the baseline's Mn d 36 pairs at E_cut/p ≈ 2 on a `BallOnly` raster (audit: density alias 1.8e-2, V_loc long 7.1e-3).  Rerun of the baseline pair at `densityEcut` 300 Ha (E_cut/p ≈ 4): U − U0 Mn d moment **+0.0229** (72 Ha: +0.0229), Mn d total −0.1281 (−0.1285), O p total +0.0831, O normal −1.682 (−1.682), Mn d band −0.851 (−0.852); energies moved ~1 mHa (U0 −93.38443 → −93.38337).  The defect is real (it collapses a basis with a tight d) but the α = 36 primitive carries too little weight here: the residual stays the GAUSSIAN BASIS, with O p and diffuse Mn d ruled out.  The sharp Mn d arm is now testable on `AliasFree` (≤ 3e-4 Ha on a d at 77 in the unit oracle).

**Sharp Mn d, Mn p, Mn f (same day; `~/Code/qchem6-runs/MnO2/basis_scan/pol/`; all four arms share basis.vet false, CholeskyPivoted 1e-7, raster AliasFree; the sharp arm no longer collapses on AliasFree).**

| U − U0 | Mn d moment | Mn d total | O p total | O p moment | O normal (eV) | Mn d band (eV) | ΔE(U) vs control |
|---|---|---|---|---|---|---|---|
| PW 80 Ry (reference) | +0.0385 | −0.1307 | +0.0825 | −0.0189 | −1.725 | −0.769 | — |
| control (baseline basis on these settings) | +0.0229 | −0.1281 | +0.0830 | −0.0216 | −1.683 | −0.851 | 0 |
| + Mn d 5.43, 11.6, 24.7, 76.7 | +0.0228 | −0.1284 | +0.0832 | −0.0216 | −1.683 | −0.851 | −9 μHa |
| + Mn p 0.4, 1.0, 2.5 | +0.0241 | −0.1269 | +0.0816 | −0.0208 | −1.704 | −0.815 | −35 mHa |
| + Mn f 1.0, 2.5 | +0.0235 | −0.1281 | +0.0829 | −0.0215 | −1.691 | −0.847 | −8 mHa |

The two residuals SEPARATE.  The BAND responses move with Mn p (the basis has none: −35 mHa, and half of both band residuals closes) and a little with Mn f; sharp Mn d does nothing (variationally 9 μHa).  The Mn SPIN response does not move with ANY variant (+0.0228…+0.0241 vs +0.0385) — after O p (both sides), diffuse/dense Mn d, sharp Mn d, Mn p, Mn f — and the read-out is not it (Mn d projector captured 0.9998).  What still differs between the Gaussian runs and the PW ones besides the basis: the Gaussian runs IMPOSE the A1-AFM Shubnikov group, the PW ones ran free; and the XC route (Becke vs FFT raster — three XC variants moved nothing, 2026-10-24).  Next: the control pair FREE (imposeSymmetry false), and a larger Mn p (+ f) set to see where the band responses saturate.

## 2026-10-26 — Mn polarization closes the BAND residuals; the imposition is ruled out; the Mn SPIN response stays short

Same settings as the polarization arms (vet off, CholeskyPivoted 1e-7, AliasFree); `~/Code/qchem6-runs/MnO2/basis_scan/pol/`.  Every arm independently from IonicSAD, every arm in the SAME basin (U: Mn d ±3.540–3.549, U'd O slot 0.960–0.961 / 0.972 = QE's seeded ground basin 0.958 / 0.969, not its metastable one; U0: Mn ±3.516–3.525, slot 0.861–0.865 / 0.835–0.840) — including the FREE control, which finds the imposed A1-AFM state unaided (Becke Mn ±3.4131 vs ±3.4130; the old "free lands at Mn +2.87/−3.79" predates the V_short/Ewald fixes).

| U − U0 | Mn d moment | Mn d total | O p total | O normal (eV) | Mn d band (eV) | A(U0) vs control |
|---|---|---|---|---|---|---|
| PW 80 Ry (reference) | +0.0385 | −0.1307 | +0.0825 | −1.725 | −0.769 | — |
| control (imposed) | +0.0229 | −0.1281 | +0.0830 | −1.683 | −0.851 | 0 |
| control (FREE) | +0.0228 | −0.1284 | +0.0832 | −1.683 | −0.851 | −27 μHa |
| + Mn p 0.4, 1.0, 2.5 | +0.0241 | −0.1269 | +0.0816 | −1.704 | −0.815 | −32.8 mHa |
| + Mn p 0.25 … 4.0 (5) | +0.0240 | −0.1258 | +0.0810 | −1.708 | −0.806 | −41.2 mHa |
| + Mn f 1.0, 2.5 | +0.0235 | −0.1281 | +0.0829 | −1.691 | −0.847 | −7.4 mHa |
| + Mn p (5) + f 0.7, 1.4, 2.8 | +0.0262 | −0.1228 | +0.0786 | −1.729 | −0.773 | −55.1 mHa |

1. **The imposition is ruled out**: free = imposed to 1e-4 in every column.
2. **The BAND residuals are Mn POLARIZATION**, and they CLOSE: p alone saturates near −0.806 / −1.708; p + f lands on PW to 4 meV (−0.773 vs −0.769, −1.729 vs −1.725).  f helps little alone (−0.847) but a lot on top of a complete p (−0.806 → −0.773).  The basis had NO Mn p or f; the U0 energy falls 55 mHa.
3. **The Mn SPIN response stays short**: +0.0262 at best (68 % of PW's +0.0385), and the CHARGE response (Mn d total, O p total) moves slightly AWAY from PW as the bands converge (−0.1228 vs −0.1307; +0.0786 vs +0.0825).  So the spin/charge response is a different deficit from the band one — still the basis by the PW verdict, but not O p, not Mn d (sharp or diffuse), not the imposition, not the grid, and only partly Mn p/f.  Next candidates: Mn d flexibility ON TOP of p5f3 (the d radial's freedom may only matter once the polarization is complete), Mn g, and the O side with Mn polarization present (O d).

**On top of p5 f3: dense Mn d, O d, Mn g (same day, same settings; `pol/p5f3_{dD,Od,g}_*`; every arm in the same basin).**

| U − U0 | Mn d moment | Mn d total | O p total | O normal (eV) | Mn d band (eV) | A(U0) vs p5f3 |
|---|---|---|---|---|---|---|
| PW 80 Ry (reference) | +0.0385 | −0.1307 | +0.0825 | −1.725 | −0.769 | — |
| p5 f3 | +0.0262 | −0.1228 | +0.0786 | −1.729 | −0.773 | 0 |
| + Mn d midpoints 0.56 … 24.7 | +0.0274 | −0.1252 | +0.0801 | −1.724 | −0.776 | −3.1 mHa |
| + O d 0.8, 2.0 | +0.0303 | −0.1255 | +0.0804 | −1.736 | −0.772 | −7.6 mHa |
| + Mn g 1.0, 2.5 | +0.0261 | −0.1229 | +0.0786 | −1.729 | −0.772 | −0.1 mHa |

**The Mn spin response climbs MONOTONICALLY as the basis completes**: +0.0229 (baseline) → +0.0241 (Mn p) → +0.0262 (Mn p5 f3) → +0.0303 (+ O d), with dense Mn d a further +0.0012 and Mn g (which cannot contribute) nothing; the charge columns track toward PW in step and the bands stay on PW.  That is the signature of orbital-basis incompleteness — each polarization / radial function that CAN carry the response closes part of the gap — and the first spin-response mover is O d, i.e. the Mn–O bonding region's polarization, which a single-element bisection on the minimal basis could not see (O p alone, without Mn polarization, moved nothing).  Remaining 0.008: plausibly the combination (O d + dense Mn d, untested together) and further O polarization.

**The U-free α·P test, first half (same day; `~/Code/qchem6-runs/MnO2/alpha/`).**  A static α·P = 0.5 eV on site 0's Mn d manifold (`alpha_eV`, QE's Hubbard_alpha), every U = 0, FREE (an imposition would star-average the two Mn sites), against the α = 0 references `control_free_U0` / `Mn_p5f3_U0`:

| α·P, U = 0 | Δm(Mn0, perturbed) | Δ\|m\|(Mn1) | Δn(Mn0 d) | Δn(Mn1 d) | ΔA |
|---|---|---|---|---|---|
| GPW control | −0.0023 | +0.0215 | −0.0501 | −0.0067 | +92.7 mHa |
| GPW p5 f3 | −0.0025 | +0.0196 | −0.0517 | −0.0060 | +92.1 mHa |

ΔA ≈ α·Tr n (0.5 eV × ~5 e = 92 mHa), as it must be.  The Gaussian response to a fixed potential is nearly BASIS-INSENSITIVE between the two (≈10 %), unlike the U response (+0.0229 vs +0.0262).  **The PW half did not converge**: `alpha_pw` (80 Ry, same α) sloshed for 59 iterations at Δρ 1–6e-2, ~0.3 Ha above its expected −92.59, under the Kerker/Pulay recipe that converged every SYMMETRIC PW run — the site-asymmetric perturbation needs a gentler recipe.  Stopped.  NEXT (OpenWork row): rerun PW α with startingRelaxRo ~0.2 (or a smaller α, 0.1–0.2 eV, read as the linear response), ideally seeded from the PW U0h state (`state.save` on the U0h deck, `restartFrom` on the α deck); the decisive comparison is Δm / Δn per unit α, GPW vs PW.

**The U-free α·P test, second half: PW (same day, after the PW XC-ball fix `c1ffd6ad`; `~/Code/qchem6-runs/MnO2/alpha/{pw_U0h_fix,alpha_pw_fix}`).**  Same α·P = 0.5 eV, 80 Ry, the U0h recipe; cold IonicSAD at `startingRelaxRo` 0.2 (the 0.45 run sloshed; a PW state cannot be saved yet, so no seed).  α = 0: −92.6795194582 Ha in 36 iterations (pre-fix −92.6795194034: the XC-ball fix is 0.05 μHa and ≤ 5e-4 e in ns on this cell at 40 Ha); α: −92.5870916146 in 41, converged without sloshing.  Same extraction as the Gaussian rows (the last `[+U]` line, Mn d manifold, m = N↑ − N↓; reproduces both rows above exactly):

| α·P, U = 0 | Δm(Mn0, perturbed) | Δ\|m\|(Mn1) | Δn(Mn0 d) | Δn(Mn1 d) | ΔA |
|---|---|---|---|---|---|
| **PW 80 Ry** | −0.0006 | +0.0206 | −0.0532 | −0.0064 | +92.4 mHa |
| GPW control | −0.0023 | +0.0215 | −0.0501 | −0.0067 | +92.7 mHa |
| GPW p5 f3 | −0.0025 | +0.0196 | −0.0517 | −0.0060 | +92.1 mHa |

**The Gaussian response to a plain potential on Mn d is PW's**: the induced moment on the unperturbed Mn is 104 % / 95 % of PW's, the perturbed site's charge 94 % / 97 %, the neighbour's charge 105 % / 94 % — against 60–79 % for the U spin response (+0.0229 … +0.0303 vs +0.0385).  (Δm on the perturbed site is a small difference of two near-equal spin channels, −0.0006 vs −0.0024; the one column that does not match, and the smallest.)  So the U deficit is NOT a generic incompleteness of the basis's response in the Mn d charge channel.  What α·P does not excite and U does: α·P is spin-INDEPENDENT and uniform over the d manifold; the U potential U(½ − n_mσ) is spin- and orbital-RESOLVED — larger on the minority / empty slots.  The deficit lives in that channel: a spin- (or slot-) resolved α (α on site 0 Mn d ↑ only) is the discriminating test, and needs a deck/HT extension (`alpha_eV` is a scalar on the manifold today).

**The spin-resolved α (same day; code `5c52dc35`: `alphaUp_eV` / `alphaDown_eV` per manifold, ADDED to `alpha_eV` in their own spin channel; decks `~/Code/qchem6-runs/MnO2/alpha/{alpha0_gpw_control,alphaUp_gpw_control,alphaUp_pw}`).**  α↑ = 0.5 eV on site-0 Mn d spin UP only (site 0 is the up-majority Mn), every U = 0; GPW control vs PW 80 Ry, each against its own α = 0 run on the merged binary (GPW −93.3837352286 in 27 iterations, re-run because the merge moved XC onto ρ[D]; PW `pw_U0h_fix`).  Same extraction:

| α↑·P, U = 0 | ΔN↑(Mn0 d) | ΔN↓(Mn0 d) | Δm(Mn0) | Δ\|m\|(Mn1) | Δn(Mn0 d) | Δn(Mn1 d) | ΔE |
|---|---|---|---|---|---|---|---|
| **PW 80 Ry** | −0.3115 | +0.3940 | −0.7053 | +0.3869 | +0.0823 | −0.0057 | +61.67 mHa |
| GPW control | −0.3217 | +0.3956 | −0.7175 | +0.3787 | +0.0741 | −0.0149 | +61.90 mHa |

**The Gaussian SPIN response to a one-spin potential is PW's too**: Δm(Mn0) 102 %, Δ\|m\|(Mn1) 98 %, Δn(Mn0 d) 90 %, both spin channels of the perturbed site within 3 %.  So neither the charge channel (spin-independent α) nor the spin channel (α↑) shows the 60–79 % U deficit; what is left is the ORBITAL (slot) resolution of U(½ − n_mσ) — or the U runs themselves.

**But 0.5 eV of α↑ is NOT linear response — both codes leave the reference state.**  The response is huge (a 0.5 eV push on the up channel moves 0.39 e INTO the minority down channel of the same site; m0 3.52 → 2.80, \|m1\| 3.52 → 3.90: ferrimagnetic, the two Mn no longer equivalent), and the energy says it is another state: Hellmann–Feynman along any continuous branch gives ΔE = ∫ N↑(α′) dα′ ≥ α·min N↑ = 0.018375 × 3.972 = 72.99 mHa (GPW), 72.80 (PW), but the runs give 61.90 / 61.67 mHa.  Equivalently the α-FREE energy E_tot − α·N↑ of the perturbed state is **11.1 mHa BELOW the α = 0 reference, in both codes** (−11.09 / −11.13 mHa).  The same signature as polarized Γ Si under a one-spin α of 1e-3 Ha (OpenWork §3).  So the 0.5 eV comparison is a comparison of the two codes' FERRI states (and they agree), not of the linear spin susceptibility; the U spin response (+0.023 … +0.030 vs +0.0385) is still unexplained by any channel measured so far.  **RESTART CHECK (GPW, `alphaUp_gpw_save` → `alpha0_from_up`, records r010 → r011): α = 0 warm-started from the α↑ state STAYS ferrimagnetic and converges in 15 iterations to −93.3954227401 Ha, 11.69 mHa BELOW the AFM reference (−93.3837352286)** — Mn0 d N↑/N↓ 4.0365 / 1.1189 (m +2.918), Mn1 0.5688 / 4.4844 (m −3.916).  So at U = 0 the β-MnO₂ A1AFM state every U / α run has been measured from is NOT the LDA ground state of this cell: a ferrimagnetic state with the two Mn inequivalent lies 0.32 eV/cell lower, and a one-spin push of 0.5 eV is enough to cross into it.  (PW agrees on the perturbed state to 0.04 mHa in its α-free energy; a PW α = 0 restart is impossible until PW states can be saved.)  Physics or a shared defect is OPEN: both codes share the polarized LDA XC, and polarized Γ Si shows the same kind of jump — an independent code (QE, unequal `starting_magnetization`) settles it.

**2026-10-27 — both follow-ups, run in parallel.**

*QE settles it: the ferrimagnetic state is LDA PHYSICS, not a qchem6 polarized-XC defect.*  `~/Code/qe-runs/MnO2/mno2b_U0_ferri{,_m0}.scf.{in,out}` — twins of `mno2b_U0_tok` (300 Ry, same GTH PPs and cell, Gaussian smearing 0.002 Ry) seeded Mn1 +0.5 / Mn2 −0.9.  With `tot_magnetization = 0` (like our multiplicity-1 deck): −186.96074763 Ry vs the AFM −186.93546944, **ΔE = −12.64 mHa** (ours −11.69 mHa, GPW control), sphere moments Mn +2.57 / −3.40, O +0.37 / +0.09 (×2 each) — the pattern of ours (Becke: Mn +2.84 / −3.78, O +0.37 / +0.10).  Net moment free: −12.70 mHa with a net −2.00 μB/cell (O moments flip sign).  So the A1AFM state of this cell is a LOCAL minimum at U = 0, and a 0.5 eV one-spin push crosses the barrier in all three codes.

*Linear spin response (α↑ = 0.05 eV, decks `alphaUp005_{gpw_control,pw}`), on the AFM branch:* ΔE = 7.8887 mHa (GPW) vs α·N̄↑ = 7.8888, and 7.8511 (PW) vs 7.8510 — first order to 0.1 μHa, so no jump.

| α↑·P 0.05 eV, U = 0 | ΔN↑(Mn0 d) | ΔN↓(Mn0 d) | Δm(Mn0) | Δ\|m\|(Mn1) | Δn(Mn0 d) | Δn(Mn1 d) |
|---|---|---|---|---|---|---|
| **PW 80 Ry** | −0.0016 | −0.0010 | −0.0006 | −0.0091 | −0.0026 | +0.0009 |
| GPW control | −0.0014 | −0.0012 | −0.0002 | −0.0097 | −0.0026 | +0.0007 |

The perturbed site's moment barely moves (both spin channels empty together: the on-site spin channel is stiff); the response is on the NEIGHBOUR's moment, Gaussian 107 % of PW, and the site charge 100 %.  **So in linear response too, the Gaussian spin channel is PW's.**  Charge (α), spin (α↑) — neither carries the 60–79 % U deficit; what U excites and neither probe does is the ORBITAL (slot) resolution of U(½ − n_mσ): next, a per-slot α.

**The slot-resolved α (same day; code `ef056880` + `63d7ea36`: `alphaSlot_eV`, one static α_k on each U slot's FIXED isotypic projector P_k, both spins; slot occupations now printed to 4 decimals).**  0.05 eV on ONE slot of site-0 Mn d at a time (slots [0],[1],[2] dim 1, [3] dim 2), U = 0, GPW control vs PW 80 Ry, decks `alphaSlot{0..3}_{gpw_control,pw}`, references re-run on the same binary for the 4-decimal print (`alpha0_gpw_control_r4` −93.3837352286, `pw_U0h_fix_r4` −92.6795194582: both reproduce exactly).  All eight on the AFM branch: ΔE = α·Tr P_k n to first order (GPW 2.051 / 1.837 / 1.703 / 3.714 mHa vs 2.053 / 1.837 / 1.706 / 3.720; PW 2.055 / 1.836 / 1.714 / 3.676 vs 2.056 / 1.837 / 1.717 / 3.682).  Responses (units 1e-4 e; the print resolution is 1e-4):

| α on slot | Δm(Mn0) GPW / PW | Δ\|m\|(Mn1) GPW / PW | ΔN(Mn0 d) GPW / PW | own slot ΔN↑ + ΔN↓ GPW / PW |
|---|---|---|---|---|
| [0] dim 1 | +3 / +6 | +17 / +17 | −7 / −4 | −15 / −16 |
| [1] dim 1 | −1 / 0 | +2 / +2 | −1 / +2 | 0 / −1 |
| [2] dim 1 | −14 / −13 | +17 / +18 | −12 / −11 | −36 / −37 |
| [3] dim 2 | +9 / +7 | −13 / −10 | −31 / −35 | −62 / −65 |

**The Gaussian ORBITAL response of Mn d is PW's too**, entry by entry within the print resolution (every per-slot, per-spin, per-site change agrees to ≤ 3e-4).  So the Mn d → Mn d response is PW's in every channel probed: charge (α), spin (α↑), orbital (α_k).

**But that was the wrong row of χ.**  Re-reading the deficit's own record (2026-10-25, "THE PLANE-WAVE VERDICT"): the U arm carries U = 4.62 eV on the **O p slot Γ2 (the plane normal)** and U = 0 on Mn d.  The "Mn d moment" deficit (+0.023 … +0.030 vs +0.0385) is therefore a CROSS response — the Mn d moment's response to a potential on O p_z — which no Mn d α probes.  The discriminating probe is α on the O p Γ2 slot (all four O, as the U arm), reading the Mn d moment; and since that U is 4.62 eV, compare its linear response before blaming a nonlinear one.

## 2026-10-25 — HubbardProjectorPlan increment 4: ACBN0 + response through the pseudo-atomic projector, re-judged

Code `b4cfc789` (`GaussianAtomCentredSet::BareCoulomb`, exact one-centre).  Decks + logs: `~/Code/qchem6-runs/TiO2/basis_scan/
b6_TiO_pseudoAtomic{.json,_20261025.log}` (record `TiO2_rutile.r003`), `~/Code/qchem6-runs/NiO/nio_q0_pseudoAtomic{.json,_20261025.log}`
(`NiO_AFM2.r009`, warm start from the saved FREE state, 11 iterations).  Each is its site-shell twin with `atomicRadial` → `"projector":
"pseudoAtomic"` (ortho-atomic over the set), same binary otherwise.

| | site-shell (10-24) | pseudo-atomic | hp.x / QE |
|---|---|---|---|
| TiO₂ b6 Γ, U(Ti 3d) 6×6 eV | 3.0504 | **2.9945** | 2.927 |
| TiO₂ b6 Γ, U(O 2p) 6×6 eV | 6.7345 | **6.8226** | 6.758 |
| TiO₂ b6 Γ, U(Ti 3d) 1×1 eV | 5.1265 | **5.0880** | 5.026 |
| TiO₂ occupation Ti 3d / O 2p | 2.211 / 4.937 | **2.210 / 4.931** | 2.208 / 4.931 |
| NiO q = 0, U 2×2 eV (Ni0 / Ni1) | 6.4677 / 6.4677 | 6.4793 / 6.4792 | (q-mesh: step 3) |
| NiO χ₀ Ni0 / E_tot | −4.02564 / −106.1384955 | −4.03786 / −106.1385368 | |

Captured norms B†S⁻¹B/O: TiO₂ 0.9999 (Ti 3d) / 1.0000 (s, O), NiO 0.9998 -- so the projectors are close to the site-shell ones on these
bases and the moves are small; TiO₂'s Ti halved its hp.x excess and the occupations reached QE's.  Si p (ACBN0, SR basis): bare U 10.6 eV
pseudo-atomic vs 23.6 eV for the site-shell contraction -- there the SR p shells could NOT hold the pseudo-atom 3p, and the contraction
was a different, compact function: the basis dependence the plan removes.

## doc/HubbardProjectorPlan.md -- retired 2026-10-25 (increments 1–5 landed), verbatim

(Source comments still cite `doc/HubbardProjectorPlan.md` D1–D7 / increment N: this is that text.)

### +U projectors independent of the orbital basis

**Born 2026-10-24, RETIRED 2026-10-25 (increments 1–5 landed)** (user: *"clean code is MUCH higher priority"* than a cheap route;
*"if we have to redo some runs that is fine"*).  The evidence is in `Records/HubbardUHistory.md` (2026-10-24
entries); the tracker row is `OpenWork.md` §2 "+U projector independent of the orbital basis".  When this
plan is done its text goes to `Records/HubbardUHistory.md` and this file is deleted.

#### 0. Status (2026-10-24, end of the first execution session)

- **Increments 1–3 LANDED** (each sweep-green, 0 of 1071 anchors moved at 1e-10): D1–D4 as written, except O
  comes from the BLOCK (`AtomCentredOverlapSource::MakeSetOverlap`) through a new engine pair kernel
  `LatticeSum1E::MakeOverlap(phase,A,g,h)` — the k convention and the lattice sums are the block's side, and the
  orbitals do not enter.  D5 = `qchem.PseudoAtomOrbital`; D6 energy/potential/slots via one synthetic
  `SphericalShellRep` shell per manifold; D7 = `"projector": "pseudoAtomic"` (default still `siteShells`).
  The set's m convention = `SphericalShell(l)[l+m]` at POSITIVE unit-sphere scale (GPW's `YlmCartesian` now
  takes |N|; bitwise-neutral for the KB sum).
- **The increment-3 gate FAILED and refuted §1's premise**: the gap does not move with the projector (VA 0.128
  → 0.128 eV; QE 0.325); the diffuse-shell gain was the orbital basis.  The projector sets E_U (QE's to 0.1 % on
  VA + O p 0.18).  Numbers: `Records/HubbardUHistory.md` 2026-10-24 "The pseudo-atomic +U projector".
- **User, 2026-10-24**: the new projector flexibility is worth keeping; the gap is to be looked for elsewhere.
  Increments 4–5 remain the plan (basis-independent U for ACBN0 / hp.x comparisons; QE's E_U).
- **Increment 4 LANDED 2026-10-25** (`b4cfc789`, fast tier 1081/1081, 0 anchors can move: the site-shell path is untouched).
  D6's ACBN0 bare integrals are the SET's own `BareCoulombSource` (one centre, exact: multipole expansion + real Gaunt
  factors + Slater integrals of Gaussian densities in closed form) -- NOT the molecular engine, which can only reach a
  pseudo-atomic radial primitive by primitive (~10^7 ERIs; CleanCode D-CONTRACTED-SHELL), and it duplicates `RkEngine`
  (CleanCode D-RK-ONCE).  The engine is the oracle: element by element to 1e-12 (d primitive, contracted d, p+d).
  `ManifoldCoefficients` and the response needed nothing (`Atomic`'s Q = T^dagger is the site-shell route's ortho convention).
  Gates: `BareCoulomb.TheAtomCentredSet*`, `GPW_Si.Γ_U_PseudoAtomic3p_Imp_Pol_eqUnpol`; physics re-judged
  (`Records/HubbardUHistory.md` 2026-10-25): **TiO₂ b6 6×6 vs hp.x -- U(Ti) 2.994 eV (+2.3 %, was +4.2 % site-shell),
  U(O 2p) 6.823 (+1.0 %, was −0.35 %), Ti 1×1 5.088 (+1.2 %, was +2.0 %); occupations Ti 3d 2.210 / O 2p 4.931 vs QE
  2.208 / 4.931.  NiO cell (q = 0, U_in 3): U 6.479 eV vs site-shell 6.468 (+11 meV; the VA site shells already
  capture 0.9998 of the 3d).**  The NiO-vs-hp.x q-mesh comparison is step (3) of OpenWork §3's ORDER.
- **Increment 5 LANDED 2026-10-25 -- THE PLAN IS DONE.**  Part A `095f0bec`: every CP2K-gated / LOWDIN manifold names
  `siteShells` (`HubbardU_SiteShells`; 14 decks); 0 of 1089 / 5 moved.  Part B `94ed536f`: `PseudoAtomic` is the default,
  `radial` / `atomicRadial` and the whole contraction path retired, 38 decks / 338 manifolds migrated, a deck asking for a
  contraction is refused, WARM `hubbard.projector` replaces `hubbard.atomicRadial`.  0 of 1084 / 5 moved vs part A (the
  first sweep caught two Si INLINE JSON decks that the part-A search missed: a LOWDIN manifold had turned pseudo-atomic).
  `GPW_MnO.Γ_U_Shub_Pol_Smear_CP2K_Long` did not move (it names `siteShells`).  REMAINDERS, all elsewhere: CleanCode
  D-RK-ONCE, D-CONTRACTED-SHELL; the campaign re-runs = OpenWork §3 ORDER step (3).

#### 1. Why

A contracted +U manifold today is χ_m = Σ_s r_s φ_{s,m} over the **site's own orbital shells**
(`HubbardManifold::radial`, `Hubbard_U::SelectOn`), and `SolidCalculation::AtomicRadial` fits the GTH
pseudo-atom orbital onto those shells.  So the projector is whatever the site's shells can hold:

- **β-MnO₂ (VA, Γ, U(O pz) 4.62 eV).**  The O 2p projector captures **0.976** of the pseudo-atom orbital
  (norm beyond 2.5 bohr: 2.1 % vs the orbital's 6.6 %).  The gap is 0.128 eV against QE's 0.325, and every
  occupied O p band sits 0.2–0.5 eV too high relative to Mn d.  U = 0 agrees with QE to 0.04 eV.  Orbital,
  energy form (Dudarev in the eigenbasis) and E_U all match QE.
- **One diffuse O p shell (0.18)** raises the captured norm to 0.9994 and the gap to 0.237 eV.  The rest is
  the projector tail's overlap with the **neighbours'** functions, which no site-only contraction can hold.
- **The diffuse shells are the ones conditioning removes**: MnO VA dropped d 0.18 (Mn 3d captures 0.991),
  and the NiO vet trims Ni s 0.06 and d 0.18.  A projector never enters S, so it need not lose them.
- **Today U changes with the basis**, which confounds the S2 basis scan and any hp.x comparison.

QE and hp.x define the projector as an atomic function (the UPF's `PP_CHI`, which our own `gth2upf` wrote
from the same pseudo-atom).  We should define it the same way.

#### 2. The design: the separable-PP pattern, applied to +U

The precedent already exists in this tree.  The nonlocal PP's projectors β are **not** basis functions:
`SpeciesProjectorSet_Gaussian` holds them as Gaussians × Y_lm, and `GPW_Evaluator` builds
b_i = ⟨χ_i^k|β⟩ per block through `LatticeSum1E::MakeOverlap(phase, cell, GaussianFunction)`
(`GPW_Evaluator.C` ~1440).  The magnitude-screened lattice sum follows pin 1 (no r-space cut), and the
symmetry-adapted and spherical blocks come for free: the block's own lattice face, including
`SphericalLatticeView`'s congruence, returns b over the block's functions.  A +U projector is the same kind
of object.

**D1. A vocabulary type: the atom-centred function set** (`BasisSet` library, beside `Orbital_PP_IBS.C`'s
vocabulary; data-free abstract face plus a Gaussian capability, the house style).  It holds, per manifold:
the centre, l, the real-harmonic m convention (the same as the spherical view's), and the radial as a
closed sum of Gaussians (`AsGaussians`).  No weights: that is the difference from `SpeciesProjectorSet`.
It is per **site** rather than per species, because manifolds are per site.  Check first whether a common
base with `SpeciesProjectorSet` falls out naturally; don't force one.

**D2. A new orbital-basis service, named for what the client consumes**: *"the overlap of every basis
function with each function of an atom-centred set"*, B (n × m) per block, Hermitian-consistent with
`MakeOverlap`.  Like `Orbital_PP_IBS`, it is an integral a basis can compute knowing nothing about
Hubbard.  Implementors: the GPW block (through `LatticeSum1E::MakeOverlap(phase, cell, g)` per Gaussian
term, as the separable PP does; refactor that loop into the shared helper rather than copy it) and the
spherical view (congruence).  A block without the face throws cleanly (pin: a call that can legitimately
fail returns an Outcome; this one is a broken invariant, so it throws).

**D3. The projectors' own overlap O (m × m), for ortho-atomic.**  QE's ortho-atomic is
φ̃ = φ O^{-1/2} with O the functions' **true** mutual overlap, lattice-summed.  Build the set's own lattice
overlap from the same analytic Gaussian lattice sums.  The clean route is to let the function set realise
a lattice overlap face over itself, which must not reach into an orbital basis.  Also report the
represented overlap B^† S^{-1} B beside it, as the "captured" diagnostic: how much of each projector the
orbital basis spans.  **Decision for the session**: confirm true-O (QE semantics) as the definition; the
represented O is a diagnostic only.

**D4. `LowdinProjector` from B and O**: T_M = B_M Õ^{-1/2}, Q_M = T_M^† (⟨φ̃|ψ⟩ = T^† c).
Forward (n = T^† D T) and Adjoint (V = T W T^†) are unchanged in form.  The column manifold
(T = S^{1/2}[:,M], CP2K's LOWDIN convention) stays as the CP2K-parity path.  **The site-contracted
radial (`radial`, `atomicRadial`) retires**: the atomic function set supersedes it.  Pin the ctor's two
contracts in unit tests: Forward/Adjoint exactly adjoint, and an idempotent density inside the projector
span giving occupations 0/1.

**D5. Who builds the function set**: the GTH pseudo-atom orbital in **its own** exponents (today's
16-exponent pool, `AtomCalculation`, LDA, unpolarized), with no projection onto the site.  It moves out of
the `SolidCalculation` facade into a builder of its own (one pseudo-atom per species, cached), reachable
from the deck.  The pool's adequacy is a gate (§4), not an assumption.

**D6. The consumers, all through the term's existing faces**:
- energy and potential (`Analyse`/Dudarev, `Forward`/`Adjoint`): no change beyond T;
- **site symmetry and slots** (`ManifoldSymmetry`): built today from the picked shells' rep matrices.
  Build it from the function set's (centre, l, m convention) instead;
- **ACBN0 bare integrals** (`ManifoldIntegrals` → `BareCoulombSource::BareCoulomb(cols)` + `Contraction`):
  the on-site ERI4 over the projector functions themselves.  The function set is a one-centre Gaussian set,
  so its 4-centre integrals come from the molecular ERI engine.  Give it the `BareCoulombSource` face rather
  than a special path;
- **ACBN0 coefficients** (`ManifoldCoefficients` = Q C) and the **response kernel** (`SlotProjectors`,
  Forward/Adjoint): no change beyond T and Q.

**D7. Deck**: a manifold names its projector kind (pseudo-atomic, the new default; site columns for CP2K
parity).  `radial` and `atomicRadial` retire.  The decks that use them migrate, found with `scripts/vet-decks`
plus a grep of `~/Code/*-runs` decks.  `doc/InputDeck.md` is updated.

#### 3. Increments (each one sweep-green and committed before the next)

1. **D1 + D2 + unit tests.**  B against a brute-force check: a basis that *contains* the projector
   function must give B = S c exactly; the separable-PP refactor must be bit-identical (gate:
   `AnalyticSeparablePPMatchesMesh` plus the full sweep with 0 moves).
2. **D3 + D4**, with `LowdinProjector` tests (adjointness, 0/1 occupations, ortho-atomic orthonormal as a
   set).  Not wired into a run yet.
3. **D5 + D6 (energy, potential, symmetry/slots) + D7**, behind the new projector kind, opt-in.
   **Gate: β-MnO₂ Γ on plain VA (no diffuse shell): gap → 0.325 eV (QE), O p bands within ~0.05 eV of QE's
   relative to Mn d** (`QCHEM_DIAGNOSTICS=frontier=16`; QE band characters in HubbardUHistory 2026-10-24).
   And basis independence: VA vs VA + O p 0.18 give the same occupations to ~1e-3.
4. **D6 (ACBN0 bare integrals + response)**.  Gates: the ACBN0 unit tests, and TiO₂ / NiO U vs hp.x
   re-judged (OpenWork §3 "NiO / SrVO₃ U vs hp.x").
5. **Make it the default; retire `radial`/`atomicRadial`**; re-judge every +U anchor (pin 10: against an
   independent route, never refreshed): `GPW_MnO.Γ_U_Shub_Pol_Smear_CP2K_Long` stays on the CP2K-parity
   column path, so it should **not** move; any anchor that does move needs a reason.
   Then the campaign runs the other session re-ran on 2026-10-24 are redone (user: fine).

#### 4. Risks and open questions (settle in the increment that meets them)

- **Is the 16-exponent pool good enough as the projector?**  It reproduces CP2K's ATOM to 3 / 0.3 / 0.1 mHa
  (Mn; `UTCalculation ValenceBasisGen.MnQ7PseudoAtomInALargePool`).  Check the O 2p tail against the UPF's
  `PP_CHI`, which this tree generated from the same pseudo-atom: normalized overlap ≈ 1.
- **Screening**: a diffuse projector widens the lattice sum.  The magnitude screen handles it, but measure
  the setup cost on MnO (Benchmark's ledger, a setup bucket of its own).
- **Real TRIM blocks / the spherical view**: B must come back real on real blocks, as `MakeOverlap` does.
- **Charge conservation of n**: with true-O ortho-atomic, the trace of n over all manifolds is no longer
  bounded by N.  This is expected (QE's is not either); don't assert it.
- **`dynamic_cast` discipline**: the new service is reached abstract→abstract (the `Orbital_PP_IBS`
  pattern); flag any concrete cast.
