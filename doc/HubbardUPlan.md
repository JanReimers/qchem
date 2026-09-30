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

⛔ **FeS₂ IN PROGRESS (2026-09-30, batch queue, `queue/02_fes2_hpx.sh`)**: structure pulled directly from
`ct3c01403.pdf` (Macke et al., the orbital-resolved-U paper itself — cubic Pa-3̄, a=5.418 Å, the paper's own
cited EXPERIMENTAL structure, not their PBE-relaxed one, since qchem6 has no GGA/relaxation yet to reproduce
that); S 8c coordinates hand-derived from Pa-3̄ symmetry at x=0.386, cross-checked against CP2K's own
`c_1_FeS2.inp` regression deck to <0.003 (that deck's coordinates are intentionally perturbed for a
symmetry-finder test).  **S converts cleanly via `gth2upf`** (E_atom −10.069183 Ha) — the "S untested"
open item from A6's scope is resolved.  SCF (nonmagnetic, `occupations='fixed'`, no 2-step) **converged in
39 iterations**, E = −330.04024036 Ry — confirms the diamagnetic low-spin Fe²⁺ d⁶ (t₂g⁶eₘg⁰) prediction, no
smearing needed.  `hp.x` (shell-averaged Fe-3d, `ortho-atomic`, q 2×2×2) running now; U value to follow.

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
