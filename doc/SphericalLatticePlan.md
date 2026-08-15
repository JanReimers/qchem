# Spherical d on the lattice path — the fixed-transform decorator plan

## ★ SESSION START HERE (updated 2026-08-14) — the V_long "defect" DISSOLVED; collapse cause reopened

**WHERE THE CAMPAIGN STANDS:** I1 (spherical lattice view) DONE + gated; I3 (valence_lowq_sph = the
EXACT CP2K-transcription span, n=136) DONE; I2 arm 1 DONE — **the ordering HEALED under spherical d
(AFM below FM by 45.5 mHa; span convicted, no code bug in the KB story)**.  I2 arm 2 (the CP2K-span
absolute A/B) is still blocked on the diffuse-re-addition collapse — but the 2026-08-13 "unit-reproduced
V_long defect" is RETRACTED (below), so the collapse cause is back under investigation.

**THE 2026-08-13 "DEFECT" WAS THE GATE'S OWN WIRING (found+fixed 2026-08-14).**
`GPW.DISABLED_DiffuseDVlongSharpFieldOracle` passed `gthMn.local` — a SINGLE-SPECIES
`HGH_LocalPotential`, whose `FormFactorLong(int /*Z*/,…)` IGNORES Z — over the Mn+O cell.  Two
consequences: (1) production integrated a Mn-q7 field (Zion 7, rloc 0.64) AT THE O SITE while the
oracle integrated the true O q6 well (Zion 6, rloc 0.2476) — different physics, guaranteed red;
(2) Mn q7 has no C terms, so `ShortRangeGaussian` returned empty ⇒ β=0 ⇒ `MakeLocalPPLong` bailed to
the kappa sweep — **the custom-G-ball path under suspicion was never even entered** (also why
`GPW_VLOC_EPS=1e-9` and `GPW_LONG_SWEEP=1` A/Bs were BIT-IDENTICAL: 0.0138439 all three).  With the
gate fixed to the `MultiSpecies_LocalPotential` router (exactly production's wiring —
`Hamiltonians.C BuildMultiSpeciesLocal`; the arm-2 runs were NEVER species-mismatched):
- fixed sharp gate, custom G-ball at true β=8.15 (its FIRST genuine unit exercise): **GREEN, 1.8e-5**;
- fixed sharp gate under `GPW_LONG_SWEEP=1` (kappa path at β=8.15): **GREEN, 1.8e-5** — same worst
  element, prod identical to 8 digits (and 477 s vs 131 s: the retirement-grade cost, as documented);
- Mn-only soft gate: still green, 6.5e-7.
**VERDICT: the V_long integrate-back (per-level restriction, harmonic routing, both paths) is
EXONERATED at unit tier.**  The 55B "−13.3 Ha in E_loc" was a density-response symptom (a collapsed
density legitimately deepens E_loc), never an operator conviction.

**THE COLLAPSE: CONDITIONING-CLASS, CONVICTED (run 58, 2026-08-14,
doc/logs/mno_probe_run58_orthotol1e-2.log).**  Variant B (Mn d 0.18) dove to −67.28@4 at
MNO_ORTHO_TOL=1e-3 — which, it turns out, dropped ZERO of the 128 functions (the pivot tolerance
never triggered; 55B solved the FULL near-dependent span, λmin 4.6e-4).  At MNO_ORTHO_TOL=1e-2 the
filter drops 6 of 128 at the door ("~reproducible by the kept set") and the dive is GONE: iters
−14.73 / −59.61 / −61.373 / **−61.3905@4**, m_stag +0.644, m_net 0.008, breakdown healthy
(Een −83.61 ≈ 55A's −83.05) — right on 55A's plateau trajectory (−61.4067@4).  The 55B −67.28 was
PROVABLY unphysical (a subset of the CP2K span cannot sit 5.8 Ha below CP2K's variational
reference), so it was numerical collapse: near-null modes × eps-screened operator inconsistencies,
exactly the −455 mechanism in milder form.  NO operator defect — the KB-nonlocal suspicion is
DROPPED for this symptom.  (Deeper cure than a stiffer pivot — tightening screens with cond(S), or
SVD-consistent F/S filtering — is a design question, not this campaign's blocker.)

**RUN 59 (AFM arm) DONE 2026-08-14: E = −61.56504, m_stag ±0.678 (m_net −0.005), NO collapse.**
Trace highlights: smeared stage ended ITSELF at iter 11 via the new ladder-exhaustion bail-out
(A=−61.5727, −TS=12.8 mHa); cold GDM stage settled −61.564575 (15 iters at |ΔE/E|~1e-9), then a
SILENT GDM fallback-diagonalize hopped the determinant at the tied metallic frontier (iter 25:
+302 mHa, cfg change, m_stag→0.59 — MNO_MOM=0 leaves the tie unpinned), GDM ground back over ~30
iters and landed 0.5 mHa BELOW the pre-hop floor with a REAL ~5 mHa aufbau gap (the pre-hop state
was a tie-cycle above the true minimum — the hop was accidentally productive).  Hit the 80-iter cap
with energy settled (relAmp 1.7e-7), Δρ tail 3.7e-4 > the 1e-5 gate ("raise NMaxIter, NOT a floor")
⇒ formal conv=0, physics done.  TWO instrument notes banked: the GDM fallback-diagonalize needs a
breadcrumb line (it caused a +302 mHa excursion with zero console trace), and the Δρ gate should be
INTENSIVE (Δρ/N — user; doc/SCFStrategyPlan.md).
**HEADLINE: we sit ~94 mHa BELOW the CP2K Γ AFM oracle (−61.4706) on a 126-of-136 SUBSET span.**
Before calling that an operator/grid residual: CP2K's printed total under SMEAR ON includes its
electronic-entropy term — the E-vs-A split and kT of the oracle runs must be pinned down first;
then the term-by-term breakdown (Ekin/Eee/Exc/E_loc/E_NL vs CP2K's energy blocks) localizes
whatever remains.

**RUN 60 (FM arm, v2 span) DONE 2026-08-14: E = −61.6011060, CONVERGED (58 cold iters, Δρ gate
passed; ends with a 3 mHa tie-hole NON-AUFBAU warning).**  v2-span table: AFM −61.5650 / FM
−61.6011 ⇒ FM below by 36.1 mHa (CP2K same file: AFM below by 8.87) — ordering reversed vs CP2K,
offsets 94/139 mHa below their AFM/FM.

**THE EXACT-SPAN SYMMETRIC MATRIX (2026-08-15, the hand-waving killer).**  qchem's rank filter
drops m-RESOLVED AO components CP2K cannot express (and asymmetrically per site — the O₁-p/O₂-s
catch), so the exact A/B was moved to SHELL-level spans both codes hold FULL RANK: variants VB
(N=128) and VA (N=118) banked as CP2K decks + basis entries (IntegrationTests/CP2K/*_v{a,b}.inp,
VALENCE-LOWQ-V{A,B}).  CP2K results (Γ, kT=5e-3, minutes each):
  full-136  AFM −61.47057 / FM −61.46170 ⇒ Δ = −8.87 mHa (AFM wins)
  VB 128    AFM −61.34454 / FM −61.33635 ⇒ Δ = −8.20 mHa (AFM wins)
  VA 118    AFM −61.30333 / FM −61.30478 ⇒ Δ = **+1.46 mHa (FM wins — CP2K's own ordering FLIPS)**
so the Mn d(0.18) shells carry ~10 mHa of AFM stabilization in CP2K — larger than the physical
6J₁+12J₂ ≈ 4 mHa scale (user, susceptibility), i.e. NONE of these Γ/LSDA orderings are
physics-grade; the matrix's purpose is CODE-vs-CODE defect localization only.
**qchem VA pair (runs 61/62: MNO_SHARED_MU=1 = CP2K's one-μ ensemble, MNO_IMPOSE=1,
MNO_ORTHO_TOL=1e-3 = full rank 118 — verified "kept 118 of 118" both arms):**
  run 61 AFM: **E = −61.4029762, CONVERGED in 14 cold iterations** (no hop, no tie grind),
  m_stag ±0.667 (still the weak-moment basin).  ⇒ **99.7 mHa BELOW CP2K's VA AFM with ZERO
  span/symmetry/ensemble excuses — the code-vs-code offset is real and operator-level.**
  run 62 FM: **E = −61.4415831, CONVERGED in 14 cold iterations** (ends NON-AUFBAU with a 12 mHa
  hole — the FM number may sit up to ~10 mHa above its own minimum; caveat carried).
**★ THE VA VERDICT TABLE (the sharpest defect coordinate this campaign has produced):**
  qchem  VA:  AFM −61.40298 / FM −61.44158 ⇒ Δ(AFM−FM) = +38.6 mHa (FM wins)
  CP2K   VA:  AFM −61.30333 / FM −61.30478 ⇒ Δ = +1.46 mHa (FM wins)
  offsets:    AFM −99.7 mHa / FM −136.8 mHa ⇒ **configuration-selective part −37.2 mHa**
Both codes order FM first on VA, but qchem carries (a) a ~100 mHa configuration-BLIND offset
(below CP2K — cannot be span, we matched it; suspects: G=0/alignment conventions, XC quadrature,
V_loc bookkeeping) and (b) a **~37 mHa configuration-SELECTIVE FM-favoring bias that is
SPAN-INDEPENDENT** (v2: Δδ ≈ 45 mHa; VA: 37 mHa; the I0 d-selective signature) — one stable
number, zero excuses, waiting for the term-by-term breakdown to name its operator.
**BREAKDOWN FORENSIC (2026-08-15).**  Exc directly (different methodologies — our site-adapted
Becke vs CP2K's uniform grid): AFM −14.4423 vs −14.4513 (+9.0 mHa), FM −14.3608 vs −14.3729
(+12.1 mHa) — real ~10 mHa methodology difference but WRONG SIGN for the blind offset, and only
3 mHa of the 37 mHa selective bias.  Δ(AFM−FM) composition: qchem Ekin +1.627 / Eloc −0.671 /
**E_NL −1.306** / Eee +0.470 / Exc −0.082 (Σ +38.6 mHa) vs CP2K CoreH +0.330 / Hartree −0.256 /
XC −0.078 / entropy +0.006 (Σ +1.5 mHa).  **Δ_NL = −1.31 Ha: the I0 signature returned on a span
with NO contaminants and NO diffuse d — the arm-1 "healing" was the trimmed s-window STARVING the
mechanism, not spherical d curing it.**

**★ THE KB ORACLE GATE (2026-08-15): the assembly is EXONERATED at oracle grade.**
`GPW.DISABLED_DiffuseDKBOracle` (~25 s): production MakeSeparablePPByL (analytic path: BetaGaussian
→ Cartesian polys → lattice-sum seam) vs a from-scratch raster quadrature (explicit Bloch-summed χ
× BetaR × independent orthonormal Y_lm; shares only the diagonalised (l,D,βr) parameter layer,
itself validated by the Mn-atom oracle) on the Mn+O diffuse-d cell, MULTI-SPECIES router:
l=0 2.4e-7 / l=1 7.7e-7 / l=2 9.0e-8, median prod/oracle = 1 in every channel (no convention
drift).  With S/T analytic, V_short gated, V_long oracle-clean (2026-08-14) and KB now
oracle-clean, **Δ_NL = −1.31 Ha is a TRUE property of our converged states under a CORRECT NL
operator** — the weak basin genuinely collects it.  The offsets' remaining habitat: the
Hartree/Poisson on the band-limited collocated ρ̃ (adjoint-exact but resolution-limited — the SPIN
density's resolution is where configuration-selectivity can enter), the ~10 mHa XC methodology
split, and the E_alphaZ ↔ CP2K core-self/overlap G=0 bookkeeping (never term-mapped precisely).
All grid/convention questions — the performance reprioritization (doc/OpenWork.md items 1–3+5)
makes probing them affordable and comes FIRST; NO further long MnO runs until it lands.

**DONE: THE ARM-2 VERDICT RUNS at MNO_ORTHO_TOL=1e-2** (the one
deviation from the 53/54 recipe, per run 58): AFM (`MNO_SKIP_FM=1`, run 59,
doc/logs/mno_afm2_run59_sph_cp2kspan_tol1e-2.log) then FM (`MNO_SKIP_AFM=1`, run 60,
doc/logs/mno_fm_run60_sph_cp2kspan_tol1e-2.log), chained sequentially, each
`MNO_ANNEAL="5e-3,0" MNO_ACC="Ladder,GDM" MNO_MOM=0 MNO_ORTHO_TOL=1e-2 GPW_SPHERICAL=1
GPW_BASIS_SPH=1` + verbose + OMP 8 + MemoryMax=12G, DEFAULT screens (never relax
GPW_SCREEN/DENSITY_EPS on this span — the eps=1e-8 collapse lesson).
TARGETS = the oracles themselves (same span, function-for-function): AFM −61.4706 / FM −61.4617;
anything far BELOW −61.5 = collapse suspicion.  Open questions the verdict answers: the honest
operator/grid residuals, and whether the deep-moment basin (m̃ → ~4.65) appears on the honest span.
Caveats to carry:
- MNO_ORTHO_TOL=1e-2 solves in a rank-filtered subspace of the 136 (the full span will drop a few
  near-dependent functions) — every dropped function is span CP2K has and we don't, so a small
  legitimate positive residual vs the oracles is expected; revisit if the numbers land within a few
  mHa.
- 55D (O diffuse, borderline −61.60@4 at the old tol) is plausibly the same near-null contamination
  in milder form — optional recheck at 1e-2 if the verdict runs leave doubt.

2026-08-12.  Born from the MnO ordering campaign (doc/SymmetryUpgradePlan.md, runs 38–44): the
FM/AFM ordering reversal traced to the KB nonlocal term rewarding the weak-moment basin with 1.3 Ha
of extra d-channel attraction that CP2K — same deck, same transcribed basis — does not see.  The
LAST structural difference between the codes at that point: **qchem runs Cartesian d (6 components),
CP2K spherical (5)**.  The d-channel NL operator matrix over the shared spherical-d span is
IDENTICAL in both codes (the contaminant is pure l=0, orthogonal to the l=2 projectors), so a
spherical-d run of OUR MnO is the exact span-matched discriminator:

- ordering heals ⇒ the reversal was **differing basis-incompleteness effects** (the 8 extra
  r²e^{-αr²} s-type functions per Mn reshaping the density toward d-hybridised character) — no code
  bug, no physics; the cure is a basis-policy option, and Cartesian-d bases are convicted of a
  QUALITATIVE failure mode (they can flip which basin wins at finite basis).
- reversal survives ⇒ the l=2 crystal KB assembly (Bloch bra, projector lattice sums) goes under
  the microscope with the span excuse eliminated.

Either way we need spherical on the lattice path — and (user) we will eventually want SALC for
lattices too.  This plan unblocks both with ONE seam.

## The blocker, named

`L3::GPWFactory` consumes the molecular basis through the **`Molecule::LatticeSum1E`** capability
(analytic periodic S/T/V, ⟨χ|g⟩ projector overlaps, `MakeLocalGaussian` 3C, `CollocateDensity` /
`IntegratePotential` streams) — implemented ONLY by `PG_Cart` / `PG_Cart_MnD`.  `PG_Spherical` is
molecular-only: no periodic capability, no collocation streams.  That is the whole blocker.

## The design: peer implementations behind the abstract capability (USER 2026-08-12)

**`Molecule::LatticeSum1E` stays the abstract interface; `PG_Cart` and `PG_Spherical` are PEER
implementations selected by virtual dispatch** — GPW cross-casts to the capability and never learns
which family answered.  No decorator, no factory plumbing: you get spherical lattice physics by
CONSTRUCTING a spherical basis (`MakeBasisLowQ` grows the choice), exactly the project's
"capabilities on the types that have them" rule.

The congruence identity is then the spherical implementation's PRIVATE strategy, not an
architectural layer: a spherical basis function is a fixed linear combination of Cartesian
components at the same exponent/centre (χ_sph = Σ T χ_cart; T block-diagonal per shell — s,p
identity, each Cartesian-d 6-block → 5 columns; C_l tables already in
`src/BasisSet/Molecule/Evaluators/PG_Spherical_MnD/SolidHarmonics.C`), so `PG_Spherical`'s
implementor COMPOSES the Cartesian lattice-sum engine internally and transforms at its own
boundary:

- matrices (S, T, V, `MakeLocalGaussian`):  M_sph = T† M_cart T
- vectors (⟨χ|g⟩ KB overlaps):              b_sph = T† b_cart      → KB assembles in spherical
  space automatically; the d channel becomes structurally identical to CP2K's
- density path: D_cart = T D_sph T† ONCE per iteration, then the existing collocation streams,
  ladder, T3 fold, LRU run UNCHANGED; the integrate-back KS block returns h_sph = T† h_cart T.

No spherical lattice sums, no new integrals.  PREREQUISITE REFACTOR: the pair-enumeration /
magnitude-screen / stream kernels currently living inside `PG_Cart_MnD`'s evaluator get HOISTED to
a shared internal module (or owned as an engine object) so both implementations call them — the
kernels are primitive-Cartesian-pair machinery either way.

**Lattice SALC later (user):** still fits — SALC is a symmetry-adaptation over WHICHEVER basis
(molecular precedent `SymmetryAdaptedBasisSet.C`); once both families answer `LatticeSum1E`, the
adapted set rides on either.  The Shubnikov ops of SymmetryUpgradePlan S1–S3 supply the magnetic
reps when that day comes.

## The second axis (USER 2026-08-12): the integral ENGINE is its own dispatch layer

The basis families must also be ENGINE-BLIND: `PG_Cart` and `PG_Spherical` may not know whether
qchem-MnD or libCint computes their primitives — a second layer of virtual dispatch beneath the
family (the bridge shape):

    capability faces (LatticeSum1E, Integrals_*)      ← what GPW/Hamiltonian consume
      basis FAMILY:  PG_Cart   | PG_Spherical         ← angular delivery, engine-blind
        engine:      qchem-MnD | libCint              ← primitive integrals, family-agnostic

Today the two axes are FUSED in the evaluator names (`PG_Cart_MnD`, `PG_Spherical_MnD`,
`PG_LibCint` is simultaneously a family and an engine).  Consequences of the split:

- The hoisted lattice kernels (pair enumeration, magnitude screen, ε-summation, stream assembly)
  are exactly the ENGINE-AGNOSTIC middle: the kernel enumerates images and asks the engine for
  FINITE per-offset primitives ⟨χᵢ|·|χⱼ(·−R)⟩.  **libCint then gains periodic support for free**
  — it never learns about lattices (THERE IS NO CUT stays kernel-side, engine-independent).
- The spherical family asks its ENGINE for spherical shell blocks: the MnD engine's default is the
  C_l congruence over its Cartesian primitives; libCint OVERRIDES with its native `sph` API — the
  parked SphericalSALCPlan S3b lands here as an engine method, not a family fork.
- The m-ordering/normalisation conventions (the known real cost) become an ENGINE-boundary
  contract: each engine must deliver blocks in the family's declared ordering.
- **Granularity contract (user): the engine interface trades in BLOCKS of integrals** — shell-pair
  blocks, whole matrices, stream batches — never per-primitive or per-point calls, so the virtual
  dispatch cost stays unmeasurable.

## Watch-outs (from the molecular SALC campaign)

- **m-ordering + normalisation conventions are THE cost** (SphericalSALCPlan S3 lesson): T must be
  expressed in PG_Cart's own component ordering and its normalised-Cartesian convention;
  SolidHarmonics.C already carries the conversion for l=0..3 — reuse, don't re-derive.
- Cache keys: the wrapped basis needs its own `BasisSetID` so `DB_Cache`/stream keys never collide
  with the Cartesian view.
- Dimensions shrink (MnO cell: 122 → 106; 8 d shells × 1 contaminant × 2 Mn) — vet/ortho and the
  EC bookkeeping just see a smaller n, nothing special.
- The T3 stream fold and Shubnikov machinery act on the INNER Cartesian streams — unchanged.

## The conditioning dividend (user)

The Mn valence basis had TRUE s exponents REMOVED to survive cond(S) — and the 8 d-shell
contaminants are near-duplicate diffuse s functions, i.e. they are PART of why the s window had to
be trimmed.  With spherical d the contaminants vanish, cond(S) improves, and the trimmed s
exponents can come back ON PURPOSE (qchem.ValenceBasisGen regen) — the s span restored where we
control it instead of by Cartesian accident.

## Increments

- **I0 DONE (2026-08-12, runs 45/46 = replays of 43/44 with `GPW_NL_PER_L=1`; commit fba8baa6).**
  Δ_NL(AFM−FM) by channel: **l=0 −0.477 / l=1 −0.00001 / l=2 −0.823 Ha** (sum −1.300 ✓; FM crystal
  per channel ≈ atomic superposition: +2.92/+0.02/−32.50 vs +3.03/0/−32.09).  THE REVISION: the
  s-contaminant was exonerated as an attraction SOURCE (h_s>0), but **37% of the weak basin's NL
  reward is an l=0 repulsion DODGE** — energy-lowering, and available ONLY to our span (the extra
  r²-Gaussians let the density rearrange s-character away from the β_s projectors; spherical CP2K
  structurally cannot).  The remaining −0.82 Ha is l=2 with a span-identical operator but a
  span-enabled density.  SHARPENED I2 PREDICTION: spherical d deletes the s-dodge freedom — if
  span is the story, the deep-moment basin wins and the ordering heals.
- **I1 DONE (2026-08-12 evening): the spherical lattice view is LIVE.**
  `MakeSphericalLatticeView(Real_BS) -> Real_BS` (`PG_Spherical/LatticeView.C` + Imp; classes
  unexported — consumers see only Real_BS): the view block answers `Molecule::Orbital_1E_IBS` +
  `LatticeSum1E` in the spherical span, composing the wrapped basis through the ABSTRACT faces only
  (engine-blind by construction; the kernel-hoist is deferred until an engine wants native-spherical
  answers).  T built per shell from `ShellRep::Monomials()` (a NEW soft capability on the rep — the
  `AoShell` "owns no monomials" pin is respected) + `Math::SphericalShell` raw harmonics, each
  column normalised against the inner shell's own overlap block (conventions measured, not
  re-derived).  Knobs: `GPW_SPHERICAL=1` (MakeBasisSR/MakeBasisLowQ) and `GPW_MN_SPHERICAL=1`
  (the Mn box A/B — previously wired to the NATIVE family and DEAD on the missing-LatticeSum1E
  cross-cast: the remembered blocker, confirmed).  GATES (suite 711→713, all green):
  `M_SphericalView.sp_null_test` (T=identity exactness) and `.d_matches_native_family` — the
  ordering-INVARIANT probes (sorted eig(S), pencil spectra eig(S⁻¹T)/eig(S⁻¹V), the span projector
  f(r)ᵀS⁻¹f(r)) against the native PG_Spherical family, since the two families ORDER functions
  differently; plus the crystal-tier A/B: Mn d⁵ box through the view −14.3498 vs the native
  spherical FACADE −14.3452 (**4.6 mHa — the same agreement class as the Cartesian pair's 12 mHa;
  the view is validated through a full GPW SCF with occupied d**).
  **FINDING that reorders the ladder: dropping the contaminants costs the Mn ATOM 281 mHa**
  (−14.626 → −14.345 vs the exact sextet −14.674) — the Cartesian contaminants were covering for
  the TRIMMED true-s window (the user's recollection, confirmed quantitatively).  So **I3 (restore
  the s exponents under spherical d) is a PREREQUISITE for a clean I2 verdict**: with the current
  basis the spherical MnO A/B would be s-starved on Mn, and the s-dodge channel was precisely the
  ordering-relevant one.
- **I2 ARM 1 DONE (2026-08-13, runs 47/48: spherical view, SR-trimmed basis, n=106): ★ THE
  ORDERING HEALED — SPAN CONVICTED, NO CODE BUG.**  AFM det −61.22677 / FM det −61.18122 ⇒ **AFM
  below FM by 45.5 mHa** (Cartesian arm: FM below by 40.1).  The per-l mechanism check: Δ_NL(AFM−FM)
  collapsed −1300 → **+108 mHa** (l=0 dodge −477 → +29: DEAD with no contaminants to dodge through;
  l=2 −823 → +80).  The user's basis-incompleteness hypothesis is CONFIRMED: the reversal was the
  Cartesian-d span differentially stabilising the weak-AFM/FM configurations — not physics, not an
  assembly defect.  NUANCE: the AFM state is STILL weak-moment (m_stag 0.656, |m̃|Ω/2 = 3.1 vs CP2K
  4.65) — the span fix flips the ordering without deepening the moment; the moment question moves to
  the span-matched absolute comparison (arm 2).  Also: cond(S) improved 2.4e3 → 505 on the cell.
  **BASIS DISCOVERY THAT REFRAMES THE RESIDUALS: the CP2K transcription (VALENCE-LOWQ-BASIS)
  predates the SR cell-trims** — the banked oracles were ALWAYS in the untrimmed spherical span (Mn
  7s+8d incl. diffuse d 0.18; O 6s incl. 0.15, 5p incl. 0.18).  So spherical-over-SR sits +273 mHa
  above CP2K FM legitimately (poorer span), and the Cartesian FM "8 mHa agreement" was
  contaminants-vs-diffuse COMPENSATION.  valence_lowq_sph v2 (ace4ec8a) now matches the
  transcription exponent-for-exponent: MnO cell n=136 == CP2K's 2×47+2×21.
- **I2 ARM 2 (in flight, runs 49/50): the CP2K-span A/B** — spherical view over valence_lowq_sph
  v2, AFM + FM.  Every residual vs the oracles (Γ AFM −61.4706 / FM −61.4617) is now an honest
  operator/grid number: span excuses eliminated.
- **I3 (basis):** restore the trimmed Mn s window under spherical d; re-vet cond(S); regen the
  library entry (ValenceBasisGen).
- **I4 (later):** lattice SALC on the same seam — symmetry-adapted T from the space-group
  (grey) or Shubnikov (magnetic) reps; blocks H per irrep at each k.  Not scheduled; the seam is
  its prerequisite and is done at I1.
