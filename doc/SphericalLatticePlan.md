# Spherical d on the lattice path — the fixed-transform decorator plan

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
