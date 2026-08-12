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

- **I0 (cheap, evidence):** per-l/per-species split of E_NL (extend the EenNL diagnostic) —
  verify the sign structure: the AFM-side −1.3 Ha must flow through l=2 (s/p slices small and
  REPULSIVE, per the Mn q7 h matrices: s diag +2.8/+2.5/+2.6, p +1.4/+0.3, d −7.995).
- **I1 (the unblock):** `TransformedLatticeSum1E` + the spherical-view wiring in `GPWFactory`
  (opt-in: `GPWParams::spherical` / test knob `MNO_SPHERICAL=1`); gates: spherical S/T/V ==
  congruence-transformed Cartesian (exact, unit tier); spherical Si Γ total == Cartesian Si Γ
  total to SCF tolerance (s/p-only basis ⇒ T=identity — the null test); Mn d⁵ pseudo-atom in-box
  spherical vs sextet ATOM oracle.
- **I2 (the physics):** MnO Γ A/B — Cartesian vs spherical, FM + AFM, imposed Shubnikov,
  {Ladder,GDM}×{5e-3,0}.  THE span-vs-bug verdict.  Then vs the banked CP2K oracles
  (doc/CP2Kresults.md): same span, same functional, same deck — residuals are now honest
  code-vs-code numbers.
- **I3 (basis):** restore the trimmed Mn s window under spherical d; re-vet cond(S); regen the
  library entry (ValenceBasisGen).
- **I4 (later):** lattice SALC on the same seam — symmetry-adapted T from the space-group
  (grey) or Shubnikov (magnetic) reps; blocks H per irrep at each k.  Not scheduled; the seam is
  its prerequisite and is done at I1.
