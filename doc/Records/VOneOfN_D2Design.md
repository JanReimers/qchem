# V-ONE-OF-N sub-step 3 (D2), redesigned: the density projects, the Hartree term applies the physics (design note, lane S 2026-11-08)

RECORD, for the user's approval before any hub code.  It supersedes `doc/Records/VOneOfNDesign.md` §5–§6 for D2, and §6's
`ProjectOnto` line (sub-step 2's ruled option 1) is revisited in §2 below.  Read-only from lane-s-h-sprint at d09e98dc1;
nothing built or run for this note.

**The user's ruling, 2026-11-08, relayed by the controller.**
1. The density answers ONE contraction of its own D with a three-index object the caller hands in.  It returns projections,
   index- or G-keyed.  No fit concept and no "Coefficients" appear in its names, and the record is renamed.
2. The periodic projection is ρ̃(G) only, over the overlap tensor.  The 4π/|G|² weight leaves the tensor and is applied by the
   Hartree term, the way the Gaussian route's repulsion lives with Vee and the fitter's metric.
3. Aim: a PW fit basis is an ordinary fit basis whose overlap and repulsion happen to be diagonal.  Inventory what retires and
   what stays PW-specific.
4. Say where numbers can move, and which gates prove 0 moves.
5. Moving `DensitySampler_*` / `FittedCD` out of qcChargeDensity is a separate row, not this step.

## 1. The finding that makes it small: the PW Hartree term already IS `Vee_DFT`, written by hand

The molecular "normal" Hartree term `Vee_DFT` (`src/Hamiltonian/Internal/Imp/Vee_DFT.C`) runs on `GridDensityFitter`
(`src/Fitting/Internal/GridDensityFitter.C`), and its shape is exactly what the user describes:

| step | `Vee_DFT` + `GridDensityFitter` (δ mesh, molecular) | `Vee_Hartree` (PW / GPW, today) |
|---|---|---|
| projection | the density contracts D into the fitter's forward: ρ(r_g) | the density contracts D into `Repulsion3C`: **k·ρ̃(G)**, kernel baked |
| "the fit" | orthogonal basis: the projection IS the fit | orthonormal {G}: the projection IS the fit (`OrthoFunctionFitter` says so) |
| the physics | `CoulombOperator::Potential(q) → v`, OWNED BY THE FITTER | inside the tensor's forward |
| matrix | the same pair's adjoint, Σ_g w_g v_g χ_iχ_j | `ContractAdjoint(Repulsion3C, V_H)` (overlap adjoint) |
| energy | `CoulombOperator::HartreeEnergy(trDV, q, v)` = ½⟨q,v⟩ | ½Ω Σ\|V_H\|²/k, a hand-written Parseval loop in the term |

So the cure is not a new PW design.  The PW route becomes this same row with a diagonal Poisson operator:

- **The density** returns ρ̃(G) through the OVERLAP tensor.  That is today's `GetFourierDensity`, which already exists and is
  already star-averaged.
- **A G-space Coulomb operator** applies k(G) = 4π/|G|², with k(0)=0, in the fitter.  That is exactly where
  `CoulombOperator` puts 1/r for the mesh.
- **The term** is `Vee_DFT`: the fitter's `Repulsion(bs)` for the matrix, and `HartreeEnergy` for the energy.

The kernel then lives in ONE place, the operator, instead of three:
- the PW `Repulsion3CTensor`'s `kernel` column array;
- the GPW collocator's `coulomb` flag;
- `SeedCD` / `FourierMixCD`'s hand-applied `k*rt`.

The Dunlap path also runs generically on an orthonormal {G} basis.  There G = diag(k), and c = G⁻¹⟨ρ|c⟩_C = ρ̃ for G≠0.  But
G=0 is singular, and Dunlap's charge constraint would be what fixes c(0).  For an orthonormal basis the first-order `Vee_DFT`
energy is exact, as the 2026-11-06 ruling for the mesh says, so **I recommend `Vee_DFT`, not Dunlap, for PW**, and leaving
Dunlap Gaussian-only.  Both give the same V_H; only the G=0 bookkeeping differs.

## 2. The density's one operation, and the record's name

```
// on tChargeDensity<T> (ChargeDensity hub) -- replaces ProjectOnto(ScalarProjector), GetRepulsion3C(DensityProjector),
// GetFourierDensity, GetRepulsion3C(fit), GetRepulsion3C_Raw, StarAverage, GetRhoOnGrid
virtual Projections Project(const BasisSet::PairProjector<T>& p) const=0;   // <c|rho> for every function c of p's basis
```

**`PairProjector<T>`** is a BasisSet face, not a Fitting one.  It is "the three-index object the caller hands in", per orbital
block.  It merges today's `Fitting::ScalarProjector` and `Fitting::DensityProjector`, which differ only in the scalar axes they
are templated on.

```
template <class TFit> class PairProjector               // qcBasisSet; implemented by the FITTERS (they build it from their fit basis)
{
    virtual const qcMesh::MatrixForward<double>& Forward(const Orbital_DFT_IBS<double,TFit>&) const=0;   // a real block
    virtual const qcMesh::MatrixForward<dcmplx>& Forward(const Orbital_DFT_IBS<dcmplx,TFit>&) const=0;   // a complex block
    virtual Projections Project(const ScalarFunction<double>& f) const=0;   // a matrix-free density projects its FIELD
};
```

- **A DM density** (`IrrepCD`, the composite) contracts its D into `p.Forward(block)` and sums blocks with `+=`.  The periodic
  composite then star-averages, as `GetFourierDensity` does today.  The only realization-specific line is that one star-average
  call, and the composite already owns the crystal point ops.
- **A matrix-free density** (`NumericCD`, `SeedCD`, `FourierMixCD`) has no D.  It answers `p.Project(*this)`, the projection face
  every fit basis already has (`OrthogonalFit`, `Overlap(ScalarFunction)`).  ⚠ **This is where numbers WOULD move**: today
  `SeedCD` returns an ANALYTIC structure-factor ρ̃ and `FourierMixCD` returns its stored map.  Sampling them and FFTing would
  change bits at the iteration-0 seed, and on every mixed density.  So `PairProjector` needs a second matrix-free entry that
  keeps today's arithmetic:
  - **Recommended:** `Projections Project(const G_StructureFactor& seedSource)` for the seed, and an identity answer
    `Project(const Projections& mine, BasisSetID)` for a density that IS a projection record on the same basis (the mixed
    density).  A mismatched basis THROWS.  It is a basis-identity check, not a realization check.
  - **The alternative**, accepting the sampled route and re-judging the anchors, is pin 10 work, not a refactor.  I would not
    take it here.

**The record's name.** It carries ⟨c|ρ⟩ from a forward and a field's expansion into an adjoint.  For an orthogonal basis (δ, PW)
those coincide.  For a Gaussian one the fitter's solve sits between them.  So the name must not claim either role:

- **`Projections`**: what a forward returns.  It is honest for the density's side, and slightly wrong for the adjoint's input.
- **`AuxComponents`**: "components on an auxiliary basis's index axis", neutral on primal vs dual.  ✔ my pick.
- `FieldComponents`: suggests real-space values.

I write `Projections` above for readability; the user picks.  The rename is mechanical: one module (`src/Math/FitCoefficients.C`)
and about 25 sites.

**Sub-step 2 (committed, unlanded) is reworked by this, not discarded.**  Its record, `GSpacePairProjectionIntegrator`, the
integrators' new signatures and every paired `.Vector()` / `.GMap()` read stand.  What changes is the rename, plus
`ProjectOnto(ScalarProjector) → Project(PairProjector)` on the hub.  So the ruling on sub-step 2's option 1 carries over, with the
new names.

## 3. Inventory: what retires, what stays PW-specific

**Retires.**

| what | where | why it can go |
|---|---|---|
| `FourierDensity` face, `FourierDensityBase<T>`, `NoFourierDensity`, the `qchem.ChargeDensity.FourierDensity` module | `src/ChargeDensity/FourierDensity.C` | `Project(PairProjector)` answers all of it |
| the 11 casts to it | `PWTerms_Hartree.C`, `DensitySampler_Pair.C`, `DensityMixer.C` ×4, `PolarizedDensityMixer.C` ×2, `FieldMixer.C`, `CompositeCD.C`, `IntegrationTests/PW/Boxes.C` | they ask the one hub method |
| `GetRepulsion3C(cFIT_CD_ABS)` on every density (IrrepCD_Fourier, Composite, SeedCD, PolarizedSeedCD, FourierMixCD, PolarizedMixCD) | ChargeDensity | the kernel moves to the operator (§1) |
| `GetRepulsion3C_Raw` and `StarAverage` on the face | `CompositeCD.C` | the average happens once, inside the composite's `Project` (it already does for ρ̃) |
| `GetRepulsion3C(DensityProjector)` + `CoulombMetric_ProjectedDensity`, `DensityProjector`, `ScalarProjector` | Fitting hub | merged into `PairProjector` (BasisSet) |
| `OverlapMetric_ProjectedDensity::GetUnconstrainedFit` (the S⁻¹ solve inside `NumericCD`) | Fitting hub + `NumericCD` | the seed projects (`p.Project(*this)`), and `ConstrainedFF` does its own solve |
| `ProjectedDensity_G`, `OrthoFunctionFitter`'s cast to it | Fitting hub | the fitter receives `Projections` |
| `OrthoFunctionFitter::Repulsion`, a throwing "obsolete" stub | Fitting | it becomes the real Hartree matrix (§1) |
| `Vee_Hartree`'s ground-state half: field memo, Parseval energy loop, `Volume()` | `PWTerms_Hartree.C` | it becomes `Vee_DFT` (see the next table for the response half) |
| the PW tensor's `kernel` array; GPW's `MakeCollocator(coulomb=true)` and `Repulsion3CTensor` | `PairProjection.C`, `PW Evaluator`, `GPW_Evaluator.C` | nothing asks for a kernel-baked forward any more |

**Stays PW-specific, because it genuinely is.**

| what | why |
|---|---|
| the diagonal Poisson operator (`ReciprocalLattice::CoulombKernel`, k(0)=0) | the periodic 1/r: a `CoulombOperator` realization beside `DirectSumCoulomb` / `RadialPoissonCoulomb` |
| the IBZ star-average (`SymmetrizeGMap`) inside the periodic composite's `Project` | crystal point ops belong to the k-composite |
| the response's δV_H(G+q) (`TransitionRepulsion`, `CoulombKernel(m,q)`, `ΔGq_Map`) | a Bloch-q field.  It stays on `Vee_Hartree`'s response half until the response route is redone.  ⚠ This means `Vee_Hartree` does not fully retire in this step: its response methods move to a small `Vee_Hartree_Response`, or stay with a ground-state body that delegates |
| `GetRhoOnGrid` / raster ρ_DM(r) | sub-step 4 (R1.0n/o), unchanged here |
| `OrthoNormalScalarFitter` (the XC raster fit, wrap-around lookup) | not a Hartree object; V-XC-ONE-TERM's |

**What has to become scalar-generic for this.**  `Vee_DFT` is a molecular term (`rsmat_t MakeMatrix`).  The PW term has both
block corners (`cobs_t` and `robs_t`).  So the merged term is `tVee_DFT<T>`, or `Vee_DFT` gains the two periodic `MakeMatrix`
overloads.  `GridDensityFitter` / `CoulombOperator` carry `rvec_t`, so they take the record (§2) and the operator gains a
G-realization.  Energy = ½⟨q,v⟩_w is a weighted inner product on the record, **the same dispatcher the mixers need in sub-step 4**,
so it is added once, here.

## 4. Where numbers can move, and the gates that say they did not

Every item below is a reordering of the same arithmetic, so it should be roundoff (≲1e-15 relative).  The gates show whether it is.

| # | change | today | after | expected |
|---|---|---|---|---|
| 1 | the kernel vs the k-sum | Σ_k (k·ρ̃_k) per block, then merged | k·(Σ_k ρ̃_k) | roundoff |
| 2 | the kernel vs the star-average | avg(k·x) over a star | k·avg(x).  ⚠ k is computed per member from G=B·m, so it is not bit-constant on a star | roundoff |
| 3 | G=0 | k(0)=0 drops it | same `CoulombKernel`, and ρ̃(0)=N/Ω is carried and multiplied by 0 | exact |
| 4 | the energy | ½Ω Σ \|kρ̃\|²/k | ½Ω Σ k\|ρ̃\|² (or ½⟨ρ̃, kρ̃⟩) | roundoff; it was already "not bit-identical" with the trace form (`eh_trace`) |
| 5 | GPW forward | a Coulomb closure: collocate, FFT, combine, ×k | the overlap closure (collocate / memo replay, FFT, combine), then ×k in the operator | ×k on the same combined value is the same IEEE op; the CollocMemo replays one closure instead of two |
| 6 | seeds / mixed density | `SeedCD`/`FourierMixCD` hand-apply k·rt | the operator applies k to their unchanged ρ̃ | exact if §2's identity and structure-factor entries are taken; **moves** if they are sampled instead |
| 7 | `Vee_DFT` charge diagnostic | absent on PW | `HartreeChargeError` = Ω ρ̃(0) − N, a new diagnostic line | a report change, not a number change |

**Gates, in order.**
1. **Before anything is deleted:** a unit gate, in `UTLattice_BS`, the PW + GPW toys, checking that k·`Contract(Overlap3C, D)`
   equals `Contract(Repulsion3C, D)` to 1e-15 relative, per G.  It proves #1, #2 and #5 on their own.
2. `eh_trace` A/B on the PW and GPW anchors: G-space energy vs Tr(D V_H), the pin-32 shape check.  Before vs after, side by side.
3. The short tier, the long tier, and `-L "panel|periodic"` (the mixers move in this step), each through
   `scripts/anchordiff --tol 1e-10`.  The PW/GPW anchors that carry the claim are listed below; the `WILL_FAIL` panel cases are
   listed apart, each with its failure reason checked.
   - `GPW_NaF.Γ_Imp_eqExactResume` / `eqColdStart`
   - the MnO AFM-II anchor
   - PW Si, NaF and CsI
   - `PW/Boxes`, including the polarized seeds
   - `GPW.RawXCConsistencyFD`
   - the TiO₂ LR CP2K/QE gates, which are the response half, unchanged
4. RelChecked, vet-decks, and the hub audit, regenerated citing the ruling.
5. **One timed MnO run**, before and after, quoting the buckets that change:
   - "scf: V_H per-block ΔG_Map merge"
   - "V_H IBZ star-average"
   - "rho ball closure"
   - "E_H G-space pairing"

   The Coulomb closure call goes away, and the merge and the average now ride ρ̃'s.  So I expect a small saving, but that is a
   claim to measure, not to assume.

## 5. Order, if approved

- **3a.** The rename of the record; `PairProjector` in BasisSet; `tChargeDensity::Project`; every density implements it.  The old
  faces stay alive one commit, forwarding to it.  0 moves.
- **3b.** A G-space `CoulombOperator`; `GridDensityFitter` on the record; `Vee_DFT` scalar-generic.  The PW / GPW Hartree term
  becomes `Vee_DFT` (with gate 1 first), and `Vee_Hartree` keeps only the response half.  0 moves, or roundoff, judged by the gates.
- **3c.** Delete everything §3 retires, the `FourierDensity` face and the 11 casts with it.  Regenerate the hub snapshot.
- Then sub-step 4 (R1.0n/o) as planned.  It is smaller now, since the mixers already ask `Project`.

Hubs touched: ChargeDensity (`Project`; the `FourierDensity` face deleted), Fitting.FunctionFitter (`ScalarProjector`,
`DensityProjector`, the two metric faces and `ProjectedDensity_G` leave), BasisSet (`PairProjector` added; `Orbital_DFT_IBS`'s
periodic `Repulsion3C` loses its last periodic client, while the molecular one stays: it is the real Coulomb ⟨ab|c⟩).

## Questions for the user

1. The PW Hartree term becomes `Vee_DFT` on a diagonal Poisson `CoulombOperator` (first order, §1), not the Dunlap path.  Agreed?
2. Matrix-free densities keep today's arithmetic through two `PairProjector` entries (a structure-factor seed; an identity for a
   projection record on the same basis, which throws on another basis) rather than being sampled (§2).  Agreed?
3. The record's name: `AuxComponents` (my pick), `Projections`, or another?
4. `Vee_Hartree` keeps only the response half, δV_H(G+q), until the response route is redone.  Agreed?
