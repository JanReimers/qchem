# +U projectors independent of the orbital basis

**Born 2026-10-24, for execution by the next session** (user: *"clean code is MUCH higher priority"* than a cheap route;
*"if we have to redo some runs that is fine"*).  The evidence is in `Records/HubbardUHistory.md` (2026-10-24
entries); the tracker row is `OpenWork.md` §2 "+U projector independent of the orbital basis".  When this
plan is done its text goes to `Records/HubbardUHistory.md` and this file is deleted.

## 0. Status (2026-10-24, end of the first execution session)

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

## 1. Why

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

## 2. The design: the separable-PP pattern, applied to +U

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

## 3. Increments (each one sweep-green and committed before the next)

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

## 4. Risks and open questions (settle in the increment that meets them)

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
