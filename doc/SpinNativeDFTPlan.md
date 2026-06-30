# Spin-Native DFT Plan — making LDA correlation + occupation spin-first

OpenWork **item B**. Companion to `doc/FacadeDFTPlan.md` (this is its deferred "D2", reframed). Living
document — *plan first, code second*.

## The reframe (the tenet, not "add the polarized case")

Per `feedback_spin_polarized_primary`: **spin-polarized (ρ↑,ρ↓ / n↑,n↓) is the native formulation;
unpolarized is the ζ=0 efficiency collapse.** So every piece below is built spin-first with unpolarized
falling out as the collapse — *not* a `_P` class bolted beside a `_U` original.

The codebase is already partway there: HF and Slater-Xα have first-class `_U`/`_P` pairs
(`Ham_HF_U`/`Ham_HF_P`, `Ham_DFT_U`/`Ham_DFT_P`, `FittedVxc`/`FittedVxcPol`). The bias lives entirely in
the **newest hand-rolled LDA-correlation pieces**: the paramagnetic-only `VWN_Correlation`, the
solitary `Ham_DFTcorr_U` (no `_P`), and `Molecule_EC`'s closed-shell-only aufbau.

## The crux: correlation does NOT separate by spin channel (exchange does)

This is the one genuinely new architectural fact, and it dictates the whole design. **Slater exchange is
channel-separable** — v_x^σ = −3α(3ρ_σ/4π)^{1/3} depends *only* on ρ_σ. That is exactly why
`FittedVxcPol` ([Terms.C:233](src/Hamiltonian/Internal/Terms.C), impl
[FittedVxcPol.C:44](src/Hamiltonian/Internal/Imp/FittedVxcPol.C)) works by holding **two independent**
`FittedVxc`, each fitting the *same* functional against a *single-channel* density (`ucd`, then `dcd`):

```cpp
rsmat_t Kab = s==Spin::Up ? itsUpVxc->GetMatrix(bs,s,ucd)   // each sees ONE channel
                          : itsDownVxc->GetMatrix(bs,s,dcd);
```

**Correlation is not separable.** v_c^σ(ρ↑,ρ↓) couples *both* channels through r_s = r_s(ρ↑+ρ↓) and the
polarization ζ = (ρ↑−ρ↓)/(ρ↑+ρ↓). So:

- The scalar `ExFunctional::GetVxc(double rho)` interface ([ExchangeFunctional.C:22](src/Hamiltonian/Internal/ExchangeFunctional.C))
  is insufficient — the functional must see **both** ρ↑(r) and ρ↓(r) at each evaluation point, plus the
  target spin σ.
- The `itsUpVxc->GetMatrix(bs,s,ucd)` pattern (feeding one channel) **breaks** — a `FittedVcorrPol`
  must hand the functional the *full* `Polarized_CD`, not a single `DM_CD`.

This is the load-bearing decision; everything in Piece 1 + Piece 2 follows from it.

---

## Piece 1 — spin-native VWN5 correlation (the physics core)  ·  ✅ DONE (B1)

**Landed** in [VWN_Correlation.C](src/Hamiltonian/Internal/VWN_Correlation.C): the para/ferro/spin-stiffness
branches share one `Gval(x,Params)`/`dGdx` closed form; `EvalRZ(rho,ζ)` returns ε_c and its (r_s,ζ)
partials; the public two-channel face is `GetEpsC(rUp,rDn)` / `GetVc(rUp,rDn,Spin)` (face B). The scalar
`GetVxc/GetEpsXc` face is the ζ=0 collapse, kept byte-identical (the SCF anchors are unmoved). Validated
in [LDA_XC_UT.C](UnitTests/LDA_XC_UT.C) vs libxc `LDA_C_VWN` polarized over an (r_s,ζ) grid to 1e-9
(`VWN5PolarizedMatchesLibxc`) + a `SpinNativeCollapsesToScalarFace` consistency check. 150/150 UTMain green.
(Gotcha caught in review: the f(ζ) denominator is 2^{4/3}−2 = 0.51984, *not* 2^{2/3}.)

The original scoping (kept for reference):

**Today:** [VWN_Correlation.C](src/Hamiltonian/Internal/VWN_Correlation.C) implements *only* the
paramagnetic ε_c^P(r_s), v_c^P(r_s) (one parameter set A,b,c,x0). It is `ExFunctional`-derived with the
scalar `GetVxc(double rho)` / `GetEpsXc(double rho)` face. Validated pointwise vs libxc `LDA_C_VWN`.

**Target:** the full VWN5 spin interpolation. ε_c(r_s,ζ) interpolates paramagnetic (ζ=0) ↔ ferromagnetic
(ζ=1) via the spin-stiffness term:

```
ε_c(r_s,ζ) = ε_c^P(r_s)
           + α_c(r_s) · [f(ζ)/f''(0)] · (1−ζ⁴)
           + [ε_c^F(r_s) − ε_c^P(r_s)] · f(ζ) · ζ⁴

f(ζ)   = [(1+ζ)^{4/3} + (1−ζ)^{4/3} − 2] / (2^{4/3} − 2)
f''(0) = 4 / (9(2^{1/3} − 1))
```

Each of ε_c^P, ε_c^F, α_c has the **same closed form** already in the file (the log + atan in X(x)=x²+bx+c,
x=√r_s) with a *different* parameter quadruple. We add two more parameter sets:

| Branch | A | b | c | x0 | note |
|---|---|---|---|---|---|
| paramagnetic ε_c^P | 0.0310907 | 3.72744 | 12.9352 | −0.10498 | ✅ present today |
| ferromagnetic ε_c^F | 0.01554535 | 7.06042 | 18.0578 | −0.32500 | new (½ of the para A) |
| spin stiffness α_c | −1/(6π²) | 1.13107 | 13.0045 | −0.0047584 | new (note negative A) |

The potentials are the two channel derivatives v_c^σ = ε_c + ρ ∂ε_c/∂ρ_σ, which expand via
∂ε_c/∂r_s and ∂ε_c/∂ζ (with ∂r_s/∂ρ_σ and ∂ζ/∂ρ_σ). The unpolarized collapse (ζ=0 ⇒ α_c, ε_c^F terms
vanish, v_c^↑=v_c^↓=v_c^P) must reproduce the *current* `Vc(rho)` **bit-identically** — that is the
regression guard (existing `LDA_XC_UT` pointwise-vs-libxc test stays green for ζ=0).

**Interface decision (Piece 1 ↔ crux).** Two options; recommend **B**:

- **A. Spin-tagged scalar (mirror `SlaterExchange(double,Spin)`).** Add a `Spin` ctor + a polarization
  channel the functional reads. *Rejected:* exchange gets away with this because it's separable; a scalar
  `GetVxc(double ρ_σ)` literally cannot compute v_c^σ — it never sees the *other* channel. Forcing it
  would mean smuggling ρ↓ through `InsertChargeDensity`, which is exactly the non-local coupling we should
  make explicit, not hide.
- **B. (recommended) A small polarized-functional face.** Give the correlation functional a two-argument
  evaluator — `GetVc(double rhoUp, double rhoDown, Spin σ)` and `GetEpsC(double rhoUp,double rhoDown)` —
  and have its `operator()(r)`/the fit read **both** channel densities at r. The unpolarized path calls it
  with ρ↑=ρ↓=ρ/2 (or routes through the existing scalar collapse). This makes the coupling a typed part of
  the interface, honouring `feedback_compile_time_over_runtime` (the capability lives on the type that has
  the two channels). Keep `VWN_Correlation` (scalar) as the ζ=0 façade if convenient, or let the new class
  subsume it via the collapse — decide when coding, guided by keeping `LDA_XC_UT` bit-identical.

Validation: extend `LDA_XC_UT` to check v_c^↑, v_c^↓, ε_c against libxc `LDA_C_VWN` (it is spin-native:
`xc_lda_vxc` with `XC_POLARIZED`) across a (r_s, ζ) grid. This is the piece's correctness anchor.

## Piece 2 — `Ham_DFTcorr_P` + `FittedVcorrPol` (carry both channels; U = collapse)  ·  ✅ DONE (B2)

**Landed:** abstract `SpinCorrelation` face (two-channel `GetEpsC`/`GetVc`) in the ExFunctional module;
`VWN_Correlation` implements it. New `FittedVcorrPol` ([Imp/FittedVcorrPol.C](src/Hamiltonian/Internal/Imp/FittedVcorrPol.C))
fits v_c^σ(ρ↑,ρ↓) against the **full** `Polarized_CD` (not the per-channel split exchange uses), with the
SAD-seed fallback (ρ↑=ρ↓=ρ/2 ⇒ v_c^P). Energy via a `FittedEpsCPol` Dynamic_CC: fit ε_c(ρ↑,ρ↓) once, the
polarized `DM_Contract` sums both channels ⇒ ∫ε_c·ρ. `Ham_DFTcorr_P` mirrors `Ham_DFTcorr_U`
(FittedVee + FittedVxcPol Dirac + FittedVcorrPol VWN5). Both Factory throws un-gated (`Model::LDA` and
`XC::DiracVWN`). Anchor `M_DFT.WaterPolarizedLDA`: closed-shell water collapses to the unpolarized LDA
anchor (-75.9324615507) to 1e-6. 151/151 UTMain green; `_U` path byte-identical. (LibXC-polarized still
throws — a separate increment; see the Factory feedback.)

The original scoping (kept for reference):

**Today:** only `Ham_DFTcorr_U` exists ([Hamiltonians.C:51](src/Hamiltonian/Internal/Hamiltonians.C),
impl [Imp/Hamiltonians.C:52](src/Hamiltonian/Internal/Imp/Hamiltonians.C)): Dirac exchange via
`FittedVxc(SlaterExchange(2/3))` + VWN5 correlation via `FittedVcorr(VWN_Correlation)`, sharing one Vxc
fit basis. The exchange half already has its polarized sibling (`FittedVxcPol`); the **correlation half
does not** — there is no `FittedVcorrPol`.

**Build:**
1. **`FittedVcorrPol`** — the polarized correlation term. Unlike `FittedVxcPol` it **cannot** delegate to
   two single-channel `FittedVxc`; it must fit the spin-native functional (Piece 1, face B) against the
   *full* `Polarized_CD` so the functional sees both channels at each mesh point. Reuse `FittedVcorr`'s
   "override `GetEnergy` to use ε_c not the ¾ exchange virial" trick ([Terms.C:283](src/Hamiltonian/Internal/Terms.C)),
   polarized: E_c = ∫ ε_c(ρ↑,ρ↓)(ρ↑+ρ↓).
2. **`Ham_DFTcorr_P`** ctor — mirror `Ham_DFT_P` ([Imp/Hamiltonians.C:164](src/Hamiltonian/Internal/Imp/Hamiltonians.C)):
   `FittedVee` (Coulomb) + `FittedVxcPol(SlaterExchange)` (exchange, already exists) + the new
   `FittedVcorrPol` (correlation). The `_U` stays as the ζ-collapse and should remain byte-identical.
3. **Seed robustness** — the SAD/spin-agnostic seed feeds an unpolarized total ρ before the first
   `Polarized_CD` exists. `FittedVxcPol::CalcMatrix` already handles this (the `cd85d13c` fix — fall back
   to the unpolarized branch when the `dynamic_cast<Polarized_CD*>` is null). `FittedVcorrPol` must do the
   same (ρ↑=ρ↓=ρ/2 on the seed), or the polarized-LDA + SAD path re-introduces the null-deref that bug
   fixed.

**Factory wiring:** delete the two "polarized LDA not yet wired" throws
([Imp/Factory.C](src/Hamiltonian/Imp/Factory.C): the `Model::LDA` branch ~L124 and the
`XC::DiracVWN`/`XC::LibXC` branches ~L141) and return `new Ham_DFTcorr_P(st,mp,bs)` for
`Pol::Polarized`. The `Model`/`XCFunctional` enums are unchanged — purely a new resolver branch.

## Piece 3 — open-shell molecular occupation (n↑, n↓)

**Today:** [Molecule_EC](src/ElectronConfigurations/Molecule_EC.C) is `Molecule_EC(int Ne)` only;
`GetN(Irrep)` does closed-shell aufbau — odd Ne splits (Ne±1)/2 across spins, even Ne splits Ne/2. There
is no way to express, say, triplet O₂ (n↑−n↓ = 2). Atoms have `Atom_EC::AssignUnpaired`
([Imp/Atom_EC.C:80](src/ElectronConfigurations/Imp/Atom_EC.C)); molecules have **no equivalent**.

**Build:** extend `Molecule_EC` to carry an explicit channel split — the spin-first formulation is
`(nUp, nDown)`; the closed-shell `Molecule_EC(Ne)` becomes the collapse nUp=nDown=Ne/2 (keep it as a
delegating ctor so existing call sites are untouched). `GetN(Irrep)` returns the stored per-spin count;
aufbau within each spin channel is already handled downstream
(`tCompositeWF::FillOrbitalsAufbau` fills each spin channel independently against `GetN(Irrep(spin,·))`,
[Imp/CompositeWF.C:178](src/WaveFunction/Internal/Imp/CompositeWF.C)) — so **no SCF-loop change is needed**;
the EC just has to report asymmetric per-spin totals. Front-end input is "multiplicity 2S+1" → 2S = nUp−nDown,
nUp+nDown = Ne (Piece 4 owns the conversion).

Spin-restricted-open vs unrestricted is a later refinement; v1 target is **unrestricted open-shell**
(matches the polarized Hamiltonian path), with the closed-shell collapse bit-identical to today.

## Piece 4 — facade multiplicity on `CalcOptions`

**Today:** neither `CalcOptions` ([Calculation.C:44](src/Calculation/Calculation.C)) nor
`AtomCalcOptions` ([AtomCalculation.C:47](src/Calculation/AtomCalculation.C)) has a multiplicity/spin
field. `Pol` selects only whether the *Hamiltonian* is spin-unrestricted; it does **not** set occupation
(today closed-shell Ne is implied by the structure).

**Build:** add `int multiplicity = 1` (= 2S+1) to `CalcOptions`. The facade converts it to (nUp,nDown)
for `Molecule_EC` (Piece 3): nUp−nDown = multiplicity−1, nUp+nDown = Ne (validate same parity; throw a
clear message otherwise — `feedback_compile_time_over_runtime` says fail loud, not assert-crash). A
multiplicity > 1 should auto-imply `Pol::Polarized` (an unrestricted-open-shell calc), or throw if the
user explicitly set `Pol::UnPolarized` with multiplicity > 1 (contradiction). `AtomCalculation` already
reaches open-shell via `Atom_EC`'s Hund's-rule machinery, so atom-side multiplicity is optional/secondary.

---

## Staging & test anchors

| Stage | What | Anchor |
|---|---|---|
| **B1** ✅ | Spin-native VWN5 (Piece 1, face B). Keep ζ=0 bit-identical. | `LDA_XC_UT` extended: v_c^↑,v_c^↓,ε_c vs libxc `LDA_C_VWN` polarized on an (r_s,ζ) grid; ζ=0 byte-identical to current `Vc` |
| **B2** ✅ | `FittedVcorrPol` + `Ham_DFTcorr_P` + Factory un-gate (Piece 2). Seed fallback. | `M_DFT.WaterPolarizedLDA`: closed-shell water collapses to the unpolarized LDA anchor (-75.9324615507) to 1e-6; `_U` path byte-identical |
| **B3** | Open-shell `Molecule_EC(nUp,nDown)` + closed-shell collapse (Piece 3). | a known open-shell molecule (e.g. triplet O₂ or the doublet from an odd-Ne radical) — "did E move" sentinel; closed-shell byte-identical |
| **B4** | `CalcOptions.multiplicity` + facade conversion + Pol auto-imply (Piece 4). | facade-level open-shell test reaching the same B3 energy through the public front door |

**Order is forced:** B1 (functional) → B2 (Hamiltonian term consuming it) → B3 (occupation) → B4 (facade).
B1+B2 are the polarized *machinery*; B3+B4 are the open-shell *occupation* that exercises it. A closed-shell
polarized run (ζ=0, multiplicity 1) is the bridge — it touches B1+B2 only and must collapse to the
unpolarized number exactly, which is the cheapest end-to-end check before open-shell occupation lands.

## Sequencing toward the north-star

`project_battery_voltage_goal` wants magnetism (spin texture is the physics) and ultimately **PBE+U**.
This plan delivers spin-native **LDA**; the rungs above it are separate library increments and are *not*
in scope here:

- **PBE / GGA** — density-gradient machinery on the mesh; the highest-value functional for batteries but a
  real increment (`doc/FacadeDFTPlan.md` flags it). The spin-native correlation interface (face B) built
  here is the right shape to extend to a spin-native GGA correlation later.
- **+U** — orbital-dependent Hubbard term; downstream of forces and multi-species PP on the battery roadmap.

So B's payoff: it removes the spin bias from the newest LDA pieces and makes the *interface* spin-native,
so PBE and +U inherit a spin-first foundation instead of a closed-shell one to retrofit.

## Open questions (resolve while coding, not before)

- **Functional face B exact signature** — `GetVc(rhoUp,rhoDown,Spin)` vs a tiny `(rho,zeta)` struct. Pick
  whichever keeps the ζ=0 collapse cleanest against `LDA_XC_UT`.
- **Does the Vxc fit basis suffice for the coupled correlation potential?** Exchange and correlation share
  one `CreateVxcFitBasisSet` today; the coupled v_c^σ is smooth, so the same basis *should* fit — verify
  against libxc-on-the-mesh, don't assume.
- **Keep `VWN_Correlation` (scalar) or subsume it?** Leaning *keep as the ζ=0 façade* to protect the
  existing test, but the spin-native class can provide the collapse itself. Decide at B1.
- **SROHF later** — v1 is unrestricted-open-shell; spin-restricted open-shell is a deliberate non-goal.
