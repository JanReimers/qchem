# Angular Math consolidation — `qchem.Math.Angular`

**Status: PARKED (spec only).** Design outline for pulling the shared angular/polynomial basis-math out of
the basis evaluators and qcSymmetry into one foundational qcMath module. Not started; recorded so it isn't
re-derived. Prompted by the `IVec3` vs `Polarization` type duality surfaced during the Spherical SALC work
(`doc/SphericalSALCPlan.md`).

## The problem: the same angular math lives in three places

Two "general angular math" concepts are currently defined inside the **basis evaluators** and partially
re-declared in **qcSymmetry** and its tests:

1. **Cartesian monomial exponents** `(nx,ny,nz)`.
   - qcSymmetry: `IVec3 = std::array<int,3>` (`src/Symmetry/CartesianRep.C`) — used as a `std::map` key and
     axis-indexed in the rep builders.
   - Evaluator: `Polarization` (`src/BasisSet/Molecule/Evaluators/PG_Cart_MnD/Polarization.C`) — the same
     `(n,l,m)` triple + ordering/arithmetic, plus evaluator-only methods.
   - The spherical extractor converts one to the other per component (`IVec3{p.n,p.l,p.m}`).

2. **Real solid harmonics as Cartesian expansions** (the c2s map).
   - Evaluator: `SolidHarmonics` / `SphericalShell(int l)` → `vector<vector<CartTerm>>`, with
     `CartTerm = {Polarization p; double c;}` (`.../PG_Spherical_MnD/SolidHarmonics.C`) — the s,p,d,f
     definitions used to build the basis.
   - qcSymmetry: `HarmonicC2S = vector<vector<pair<IVec3,double>>>` (`src/Symmetry/SphericalRep.C`) — the
     **same shape**; qcSymmetry holds the *type*, the evaluator holds the *data*.
   - The qcSymmetry tests (`M_SphericalRep.C`, `M_CartesianRep.C`) hardcode the d-shell harmonics.

So `HarmonicC2S` and `vector<vector<CartTerm>>` are literally the same structure once `IVec3 == Monomial`.

## Why qcMath (not qcSymmetry)

The **layering is decisive**: qcMath is *already* a shared dependency of both qcSymmetry
(`target_link_libraries(qcSymmetry PUBLIC qcMath qcCommon)`) and the basis evaluators. Moving these shared
types into qcMath therefore adds **zero new library edges**.

Putting them in qcSymmetry (the tempting "it's symmetry-ish" choice) would instead force the low-level
Cartesian/spherical **integral evaluators to depend on qcSymmetry** — integrals depending on group theory,
the wrong direction (symmetry sits *above* the basis; the extractor is a deliberate one-way bridge up).
`PG_Cart_MnD` imports no `qchem.Symmetry.*` today, and it should stay that way.

Angular math (monomials, solid harmonics, and later Wigner/Clebsch-Gordan) is genuinely *foundational* — it
belongs below both.

## Proposed module: `qchem.Math.Angular`

One qcMath module holding:

- **`Monomial`** — the `(n,l,m)` Cartesian exponents. Aggregate (so `Monomial{2,0,0}` works), plus the pure
  index operations both consumers need:
  - `int operator[](int)` (read) + `int& operator[](int)` (write) — for the rep builders' `e[i]` / `te[j]+=1`.
  - `operator<` — **lexicographic** (n, then l, then m). This is IDENTICAL to `Polarization`'s current
    `n*LMax^2+l*LMax+m` ordering for valid exponents (n,l,m ∈ [0,64)), so `std::map<Polarization,…>` iteration
    order is preserved (the one live map is `Hermite.C:87` `indexCache`). No `LMax` needed.
  - `operator==` / `operator!=`.
- **`CartTerm`** = `{Monomial p; double c;}` — a monomial with a coefficient.
- **the c2s type** — `HarmonicC2S = vector<vector<CartTerm>>` (unifies qcSymmetry's `HarmonicC2S` and the
  evaluator's `vector<vector<CartTerm>>`).
- **`SphericalShell(int l)`** — the real solid harmonic definitions (s,p,d,f today; l≥4 asserts, as now).
- **Future:** Wigner 3j/6j, Clebsch-Gordan, Ω/Yl helpers that are currently scattered — the reason for the
  general name `Angular` rather than `SolidHarmonics`.

## What stays where (the split of `Polarization`)

`Polarization` is NOT moved wholesale — it keeps its **evaluator-coupled** parts and derives from `Monomial`:

```
class Polarization : public qchem::Math::Monomial   // n,l,m + operator[]/< /== inherited
{
    // KEEP here (evaluator/basis-coupled):
    operator()(rvec3_t) / Gradient(rvec3_t);         // monomial evaluated in real space (uses IntPow)
    operator Index3() const;                          // Cartesian -> MnD Hermite index seam (Index3 is a
                                                      //   BasisSet-evaluator type; CANNOT go to qcMath)
    Polarization operator+/-(const Polarization&);    // MUST return Polarization, not Monomial, so the
                                                      //   rnlm(pa+pb) -> Index3 chain still works
                                                      //   (GaussianRF.C:370/398/434). Do NOT move to base.
    GetTotalL / GetSign / GetMaximumL / operator>;
};
```

The gotcha to respect: if `operator+`/`operator-` moved to `Monomial` they'd return `Monomial`, and
`const Polarization Pab = pa+pb;` would break (no `Monomial→Polarization` conversion). So arithmetic stays on
`Polarization`. Only the pure index skeleton (`n,l,m`, `operator[]`, `operator<`, `==`) lives in `Monomial`.

## Migration — two full-suite-verified stages

**Stage 1 — `Monomial` → qcMath (the clear win).**
- Add `Monomial` to `qchem.Math.Angular`.
- qcSymmetry: `IVec3` becomes `Monomial` (rep builders, `HarmonicC2S`, extractors, tests).
- Evaluator: `Polarization : Math::Monomial`, dropping its duplicated `n,l,m` / `operator<` (inherit) while
  keeping arithmetic + evaluator methods; drop the now-unused `LMax` (used only by the old `operator<`).
- Extractors: the `IVec3{p.n,p.l,p.m}` conversions vanish (a `Polarization` *is* a `Monomial`).
- Verify: **full** UTMain (the integral core is a `Polarization` consumer — energies must be byte-identical),
  plus UTSymmetry.

**Stage 2 — `SolidHarmonics` / `CartTerm` / c2s → qcMath.**
- Move `SphericalShell(l)` + `CartTerm` + the c2s type into `qchem.Math.Angular`; unify with qcSymmetry's
  `HarmonicC2S`.
- `PG_Spherical_MnD` consumes qcMath's `SolidHarmonics` to build `SphData::comps[].terms`.
- qcSymmetry's `SphericalShellRep` takes the qcMath c2s type directly; its tests use the canonical
  `SphericalShell(2)` d-harmonics instead of hardcoding.
- Verify: full UTMain + UTSymmetry (spherical SCF energies byte-identical).

## Risks / notes

- **Core-type surgery.** `Polarization` threads through the entire Cartesian integral engine
  (`GaussianRF`, `PGData.pols`, `CartTerm`, `Hermite`), so Stage 1's regression bar is the **whole** UTMain
  energy suite, not just the symmetry tests.
- **Map-order preservation.** Lexicographic `Monomial::operator<` ≡ the old `LMax`-radix order for valid
  exponents, so `map<Polarization,size_t>` behaviour is unchanged. Confirmed the only live such map is
  `PG_Cart_MnD/Internal/Hermite.C:87`.
- **No new library edges** — both consumers already link qcMath. This is the whole point of choosing qcMath.
- **Naming.** `qchem.Math.Angular` (general) over `SolidHarmonics` (too narrow) so Wigner/CG helpers can join.
- **Payoff is asymmetric.** `Monomial` removes a real type duplication; `SolidHarmonics` removes a softer one
  (type-in-symmetry vs data-in-evaluator + test hardcoding). Both cohere into one clean unit, but Stage 1 is
  the higher-value half — fine to ship it alone and leave Stage 2 for when convenient.

## Relationship to other work

- Surfaced by the ShellRep DIP refactor (`2cb1283f`) and the `rmat3d_t`/`rvec_t` cleanup (`96825073`);
  see `doc/SphericalSALCPlan.md`.
- Independent of Spherical SALC S3b (libcint-spherical), which is separately parked.
