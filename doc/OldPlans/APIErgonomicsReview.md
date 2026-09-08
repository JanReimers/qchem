# qchem API ergonomics — an integrator's review

Feedback gathered while building the Python/nanobind binding (`pybind/`,
`viz-demo/`) in 2026-06. Perspective: someone *consuming* the public interfaces to
do ordinary things (build a molecule, run an SCF, sample the density/orbitals),
not developing the internals. Grounded in friction actually hit, not theory.

Each item: **what I hit → impact → "all I should have to write" → suggestion**.
Ordered by payoff.

---

## ✅ The standout win (keep doing this)

`ScalarFunction<double>` (`operator()` + `Gradient`), and the fact that **both
`ChargeDensity` and `Orbital` ARE `ScalarFunction`s**. One sampling loop works for
density, any MO, and the gradient field — zero special-casing. This is the
"expose high-level questions, not getters" philosophy paying off perfectly. The
whole grid-sampling bridge is ~15 lines because of it. The per-module `Factory`
consistency is also good.

---

## 1. No production "run an SCF" entry point  ⟵ biggest gap

**What I hit:** the only place the full assembly recipe exists is
`QchemTester`/`TestMolecule` — in `UnitTests/`, pulling `gtest`. To do the most
basic thing from non-test code I had to reverse-engineer the test harness and
re-implement its 9-step assembly in the bridge:

```cpp
// what I actually had to write (paraphrased from pybind/qchem_bridge.cpp):
Molecule* mol = new Molecule();
mol->Insert(new Atom(8,0,{0,0,0}));  mol->Insert(new Atom(1,0,{0,1.431,1.107})); ...
std::shared_ptr<const Structure> st(mol);
auto* ec    = new Molecule_EC(int(mol->GetNumElectrons()));
auto* basis = BasisSet::Molecule::Factory(json{{"basis","dzvp"}}, mol);
auto* ham   = Factory(Model::HF, Pol::UnPolarized, st);
auto* acc   = SCFAccelerators::Factory(Type::DIIS, {{"NProj",4},{"EMax",...},...});
auto* scf   = new SCFIterator(basis, ec, ham, acc, SeedStrategy::Default, mol);
scf->Iterate({60,1e-7,1e-9,1e2,1e-7,0.5,1e-4,false});
auto* rho   = scf->GetWaveFunction()->GetChargeDensity();   // caller owns
// ...and I own basis/ec/structure but NOT ham/acc (see item 2)
```

**Impact:** the canonical "how do I use this library" lives only in test
scaffolding. Every caller re-derives it. A chemist (or I) can't self-serve.

**Here's all I should have to write:**

```cpp
import qchem;                              // umbrella (item 7) — or: import qchem.Calculation;

qchem::Molecule h2o = {                    // Z + position (bohr)
    {8, {0, 0,      0}},
    {1, {0, 1.431,  1.107}},
    {1, {0,-1.431,  1.107}},
};
qchem::Calculation calc(h2o, {.basis="dzvp", .model=HF, .symmetry=true});  // build + converge

double E                          = calc.Energy();
const ScalarFunction<double>& rho = calc.Density();   // sample rho(r) directly
const ScalarFunction<double>& mo  = calc.HOMO();      // or calc.Orbital(i)
calc.OnIteration([](const SCFProgress& p){ ... });    // opt-in live trace
```

**Suggestion:** promote a `qchem::Calculation` (a.k.a. `MolecularSCF`) facade into
`src/`. It owns the whole graph, exposes `Energy()/Density()/HOMO()/Orbital(i)/
Orbitals(irrep)`, and becomes the *single tested recipe*. I already wrote ~80% of
it in `qchem_bridge.cpp` (`struct Calc`) — it can be lifted almost verbatim.

---

## 2. Ownership is asymmetric and invisible at the call site

**What I hit:** `SCFIterator` **owns** `ham` and `acc` (deletes them) but **not**
`basis`, `ec`, or `structure`. I only learned this by reading the destructor. Of
the things I pass in, two are silently adopted and three remain mine to free —
with nothing in the signature saying so.

**Impact:** a leak-or-double-free trap that requires reading `.C` impls to avoid.

**Here's the contract made visible:**

```cpp
SCFIterator(std::unique_ptr<Hamiltonian>  ham,    // adopted — type says so
            std::unique_ptr<Accelerator>  acc,    // adopted
            const BasisSet<double>&       basis,  // borrowed
            const ElectronConfiguration&  ec,     // borrowed
            const Structure&              st);
```

**Suggestion:** `unique_ptr` parameters for what's adopted, `const&`/raw for what's
borrowed. The compiler then documents (and enforces) ownership. (A `Calculation`
facade sidesteps this for the common path entirely.)

---

## 3. Symmetry is a deep, evaluator-specific decorator (should be a flag)

**What I hit:** to get point-group blocking + aufbau I must
`import qchem.BasisSet.Molecule.PG_Cart.SymmetryAdapt` and call
`PG::SymmetryAdapt(rawBasis, *mol, 1e-4)` — a `PG_Cart`-specific free function
four namespaces deep — then keep `Molecule_EC` for global aufbau.

**Impact:** symmetry is a *conceptual* property of the calculation, but it's
exposed as an artifact of one Cartesian evaluator. Least-discoverable part of the
API — unfortunate, since it's what you most want people to use.

**Plus a naming trap:** `Molecule_EC` (aufbau, what you want) vs `MolecularSym_EC`
(fixed per-irrep occupation stopgap, e.g. `{A1:6,B1:2,B2:2}`). The names are so
similar that even the library author reached for the wrong one in conversation.

**Here's all I should have to write:**

```cpp
qchem::Calculation calc(mol, {.basis="dzvp", .symmetry=true});  // detects point group, blocks H, aufbau
// trace then shows the irrep occupation per level and the detected point group
```

**Suggestion:** (a) make symmetry a basis-factory / Calculation option; the
factory detects the point group and returns the blocked basis. (b) Rename the
fixed-occupation EC to telegraph it's a special case, e.g.
`FixedIrrepOcc_EC`, leaving `Molecule_EC` as the obvious aufbau default.

---

## 4. Namespace placement is unpredictable

**What I hit:** `qchem::ChargeDensity`, `qchem::Orbitals`, `qchem::SCFIterator`,
`qchem::Hamiltonian` are under `qchem::` — but `ScalarFunction`, `Spin`,
`Vector3D`, and the entire `BasisSet::` tree (`BasisSet::Real_BS`,
`BasisSet::Molecule::Factory`) are **global**. I wrote `qchem::ScalarFunction` and
`qchem::BasisSet::Real_BS`; both failed; I had to grep for the real homes.

**Impact:** a consumer can't predict where a public symbol lives.

**Suggestion:** put all public types under `qchem::` (with sub-namespaces).
Mechanical, but removes a whole class of "where does this live" lookups.

---

## 5. `json` config is stringly-typed and undiscoverable

**What I hit:** `{"NProj":4,"EMax":0.1,"EMin":1e-7,"SVTol":5e-9,"type":"DIIS"}` and
`{"basis":"dzvp"}`. I learned every key by grepping tests — no header enumerates
the valid keys, types, defaults, or units, and nothing is checked at compile time.

**Here's all I should have to write:**

```cpp
AcceleratorConfig acc{.kind=DIIS, .nProj=4};   // documented fields, sane defaults
```

**Suggestion:** a documented options struct (with defaults) for the common path;
keep `json` as the escape hatch for the long tail.

---

## 6. `SCFParams` is an 8-field positional aggregate

**What I hit:** `{60, 1e-7, 1e-9, 1e2, 1e-7, 0.5, 1e-4, false}` — eight unlabeled
values. The tell: **every call site in the codebase writes a decoder comment
above it** (`// NMaxIter MinΔρ MinΔFD MinVirial MinFD relax MergeTol verbose`).
When the code always annotates a struct, the struct is too hard to use. (The
unicode field names `MinΔFD`/`MinΔρ` also read nicely but are painful to type.)

**Here's all I should have to write:**

```cpp
calc.Converge({.maxIter=60, .minFD=1e-7});     // everything else defaulted
```

**Suggestion:** ASCII field names + defaults so designated initializers carry the
labels and callers omit what they don't set.

---

## 7. One umbrella `import qchem;`

**What I hit:** the bridge imports ~16 modules by hand
(`qchem.Structure`, `qchem.ScalarFunction`, `qchem.Hamiltonian.Factory`,
`qchem.SCFIterator`, …) and a consumer must know each module's name *and* which
`using`/alias it provides.

**Here's all I should have to write:**

```cpp
import qchem;     // re-exports the public modules AND the public usings/aliases
```

**Suggestion (lib team):** a thin umbrella module `qchem` that `export import`s the
public-facing modules and re-exports the common `using`s. (Mirror the existing
`qchem.Math` umbrella convention — and the CMake gotcha noted in
`project_qcmath_library_split`: an umbrella must not re-export same-library
siblings, which CMake reads as a bogus cycle.) Internals (`*.Internal.*`) stay out,
per the re-export rule in CLAUDE.md.

**Relationship to the snippets above:** all the "all I should have to write"
examples assume this umbrella — `import qchem;` makes `qchem::Calculation`,
`qchem::Molecule`, `ScalarFunction`, … all visible at once. The granular form
(`import qchem.Calculation;`) is the equivalent narrow alternative: same symbols,
but it only depends on (and recompiles for) that one module. Rule of thumb:
end-user/binding code → `import qchem;` (convenience); library-internal code that
wants minimal recompile coupling → granular per-module imports.

---

## Division of labor: lib hands out containers, the binding marshals to raw

A boundary rule worth making explicit (from the lib team):

- **`void*` and raw `double*` are the BINDING's job, not the lib's.** The C ABI
  handle (`void*`) and the `double*` grid buffers in `qcb_api.h` / `qchem_py.cpp`
  are correct *there* — that's the language-interop layer. They must never leak
  back into `src/`.
- **The lib should hand out blaze / `std::` containers** (`rvec_t`, `rmat_t`,
  `ScalarFunction`, …). The transition container → `double*` → NumPy stays on the
  binding side.

**Things currently in `pybind/qchem_bridge.cpp` that are lib-side candidates**
(they return/operate on containers, no raw pointers needed — migrate them to
`src/` returning blaze/std, and let the binding do the final `double*` wrap):

- `Calc` (the whole assembly + converge) → the `qchem::Calculation` facade (item 1).
- `bbox(positions, pad)` → a geometry helper on `Structure` (e.g.
  `Structure::BoundingBox(pad)` returning two `rvec3_t`).
- `sample_scalar` / `sample_gradient`: sampling a `ScalarFunction<double>` onto a
  regular grid → a lib utility returning an `rvec_t` (flat) or `rmat_t`; the
  binding then wraps that buffer as NumPy. (The grid *definition* — origin/
  spacing/dims — is a natural small lib struct too.)

What stays binding-side: the `extern "C"` ABI, the `void*` handle, the
`double*`→`nb::ndarray` wrapping, the Python callback trampoline.

---

## Net

The low-level abstractions (especially `ScalarFunction`, and `Factory`
consistency) are genuinely good. The gap is the **missing production facade**: the
"set up and run a calculation" knowledge lives only in test code, ownership rules
live only in destructors, and config keys live only in grep results. Close that —
a `qchem::Calculation` front door with options structs and visible ownership — and
the library becomes self-service for chemists and bindings alike.

Most of the facade already exists as `struct Calc` in `pybind/qchem_bridge.cpp`;
lifting it into `src/` is the highest-leverage single change.

---
---

# Lib-team response — review of the review (2026-06-29)

Read against the codebase by someone who'd just been working the high-level layer
(SCFIterator / Hamiltonian / ChargeDensity). Verdict: **on the mark.** The
diagnostics are sharp — the "every call site writes a decoder comment → the struct
is too hard to use" tell is exactly right (those decoder comments are sitting above
the `SCFParams` literal in `M_HF_U.C` and `SCFIterator.C` today), and the central
claim holds: the production recipe lives only in `QchemTester`, ownership lives only
in the `SCFIterator` destructor, config keys live only in grep results.

Below: one reframe of the *sequencing*, per-item verdicts with the few caveats, and
a verified scoping result for the symmetry flag (item 3).

## The reframe: build one additive layer, don't refactor internals

The seven items are two different kinds of work, and they should be split hard:

- **Additive consumer layer — low risk, ~80% of the win, do now:** item 1 (facade),
  5 (options struct), 6 (`SCFParams` struct), 7 (umbrella), and the
  container-returning utilities. None of these *change* existing code — they sit on
  top. Shippable without touching a single internal call site.
- **Invasive internal churn — high effort, defer or sidestep:** item 2 (rewrite
  `SCFIterator`'s ownership signature — ripples to every call site, both
  `double`/`dcmplx` instantiations, and the binding) and item 4 (move
  `ScalarFunction`/`Spin`/`Vector3D`/`BasisSet::` under `qchem::` — touches nearly
  every file).

Strategic point: **the facade lets you *sidestep* the invasive items, not do them.**
Once `qchem::Calculation` owns the whole graph, `SCFIterator`'s asymmetric ownership
(item 2) becomes an internal detail behind the front door — so don't churn the
iterator's ctor; make the facade the sole owner and apply the `unique_ptr`-adopted /
`const&`-borrowed idiom *there*, where there are no existing callers to break.
Likewise the umbrella (7) gives consumers one predictable surface, which removes most
of item 4's pain without the file-wide namespace sweep.

## Per-item verdicts

- **1 (facade) — yes, first.** Highest leverage. **Non-negotiable caveat:** it must
  become *the single recipe* — refactor `QchemTester` to sit on it. Two parallel
  assembly paths recreate the exact drift problem being solved.
- **2 (ownership) — real trap, but fix by encapsulation.** (Confirmed concretely: the
  HF-bootstrap work in `SCFIterator::Initialize` had to pass `st` as a raw
  `const Structure*` and fake a null-deleter `shared_ptr` *because* the iterator's
  ownership story is muddy.) Let the facade own everything; leave the iterator's
  signature alone for now.
- **3 (symmetry as a flag) — right UX, scope it honestly.** See the verified result
  below: it's a ~3-line facade branch, safe. Wire it "real for PG_Cart, throw/no-op
  otherwise." The `Molecule_EC` vs `MolecularSym_EC` → `FixedIrrepOcc_EC` rename is
  good and dovetails with the queued symmetry-naming cleanup.
- **4 (namespaces) — real wart, lowest priority.** Highest churn-to-benefit ratio
  here; the umbrella buys most of the discoverability. If ever done, do it as one
  deliberate sweep, not piecemeal.
- **5 + 6 (options/params structs) — easy wins.** Additive, contained; keep `json` as
  the escape hatch. Offer ASCII designated-init names (keep the pretty `MinΔρ` in the
  struct if you like it).
- **7 (umbrella) — yes, already in motion.** The `qchem.Math` precedent and the
  CMake no-sibling-cycle gotcha are correctly recalled.
- **Container-vs-raw division — fully agree.** `void*`/`double*` staying binding-side
  matches "void\* is banished" in CLAUDE.md. `sample_scalar`/`sample_gradient` belong
  lib-side anyway (cube-file export and the viz both want them, not just Python).

## Suggested sequence

1. `qchem::Calculation` facade + `Converge(SCFParams{...})` + options structs
   (items 1, 5, 6) — lifted from `struct Calc`, made *the* recipe.
2. Umbrella `import qchem;` (7).
3. Container utilities into `src/` (`sample_*`, `Structure::BoundingBox`).
4. `FixedIrrepOcc_EC` rename + symmetry flag (3), coordinated with the symmetry-naming
   cleanup.
5. Ownership (2) absorbed by the facade; namespace unification (4) deferred to a
   standalone sweep if ever.

## Decision: molecule-only facade first

`qchem::Calculation` is a **concrete `double`** class for now (not `tCalculation<T>`),
matching the molecule GUI. It can lean directly on `Molecule_EC` aufbau, the PG_Cart
symmetry path, and `BasisSet::Molecule::Factory`. The plane-wave (`dcmplx`)
generalization stays a clean future move: extract `tCalculation<T>` and make today's
class the `<double>` alias — the same rX/cX pattern already used for `tHamiltonian` /
`tSCFIterator`. No second front door, just a templatization when a crystal GUI needs
one.

## Verified: the symmetry flag is a facade one-liner, not a basis-library change

Traced against the `M_Sym` wiring. The entire symmetry machinery already exists and is
tested end-to-end; the symmetry-specific part of a calculation is **one call** between
the basis factory and the iterator:

```cpp
Real_BS*   raw     = Molecule::Factory(js, mol);          // unchanged
auto       blocked = PG::SymmetryAdapt(raw, *mol, 1e-4);  // <-- the ONLY symmetry step
Molecule_EC ec(Ne);                                       // unchanged: UsesAufbau()=true, global across irreps
SCFIterator scf(blocked, &ec, ham, acc);                  // unchanged
```

Everything downstream of the basis is **symmetry-agnostic**: `Molecule_EC(Ne)` does
global aufbau across point-group irreps with zero extra config (it just sees a basis
with N IBS blocks), and `Hamiltonian` / `SCFIterator` / `WaveFunction` iterate IBS
blocks whether there's 1 or 5. So the SALC builder, the `SymmetryAdapted_IBS`
decorator, the per-irrep SCF, and global aufbau are all **done**.

Therefore `{.symmetry=true}` is purely front-end plumbing — a ~3-line branch in the
facade that calls the existing `PG::SymmetryAdapt`. **Placement matters:**

- **Facade (do this) — zero `qcMolecule_BS` change.** The facade has the structure and
  owns the basis; it just calls the existing builder.
- **`Molecule::Factory` flag (avoid) — a layering downgrade.** `SymmetryAdapt` pulls in
  the `Symmetry` / group-theory library and `PG_Cart_MnD`; making the *generic* basis
  factory call it drags group theory into raw basis construction. Keep that dependency
  at the calculation layer, where "this molecule is C2v" actually belongs.

Net: this is a *narrower* change than feared — no new machinery, no basis-library edit,
guarded to PG_Cart, with `M_Sym` already proving the path. A safe upgrade.
