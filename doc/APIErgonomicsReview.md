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
import qchem.Calculation;                 // the missing front door

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

## Net

The low-level abstractions (especially `ScalarFunction`, and `Factory`
consistency) are genuinely good. The gap is the **missing production facade**: the
"set up and run a calculation" knowledge lives only in test code, ownership rules
live only in destructors, and config keys live only in grep results. Close that —
a `qchem::Calculation` front door with options structs and visible ownership — and
the library becomes self-service for chemists and bindings alike.

Most of the facade already exists as `struct Calc` in `pybind/qchem_bridge.cpp`;
lifting it into `src/` is the highest-leverage single change.
