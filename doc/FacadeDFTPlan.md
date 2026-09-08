# Facade DFT Plan — adding DFT to `qchem::Calculation`

Plan for the one remaining lib-side facade follow-up: make `qchem::Calculation` run **DFT**, not
just HF. Companion to the API-ergonomics work (`doc/OldPlans/APIErgonomicsReview.md`). Living document.

## First, untangle the two "item 4"s — they are unrelated

- **"Consistent namespaces"** = the *original review's item 4* (put `ScalarFunction`/`Spin`/`Vector3D`/
  `BasisSet::` under `qchem::`). That is the **invasive, whole-tree** change deliberately deferred — it
  touches nearly every file and is *not* part of DFT. The umbrella `import qchem;` already bought most of
  its discoverability benefit. Leave it parked; it's a standalone mechanical sweep if ever done.
- **DFT support** = the *facade follow-up* this doc plans. Additive — a new code path in the facade, no
  existing call sites change.

They share a number by accident; nothing about DFT requires the namespace sweep.

## What HF already does (the baseline)

Today `Calculation` (HF) assembles: orbital basis → `Molecule_EC` → `Factory(Model::HF, Pol, st)` →
DIIS → `SCFIterator` with `SeedStrategy::Default` (= core guess). The Hamiltonian is the *only* piece
that changes for DFT.

## The DFT delta — what actually changes (the "stuff")

DFT is the **same SCF loop** with a different Hamiltonian. Tracing `M_DFT` / `A_DFT`, exactly five
things differ from HF, and three of them are nearly free:

| # | Delta | Cost | Notes |
|---|---|---|---|
| 1 | **Functional** (the Hamiltonian) | the real work | two sub-paths — see below |
| 2 | **Numerical mesh** (`MeshParams`) | small | HF needs none; DFT integrates XC on a grid |
| 3 | **SAD seed** (`SeedStrategy::SAD`) | trivial | DFT seeds from superposition-of-atomic-densities, not core guess |
| 4 | **Fit/auxiliary basis** | **free** | derived from the orbital basis *inside* the Factory — facade just passes `itsBasis` |
| 5 | **Polarized (spin) DFT** | free-ish | `Pol` is already a `CalcOptions` field; DFT-polarized = `Factory(Pol::Polarized, …)` |

So #4 (the scary-sounding "fit basis") is a non-issue — the orbital IBS's `CreateCDFitBasisSet` /
`CreateVxcFitBasisSet` build it internally. The facade owns no new objects.

### #1 in detail — the functional, and the resolution principle

**Maintainer constraint (decided 2026-06-30):** the `Method{HF,Xalpha,LDA}` enum is a *thin front-end
token*. It must resolve into the polymorphic Hamiltonian types **quickly — 1 or 2 function calls, not
5–10 — and that resolution lives in `Factory` (or at the latest `SCFIterator`), never threaded deeper.**
The token is for test-writer/facade ergonomics; the moment it crosses the Factory boundary it becomes a
concrete `Hamiltonian*` and the enum disappears.

→ Concretely: add a single resolving overload, e.g.
`Factory(Method, Pol, st, MeshParams, orbitalBasis, /*xalpha*/)`, that switches once on `Method` and
returns the concrete polymorphic Hamiltonian (`Ham_DFT…` for the DFT methods, the HF term for HF). This
overload is **also the "compact default Hamiltonian" the tests want** — one call, no manual
mesh/functional/fit assembly. The `ExFunctional` machinery (`SlaterExchange`, `Libxc_LDA_Exchange`,
`VWN_Correlation`, under `Hamiltonian/Internal/`) is constructed *inside* this switch and never leaks
out — the enum is the only public handle, so the facade never imports `Internal/*`.

The two DFT functionals behind the enum:

- **LDA — the real molecular default (parameter-free: Dirac exchange + VWN correlation).** This is what
  the **molecular DFT regression test should anchor to**, not Xα=0.7. Widely accepted, parameter-free,
  not-necessarily-modern is fine. Today `A_LDA_U` reaches it via `Ham_DFTcorr_U(st, mesh, basis)`
  directly; the resolver wraps that construction behind the enum's `LDA`.
- **Xα — kept, but demoted to the tuned-`alpha` use (atoms / Slater-α-per-Z), not the molecular anchor.**
  `Xalpha` + a `double xalpha`. Still the simplest path and useful, just no longer what a "did molecular
  DFT work" test pins.

### ONE unified enum (decided 2026-06-30) + the friendly-libxc interface

Unify HF/DFT/relativistic into a **single** `Hamiltonian::Model` enum (extend the existing one; no second
sibling). The 1-electron / Dirac members aren't for the GUI but **are important for unit testing**, so
they stay:

```
enum class Model { E1, HF, DE1, DHF,        // existing: 1-e, HF, Dirac-1e, Dirac-HF  (unit tests)
                   Xalpha, LDA, /*PBE,…*/ }; // DFT methods (the friendly-functional vocabulary)
```

The resolver `Factory(Model, Pol, st, MeshParams, orbitalBasis, xalpha)` switches once: for `HF`/`E1`/
Dirac it ignores mesh/basis/alpha and builds as today; for `Xalpha`/`LDA`/… it builds the DFT
Hamiltonian. One enum, one switch, ≤2 calls — the maintainer constraint holds.

**Friendly-libxc TODO, folded in:** clients should never pass a libxc integer id — they pick a named
functional and we map it to the id(s) + handle params under the hood (`Libxc_LDA_Exchange` already takes
an `int id`, so the LDA rung is a pure enum→id mapping, zero new machinery). The honest catch: **what's
offerable is bounded by the XC rung the code supports, which is LDA-only today.** So the vocabulary splits:

| Functional (enum) | Rung | Relevance here | Status |
|---|---|---|---|
| `LDA` (SVWN5: Dirac X + VWN5 C) | LDA | the parameter-free baseline / molecular anchor | ✅ wired now |
| `Xalpha` (Slater X, tuned α) | LDA | atomic/tuned tests | ✅ wired now |
| `LDA_PW92` (Dirac X + PW92 C) | LDA | the LDA flavour standard in solid-state codes (VASP/QE) | ◑ small: a `Libxc_LDA_Correlation(id)` wrapper (mirror the exchange one) |
| **`PBE`** (GGA_X_PBE + GGA_C_PBE) | **GGA** | **THE materials/solid-state workhorse — #1 for the battery-voltage north-star** | ✗ needs GGA machinery (density gradient on the mesh + gradient-dependent potential) |
| `BLYP` (B88 + LYP) | GGA | the molecular-chemistry standard GGA | ✗ needs GGA machinery |
| `PBEsol` | GGA | PBE retuned for solids (batteries) | ✗ needs GGA machinery |
| `B3LYP` / `PBE0` / `HSE06` | hybrid | molecular gold-standard / solid band gaps | ✗✗ needs exact (HF) exchange mixing |

**Reading of this for the plan:** the friendly-functional *interface* (named enum, hide libxc ids, handle
params) lands now and covers the **LDA family** with zero/▱tiny new machinery. But the functional that
actually matters most for the battery goal — **PBE** — is a GGA, and **GGA support is a genuine library
increment** (gradients on the mesh), not an enum entry. So: ship the named-functional API + the LDA family
in this facade work, and treat **PBE-GGA as the clearly-flagged next functional milestone** (highest value,
own track, coordinate with `project_dft_upgrade_plan` and the battery roadmap). Hybrids are further out
(exact exchange). The enum can *list* `PBE` from day one with the resolver throwing "GGA not yet wired"
until that milestone lands — so the API vocabulary is forward-looking but honest.

## Design decisions

- **A. ONE unified enum (decided).** Extend `Hamiltonian::Model` to `{E1, HF, DE1, DHF, Xalpha, LDA, …}`
  — not a second `Method` sibling. 1-e/Dirac stay (unit-test-only, not GUI). `CalcOptions.model` (default
  `HF`) selects it; the facade passes it straight to the resolver. Keeps the existing name → no rename
  ripple.
- **B. Xα alpha source.** A `CalcOptions` field `double xalpha = 0.7` (molecule-wide), used only when
  `model==Xalpha`. The per-Z Slater alpha the atom tests use is an atomic-physics detail, not a GUI default.
- **C. Mesh — `MeshParams` (designated init) + `Structure*` → `Mesh`, resolved internally.** The *only*
  extra input the mesh needs is the structure; the resolver (Factory side) turns a designated-initialized
  `MeshParams{.nRadial=…, .nAngular=…}` into the concrete `Mesh`. Keep `MeshParams` as the ergonomic
  front-end knob (it already is a designated-init struct with sane defaults), with a molecular default
  mirroring `TestMolecule`: `{MHL, nRadial=30, mhl_m=3, mhl_alpha=2.0, Gauss, nAngular=12, beckeOrder=2}`.
  On `CalcOptions`, surface it as a nested `MeshParams` (or a thin `GridOptions`) rather than flattening
  its ~7 fields — designated init makes overriding just `nRadial`/`nAngular` clean.
- **D. Seed.** Auto-select by method: DFT → `SeedStrategy::SAD`, HF → `Default` (core guess). SAD is
  DFT-only (must not be used with HF). Keep it automatic; no user knob needed in v1.
- **E. The `Method` resolver is the public functional handle.** One switch in `Factory`, returning the
  concrete Hamiltonian; the `ExFunctional` internals stay inside it. This single overload serves the
  facade, the test "default Hamiltonian" need, and keeps the enum from leaking deeper than `Factory`.

## Staging

| Stage | What | Test anchor |
|---|---|---|
| **D1** | The `Factory(Method, Pol, st, MeshParams, basis, xalpha)` resolver (enum→polymorphic, ≤2 calls) + mesh-from-`MeshParams`+`Structure*` + auto SAD seed. Wire `Method::LDA` **and** `Method::Xalpha`; facade `CalcOptions.method`. | **parameter-free LDA** molecular (`Ham_DFTcorr_U`-equivalent) through the front door — the new molecular DFT sentinel |
| **D2** | Polarized (spin) DFT (`Pol::Polarized`) | `A_DFT_U`-style |
| — | GGA / meta-functionals | future, out of scope |

LDA is in D1 now (not a deferred increment): it's the molecular anchor the maintainer wants, it already
exists via `Ham_DFTcorr_U`, and the `Method` resolver *is* the public handle that exposes it cleanly.
Xα rides along in the same switch for the tuned-α atom use. So D1 lands real, parameter-free molecular
DFT through one Factory overload.

## Risks / open questions

- **Enum unification (DECIDED): one `Model`.** Extend the existing `Hamiltonian::Model` with the DFT
  members rather than adding a parallel `Method`. The resolver overload takes the DFT extras (mesh, basis,
  xalpha) and ignores them for the HF/1-e/Dirac members.
- **XC rung is LDA-only today.** The named-functional API + LDA family ship now; GGA (PBE — the
  battery-critical one) and hybrids are separate library increments. The enum may list `PBE` early with a
  "GGA not yet wired" throw, so the vocabulary is forward-looking without over-promising.
- **Keep the enum shallow.** The whole point of the maintainer constraint: resolve in `Factory`, return
  `Hamiltonian*`. If the switch needs the structure/mesh/basis, that's fine (it has them); it must not pass
  `Method` onward.
- **Anchors are "did E move" sentinels** (per the test convention), not physical-accuracy claims — pin
  the converged LDA value, same as HF pins `-76.022903`.
- **Mesh defaults** are a quality/perf tradeoff; start with `TestMolecule`'s proven molecular values.
