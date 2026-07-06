# scfrun: add a molecule mode

Spec for a future session.  Goal: extend `UnitTests/scfrun.C` (today an **atom-only** command-line SCF
driver) to also run **molecules**, chosen from a small hard-coded menu on the command line — no xyz files.

Status: NOT STARTED. This is a design note; nothing below is built yet.

---

## 1. Why it's easy now

Everything the molecule path needs already exists and is exercised by the unit tests:

- **Fixture:** `TestMolecule` (`UnitTests/QchemTester.C:108`, ctor `TestMolecule(Cluster*)`) is the molecular
  analogue of `TestAtom`. It already drives `Iterate(SCFParams)`, `TotalEnergy()`, `GetIterationCount()`,
  `Converged()`, accelerator config, and holds the built basis (`itsBasisSet`) + mesh (`GetMeshParams()`).
  See `UnitTests/M_HF_U.C` (HF) and `UnitTests/M_DFT.C` (DFT) for the exact usage.
- **Basis selection:** `BasisSet::Molecule::Factory(json, cluster)` now takes the three orthogonal axes we
  built — pass them straight through from the CLI:
  `{ "basis": "dzvp", "engine": "mnd"|"libcint", "angular": "cartesian"|"spherical" }`
  (`src/BasisSet/Molecule/Factory.C` documents the axes; `QchemTester::Init(js)` forwards `js` to it.)
- **Hamiltonian:** `qchem::Hamiltonian::Factory(Model::HF, Pol, cluster)` for HF; the DFT-Xα overload is
  `Factory(Pol, cluster, alpha, GetMeshParams(), itsBasisSet)` (see `M_DFT.C:43`).

So the upgrade is mostly: a molecule menu + a `CliMolecule` fixture + wiring the three basis axes — almost all
of `scfrun.C`'s existing machinery (arg parsing, accelerator JSON, convergence flags, `theGlobalCache`
setup, reporting) is reused verbatim.

---

## 2. CLI design

Keep the atom path 100% backward-compatible.  Add a molecule path selected by a new `--mol` flag (presence of
`--mol` ⇒ molecule mode; otherwise the existing `--Z` atom mode runs unchanged).

New / reused flags in molecule mode:

```
  --mol <name>       N2 | H2O | CH4 | TiCl4   (selects mode = molecule)
  --basis <name>     dzvp | dzvp2 | tzvp | orb | orb1      (BasisSetData; default dzvp)
  --engine <name>    mnd | libcint                          (default mnd)
  --angular <name>   cartesian | spherical                  (default cartesian)
  --model <name>     HF | DFT                                (default HF)
  --alpha <float>    Xα exchange parameter (DFT only)        (default 0.7)
  --pol <U|P>        UnPolarized | Polarized                 (default U)
  --maxiter, --accel, --nproj/--emax/..., --minfd/--minde/--virial/--minro/--relax  (reused as-is)
```

Atom-only flags (`--Z --q --acc`) are ignored in molecule mode (and vice-versa).  Note `--basis` is overloaded
between the two modes (atom basis Type vs molecular BasisSetData name) — that's fine since the mode is known
before it's interpreted.

---

## 3. Molecule menu (geometries in BOHR)

`Atom`/`Molecule` coordinates are in **Bohr** (1 Å = 1.8897 Bohr).  A `std::map<string, Molecule*(*)()>` of
builders, each `new Molecule(); m->Insert(new Atom(Z, q=0, rvec3_t(x,y,z)))`.  Suggested geometries
(experimental-ish; exact values don't matter for a tuning driver):

- **N2** (14 e, closed shell): N(±1.037, 0, 0).
- **H2O** (10 e, closed shell): O(0,0,0), H(0, ±1.431, 1.107).   *(matches M_HF_U)*
- **CH4** (10 e, Td, closed shell): C(0,0,0); C–H = 2.067 Bohr, H at the 4 tetrahedral dirs ×(2.067/√3=1.193):
  H(1.193,1.193,1.193), (1.193,−1.193,−1.193), (−1.193,1.193,−1.193), (−1.193,−1.193,1.193).
- **Transition-metal** — two options, pick in the session:
  - **TiCl4** (Td, closed-shell singlet, 90 e — the textbook closed-shell TM molecule, but a *big* run):
    Ti(0,0,0); Ti–Cl = 2.17 Å = 4.10 Bohr, Cl at the 4 tetrahedral dirs ×(4.10/√3=2.367):
    Cl(2.367,2.367,2.367), (2.367,−2.367,−2.367), (−2.367,2.367,−2.367), (−2.367,−2.367,2.367).
  - **TiH4** (Td, 26 e) or **ScH3** (D3h, 24 e) — much lighter closed-shell hydrides if TiCl4 is too slow to
    iterate on; H basis is trivial so they exercise the TM orbital basis with far fewer functions.

  **All TM molecules require `--basis dzvp2`** (dzvp2 covers Sc–Zn; dzvp's TM coverage is partial). Confirm the
  chosen molecule is closed-shell and run `--pol U` to be safe; this is a convergence/tuning toy, not a
  benchmark, so spin state need only be "good enough to converge."

---

## 4. The CliMolecule fixture (sketch)

Mirror `CliAtom`; the only molecule-specific parts are the cluster (from the menu), the basis JSON (three
axes), and the HF-vs-DFT Hamiltonian branch.

```cpp
class CliMolecule : public TestMolecule
{
    Model m; Pol p; bool dft; double alpha;
public:
    CliMolecule(Molecule* mol, Model _m, Pol _p, bool _dft, double _a)
        : TestMolecule(mol), m(_m), p(_p), dft(_dft), alpha(_a) {}
    Hamiltonian* GetHamiltonian(cl_t& c) const override
    {
        if (dft) return Hamiltonian::Factory(p, c, alpha, GetMeshParams(), itsBasisSet); // M_DFT.C:43
        return Hamiltonian::Factory(m, p, c);                                            // HF
    }
};
```

Driver flow (molecule mode), reusing the atom path's tail:

```cpp
Molecule* mol = theMolecules.at(name)();                 // menu builder
auto* t = new CliMolecule(mol, Model::HF, pp, dft, alpha);
t->SetAcceleratorConfig(accj);
nlohmann::json bjs = {{"basis",basis}, {"engine",engine}, {"angular",angular}};
t->Init(bjs);                                            // QchemTester::Init(json) -> Molecule::Factory
t->Iterate({maxiter, minro, minde, virial, minfd, relax, 1e-7, true});
cout << "RESULT E=" << t->TotalEnergy() << " iters=" << t->GetIterationCount()
     << " converged=" << (t->Converged()?"yes":"no") << endl;
```

`theGlobalCache = new IntegralsCache_RAM<double>(true);` once at startup (already done in `scfrun.C:125`).

---

## 5. Constraints / gotchas to honour

- **libcint is HF-only** (no DFT 3-centre fit yet). Reject `--engine libcint --model DFT` with a clear message.
- **libcint spherical** is an HF *oracle* in libcint's native harmonic order (`--engine libcint --angular
  spherical`); fine for energies (basis-ordering invariant), but it has no grid eval, so it's HF-only too.
- **DFT** needs `--alpha`, `GetMeshParams()`, and `itsBasisSet` — the separate Factory overload, not a `Model`.
- **No reference energies for molecules** (unlike atoms' `RelativeHFError`). The natural correctness check is
  **cross-engine agreement**: same molecule/basis via `--engine mnd` vs `--engine libcint` should give the
  same HF energy (they match to <1e-10 element-wise — see `M_LibCint`). Consider a `--engine both` convenience
  that runs both and prints the difference.
- **Geometry units are Bohr.** Keep a `1 Å = 1.8897 Bohr` note by the menu.
- `--basis` is mode-dependent (atom Type vs molecular BasisSetData). Document it; don't try to unify.

---

## 6. Suggested increments

1. Menu + `CliMolecule` + HF only, `--engine mnd`, `--angular cartesian`. Get N2/H2O/CH4 converging.
2. Wire `--engine` (mnd|libcint) and `--angular` (cartesian|spherical) through the basis JSON; add
   `--engine both` cross-check.
3. Add the TM molecule (dzvp2) and confirm it converges (may need `--accel GDM` / more `--maxiter`).
4. Add `--model DFT --alpha` (Xα) via the DFT Factory overload. (mnd engine only.)

---

## 7. Related

- Basis axes + engines: `src/BasisSet/Molecule/Factory.C` (BasisSetData × Engine × Angular).
- The integral engines: `PG_Cart_MnD` (native McMurchie-Davidson), `PG_Spherical_MnD`, `PG_LibCint`
  (libcint, cart + spherical).
- Patterns to copy: `UnitTests/M_HF_U.C`, `UnitTests/M_DFT.C` (fixtures), `UnitTests/scfrun.C` (CLI scaffold).
- Perf context: `MnDBench` (`src/BasisSet/Molecule/bench/`) if you want timings while tuning.
