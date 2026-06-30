# Test → Facade Migration — make `qchem::Calculation` the single molecular recipe, slim `QchemTester`

The "4th topic". Goal: the molecular unit tests drive `qchem::Calculation` (the public front door)
instead of the `QchemTester`/`TestMolecule` scaffold, so the facade is **the** molecular recipe and the
scaffold shrinks to what only it covers. This is the litmus test we flagged: if the facade is good, the
molecular half of `QchemTester` evaporates. Sequenced **right after the namespace sweep** (so migrated
tests are written in the final namespaces). Goal is "thumbs-up when ready," not speed.

## Why

- One recipe, not two. Today the production assembly lives in BOTH `QchemTester` and the facade — the
  exact drift risk the facade was meant to kill. Migrating the tests makes the facade the sole molecular
  path, continuously exercised by the suite.
- Molecular tests start reading like end-user code (`Calculation calc(water, {.basis="dzvp"}); …`), which
  is itself documentation for the GUI/binding consumers.
- Measurable payoff: `TestMolecule` deleted, `QchemTester` slimmed. (Track the LOC drop as the litmus.)

## Scope

**Migrate (molecule-only path):** the `TestMolecule`-derived fixtures —
- `M_HF_U.C` (water/N2 HF, + the libcint / spherical / libcint-spherical variants),
- `M_DFT.C` (Xα water/N2, + spherical variant),
- the molecular fixtures in `A_DFT.C`,
- `M_Sym.C` — already bypasses `QchemTester` (drives `SCFIterator` directly via `PG::SymmetryAdapt`);
  fold it onto `Calculation{.symmetry=true}`, which now covers exactly that path.

**Leave on the scaffold (out of scope — facade is molecule-only by decision):** `TestAtom` / `TestDiracAtom`
fixtures (`A_HF_U/P`, `A_DHF`, `A_DFT` atomic, `A_PP`, `A_DFT_U`) and `scfrun.C`. These use atomic bases,
atomic ECs, and the NIST/Dirac **oracle** comparisons the molecular facade has no business with.

## Prerequisite — fill the facade gaps (additive) so molecular tests CAN fully migrate

Audited what the molecular tests call on `QchemTester`; the facade already covers `Energy()`,
`IsConverged()`, `IterationCount()`, `EnergyTerms()`, and auto-picks SAD seed + DIIS-from-start for DFT.
The gaps:

1. **Basis engine/angular selection (the real one).** `CalcOptions.basis` is just a name → the facade
   always builds the default Cartesian MnD basis. The variant tests need `engine=libcint` and
   `angular=spherical`. Add to `CalcOptions` (e.g. `enum class Engine{MnD,LibCint}; enum class Angular{
   Cartesian,Spherical};`) and thread into the `Molecule::Factory` json. **Interaction:** `angular=spherical`
   + `{.symmetry=true}` must stay guarded until the Spherical SALC track (A) lands — throw clearly meanwhile.
2. **Per-molecule Xα alpha** — already covered by `CalcOptions.xalpha`.
3. **Irrep / orbital enumeration for assertions** — facade has `Orbitals(irrep)` but check whether any
   migrated test needs to *enumerate* irreps or read per-irrep occupations; if so, expose a small accessor
   (the wave function already has `GetQNs()`).
4. **Virial / breakdown asserts** — `EnergyTerms()` returns `EnergyBreakdown`; confirm it carries what the
   tests assert (virial ratio etc.); expose if a gap.

## Extract the assertion helpers

`QchemTester::RelativeError(Eref)` (the (E_ref−E)/E_ref check with the ppt/ppb/ppm reporter) is a *test*
helper, not calculation state. Pull it into a tiny test-utility (free function `RelativeError(double E,
double Eref, bool quiet=false)`), usable by facade-based molecular tests. The **oracle** helpers
(`RelativeHFError/DFTError/DHFError` vs NIST atomic energies) are Z-keyed and stay with the atom scaffold.

## Migration mechanics (per file, pure refactor)

For each molecular fixture: replace `TestMolecule(mol)` + `GetHamiltonian` override + `Init(js)` +
`Iterate(scf)` with `Calculation calc(mol, {.basis=…, .model=…, …});` and assert via the free
`RelativeError(calc.Energy(), anchor)`. **Anchors stay byte-for-byte identical** — this is a refactor, not
a re-tuning; any energy change is a bug. Build + run after each file.

## Strip `QchemTester`

Once no test uses `TestMolecule`: delete `TestMolecule` and the molecule-specific virtual hooks/branches
in `QchemTester`/`QchemTesterImp`. What remains is the atom/Dirac harness + oracle asserts. Then consider
renaming `QchemTester` → `AtomTester` (or similar) to telegraph its now-atomic scope, and prune the umbrella
imports it no longer needs.

## Risks / watch-list

- **`scfrun.C`** (the CLI accelerator-tuning driver) depends on `QchemTester` — keep it working (it's
  atom-oriented; the slim atom harness must still serve it).
- **Exact-anchor reproduction** — verified for HF (−76.022903) and LDA (−75.9324615507) through the facade;
  confirm Xα and the libcint/spherical variants reproduce their anchors once the engine/angular gap is filled.
- **Spherical+symmetry** stays guarded until track A.
- Don't migrate atom/Dirac/PP tests — the molecule-only facade can't express them, by decision.

## Stages

1. Facade gaps: `CalcOptions` engine/angular (+ any small accessor a test needs).
2. Extract the `RelativeError` free helper.
3. Migrate molecular fixtures file-by-file (`M_HF_U`, `M_DFT`, `A_DFT` molecular, `M_Sym`), anchors unchanged.
4. Delete `TestMolecule`; slim `QchemTester`; optional `→ AtomTester` rename; keep `scfrun` working.
5. Measure the LOC drop (the litmus).
