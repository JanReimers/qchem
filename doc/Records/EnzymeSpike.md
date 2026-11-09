# The Enzyme AD spike (lane B, 2026-11-08)

RECORD.  The question was set in `doc/GeometryOptimizationPlan.md` ("Automatic differentiation spike", user 2026-11-07): before any term
gradient is hand-derived, try Enzyme on the Becke-mesh δ XC energy at fixed D on H₂O.  Check it against the per-term finite difference,
check Σ_A F_A = 0, measure its cost against one energy and its tape memory, and record the friction with modules, Blaze and OpenMP.  Then
take dE/dD by AD and check it against the term's H (pin 32).  The verdict is in the plan; this file holds the evidence.

## 1. Toolchain

- **Enzyme** `~/Code/Enzyme`, commit `1a87ccc`, built against `/opt/LLVM-21.1.6-Linux-X64` (Release).  It produces three plugins in
  `~/Code/Enzyme/build/Enzyme`: `ClangEnzyme-21.so` (the one we use), `LLVMEnzyme-21.so` and `LLDEnzyme-21.so`.
- **Two fixes were needed to link it.**
  1. The `/opt` LLVM archives are LLVM bitcode (that release is ThinLTO-built), and GNU `ld.bfd` cannot read them.  Link with lld.
  2. The `/opt` `ld.lld` needs `libxml2.so.2`, which this box does not have (it has `.so.16`).  A private shim in `~/Code/Enzyme/build/lldshim/`
     fixes it: an `ld.lld` wrapper that sets `LD_LIBRARY_PATH`, plus `libxml2.so.2 -> /usr/lib/x86_64-linux-gnu/libxml2.so.16`.  lld uses
     libxml2 only for COFF manifests.  The build passes `--ld-path=<shim>/ld.lld` in the three `CMAKE_*_LINKER_FLAGS`.  Nothing system-wide changed.
- **In qchem** the path was `-DQCHEM_ENZYME=ON -DQCHEM_ENZYME_PLUGIN=…/ClangEnzyme-21.so -DQCHEM_ENZYME_INCLUDE=…/enzyme/include`, default
  OFF, and then every AD gradient refused with `NoGradient`.  **Superseded 2026-11-09: Enzyme is REQUIRED** (user ruling).  `QCHEM_ENZYME`
  is retired; configure finds `~/Code/Enzyme-main` (871de19f, ⊇ v0.0.308, the same build recipe) by default, checks the pin and runs a
  smoke derivative (root `CMakeLists.txt`, `cmake/EnzymeSmoke.cpp`).  The plugin flags ride `qcAD` and the files named by
  `qchem_enzyme_sources` (e.g. `BeckePartition.C`), nothing else.  Toggling the old option rebuilt `qcAD` and
  nothing else: `ninja -n` claims ~600 objects, but that is a dry-run artifact of the dyndep step.

## 2. What was built

- **`qcAD`** (`src/AD`): `qchem.AD` (`NoGradient`, `Available()`) and `qchem.AD.XCKernel`.  The XC kernel is a pure function E_xc(R, D) on
  plain arrays: one mesh template shifted to every atom, the Becke partition, Cartesian Gaussians riding with their atoms, and Dirac +
  VWN5 paramagnetic.  Its gradients are dE/dR and dE/dD.
- **Tests:**
  - `UTAD`: the kernel's Becke weight equals `qcMesh::BeckePartitionWeight`, and its AD gradient equals the hand-derived
    `BeckePartitionGradient`.  A synthetic system checks both gradients against the kernel's own FD.
  - `UTCalculation` `ADSpike.*` (in `src/Calculation/tests/ForceFacade.C`, through the `ForceFacadeTests` friend): the checks against the
    library's term on converged H₂O LDA / dzvp.

## 3. Results (H₂O, LDA, dzvp: 19 Cartesian functions; mesh (80, 29): 72,480 points; Becke order 3)

| check | result |
|---|---|
| kernel E_xc vs the library's δ XC term at the converged D | −8.7436664411944 vs −8.7436664411952 Ha (8e-13) |
| kernel basis values vs the library's IBS at 4 points | ≤ 1e-13 relative |
| AD dE/dR vs the **per-term FD** (library term at R ± 1e-4, D frozen; basis, points and partition all move) | 9 components agree to ≤ 2.3e-9, the FD's own O(h²) error |
| Σ_A dE_xc/dR_A | ~1e-16, machine zero (the grid response is in) |
| AD dE/dD vs the term's H (`MakeTermMatrices`, pin 32) | max \|ΔH\| 2e-12, max\|H\| 3 |
| AD Becke-weight gradient vs the hand-derived one | ≤ 1e-13, and it sums to zero |
| determinism | gradients bit-identical run to run, and at 1 and 4 threads |

Cost, against one energy evaluated by the same kernel at `-O2`:

| | 1 thread | 4 threads |
|---|---|---|
| E_xc | 0.045 s | 0.012 s |
| dE/dR (9 values) | 0.75 s, **17×** | 0.26 s, **21×** |
| dE/dD (361 values) | 0.29 s, **6.5×** | 0.095 s, **7.7×** |
| tape (peak RSS over the call) | < 1 MB | < 1 MB |

- **Tape.** Each call differentiates one block of ~280 points, so the tape is per block.  The whole kernel differentiated in one call
  (before the restructure) took 17 MB.
- **Cost against the alternatives.** The per-term FD of the same 9 components took 7.8 s: 18 rebuilds of the point.
- **Cost against hand-coded.** A hand-coded XC gradient (Johnson, Gill & Pople; it needs ∇φ) is typically 3–4× an energy, so AD's dE/dR is
  about 5× a hand-coded one here.  Nothing was tuned: no caching of φ between blocks, and `-O2` instead of the tree's `-O3`.

## 4. Friction (each one was hit, not predicted)

1. **Enzyme sees only bodies in the same TU, or in one LTO link.**
   - A call into another TU is refused at compile time ("No augmented forward pass found").
   - A call into another shared library is refused the same way.  Every qchem library is a `.so`, so LTO cannot help across them.
   - LTO plus `LLDEnzyme` worked across two TUs (probe p3).
   - The spike's kernel therefore restated the library formulas it needed (the Becke cutoff, the functional), each held equal to its
     library home by a test.  That was the demonstration, not the architecture: the user's ruling (2026-11-08, the plan's verdict) is that
     each Imp TU differentiates its OWN leaves in place, and the derivatives compose by the chain rule across the faces (item 9).
2. **No virtual dispatch.** A virtual call through a `Const` object is refused ("No create nofree of unknown value").  So Enzyme never
   differentiates THROUGH the OOD term hierarchy (`FittedVxc` → `ExFunctional` → …).  Each leaf behind a face is differentiated where it
   is defined, and the faces compose those derivatives.
3. **OpenMP.**
   - A `reduction` clause is refused (`__kmpc_reduce_nowait`).
   - A simple `parallel for` with per-block partials works and threads its reverse pass (probe p2b).
   - The XC kernel's parallel region failed: a vector allocated inside it segfaulted, and with the partials passed in, the heap was
     corrupted ("malloc(): unaligned tcache chunk") whenever the region captured a constant pointer beside an active one.
   - The fix keeps OpenMP outside the differentiated code: each serial block is differentiated, the calls run in parallel, and their
     gradients are summed in block order.  This is deterministic, and it keeps the tape per block.
4. **A WRONG gradient, silently: Enzyme's loop caching (the dangerous one).**
   - Minimised after the spike to `doc/Records/enzyme_loopcache_repro.cpp` (§6, issue A), about 30 lines.  A callee fills a loop-local
     stack array through an inner loop whose trip count is read from data, and the reverse pass returns a wrong dE/dD at EVERY `-O` level:
     0.851 against the FD's 4.571.
   - The kernel met it in a shape-dependent way.  The spike first saw it as "wrong at `-O3`, right at `-O2`" on the whole kernel, and as
     dE/dD ≡ 0 in the tree.  The tree's present kernel shape is exact at the flags it builds with, which is why **every AD gradient is
     gated by a finite-difference test**: no flag set makes Enzyme safe.
   - **A segfault, separately.** At `-O2 -march=native` (AVX2, i7-10700) the reverse pass of the spike's three-block kernel crashes:
     `doc/Records/enzyme_flags_repro.cpp` (§6, issue B, not minimised below that file).  Without `-march`, with
     `-mprefer-vector-width=128`, or with one block, it is exact.
   - The tree compiles `-O3 -march=native`, and `-march` must stay (a .pcm and its importers agree on target features).  So `qcAD` appends
     `-O2 -mprefer-vector-width=128`.
5. **RETRACTED (2026-11-08, after minimising): "Enzyme's tuple wrapper `enzyme::autodiff` zeroes the gradient".**  On a correct case the
   wrapper gives the right gradient.  On the loop-cache case it gives the same wrong one as the marker entry.  The zero seen in the tree
   was defect 4 in another shape.  `qcAD` stays on the classic marker entry (`__enzyme_autodiff(fn, enzyme_dup, x, dx, enzyme_const, ctx)`),
   because it states each argument's activity explicitly; it was never a fix.
6. **C++ details.**
   - Enzyme's entry is a declared-never-defined template, so an argument type in an anonymous namespace is ill-formed ("used but not
     defined … its type does not have linkage"); the context structs need module linkage.
   - A mesh node at r = 0 (MHL's zero-weight node) makes d√x/dx infinite, and 0·∞ is a NaN in the reverse pass, so it is skipped before any
     distance is taken.
   - The active array must not also be reachable through the constant context (the D gradient differentiates a copy).
7. **Modules: no friction beyond the above.** A module interface plus implementation unit compiled with the plugin works.  It imports
   `qchem.Types` and `qchem.Outcome` .pcm's built without the plugin, and their importers are untouched.
8. **Blaze: worked in a probe** (`x^T A x + Σ exp(x)` through Blaze expressions, BLAS mode on and off).  The kernel itself uses plain
   arrays, with Blaze only at its boundary.

9. **What Enzyme can see from an imported module** (scratch probes, 2026-11-08).
   - An exported `inline` function and an exported template differentiate in the importer.
   - An exported non-inline function is "No reverse pass found" in the importer, and so is an in-class member of an exported class (not
     implicitly inline in a module).
   - Both differentiate when the `__enzyme_autodiff` call sits in their own defining unit, and the importer of that derivative needs no
     plugin: the in-place architecture (the plan's verdict).
10. **The AD flags' cost on ENERGY code** (no plugin calls, serial).
    - The XC kernel energy (72,480 points, 16 functions) is 0.0290 s at `-O3 -march=native`, 0.0302 s at `-O2 -march=native`, and 0.0288 s
      at `-O2 -march=native -mprefer-vector-width=128`.
    - A strict-order φᵀDφ contraction (nbf 200, 2e5 points) is 7.84 s and 7.91 s, and 7.94 s with no AVX at all: it is not
      autovectorised.
    - A .pcm built at either flag set imports into a TU at the other.

## 5. What the spike did not cover

- Forward mode and Hessian-vector products (level 2).
- A polarized density, spherical shells, GGA, and the δ Hartree.
- A larger system's cost scaling: the block structure is linear in points, as the energy is.
- Upstream bug reports for friction items 3 and 4 (OpenWork row).

## 6. Draft upstream issues (NOT filed: posting to GitHub is the user's call)

Versions:
- Enzyme `main` at `1a87ccc` (2026-11-07, after v0.0.307), the clang plugin `ClangEnzyme-21.so`.
- LLVM / clang 21.1.6, the official Linux x86-64 release binary.
- Ubuntu, Intel i7-10700 (AVX2, no AVX-512).

**Re-run on Enzyme `main` 871de19f (2026-11-08), which contains the latest release v0.0.308.**  It was built by the controller's session on
the user's approval, with the same `/opt` LLVM 21.1.6 and lld shim; no newer clang is needed.  Both defects are UNCHANGED, value for value:

| repro | flags | 1a87ccc | 871de19f (⊇ v0.0.308) |
|---|---|---|---|
| A `enzyme_loopcache_repro.cpp` | `-O1`, `-O2`, `-O3`, `-O2 -march=native` | AD 0.8510516320, FD 4.5709930… | the same |
| B `enzyme_flags_repro.cpp` | `-O2` | exact (−1.1419468942) | exact |
| B | `-O3` | AD −0.0130734452 (wrong) | the same |
| B | `-O2 -march=native` | SEGFAULT | SEGFAULT |
| B | `-O2 -march=native -mprefer-vector-width=128` | exact | exact |

So v0.0.308's 562c47b4 (the CacheUtility loop-sublimit clamp) does not cure A.  The issues below can say "reproduced on main 871de19f".

The spike's own checks also pass on 871de19f (`qcAD` and `BeckePartition.C` built against its plugin, 2026-11-08), with the same numbers:
- `ADSpike.*` 3/3 at 1 and 4 threads: the kernel equals the term to 8e-13, Σ_A F_A ≈ 1e-16, and max|dE/dD − H| = 2e-12.
- Cost: dE/dR 15× / 22×, dE/dD 6.0× / 7.6×.
- `UTAD` 3/3, and `UTMesh BeckePartition*` 3/3: the in-place Becke leaf equals the hand oracle to 1e-13.

The pin can move to 871de19f whenever useful.  Nothing in the tree depends on the difference.

**Issue A.**  Title: "Reverse mode returns a wrong gradient when a callee fills a loop-local array through an inner loop with a
data-dependent trip count".

> Enzyme main 1a87ccc and 871de19f (after v0.0.308), clang/LLVM 21.1.6 (release binary), Linux x86-64, Intel i7-10700.  The attached file (30 lines, no headers beyond <cmath>/<cstdio>)
> computes E(D) = Σ_p (φ_pᵀ D φ_p)², where φ_p is filled by a callee `Fill` through an inner loop whose trip count `nk[i]` is read from a
> constant data struct (all 1 here).  `__enzyme_autodiff(E, enzyme_dup, D, dD, enzyme_const, &data)`, then dD[0] is compared with a
> central finite difference.
> Expected: `AD 4.5709930…  FD 4.5709930…`.  Actual: `AD 0.8510516320  FD 4.5709930507`, at -O1, -O2 and -O3, with or without
> -march=native.  No warning or error is printed.  Replacing the inner loop by its single iteration gives the correct value.
> Forward mode (`__enzyme_fwddiff` with the tangent on D_00) gives the correct 4.5709930494 on both commits, so the defect is in reverse mode.
> Build: `clang++ -std=c++17 -O2 -fplugin=ClangEnzyme-21.so enzyme_loopcache_repro.cpp`.

**Issue B.**  Title: "Reverse-mode segfault at -O2 -march=native (AVX2) for a blocked reduction into a duplicated scratch buffer".

> Same versions (both commits).  The attached file (~150 lines, self-contained) differentiates a blocked sum: 3 blocks write partials into a buffer passed
> `enzyme_dup`, and are summed after the loop.  At `-O2` it is correct.  At `-O2 -march=native` on an AVX2 machine the reverse pass
> segfaults.  `-mprefer-vector-width=128`, or a single block, avoids it.  (The same file at `-O3` also shows issue A's wrong value.)
