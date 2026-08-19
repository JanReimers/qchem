# Building CP2K 2026.1 from source (the oracle for doc/Benchmark.md)

Source-built on purpose (user, 2026-08-19): we get control over output and instrumentation, and the user
has done this build before.  CP2K is a reference ORACLE we never modify — but a build we own is still the
better tool when the alternative is a packaged version that does not match the banked 2026.1 numbers.

## Survey (2026-08-19, this box)

| | state |
|---|---|
| CP2K source | `~/Code/cp2k`, **2026.1**, `git:5e54ba2`, CMake build system |
| DBCSR source | `~/Code/dbcsr` — cloned `--recursive` (CP2K wants `find_package(DBCSR 2.6 REQUIRED)`; NOT vendored) |
| fypp | **vendored** at `cp2k/tools/build_utils/fypp`, and DBCSR carries it as a submodule.  `fypp-sources.cmake` prefers a system fypp and falls back to the bundled one, so **no pip needed** (there is no pip on this box) |
| CMake | 4.2.3 — **fine**: CP2K declares `cmake_minimum_required(3.24)`, DBCSR `3.22`, both clear of CMake 4's 3.5 cutoff |
| BLAS / LAPACK | `libblas-dev`, `liblapack-dev` already installed (netlib) |
| Fortran | ⛔ **none** — only LLVM flang at `/opt/LLVM-21.1.6-Linux-X64/bin`, and CP2K's `CP2K_Fortran_COMPILER_LIST` is `GNU;Intel;IntelLLVM;NAG;Cray;PGI` (`IntelLLVM` = `ifx`, not flang).  **gfortran is required.** |
| FFTW3 | runtime `.so.3` present, headers missing → `libfftw3-dev` |

## The one privileged step

```bash
sudo apt install -y gfortran libfftw3-dev
```

## Then, no sudo needed

```bash
# 1. DBCSR (serial + OpenMP, no MPI -- matches the banked 'ssmp' oracles)
cmake -S ~/Code/dbcsr -B ~/Code/dbcsr/build -DCMAKE_BUILD_TYPE=Release \
      -DUSE_MPI=OFF -DUSE_OPENMP=ON -DWITH_EXAMPLES=OFF \
      -DCMAKE_INSTALL_PREFIX=$HOME/.local
cmake --build ~/Code/dbcsr/build -j8 && cmake --install ~/Code/dbcsr/build

# 2. CP2K (minimal GPW build: everything optional defaults OFF via CP2K_USE_EVERYTHING)
cmake -S ~/Code/cp2k -B ~/Code/cp2k/build -DCMAKE_BUILD_TYPE=Release \
      -DCP2K_USE_MPI=OFF -DCP2K_USE_FFTW3=ON -DCP2K_USE_LIBXC=OFF -DCP2K_USE_LIBINT2=OFF \
      -DCMAKE_PREFIX_PATH=$HOME/.local
cmake --build ~/Code/cp2k/build -j8
```

Notes: `CP2K_USE_LIBXC=OFF` is fine for the oracle decks — they use CP2K's built-in `LDA_X + LDA_C_VWN`
(== our Slater/Dirac + VWN5).  Turn libxc on later only if a deck needs a functional CP2K lacks natively.

## VALIDATE BEFORE TRUSTING ANY ROW

A freshly built oracle is not an oracle until it reproduces the banked numbers.

```bash
scripts/bench "Si Gamma cp2k" -- <cp2k> -i IntegrationTests/CP2K/si_fcc_gpw.inp
```

Banked (`doc/CP2Kresults.md`): Si Γ **−7.11506** (re-validated −7.11505788), Si 2×2×2 **−7.86744**.
If those reproduce, the build is sound and every row in `doc/Benchmark.md` can be filled.
**Do not quote a CP2K row before this passes.**

## Why RAM needs no CP2K-side hook

Peak RSS is observable from outside: `scripts/bench` wraps both codes in `/usr/bin/time -v` and reads
"Maximum resident set size" — the same kernel `VmHWM` high-water mark qchem now reports for itself
(cross-checked on Si Γ: external 266 MB vs internal 266.7 MB).  One wrapper, both codes, comparable columns.
