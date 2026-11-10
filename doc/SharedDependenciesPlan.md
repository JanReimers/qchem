# Shared dependencies across worktrees — DESIGN NOTE (for review, no code yet)

OpenWork §4 row "Dependencies are vendored per worktree".  Lane G, 2026-11-10.  Status: **design for the user's and the
controller's look**; it touches the root `CMakeLists.txt` and every tree's build, so nothing is built until it is ruled.
The user: *"Each time we add repo tree we should upgrade the priority of the submodule sharing work item"* … *"no huge rush."*

## 1. What is duplicated today (measured 2026-11-10)

| item | per tree | trees | note |
|---|---|---|---|
| `submodules/` sources (8 submodules) | 620 MB (libxc 291, blaze 174, libcint 127, json 20, …) | 9 worktrees + their `.git/modules` | every tree pins the SAME eight commits (`git submodule status`, all trees identical) |
| libxc + libcint BUILT (`add_subdirectory`) | ~51 MB per build tree (`build/*/submodules`) | ~15 build trees (Release, RelChecked, Dev, TSan) | rebuilt in every tree and every flavour |
| nanobind (pybind only, `FetchContent` v2.4.0) | 11 MB | lane-g only | already FetchContent |
| Enzyme | — | one, `~/Code/Enzyme/build` | RULED toolchain (2026-11-08), not a source dependency |

Disk: 45 % of 468 GB used.  The cost that grows is not the bytes but the REBUILDS: a fresh tree, or a fresh flavour, recompiles
libxc (≈ the longest single C target) before the first test can run.

## 2. Proposal

**(1) One shared install prefix per (dependency pin, compiler, flags)**, outside every tree:

```
~/Code/deps/<compiler-id>/<dep>-<pin>/        e.g. ~/Code/deps/clang-21.1.6/libxc-d00e5d0e2/{include,lib,lib/cmake/Libxc}
```

* Each dependency is built ONCE per key and installed (`cmake --install`), never edited in place.
* The root `CMakeLists.txt` asks `find_package(Libxc <pin> CONFIG)` / `find_package(cint …)` / `find_package(blaze …)` /
  `nlohmann_json` / `GTest` with `CMAKE_PREFIX_PATH` (or `<Dep>_ROOT`) pointing at the prefix -- a configure-time choice, not
  an environment variable the library reads (the library itself reads none).
* **Outsider fallback, at the SAME pins:** when `find_package` misses, `FetchContent_Declare(… GIT_TAG <pin>)` +
  `FetchContent_MakeAvailable` builds it in-tree exactly as the submodules do now.  So a fresh clone still builds with no
  prefix, and the pin lives in ONE place (a `cmake/Dependencies.cmake` table: name, repo, tag, options), which both paths read.
* The submodules are then RETIRED (the FetchContent table replaces `.gitmodules`); a tree no longer carries 620 MB.
* A small script (`scripts/deps-install`) builds a missing key into the prefix; a key is immutable (a new pin is a new
  directory), so trees on different pins coexist and no tree can break another by rebuilding.

**(2) The user's forks (blaze, libxc): the patch flow.**  Both are forks on GitHub (JanReimers/blaze, JanReimers/libxc); the
pins are commits ON THE FORKS (blaze 7478a70b9, libxc d00e5d0e2 -- the libxc one carries the "write xc.h only when its content
changes" CMake patch, also kept as `~/Code/libxc.patch`).  The flow stays: change the fork, commit, push, bump the pin in the
table; the prefix gains a new key; trees pick it up when they bump.  Nothing patches a dependency inside a tree any more
(today `submodules/libxc` shows "modified content" in some trees -- exactly the drift this removes).  Blaze is header-only
here (only its CMake module is used), so its "install" is its headers + `Blaze_Import` module at the pin.

**(3) NATIVE flags across trees.**  Every tree's Release is `QCHEM_ARCH_EXPERIMENT=NATIVE` (`-O3 -march=native`, 15 of 15
cache entries), but the root CMake applies that codegen to the C++ TUs ONLY: the C and Fortran TUs (libxc, libcint) keep plain
`-O3` (root `CMakeLists.txt`, the NATIVE comment).  So the C dependencies are flag-independent of the tree's arch choice and a
key needs only (compiler, build type, PIC).  Two cautions to keep in the key: (a) RelChecked / TSan / Debug trees -- TSan
needs an instrumented libxc only if libxc's own threads matter (they do not: libxc is called from our threads), so the Release
prefix serves them; ASan/MSan would need their own key; (b) janr1 (the second machine, no FMA) gets its own prefix -- the
prefix is per machine, never on NFS.

**(4) C++ module BMIs are never shared.**  `qchem.Blaze`, `qchem.Math` and every project module compile their BMIs per tree
and per flag set (a BMI encodes the exact flags: `-march=native` once broke `qchem.Math.pcm`, root CMake comment).  Nothing
in the prefix is a BMI: blaze and json are consumed as HEADERS through the tree's own module wrappers, so the wrappers keep
building per tree and the prefix holds only headers and C libraries.  (If a dependency ever ships C++ modules, it gets built
per tree, not installed.)

**(5) Outsiders and CI: an image.**  A container image (or devcontainer) = the toolchain (LLVM 21.1.6 + flang, libomp from
apt, Enzyme built against that LLVM, the lld libxml2 shim) + the prefix for that compiler, so a CI job or a new contributor
configures with `CMAKE_PREFIX_PATH` already set.  **Spack** is the HPC-native alternative (it models exactly the
(pin, compiler, flags) key as a spec hash); worth it if the code is ever installed on a cluster, heavier than an image for the
desktop.  Recommendation: the FetchContent table first (it is useful alone and costs nothing to outsiders), the shared prefix
second, the image third; Spack only on a cluster request.

## 3. What changes, and the order

1. `cmake/Dependencies.cmake`: the pin table + `find_package`-else-`FetchContent` per dependency.  Root CMake includes it in
   place of the four `add_subdirectory(submodules/…)` / `CMAKE_MODULE_PATH` lines.  Gate: a fresh tree with NO prefix builds and
   sweeps identically (0 anchors moved) -- the outsider path.
2. `scripts/deps-install <compiler>`: build + install every pin into the prefix.  Gate: a tree configured against the prefix
   sweeps identically, and its first build skips libxc/libcint entirely.
3. Retire the submodules (`.gitmodules`, `submodules/`) once every tree has moved; the lanes move at their next main merge.
4. The image (separate, later).

## 4. Open questions for the user

* The prefix location (`~/Code/deps/` proposed) and whether janr1 should mirror it.
* Keep googletest as a FetchContent pin too (proposed), or the distro's `libgtest-dev`?
* Retire the submodules in one landing, or keep `.gitmodules` for a transition period?
