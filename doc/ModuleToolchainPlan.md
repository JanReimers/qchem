# Module / Toolchain Modernization Plan

Goal: finish the job of banishing the C preprocessor from the project by turning every textual `#include`
into a real module import — stdlib via `import std;`, our Blaze fork and nlohmann::json via proper modules.
Triggered by a recurring papercut: a module TU that only *imports* `qchem.Types` (which `#include <complex>`
in its global-module fragment and `export using dcmplx = std::complex<double>;`) sees the **type** but not
`std::operator*(complex,complex)` — because a textual include in module A is not visible to module B. Today's
workaround is `#include <complex>` in each such TU's global-module fragment (e.g. `PlaneWaveFit_IBS.C`,
`Imp/PlaneWave_IBS.C`). `import std;` is the *root* fix.

## Toolchain reality (2026-07)

- **Build compiler: Clang** (`.pcm` BMIs; currently `clang++ 21.1.6`, `/opt/LLVM-21.1.6-Linux-X64`).
- **GCC is OUT of the picture for now** — its C++20-module + `import std` support is too far behind Clang's.
  Do not target GCC for the module build; do not add GCC-compatibility workarounds. (`g++` may still exist on
  the box for one-off checks, but it is not a supported build path.)
- CMake 4.2.3, ninja 1.13.2 (bleeding-edge modules stack; see `reference_ninja_dyndep_recovery`).
- CMakeLists.txt already has the `import std` toggle **scaffolded and commented out**:
  `# set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "0e5b6991-d74f-4b3d-a41c-cf096e0b2508")` and
  `# set(CMAKE_CXX_MODULE_STD cxx_std_23)`. So step 1 is a toggle, not a build-system research project.

## The one structural caveat — gtest & nanobind are `#include` islands

gtest (`#include "gtest/gtest.h"`) and nanobind (`pybind/`) are preprocessor-heavy and non-modular; they pull
stdlib in textually. A single TU that does **both** `import std;` and `#include` of those can hit duplicate-
declaration / ambiguity. Resolution: **`import std;` is for the LIBRARY TUs; test TUs and `pybind/` glue stay
on `#include`.** That is where the preprocessor pain and BMI bloat actually live anyway, so nothing is lost.

## What each lever actually buys (don't conflate them)

- **`import std;` → ergonomics.** Kills stdlib `#include`s and fixes the operator-visibility class of bug.
  It does **NOT** fix the big BMIs.
- **Modular Blaze → BMI size + compile time.** The 50-80 MB BMIs are dominated by Blaze's expression-template
  headers being absorbed *textually* into `qchem.Blaze`'s BMI (and re-absorbed anywhere it is re-exported —
  hence the standing "never umbrella Blaze" rule). Making Blaze a real module (`import blaze;`, built once,
  referenced) is what collapses those numbers and speeds incremental builds.

## Sequenced steps (independent, each individually verified — NOT big-bang)

Stacking Clang-modules + experimental-CMake-import-std + modular-Blaze-fork + dev-branch-json all at once
multiplies "which layer broke?" debugging (the same tax as the ninja dyndep bug). Keep each step separate and
green before the next.

- **Step 0 — Upgrade Clang 21.1.6 → 22.1.6 (latest).** Newest module + `import std` fixes; do this first so the
  rest is on the best-supported base. Verify a clean `UTMain` + `allTests` build on 22.1.6 before changing any
  code. (GCC remains out — Clang-only.)
- **Step 1 — `import std;` spike on ONE leaf lib** (qcMath or qcCommon). Uncomment the two CMake lines, bump
  that lib to C++23 + `import std;`, delete its stdlib `#include`s, confirm it builds and the `<complex>`-style
  problem is gone. Cheapest step; directly validates the root fix.
- **Step 2 — Roll `import std;` across the library TUs.** Leave gtest/`pybind/` TUs on `#include` (the island
  rule above). This is the preprocessor-elimination win for the library half.
- **Step 3 — nlohmann::json module** (used in Atom/Molecule factories + several tests). A contained proof-of-
  concept for "modular 3rd-party dep alongside `import std;`". Asterisk: json's module is on `develop`
  (unreleased) — pinning a dependency to a dev branch is a maintenance smell; revisit before making permanent.
  Optional / skippable.
- **Step 4 — Modularize the Blaze fork.** The big BMI/compile-time lever, uniquely feasible because we own the
  fork. Do it LAST: on top of C++23 (so modular Blaze can `import std;` itself) and after json has de-risked the
  modular-3rd-party path. Expect expression-template/ADL visibility quirks — Clang 22 handles these best.

## Acceptance per step

Standard: clean `UTMain` (Release) + `allTests` build, `-A_*` fast suite + PW/DFT anchors green. For Step 4,
also watch BMI sizes on disk (the whole point) and incremental-build wall time.
