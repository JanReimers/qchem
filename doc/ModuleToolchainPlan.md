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

- **Step 5 — REPLACE googletest (action item, user, 2026-09-21).**  Two reasons, one of them disqualifying:
  **(a) gtest silently passes an empty selection.**  On 2026-09-20 we found that every ctest sweep since the
  2026-09-15 renaming had run NONE of the 48 `Γ`-named SCF tests: CMake's JSON discovery double-encoded the
  name, the filter matched nothing, gtest printed a WARNING and **exited 0**, ctest said "Passed".  Verified
  on gtest 1.16: `--gtest_filter=Nope` → rc 0; the new `--gtest_fail_if_no_test_linked` does not cover it.  A
  runner that can report success while running nothing defeats the purpose of a runner (user).  The interim
  guard is `qchem_discover_tests` in the root CMakeLists (text-listing discovery + `FAIL_REGULAR_EXPRESSION
  "no tests were run"`), which makes the NEXT such accident red — it does not make gtest honest.  **(b) gtest
  is the preprocessor island** of this plan: `TEST`, `EXPECT_*`, `ASSERT_*` are macros, `gtest.h` pulls the
  stdlib in textually, and it is the reason test TUs cannot take `import std;` (the caveat above).
  **Candidate: Boost.UT (μt, boost-ext/ut)** — C++20, single header, **macro-free** (`"name"_test = []{
  expect(x == 1_i); };`, `expect(approx(a, b, tol))`), ships as a module (`import boost.ut;`), no
  dependency on Boost proper.  Catch2 v3, doctest and snitch are all macro-based (`TEST_CASE`, `REQUIRE`) and
  fail criterion (b), whatever their runners do.  **Evaluation gates, in order:** (1) does μt's runner FAIL
  (non-zero exit) on an empty filter / zero registered tests — test it first, it is the reason for the item;
  (2) a UTF-8 test name survives its listing and filtering; (3) the ctest bridge — there is no upstream
  `ut_discover_tests`; we write our own listing→`add_test` script, which we effectively own already
  (`qchem_discover_tests`), with the name carried as bytes; (4) the IDE cost, stated honestly: the user drives
  tests through the VSCode C++ TestMate tree, which speaks gtest/Catch2/doctest, not μt — a TestMate-compatible
  listing/reporting mode or a different explorer is part of the price; (5) the migration: ~900 cases,
  `TEST(Suite, Name)` → `suite`/`"Name"_test`, `EXPECT_NEAR` → `expect(approx(...))`, `ASSERT_TRUE(o) << Why(o)`
  → μt's `expect(...) << msg` with `fatal` — mechanical, scriptable, and the test-name grammar
  (`scripts/testgrid`) must keep parsing.  Do it as a per-exe migration (one `UT*` exe at a time, both runners
  in the tree meanwhile), ITMain last.  Sequence: after Step 2 (so the test TUs are the only `#include`
  islands left and the win is measurable), before Step 4.

## Acceptance per step

Standard: clean `UTMain` (Release) + `allTests` build, `-A_*` fast suite + PW/DFT anchors green. For Step 4,
also watch BMI sizes on disk (the whole point) and incremental-build wall time.
