# CLAUDE — Project Conventions

Brief notes about module/library conventions, naming, and includes.

## My high level biases

- I used to work on frustrated magnets so one of my biases is that I don't like the "closed shell only, polarized is a special case" mind set. Ultimately I want everything available Pol and UnPol ... and for me UnPol is the special case for efficiency.

## Build & test

- Build & test: `cd build/Release && ninja UTMain`, then `./UnitTests/UTMain` (filter `-A_*` for
    fast runs; a full run including the `A_*` tests is the regression anchor).
- Calling ninja directly is fine. (I earlier suggested cmake only because ninja rebuilds were flaky —
    needing file-touches to trigger them. If ninja works well for you, use it.)
- There are interation tests in the (now mis-named) Unitests folder.  The exe target is UTMain
- But there are also lot of actual unit tests in most of project-module (library) folders, under tests in each one.
- allTests is the CMake target to include them.  It if fine to just focuse on UTMain or one pertinent unit test while devloping, but please build and pass every thing before any big commits.

## pybind/ — do not modify (GUI/binding-owned)
The `pybind/` directory is the Python binding (nanobind C++ glue that compiles
`qchem_py.so`). It lives here only because it must compile against the C++20
modules; it is owned by the GUI/binding side, not the lib team. Lib-side work
MUST NOT edit `pybind/`. It links `libqchem` and rides the module build, so a lib
API change CAN break its compile — if that happens, FLAG it for the binding owner;
do not "fix" pybind/ yourself. (Built only under `-DQCHEM_PYBIND=ON`.)

## Modules & libraries

- There is a hierarchy: a library contains many modules, which contain many classes. The compiler
    enforces no cyclic dependencies between modules. The linker enforces no cyclic dependencies between
    libraries. Cyclic dependencies between tightly coupled classes within a single module are allowed.
    I may also refer to the libraries as project-modules, in order to distinguish them from c++ modules.
- Only re-export modules that don't have `.Internal.` in the module name; the goal is to avoid
    importing internals across library boundaries. Unit tests are allowed to cheat and import Internal stuff.
- The only header in the project is `forward.H`. It is used to forward-declare unit test classes that
    need access to internals of certain library classes. Header files for submodules are used where appropriate.

## Naming

- Use CamelCase / PascalCase for names, breaking to underscores at acronyms — for example
    `HF_Evaluator` or `HF_DFT_Evaluator`.
- Cosmetic: I like using `dcmplx = std::complex<double>;` (already in `src/Common/Types.C`) — it has the
    same number of chars as `double`, so things line up nicely in a class interface.

## Includes & types

- Prefer `import qchem.Math` over `#include <cmath>` — feel free to add symbols you need to `src/Common/Math.C`.
- Prefer `import qchem.Blaze` over `#include <blaze/Math.h>` — feel free to add symbols you need to
    `src/Common/Blaze.C`. There are a couple of exceptions, for example `std::sort` cannot see exported
    Blaze `op==` or `op!=` for iterators.
- Prefer `rvec_t` over `std::vector<double>`. We even have a push_back builder class (`VecBuilder<T>` in
    `src/Common/Blaze.C`) for `rvec_t` if you need to accumulate when the final length isn't known up front.
- Prefer `rmat_t` over `rvec_t<rvec_t>`. Sometimes there can be a very good reason for `rvec_t<rvec_t>`
    (push_back) — if so, let us discuss.

## Memory & ownership

- Raw `new` ops are fine if the pointer quickly goes into a `std::unique_ptr` or `std::shared_ptr`
    (within a few lines). As a result `delete` should be rare or non-existent.
- void* is banashed from this project. It has no place in modern c++.
- Do not key map/set etc off pointers.

## Design

- The codebase makes liberal use of abstract interface base classes with no data (default constructors).
    These are derived from using virtual inheritance. This design creates many diamond inheritance patterns,
    which are considered harmless in this project. These are not "diamonds of death". Editorial: there is a lot
    of misunderstanding around this, which is why multiple inheritance was incorrectly banished from Java and C#.
- Prefer classes to do/answer high-level operations/questions rather than expose internal data with lots of
    Get functions. Not a rule, just a "prefer" — `src/ChargeDensity/Internal/IrrepCD.C` is a great example: it
    has no `GetDensityMatrix()`. And from the other direction, the outer ChargeDensity classes have no `GetIrrepCD()` methods.
- There is a lot of `dynamic_cast` in the code. The intention is that we only cast between (almost) pure
    abstract base classes, never from an abstract base to a concrete implementation. Please flag violations as
    you work. I have a TODO item to do a system-wide survey of these casts and throw custom exceptions full of
    relevant information in the event they fail. Again, unit tests are allowed to cheat.

## Style

- Avoid
    ```
    for (...) {
        ...
    }
    ```
    Prefer one line, or
    ```
    for (...)
    {
        ...
    }
    ```
- Liberal use of asserts is encouraged — you know the reasons. And if you ever have to fire up the debugger,
    consider that this is usually an indication that there is a missing assert somewhere!
- Avoid the preprocessor as much as possible. An example of good use is getting info from the compile command
    line (`-DXXX`); as far as I know c++ does not have a modern way of doing this.

## Documentation

- I intend to document this project with Doxygen (if you think a different tool is better, let me know). If you
    can document interfaces as you go, that should result in better docs, since you have full context while
    writing the code.

## Heads-up

- FYI: I am constantly editing TODO, CLAUDE.md, and NOTES as you work. Just so you're not surprised when you
    do a `git status`.
