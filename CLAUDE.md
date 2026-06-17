# CLAUDE — Project Conventions

Brief notes about module/library conventions, naming, and includes.

- There is a hierarchy: a library contains many modules, which contain many classes. The compiler enforces no cyclic
    dependencies between modules. The linker enforces no cyclic dependencies between libraries. Cyclic dependencies
    between tightly coupled classes within a single module are allowed.
- Only re-export modules that don't have `.Internal.` in the module name; the goal is to avoid importing internals
    across library boundaries.
- The only header in the project is `forward.H`. It is used to forward-declare unit test classes that need access
    to internals of certain library classes. Header files for submodules are used where appropriate.
- Use CamelBack for names, but transition to acronyms with underscores, for example `HF_Evaluator` or
    `HF_DFT_Evaluator`.
- The codebase makes liberal use of abstract interface base classes with no data (default constructors). These are
    derived from using virtual inheritance. This design creates many diamond inheritance patterns, which are expected
    and considered harmless in this project.
- Prefer `import qchem.Math` over `#include <cmath>` — feel free to add symbols you need to `src/Common/Math.C`.
- Prefer `import qchem.Blaze` over `#include <blaze/Math.h>` — feel free to add symbols you need to
    `src/Common/Blaze.C`. There are a couple of exceptions, for example `std::sort` cannot see exported Blaze
    `op==` or `op!=` for iterators.