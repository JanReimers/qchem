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
- Miss understanding: " you've told me never to call ninja directly (cmake handles it)" I only recommended cmake because you were
  struggling to get the ninja to build.  You had to iterate by touching certain files to trigger rebuilds.  It looked inefficient at the time.  Any if ninja is working well and prefer, then go ahead.
- Prefer classes to do/answer high level operations/questions rather expose internal data with lots of Get functions.  Not a rule just a "prefer" ... src/ChargeDensity/Internal/IrrepCD.C is a great example it has no GetDensityMatrix().  And from other direction, the outer ChargeDensity classes have no GetIrrepCD() methods. 
- Raw new ops are fine if the pointer quickly (withint a few lines) goes into a std::unique_ptr or std::shared_ptr.  As a result delete should be rare or non-existent.  
- I like code consistency: So for example if we use the Cache4 looping mechanism for HF performance and the Cache3 looping mechanism for DFT does not yield any great performance improvement, we should still use Cache3 for DFT just so the code looks the same.
- Avoid the preprocessor as much as possible.  An example of good use is getting info from the compile command line -DXXX, as far as I know c++ does not have a modern way of doing this.   

