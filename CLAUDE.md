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
- Prefer rvec_t over std::vector<double>.  We even have a push_back builder class for rvec_t if you need that.
- Prefer rmat_t over rvec_t<rvec_t> Sometimes there can be a very good reason for rvec_t<rvec_t> (push_back) if so let us discuss.
- Avoid
    for (...) {
        ...
    }
    Prefer one line or
    for (...)
    {
        ...
    }

- Liberal use of asserts is encouraged, you know the reasons.  And if you ever have to fire up the debugger consider that
  this is usually and indication that there is a missing assert somewhere!
- Miss understanding: " you've told me never to call ninja directly (cmake handles it)" I only recommended cmake because you were
  struggling to get the ninja to build.  You had to iterate by touching certain files to trigger rebuilds.  It looked inefficient at the time.  Any if ninja is working well and prefer, then go ahead.
- Prefer classes to do/answer high level operations/questions rather expose internal data with lots of Get functions.  Not a rule just a "prefer" ... src/ChargeDensity/Internal/IrrepCD.C is a great example it has no GetDensityMatrix().  And from other direction, the outer ChargeDensity classes have no GetIrrepCD() methods.   
- Raw new ops are fine if the pointer quickly  goes into a std::unique_ptr or std::shared_ptr (withint a few lines).  As a result delete should be rare or non-existent.  
- Avoid the preprocessor as much as possible.  An example of good use is getting info from the compile command line -DXXX, as far as I know c++ does not have a modern way of doing this.
- Cosmetic: I like using dcmplx=std::complex<double>; (already in src/Common/Types.C)... it has the same number of chars as double so things line up nicely in class interface.
- There is a lot of dynamic_cast in the code ... the intention is that we are only casting between (almost) pure abstract base classes, never from abstract base to concrete implementation.   Please flag violations as you work.  I have a TODO item to do a system wide survey of these casts and throw custom exceptions full of relevent information in the even they fail.  
- FYI: I am constantly editing TODO, Claude.md and NOTES as you work.

