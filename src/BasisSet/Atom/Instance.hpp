// File: Instance.hpp  use preprocessor tricks to control instancing of BSpline basis sets.

// Usage:
//> #define INSTANCEk(k) template class MyTemplateClass<k>;
//> #include "Instance.hpp"


#define INSTANCE(seq) END(A seq)
#define A(k) INSTANCEk(k) B
#define B(k) INSTANCEk(k) A
#define A_END
#define B_END
#define END(...) END_(__VA_ARGS__)
#define END_(...) __VA_ARGS__##_END

// Instance with k={5,6,7}
INSTANCE((5)(6)(7)) 

#undef END
#undef END_
#undef B_END
#undef A_END
#undef A
#undef B
#undef INSTANCE
#undef INSTANCEk

