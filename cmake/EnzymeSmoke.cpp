// cmake/EnzymeSmoke.cpp  Configure's proof that the Enzyme plugin loads into THIS clang and differentiates: d(x^2)/dx at 3 is 6.
// Compiled and run by the root CMakeLists.txt with the same plugin flags the differentiating TUs use.
#include <cstdio>
#include <enzyme/enzyme>
static double Square(double x) {return x*x;}
int main()
{
    const double g=__enzyme_autodiff<double>(reinterpret_cast<void*>(Square), 3.0);
    std::printf("%.17g\n", g);
    return g==6.0 ? 0 : 1;
}
