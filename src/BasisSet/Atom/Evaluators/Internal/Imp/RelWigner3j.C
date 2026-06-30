// File: BasisSet/Atom/Evaluators/Internal/Imp/RelWigner3j.C
module;
#include <cstdlib>   // std::abs(int)
#include <cassert>
module qchem.BasisSet.Atom.Evaluators.Internal.RelWigner3j;
import qchem.BasisSet.Atom.Evaluators.Internal.Wigner3j; //Wigner::wigner3j (home-grown half-integer 3j)
import qchem.Symmetry.Spherical;                          //SphericalSpinor::j(κ)

namespace qchem {

const RelWigner3j RelWigner3j::w3j{};

double RelWigner3j::operator()(int κa, int κb, int k) const
{
    assert(κa!=0 && std::abs(κa)<=KMax);
    assert(κb!=0 && std::abs(κb)<=KMax);
    assert(k>=0 && k<=KkMax);
    double ja=Symmetry::SphericalSpinor::j(κa), jb=Symmetry::SphericalSpinor::j(κb);
    return Wigner::wigner3j(ja, (double)k, jb, 0.5, 0.0, -0.5);
}

double RelWigner3j::operator()(int κa, int κb, int k, double mja, double mjb) const
{
    assert(κa!=0 && std::abs(κa)<=KMax);
    assert(κb!=0 && std::abs(κb)<=KMax);
    assert(k>=0 && k<=KkMax);
    double ja=Symmetry::SphericalSpinor::j(κa), jb=Symmetry::SphericalSpinor::j(κb);
    assert(mja>=-ja && mja<=ja);
    assert(mjb>=-jb && mjb<=jb);
    return Wigner::wigner3j(ja, (double)k, jb, -mja, mja-mjb, mjb);
}

} // namespace qchem