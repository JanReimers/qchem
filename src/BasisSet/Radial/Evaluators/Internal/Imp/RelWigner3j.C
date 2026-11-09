// File: BasisSet/Radial/Evaluators/Internal/Imp/RelWigner3j.C
module;
#include <cstdlib>   // std::abs(int)
#include <cassert>
#include <stdexcept>
#include <string>
module qchem.BasisSet.Radial.Evaluators.Internal.RelWigner3j;
import qchem.BasisSet.Radial.Evaluators.Internal.Wigner3j; //Wigner::wigner3j (home-grown half-integer 3j)
import qchem.Symmetry.Atom.Spherical;
import qchem.Factorials;                          //SphericalSpinor::j(κ)

namespace qchem {

const RelWigner3j RelWigner3j::w3j{};

void RelWigner3j::Check(int κa, int κb, int k)
{
    // The Racah core indexes the Fact table up to n = ja+k+jb+1 = (2*KMax-1)+KkMax+1 at the bounds below, so
    // these bounds ARE the Fact table's.  Checked at this (non-hot) entry and THROWN, never assert-only.
    static_assert(2*KMax-1+KkMax+1 < int(qchem::Fact.size()), "RelWigner3j bounds outgrew the Fact table");
    if (κa==0 || std::abs(κa)>KMax || κb==0 || std::abs(κb)>KMax || k<0 || k>KkMax)
        throw std::out_of_range("RelWigner3j: (κa,κb,k)=("+std::to_string(κa)+","+std::to_string(κb)+","+std::to_string(k)
                               +") outside 0<|κ|<="+std::to_string(KMax)+", 0<=k<="+std::to_string(KkMax));
}

double RelWigner3j::operator()(int κa, int κb, int k) const
{
    Check(κa,κb,k);
    double ja=Symmetry::Atom::SphericalSpinor::j(κa), jb=Symmetry::Atom::SphericalSpinor::j(κb);
    return Wigner::wigner3j(ja, (double)k, jb, 0.5, 0.0, -0.5);
}

double RelWigner3j::operator()(int κa, int κb, int k, double mja, double mjb) const
{
    Check(κa,κb,k);
    double ja=Symmetry::Atom::SphericalSpinor::j(κa), jb=Symmetry::Atom::SphericalSpinor::j(κb);
    assert(mja>=-ja && mja<=ja);
    assert(mjb>=-jb && mjb<=jb);
    return Wigner::wigner3j(ja, (double)k, jb, -mja, mja-mjb, mjb);
}

} // namespace qchem