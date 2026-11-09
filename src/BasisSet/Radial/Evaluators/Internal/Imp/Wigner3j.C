// File: src/BasisSet/Radial/Evaluators/Internal/Imp/Wigner3j.C
module;
#include <cmath>      // std::sqrt (runtime reconstruction of the signed value)
#include <cassert>
#include <stdexcept>
#include <string>
module qchem.BasisSet.Radial.Evaluators.Internal.Wigner3j;

namespace qchem {

using namespace w3j_detail;

const Wigner3j Wigner3j::w3j{};

void Wigner3j::CheckL(int l)
{
    if (l<0 || l>LMax)
        throw std::out_of_range("Wigner3j: l="+std::to_string(l)+" exceeds the 3j table (LMax="+std::to_string(LMax)+")");
}

double Wigner3j::operator()(int la, int lb, int k, int ma, int mb) const
{
    assert(la>=0); assert(la<=LMax);
    assert(lb>=0); assert(lb<=LMax);
    assert(k >=0); assert(k <=2*LMax);
    assert(ma>=-la); assert(ma<=la);
    assert(mb>=-lb); assert(mb<=lb);
    const W& w=Table[Index(la,lb,k,ma,mb)];
    return w.sgn*std::sqrt(w.sq);
}

double Wigner3j::sq(int la, int lb, int k, int ma, int mb) const
{
    assert(la>=0); assert(la<=LMax);
    assert(lb>=0); assert(lb<=LMax);
    assert(k >=0); assert(k <=2*LMax);
    assert(ma>=-la); assert(ma<=la);
    assert(mb>=-lb); assert(mb<=lb);
    return Table[Index(la,lb,k,ma,mb)].sq;
}

} // namespace qchem