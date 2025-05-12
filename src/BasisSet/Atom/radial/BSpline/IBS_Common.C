// File: Atom/radial/BSpline/IBS_Common.H  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom BSpline functions.

#include "Imp/BasisSet/Atom/radial/BSpline/IBS_Common.H"
//#include "Imp/BasisSet/Atom/radial/BSpline/Integrals.H"
#include <BasisFunction.H>
#include <Symmetry.H>
//#include "oml/vector.h"


namespace BSpline
{
//----------------------------------------------------------------
//
//  Common implementation for orbital and fit basis sets.
//

template <size_t K> IrrepBasisSet<K>::IrrepBasisSet(size_t Ngrid, double rmin, double rmax, Symmetry* sym, size_t L, int m)
    : IBS_Common(sym)
    , BSpline::IrrepIEClient<K>(Ngrid,rmin, rmax,L,m)
{
    
    
};

template <size_t K> std::ostream&  IrrepBasisSet<K>::Write(std::ostream& os) const
{
    if (Pretty())
    {
        os << "Spherical BSpline L=" << GetSymmetry()
        << " with " << GetNumFunctions() << " basis functions, {";
        for (auto b:*this) os << *b;
        os << "}" << std::endl;
    }
    return os;
}

template class IrrepBasisSet<6>;

} //namespace 