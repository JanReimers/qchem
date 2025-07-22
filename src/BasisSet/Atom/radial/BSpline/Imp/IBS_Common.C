// File: Atom/radial/BSpline/IBS_Common.H  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom BSpline functions.
module;
#include <iosfwd>
#include <iostream>
#include <vector>
module qchem.Basisset.Atom.radial.BSpline.BS_Common; 
import qchem.BasisFunction;
import qchem.Symmetry;


namespace BSpline
{
//----------------------------------------------------------------
//
//  Common implementation for orbital and fit basis sets.
//

template <size_t K> IrrepBasisSet<K>::IrrepBasisSet(size_t Ngrid, double rmin, double rmax, Symmetry* sym, size_t L)
    : IBS_Common(sym)
    , BSpline::IrrepIEClient<K>(Ngrid,rmin, rmax,L)
{};
template <size_t K> IrrepBasisSet<K>::IrrepBasisSet(size_t Ngrid, double rmin, double rmax, Symmetry* sym, size_t L, const std::vector<int>& ml)
    : IBS_Common(sym)
    , BSpline::IrrepIEClient<K>(Ngrid,rmin, rmax,L,ml)
{};

template <size_t K> std::ostream&  IrrepBasisSet<K>::Write(std::ostream& os) const
{
    os << "Spherical BSpline L=" << *GetSymmetry()
    << " with " << GetNumFunctions() << " basis functions, {";
    // for (auto b:*this) os << *b;
    os << *front() << " ... " << *back();
    os << "}" << std::endl;
    return os;
}

#define INSTANCEk(k) template class IrrepBasisSet<k>;
#include "../../../Instance.hpp"

} //namespace 