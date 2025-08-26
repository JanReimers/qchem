// File: Atom/radial/BSpline/IBS_Common.H  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom BSpline functions.
module;
#include <iosfwd>
#include <iostream>
#include <vector>
#include <cassert>
#include "bspline/operators/Derivative.h"
module qchem.Basisset.Atom.radial.BSpline.BS_Common; 
import qchem.Symmetry;


namespace BSpline
{
//----------------------------------------------------------------
//
//  Common implementation for orbital and fit basis sets.
//

// template <size_t K> IrrepBasisSet<K>::IrrepBasisSet(size_t Ngrid, double rmin, double rmax, Symmetry* sym, size_t L)
//     : IrrepBasisSet_Common<double>(sym)
//     , BSpline::IrrepIEClient<K>(Ngrid,rmin, rmax,L)
// {};
// template <size_t K> IrrepBasisSet<K>::IrrepBasisSet(size_t Ngrid, double rmin, double rmax, Symmetry* sym, size_t L, const std::vector<int>& ml)
//     : IrrepBasisSet_Common<double>(sym)
//     , BSpline::IrrepIEClient<K>(Ngrid,rmin, rmax,L,ml)
// {};

// template <size_t K> IrrepBasisSet<K>::Vec     IrrepBasisSet<K>::operator() (const RVec3& r) const
// {
//     Vec ret(size());
//     double mr=norm(r);
//     size_t i=0;
//     for (auto s:IEC::splines) 
//     {
//         ++i;
//         ret(i)=ns(i)*s(mr);
//     }
//     return ret;
// }
// template <size_t K> IrrepBasisSet<K>::Vec3Vec IrrepBasisSet<K>::Gradient   (const RVec3& r) const
// {
//     Vec3Vec ret(size());
//     double mr=norm(r);
//     if (mr==0.0) 
//     {
        
//         Fill(ret,RVec3(0,0,0));
//         return ret; //Cusp at the origin so grad is undefined.
//     }
//     assert(mr>0);
//     Fill(ret,r/mr);
//     size_t i=0;
//     for (auto s:splines) 
//     {
//         auto dsdx=transformSpline(bspline::operators::Dx<1>{},s);
//         ++i;
//         ret(i)*=ns(i)*dsdx(mr);
//     }
//     return ret;
// }

template <size_t K> std::ostream&  IrrepBasisSet<K>::Write(std::ostream& os) const
{
    os << "Spherical BSpline L=" << *GetSymmetry();
    itsEval->Write(os);
    return os;
}

#define INSTANCEk(k) template class IrrepBasisSet<k>;
#include "../Instance.hpp"

} //namespace 