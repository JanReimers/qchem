// File: Atom/radial/BSpline/IBS_Common.H  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom BSpline functions.
#ifndef _Atom_radial_BSpline_IBS_Common_H_
#define _Atom_radial_BSpline_IBS_Common_H_

#include "IEC.H"
import qchem.BasisSet.qchem.BasisSet.IBS_Common;

namespace BSpline
{
    // Common implementation for orbital and fit basis sets.
template <size_t K> class IrrepBasisSet
    : public virtual ::IrrepBasisSet
    , public         IBS_Common
    , public         IrrepIEClient<K>
{
public:
    IrrepBasisSet(size_t Ngrid,double rmin, double rmax, Symmetry*,size_t L);
    IrrepBasisSet(size_t Ngrid,double rmin, double rmax, Symmetry*,size_t L, const std::vector<int>& ml);
    virtual std::ostream&  Write(std::ostream&    ) const;

};


}

#endif //_Atom_radial_BSpline_IBS_Common_H_