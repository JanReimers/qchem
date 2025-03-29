// File: Atom/radial/Slater/IBS_Common.H  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom Slater functions.
#ifndef _Atom_radial_Slater_IBS_Common_H_
#define _Atom_radial_Slater_IBS_Common_H_

#include "../../IEC.H"
#include "../../../IBS_Common.H"

namespace Slater
{
    // Common implementation for orbital and fit basis sets.
class IrrepBasisSet
    : public virtual ::IrrepBasisSet
    , public         IBS_Common
    , public         AtomIrrepIEClient
{
public:
    IrrepBasisSet(const Vector<double>& exponents, size_t L);
    IrrepBasisSet(const Vector<double>& exponents, size_t L, int m);
    virtual std::ostream&  Write(std::ostream&    ) const;
private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;
};


}

#endif //_Atom_radial_Slater_IBS_Common_H_