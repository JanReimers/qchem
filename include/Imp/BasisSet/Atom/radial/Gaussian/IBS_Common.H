// File: Atom/radial/Gaussian/IBS_Common.H  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom Gaussians.
#ifndef _Atom_radial_Gaussian_IBS_Common_H_
#define _Atom_radial_Gaussian_IBS_Common_H_
#include <vector>

import qchem.BasisSet.qchem.BasisSet.IBS_Common;
import qchem.BasisSet.Atom.IEClient;


namespace Gaussian
{
    // Common implementation for orbital and fit basis sets.
class IrrepBasisSet
: public virtual ::IrrepBasisSet
, public         IBS_Common
, public         AtomIrrepIEClient
{
public:
    IrrepBasisSet(const Vector<double>& exponents, Symmetry*, size_t l);
    IrrepBasisSet(const Vector<double>& exponents, Symmetry*, size_t l, const std::vector<int>& ml);
    virtual std::ostream&  Write(std::ostream&    ) const;
private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;
};


} //namespace

#endif //_Atom_radial_Gaussian_IBS_Common_H_
