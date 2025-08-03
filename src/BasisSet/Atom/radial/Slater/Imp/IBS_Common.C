// File: Atom/radial/Slater/IBS_Common.H  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom Slater functions.
module;
#include <iostream>
#include <cmath>
#include <vector>
module qchem.BasisSet.Atom.Internal.radial.SlaterBS;
import qchem.BasisSet.Atom.Internal.radial.Slater.Integrals;
import qchem.BasisFunction;
import qchem.Symmetry;
import oml;

namespace Slater
{
//----------------------------------------------------------------
//
//  Common implementation for orbital and fit basis sets.
//
IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents,Symmetry* sym,size_t L)
    : TIBS_Common1<double>(sym)
    , AtomIrrepIEClient(exponents.size())
{
    Init(exponents,Norms(exponents,L),L);

};

IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents,Symmetry* sym, size_t L, const std::vector<int>& ml)
    : TIBS_Common1<double>(sym)
    , AtomIrrepIEClient(exponents.size())
{
    Init(exponents,Norms(exponents,L),L,ml);
};

Vector<double> IrrepBasisSet::Norms(const Vector<double>& es, size_t l) const
{
    Vector<double> ns(es.size());
    int i=0;
    for (auto e:es) ns(++i)=::Slater::Norm(e,l+1);
    return ns;
}

std::ostream&  IrrepBasisSet::Write(std::ostream& os) const
{
    os << "Spherical Slater L=" << *GetSymmetry()
    << " with " << GetNumFunctions() << " basis functions, alpha={";
    // for (auto b:*this) os << *b;
    os << *front() << " ... " << *back();
    os << "}" << std::endl;
    return os;
}

} //namespace 