// File: Atom/radial/Gaussian/IBS_Common.C  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom Gaussians.
module;
#include <iostream>
#include <iomanip>
#include <iosfwd>
#include <vector>
module qchem.BasisSet.Atom.Internal.radial.GaussianBS;
import qchem.BasisSet.Atom.Internal.radial.GaussianIntegrals;
import qchem.BasisFunction;
import qchem.Symmetry;
import oml;



namespace Gaussian
{
//----------------------------------------------------------------
//
//  Common implementation for orbital and fit basis sets.
//
IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents, Symmetry* sym,size_t l)
    : TIBS_Common1<double>(sym)
    , AtomIrrepIEClient(exponents.size())
{
    Init(exponents,Norms(exponents,l),l);
};

IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents, Symmetry* sym,size_t l, const std::vector<int>& ml)
    : TIBS_Common1<double>(sym)
    , AtomIrrepIEClient(exponents.size())
{
    Init(exponents,Norms(exponents,l),l,ml);
};

Vector<double> IrrepBasisSet::Norms(const Vector<double>& es, size_t l) const
{
    Vector<double> ns(es.size());
    int i=0;
    for (auto e:es) ns(++i)=::Gaussian::Norm(e,l);
    return ns;
}

std::ostream&  IrrepBasisSet::Write(std::ostream& os) const
{
    os << "Spherical Gaussian L=" << *GetSymmetry()
    << " with " << GetNumFunctions() << " basis functions, alpha={";
    // for (auto b:*this) os << *b;
    os << *front() << " ... " << *back();
    os << "}" << std::endl;
    return os;
}




} //namespace
