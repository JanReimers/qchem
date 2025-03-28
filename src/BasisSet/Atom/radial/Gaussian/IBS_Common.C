// File: Atom/radial/Gaussian/IBS_Common.C  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom Gaussians.

#include "Imp/BasisSet/Atom/radial/Gaussian/IBS_Common.H"
#include "Imp/BasisSet/Atom/radial/Gaussian/Integrals.H"
#include "Imp/Symmetry/YlQN.H"
#include "Imp/Symmetry/YlmQN.H"
#include <BasisFunction.H>


namespace Gaussian
{
//----------------------------------------------------------------
//
//  Common implementation for orbital and fit basis sets.
//
IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents,size_t l)
    : IrrepBasisSetCommon(new YlQN(l))
    , AtomIrrepIEClient(exponents.size())
{
    Init(exponents,Norms(exponents,l),l);
};

IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents,size_t l, int m)
    : IrrepBasisSetCommon(new YlmQN(l,m))
    , AtomIrrepIEClient(exponents.size())
{
    Init(exponents,Norms(exponents,l),l,m);
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
    if (Pretty())
    {
        os << "Spherical Gaussian L=" << GetQuantumNumber()
        << " with " << GetNumFunctions() << " basis functions, alpha={";
        for (auto b:*this) os << *b;
        os << "}" << std::endl;
    }
    return os;
}




} //namespace
