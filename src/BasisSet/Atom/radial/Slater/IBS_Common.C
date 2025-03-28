// File: Atom/radial/Slater/IBS_Common.H  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom Slater functions.

#include "Imp/BasisSet/Atom/l/Slater_BF.H"
#include "Imp/BasisSet/Atom/ml/Slater_BF.H"
#include "Imp/BasisSet/Atom/radial/Slater/IBS_Common.H"
#include "Imp/BasisSet/Atom/radial/Slater/Integrals.H"
#include "Imp/Symmetry/YlQN.H"
#include "Imp/Symmetry/YlmQN.H"
#include <BasisFunction.H>

namespace Slater
{
//----------------------------------------------------------------
//
//  Common implementation for orbital and fit basis sets.
//
IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents,size_t L)
    : IrrepBasisSetCommon(new YlQN(L))
    , AtomIrrepIEClient(exponents.size())
{
    Init(exponents,Norms(exponents,L),L);
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new Atoml::Slater::BasisFunction(e,L+1,L,ns(i++))); //ns from SlaterIEClient

};

IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents,size_t L, int m)
    : IrrepBasisSetCommon(new YlmQN(L,m))
    , AtomIrrepIEClient(exponents.size())
{
    Init(exponents,Norms(exponents,L),L,m);
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new Atom_ml::Slater::BasisFunction(e,L+1,L,m,ns(i++))); //ns from SlaterIEClient

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
    if (Pretty())
    {
        os << "Spherical Slater L=" << GetQuantumNumber()
        << " with " << GetNumFunctions() << " basis functions, alpha={";
        for (auto b:*this) os << *b;
        os << "}" << std::endl;
    }
    return os;
}

} //namespace 