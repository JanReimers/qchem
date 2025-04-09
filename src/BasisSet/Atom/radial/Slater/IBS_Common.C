// File: Atom/radial/Slater/IBS_Common.H  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom Slater functions.

#include "Imp/BasisSet/Atom/radial/Slater/IBS_Common.H"
#include "Imp/BasisSet/Atom/radial/Slater/Integrals.H"
#include "Imp/Symmetry/YlQN.H"
#include "Imp/Symmetry/YlmQN.H"
#include <BasisFunction.H>
#include "oml/vector.h"


namespace Slater
{
//----------------------------------------------------------------
//
//  Common implementation for orbital and fit basis sets.
//
IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents,size_t L)
    : IBS_Common(new Yl_Sym(L))
    , AtomIrrepIEClient(exponents.size())
{
    Init(exponents,Norms(exponents,L),L);

};

IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents,size_t L, int m)
    : IBS_Common(new Ylm_Sym(L,m))
    , AtomIrrepIEClient(exponents.size())
{
    Init(exponents,Norms(exponents,L),L,m);
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
        os << "Spherical Slater L=" << GetSymmetry()
        << " with " << GetNumFunctions() << " basis functions, alpha={";
        for (auto b:*this) os << *b;
        os << "}" << std::endl;
    }
    return os;
}

} //namespace 