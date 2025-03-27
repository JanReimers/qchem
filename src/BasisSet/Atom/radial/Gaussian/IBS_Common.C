// File: Atom/radial/Gaussian/IBS_Common.C  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom Gaussians.

#include "Imp/BasisSet/Atom/radial/Gaussian/IBS_Common.H"
#include "Imp/BasisSet/Atom/radial/Gaussian/BasisFunction.H"
#include "Imp/BasisSet/Atom/radial/Gaussian/Integrals.H"
#include "Imp/Symmetry/YlQN.H"
#include "Imp/Symmetry/YlmQN.H"


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
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new ::Gaussian::BasisFunction(e,l,ns(i++)));
};

IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents,size_t l, int m)
    : IrrepBasisSetCommon(new YlmQN(l,m))
    , AtomIrrepIEClient(exponents.size())
{
    Init(exponents,Norms(exponents,l),l,m);
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new BasisFunction_ml(e,l+1,l,m,ns(i++))); 
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

::Fit_IBS* Fit_IBS::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a spherical Gaussian basis set?!" << std::endl;
    return 0;
}


} //namespace
