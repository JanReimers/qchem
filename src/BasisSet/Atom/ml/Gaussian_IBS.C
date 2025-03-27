// File: Atom/ml/Gaussian_IBS.H  r^l exp(-ar^2)*Y_lm type Irrep Basis set (IBS).

#include "Imp/BasisSet/Atom/ml/Gaussian_IBS.H"
#include "Imp/BasisSet/Atom/ml/Gaussian_BF.H"
#include "Imp/Symmetry/YlmQN.H"
#include <iostream>
#include <cassert>

namespace Atom_ml
{
namespace Gaussian
{

IrrepBasisSet::IrrepBasisSet(const std::set<double>& exponents, size_t L, int m)
    : IrrepBasisSetCommon(new YlmQN(L,m))
    , IrrepIEClient(exponents.size())
{
    IrrepIEClient::Init(exponents,L,m);
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new BasisFunction(e,L+1,L,m,ns(i++))); //ns from SlaterIEClient

};

IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents,size_t L, int m)
    : IrrepBasisSetCommon(new YlmQN(L,m))
    , IrrepIEClient(exponents.size())
{
    IrrepIEClient::Init(exponents,L,m);
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new BasisFunction(e,L+1,L,m,ns(i++))); //ns from SlaterIEClient

};
std::ostream&  IrrepBasisSet::Write(std::ostream& os) const
{
    if (Pretty())
    {
        os << "Gaussian functions l,m=" << GetQuantumNumber()
        << " with " << GetNumFunctions() << " basis functions, alpha={";
        for (auto b:*this) os << *b;
        os << "}" << std::endl;
    }
    return os;
}
::Fit_IBS* Orbital_IBS::CreateCDFitBasisSet(const Cluster*) const
{
    // return new IrrepBasisSet(itsLAParams,GetDataBase(),0,es*2,0,0);
    assert(false);
    return 0;
}

::Fit_IBS* Orbital_IBS::CreateVxcFitBasisSet(const Cluster*) const
{
    // return new IrrepBasisSet(itsLAParams,GetDataBase(),0,es*2.0/3.0,0,0);    
    assert(false);
    return 0;
}

::IrrepBasisSet* Orbital_IBS::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a Gaussian atomic basis set?!" << std::endl;
    return 0;
}


} //namespace
} //namespace
