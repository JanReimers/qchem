// File: Atom/l/Slater_IBS.H  Slater Irrep Basis Set (IBS) with orbital angular momentum l.

#include "Imp/BasisSet/Atom/l/Slater_IBS.H"
#include "Imp/BasisSet/Atom/l/Slater_IE.H"
#include "Imp/BasisSet/Atom/radial/Slater/BasisFunction.H"
#include "Imp/BasisSet/Atom/radial/Slater/Integrals.H"
#include "Imp/Symmetry/YlQN.H"
#include <BasisSet.H>
#include <iostream>
#include <cassert>

namespace Atoml
{
namespace Slater
{
//----------------------------------------------------------------
//
//  Common implementation for orbital and fit basis sets.
//
IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents,size_t L)
    : IrrepBasisSetCommon(new YlQN(L))
    , IrrepIEClient(exponents.size())
{
    IrrepIEClient::Init(exponents,Norms(exponents,L),L);
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new ::Slater::BasisFunction(e,L+1,L,ns(i++))); //ns from SlaterIEClient

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

//----------------------------------------------------------------
//
// Orbital SL basis set.
//

::Fit_IBS* Orbital_IBS::CreateCDFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    return new Fit_IBS(db,es*2,0);
}

::Fit_IBS* Orbital_IBS::CreateVxcFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    return new Fit_IBS(db,es*2.0/3.0,0);    
}

::IrrepBasisSet* Orbital_IBS::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a spherical Slater basis set?!" << std::endl;
    return 0;
}
//----------------------------------------------------------------
//
//  Fit PG basis set.
//
::Fit_IBS* Fit_IBS::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a spherical Slater basis set?!" << std::endl;
    return 0;
}

}} //namespace
