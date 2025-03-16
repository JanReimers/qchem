// File: SlaterBS.C  Spherical Slater basis set.

#include "Imp/BasisSet/Slater/IrrepBasisSet.H"
#include "Imp/BasisSet/Slater/BasisFunction.H"
#include "Imp/BasisSet/Slater/IntegralEngine.H"
#include "Imp/Symmetry/YlQN.H"
#include <BasisSet.H>
#include <iostream>
#include <cassert>

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
    IrrepIEClient::Init(exponents,L);
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new BasisFunction(e,L+1,L,ns(i++))); //ns from SlaterIEClient

};

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
    return new Fit_IBS(itsLAParams,db,es*2,0);
}

::Fit_IBS* Orbital_IBS::CreateVxcFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    return new Fit_IBS(itsLAParams,db,es*2.0/3.0,0);    
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

} //namespace
