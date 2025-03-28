// File: Atom/l/Slater_IBS.H  Slater Irrep Basis Set (IBS) with orbital angular momentum l.

#include "Imp/BasisSet/Atom/l/Slater_IBS.H"
#include "Imp/BasisSet/Atom/l/Slater_IE.H"
#include "Imp/BasisSet/Atom/l/Slater_BF.H"
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
// Orbital SL basis set.
//
void Orbital_IBS::InsertBasisFunctions()
{
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new BasisFunction(e,l+1,l,ns(i++))); //ns from SlaterIEClient
}

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
void Fit_IBS::InsertBasisFunctions()
{
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new BasisFunction(e,l+1,l,ns(i++))); //ns from SlaterIEClient
}

::Fit_IBS* Fit_IBS::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a spherical Slater basis set?!" << std::endl;
    return 0;
}

}} //namespace
