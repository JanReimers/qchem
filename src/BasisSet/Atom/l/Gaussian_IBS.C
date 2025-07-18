// File: Atom/l/Gaussian_IBS.H  Gaussian Irrep Basis Set (IBS) with orbital angular momentum l.

#include <iostream>
#include <cassert>
#include "l/Gaussian_IBS.H"
#include "l/Gaussian_BF.H"
#include "radial/Gaussian/Integrals.H"
#include <BasisSet/BasisSet.H>

import qchem.Symmetry.Yl;

namespace Atoml
{
namespace Gaussian
{
  

//----------------------------------------------------------------
//
// Orbital SG basis set.
//
Orbital_IBS::Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L)
    : ::Gaussian::IrrepBasisSet(exponents,new Yl_Sym(L),L)
    , Orbital_IBS_Common<double>()
    , Orbital_IE(db)
    {
        InsertBasisFunctions();   
    };

void Orbital_IBS::InsertBasisFunctions()
{
    size_t i=1;
    for (auto e:es) 
        IBS_Common::Insert(new BasisFunction(e,l,ns(i++)));
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
    std::cerr << "Why are you relocating a spherical Gaussian basis set?!" << std::endl;
    return 0;
}

Fit_IBS::Fit_IBS(const DB_cache<double>* db,const Vector<double>& exponents, size_t L)
: ::Gaussian::IrrepBasisSet(exponents,new Yl_Sym(L), L)
, TIBS_Common<double>()
, Fit_IE(db)
{
    InsertBasisFunctions();   
};

void Fit_IBS::InsertBasisFunctions()
{
    size_t i=1;
    for (auto e:es) 
        IBS_Common::Insert(new BasisFunction(e,l,ns(i++)));
}

::Fit_IBS* Fit_IBS::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a spherical Gaussian basis set?!" << std::endl;
    return 0;
}


} //namespace
} //namespace
