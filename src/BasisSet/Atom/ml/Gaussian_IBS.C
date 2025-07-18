// File: Atom/ml/Gaussian_IBS.H  r^l exp(-ar^2)*Y_lm type Irrep Basis set (IBS).

#include <iostream>
#include "ml/Gaussian_IBS.H"
#include "ml/Gaussian_BF.H"

#include <cassert>

import qchem.Symmetry.Ylm;

namespace Atom_ml
{
namespace Gaussian
{

Orbital_IBS::Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L, const std::vector<int>& ml)
    : ::Gaussian::IrrepBasisSet(exponents,new Ylm_Sym(L,ml),L,ml)
    , Orbital_IBS_Common<double>()
    , Atoml::Gaussian::Orbital_IE(db)
{
    InsertBasisFunctions();   

};

void Orbital_IBS::InsertBasisFunctions()
{
    size_t i=1;
    for (auto e:es) 
        IBS_Common::Insert(new BasisFunction(e,l+1,l,ml[0],ns(i++))); 
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

::IrrepBasisSet* Orbital_IBS::Clone(const ::IrrepBasisSet::RVec3&) const
{
    std::cerr << "Why are you relocating a Gaussian atomic basis set?!" << std::endl;
    return 0;
}



} //namespace
} //namespace
