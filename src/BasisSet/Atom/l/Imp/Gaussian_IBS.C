// File: Atom/l/Gaussian_IBS.H  Gaussian Irrep Basis Set (IBS) with orbital angular momentum l.
module;
#include <iostream>
#include <cassert>
#include <vector>
module qchem.BasisSet.Atom.Internal.l.GaussianBS;
import qchem.BasisSet.Atom.Internal.radial.GaussianIntegrals;
import qchem.BasisSet.Atom.Internal.radial.GaussianBS;
import qchem.BasisSet;
import qchem.Symmetry.Yl;
import qchem.Symmetry.Ylm;

namespace Atoml
{
namespace Gaussian
{
  

//----------------------------------------------------------------
//
// Orbital SG basis set.
//
Orbital_IBS::Orbital_IBS(const DB_BS_2E<double>* db,const IE_Primatives* pie,const Vector<double>& exponents, size_t L)
    : ::Gaussian::IrrepBasisSet(exponents,new Yl_Sym(L),L)
    , Orbital_IBS_Common<double>()
    , Orbital_HF_IBS_Common<double>(db)
    , Orbital_IE(db,pie)
    {
    };

Orbital_IBS::Orbital_IBS(const DB_BS_2E<double>* db,const IE_Primatives* pie,const Vector<double>& exponents, size_t L, const std::vector<int>& ml)
    : ::Gaussian::IrrepBasisSet(exponents,new Ylm_Sym(L,ml),L,ml)
    , Orbital_IBS_Common<double>()
    , Orbital_HF_IBS_Common<double>(db)
    , Orbital_IE(db,pie)
{

};



::Fit_IBS* Orbital_IBS::CreateCDFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    auto pie=dynamic_cast<const IE_Primatives*>(bs);
    return new Fit_IBS(db,pie,es*2,0);
}

::Fit_IBS* Orbital_IBS::CreateVxcFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    auto pie=dynamic_cast<const IE_Primatives*>(bs);
    return new Fit_IBS(db,pie,es*2.0/3.0,0);    
}

Fit_IBS::Fit_IBS(const DB_cache<double>* db,const IE_Primatives* pie,const Vector<double>& exponents, size_t L)
: ::Gaussian::IrrepBasisSet(exponents,new Yl_Sym(L), L)
, Fit_IE(db,pie)
{
};


} //namespace
} //namespace
