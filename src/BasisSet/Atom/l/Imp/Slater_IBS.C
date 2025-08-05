// File: Atom/l/Slater_IBS.H  Slater Irrep Basis Set (IBS) with orbital angular momentum l.
module;
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
module qchem.BasisSet.Atom.Internal.l.SlaterBS;
import qchem.BasisSet.Atom.Internal.radial.SlaterBS;
import qchem.BasisSet.Atom.Internal.radial.Slater.Integrals;
import qchem.BasisSet;
import qchem.Symmetry.Yl;
import qchem.Symmetry.Ylm;

namespace Atoml
{
namespace Slater
{
//----------------------------------------------------------------
//
// Orbital SL basis set.
//
Orbital_IBS::Orbital_IBS(const DB_BS_2E<double>* db,const IE_Primatives* pie,const Vector<double>& exponents, size_t L)
: ::Slater::IrrepBasisSet(exponents,new Yl_Sym(L),L)
    , Orbital_IBS_Common<double>()
    , Orbital_HF_IBS_Common<double>(db)
    , IE_Common(db,pie)
    , AtomIE_DFT<double>(db,pie)
{
};

Orbital_IBS::Orbital_IBS(const DB_BS_2E<double>* db,const IE_Primatives* pie,const Vector<double>& exponents, size_t L, const std::vector<int>& ml)
    : IrrepBasisSet(exponents,new Ylm_Sym(L,ml),L,ml)
    , Orbital_IBS_Common<double>()
    , Orbital_HF_IBS_Common<double>(db)
    , IE_Common(db,pie)
    , AtomIE_DFT<double>(db,pie)
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

//----------------------------------------------------------------
//
//  Fit with Slater_l  basis set.
//
Fit_IBS::Fit_IBS(const DB_cache<double>* db,const IE_Primatives* pie,const Vector<double>& exponents, size_t L)
    : ::Slater::IrrepBasisSet(exponents,new Yl_Sym(L),L)
    , Fit_IE(db,pie)
    {
    };



}} //namespace
