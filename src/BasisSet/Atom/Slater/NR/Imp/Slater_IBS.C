// File: BasisSet/Atom/l/Imp/Slater_IBS.C  Slater Irrep Basis Set (IBS) with orbital angular momentum l.
module;
#include <vector>
#include <valarray>
module qchem.BasisSet.Atom.Slater.NR.BS;
import qchem.Symmetry.Yl;
import qchem.Symmetry.Ylm;

namespace AtomBS
{
namespace Slater
{
//----------------------------------------------------------------
//
// Orbital SL basis set.
//
Orbital_IBS::Orbital_IBS(const DB_BS_2E<double>* db,const ds_t& exponents, size_t L)
    : Slater_IBS(exponents,L)
    , AtomBS::IrrepBasisSet(this,new Yl_Sym(L))
    , AtomBS::Orbital_HF_IBS <double>(db)
    , AtomBS::Orbital_IBS    <double>(db,this)
    , AtomBS::Orbital_DFT_IBS<double>(db,this)
{};

Orbital_IBS::Orbital_IBS(const DB_BS_2E<double>* db,const ds_t& exponents, size_t L, const std::vector<int>& ml)
    : Slater_IBS(exponents,L,ml)
    , AtomBS::IrrepBasisSet(this,new Ylm_Sym(L,ml))
    , AtomBS::Orbital_HF_IBS <double>(db)
    , AtomBS::Orbital_IBS    <double>(db,this)
    , AtomBS::Orbital_DFT_IBS<double>(db,this)
{};



::Fit_IBS* Orbital_IBS::CreateCDFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    return new Fit_IBS(db,Slater_IBS::es*2.0,0);
}

::Fit_IBS* Orbital_IBS::CreateVxcFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    return new Fit_IBS(db,Slater_IBS::es*2.0/3.0,0);    
}

//----------------------------------------------------------------
//
//  Fit with Slater_l  basis set.
//
Fit_IBS::Fit_IBS(const DB_cache<double>* db,const ds_t& exponents, size_t L)
    : Slater_IBS(exponents,L)
    , AtomBS::IrrepBasisSet(this,new Yl_Sym(L))
    , AtomBS::Fit_IBS(db,this)
    {};

}} //namespace
