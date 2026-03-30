// File: Atom/l/Gaussian_IBS.H  Gaussian Irrep Basis Set (IBS) with orbital angular momentum l.
module;
#include <blaze/math/DynamicVector.h>

module qchem.BasisSet.Atom.Gaussian.NR.BS;
import qchem.Symmetry.Yl;
import qchem.Symmetry.Ylm;

namespace AtomBS
{
namespace Gaussian
{
  
Orbital_IBS::Orbital_IBS(const DB_BS_2E<double>* db,const rvec_t& exponents, size_t L)
    : Gaussian_IBS(exponents,L)
    , AtomBS::IrrepBasisSet(this,new Yl_Sym(L))
    , AtomBS::Orbital_HF_IBS <double>(db)
    , AtomBS::Orbital_IBS    <double>(db,this)
    , AtomBS::Orbital_DFT_IBS<double>(db,this)
{};

Orbital_IBS::Orbital_IBS(const DB_BS_2E<double>* db,const rvec_t& exponents, size_t L, const std::vector<int>& ml)
    : Gaussian_IBS(exponents,L,ml)
    , AtomBS::IrrepBasisSet(this,new Ylm_Sym(L,ml))
    , AtomBS::Orbital_HF_IBS <double>(db)
    , AtomBS::Orbital_IBS    <double>(db,this)
    , AtomBS::Orbital_DFT_IBS<double>(db,this)
{};

::Fit_IBS* Orbital_IBS::CreateCDFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    return new Fit_IBS(db,Gaussian_IBS::es*2.0,0);
}

::Fit_IBS* Orbital_IBS::CreateVxcFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    return new Fit_IBS(db,Gaussian_IBS::es*2.0/3.0,0);    
}

Fit_IBS::Fit_IBS(const DB_cache<double>* db,const rvec_t& exponents, size_t L)
    : Gaussian_IBS(exponents,L)
    , AtomBS::IrrepBasisSet(this,new Yl_Sym(L))
    , AtomBS::Fit_IBS(db,this)
    {};


} //namespace
} //namespace
