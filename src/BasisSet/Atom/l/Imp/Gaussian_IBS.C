// File: Atom/l/Gaussian_IBS.H  Gaussian Irrep Basis Set (IBS) with orbital angular momentum l.
module;
#include <vector>
module qchem.BasisSet.Atom.Internal.l.GaussianBS;
import BasisSet.Atom.Gaussian_IBS;
import qchem.Symmetry.Yl;
import qchem.Symmetry.Ylm;

namespace Atoml
{
namespace Gaussian
{
  
Orbital_IBS::Orbital_IBS(const DB_BS_2E<double>* db,const IBS_Evaluator* eval,const Vector<double>& exponents, size_t L)
    : ::Gaussian::IrrepBasisSet(exponents,new Yl_Sym(L),L)
    , Atom::Orbital_HF_IBS <double>(db)
    , Atom::Orbital_IBS    <double>(db,eval)
    , Atom::Orbital_DFT_IBS<double>(db,eval)
{};

Orbital_IBS::Orbital_IBS(const DB_BS_2E<double>* db,const IBS_Evaluator* eval,const Vector<double>& exponents, size_t L, const std::vector<int>& ml)
    : ::Gaussian::IrrepBasisSet(exponents,new Ylm_Sym(L,ml),L,ml)
    , Atom::Orbital_HF_IBS <double>(db)
    , Atom::Orbital_IBS    <double>(db,eval)
    , Atom::Orbital_DFT_IBS<double>(db,eval)
{};

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

Fit_IBS::Fit_IBS(const DB_cache<double>* db,const Vector<double>& exponents, size_t L)
    : ::Gaussian::IrrepBasisSet(exponents,new Yl_Sym(L), L)
    , Gaussian_IBS(exponents,L,{})
    , Atom::Fit_IBS(db,this)
{};


} //namespace
} //namespace
