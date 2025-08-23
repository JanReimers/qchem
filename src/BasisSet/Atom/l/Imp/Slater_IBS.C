// File: BasisSet/Atom/l/Imp/Slater_IBS.C  Slater Irrep Basis Set (IBS) with orbital angular momentum l.
module;
#include <vector>
module qchem.BasisSet.Atom.Internal.l.SlaterBS;
import BasisSet.Atom.Slater_IBS;
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
Orbital_IBS::Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L)
: ::Slater::IrrepBasisSet(exponents,new Yl_Sym(L),L)
    , Slater_IBS(exponents,L,{})
    , Atom::Orbital_HF_IBS <double>(db)
    , Atom::Orbital_IBS    <double>(db,this)
    , Atom::Orbital_DFT_IBS<double>(db,this)
{};

Orbital_IBS::Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L, const std::vector<int>& ml)
    : IrrepBasisSet(exponents,new Ylm_Sym(L,ml),L,ml)
    , Slater_IBS(exponents,L,ml)
    , Atom::Orbital_HF_IBS <double>(db)
    , Atom::Orbital_IBS    <double>(db,this)
    , Atom::Orbital_DFT_IBS<double>(db,this)
{};



::Fit_IBS* Orbital_IBS::CreateCDFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    return new Fit_IBS(db,AtomIrrepIEClient::es*2,0);
}

::Fit_IBS* Orbital_IBS::CreateVxcFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    return new Fit_IBS(db,AtomIrrepIEClient::es*2.0/3.0,0);    
}

//----------------------------------------------------------------
//
//  Fit with Slater_l  basis set.
//
Fit_IBS::Fit_IBS(const DB_cache<double>* db,const Vector<double>& exponents, size_t L)
    : ::Slater::IrrepBasisSet(exponents,new Yl_Sym(L),L)
    , Slater_IBS(exponents,L,{})
    , Atom::Fit_IBS(db,this)
    {};

}} //namespace
