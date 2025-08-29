// File: Atom/kappa/Gaussian_IBS.C  Restricted Kinetic Balance (RKB) Irrep Basis Set (IBS).
module;
#include <memory>
#include <iostream>
#include <cassert>
#include <cmath>
module qchem.BasisSet.Atom.Internal.kappa.GaussianBS;
import BasisSet.Atom.Gaussian_IBS;
import qchem.Symmetry.Okmj;

using std::endl;

namespace Atom_kappa
{
namespace Gaussian
{

Orbital_RKB_IBS::Orbital_RKB_IBS
    (const DB_cache<double>* db, const ds_t& exponents, int kappa)
    : IrrepBasisSet_Common<double>(new Omega_k_Sym(kappa))
    , Orbital_RKB_IBS_Common<double>
        (db, kappa
        , new Orbital_RKBL_IBS<double>(db,exponents,kappa)
        , new Orbital_RKBS_IBS<double>(db,exponents,kappa) //Known memory leak, need redesign
        )
{
    
};

std::ostream&  Orbital_RKB_IBS::Write(std::ostream& os) const
{
    os << "Dirac basis set." << endl << "    Large: " << *itsRKBL << endl << "    Small: " << *itsRKBS << endl;
    return os;
}

//-----------------------------------------------------------------------------------------------
//
//  Large sector
//
template <class T> Orbital_RKBL_IBS<T>::Orbital_RKBL_IBS
(const DB_cache<T>* db,const ds_t& exponents,int kappa)
    : Gaussian_IBS(exponents,Omega_k_Sym::l(kappa))
    , Atom::IrrepBasisSet(this,new Omega_k_Sym(kappa))
    , Atom::Orbital_RKBL_IBS<T>(db,this,kappa)
{};



//-----------------------------------------------------------------------------------------------
//
//  Small sector
//
template <class T> Orbital_RKBS_IBS<T>::Orbital_RKBS_IBS(const DB_cache<T>* db, const ds_t& exponents,int kappa)
    : Gaussian_RKBS_IBS(exponents,kappa,Omega_k_Sym::l(kappa))
    , Atom::IrrepBasisSet(this,new Omega_k_Sym(-kappa))
    , Atom::Orbital_RKBS_IBS<T>(db,this,kappa)
{
    
}

template <class T> std::ostream&  Orbital_RKBS_IBS<T>::Write(std::ostream& os) const
{
    os << "Gaussian     " << this->GetSymmetry()
    << "               r^" << Gaussian_IBS::l << "*exp(-e*r^2), e={";
    for (auto e:Gaussian_IBS::es) os << e << ",";
    os << "}";
    return os;
}



}} //namespace
