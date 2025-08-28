// File: Atom/kappa/Slater_IBS.C  Slater Irrep Basis Set (IBS) with Restricted Kinetic Balance (RKB).
module;
#include <iostream>
#include <cassert>
#include <iomanip>
#include <memory>
module qchem.BasisSet.Atom.Internal.kappa.SlaterBS;
import BasisSet.Atom.Slater_IBS;
import qchem.Symmetry.Okmj;

namespace Atom_kappa
{
namespace Slater
{

//
//  Concrete  Slater basis set.
//

Orbital_RKB_IBS::Orbital_RKB_IBS
    (const DB_cache<double>* db,const Vector<double>& exponents, int kappa)
    : IrrepBasisSet_Common<double>(new Omega_k_Sym(kappa))
    , Orbital_RKB_IBS_Common<double>(db, kappa
        , new Orbital_RKBL_IBS<double>(db,exponents, kappa)
        , new Orbital_RKBS_IBS<double>(db,exponents,kappa) 
        )
{
    
};


std::ostream&  Orbital_RKB_IBS::Write(std::ostream& os) const
{
    os << "Dirac basis set." << std::endl << "    Large: " << *itsRKBL << std::endl << "    Small: " << *itsRKBS << std::endl;
    return os;
}


//-----------------------------------------------------------------------------------------------
//
//  Large sector
//
template <class T> Orbital_RKBL_IBS<T>::Orbital_RKBL_IBS
(const DB_cache<T>* db,const Vector<T>& exponents,int kappa)
    : Slater_IBS(exponents,Omega_kmj_Sym::l(kappa),{})
    , Atom::IrrepBasisSet(this,new Omega_k_Sym(kappa))
    , Atom::Orbital_RKBL_IBS<T>(db,this,kappa)
{
};




//-----------------------------------------------------------------------------------------------
//
//  Small sector
//
template <class T> Orbital_RKBS_IBS<T>::Orbital_RKBS_IBS
    (const DB_cache<double>* db,const Vector<T>& exponents, int kappa)
    : IrrepBasisSet_Common<T> (new Omega_k_Sym(-kappa))
    , Slater_RKBS_IBS(exponents,kappa,Omega_k_Sym::l(kappa),{})
    , Atom::Orbital_RKBS_IBS<T>(db,this,kappa)
{
};

template <class T> std::ostream&  Orbital_RKBS_IBS<T>::Write(std::ostream& os) const
{
    os << "Slater RKB " << this->GetSymmetry();
    if (kappa>0)
        os << "[ " << std::setw(2) << 2*kappa+1 << "/r - e ]";
    else
        os << "[       -e ]";
    os << "*r^" << Slater_IBS::l << "*exp(-e*r), e={";
    for (auto e: Slater_IBS::es) os << e << ",";
    os << "}";
    return os;
}


}} //namespace
