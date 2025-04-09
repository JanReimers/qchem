// File: Atom/kappa/Gaussian_IBS.C  Restricted Kinetic Balance (RKB) Irrep Basis Set (IBS).


#include "Imp/BasisSet/Atom/kappa/Gaussian_IBS.H"
#include "Imp/BasisSet/Atom/kappa/Gaussian_BF.H"
#include "Imp/BasisSet/Atom/radial/Gaussian/Integrals.H"

#include "Imp/Symmetry/Okmj.H"
#include <Symmetry.H>
#include <iostream>
#include <cassert>

using std::endl;

namespace Atom_kappa
{
namespace Gaussian
{
  
Orbital_IBS::Orbital_IBS
    (const DB_cache<double>* db
        , const Vector<double>& exponents
        , int kappa)
    : Orbital_RKB_IBS_Common<double>
        (db
        , new Omega_k_Sym(kappa)
        , kappa
        , new Large_Orbital_IBS<double>(db,exponents, kappa)
        , new Small_Orbital_IBS<double>(db,exponents,kappa)
        )
{
    for (auto b:itsRKBL->Iterate<BasisFunction>()) Insert(b);
    for (auto b:itsRKBS->Iterate<BasisFunction>()) Insert(b);
};

std::ostream&  Orbital_IBS::Write(std::ostream& os) const
{
    if (!Pretty())
    {
        WriteBasisFunctions(os);
        IBS_Common::Write(os);
    }
    else
    {
        os << "Dirac basis set." << endl << "    Large: " << *itsRKBL << endl << "    Small: " << *itsRKBS << endl;
    }
    return os;
}
::IrrepBasisSet* Orbital_IBS::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a Slater atomic basis set?!" << std::endl;
    return 0;
}


template <class T> Large_Orbital_IBS<T>::Large_Orbital_IBS(const DB_cache<T>* db,
    const Vector<T>& exponents,int kappa)
    : Orbital_RKBL_IBS_Common<T>(new Omega_k_Sym(kappa),kappa)
    , Orbital_RKBL_IE<T>(db)
    , AtomIrrepIEClient(exponents.size())
{
    size_t l=Omega_kmj_Sym::l(kappa);
    AtomIrrepIEClient::Init(exponents,Norms(exponents,l),l);    
    size_t i=1;
    for (auto e:es) 
        IBS_Common::Insert(new Large_BasisFunction(e,kappa,ns(i++))); //ns from Slater_mj::IEClient
};
template <class T> Vector<double> Large_Orbital_IBS<T>::Norms(const Vector<double>& es, size_t l) const
{
    Vector<double> ns(es.size());
    int i=0;
    for (auto e:es) ns(++i)=::Gaussian::Norm(e,l);
    return ns;
}
template <class T> std::ostream&  Large_Orbital_IBS<T>::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
    {
        os << "Gaussian     " << this->GetSymmetry()
        << "               r^" << l << "*exp(-e*r^2), e={";
        for (auto b:*this) os << *b;
        os << "}";
    }
    return os;
}
template <class T> ::IrrepBasisSet* Large_Orbital_IBS<T>::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a Slater atomic basis set?!" << std::endl;
    assert(false);
    return 0;
}

template <class T> Small_Orbital_IBS<T>::Small_Orbital_IBS(const DB_cache<T>* db,
    const Vector<T>& exponents,int kappa)
    : Orbital_RKBS_IBS_Common<T>(new Omega_k_Sym(-kappa),kappa)
    , Orbital_RKBS_IE<T>(db)
    , AtomIrrepIEClient(exponents.size())
{
    size_t l=Omega_kmj_Sym::l(kappa);
    AtomIrrepIEClient::Init(exponents,Norms(exponents,l),l);  
    
}

template <class T> void Small_Orbital_IBS<T>::InsertBasisFunctions(const Orbital_RKBL_IBS<T>* lbs)
{
    size_t i=1;
    for (auto lb:lbs->template Iterate<Large_BasisFunction>()) 
        IBS_Common::Insert(new Small_BasisFunction(lb,ns(i++))); 
}

template <class T> Vector<double> Small_Orbital_IBS<T>::Norms(const Vector<double>& es, size_t l) const
{
    Vector<double> ns(es.size());
    int i=0;
    for (auto e:es) ns(++i)=1.0/sqrt(::Gaussian::IE_Primatives::Grad2(e,e,l,l));
    return ns;
}
template <class T> std::ostream&  Small_Orbital_IBS<T>::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
    {
        os << "Gaussian     " << this->GetSymmetry()
        << "               r^" << l << "*exp(-e*r^2), e={";
        for (auto b:*this) os << *b;
        os << "}";
    }
    return os;
}
template <class T> ::IrrepBasisSet* Small_Orbital_IBS<T>::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a Slater atomic basis set?!" << std::endl;
    return 0;
}



}} //namespace
