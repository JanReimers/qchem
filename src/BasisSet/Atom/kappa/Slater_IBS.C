// File: Atom/kappa/Slater_IBS.C  Slater Irrep Basis Set (IBS) with Restricted Kinetic Balance (RKB).

#include "Imp/BasisSet/Atom/kappa/Slater_IBS.H"
#include "Imp/BasisSet/Atom/kappa/Slater_BF.H"
#include "Imp/BasisSet/Atom/radial/Slater/Integrals.H"

#include "Imp/Symmetry/OkmjQN.H"
#include <QuantumNumber.H>
#include <iostream>
#include <cassert>

using std::endl;

namespace Atom_kappa
{
namespace Slater
{
//
//  Concrete  Slater basis set.
//

Dirac_IrrepBasisSet::Dirac_IrrepBasisSet(const DB_cache<double>* db,
    const Vector<double>& exponents,int kappa)
: Orbital_RKB_IBS_Common<double>(db,kappa )
{
    auto rkbl=new Large_Orbital_IBS<double>(db,exponents, kappa);
    auto rkbs=new Small_Orbital_IBS<double>(db,rkbl);
    Orbital_RKB_IBS_Common<double>::Init(rkbl,rkbs);
    Dirac_IrrepIEClient::Init(rkbl,rkbs);
    for (auto b:itsRKBL->Iterate<BasisFunction>()) Insert(b);
    for (auto b:itsRKBS->Iterate<BasisFunction>()) Insert(b);
};

std::ostream&  Dirac_IrrepBasisSet::Write(std::ostream& os) const
{
    if (!Pretty())
    {
        WriteBasisFunctions(os);
        IrrepBasisSetCommon::Write(os);
    }
    else
    {
        os << "Dirac basis set." << endl << "    Large: " << *itsRKBL << endl << "    Small: " << *itsRKBS << endl;
    }
    return os;
}

::IrrepBasisSet* Dirac_IrrepBasisSet::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a Slater atomic basis set?!" << std::endl;
    return 0;
}

//-----------------------------------------------------------------------------------------------
//
//  Large sector
//
template <class T> Large_Orbital_IBS<T>::Large_Orbital_IBS(const DB_cache<T>* db,
        const Vector<T>& exponents,int kappa)
    : Orbital_RKBL_IBS_Common<T>(kappa)
    , Orbital_RKBL_IE<T>(db)
    , IrrepIEClient(exponents.size(),kappa)
{
    size_t l=Omega_kmjQN::l(kappa);
    Init(exponents,Norms(exponents,l),l);
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new Large_BasisFunction(e,kappa,0.5,ns(i++))); //ns from Slater_mj::IEClient

};

template <class T> Vector<double> Large_Orbital_IBS<T>::Norms(const Vector<double>& es, size_t l) const
{
    Vector<double> ns(es.size());
    int i=0;
    for (auto e:es) ns(++i)=::Slater::Norm(e,l+1);
    return ns;
}


template <class T> std::ostream&  Large_Orbital_IBS<T>::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
    {
        os << "Slater     " << this->GetQuantumNumber()
        << "             r^" << l << "*exp(-e*r), e={";
        for (auto b:*this) os << *b;
        os << "}";
    }
    return os;
}

template <class T> ::IrrepBasisSet* Large_Orbital_IBS<T>::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a Slater atomic basis set?!" << std::endl;
    return 0;
}

//-----------------------------------------------------------------------------------------------
//
//  Small sector
//
template <class T> Small_Orbital_IBS<T>::Small_Orbital_IBS(const DB_cache<double>* db,const Large_Orbital_IBS<T>* lbs)
    : Orbital_RKBS_IBS_Common<T>(lbs->kappa)
    , Orbital_RKBS_IE<T>(db)
    , Small_IrrepIEClient(lbs->size(),lbs->kappa)
{
    size_t l=Omega_kmjQN::l(kappa);
    Init(lbs->es,Norms(lbs->es,l),l);
    size_t i=1;
    for (auto b:*lbs) 
    {
        const Large_BasisFunction* lb=dynamic_cast<const Large_BasisFunction*>(b);
        IrrepBasisSetCommon::Insert(new Small_BasisFunction(lb,ns(i++))); 
    }

};

template <class T> Vector<double> Small_Orbital_IBS<T>::Norms(const Vector<double>& es, size_t l) const
{
    Vector<double> ns(es.size());
    int i=0;
    for (auto e:es) ns(++i)=1.0/sqrt(::Slater::IE_Primatives::Grad2(e,e,l,l));
    return ns;
}

template <class T> std::ostream&  Small_Orbital_IBS<T>::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
    {
        os << "Slater RKB " << kappa << endl << this->GetQuantumNumber();
        if (kappa>0)
            os << "[ " << std::setw(2) << 2*kappa+1 << "/r - e ]";
        else
            os << "[       -e ]";
        os << "*r^" << l << "*exp(-e*r), e={";
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
