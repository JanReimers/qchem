// File: Slater_m/IrrepBasisSet.C  Spherical Slater basis set with orbital angular momentum l,m.

#include "Imp/BasisSet/Slater_mj/IrrepBasisSet.H"
#include "Imp/BasisSet/Slater_mj/BasisFunction.H"
#include "Imp/BasisSet/Slater_mj/IntegralEngine.H"
#include "Imp/Symmetry/OkmjQN.H"
#include "Imp/Integrals/SlaterIntegrals.H"

#include <iostream>
#include <cassert>

using std::endl;

namespace Slater_mj
{
//
//  Concrete  Slater basis set.
//

Dirac_IrrepBasisSet::Dirac_IrrepBasisSet(const LAParams& lap,const DB_cache<double>* db,
    const Vector<double>& exponents,int kappa)
: Dirac::IrrepBasisSet<double>(lap,db,kappa )
{
    auto rkbl=new Large_Orbital_IBS<double>(lap,db,exponents, kappa);
    auto rkbs=new Small_Orbital_IBS<double>(lap,db,rkbl);
    Dirac::IrrepBasisSet<double>::Init(rkbl,rkbs);
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
        TIrrepBasisSetCommon<double>::Write(os);
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
template <class T> Large_Orbital_IBS<T>::Large_Orbital_IBS(const LAParams& lap,const DB_cache<T>* db,
        const Vector<T>& exponents,int kappa)
    : IrrepBasisSetCommon(new Omega_kQN(kappa))
    , TIrrepBasisSetCommon<T>(lap)
    , Orbital_RKBL_IE<T>(db)
    , IrrepIEClient(exponents.size(),kappa)
{
    IrrepIEClient::Init(exponents);
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new Large_BasisFunction(e,kappa,0.5,ns(i++))); //ns from Slater_mj::IEClient

};


template <class T> std::ostream&  Large_Orbital_IBS<T>::Write(std::ostream& os) const
{
    if (Pretty())
    {
        os << "Slater     " << GetQuantumNumber()
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
template <class T> Small_Orbital_IBS<T>::Small_Orbital_IBS(const LAParams& lap,const DB_cache<double>* db,const Large_Orbital_IBS<T>* lbs)
    : IrrepBasisSetCommon(new Omega_kQN(-lbs->kappa))
    , TIrrepBasisSetCommon<T>(lap)
    , AtomIE_RKBS<T>(db)
    , Small_IrrepIEClient(lbs->size(),lbs->kappa)
{
  Small_IrrepIEClient::Init(lbs->es);
  size_t i=1;
  for (auto b:*lbs) 
    {
        const Large_BasisFunction* lb=dynamic_cast<const Large_BasisFunction*>(b);
        IrrepBasisSetCommon::Insert(new Small_BasisFunction(lb,ns(i++))); 
    }

};

template <class T> std::ostream&  Small_Orbital_IBS<T>::Write(std::ostream& os) const
{
    if (Pretty())
    {
        os << "Slater RKB " << GetQuantumNumber();
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


template <class T> double Small_Orbital_IBS<T>::Overlap(double ea , double eb,size_t l_total) const
{
    assert(false);
    return 0.0;
}
template <class T> double Small_Orbital_IBS<T>::Kinetic(double ea , double eb,size_t la, size_t lb) const
{
    double ab=ea+eb;
    int na=l+1,nb=l+1;
    size_t ll=(l*(l+1)+l*(l+1))/2;
    int n=na+nb;
    double Term1=0.5*(na*nb+ll)*SlaterIntegral(ab,n-2); //SlaterIntegral already has 4*Pi
    double Term2=-0.5*(na*eb+nb*ea)* SlaterIntegral(ab,n-1);
    double Term3=0.5*ea*eb*SlaterIntegral(ab,n);
    //cout << "Slater::IntegralEngine::Kinetic Terms 1,2,3=" << Term1 << " " << Term2 << " " << Term3 << endl;

    return 2.0*(Term1+Term2+Term3);
}
template <class T> double Small_Orbital_IBS<T>::Nuclear(double ea , double eb,size_t l_total) const
{
    return ea*eb*SlaterIntegral(ea+eb,l_total+1);
}

} //namespace
