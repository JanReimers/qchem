// File: SphericalGaussian_RKB_IrrepBasisSet.C  Spherical gaussian basis set.


#include "Imp/BasisSet/SG_RKB/IrrepBasisSet.H"
#include "Imp/BasisSet/SG_RKB/BasisFunction.H"
#include "Imp/Symmetry/OkmjQN.H"
#include <iostream>
#include <cassert>

using std::endl;

namespace SphericalGaussian_RKB
{
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
            IrrepBasisSetCommon::Insert(new Large_BasisFunction(e,kappa,ns(i++))); //ns from Slater_mj::IEClient
    };

template <class T> Small_Orbital_IBS<T>::Small_Orbital_IBS(const LAParams& lap,const DB_cache<T>* db,
    const Large_Orbital_IBS<T>* lbs)
    : IrrepBasisSetCommon(new Omega_kQN(-lbs->kappa))
    , TIrrepBasisSetCommon<T>(lap)
    , Orbital_RKBS_IE<T>(db)
    , Small_IrrepIEClient(lbs->size(),lbs->kappa)
{
    IrrepIEClient::Init(lbs->es);
    size_t i=1;
    for (auto b:*lbs) 
    {
        const Large_BasisFunction* lb=dynamic_cast<const Large_BasisFunction*>(b);
        IrrepBasisSetCommon::Insert(new Small_BasisFunction(lb,ns(i++))); 
    };
}


template <class T> std::ostream&  Large_Orbital_IBS<T>::Write(std::ostream& os) const
{
    if (Pretty())
    {
        os << "Gaussian     " << GetQuantumNumber()
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
template <class T> std::ostream&  Small_Orbital_IBS<T>::Write(std::ostream& os) const
{
    if (Pretty())
    {
        os << "Gaussian     " << GetQuantumNumber()
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


  
Dirac_IrrepBasisSet::Dirac_IrrepBasisSet(const LAParams& lap,const DB_cache<double>* db, const Vector<double>& exponents, int kappa)
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


} //namespace
