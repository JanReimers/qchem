// File: SphericalGaussian_RKB_IrrepBasisSet.C  Spherical gaussian basis set.


#include "Imp/BasisSet/SG_RKB/IrrepBasisSet.H"
#include "Imp/BasisSet/SG_RKB/BasisFunction.H"
//#include "Imp/BasisSet/SG_RKB/IntegralEngine.H"
#include "Imp/Integrals/GaussianIntegrals.H"
#include "Imp/Symmetry/OkmjQN.H"
#include <iostream>
#include <cassert>

using std::endl;

namespace SphericalGaussian_RKB
{
    template <class T> Large_Orbital_IBS<T>::Large_Orbital_IBS(const LAParams& lap,
        const Vector<T>& exponents,int kappa)
        : IrrepBasisSetCommon(new Omega_kQN(kappa))
        , TIrrepBasisSetCommon<T>(lap)
        , IrrepIEClient(exponents.size(),kappa)
    {
        IrrepIEClient::Init(exponents);
        size_t i=1;
        for (auto e:es) 
            IrrepBasisSetCommon::Insert(new Large_BasisFunction(e,kappa,ns(i++))); //ns from Slater_mj::IEClient
    };

template <class T> Small_Orbital_IBS<T>::Small_Orbital_IBS(const LAParams& lap
    ,const Large_Orbital_IBS<T>* lbs)
    : IrrepBasisSetCommon(new Omega_kQN(-lbs->kappa))
    , TIrrepBasisSetCommon<T>(lap)
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

template <class T> T Large_Orbital_IBS<T>::Integral(qchem::IType it,double ea , double eb,size_t l) const
{
    if (it==qchem::Overlap1)
    {
        return GaussianIntegral(ea+eb,2*l);
    }
    else if(it==qchem::Kinetic1)
    {
        double t=ea+eb;
        size_t l1=l+1;
        return 0.5*(
                (l1*l1 + l*l1) * GaussianIntegral(t,2*l-2)
                -2*l1 * t      * GaussianIntegral(t,2*l  )
                +4*ea*eb       * GaussianIntegral(t,2*l+2)
            );
    }
    else if (it==qchem::Nuclear1)
    {
        return GaussianIntegral(ea+eb,2*l-1);
    }
    assert(false);
    return 0.0;
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


template <class T> T Small_Orbital_IBS<T>::Integral(qchem::IType it,double ea , double eb,size_t l) const
{
    if (it==qchem::Overlap1)
    {
        assert(l==0);
        double eab=ea+eb;
        size_t l1=l+1;
        return 1.0*(
                (l1*l1 + l*l1) * GaussianIntegral(eab,2*l-2)
                -2*l1 * eab      * GaussianIntegral(eab,2*l  )
                +4*ea*eb       * GaussianIntegral(eab,2*l+2)
            );
    }
    else if (it==qchem::Nuclear1)
    {
        assert(l==0);
        //int kappa = -l -1;
        return 4*ea*eb*GaussianIntegral(ea+eb,l+1); //Don't count the r^2 in dr^3
    }
    assert(false);
    return 0.0;
}

  
Dirac_IrrepBasisSet::Dirac_IrrepBasisSet(const LAParams& lap,const DB_cache<double>* db, const Vector<double>& exponents, int kappa)
    : Dirac::IrrepBasisSet<double>(lap,db,new Large_Orbital_IBS<double>(lap,exponents, kappa),kappa )
{
    auto rkbl=dynamic_cast<Large_Orbital_IBS<double>*>(itsRKBL);
    assert(rkbl);
    auto rkbs=new Small_Orbital_IBS<double>(lap,rkbl);
    itsRKBS=rkbs;
    assert(itsRKBS);
    RKB_IE::itsRKBS=rkbs;
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


} //namespace
