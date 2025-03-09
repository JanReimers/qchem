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
    template <class T> Large_Orbital_IBS<T>::Large_Orbital_IBS(const LAParams& lap,IntegralDataBase<T>* db,
        const Vector<T>& exponents,int kappa)
        : IrrepBasisSetCommon(new Omega_kQN(kappa))
        , TIrrepBasisSetCommon<T>(lap,db)
        , IrrepIEClient(exponents.size(),kappa)
    {
        IrrepIEClient::Init(exponents);
        size_t i=1;
        for (auto e:es) 
            IrrepBasisSetCommon::Insert(new Large_BasisFunction(e,kappa,ns(i++))); //ns from Slater_mj::IEClient
    };

template <class T> Small_Orbital_IBS<T>::Small_Orbital_IBS(const LAParams& lap,IntegralDataBase<T>* db
    ,const Large_Orbital_IBS<T>* lbs)
    : IrrepBasisSetCommon(new Omega_kQN(-lbs->kappa))
    , TIrrepBasisSetCommon<T>(lap,db)
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

  
Dirac_IrrepBasisSet::Dirac_IrrepBasisSet(const LAParams& lap,IntegralDataBase<double>* theDB,
        const Vector<double>& exponents,int kappa)
    : Dirac::IrrepBasisSet<double>(lap,theDB,new Large_Orbital_IBS<double>(lap,theDB,exponents, kappa),kappa )
{
    auto rkbl=dynamic_cast<Large_Orbital_IBS<double>*>(itsRKBL);
    assert(rkbl);
    auto rkbs=new Small_Orbital_IBS<double>(lap,theDB,rkbl);
    itsRKBS=rkbs;
    assert(itsRKBS);
    IntegralEngine::itsRKBS=rkbs;
    Dirac_IrrepIEClient::Init(rkbl,rkbs);
    for (auto b:*itsRKBL) Insert(b);
    for (auto b:*itsRKBS) Insert(b);
};

::Fit_IBS* Dirac_IrrepBasisSet::CreateCDFitBasisSet(const Cluster*) const
{
    assert(false);
    return 0;
//    return new Dirac_IrrepBasisSet(itsLAParams,GetDataBase(),es*2,-1,0.5);
}

::Fit_IBS* Dirac_IrrepBasisSet::CreateVxcFitBasisSet(const Cluster*) const
{
    assert(false);
    return 0;
//    return new Dirac_IrrepBasisSet(itsLAParams,GetDataBase(),es*2.0/3.0,-1,0.5);    
}

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

::IrrepBasisSet* Dirac_IrrepBasisSet::Clone() const
{
    return new Dirac_IrrepBasisSet(*this);
}

::IrrepBasisSet* Dirac_IrrepBasisSet::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a Slater atomic basis set?!" << std::endl;
    return Clone();
}

// Dirac_IrrepBasisSet::SMat Dirac_IrrepBasisSet::MakeOverlap() const
// {
//     SMat ol=itsLargeBS->Overlap();
//     SMat os=itsSmallBS->Overlap();
//     return DiracIntegralEngine::merge_diag(ol,os);
// }

// Dirac_IrrepBasisSet::SMat Dirac_IrrepBasisSet::MakeKinetic() const
// {
//     const Integrals_RKB<double>* irkb=itsSmallBS;
//     Matrix<double> k=-irkb->Kinetic(itsLargeBS); 
//     //std::cout << "k=" << k << std::endl;
//     return DiracIntegralEngine::merge_off_diag(k);
// }

// Dirac_IrrepBasisSet::SMat Dirac_IrrepBasisSet::MakeNuclear(const Cluster* cl) const
// {
//     SMat ol=itsLargeBS->Nuclear(cl);
//     SMat os=itsSmallBS->Nuclear(cl);
//     return DiracIntegralEngine::merge_diag(ol,os);
// }




//         assert(false);
//     }
//     return SMat();
// }


} //namespace
