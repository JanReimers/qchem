// File: SphericalGaussian_RKB_IrrepBasisSet.C  Spherical gaussian basis set.


#include "Imp/BasisSet/SG_RKB/IrrepBasisSet.H"
#include "Imp/BasisSet/SG_RKB/BasisFunction.H"
#include "Imp/BasisSet/SG_RKB/IntegralEngine.H"
#include "Imp/Symmetry/OkmjQN.H"
#include <iostream>
#include <cassert>

using std::endl;

namespace SphericalGaussian_RKB
{
  
Dirac_IrrepBasisSet::Dirac_IrrepBasisSet(const LAParams& lap,IntegralDataBase<double>* theDB,
        const Vector<double>& exponents,int kappa)
    : IrrepBasisSetCommon(new Omega_kQN(kappa))
    , TIrrepBasisSetCommon<double>(lap,theDB)
    , itsLargeBS(new Large_IrrepBasisSet(lap,theDB,exponents,kappa))
    , itsSmallBS(new Small_IrrepBasisSet(lap,theDB,itsLargeBS))
{
    Dirac_IrrepIEClient::Init(itsLargeBS,itsSmallBS);
    for (auto b:*itsLargeBS) Insert(b);
    for (auto b:*itsSmallBS) Insert(b);
};

IrrepBasisSet* Dirac_IrrepBasisSet::CreateCDFitBasisSet(const Cluster*) const
{
    assert(false);
    return 0;
//    return new Dirac_IrrepBasisSet(itsLAParams,GetDataBase(),es*2,-1,0.5);
}

IrrepBasisSet* Dirac_IrrepBasisSet::CreateVxcFitBasisSet(const Cluster*) const
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
        os << "Dirac basis set." << endl << "    Large: " << *itsLargeBS << endl << "    Small: " << *itsSmallBS << endl;
    }
    return os;
}

std::istream&  Dirac_IrrepBasisSet::Read (std::istream& is)
{
    ReadBasisFunctions(is);
    IrrepBasisSetCommon::Read(is);
    TIrrepBasisSetCommon<double>::Read(is);
    return is;
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

//-----------------------------------------------------------------------------------------------
//
//  Large sector
//
Large_IrrepBasisSet::Large_IrrepBasisSet(const LAParams& lap,IntegralDataBase<double>* theDB,
        const Vector<double>& exponents,int kappa)
    : IrrepBasisSetCommon(new Omega_kQN(kappa))
    , TIrrepBasisSetCommon<double>(lap,theDB)
    , IrrepIEClient(exponents.size(),kappa)
{
    IrrepIEClient::Init(exponents);
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new Large_BasisFunction(e,kappa,ns(i++))); //ns from Slater_mj::IEClient

};

::IrrepBasisSet* Large_IrrepBasisSet::CreateCDFitBasisSet(const Cluster*) const
{
    return new Large_IrrepBasisSet(itsLAParams,GetDataBase(),es*2,0);
}

::IrrepBasisSet* Large_IrrepBasisSet::CreateVxcFitBasisSet(const Cluster*) const
{
    return new Large_IrrepBasisSet(itsLAParams,GetDataBase(),es*2.0/3.0,0);    
}

std::ostream&  Large_IrrepBasisSet::Write(std::ostream& os) const
{
    if (!Pretty())
    {
        WriteBasisFunctions(os);
        IrrepBasisSetCommon::Write(os);
        TIrrepBasisSetCommon<double>::Write(os);
    }
    else
    {
        os << "Gaussian     " << GetQuantumNumber()
        << "               r^" << l << "*exp(-e*r^2), e={";
        for (auto b:*this) os << *b;
        os << "}";
    }
    return os;
}

std::istream&  Large_IrrepBasisSet::Read (std::istream& is)
{
    ReadBasisFunctions(is);
    IrrepBasisSetCommon::Read(is);
    TIrrepBasisSetCommon<double>::Read(is);
    return is;
}

::IrrepBasisSet* Large_IrrepBasisSet::Clone() const
{
    return new Large_IrrepBasisSet(*this);
}

::IrrepBasisSet* Large_IrrepBasisSet::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a Slater atomic basis set?!" << std::endl;
    return Clone();
}

//-----------------------------------------------------------------------------------------------
//
//  Small sector
//
Small_IrrepBasisSet::Small_IrrepBasisSet(const LAParams& lap,IntegralDataBase<double>* db,const Large_IrrepBasisSet* lbs)
    : IrrepBasisSetCommon(new Omega_kQN(-lbs->kappa))
    , TIrrepBasisSetCommon<double>(lap,db)
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

::IrrepBasisSet* Small_IrrepBasisSet::CreateCDFitBasisSet(const Cluster*) const
{
    assert(false);
    return 0;
//    return new Small_IrrepBasisSet(itsLAParams,GetDataBase(),es*2,0,0);
}

::IrrepBasisSet* Small_IrrepBasisSet::CreateVxcFitBasisSet(const Cluster*) const
{
    assert(false);
    return 0;
//    return new Small_IrrepBasisSet(itsLAParams,GetDataBase(),es*2.0/3.0,0,0);    
}

std::ostream&  Small_IrrepBasisSet::Write(std::ostream& os) const
{
    if (!Pretty())
    {
        WriteBasisFunctions(os);
        IrrepBasisSetCommon::Write(os);
        TIrrepBasisSetCommon<double>::Write(os);
    }
    else
    {
        os << "Gaussian RKB " << GetQuantumNumber();
        if (kappa>0)
            os << "[ " << std::setw(2) << 2*kappa+1 << "/r - 2er ]";
        else
            os << "[       -2er ]";
        os << "*r^" << l << "*exp(-e*r^2), e={";
        for (auto b:*this) os << *b;
        os << "}";
        
    }
    return os;
}

std::istream&  Small_IrrepBasisSet::Read (std::istream& is)
{
    ReadBasisFunctions(is);
    IrrepBasisSetCommon::Read(is);
    TIrrepBasisSetCommon<double>::Read(is);
    return is;
}

::IrrepBasisSet* Small_IrrepBasisSet::Clone() const
{
    return new Small_IrrepBasisSet(*this);
}

::IrrepBasisSet* Small_IrrepBasisSet::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a Slater atomic basis set?!" << std::endl;
    return Clone();
}

} //namespace
