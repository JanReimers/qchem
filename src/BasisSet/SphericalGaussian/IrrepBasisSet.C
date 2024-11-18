// File: SphericalGaussianBS.C  Spherical gaussian basis set.



#include "Imp/BasisSet/SphericalGaussian/IrrepBasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/BasisFunction.H"
#include "Imp/BasisSet/SphericalGaussian/IntegralEngine.H"
#include "Imp/BasisSet/SphericalGaussian/QuantumNumber.H"
#include <iostream>
#include <cassert>

namespace SphericalGaussian
{
  
//#######################################################################
//
//  Concrete  gaussian basis set.
//
IrrepBasisSet::IrrepBasisSet()
    :  IrrepBasisSetCommon        ()
    , TIrrepBasisSetCommon<double>()
{};

//
//  We need three constructors type here.  They all need DB, size, exponents,L
//    1) For HF orbitals also need LAParams for secular eq. solving
//    2) For DFT orbitals need 1+LAParams for overlap inversion.
//    3) For DFT Vxc, and ro fitting we need defulat + mesh.
//

IrrepBasisSet::IrrepBasisSet(
        const LAParams& lap,
        IntegralDataBase<double>* theDB,
        size_t size,
        double minexp,
        double maxexp,
        size_t L)
    : IrrepBasisSetCommon(new SphericalSymmetryQN(L))
    , TIrrepBasisSetCommon<double>(lap,theDB)
    , IrrepIEClient(size)
{
    IrrepIEClient::Init(minexp,maxexp,L);
    TIrrepBasisSetCommon<double>::Insert(new IntegralEngine());  
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new BasisFunction(e,L,ns(i++))); //ns from SphericalGaussianIEClient

};

IrrepBasisSet::IrrepBasisSet(const LAParams& lap,IntegralDataBase<double>* theDB,
        const Vector<double>& exponents,size_t L)
    : IrrepBasisSetCommon(new SphericalSymmetryQN(L))
    , TIrrepBasisSetCommon<double>(lap,theDB)
    , IrrepIEClient(exponents.size())
{
    IrrepIEClient::Init(exponents,L);
    TIrrepBasisSetCommon<double>::Insert(new IntegralEngine());  
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new BasisFunction(e,L,ns(i++))); //ns from SlaterIEClient

};


IrrepBasisSet* IrrepBasisSet::CreateCDFitBasisSet(const Cluster*) const
{
    double emin=es(1), emax=es(size());
    return new IrrepBasisSet(itsLAParams,GetDataBase(),size(),emin*2,emax*2,0);
}

IrrepBasisSet* IrrepBasisSet::CreateVxcFitBasisSet(const Cluster*) const
{
    double emin=es(1), emax=es(size());
    return new IrrepBasisSet(itsLAParams,GetDataBase(),size(),emin*2.0/3.0,emax*2.0/3.0,0);    
}


std::ostream&  IrrepBasisSet::Write(std::ostream& os) const
{
    if (!Pretty())
    {
        WriteBasisFunctions(os);
        IrrepBasisSetCommon::Write(os);
        TIrrepBasisSetCommon<double>::Write(os);
    }
    else
    {
        os << "Spherical Gaussian L=" << GetQuantumNumber()
        << " with " << GetNumFunctions() << " basis functions, alpha={";
        for (auto b:*this) os << *b;
        os << "}" << std::endl;
    }
    return os;
}

std::istream&  IrrepBasisSet::Read (std::istream& is)
{
    ReadBasisFunctions(is);
    IrrepBasisSetCommon::Read(is);
    TIrrepBasisSetCommon<double>::Read(is);
    return is;
}

::IrrepBasisSet* IrrepBasisSet::Clone() const
{
    return new IrrepBasisSet(*this);
}

::IrrepBasisSet* IrrepBasisSet::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a spherical Gaussian basis set?!" << std::endl;
    return Clone();
}


} //namespace
