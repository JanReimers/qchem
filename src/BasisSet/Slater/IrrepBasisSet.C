// File: SlaterBS.C  Spherical Slater basis set.

#include "Imp/BasisSet/Slater/IrrepBasisSet.H"
#include "Imp/BasisSet/Slater/BasisFunction.H"
#include "Imp/BasisSet/Slater/IntegralEngine.H"
#include "Imp/BasisSet/SphericalGaussian/QuantumNumber.H"
#include <iostream>
#include <cassert>

namespace Slater
{
//
//  Concrete  Slater basis set.
//
IrrepBasisSet::IrrepBasisSet()
    :  IrrepBasisSetCommon        ()
    , TIrrepBasisSetCommon<double>()
{};


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
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new BasisFunction(e,L+1,L,ns(i++))); //ns from SlaterIEClient

};

IrrepBasisSet::IrrepBasisSet(const LAParams& lap,IntegralDataBase<double>* theDB,
        const Vector<double>& exponents,size_t L)
    : IrrepBasisSetCommon(new SphericalSymmetryQN(L))
    , TIrrepBasisSetCommon<double>(lap,theDB)
    , IrrepIEClient(exponents.size())
{
    IrrepIEClient::Init(exponents,L);
     size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new BasisFunction(e,L+1,L,ns(i++))); //ns from SlaterIEClient

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
        os << "Slater functions L=" << GetQuantumNumber()
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
    std::cerr << "Why are you relocating a Slater atomic basis set?!" << std::endl;
    return Clone();
}


} //namespace
