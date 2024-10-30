// File: Slater_m/IrrepBasisSet.C  Spherical Slater basis set with orbital angular momentum l,m.

#include "Imp/BasisSet/Slater_m/IrrepBasisSet.H"
#include "Imp/BasisSet/Slater_m/BasisFunction.H"
#include "Imp/BasisSet/Slater_m/IntegralEngine.H"
#include "Imp/BasisSet/Slater_m/QuantumNumber.H"
#include <iostream>
#include <cassert>

namespace Slater_m
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
        size_t L, int m)
    : IrrepBasisSetCommon(new YlmQN(L,m))
    , TIrrepBasisSetCommon<double>(lap,theDB)
    , IrrepIEClient(size)
{
    IrrepIEClient::Init(minexp,maxexp,L,m);
    TIrrepBasisSetCommon<double>::Insert(new IntegralEngine());  
    size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new BasisFunction(e,L+1,L,m,ns(i++))); //ns from SlaterIEClient

};

IrrepBasisSet* IrrepBasisSet::CreateCDFitBasisSet(const Cluster*) const
{
    double emin=es(1), emax=es(size());
    return new IrrepBasisSet(itsLAParams,GetDataBase(),size(),emin*2,emax*2,0,0);
}

IrrepBasisSet* IrrepBasisSet::CreateVxcFitBasisSet(const Cluster*) const
{
    double emin=es(1), emax=es(size());
    return new IrrepBasisSet(itsLAParams,GetDataBase(),size(),emin*2.0/3.0,emax*2.0/3.0,0,0);    
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
        os << "Slater functions l,m=" << GetQuantumNumber()
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
