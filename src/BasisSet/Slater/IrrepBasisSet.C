// File: SlaterBS.C  Spherical Slater basis set.

#include "Imp/BasisSet/Slater/IrrepBasisSet.H"
#include "Imp/BasisSet/Slater/BasisFunction.H"
#include "Imp/BasisSet/Slater/IntegralEngine.H"
#include "Imp/Symmetry/YlQN.H"
#include <iostream>
#include <cassert>

namespace Slater
{
//
//  Concrete  Slater basis set.
//
// IrrepBasisSet::IrrepBasisSet()
//     :  IrrepBasisSetCommon        ()
//     , Orbital_IBS_Common<double>()
// {};


IrrepBasisSet::IrrepBasisSet(const LAParams& lap,IntegralDataBase<double>* theDB,const DB_BS_2E<double>* db,
        const Vector<double>& exponents,size_t L)
    : IrrepBasisSetCommon(new YlQN(L))
    , Orbital_IBS_Common<double>(lap,theDB)
    , IntegralEngine1(db)
    , IrrepIEClient(exponents.size())
{
    IrrepIEClient::Init(exponents,L);
     size_t i=1;
    for (auto e:es) 
        IrrepBasisSetCommon::Insert(new BasisFunction(e,L+1,L,ns(i++))); //ns from SlaterIEClient

};

::Fit_IBS* IrrepBasisSet::CreateCDFitBasisSet(const Cluster*) const
{
    // return new IrrepBasisSet(itsLAParams,GetDataBase(),0,es*2,0);
    assert(false);
    return 0;
}

::Fit_IBS* IrrepBasisSet::CreateVxcFitBasisSet(const Cluster*) const
{
    // return new IrrepBasisSet(itsLAParams,GetDataBase(),0,es*2.0/3.0,0);    
    assert(false);
    return 0;
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
