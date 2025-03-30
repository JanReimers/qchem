// File: BasisSetImplementation.C  General basis set implementation as a list.



#include "Imp/BasisSet/IBS_Common.H"
#include <QuantumNumber.H>
#include <BasisFunction.H>
#include "Imp/Containers/ptr_vector_io.h"

#include <iostream>
#include <iomanip>
#include <cassert>

//-----------------------------------------------------------------------------
//
//  Construction zone
//
IBS_Common::IBS_Common()
    : itsQuantumNumber(0)
{
};

IBS_Common::IBS_Common(QNs* theQN)
    : itsQuantumNumber(theQN)
{
    assert(itsQuantumNumber);
};

IBS_Common::IBS_Common(const IBS_Common& bs)
  : itsQuantumNumber (bs.itsQuantumNumber->Clone())
  , itsBasisFunctions(bs.itsBasisFunctions)
  {
    assert(itsQuantumNumber);
  };

IBS_Common::~IBS_Common()
{
    delete itsQuantumNumber;
}

//-----------------------------------------------------------------------------
//
//  Post construction initializations called by dervied classes.
//
void IBS_Common::Insert(const BasisFunction* bf)
{
    assert(bf);
    itsBasisFunctions.push_back(bf);
}

void IBS_Common::EmptyBasisFunctions()
{
    itsBasisFunctions.clear();
}


//-----------------------------------------------------------------------------
//
//  Basis Set Stuff.
//


size_t IBS_Common::GetNumFunctions() const
{
    return itsBasisFunctions.size();
}

bool IBS_Common::operator==(const IrrepBasisSet& bs) const
{
    // No UT coverage
    if (GetNumFunctions() != bs.GetNumFunctions()) return false;
    bool ret=true;
    auto b2=bs.Iterate<BasisFunction>().begin(); 
    for (auto b1:Iterate<BasisFunction>()) 
    {
        ret=ret && (*b1)==(**b2);
        ++b2;
    }
    return ret;
}

QNs*  IBS_Common::GetQuantumNumber(int index) const
{
    return itsQuantumNumber->AddPrincipleQN(index);
}


//-----------------------------------------------------------------------------
//
//  Streamable stuff.
//
std::ostream& IBS_Common::Write(std::ostream& os) const
{
    assert(itsQuantumNumber);
    if(!Pretty())
    {
        UniqueIDImp::Write(os);
        if(!Binary()) os << std::endl;
        os << *itsQuantumNumber;
        if(!Binary()) os << std::endl;
    }
    else
    {
        os << "IrrepBasisSet " << " with " << GetNumFunctions() << " basis functions"
        << ", Quantum number=" << *itsQuantumNumber <<  std::endl;
        os << "   #          Center      Pol     Radial     Exponents" << std::endl;
        // No UT coverage
        int i=1;
        for(auto b:*this) os << std::setw(4) << i++ << " " << *b << std::endl;
        
        os << "-------------------------------------------------------------------------------" << std::endl;
    }

    return os;
}

// std::istream& IrrepBasisSetCommon::Read(std::istream& is)
// {
//     UniqueID::Read(is);

//     delete itsQuantumNumber;
//     itsQuantumNumber=QuantumNumber::Factory(is);
//     is >> *itsQuantumNumber;

//     return is;
// };

//-----------------------------------------------------------------------------
//
//  Controlled by derived class for IO.
//
std::ostream& IBS_Common::WriteBasisFunctions(std::ostream& os) const
{
    os << itsBasisFunctions;
    return os;
}

std::istream& IBS_Common::ReadBasisFunctions(std::istream& is)
{
    // is >> itsBasisFunctions;
    return is;
}


