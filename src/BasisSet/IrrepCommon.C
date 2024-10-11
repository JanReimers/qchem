// File: BasisSetImplementation.C  General basis set implementation as a list.



#include "Imp/BasisSet/IrrepCommon.H"
#include <QuantumNumber.H>
#include "Imp/Containers/ptr_vector_io.h"

#include <iostream>
#include <iomanip>
#include <cassert>

//-----------------------------------------------------------------------------
//
//  Construction zone
//
IrrepBasisSetCommon::IrrepBasisSetCommon()
    : itsQuantumNumber(0)
    , itsStartIndex(0)
{
};

IrrepBasisSetCommon::IrrepBasisSetCommon(QuantumNumber* theQN)
    : itsQuantumNumber(theQN)
    , itsStartIndex(0)
{
};

IrrepBasisSetCommon::IrrepBasisSetCommon(const IrrepBasisSetCommon& bs)
  : itsQuantumNumber (bs.itsQuantumNumber->Clone())
  , itsBasisFunctions(bs.itsBasisFunctions)
  , itsStartIndex(bs.itsStartIndex)
  {
    assert(itsQuantumNumber);
  };

IrrepBasisSetCommon::~IrrepBasisSetCommon()
{
    delete itsQuantumNumber;
}

//-----------------------------------------------------------------------------
//
//  Post construction initializations called by dervied classes.
//
void IrrepBasisSetCommon::Insert(BasisFunction* bf)
{
    assert(bf);
    itsBasisFunctions.push_back(bf);
}

void IrrepBasisSetCommon::EmptyBasisFunctions()
{
    itsBasisFunctions.clear();
}


//-----------------------------------------------------------------------------
//
//  Basis Set Stuff.
//
void IrrepBasisSetCommon::SetStartIndex(int si)
{
    assert(si>0);
    itsStartIndex=si;
}

int  IrrepBasisSetCommon::GetStartIndex() const
{
    assert(itsStartIndex>0);
    return itsStartIndex;
}


size_t IrrepBasisSetCommon::GetNumFunctions() const
{
    return itsBasisFunctions.size();
}

bool IrrepBasisSetCommon::operator==(const IrrepBasisSet& bs) const
{
    // No UT coverage
    if (GetNumFunctions() != bs.GetNumFunctions()) return false;
    bool ret=true;
    auto b2=bs.begin(); 
    for (auto b1:*this) 
    {
        ret=ret && (*b1)==(**b2);
        b2++;
    }
    return ret;
}



//-----------------------------------------------------------------------------
//
//  Streamable stuff.
//
std::ostream& IrrepBasisSetCommon::Write(std::ostream& os) const
{
    assert(itsQuantumNumber);
    if(!Pretty())
    {
        UniqueID::Write(os);
        if(!Binary()) os << std::endl;
        os << *itsQuantumNumber;
        if(!Binary()) os << std::endl;
    }
    else
    {
        os << "Type" << " with " << GetNumFunctions() << " basis functions"
        << ", Quantum number=" << "QN Type" << std::endl;
        os << "   #          Center      Pol     Radial     Exponents" << std::endl;
        // No UT coverage
        int i=1;
        for(auto b:*this) os << std::setw(4) << i++ << " " << *b << std::endl;
    }

    return os;
}

std::istream& IrrepBasisSetCommon::Read(std::istream& is)
{
    UniqueID::Read(is);

    delete itsQuantumNumber;
    itsQuantumNumber=QuantumNumber::Factory(is);
    is >> *itsQuantumNumber;

    return is;
};

//-----------------------------------------------------------------------------
//
//  Controlled by derived class for IO.
//
std::ostream& IrrepBasisSetCommon::WriteBasisFunctions(std::ostream& os) const
{
    os << itsBasisFunctions;
    return os;
}

std::istream& IrrepBasisSetCommon::ReadBasisFunctions(std::istream& is)
{
    is >> itsBasisFunctions;
    return is;
}


