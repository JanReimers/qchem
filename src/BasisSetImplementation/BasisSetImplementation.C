// File: BasisSetImplementation.C  General basis set implementation as a list.



#include "BasisSetImplementation/BasisSetImplementation.H"
#include "BasisSet/QuantumNumber.H"
#include "Misc/ptr_vector1_io.h"

#include <iostream>
#include <iomanip>
#include <cassert>

//-----------------------------------------------------------------------------
//
//  Construction zone
//
BasisSetImplementation::BasisSetImplementation()
    : itsQuantumNumber(0)
    , itsStartIndex(0)
{
};

BasisSetImplementation::BasisSetImplementation(QuantumNumber* theQN)
    : itsQuantumNumber(theQN)
    , itsStartIndex(0)
{
};

BasisSetImplementation::BasisSetImplementation(const BasisSetImplementation& bs)
  : itsQuantumNumber (bs.itsQuantumNumber->Clone())
  , itsBasisFunctions(bs.itsBasisFunctions)
  , itsStartIndex(bs.itsStartIndex)
  {
    assert(itsQuantumNumber);
  };

BasisSetImplementation::~BasisSetImplementation()
{
    delete itsQuantumNumber;
}

//-----------------------------------------------------------------------------
//
//  Post construction initializations called by dervied classes.
//
void BasisSetImplementation::Insert(BasisFunction* bf)
{
    assert(bf);
    itsBasisFunctions.push_back(bf);
}

void BasisSetImplementation::EmptyBasisFunctions()
{
    itsBasisFunctions.clear();
}


//-----------------------------------------------------------------------------
//
//  Basis Set Stuff.
//
void BasisSetImplementation::SetStartIndex(int si)
{
    assert(si>0);
    itsStartIndex=si;
}

int  BasisSetImplementation::GetStartIndex() const
{
    assert(itsStartIndex>0);
    return itsStartIndex;
}


size_t BasisSetImplementation::GetNumFunctions() const
{
    return itsBasisFunctions.size();
}

bool BasisSetImplementation::operator==(const BasisSet& bs) const
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
std::ostream& BasisSetImplementation::Write(std::ostream& os) const
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

std::istream& BasisSetImplementation::Read(std::istream& is)
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
std::ostream& BasisSetImplementation::WriteBasisFunctions(std::ostream& os) const
{
    os << itsBasisFunctions;
    return os;
}

std::istream& BasisSetImplementation::ReadBasisFunctions(std::istream& is)
{
    is >> itsBasisFunctions;
    return is;
}


