// File: BasisSetImplementation.C  General basis set implementation as a list.



#include "Imp/BasisSet/IBS_Common.H"
#include <Symmetry.H>
#include <BasisFunction.H>
#include "Imp/Containers/stl_io.h"

#include <iomanip>

//-----------------------------------------------------------------------------
//
//  Construction zone
//
IBS_Common::IBS_Common()
    : itsSymmetry(0)
{
};

IBS_Common::IBS_Common(Symmetry* theQN)
    : itsSymmetry(theQN)
{
    assert(itsSymmetry);
};

IBS_Common::IBS_Common(const IBS_Common& bs)
  : itsSymmetry (bs.itsSymmetry->Clone())
  , itsBasisFunctions(bs.itsBasisFunctions)
  {
    assert(itsSymmetry);
  };

IBS_Common::~IBS_Common()
{
    delete itsSymmetry;
}

//-----------------------------------------------------------------------------
//
//  Post construction initializations called by dervied classes.
//
void IBS_Common::Insert(bf_t* bf)
{
    assert(bf);
    itsBasisFunctions.push_back(std::shared_ptr<bf_t>(bf));
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

//-----------------------------------------------------------------------------
//
//  Streamable stuff.
//
std::ostream& IBS_Common::Write(std::ostream& os) const
{
    assert(itsSymmetry);
    if(!Pretty())
    {
        UniqueIDImp::Write(os);
        if(!Binary()) os << std::endl;
        os << *itsSymmetry;
        if(!Binary()) os << std::endl;
    }
    else
    {
        os << "IrrepBasisSet " << " with " << GetNumFunctions() << " basis functions"
        << ", Quantum number=" << *itsSymmetry <<  std::endl;
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


