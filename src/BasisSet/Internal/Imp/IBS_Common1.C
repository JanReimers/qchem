// File: Imp/IBS_Common.C  Irrep Basis set common implementation.
module;
#include <iosfwd>
#include <iomanip>
#include <memory>
#include <cassert>
module qchem.BasisSet.Internal.IBS_Common1;
import qchem.BasisFunction;
import qchem.Symmetry;
import qchem.stl_io;

//-----------------------------------------------------------------------------
//
//  Construction zone
//
IBS_Common1::IBS_Common1()
    : itsSymmetry(0)
{
};

IBS_Common1::IBS_Common1(Symmetry* theQN)
    : itsSymmetry(theQN)
{
    assert(itsSymmetry);
};

IBS_Common1::IBS_Common1(const IBS_Common1& bs)
  : itsSymmetry (bs.itsSymmetry)
  , itsBasisFunctions(bs.itsBasisFunctions)
  {
    assert(itsSymmetry);
  };

IBS_Common1::~IBS_Common1()
{
}

//-----------------------------------------------------------------------------
//
//  Post construction initializations called by dervied classes.
//
void IBS_Common1::Insert(bf_t* bf)
{
    assert(bf);
    itsBasisFunctions.push_back(std::shared_ptr<bf_t>(bf));
}

void IBS_Common1::EmptyBasisFunctions()
{
    itsBasisFunctions.clear();
}


//-----------------------------------------------------------------------------
//
//  Basis Set Stuff.
//


size_t IBS_Common1::GetNumFunctions() const
{
    return itsBasisFunctions.size();
}

//-----------------------------------------------------------------------------
//
//  Streamable stuff.
//
std::ostream& IBS_Common1::Write(std::ostream& os) const
{
    assert(itsSymmetry);
    os << "IrrepBasisSet " << " with " << GetNumFunctions() << " basis functions"
    << ", Quantum number=" << *itsSymmetry <<  std::endl;
    os << "   #          Center      Pol     Radial     Exponents" << std::endl;
    // No UT coverage
    int i=1;
    for(auto b:*this) os << std::setw(4) << i++ << " " << *b << std::endl;
    
    os << "-------------------------------------------------------------------------------" << std::endl;

    return os;
}

//-----------------------------------------------------------------------------
//
//  Controlled by derived class for IO.
//
std::ostream& IBS_Common1::WriteBasisFunctions(std::ostream& os) const
{
    os << itsBasisFunctions;
    return os;
}

std::istream& IBS_Common1::ReadBasisFunctions(std::istream& is)
{
    // is >> itsBasisFunctions;
    return is;
}


