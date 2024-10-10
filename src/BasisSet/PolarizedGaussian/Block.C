// File: BasisFunctionBlock.C  A block of basis functions with the same radial function.



#include "Imp/BasisSet/PolarizedGaussian/Block.H"
#include "oml/imp/binio.h"
#include "oml/io3d.h"
#include "Misc/stl_io.h"
#include <iomanip>

BasisFunctionBlock::BasisFunctionBlock()
    : itsRadial(0)
    , itsN     (0)
{};

BasisFunctionBlock::BasisFunctionBlock(RadialFunction* rf, index_t N)
    : itsRadial(rf)
    , itsN     (N)
{};

BasisFunctionBlock::BasisFunctionBlock(const BasisFunctionBlock& bfb)
    : itsRadial(bfb.itsRadial->Clone())
    , itsPols  (bfb.itsPols)
    , itsN     (bfb.itsN)
{};

BasisFunctionBlock::~BasisFunctionBlock()
{
    delete itsRadial;
};

std::ostream& BasisFunctionBlock::Write(std::ostream& os) const
{
    if (!Pretty())
    {
        os << itsRadial << itsPols;
        if(Binary()) BinaryWrite(itsN,os);
        if(Ascii ()) os << itsN << " ";
    }
    else
    {
        os.setf(std::ios::fixed,std::ios::floatfield);
        os << std::setw(5) << std::setprecision(2) << itsRadial->GetCenter() << " "
        << *itsRadial << " " ;
        for (std::vector<Polarization>::const_iterator b(itsPols.begin()); b!=itsPols.end(); b++) os << *b;
        os << std::endl;
    }
    return os;
}

std::istream& BasisFunctionBlock::Read(std::istream& is)
{
    delete itsRadial;
    itsRadial=RadialFunction::Factory(is);
    is >> itsRadial >> itsPols;
    if(Binary()) BinaryRead(itsN,is);
    if(Ascii ()) is >> itsN;
    return is;
}

BasisFunctionBlock* BasisFunctionBlock::Clone() const
{
    return new BasisFunctionBlock(*this);
}

BasisFunctionBlock* BasisFunctionBlock::Clone(const RVec3& newCenter) const
{
    RadialFunction* newRF=itsRadial->Clone(newCenter);
    BasisFunctionBlock* ret= new BasisFunctionBlock(newRF,itsN);
    for (std::vector<Polarization>::const_iterator b(itsPols.begin()); b!=itsPols.end(); b++) ret->Add(*b);
    return ret;
}

