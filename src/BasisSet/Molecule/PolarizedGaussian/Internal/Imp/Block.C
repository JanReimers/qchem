// File: Block.C  A block of basis functions with the same radial function.
module;
#include <iomanip>
#include "RadialFunction.H"
module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Block;
import qchem.stl_io;

namespace PolarizedGaussian
{

Block::Block()
    : itsRadial(0)
    , itsN     (0)
{};

Block::Block(RadialFunction* rf, size_t N)
    : itsRadial(rf)
    , itsN     (N)
{};

Block::Block(const Block& bfb)
    : itsRadial(bfb.itsRadial->Clone())
    , itsPols  (bfb.itsPols)
    , itsN     (bfb.itsN)
{};

Block::~Block()
{
    delete itsRadial;
};

size_t Block::LMax() const
{
    size_t lmax=0;
    for (auto p:itsPols) lmax=std::max(lmax,(size_t)p.GetTotalL());
    return lmax;
}

std::ostream& Block::Write(std::ostream& os) const
{
    os << "Block start:" << std::endl;
    os.setf(std::ios::fixed,std::ios::floatfield);
    os << std::setw(5) << std::setprecision(2) << itsRadial->GetCenter() << " "
    << *itsRadial << " " ;
    for (std::vector<Polarization>::const_iterator b(itsPols.begin()); b!=itsPols.end(); b++) os << *b;
    os << std::endl << "-------------------------------------------------------------------" << std::endl;
    return os;
}


Block* Block::Clone() const
{
    return new Block(*this);
}

Block* Block::Clone(const RVec3& newCenter) const
{
    RadialFunction* newRF=itsRadial->Clone(newCenter);
    Block* ret= new Block(newRF,itsN);
    for (std::vector<Polarization>::const_iterator b(itsPols.begin()); b!=itsPols.end(); b++) ret->Add(*b);
    return ret;
}

} //namespace PolarizedGaussian
