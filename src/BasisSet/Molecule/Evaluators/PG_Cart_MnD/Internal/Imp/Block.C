// File: Block.C  A block of basis functions with the same radial function.
module;
#include <iomanip>
#include <vector>
module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Internal.Block;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.stl_io;

namespace qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD
{

Block::Block()
    : itsRadial(0)
{};

Block::Block(GaussianRF* rf)
    : itsRadial(rf)
{};

Block::Block(const Block& bfb)
    : itsRadial(new GaussianRF(*bfb.itsRadial))
    , itsPols  (bfb.itsPols)
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

Block* Block::Clone(const rvec3_t& newCenter) const
{
    GaussianRF* newRF=new GaussianRF(itsRadial->AtCenter(newCenter));
    Block* ret= new Block(newRF);
    for (std::vector<Polarization>::const_iterator b(itsPols.begin()); b!=itsPols.end(); b++) ret->Add(*b);
    return ret;
}

} //namespace qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD
