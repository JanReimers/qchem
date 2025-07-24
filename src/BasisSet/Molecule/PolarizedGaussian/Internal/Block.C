// File: Block.C  A block of basis functions with the same radial function.
module;
#include <iosfwd>
#include "RadialFunction.H"
#include <cassert>
export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Block;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;
import qchem.Streamable;

export namespace PolarizedGaussian
{

//-----------------------------------------------------------------------
//
//  This class or structure represents a group of basis functions
//  with differing polarizations and the same radial part.
//
class Block
    : public virtual Streamable
{
public:
    Block(                         );
    Block(RadialFunction*, size_t  );
    Block(const Block&);
    ~Block(                         );

    void Add(const Polarization& p)
    {
        itsPols.push_back(p);
    }
    size_t  size() const
    {
        return itsPols.size();
    }
    size_t LMax() const;

    virtual std::ostream&       Write  (std::ostream&) const;
    virtual Block* Clone  (             ) const;
    virtual Block* Clone  (const RVec3& ) const;

    RadialFunction*           itsRadial; //Common radial function.
    std::vector<Polarization> itsPols;   //All polarizations for this block.
    size_t                    itsN;      //Index of first basis function in block.
};

} //namespace PolarizedGaussian

