// File: Block.H  A block of basis functions with the same radial function.
#ifndef _Block_H_
#define _Block_H_

#include "RadialFunction.H"
#include "Polarization.H"
#include "Common/pmstream.h"
#include <cassert>

namespace PolarizedGaussian
{

//-----------------------------------------------------------------------
//
//  This class or structure represents a group of basis functions
//  with differing polarizations and the same radial part.
//
class Block
    : public virtual PMStreamableObject
{
public:
    Block(                         );
    Block(RadialFunction*, index_t );
    Block(const Block&);
    ~Block(                         );

    void Add(const Polarization& p)
    {
        itsPols.push_back(p);
    }
    index_t size() const
    {
        return itsPols.size();
    }
    size_t LMax() const;

    virtual std::ostream&       Write  (std::ostream&) const;
    virtual Block* Clone  (             ) const;
    virtual Block* Clone  (const RVec3& ) const;
    static  Block* Factory(std::istream&)      ;

    RadialFunction*           itsRadial; //Common radial function.
    std::vector<Polarization> itsPols;   //All polarizations for this block.
    index_t                   itsN;      //Index of first basis function in block.
};

} //namespace PolarizedGaussian

#endif //_Block_H_
