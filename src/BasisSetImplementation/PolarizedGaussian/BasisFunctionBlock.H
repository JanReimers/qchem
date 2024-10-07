// File: BasisFunctionBlock.H  A block of basis functions with the same radial function.
#ifndef _BasisFunctionBlock_H_
#define _BasisFunctionBlock_H_



#include "Misc/pmstream.h"
#include "BasisSetImplementation/PolarizedGaussian/RadialFunction.H"
#include "Misc/Polarization.H"
#include <cassert>

//-----------------------------------------------------------------------
//
//  This class or structure represents a group of basis functions
//  with differing polarizations and the same radial part.
//
class BasisFunctionBlock
    : public virtual PMStreamableObject
{
public:
    BasisFunctionBlock(                         );
    BasisFunctionBlock(RadialFunction*, index_t );
    BasisFunctionBlock(const BasisFunctionBlock&);
    ~BasisFunctionBlock(                         );

    void Add(const Polarization& p)
    {
        itsPols.push_back(p);
    }
    index_t size() const
    {
        return itsPols.size();
    }

    virtual std::ostream&       Write  (std::ostream&) const;
    virtual std::istream&       Read   (std::istream&)      ;
    virtual BasisFunctionBlock* Clone  (             ) const;
    virtual BasisFunctionBlock* Clone  (const RVec3& ) const;
    static  BasisFunctionBlock* Factory(std::istream&)      ;

    RadialFunction*           itsRadial; //Common radial function.
    std::vector<Polarization> itsPols;   //All polarizations for this block.
    index_t                   itsN;      //Index of first basis function in block.
};

#endif //_BasisFunctionBlock_H_
