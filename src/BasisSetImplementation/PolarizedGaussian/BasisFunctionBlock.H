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

class GaussianRF;
struct RadialFunctionPair
{
    RadialFunctionPair() : itsIndex(-1){};
    void push_back(const GaussianRF* rf)
    {
        assert(rf);
        itsIndex++;
        assert(itsIndex<2);
        itsRFs[itsIndex]=rf;
    }

    void pop_back()
    {
        itsRFs[itsIndex]=0;
        itsIndex--;
    }
    const GaussianRF* a() const {return itsRFs[0];}
    const GaussianRF* b() const {return itsRFs[1];}

    const GaussianRF* itsRFs[2];
    int itsIndex;
};

struct RadialFunctionQuartet
{
    RadialFunctionQuartet() : itsIndex(-1){};
    void push_back(const GaussianRF* rf)
    {
        assert(rf);
        itsIndex++;
        assert(itsIndex<4);
        itsRFs[itsIndex]=rf;
    }

    void pop_back()
    {
        itsRFs[itsIndex]=0;
        itsIndex--;
    }
    const GaussianRF* a() const {return itsRFs[0];}
    const GaussianRF* b() const {return itsRFs[1];}
    const GaussianRF* c() const {return itsRFs[2];}
    const GaussianRF* d() const {return itsRFs[3];}

    const GaussianRF* itsRFs[4];
    int itsIndex;
};

struct RadialFunctionTriplet
{
    RadialFunctionTriplet() : itsIndex(-1){};
    void push_back(const GaussianRF* rf)
    {
        assert(rf);
        itsIndex++;
        assert(itsIndex<3);
        itsRFs[itsIndex]=rf;
    }

    void pop_back()
    {
        itsRFs[itsIndex]=0;
        itsIndex--;
    }
    const GaussianRF* a() const {return itsRFs[0];}
    const GaussianRF* b() const {return itsRFs[1];}
    const GaussianRF* c() const {return itsRFs[2];}

    const GaussianRF* itsRFs[3];
    int itsIndex;
};

// Syntactic sugar for function arguments.
struct BasisFunctionBlockPair
{
    typedef const BasisFunctionBlock BFB;
    BasisFunctionBlockPair(BFB* _a,BFB* _b) : a(_a), b(_b) {};
    BFB* a;
    BFB* b;
    RadialFunctionPair radials; //These are sorted and guaranteed to be primative.
};

struct BasisFunctionBlockTriplet
{
    typedef const BasisFunctionBlock BFB;
    BasisFunctionBlockTriplet(BFB* _a,BFB* _b, BFB* _c) : a(_a), b(_b), c(_c) {};
    BFB* a;
    BFB* b;
    BFB* c;
    RadialFunctionTriplet radials; //These are sorted and guaranteed to be primative.
};

struct BasisFunctionBlockQuartet
{
    typedef const BasisFunctionBlock BFB;
    BasisFunctionBlockQuartet(BFB* _a,BFB* _b, BFB* _c, BFB* _d) : a(_a), b(_b), c(_c), d(_d) {};
    BFB* a;
    BFB* b;
    BFB* c;
    BFB* d;
    RadialFunctionQuartet radials; //These are sorted and guaranteed to be primative.
};


#endif //_BasisFunctionBlock_H_
