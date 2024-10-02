// File: BasisGroup.H  Abstract interface for a basis set.
#ifndef _BasisGroup_H_
#define _BasisGroup_H_

#include "BasisSet.H"
#include "Misc/ptr_vector1.h"
#include "Misc/ERIList.H"

//----------------------------------------------------------------------------
//
//  Interface for all basis sets.  All matrices of integrals are supplied by the
//  IntegralDataBase, which only calculates the integrals once and buffers them.
//  The Factory allows one to read basis sets from a stream polymorphically.
//  Clone(RVec&) creates a similar set with a new center.
//  The quantum number could be L for atoms, Irreducable rep for molecules, or
//  the wave vector k for solids.
//  Persistent is derived from StreamableObject which has a virtual destructor
//
class BasisGroup
    : public UniqueID
{
public:
    typedef optr_vector1<BasisSet*>::const_iterator const_iterator;

    BasisGroup();
    virtual ~BasisGroup();
    virtual size_t GetNumFunctions() const;
    virtual size_t GetNumBasisSets() const;
    virtual void Insert(BasisSet*);
    // Share ERI tables with other DBs
    virtual void Insert(const ERIList& C, const ERIList& X) const;

    virtual const_iterator begin() const {return itsBasisSets.begin();}
    virtual const_iterator end  () const {return itsBasisSets.end  ();}
//
//  StreamableObject stuff.
//
    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    static  BasisGroup*  Factory(std::istream&    )      ;

private:
    optr_vector1<BasisSet*> itsBasisSets;
};

#endif //_BasisGroup_H_
