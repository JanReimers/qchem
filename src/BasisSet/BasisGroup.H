// File: BasisGroup.H  Abstract interface for a basis set.
#ifndef _BasisGroup_H_
#define _BasisGroup_H_

#include "BasisSet/BasisSet.H"
#include "Misc/ptr_vector.h"
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
    BasisGroup();
    virtual ~BasisGroup();
    virtual index_t GetNumFunctions() const;
    virtual index_t GetNumBasisSets() const;
    virtual void Insert(BasisSet*);
    // Share ERI tables with other DBs
    virtual void Insert(const ERIList& C, const ERIList& X) const;
//
//  StreamableObject stuff.
//
    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    static  BasisGroup*  Factory(std::istream&    )      ;
private:
    typedef optr_vector<BasisSet*>::const_iterator CITER;
    typedef optr_vector<BasisSet*>::iterator ITER;
    friend class BasisGroupBrowser;
    friend class BasisGroupIterator;
    virtual CITER begin() const {return itsBasisSets.begin();}
    virtual CITER end  () const {return itsBasisSets.end  ();}
    virtual ITER  begin()       {return itsBasisSets.begin();}
    virtual ITER  end  ()       {return itsBasisSets.end  ();}

    optr_vector<BasisSet*> itsBasisSets;
};

#endif //_BasisGroup_H_
