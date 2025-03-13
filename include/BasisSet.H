// File: BasisSet.H  Abstract interface for a basis set.
#ifndef _BasisSet_H_
#define _BasisSet_H_

#include <types.H>
#include "Imp/Misc/UniqueID/UniqueID.H"

#include <vector>
#include "IEClient.H"
#include "Imp/Containers/ptr_vector.h"
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
class BasisSet
    : public UniqueID
{
public:
    typedef Orbital_IBS bs_t;
    typedef TOrbital_IBS<double> tbs_t;
    typedef optr_vector1<bs_t*>::const_iterator const_iterator;
    typedef dynamic_cast_iterator<bs_t*,tbs_t* > const_Titerator;

    BasisSet() {};

    virtual ~BasisSet() {}; 
    virtual size_t GetNumFunctions() const=0;
//    virtual size_t GetNumIrreps() const=0;

    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;

    virtual UniqueID::IDtype GetID() const {return UniqueID::GetID();}
    const_Titerator beginT() const {return const_Titerator(begin());};
    
    virtual Fit_IBS* CreateCDFitBasisSet(const Cluster* cl) const;
    virtual Fit_IBS* CreateVxcFitBasisSet(const Cluster* cl) const;
    static  BasisSet*  Factory(std::istream&    )      ;
};


#endif //_BasisSet_H_
