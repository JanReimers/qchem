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
    //typedef dynamic_cast_iterator<bs_t*,tbs_t* > const_Titerator;

    BasisSet() {};

    virtual ~BasisSet() {}; 
    virtual size_t GetNumFunctions() const=0;
//    virtual size_t GetNumIrreps() const=0;


    virtual UniqueID::IDtype GetID() const {return UniqueID::GetID();}
    
    virtual Fit_IBS* CreateCDFitBasisSet(const Cluster* cl) const;
    virtual Fit_IBS* CreateVxcFitBasisSet(const Cluster* cl) const;
    static  BasisSet*  Factory(std::istream&    )      ;

private:
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
public:

    template <class T> class BS_iterator
    {
        typedef BasisSet::const_iterator it_t;
    public:
        BS_iterator(const it_t& b) : current(b) {};
        BS_iterator(const it_t& b, const bs_t* c) : current(b) //STL won't let us construct current(c).
        {
            while (*current!=c) 
                ++current; //Clunky
        };
        it_t operator++() {return ++current;} //Prefix only.
        const T* operator*() const
        {
            const T* ret(dynamic_cast<const T*>(*current));
            assert(ret);
            return ret;
        }
        friend bool operator!=(const BS_iterator& a, const BS_iterator& b)
        {
            return a.current!=b.current;
        }
    private:
        it_t current;
    };
    template <class T> class BS_iterator_proxy
    {
        typedef BS_iterator<T> it_t;
        typedef BasisSet::const_iterator cit_t;
    public:
        BS_iterator_proxy(const cit_t& b, const cit_t& e) : ib(b), ie(e) {};
        BS_iterator_proxy(const cit_t& b, const cit_t& e, const T* start) : ib(b,start), ie(e) {};
        it_t begin() const {return ib;}
        it_t end  () const {return ie;}
    private:
        it_t ib,ie;
    };

    template <class T> BS_iterator_proxy<T> Iterate() const
    {
        return BS_iterator_proxy<T>(begin(),this->end());
    }
    template <class T> BS_iterator_proxy<T> Iterate(const T* start) const
    {
        return BS_iterator_proxy<T>(begin(),this->end(),start);
    }

};


#endif //_BasisSet_H_
