// File: BasisSet.H  Abstract interface for a basis set.
#ifndef _BasisSet_H_
#define _BasisSet_H_

#include <types.H>
#include <UniqueID.H>
#include "Imp/Misc/pmstream.h"

#include <vector>
#include "IEClient.H"
#include "Imp/Containers/ptr_vector.h"
//----------------------------------------------------------------------------
//
//  Interface for a BasisSet which is assumed to a list of Irrep Basis Sets.
//
class BasisSet
    : public virtual UniqueID
    , public virtual PMStreamableObject
{
public:
    typedef Orbital_IBS bs_t;
    typedef TOrbital_IBS<double> tbs_t;
    typedef optr_vector1<bs_t*>::const_iterator const_iterator;
    
    BasisSet() {};
    virtual ~BasisSet() {}; 
    virtual void Set(const LAParams&)=0;
    virtual size_t GetNumFunctions() const=0;
    
    virtual Fit_IBS* CreateCDFitBasisSet(const Cluster* cl) const;
    virtual Fit_IBS* CreateVxcFitBasisSet(const Cluster* cl) const;
    static  BasisSet*  Factory(std::istream&    )      ;

private:
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;

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
    
public:

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
