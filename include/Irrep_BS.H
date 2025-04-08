// File: Irrep_BS.H
#ifndef _Irrep_BS_H_
#define _Irrep_BS_H_

#include <IEClient.H>
#include <types.H>
#include <UniqueID.H>
#include "Imp/Misc/pmstream.h" 
#include <vector>
#include <cassert>
#include <iosfwd>

//----------------------------------------------------------------------------
//
//  Interface for an irreducible representation basis sets.  H is block diagonal with one
//  block for  IrrepBasisSet,  For atoms each L get.s an IrrepBasisSet and an H  block. 
//  
//  The quantum number could be L for atoms, Irreducable rep for molecules, or
//  the wave vector k for solids.
//

class IrrepBasisSet
    : public virtual IrrepIEClient
    , public virtual UniqueID
    , public virtual PMStreamableObject
{
public:
    typedef std::vector<const BasisFunction*>::      iterator       iterator;
    typedef std::vector<const BasisFunction*>::const_iterator const_iterator;

    virtual void Set(const LAParams&)=0;
    virtual bool                 operator==      (const IrrepBasisSet&) const=0;
    virtual size_t               GetNumFunctions (               ) const=0;
    virtual const Symmetry& GetSymmetry(               ) const=0;
private:
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
    
    template <class T> class IBS_iterator
    {
        typedef IrrepBasisSet::const_iterator it_t;
    public:
        IBS_iterator(const it_t& b) : current(b) {};
        IBS_iterator(const it_t& b, const BasisFunction* c) : current(b) //STL won't let us construct current(c).
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
        friend bool operator!=(const IBS_iterator& a, const IBS_iterator& b)
        {
            return a.current!=b.current;
        }
    private:
        it_t current;
    };
    template <class T> class IBS_iterator_proxy
    {
        typedef IBS_iterator<T> it_t;
        typedef IrrepBasisSet::const_iterator cit_t;
    public:
        IBS_iterator_proxy(const cit_t& b, const cit_t& e) : ib(b), ie(e) {};
        IBS_iterator_proxy(const cit_t& b, const cit_t& e, const T* start) : ib(b,start), ie(e) {};
        it_t begin() const {return ib;}
        it_t end  () const {return ie;}
    private:
        it_t ib,ie;
    };
    
public:
    template <class T> IBS_iterator_proxy<T> Iterate() const
    {
        return IBS_iterator_proxy<T>(begin(),this->end());
    }
    template <class T> IBS_iterator_proxy<T> Iterate(const T* start) const
    {
        return IBS_iterator_proxy<T>(begin(),this->end(),start);
    }
    
//
//  StreamableObject stuff.
//
    virtual IrrepBasisSet* Clone  (const RVec3&) const=0;
    static  IrrepBasisSet* Factory(std::istream&    )        ;
private:
    friend class BS_Common;
};

#include <VectorFunction.H>
#include "Imp/Containers/ptr_vector.h"
//----------------------------------------------------------------------------
//
//  Extend basis to be a set of real or complex valued functions
//
template <class T> class TIrrepBasisSet
    : public virtual IrrepBasisSet
    , public virtual VectorFunction<T>
{
public:
    size_t GetVectorSize() const {return GetNumFunctions();}
};

#include "AnalyticIE.H"
//
// Define an orbital irrep basis set which supports integrals for SCF orbital calculations.
// Mix-in the integral interfaces required for an orbital basis. 
//
class Orbital_IBS
    : public virtual IrrepBasisSet
{
    public:
    virtual LASolver<double>* CreateSolver() const=0;
    
};

template <class T> class TOrbital_IBS
    : public virtual Orbital_IBS
    , public virtual TIrrepBasisSet<T>
    , public virtual Integrals_Overlap<T> 
    , public virtual Integrals_Grad2<T> 
    , public virtual Integrals_Nuclear<T> 
{
    public:
    
};

#endif //_Irrep_BS_H_