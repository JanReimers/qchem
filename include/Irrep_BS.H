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
//  All matrices of integrals are supplied by the
//  IntegralDataBase, which only calculates the integrals once and buffers them.
//  The Factory allows one to read basis sets from a stream polymorphically.
//  Clone(RVec&) creates a similar set with a new center.
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
    virtual const QuantumNumber& GetQuantumNumber(               ) const=0;
    virtual       QuantumNumber* GetQuantumNumber(int index      ) const=0;
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
    friend class BasisSetImp;
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
class Fit_IBS;
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

template <class T> class TOrbital_DFT_IBS
    : public virtual TOrbital_IBS<T>
    , public virtual Integrals_DFT<T> //DFT integrals
    
{
    typedef typename Integrals_Base<T>::SMat SMat;
    typedef typename Integrals_Base<T>::fbs_t fbs_t;
    typedef Vector<T> Vec;
public:
    virtual fbs_t*    CreateCDFitBasisSet(const BasisSet*,const Cluster*) const=0;
    virtual fbs_t*    CreateVxcFitBasisSet(const BasisSet*,const Cluster*) const=0;
    using Integrals_DFT<T>::Overlap3C; //Unhide
    using Integrals_DFT<T>::Repulsion3C; //Unhide
    virtual Vec Overlap3C  (const SMat& Dcd, const fbs_t* ff) const=0;
    virtual Vec Repulsion3C(const SMat& Dcd, const fbs_t* ff) const=0;
};

template <class T> class TOrbital_HF_IBS
    : public virtual TOrbital_IBS<T>
    , public virtual Integrals_2E<T> //Two electron integrals used for HF
{
    typedef typename Integrals_Base<T>::SMat SMat;
    typedef typename Integrals_Base<T>::obs_t obs_t;
public:
    using Integrals_2E<T>::Direct;
    using Integrals_2E<T>::Exchange;
    virtual SMat Direct  (const SMat& Dcd, const obs_t* bs_cd) const=0;
    virtual SMat Exchange(const SMat& Dcd, const obs_t* bs_cd) const=0;

};

template <class T> class Orbital_RKB_IBS
    : public virtual TIrrepBasisSet<T>
    , public virtual Integrals_RKB<T> 
{
public:
    
};

template <class T> class Orbital_RKBL_IBS
    : public virtual TIrrepBasisSet<T>
    , public virtual Integrals_RKBL<T> //One electron integrals used for everything
{
public:
    //int GetKappa() const;
};

template <class T> class Orbital_RKBS_IBS
    : public virtual TIrrepBasisSet<T>
    , public virtual Integrals_RKBS<T> //One electron integrals used for everything
{
public:
    //int GetKappa() const;
};

class Fit_IBS
    : public virtual TIrrepBasisSet<double>
    , public virtual FitIntegrals_1E 
{
    public:
    virtual Fit_IBS* Clone  (const RVec3&) const=0;  

};




#endif //_Irrep_BS_H_