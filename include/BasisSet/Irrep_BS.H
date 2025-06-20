// File: Irrep_BS.H
#ifndef _Irrep_BS_H_
#define _Irrep_BS_H_

#include <omlfwd.H>
#include <Symmetry/fwd.H>
#include <LASolver/fwd.H>
#include <BasisSet/fwd.H>
#include <Common/UniqueID.H>
#include <Common/Iterators.H>
#include "Common/pmstream.h" 

#include <vector>
#include <memory>

//----------------------------------------------------------------------------
//
//  Interface for an irreducible representation basis sets.  H is block diagonal with one
//  block for  IrrepBasisSet,  For atoms each L get.s an IrrepBasisSet and an H  block. 
//  
//  The quantum number could be L for atoms, Irreducable rep for molecules, or
//  the wave vector k for solids.
//

class IrrepBasisSet
    : public virtual UniqueID
    , public virtual PMStreamableObject
{
public:
    typedef BasisFunction bf_t;
    typedef std::vector<std::shared_ptr<bf_t>> bfv_t;
    typedef bfv_t::      iterator       iterator;
    typedef bfv_t::const_iterator const_iterator;
    typedef std::shared_ptr<const Symmetry> sym_t;

    virtual void   Set(const LAParams&)=0;
    virtual size_t GetNumFunctions() const=0;
    virtual size_t size() const {return GetNumFunctions();}
    virtual sym_t  GetSymmetry() const=0;
//
//  PMStreamableObject stuff.
//
    virtual IrrepBasisSet* Clone  (const RVec3&) const=0;
    static  IrrepBasisSet* Factory(std::istream&    )        ;

private:
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
    virtual       iterator begin()      =0;
    virtual       iterator end  ()      =0;
    
public:
    // These facilutate iteration but return pointers dyn-cast'ed to derived types (D);
    template <class D> auto Iterate() const
    {
        return D_iterator_proxy<const D,const_iterator>(begin(),end());
    }
    template <class D> auto Iterate(const D* start) const
    {
        return D_iterator_proxy<const D,const_iterator>(begin(),end(),start);
    }
    template <class D> auto Iterate() 
    {
        return D_iterator_proxy<D,iterator>(begin(),end());
    }
    template <class D> auto Iterate(D* start) 
    {
        return D_iterator_proxy<D,iterator>(begin(),end(),start);
    }
    
    friend class BS_Common;
};

#include <Mesh/VectorFunction.H>
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

#include <BasisSet/Integrals.H>
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
    , public virtual Integrals_Kinetic<T> 
    , public virtual Integrals_Nuclear<T> 
{
    public:
    
};

#endif //_Irrep_BS_H_