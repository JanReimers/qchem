// File: BasisSet.H  Abstract interface for a basis set.
#ifndef _BasisSet_H_
#define _BasisSet_H_

#include <vector>
#include <memory>

#include <Symmetry/fwd.H>
#include <LASolver/fwd.H>
#include <BasisSet/fwd.H>
#include "Common/pmstream.h"
#include <LASolver/LAParams.H>


import Common.UniqueID;
import Common.Iterators;
import qchem.Cluster;
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
    typedef std::vector<std::unique_ptr<bs_t>> bsv_t;
    typedef bsv_t::const_iterator const_iterator;
    typedef std::shared_ptr<const Symmetry> sym_t;
    typedef std::vector<sym_t> symv_t;
    
    BasisSet() {};
    virtual ~BasisSet() {}; 
    virtual void Set(const LAParams&)=0;
    virtual size_t GetNumFunctions() const=0;
    virtual symv_t GetSymmetries() const=0;
    
    virtual Fit_IBS* CreateCDFitBasisSet(const Cluster* cl) const;
    virtual Fit_IBS* CreateVxcFitBasisSet(const Cluster* cl) const;
    static  BasisSet*  Factory(std::istream&    )      ;

private:
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;

    
public:
    // These facilutate iteration but return pointers dyn-cast'ed to derived types (T);
    template <class D> auto Iterate() const
    {
        return D_iterator_proxy<const D,const_iterator>(begin(),end());
    }
    template <class D> auto Iterate(const D* start) const
    {
        return D_iterator_proxy<const D,const_iterator>(begin(),end(),start);
    }

};


#endif //_BasisSet_H_
