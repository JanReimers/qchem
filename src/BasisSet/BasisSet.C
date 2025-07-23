// File: BasisSet.C Quantum Chemistry basis set expressed as a sequence of Irrep basis sets.
module;
#include <vector>
#include <memory>

export module qchem.BasisSet;
export import qchem.Irrep_BS;
export import qchem.Fit_IBS;
export import qchem.Cluster;
export import qchem.Symmetry;
export import qchem.Symmetry.ElectronConfiguration;

import qchem.LAParams;
import Common.UniqueID;
import Common.Iterators;
import qchem.Streamable;

//----------------------------------------------------------------------------
//
//  Interface for a BasisSet which is assumed to a list of Irrep Basis Sets.
//
export class BasisSet
    : public virtual UniqueID
    , public virtual Streamable
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
