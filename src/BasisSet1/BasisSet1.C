// File: BasisSet.C Quantum Chemistry basis set expressed as a sequence of Irrep basis sets.
module;
#include <vector>
#include <memory>

export module qchem.BasisSet1;
export import qchem.Orbital_1E_IBS1;
export import qchem.Fit_IBS;
export import qchem.Cluster;
export import qchem.Symmetry;
export import qchem.Symmetry.ElectronConfiguration;

import Common.UniqueID;
import Common.Iterators;
import qchem.Streamable;

//----------------------------------------------------------------------------
//
//  Interface for a BasisSet which is assumed to a list of Irrep Basis Sets.
//
export class BasisSet1
    : public virtual UniqueID
    , public virtual Streamable
{
public:
    typedef Orbital_1E_IBS1<double> bs_t; 
    typedef std::vector<std::unique_ptr<bs_t>> bsv_t;
    typedef bsv_t::const_iterator const_iterator;
    typedef std::vector<Irrep_QNs> irrepv_t; 
    
    BasisSet1() {};
    virtual ~BasisSet1() {}; 
    virtual size_t GetNumFunctions() const=0;
    virtual irrepv_t GetIrreps(const Spin& ms) const=0;
    
    // virtual Fit_IBS* CreateCDFitBasisSet(const Cluster* cl) const;
    // virtual Fit_IBS* CreateVxcFitBasisSet(const Cluster* cl) const;

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
