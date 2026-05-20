// File: BasisSet.C Quantum Chemistry basis set expressed as a sequence of Irrep basis sets.
module;
#include <vector>
#include <memory>

export module qchem.BasisSet;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Fit_IBS;
export import qchem.Cluster;
export import qchem.Symmetry;
export import qchem.Symmetry.ElectronConfiguration;

import Common.Iterators;
export import qchem.Streamable;

export namespace BasisSet
{
typedef std::vector<Irrep_QNs> irrepv_t; 

//----------------------------------------------------------------------------
//
//  Interface for a BasisSet which is assumed to a list of Irrep Basis Sets.
//
template <class T> class BasisSet
    : public virtual Streamable
{
public:
    typedef Orbital_1E_IBS<T> bs_t; 
    typedef std::vector<std::unique_ptr<bs_t>> bsv_t;
    typedef bsv_t::const_iterator const_iterator;
    
    virtual ~BasisSet() {}; 
    virtual size_t   GetNumFunctions() const=0;
    virtual irrepv_t GetIrreps(const Spin& ms) const=0;
    
    virtual Fit_IBS* CreateCDFitBasisSet(const Cluster* cl) const;
    virtual Fit_IBS* CreateVxcFitBasisSet(const Cluster* cl) const;

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

typedef BasisSet<double>    Real_BS;
typedef BasisSet<dcmplx> Complex_BS;

}//namespace