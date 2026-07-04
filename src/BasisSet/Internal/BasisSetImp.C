// File: src/BasisSet/Internal/BasisSetImp.C  Common implementation for a full basis set.
module;
#include <memory>
#include <vector>
#include <cassert>
#include "forward.H"
export module qchem.BasisSet.Internal.BasisSetImp;
export import qchem.Types;
export import qchem.BasisSet;
import qchem.stl_io;

export namespace qchem::BasisSet
{

template <class T> class BasisSetImp
    : public virtual tBasisSet<T>
{
public:
    
    using obs_t = typename tBasisSet<T>::obs_t;

    virtual size_t GetNumFunctions() const
    {
        size_t ret=0;
        for (auto& bs:itsBasisSets)
            ret+=bs->GetNumFunctions();
        return ret;
    }
    virtual irrepv_t GetIrreps(const Spin& ms) const
    {
        irrepv_t irrepv;
        for (auto& b:itsBasisSets) irrepv.push_back(b->GetIrrep(ms));
        return irrepv;
    }
//
//  StreamableObject stuff.
//
    virtual std::ostream&  Write(std::ostream& os) const
    {
        return os << itsBasisSets;
    }
protected:
    void Insert(Orbital_1E_IBS<T>* bs)
    {
        assert(bs);
        itsBasisSets.push_back(std::unique_ptr<obs_t>(bs));
    }

    virtual size_t      GetNumIBS()      const {return itsBasisSets.size();}
    virtual const obs_t* GetIBS(size_t i) const {return itsBasisSets[i].get();}

    friend class ::SlaterRadialIntegralTests;
    friend class ::DiracIntegralTests;

    std::vector<std::unique_ptr<obs_t>> itsBasisSets;
};

} //namespace
