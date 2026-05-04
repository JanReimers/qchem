// File: IBS_Common.C  Irrep Basis set common implementation.
module;
#include <iosfwd>
#include <memory>
#include <cassert>
#include "forward.H"
export module qchem.BasisSet1.Internal.Common;
export import qchem.Types;
export import qchem.BasisSet1;
import qchem.stl_io;

export namespace BasisSet1
{

template <class T> class BS_Common
    : public virtual BasisSet<T>
{
public:
    
    virtual size_t GetNumFunctions() const
    {
        size_t ret=0;
        for (auto& bs:*this) 
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
        itsBasisSets.push_back(std::unique_ptr<typename BasisSet<T>::bs_t>(bs));
    }

    using const_iterator=BasisSet<T>::const_iterator;
    virtual const_iterator begin() const {return itsBasisSets.begin();}
    virtual const_iterator end  () const {return itsBasisSets.end  ();}
    friend class SlaterRadialIntegralTests;
    friend class DiracIntegralTests;

    BasisSet<T>::bsv_t itsBasisSets;
};

} //namespace
