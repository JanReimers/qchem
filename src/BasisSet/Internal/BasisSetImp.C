// File: src/BasisSet/Internal/BasisSetImp.C  Common implementation for a full basis set.
module;
#include <memory>
#include <vector>
#include <cassert>
#include <variant>
#include <stdexcept>
#include "forward.H"
export module qchem.BasisSet.Internal.BasisSetImp;
export import qchem.Types;
export import qchem.BasisSet;
import qchem.stl_io;

export namespace qchem::BasisSet
{

//! \brief THE CHILD SLOT AT THE DECISION POINT (doc/RealComplexPlan.md §4 item 1, Step 3b): one irrep
//! basis block, typed by ITS OWN orbital scalar rather than the whole set's face -- the basis is where
//! the realness decision is COMPUTED (`irrep.IsReal()`), so its container is where per-block variation
//! originates; the composite WF/CD slots (2a/2b) downstream inherit it.  Same idiom as those: the
//! alternatives are the ABSTRACT per-block face, ownership included; the aggregation that is
//! scalar-independent visits with one generic lambda; the T-typed face (`GetIBS`, feeding `Iterate<D>`)
//! sees only the same-scalar children, with the cross arm throwing until Step 3c re-threads the
//! consumers (the Hamiltonian terms, MakeIrrepWFs) per block.
using ibs_child_t = std::variant<std::unique_ptr<Orbital_1E_IBS<double>>, std::unique_ptr<Orbital_1E_IBS<dcmplx>>>;

template <class T> class BasisSetImp
    : public virtual tBasisSet<T>
{
public:

    using obs_t = typename tBasisSet<T>::obs_t;

    virtual size_t GetNumFunctions() const
    {
        size_t ret=0;
        for (auto& c:itsBasisSets)
            ret+=std::visit([](const auto& b){return b->GetNumFunctions();}, c);
        return ret;
    }
    virtual irrepv_t GetIrreps(const Spin& ms) const
    {
        irrepv_t irrepv;
        for (auto& c:itsBasisSets)
            irrepv.push_back(std::visit([&](const auto& b){return b->GetIrrep(ms);}, c));
        return irrepv;
    }
//
//  StreamableObject stuff.
//
    virtual std::ostream&  Write(std::ostream& os) const
    {
        for (auto& c:itsBasisSets) std::visit([&](const auto& b){os << *b;}, c);
        return os;
    }

    //! Typed child view for MIXED-AWARE consumers (Step 3c) and the unit tests: the raw slot, both
    //! alternatives.  The tBasisSet<T> face's GetIBS below stays the same-scalar view.
    const ibs_child_t& GetChild(size_t i) const {return itsBasisSets[i];}

protected:
    //! TAKES OWNERSHIP.  TWO overloads, one per block scalar (Step 3b): the child slot is typed by the
    //! BLOCK, not by this set's face T, so a complex-faced periodic set may hold a real TRIM block.
    void Insert(Orbital_1E_IBS<double>* bs)
    {
        assert(bs);
        itsBasisSets.push_back(ibs_child_t(std::unique_ptr<Orbital_1E_IBS<double>>(bs)));
    }
    void Insert(Orbital_1E_IBS<dcmplx>* bs)
    {
        assert(bs);
        itsBasisSets.push_back(ibs_child_t(std::unique_ptr<Orbital_1E_IBS<dcmplx>>(bs)));
    }

    virtual size_t      GetNumIBS()      const {return itsBasisSets.size();}
    //! The tBasisSet<T> face's storage primitive: the SAME-SCALAR view (Iterate<D> casts from it).  A
    //! cross-scalar child THROWS -- loud in Release -- until Step 3c makes the consumers per-block-typed;
    //! until then nothing builds a mixed set on a shipped path.
    virtual const obs_t* GetIBS(size_t i) const
    {
        if (auto* p=std::get_if<std::unique_ptr<Orbital_1E_IBS<T>>>(&itsBasisSets[i])) return p->get();
        throw std::logic_error("tBasisSet: a basis block's scalar differs from the set's face -- the typed "
                               "Iterate view cannot reach it until RealComplexPlan Step 3c lands");
    }

    friend class ::SlaterRadialIntegralTests;
    friend class ::DiracIntegralTests;

    std::vector<ibs_child_t> itsBasisSets;
};

} //namespace
