// File: HamiltonianTerm.C  General implementation of a HamiltonianTerm term in the Hamiltonian.
//
// Templated on the matrix element type T (double for atoms/molecules; dcmplx for the plane-wave
// lattice lineage).  hmat_t<double> IS rsmat_t and tobs_t<double> IS robs_t, so the <double> aliases
// at the bottom leave existing real code unchanged.  The cache-lookup bodies are inline here (not in a
// separate Imp unit) because the dcmplx terms instantiate these templates in other translation units
// and so need the definitions visible.
module;
#include <stdexcept>
#include <map>
#include <cassert>
#include <cstddef>
#include <type_traits>   // std::conditional_t/is_same_v (StaticRealBlockBase -- Step 3c)
export module qchem.Hamiltonian.Internal.Term;
export import qchem.Hamiltonian;
import qchem.Hamiltonian.Types;

import qchem.Symmetry.Irrep;

export namespace qchem::Hamiltonian
{

//! \brief WHICH quantum numbers a term's matrix actually DEPENDS ON -- and therefore how many times the
//! term is EVALUATED, which is the point of it.
//!
//! Every cache below keys on \c Irrep, i.e. on (spatial symmetry, SPIN).  For a spin-DEPENDENT term that
//! is exactly right.  For a term whose matrix does not depend on the channel it is a doubled evaluation:
//! the up and down calls miss each other and build the identical matrix twice.
//!
//! ⛔ THAT IS NOT A CACHE INEFFICIENCY, IT IS A DOUBLED BUILD OF THE MOST EXPENSIVE OBJECT IN THE RUN.
//! Measured 2026-09-04 on the MnO parity row: \c Vee_Hartree -- whose \c MakeMatrixT takes an UNNAMED
//! \c Spin&, i.e. provably ignores it -- was gathering ~2.1 duplicate KS matrices per SCF iteration, and
//! the integrate-back is 71% of that run (doc/OpenWork.md bin 1).
//!
//! \note IT ASKS WHAT THE ANSWER DEPENDS ON, NOT WHAT KIND OF TERM THIS IS.  A predicate spelled
//! \c IsSpinIndependent() would be the IDENTITY question this codebase rejects -- the caller would branch
//! on the term's type, and a third answer would mean a third branch (doc/OldPlans/ScreeningPlan.md §5).  Overriding
//! the KEY instead changes nothing at any call site: the caches go on keying, on one fewer quantum number.
class HT_SpinDependence
{
public:
    virtual ~HT_SpinDependence() = default;
    //! The spin this term's matrix depends on.  Default \a s (spin-resolved, the safe answer).
    //! Override to \c Spin::None when the matrix is the SAME for every channel.
    virtual Spin CacheSpin(const Spin& s) const {return s;}
};

template <class T> class tHT_Common
    : public virtual HT_SpinDependence
{
protected:
    typedef std::map<Irrep,hmat_t<T>> CacheMap;
    mutable CacheMap   itsCache;       //Cache the H matrices for total energy calculations.

    //! \brief CREATE (empty) one slot per irrep of \a bs, so the block loop only fills existing nodes.
    //! R1.0h; the rationale is on \c HT_SlotOwner::PrepareSlots.
    //!
    //! ★ AN EMPTY MATRIX IS THE "NOT FILLED YET" SENTINEL, and it is a safe one: a Fock/KS block always has
    //! rows, so 0x0 cannot be a legitimate cached value.  That is what lets a pre-created slot be
    //! distinguished from a filled one WITHOUT changing the map's value type (an `optional` would have
    //! rippled through five cache holders for the same information).
    //!
    //! All three spins, because \c CacheSpin decides which of them actually key distinct slots and this
    //! layer must not second-guess it: a spin-independent term folds all three onto one key (duplicates
    //! collapse), a polarized one keys Up and Down apart.  A few empty map nodes is the whole cost.
    void PrepareCacheSlots(const tbs_t<T>* bs) const
    {
        if (!bs) return;
        for (size_t i=0;i<bs->GetNumIBS();++i)
        {
            // ⛔ SKIP THE REAL BLOCKS OF A COMPLEX-FACED SET -- "walk the blocks" is NOT one loop on a
            // MIXED set (doc/RealComplexPlan.md 3c-3), and `operator[]` says so by THROWING: *"a basis
            // block's scalar differs from the set's face"*.  Those blocks are served by the real cache and
            // its own PrepareRealSlots; reaching them through the typed view is the error, not the skip.
            // (Measured the hard way 2026-09-09: 32 integration tests, every one a mixed-mesh run.)
            if (bs->GetRealIBS(i)) continue;
            if (const tobs_t<T>* b=(*bs)[i])
                for (const Spin& s : {Spin::None, Spin::Up, Spin::Down})
                    itsCache[b->GetIrrep(this->CacheSpin(s))];   // default-construct: 0x0 == unfilled
        }
    }
};


template <class T> class tStatic_HT_Imp
    : public virtual tStatic_HT<T>
    , protected tHT_Common<T>
{
public:
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>* bs,const Spin& s) const
    {
        assert(bs);
        Irrep qns(bs->GetIrrep(this->CacheSpin(s)));
        // FILL-IF-EMPTY (R1.0h).  operator[] mutates the tree ONLY when the key is absent, which the
        // pre-create phase is there to prevent; a 0x0 matrix is the "slot exists, not filled" sentinel.
        // A static term's cache is never cleared, so after the first iteration this is pure read.
        hmat_t<T>& slot=this->itsCache[qns];
        if (slot.rows()==0) slot=MakeMatrix(bs,s);
        return slot;
    }

    //! \copydoc HT_SlotOwner::PrepareSlots
    //! No clear: a static term's blocks are geometry-fixed, so a slot once filled stays valid for the run.
    virtual void PrepareSlots(const tbs_t<T>* bs, const tChargeDensity<T>*) const override {this->PrepareCacheSlots(bs);}

protected:
    // Unconditional calculation, does not use cache.
    virtual hmat_t<T> MakeMatrix(const tobs_t<T>*,const Spin&) const=0;
};

template <class T> class tDynamic_HT_Imp
    : public virtual tDynamic_HT<T>
    , protected tHT_Common<T>
{
public:
    tDynamic_HT_Imp() : itsCacheVersion(size_t(-1)), itsFitVersion(size_t(-1)) {}
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>* bs,const Spin& s,DensityID id) const
    {
        assert(bs);
        // The cache keys on Irrep, not the density.  If the density has changed (its logical-clock serial
        // differs from the one the cache was built for), the cached matrices are stale -- drop them and
        // recompute for the new cd.  This makes GetMatrix self-correcting for ANY caller (the serial is
        // intrinsic to the density, stamped at construction/mutation), not just densities the SCFIterator
        // hands us -- which is why a plain density-pointer key, or an SCFIterator-only stamp, both failed.
        // This is now the SOLE owner of cache-freshness; newCD() only signals the concrete term to refit.
        // THE STAMP IS THE CONTRACT (doc/KAxisPlan.md P2): PrepareSlots(bs, cd) -- the refresh -- recorded this density;
        // GetMatrix never resets the cache itself (it used to, on first use, which the concurrent block loop cannot have).
        if (id.IsNone() || !(id==itsCacheID))
            throw std::logic_error("Hamiltonian term GetMatrix for a density it was not refreshed for -- call "
                                   "tHamiltonian::RefreshForDensity(bs, cd) first (doc/KAxisPlan.md P2)");
        Irrep qns(bs->GetIrrep(this->CacheSpin(s)));
        // FILL-IF-EMPTY (R1.0h): the refresh pre-created this node, so operator[] finds it and mutates nothing; a 0x0
        // matrix means "slot exists, not filled for this density".  A serial caller that refreshed with bs = null (energy
        // evaluation) still creates the node here -- an insertion, so never from inside a concurrent loop.
        hmat_t<T>& slot=this->itsCache[qns];
        if (slot.rows()==0) slot=MakeMatrix(bs,s,id);
        return slot;   //Cache hit (same density serial, already computed for this Irrep)
    }

    //! \copydoc HT_SlotOwner::PrepareSlots
    //! Drops the previous iteration's blocks first (the same density-serial test \c GetMatrix makes, and
    //! it must stay in BOTH: this phase is optional, that one is the correctness mechanism), then creates
    //! this iteration's empty slots.
    virtual void PrepareSlots(const tbs_t<T>* bs, const tChargeDensity<T>* cd) const override
    {
        if (cd && cd->Version()!=itsCacheVersion)   // a NEW density: last one's values are stale -- drop them, stamp this one
        {
            this->itsCache.clear();
            itsCacheVersion=cd->Version();
            itsCacheID=DensityID(cd);
        }
        this->PrepareCacheSlots(bs);
    }

protected:
    // Unconditional calculation, does not use cache: a READ of what the term's refresh built for the density \a id names
    // (V1.41 stage 2), never a compute from rho.
    virtual hmat_t<T> MakeMatrix(const tobs_t<T>*,const Spin&,DensityID id) const=0;
    // Refit trigger only (the OTHER half of the former itsCD double-duty): true exactly once per new
    // density serial, so a concrete term (FittedVee/FittedVxc) refits its fitted potential exactly once.
    // Cache-freshness is GetMatrix's job (itsCacheVersion) -- newCD no longer touches the cache.
    bool newCD(const tChargeDensity<T>* cd) const
    {
        assert(cd);
        if (cd->Version()==itsFitVersion) return false;
        itsFitVersion=cd->Version();
        itsFitID=DensityID(cd);
        return true;
    }

    //! A term whose PARAMETERS changed under the SAME density (DFT+U's U rewritten by the ACBN0 outer
    //! loop) must say so: the density serial cannot know, so the next GetMatrix would serve the old matrix.
    //! The block path's check (V1.41 stage 2): the term's refresh fitted \a cd (\c newCD stamped it).  Throws otherwise.
    void RequireFitFor(DensityID id) const
    {
        if (id.IsNone() || !(id==itsFitID))
            throw std::logic_error("Hamiltonian term: block asked for a density its refresh did not fit (RefreshForDensity first)");
    }

    void InvalidateCache() const {itsCacheVersion=size_t(-1); itsFitVersion=size_t(-1); itsCacheID=itsFitID=NoDensity();}
    static DensityID NoDensity() {return DensityID(static_cast<const tChargeDensity<T>*>(nullptr));}

    mutable size_t itsCacheVersion;   //!< density serial the Irrep cache was built for
    mutable size_t itsFitVersion;     //!< density serial the concrete term last refit for
    mutable DensityID itsCacheID=NoDensity();   //!< the identity of the density the Irrep cache was built for
    mutable DensityID itsFitID  =NoDensity();   //!< ...and of the one the concrete term last refit for
};

//! The REAL-BLOCK caching mixin of a complex run's static term (V1.35): the term's ordinary \c tStatic_CC<double> face for a real
//! TRIM block.  Its slot hook is the run's (HT_SlotOwner<dcmplx>).
class Static_HT_RealBlock_Imp : public virtual ChargeDensity::tStatic_CC<double>
                              , public virtual HT_SlotOwner<dcmplx>   //!< a SECOND Irrep-keyed cache -> the SAME hook (R1.0h)
                              , public virtual HT_SpinDependence
{
public:
    virtual const hmat_t<double>& GetMatrix(const tobs_t<double>* bs,const Spin& s) const
    {
        assert(bs);
        Irrep qns(bs->GetIrrep(this->CacheSpin(s)));
        hmat_t<double>& slot=itsRealCache[qns];        // FILL-IF-EMPTY (R1.0h) -- see tDynamic_HT_Imp
        if (slot.rows()==0) slot=MakeMatrix(bs,s);
        return slot;
    }
    //! \copydoc HT_SlotOwner::PrepareSlots
    //! MY (real) cache's slots.  A term that also has a scalar cache overrides this once and calls BOTH
    //! mixins -- see HT_SlotOwner for why the ambiguity until it does is the point.
    virtual void PrepareSlots(const tbs_t<dcmplx>* bs, const tChargeDensity<dcmplx>*) const override
    {
        if (!bs) return;                                // no clear: geometry-fixed, like the scalar sibling
        for (size_t i=0;i<bs->GetNumIBS();++i)
            if (const tobs_t<double>* b=bs->GetRealIBS(i))
                for (const Spin& s : {Spin::None, Spin::Up, Spin::Down})
                    itsRealCache[b->GetIrrep(this->CacheSpin(s))];
    }
protected:
    virtual hmat_t<double> MakeMatrix(const tobs_t<double>*,const Spin&) const=0;
    mutable std::map<Irrep,hmat_t<double>> itsRealCache;
};

//! \brief THE BLOCK HELPER of a T-templated static term (V1.35 wrinkle 2): a virtual cannot carry a \c requires clause, so on
//! \c T=double a term's real-block \c MakeMatrix(const tobs_t<double>*) would collide with its native one.  The term instead writes
//! ONE \c template<class U> MakeBlock body, and this mixin (present on the \c dcmplx instantiation only) routes the real block to it.
template <class D> class Static_HT_RealBlock_Of : public Static_HT_RealBlock_Imp
{
protected:
    virtual hmat_t<double> MakeMatrix(const tobs_t<double>* bs,const Spin& s) const override
    {return static_cast<const D&>(*this).template MakeBlock<double>(bs,s);}
};

//! The REAL-BLOCK caching mixin of a complex run's dynamic term (V1.35): the term's ordinary \c tDynamic_CC<double> face, for a
//! real TRIM block, asked with the RUN's density identity.  Its slot hook is the run's (HT_SlotOwner<dcmplx>).
class Dynamic_HT_RealBlock_Imp : public virtual ChargeDensity::tDynamic_CC<double>
                               , public virtual HT_SlotOwner<dcmplx>   //!< a SECOND Irrep-keyed cache -> the SAME hook (R1.0h)
                               , public virtual HT_SpinDependence
{
public:
    virtual const hmat_t<double>& GetMatrix(const tobs_t<double>* bs,const Spin& s,DensityID id) const
    {
        assert(bs);
        if (!(id==itsRealCacheID))
            throw std::logic_error("Hamiltonian term (real block) GetMatrix for a density it was not refreshed for -- call "
                                   "tHamiltonian::RefreshForDensity(bs, cd) first (doc/KAxisPlan.md P2)");
        Irrep qns(bs->GetIrrep(this->CacheSpin(s)));
        hmat_t<double>& slot=itsRealCache[qns];        // FILL-IF-EMPTY (R1.0h) -- see tDynamic_HT_Imp
        if (slot.rows()==0) slot=MakeMatrix(bs,s,id);
        return slot;
    }
    virtual void PrepareSlots(const tbs_t<dcmplx>* bs, const tChargeDensity<dcmplx>* cd) const override
    {
        if (cd && !(DensityID(cd)==itsRealCacheID))   // a NEW density: drop the old values, stamp this one
        {
            itsRealCache.clear();
            itsRealCacheID=DensityID(cd);
        }
        if (!bs) return;
        for (size_t i=0;i<bs->GetNumIBS();++i)
            if (const tobs_t<double>* b=bs->GetRealIBS(i))
                for (const Spin& s : {Spin::None, Spin::Up, Spin::Down})
                    itsRealCache[b->GetIrrep(this->CacheSpin(s))];
    }
    //! E = D.V for the real block too (the periodic xc terms never take the DM_Contract route, so no eps override is needed).
    virtual const hmat_t<double>& GetEMatrix(const tobs_t<double>* bs,const Spin& s,DensityID id) const
    { return GetMatrix(bs,s,id); }
protected:
    virtual hmat_t<double> MakeMatrix(const tobs_t<double>*,const Spin&,DensityID) const=0;
    mutable std::map<Irrep,hmat_t<double>> itsRealCache;
    void InvalidateRealCache() const {itsRealCacheID=DensityID(static_cast<const tChargeDensity<dcmplx>*>(nullptr));}
    mutable DensityID itsRealCacheID=DensityID(static_cast<const tChargeDensity<dcmplx>*>(nullptr));
};

struct NoRealBlock {};
template <class T, class D> using StaticRealBlockBase =
    std::conditional_t<std::is_same_v<T,dcmplx>, Static_HT_RealBlock_Of<D>, NoRealBlock>;

// tFittablePotential is GONE (R2.6).  It existed for exactly one derivation, LDAVxc, whose only real job
// was to answer the fitter's "what field am I fitting?" question -- and which had to fake MakeMatrix and
// GetEnergy (both exit(-1)) to pay for being a tDynamic_HT at all.  A fittable field is now expressed as
// what it IS, a Fitting::ProjectedScalar_R adapter holding (functional, density), built where it is used
// (Imp/FittedVxc.C's VxcDensity / EpsXcDensity pair).  No term interface is involved.

// r* = <double>, c* = <dcmplx> (mirrors rsmat_t/chmat_t).
using rHT_Common              = tHT_Common<double>;              using cHT_Common              = tHT_Common<dcmplx>;
using rStatic_HT_Imp          = tStatic_HT_Imp<double>;          using cStatic_HT_Imp          = tStatic_HT_Imp<dcmplx>;
using rDynamic_HT_Imp         = tDynamic_HT_Imp<double>;         using cDynamic_HT_Imp         = tDynamic_HT_Imp<dcmplx>;

} //namespace
