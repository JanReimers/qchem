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
#include <vector>
#include <cassert>
#include <cstddef>
#include <type_traits>   // std::conditional_t/is_same_v (StaticRealBlockBase -- Step 3c)
export module qchem.Hamiltonian.Internal.Term;
export import qchem.Hamiltonian;
import qchem.Hamiltonian.Types;

import qchem.Symmetry.Irrep;
import qchem.SlotMap;      // every block cache below (doc/IntegralCacheOwnershipPlan.md S1)

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
    //! The spins the run will ask this term's blocks for, which \c BuildBlocks builds (each folded through \c CacheSpin, so a
    //! spin-blind term builds one block per irrep).  Default all three; a term built for an imposed spin subgroup answers
    //! \c SpinIrreps(g), so it builds exactly the channels the Fock loop reads and no more.
    virtual std::vector<Spin> BlockSpins() const {return {Spin::None, Spin::Up, Spin::Down};}
};

//! \brief A term's block cache: one \c hmat_t<T> per \c Irrep, held in a \c SlotMap (doc/IntegralCacheOwnershipPlan.md S1).  \a D is
//! the discipline: a STATIC term's blocks are a pure function of the geometry (\c FillOnce: built on first read, never cleared); a
//! DYNAMIC term's depend on the density (\c Refresh: built only in \c BuildBlocks, the refresh, and read purely after it -- V1.41).
//! The map's lock and per-slot once-fill replace the old pre-created-slot dance (R1.0h): the block loop may ask concurrently, and two
//! spins a term folds onto one key (\c CacheSpin) build that block once.
template <class T, SlotDiscipline D> class tHT_Common
    : public virtual HT_SpinDependence
{
protected:
    mutable SlotMap<Irrep,hmat_t<T>,D> itsCache;
};


template <class T> class tStatic_HT_Imp
    : public virtual tStatic_HT<T>
    , protected tHT_Common<T,SlotDiscipline::FillOnce>
{
public:
    //! Built on the first read of each block (once, whichever thread asks first), then served forever: a static term's cache is never
    //! cleared.
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>* bs,const Spin& s) const
    {
        assert(bs);
        return this->itsCache.GetOrFill(bs->GetIrrep(this->CacheSpin(s)), [&]{return MakeMatrix(bs,s);});
    }

    //! \copydoc HT_BlockOwner::BuildBlocks
    //! A static term's blocks NEVER change (geometry-fixed), the one case where lazy caching is allowed (CLAUDE.md): nothing to do
    //! here, \c GetMatrix fills each block on its first read -- the same value whenever it is filled.
    virtual void BuildBlocks(const tbs_t<T>*, const tChargeDensity<T>*) const override {}

protected:
    // Unconditional calculation, does not use cache.
    virtual hmat_t<T> MakeMatrix(const tobs_t<T>*,const Spin&) const=0;
};

template <class T> class tDynamic_HT_Imp
    : public virtual tDynamic_HT<T>
    , protected tHT_Common<T,SlotDiscipline::Refresh>
{
public:
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>* bs,const Spin& s,DensityID id) const
    {
        assert(bs);
        // A PURE READ (V1.41 stage 3, user rulings 2026-11-07): BuildBlocks -- the refresh -- built every block the run asks for,
        // for the density it stamped.  No fill, no insertion, so the block loop writes nothing and can run concurrently.
        if (id.IsNone() || !(id==itsCacheID))
            throw std::logic_error("Hamiltonian term GetMatrix for a density it was not refreshed for -- call "
                                   "tHamiltonian::RefreshForDensity(bs, cd) first (doc/KAxisPlan.md P2)");
        const hmat_t<T>* block=this->itsCache.Read(bs->GetIrrep(this->CacheSpin(s)));
        if (!block)
            throw std::logic_error("Hamiltonian term GetMatrix: a block the refresh did not build (its BlockSpins, or the "
                                   "basis BuildBlocks was handed, does not cover it)");
        return *block;
    }

    //! \copydoc HT_BlockOwner::BuildBlocks
    //! ALWAYS rebuilds (no "same density, skip"): no Fock refresh runs twice at one density, and a parameter write under the
    //! same density (DFT+U's U) is then visible as the next refresh, with nothing to invalidate.  The term's own
    //! RefreshForDensity(cd) has run first, so MakeMatrix reads state built for \a cd.  \a bs null (an energy evaluation):
    //! nothing to build, and the memo keeps the Fock density's blocks.
    virtual void BuildBlocks(const tbs_t<T>* bs, const tChargeDensity<T>* cd) const override
    {
        if (!cd || !bs) return;
        const auto refresh=this->itsCache.BeginRefresh();   // clears the previous density's blocks
        itsCacheID=DensityID(cd);
        const DensityID id(cd);
        for (size_t i=0;i<bs->GetNumIBS();++i)
        {
            if (bs->GetRealIBS(i)) continue;   // the real cache's (Dynamic_HT_RealBlock_Imp)
            if (const tobs_t<T>* b=(*bs)[i])
                for (const Spin& s : this->BlockSpins())   // a spin-blind term folds the spins onto one block: built once
                    refresh.Fill(b->GetIrrep(this->CacheSpin(s)), [&]{return MakeMatrix(b,s,id);});
        }
    }

protected:
    // Unconditional calculation, does not use cache: a READ of what the term's refresh built for the density \a id names
    // (V1.41 stage 2), never a compute from rho.
    virtual hmat_t<T> MakeMatrix(const tobs_t<T>*,const Spin&,DensityID id) const=0;
    // Refit trigger only (the OTHER half of the former itsCD double-duty): true exactly once per new
    // density serial, so a concrete term (FittedVee/FittedVxc) refits its fitted potential exactly once.
    // Cache-freshness is GetMatrix's job (itsCacheID) -- newCD no longer touches the cache.
    bool newCD(const tChargeDensity<T>* cd) const
    {
        assert(cd);
        if (DensityID(cd)==itsFitID) return false;
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

    static DensityID NoDensity() {return DensityID(static_cast<const tChargeDensity<T>*>(nullptr));}

    mutable DensityID itsCacheID=NoDensity();   //!< the identity of the density the Irrep cache was built for
    mutable DensityID itsFitID  =NoDensity();   //!< ...and of the one the concrete term last refit for
};

//! The REAL-BLOCK caching mixin of a complex run's static term (V1.35): the term's ordinary \c tStatic_CC<double> face for a real
//! TRIM block.  Its BuildBlocks hook is the run's (HT_BlockOwner<dcmplx>).
class Static_HT_RealBlock_Imp : public virtual ChargeDensity::tStatic_CC<double>
                              , public virtual HT_BlockOwner<dcmplx>   //!< a SECOND Irrep-keyed block cache -> the SAME BuildBlocks (R1.0h)
                              , public virtual HT_SpinDependence
{
public:
    virtual const hmat_t<double>& GetMatrix(const tobs_t<double>* bs,const Spin& s) const
    {
        assert(bs);
        return itsRealCache.GetOrFill(bs->GetIrrep(this->CacheSpin(s)), [&]{return MakeMatrix(bs,s);});   // see tStatic_HT_Imp
    }
    //! \copydoc HT_BlockOwner::BuildBlocks
    //! MY (real) cache is filled on first read (static: never changes), so nothing to do.  A term that also has a scalar cache
    //! overrides this once and calls BOTH mixins -- see HT_BlockOwner for why the ambiguity until it does is the point.
    virtual void BuildBlocks(const tbs_t<dcmplx>*, const tChargeDensity<dcmplx>*) const override {}
protected:
    virtual hmat_t<double> MakeMatrix(const tobs_t<double>*,const Spin&) const=0;
    mutable SlotMap<Irrep,hmat_t<double>,SlotDiscipline::FillOnce> itsRealCache;
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
//! real TRIM block, asked with the RUN's density identity.  Its BuildBlocks hook is the run's (HT_BlockOwner<dcmplx>).
class Dynamic_HT_RealBlock_Imp : public virtual ChargeDensity::tDynamic_CC<double>
                               , public virtual HT_BlockOwner<dcmplx>   //!< a SECOND Irrep-keyed block cache -> the SAME BuildBlocks (R1.0h)
                               , public virtual HT_SpinDependence
{
public:
    virtual const hmat_t<double>& GetMatrix(const tobs_t<double>* bs,const Spin& s,DensityID id) const
    {
        assert(bs);
        if (!(id==itsRealCacheID))
            throw std::logic_error("Hamiltonian term (real block) GetMatrix for a density it was not refreshed for -- call "
                                   "tHamiltonian::RefreshForDensity(bs, cd) first (doc/KAxisPlan.md P2)");
        const hmat_t<double>* block=itsRealCache.Read(bs->GetIrrep(this->CacheSpin(s)));   // a PURE READ -- see tDynamic_HT_Imp
        if (!block)
            throw std::logic_error("Hamiltonian term (real block) GetMatrix: a block the refresh did not build");
        return *block;
    }
    //! \copydoc HT_BlockOwner::BuildBlocks
    //! The real TRIM blocks of \a bs, built like tDynamic_HT_Imp's (always rebuilt; \a bs null builds nothing).
    virtual void BuildBlocks(const tbs_t<dcmplx>* bs, const tChargeDensity<dcmplx>* cd) const override
    {
        if (!cd || !bs) return;
        const auto refresh=itsRealCache.BeginRefresh();
        itsRealCacheID=DensityID(cd);
        const DensityID id(cd);
        for (size_t i=0;i<bs->GetNumIBS();++i)
            if (const tobs_t<double>* b=bs->GetRealIBS(i))
                for (const Spin& s : this->BlockSpins())
                    refresh.Fill(b->GetIrrep(this->CacheSpin(s)), [&]{return MakeMatrix(b,s,id);});
    }
    //! E = D.V for the real block too (the periodic xc terms never take the DM_Contract route, so no eps override is needed).
    virtual const hmat_t<double>& GetEMatrix(const tobs_t<double>* bs,const Spin& s,DensityID id) const
    { return GetMatrix(bs,s,id); }
protected:
    virtual hmat_t<double> MakeMatrix(const tobs_t<double>*,const Spin&,DensityID) const=0;
    mutable SlotMap<Irrep,hmat_t<double>,SlotDiscipline::Refresh> itsRealCache;
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
using rStatic_HT_Imp          = tStatic_HT_Imp<double>;          using cStatic_HT_Imp          = tStatic_HT_Imp<dcmplx>;
using rDynamic_HT_Imp         = tDynamic_HT_Imp<double>;         using cDynamic_HT_Imp         = tDynamic_HT_Imp<dcmplx>;

} //namespace
