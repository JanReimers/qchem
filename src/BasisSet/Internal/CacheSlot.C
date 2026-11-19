// File: BasisSet/Internal/CacheSlot.C  A block's way into its PARENT's cache, and the cached overlap that reads through it.
//
// D-BASISCACHE-FANIN (the user's ruling 2026-11-18, option (a): "IrrepBasisSetImp has historically been (almost) everything except
// integrals").  The public faces used to declare Cache() / CacheIndex() on IrrepBasisSet_IDs, and that one protected pair forced
// `import qchem.BasisSet.BasisCache` (with ERI4, Cache4, PrimitiveStore and Cache2 behind it) into a face 381 TUs import.  The faces now
// declare only the cached Gets (Overlap, Kinetic, Overlap3C, Repulsion3C, RestMass), pure virtual; the cache access lives here, Internal,
// and the Cached* mixins implement the Gets through it.  A family lists the mixin beside the family mixin that supplies the matching Make*.
module;
#include <cstddef>
export module qchem.BasisSet.Internal.CacheSlot;
export import qchem.BasisSet.BasisCache;
export import qchem.BasisSet.IrrepBasisSet;

export namespace qchem::BasisSet
{

//! \brief A block's slot in its parent's cache: the cache itself and this block's index in it.  Data free; IrrepBasisSetImp is the one
//! implementer (it holds the cache handed to every block at construction, doc/IntegralCacheOwnershipPlan.md S3a).
//!
//! Every member is PROTECTED: only a block's own cached accessors read through it.  Inherit it `protected virtual`, listed BEFORE the
//! public faces (CLAUDE.md, the cross-cast trap: a non-public base listed after a public face can make a cross-cast return null).
class CacheSlot
{
protected:
    virtual ~CacheSlot() = default;
    virtual const BasisCache& Cache     () const=0;
    virtual size_t            CacheIndex() const=0;
    //! Is \a other a block of the SAME parent (one cache)?  CLIENT: an ERI4 pair, which a cache stores only between two blocks of one
    //! parent and refuses otherwise.
    bool SharesCacheWith(const CacheSlot& other) const;
};

//! The cached overlap: MakeOverlap() built once, in the parent's cache under this block's index.  Hermitian-correct for both lineages:
//! hmat_t<double> is symmetric (rsmat_t), hmat_t<dcmplx> is chmat_t.  No string key, so no collision class (the NaF 113x113 overlap once
//! served to a later CsI basis at the same address cannot recur: the cache dies with its basis).
template <class T> class CachedOverlap
    : protected virtual CacheSlot
    , public    virtual Integrals_Overlap<T>
{
public:
    const hmat_t<T>& Overlap() const override;
};

} //namespace
