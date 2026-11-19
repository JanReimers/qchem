// File: BasisSet/Internal/CachedIntegrals.C  The cached kinetic and 3-centre accessors, read through a block's CacheSlot.
//
// D-BASISCACHE-FANIN (the user's ruling 2026-11-18, option (a)): the faces declare these Gets pure virtual and never see the cache; a
// family lists the mixin beside the family mixin that supplies the matching Make* (see Internal/CacheSlot.C).
module;
export module qchem.BasisSet.Internal.CachedIntegrals;
export import qchem.BasisSet.Internal.CacheSlot;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Orbital_DFT_IBS;

export namespace qchem::BasisSet
{

//! The cached \f$\langle p^2\rangle=\langle-\nabla^2\rangle\f$ block (NOT the kinetic energy: no 1/2; see BasisSet/Orbital_1E_IBS.C):
//! MakeKinetic() built once, in the parent's cache under this block's index.
template <class T> class CachedKinetic
    : protected virtual CacheSlot
    , public    virtual Integrals_Kinetic<T>
{
public:
    const hmat_t<T>& Kinetic() const override;
};

//! The cached 3-centre tensors \f$\langle ab|c\rangle\f$ (overlap and repulsion metric), keyed by (this block's index, the fit's INSTANCE
//! serial) (doc/IntegralCacheOwnershipPlan.md S3a, §3b): never the fit's address (reused after it dies) nor its content ID.  The element
//! type follows the FIT axis (TFit, NOT T), so a real TRIM block on the run's complex fit lands in the cache's TFit map.
template <class T, class TFit = T> class CachedFit3C
    : protected virtual CacheSlot
    , public    virtual Orbital_DFT_IBS<T,TFit>
{
public:
    const PairProjection<TFit>& Overlap3C  (const FIT_SF_ABS<TFit>& c) const override;
    const PairProjection<TFit>& Repulsion3C(const FIT_CD_ABS<TFit>& c) const override;
};

} //namespace
