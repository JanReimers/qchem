// File: BasisSet/PointCoulomb.C  A block's point-charge Coulomb integrals (ij|δ_g): what the δ-projection Hartree asks of an orbital block
// (lane S 2026-11-25; design note lane-s-design-ij-delta.md, the user's ruling on its Q1 2026-11-25).
//
// The δ fit basis keeps its points private (DeltaFit_IBS, ruling 2026-08-22: it answers, it does not hand out), so the CLIENT that holds the
// points (the fitter, handed the mesh by its factory) passes them to the block, exactly as it does for PointSetPullback.  Not a tier-1 face:
// a capability a block may have, reached by an abstract -> abstract cross-cast from the block the fitter is handed.
module;
export module qchem.BasisSet.PointCoulomb;
export import qchem.Types;                    // rvec3vec_t
export import qchem.BasisSet.PairProjection;  // PairProjection<T>: the raw forward/adjoint pair

export namespace qchem::BasisSet
{
//! \brief THE POINT-CHARGE COULOMB INTEGRALS of this block's functions at a point set,
//! \f$(ij|\delta_g)=\int\frac{\chi_i(r)\chi_j(r)}{|r-r_g|}\,dr\f$: analytic projection integrals of the pair products onto the δ functions in
//! the Coulomb metric (no quadrature anywhere).  Returned as the RAW pair of a \c PairProjection, matrix-free (dense \f$n^2N_{pts}\f$ is tens
//! of GB on a production mesh):
//!   - \c applyRaw (forward): \f$D\mapsto\tilde v_g=\sum_{ij}D_{ij}(ij|\delta_g)\f$, the exact potential of \f$\rho_D\f$ at each point;
//!   - \c applyRawAdjoint: \f$q\mapsto\sum_g q_g(ij|\delta_g)\f$, the matrix of the point charges' potential.
//! D-free, as every \c PairProjection: a closure over the block's own integrals, recomputed per call (nothing is stored).
//! CLIENT: the δ-projection density fitter (\c Fitting::DeltaProjectionDensityFitter).
template <class T> class PointCoulomb
{
public:
    virtual ~PointCoulomb() = default;
    //! Can \c MakePointCoulomb answer?  A block that implements it says yes; a DECORATOR asks what it wraps (the SALC block).
    virtual bool CanPointCoulomb() const {return true;}
    //! The pair at the points \a pts (copied into the closures).
    virtual PairProjection<T> MakePointCoulomb(const rvec3vec_t& pts) const=0;
};
} // namespace qchem::BasisSet
