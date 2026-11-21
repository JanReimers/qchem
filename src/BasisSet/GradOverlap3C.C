// File: BasisSet/GradOverlap3C.C  The GRADIENT 3-centre overlap of a point (δ) fit basis: what a GGA's density gradient and its
// gradient term of H are contracted through (GGA stage 2, the molecular point route; GGA design note §3b).
//
// A GGA's energy reads \f$\sigma_{\sigma\tau}=\nabla\rho_\sigma\cdot\nabla\rho_\tau\f$ at the XC points, and its KS matrix gains
// \f$\sum_g w_g\,\mathbf{w}_g\cdot\nabla(\chi_i\chi_j)(r_g)\f$.  Both are ONE three-index object, \f$\nabla(\chi_i\chi_j)(r_g)\f$, contracted the
// two ways a \c PairProjection is: FORWARD with a density matrix (\f$\nabla\rho(r_g)=\sum_{ij}D_{ij}\nabla(\chi_i\chi_j)(r_g)\f$) and ADJOINT
// with a field (\f$\sum_g w_g\,c_g\cdot\nabla(\chi_i\chi_j)(r_g)\f$).  So it is a \c PairProjection, with the three Cartesian components laid
// out as ONE coefficient axis of length 3N -- x for every point, then y, then z -- and the density's ONE operation (\c Project) gives
// \f$\nabla\rho\f$ through it unchanged: no density-hub change, no second contraction protocol.
//
// THE GRADIENT IS A TYPE (user, 2026-11-20): a fit basis that can give it derives this face, and one that cannot does not have it.  Only the
// REAL (molecular) δ basis has it in stage 2; the Bloch one gains it in stage 3, on the shared-Φ images.  Not a hub (V-INTERFACE-TIERS).
module;
export module qchem.BasisSet.GradOverlap3C;
export import qchem.BasisSet.Orbital_DFT_IBS;   // Orbital_DFT_IBS<U,TFit>, PairProjection<U>

export namespace qchem::BasisSet
{

//! \brief \f$\nabla(\chi_i\chi_j)(r_g)\f$ over a point fit basis's functions, per orbital block, as a \c PairProjection whose coefficient
//! axis is 3N long (component-major: \f$c_{aN+g}\f$ is Cartesian component a at point g).
//!   - \c applyRaw(D): \f$[\nabla\rho]_{aN+g}=\sum_{ij}D_{ij}\,\partial_a(\chi_i\chi_j)(r_g)\f$, for a Hermitian D;
//!   - \c applyRawAdjoint(c): \f$H_{ij}=\sum_g w_g\sum_a c_{aN+g}\,\partial_a(\chi_i\chi_j)(r_g)\f$, the basis's own weights inside, exactly as
//!     \c Overlap3C's adjoint carries them -- so the pair is the exact transpose and \f$H=\partial E/\partial D\f$ holds for the gradient part.
//! Cached per block for the run (geometry-fixed: the lazy-caching rule allows it), keyed as the \f$\Phi\f$ tables are.
template <class U, class TFit> class Integrals_GradOverlap3C
{
public:
    virtual ~Integrals_GradOverlap3C() = default;
    virtual const PairProjection<U>& GradOverlap3C(const Orbital_DFT_IBS<U,TFit>& orb) const=0;
};

} //namespace
