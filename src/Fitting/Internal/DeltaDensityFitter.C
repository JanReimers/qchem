// File: Fitting/Internal/DeltaDensityFitter.C  A density fit on a δ fit basis: the projection is the fit, the operator finds v.
//
// The density-side sibling of DeltaScalarFitter, on the same δ fit basis and the same forward/adjoint pair:
//     forward    rho(r_g)                          -- the density contracts its own D, per block, as it does for a Gaussian fit
//     potential  v(r_g) = op(w_g rho(r_g))         -- the BASIS's Coulomb integrals (PointCoulombIntegrals, V-ONE-OF-N 3b-ii)
//     adjoint    V_ij = sum_g c_g (ij|δ_g)         -- c_g = w_g v(r_g); the δ basis's exact projection integrals, contracted
//
// NAMING (user ruling 2026-11-24): this is the Hartree FIT with a δ fit basis ("δ-fit Hartree"), never "grid Hartree".  A fit basis is
// named as a fit basis whatever its functions are; "mesh" and "grid" are kept for quadrature.  The XC mesh only SUPPLIES the δ points
// (and their weights, carried inside c_g); (ij|δ_g) are analytic projection integrals that get summed, not a quadrature.
//
// ENERGY IS FIRST ORDER (user ruling 2026-11-06).  Tr(DV) = sum_g w_g rho_g v_g = <rho~|rho~> here, so the energy is the plain
// 1/2 <rho~|rho~>; Dunlap's 2 EeeFit - EeeFitFit is second-order only where the cross term is exact, which needs a real rho~ with a
// solved potential.  The operator answers (HartreeEnergy), so a multipole fit keeps its two legs and Vee_DFT never asks which it got.
//
// THE FIXED-D FORCE (lane S 2026-11-23; design note lane-s-DeltaHartree-force-design.md).  E_H = ½ qᵀ K q with q_g = w_g ρ_D(r_g):
//     dE_H/dR = Σ_g v̄_g dq_g + ½ qᵀ dK q,   v̄ = ½(K + Kᵀ) q,   dq_g = ρ_g dw_g + w_g dρ_D(r_g)   (the points ride their owners)
// in three legs: (1) the weights, Σ_g v_g ρ_g dw_g, the mesh's weight response (lane B's piece A); (2) ρ at the points, Σ_g w_g v_g dρ_g,
// per block, the block's point-set pullback (lane B's piece B) -- MakeRepulsionPullback; (3) the operator, ½ qᵀ dK q, its cross-site
// reads moving (CoulombOperator::MakeSelfRepulsionPullback).  (1) + (3) are the fit side, MakeHartreePullback, taken once.
// Legs 1-2 read v̄, and so does the Fock: the δ operator answers the symmetrized ½(K + Kᵀ) itself (user ruling 2026-11-23; why, and the
// measurement that forced it, at RadialPoissonCoulomb), so v = Kq IS v̄ here and H = ∂E_H/∂D exactly.
//
// NO CHARGE CONSTRAINT (user ruling 2026-11-06).  Dunlap's shift is a Coulomb-metric object, and a uniform rescale was declined
// for now: the mesh's charge error is reported instead (Vee_DFT, as a diagnostic).  So Integral() is the mesh's charge, not N.
//
// The bodies are in Internal/Imp/DeltaDensityFitter.C.
module;
#include <functional>
#include <memory>
#include <string>
export module qchem.Fitting.Internal.DeltaDensityFitter;
export import qchem.Fitting.CoulombOperator;
import qchem.BasisSet.CoulombDeltaFit_IBS;   // PointCoulombIntegrals (the basis answers its Coulomb integrals)
import qchem.Fitting.Types;           // robs_t
import qchem.Fitting.FitOperations;   // OrthogonalFit: a matrix-free seed samples itself
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.Mesh;                        // the δ basis's point set: owners, weights, the weight response
import qchem.Structure;                   // the force's atoms
import qchem.Blaze;
import qchem.Parallel;                // LazyLock -- the per-block cache's guard
import qchem.SlotMap;                 // InstanceSlotMap

export namespace qchem::Fitting
{

class DeltaDensityFitter
    : public virtual FirstOrderDensityFit<double>
{
public:
    typedef std::shared_ptr<const BasisSet::rFIT_SF_ABS> fbs_t;
    //! The constructor DEMANDS the basis's Coulomb-integral face \a coulomb (the same object as \a fbs, gathered by the factory): the
    //! physics belongs to the basis (ruling (f)), and this fitter only asks it.
    //! And it demands the δ basis's own point set \a mesh (its quadrature's mesh, gathered by the factory too): the force reads its owners
    //! and its weight response, which the basis keeps private.
    DeltaDensityFitter(fbs_t fbs, const BasisSet::PointCoulombIntegrals& coulomb, std::shared_ptr<const qcMesh::Mesh> mesh);

    double HartreeEnergy(const std::function<double()>& trDV) const override;
    //! \copydoc FirstOrderDensityFit<double>::DoFit -- on a δ basis the projections ARE the fit (values at the points).
    void DoFit(const Projections& rho) override;
    //! \copydoc FirstOrderDensityFit<double>::Repulsion
    hmat_t<double> Repulsion(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const override;
    //! Leg 2 of the force for one block: \f$\sum_g w_gv_g\,\partial\rho_D(r_g)\f$, the block's functions moving and the points riding their
    //! owners, the weights held (\c BasisSet::PointSetPullback with \f$c_g=w_gv_g\f$).
    bool CanPullbackRepulsion(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const override;
    Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsionPullback(const BasisSet::Orbital_DFT_IBS<double,double>& orb, const hmat_t<double>& D,
                                                             const Structure* st) const override;
    //! Legs 1 and 3, the fit side: \f$\sum_gv_g\rho_g\,\partial w_g+\tfrac12q^T\partial Kq\f$.
    bool CanPullbackHartree() const override {return ForceGap().empty();}
    Outcome<rvec3vec_t,AD::NoGradient> MakeHartreePullback(const Structure* st) const override;
    double Integral() const override;
    // PairProjector<double>
    const qcMesh::MatrixForward<double>& Forward(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const override {return Integrator(orb);}
    size_t NumCoefficients() const override {return itsW->size();}
    //! A field's projection on the δ basis: its values at the points (\c OrthogonalFit), the axis \c Forward uses.
    Projections ProjectField(const ScalarFunction<double>& f) const override {return Projections(OrthogonalFit(*itsFitBasis,f));}

private:
    //! Why this fit gives no fixed-D force, or "" when it does.
    std::string ForceGap() const;
    //! \c ForceGap, plus a fit to differentiate and the structure's atoms being the mesh's sites.
    std::string AtomsGap(const Structure* st) const;
    rvec_t Projection() const {return (*itsW)*itsC;}   // <f_a|rho>.  Blaze: element-wise for two column vectors
    //! The forward+adjoint pair for one orbital block, off the basis's cached tensor.  Locked for the same reason
    //! DeltaScalarFitter::Integrator is: a parallel Fock pass's blocks share this fitter.
    const ScreenedMatrixIntegrator<double>& Integrator(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const;

    mutable InstanceSlotMap<ScreenedMatrixIntegrator<double>> itsInt;   //!< per block, by the block's instance
    qchem::LazyLock itsCacheLock;
    fbs_t  itsFitBasis;                             //!< co-owned: the integrators' closures capture it
    const BasisSet::PointCoulombIntegrals* itsCoulomb;   //!< the basis's own Coulomb integrals (a face of itsFitBasis, co-owned there)
    std::shared_ptr<const qcMesh::Mesh> itsMesh;   //!< the δ basis's points (the force's owners and weight response)
    std::shared_ptr<const rvec_t> itsW;   //!< \f$\langle f_a|f_a\rangle\f$: the mesh weights on a \f$\delta\f$ basis (shared with the views)
    rvec_t itsCharge;   //!< \f$\langle f_a|1\rangle\f$: the same numbers there, not on a multipole basis
    rvec_t itsC;   //!< the fit: on a \f$\delta\f$ basis, the density's values at the points
    rvec_t itsV;   //!< the potential's coefficients; empty until the first DoFit
};

} //namespace
