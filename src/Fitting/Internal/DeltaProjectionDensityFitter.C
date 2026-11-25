// File: Fitting/Internal/DeltaProjectionDensityFitter.C  The δ-projection Hartree: a density fit on a δ fit basis whose Coulomb side is the
// analytic point-charge integrals (ij|δ_g), contracted -- no Poisson operator (lane S 2026-11-25; design note lane-s-design-ij-delta.md).
//
// The sibling of DeltaDensityFitter on the SAME δ fit basis and the same overlap projection; it differs in where the potential comes from:
//     DeltaDensityFitter              v = K q            K a CoulombOperator over the points (radial Poisson, multipole, direct sum)
//     DeltaProjectionDensityFitter    ṽ_g = (ρ|δ_g)      the exact potential of ρ_D at each point, from the block's PointCoulomb integrals
// So nothing here is a quadrature of a potential: the projection integrals ⟨ij|δ_g⟩ (overlap) and (ij|δ_g) (Coulomb) are analytic, and the
// only mesh error left is the δ fit's own (its charges q_g = w_g ρ(r_g), the mesh's charge error; reported by Vee_Fit as before).
// USER RULING 2026-11-25 (Q1 of the design note): "OK I agree on Q1, it is the analytically correct expression."  The derivation is at
// HartreeEnergy and Repulsion below.
//
// THE PROJECTIONS ARE TWO PER POINT.  ṽ needs D, and D stays in the density (a fitter never sees it), so the density contracts its D through
// this fitter's Forward(orb) into BOTH projections at once: [ρ_g ; ṽ_g], 2 N_pts coefficients.  The record's only consumer is this fitter's
// own DoFit (V-ONE-OF-N: its own client reads the realization), and the sum over blocks stays a plain vector sum, since both are linear in D.
// A density with no D (a field seed) has no (ρ|δ_g): ProjectField REFUSES and CanProjectField says so up front, so the Hartree term
// RequiresDensityMatrix and the SCF's matrix-free-seed bootstrap (SCFIterator, the HF path) hands it a D0 instead.
//
// The bodies are in Internal/Imp/DeltaProjectionDensityFitter.C.
module;
#include <functional>
#include <memory>
#include <string>
export module qchem.Fitting.Internal.DeltaProjectionDensityFitter;
export import qchem.Fitting.CoulombOperator;      // FirstOrderDensityFit
import qchem.Fitting.Types;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.PointCoulomb;               // (ij|δ_g): the block's face
import qchem.Mesh;                                 // the δ basis's point set (its points; the basis keeps them private)
import qchem.Blaze;
import qchem.Parallel;                             // LazyLock -- the per-block cache's guard
import qchem.SlotMap;                              // InstanceSlotMap

export namespace qchem::Fitting
{

class DeltaProjectionDensityFitter
    : public virtual FirstOrderDensityFit<double>
{
public:
    typedef std::shared_ptr<const BasisSet::rFIT_SF_ABS> fbs_t;
    //! \a fbs is the δ fit basis (orthogonal: its projection is its fit; its \c Overlap3C gives \f$\langle ij|\delta_g\rangle\f$), and \a mesh
    //! its own point set (the factory gathers it: the basis keeps its points private), whose points the blocks' \c PointCoulomb integrals are
    //! taken at.  THROWS when the basis is not orthogonal or the mesh is not its points.
    DeltaProjectionDensityFitter(fbs_t fbs, std::shared_ptr<const qcMesh::Mesh> mesh);
    ~DeltaProjectionDensityFitter();   //!< out of line: \c Block is complete only in the Imp unit

    //! \brief THE ENERGY, \f$E_H=\tfrac12(\rho|\tilde\rho)\f$ with \f$\tilde\rho=\sum_g q_g\delta_g\f$ the δ fit of the density:
    //! \f[ E_H=\tfrac12\sum_g q_g\,\tilde v_g,\qquad q_g=\sum_{ij}D_{ij}\,w_g\chi_i\chi_j(r_g)=w_g\rho(r_g),\qquad
    //!     \tilde v_g=(\rho|\delta_g)=\sum_{kl}D_{kl}\,(kl|\delta_g)=\int\frac{\rho(r)}{|r-r_g|}dr. \f]
    //! The cross term \f$(\rho|\tilde\rho)\f$ is exact (analytic \f$(kl|\delta_g)\f$); the self term \f$(\tilde\rho|\tilde\rho)\f$ that
    //! Dunlap's robust form also needs is NOT available (\f$(\delta_g|\delta_h)\f$ is singular at \f$g=h\f$), so the energy is first order
    //! in the fit error: \f$E_H-E_H^{exact}=\tfrac12(\rho|\tilde\rho-\rho)\f$, the mesh's charge error weighted by the exact potential.
    //! \a trDV is not read.
    double HartreeEnergy(const std::function<double()>& trDV) const override;
    //! \copydoc FirstOrderDensityFit<double>::DoFit -- the projections are \f$[\rho(r_g);\ \tilde v_g]\f$ (two per point, see the header).
    void DoFit(const Projections& proj) override;
    //! \brief THE FOCK MATRIX, \f$V_{ij}=\partial E_H/\partial D_{ij}\f$ exactly.  \f$E_H=\tfrac12\sum_gq_g\tilde v_g\f$ is bilinear in D
    //! through BOTH factors, so both carry a leg:
    //! \f[ V_{ij}=\tfrac12\sum_g\frac{\partial q_g}{\partial D_{ij}}\tilde v_g+\tfrac12\sum_g q_g\frac{\partial\tilde v_g}{\partial D_{ij}}
    //!          =\tfrac12\sum_g\big(w_g\tilde v_g\big)\langle ij|\delta_g\rangle+\tfrac12\sum_g q_g\,(ij|\delta_g), \f]
    //! with \f$\langle ij|\delta_g\rangle=\chi_i\chi_j(r_g)\f$ the overlap projection and \f$(ij|\delta_g)\f$ the Coulomb one: two
    //! contractions of fit coefficients with analytic projection integrals.  (By symmetry the two legs have the same trace against D:
    //! \f$\mathrm{Tr}(DV)=2E_H\f$.)
    //!
    //! WHY THE q-LEG SURVIVES HERE AND NOT IN DUNLAP'S FIT.  Dunlap's \f$E=(\rho|\tilde\rho)-\tfrac12(\tilde\rho|\tilde\rho)\f$ chooses its
    //! coefficients to make E STATIONARY in them (\f$(\rho|f_k)=\sum_l c_l(f_k|f_l)\f$), so \f$\partial E/\partial c=0\f$ and the chain
    //! through \f$dc/dD\f$ vanishes: its Fock is \f$\sum_kc_k(ij|f_k)\f$ alone.  Here the coefficients \f$q_g\f$ are fixed by the OVERLAP
    //! projection, not by stationarity of \f$E_H\f$ (which would need the singular \f$(\delta_g|\delta_h)\f$), so \f$\partial E_H/\partial q
    //! =\tfrac12\tilde v\neq0\f$ and the leg through \f$dq/dD\f$ stays.  Dropping it (\f$V=\sum_gq_g(ij|\delta_g)\f$ alone) is the gradient
    //! of no energy: an SCF on it would not stop where its reported \f$E_H\f$ is stationary (measured on TeF\f$_6\f$ PP 80/29: the two legs
    //! differ by 4.5e-4 Ha, lane-s-ij-delta-measurement.md §2).
    hmat_t<double> Repulsion(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const override;
    //! The fixed-D force is not built yet (design note §4: legs 1-2 of the Poisson route at ½ṽ, leg 3' = the blocks' PointCoulomb
    //! pullback): both say no, so a forces run names this term in MissingGradients.
    bool CanPullbackRepulsion(const BasisSet::Orbital_DFT_IBS<double,double>&) const override {return false;}
    Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsionPullback(const BasisSet::Orbital_DFT_IBS<double,double>&, const hmat_t<double>&,
                                                             const Structure*) const override;
    //! The mesh's charge \f$\sum_g w_g\rho(r_g)\f$ (no charge constraint, as \c DeltaDensityFitter).
    double Integral() const override;
    // PairProjector<double>
    //! The density's forward for one block: \f$[\rho(r_g);\ \tilde v_g]\f$ off the block's overlap and Coulomb projections.
    const qcMesh::MatrixForward<double>& Forward(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const override;
    size_t NumCoefficients() const override {return 2*itsW->size();}
    //! REFUSES: a field has no analytic \f$(\rho|\delta_g)\f$ (the Gaussian Coulomb fit refuses the same way).
    Projections ProjectField(const ScalarFunction<double>&) const override;
    bool CanProjectField() const override {return false;}   //!< a field has no D: the SCF's seed bootstrap must give one

private:
    struct Block;   //!< one orbital block's two projections and their joint forward
    const Block& BlockFor(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const;

    mutable InstanceSlotMap<std::unique_ptr<Block>> itsBlocks;   //!< per block, by the block's instance
    qchem::LazyLock itsCacheLock;
    fbs_t  itsFitBasis;                             //!< co-owned: the overlap integrators borrow its cached tensors
    std::shared_ptr<const qcMesh::Mesh> itsMesh;   //!< the δ basis's points (the Coulomb integrals are taken there)
    std::shared_ptr<const rvec_t> itsW;   //!< the mesh weights \f$w_g\f$ (the δ basis's \f$\langle\delta_g|\delta_g\rangle\f$)
    rvec_t itsCharge;   //!< \f$\langle\delta_g|1\rangle\f$ (the same numbers on a δ basis)
    rvec_t itsRho;      //!< \f$\rho(r_g)\f$ of the last fit
    rvec_t itsQ;        //!< \f$q_g=w_g\rho(r_g)\f$ of the last fit
    rvec_t itsV;        //!< \f$\tilde v_g\f$ of the last fit; empty until the first DoFit
};

} //namespace
