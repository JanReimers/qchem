// File: Fitting/Internal/GridDensityFitter.C  A density fit on a mesh: the projection is the fit, the operator finds v.
//
// The density-side sibling of DeltaScalarFitter, on the same δ fit basis and the same forward/adjoint pair:
//     forward    rho(r_g)                          -- the density contracts its own D, per block, as it does for a Gaussian fit
//     potential  v(r_g) = op(w_g rho(r_g))         -- the BASIS's Coulomb integrals (PointCoulombIntegrals, V-ONE-OF-N 3b-ii)
//     adjoint    sum_g w_g v(r_g) chi_i chi_j      -- the same quadrature XC uses
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
module;
#include <cassert>
#include <functional>
#include <memory>
#include <mutex>
#include <ostream>
#include <stdexcept>
#include <string>
export module qchem.Fitting.Internal.GridDensityFitter;
export import qchem.Fitting.CoulombOperator;
import qchem.BasisSet.CoulombDeltaFit_IBS;   // PointCoulombIntegrals (the basis answers its Coulomb integrals)
import qchem.Fitting.Types;           // robs_t
import qchem.Fitting.FitOperations;   // OrthogonalFit: a matrix-free seed samples itself
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.PointSetPullback;   // the block face leg 2 asks (lane B, READY 4)
import qchem.Mesh;                        // the δ basis's point set: owners, weights, the weight response
import qchem.Structure;                   // the force's atoms
import qchem.Blaze;
import qchem.Parallel;                // LazyLock -- the per-block cache's guard
import qchem.SlotMap;                 // InstanceSlotMap

export namespace qchem::Fitting
{

class GridDensityFitter
    : public virtual FirstOrderDensityFit<double>
{
public:
    typedef std::shared_ptr<const BasisSet::rFIT_SF_ABS> fbs_t;
    //! The constructor DEMANDS the basis's Coulomb-integral face \a coulomb (the same object as \a fbs, gathered by the factory): the
    //! physics belongs to the basis (ruling (f)), and this fitter only asks it.
    //! And it demands the δ basis's own point set \a mesh (its quadrature's mesh, gathered by the factory too): the force reads its owners
    //! and its weight response, which the basis keeps private.
    GridDensityFitter(fbs_t fbs, const BasisSet::PointCoulombIntegrals& coulomb, std::shared_ptr<const qcMesh::Mesh> mesh)
        : itsFitBasis(std::move(fbs)), itsCoulomb(&coulomb), itsMesh(std::move(mesh))
    {
        if (!itsFitBasis) throw std::logic_error("GridDensityFitter: needs a fit basis");
        if (!itsMesh) throw std::logic_error("GridDensityFitter: needs the fit basis's mesh");
        if (!itsFitBasis->isOrtho())
            throw std::logic_error("GridDensityFitter: the fit basis '"+itsFitBasis->Name()+"' is not orthogonal, so its projection "
                                   "is not its fit -- a grid density fit needs a diagonal metric");
        itsW=std::make_shared<const rvec_t>(itsFitBasis->OverlapDiagonal());
        itsCharge=itsFitBasis->Charge();
        double total=0.0;
        for (size_t g=0; g<itsCharge.size(); g++) total+=itsCharge[g];
        itsC=rvec_t(itsW->size(),1.0/total);   // the initial guess: a uniform density of unit charge, as ConstrainedFF's guess is
    }

    double HartreeEnergy(const std::function<double()>& trDV) const override
    {
        if (itsV.size()==0) throw std::logic_error("GridDensityFitter: asked for an energy before any DoFit");
        return itsCoulomb->HartreeEnergy(itsCoulomb->UsesTrDV() ? trDV() : 0.0, Projection(), itsV);
    }
    //! \copydoc FirstOrderDensityFit<double>::DoFit -- on a δ basis the projections ARE the fit (values at the points).
    void DoFit(const Projections& rho) override
    {
        itsC=rho.Vector();
        itsV=itsCoulomb->Potential(Projection());
    }
    //! \copydoc FirstOrderDensityFit<double>::Repulsion
    hmat_t<double> Repulsion(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const override
    {
        if (itsV.size()==0) throw std::logic_error("GridDensityFitter: asked for a repulsion before any DoFit");
        return Integrator(orb).Adjoint(itsV);   // sum_g w_g v_g chi_i chi_j
    }
    //! Leg 2 of the force for one block: \f$\sum_g w_gv_g\,\partial\rho_D(r_g)\f$, the block's functions moving and the points riding their
    //! owners, the weights held (\c BasisSet::PointSetPullback with \f$c_g=w_gv_g\f$).
    bool CanPullbackRepulsion(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const override
    {
        if (!ForceGap().empty()) return false;
        const auto* face=dynamic_cast<const BasisSet::PointSetPullback<double>*>(&orb);   // abstract -> abstract, as the XC sampler's
        return face && face->CanPointSetPullback();
    }
    Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsionPullback(const BasisSet::Orbital_DFT_IBS<double,double>& orb, const hmat_t<double>& D,
                                                             const Structure* st) const override
    {
        using R=Outcome<rvec3vec_t,AD::NoGradient>;
        if (const std::string gap=AtomsGap(st); !gap.empty()) return R::Fail({gap});
        const auto* face=dynamic_cast<const BasisSet::PointSetPullback<double>*>(&orb);
        if (!face || !face->CanPointSetPullback()) return R::Fail({"the orbital block "+orb.Name()+" gives no point-set pullback"});
        const qcMesh::Mesh& m=*itsMesh;
        rvec_t c(m.size());
        for (size_t g=0; g<c.size(); g++) c[g]=(*itsW)[g]*itsV[g];
        const auto p=face->MakePointSetPullback(m.Points(), c, D, st);
        if (!p) return R::Fail(p.Error());
        rvec3vec_t G=p->gradient;
        for (size_t a=0; a<m.NSites(); a++)
            for (size_t g=m.SiteBegin(a); g<m.SiteEnd(a); g++) G[a]+=p->perPoint[g];   // the point rides its owner
        return R::Ok(std::move(G));
    }
    //! Legs 1 and 3, the fit side: \f$\sum_gv_g\rho_g\,\partial w_g+\tfrac12q^T\partial Kq\f$.
    bool CanPullbackHartree() const override {return ForceGap().empty();}
    Outcome<rvec3vec_t,AD::NoGradient> MakeHartreePullback(const Structure* st) const override
    {
        using R=Outcome<rvec3vec_t,AD::NoGradient>;
        if (const std::string gap=AtomsGap(st); !gap.empty()) return R::Fail({gap});
        const qcMesh::Mesh& m=*itsMesh;
        rvec_t f(m.size());
        for (size_t g=0; g<f.size(); g++) f[g]=itsV[g]*itsC[g];
        const auto w=m.WeightResponse()->Pullback(m, f);                       // leg 1
        const rvec3vec_t k=itsCoulomb->MakeSelfRepulsionPullback(Projection());  // leg 3, scaled as qᵀKq
        rvec3vec_t G(m.NSites());
        for (size_t a=0; a<G.size(); a++) G[a]=w.gradient[a]+0.5*k[a];
        return R::Ok(std::move(G));
    }
    double Integral() const override
    {
        double n=0.0;
        for (size_t g=0; g<itsC.size(); g++) n+=itsCharge[g]*itsC[g];
        return n;
    }
    // PairProjector<double>
    const qcMesh::MatrixForward<double>& Forward(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const override {return Integrator(orb);}
    size_t NumCoefficients() const override {return itsW->size();}
    //! A field's projection on the δ basis: its values at the points (\c OrthogonalFit), the axis \c Forward uses.
    Projections ProjectField(const ScalarFunction<double>& f) const override {return Projections(OrthogonalFit(*itsFitBasis,f));}

private:
    //! Why this fit gives no fixed-D force, or "" when it does.
    std::string ForceGap() const
    {
        if (itsCoulomb->UsesTrDV()) return "a multipole fit's two-leg Hartree energy gives no force yet";
        if (!itsCoulomb->CanPullbackSelfRepulsion()) return "the δ fit's Coulomb operator gives no pullback (the direct sum)";
        const qcMesh::Mesh& m=*itsMesh;
        if (m.size()!=itsW->size()) return "the δ fit's functions are not its mesh's points";
        if (m.NSites()==0 || m.OwnerOffsets().size()!=m.size()) return "the δ fit's mesh records no point owners (not an atom-centred build)";
        if (!m.WeightResponse()) return "the δ fit's mesh builder installed no weight response";
        return "";
    }
    //! \c ForceGap, plus a fit to differentiate and the structure's atoms being the mesh's sites.
    std::string AtomsGap(const Structure* st) const
    {
        if (std::string gap=ForceGap(); !gap.empty()) return gap;
        if (itsV.size()==0) return "a δ fit asked for a force before any DoFit";
        if (!st || st->GetNumAtoms()!=itsMesh->NSites()) return "the δ fit's mesh site blocks are not the structure's atoms";
        return "";
    }
    rvec_t Projection() const {return (*itsW)*itsC;}   // <f_a|rho>.  Blaze: element-wise for two column vectors
    //! The forward+adjoint pair for one orbital block, off the basis's cached tensor.  Locked for the same reason
    //! DeltaScalarFitter::Integrator is: a parallel Fock pass's blocks share this fitter.
    const ScreenedMatrixIntegrator<double>& Integrator(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const
    {
        std::lock_guard<std::mutex> lock(itsCacheLock.m);
        return itsInt.GetOrFill(orb.Instance(), [&]   // by the block's INSTANCE, not its symmetry (lane H 2026-11-15)
        {
            const auto& face=dynamic_cast<const BasisSet::Integrals_Overlap3C<double,double>&>(*itsFitBasis);
            return ScreenedMatrixIntegrator<double>(face.Overlap3C(orb), itsW);   // the weights are shared, not copied per block
        });
    }

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
