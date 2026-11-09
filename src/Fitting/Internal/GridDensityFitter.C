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
// NO CHARGE CONSTRAINT (user ruling 2026-11-06).  Dunlap's shift is a Coulomb-metric object, and a uniform rescale was declined
// for now: the mesh's charge error is reported instead (Vee_DFT, as a diagnostic).  So Integral() is the mesh's charge, not N.
module;
#include <cassert>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <ostream>
#include <stdexcept>
export module qchem.Fitting.Internal.GridDensityFitter;
export import qchem.Fitting.CoulombOperator;
import qchem.BasisSet.CoulombDeltaFit_IBS;   // PointCoulombIntegrals (the basis answers its Coulomb integrals)
import qchem.Fitting.Types;           // robs_t
import qchem.Fitting.FitOperations;   // OrthogonalFit: a matrix-free seed samples itself
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.Blaze;
import qchem.Parallel;                // LazyLock -- the per-block cache's guard

export namespace qchem::Fitting
{

class GridDensityFitter
    : public virtual FirstOrderDensityFit<double>
{
public:
    typedef std::shared_ptr<const BasisSet::rFIT_SF_ABS> fbs_t;
    //! The constructor DEMANDS the basis's Coulomb-integral face \a coulomb (the same object as \a fbs, gathered by the factory): the
    //! physics belongs to the basis (ruling (f)), and this fitter only asks it.
    GridDensityFitter(fbs_t fbs, const BasisSet::PointCoulombIntegrals& coulomb)
        : itsFitBasis(std::move(fbs)), itsCoulomb(&coulomb)
    {
        if (!itsFitBasis) throw std::logic_error("GridDensityFitter: needs a fit basis");
        if (!itsFitBasis->isOrtho())
            throw std::logic_error("GridDensityFitter: the fit basis '"+itsFitBasis->Name()+"' is not orthogonal, so its projection "
                                   "is not its fit -- a grid density fit needs a diagonal metric");
        itsW=itsFitBasis->OverlapDiagonal();
        itsCharge=itsFitBasis->Charge();
        double total=0.0;
        for (size_t g=0; g<itsCharge.size(); g++) total+=itsCharge[g];
        itsC=rvec_t(itsW.size(),1.0/total);   // the initial guess: a uniform density of unit charge, as ConstrainedFF's guess is
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
    double Integral() const override
    {
        double n=0.0;
        for (size_t g=0; g<itsC.size(); g++) n+=itsCharge[g]*itsC[g];
        return n;
    }
    // PairProjector<double>
    const qcMesh::MatrixForward<double>& Forward(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const override {return Integrator(orb);}
    size_t NumCoefficients() const override {return itsW.size();}
    //! A field's projection on the δ basis: its values at the points (\c OrthogonalFit), the axis \c Forward uses.
    Projections ProjectField(const ScalarFunction<double>& f) const override {return Projections(OrthogonalFit(*itsFitBasis,f));}

private:
    rvec_t Projection() const {return itsW*itsC;}   // <f_a|rho>.  Blaze: element-wise for two column vectors
    //! The forward+adjoint pair for one orbital block, off the basis's cached tensor.  Locked for the same reason
    //! DeltaScalarFitter::Integrator is: a parallel Fock pass's blocks share this fitter.
    const ScreenedMatrixIntegrator<double>& Integrator(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const
    {
        std::lock_guard<std::mutex> lock(itsCacheLock.m);
        const sym_t& id=orb.GetSymt();
        auto it=itsInt.find(id);
        if (it!=itsInt.end()) return it->second;
        const auto& face=dynamic_cast<const BasisSet::Integrals_Overlap3C<double,double>&>(*itsFitBasis);
        return itsInt.emplace(id, ScreenedMatrixIntegrator<double>(face.Overlap3C(orb), itsW)).first->second;
    }

    mutable SymMap<ScreenedMatrixIntegrator<double>> itsInt;   //!< per block, by SPATIAL symmetry
    qchem::LazyLock itsCacheLock;
    fbs_t  itsFitBasis;                             //!< co-owned: the integrators' closures capture it
    const BasisSet::PointCoulombIntegrals* itsCoulomb;   //!< the basis's own Coulomb integrals (a face of itsFitBasis, co-owned there)
    rvec_t itsW;        //!< \f$\langle f_a|f_a\rangle\f$: the mesh weights on a \f$\delta\f$ basis
    rvec_t itsCharge;   //!< \f$\langle f_a|1\rangle\f$: the same numbers there, not on a multipole basis
    rvec_t itsC;   //!< the fit: on a \f$\delta\f$ basis, the density's values at the points
    rvec_t itsV;   //!< the potential's coefficients; empty until the first DoFit
};

} //namespace
