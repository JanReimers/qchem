// File: Fitting/Internal/GridDensityFitter.C  A density fit on a mesh: the projection is the fit, the operator finds v.
//
// The density-side sibling of DeltaScalarFitter, on the same δ fit basis and the same forward/adjoint pair:
//     forward    rho(r_g)                          -- the density contracts its own D, per block, as it does for a Gaussian fit
//     potential  v(r_g) = op(w_g rho(r_g))         -- the CoulombOperator's whole job
//     adjoint    sum_g w_g v(r_g) chi_i chi_j      -- the same quadrature XC uses
//
// ENERGY IS FIRST ORDER (user ruling 2026-11-06).  Tr(DV) = sum_g w_g rho_g v_g = <rho~|rho~> here, so FittedVee's Dunlap
// combination 2 EeeFit - EeeFitFit collapses to the plain 1/2 <rho~|rho~>.  It is second-order only where the cross term is
// exact, which on the Gaussian route it is.  Nothing needs changing in FittedVee; the term reports which route it is.
//
// NO CHARGE CONSTRAINT (user ruling 2026-11-06).  Dunlap's shift is a Coulomb-metric object, and a uniform rescale was declined
// for now: the mesh's charge error is reported instead (FittedVee, as a diagnostic).  So Integral() is the mesh's charge, not N.
module;
#include <cassert>
#include <map>
#include <memory>
#include <mutex>
#include <ostream>
#include <stdexcept>
export module qchem.Fitting.Internal.GridDensityFitter;
export import qchem.Fitting.CoulombOperator;
import qchem.Fitting.Types;           // robs_t
import qchem.Fitting.FitOperations;   // OrthogonalFit: a matrix-free seed samples itself
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.Blaze;
import qchem.Parallel;                // LazyLock -- the per-block cache's guard

export namespace qchem::Fitting
{

class GridDensityFitter
    : public virtual FunctionFitter_Density_NonOrtho<double>
    , public virtual DensityProjector
{
public:
    typedef std::shared_ptr<const BasisSet::rFIT_SF_ABS> fbs_t;
    GridDensityFitter(fbs_t fbs, std::unique_ptr<CoulombOperator> op)
        : itsFitBasis(std::move(fbs)), itsOperator(std::move(op))
    {
        if (!itsFitBasis || !itsOperator) throw std::logic_error("GridDensityFitter: needs a fit basis and a Coulomb operator");
        if (!itsFitBasis->isOrtho())
            throw std::logic_error("GridDensityFitter: the fit basis '"+itsFitBasis->Name()+"' is not orthogonal, so its projection "
                                   "is not its fit -- a grid density fit needs a diagonal metric");
        itsW=itsFitBasis->Charge();
        double total=0.0;
        for (size_t g=0; g<itsW.size(); g++) total+=itsW[g];
        itsC=rvec_t(itsW.size(),1.0/total);   // the initial guess: a uniform density of unit charge, as ConstrainedFF's guess is
    }

    //! A density with a matrix contracts its D into the forward I vend; a matrix-free seed has only its own field to offer.
    void DoFit(const ProjectedDensity<double>& pd) override
    {
        if (auto* cm=dynamic_cast<const CoulombMetric_ProjectedDensity*>(&pd))
            itsC=cm->GetRepulsion3C(*this);
        else if (auto* f=dynamic_cast<const ScalarFunction<double>*>(&pd))
            itsC=OrthogonalFit(*itsFitBasis,*f);
        else
            throw std::logic_error("GridDensityFitter: this density offers neither a matrix (the Coulomb face) nor a field");
        itsV=itsOperator->Potential(Charges());
    }

    hmat_t<double> Repulsion(const robs_t<double>* bs) const override
    {
        if (itsV.size()==0) throw std::logic_error("GridDensityFitter: asked for a repulsion before any DoFit");
        const auto& dftbs=dynamic_cast<const BasisSet::Orbital_DFT_IBS<double,double>&>(*bs);
        return Integrator(dftbs).Adjoint(itsV);   // sum_g w_g v_g chi_i chi_j
    }

    double FitGetSelfRepulsion() const override
    {
        const rvec_t q=Charges();
        double e=0.0;
        for (size_t g=0; g<q.size(); g++) e+=q[g]*itsV[g];
        return e;
    }
    double Integral() const override
    {
        const rvec_t q=Charges();
        double n=0.0;
        for (size_t g=0; g<q.size(); g++) n+=q[g];
        return n;
    }
    void ReScale(double factor) override {itsC*=factor; itsV*=factor;}   // v is linear in rho

    std::ostream& Write(std::ostream& os) const override
        {return os << "GridDensityFitter (a density on " << itsFitBasis->Name() << ", first-order energy)" << std::endl;}

    // DensityProjector
    const qcMesh::MatrixForward<double>& Forward(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const override {return Integrator(orb);}
    size_t NumCoefficients() const override {return itsW.size();}

private:
    rvec_t Charges() const {return itsW*itsC;}   // Blaze: element-wise for two column vectors
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
    std::unique_ptr<CoulombOperator> itsOperator;
    rvec_t itsW;   //!< \f$\langle f_g|1\rangle\f$: the mesh weights on a \f$\delta\f$ basis
    rvec_t itsC;   //!< the fit: the density's values at the points
    rvec_t itsV;   //!< the potential at the points; empty until the first DoFit
};

} //namespace
