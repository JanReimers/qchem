// File: ChargeDensity/NumericCD.C  A DFT-only charge density built as a sum of (fitted) atomic
// densities -- the molecular SAD seed (see doc/SCFSeedingPlan.md section 9).
//
// It is a superposition of per-atom real-space densities rho(r) (e.g. RecentredAtomicDensity).  It carries
// NO density matrix, so it is a DFT seed only: the Hamiltonian's DFT terms consume it purely through rho(r)
// (FittedVxc samples op(r); FittedVee overlap-fits op(r) onto its fit basis -- the seed path of FittedCD).
// op(r) (plus its total charge) is the ONLY state this object holds.
//
// It is a slim tChargeDensity (rho(r) + total charge), NOT a tDM_CD: it has no density matrix, so it cannot
// (and need not) provide the matrix-only operations (DM contraction, SCF mixing, HF J/K) -- the DFT Fock build
// consumes it through the tChargeDensity/ScalarFunction face, and the HF terms (which DO need a matrix) cross-
// cast to tDM_CD and so simply never see this object.  It DOES, however, own its OWN density-fit projection:
// like every finite density it IS-A Fitting::ProjectedDensity_AO, answering the fitter's <rho|c> query by
// overlap-fitting its own rho(r) (a numeric seed has no matrix, so <rho|c> is the S-metric fit e=S^-1<f|rho>,
// Coulomb-projected).  This mirrors how the plane-wave seed (FourierSeedCD) owns its GetFourierDensity, and
// keeps FittedCD ignorant that seeding exists.  (Reverses the stage-3 note: the overlap-fit is back on the
// seed, as a genuine projection, not a FittedCD stub.  The CD->SF cross-cast inside GetRepulsion3C is item F.)
module;
#include <vector>
#include <memory>
#include <cstddef>
export module qchem.ChargeDensity.NumericCD;
export import qchem.ChargeDensity;             // tChargeDensity<double>
import qchem.ScalarFunction;                   // ScalarFunction<double>
import qchem.Fitting.FunctionFitter;           // ProjectedDensity_AO (the seed's own density-fit projection)
import qchem.BasisSet.Fit_IBS;                 // rFIT_CD_ABS (the fitter's narrow density-fit-basis face)

export namespace qchem::ChargeDensity
{

class NumericCD
    : public virtual tChargeDensity<double>
    , public virtual Fitting::ProjectedDensity_AO   // the seed owns its own density-fit projection <rho|c>
{
public:
    //! \a totalCharge is the seed's electron count N (the fit charge constraint).
    explicit NumericCD(double totalCharge);

    //! Add one atom's (recentred) real-space density to the superposition.
    void Insert(std::shared_ptr<const ScalarFunction<double>>);

    // ScalarFunction -- the superposed rho(r) (FittedVxc samples it; FittedVee overlap-fits it) and gradient.
    virtual double  operator()(const rvec3_t&) const override;
    virtual rvec3_t Gradient  (const rvec3_t&) const override;

    // ProjectedDensity_AO -- the seed's own density-fit projection: no matrix, so <rho|c> is an overlap-metric
    // fit of rho(r), Coulomb-projected; the charge constraint is the seed's total charge.
    virtual double FitGetConstraint() const override {return GetTotalCharge();}          //!< the fit RHS charge N
    virtual rvec_t GetRepulsion3C(const BasisSet::rFIT_CD_ABS*) const override;          //!< <rho|c> (overlap-fit)

    // tChargeDensity -- the matrix-free face (everything a DFT seed needs).
    virtual double GetTotalCharge() const override {return itsScale*itsCharge;}
    virtual size_t Version       () const override {return itsVersion;}
    virtual void   ReScale(double factor) override;

private:
    std::vector<std::shared_ptr<const ScalarFunction<double>>> itsDensities;  //!< per-atom rho(r)
    double itsCharge;        //!< total electron count N (fit constraint), pre-scale
    double itsScale=1.0;     //!< uniform scale applied by ReScale
    size_t itsVersion;       //!< transient freshness serial (newCD cache check)
};

} //namespace
