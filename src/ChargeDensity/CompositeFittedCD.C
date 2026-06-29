// File: ChargeDensity/CompositeFittedCD.C  A DFT-only charge density built as a sum of (fitted) atomic
// densities -- the molecular SAD seed (see doc/SCFSeedingPlan.md section 9).
//
// It is a superposition of per-atom real-space densities rho(r) (e.g. RecentredAtomicDensity).  It carries
// NO density matrix, so it is a DFT seed only: the Hamiltonian's DFT terms consume it through rho(r)
// (FittedVxc samples op(r)) and the Coulomb projection <rho|c> (FittedVee re-fits via ProjectedDensity_AO).
// The latter is derived ON DEMAND from op(r): overlap-fit onto the term's own fit basis, then Coulomb-
// project -- so op(r) is the only state this object holds.
//
// It is a tChargeDensity (rho(r) + total charge), NOT a tDM_CD: a fit has no density matrix, so it cannot
// (and need not) provide the matrix-only operations (DM contraction, SCF mixing, HF J/K).  That is exactly
// the ISP split that lets it carry no assert(false) stubs -- the DFT Fock build consumes it through the
// tChargeDensity face (op(r) + the cross-cast ProjectedDensity_AO for the Hartree projection), and the HF
// terms (which DO need a matrix) cross-cast to tDM_CD and so simply never see this object.
module;
#include <vector>
#include <memory>
#include <cstddef>
export module qchem.ChargeDensity.CompositeFittedCD;
export import qchem.ChargeDensity;             // tChargeDensity<double>
import qchem.ChargeDensity.Types;              // FIT_CD_ABS (via Fit_IBS)
import qchem.Fitting.FunctionFitter;           // ProjectedDensity_AO
import qchem.ScalarFunction;                   // ScalarFunction<double>

export namespace qchem::ChargeDensity
{

class CompositeFittedCD
    : public virtual tChargeDensity<double>
    , public virtual Fitting::ProjectedDensity_AO   // the AO projection face FittedVee::DoFit cross-casts to
{
public:
    //! \a totalCharge is the seed's electron count N (the AO-fit charge constraint).
    explicit CompositeFittedCD(double totalCharge);

    //! Add one atom's (recentred) real-space density to the superposition.
    void Insert(std::shared_ptr<const ScalarFunction<double>>);

    // ScalarFunction -- the superposed rho(r) (FittedVxc) and its gradient.
    virtual double  operator()(const rvec3_t&) const;
    virtual rvec3_t Gradient  (const rvec3_t&) const;

    // ProjectedDensity_AO -- the Coulomb projection <rho|c> (FittedVee), derived on demand from op(r).
    virtual double FitGetConstraint() const {return itsScale*itsCharge;}   // the AO fit RHS charge N
    virtual rvec_t GetRepulsion3C(const BasisSet::FIT_CD_ABS*) const;

    // tChargeDensity -- the matrix-free face (everything a DFT seed needs).
    virtual double GetTotalCharge() const {return itsScale*itsCharge;}
    virtual size_t Version       () const {return itsVersion;}
    virtual void   ReScale(double factor);

private:
    std::vector<std::shared_ptr<const ScalarFunction<double>>> itsDensities;  //!< per-atom rho(r)
    double itsCharge;        //!< total electron count N (fit constraint), pre-scale
    double itsScale=1.0;     //!< uniform scale applied by ReScale
    size_t itsVersion;       //!< transient freshness serial (newCD cache check)
};

} //namespace
