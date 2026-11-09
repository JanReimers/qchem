// File: ChargeDensity/NumericCD.C  A DFT-only charge density built as a sum of (fitted) atomic
// densities -- the molecular SAD seed (see doc/SCFSeedingPlan.md section 9).
//
// It is a superposition of per-atom real-space densities rho(r) (e.g. RecentredAtomicDensity).  It carries
// NO density matrix, so it is a DFT seed only: the Hamiltonian's DFT terms consume it purely through rho(r)
// (FittedVxc samples op(r); a density fitter projects it through its own ProjectField).
// op(r) (plus its total charge) is the ONLY state this object holds.
//
// It is a slim tChargeDensity (rho(r) + total charge), NOT a tDM_CD: it has no density matrix, so it cannot
// (and need not) provide the matrix-only operations (DM contraction, SCF mixing, HF J/K) -- the DFT Fock build
// consumes it through the tChargeDensity/ScalarFunction face, and the HF terms (which DO need a matrix) cross-
// cast to tDM_CD and so simply never see this object.  A density fit reaches it through tChargeDensity::Project's
// default, p.ProjectField(*this): the fitter samples the field in its own metric (V-ONE-OF-N 3c; the old
// overlap-metric face and GetUnconstrainedFit are gone, and the Dunlap fitter's G^-1 G c0 round trip is the price).
module;
#include <vector>
#include <memory>
#include <cstddef>
export module qchem.ChargeDensity.NumericCD;
export import qchem.ChargeDensity;             // tChargeDensity<double>
import qchem.ScalarFunction;                   // ScalarFunction<double>
import qchem.BasisSet.Orbital_DFT_IBS;                 // rFIT_CD_ABS (the fitter's narrow density-fit-basis face)

export namespace qchem::ChargeDensity
{

class NumericCD
    : public virtual tChargeDensity<double>
{
public:
    //! \a totalCharge is the seed's electron count N (the fit charge constraint).
    explicit NumericCD(double totalCharge);

    //! Add one atom's (recentred) real-space density to the superposition.
    void Insert(std::shared_ptr<const ScalarFunction<double>>);

    // ScalarFunction -- the superposed rho(r) (FittedVxc samples it; FittedVee overlap-fits it) and gradient.
    virtual double  operator()(const rvec3_t&) const override;
    virtual rvec3_t Gradient  (const rvec3_t&) const override;

    // No matrix, so a fitter asks for the seed's FIELD: tChargeDensity::Project's default, p.ProjectField(*this) (V-ONE-OF-N 3c).

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
