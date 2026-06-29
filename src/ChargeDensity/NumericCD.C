// File: ChargeDensity/NumericCD.C  A DFT-only charge density built as a sum of (fitted) atomic
// densities -- the molecular SAD seed (see doc/SCFSeedingPlan.md section 9).
//
// It is a superposition of per-atom real-space densities rho(r) (e.g. RecentredAtomicDensity).  It carries
// NO density matrix, so it is a DFT seed only: the Hamiltonian's DFT terms consume it purely through rho(r)
// (FittedVxc samples op(r); FittedVee overlap-fits op(r) onto its fit basis -- the seed path of FittedCD).
// op(r) (plus its total charge) is the ONLY state this object holds.
//
// It is therefore a PURE ScalarFunction + charge: a slim tChargeDensity (rho(r) + total charge), NOT a
// tDM_CD and NOT a ProjectedDensity_AO.  A fit has no density matrix, so it cannot (and need not) provide
// the matrix-only operations (DM contraction, SCF mixing, HF J/K) NOR its own Coulomb projection <rho|c>:
// the latter overlap-fit now lives on the FittedCD side (DoFit(ScalarFunction,charge), stage 3).  That is
// the ISP split that lets it carry no assert(false) stubs -- the DFT Fock build consumes it through the
// tChargeDensity/ScalarFunction face, and the HF terms (which DO need a matrix) cross-cast to tDM_CD and so
// simply never see this object.
module;
#include <vector>
#include <memory>
#include <cstddef>
export module qchem.ChargeDensity.NumericCD;
export import qchem.ChargeDensity;             // tChargeDensity<double>
import qchem.ScalarFunction;                   // ScalarFunction<double>

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
    virtual double  operator()(const rvec3_t&) const;
    virtual rvec3_t Gradient  (const rvec3_t&) const;

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
