// File: ChargeDensity/CompositeFittedCD.C  A DFT-only charge density built as a sum of (fitted) atomic
// densities -- the molecular SAD seed (see doc/SCFSeedingPlan.md section 9).
//
// It is a superposition of per-atom real-space densities rho(r) (e.g. RecentredAtomicDensity).  It carries
// NO density matrix, so it is a DFT seed only: the Hamiltonian's DFT terms consume it through rho(r)
// (FittedVxc samples op(r)) and the Coulomb projection <rho|c> (FittedVee re-fits via ProjectedDensity_AO).
// The latter is derived ON DEMAND from op(r): overlap-fit onto the term's own fit basis, then Coulomb-
// project -- so op(r) is the only state this object holds.
//
// It nominally satisfies tDM_CD (the type the whole Hamiltonian framework is keyed on), but the genuine
// density-MATRIX capabilities (HF exchange/direct, DM contraction, SCF mixing) assert(false): a fit has no
// DM.  The clean fix is the deferred tChargeDensity/tDM_CD ISP split (option 1) -- a DFT-only base the
// framework takes, with the HF methods on the tDM_CD derived face.  Until then these stubs are unreachable
// on a seed (the seed is consumed once, at iteration 0, by the DFT Fock build only).
module;
#include <vector>
#include <memory>
#include <cstddef>
export module qchem.ChargeDensity.CompositeFittedCD;
export import qchem.ChargeDensity;             // tDM_CD<double>, tStatic_CC, tDynamic_CC
import qchem.ChargeDensity.Types;              // ohfbs_t, FIT_CD_ABS (via Fit_IBS)
import qchem.Fitting.FunctionFitter;           // ProjectedDensity_AO
import qchem.ScalarFunction;                   // ScalarFunction<double>

export namespace qchem::ChargeDensity
{

class CompositeFittedCD
    : public virtual tDM_CD<double>
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
    virtual double FitGetConstraint() const {return itsCharge;}      // the AO fit RHS charge N
    virtual rvec_t GetRepulsion3C(const BasisSet::FIT_CD_ABS*) const;

    // tDM_CD scalars a DFT seed legitimately answers.
    virtual double GetTotalCharge() const {return itsCharge;}
    virtual size_t Version       () const {return itsVersion;}

    // --- density-MATRIX capabilities: a fitted density has no DM (see the tChargeDensity ISP TODO). ---
    virtual double DM_Contract(const tStatic_CC<double>*) const;
    virtual double DM_Contract(const tDynamic_CC<double>*,const tDM_CD<double>*) const;
    virtual void   ReScale      (double)                         ;
    virtual void   MixIn        (const tDM_CD<double>&,double)   ;
    virtual double GetChangeFrom(const tDM_CD<double>&) const    ;
    virtual void   AccumulateDirect  (hmat_t<double>&, const ohfbs_t*) const;
    virtual void   AccumulateExchange(hmat_t<double>&, const ohfbs_t*) const;

private:
    std::vector<std::shared_ptr<const ScalarFunction<double>>> itsDensities;  //!< per-atom rho(r)
    double itsCharge;     //!< total electron count N (fit constraint)
    size_t itsVersion;    //!< transient freshness serial (newCD cache check)
};

} //namespace
