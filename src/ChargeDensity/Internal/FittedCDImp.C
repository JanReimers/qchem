// File: FittedCDImp.C  Fitted charge density.
module;
#include <memory>
#include <cassert>
export module qchem.ChargeDensity.Imp.FittedCD;
import qchem.ChargeDensity.Types;
import qchem.ChargeDensity;            // rDM_CD (cross-cast to its AO projection face)
import qchem.FittedCD;
import qchem.Fitting.FunctionFitter;   // FunctionFitter_Density (composed) + ProjectedDensity_AO (the AO face)
import qchem.Outcome;
import qchem.AD;
import qchem.Structure;

export namespace qchem::ChargeDensity
{

//---------------------------------------------------------------------------------
//
//  A charge density implemented by fitting the real charge density to an
//  auxillary basis set.
//
template <class T> class FittedCDImp
    : public virtual FittedCD
{
    typedef std::shared_ptr<const BasisSet::rFIT_CD_ABS> fbs_t;   //!< Coulomb-metric aux basis (narrow face)
public:
    typedef std::unique_ptr<Fitting::FunctionFitter_Density_NonOrtho<T>> fitter_t;
    FittedCDImp(fbs_t&, double totalCharge);
    FittedCDImp(fitter_t, double totalCharge);
    FittedCDImp(const FittedCDImp&) = delete;   //!< copying would slice the fitter's constraint

    // FittedCD  (every finite density -- matrix-backed or matrix-free seed -- fits through its own AO face)
    virtual void      DoFit           (const rChargeDensity&             );
    virtual smat_t<T> GetRepulsion    (const odftbs_t*) const;
    virtual double    GetSelfRepulsion(               ) const;  //Does GetRepulsion(*this);
    virtual double    GetChargeError  (               ) const;
    virtual bool CanPullbackRepulsion    (const odftbs_t*) const;
    virtual bool CanPullbackSelfRepulsion(               ) const;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsionPullback    (const odftbs_t*, const hmat_t<double>& D, const Structure*) const;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeSelfRepulsionPullback(const Structure*) const;

private:
    std::unique_ptr<Fitting::FunctionFitter_Density_NonOrtho<T>> itsFitter;   //!< COMPOSED non-ortho fit (was inherited)
    double itsTotalCharge;
};

} //namespace