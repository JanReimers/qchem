// File: FittedCDImp.C  Fitted charge density.
module;
#include <memory>
#include <cassert>
export module qchem.ChargeDensity.Imp.FittedCD;
import qchem.ChargeDensity.Types;
import qchem.ChargeDensity;            // DM_CD (cross-cast to its AO projection face)
import qchem.FittedCD;
import qchem.Fitting.FunctionFitter;   // FunctionFitter_Density (composed) + ProjectedDensity_AO (the AO face)

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
     typedef typename Fitting::FunctionFitter_Density<T>::bs_t   bs_t;
public:
    FittedCDImp(bs_t&, double totalCharge);
    FittedCDImp(const FittedCDImp&) = delete;   //!< copying would slice the fitter's constraint

    // FittedCD  (DoFit cross-casts the density to its AO face, then delegates to the COMPOSED fitter)
    virtual void      DoFit           (const DM_CD& cd)
    {
        auto* ao=dynamic_cast<const Fitting::ProjectedDensity_AO*>(&cd);
        assert(ao && "FittedCD::DoFit: a fitted (molecular) density must be a ProjectedDensity_AO");
        itsFitter->DoFit(*ao);
    }
    virtual smat_t<T> GetRepulsion    (const odftbs_t*) const;
    virtual double    GetSelfRepulsion(               ) const;  //Does GetRepulsion(*this);
    virtual FittedCD* Clone           (               ) const;

    // ScalarFunction
    virtual double  operator()(const rvec3_t& r) const {return (*itsFitter)(r);}      // No UT coverage
    virtual rvec3_t Gradient  (const rvec3_t& r) const {return itsFitter->Gradient(r);} // No UT coverage

private:
    std::unique_ptr<Fitting::FunctionFitter_Density<T>> itsFitter;   //!< COMPOSED fit (was inherited)
};

} //namespace