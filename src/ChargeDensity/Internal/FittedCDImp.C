// File: FittedCDImp.C  Fitted charge density.
module;
#include <memory>
export module qchem.ChargeDensity.Imp.FittedCD;
import qchem.ChargeDensity.Types;
import qchem.FittedCD;
import qchem.Fitting.FunctionFitter;   // Fitting::FunctionFitter (composed, via the Factory)

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
     typedef typename Fitting::FunctionFitter<T>::bs_t   bs_t;
public:
    FittedCDImp(bs_t&, double totalCharge);
    FittedCDImp(const FittedCDImp&) = delete;   //!< copying would slice the fitter's constraint

    // FittedCD  (DoFit delegates to the COMPOSED fitter; the fitter answers the energy queries)
    virtual void      DoFit           (const Fitting::ProjectedDensity_AO& c)      {itsFitter->DoFit(c);}
    virtual smat_t<T> GetRepulsion    (const odftbs_t*) const;
    virtual double    GetSelfRepulsion(               ) const;  //Does GetRepulsion(*this);
    virtual FittedCD* Clone           (               ) const;

    // ScalarFunction
    virtual double  operator()(const rvec3_t& r) const {return (*itsFitter)(r);}      // No UT coverage
    virtual rvec3_t Gradient  (const rvec3_t& r) const {return itsFitter->Gradient(r);} // No UT coverage

private:
    std::unique_ptr<Fitting::FunctionFitter<T>> itsFitter;   //!< COMPOSED fit (was inherited)
};

} //namespace