// File: FittedCDImp.C  Fitted charge density.
module;
#include <memory>
export module qchem.ChargeDensity.Imp.FittedCD;
import qchem.ChargeDensity.Types;
import qchem.FittedCD;
import qchem.FittedFunctionImp;

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
     typedef typename Fitting::FunctionFitter<T>::mesh_t mesh_t;
     typedef typename Fitting::FunctionFitter<T>::bs_t   bs_t;
public:
    FittedCDImp(bs_t&, mesh_t&, double totalCharge);
    FittedCDImp(const FittedCDImp&) = delete;   //!< copying would slice the fitter's constraint

    // FittedCD
    virtual smat_t<T> GetRepulsion    (const odftbs_t*) const;
    virtual double    GetSelfRepulsion(               ) const;  //Does GetRepulsion(*this);
    virtual FittedCD* Clone           (               ) const;

    // FittedFunction -- delegate to the COMPOSED fitter (was inherited from IntegralConstrainedFF).
    virtual void   DoFit           (const Fitting::ScalarFFClient&  c)      {itsFitter->DoFit(c);}
    virtual void   DoFit           (const Fitting::DensityFFClient& c)      {itsFitter->DoFit(c);}
    virtual void   ReScale         (double factor)                         {itsFitter->ReScale(factor);}
    virtual void   FitMixIn        (const Fitting::FittedFunction& g,double f);
    virtual double FitGetChangeFrom(const Fitting::FittedFunction& g) const;

    // ScalarFunction
    virtual double  operator()(const rvec3_t& r) const {return (*itsFitter)(r);}      // No UT coverage
    virtual rvec3_t Gradient  (const rvec3_t& r) const {return itsFitter->Gradient(r);} // No UT coverage

private:
    std::unique_ptr<Fitting::FunctionFitter<T>> itsFitter;   //!< COMPOSED fit (was inherited)
};

} //namespace