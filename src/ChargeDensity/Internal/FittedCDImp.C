// File: FittedCDImp.C  Fitted charge density.
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
    , public         Fitting::IntegralConstrainedFF<T> //Pick up the DoFit method.
{
     typedef typename Fitting::IntegralConstrainedFF<T>::mesh_t mesh_t;
     typedef typename Fitting::IntegralConstrainedFF<T>::bs_t   bs_t;
public:
    FittedCDImp(bs_t&, mesh_t&, double totalCharge);

    virtual smat_t<T> GetRepulsion    (const odftbs_t*) const;
    virtual double    GetSelfRepulsion(                      ) const;  //Does GetRepulsion(*this);
    
    virtual double  operator()(const rvec3_t&) const; // No UT coverage
    virtual rvec3_t Gradient  (const rvec3_t&) const; // No UT coverage

    virtual FittedCD* Clone(        ) const;

private:
    using Fitting::FittedFunctionImp<T>::FitGetCharge;
    using Fitting::FittedFunctionImp<T>::FitGet2CenterOverlap;
    using Fitting::FittedFunctionImp<T>::FitGet2CenterRepulsion;
    using Fitting::FittedFunctionImp<T>::itsFitCoeff;
    using Fitting::FittedFunctionImp<T>::itsBasisSet;
    
};

} //namespace