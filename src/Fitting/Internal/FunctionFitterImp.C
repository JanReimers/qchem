// File: FunctionFitterImp.C  Concrete least-squares fitter implementing Fitting::FunctionFitter.
module;
#include <memory> // for std::shared_ptr
export module qchem.Fitting.Internal.FunctionFitterImp;
export import qchem.Fitting.FunctionFitter;  // FunctionFitter<T>, ScalarFFClient, DensityFFClient, ScalarFunction
import qchem.Fitting.Types;
import qchem.Blaze;
//--------------------------------------------------------------------------
//
//  The fit function and fit basis set are assumed to be real valued.
//  But the orbital basis set and coefficients can be complex (T).
//
export namespace qchem::Fitting
{

template <class T> class FunctionFitterImp
    : public virtual FunctionFitter<T>
{
public:
    typedef std::shared_ptr<const fbs_t> bs_t;

    FunctionFitterImp(                                         );
    FunctionFitterImp(bs_t&);
    ~FunctionFitterImp();

    virtual void   DoFit           (const ScalarFFClient& )      ;
    virtual void   DoFit           (const DensityFFClient& )      ;
    virtual void   ReScale         (double factor               )      ; //Fit *= factor
    virtual void   FitMixIn        (const FunctionFitter<T>&,double)      ; // this = this*(1-c) + that*c.
    virtual double FitGetChangeFrom(const FunctionFitter<T>&       ) const;

    virtual double  operator()(const rvec3_t&) const;
    virtual rvec3_t Gradient  (const rvec3_t&) const;
    virtual std::ostream& Write(std::ostream&) const;

    virtual smat_t<T> FitGet3CenterOverlap  (const obs_t<T>*) const;
    virtual smat_t<T> FitGet3CenterRepulsion(const obs_t<T>*) const;
    virtual double    FitGetSelfRepulsion   ()                const;  // <fit|1/r12|fit>
    virtual double    FitGetCharge          ()                const;
protected:
    virtual void   DoFitInternal(const ScalarFFClient&,double constraint=0);
    virtual void   DoFitInternal(const DensityFFClient&,double constraint=0);
    //! Coulomb repulsion energy with another fit (self-repulsion uses *this).
    double FitGetRepulsion(const FunctionFitterImp*) const;

public: //Client code needs read access to this data.
    bs_t     itsBasisSet;
    vec_t<T> itsFitCoeff;
};

template <class T> class ConstrainedFF
    : public FunctionFitterImp<T>
{
    typedef FunctionFitterImp<T> Base;
public:
    typedef typename Base::bs_t   bs_t;

    ConstrainedFF();
    ConstrainedFF(bs_t&, const vec_t<T>& g);

    virtual void   DoFit(const ScalarFFClient&);
    virtual void   DoFit(const DensityFFClient&);

    virtual std::ostream& Write    (std::ostream&) const;
private:
    vec_t<T> g;
    row_t<T> gS;
    T        gSg;
};

template <class T> class IntegralConstrainedFF
    : public ConstrainedFF<T>
{
public:
    typedef typename ConstrainedFF<T>::bs_t   bs_t;

    IntegralConstrainedFF(              );
    IntegralConstrainedFF(bs_t&);
};

} //namespace
