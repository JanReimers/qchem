// File: FunctionFitterImp.C  Concrete least-squares fitters implementing the two Fitting faces.
module;
#include <memory> // for std::shared_ptr
export module qchem.Fitting.Internal.FunctionFitterImp;
export import qchem.Fitting.FunctionFitter;  // FunctionFitter_Scalar/_Density, the *FFClients, ScalarFunction
import qchem.Fitting.Types;
import qchem.Blaze;
//--------------------------------------------------------------------------
//
//  The fit function and fit basis set are assumed to be real valued.
//  But the orbital basis set and coefficients can be complex (T).
//
export namespace qchem::Fitting
{

//! \brief Shared implementation of the coefficient + real-space machinery common to BOTH fitter faces.
//! Parametrised on the public face it implements (Face = FunctionFitter_Scalar/_Density<T>) and the narrow
//! fit-basis face it holds (FBS = FIT_SF_ABS/FIT_CD_ABS).  Carries the fit coefficients and the shared
//! ScalarFunction (operator()/Gradient), ReScale, FitMixIn/FitGetChangeFrom (coefficient-only) and Write;
//! the metric-specific DoFit + contraction live in the leaf impls below.
template <class T, class Face, class FBS> class FitImpBase
    : public virtual Face
{
public:
    typedef std::shared_ptr<const FBS> fbs_t;

    FitImpBase(     ) : itsBasisSet( ), itsFitCoeff( ) {}
    FitImpBase(fbs_t& fbs) : itsBasisSet(fbs), itsFitCoeff(fbs->GetNumFunctions(),0.0) {}

    virtual void   ReScale         (double factor)            override; // Fit *= factor
    virtual void   FitMixIn        (const Face&,double)        override; // this = this*(1-c) + that*c
    virtual double FitGetChangeFrom(const Face&) const         override;

    virtual double  operator()(const rvec3_t&) const          override;
    virtual rvec3_t Gradient  (const rvec3_t&) const          override;
    virtual std::ostream& Write(std::ostream&) const          override;

public: // Client code needs read access to this data.
    fbs_t     itsBasisSet;
    vec_t<T> itsFitCoeff;
};

//---------------------------------------------------------------------- Scalar (overlap-metric) impl
template <class T> class FunctionFitterImp
    : public FitImpBase<T, FunctionFitter_Scalar<T>, BasisSet::FIT_SF_ABS>
{
    typedef FitImpBase<T, FunctionFitter_Scalar<T>, BasisSet::FIT_SF_ABS> Base;
public:
    typedef typename Base::fbs_t fbs_t;

    FunctionFitterImp(     ) : Base( ) {}
    FunctionFitterImp(fbs_t& fbs) : Base(fbs) {}

    virtual void      DoFit  (const ScalarFFClient&)      override;  // overlap-metric projection
    virtual hmat_t<T> Overlap(const obs_t<T>*) const      override;  // Sum_a c_a <Oi|f_a|Oj>
};

//---------------------------------------------------------------- Density (Coulomb-metric) impl
template <class T> class ConstrainedFF
    : public FitImpBase<T, FunctionFitter_Density<T>, BasisSet::FIT_CD_ABS>
{
    typedef FitImpBase<T, FunctionFitter_Density<T>, BasisSet::FIT_CD_ABS> Base;
public:
    typedef typename Base::fbs_t fbs_t;

    ConstrainedFF();
    ConstrainedFF(fbs_t&, const vec_t<T>& g);

    virtual void      DoFit    (const ProjectedDensity_AO&)  override; // Dunlap charge-constrained fit
    virtual hmat_t<T> Repulsion(const obs_t<T>*) const       override; // Sum_a c_a <Oi|f_a/r12|Oj>
    virtual double    FitGetSelfRepulsion() const            override; // <fit|1/r12|fit>
    virtual double    Integral () const                      override;

    virtual std::ostream& Write(std::ostream&) const         override;
protected:
    //! Unconstrained Coulomb-metric solve c0 = J^-1 <rho|f> (the constraint is applied on top in DoFit).
    void   DoFitUnconstrained(const ProjectedDensity_AO&);
    //! Coulomb repulsion energy with another fit (self-repulsion uses *this).
    double FitGetRepulsion(const ConstrainedFF*) const;
private:
    vec_t<T> g;
    row_t<T> gS;
    T        gSg;
};

template <class T> class IntegralConstrainedFF
    : public ConstrainedFF<T>
{
public:
    typedef typename ConstrainedFF<T>::fbs_t   fbs_t;

    IntegralConstrainedFF(              );
    IntegralConstrainedFF(fbs_t&);
};

} //namespace
