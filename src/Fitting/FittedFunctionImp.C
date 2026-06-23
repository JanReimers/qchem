// File: FittedFunctionImp.C  Implementation for Fitted Functions.
module;
#include <memory> // for std::shared_ptr
export module qchem.FittedFunctionImp;
export import qchem.FittedFunctionClient;
export import qchem.FittedFunction;
export import qchem.Mesh;
import qchem.Fitting.Types;
import qchem.Blaze;
//--------------------------------------------------------------------------
//
//  The fit function and fit basis set are assumed to be real valued. 
//  But the orbtial basis set and coefficients can be complex. 
//
export namespace qchem::Fitting
{

template <class T> class FittedFunctionImp
    : public virtual FittedFunction
{
public:
    typedef std::shared_ptr<const Mesh>  mesh_t;
    typedef std::shared_ptr<const fbs_t> bs_t;
    
    FittedFunctionImp(                                         );
    FittedFunctionImp(bs_t&, mesh_t&);
    ~FittedFunctionImp();
    
    virtual void   DoFit           (const ScalarFFClient& )      ;
    virtual void   DoFit           (const DensityFFClient& )      ;
    virtual void   ReScale         (double factor               )      ; //Fit *= factor
    virtual void   FitMixIn        (const FittedFunction&,double)      ; // this = this*(1-c) + that*c.
    virtual double FitGetChangeFrom(const FittedFunction&       ) const;

    virtual double operator()(const rvec3_t&) const;
    virtual rvec3_t  Gradient  (const rvec3_t&) const;

    virtual std::ostream& Write(std::ostream&) const;

    // Fit-derived quantities.  PUBLIC so client code can COMPOSE a fitter (hold one and call these)
    // rather than INHERIT FittedFunctionImp -- see FunctionFitter alias below.
    virtual vec_t<T>    FitGet2CenterOverlap  (const fbs_t*) const;
    virtual vec_t<T>    FitGet2CenterRepulsion(const fbs_t*) const;
    virtual smat_t<T>   FitGet3CenterOverlap  (const obs_t<T>*) const;
    virtual double FitGetCharge   (                    ) const;
    virtual double FitGetRepulsion(const FittedFunctionImp*) const;
    virtual double FitGetOverlap  (const FittedFunctionImp*) const;
protected:
    virtual void   DoFitInternal(const ScalarFFClient&,double constraint=0);
    virtual void   DoFitInternal(const DensityFFClient&,double constraint=0);

public: //Client code needs read access to this data.
    bs_t     itsBasisSet;
    vec_t<T> itsFitCoeff;
    mesh_t   itsMesh;
};

//! \brief A least-squares function fitter that client code COMPOSES (holds and delegates to) rather than
//! inheriting.  Currently an alias for the unconstrained FittedFunctionImp; the constraint variants
//! (ConstrainedFF/IntegralConstrainedFF) are dormant -- their constraint is computed but never applied
//! (Imp/ConstrainedFF.C).  Step 2 reintroduces a real charge constraint per Dunlap, Connolly & Sabin,
//! J. Chem. Phys. 71, 3396 (1979).
template <class T> using FunctionFitter = FittedFunctionImp<T>;

template <class T> class ConstrainedFF
    : public FittedFunctionImp<T>
{
    typedef FittedFunctionImp<T> Base;
public:
    typedef typename Base::mesh_t mesh_t;
    typedef typename Base::bs_t   bs_t;

    ConstrainedFF();
    ConstrainedFF(bs_t&, const vec_t<T>& g, mesh_t&  m);

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
    typedef typename ConstrainedFF<T>::mesh_t mesh_t;
    typedef typename ConstrainedFF<T>::bs_t   bs_t;

    IntegralConstrainedFF(              );
    IntegralConstrainedFF(bs_t&, mesh_t&);
};

} //namespace