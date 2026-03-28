// File: FittedFunctionImp.C  Implementation for Fitted Functions.
module;
#include <memory> // for std::shared_ptr
#include "blaze/Math.h"

export module qchem.FittedFunctionImp;
export import qchem.FittedFunctionClient;
export import qchem.FittedFunction;
import qchem.Orbital_DFT_IBS;
import qchem.Fit_IBS;
export import qchem.Mesh;


//--------------------------------------------------------------------------
//
//  The fit function is assumed to be real valued, but the basis set and
//  coefficients can be complex.
//
export template <class T> class FittedFunctionImp
    : public virtual FittedFunction
{
public:
    typedef std::shared_ptr<const Mesh>    mesh_t;
    typedef std::shared_ptr<const Fit_IBS> bs_t;
    
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
protected:
    virtual void   DoFitInternal(const ScalarFFClient&,double constraint=0);
    virtual void   DoFitInternal(const DensityFFClient&,double constraint=0);

    virtual vec_t<T>    FitGet2CenterOverlap  (const Fit_IBS*) const;
    virtual vec_t<T>    FitGet2CenterRepulsion(const Fit_IBS*) const;
    virtual smat_t<T>   FitGet3CenterOverlap  (const Orbital_DFT_IBS<double>*) const;
    virtual double FitGetCharge   (                    ) const;
    virtual double FitGetRepulsion(const FittedFunctionImp*) const;
    virtual double FitGetOverlap  (const FittedFunctionImp*) const;

public: //Client code needs read access to this data.
    bs_t     itsBasisSet;
    vec_t<T> itsFitCoeff;
    mesh_t   itsMesh;
};

export template <class T> class ConstrainedFF
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

export template <class T> class IntegralConstrainedFF
    : public ConstrainedFF<T>
{
public:
    typedef typename ConstrainedFF<T>::mesh_t mesh_t;
    typedef typename ConstrainedFF<T>::bs_t   bs_t;

    IntegralConstrainedFF(              );
    IntegralConstrainedFF(bs_t&, mesh_t&);
};
