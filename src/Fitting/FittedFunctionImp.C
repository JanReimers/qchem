// File: FittedFunctionImp.C  Implementation for Fitted Functions.
module;
#include <memory> // for std::shared_ptr
import qchem.LAParams;
export module qchem.FittedFunctionImp;
export import qchem.FittedFunctionClient;
export import qchem.FittedFunction;
import qchem.DFT_IBS;
import qchem.Fit_IBS;
export import qchem.Mesh;
import oml;


//--------------------------------------------------------------------------
//
//  The fit function is assumed to be real valued, but the basis set and
//  coefficients can be complex.
//
export template <class T> class FittedFunctionImp
    : public virtual FittedFunction
{
public:
    typedef Vector<T> Vec;
    typedef std::shared_ptr<const Mesh>    mesh_t;
    typedef std::shared_ptr<const Fit_IBS> bs_t;
    
    FittedFunctionImp(                                         );
    FittedFunctionImp(bs_t&, mesh_t&);
    ~FittedFunctionImp();
    
    virtual double DoFit           (const ScalarFFClient& )      ;
    virtual double DoFit           (const DensityFFClient& )      ;
    virtual void   ReScale         (double factor               )      ; //Fit *= factor
    virtual void   ShiftOrigin     (const RVec3& newCenter      )      ;
    virtual void   FitMixIn        (const FittedFunction&,double)      ; // this = this*(1-c) + that*c.
    virtual double FitGetChangeFrom(const FittedFunction&       ) const;

    virtual double operator()(const RVec3&) const;
    virtual RVec3  Gradient  (const RVec3&) const;

    virtual std::ostream& Write(std::ostream&) const;
protected:
    virtual double DoFitInternal(const ScalarFFClient&,double constraint=0);
    virtual double DoFitInternal(const DensityFFClient&,double constraint=0);


    virtual Vec    FitGet2CenterOverlap  (const Fit_IBS*) const;
    virtual Vec    FitGet2CenterRepulsion(const Fit_IBS*) const;
    virtual SMatrix<T>   FitGet3CenterOverlap  (const TOrbital_DFT_IBS<double>*) const;
    virtual double FitGetCharge   (                    ) const;
    virtual double FitGetRepulsion(const FittedFunctionImp*) const;
    virtual double FitGetOverlap  (const FittedFunctionImp*) const;

    virtual void  Eval(const Mesh&, Vec&) const;

public: //Client code needs read access to this data.
    bs_t     itsBasisSet;
    Vec      itsFitCoeff;
    mesh_t   itsMesh;
protected:
    LAParams itsLAParams; //Decides about matrix inversion.
};

export template <class T> class ConstrainedFF
    : public FittedFunctionImp<T>
{
    typedef FittedFunctionImp<T> Base;
    typedef typename Base::Vec Vec;
public:
    typedef typename Base::mesh_t mesh_t;
    typedef typename Base::bs_t   bs_t;

    ConstrainedFF();
    ConstrainedFF(bs_t&, const Vec&, mesh_t&  m);

    virtual double DoFit(const ScalarFFClient&);
    virtual double DoFit(const DensityFFClient&);

    virtual std::ostream& Write    (std::ostream&) const;
private:
    using Base::itsLAParams;
    Vec g,gS;
    T   gSg;
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
