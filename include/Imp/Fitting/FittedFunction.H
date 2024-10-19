// File: FittedFunctionImp.H  Common imp for Fitted Functions.
#ifndef _FittedFunctionImp_H_
#define _FittedFunctionImp_H_



#include <BasisSet.H>
#include "Misc/rc_ptr.H"
#include <FittedFunction.H>
#include "Mesh/Mesh.H"
#include "Imp/LAParams.H"
#include "oml/vector.h"
#include "oml/smatrix.h"

template <class T> class FittedCDImp;
template <class T> class ConstrainedFF;

//--------------------------------------------------------------------------
//
//  The fit function is assumed to be real valued, but the basis set and
//  coefficients can be complex.
//
template <class T> class FittedFunctionImp
    : public virtual FittedFunction
{
public:
    typedef Vector<T> Vec;

    FittedFunctionImp(                                         );
    FittedFunctionImp(const rc_ptr<IrrepBasisSet>&, const rc_ptr<Mesh>&);
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
    virtual std::istream& Read (std::istream&)      ;

protected:
    virtual double DoFitInternal(const ScalarFFClient&,double constraint=0);
    virtual double DoFitInternal(const DensityFFClient&,double constraint=0);


    virtual Vec    FitGet2CenterOverlap  (const IrrepBasisSet*) const;
    virtual Vec    FitGet2CenterRepulsion(const IrrepBasisSet*) const;
    virtual SMat   FitGet3CenterOverlap  (const IrrepBasisSet*) const;
    virtual double FitGetCharge   (                    ) const;
    virtual double FitGetRepulsion(const FittedFunctionImp*) const;
    virtual double FitGetOverlap  (const FittedFunctionImp*) const;

    virtual void  Eval(const Mesh&, Vec&) const;

    //SMat    GetInverseOverlap() const; //Inverse of the normalized fit basis overlap or repulsion.

private:
    friend class FittedCDImp<T>;
    friend class ConstrainedFF<T>;
    const TIrrepBasisSet<T>* CastBasisSet() const
    {
        return dynamic_cast<TIrrepBasisSet<T>*>(itsBasisSet.get());
    }

    rc_ptr<IrrepBasisSet> itsBasisSet;
    Vec   itsFitCoeff;
    SMat  itsInvOvlp;
    SMat  itsInvRepl;
    rc_ptr<Mesh> itsMesh;
    LAParams itsLAParams; //Decides about matrix inversion.
};

#endif //_FittedFunctionImp_H_

