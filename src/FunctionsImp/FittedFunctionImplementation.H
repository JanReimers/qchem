// File: FittedFunctionImplementation.H  General Fitted Function.
#ifndef _FittedFunctionImplementation_H_
#define _FittedFunctionImplementation_H_



#include "BasisSet/TBasisSet.H"
#include "Misc/rc_ptr.H"
#include "Functions/FittedFunction.H"
#include "oml/vector.h"

template <class T> class TBasisSetImplementation;

//--------------------------------------------------------------------------
//
//  The fit function is assumed to be real valued, but the basis set and
//  coefficients can be complex.
//
template <class T> class FittedFunctionImplementation
    : public virtual FittedFunction
{
public:
    typedef Vector<T> Vec;

    FittedFunctionImplementation(                                         );
    FittedFunctionImplementation(const rc_ptr<BasisSet>&, bool CDfit=false);

    virtual double DoFit           (const FittedFunctionClient& )      ;
    virtual void   ReScale         (double factor               )      ; //Fit *= factor
    virtual void   ShiftOrigin     (const RVec3& newCenter      )      ;
    virtual void   FitMixIn        (const FittedFunction&,double)      ; // this = this*(1-c) + that*c.
    virtual double FitGetChangeFrom(const FittedFunction&       ) const;

    virtual double operator()(const RVec3&) const;
    virtual RVec3  Gradient  (const RVec3&) const;

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

    Vec&    GetFitCoeff()
    {
        return itsFitCoeff;
    }

protected:
    virtual double DoFit(double, const FittedFunctionClient&);
    virtual double DoFit(double, const Vec&);

    virtual Vec    FitGet2CenterOverlap  (const BasisSet*) const;
    virtual Vec    FitGet2CenterRepulsion(const BasisSet*) const;
    virtual double FitGetCharge   (                    ) const;
    virtual double FitGetRepulsion(const FittedFunctionImplementation*) const;
    virtual double FitGetOverlap  (const FittedFunctionImplementation*) const;

    virtual void  Eval(const Mesh&, Vec&) const;

    void    SetFitCoeff(const Vec& fc);
    SMat    GetInverseOverlap() const; //Inverse of the normalized fit basis overlap or repulsion.

private:
    friend class TBasisSetImplementation<T>;
    const TBasisSet<T>* CastBasisSet() const
    {
        return dynamic_cast<TBasisSet<T>*>(itsBasisSet.get());
    }

    rc_ptr<BasisSet> itsBasisSet;
    Vec              itsFitCoeff;
    bool             itsCDFitFlag;
};

#endif //_FittedFunctionImplementation_H_

