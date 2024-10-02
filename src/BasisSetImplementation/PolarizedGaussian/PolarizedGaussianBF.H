// File: PolarizedGaussianBF.H  Polarized Gaussian in 3D space.
#ifndef _PolarizedGaussianBF_H_
#define _PolarizedGaussianBF_H_



#include "BasisSet.H"
#include "BasisSetImplementation/PolarizedGaussian/RadialFunction.H"
#include "Misc/Polarization.H"
#include "oml/vector3d.h"

class PolarizedGaussianBF
    : public TBasisFunction<double>
{
public:
    PolarizedGaussianBF(                                          );
    PolarizedGaussianBF(const RadialFunction*, const Polarization&);

    virtual bool   operator==(const BasisFunction&) const;
    virtual double GetNormalization() const;
    virtual double GetCharge       () const;

    virtual std::ostream&       Write(std::ostream&   ) const;
    virtual std::istream&       Read (std::istream&   )      ;
    virtual BasisFunction* Clone(           ) const;

    virtual double operator()(const RVec3&) const;
    virtual RVec3  Gradient  (const RVec3&) const;

private:
    friend class PolarizedGaussianIE;
    friend class PolarizedGaussianBS;

    virtual void Eval    (const Mesh&,     Vec&) const;
    virtual void EvalGrad(const Mesh&, Vec3Vec&) const;
    void Insert  (const RadialFunction* theRF,const Polarization& thePol);

    const RadialFunction*  itsRadial;
    Polarization           itsPol;
    double                 itsNormalization;
};

#endif //_PolarizedGaussianBF_H_
