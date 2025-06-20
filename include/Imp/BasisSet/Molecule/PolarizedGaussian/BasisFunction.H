// File: PolarizedGaussianBF.H  Polarized Gaussian in 3D space.
#ifndef _PolarizedGaussianBF_H_
#define _PolarizedGaussianBF_H_


#include <BasisSet/BasisFunction.H>
#include "RadialFunction.H"
#include "Polarization.H"
#include "oml/vector3d.h"

namespace PolarizedGaussian
{

class BasisFunction
    : public TBasisFunction<double>
{
public:
    BasisFunction(                                          );
    BasisFunction(const RadialFunction*, const Polarization&, double norm);

    virtual bool   operator==(const ::BasisFunction&) const;
 
    virtual std::ostream&       Write(std::ostream&   ) const;
    virtual std::istream&       Read (std::istream&   )      ;
    virtual BasisFunction* Clone(           ) const;

    virtual double operator()(const RVec3&) const;
    virtual RVec3  Gradient  (const RVec3&) const;

private:
    void Insert  (const RadialFunction* theRF,const Polarization& thePol);

    const RadialFunction*  itsRadial;
    Polarization           itsPol;
    double                 itsNormalization;
};

} //namespace PolarizedGaussian

#endif //_PolarizedGaussianBF_H_
