// File: SphericalGaussianBF.H  r^l exp(-ar^2) type basis function for an Spherical.
#ifndef _SphericalGaussianBF_H_
#define _SphericalGaussianBF_H_



#include "BasisSet.H"

class SphericalGaussianBF
    : public virtual TBasisFunction<double>
{
public:
    SphericalGaussianBF(                             );
    SphericalGaussianBF(double theExponent, int theL);
    virtual void Init(double norm, double charge); //Use an integral engine to calculate normalization and charge constants.
    
    virtual bool   operator==(const BasisFunction&) const;
    virtual double GetNormalization() const;
    virtual double GetCharge       () const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;
    virtual BasisFunction* Clone(        ) const;

    virtual double operator()(const Vec3&) const;
    virtual Vec3   Gradient  (const Vec3&) const;

private:
    friend class SphericalGaussianIE;
    friend class SphericalGaussianIE1;
    friend class SphericalSymmetryBS;

    double itsExponent;
    int    itsL;
    double itsNormalization;
    double itsCharge;
};

#endif //_SphericalGaussianBF_H_
