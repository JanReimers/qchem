// File: PlaneWaveBF.H  Polarized Gaussian in 3D space.
#ifndef _PlaneWaveBF_H_
#define _PlaneWaveBF_H_


#include <BasisFunction.H>
#include "oml/vector3d.h"

class PlaneWaveBF
    : public TBasisFunction<std::complex<double> >
{
public:
    PlaneWaveBF(                             );
    PlaneWaveBF(const RVec3& G, double Volume);
    virtual void Init(double norm, double charge); //Use an integral engine to calculate normalization and charge constants.

    virtual bool   operator==(const BasisFunction&) const;
    virtual double GetNormalization() const;
    virtual double GetCharge       () const; //Zero, unless G=0.

    
    virtual std::ostream&       Write(std::ostream&   ) const;
    virtual std::istream&       Read (std::istream&   )      ;
    virtual BasisFunction* Clone(           ) const;

    virtual std::complex<double>  operator()(const RVec3&) const;
    virtual Vec3             Gradient  (const RVec3&) const;

private:
    friend class PlaneWaveIE;
    friend class PlaneWaveBS;

    virtual void Eval    (const Mesh&,     Vec&) const;
    virtual void EvalGrad(const Mesh&, Vec3Vec&) const;

    RVec3   itsG;             //RL vector.
    double  itsNormalization; //Born Von-Karmen lattice volume
    double  itsCharge;
};

#endif //_PlaneWaveBF_H_
