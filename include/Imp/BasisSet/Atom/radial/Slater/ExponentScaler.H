// File: ExponentScaler.H  Rescale Slater exponents based in angular momentum L.
#ifndef _SlaterScaler_H_
#define _SlaterScaler_H_

#include "oml/vector.h"
namespace Slater
{

class ExponentScaler
{
public:
    typedef Vector<double> RVec;
    ExponentScaler(size_t N, double emin, double emax, size_t LMax);
    RVec   Get_es (size_t L) const;
private:        
    size_t itsN,itsLMax;
    double itsemin,itsemax;
    RVec es;
};

} //namespace

#endif // _SlaterScaler_H_
