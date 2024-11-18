// File: GaussianScaler.H  Rescale Gaussian exponents based in angular momentum L.
#ifndef _GaussianScaler_H_
#define _GaussianScaler_H_

#include "oml/vector.h"

class GaussianScaler
{
public:
    typedef Vector<double> RVec;
    GaussianScaler(size_t N, double emin, double emax, size_t LMax);
    RVec   Get_es (size_t L) const;
private:        
    size_t itsN,itsLMax;
    double itsemin,itsemax;
    RVec es;
};
#endif // _GaussianScaler_H_
