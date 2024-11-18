// File: SlaterScaler.H  Rescale Slater exponents based in angular momentum L.
#ifndef _SlaterScaler_H_
#define _SlaterScaler_H_

#include "oml/vector.h"

class SlaterScaler
{
public:
    typedef Vector<double> RVec;
    SlaterScaler(size_t N, double emin, double emax, size_t LMax);
    RVec   Get_es (size_t L) const;
private:        
    size_t itsN,itsLMax;
    double itsemin,itsemax;
    RVec es;
};
#endif // _SlaterScaler_H_
