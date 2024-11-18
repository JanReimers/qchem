// File: GaussianScaler.H  Rescale Gaussian exponents based in angular momentum L.
#ifndef _GaussianScaler_H_
#define _GaussianScaler_H_

#include "oml/vector.h"

class GaussianScaler
{
public:
    typedef Vector<double> RVec;
    GaussianScaler(size_t N, double emin, double emax, size_t LMax);
    size_t N   (size_t L) const {return (itsN*(itsLMax+1-L))/3;}
    double emin(size_t L) const {return itsemin;}
    double emax(size_t L) const {return itsemax/pow(50,L);}
    RVec   Get_es (size_t L) const;
private:        
    size_t itsN,itsLMax;
    double itsemin,itsemax;
    RVec es;
};
#endif // _GaussianScaler_H_
