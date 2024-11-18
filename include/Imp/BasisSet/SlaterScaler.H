// File: SlaterScaler.H  Rescale Slater exponents based in angular momentum L.
#ifndef _SlaterScaler_H_
#define _SlaterScaler_H_

#include "oml/vector.h"

class SlaterScaler
{
public:
    typedef Vector<double> RVec;
    SlaterScaler(size_t N, double emin, double emax, size_t LMax);
    size_t N   (size_t L) const 
    {
        size_t Nt=(itsN*(itsLMax+1))/(itsLMax+1);
        return Nt>0 ? Nt : 1;
    }
    double emin(size_t L) const {return itsemin/pow(1.0,L);}
    double emax(size_t L) const {return itsemax/pow(3,L);}    
    RVec   Get_es (size_t L) const;
private:        
    size_t itsN,itsLMax;
    double itsemin,itsemax;
    RVec es;
};
#endif // _SlaterScaler_H_
