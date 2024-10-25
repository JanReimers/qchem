// File: SlaterScaler.H  Rescale Slater exponents based in angular momentum L.
#ifndef _SlaterScaler_H_
#define _SlaterScaler_H_

class SlaterScaler
{
public:
    SlaterScaler(size_t N, double emin, double emax, size_t LMax) :
        itsN(N), itsLMax(LMax), itsemin(emin), itsemax(emax) {};
    size_t N   (size_t L) const {return (itsN*(itsLMax+1-L))/(itsLMax+1);}
    double emin(size_t L) const {return itsemin/pow(.7,L);}
    double emax(size_t L) const {return itsemax/pow(3,L);}    
private:        
    size_t itsN,itsLMax;
    double itsemin,itsemax;
};
#endif // _SlaterScaler_H_
