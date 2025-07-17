// File: Hermite2.C


#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>

#include "PolarizedGaussian/MnD/Hermite2.H"
#include "PolarizedGaussian/Polarization.H"

import Common.IntPow;

std::ostream& operator << (std::ostream& os, const std::vector<double>& v)
{
    for (auto i:v) os << i << " ";
    return os;
}

namespace PolarizedGaussian
{
//----------------------------------------------------------------------------------------
//
//  None of this will make any sense unless you read McMurchie and Davidson,
//  J. Comp. Phys. 26, 218 (1978).
//

/*
 The storage here is non-trivial:
 I will label the elements as d(N,n1,n2) with 0 <= n1+n2 = n < LMax, and 0 < N <= n.
 d(0,n1,n2) looks like this up to n=2:

    1         B       B**2+1/2a .
    A         AB+1/2a .         .
    A**2+1/2a .       .         .
    .         .       .         .


 d(1,n1,n2) looks like

    0     1/2a     2B/2a .
    1/2a  (A+B)/2a .     .
    2A/2a .        .     .
    .     .        .     .

 d(2,n1,n2) looks like

    0         0         1/(2a)**2 .
    0         1/(2a)**2 .         .
    1/(2a)**2 .         .         .
    .         .         .         .

 You can see the all the 0 and 1/(2a)**N elements need not be stored.   And all the zeros
also need not be stored.  The elements are stored using the fancy lookup table below to transform
(N,na,nb) ---> array index.  ANother lookup table is used to get the required size for a given
LA,LB pair.
   The N arrays are filled in the following order:  As you move from top left to bottom right
 n=n1+n2 increases.  Start by loading the n=1 diagonal for n2=0 to 1, in the N=0 array.  N>0
 does not need this row filled.  The move on to n=2, and do n2=0 to 2 for N=0 and 1.  Then n=3
 for N=1 to n-1, and n2=0 to n+1.
*/

using std::cout;
using std::endl;
//----------------------------------------------------------------------------------------




inline bool near(double a, double b) 
{
    static double eps=1e-13;
    double err=fabs(a)<1 ? fabs(a-b) : fabs(a-b)/fabs(a) ;
    if (err>=eps)
        cout << std::setprecision(14) << a << " " << b << " " << err << endl;
    return err<eps;
}

inline bool near(const RVec3& a, const RVec3& b) 
{
    return near(a.x,b.x) && near(a.y,b.y) &&  near(a.z,b.z);
}

inline RVec3 DirectMultiply(const RVec3& a, const RVec3& b)
{
    return RVec3(a.x*b.x, a.y*b.y, a.z*b.z);
}

inline RVec3 operator*(const RVec3& a, const int& b)
{
    return RVec3(a.x*b,a.y*b,a.z*b);
}

//using std::cout;
//using std::endl;
Hermite2::Hermite2(double AlphaP, const RVec3& PA, const RVec3& PB, int _LA, int _LB)
    : LA(_LA)
    , LB(_LB)
    , LAB((LA+1)*(LB+1))
    , LB1(LB+1)    
    , d(GetSize())
    , e(GetSize())
    , f(GetSize())
{
    assert(_LA>=0);
    assert(_LB>=0);
   // cout << "Hermite 2 constructor" << endl;
    
    double a12=1.0/(2*AlphaP);
    Assign(0,0,0,RVec3(1,1,1)); //d_0^00 
    //
    // Layered approach:  na+nb = layer.  
    //
    for (int layer=0;layer<=LA+LB-1;layer++)
        for (int na=0;na<=layer;na++)
        {
            int nb=layer-na;
            if (na>LA || nb>LB) continue;
            for (int N=0;N<=layer+1;N++)
            {
                assert(N<=na+nb+1);
                RVec3 t1(0,0,0),ta(0,0,0),tb(0,0,0),t3(0,0,0);
                if (N-1>=0)
                    t1=Get(N-1,na,nb)*a12;
                if (N<=na+nb)
                {
                    ta=DirectMultiply(Get(N,na,nb),PA);
                    tb=DirectMultiply(Get(N,na,nb),PB);
                }
                if (N+1<=na+nb)
                    t3=Get(N+1,na,nb)*(N+1);
                
                if (na+1<=LA)
                    Assign(N,na+1,nb,t1+ta+t3);                 
                
                if (nb+1<=LB)
                    Assign(N,na,nb+1,t1+tb+t3);    
            }
        }
}

Hermite2::~Hermite2()
{
    //cout << (void*)this << " destructor! " << endl;
}
void Hermite2::Assign(int N,int na,int nb,const RVec3& a)
{
    size_t index=GetIndex(N,na,nb);
    #ifdef USE_CACHE
    Polarization key(N,na,nb);
    if (auto i=indexCache.find(key);i==indexCache.end())
        indexCache[key]=index;
    #endif
    d[index]=a.x;
    e[index]=a.y;
    f[index]=a.z;
}

double Hermite2::operator()(const Polarization& P,const Polarization& Pa,const Polarization& Pb) const
{
#if DEBUG
    if (P.GetTotalL() > LA+LB)
    {
        std::cerr << "Hermite 2 overflow for P (P,Pa,Pb)=(" << P << "," << Pa << "," << Pb << ")  (La,Lb)=(" << LA << "," << LB << ")" << std::endl;
    }
    if (Pa.GetTotalL() >LA)
    {
        std::cerr << "Hermite 2 overflow for Pa (P,Pa,Pb)=(" << P << "," << Pa << "," << Pb << ")  (La,Lb)=(" << LA << "," << LB << ")" << std::endl;
    }
    if (Pb.GetTotalL() >LB)
    {
        std::cerr << "Hermite 2 overflow for Pb (P,Pa,Pb)=(" << P << "," << Pa << "," << Pb << ")  (La,Lb)=(" << LA << "," << LB << ")" << std::endl;
    }
#endif

    assert(P .GetTotalL()<=LA+LB);
    assert(Pa.GetTotalL()<=LA      );
    assert(Pb.GetTotalL()<=LB      );
//
//  Take care of negative n's generated by derivatives of distributions.
//
    if ( Pa.n < 0 || Pb.n < 0 || P.n < 0) return 0.0;
    if ( Pa.l < 0 || Pb.l < 0 || P.l < 0) return 0.0;
    if ( Pa.m < 0 || Pb.m < 0 || P.m < 0) return 0.0;
    
    #ifdef USE_CACHE
    auto dindex=indexCache.find(Polarization(P.n,Pa.n,Pb.n));
    auto eindex=indexCache.find(Polarization(P.l,Pa.l,Pb.l));
    auto findex=indexCache.find(Polarization(P.m,Pa.m,Pb.m));
//    cout << (void*)this << "  Looking up index at " << Polarization(P.n,Pa.n,Pb.n) << endl;
//    if (!(dindex!=indexCache.end()))
//        for (auto i:indexCache) cout << i.first << " ";
    assert(dindex!=indexCache.end());
    assert(eindex!=indexCache.end());
    assert(findex!=indexCache.end());
    
    return d[dindex->second]*e[eindex->second]*f[findex->second];
    #endif
    return d[GetIndex(P.n,Pa.n,Pb.n)]*e[GetIndex(P.l,Pa.l,Pb.l)]*f[GetIndex(P.m,Pa.m,Pb.m)];
}

std::ostream& Hermite2::Write(std::ostream& os) const
{
    os << ": LA=" << LA << ", LB = " << LB << std::endl;
    os.precision(4);
    os.width(8);
    os.setf(std::ios::fixed,std::ios::floatfield);
    for (int N=0; N<=LA+LB; N++)
    {
        os << "---------------------------------------------------------------------" << std::endl;
        os << "N = " << N << std::endl;
        for (int na=0; na<=LA; na++)
        {
            for (int nb=0; nb<=LB; nb++)  os << d[GetIndex(N,na,nb)] << " ";
            os << std::endl;
        }
    }
    for (int N=0; N<=LA+LB; N++)
    {
        os << "---------------------------------------------------------------------" << std::endl;
        os << "L = " << N << std::endl;
        for (int na=0; na<=LA; na++)
        {
            for (int nb=0; nb<=LB; nb++)  os << e[GetIndex(N,na,nb)] << " ";
            os << std::endl;
        }
    }
    for (int N=0; N<=LA+LB; N++)
    {
        os << "---------------------------------------------------------------------" << std::endl;
        os << "M = " << N << std::endl;
        for (int na=0; na<=LA; na++)
        {
            for (int nb=0; nb<=LB; nb++)  os << f[GetIndex(N,na,nb)] << " ";
            os << std::endl;
        }
    }
    return os;
}


Hermite2* Hermite2::Clone() const
{
    return new Hermite2(*this);
}


} //namespace PolarizedGaussian
