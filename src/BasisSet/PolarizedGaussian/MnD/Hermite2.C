// File: Hermite2.C



#include "Imp/BasisSet/PolarizedGaussian/MnD/Hermite2.H"
#include "Imp/Misc/IntPower.H"
#include "oml/imp/binio.h"
#include <iostream>
#include <cassert>

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


Hermite2::Hermite2()
    : LA(0)
    , LB(0)
{}


inline bool near(double a, double b) {return fabs(a-b)<1e-14;}

//using std::cout;
//using std::endl;
Hermite2::Hermite2(double AlphaP, const RVec3& PA, const RVec3& PB, int _LA, int _LB)
    : LA(_LA)
    , LB(_LB)
    , d1(GetSize())
    , e1(GetSize())
    , f1(GetSize())
{
    assert(_LA>=0);
    assert(_LB>=0);
    //cout << "Constructor: " << LA << " " << LB << " " <<  d1.size()  << " " << AlphaP << endl;
    assert(LA<LMAX+1);
    assert(LB<LMAX+1);

    double a12=1.0/(2*AlphaP);

    size_t index=GetIndex(0,0,0);
    d1[index]=e1[index]=f1[index]=1.0;
    for(size_t N=1;N<=LA+LB;N++)
        for (int na=0;na<=LA;na++)
        {
            int nb=N-na;
            if (nb>LB) continue;
            //cout << "N,na,nb = " << N << " " << na << " " << nb << endl;
            size_t index=GetIndex(N,na,nb);
            //cout << index << " " << d1.size() << endl;
            d1[index]=e1[index]=f1[index]=intpow(a12,N);        
        }

    for (int n=1; n<=LA+LB; n++)
    {
        for (int N=0; N<n ; N++)
        {
            if (LA>=LB)
            {
                int nb_low=n-LA;
                if (nb_low<0) nb_low=0;
                for (int nb=nb_low; nb<n && nb<=LB; nb++)
                {
                    int na=n-nb;
                    //cout << "N,na,nb = " << N << " " << na << " " << nb << endl;
                    //cout << index << " " << Getd(N-1,na-1,nb) << endl;
                    size_t index1=GetIndex(N,na,nb);
                    d1[index1]=Getd1(N-1,na-1,nb)*a12 + Getd1(N,na-1,nb)*PA.x + Getd1(N+1,na-1,nb)*(N+1);
                    e1[index1]=Gete1(N-1,na-1,nb)*a12 + Gete1(N,na-1,nb)*PA.y + Gete1(N+1,na-1,nb)*(N+1);
                    f1[index1]=Getf1(N-1,na-1,nb)*a12 + Getf1(N,na-1,nb)*PA.z + Getf1(N+1,na-1,nb)*(N+1);
                    
                   
                }
                if (n<=LB)
                {
                    //cout << "N,na,nb = " << N << " " << 0 << " " << n << endl;
                    size_t index1=GetIndex(N,0,n);
                    d1[index1]=Getd1(N-1,0,n-1)*a12 + Getd1(N,0,n-1)*PB.x + Getd1(N+1,0,n-1)*(N+1);
                    e1[index1]=Gete1(N-1,0,n-1)*a12 + Gete1(N,0,n-1)*PB.y + Gete1(N+1,0,n-1)*(N+1);
                    f1[index1]=Getf1(N-1,0,n-1)*a12 + Getf1(N,0,n-1)*PB.z + Getf1(N+1,0,n-1)*(N+1);
                }
            }
            else
            {
                int na_low=n-LB;
                if (na_low<0) na_low=0;
                for (int na=na_low; na<n && na<=LA; na++)
                {
                    int nb=n-na;
                    //cout << "N,na,nb = " << N << " " << na << " " << nb << endl;
                    size_t index1=GetIndex(N,na,nb);
                    d1[index1]=Getd1(N-1,na,nb-1)*a12 + Getd1(N,na,nb-1)*PB.x + Getd1(N+1,na,nb-1)*(N+1);
                    e1[index1]=Gete1(N-1,na,nb-1)*a12 + Gete1(N,na,nb-1)*PB.y + Gete1(N+1,na,nb-1)*(N+1);
                    f1[index1]=Getf1(N-1,na,nb-1)*a12 + Getf1(N,na,nb-1)*PB.z + Getf1(N+1,na,nb-1)*(N+1);
//                    cout << Getd (N-1,na-1,nb) << " " << Getd (N,na-1,nb) << " " << Getd (N+1,na-1,nb) << " " << endl;
//                    cout << Getd1(N-1,na-1,nb) << " " << Getd1(N,na-1,nb) << " " << Getd1(N+1,na-1,nb) << " " << endl << endl;
                  
                }
                if (n<=LA)
                {
                    size_t index1=GetIndex(N,n,0);
                    d1[index1]=Getd1(N-1,n-1,0)*a12 + Getd1(N,n-1,0)*PA.x + Getd1(N+1,n-1,0)*(N+1);
                    e1[index1]=Gete1(N-1,n-1,0)*a12 + Gete1(N,n-1,0)*PA.y + Gete1(N+1,n-1,0)*(N+1);
                    f1[index1]=Getf1(N-1,n-1,0)*a12 + Getf1(N,n-1,0)*PA.z + Getf1(N+1,n-1,0)*(N+1);
                }

            }
        }
    }
//    if (LB>=2)
//        cout << "after  d1[0,0,1]=" << d1[GetIndex(0,0,1)] << " Getd(0,0,1)=" << Getd(0,0,1) << " Getd1(0,0,1)=" << Getd1(0,0,1) << endl;;

}

double Hermite2::GetAny(const std::vector<double>& def, int N, int na, int nb) const
{
    assert(N <=LA+LB);
    assert(na<=LA      );
    assert(nb<=LB      );

//
//  Take care of negative n's generated by derivatives of distributions.
//
    if ( na < 0 || nb < 0 || N < 0) return 0.0;

    int n=na+nb; // we now know that n1 and n2 are >=0.
//
//  Take care of s-s charge distributions.  We know n1 and n2 are >=0.
//
    if (n==0 && N==0) return 1.0;
//
//  All other 0.0 elements.
//
    if (N >  n) return 0.0;
//
//  Use std indexing to get at the numbers
//
    int index = GetIndex(N,na,nb);
    assert(index>=0        );
    return def[index];
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

    int n=Pa.n+Pb.n; // we now know that na and nb are >=0.
    int l=Pa.l+Pb.l; // we now know that na and nb are >=0.
    int m=Pa.m+Pb.m; // we now know that na and nb are >=0.
//
//  All other 0.0 elements.
//
    if (P.n>n || P.l>l || P.m>m ) return 0.0;
//
//  Take care of s-s charge distributions.  We know n1 and n2 are >=0.
//
    double retd=0.0, rete=0.0, retf=0.0;
    if (n==0 && P.n==0) retd=1.0;
    if (l==0 && P.l==0) rete=1.0;
    if (m==0 && P.m==0) retf=1.0;
//
//  Use fancy indexing to get at the numbers
//
    if (retd==0.0)
    {
        retd=Getd1(P.n,Pa.n,Pb.n);
    }
    if (rete==0.0)
    {
        rete=Gete1(P.l,Pa.l,Pb.l);
    }
    if (retf==0.0)
    {
        retf=Getf1(P.m,Pa.m,Pb.m);
    }
    return retd * rete * retf;
}

std::ostream& Hermite2::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
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
                for (int nb=0; nb<=LB; nb++)  os << Getd1(N,na,nb) << " ";
                os << std::endl;
            }
        }
        for (int N=0; N<=LA+LB; N++)
        {
            os << "---------------------------------------------------------------------" << std::endl;
            os << "L = " << N << std::endl;
            for (int na=0; na<=LA; na++)
            {
                for (int nb=0; nb<=LB; nb++)  os << Gete1(N,na,nb) << " ";
                os << std::endl;
            }
        }
        for (int N=0; N<=LA+LB; N++)
        {
            os << "---------------------------------------------------------------------" << std::endl;
            os << "M = " << N << std::endl;
            for (int na=0; na<=LA; na++)
            {
                for (int nb=0; nb<=LB; nb++)  os << Getf1(N,na,nb) << " ";
                os << std::endl;
            }
        }
    }
    return os;
}

std::istream& Hermite2::Read (std::istream& is)
{
    return is;
}

Hermite2* Hermite2::Clone() const
{
    return new Hermite2(*this);
}


} //namespace PolarizedGaussian
