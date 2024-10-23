// File: Hermite2.C



#include "Imp/BasisSet/PolarizedGaussian/MnD/Hermite2.H"
#include "Imp/Misc/IntPower.H"
#include "oml/imp/binio.h"

#include <iostream>
#include <cassert>

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


Hermite2::Hermite2()
    : LA(0)
    , LB(0)
{}


inline bool near(double a, double b) 
{
    static double eps=1e-13;
    double err=fabs(a)<1 ? fabs(a-b) : fabs(a-b)/fabs(a) ;
    if (err>=eps)
        cout << std::setprecision(14) << a << " " << b << " " << err << endl;
    return err<eps;
}

//using std::cout;
//using std::endl;
Hermite2::Hermite2(double AlphaP, const RVec3& PA, const RVec3& PB, int _LA, int _LB)
    : LA(_LA)
    , LB(_LB)
    , d(GetSize())
    , e(GetSize())
    , f(GetSize())
    , d1(GetSize())
    , e1(GetSize())
    , f1(GetSize())
{
    assert(_LA>=0);
    assert(_LB>=0);
    //cout << "Constructor: " << LA << " " << LB << " " <<  d1.size()  << " " << AlphaP << endl;

    double a12=1.0/(2*AlphaP);

    size_t index=GetIndex(0,0,0);
    d[index]=e[index]=f[index]=1.0;
    d1[index]=e1[index]=f1[index]=1.0;
    for(size_t N=1;N<=LA+LB;N++)
        for (int na=0;na<=LA;na++)
        {
            int nb=N-na;
            if (nb>LB) continue;
            //cout << "old N,na,nb = " << N << " " << na << " " << nb << endl;
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
                    //cout << "old N,na,nb = " << N << " " << na << " " << nb << endl;
                    //cout << index << " " << Getd(N-1,na-1,nb) << endl;
                    size_t index1=GetIndex(N,na,nb);
                    d1[index1]=Getd1(N-1,na-1,nb)*a12 + Getd1(N,na-1,nb)*PA.x + Getd1(N+1,na-1,nb)*(N+1);
                    e1[index1]=Gete1(N-1,na-1,nb)*a12 + Gete1(N,na-1,nb)*PA.y + Gete1(N+1,na-1,nb)*(N+1);
                    f1[index1]=Getf1(N-1,na-1,nb)*a12 + Getf1(N,na-1,nb)*PA.z + Getf1(N+1,na-1,nb)*(N+1);
                    
                   
                }
                if (n<=LB)
                {
                    //cout << "old N,na,nb = " << N << " " << 0 << " " << n << endl;
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
                    //cout << "old N,na,nb = " << N << " " << na << " " << nb << endl;
                    size_t index1=GetIndex(N,na,nb);
                    d1[index1]=Getd1(N-1,na,nb-1)*a12 + Getd1(N,na,nb-1)*PB.x + Getd1(N+1,na,nb-1)*(N+1);
                    e1[index1]=Gete1(N-1,na,nb-1)*a12 + Gete1(N,na,nb-1)*PB.y + Gete1(N+1,na,nb-1)*(N+1);
                    f1[index1]=Getf1(N-1,na,nb-1)*a12 + Getf1(N,na,nb-1)*PB.z + Getf1(N+1,na,nb-1)*(N+1);
//                    cout << Getd (N-1,na-1,nb) << " " << Getd (N,na-1,nb) << " " << Getd (N+1,na-1,nb) << " " << endl;
//                    cout << Getd1(N-1,na-1,nb) << " " << Getd1(N,na-1,nb) << " " << Getd1(N+1,na-1,nb) << " " << endl << endl;
                  
                }
                if (n<=LA)
                {
                    //cout << "old N,na,nb = " << N << " " << n << " " << 0 << endl;
                    size_t indexb=GetIndex(N,n,0);
                    d1[indexb]=Getd1(N-1,n-1,0)*a12 + Getd1(N,n-1,0)*PA.x + Getd1(N+1,n-1,0)*(N+1);
                    e1[indexb]=Gete1(N-1,n-1,0)*a12 + Gete1(N,n-1,0)*PA.y + Gete1(N+1,n-1,0)*(N+1);
                    f1[indexb]=Getf1(N-1,n-1,0)*a12 + Getf1(N,n-1,0)*PA.z + Getf1(N+1,n-1,0)*(N+1);
                }

            }
        }
    }
//    if (LB>=2)
//        cout << "after  d1[0,0,1]=" << d1[GetIndex(0,0,1)] << " Getd(0,0,1)=" << Getd(0,0,1) << " Getd1(0,0,1)=" << Getd1(0,0,1) << endl;;


    //
    // Layered approach:  na+nb = layer.  
    //
    for (int layer=0;layer<=LA+LB-1;layer++)
    {
        for (int na=0;na<=layer;na++)
        {
            int nb=layer-na;
            if (na>LA || nb>LB) continue;
            for (int N=0;N<=layer+1;N++)
            {
                //cout << "new LA,LB,layer,N,na,nb = " << LA << " " << LB << " " << layer << " " << N << " " << na << " " << nb << endl;
                assert(near(Getd(N-1,na,nb),Getd1(N-1,na,nb)));
                assert(near(Getd(N  ,na,nb),Getd1(N  ,na,nb)));
                assert(near(Getd(N+1,na,nb),Getd1(N+1,na,nb)));
                assert(near(Gete(N-1,na,nb),Gete1(N-1,na,nb)));
                assert(near(Gete(N  ,na,nb),Gete1(N  ,na,nb)));
                assert(near(Gete(N+1,na,nb),Gete1(N+1,na,nb)));
                assert(near(Getf(N-1,na,nb),Getf1(N-1,na,nb)));
                assert(near(Getf(N  ,na,nb),Getf1(N  ,na,nb)));
                assert(near(Getf(N+1,na,nb),Getf1(N+1,na,nb)));

                assert(N<=na+nb+1);
                //if (na<LA) //is na+1 in range?
                double dt1=0.0,et1=0.0,ft1=0.0;
                double dta=0.0,eta=0.0,fta=0.0;
                double dtb=0.0,etb=0.0,ftb=0.0;
                double dt3=0.0,et3=0.0,ft3=0.0;
                if (N-1>=0)
                {
                    dt1=Getd(N-1,na,nb)*a12;
                    et1=Gete(N-1,na,nb)*a12;
                    ft1=Getf(N-1,na,nb)*a12;
                }
                if (N<=na+nb)
                {
                    dta=Getd(N,na,nb)*PA.x;
                    eta=Gete(N,na,nb)*PA.y;
                    fta=Getf(N,na,nb)*PA.z;
                    dtb=Getd(N,na,nb)*PB.x;
                    etb=Gete(N,na,nb)*PB.y;
                    ftb=Getf(N,na,nb)*PB.z;
                }
                if (N+1<=na+nb)
                {
                    dt3=(N+1)*Getd(N+1,na,nb);
                    et3=(N+1)*Gete(N+1,na,nb);
                    ft3=(N+1)*Getf(N+1,na,nb);
                }
                
                if (na+1<=LA)
                {
                    //cout << "new LA,LB,layer,N,na+1,nb = " << LA << " " << LB << " " << layer << " " << N << " " << na+1 << " " << nb << endl;
                    size_t indexa=GetIndex(N,na+1,nb);
                    d[indexa]=dt1+dta+dt3;
                    e[indexa]=et1+eta+et3;
                    f[indexa]=ft1+fta+ft3;                    
                    assert(near(Getd(N,na+1,nb),Getd1(N,na+1,nb)));
                    assert(near(Gete(N,na+1,nb),Gete1(N,na+1,nb)));
                    assert(near(Getf(N,na+1,nb),Getf1(N,na+1,nb)));
    //                    cout << Getd (N-1,na,nb) << " " << Getd (N,na,nb) << " " << Getd (N+1,na,nb) << " " << endl;
    //                    cout << Getd1(N-1,na,nb) << " " << Getd1(N,na,nb) << " " << Getd1(N+1,na,nb) << " " << endl ;
                    assert(near(d[indexa],d1[indexa]));                    
                    assert(near(e[indexa],e1[indexa]));                    
                    assert(near(f[indexa],f1[indexa]));                    
                }
                
                if (nb+1<=LB)
                {
                    //cout << "new LA,LB,layer,N,na,nb+1 = " << LA << " " << LB << " " << layer << " " << N << " " << na << " " << nb+1 << endl;
                    size_t indexb=GetIndex(N,na,nb+1);
                    d[indexb]=dt1+dtb+dt3;
                    e[indexb]=et1+etb+et3;
                    f[indexb]=ft1+ftb+ft3;
                    
                    assert(near(Getd(N,na,nb+1),Getd1(N,na,nb+1))); 
                    assert(near(Gete(N,na,nb+1),Gete1(N,na,nb+1)));
                    assert(near(Getf(N,na,nb+1),Getf1(N,na,nb+1)));

                    //cout << "Getd(N,na+1,nb)=" << Getd(N,na+1,nb) << " Getd(N,na,nb+1)=" << Getd(N,na,nb+1) << endl<< endl;
                    assert(near(d[indexb],d1[indexb]));
                    assert(near(e[indexb],e1[indexb]));
                    assert(near(f[indexb],f1[indexb]));
                }
                

            }
        }
    }
    //cout << d << endl << d1 << endl << endl;
}

double Hermite2::GetAny(const std::vector<double>& def, int N, int na, int nb) const
{
    if (nb>LB)
        cout << nb << endl;
    
    //assert(N <=LA+LB);
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
        retd=Getd(P.n,Pa.n,Pb.n);
    }
    if (rete==0.0)
    {
        rete=Gete(P.l,Pa.l,Pb.l);
    }
    if (retf==0.0)
    {
        retf=Getf(P.m,Pa.m,Pb.m);
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
