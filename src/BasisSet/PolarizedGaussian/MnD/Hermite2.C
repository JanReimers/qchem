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

#define Z -2*LMAX-1
//
//  Indexing lookup table for LA>LB
//
int Hermite2::theIndex[2*LMAX+1][LMAX+1][LMAX+1] =
{
    {
        //N=0
        {0 ,2 ,10,31},
        {1 ,3 ,11,32},
        {5 ,7 ,12,33},
        {19,22,26,34}
    },
    {
        //N=1
        {Z ,-1,13,35},
        {-1,4 ,14,36},
        {6 ,8 ,15,37},
        {20,23,27,38}
    },
    {
        //N=2
        {Z ,Z ,-2,39},
        {Z ,-2,16,40},
        {-2, 9,17,41},
        {21,24,28,42}
    },
    {
        //N=3
        {Z ,Z ,Z ,-3},
        {Z ,Z ,-3,43},
        {Z ,-3,18,44},
        {-3,25,29,45}
    },
    {
        //N=4
        {Z ,Z ,Z ,Z },
        {Z ,Z ,Z ,-4},
        {Z ,Z ,-4,46},
        {Z ,-4,30,47}
    },
    {
        //N=5
        {Z ,Z ,Z ,Z },
        {Z ,Z ,Z ,Z },
        {Z ,Z ,Z ,-5},
        {Z ,Z ,-5,48}
    },
    {
        //N=6
        {Z ,Z ,Z ,Z },
        {Z ,Z ,Z ,Z },
        {Z ,Z ,Z ,Z },
        {Z ,Z ,Z ,-6}
    }
};

int Hermite2::theArraySizes[LMAX+1][LMAX+1] =
{
    {0 ,1 ,6 ,21},
    {1 ,4 ,9 ,25},
    {6 ,9 ,18,30},
    {21,25,30,48}
};
#undef Z

//----------------------------------------------------------------------------------------


Hermite2::Hermite2()
    : itsLA(0)
    , itsLB(0)
{}


//using std::cout;
//using std::endl;
Hermite2::Hermite2(double AlphaP, const RVec3& PA, const RVec3& PB, int LA, int LB)
    : itsLA(LA)
    , itsLB(LB)
    , d(-(itsLA+itsLB)+1,theArraySizes[itsLA][itsLB],0.0)
    , e(-(itsLA+itsLB)+1,theArraySizes[itsLA][itsLB],0.0)
    , f(-(itsLA+itsLB)+1,theArraySizes[itsLA][itsLB],0.0)
{
    //cout << LA << " " << LB << " " <<  d.GetLimits() << " " << e.GetLimits() << " " <<  f.GetLimits() << endl;
    assert(LA<LMAX+1);
    assert(LB<LMAX+1);
    Vector<double>::Subscriptor sd(d);
    Vector<double>::Subscriptor se(e);
    Vector<double>::Subscriptor sf(f);

    double a12=1.0/(2*AlphaP);
    for (int N=1; N<=itsLA+itsLB ; N++)
    {
        sd(-N+1)=se(-N+1)=sf(-N+1)=intpow(a12,N);
    }


    for (int n=1; n<=itsLA+itsLB; n++)
    {
        for (int N=0; N<n ; N++)
        {
            if (itsLA>=itsLB)
            {
                int nb_low=n-itsLA;
                if (nb_low<0) nb_low=0;
                for (int nb=nb_low; nb<n && nb<=itsLB; nb++)
                {
                    int na=n-nb;
                    int index=theIndex[N][na][nb];
                    sd(index) = Getd(N-1,na-1,nb)*a12 + Getd(N,na-1,nb)*PA.x + Getd(N+1,na-1,nb)*(N+1);
                    se(index) = Gete(N-1,na-1,nb)*a12 + Gete(N,na-1,nb)*PA.y + Gete(N+1,na-1,nb)*(N+1);
                    sf(index) = Getf(N-1,na-1,nb)*a12 + Getf(N,na-1,nb)*PA.z + Getf(N+1,na-1,nb)*(N+1);
                }
                if (n<=itsLB)
                {
                    int index=theIndex[N][0][n];
                    sd(index) = Getd(N-1,0,n-1)*a12 + Getd(N,0,n-1)*PB.x + Getd(N+1,0,n-1)*(N+1);
                    se(index) = Gete(N-1,0,n-1)*a12 + Gete(N,0,n-1)*PB.y + Gete(N+1,0,n-1)*(N+1);
                    sf(index) = Getf(N-1,0,n-1)*a12 + Getf(N,0,n-1)*PB.z + Getf(N+1,0,n-1)*(N+1);
                }

            }
            else
            {
                int na_low=n-itsLB;
                if (na_low<0) na_low=0;
                for (int na=na_low; na<n && na<=itsLA; na++)
                {
                    int nb=n-na;
                    int index=theIndex[N][nb][na];
                    sd(index) = Getd(N-1,na,nb-1)*a12 + Getd(N,na,nb-1)*PB.x + Getd(N+1,na,nb-1)*(N+1);
                    se(index) = Gete(N-1,na,nb-1)*a12 + Gete(N,na,nb-1)*PB.y + Gete(N+1,na,nb-1)*(N+1);
                    sf(index) = Getf(N-1,na,nb-1)*a12 + Getf(N,na,nb-1)*PB.z + Getf(N+1,na,nb-1)*(N+1);
                }
                if (n<=itsLA)
                {
                    int index=theIndex[N][0][n];
                    sd(index) = Getd(N-1,n-1,0)*a12 + Getd(N,n-1,0)*PA.x + Getd(N+1,n-1,0)*(N+1);
                    se(index) = Gete(N-1,n-1,0)*a12 + Gete(N,n-1,0)*PA.y + Gete(N+1,n-1,0)*(N+1);
                    sf(index) = Getf(N-1,n-1,0)*a12 + Getf(N,n-1,0)*PA.z + Getf(N+1,n-1,0)*(N+1);
                }

            }
        }
    }
}


double Hermite2::GetAny(const Vector<double>& def, int N, int na, int nb) const
{
    assert(N <=itsLA+itsLB);
    assert(na<=itsLA      );
    assert(nb<=itsLB      );

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
//  Use fancy indexing to get at the numbers
//
    int index = (itsLA>=itsLB) ? theIndex[N][na][nb] : theIndex[N][nb][na];
    assert(index!=0        );
    assert(index!=-2*LMAX-1);
    if (index <  0      )  return def(index+1);
    else return def(index);
}


double Hermite2::operator()(const Polarization& P,const Polarization& Pa,const Polarization& Pb) const
{
#if DEBUG
    if (P.GetTotalL() > itsLA+itsLB)
    {
        std::cerr << "Hermite 2 overflow for P (P,Pa,Pb)=(" << P << "," << Pa << "," << Pb << ")  (La,Lb)=(" << itsLA << "," << itsLB << ")" << std::endl;
    }
    if (Pa.GetTotalL() >itsLA)
    {
        std::cerr << "Hermite 2 overflow for Pa (P,Pa,Pb)=(" << P << "," << Pa << "," << Pb << ")  (La,Lb)=(" << itsLA << "," << itsLB << ")" << std::endl;
    }
    if (Pb.GetTotalL() >itsLB)
    {
        std::cerr << "Hermite 2 overflow for Pb (P,Pa,Pb)=(" << P << "," << Pa << "," << Pb << ")  (La,Lb)=(" << itsLA << "," << itsLB << ")" << std::endl;
    }
#endif

    assert(P .GetTotalL()<=itsLA+itsLB);
    assert(Pa.GetTotalL()<=itsLA      );
    assert(Pb.GetTotalL()<=itsLB      );
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
        int index = (itsLA>=itsLB) ? theIndex[P.n][Pa.n][Pb.n] : theIndex[P.n][Pb.n][Pa.n];
        assert(index!=0        );
        assert(index!=-2*LMAX-1);
        retd = index<0 ? d(index+1) : d(index);
    }
    if (rete==0.0)
    {
        int index = (itsLA>=itsLB) ? theIndex[P.l][Pa.l][Pb.l] : theIndex[P.l][Pb.l][Pa.l];
        assert(index!=0        );
        assert(index!=-2*LMAX-1);
        rete = index<0 ? e(index+1) : e(index);
    }
    if (retf==0.0)
    {
        int index = (itsLA>=itsLB) ? theIndex[P.m][Pa.m][Pb.m] : theIndex[P.m][Pb.m][Pa.m];
        assert(index!=0        );
        assert(index!=-2*LMAX-1);
        retf = index<0 ? f(index+1) : f(index);
    }
    return retd * rete * retf;
}

std::ostream& Hermite2::Write(std::ostream& os) const
{
    if (StreamableObject::Binary())
    {
        BinaryWrite(itsLA,os);
        BinaryWrite(itsLB,os);
        os << d << e << f;
    }
    if (StreamableObject::Ascii())
    {
        os << itsLA << " " << itsLB << " " << d << " " << e << " " << f;
    }
    if (StreamableObject::Pretty())
    {
        os << ": LA=" << itsLA << ", LB = " << itsLB << std::endl;
        os << "d Data = " << d << "e Data = " << e << "f Data = " << f << std::endl;
        os.precision(4);
        os.width(8);
        os.setf(std::ios::fixed,std::ios::floatfield);
        for (int N=0; N<=itsLA+itsLB; N++)
        {
            os << "---------------------------------------------------------------------" << std::endl;
            os << "N = " << N << std::endl;
            for (int na=0; na<=itsLA; na++)
            {
                for (int nb=0; nb<=itsLB; nb++)  os << Getd(N,na,nb) << " ";
                os << std::endl;
            }
        }
        for (int N=0; N<=itsLA+itsLB; N++)
        {
            os << "---------------------------------------------------------------------" << std::endl;
            os << "L = " << N << std::endl;
            for (int na=0; na<=itsLA; na++)
            {
                for (int nb=0; nb<=itsLB; nb++)  os << Gete(N,na,nb) << " ";
                os << std::endl;
            }
        }
        for (int N=0; N<=itsLA+itsLB; N++)
        {
            os << "---------------------------------------------------------------------" << std::endl;
            os << "M = " << N << std::endl;
            for (int na=0; na<=itsLA; na++)
            {
                for (int nb=0; nb<=itsLB; nb++)  os << Getf(N,na,nb) << " ";
                os << std::endl;
            }
        }
    }
    return os;
}

std::istream& Hermite2::Read (std::istream& is)
{
    if (StreamableObject::Binary())
    {
        BinaryRead(itsLA,is);
        BinaryRead(itsLB,is);
        is >> d >> e >> f;
    }
    if (StreamableObject::Ascii())
    {
        is >> itsLA >> itsLB >> d >> e >> f;
    }
    return is;
}

Hermite2* Hermite2::Clone() const
{
    return new Hermite2(*this);
}


} //namespace PolarizedGaussian
