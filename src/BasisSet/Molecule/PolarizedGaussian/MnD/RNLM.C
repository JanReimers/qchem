// File: RNLM.C  Manager the RNLM Auxillary functions.



#include "Imp/BasisSet/Molecule/PolarizedGaussian/MnD/RNLM.H"
#include "Imp/Integrals/AuxillaryFJ.H"
#include "Base/IntPower.H"
#include "oml/vector.h"
#include <iostream>
#include <cassert>

#define MAX 16

namespace PolarizedGaussian
{


RNLM::RNLM() : itsData() {};

RNLM::RNLM(int Max, double Alpha, const RVec3& dR)
    : itsLMax(Max)
    , itsData()
{
    static Triangle3D theRjs[MAX+1];
    static bool     theRjsInitialized = false;
    if(!theRjsInitialized)
    {
        assert(MAX==AuxillaryFJ::thejMax);
        int GlobalMax=AuxillaryFJ::thejMax;
        for (int j=0; j<=GlobalMax; j++) theRjs[j] = Triangle3D(GlobalMax-j);
        theRjsInitialized=true;
    }

    assert(Max<=AuxillaryFJ::thejMax);
    double T=Alpha*dR*dR;
    Vector<double> Fj(0,Max,0.0);
    {
        AuxillaryFJ FjCalculator;
        FjCalculator.GetFjAt(T,Fj);
    }

//  cout << "T=" << T << " Max = " << Max << "  Fj=" << Fj << std::endl;
    Triangle3D* Rjs=&theRjs[AuxillaryFJ::thejMax-Max];
//  for (int j=0;j<=Max;j++) Rjs.push_back(new Triangle(Max-j));

    Rjs[Max](0,0,0)=uintpow(-2.0*Alpha,Max)*Fj(Max);
    for (int j=Max-1; j>=0; j--)
    {
        Triangle3D& Rj  =Rjs[j  ];
        const Triangle3D& Rjp1=Rjs[j+1];

        Rj(0,0,0)=intpow(-2.0*Alpha,j)*Fj(j);
        int n=Max-j;
        for (int M=0; M<n; M++)
        {
            Rj(0,0,M+1)= M>0 ? dR.z*Rjp1(0,0,M)+M*Rjp1(0,0,M-1) : dR.z*Rjp1(0,0,M);
            for (int L=0; L<n-M; L++)
            {
                Rj(0,L+1,M)= L>0 ? dR.y*Rjp1(0,L,M)+L*Rjp1(0,L-1,M) : dR.y*Rjp1(0,L,M);
                for (int N=0; N<n-M-L; N++)
                {
                    Rj(N+1,L,M)= N>0 ? dR.x*Rjp1(N,L,M)+N*Rjp1(N-1,L,M) : dR.x*Rjp1(N,L,M);
                }
            }
        }
    }

    itsData = Rjs[0];
}

void RNLM::Add(const RNLM& theR, double theScale)
{
    itsData.Add(theR.itsData,theScale);
}


std::ostream& RNLM::Write(std::ostream& os) const
{
    return os << itsData;
}

std::istream& RNLM::Read (std::istream& is)
{
    return is >> itsData;
}

RNLM* RNLM::Clone() const
{
    return new RNLM(*this);
}

} //namespace PolarizedGaussian
