// File: GaussianCD.H  Charge distribution for two primative gaussians.
#ifndef _GaussianCD_H_
#define _GaussianCD_H_



#include "BasisSetImplementation/PolarizedGaussian/RadialFunction.H"
#include "BasisSetImplementation/PolarizedGaussian/Hermite/Hermite2.H"
#include "oml/matrix.h"
#include "Misc/ptr_vector1.h"
#include "oml/vector3d.h"


class RNLM;
class GaussianRF;

class GaussianCD
{
public:
    GaussianCD(const GaussianRF&,const GaussianRF&);


private:
    friend class GaussianRF;
//
//  Raw data defining the charge distribution, all pickeled.
//
    const RadialFunction&  r1;           //All this stuff gets pickled.
    const RadialFunction&  r2;           //All this stuff gets pickled.
    int                    Ltotal;       //Total angular momentum.
    double                 a,b;          //exponents.
    double                 ab,AlphaP;    //a*b, a+b.
    RVec3                  AB,P;         //A-B, new center.
    double                 Eij;          //scale factor.
    Hermite2               H2;           //Hermite coefficients.
//
//  Polarizations lists used for some summations.
//
    static const std::vector<Polarization>&  GetNMLs(int LMax)
    {
        return *theNMLs[LMax];
    }
    static void MakeNMLs();
    static optr_vector1<std::vector<Polarization>* > theNMLs; //A list of all NMLs for each LMax.
};

#endif //_GaussianCD_H_
