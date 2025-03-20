// File: GaussianCD.H  Charge distribution for two primative gaussians.
#ifndef _GaussianCD_H_
#define _GaussianCD_H_


#include "Imp/BasisSet/PolarizedGaussian/MnD/Hermite2.H"
#include "types.H"
#include "Imp/Containers/ptr_vector.h"
#include "Imp/Misc/UniqueID/UniqueIDImp.H"
#include "oml/matrix.h"
#include "oml/vector.h"
#include "oml/vector3d.h"
#include <map>

namespace PolarizedGaussian
{

class RNLM;
class GaussianRF;
class RadialFunction;
class Hermite1;

class GaussianCD : public UniqueIDImp
{
public:
    GaussianCD(const GaussianRF&,const GaussianRF&);
    ~GaussianCD();
    
    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;


private:
    friend class GaussianRF;
    friend class CDCache;
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

} //namespace PolarizedGaussian

#endif //_GaussianCD_H_
