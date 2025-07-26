// File: GaussianCD.H  Charge distribution for two primative gaussians.
#ifndef _GaussianCD_H_
#define _GaussianCD_H_


#include "PolarizedGaussian/MnD/Hermite2.H"

import Common.UniqueIDImp;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.GData;

namespace PolarizedGaussian
{

struct GaussianCD : public UniqueIDImp
{
    GaussianCD(const GData&,const GData&);
    GData GetGData() const {return GData{GetID(),AlphaP,P,Ltotal};};

//
//  Raw data defining the charge distribution, all pickeled.
//
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
        return theNMLs[LMax];
    }
    static void MakeNMLs();
    static std::vector<std::vector<Polarization>> theNMLs; //A list of all NMLs for each LMax.
};

} //namespace PolarizedGaussian

#endif //_GaussianCD_H_
