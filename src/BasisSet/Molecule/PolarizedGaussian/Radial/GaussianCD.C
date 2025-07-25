// File: GaussianCD.C  Charge distribution for two  gaussians.

#include <cmath>
#include "PolarizedGaussian/Radial/GaussianCD.H"
#include "PolarizedGaussian/Radial/GaussianRF.H"

namespace PolarizedGaussian
{

std::vector<std::vector<Polarization>> GaussianCD::theNMLs;

std::vector<Polarization> MakeAllPolarizations(int Lmax)
{
    std::vector<Polarization> list;
    for (int n=0; n<=Lmax; n++)
        for (int l=0; l<=Lmax-n; l++)
            for (int m=0; m<=Lmax-n-l; m++)
                list.push_back(Polarization(n,l,m));
    return list;
}

void GaussianCD::MakeNMLs()
{
    for (int L=0; L<=10; L++) theNMLs.push_back(MakeAllPolarizations(L));
}

//--------------------------------------------------------------------------------------
//
//  Construction zone
//

GaussianCD::GaussianCD(const GaussianRF& g1,const GaussianRF& g2)
    // : r1    (the_g1)
    // , r2    (the_g2)
    : Ltotal(g1.GetL() + g2.GetL())
    , a     (g1.itsExponent)
    , b     (g2.itsExponent)
    , ab    (a * b)
    , AlphaP(a + b)
    , AB    (g1.GetCenter() - g2.GetCenter())
    , P     ( (a*g1.GetCenter() + b*g2.GetCenter()) / AlphaP)
    , Eij   ( exp(-ab / AlphaP * (AB*AB)) )
    , H2    (AlphaP, P - g1.GetCenter(), P - g2.GetCenter(), g1.GetL()+1, g2.GetL()+1)
{
    if (theNMLs.size()==0) MakeNMLs();
};

GaussianCD::~GaussianCD()
{
    
}


std::ostream& GaussianCD::Write(std::ostream& os) const {return os;}

} //namespace PolarizedGaussian
