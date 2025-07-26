// File: GaussianCD.C  Charge distribution for two  gaussians.

#include <cmath>
#include "PolarizedGaussian/Radial/GaussianCD.H"

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
GaussianCD::GaussianCD(const GData& g1,const GData& g2)
    : Ltotal(g1.L + g2.L)
    , a     (g1.Alpha)
    , b     (g2.Alpha)
    , ab    (a * b)
    , AlphaP(a + b)
    , AB    (g1.R - g2.R)
    , P     ( (a*g1.R + b*g2.R) / AlphaP)
    , Eij   ( exp(-ab / AlphaP * (AB*AB)) )
    , H2    (AlphaP, P - g1.R, P - g2.R, g1.L+1, g2.L+1)
{
    if (theNMLs.size()==0) MakeNMLs();
};



} //namespace PolarizedGaussian
