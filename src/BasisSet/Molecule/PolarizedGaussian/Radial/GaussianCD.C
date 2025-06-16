// File: GaussianCD.C  Charge distribution for two  gaussians.


#include "Molecule/PolarizedGaussian/Radial/GaussianCD.H"
#include "Molecule/PolarizedGaussian/Radial/GaussianRF.H"

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

GaussianCD::GaussianCD(const GaussianRF& the_g1,const GaussianRF& the_g2)
    : r1    (the_g1)
    , r2    (the_g2)
    , Ltotal(r1.GetL() + r2.GetL())
    , a     (the_g1.itsExponent)
    , b     (the_g2.itsExponent)
    , ab    (a * b)
    , AlphaP(a + b)
    , AB    (r1.GetCenter() - r2.GetCenter())
    , P     ( (a*r1.GetCenter() + b*r2.GetCenter()) / AlphaP)
    , Eij   ( exp(-ab / AlphaP * (AB*AB)) )
    , H2    (AlphaP, P - r1.GetCenter(), P - r2.GetCenter(), r1.GetL()+1, r2.GetL()+1)
{
    if (theNMLs.size()==0) MakeNMLs();
};

GaussianCD::~GaussianCD()
{
    
}


std::ostream& GaussianCD::Write(std::ostream& os) const {return os;}
std::istream& GaussianCD::Read (std::istream& is)       {return is;}

} //namespace PolarizedGaussian
