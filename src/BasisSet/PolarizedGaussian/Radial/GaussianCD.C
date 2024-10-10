// File: GaussianCD.C  Charge distribution for two  gaussians.


#include "oml/smatrix.h"
#include "Imp/BasisSet/PolarizedGaussian/MnD/Hermite1.H"
#include "Imp/BasisSet/PolarizedGaussian/MnD/Hermite3.H"
#include "Imp/BasisSet/PolarizedGaussian/MnD/RNLM.H"
#include "Imp/BasisSet/PolarizedGaussian/Radial/GaussianRF.H"
#include "Imp/BasisSet/PolarizedGaussian/Radial/GaussianCD.H"
#include "Imp/BasisSet/PolarizedGaussian/Block.H"
#include <Cluster.H>
#include "Misc/DFTDefines.H"
#include "oml/matrix.h"
#include "oml/io3d.h"
#include <cmath>
#include <cassert>

optr_vector1<std::vector<Polarization>* > GaussianCD::theNMLs;

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
    for (int L=0; L<=10; L++) theNMLs.push_back(new std::vector<Polarization>(MakeAllPolarizations(L)));
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

