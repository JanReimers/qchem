// File: GaussianCD.C  Charge distribution for two  gaussians.


#include "oml/smatrix.h"
#include "BasisSetImplementation/PolarizedGaussian/Hermite/Hermite1.H"
#include "BasisSetImplementation/PolarizedGaussian/Gaussian/GaussianRF.H"
#include "BasisSetImplementation/PolarizedGaussian/Gaussian/GaussianCD.H"
#include "BasisSetImplementation/PolarizedGaussian/BasisFunctionBlock.H"
#include "BasisSetImplementation/PolarizedGaussian/Hermite/Hermite1.H"
#include "BasisSetImplementation/PolarizedGaussian/Hermite/Hermite3.H"
#include "BasisSetImplementation/PolarizedGaussian/Auxillary/RNLM.H"
#include "Cluster.H"
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
    , lookups(0)
{
    if (theNMLs.size()==0) MakeNMLs();
};

GaussianCD::~GaussianCD()
{
    
    for (auto r:RNLMs) delete r.second;
}

void GaussianCD::Report(std::ostream& os) const
{
    if (lookups>0)
    {
        double eff=(100.0*lookups)/(RNLMs.size()+lookups-1);
        os << "   RNLMs cache N=" << RNLMs.size() << " lookups=" << lookups << " efficiencty=" << eff << "%" << std::endl;
    }

}

// third Center at c.
const RNLM& GaussianCD::GetRNLM(const RVec3& c) const
{
    auto i=RNLMs.find(c);
    if (i==RNLMs.end())
    {
        RNLMs[c]=new RNLM(Ltotal,ab/AlphaP,c);
        i=RNLMs.find(c);
    }
    else
        assert(i->second->CheckLMax(Ltotal));
    
    lookups++;
    return *i->second;
}

Vector<double>  GaussianCD::GetRNLMs(const Polarization& p, const Hermite1& H1) const
{
     const RNLM& R=GetRNLM(AB);
     //RNLM R(Ltotal,alphaR,AB);
    // std::cout << "GaussianCD::GetRNLMs " << Ltotal << " " << alphaR << " " << AB << std::endl;
     auto NLMs=GetNMLs(r1.GetL());
     Vector<double> ret(NLMs.size());
     Fill(ret,0.0);
     double factor=1.0/(ab*sqrt(AlphaP));
     double sign = (p.GetTotalL()%2) ? -factor : factor;

    for (int n=0; n<=p.n; n++)
        for (int l=0; l<=p.l; l++)
            for (int m=0; m<=p.m; m++)
            {
                Polarization NLMp(n,l,m);
                double h=H1(NLMp,p);
               // std::cout << "GaussianCD::GetRNLMs h,sign=" << h << " " << sign << std::endl;
                if (h!=0.0)
                {
                    h*=sign;
                    std::vector<Polarization>::const_iterator bNLM(NLMs.begin());
                    for (int nNLM=1; bNLM!=NLMs.end(); bNLM++,nNLM++) 
                    {
//                        std::cout << "GaussianCD::GetRNLMs " << p << *bNLM << " " << NLMp << std::endl;
                        ret(nNLM)+=h*R(*bNLM+NLMp);
                    }
                }
            }
     return ret;
}

std::ostream& GaussianCD::Write(std::ostream& os) const {return os;}
std::istream& GaussianCD::Read (std::istream& is)       {return is;}

