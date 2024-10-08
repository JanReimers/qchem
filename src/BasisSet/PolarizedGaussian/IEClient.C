// File: PolarizedGaussian/IEClient.C

#include "Imp/BasisSet/PolarizedGaussian/IEClient.H"
#include "Imp/BasisSet/PolarizedGaussian/CDCache.H"
#include "Imp/BasisSet/PolarizedGaussian/Block.H"

namespace PolarizedGaussian
{
    
void IEData::Init(std::vector<const Block*>& blocks)
{
     for (auto bl:blocks)
        for (auto p:bl->itsPols)
        {
            radials.push_back(bl->itsRadial);
            pols.push_back(p);
        }
   
    
    size_t N=radials.size();
    ns.SetLimits(N);
    CDCache cache;
    for (size_t i=0;i<N;i++)
        ns(i+1)=radials[i]->Integrate(RadialFunction::Overlap2C,radials[i],pols[i],pols[i],cache);
    ns=1.0/sqrt(ns);
}
} //namespace PolarizedGaussian
