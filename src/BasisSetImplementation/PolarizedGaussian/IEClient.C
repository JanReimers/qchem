// File: PolarizedGaussian/IEClient.C

#include "BasisSetImplementation/PolarizedGaussian/IEClient.H"

void PolarizedGaussianIEClient::Init(std::vector<const BasisFunctionBlock*>& blocks)
{
     for (auto bl:blocks)
        for (auto p:bl->itsPols)
        {
            radials.push_back(bl->itsRadial);
            pols.push_back(p);
        }
   
    size_t N=size();
    ns.SetLimits(N);
    CDcache cache;
    for (size_t i=0;i<N;i++)
        ns(i+1)=radials[i]->Integrate(RadialFunction::Overlap2C,radials[i],pols[i],pols[i],cache);
    ns=1.0/sqrt(ns);
}
