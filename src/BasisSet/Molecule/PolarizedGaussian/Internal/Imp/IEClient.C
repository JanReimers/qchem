// File: PolarizedGaussian/IEClient.C
module;
#include <vector>
#include <blaze/Math.h>
module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.IEClient;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.CDCache;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Block;

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
    ns.resize(N);
    CDCache cache;
    for (size_t i=0;i<N;i++)
        ns[i]=radials[i]->Integrate(Overlap2C,radials[i],pols[i],pols[i],cache);
    ns=1.0/sqrt(ns);
}

} //namespace PolarizedGaussian
