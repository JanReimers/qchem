// File: BasisSet/Molecule/PolarizedGaussian/Internal/Imp/PGData.C
module;
#include <vector>
#include <blaze/Math.h>
module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.PGData;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.CDCache;

namespace PolarizedGaussian
{
    void PGData::Init(std::vector<const Block*>& blocks)
{
     for (auto bl:blocks)
        for (auto p:bl->itsPols)
        {
            radials1.push_back(bl->itsRadial);
            pols1.push_back(p);
        }
   
    
    size_t N=radials1.size();
    ns1.resize(N);
    CDCache cache;
    for (size_t i=0;i<N;i++)
        ns1[i]=radials1[i]->Integrate(Overlap2C,radials1[i],pols1[i],pols1[i],cache);
    ns1=1.0/sqrt(ns1);
}

}