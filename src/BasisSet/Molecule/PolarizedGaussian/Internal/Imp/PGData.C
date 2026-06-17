// File: BasisSet/Molecule/PolarizedGaussian/Internal/Imp/PGData.C
module;
#include <vector>
#include <blaze/Math.h>
module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.PGData;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.CDCache;

namespace BasisSet::Molecule::PolarizedGaussian
{
    void PGData::Init(std::vector<const Block*>& blocks)
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

std::string PGData::RadialID () const
{
    std::ostringstream os;
    os << " PG { ";
    for (auto r:radials) os << *r << " ";
    os << "}";
    return os.str();
}
std::string PGData::AngularID() const
{
    std::ostringstream os;
    os << "{ ";
    for (auto p:pols) os << p << " ";
    os << "}";
    return os.str();
}

}