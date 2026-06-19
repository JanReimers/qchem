// File: BasisSet/Molecule/PolarizedGaussian/Internal/Imp/PGData.C
module;
#include <vector>
#include <string>
module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.PGData;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.CDCache;
import qchem.Blaze;

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
    ns=1.0/blazem::sqrt(ns);
}

// Cache identity of the whole molecular basis.  The radial/angular split is atomic bias (for an
// atom the centre is pinned, so radial x angular is a complete key); a molecule also needs the
// centres -- the 2-electron ERIs, overlap and kinetic are all orientation-dependent, so two
// geometries with the same basis must not collide in the cache.  We therefore fold everything --
// radial, centre and polarization per function -- into RadialID and leave AngularID empty.
std::string PGData::RadialID () const
{
    std::ostringstream os;
    os << " PG { ";
    for (size_t i=0;i<radials.size();++i)
        os << *radials[i] << "@" << radials[i]->GetCenter() << ":" << pols[i] << " ";
    os << "}";
    return os.str();
}
std::string PGData::AngularID() const { return ""; }

}