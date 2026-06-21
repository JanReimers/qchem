// File: BasisSet/Molecule/Evaluators/PG_Cart_MnD/Imp/PGData.C
module;
#include <vector>
#include <string>
module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;
import qchem.Blaze;

namespace BasisSet::Molecule::Evaluators::PG_Cart_MnD
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
    for (size_t i=0;i<N;i++)
        ns[i]=radials[i]->Overlap2C(*radials[i],pols[i],pols[i]);
    ns=1.0/blazem::sqrt(ns);
}

// RadialID / AngularID are the basis-only pieces (exponents/contraction and polarizations); they
// carry no geometry.  BasisSetID below is the cache identity and folds in the centres.
std::string PGData::RadialID () const
{
    std::ostringstream os;
    os << " PG { ";   // distinct from PG so the global integral cache never serves PG's
                       // matrices to PG (and vice-versa): each tree computes independently.
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

// Geometry-aware cache identity of the whole molecular basis.  An atom's radial x angular is a
// complete key (centre pinned at the nucleus), but a molecule also needs the centres: overlap,
// kinetic and the 2-electron ERIs are all orientation-dependent, so two geometries with the same
// basis must not collide.  We fold radial, centre and polarization per function into one string.
std::string PGData::BasisSetID() const
{
    std::ostringstream os;
    os << " PG { ";   // distinct from PG so the global integral cache never serves PG's
                       // matrices to PG (and vice-versa): each tree computes independently.
    for (size_t i=0;i<radials.size();++i)
        os << *radials[i] << "@" << radials[i]->GetCenter() << ":" << pols[i] << " ";
    os << "}";
    return os.str();
}

}