// File: BasisSet/Molecule/Evaluators/PG_Spherical_MnD/Imp/Evaluator.C
module;
#include <string>
#include <sstream>
module qchem.BasisSet.Molecule.Evaluators.PG_Spherical_MnD;
import qchem.Math;   // sqrt

namespace BasisSet::Molecule::Evaluators::PG_Spherical_MnD
{

// Per-component normalisation: 1/sqrt(<chi_i|chi_i>), the raw spherical self-overlap assembled from the
// Cartesian overlap of the component's monomials.
void SphData::Init()
{
    ns.resize(comps.size());
    for (size_t i=0;i<comps.size();++i)
    {
        double raw = 0.0;
        for (const auto& ta : comps[i].terms)
            for (const auto& tb : comps[i].terms)
                raw += ta.c * tb.c * comps[i].radial->Overlap2C(*comps[i].radial, ta.p, tb.p);
        ns[i] = 1.0/sqrt(raw);
    }
}

// Cache identity.  Streams each component's radial (exponents/centre) and its monomial expansion; the
// "PGSph" tag keeps it distinct from the Cartesian PG basis so the global matrix cache never confuses
// the two (the underlying primitive/Omega caches ARE shared -- that is correct and desirable).
std::string SphData::RadialID() const
{
    std::ostringstream os;
    os << " PGSph { ";
    for (const auto& c : comps)
    {
        os << *c.radial << " [";
        for (const auto& t : c.terms) os << t.c << "*" << t.p << " ";
        os << "] ";
    }
    os << "}";
    return os.str();
}

} //namespace
