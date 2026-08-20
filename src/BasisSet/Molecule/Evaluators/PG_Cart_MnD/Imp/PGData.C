// File: BasisSet/Molecule/Evaluators/PG_Cart_MnD/Imp/PGData.C
module;
#include <vector>
#include <string>
// R1.9: <sstream> (and the <ostream> it pulls in) MUST be included here even though the code below
// compiled without it.  ostream::operator<<(const void*) is a MEMBER and is always found; the one for
// const char* is a FREE FUNCTION TEMPLATE in <ostream>, and a C++20 module does not export the free
// operators a header-including TU happens to have seen.  Without this include every string literal
// below silently printed as a hex ADDRESS -- and BasisSetID() is a cache key.
#include <sstream>
#include <cmath>    // std::pow/std::fabs (the Reaches() scan)
module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;
import qchem.Blaze;

namespace qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD
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
    BuildReaches();          // eager: the pointwise screen must be ready before any threaded sweep
}

// The per-function magnitude-screen radius (see the interface doc).  Scanned, not solved: a CONTRACTED
// radial is a sum of Gaussians and need not be monotone, so we walk out and keep the LAST radius whose
// bound still exceeds eps -- never the first one below it.  One-time and cheap (~1200 radial evaluations
// per function, against the millions of image evaluations it screens).
void PGData::BuildReaches()
{
    const double eps=1e-10;    // the house pointwise magnitude floor (BetaSupportRadius uses the same)
    const double h=0.05, rmax=60.0;
    itsReach.resize(size());
    for (size_t i=0;i<size();i++)
    {
        const GaussianRF& rf=*radials[i];
        const rvec3_t     c =rf.GetCenter();
        const int         L =pols[i].GetTotalL();
        double last=0.0;
        for (double d=0.0; d<=rmax; d+=h)
        {
            // |chi_i| <= ns * d^L * |R(d)|: R is spherically symmetric and |x^a y^b z^c| <= d^L over
            // every direction, so one radial ray bounds the whole sphere.
            const double bound=std::fabs(ns[i])*std::pow(d,double(L))*std::fabs(rf(c+rvec3_t(d,0,0)));
            if (bound>eps) last=d;
        }
        itsReach[i]=last+2.0*h;   // a small margin past the last significant radius
    }
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