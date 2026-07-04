// File: BasisSet/Molecule/PG_Spherical/Imp/Symmetry.C
module;
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <utility>
module qchem.BasisSet.Molecule.PG_Spherical.Symmetry;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;           // Cart::GaussianRF (TypeID/GetCenter)
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;         // Cart::Polarization (the monomial exps)
import qchem.BasisSet.Molecule.Evaluators.PG_Spherical_MnD.SolidHarmonics;  // CartTerm
import qchem.Symmetry.Molecule.SphericalRep;   // HarmonicC2S, SphericalShellRep (the concrete ShellRep this basis produces)
import qchem.Blaze;                   // blazem::VecBuilder (accumulate the per-shell norms into an rvec_t)

namespace qchem::BasisSet::Molecule::PG_Spherical
{
namespace Cart = ::qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD;
namespace Sph  = ::qchem::BasisSet::Molecule::Evaluators::PG_Spherical_MnD;
using Symmetry::Molecule::AoShell;
using Symmetry::Molecule::IVec3;
using Symmetry::Molecule::HarmonicC2S;
using Symmetry::Molecule::SphericalShellRep;

// A center-independent id for a radial shell (L + exponents + coefficients): symmetry-equivalent shells on
// different atoms share it, so the center permutation can match them.  (Same logic as PG_Cart's helper; the
// radial is the same Cart::GaussianRF type -- kept file-local, a shared bridge helper can wait for S3b.)
static int ShellTypeId(const Cart::GaussianRF* r, std::map<std::string,int>& types)
{
    const std::string& key = r->TypeID();
    auto it = types.find(key);
    if (it != types.end()) return it->second;
    int id = (int)types.size();
    types[key] = id;
    return id;
}

std::vector<AoShell> ExtractAoShells(const Sph::SphData& sph)
{
    std::vector<AoShell> shells;
    std::map<std::string,int> types;
    const size_t n = sph.size();
    size_t i = 0;
    while (i < n)
    {
        const Cart::GaussianRF* r = sph.comps[i].radial;   // all harmonics of one shell share this radial
        AoShell sh;
        sh.center    = r->GetCenter();
        sh.offset    = i;
        sh.shellType = ShellTypeId(r, types);
        HarmonicC2S                c2s;
        blazem::VecBuilder<double> norm;
        size_t j = i;
        for (; j < n && sph.comps[j].radial == r; ++j)
        {
            // This harmonic's Cartesian expansion IS the basis's own c2s (comps[j].terms), m-ordering
            // preserved -> the rep is built in the basis's exact convention (no re-derivation, no mismatch).
            std::vector<std::pair<IVec3,double>> harmonic;
            for (const auto& t : sph.comps[j].terms)
                harmonic.push_back({IVec3{t.p.n, t.p.l, t.p.m}, t.c});
            c2s.push_back(std::move(harmonic));
            norm.Append(sph.ns[j]);                  // per-harmonic normalization (folded into the rep)
        }
        sh.rep  = std::make_shared<SphericalShellRep>(std::move(c2s));   // this shell's spherical rep
        sh.norm = norm.take();
        shells.push_back(std::move(sh));
        i = j;
    }
    return shells;
}

} //namespace
