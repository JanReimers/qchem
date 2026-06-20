// File: BasisSet/Molecule/PolarizedGaussian/Imp/Symmetry.C
module;
#include <vector>
#include <string>
#include <map>
module qchem.BasisSet.Molecule.PolarizedGaussian.Symmetry;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.GaussianRF;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;

namespace BasisSet::Molecule::PolarizedGaussian
{
using Symmetry::AoShell;
using Symmetry::SymPoint;
using Symmetry::IVec3;

// A center-independent id for a radial shell (L + exponents + coefficients): symmetry-
// equivalent shells on different atoms share it, so the center permutation can match them.
static int ShellTypeId(const GaussianRF* r, std::map<std::string,int>& types)
{
    const std::string& key = r->TypeID();
    auto it = types.find(key);
    if (it != types.end()) return it->second;
    int id = (int)types.size();
    types[key] = id;
    return id;
}

std::vector<AoShell> ExtractAoShells(const PGData& pg)
{
    std::vector<AoShell> shells;
    std::map<std::string,int> types;
    const size_t n = pg.size();
    size_t i = 0;
    while (i < n)
    {
        const GaussianRF* r = pg.radials[i];   // all functions of one block share this radial
        AoShell sh;
        sh.center    = r->GetCenter();
        sh.offset    = i;
        sh.shellType = ShellTypeId(r, types);
        size_t j = i;
        for (; j < n && pg.radials[j] == r; ++j)
        {
            const Polarization& p = pg.pols[j];
            sh.monomials.push_back(IVec3{p.n, p.l, p.m});
            sh.norm.push_back(pg.ns[j]);
        }
        shells.push_back(std::move(sh));
        i = j;
    }
    return shells;
}

std::vector<SymPoint> ClusterToSymPoints(const Cluster& cl)
{
    std::vector<SymPoint> pts;
    for (auto atom : cl) pts.push_back(SymPoint{atom->itsZ, atom->itsR});
    return pts;
}

} //namespace
