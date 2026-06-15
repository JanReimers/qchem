// File: UnitTests/M_PGSymmetry.C  PG basis -> symmetry bridge (stage 5 wiring) on a real basis.
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include "blaze/Math.h"
import qchem.BasisSet.Molecule.PolarizedGaussian.Symmetry;  // ExtractAoShells, ClusterToSymPoints
import qchem.BasisSet.Molecule.PolarizedGaussian;            // Orbital_IBS
import qchem.Symmetry.SALC;                                   // BuildAbelianGroup, BuildSALCs, BuildOperationRep
import qchem.Cluster;                                         // Molecule, Atom
import qchem.Types;

using namespace BasisSet::Molecule::PolarizedGaussian;
using namespace Symmetry;

// Build a real Cartesian-Gaussian basis on H2O, extract its shells, and run the full
// symmetry pipeline (detect -> abelian group -> SALCs).  Validates the bridge end to end.
TEST(PGSymmetry, water_extract_and_SALCs)
{
    Molecule h2o;
    h2o.Insert(new Atom(8, 0, rvec3_t(0, 0,      0.117)));   // O on the C2 (z) axis
    h2o.Insert(new Atom(1, 0, rvec3_t(0, 0.757, -0.467)));   // H
    h2o.Insert(new Atom(1, 0, rvec3_t(0,-0.757, -0.467)));   // H

    rvec_t exps{1.0, 0.25};
    Orbital_IBS ibs(exps, 1, &h2o);                          // s + p shells, 2 exponents each

    auto shells = ExtractAoShells(ibs);
    auto pts    = ClusterToSymPoints(h2o);

    // nAO = 3 atoms * (2 s-functions + 2 p-shells * 3) = 3 * 8 = 24.
    size_t nAO = 0;
    for (const auto& s : shells) nAO = std::max(nAO, s.offset + s.monomials.size());
    EXPECT_EQ(nAO, 24u);
    EXPECT_EQ(pts.size(), 3u);

    auto g = BuildAbelianGroup(pts, 1e-4);
    EXPECT_EQ(g.table.symbol, "C2v");

    rvec3_t o = Centroid(pts);
    auto salc = BuildSALCs(shells, g, o, 1e-4);
    EXPECT_EQ(salc.O.columns(), nAO);     // SALCs span the full AO space

    // Defining property: each SALC column v of irrep Gamma satisfies M(g) v = chi^Gamma(g) v.
    std::map<std::string,size_t> irow;
    for (size_t r=0;r<g.table.irreps.size();++r) irow[g.table.irreps[r]] = r;
    std::vector<rmat_t> M;
    for (const auto& op : g.ops) M.push_back(BuildOperationRep(shells, op.Matrix(), o, 1e-4));

    for (size_t c=0;c<salc.O.columns();++c)
    {
        blaze::DynamicVector<double> v = column(salc.O, c);
        size_t r = irow[salc.irrep[c]];
        for (size_t k=0;k<M.size();++k)
        {
            blaze::DynamicVector<double> Mv  = M[k]*v;
            blaze::DynamicVector<double> chv = double(g.table.chi[r][k]) * v;
            EXPECT_LT(norm(Mv - chv), 1e-8) << "irrep " << salc.irrep[c] << " op " << g.table.opTags[k];
        }
    }
}
