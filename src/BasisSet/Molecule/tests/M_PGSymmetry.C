// File: BasisSet/Molecule/tests/M_PGSymmetry.C  PG basis -> symmetry bridge (stage 5 wiring) on a real basis.
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <map>
#include <algorithm>
import qchem.BasisSet.Molecule.PolarizedGaussian.Symmetry;  // ExtractAoShells, ClusterToSymPoints
import qchem.BasisSet.Molecule.PolarizedGaussian;            // Orbital_IBS
import qchem.BasisSet.SymmetryAdapted_IBS;                    // SymmetryAdapted_IBS (1-e decorator)
import qchem.BasisSet.Molecule.SymmetryAdaptedBasisSet;       // SymmetryAdaptedBasisSet (per-irrep)
import qchem.Symmetry.SALC;                                   // BuildAbelianGroup, BuildSALCs, BuildOperationRep
import qchem.Cluster;                                         // Molecule, Atom
import qchem.Types;
import qchem.Blaze;
import qchem.Math;                                            // fabs

using namespace BasisSet::Molecule::PolarizedGaussian;
using namespace Symmetry;
using SymmetryAdapted_IBS      = ::BasisSet::SymmetryAdapted_IBS;          // ::BasisSet (the class clashes)
using SymmetryAdaptedBasisSet  = ::BasisSet::Molecule::SymmetryAdaptedBasisSet;

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
        rvec_t v = blazem::column(salc.O, c);
        size_t r = irow[salc.irrep[c]];
        for (size_t k=0;k<M.size();++k)
        {
            rvec_t Mv  = M[k]*v;
            rvec_t chv = double(g.table.chi[r][k]) * v;
            EXPECT_LT(blazem::norm(Mv - chv), 1e-8) << "irrep " << salc.irrep[c] << " op " << g.table.opTags[k];
        }
    }
}

// The 1-electron decorator: each irrep's SymmetryAdapted_IBS must return the matching block
// of O^T S_raw O, computed from the REAL overlap integrals -- and O must block-diagonalize
// that real overlap (cross-irrep elements vanish).
TEST(PGSymmetry, decorator_blocks_real_overlap)
{
    Molecule h2o;
    h2o.Insert(new Atom(8, 0, rvec3_t(0, 0,      0.117)));
    h2o.Insert(new Atom(1, 0, rvec3_t(0, 0.757, -0.467)));
    h2o.Insert(new Atom(1, 0, rvec3_t(0,-0.757, -0.467)));
    rvec_t exps{1.0, 0.25};
    Orbital_IBS ibs(exps, 1, &h2o);

    auto shells = ExtractAoShells(ibs);
    auto pts    = ClusterToSymPoints(h2o);
    auto g      = BuildAbelianGroup(pts, 1e-4);
    rvec3_t o   = Centroid(pts);
    auto salc   = BuildSALCs(shells, g, o, 1e-4);

    const rsmat_t& Sraw = ibs.Overlap();        // REAL overlap integrals (normalized, 24x24)
    size_t nAO = Sraw.rows();
    ASSERT_EQ(nAO, salc.O.columns());

    // O block-diagonalizes the real overlap: cross-irrep elements vanish.
    rmat_t OtSO = blazem::trans(salc.O) * Sraw * salc.O;
    for (size_t i=0;i<nAO;i++) for (size_t j=0;j<nAO;j++)
        if (salc.irrep[i] != salc.irrep[j])
            EXPECT_LT(fabs(OtSO(i,j)), 1e-8) << "cross-irrep " << salc.irrep[i] << "/" << salc.irrep[j];

    // Each irrep's decorator returns its diagonal block of O^T S_raw O.
    for (size_t r=0; r+1<salc.blockStart.size(); ++r)
    {
        size_t start = salc.blockStart[r], dG = salc.blockStart[r+1]-start;
        if (dG==0) continue;
        rmat_t Or = blazem::submatrix(salc.O, 0, start, nAO, dG);
        SymmetryAdapted_IBS sa(&ibs, Or, salc.irrep[start], ibs.GetSymt());
        EXPECT_EQ(sa.GetNumFunctions(), dG);
        rsmat_t Sblk = sa.MakeOverlap();
        for (size_t a=0;a<dG;a++) for (size_t b=0;b<dG;b++)
            EXPECT_LT(fabs(Sblk(a,b) - OtSO(start+a,start+b)), 1e-9);
    }
}

// The top-level SymmetryAdaptedBasisSet: one labelled IrrepBasisSet per non-empty irrep,
// iterable exactly like an atomic basis, the irrep blocks summing to the full AO space.
TEST(PGSymmetry, symmetry_adapted_basis_set)
{
    Molecule h2o;
    h2o.Insert(new Atom(8, 0, rvec3_t(0, 0,      0.117)));
    h2o.Insert(new Atom(1, 0, rvec3_t(0, 0.757, -0.467)));
    h2o.Insert(new Atom(1, 0, rvec3_t(0,-0.757, -0.467)));
    rvec_t exps{1.0, 0.25};
    Orbital_IBS ibs(exps, 1, &h2o);

    auto shells = ExtractAoShells(ibs);
    auto pts    = ClusterToSymPoints(h2o);
    auto g      = BuildAbelianGroup(pts, 1e-4);
    auto salc   = BuildSALCs(shells, g, Centroid(pts), 1e-4);

    SymmetryAdaptedBasisSet sab(&ibs, salc);

    // Iterate the irrep IBSs: blocks sum to the full AO count, each carries a valid C2v label.
    size_t nIBS=0, nfunc=0;
    for (auto oi : sab.Iterate<::BasisSet::Real_OIBS>())
    {
        ++nIBS;
        nfunc += oi->GetNumFunctions();
        std::string lab = oi->GetSymmetry().GetLabel();
        EXPECT_TRUE(lab=="A1"||lab=="A2"||lab=="B1"||lab=="B2") << "unexpected irrep label: " << lab;
    }
    EXPECT_EQ(nfunc, ibs.GetNumFunctions());        // SALC blocks span the AO space
    EXPECT_EQ(sab.GetNumFunctions(), ibs.GetNumFunctions());
    EXPECT_GE(nIBS, 2u);                            // water (s+p) populates A1, B1, B2
}
