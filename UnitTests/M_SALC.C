// File: UnitTests/M_SALC.C  Symmetry-adapted linear combinations (stage 3b).
#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <map>
#include "blaze/Math.h"
import qchem.Symmetry.SALC;          // BuildSALCs, SALCs
import qchem.Symmetry.AbelianGroup;  // BuildAbelianGroup
import qchem.Symmetry.PointGroup;    // SymPoint
import qchem.Math;

using namespace Symmetry;

// Water, C2 along z, molecule in the yz-plane.  AO basis: O s [0], O p [1,2,3], H1 s [4],
// H2 s [5].  Centers shared between the nuclear point set and the basis shells.
static rvec3_t O_(0,0,0.117), H1_(0,0.757,-0.467), H2_(0,-0.757,-0.467);

static std::vector<SymPoint> WaterPts()
{
    return { {8,O_}, {1,H1_}, {1,H2_} };
}
static std::vector<AoShell> WaterShells()
{
    return {
        {0, O_,  {{0,0,0}},                  {1.0},         0},
        {1, O_,  {{1,0,0},{0,1,0},{0,0,1}},  {1.0,1.0,1.0}, 1},
        {2, H1_, {{0,0,0}},                  {1.0},         4},
        {2, H2_, {{0,0,0}},                  {1.0},         5},
    };
}
static rvec3_t WaterOrigin() { return (O_+H1_+H2_)/3.0; }

TEST(SALC, water_spans_and_block_dims)
{
    auto g = BuildAbelianGroup(WaterPts(), 1e-6);
    auto s = BuildSALCs(WaterShells(), g, WaterOrigin(), 1e-6);

    EXPECT_EQ(s.O.columns(), 6u);     // SALCs span the full AO space
    EXPECT_EQ(s.O.rows(),    6u);

    std::map<std::string,int> dim;
    for (const auto& lab : s.irrep) dim[lab]++;
    EXPECT_EQ(dim["A1"], 3);          // O s, O pz, (H1+H2)
    EXPECT_EQ(dim["A2"], 0);
    EXPECT_EQ(dim["B1"] + dim["B2"], 3);   // O px, O py, (H1-H2) split across b1/b2
    EXPECT_EQ(dim["B1"] * dim["B2"], 2);   // {1,2} in some order (convention-dependent)
}

// The defining property: each SALC column v of irrep Gamma satisfies M(g) v = chi^Gamma(g) v.
TEST(SALC, columns_are_irrep_eigenvectors)
{
    auto shells = WaterShells();
    rvec3_t o = WaterOrigin();
    auto g = BuildAbelianGroup(WaterPts(), 1e-6);
    auto s = BuildSALCs(shells, g, o, 1e-6);

    // index of each irrep label in the character table
    std::map<std::string,size_t> irow;
    for (size_t r=0;r<g.table.irreps.size();++r) irow[g.table.irreps[r]] = r;

    std::vector<rmat_t> M;
    for (const auto& op : g.ops) M.push_back(BuildOperationRep(shells, op.Matrix(), o, 1e-6));

    for (size_t c=0;c<s.O.columns();++c)
    {
        blaze::DynamicVector<double> v = column(s.O, c);
        size_t r = irow[s.irrep[c]];
        for (size_t k=0;k<M.size();++k)
        {
            blaze::DynamicVector<double> Mv  = M[k]*v;
            blaze::DynamicVector<double> chv = double(g.table.chi[r][k]) * v;
            EXPECT_LT(norm(Mv - chv), 1e-9) << "irrep " << s.irrep[c] << " col " << c
                                            << " op " << g.table.opTags[k];
        }
    }
}
