// File: src/Structure/tests/SiteGroupsUT.C  The two site groups a Hubbard manifold is resolved under
// (doc/OpenWork.md step 5 increment 2): Lattice_3D::SiteRotations (the decoration's Shubnikov stabiliser,
// a SPACE-GROUP property) and Lattice_3D::SiteEnvironmentRotations (the coordination polyhedron's point
// group, a CHEMICAL property).  The MnO AFM-II supercell is the case that separates them: rhombohedral, so
// its space-group stabiliser at Mn is D_3d (12) with OR without decoration, while the MnO6 octahedron is
// O_h (48) -- the t2g/e_g parentage a d level is named by lives in the environment, not the cell
// (measured 2026-09-21: the tree's "grey stabiliser = O_h" claim was false for the supercell).
#include "gtest/gtest.h"
#include <vector>
#include <cmath>

import qchem.Lattice_3D;    // Lattice_3D, UnitCell, FCCUnitCell, BravaisCell
import qchem.Types;

using namespace qchem;

namespace
{
bool Same(const rmat3d_t& A, const rmat3d_t& B)
{
    for (int i=1;i<=3;i++) for (int j=1;j<=3;j++) if (std::fabs(A(i,j)-B(i,j))>1e-8) return false;
    return true;
}
bool Contains(const std::vector<rmat3d_t>& G, const rmat3d_t& R) { for (const auto& g : G) if (Same(g,R)) return true; return false; }
bool IsOrthogonal(const rmat3d_t& R)
{
    const rmat3d_t G=Transpose(R)*R;
    for (int i=1;i<=3;i++) for (int j=1;j<=3;j++) if (std::fabs(G(i,j)-(i==j?1.0:0.0))>1e-8) return false;
    return true;
}
//! MnO AFM-II: the FCC cell doubled along [111], Mn at 0 and (1/2,1/2,1/2), O at (1/4,1/4,1/4) and (3/4,3/4,3/4)
//! in the doubled cell's fractional coordinates.
UnitCell MnO_AFM2(double a=8.40)
{
    UnitCell cell=BravaisCell(Bravais::CubicF, {.a=a}, Matrix3D<int>(0,1,1, 1,0,1, 1,1,0));
    cell.AddAtom(25, {0,0,0});
    cell.AddAtom(25, {0.5,0.5,0.5}, true);
    cell.AddAtom(8,  {0.25,0.25,0.25});
    cell.AddAtom(8,  {0.75,0.75,0.75});
    return cell;
}
}

TEST(SiteGroups, MnOAFM2_TheCellStabiliserIsD3dEitherWayAndTheEnvironmentIsOh)
{
    Lattice_3D lat(MnO_AFM2(), ivec3_t(1,1,1));
    const std::vector<int> afm{+1,-1,0,0};
    const auto site=lat.SiteRotations(0, afm);
    const auto grey=lat.SiteRotations(0, {});
    const auto env =lat.SiteEnvironmentRotations(0);
    EXPECT_EQ(site.size(), 12u) << "D_3d: the Shubnikov stabiliser of the ordered Mn";
    EXPECT_EQ(grey.size(), 12u) << "the SUPERCELL's grey stabiliser is D_3d too -- it cannot name t2g/e_g";
    EXPECT_EQ(env.size(),  48u) << "the MnO6 octahedron + fcc Mn12 shell: O_h";
    for (const auto& R : env) EXPECT_TRUE(IsOrthogonal(R));
    // The site group is a subgroup of the environment group (the slot table needs Tr(P_site P_grey) integral).
    for (const auto& R : site) EXPECT_TRUE(Contains(env, R));
    for (const auto& R : grey) EXPECT_TRUE(Contains(env, R));
    // The second Mn sees the same chemistry.
    EXPECT_EQ(lat.SiteEnvironmentRotations(1).size(), 48u);
    // O sits in an Mn6 octahedron: O_h as well.
    EXPECT_EQ(lat.SiteEnvironmentRotations(2).size(), 48u);
}

TEST(SiteGroups, DiamondSi_TheEnvironmentIsTd_AndAgreesWithTheCellStabiliser)
{
    FCCUnitCell cell(10.26);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    const auto site=lat.SiteRotations(0, {});
    const auto env =lat.SiteEnvironmentRotations(0);
    EXPECT_EQ(site.size(), 24u);
    EXPECT_EQ(env.size(),  24u) << "the primitive cell already carries the full site symmetry: T_d both ways";
    for (const auto& R : site) EXPECT_TRUE(Contains(env, R));
    // Exactly one shell (the 4 nearest neighbours alone) is a regular tetrahedron: still T_d, and no more.
    EXPECT_EQ(lat.SiteEnvironmentRotations(0, 1).size(), 24u);
}
