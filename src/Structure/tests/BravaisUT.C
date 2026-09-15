// File: src/Structure/tests/BravaisUT.C  The 14 Bravais lattices (doc/OpenWork.md row BL).
//
// Each type's PRIMITIVE cell must carry exactly its lattice system's point symmetry: one atom at the origin
// and SpaceGroup::Detect must find the holohedry -- 48 cubic, 24 hexagonal, 16 tetragonal, 12 rhombohedral,
// 8 orthorhombic, 4 monoclinic, 2 triclinic -- and nothing more (a wrong primitive vector lowers the count)
// or less (a parameter choice that accidentally raises the symmetry is a test defect, hence the ugly numbers).
// The detector is the oracle: it enumerates integer ops preserving the metric, so a correct cell cannot fail
// it and an incorrect one cannot pass it.
#include "gtest/gtest.h"
#include <vector>
#include <cmath>
#include <stdexcept>

import qchem.Symmetry.Lattice_3D.SpaceGroup;   // SpaceGroup::Detect, AtomSite
import qchem.UnitCell;                         // Bravais, LatticeParams, BravaisCellMatrix, BravaisCell, FCCUnitCell
import qchem.Types;

using namespace qchem;
namespace SL=qchem::Symmetry::Lattice_3D;

namespace
{
size_t Holohedry(Bravais type, const LatticeParams& p)
{
    std::vector<SL::AtomSite> one = {{14, rvec3_t(0,0,0)}};
    return SL::SpaceGroup::Detect(BravaisCellMatrix(type,p), one).Order();
}
double Volume(Bravais type, const LatticeParams& p) { return std::fabs(Determinant(BravaisCellMatrix(type,p))); }
}

TEST(Bravais, EveryTypeDetectsItsHolohedry)
{
    // Generic parameters on purpose: nothing commensurate, no accidental extra symmetry.
    EXPECT_EQ(Holohedry(Bravais::CubicP,        {.a=5.1}),                        48u);
    EXPECT_EQ(Holohedry(Bravais::CubicI,        {.a=5.1}),                        48u);
    EXPECT_EQ(Holohedry(Bravais::CubicF,        {.a=5.1}),                        48u);
    EXPECT_EQ(Holohedry(Bravais::TetragonalP,   {.a=5.1, .c=7.3}),                16u);
    EXPECT_EQ(Holohedry(Bravais::TetragonalI,   {.a=5.1, .c=7.3}),                16u);
    EXPECT_EQ(Holohedry(Bravais::OrthorhombicP, {.a=5.1, .b=6.2, .c=7.3}),         8u);
    EXPECT_EQ(Holohedry(Bravais::OrthorhombicC, {.a=5.1, .b=6.2, .c=7.3}),         8u);
    EXPECT_EQ(Holohedry(Bravais::OrthorhombicI, {.a=5.1, .b=6.2, .c=7.3}),         8u);
    EXPECT_EQ(Holohedry(Bravais::OrthorhombicF, {.a=5.1, .b=6.2, .c=7.3}),         8u);
    EXPECT_EQ(Holohedry(Bravais::HexagonalP,    {.a=5.1, .c=7.3}),                24u);
    EXPECT_EQ(Holohedry(Bravais::RhombohedralR, {.a=5.1, .α=67.0}),               12u);
    EXPECT_EQ(Holohedry(Bravais::MonoclinicP,   {.a=5.1, .b=6.2, .c=7.3, .β=101.0}), 4u);
    EXPECT_EQ(Holohedry(Bravais::MonoclinicC,   {.a=5.1, .b=6.2, .c=7.3, .β=101.0}), 4u);
    EXPECT_EQ(Holohedry(Bravais::TriclinicP,    {.a=5.1, .b=6.2, .c=7.3, .α=83.0, .β=101.0, .γ=97.0}), 2u);
}

// A centred type's primitive cell holds 1/2 (I, C) or 1/4 (F) of the conventional cell.
TEST(Bravais, CentredPrimitiveCellsHaveTheRightVolume)
{
    const double a=5.1, b=6.2, c=7.3;
    EXPECT_NEAR(Volume(Bravais::CubicI,        {.a=a}),             a*a*a/2, 1e-12);
    EXPECT_NEAR(Volume(Bravais::CubicF,        {.a=a}),             a*a*a/4, 1e-12);
    EXPECT_NEAR(Volume(Bravais::TetragonalI,   {.a=a, .c=c}),       a*a*c/2, 1e-12);
    EXPECT_NEAR(Volume(Bravais::OrthorhombicC, {.a=a, .b=b, .c=c}), a*b*c/2, 1e-12);
    EXPECT_NEAR(Volume(Bravais::OrthorhombicI, {.a=a, .b=b, .c=c}), a*b*c/2, 1e-12);
    EXPECT_NEAR(Volume(Bravais::OrthorhombicF, {.a=a, .b=b, .c=c}), a*b*c/4, 1e-12);
    EXPECT_NEAR(Volume(Bravais::HexagonalP,    {.a=a, .c=c}),       a*a*c*std::sqrt(3.0)/2, 1e-12);
    EXPECT_NEAR(Volume(Bravais::MonoclinicC,   {.a=a, .b=b, .c=c, .β=101.0}), a*b*c*std::sin(101.0*M_PI/180)/2, 1e-12);
}

// The one type the suite already has a class for: BravaisCell(CubicF) IS FCCUnitCell, bitwise.
TEST(Bravais, CubicFIsTheFCCUnitCell)
{
    const double a=10.26;
    const Matrix3D<double> A=FCCUnitCell(a).GetCellMatrix(), B=BravaisCellMatrix(Bravais::CubicF, {.a=a});
    for (int i=1;i<=3;i++) for (int j=1;j<=3;j++) EXPECT_EQ(A(i,j), B(i,j)) << i << "," << j;
}

// The superlattice re-basing: MnO's AFM-II cell is the FCC cell doubled along [111], i.e. CubicF re-based by
// T = [[0,1,1],[1,0,1],[1,1,0]] -- which must reproduce the nine literals the campaign wrote by hand
// (GPW_SCF_UT.C RunMnO: A = (a, a/2, a/2; a/2, a, a/2; a/2, a/2, a)), and carry the rhombohedral holohedry.
TEST(Bravais, SuperlatticeRebasingNamesTheMnOCell)
{
    const double a=8.40;
    const Matrix3D<int> T(0,1,1, 1,0,1, 1,1,0);
    UnitCell mno=BravaisCell(Bravais::CubicF, {.a=a}, T);
    const Matrix3D<double> hand(a, a/2, a/2,  a/2, a, a/2,  a/2, a/2, a);
    for (int i=1;i<=3;i++) for (int j=1;j<=3;j++) EXPECT_NEAR(mno.GetCellMatrix()(i,j), hand(i,j), 1e-14) << i << "," << j;
    EXPECT_NEAR(mno.GetCellVolume(), 2*a*a*a/4, 1e-9) << "twice the FCC primitive cell";
    std::vector<SL::AtomSite> one = {{25, rvec3_t(0,0,0)}};
    EXPECT_EQ(SL::SpaceGroup::Detect(mno.GetCellMatrix(), one).Order(), 12u) << "the [111]-doubled FCC cell is rhombohedral";
    EXPECT_THROW(BravaisCell(Bravais::CubicF, {.a=a}, Matrix3D<int>(1,1,0, 1,1,0, 0,0,1)), std::invalid_argument);   // singular T
}

// A parameter the type does not read is a caller error, not a silently different lattice.
TEST(Bravais, StrayParametersThrow)
{
    EXPECT_THROW(BravaisCellMatrix(Bravais::CubicF,      {.a=5.0, .c=7.0}),  std::invalid_argument);
    EXPECT_THROW(BravaisCellMatrix(Bravais::CubicF,      {.a=5.0, .α=60.0}), std::invalid_argument);
    EXPECT_THROW(BravaisCellMatrix(Bravais::HexagonalP,  {.a=5.0, .b=5.0, .c=7.0}), std::invalid_argument);
    EXPECT_THROW(BravaisCellMatrix(Bravais::TetragonalP, {.a=5.0}),          std::invalid_argument);   // c missing
    EXPECT_THROW(BravaisCellMatrix(Bravais::MonoclinicP, {.a=5.0, .b=6.0, .c=7.0, .γ=100.0}), std::invalid_argument);
    EXPECT_NO_THROW(BravaisCellMatrix(Bravais::RhombohedralR, {.a=5.0, .α=60.0}));
}
