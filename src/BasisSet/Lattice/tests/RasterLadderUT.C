// File: RasterLadderUT.C  THE RASTER LADDER EQUIVALENCE (V-BS-RASTER step 1): the ONE ladder MakeRasterLadder builds per geometry point
// against the ladder every GPW block builds for itself today, level by level, to the bit.
//
// The blocks' ladder is read from what they FILE (grids.ladder, the run report EmitGpwGrids writes) -- the only face it has -- so this
// gate needs no access to the evaluator's internals.  Once the evaluator is wired to the module and its own copy is gone, this becomes
// the regression test of the one implementation.  Cases span the ladder's branches: Si Γ and k222 (narrow exponent range, the
// completion rung), a MnO-scale rocksalt (a wide range: several coarse levels), an O2 box (a big cubic cell), the HartreeOnly routing,
// the explicit mgridEcuts list and a non-default progression factor.
#include "gtest/gtest.h"
#include <memory>
#include <stdexcept>
#include <string>
import qchem.UnitCell;
import qchem.Lattice_3D;
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Point.Factory;      // Gaussian::Factory (the molecular basis)
import qchem.BasisSet.Lattice.BasisSet;            // GPWFactory, GPWParams, EmitGpwGrids, MakeRasterLadder
import qchem.Reporting;                            // the run the blocks file their grids into (a test may read it)
import qchem.Types;

using namespace qchem;
using BasisSet::Real_BS;
using BasisSet::Gaussian::BasisSetData;

namespace
{
std::shared_ptr<const Real_BS> Mol(BasisSetData d, const Structure& st)
{
    return std::shared_ptr<const Real_BS>(BasisSet::Gaussian::Factory(d, &st, BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
}

//! Every level of the blocks' filed ladder equals the shared ladder's: FFT divisions, cutoff (bitwise), |{G}| and role.
void ExpectSameLadder(const Lattice_3D& lat, std::shared_ptr<const Real_BS> mol, const BasisSet::Lattice::GPWParams& p, const std::string& name)
{
    report::json filed;
    {
        report::Run run("raster ladder "+name);
        std::unique_ptr<BasisSet::Complex_BS> bs(BasisSet::Lattice::GPWFactory(lat, mol, p));
        BasisSet::Lattice::EmitGpwGrids(*bs);
        filed=report::CurrentRunReport()["grids"]["ladder"];
    }
    const auto ladder=BasisSet::Lattice::MakeRasterLadder(lat, mol, p);
    ASSERT_TRUE(ladder->DensityGrid()) << name;
    const auto& L=ladder->Levels();
    EXPECT_EQ(L.levels.front(), ladder->DensityGrid()) << name << ": level 0 IS the density grid";
    ASSERT_EQ(filed.size(), L.levels.size()) << name << ": level count";
    for (size_t l=0; l<L.levels.size(); l++)
    {
        const auto& row=filed[l];
        const ivec3_t N=L.N[l];
        EXPECT_EQ(row["N"][0].get<int>(), N.x) << name << " L" << l;
        EXPECT_EQ(row["N"][1].get<int>(), N.y) << name << " L" << l;
        EXPECT_EQ(row["N"][2].get<int>(), N.z) << name << " L" << l;
        EXPECT_EQ(row["ecut"].get<double>(), L.ecut[l]) << name << " L" << l << ": the cutoff, bitwise";
        EXPECT_EQ(row["nG"].get<long>(), long(L.levels[l]->size())) << name << " L" << l;
        const std::string role = l==0 ? "reference" : (l>=L.nBase ? "rung" : "coarse");
        EXPECT_EQ(row["role"].get<std::string>(), role) << name << " L" << l;
    }
}

FCCUnitCell SiDiamond()
{
    FCCUnitCell cell(10.26);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    return cell;
}
FCCUnitCell MnORocksalt()
{
    FCCUnitCell cell(8.40);
    cell.AddAtom(25, {0,0,0});
    cell.AddAtom(8,  {0.5,0.5,0.5});
    return cell;
}
UnitCell O2Box()
{
    const double a=16.0, d=2.282;
    UnitCell cell(a);
    cell.AddAtom(8, {0.5-0.5*d/a,0.5,0.5});
    cell.AddAtom(8, {0.5+0.5*d/a,0.5,0.5});
    return cell;
}
}

TEST(RasterLadder, SiGammaAndK222MatchTheBlocksLadder)
{
    const auto cell=SiDiamond();
    const Lattice_3D g(cell, ivec3_t(1,1,1)), k(cell, ivec3_t(2,2,2));
    const auto mol=Mol(BasisSetData::SIPP, *g.GetStructure());
    ExpectSameLadder(g, mol, {.densityEcut=20.0}, "Si Γ Ecut 20");
    ExpectSameLadder(g, mol, {},                  "Si Γ auto floor");
    ExpectSameLadder(k, mol, {.densityEcut=20.0}, "Si k222 Ecut 20");
}

TEST(RasterLadder, AWideExponentRangeMatchesLevelForLevel)
{
    const auto cell=MnORocksalt();
    const Lattice_3D lat(cell, ivec3_t(1,1,1));
    const auto mol=Mol(BasisSetData::VALENCE_LOWQ_SR, *lat.GetStructure());
    ExpectSameLadder(lat, mol, {}, "MnO auto floor");
    ExpectSameLadder(lat, mol, {.rasterFields=BasisSet::PlaneWave::RasterFields::HartreeOnly}, "MnO HartreeOnly");
}

TEST(RasterLadder, ABoxedMoleculeAndTheLadderKnobsMatch)
{
    const auto cell=O2Box();
    const Lattice_3D lat(cell, ivec3_t(1,1,1));
    const auto mol=Mol(BasisSetData::VALENCE_LOWQ_SR, *lat.GetStructure());
    ExpectSameLadder(lat, mol, {.densityEcut=40.0}, "O2 box Ecut 40");
    ExpectSameLadder(lat, mol, {.densityEcut=40.0, .ladderFactor=3.0}, "O2 box progression 3");
    BasisSet::Lattice::GPWParams mg{.densityEcut=40.0};
    mg.tol.mgridEcuts={13.33, 4.444};
    ExpectSameLadder(lat, mol, mg, "O2 box mgridEcuts");
}

//! densityEcut 0 is the 1E-only tier: no grid, and asking for a ladder is a broken invariant (a throw, not a silent empty).
TEST(RasterLadder, NoDensityGridMeansNoLadder)
{
    const auto cell=SiDiamond();
    const Lattice_3D lat(cell, ivec3_t(1,1,1));
    const auto ladder=BasisSet::Lattice::MakeRasterLadder(lat, Mol(BasisSetData::SIPP, *lat.GetStructure()), {.densityEcut=0.0});
    EXPECT_FALSE(ladder->DensityGrid());
    EXPECT_THROW((void)ladder->Levels(), std::logic_error);
}
