// File: RasterLadderUT.C  THE RASTER LADDER (V-BS-RASTER step 1): the ONE ladder MakeRasterLadder builds per geometry point is the one
// every GPW block of that point runs on, and the one that reports the grids.
//
// What the ladder FILES (grids.ladder, written into the caller's "grids" scope by EmitGpwGrids) is its own levels, bit for bit; a
// factory handed a ladder shares it, and refuses one built from other settings; V_loc long's audit row is a pure function of the
// ladder and the PP sharpness it was sized with.  Cases span the ladder's branches: Si Γ and k222 (narrow exponent range, the
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
import qchem.BasisSet.Lattice.BasisSet;            // GPWFactory, GPWParams, EmitGpwGrids, MakeRasterLadder, ResolveImposition
import qchem.Pseudopotential.GTH_Potentials;       // GetGTH (a local PP to size V_long's ladder for)
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

//! Every level the ladder files (grids.ladder) is its own: FFT divisions, cutoff (bitwise), |{G}| and role.
void ExpectSameLadder(const Lattice_3D& lat, std::shared_ptr<const Real_BS> mol, const BasisSet::Lattice::GPWParams& p, const std::string& name)
{
    const auto ladder=BasisSet::Lattice::MakeRasterLadder(lat, mol, p);
    report::json filed;
    {
        report::Run run("raster ladder "+name);
        {
            report::Section grids("grids");
            BasisSet::Lattice::EmitGpwGrids(*ladder, *mol, ResolveImposition(lat, {}));
        }
        filed=report::CurrentRunReport()["grids"]["ladder"];
    }
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

TEST(RasterLadder, SiGammaAndK222FileTheirOwnLevels)
{
    const auto cell=SiDiamond();
    const Lattice_3D g(cell, ivec3_t(1,1,1)), k(cell, ivec3_t(2,2,2));
    const auto mol=Mol(BasisSetData::SIPP, *g.GetStructure());
    ExpectSameLadder(g, mol, {.densityEcut=20.0}, "Si Γ Ecut 20");
    ExpectSameLadder(g, mol, {},                  "Si Γ auto floor");
    ExpectSameLadder(k, mol, {.densityEcut=20.0}, "Si k222 Ecut 20");
}

TEST(RasterLadder, AWideExponentRangeFilesLevelForLevel)
{
    const auto cell=MnORocksalt();
    const Lattice_3D lat(cell, ivec3_t(1,1,1));
    const auto mol=Mol(BasisSetData::VALENCE_LOWQ_SR, *lat.GetStructure());
    ExpectSameLadder(lat, mol, {}, "MnO auto floor");
    ExpectSameLadder(lat, mol, {.rasterFields=BasisSet::PlaneWave::RasterFields::HartreeOnly}, "MnO HartreeOnly");
}

TEST(RasterLadder, ABoxedMoleculeAndTheLadderKnobsFileTheirLevels)
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

//! A factory handed a ladder runs every block on it (the same density grid object), and refuses a ladder built from other settings
//! rather than run on it: two statements of the ladder must agree.
TEST(RasterLadder, AHandedInLadderIsSharedAndAMismatchIsRefused)
{
    const auto cell=SiDiamond();
    const Lattice_3D k(cell, ivec3_t(2,2,2));
    const auto mol=Mol(BasisSetData::SIPP, *k.GetStructure());
    BasisSet::Lattice::GPWParams p{.densityEcut=20.0};
    p.ladder=BasisSet::Lattice::MakeRasterLadder(k, mol, p);
    std::unique_ptr<BasisSet::Complex_BS> bs(BasisSet::Lattice::GPWFactory(k, mol, p));
    EXPECT_GT(bs->GetNumIBS(), 1u);
    BasisSet::Lattice::GPWParams q=p;
    q.densityEcut=30.0;                           // the params now state another ball than the ladder was built for
    EXPECT_THROW(std::unique_ptr<BasisSet::Complex_BS>(BasisSet::Lattice::GPWFactory(k, mol, q)), std::invalid_argument);
}

//! V_loc long's audit row is the LADDER's: present exactly when it was sized for a local PP with a closed Gaussian sharpness, and a
//! pure function of that -- so it is the same whether or not the V_long sweep has run (a cached matrix used to hide it).
TEST(RasterLadder, VlocLongsAuditIsAFunctionOfTheLadderAndThePP)
{
    const auto cell=SiDiamond();
    const Lattice_3D lat(cell, ivec3_t(1,1,1));
    const auto mol=Mol(BasisSetData::SIPP, *lat.GetStructure());
    const auto gth=Pseudopotential::GetGTH("Si","LDA",4);
    const double beta=BasisSet::Gaussian::LocalPPSharpness(*lat.GetStructure(), gth.local);
    ASSERT_GT(beta, 0.0) << "the Si GTH local has a closed Gaussian short part";
    const auto bare =BasisSet::Lattice::MakeRasterLadder(lat, mol, {.densityEcut=20.0});
    const auto sized=BasisSet::Lattice::MakeRasterLadder(lat, mol, {.densityEcut=20.0}, beta);
    EXPECT_EQ(bare ->DeliveredRasterVerdicts().size(), 1u) << "no PP sharpness: the density row only";
    const auto v=sized->DeliveredRasterVerdicts();
    ASSERT_EQ(v.size(), 2u);
    EXPECT_EQ(v[1].audit, "auditVlocLong");
    EXPECT_EQ(bare->DeliveredRasterExponents().vlocLong, 0.0);
    EXPECT_GT(sized->DeliveredRasterExponents().vlocLong, 0.0);
    EXPECT_EQ(sized->DeliveredRasterExponents().density, bare->DeliveredRasterExponents().density) << "the PP touches only V_long's row";
}
