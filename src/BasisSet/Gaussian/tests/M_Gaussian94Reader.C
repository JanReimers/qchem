// File: M_Gaussian94Reader.C  The Gaussian-94 basis reader names an element its file does not hold.
//
// Found 2026-10-21: TiO2 on VALENCE_LOWQ_VA (a file with no Ti).  FindAtom ran to EOF and reported "not found" through a bool every
// caller ignored, so the next ReadNext parsed whatever bytes it stood on and died as "'\xef' is not an angular-momentum label".
// FindAtom now THROWS, naming the element and the file.  The file here is written by the test, so no data-file edit can move it.

#include "gtest/gtest.h"
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

import qchem.BasisSet.Gaussian.Point.Readers.Gaussian94;
import qchem.BasisSet.Gaussian.Point.Reader;    // ReadShells / ShellMerge
import qchem.BasisSet.Gaussian.Point.Factory;   // FactoryFromFile -- the orbital bases end to end
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.Structure;
import qchem.Types;

using namespace qchem;
using namespace qchem::BasisSet::Gaussian;
using namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD;

namespace
{
//! A one-element (H) Gaussian-94 file: the "!" title block the reader seeks, one s shell, the "****" terminator.
std::filesystem::path HydrogenOnlyFile()
{
    const auto p=std::filesystem::temp_directory_path() / "qchem_UT_Gaussian94Reader_H_only.bsd";
    std::ofstream f(p);
    f << "header text the reader skips\n!\nBASIS=\"UT: hydrogen only\"\n H   0\n S   1  1.00\n      1.0   1.0\n****\n";
    return p;
}
}

TEST(Gaussian94Reader, AnElementThatIsInTheFileIsFoundAndRead)
{
    Gaussian94Reader reader(HydrogenOnlyFile().string());
    const Atom H(1, 0.0, Vector3D(0,0,0));
    reader.FindAtom(H);
    std::unique_ptr<GaussianRF> rf{reader.ReadNext(H)};
    EXPECT_TRUE(rf) << "the H s shell";
}

TEST(Gaussian94Reader, AnElementThatIsNotInTheFileThrowsNamingTheElementAndTheFile)
{
    const auto file=HydrogenOnlyFile();
    Gaussian94Reader reader(file.string());
    const Atom Ti(22, 0.0, Vector3D(0,0,0));
    try
    {
        reader.FindAtom(Ti);
        FAIL() << "FindAtom must throw for an element the file does not hold";
    }
    catch (const std::runtime_error& e)
    {
        const std::string what=e.what();
        EXPECT_NE(what.find("Ti"), std::string::npos) << what;
        EXPECT_NE(what.find(file.string()), std::string::npos) << what;
    }
}


// TWO SHELLS OF DIFFERENT l SHARING ONE EXPONENT ON ONE ATOM (OpenWork §3, fixed 2026-10-26).  GaussianRF::operator== ignores L, so
// the s and the p below are "the same radial"; the historical merge (kept for the density-fit bases) never added the new l and
// loaded ONE function where an orbital basis needs four.  Bit NiO's scan the same day: Ni p 0.30 beside the Ni d window starting
// at 0.30 loaded +2 functions, not +12.
namespace
{
std::filesystem::path SharedExponentFile()
{
    const auto p=std::filesystem::temp_directory_path() / "qchem_UT_Gaussian94Reader_H_sp_shared.bsd";
    std::ofstream f(p);
    f << "header\n!\nBASIS=\"UT: H s + p at ONE exponent\"\n H   0\n S   1  1.00\n      1.0   1.0\n P   1  1.00\n      1.0   1.0\n****\n";
    return p;
}
}

TEST(Gaussian94Reader, ReadShellsKeepsEachLForAnOrbitalBasisAndMergesForAFitBasis)
{
    Molecule m; m.Insert(new Atom(1, 0.0, Vector3D(0,0,0)));
    {
        Gaussian94Reader reader(SharedExponentFile().string());
        const auto sh=ReadShells(reader, m, ShellMerge::KeepEachL);
        ASSERT_EQ(sh.size(), 2u) << "the s and the p are two blocks";
        EXPECT_EQ(sh[0].Ls, std::vector<int>{0});
        EXPECT_EQ(sh[1].Ls, std::vector<int>{1});
    }
    {
        Gaussian94Reader reader(SharedExponentFile().string());
        const auto sh=ReadShells(reader, m, ShellMerge::LegacyMerge);
        ASSERT_EQ(sh.size(), 1u) << "the fit bases' historical merge, kept bit-identical";
        EXPECT_EQ(sh[0].Ls, std::vector<int>{0}) << "...which never adds the new l (the documented quirk)";
    }
}

TEST(Gaussian94Reader, AnOrbitalBasisLoadsEveryShellOfASharedExponent)
{
    Molecule m; m.Insert(new Atom(1, 0.0, Vector3D(0,0,0)));
    const std::string file=SharedExponentFile().string();
    auto count=[&](Engine e, Angular a)
    {
        std::unique_ptr<qchem::BasisSet::tBasisSet<double>> bs(FactoryFromFile(file, &m, e, a));
        size_t n=0;
        for (auto ibs : bs->template Iterate<qchem::BasisSet::Real_OIBS>()) n+=ibs->GetNumFunctions();
        return n;
    };
    EXPECT_EQ(count(Engine::MnD,     Angular::Cartesian), 4u) << "1 s + 3 p";
    EXPECT_EQ(count(Engine::MnD,     Angular::Spherical), 4u);
    EXPECT_EQ(count(Engine::LibCint, Angular::Cartesian), 4u);
}
