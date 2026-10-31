// File: M_Gaussian94Reader.C  The Gaussian-94 basis reader names an element its file does not hold.
//
// Found 2026-10-21: TiO2 on VALENCE_LOWQ_VA (a file with no Ti).  FindAtom ran to EOF and reported "not found" through a bool every
// caller ignored, so the next ReadNext parsed whatever bytes it stood on and died as "'\xef' is not an angular-momentum label".
// FindAtom now THROWS, naming the element and the file.  The file here is written by the test, so no data-file edit can move it.

#include "gtest/gtest.h"
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>
#include <nlohmann/json.hpp>

import qchem.BasisSet.Gaussian.Point.Readers.Gaussian94;
import qchem.BasisSet.Gaussian.Point.Reader;    // ReadShells
import qchem.BasisSet.Gaussian.Point.Factory;   // FactoryFromFile -- the orbital bases end to end
import qchem.BasisSet.Gaussian.Point.BasisFiles; // BasisNames / BasisFileOfName / BasisFileChecksum (D-BASIS-REGISTRY)
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.Structure;
import qchem.Types;
import qchem.Mesh;                               // qcMesh::MeshParams (the fit bases' factories take one)

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

TEST(Gaussian94Reader, ReadShellsKeepsEachLOfASharedExponent)
{
    Molecule m; m.Insert(new Atom(1, 0.0, Vector3D(0,0,0)));
    Gaussian94Reader reader(SharedExponentFile().string());
    const auto sh=ReadShells(reader, m);
    ASSERT_EQ(sh.size(), 2u) << "the s and the p are two blocks";
    EXPECT_EQ(sh[0].Ls, std::vector<int>{0});
    EXPECT_EQ(sh[1].Ls, std::vector<int>{1});
}

// THE FIT BASES READ EVERY SHELL (2026-10-31, OpenWork "The molecular FIT bases drop their D shells"): every A1 file puts SP and D
// shells on the SAME exponents, and the fit bases' old merge never added the new l -- N2's A1_coul loaded 32 functions, not 68, and
// the Coulomb fit sat 4-62 mHa low against PySCF's own A1_coul density fit.  Invariant: the fit basis has exactly the Cartesian
// functions of the file's shells for each atom, counted here from the file through the reader, independently of the basis.
TEST(Gaussian94Reader, TheFitBasesCarryEveryShellOfTheirFile)
{
    Molecule n2;
    n2.Insert(new Atom(7, 0.0, Vector3D(-1.03,0.0,0.0)));
    n2.Insert(new Atom(7, 0.0, Vector3D( 1.04,0.0,0.0)));
    auto cartesianCount=[&](const std::string& file)
    {
        Gaussian94Reader reader(BasisFile(file));
        size_t n=0;
        for (const auto& sh : ReadShells(reader, n2)) for (int l : sh.Ls) n+=size_t((l+1)*(l+2)/2);
        return n;
    };
    std::unique_ptr<BasisSet::Real_BS> bs(Factory(BasisSetData::DZVP, &n2));
    const qcMesh::MeshParams mp{};
    std::unique_ptr<BasisSet::rFIT_CD_ABS> cd(bs->CreateCDFitBasisSet(&n2, mp));
    std::unique_ptr<BasisSet::rFIT_SF_ABS> sf(bs->CreateVxcFitBasisSet(&n2, mp));
    EXPECT_EQ(cd->GetNumFunctions(), cartesianCount("A1_coul.bsd")) << "the Coulomb fit basis dropped shells of its file";
    EXPECT_EQ(sf->GetNumFunctions(), cartesianCount("A1_exch.bsd")) << "the Vxc fit basis dropped shells of its file";
    // ...and the d shells are really there: a d shell on a shared exponent adds 6 Cartesian functions per atom
    Gaussian94Reader reader(BasisFile("A1_coul.bsd"));
    size_t nd=0; for (const auto& sh : ReadShells(reader, n2)) for (int l : sh.Ls) if (l==2) ++nd;
    EXPECT_GT(nd, 0u) << "A1_coul N has d shells (on the SP exponents) -- else this test tests nothing";
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

// C4 step 4 (D-BASIS-PP): a file's PP-PROVENANCE preamble line is checked against the run's {element, q}.  A q3 run on a q1
// block THROWS naming file, element and both q; an element or a file that declares nothing is accepted; the tag is invisible
// to the reader (it sits before the lone `!`).
TEST(PPProvenance, AQ3RunOnAQ1BlockThrowsNamingBoth)
{
    const std::string path=(std::filesystem::path(testing::TempDir())/(std::to_string(::getpid())+"_provenance.bsd")).string();
    {
        std::ofstream f(path);
        f<<"title\nPP-PROVENANCE LI=1 F=7\n\n!\nBASIS=\"t\"\n";
    }
    EXPECT_NO_THROW(CheckPPProvenance(path, {{"Li",1},{"F",7}}));
    EXPECT_NO_THROW(CheckPPProvenance(path, {{"Na",9}})) << "an element the file does not declare is accepted";
    try { CheckPPProvenance(path, {{"Li",3}}); FAIL() << "q3 on a q1 block must throw"; }
    catch (const std::runtime_error& e)
    {
        const std::string m=e.what();
        EXPECT_NE(m.find("Li"),std::string::npos) << m;
        EXPECT_NE(m.find("q1"),std::string::npos) << m;
        EXPECT_NE(m.find("q3"),std::string::npos) << m;
    }
    const std::string bare=(std::filesystem::path(testing::TempDir())/(std::to_string(::getpid())+"_noprov.bsd")).string();
    { std::ofstream f(bare); f<<"title\n\n!\nBASIS=\"t\"\n"; }
    EXPECT_NO_THROW(CheckPPProvenance(bare, {{"Li",3}})) << "an untagged (legacy) file declares nothing";
}

// D-BASIS-REGISTRY (CleanCode, 2026-10-30): a basis is NAMED by its file stem in the data directory -- every .bsd there is
// usable by name with no registration; the enum is only the C++ handle of the files code names directly, and each handle's
// name is a stem that is there.  A name that is not there throws, listing the names that are.
TEST(BasisNames, ABasisIsAFileStemAndEveryHandleNamesOneThatIsThere)
{
    const std::vector<std::string> names=BasisNames();
    auto has=[&](const std::string& n){ return std::find(names.begin(), names.end(), n)!=names.end(); };
    EXPECT_TRUE(has("valence_lowq_tm")) << "a file with no handle is still a name";
    for (BasisSetData d : {BasisSetData::DZVP, BasisSetData::DZVP2, BasisSetData::TZVP, BasisSetData::ORB, BasisSetData::ORB1,
                           BasisSetData::SIPP, BasisSetData::SIPP_SR, BasisSetData::VALENCE_LOWQ, BasisSetData::VALENCE_LOWQ_SR,
                           BasisSetData::VALENCE_LOWQ_SR2, BasisSetData::VALENCE_LOWQ_SPH, BasisSetData::VALENCE_LOWQ_VA,
                           BasisSetData::VALENCE_LOWQ_VB})
    {
        const std::string stem=BasisName(d);
        EXPECT_TRUE(has(stem)) << stem << ": a handle names a file that is not in the data directory";
        std::string upper=stem; for (char& c : upper) c=(char)std::toupper((unsigned char)c);
        BasisSetData back{};
        ASSERT_TRUE(RegisteredBasis(upper, back)) << upper;
        EXPECT_EQ(back, d) << "the handle's name, in any case, gives the handle back";
    }
    BasisSetData unused{};
    EXPECT_FALSE(RegisteredBasis("valence_lowq_tm", unused));
    EXPECT_EQ(BasisFileOfName("Valence_LowQ_TM"), "valence_lowq_tm.bsd") << "any case";
    try { (void)BasisFileOfName("valence_lowq_typo"); FAIL() << "an unknown name must throw"; }
    catch (const std::runtime_error& e)
    {
        const std::string m=e.what();
        EXPECT_NE(m.find("valence_lowq_typo"), std::string::npos) << m;
        EXPECT_NE(m.find("valence_lowq_tm"),   std::string::npos) << "the names that ARE there are listed: " << m;
    }
    // The json door resolves the same names: one H atom on the dzvp file, by stem.
    Molecule m; m.Insert(new Atom(1, 0.0, Vector3D(0,0,0)));
    std::unique_ptr<qchem::BasisSet::tBasisSet<double>> bs(Factory(nlohmann::json{{"basis","DZVP"}}, &m));
    EXPECT_GT(bs->GetNumFunctions(), 0u);
    EXPECT_THROW((void)Factory(nlohmann::json{{"basis","nope"}}, &m), std::runtime_error);
    // The checksum is of the CONTENTS: stable, and different for a different file.
    EXPECT_EQ(BasisFileChecksum("dzvp.bsd"), BasisFileChecksum("dzvp.bsd"));
    EXPECT_NE(BasisFileChecksum("dzvp.bsd"), BasisFileChecksum("dzvp2.bsd"));
    EXPECT_EQ(BasisFileChecksum("dzvp.bsd").size(), 16u);
}
