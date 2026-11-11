// File: M_Gaussian94Reader.C  The Gaussian-94 basis reader names an element its file does not hold.
//
// Found 2026-10-21: TiO2 on VALENCE_LOWQ_VA (a file with no Ti).  FindAtom ran to EOF and reported "not found" through a bool every
// caller ignored, so the next ReadNext parsed whatever bytes it stood on and died as "'\xef' is not an angular-momentum label".
// FindAtom now THROWS, naming the element and the file.  The file here is written by the test, so no data-file edit can move it.

#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
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
import qchem.PeriodicTable;                      // the element of each header line
import qchem.Mesh;                               // qcMesh::MeshParams (the Vxc fit's factory takes one)

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
    const Site H{1, Vector3D(0,0,0)};
    reader.FindAtom(H);
    std::unique_ptr<GaussianRF> rf{reader.ReadNext(H)};
    EXPECT_TRUE(rf) << "the H s shell";
}

TEST(Gaussian94Reader, AnElementThatIsNotInTheFileThrowsNamingTheElementAndTheFile)
{
    const auto file=HydrogenOnlyFile();
    Gaussian94Reader reader(file.string());
    const Site Ti{22, Vector3D(0,0,0)};
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
    Molecule m; m.Insert(Site{1, Vector3D(0,0,0)});
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
    n2.Insert(Site{7, Vector3D(-1.03,0.0,0.0)});
    n2.Insert(Site{7, Vector3D( 1.04,0.0,0.0)});
    auto cartesianCount=[&](const std::string& file)
    {
        Gaussian94Reader reader(BasisFile(file));
        size_t n=0;
        for (const auto& sh : ReadShells(reader, n2)) for (int l : sh.Ls) n+=size_t((l+1)*(l+2)/2);
        return n;
    };
    std::unique_ptr<BasisSet::Real_BS> bs(Factory(BasisSetData::DZVP, &n2));
    const qcMesh::MeshParams mp{};
    std::unique_ptr<BasisSet::rFIT_CD_ABS> cd(bs->CreateCDFitBasisSet(&n2));
    std::unique_ptr<BasisSet::rFIT_SF_ABS> sf(bs->CreateVxcFitBasisSet(&n2, mp));
    EXPECT_EQ(cd->GetNumFunctions(), cartesianCount("A1_coul.bsd")) << "the Coulomb fit basis dropped shells of its file";
    EXPECT_EQ(sf->GetNumFunctions(), cartesianCount("A1_exch.bsd")) << "the Vxc fit basis dropped shells of its file";
    // ...and the d shells are really there: a d shell on a shared exponent adds 6 Cartesian functions per atom
    Gaussian94Reader reader(BasisFile("A1_coul.bsd"));
    size_t nd=0; for (const auto& sh : ReadShells(reader, n2)) for (int l : sh.Ls) if (l==2) ++nd;
    EXPECT_GT(nd, 0u) << "A1_coul N has d shells (on the SP exponents) -- else this test tests nothing";
}

// THE COULOMB (ρ) FIT IS 100% ANALYTIC (user 2026-11-11, "The rho fit should be 100% analytic, no grid"): its factory takes no mesh,
// its Norm is the analytic metric diagonal's (1 to roundoff: the shells are analytically normalised), and it has no points to project a
// field on, so that refuses with a reason.  OverlapDiagonal refuses on any Gaussian fit (R1.0: delete, don't correct).  The Vxc fit
// keeps its projection mesh and answers <f_a|f> there.
TEST(Gaussian94Reader, TheCoulombFitIsAnalyticAndOwnsNoMesh)
{
    struct One : ScalarFunction<double>
    {
        double  operator()(const rvec3_t&) const override {return 1.0;}
        rvec3_t Gradient  (const rvec3_t&) const override {return rvec3_t(0,0,0);}
    } one;
    Molecule n2;
    n2.Insert(Site{7, Vector3D(-1.03,0.0,0.0)});
    n2.Insert(Site{7, Vector3D( 1.04,0.0,0.0)});
    std::unique_ptr<BasisSet::Real_BS> bs(Factory(BasisSetData::DZVP, &n2));
    std::unique_ptr<BasisSet::rFIT_CD_ABS> cd(bs->CreateCDFitBasisSet(&n2));
    const auto& asSF=dynamic_cast<const BasisSet::FIT_SF_NonOrtho&>(*cd);   // the Dunlap fitter's overlap-solve face
    const rvec_t& n=asSF.Norm();
    const rsmat_t& S=asSF.Overlap();
    ASSERT_EQ(n.size(), cd->GetNumFunctions());
    for (size_t a=0; a<n.size(); a++)
    {
        EXPECT_EQ(n[a], 1.0/std::sqrt(S(a,a))) << "function " << a << ": the norm is the analytic metric's";
        EXPECT_NEAR(n[a], 1.0, 1e-12) << "function " << a << ": the shells are analytically normalised";
    }
    EXPECT_THROW((void)asSF.Overlap(one), std::logic_error) << "a Coulomb fit owns no mesh to project a field on";
    EXPECT_THROW((void)asSF.OverlapDiagonal(), std::logic_error) << "R1.0: a non-orthogonal fit refuses a diagonal-only metric";
    const qcMesh::MeshParams mesh{.radial=qcMesh::RadialKind::MHL, .nRadial=80, .mhl_m=3, .mhl_alpha=2.0, .angular=qcMesh::AngularKind::Lebedev,
                                 .angularDegree=29, .beckeOrder=2};   // the molecular entry rung
    std::unique_ptr<BasisSet::rFIT_SF_ABS> sf(bs->CreateVxcFitBasisSet(&n2, mesh));
    const rvec_t p=sf->Overlap(one);   // the Vxc fit projects on its own mesh: <f_a|1> by quadrature, against the analytic Charge
    ASSERT_EQ(p.size(), sf->GetNumFunctions());
    const rvec_t q=sf->Charge();
    double worst=0.0; for (size_t a=0; a<p.size(); a++) worst=std::max(worst, std::fabs(p[a]-q[a])/std::max(1.0,std::fabs(q[a])));
    EXPECT_LT(worst, 1e-4) << "the Vxc fit's projection of 1 is its charge, to the mesh's quadrature error";
}

TEST(Gaussian94Reader, AnOrbitalBasisLoadsEveryShellOfASharedExponent)
{
    Molecule m; m.Insert(Site{1, Vector3D(0,0,0)});
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
    Molecule m; m.Insert(Site{1, Vector3D(0,0,0)});
    std::unique_ptr<qchem::BasisSet::tBasisSet<double>> bs(Factory(nlohmann::json{{"basis","DZVP"}}, &m));
    EXPECT_GT(bs->GetNumFunctions(), 0u);
    EXPECT_THROW((void)Factory(nlohmann::json{{"basis","nope"}}, &m), std::runtime_error);
    // The checksum is of the CONTENTS: stable, and different for a different file.
    EXPECT_EQ(BasisFileChecksum("dzvp.bsd"), BasisFileChecksum("dzvp.bsd"));
    EXPECT_NE(BasisFileChecksum("dzvp.bsd"), BasisFileChecksum("dzvp2.bsd"));
    EXPECT_EQ(BasisFileChecksum("dzvp.bsd").size(), 16u);
}

// THE BODY IS LOCATED BY ITS HEADER (lane G 2026-11-10, OpenWork §3): the reader once stopped at the first "!" TOKEN, skipped two lines
// and then to the next "****" -- a preamble opening with its own "!" lines silently cost the file its FIRST block (decks/panel/
// MnO2_control.bsd lost F).  Now the body starts after BASIS=, every block is found, the historical "!"+title layout still reads, and
// a file with neither is refused by name.
namespace
{
std::filesystem::path Wrote(const std::string& name, const std::string& text)
{
    const auto p=std::filesystem::temp_directory_path()/("qchem_UT_Gaussian94Reader_"+name+".bsd");
    std::ofstream(p) << text;
    return p;
}
size_t Count(const std::filesystem::path& f, int Z)
{
    Molecule m; m.Insert(Site{Z, Vector3D(0,0,0)});
    std::unique_ptr<qchem::BasisSet::tBasisSet<double>> bs(FactoryFromFile(f.string(), &m, Engine::MnD, Angular::Cartesian));
    return bs->GetNumFunctions();
}
const char* kTwoBlocks=" F   0\n S   1  1.00\n   1.0   1.0\n ****\n O   0\n S   1  1.00\n   2.0   1.0\n P   1  1.00\n   1.0   1.0\n ****\n";
}

TEST(Gaussian94Reader, APreambleThatOpensWithBangLinesLosesNoBlock)
{
    const auto f=Wrote("bang_preamble", std::string("!\n! a note\nfree text\n!\nBASIS=\"t\"\n")+kTwoBlocks);
    EXPECT_EQ(Count(f, 9), 1u) << "the FIRST block (F) is found";
    EXPECT_EQ(Count(f, 8), 4u);
}

TEST(Gaussian94Reader, TheHistoricalBangTitleLayoutStillReads)
{
    const auto f=Wrote("bang_title", std::string("free text\n!\nA title, no BASIS= line\n")+kTwoBlocks);
    EXPECT_EQ(Count(f, 9), 1u);
    EXPECT_EQ(Count(f, 8), 4u);
}

TEST(Gaussian94Reader, AFileWithNoLocatableBodyIsRefusedByName)
{
    const auto f=Wrote("no_body", std::string("<pre>! not a marker\n<p>BASIS=\"t\"\n")+kTwoBlocks);
    try { (void)Count(f, 8); FAIL() << "a file with no BASIS= line and no \"!\" line was read"; }
    catch (const std::runtime_error& e) { EXPECT_NE(std::string(e.what()).find("cannot be located"), std::string::npos) << e.what(); }
}

TEST(Gaussian94Reader, EveryShippedFileReadsEveryBlockItHolds)
{
    for (const std::string& name : BasisNames())
    {
        std::ifstream in(BasisFile(name+".bsd"));
        std::string line; bool body=false;
        for (; std::getline(in, line);)
        {
            std::istringstream is(line); std::string a, b, extra;
            if (line.find("BASIS=")!=std::string::npos) { body=true; continue; }
            if (!body || !(is >> a >> b) || (is >> extra) || b!="0") continue;
            std::string sym=a; for (char& c : sym) c=(char)std::tolower((unsigned char)c);
            sym[0]=(char)std::toupper((unsigned char)sym[0]);
            const int Z=int(thePeriodicTable().GetZ(sym));
            if (Z==0) continue;
            Gaussian94Reader reader(BasisFile(name+".bsd"));
            EXPECT_NO_THROW(reader.FindAtom(Site{Z, Vector3D(0,0,0)})) << name << ": " << a;
        }
    }
}
