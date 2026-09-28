// Saved periodic SCF states (CK-1; qchem.SolidState + SolidCalculation::SaveState/Restart).
//
// The gates, in the order the design argues them:
//   - an EXACT RESUME reconverges at once to the energy it was saved at (the restored D IS the state);
//   - the file is what the layout says (h5py's reader sees the same): one group per (k,σ) block, f summing
//     to the block's electrons, D = C f C^† exactly;
//   - a WARM-class difference restarts and says what differs; a REFUSE-class one is a failed Outcome built
//     from nothing more than the basis (no SCF paid);
//   - a complex-saved state restarts onto REAL TRIM blocks (D real at a TRIM point) and back.
// Si diamond, GTH SIPP_SR, densityEcut 20 -- the Response kernel gates' cell, on a 2x2x2 mesh so every block
// is TRIM (real by default) and there is more than one block to get in the wrong order.
#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <map>
#include <utility>
#include <memory>
#include <string>
#include <vector>
import qchem.SolidCalculation;
import qchem.SolidState;
import qchem.Lattice_3D;
import qchem.Structure;
import qchem.UnitCell;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.SCFParams;
import qchem.SCFIterator;   // SCFProgress (the first-iterate probe)
import qchem.HDF5;
import qchem.Types;

using namespace qchem;
namespace {

std::string TempPath(const std::string& leaf) {return (std::filesystem::path(testing::TempDir())/leaf).string();}

struct SiRun
{
    std::unique_ptr<Lattice_3D>                   lat;
    std::shared_ptr<const BasisSet::Real_BS>      mol;
};
SiRun Si(ivec3_t N)
{
    FCCUnitCell cell(10.26);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    SiRun r;
    r.lat=std::make_unique<Lattice_3D>(cell, N);
    r.mol=std::shared_ptr<const BasisSet::Real_BS>(
        BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::SIPP_SR, &cell,
                                    BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
    return r;
}
SCFParams Params()
{
    SCFParams p;
    p.NMaxIter=80; p.MinΔρ=1e-7; p.MinΔE=1e-10; p.MinΔFD=1e-7; p.MinVirial=1e30; p.MinFD=1e30;
    return p;
}
SolidCalcOptions Options(int multiplicity, bool forceComplex, const std::string& save)
{
    return SolidCalcOptions{.Nelec=8, .multiplicity=multiplicity, .species={{"Si",4}}, .densityEcut=20.0,
                            .forceComplex=forceComplex, .label="ck1", .saveStateTo=save};
}

//! One converged-and-saved state per (spin, ansatz), built once for the whole suite.
struct Saved
{
    std::string path;
    double      energy=0.0;
    size_t      iterations=0;
};
Saved& SavedState(int multiplicity, bool forceComplex)
{
    static std::map<std::pair<int,bool>, Saved> cache;
    auto key=std::make_pair(multiplicity, forceComplex);
    if (auto i=cache.find(key); i!=cache.end()) return i->second;
    Saved s;
    s.path=TempPath("ck1_si_m"+std::to_string(multiplicity)+(forceComplex?"_c":"_r")+".h5");
    SiRun si=Si(ivec3_t(2,2,2));
    SolidCalculation calc(*si.lat, si.mol, Options(multiplicity, forceComplex, s.path), Params());
    auto r=calc.Result();
    EXPECT_TRUE(r) << "the Si reference run must converge";
    if (r) { s.energy=r->Energy(); s.iterations=r->IterationCount(); }
    return cache[key]=s;
}

void ExpectExactResume(int multiplicity, bool forceComplex)
{
    const Saved& s=SavedState(multiplicity, forceComplex);
    ASSERT_TRUE(std::filesystem::exists(s.path));
    SiRun si=Si(ivec3_t(2,2,2));
    auto c=SolidCalculation::Restart(s.path, *si.lat, si.mol, Options(multiplicity, forceComplex, ""), Params());
    ASSERT_TRUE(c) << c.Error().details;
    auto r=(*c)->Result();
    ASSERT_TRUE(r);
    EXPECT_NEAR(r->Energy(), s.energy, 1e-8);
    // The restored density IS the converged one: the SCF only has to confirm it (ΔE needs two iterates).
    EXPECT_LE(r->IterationCount(), 3u) << "the reference run took " << s.iterations;
    EXPECT_LT(r->IterationCount(), s.iterations);
}

TEST(SolidState, ExactResume_UnPol_Real)  {ExpectExactResume(0, false);}
TEST(SolidState, ExactResume_Pol_Complex) {ExpectExactResume(1, true);}

TEST(SolidState, FileLayoutIsWhatTheDocSays)
{
    const Saved& s=SavedState(0, false);
    auto o=H5::File::Open(s.path);
    ASSERT_TRUE(o) << o.Error();
    H5::File f=o.TakeValue();
    EXPECT_EQ(f.AttrString("format"), "qchem-solid-state 1");
    EXPECT_EQ(f.AttrInt("converged"), 1);
    EXPECT_NEAR(f.AttrReal("energy"), s.energy, 1e-12);
    H5::Group blocks=f.OpenGroup("blocks");
    ASSERT_EQ(blocks.Children().size(), 8u);                        // 2x2x2, unpolarized, nothing folded
    double wsum=0.0;
    for (size_t b=0;b<8;b++)
    {
        H5::Group g=blocks.OpenGroup(std::to_string(b));
        const size_t n=size_t(g.AttrInt("n")), nmo=size_t(g.AttrInt("nmo"));
        EXPECT_EQ(g.AttrInt("real"), 1) << "every 2x2x2 block is TRIM, so real by default";
        EXPECT_FALSE(g.IsComplex("D"));
        const std::vector<double> D=g.ReadReal("D"), C=g.ReadReal("C"), fv=g.ReadReal("f");
        ASSERT_EQ(C.size(), n*nmo);
        double ne=0.0;
        for (double x : fv) ne+=x;
        EXPECT_NEAR(ne, 8.0, 1e-8) << "block " << b << ": f is the PHYSICAL occupation (unweighted)";
        // D == C f C^T, element by element (the writer's own definition, checked from the file side).
        double dmax=0.0;
        for (size_t a=0;a<n;a++)
            for (size_t c=0;c<n;c++)
            {
                double x=0.0;
                for (size_t i=0;i<nmo;i++) x+=C[a*nmo+i]*fv[i]*C[c*nmo+i];
                dmax=std::max(dmax, std::fabs(x-D[a*n+c]));
            }
        EXPECT_LT(dmax, 1e-12);
        wsum+=g.AttrReal("weight");
    }
    EXPECT_NEAR(wsum, 1.0, 1e-12);
    EXPECT_TRUE(f.OpenGroup("fingerprint").OpenGroup("refuse").Has("basis.exponents"));
}

TEST(SolidState, WarmClassDifferencesRestartAndAreNamed)
{
    auto st=ReadSolidState(SavedState(0, false).path);
    ASSERT_TRUE(st) << st.Error().details;
    StateFingerprint now=st->fingerprint;
    auto same=CompareFingerprints(st->fingerprint, now);
    ASSERT_TRUE(same);
    EXPECT_TRUE(same->empty()) << "a state compared with itself is an EXACT resume";
    now.warmNum["scf.kT"]={0.01};
    now.warmNum["hubbard.U"]={0.147};
    auto warm=CompareFingerprints(st->fingerprint, now);
    ASSERT_TRUE(warm);
    ASSERT_EQ(warm->size(), 2u);
    EXPECT_NE((*warm)[0].find("hubbard.U"), std::string::npos);
    EXPECT_NE((*warm)[1].find("scf.kT"), std::string::npos);
    now.refuseNum["structure.positions"][0]+=0.01;
    now.refuseStr["functional"]="PBE";
    auto refused=CompareFingerprints(st->fingerprint, now);
    ASSERT_FALSE(refused);
    EXPECT_EQ(refused.Error().why, RestartRefusal::Why::Mismatch);
    EXPECT_NE(refused.Error().details.find("structure.positions"), std::string::npos);
    EXPECT_NE(refused.Error().details.find("functional"), std::string::npos);
}

TEST(SolidState, WarmStart_OtherGrid_StartsAtTheAnswer)
{
    // Grid continuation across processes: the same orbital basis on a finer density grid.  THE CLAIM IS WHERE THE
    // SCF STARTS, not how many iterations it takes: measured 2026-09-28, the warm start's FIRST iterate is 2e-10
    // Ha from the converged Ecut=30 energy (a fresh IonicSAD start: 6.5e-4), yet it still takes 9 iterations to
    // the fresh run's 8 -- the recipe's near-convergence tail (rho_mix halves to 0.5/0.4 and DIIS resets on
    // 1e-12 energy noise), a property of the accelerator, not of the restart.
    const Saved& s=SavedState(0, false);
    SiRun si=Si(ivec3_t(2,2,2));
    SolidCalcOptions o=Options(0, false, "");
    o.densityEcut=30.0;
    std::vector<double> E;
    o.onIteration=[&E](const SCFIterator::SCFProgress& p){ E.push_back(p.energy); };
    auto c=SolidCalculation::Restart(s.path, *si.lat, si.mol, o, Params());
    ASSERT_TRUE(c) << c.Error().details;
    auto r=(*c)->Result();
    ASSERT_TRUE(r);
    ASSERT_FALSE(E.empty());
    EXPECT_NEAR(r->Energy(), s.energy, 1e-3) << "a finer grid moves E a little, not a lot";
    EXPECT_NEAR(E.front(), r->Energy(), 1e-8) << "the first iterate on the new grid is already the answer";
}

TEST(SolidState, RealnessCrossesBothWays)
{
    // Complex-saved -> real run (D at a TRIM point is real) and real-saved -> forced-complex run: both WARM
    // (blocks.real differs), both land on the saved energy.
    for (bool savedComplex : {true, false})
    {
        const Saved& s=SavedState(savedComplex ? 1 : 0, savedComplex);
        SiRun si=Si(ivec3_t(2,2,2));
        auto c=SolidCalculation::Restart(s.path, *si.lat, si.mol, Options(savedComplex ? 1 : 0, !savedComplex, ""), Params());
        ASSERT_TRUE(c) << c.Error().details;
        auto r=(*c)->Result();
        ASSERT_TRUE(r);
        EXPECT_NEAR(r->Energy(), s.energy, 1e-8);
        EXPECT_LE(r->IterationCount(), 3u);
    }
}

TEST(SolidState, RefusedOnAnotherKMesh_OrSpinGroup_OrMissingFile)
{
    const Saved& s=SavedState(0, false);
    {
        SiRun si=Si(ivec3_t(1,1,1));
        auto c=SolidCalculation::Restart(s.path, *si.lat, si.mol, Options(0, false, ""), Params());
        ASSERT_FALSE(c);
        EXPECT_EQ(c.Error().why, RestartRefusal::Why::Mismatch);
        EXPECT_NE(c.Error().details.find("blocks.k"), std::string::npos);
    }
    {
        SiRun si=Si(ivec3_t(2,2,2));
        auto c=SolidCalculation::Restart(s.path, *si.lat, si.mol, Options(1, false, ""), Params());
        ASSERT_FALSE(c);
        EXPECT_NE(c.Error().details.find("spinGroup"), std::string::npos);
    }
    {
        SiRun si=Si(ivec3_t(2,2,2));
        auto c=SolidCalculation::Restart(TempPath("ck1_no_such_state.h5"), *si.lat, si.mol, Options(0, false, ""), Params());
        ASSERT_FALSE(c);
        EXPECT_EQ(c.Error().why, RestartRefusal::Why::Unreadable);
    }
}

} //namespace
