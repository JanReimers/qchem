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
#include <unistd.h>    // getpid: per-PROCESS temp names (ctest runs each case as its own process, in parallel)
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
import qchem.SCFCap;                             // IterationBudget: the convergence PANEL's budget (doc/SCFRobustnessPlan.md §4)
import qchem.SolidState;
import qchem.Lattice_3D;
import qchem.Structure;
import qchem.UnitCell;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.SCFParams;
import qchem.SCFAccelerator.Factory;   // SCFAccelerators::Type::Ladder (the restart-rung gate)
import qchem.SCFIterator;   // SCFProgress (the first-iterate probe)
import qchem.HDF5;
import qchem.Types;
import qchem.Mesh;   // qcMesh::UnitCellKind (the field-mixer restart gate's XC mesh)

using namespace qchem;
namespace {

//! Per PROCESS (ctest runs each case as its own process, in parallel, and the cached saved states are rebuilt in each).
std::string TempPath(const std::string& leaf) {return (std::filesystem::path(testing::TempDir())/(std::to_string(::getpid())+"_"+leaf)).string();}

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
    p.NMaxIter=80; p.MinΔρPerElectron=1e-7; p.MinΔE=1e-10; p.MinΔFD=1e-7; p.MinVirial=1e30; p.MinFD=1e30;
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
    SolidCalculation calc(SolidSetup::Resolve(*si.lat, si.mol, Options(multiplicity, forceComplex, s.path)), Params());
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
    auto c=SolidCalculation::Restart(s.path, SolidSetup::Resolve(*si.lat, si.mol, Options(multiplicity, forceComplex, "")), Params());
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
    EXPECT_EQ(f.AttrString("format"), "qchem-solid-state 1.1");
    EXPECT_EQ(f.AttrInt("converged"), 1);
    EXPECT_TRUE(std::isfinite(f.AttrReal("densityChange"))) << "1.1: the residual the state was saved at";
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
    now.refuseNum["structure.cell"][0]+=0.01;
    now.refuseStr["pp.set"]="some-other-family";
    auto refused=CompareFingerprints(st->fingerprint, now);
    ASSERT_FALSE(refused);
    EXPECT_EQ(refused.Error().why, RestartRefusal::Why::Mismatch);
    EXPECT_NE(refused.Error().details.find("structure.cell"), std::string::npos);
    EXPECT_NE(refused.Error().details.find("pp.set"), std::string::npos);
}

// C4 step 5: the functional and the atom positions are WARM (an LDA->PBE restart, a geometry step), and a moved atom is
// reported as ONE line that says how far.
TEST(SolidState, AFunctionalChangeAndAMovedAtomAreWarmAndSayHowFar)
{
    auto st=ReadSolidState(SavedState(0, false).path);
    ASSERT_TRUE(st) << st.Error().details;
    StateFingerprint now=st->fingerprint;
    now.warmStr["functional"]="PBE";
    now.warmNum["structure.positions"][0]+=0.03;
    now.warmNum["structure.positions"][1]+=0.04;
    now.warmNum["basis.shellCenters"][0]+=0.03;
    auto warm=CompareFingerprints(st->fingerprint, now);
    ASSERT_TRUE(warm) << warm.Error().details;
    ASSERT_EQ(warm->size(), 2u);
    EXPECT_NE((*warm)[0].find("functional"), std::string::npos);
    EXPECT_NE((*warm)[1].find("max displacement 0.05"), std::string::npos) << (*warm)[1];
}

// ...and a state written under the OLD schema (functional, positions and shell centres filed under REFUSE, no pp.set) still loads.
TEST(SolidState, AStateFromBeforeStep5StillLoads)
{
    auto st=ReadSolidState(SavedState(0, false).path);
    ASSERT_TRUE(st) << st.Error().details;
    StateFingerprint old=st->fingerprint;
    for (const char* k : {"structure.positions","basis.shellCenters"}) { old.refuseNum[k]=old.warmNum.at(k); old.warmNum.erase(k); }
    old.refuseStr["functional"]=old.warmStr.at("functional"); old.warmStr.erase("functional");
    old.refuseStr.erase("pp.set");
    auto same=CompareFingerprints(old, st->fingerprint);
    ASSERT_TRUE(same) << same.Error().details;
    EXPECT_TRUE(same->empty()) << "an old file against the same run is an EXACT resume";
}

// Gap 2 (2026-10-30): the numerics that move the ground state are recorded WARM -- and a state from before them is still an
// exact resume (none of them compared), while a state that HAS them is held to them.
TEST(SolidState, TheNumericsKeysAreWarmAndAnOlderStateIsStillExact)
{
    auto st=ReadSolidState(SavedState(0, false).path);
    ASSERT_TRUE(st) << st.Error().details;
    const StateFingerprint& fp=st->fingerprint;
    EXPECT_EQ(fp.warmStr.at("xc.fit"), "PlaneWave") << "Auto on a uniform mesh resolves to the plane-wave fit";
    EXPECT_EQ(fp.warmStr.at("grids.raster"), "BallOnly");
    for (const char* k : {"grids.ladderFactor","gpw.vlocEps","gpw.screenEps","gpw.relCutoff","gpw.mgridEcuts"})
        EXPECT_TRUE(fp.warmNum.count(k)) << k;
    EXPECT_FALSE(fp.warmNum.count("grids.beckeEps")) << "no Becke mesh, no Becke eps";

    StateFingerprint old=fp;
    old.warmStr.erase("xc.fit"); old.warmStr.erase("grids.raster"); old.warmNum.erase("grids.ladderFactor");
    std::erase_if(old.warmNum, [](const auto& e){ return e.first.rfind("gpw.",0)==0; });
    auto same=CompareFingerprints(old, fp);
    ASSERT_TRUE(same) << same.Error().details;
    EXPECT_TRUE(same->empty()) << "a state written before the numerics keys is an EXACT resume of the same run";

    StateFingerprint now=fp;
    now.warmStr["xc.fit"]="Delta";
    now.warmNum["gpw.screenEps"]={1e-12};
    auto warm=CompareFingerprints(fp, now);
    ASSERT_TRUE(warm) << warm.Error().details;
    ASSERT_EQ(warm->size(), 2u);
    // Gap 3: each line says what the difference MEANS, not only the two values.
    EXPECT_NE((*warm)[0].find("gpw.screenEps"), std::string::npos);
    EXPECT_NE((*warm)[0].find("GPW numerical tolerance changed"), std::string::npos) << (*warm)[0];
    EXPECT_NE((*warm)[1].find("saved 'PlaneWave', now 'Delta'"), std::string::npos) << (*warm)[1];
    EXPECT_NE((*warm)[1].find("v_xc fit basis changed: the converged ground state moves"), std::string::npos) << (*warm)[1];
}

// A restart ACROSS one of the new keys (the multigrid ladder) is a warm start that says so, and converges to the new run's answer.
// (The raster was the first choice: AliasFree on this cell costs ~230 s, far over a unit test's budget.)
TEST(SolidState, ARestartAcrossTheLadderPrintsItsLineAndConverges)
{
    const Saved& s=SavedState(0, false);
    SiRun si=Si(ivec3_t(2,2,2));
    SolidCalcOptions o=Options(0, false, "");
    o.ladderFactor=3.0;
    testing::internal::CaptureStdout();
    auto c=SolidCalculation::Restart(s.path, SolidSetup::Resolve(*si.lat, si.mol, o), Params());
    const std::string out=testing::internal::GetCapturedStdout();
    ASSERT_TRUE(c) << c.Error().details;
    EXPECT_NE(out.find("WARM START"), std::string::npos);
    EXPECT_NE(out.find("differs: grids.ladderFactor: saved [4], now [3] -- the multigrid ladder changed"), std::string::npos) << out;
    EXPECT_NE(out.find("saved at Δρ"), std::string::npos) << "the saved state's own residual is reported";
    auto r=(*c)->Result();
    ASSERT_TRUE(r);
    EXPECT_NEAR(r->Energy(), s.energy, 1e-3) << "the ladder moves the answer only slightly";
}

//! The same state as the format-1.0 writer made it: every attribute, group and dataset EXCEPT the 1.1 fields, stamped
//! "qchem-solid-state 1".  Built by copying, so the test needs no 1.0 writer kept alive in the library.
void CopyAsFormat10(const H5::Group& from, H5::Group& to, bool root)
{
    for (const std::string& a : from.AttrNames())
    {
        if (root && a=="densityChange") continue;
        if (root && a=="format") { to.SetAttr("format", "qchem-solid-state 1"); continue; }
        if (from.AttrIsString(a)) { to.SetAttr(a, from.AttrString(a)); continue; }
        const double x=from.AttrReal(a);
        if (std::isfinite(x) && x==std::floor(x) && std::fabs(x)<1e15) to.SetAttr(a, std::int64_t(x)); else to.SetAttr(a, x);
    }
    for (const std::string& c : from.Children())
    {
        bool isGroup=true;   // the wrapper has no "is a group" ask: a dataset refuses OpenGroup
        try { (void)from.OpenGroup(c); } catch (...) { isGroup=false; }
        if (isGroup) { H5::Group g=to.CreateGroup(c); CopyAsFormat10(from.OpenGroup(c), g, false); continue; }
        if (from.IsComplex(c)) to.Write(c, from.ReadComplex(c), from.Shape(c));
        else                   to.Write(c, from.ReadReal(c),    from.Shape(c));
    }
}

// D-STATE-SCHEMA (CleanCode): a minor version only ADDS optional fields, so a state written BEFORE one (format 1.0, no
// densityChange) still reads and restores, and the restart says the field was not recorded.
TEST(SolidState, AFormat10StateStillRestores)
{
    const Saved& s=SavedState(0, false);
    const std::string old=TempPath("ck1_si_format10.h5");
    {
        auto o=H5::File::Open(s.path);
        ASSERT_TRUE(o) << o.Error();
        H5::File from=o.TakeValue();
        H5::File to=H5::File::Create(old);
        CopyAsFormat10(from, to, true);
    }
    auto st=ReadSolidState(old);
    ASSERT_TRUE(st) << st.Error().details;
    EXPECT_TRUE(std::isnan(st->summary.densityChange)) << "absent in 1.0 = not recorded";
    SiRun si=Si(ivec3_t(2,2,2));
    testing::internal::CaptureStdout();
    auto c=SolidCalculation::Restart(old, SolidSetup::Resolve(*si.lat, si.mol, Options(0, false, "")), Params());
    const std::string out=testing::internal::GetCapturedStdout();
    ASSERT_TRUE(c) << c.Error().details;
    EXPECT_NE(out.find("EXACT RESUME"), std::string::npos) << out;
    EXPECT_NE(out.find("saved residual not recorded (a format-1.0 state"), std::string::npos) << out;
    auto r=(*c)->Result();
    ASSERT_TRUE(r);
    EXPECT_NEAR(r->Energy(), s.energy, 1e-8);
    EXPECT_LE(r->IterationCount(), 3u);
    std::filesystem::remove(old);
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
    auto c=SolidCalculation::Restart(s.path, SolidSetup::Resolve(*si.lat, si.mol, o), Params(),
                                     [&E](const SCFIterator::SCFProgress& p){ E.push_back(p.energy); });
    ASSERT_TRUE(c) << c.Error().details;
    auto r=(*c)->Result();
    ASSERT_TRUE(r);
    ASSERT_FALSE(E.empty());
    EXPECT_NEAR(r->Energy(), s.energy, 1e-3) << "a finer grid moves E a little, not a lot";
    EXPECT_NEAR(E.front(), r->Energy(), 1e-8) << "the first iterate on the new grid is already the answer";
}

// A restart enters the Ladder at its GDM rung at kT = 0 only (user 2026-10-21; the NaF gates GPW_NaF.Γ_Imp_eq{ExactResume,
// ColdStart} pin the integer-occupation win).  A SMEARED restart from a kT = 0 state must start on DIIS: the restored D' is
// idempotent, so GDM would ENGAGE and minimise with integer occupations, never applying the smearing (measured 2026-10-21
// once the carried density became the working density: 6 mHa off the cold smeared run).  It must land where a cold smeared
// Ladder run does.
TEST(SolidState, LadderRestart_Smeared_StartsOnDIIS)
{
    const Saved& s=SavedState(0, false);                     // a kT=0 state: the restart is a WARM start (kT differs)
    SiRun si=Si(ivec3_t(2,2,2));
    SolidCalcOptions o=Options(0, false, "");
    o.accelerator=SCFAccelerators::Type::Ladder;
    SCFParams p=Params();
    p.SmearingkT=0.02;
    double cold=0.0;
    {
        SolidCalculation calc(SolidSetup::Resolve(*si.lat, si.mol, o), p);
        auto r=calc.Result();
        ASSERT_TRUE(r) << "the cold smeared Ladder reference must converge";
        cold=r->Energy();
    }
    testing::internal::CaptureStdout();
    auto c=SolidCalculation::Restart(s.path, SolidSetup::Resolve(*si.lat, si.mol, o), p);
    const std::string out=testing::internal::GetCapturedStdout();
    ASSERT_TRUE(c) << c.Error().details;
    auto r=(*c)->Result();
    ASSERT_TRUE(r) << "a smeared Ladder restart must still converge";
    EXPECT_EQ(out.find("starting on rung 1 (GDM)"), std::string::npos) << "a smeared restart must not start on the GDM rung";
    EXPECT_NEAR(r->Energy(), cold, 1e-6);
}

// A CARRIED density is the working density (OpenWork §3 "A carried density is re-filled at kT = 0", 2026-10-21).  A SMEARED
// converged state restarted EXACTLY (same kT) must start AT the saved energy: before the fix, iteration 0 re-filled the
// restored density at kT = 0 (the occupation policy's default) and iteration 1 started from that integer re-fill.
TEST(SolidState, SmearedExactResume_StartsAtTheSavedState)
{
    SiRun si=Si(ivec3_t(2,2,2));
    SCFParams p=Params();
    p.SmearingkT=0.02;                                       // Si's gap is ~0.02 Ha: the frontier is genuinely fractional
    const std::string path=TempPath("ck1_si_smeared.h5");
    double saved=0.0;
    {
        SolidCalculation calc(SolidSetup::Resolve(*si.lat, si.mol, Options(0, false, path)), p);
        auto r=calc.Result();
        ASSERT_TRUE(r) << "the smeared Si reference must converge";
        saved=r->Energy();
    }
    SolidCalcOptions o=Options(0, false, "");
    std::vector<double> E;
    auto c=SolidCalculation::Restart(path, SolidSetup::Resolve(*si.lat, si.mol, o), p,
                                     [&E](const SCFIterator::SCFProgress& pr){ E.push_back(pr.energy); });
    ASSERT_TRUE(c) << c.Error().details;
    auto r=(*c)->Result();
    ASSERT_TRUE(r);
    ASSERT_FALSE(E.empty());
    EXPECT_NEAR(E.front(), saved, 1e-7) << "the first iterate of an exact smeared resume is the saved state, not a kT = 0 re-fill";
    EXPECT_NEAR(r->Energy(), saved, 1e-8);
    EXPECT_LE(r->IterationCount(), 3u);
}

// A FIELD-MIXER restart starts at its own fixed point (OpenWork §3 "A CARRIED density's first iteration feeds XC the FITTED
// rho", 2026-10-28).  With the Null accelerator the Fock is built from the mixer's PRESENTED field; under a Kerker filter that
// field's first iterate had no DM source seated, so iteration 1's XC took the fitted rho and the exact resume re-opened off its
// fixed point (TiO2 k222: drho 6e-4, 15 iterations back).  Pulay with and without the Kerker filter: the first iterate is the
// saved energy, XC sees rho[D] from iteration 1, and the resume confirms in <= 3 iterations.  Polarized on the Becke mesh: the
// singles sampler reaches rho[D] only through the field's N4 correction, which a mixer formed only inside Mix -- so a field
// taken straight from a D (the seed, a re-seat) had none.
namespace {
void ExpectFieldMixerExactResume(double kerkerG0)
{
    SiRun si=Si(ivec3_t(2,2,2));
    SCFParams p=Params();
    p.Mixer=SCFParams::Mixing::Pulay; p.PulayDepth=8; p.PulayStart=2; p.KerkerG0=kerkerG0;
    const std::string path=TempPath("ck1_si_fieldmix_m1_"+std::to_string(int(kerkerG0*10))+".h5");
    SolidCalcOptions o=Options(1, true, path);         // POLARIZED: the (ρ,m) mixer, as on the TiO2 / NiO campaign decks
    o.accelerator=SCFAccelerators::Type::Null;
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Becke;     // the singles route: rho[D] only through the N4 correction
    double saved=0.0;
    {
        SolidCalculation calc(SolidSetup::Resolve(*si.lat, si.mol, o), p);
        auto r=calc.Result();
        ASSERT_TRUE(r) << "the Si field-mixer reference must converge";
        saved=r->Energy();
    }
    o.saveStateTo="";
    std::vector<SCFIterator::SCFProgress> it;
    qchem::IterationBudget budget(3, "field-mixer exact resume (regime 4)");   // the PANEL case: the resume only
    testing::internal::CaptureStdout();
    auto c=SolidCalculation::Restart(path, SolidSetup::Resolve(*si.lat, si.mol, o), p,
                                     [&it](const SCFIterator::SCFProgress& pr){ it.push_back(pr); });
    const std::string out=testing::internal::GetCapturedStdout();
    ASSERT_TRUE(c) << c.Error().details;
    auto r=(*c)->Result();
    ASSERT_TRUE(r);
    ASSERT_FALSE(it.empty());
    // R1.3: the exact resume is a HAND-OFF, printed in the one format: "[handoff] iteration 0: restored state -> Pulay history from iteration 1 (exact resume: ...)".
    EXPECT_NE(out.find("-> Pulay history from iteration 1 (exact resume"), std::string::npos)
        << "an exact resume starts AT the fixed point: no plain-Kerker priming (beta-MnO2 r015: 23 -> 4 iterations)";
    EXPECT_NE(it.front().eb.charge.xcRho, "fitted") << "iteration 1's XC must see rho[D] (pin 28), not the fitted field";
    EXPECT_NEAR(it.front().energy, saved, 1e-7) << "the first iterate of an exact resume is the saved state";
    EXPECT_NEAR(r->Energy(), saved, 1e-7) << "the reference stopped at drho/N 1e-7: its energy is that far from the fixed point";
    EXPECT_LE(r->IterationCount(), 3u);
}
} // namespace
TEST(SolidState, ExactResume_Pulay_Null_StartsAtTheSavedState_Panel)       {ExpectFieldMixerExactResume(0.0);}   // PANEL case 3
TEST(SolidState, ExactResume_KerkerPulay_Null_StartsAtTheSavedState_Panel) {ExpectFieldMixerExactResume(1.0);}

TEST(SolidState, RealnessCrossesBothWays)
{
    // Complex-saved -> real run (D at a TRIM point is real) and real-saved -> forced-complex run: both WARM
    // (blocks.real differs), both land on the saved energy.
    for (bool savedComplex : {true, false})
    {
        const Saved& s=SavedState(savedComplex ? 1 : 0, savedComplex);
        SiRun si=Si(ivec3_t(2,2,2));
        auto c=SolidCalculation::Restart(s.path, SolidSetup::Resolve(*si.lat, si.mol, Options(savedComplex ? 1 : 0, !savedComplex, "")), Params());
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
        auto c=SolidCalculation::Restart(s.path, SolidSetup::Resolve(*si.lat, si.mol, Options(0, false, "")), Params());
        ASSERT_FALSE(c);
        EXPECT_EQ(c.Error().why, RestartRefusal::Why::Mismatch);
        EXPECT_NE(c.Error().details.find("blocks.k"), std::string::npos);
    }
    {
        SiRun si=Si(ivec3_t(2,2,2));
        auto c=SolidCalculation::Restart(s.path, SolidSetup::Resolve(*si.lat, si.mol, Options(1, false, "")), Params());
        ASSERT_FALSE(c);
        EXPECT_NE(c.Error().details.find("spinGroup"), std::string::npos);
    }
    {
        SiRun si=Si(ivec3_t(2,2,2));
        auto c=SolidCalculation::Restart(TempPath("ck1_no_such_state.h5"), SolidSetup::Resolve(*si.lat, si.mol, Options(0, false, "")), Params());
        ASSERT_FALSE(c);
        EXPECT_EQ(c.Error().why, RestartRefusal::Why::Unreadable);
    }
}

//! A PLANE-WAVE state (OpenWork §3 "PW states cannot be saved (or restarted)", 2026-10-30): saved through the same fingerprint and
//! block path as a Gaussian one, the cutoff standing in for the shells.  \a numBands > 0 is the partial eigensolve a production
//! PW run uses, whose restored stage starts on DIIS (GDM needs every virtual).
SolidCalcOptions PWOptions(size_t numBands, const std::string& save)
{
    SolidCalcOptions o=Options(0, false, save);
    o.orbitalBasis=SolidCalcOptions::OrbitalBasis::PlaneWave; o.planeWaveEcut=4.0; o.numBands=numBands; o.label="ck1pw";
    return o;
}
void ExpectPlaneWaveExactResume(size_t numBands)
{
    const std::string path=TempPath("ck1_si_pw_"+std::to_string(numBands)+".h5");
    SiRun si=Si(ivec3_t(1,1,1));
    SolidCalculation calc(SolidSetup::Resolve(*si.lat, si.mol, PWOptions(numBands, path)), Params());
    auto r=calc.Result();
    ASSERT_TRUE(r) << "the Si plane-wave reference run must converge";
    ASSERT_TRUE(std::filesystem::exists(path)) << "a plane-wave state must be SAVED";

    auto c=SolidCalculation::Restart(path, SolidSetup::Resolve(*si.lat, si.mol, PWOptions(numBands, "")), Params());
    ASSERT_TRUE(c) << c.Error().details;
    auto rr=(*c)->Result();
    ASSERT_TRUE(rr);
    // Full eigensolve: the restored stage starts on the Ladder's GDM rung, AT the saved state.  Partial: it starts on DIIS, which
    // takes one real step from an iterate the saved run stopped at Δρ/N 1e-7 -- measured 1.2e-8 Ha (2026-10-30), inside that stop.
    EXPECT_NEAR(rr->Energy(), r->Energy(), numBands==0 ? 1e-8 : 5e-8);
    EXPECT_LE(rr->IterationCount(), 3u) << "the reference run took " << r->IterationCount();

    // The family is REFUSE: a plane-wave state cannot seed a Gaussian run, nor the other way round.
    auto g=SolidCalculation::Restart(path, SolidSetup::Resolve(*si.lat, si.mol, Options(0, false, "")), Params());
    ASSERT_FALSE(g);
    EXPECT_NE(g.Error().details.find("basis.family"), std::string::npos) << g.Error().details;
}
TEST(SolidState, ExactResume_PlaneWave)                 {ExpectPlaneWaveExactResume(0);}
TEST(SolidState, ExactResume_PlaneWave_PartialEigensolve) {ExpectPlaneWaveExactResume(8);}

} //namespace
