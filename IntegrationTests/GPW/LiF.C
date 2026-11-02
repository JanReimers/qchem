// File: IntegrationTests/GPW/LiF.C  GPW LiF (rocksalt, a = 7.6 bohr): the Li-site curvature's XC-MESH error (doc/OpenWork.md §3 "A curvature
// needs an XC-MESH ladder as well as the egg-box").  An atom-centred Becke mesh moves WITH the atoms, so the egg-box (a rigid translation) is
// blind to it; a ONE-ion displacement moves the Becke partition between dissimilar neighbours, and what that costs a curvature is read here
// at a FROZEN D: each geometry's converged D re-evaluated on a finer reference mesh (the postSCF `xcMeshAudit`; a restart with NMaxIter 1 is
// NOT a frozen D -- it rebuilds and re-diagonalises).
#include "gtest/gtest.h"
#include <unistd.h>    // getpid: per-PROCESS temp names (ctest runs each case as its own process, in parallel)
#include <cmath>
#include <fstream>
#include <string>
#include <filesystem>
#include <nlohmann/json.hpp>
import qchem.Deck;     // deck::Run -- the gate is run FROM A DECK, so it also gates the xcMeshAudit action

namespace
{
namespace fs=std::filesystem;
// Li q1 (lane B's light-PP config) at Γ, the plain atom-centred mesh (a displacement path needs it), cellKind Becke with NO recipe fields =
// BeckeXCParams() (40, 2.0, Lebedev 29) since 2026-11-01.  The vet-stage trim (pin 22) is STATED so both geometries carry the same basis.
nlohmann::json LiFDeck(double dx)
{
    const double a=3.8, f=dx/(2.0*a);   // frac (-f, f, f) of the fcc vectors (0,a,a),(a,0,a),(a,a,0) = Cartesian (dx, 0, 0)
    return {
        {"structure",{{"name","LiF_xcaudit"},{"lattice",{{"vectors_bohr",{{0.0,a,a},{a,0.0,a},{a,a,0.0}}}}},
                      {"atoms",{{{"el","Li"},{"frac",{-f,f,f}}},{{"el","F"},{"frac",{0.5,0.5,0.5}}}}}}},
        {"kmesh",{1,1,1}}, {"basis",{{"data","VALENCE_LOWQ_VA"},{"trim",{{{"Z",3},{"l",0},{"alpha",0.03}}}}}},
        {"solid",{{"seed","IonicSAD"},{"imposeSymmetry",false},{"xcMeshSymmetric",false},{"species",nlohmann::json::array({nlohmann::json::array({"Li",1}),nlohmann::json::array({"F",7})})},
                  {"xcMesh",{{"cellKind","Becke"}}},{"ortho","CholeskyPivoted"},{"orthoTol",1e-4}}},
        {"scf",{{"NMaxIter",100},{"minDeltaRhoPerElectron",1e-7},{"minDeltaE",1e-10}}},
        {"postSCF",{{{"xcMeshAudit",nlohmann::json::object()}}}}};   // default reference (80, 48)
}
struct Audited { double energy, dExc, dHxc; };
Audited RunAudited(double dx, const std::string& tag)
{
    qchem::deck::RunSpec spec; qchem::deck::FromJson(LiFDeck(dx), spec);
    const fs::path dir=fs::temp_directory_path()/("qchem_it_lif_xcaudit_"+tag+"_"+std::to_string(::getpid()));
    fs::remove_all(dir);
    qchem::deck::Provenance pv; pv.codeVersion="it";
    const auto out=qchem::deck::Run(spec, pv, dir);
    EXPECT_TRUE(out.converged) << out.summary;
    std::ifstream in(out.revision); const nlohmann::json rec=nlohmann::json::parse(in);
    const auto& a=rec["results"]["postSCF"][0];
    EXPECT_EQ(a["action"], "xcMeshAudit");
    EXPECT_NE(a["reference"].get<std::string>().find("nR=80"), std::string::npos) << a["reference"];
    Audited r{out.energy.value_or(0.0), 0.0, 0.0};
    for (const auto& t : a["terms"])
    {
        if (t["entries"]=="Exc") { r.dExc=t["dE"].get<double>(); r.dHxc=t["maxAbsDH"].get<double>(); }
        else EXPECT_EQ(t["maxAbsDH"].get<double>(), 0.0) << t["term"] << ": only the XC term reads the XC mesh";
    }
    EXPECT_NEAR(r.dExc, a["dE"].get<double>(), 1e-14) << "the audit's total is its XC term";
    fs::remove_all(dir);
    return r;
}
} // namespace

// Measured 2026-11-01 (lane A; ~30 s of SCF per geometry): Etot(0) -24.3245668693, Etot(0.1) -24.3245349809 -> the curvature dE(0.1) = +3.19e-5 Ha.
// dE_xc vs (80, 48): -1.330e-4 at 0, -1.134e-4 at 0.1, max|dH_xc| 2.1e-3 / 1.8e-3 -> the default mesh's error ON THE CURVATURE is +1.96e-5, 61% of
// it.  The calibrated default (an ENERGY calibration) cannot carry a Li-site curvature: a geometry run climbs the nested ladder instead
// (SolidCalculation's pilot, pin 31).  Lane B's k222 calibration (CP2K -3.78e-5): (20,24) -1.34e-4, (40,36) -2.08e-5, (60,48) -2.63e-5.
// Thresholds: each audit to 10% of its measured value, the curvature error to 15% -- a mesh change, or an audit that stops reading the
// frozen D, moves them by far more.
TEST(GPW_LiF, Γ_Becke_XCMeshAuditReadsTheCurvatureError_Long)
{
    const Audited e0=RunAudited(0.0, "0"), e1=RunAudited(0.1, "p10");
    const double curvature=e1.energy-e0.energy, meshError=e1.dExc-e0.dExc;
    EXPECT_NEAR(e0.dExc, -1.330e-4, 1.3e-5);
    EXPECT_NEAR(e1.dExc, -1.134e-4, 1.1e-5);
    EXPECT_GT(e0.dHxc, 0.0); EXPECT_GT(e1.dHxc, 0.0);
    EXPECT_NEAR(curvature, 3.19e-5, 2e-6) << "Li q1 at Γ: the 0.1 bohr step";
    EXPECT_NEAR(meshError, 1.96e-5, 3e-6) << "the default Becke mesh's error on the curvature, at frozen D";
}
