// File: IntegrationTests/GPW/Li.C  GPW bcc Li (q1, imposed Oh, Becke XC, smeared metal): the THREADED multi-k run is a run, not a guess.
//
// OpenWork_History5 "THREADED multi-k run HANGS FOREVER" (closed 2026-11-02, lane A): before the delta fitter's per-block caches locked, a
// threaded multi-k Becke run either spun forever in a corrupted std::map or -- what this gate is for -- CONVERGED to a silently WRONG
// energy: lane H's Li_bcc_q1_k4_u at 8 threads read Etot -0.272725 vs -0.252040 serial (-20.7 mHa, almost all E_xc, no warning).  The
// Si k311 Fock bit-identity gates (UTResponse) passed through that defect, so this one runs the real thing: the deck, to convergence.
#include "gtest/gtest.h"
#include <unistd.h>    // getpid: per-PROCESS temp names (ctest runs each case as its own process, in parallel)
#include <cstdio>
#include <filesystem>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <nlohmann/json.hpp>
import qchem.Deck;
import qchem.Parallel;     // ScopedWorkerThreads -- one process runs the same deck at two thread counts
import qchem.Tests.Panel;  // panel::BasisDataDir (the in-tree VALENCE_LOWQ_VA file)
import qchem.Reporting;    // the PP / functional policy rows the Li box runs file

namespace
{
namespace fs=std::filesystem;
// lane H's ~/Code/qchem6-runs/Li/mu/Li_bcc_q1_k4_u.json, inline: bcc a = 6.52145 bohr, Li q1, k444 (8 irreducible of 64), stated x 8 rasters.
nlohmann::json LiDeck()
{
    const double h=3.260725;
    return {
        {"structure",{{"name","Li_bcc"},{"lattice",{{"vectors_bohr",{{-h,h,h},{h,-h,h},{h,h,-h}}}}},{"atoms",{{{"el","Li"},{"frac",{0.0,0.0,0.0}}}}}}},
        {"kmesh",{4,4,4}},
        {"basis",{{"data","VALENCE_LOWQ_VA"},{"file",(qchem::panel::BasisDataDir()/"valence_lowq_va.bsd").string()},{"autoTrim",false}}},
        {"solid",{{"imposeSymmetry",true},{"xcMesh",{{"cellKind","Becke"}}},{"species",nlohmann::json::array({nlohmann::json::array({"Li",1})})},
                  {"ortho","CholeskyPivoted"},{"orthoTol",1e-7},{"tolerances",{{"rasterX",8.0},{"vlocRasterX",8.0},{"xcRasterX",8.0}}}}},
        {"accuracy",0},
        {"scf",{{"NMaxIter",100},{"minDeltaRhoPerElectron",1e-7},{"smearingkT",0.001}}}};
}
double RunAt(int workers, const std::string& tag)
{
    qchem::ScopedWorkerThreads scope(workers);
    qchem::deck::RunSpec spec; qchem::deck::FromJson(LiDeck(), spec);
    const fs::path dir=fs::temp_directory_path()/("qchem_it_li_thr_"+tag+"_"+std::to_string(::getpid()));
    fs::remove_all(dir);
    qchem::deck::Provenance pv; pv.codeVersion="it";
    const auto out=qchem::deck::Run(spec, pv, dir);
    fs::remove_all(dir);
    EXPECT_TRUE(out.converged) << out.summary;
    return out.energy.value_or(0.0);
}
} // namespace

// Measured 2026-11-02 (lane A, rundeck): 1 / 2 / 8 / 8 threads -0.2520403163 / -0.2520403162 / -0.2520403180 / -0.2520403180, 10 iterations
// each.  The SAME thread count is bit-identical (the run-to-run fix holds); a DIFFERENT count is a different collocation partition, so it
// legitimately moves the last digits (OpenWork row 160's thread-count remainder) -- 1e-8 holds that and fails the 20.7 mHa defect 2e6 times over.
TEST(GPW_Li, k444_Becke_Imp_Smear_ThreadedRunIsReproducibleAndMatchesSerial_Long)
{
    const double e8a=RunAt(8, "8a"), e8b=RunAt(8, "8b"), e1=RunAt(1, "1");
    EXPECT_EQ(e8a, e8b) << "two runs at the SAME thread count must be bit-identical";
    EXPECT_NEAR(e8a, e1, 1e-8) << "threaded vs serial: only the partition's rounding (row 160) may separate them";
    EXPECT_NEAR(e1, -0.2520403163, 1e-8) << "the serial bcc Li q1 k444 state (lane H's serial run reads the same)";
}

// ---------------------------------------------------------------------------------------------------------------------------------------
// GPW_LiBox: the Li ATOM with GTH-PADE-q3 (1s² 2s¹ in the valence) in a 12 bohr box -- the only SCF gate on a C3/C4 local PP.
//
// OpenWork §3 "Li GTH-PADE-q3 is NON-VARIATIONAL in GPW" (closed 2026-11-04, lane H; doc/Records/HubbardUHistory.md): Hermite3 sized its
// McMurchie-Davidson table by a fixed LMAX=3, and the C3/C4 short-range local terms (r⁴, r⁶ ⇒ per-axis powers 4 and 6) indexed past it.
// The analytic V_short matrix was 38% wrong and the atom fell 9-30 mHa BELOW CP2K ATOM's complete-basis value.  No suite test used a
// C3/C4 GTH (only Li q3 and Be q4 carry them), so a 1263-test sweep never saw it; GPW_UT LocalPPShortAnalyticMatchesGrid_C3C4 guards the
// operator, this guards the energy.
//
// THE ORACLES (2026-11-04; the like-for-like one re-run 2026-11-11 with OUR functional, lane H):
//   CP2K ATOM, complete basis, PADE   -7.3007906   a VARIATIONAL floor for any finite basis: falling below it is the defect's
//                                                  signature (PADE, not VWN5: the functionals differ by ~0.6 mHa here, the margin is 15 mHa).
//   CP2K QS GPW, this basis, L=12,    -7.2857130   XALPHA 2/3 + VWN5, i.e. OUR LDA (Dirac + VWN5); EPS_DEFAULT 1e-14, CUTOFF 2000/150.
//     Dirac + VWN5
//   ours, default grids               -7.2857262   13 μHa below it (grid-converged Becke 100/35: -7.2857129, 0.1 μHa).
// The 0.64 mHa "gap" this gate once sat above was the FUNCTIONAL, not the code: the old oracle ran CP2K's XC_FUNCTIONAL PADE (the
// Goedecker-Teter-Hutter Padé LDA the GTH-PADE PPs are fitted with), -7.2850821, unchanged to 2e-10 at EPS_PGF_ORB 1e-30 / CUTOFF
// 3000/300.  With VWN5 on both sides q1 on this basis agrees to 15 μHa too (-0.1840431 vs CP2K -0.1840584).  OpenWork §3 row closed
// 2026-11-11; decks + outputs ~/Code/cp2k-runs/Li/atom/qs/li_q3_qs_L12{,_tight,_vwn5}.{inp,out}.
namespace
{
//! The Li q3 box deck under \a model ("LDA" | "PADE"): the same basis, box and grids for both oracles.
nlohmann::json LiBoxDeck(const std::string& model)
{
    return {
        {"structure",{{"name","Li_atom"},{"lattice",{{"vectors_bohr",{{12.0,0.0,0.0},{0.0,12.0,0.0},{0.0,0.0,12.0}}}}},
                      {"atoms",{{{"el","Li"},{"frac",{0.0,0.0,0.0}}}}}}},
        {"kmesh",{1,1,1}},
        {"basis",{{"data","VALENCE_LOWQ_VA"},{"file",(qchem::panel::BasisDataDir()/"valence_lowq_va_liq3.bsd").string()},{"autoTrim",true}}},
        {"solid",{{"species",nlohmann::json::array({nlohmann::json::array({"Li",3})})},{"multiplicity",2},{"model",model},
                  {"ortho","CholeskyPivoted"},{"orthoTol",1e-7}}},
        {"scf",{{"NMaxIter",100},{"minDeltaRhoPerElectron",1e-7}}}};
}
//! The PPFittedFunctional rows of the LAST filed run whose name starts with \a prefix (the deck files its run under its stem).
size_t PPFunctionalRows(const std::string& prefix)
{
    const auto& G=qchem::report::GlobalReport();
    const qchem::report::json* last=nullptr;
    for (auto it=G.begin(); it!=G.end(); ++it) if (it.key().rfind(prefix,0)==0) last=&it.value();
    if (!last) throw std::logic_error("no filed run "+prefix);
    size_t n=0;
    std::function<void(const qchem::report::json&)> walk=[&](const qchem::report::json& x)
    {
        if (x.is_object())
            for (auto it=x.begin(); it!=x.end(); ++it)
            {
                if (it.key()=="policy" && it.value().is_array())
                    for (const auto& r : it.value()) if (r.contains("name") && r["name"]=="PPFittedFunctional") n++;
                walk(it.value());
            }
        else if (x.is_array()) for (const auto& y : x) walk(y);
    };
    walk(*last);
    return n;
}
//! Run the Li q3 box under \a model; returns its outcome and how many PP / functional warnings its run filed.
std::pair<qchem::deck::RunOutcome,size_t> RunLiBox(const std::string& model)
{
    qchem::deck::RunSpec spec; qchem::deck::FromJson(LiBoxDeck(model), spec);
    const fs::path dir=fs::temp_directory_path()/("qchem_it_libox_q3_"+model+"_"+std::to_string(::getpid()));
    fs::remove_all(dir);
    qchem::deck::Provenance pv; pv.codeVersion="it";
    auto out=qchem::deck::Run(spec, pv, dir);
    fs::remove_all(dir);
    return {out, PPFunctionalRows("Li_atom")};
}
} // namespace

TEST(GPW_LiBox, Γ_M2_CP2K_Long)
{
    const auto [out, warned]=RunLiBox("LDA");
    ASSERT_TRUE(out.converged) << out.summary;
    const double e=out.energy.value_or(0.0);
    EXPECT_GT(e, -7.3007906) << "below CP2K ATOM's complete-basis energy: the q3 local PP is non-variational again (the Hermite3 overrun class)";
    EXPECT_NEAR(e, -7.2857130, 5e-5) << "CP2K QS GPW on the same basis + box + functional (Dirac + VWN5); default grids sit 13 μHa off";
    // The policy fires once per Hamiltonian BUILD (the solid schedule builds one per SCF stage: two here), so the claim is that the
    // mismatch is REPORTED, not how many builds the schedule makes (no count pins).
    EXPECT_GE(warned, 1u) << "a GTH-PADE PP under LDA (Dirac + VWN5) files a PPFittedFunctional warn row (user ruling 2026-11-12)";
}

// THE PADE TWIN (user ruling 2026-11-12, OpenWork PADE row): the functional the GTH-PADE PPs were fitted with, stated, against the OLD
// oracle -- CP2K QS GPW, this basis, L=12, XC_FUNCTIONAL PADE: -7.2850821 (unchanged to 2e-10 at EPS_PGF_ORB 1e-30 / CUTOFF 3000/300;
// ~/Code/cp2k-runs/Li/atom/qs/li_q3_qs_L12{,_tight}.out).  Same grids as the VWN5 gate, so the same ~13 μHa grid offset is expected.
// And the PP now MATCHES the run's functional, so the run files no PPFittedFunctional row.
TEST(GPW_LiBox, Γ_PADE_M2_CP2K_Long)
{
    const auto [out, warned]=RunLiBox("PADE");
    ASSERT_TRUE(out.converged) << out.summary;
    const double e=out.energy.value_or(0.0);
    std::printf("  [Li q3 PADE] E %.10f  CP2K PADE -7.2850821  diff %+.2e Ha\n", e, e+7.2850821);
    EXPECT_GT(e, -7.3007906) << "below CP2K ATOM's complete-basis PADE energy: non-variational";
    EXPECT_NEAR(e, -7.2850821, 5e-5) << "CP2K QS GPW, XC_FUNCTIONAL PADE, on the same basis + box";
    EXPECT_EQ(warned, 0u) << "a GTH-PADE PP under PADE matches: no PPFittedFunctional row";
}
