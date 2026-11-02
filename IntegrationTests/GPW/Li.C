// File: IntegrationTests/GPW/Li.C  GPW bcc Li (q1, imposed Oh, Becke XC, smeared metal): the THREADED multi-k run is a run, not a guess.
//
// OpenWork_History5 "THREADED multi-k run HANGS FOREVER" (closed 2026-11-02, lane A): before the delta fitter's per-block caches locked, a
// threaded multi-k Becke run either spun forever in a corrupted std::map or -- what this gate is for -- CONVERGED to a silently WRONG
// energy: lane H's Li_bcc_q1_k4_u at 8 threads read Etot -0.272725 vs -0.252040 serial (-20.7 mHa, almost all E_xc, no warning).  The
// Si k311 Fock bit-identity gates (UTResponse) passed through that defect, so this one runs the real thing: the deck, to convergence.
#include "gtest/gtest.h"
#include <unistd.h>    // getpid: per-PROCESS temp names (ctest runs each case as its own process, in parallel)
#include <filesystem>
#include <string>
#include <nlohmann/json.hpp>
import qchem.Deck;
import qchem.Parallel;     // ScopedWorkerThreads -- one process runs the same deck at two thread counts
import qchem.Tests.Panel;  // panel::BasisDataDir (the in-tree VALENCE_LOWQ_VA file)

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
        {"basis",{{"data","VALENCE_LOWQ_VA"},{"file",(qchem::panel::BasisDataDir()/"valence_lowq_va.bsd").string()},{"vet",false}}},
        {"solid",{{"imposeSymmetry",true},{"xcMesh",{{"cellKind","Becke"}}},{"species",nlohmann::json::array({nlohmann::json::array({"Li",1})})},
                  {"ortho","CholeskyPivoted"},{"orthoTol",1e-7},{"accuracy",0},{"tolerances",{{"rasterX",8.0},{"vlocRasterX",8.0},{"xcRasterX",8.0}}}}},
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
