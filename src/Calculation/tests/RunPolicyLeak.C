// File: Calculation/tests/RunPolicyLeak.C  A run's density route is ITS OWN, not the last solid's (D-MUTABLE-STATICS-RUNPOLICY, 2026-11-23).
//
// THE LEAK.  Only SolidSetup::Resolve installs a run policy (SetRunPolicy, process-wide), and the orbitals build every iteration's density
// through ChargeDensity::DefaultRhoRoute(), which reads that process-wide policy.  So a molecule (or atom) run takes the density route the
// last SOLID in the process installed: after a cp2kCompat solid, a molecule runs the Direct route where its own default is the factored
// (pivoted-Cholesky) one.  A GUI session or a test binary that runs a CP2K-parity solid and then a molecule is exposed.
//
// KNOWN RED (root CMakeLists QCHEM_KNOWN_RED, OpenWork §3 "D-RUNPOLICY-LEAK") until group E of doc/Records/RunPolicyOnTheRun.md lands (the
// user's ruling 2026-11-23, E2': the setup builds a density maker and the orbitals call it), which makes the route the molecule setup's own.
#include "gtest/gtest.h"
#include "forward.H"
#include <functional>
#include <memory>
#include <string>
#include <variant>
#include <vector>

import qchem.MoleculeCalculation;
import qchem.SolidCalculation;
import qchem.Hamiltonian.Factory;                 // Hamiltonian::Model
import qchem.Materials;
import qchem.Lattice_3D;
import qchem.Structure;
import qchem.UnitCell;
import qchem.BasisSet.Gaussian.Point.Factory;     // BasisSetData
import qchem.SCFParams;
import qchem.Mesh;
import qchem.Types;
import qchem.RunPolicy;                           // ScopedRunPolicy: the test leaves the process policy as it found it
import qchem.CompositeCD;
import qchem.Reporting;                           // report::Run / GlobalReport: the setup's DensityRoute row
import qchem.ChargeDensity.Imp.IrrepCD;          // FactoredRho / FiniteIrrepCD (tests may cheat: abstract -> concrete)
using namespace qchem;

//! THE FRIEND of tComposite_CD (src/forward.H): which concrete block densities a converged composite holds.  The composite's face never
//! names its blocks' realization; a test asking WHICH ρ route was built has to look.
class CompositeDensityTests
{
public:
    //! The number of blocks of \a cd, and how many of them are the factored (pivoted-Cholesky) finite leaf.
    static std::pair<size_t,size_t> FactoredBlocks(const ChargeDensity::tDM_CD<double>& cd)
    {
        auto* comp=dynamic_cast<const ChargeDensity::tComposite_CD<double>*>(&cd);
        if (!comp) throw std::logic_error("CompositeDensityTests: a molecule's density is a composite");
        size_t n=0, factored=0;
        for (const auto& b : comp->itsBlocks)
            std::visit([&](const auto* d)
            {
                ++n;
                if (dynamic_cast<const ChargeDensity::FactoredRho<ChargeDensity::FiniteIrrepCD>*>(d)) ++factored;
            }, b.cd);
        return {n, factored};
    }
};

namespace
{
//! H₂ HF/DZVP, converged: how many of its density blocks the factored route built.
std::pair<size_t,size_t> MoleculeBlocks()
{
    Molecule m;
    m.Insert(Site{1, rvec3_t(0.0, 0.0,-0.7)});
    m.Insert(Site{1, rvec3_t(0.0, 0.0, 0.7)});
    SCFParams p; p.NMaxIter=40; p.MinΔρPerElectron=1e-6;
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(m, {.model=Hamiltonian::Model::HF, .scf=p, .molecule={.basis="dzvp"}})).TakeValue();
    auto proof=calc->Converge();
    if (!proof) {ADD_FAILURE()<<"H2 HF/DZVP did not converge"; return {0,0};}
    return CompositeDensityTests::FactoredBlocks(proof.Value().Density());
}
//! The DensityRoute policy rows of the last filed run named \a runName (as GGARoute.C reads its own rows).
std::vector<report::json> DensityRouteRows(const std::string& runName)
{
    const auto& G=report::GlobalReport();
    const report::json* last=nullptr;
    for (auto it=G.begin(); it!=G.end(); ++it) if (it.key().rfind(runName+"@",0)==0) last=&it.value();
    if (!last) throw std::logic_error("no filed run "+runName);
    std::vector<report::json> out;
    std::function<void(const report::json&)> walk=[&](const report::json& x)
    {
        if (x.is_object())
            for (auto it=x.begin(); it!=x.end(); ++it)
            {
                if (it.key()=="policy" && it.value().is_array())
                    for (const auto& r : it.value()) if (r.contains("name") && r["name"]=="DensityRoute") out.push_back(r);
                walk(it.value());
            }
        else if (x.is_array()) for (const auto& y : x) walk(y);
    };
    walk(*last);
    return out;
}
} // namespace

// THE CONTROL first (no solid has run): the molecule's own default is the factored route, on every block.  THEN a cp2kCompat solid setup
// (its Resolve is what installs the policy -- no SCF needed), and the SAME molecule again: its route must not have changed.
TEST(RunPolicyLeak, AMoleculeAfterACP2KCompatSolidKeepsItsOwnDensityRoute)
{
    ScopedRunPolicy restore(RunPolicySpec{});   // whatever this test leaves installed, the next test in the process starts from the default
    const auto [n0, f0]=MoleculeBlocks();
    ASSERT_GT(n0, 0u);
    ASSERT_EQ(f0, n0) << "the control: with no solid run before it, every block of the molecule's density is the factored route";

    const auto mat=Materials::Get("Si_diamond");
    RunOptions o; o.Nelec=8; o.solid.densityEcut=20.0; o.label="cp2kCompat solid"; o.accuracy.eps=0.0;
    SolidXCMesh(o).cellKind=qcMesh::UnitCellKind::Uniform;
    o.solid.policy.cp2kCompat=true;
    const std::string run="runpolicyleak cp2kCompat solid";
    {
        report::Run r(run);
        auto solid=SolidSetup::Resolve(Lattice_3D(*mat.cell, ivec3_t(1,1,1)), GaussianBasisRecipe{.data=BasisSet::Gaussian::BasisSetData::SIPP_SR}, o);
        ASSERT_TRUE(solid) << "the cp2kCompat solid setup resolves";
    }
    // THE SETUP SAYS WHICH ROUTE IT BUILT: one DensityRoute row in the filed run's policy scope, direct, for cp2kCompat's reason.
    const auto rows=DensityRouteRows(run);
    ASSERT_EQ(rows.size(), 1u) << "one DensityRoute policy row in the cp2kCompat solid's setup";
    EXPECT_EQ(rows[0]["verdict"], "choose");
    EXPECT_NE(std::string(rows[0]["decided"]).find("direct"), std::string::npos) << std::string(rows[0]["decided"]);
    EXPECT_NE(std::string(rows[0]["reason"]).find("cp2kCompat"), std::string::npos) << std::string(rows[0]["reason"]);

    const auto [n1, f1]=MoleculeBlocks();
    ASSERT_EQ(n1, n0);
    EXPECT_EQ(f1, n1) << "the LEAK: a molecule run after a cp2kCompat solid took that solid's Direct density route (" << n1-f1 << " of "
                      << n1 << " blocks unfactored) -- its route must be its own setup's (D-MUTABLE-STATICS-RUNPOLICY group E)";
}
