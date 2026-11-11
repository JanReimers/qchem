// File: Calculation/tests/MeshAudit.C  The molecular XC-MESH AUDIT at a fixed D and the MESH PILOT that climbs on it
// (doc/GeometryOptimizationPlan.md §14.2-14.3; pin 31).  H2O LDA/dzvp and F2 GTH q7.  Invariants, not pinned numbers: the run against its
// OWN mesh reads zero everywhere; against a finer mesh only the mesh-borne terms move (XC; and Hartree once it rides the mesh,
// CoulombFit::Delta); the total is the sum of the terms; an audit before a Converge throws.  The pilot: an unstated mesh settles where
// its audit meets eps (water at the entry rung, F2 above it), a stated one is kept and reported, HF has none, and a session freezes it.
#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>
import qchem.MoleculeCalculation;
import qchem.StructureData;
import qchem.Structure;
import qchem.SCFParams;
import qchem.Types;
import qchem.Mesh;                // MeshParams (the reference mesh)
import qchem.Hamiltonian.Types;   // CoulombFit (Hartree on the mesh)

using namespace qchem;
using StructureData::GetMolecule;

namespace
{
SCFParams Tight()
{
    SCFParams p; p.NMaxIter=100; p.MinΔρPerElectron=1e-9;
    return p;
}
//! A mesh finer than \a m on both axes (the ladder's next radial count and an angular degree above it).
qcMesh::MeshParams Finer(qcMesh::MeshParams m) { m.nRadial+=20; m.angularDegree = m.angularDegree<35 ? 35 : 53; return m; }
Molecule MakeF2()   // ~/Code/cp2k-runs/F2 geometry (bohr), as IntegrationTests/M_DFT.C
{
    Molecule m;
    m.Insert(Site{9, Vector3D<double>(0, 0, -1.3343)});
    m.Insert(Site{9, Vector3D<double>(0, 0,  1.3343)});
    return m;
}
CalcOptions F2Options() {return {.basis="valence_lowq_vb", .model=Model::LDA, .pp="gth-lda-valence"};}
}

TEST(MeshAudit, AgainstItsOwnMeshEveryTermReadsZero)
{
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA});
    MoleculeCalculation c(setup, Tight());
    ASSERT_TRUE(c.Result());
    const GridAuditResult R=c.MeshAudit(c.SizedSetup()->mesh);
    ASSERT_FALSE(R.terms.empty());
    for (const auto& t : R.terms)
    {
        EXPECT_EQ(t.ΔE, 0.0)       << t.label;
        EXPECT_EQ(t.maxAbsΔH, 0.0) << t.label;
    }
    EXPECT_EQ(R.ΔE, 0.0);
    EXPECT_EQ(R.maxAbsΔH, 0.0);
}

TEST(MeshAudit, AFinerMeshMovesOnlyTheMeshBorneTerms)
{
    for (auto fit : {Hamiltonian::CoulombFit::Auto, Hamiltonian::CoulombFit::Delta})
    {
        const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA, .coulomb={.fit=fit}});
        MoleculeCalculation c(setup, Tight());
        ASSERT_TRUE(c.Result());
        const double E0=c.Energy();
        const GridAuditResult R=c.MeshAudit(Finer(c.SizedSetup()->mesh));
        const bool hartreeOnMesh = fit==Hamiltonian::CoulombFit::Delta;
        double sum=0.0, largest=0.0;
        bool sawExc=false, sawEee=false;
        for (const auto& t : R.terms)
        {
            sum+=t.ΔE; largest=std::max(largest,t.maxAbsΔH);
            const bool meshBorne = t.entries=="Exc" || (hartreeOnMesh && t.entries=="Eee");
            sawExc = sawExc || t.entries=="Exc";
            sawEee = sawEee || t.entries=="Eee";
            if (meshBorne)
            {
                EXPECT_GT(t.maxAbsΔH, 0.0) << t.label << ": a finer mesh must change a mesh-borne term";
                EXPECT_EQ(t.trDdH.size(), 1u) << "one channel: an unpolarized run without symmetry has one block";
            }
            else
            {
                EXPECT_EQ(t.ΔE, 0.0)       << t.label << ": the mesh does not touch this term";
                EXPECT_EQ(t.maxAbsΔH, 0.0) << t.label;
            }
            if (t.maxAbsΔH>0.0) EXPECT_FALSE(t.maxElement.functionI.empty()) << t.label << ": an AO block names its largest element";
        }
        EXPECT_TRUE(sawExc && sawEee);
        EXPECT_NEAR(R.ΔE, sum, 1e-14*std::fabs(E0));
        EXPECT_EQ(R.maxAbsΔH, largest);
        EXPECT_EQ(c.Energy(), E0) << "the audit is fixed-D: the run's answer is untouched";
        // The pilot sized the mesh to the energy run's eps x N_atoms, so a finer mesh moves the energy by less than that.
        EXPECT_LT(std::fabs(R.ΔE), setup->opts.accuracy.Target(ClientConsumes::Energy, 3));
    }
}

// THE COULOMB (ρ) FIT IS ANALYTIC: no mesh reaches the Dunlap Hartree term (user 2026-11-11, "The rho fit should be 100% analytic, no
// grid").  At a fixed D its energy and matrix are BIT-identical on a finer XC mesh, on methane as on water above.  Until 2026-11-11
// the A1_coul fit took the XC MeshParams and built a private mesh for its normalisation and for projecting the SAD seed; that never
// reached the Hartree term at a fixed D (this read 0 then too), only the seed's first iterate, which now names its own quadrature.
TEST(MeshAudit, TheDunlapHartreeIsMeshFreeAtAFixedD)
{
    Molecule ch4;
    const double a=1.186;   // C-H 2.054 bohr
    ch4.Insert(Site{6, Vector3D<double>(0,0,0)});
    for (const rvec3_t& h : {rvec3_t(a,a,a), rvec3_t(-a,-a,a), rvec3_t(-a,a,-a), rvec3_t(a,-a,-a)}) ch4.Insert(Site{1, Vector3D<double>(h.x,h.y,h.z)});
    const auto setup=MoleculeSetup::Resolve(ch4, {.basis="dzvp", .model=Model::LDA});
    MoleculeCalculation c(setup, Tight());
    ASSERT_TRUE(c.Result());
    const GridAuditResult R=c.MeshAudit(Finer(c.SizedSetup()->mesh));
    bool sawEee=false, sawExc=false;
    for (const auto& t : R.terms)
    {
        if (t.entries=="Eee") {sawEee=true; EXPECT_EQ(t.ΔE, 0.0) << t.label; EXPECT_EQ(t.maxAbsΔH, 0.0) << t.label;}
        if (t.entries=="Exc") {sawExc=true; EXPECT_GT(t.maxAbsΔH, 0.0) << t.label << ": the finer mesh is a different mesh";}
    }
    EXPECT_TRUE(sawEee && sawExc);
}

TEST(MeshAudit, BeforeAConvergeItThrows)
{
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA});
    const auto c=MoleculeCalculation::Prepare(setup);
    EXPECT_THROW((void)c->MeshAudit(Finer(setup->mesh)), std::logic_error);
}

// MEASURED 2026-11-08 at the converged D, against (300, 65): water at the entry rung (80, 29) |ΔE| 2e-7, max|ΔH| 1e-6; F2 GTH q7 |ΔE|
// 2.3e-6, max|ΔH| 2.4e-5.  So both meet an ENERGY run's target (eps x N_atoms), and F2 misses a GEOMETRY run's (max|ΔH_ij| <= eps).
TEST(MeshPilot, AnEnergyRunSettlesWhereItsAuditMeetsEps)
{
    for (bool f2 : {false, true})
    {
        const auto setup = f2 ? MoleculeSetup::Resolve(MakeF2(), F2Options())
                              : MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA});
        EXPECT_FALSE(setup->meshStated);
        EXPECT_FALSE(setup->meshResolved);
        EXPECT_EQ(setup->mesh.nRadial,       kMolecularMeshEntryRadial);
        EXPECT_EQ(setup->mesh.angularDegree, kMolecularMeshEntryDegree);
        MoleculeCalculation c(setup, Tight());
        ASSERT_TRUE(c.Result());
        EXPECT_EQ(c.Consumes(), ClientConsumes::Energy);
        const auto cost=c.PilotCost();
        EXPECT_TRUE(cost.settled);
        EXPECT_FALSE(cost.statedOver);
        EXPECT_GE(cost.rounds, 1);
        EXPECT_GE(cost.referencesBuilt, 2) << "each axis's rung above, at least";
        const auto sized=c.SizedSetup();
        EXPECT_TRUE(sized->meshResolved);
        EXPECT_GE(sized->mesh.nRadial,       setup->mesh.nRadial)       << "the climb only goes up";
        EXPECT_GE(sized->mesh.angularDegree, setup->mesh.angularDegree) << "the climb only goes up";
        const bool raised = sized->mesh.nRadial!=setup->mesh.nRadial || sized->mesh.angularDegree!=setup->mesh.angularDegree;
        EXPECT_EQ(cost.reconvergeIterations>0, raised) << "a raise re-converges once from the carried D";
        // The settled mesh meets the energy run's target against a finer one, at the converged D.
        const GridAuditResult R=c.MeshAudit(Finer(sized->mesh));
        EXPECT_LT(std::fabs(R.ΔE), setup->opts.accuracy.Target(ClientConsumes::Energy, c.GetStructure().GetNumAtoms()));
    }
}

TEST(MeshPilot, AGeometryRunRaisesAMeshThatMissesEpsOnMatrixElements)
{
    const auto setup=MoleculeSetup::Resolve(MakeF2(), F2Options());
    MoleculeSession pes(setup, Tight());
    ASSERT_TRUE(pes.ConvergeAt(std::nullopt, pes.Reference().R));
    const MoleculeCalculation& c=*pes.Current();
    EXPECT_EQ(c.Consumes(), ClientConsumes::MatrixElements);
    const auto cost=c.PilotCost();
    EXPECT_TRUE(cost.settled);
    const auto sized=c.SizedSetup();
    EXPECT_TRUE(sized->meshResolved);
    const bool raised = sized->mesh.nRadial>setup->mesh.nRadial || sized->mesh.angularDegree>setup->mesh.angularDegree;
    EXPECT_TRUE(raised) << "F2 GTH misses max|ΔH| <= eps at the entry rung";
    EXPECT_GT(cost.reconvergeIterations, 0u) << "a raise re-converges once from the carried D";
    const GridAuditResult R=c.MeshAudit(Finer(sized->mesh));
    EXPECT_LT(R.maxAbsΔH, setup->opts.accuracy.Target(ClientConsumes::MatrixElements, 2));
}

TEST(MeshPilot, AStatedMeshIsKeptAndAMissIsReported)
{
    CalcOptions o=F2Options();
    o.mesh.nRadial=kMolecularMeshEntryRadial; o.mesh.angularDegree=kMolecularMeshEntryDegree;   // the entry rung, STATED
    const auto setup=MoleculeSetup::Resolve(MakeF2(), o);
    EXPECT_TRUE(setup->meshStated);
    MoleculeSession pes(setup, Tight());   // a geometry run: matrix elements, which F2 misses at this rung
    ASSERT_TRUE(pes.ConvergeAt(std::nullopt, pes.Reference().R));
    const MoleculeCalculation& c=*pes.Current();
    const auto cost=c.PilotCost();
    EXPECT_TRUE(cost.statedOver);
    EXPECT_FALSE(cost.settled);
    EXPECT_EQ(cost.reconvergeIterations, 0u) << "a stated mesh is never raised";
    EXPECT_EQ(c.SizedSetup()->mesh.nRadial,       o.mesh.nRadial);
    EXPECT_EQ(c.SizedSetup()->mesh.angularDegree, o.mesh.angularDegree);
}

TEST(MeshPilot, HartreeFockHasNoMeshToPilot)
{
    MoleculeCalculation c(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::HF}), Tight());
    ASSERT_TRUE(c.Result());
    EXPECT_EQ(c.PilotCost().rounds, 0);
    EXPECT_FALSE(c.SizedSetup()->meshResolved);
}

TEST(MeshPilot, ASessionPilotsItsFirstPointOnMatrixElementsAndStepsTheRestOnIt)
{
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA});
    MoleculeSession pes(setup, Tight());
    const PESGeometry g=pes.Reference();
    ASSERT_TRUE(pes.ConvergeAt(std::nullopt, g.R));
    const MoleculeCalculation* first=pes.Current();
    EXPECT_EQ(first->Consumes(), ClientConsumes::MatrixElements);
    EXPECT_GE(first->PilotCost().rounds, 1);
    EXPECT_TRUE(first->SizedSetup()->meshResolved);
    const qcMesh::MeshParams frozen=first->SizedSetup()->mesh;
    std::vector<rvec3_t> R=g.R;
    R[0]+=rvec3_t(0.0, 0.0, 0.02);
    ASSERT_TRUE(pes.ConvergeAt(std::nullopt, R));
    EXPECT_EQ(pes.Current()->PilotCost().rounds, 0) << "a later point steps on the frozen mesh: no pilot";
    EXPECT_EQ(pes.Current()->SizedSetup()->mesh.nRadial,       frozen.nRadial);
    EXPECT_EQ(pes.Current()->SizedSetup()->mesh.angularDegree, frozen.angularDegree);
}
