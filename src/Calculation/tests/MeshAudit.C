// File: Calculation/tests/MeshAudit.C  The molecular XC-MESH AUDIT at a fixed D (doc/GeometryOptimizationPlan.md §14.3; pin 31).
// H2O LDA/dzvp.  Invariants, not pinned numbers: the run against its OWN mesh reads zero everywhere; against a finer mesh only the
// mesh-borne terms move (XC; and Hartree once it rides the mesh, CoulombFit::Delta); the total is the sum of the terms; the verdict is
// the run's measure against eps x N_atoms; and an audit before a Converge throws.
#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <stdexcept>
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
qcMesh::MeshParams Finer(qcMesh::MeshParams m) { m.nRadial=100; m.angularDegree=35; return m; }
}

TEST(MeshAudit, AgainstItsOwnMeshEveryTermReadsZero)
{
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA});
    MoleculeCalculation c(setup, Tight());
    ASSERT_TRUE(c.Result());
    const GridAuditResult R=c.MeshAudit(setup->opts.mesh);
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
        const GridAuditResult R=c.MeshAudit(Finer(setup->opts.mesh));
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
        // The default (80, 29) is converged on water to about a microhartree per atom (MoleculeSetup.C's measurement), well inside the
        // energy run's eps x N_atoms at the default eps 1e-5.
        EXPECT_LT(std::fabs(R.ΔE), setup->opts.accuracy.Target(ClientConsumes::Energy, 3));
    }
}

TEST(MeshAudit, BeforeAConvergeItThrows)
{
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA});
    const auto c=MoleculeCalculation::Prepare(setup);
    EXPECT_THROW((void)c->MeshAudit(Finer(setup->opts.mesh)), std::logic_error);
}
