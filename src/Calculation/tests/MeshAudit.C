// File: Calculation/tests/MeshAudit.C  The molecular XC-MESH AUDIT at a fixed D and the MESH PILOT that climbs on it
// (doc/GeometryOptimizationPlan.md §14.2-14.3; pin 31).  H2O LDA/dzvp and F2 GTH q7.  Invariants, not pinned numbers: the run against its
// OWN mesh reads zero everywhere; against a finer mesh only the mesh-borne terms move (XC; and Hartree once it rides the mesh,
// CoulombFit::Delta); the total is the sum of the terms; an audit before a Converge throws.  The pilot: an unstated mesh settles where
// its audit meets eps (water at the entry rung, F2 above it), a stated one is kept and reported, HF has none, and a session freezes it.
#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
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
RunOptions F2Options() {return {.model=Model::LDA, .pp="gth-lda-valence", .molecule={.basis="valence_lowq_vb"}};}
}

TEST(MeshAudit, AgainstItsOwnMeshEveryTermReadsZero)
{
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.model=Model::LDA, .scf=Tight(), .molecule={.basis="dzvp"}});
    auto c=MoleculeCalculation::Prepare(setup).TakeValue(); (void)c->Converge();
    ASSERT_TRUE(c->Result());
    const GridAuditResult R=c->MeshAudit(c->SizedSetup()->xcMesh);
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
    // Both Hartree routes, STATED: the Dunlap fit (the unnamed Auxiliary file, A2_coul) is mesh-free, the δ fit (the default since the user
    // ruling 2026-11-22) is mesh-borne.
    for (auto fit : {Hamiltonian::CoulombFit::Auxiliary, Hamiltonian::CoulombFit::Delta})
    {
        const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.model=Model::LDA, .coulomb={.fit=fit}, .scf=Tight(), .molecule={.basis="dzvp"}});
        auto c=MoleculeCalculation::Prepare(setup).TakeValue(); (void)c->Converge();
        ASSERT_TRUE(c->Result());
        const double E0=c->Energy();
        const GridAuditResult R=c->MeshAudit(Finer(c->SizedSetup()->xcMesh));
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
        EXPECT_EQ(c->Energy(), E0) << "the audit is fixed-D: the run's answer is untouched";
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
    const auto setup=MoleculeSetup::Resolve(ch4, {.model=Model::LDA, .coulomb={.fit=Hamiltonian::CoulombFit::Auxiliary}, .scf=Tight(), .molecule={.basis="dzvp"}});   // the Dunlap route STATED
    auto c=MoleculeCalculation::Prepare(setup).TakeValue(); (void)c->Converge();
    ASSERT_TRUE(c->Result());
    const GridAuditResult R=c->MeshAudit(Finer(c->SizedSetup()->xcMesh));
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
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.model=Model::LDA, .molecule={.basis="dzvp"}});
    const auto c=MoleculeCalculation::Prepare(setup).TakeValue();
    EXPECT_THROW((void)c->MeshAudit(Finer(setup->xcMesh)), std::logic_error);
}

// MEASURED 2026-11-08 at the converged D, against (300, 65): water at the entry rung (80, 29) |ΔE| 2e-7, max|ΔH| 1e-6; F2 GTH q7 |ΔE|
// 2.3e-6, max|ΔH| 2.4e-5.  So both meet an ENERGY run's target (eps x N_atoms), and F2 misses a GEOMETRY run's (max|ΔH_ij| <= eps).
TEST(MeshPilot, AnEnergyRunSettlesWhereItsAuditMeetsEps)
{
    for (bool f2 : {false, true})
    {
        RunOptions F2o=F2Options(); F2o.scf=Tight();
        const auto setup = f2 ? MoleculeSetup::Resolve(MakeF2(), F2o)
                              : MoleculeSetup::Resolve(GetMolecule("H2O"), {.model=Model::LDA, .scf=Tight(), .molecule={.basis="dzvp"}});
        EXPECT_FALSE(setup->meshStated);
        EXPECT_FALSE(setup->meshResolved);
        EXPECT_EQ(setup->xcMesh.nRadial,       kMolecularMeshEntryRadial);
        EXPECT_EQ(setup->xcMesh.angularDegree, kMolecularMeshEntryDegree);
        auto c=MoleculeCalculation::Prepare(setup).TakeValue(); (void)c->Converge();
        ASSERT_TRUE(c->Result());
        EXPECT_EQ(c->Consumes(), ClientConsumes::Energy);
        const auto cost=c->PilotCost();
        EXPECT_TRUE(cost.settled);
        EXPECT_FALSE(cost.statedOver);
        EXPECT_GE(cost.rounds, 1);
        EXPECT_GE(cost.referencesBuilt, 2) << "each axis's rung above, at least";
        const auto sized=c->SizedSetup();
        EXPECT_TRUE(sized->meshResolved);
        EXPECT_GE(sized->xcMesh.nRadial,       setup->xcMesh.nRadial)       << "the climb only goes up";
        EXPECT_GE(sized->xcMesh.angularDegree, setup->xcMesh.angularDegree) << "the climb only goes up";
        const bool raised = sized->xcMesh.nRadial!=setup->xcMesh.nRadial || sized->xcMesh.angularDegree!=setup->xcMesh.angularDegree;
        EXPECT_EQ(cost.reconvergeIterations>0, raised) << "a raise re-converges once from the carried D";
        // The settled mesh meets the energy run's target against a finer one, at the converged D.
        const GridAuditResult R=c->MeshAudit(Finer(sized->xcMesh));
        EXPECT_LT(std::fabs(R.ΔE), setup->opts.accuracy.Target(ClientConsumes::Energy, c->GetStructure().GetNumAtoms()));
    }
}

// The Hartree is STATED as the Dunlap fit, so the climb is the XC term's.  On the δ default (user ruling 2026-11-22) the δ Hartree's matrix
// elements still miss max|ΔH| <= eps at (140, 53), the ladder's top, so the pilot never settles (lane S 2026-11-22: OpenWork row).
TEST(MeshPilot, AGeometryRunRaisesAMeshThatMissesEpsOnMatrixElements)
{
    RunOptions o=F2Options(); o.coulomb.fit=Hamiltonian::CoulombFit::Auxiliary;
    const auto setup=MoleculeSetup::Resolve(MakeF2(), o);
    MoleculeSession pes(setup, Tight());
    ASSERT_TRUE(pes.ConvergeAt(std::nullopt, pes.Reference().R));
    const MoleculeCalculation& c=*pes.Current();
    EXPECT_EQ(c.Consumes(), ClientConsumes::MatrixElements);
    const auto cost=c.PilotCost();
    EXPECT_TRUE(cost.settled);
    const auto sized=c.SizedSetup();
    EXPECT_TRUE(sized->meshResolved);
    const bool raised = sized->xcMesh.nRadial>setup->xcMesh.nRadial || sized->xcMesh.angularDegree>setup->xcMesh.angularDegree;
    EXPECT_TRUE(raised) << "F2 GTH misses max|ΔH| <= eps at the entry rung";
    EXPECT_GT(cost.reconvergeIterations, 0u) << "a raise re-converges once from the carried D";
    // The run's iteration count is its WHOLE run (OpenWork row, the atom's twin 00dcde5b0): the pilot rebuilds the iterator per raised
    // mesh, so the last iterator alone counted only the last re-convergence (equal to or below the pilot's share).
    EXPECT_GT(c.IterationCount(), cost.reconvergeIterations) << "the entry SCF's iterations are part of the run";
    const GridAuditResult R=c.MeshAudit(Finer(sized->xcMesh));
    EXPECT_LT(R.maxAbsΔH, setup->opts.accuracy.Target(ClientConsumes::MatrixElements, 2));
}

TEST(MeshPilot, AStatedMeshIsKeptAndAMissIsReported)
{
    RunOptions o=F2Options();
    MoleculeXCMesh(o).nRadial=kMolecularMeshEntryRadial; MoleculeXCMesh(o).angularDegree=kMolecularMeshEntryDegree;   // the entry rung, STATED
    const auto setup=MoleculeSetup::Resolve(MakeF2(), o);
    EXPECT_TRUE(setup->meshStated);
    MoleculeSession pes(setup, Tight());   // a geometry run: matrix elements, which F2 misses at this rung
    ASSERT_TRUE(pes.ConvergeAt(std::nullopt, pes.Reference().R));
    const MoleculeCalculation& c=*pes.Current();
    const auto cost=c.PilotCost();
    EXPECT_TRUE(cost.statedOver);
    EXPECT_FALSE(cost.settled);
    EXPECT_EQ(cost.reconvergeIterations, 0u) << "a stated mesh is never raised";
    EXPECT_EQ(c.SizedSetup()->xcMesh.nRadial,       MoleculeXCMesh(o).nRadial);
    EXPECT_EQ(c.SizedSetup()->xcMesh.angularDegree, MoleculeXCMesh(o).angularDegree);
}

TEST(MeshPilot, HartreeFockHasNoMeshToPilot)
{
    auto c=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), {.model=Model::HF, .scf=Tight(), .molecule={.basis="dzvp"}})).TakeValue(); (void)c->Converge();
    ASSERT_TRUE(c->Result());
    EXPECT_EQ(c->PilotCost().rounds, 0);
    EXPECT_FALSE(c->SizedSetup()->meshResolved);
}

TEST(MeshPilot, ASessionPilotsItsFirstPointOnMatrixElementsAndStepsTheRestOnIt)
{
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.model=Model::LDA, .molecule={.basis="dzvp"}});
    MoleculeSession pes(setup, Tight());
    const PESGeometry g=pes.Reference();
    ASSERT_TRUE(pes.ConvergeAt(std::nullopt, g.R));
    const MoleculeCalculation* first=pes.Current();
    EXPECT_EQ(first->Consumes(), ClientConsumes::MatrixElements);
    EXPECT_GE(first->PilotCost().rounds, 1);
    EXPECT_TRUE(first->SizedSetup()->meshResolved);
    const qcMesh::MeshParams frozen=first->SizedSetup()->xcMesh;
    std::vector<rvec3_t> R=g.R;
    R[0]+=rvec3_t(0.0, 0.0, 0.02);
    ASSERT_TRUE(pes.ConvergeAt(std::nullopt, R));
    EXPECT_EQ(pes.Current()->PilotCost().rounds, 0) << "a later point steps on the frozen mesh: no pilot";
    EXPECT_EQ(pes.Current()->SizedSetup()->xcMesh.nRadial,       frozen.nRadial);
    EXPECT_EQ(pes.Current()->SizedSetup()->xcMesh.angularDegree, frozen.angularDegree);
}

// ---- THE FOLD AUDIT (the per-orbit consumer of the molecular mesh fold, OpenWork §3; lane H 2026-11-24) -----------------------------------
// On a SALC-imposed run every orbital block is a 1-D real irrep of the mesh's group H, so every integrand the δ XC route forms in a block is
// H-invariant pointwise, and one point per orbit at the orbit's summed weight integrates it EXACTLY.  The claim: at the converged D the
// folded Hamiltonian (one point per orbit) and the unfolded one (every point) agree term by term to roundoff -- E to 1e-10, every
// matrix element to 1e-12 -- LDA and PBE, closed and open shell.
namespace
{
Molecule MakeO2()   // as IntegrationTests/M_DFT.C, along z: D∞h, its D2h subgroup
{
    Molecule m;
    m.Insert(Site{8, Vector3D<double>(0, 0,  1.141)});
    m.Insert(Site{8, Vector3D<double>(0, 0, -1.141)});
    return m;
}
Molecule MakeCH4()   // Td; H's D2h-type subgroup is the C2v BuildAbelianGroup picks (its frame 45° from these axes)
{
    Molecule ch4;
    const double a=1.186;
    ch4.Insert(Site{6, Vector3D<double>(0,0,0)});
    for (const rvec3_t& h : {rvec3_t(a,a,a), rvec3_t(-a,-a,a), rvec3_t(-a,a,-a), rvec3_t(a,-a,-a)}) ch4.Insert(Site{1, Vector3D<double>(h.x,h.y,h.z)});
    return ch4;
}
void ExpectTheFoldIsExact(const Molecule& mol, RunOptions o, const char* what)
{
    o.imposeSymmetry=true; o.scf=Tight(); o.molecule.basis="dzvp";
    auto c=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(mol, o)).TakeValue(); (void)c->Converge();
    ASSERT_TRUE(c->Result()) << what;
    const GridAuditResult R=c->FoldAudit();
    bool sawExc=false;
    for (const auto& t : R.terms)
    {
        sawExc = sawExc || t.entries=="Exc";
        EXPECT_LT(std::fabs(t.ΔE), 1e-10)   << what << ", " << t.label;
        EXPECT_LT(t.maxAbsΔH,      1e-12)   << what << ", " << t.label;
    }
    EXPECT_TRUE(sawExc) << what;
    std::printf("[fold audit] %s: |ΔE| %.2e, max|ΔH| %.2e\n", what, std::fabs(R.ΔE), R.maxAbsΔH);
}
} //anon

TEST(FoldAudit, OnePointPerOrbitIsExactOnAnImposedRun_H2O_LDA)    {ExpectTheFoldIsExact(GetMolecule("H2O"), {.model=Model::LDA}, "H2O C2v LDA");}
TEST(FoldAudit, OnePointPerOrbitIsExactOnAnImposedRun_H2O_PBE)    {ExpectTheFoldIsExact(GetMolecule("H2O"), {.model=Model::PBE}, "H2O C2v PBE");}
TEST(FoldAudit, OnePointPerOrbitIsExactOnAnImposedRun_CH4_LDA)    {ExpectTheFoldIsExact(MakeCH4(), {.model=Model::LDA}, "CH4 Td->C2v LDA");}
TEST(FoldAudit, OnePointPerOrbitIsExactOnAnImposedRun_O2_LSDA_Long)    {ExpectTheFoldIsExact(MakeO2(), {.multiplicity=3, .model=Model::LDA}, "O2 triplet D2h LSDA");}
TEST(FoldAudit, OnePointPerOrbitIsExactOnAnImposedRun_O2_PBE_Long)    {ExpectTheFoldIsExact(MakeO2(), {.multiplicity=3, .model=Model::PBE}, "O2 triplet D2h PBE");}

TEST(FoldAudit, AFreeRunIsRefused)
{
    auto c=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), {.model=Model::LDA, .scf=Tight(), .molecule={.basis="dzvp"}})).TakeValue();
    (void)c->Converge();
    ASSERT_TRUE(c->Result());
    EXPECT_THROW((void)c->FoldAudit(), std::invalid_argument) << "a free run's density is not invariant: folding would change the answer";
}
