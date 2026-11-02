// File: Calculation/tests/SolidSession.C  The SOLID potential energy surface (doc/CalcNetworkReview.md §5 step 8, §7): a
// session over one setup, sampled at positions the caller chooses.  Si Γ (two atoms, eight electrons, ~1 s an SCF) carries
// the gates: (a) ConvergeAt the reference == the single point, bit for bit; (b) a step from the carried density == a fresh run
// at the new positions, in fewer iterations; and the refusals that are outcomes, not crashes (Q3 BasisConditioning included).
#include <gtest/gtest.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>
import qchem.BasisSet.Gaussian.Lattice.GPWTolerances;   // ClientConsumes (what a run's client consumes, C2 step 3)
import qchem.SolidCalculation;
import qchem.Mesh.XCPolicy;   // BeckeXCParams / BeckeLadderAbove (the nested Becke ladder)
import qchem.Materials;
import qchem.Lattice_3D;
import qchem.Structure;
import qchem.UnitCell;
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.SCFParams;
import qchem.SCFAccelerator.Factory;   // SCFAccelerators::Type
import qchem.LASolver;                 // qchem::CholeskyPivoted (the Q3 conditioning window)
import qchem.Types;

using namespace qchem;

namespace
{
SCFParams Tight()
{
    SCFParams p;
    // Δρ/N 1e-8: above Si Γ's Δρ quadrature floor (ρ_lost/N ~1.8e-9 at densityEcut 20), where a 1e-9 gate wanders.
    p.NMaxIter=100; p.MinΔρPerElectron=1e-8; p.MinΔE=std::nullopt; p.MinΔFD=std::nullopt; p.MinVirial=std::nullopt; p.MinFD=std::nullopt;
    return p;
}
SolidCalcOptions SiGamma(bool imposed=false, bool statedRasters=true)
{
    SolidCalcOptions o; o.Nelec=8; o.densityEcut=20.0; o.imposeSymmetry=imposed; o.label="pes";
    o.xcMeshSymmetric=false;   // a session's free setup: the plain XC mesh, the same at every geometry (SolidSession requires it)
    // The SESSION machinery's gates state the rasters (accuracy 0: the ball's own), so a session point and a single point sample on the
    // same grids -- C2 step 3 sizes them by what the client consumes (a session: matrix elements; a single point: the energy), which
    // SessionSizesItsRastersOnMatrixElements gates on its own.
    if (statedRasters) o.accuracy.eps=0.0;
    // ... and the XC mesh: a session point is a geometry run, which pin 31 sends to Becke while the single point's Auto picks uniform on
    // Si Γ -- the machinery gates compare the two on ONE stated mesh (the uniform one, as before pin 31).
    if (statedRasters) o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;
    return o;
}
const GaussianBasisRecipe kSiBasis{.data=BasisSet::Gaussian::BasisSetData::SIPP_SR};

std::shared_ptr<const SolidSetup> SiSetup(const std::vector<rvec3_t>* at=nullptr, bool imposed=false, bool statedRasters=true)
{
    const auto mat=Materials::Get("Si_diamond");
    UnitCell cell(*mat.cell);
    if (at) for (size_t a=0;a<at->size();a++) cell[a]->itsR=(*at)[a];
    return SolidSetup::Resolve(Lattice_3D(cell, ivec3_t(1,1,1)), kSiBasis, SiGamma(imposed, statedRasters));
}
} // namespace

TEST(SolidSession, ConvergeAtTheReferenceIsTheSinglePointBitForBit)
{
    const auto setup=SiSetup();
    SolidCalculation single(setup, Tight());
    auto S=single.Result(); ASSERT_TRUE(S);

    SolidSession pes(setup, {SCFStage{Tight(), SCFAccelerators::Type::Auto}});
    const PESGeometry ref=pes.Reference();
    ASSERT_TRUE(ref.cell) << "a solid's reference has a cell";
    auto P=pes.ConvergeAt(ref.cell, ref.R);
    ASSERT_TRUE(P) << P.Error().details;
    EXPECT_DOUBLE_EQ(P->energy, S->Energy());
    EXPECT_EQ(P->iterations, S->IterationCount());
    EXPECT_FALSE(P->forces);
}

TEST(SolidSession, ACarriedStepEqualsAFreshRunAtTheNewPositionsInFewerIterations)
{
    const auto setup=SiSetup();
    SolidSession pes(setup, {SCFStage{Tight(), SCFAccelerators::Type::Auto}});
    const PESGeometry ref=pes.Reference();
    std::vector<rvec3_t> R=ref.R;
    auto P0=pes.ConvergeAt(ref.cell, R); ASSERT_TRUE(P0) << P0.Error().details;
    R[1].x+=0.05;   // bohr, Cartesian: one atom off its site
    auto P1=pes.ConvergeAt(ref.cell, R); ASSERT_TRUE(P1) << P1.Error().details;

    SolidCalculation fresh(SiSetup(&R), Tight());   // the same recipe resolved AT the new positions, seeded from scratch
    auto F=fresh.Result(); ASSERT_TRUE(F);
    // Each run's stop at Δρ/N 1e-8 is good to ~1e-8 Ha, so two runs agree to their SUM.  Measured 2026-10-30: carried
    // -7.11461659840 (equal to the same run taken to the floor, 29 iterations), fresh 1.1e-8 above it.
    EXPECT_NEAR(P1->energy, F->Energy(), 2e-8);
    EXPECT_GT(P1->energy, P0->energy) << "off the equilibrium site the energy rises";
    EXPECT_LT(P1->iterations, F->IterationCount()) << "the carried density is the point of the session";
}

// C2 step 3 (user ruling 2026-10-31): the run answers what its client consumes, and the pilot sizes on THAT -- a single point on its
// energy (eps per atom), a session point on max|dH_ij| (eps per element: forces and paths read H), so the session's rasters are at least
// the single point's; the session FREEZES them at its first point, so a later step samples the same grids (no surface jump).
TEST(SolidSession, SessionSizesItsRastersOnMatrixElements)
{
    using BasisSet::Gaussian::ClientConsumes;
    // The AUTO density floor (the ball's own raster delivers x ~ 4) and the Becke XC mesh: the LINEAR raster terms (Hartree, V_long) are
    // what an H-graded pilot raises here, in one round.  (A uniform XC raster cannot reach 1e-5 per element on Si at any affordable x --
    // 3e-5 at x 137, C2 step 3 -- which the pilot reports as NOT SETTLED with its advice; not a unit test's cost.)
    const auto mat=Materials::Get("Si_diamond");
    SolidCalcOptions o=SiGamma(false, /*statedRasters*/false);
    o.densityEcut=-1.0;
    // A STATED Becke mesh: the pilot audits it on the geometry run and reports it, never raises it (the grid is a user decision) -- the
    // facade's own choice, sized on the nested ladder, is GeometryRunPicksAndSizesTheBeckeMesh_Long.  (A bare cellKind Becke is
    // BeckeXCParams(), (40, 2.0, 29); this gate STATES a deliberately coarse (30, 1.0, Lebedev 5) one -- its subject is the raster terms,
    // and the calibrated mesh costs 107 s here against 17 s.  The audit reads 1.8e-2 per element on it: OVER, and honoured.)
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Becke; o.xcMesh.nRadial=30; o.xcMesh.mhl_alpha=1.0; o.xcMesh.angularDegree=5;
    const auto setup=SolidSetup::Resolve(Lattice_3D(*mat.cell, ivec3_t(1,1,1)), kSiBasis, o);
    ASSERT_FALSE(setup->rastersResolved) << "the default setup self-calibrates";
    SolidCalculation single(setup, Tight());
    ASSERT_TRUE(single.Result());
    EXPECT_EQ(single.Consumes(), ClientConsumes::Energy);
    const auto sizedE=single.SizedSetup();
    EXPECT_TRUE(sizedE->rastersResolved);
    // What the pilot costs (reported and gated, not hidden): the energy run raises the density raster once from the ball's.
    EXPECT_TRUE(single.PilotCost().settled);
    EXPECT_EQ(single.PilotCost().rounds, 2);
    // 39 measured 2026-10-31: the re-converged energy is settled at its FIRST iterate (dE 1e-9), the rest is GDM's linear Δρ tail down
    // to this gate's Δρ/N 1e-8 (x0.81 per iteration) -- a cost a user pays at a tight Δρ gate, so it is pinned here.
    EXPECT_LE(single.PilotCost().reconvergeIterations, 45u);

    SolidSession pes(setup, {SCFStage{Tight(), SCFAccelerators::Type::Auto}});
    const PESGeometry ref=pes.Reference();
    std::vector<rvec3_t> R=ref.R;
    ASSERT_TRUE(pes.ConvergeAt(ref.cell, R));
    ASSERT_TRUE(pes.Current());
    EXPECT_EQ(pes.Current()->Consumes(), ClientConsumes::MatrixElements);
    // A geometry pilot STARTS AT THE BOUND (rigorous L(eps), eta = 1): the audit confirms it, no re-convergence.
    EXPECT_EQ(pes.Current()->PilotCost().rounds, 1) << "seeded at the bound, the first audit settles";
    EXPECT_EQ(pes.Current()->PilotCost().reconvergeIterations, 0u);
    const auto sizedH=pes.Current()->SizedSetup();
    const auto& tE=sizedE->opts.tolerances; const auto& tH=sizedH->opts.tolerances;
    EXPECT_GE(tH.rasterX, tE.rasterX);  EXPECT_GE(tH.vlocRasterX, tE.vlocRasterX);  EXPECT_GE(tH.xcRasterX, tE.xcRasterX);
    std::cout<<"[raster sizing] energy run x=("<<tE.rasterX<<", "<<tE.vlocRasterX<<", "<<tE.xcRasterX<<")  geometry run x=("
             <<tH.rasterX<<", "<<tH.vlocRasterX<<", "<<tH.xcRasterX<<")"<<std::endl;
    R[1].x+=0.05;
    ASSERT_TRUE(pes.ConvergeAt(ref.cell, R));
    EXPECT_TRUE(pes.Current()->SizedSetup()->opts.tolerances==tH) << "a later point samples on the rasters the first point settled";
}

// THE XC ROUTE FOLLOWS THE CONSUMER (user ruling 2026-11-01, pin 31): with NO deck change -- XC Auto, which routes a Si energy run to the
// uniform raster -- a geometry point takes the Becke mesh and its pilot sizes that mesh on the nested ladder to max|dH_xc| <= eps per element
// (BeckeXCParams() is calibrated on the ENERGY).  Gated by an INDEPENDENT audit two rungs above the settled mesh, not the pilot's own
// reference.  MEASURED 2026-11-01 (Si Gamma, auto floor): the rung differences 1.8e-5 (40/29 vs 60/35), 8.9e-6 (60/35 vs 80/41 GL), 4.7e-7
// (80/41 vs 100/47); the pilot's geometric tail estimate puts 60/35 at 1.8e-5, so it settles on 80/41 in 3 rounds, 9 re-convergence iterations;
// the gate's audit against 120/53 passes.  ~520 s serial (the GL references), hence _Long.  The NaF twin: GPW_NaF.Γ_GeometryRunSizesItsBeckeMesh_Long.
TEST(SolidSession, GeometryRunPicksAndSizesTheBeckeMesh_Long)
{
    const auto mat=Materials::Get("Si_diamond");
    SolidCalcOptions o=SiGamma(false, /*statedRasters*/false);
    o.densityEcut=-1.0;
    const auto setup=SolidSetup::Resolve(Lattice_3D(*mat.cell, ivec3_t(1,1,1)), kSiBasis, o);
    ASSERT_EQ(setup->xcMesh.cellKind, qcMesh::UnitCellKind::Uniform) << "the premise: Auto routes this system's energy run to the uniform raster";

    SolidSession pes(setup, {SCFStage{Tight(), SCFAccelerators::Type::Auto}});
    const PESGeometry ref=pes.Reference();
    std::vector<rvec3_t> R=ref.R;
    ASSERT_TRUE(pes.ConvergeAt(ref.cell, R));
    const SolidCalculation& c=*pes.Current();
    const qcMesh::MeshParams xc=c.ResolvedXCMesh();
    ASSERT_EQ(xc.cellKind, qcMesh::UnitCellKind::Becke) << "a geometry run's XC follows the atoms";
    EXPECT_GT(xc.nRadial, qcMesh::BeckeXCParams().nRadial) << "the energy-calibrated default does not meet 1e-5 per element on Si Gamma";
    EXPECT_TRUE(c.PilotCost().settled);
    EXPECT_LE(c.PilotCost().rounds, 3);
    EXPECT_LE(c.PilotCost().reconvergeIterations, 20u);

    const double eps=o.accuracy.eps;   // per matrix element on a geometry run (GridAccuracy::Target)
    const auto A=c.XCMeshAudit(qcMesh::BeckeLadderAbove(qcMesh::BeckeLadderAbove(xc)));
    double dXC=0.0;
    for (const auto& t : A.terms) if (t.entries=="Exc") dXC=std::max(dXC, t.maxAbsΔH);
    EXPECT_GT(dXC, 0.0) << "the audit compared two different meshes";
    EXPECT_LE(dXC, eps) << "max|dH_xc| against a mesh two rungs finer";

    R[1].x+=0.05;
    ASSERT_TRUE(pes.ConvergeAt(ref.cell, R));
    EXPECT_EQ(pes.Current()->ResolvedXCMesh().nRadial, xc.nRadial) << "a later point samples on the XC mesh the first point settled";
    EXPECT_EQ(pes.Current()->ResolvedXCMesh().angularDegree, xc.angularDegree);
}

TEST(SolidSession, RefusalsAreOutcomesAndAWrongCountThrows)
{
    const auto setup=SiSetup(nullptr, /*imposed*/true);
    SolidSession pes(setup, {SCFStage{Tight(), SCFAccelerators::Type::Auto}});
    const PESGeometry ref=pes.Reference();
    const std::vector<rvec3_t>& R=ref.R;

    auto F=pes.ConvergeAt(ref.cell, R, {.forces=true});
    ASSERT_FALSE(F); EXPECT_EQ(F.Error().why, PESFailure::Why::NoGradient);

    std::vector<rvec3_t> broken=R; broken[1].x+=0.05;
    auto B=pes.ConvergeAt(ref.cell, broken);
    ASSERT_FALSE(B); EXPECT_EQ(B.Error().why, PESFailure::Why::SymmetryBroken) << B.Error().details;

    rmat3d_t strained=*ref.cell; for (size_t i=1;i<=3;i++) for (size_t j=1;j<=3;j++) strained(i,j)*=1.01;   // a 1% isotropic strain
    auto C=pes.ConvergeAt(strained, R);
    ASSERT_FALSE(C); EXPECT_EQ(C.Error().why, PESFailure::Why::CellChange) << C.Error().details;

    EXPECT_THROW((void)pes.ConvergeAt(ref.cell, {R[0]}), std::invalid_argument);
    EXPECT_THROW((void)pes.ConvergeAt(std::nullopt, R), std::invalid_argument) << "a solid's geometry has a cell";

    // A setup over a PRE-BUILT basis cannot follow moved atoms: a programming error, said at once.
    const auto mat=Materials::Get("Si_diamond");
    std::shared_ptr<const BasisSet::Real_BS> mol(BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::SIPP_SR, mat.cell.get()));
    SolidSession fixed(SolidSetup::Resolve(Lattice_3D(*mat.cell, ivec3_t(1,1,1)), mol, SiGamma()), {SCFStage{Tight(), SCFAccelerators::Type::Auto}});
    const PESGeometry fr=fixed.Reference();
    std::vector<rvec3_t> moved=fr.R; moved[1].x+=0.05;
    EXPECT_THROW((void)fixed.ConvergeAt(fr.cell, moved), std::invalid_argument);

    // A FREE setup whose XC mesh follows the group detected at each geometry would switch meshes along a path: refused at once.
    SolidCalcOptions sym=SiGamma(); sym.xcMeshSymmetric=true;
    const auto symSetup=SolidSetup::Resolve(Lattice_3D(*mat.cell, ivec3_t(1,1,1)), kSiBasis, sym);
    EXPECT_THROW(SolidSession(symSetup, {SCFStage{Tight(), SCFAccelerators::Type::Auto}}), std::invalid_argument);
}

// Q3 (doc/CalcNetworkReview.md §7): the trim is decided ONCE, at the reference, and frozen.  A step that pushes the frozen basis
// out of the conditioning window it was trimmed in -- here one Si squeezed to a quarter of the bond, where its diffuse shells and
// its neighbour's become near-dependent -- is REFUSED, never run on a silently re-trimmed (different) basis: a per-step re-trim
// would make the surface discontinuous.  A modest step inside the window still converges.
TEST(SolidSession, AStepOutOfTheTrimWindowIsRefusedAsBasisConditioning)
{
    const auto mat=Materials::Get("Si_diamond");
    SolidCalcOptions o=SiGamma();
    // The window: measured 2026-10-30, min eig S at Γ = 0.027 at the reference and 4.5e-4 squeezed to a quarter bond, so 1e-3
    // sits between them -- the reference vets clean, the squeezed point does not.
    o.ortho=qchem::CholeskyPivoted; o.orthoTol=1e-3;
    GaussianBasisRecipe vetted=kSiBasis; vetted.vet=true;
    SolidSession pes(SolidSetup::Resolve(Lattice_3D(*mat.cell, ivec3_t(1,1,1)), vetted, o), {SCFStage{Tight(), SCFAccelerators::Type::Auto}});
    const PESGeometry ref=pes.Reference();

    std::vector<rvec3_t> squeezed=ref.R; squeezed[1]=ref.R[0]+0.25*(ref.R[1]-ref.R[0]);
    auto B=pes.ConvergeAt(ref.cell, squeezed);
    ASSERT_FALSE(B) << "a near-dependent frozen basis must not run";
    EXPECT_EQ(B.Error().why, PESFailure::Why::BasisConditioning) << B.Error().details;

    std::vector<rvec3_t> step=ref.R; step[1].x+=0.05;
    auto P=pes.ConvergeAt(ref.cell, step);
    ASSERT_TRUE(P) << P.Error().details;
}

// THE LÖWDIN CARRY (OpenWork "PES / run-setup remainders" item 3): along one displacement path, a session that carries D through
// the Löwdin frame lands on the SAME energies as one that carries the AO coefficients verbatim -- the carry is a seed, never the
// answer -- and is measured, not assumed, to save iterations.  Si Γ, one atom walked 0.05 bohr at a time.
TEST(SolidSession, TheLowdinCarryLandsOnTheSameSurface)
{
    const auto setup=SiSetup();
    const PESGeometry ref=SolidSession(setup, {SCFStage{Tight(), SCFAccelerators::Type::Auto}}).Reference();
    auto walk=[&](SolidSession::Carry how)
    {
        SolidSession pes(setup, {SCFStage{Tight(), SCFAccelerators::Type::Auto}}, nullptr, how);
        std::vector<std::pair<double,size_t>> out;
        std::vector<rvec3_t> R=ref.R;
        for (int i=0;i<5;i++)
        {
            auto P=pes.ConvergeAt(ref.cell, R);
            EXPECT_TRUE(P) << P.Error().details;
            out.push_back(P ? std::make_pair(P->energy, P->iterations) : std::make_pair(0.0, size_t(0)));
            R[1].x+=0.05;
        }
        return out;
    };
    const auto ao=walk(SolidSession::Carry::AO), lw=walk(SolidSession::Carry::Lowdin);
    size_t nAO=0, nLW=0;
    for (size_t i=0;i<ao.size();i++)
    {
        std::cout << "[Si Γ carry] step " << i << "  E(AO)=" << std::setprecision(12) << ao[i].first << " (" << ao[i].second
                  << " it)  E(Löwdin)=" << lw[i].first << " (" << lw[i].second << " it)" << std::endl;
        EXPECT_NEAR(lw[i].first, ao[i].first, 2e-8) << "step " << i << ": each stop is good to ~1e-8 Ha, two agree to their sum";
        if (i>0) { nAO+=ao[i].second; nLW+=lw[i].second; }
    }
    std::cout << "[Si Γ carry] carried steps: AO " << nAO << " iterations, Löwdin " << nLW << std::endl;
}
