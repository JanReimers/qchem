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
import qchem.Mesh.XCPolicy;   // BeckeXCParams / the Becke ladder's rungs
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
import qchem.Blaze;                    // rvec_t arithmetic (q0 + step)
import qchem.SCFIterator;              // IterationRow (the observer counts SCF iterations)

using namespace qchem;

namespace
{
SCFParams Tight(double minΔρPerElectron=1e-8)
{
    SCFParams p;
    // Δρ/N 1e-8: above Si Γ's Δρ quadrature floor (ρ_lost/N ~1.8e-9 at densityEcut 20), where a 1e-9 gate wanders.
    p.NMaxIter=100; p.MinΔρPerElectron=minΔρPerElectron; p.MinΔE=std::nullopt; p.MinΔFD=std::nullopt; p.MinVirial=std::nullopt; p.MinFD=std::nullopt;
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

//! Si at the reference with a cell STRAIN ALLOWANCE: every size decided on the reference grown by \a allowance (the frozen rule).
std::shared_ptr<const SolidSetup> SiStrainSetup(double allowance)
{
    SolidCalcOptions o=SiGamma();
    o.cellStrainAllowance=allowance;
    return SolidSetup::Resolve(Lattice_3D(UnitCell(*Materials::Get("Si_diamond").cell), ivec3_t(1,1,1)), kSiBasis, o);
}
//! The Si reference under the homogeneous strain F: the cell and the atoms ride it.
std::pair<rmat3d_t, std::vector<rvec3_t>> Strained(const PESGeometry& ref, const rmat3d_t& F)
{
    std::vector<rvec3_t> R; for (const auto& r : ref.R) R.push_back(F*r);
    return {F*(*ref.cell), R};
}
rmat3d_t Hydrostatic(double e) {return rmat3d_t(1+e,0,0, 0,1+e,0, 0,0,1+e);}

std::shared_ptr<const SolidSetup> SiSetup(const std::vector<rvec3_t>* at=nullptr, bool imposed=false, bool statedRasters=true)
{
    const auto mat=Materials::Get("Si_diamond");
    UnitCell cell(*mat.cell);
    if (at) cell.MoveSites(*at);
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
    o.xcMesh.cellKind=qcMesh::CellKind::RecipeStated(qcMesh::UnitCellKind::Becke); o.xcMesh.nRadial=30; o.xcMesh.mhl_alpha=1.0; o.xcMesh.angularDegree=5; o.xcMesh.angRot=0.0;
    const auto setup=SolidSetup::Resolve(Lattice_3D(*mat.cell, ivec3_t(1,1,1)), kSiBasis, o);
    ASSERT_FALSE(setup->rastersResolved) << "the default setup self-calibrates";
    // Δρ/N 1e-7, not Tight()'s 1e-8 (user, 2026-11-10; see GeometryRunPicksAndSizesTheBeckeMesh_Long): a pilot re-converges from an
    // already-converged D, and at 1e-8 GDM's line search works at the energy's roundoff floor, so the counts below went chaotic (8 stalls,
    // 41 -> 54 and 15 -> 35 under a 4-thread roundoff change, 2026-11-10).  Above the floor the caps are regression detectors again.
    SolidCalculation single(setup, Tight(1e-7));
    ASSERT_TRUE(single.Result());
    EXPECT_EQ(single.Consumes(), ClientConsumes::Energy);
    const auto sizedE=single.SizedSetup();
    EXPECT_TRUE(sizedE->rastersResolved);
    // What the pilot costs (reported and gated, not hidden): the energy run raises the density raster once from the ball's.
    EXPECT_TRUE(single.PilotCost().settled);
    EXPECT_EQ(single.PilotCost().rounds, 2);
    // 2 measured 2026-11-10 at Δρ/N 1e-7, serial and 4 threads alike (39-54 at 1e-8, GDM's tail into the noise floor).
    EXPECT_LE(single.PilotCost().reconvergeIterations, 4u);

    SolidSession pes(setup, {SCFStage{Tight(1e-7), SCFAccelerators::Type::Auto}});
    const PESGeometry ref=pes.Reference();
    std::vector<rvec3_t> R=ref.R;
    ASSERT_TRUE(pes.ConvergeAt(ref.cell, R));
    ASSERT_TRUE(pes.Current());
    EXPECT_EQ(pes.Current()->Consumes(), ClientConsumes::MatrixElements);
    // A geometry pilot STARTS ITS RASTERS AT THE BOUND (rigorous L(eps), eta = 1): they meet it at once.  The DENSITY BALL does not
    // (2026-11-02, OpenWork F "density ball unjudged"): measured on this cell, max|dH| of Hartree against twice the ball 1.3e-2 / 6.7e-5 /
    // 1.4e-5 / 6.4e-7 at cutoffFactor 2 / 4 / 8 / 16, so it settles on 16 in 4 rounds, 24 re-convergence iterations.
    // ... but its STATED coarse Becke mesh reads 1.8e-2 per element, OVER: honoured and reported, so the pilot is NOT settled (pin 31 rule 3;
    // it read settled until lane G's Si Γ catch, 2026-11-07) and says why: statedOver, with a StatedXCMeshMissesEps warn row.
    EXPECT_FALSE(pes.Current()->PilotCost().settled) << "a stated grid over its target is not settled";
    EXPECT_TRUE(pes.Current()->PilotCost().statedOver);
    EXPECT_FALSE(single.PilotCost().statedOver) << "the energy run does not audit the Becke mesh (an H-run audit)";
    EXPECT_EQ(pes.Current()->PilotCost().rounds, 4) << "the rasters at the bound; the ball climbs cutoffFactor 2 -> 16";
    EXPECT_LE(pes.Current()->PilotCost().reconvergeIterations, 16u) << "9 serial / 13 at 4 threads, 0 stalls (2026-11-10, Δρ/N 1e-7)";
    EXPECT_EQ(pes.Current()->SizedSetup()->opts.cutoffFactor, 16.0);
    EXPECT_EQ(single.SizedSetup()->opts.cutoffFactor, 4.0) << "the energy run's ball: 1.2e-3 Ha at 2, 3.3e-6 at 4";
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
// reference.  MEASURED 2026-11-05 (Si Gamma, auto floor): the run enters at the calibrated rung, nR 100 / rotated Lebedev 53, and settles there in
// 4 rounds, 30 re-convergence iterations, the density ball setting the rounds: angular 1.4e-7, radial 5.1e-8, so Si sits far under its target on
// a rung calibrated on MnO.  ~1600 s at -j3, hence _Long.  The NaF twin: GPW_NaF.Γ_GeometryRunSizesItsBeckeMesh_Long.
TEST(SolidSession, GeometryRunPicksAndSizesTheBeckeMesh_Long)
{
    const auto mat=Materials::Get("Si_diamond");
    SolidCalcOptions o=SiGamma(false, /*statedRasters*/false);
    o.densityEcut=-1.0;
    const auto setup=SolidSetup::Resolve(Lattice_3D(*mat.cell, ivec3_t(1,1,1)), kSiBasis, o);
    ASSERT_EQ(setup->xcMesh.cellKind, qcMesh::UnitCellKind::Uniform) << "the premise: Auto routes this system's energy run to the uniform raster";

    // Δρ/N 1e-7, not Tight()'s 1e-8 (user, 2026-11-10): the pilot RE-converges from an already-converged D, and GDM descends E, which
    // resolves a density change δ only as δ²·E -- at 1e-8 that is ~1e-15 Ha on Si's 7.1 Ha, the energy's roundoff floor.  There GDM's
    // line search compared energies 0-16 ulps apart, stalled, reset its history, and the count went 20/24/30/52 across bit-only changes
    // (lane S bisect, 2026-11-10).  Above the floor the iteration cap below is a regression detector again.
    SolidSession pes(setup, {SCFStage{Tight(1e-7), SCFAccelerators::Type::Auto}});
    const PESGeometry ref=pes.Reference();
    std::vector<rvec3_t> R=ref.R;
    ASSERT_TRUE(pes.ConvergeAt(ref.cell, R));
    const SolidCalculation& c=*pes.Current();
    const qcMesh::MeshParams xc=c.ResolvedXCMesh();
    ASSERT_EQ(xc.cellKind, qcMesh::UnitCellKind::Becke) << "a geometry run's XC follows the atoms";
    EXPECT_GT(xc.nRadial, qcMesh::BeckeXCParams().nRadial) << "the energy-calibrated default does not meet 1e-5 per element on Si Gamma";
    EXPECT_TRUE(c.PilotCost().settled);
    EXPECT_LE(c.PilotCost().rounds, 4) << "the density ball doubles three times beside the mesh's climb (cutoffFactor 2 -> 16, 2026-11-02)";
    EXPECT_LE(c.PilotCost().reconvergeIterations, 28u) << "three re-convergences (24 measured 2026-11-10 at Δρ/N 1e-7, 0 GDM stalls)";
    // The pilot reads the mesh's error off its raster reference, so the mesh's XC term must not feel a raster at all.
    const auto G=c.GridAudit();
    ASSERT_TRUE(G) << G.Error().Summary();
    for (const auto& t : G->terms) if (t.entries=="Exc") EXPECT_EQ(t.maxAbsΔH, 0.0) << "a raster reached the Becke XC term: " << t.label;

    const double eps=o.accuracy.eps;   // per matrix element on a geometry run (GridAccuracy::Target)
    const qcMesh::BeckeLadderCalibration cal;   // two rungs up on the angular axis, one on the radial
    const auto A=c.XCMeshAudit(*qcMesh::BeckeRadialRungAbove(*qcMesh::BeckeAngularRungAbove(*qcMesh::BeckeAngularRungAbove(xc, cal, false), cal, false)));
    ASSERT_TRUE(A) << A.Error().Summary();
    double dXC=0.0;
    for (const auto& t : A->terms) if (t.entries=="Exc") dXC=std::max(dXC, t.maxAbsΔH);
    EXPECT_GT(dXC, 0.0) << "the audit compared two different meshes";
    EXPECT_LE(dXC, eps) << "max|dH_xc| against a mesh two rungs finer";
    // Cauchy-Schwarz in the S metric ties the D-weighted measure to the norm: |Tr(D dH)| <= |dH|_S |D|_S, with |D|_S = sqrt(sum f^2).
    const double normD=4.0;   // Si2, 8 electrons: at most four orbitals at f = 2 in a channel
    for (const auto& t : A->terms)
    {
        ASSERT_FALSE(t.trDdH.empty()) << t.label;
        for (double c : t.trDdH) EXPECT_LE(std::fabs(c), t.sNormDH*normD*(1.0+1e-10)) << t.label;
        if (t.maxAbsΔH>0.0) { EXPECT_GT(t.sNormDH, 0.0) << t.label; EXPECT_FALSE(t.maxElement.functionI.empty()) << t.label; }
    }

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
    EXPECT_THROW((void)pes.ConvergeAt(ref.cell, broken), std::invalid_argument) << "an off-group R on an imposed setup is a wiring defect";

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
    GaussianBasisRecipe vetted=kSiBasis; vetted.autoTrim=true;
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

// THE GEOMETRY RUN FACTORY (OpenWork §2 geometry row, (a)): the session and its coordinate map are built from ONE group, so every
// geometry the map can produce is one the session accepts.  Si Γ carries it: perfect diamond imposed has no free parameter (both atoms
// on 8a); the bond stretched along [111] lowers Oh to D3d -- the inversion through the bond centre still swaps the two atoms, so ONE
// orbit with one free z, the bond length; a free setup gets P1, 3N coordinates with the rigid translations held.
TEST(SolidSession, AGeometryRunPairsTheSessionWithAMapThatKeepsItsGroup)
{
    const std::vector<SCFStage> stages{SCFStage{Tight(), SCFAccelerators::Type::Auto}};
    {
        auto run=MakeSolidGeometryRun(SiSetup(nullptr, /*imposed*/true), stages);
        EXPECT_EQ(run.map->Dimension(), 0u) << "diamond Si: both atoms on the 8a sites, no free Wyckoff parameter";
        // A CELL relaxation needs the frozen rule's sizing box: this setup states no strain allowance, so the factory refuses it.
        EXPECT_THROW((void)MakeSolidGeometryRun(SiSetup(nullptr, /*imposed*/true), stages, nullptr, SolidSession::Carry::AO, /*relaxCell*/true),
                     std::invalid_argument);
        // With one, the map is the cell map over the SAME group: cubic diamond keeps one strain (the breathing mode) and no Wyckoff
        // parameter, and every point asks for the stress as well as the forces.
        SolidCalcOptions o=SiGamma(/*imposed*/true, /*statedRasters*/true); o.cellStrainAllowance=0.02;
        const auto mat=Materials::Get("Si_diamond");
        auto cellRun=MakeSolidGeometryRun(SolidSetup::Resolve(Lattice_3D(UnitCell(*mat.cell), ivec3_t(1,1,1)), kSiBasis, o), stages, nullptr,
                                          SolidSession::Carry::AO, /*relaxCell*/true);
        EXPECT_EQ(cellRun.map->Dimension(), 1u) << "cubic: the breathing strain only";
        EXPECT_TRUE(cellRun.map->PointRequest().forces);
        EXPECT_TRUE(cellRun.map->PointRequest().stress);
    }
    const PESGeometry perfect=SolidSession(SiSetup(), stages).Reference();
    std::vector<rvec3_t> stretched=perfect.R;
    stretched[1]+=rvec3_t(1,1,1)*(0.05/std::sqrt(3.0));   // along the bond, the 3-fold axis through atom 0
    auto run=MakeSolidGeometryRun(SiSetup(&stretched, /*imposed*/true), stages);
    ASSERT_EQ(run.map->Dimension(), 1u) << "D3d: the two atoms are one orbit, one free z along the 3-fold axis";
    const rvec_t q0=run.map->Start();
    rvec_t v(1); v[0]=0.04;
    const rvec_t s=run.map->Project(q0, v);
    EXPECT_NEAR(s[0], v[0], 1e-12) << "a non-polar group holds nothing";
    const PESGeometry g=run.map->At(q0+s);
    const rvec3_t b0=stretched[1]-stretched[0], b1=g.R[1]-g.R[0];
    EXPECT_NEAR(norm(b1-b0*(norm(b1)/norm(b0))), 0.0, 1e-12) << "the step stretches the bond along [111], nothing else";
    EXPECT_GT(std::fabs(norm(b1)-norm(b0)), 1e-3) << "and it does move";
    for (size_t a=0;a<g.R.size();a++) EXPECT_LT(norm(run.map->At(q0).R[a]-stretched[a]), 1e-12) << "q0 is the reference";
    // The real wiring proof: the map's step is a geometry the IMPOSED session accepts (an off-group one would throw).
    auto P=run.surface->ConvergeAt(g.cell, g.R);
    ASSERT_TRUE(P) << P.Error().details;

    // FREE: P1, all 3N coordinates; a step never moves the cell's contents as a whole.
    auto free=MakeSolidGeometryRun(SiSetup(), stages);
    ASSERT_EQ(free.map->Dimension(), 6u);
    rvec_t w(6); for (size_t i=0;i<6;i++) w[i]=0.1*double(i+1)-0.25;
    const PESGeometry f=free.map->At(free.map->Start()+free.map->Project(free.map->Start(), w));
    const rvec3_t shift=(f.R[0]+f.R[1])-(perfect.R[0]+perfect.R[1]);
    EXPECT_LT(norm(shift), 1e-12) << "the rigid translations are held";
}

// FORCES: the session names what is missing, from the point's own Hamiltonian, instead of a blanket refusal -- and spends NO SCF finding
// out: MissingGradients needs the built terms and blocks, not a converged density.  The first point is built to its Hamiltonian, asked and
// discarded; a later ask reads the point the session holds.  The optimizer over a factory-built run stops on that, named.
TEST(SolidSession, AForceRequestNamesTheMissingGradientsBeforeAnySCF)
{
    size_t rows=0;
    auto run=MakeSolidGeometryRun(SiSetup(), {SCFStage{Tight(), SCFAccelerators::Type::Auto}}, [&](const SCFIterator::IterationRow&){ ++rows; });
    const PESGeometry ref=run.surface->Reference();
    auto F=run.surface->ConvergeAt(ref.cell, ref.R, {.forces=true});
    ASSERT_FALSE(F);
    EXPECT_EQ(F.Error().why, PESFailure::Why::NoGradient);
    EXPECT_NE(F.Error().details.find("MissingGradients"), std::string::npos) << F.Error().details;
    EXPECT_NE(F.Error().details.find("\n    "), std::string::npos) << "at least one term named: " << F.Error().details;
    std::cout << "[Si Γ forces] " << F.Error().details << std::endl;
    EXPECT_EQ(rows, 0u) << "the refusal spent no SCF iteration";
    EXPECT_EQ(run.surface->Current(), nullptr) << "the unconverged point is discarded, never held as a carry";
    // STRESS, the same way: the GPW terms without a virial (MissingVirials) are named, and no SCF is spent.
    auto S=run.surface->ConvergeAt(ref.cell, ref.R, {.stress=true});
    ASSERT_FALSE(S);
    EXPECT_EQ(S.Error().why, PESFailure::Why::NoGradient);
    EXPECT_NE(S.Error().details.find("MissingVirials"), std::string::npos) << S.Error().details;
    std::cout << "[Si Γ stress] " << S.Error().details << std::endl;
    EXPECT_EQ(rows, 0u) << "the stress refusal spent no SCF iteration";
    EXPECT_EQ(run.surface->Current(), nullptr);

    ASSERT_TRUE(run.surface->ConvergeAt(ref.cell, ref.R));   // an ENERGY point: held
    const SolidCalculation* held=run.surface->Current();
    ASSERT_NE(held, nullptr);
    const size_t afterEnergy=rows;
    std::vector<rvec3_t> R=ref.R; R[1].x+=0.05;
    auto G=run.surface->ConvergeAt(ref.cell, R, {.forces=true});
    ASSERT_FALSE(G); EXPECT_EQ(G.Error().why, PESFailure::Why::NoGradient);
    EXPECT_EQ(run.surface->Current(), held) << "a later ask reads the held point: no new point";
    EXPECT_EQ(rows, afterEnergy) << "and no SCF";

    OptimizerOptions o; o.maxSteps=3;
    auto opt=GeometryOptimizer(*run.surface, *run.map, o).Run();
    ASSERT_FALSE(opt);
    EXPECT_EQ(opt.Error().why, OptimizeFailure::Why::Cost);
    ASSERT_TRUE(opt.Error().cost);
    EXPECT_EQ(opt.Error().cost->why, PESFailure::Why::NoGradient);
    EXPECT_EQ(rows, afterEnergy) << "the optimizer's refusal spent no SCF either";
}

// THE CELL HALF, step 2 (doc/Records/CellHalf_Design.md §3 and §4 step 2, ruled 2026-11-10): the session FOLLOWS a strained cell under
// the frozen rule -- every size decided once on the reference grown by the allowance -- so E(cell) is one smooth surface.
TEST(SolidSession, ACellWithinTheAllowanceIsFollowedAndARoundTripReturns)
{
    SolidSession pes(SiStrainSetup(0.03), {SCFStage{Tight(), SCFAccelerators::Type::Auto}});
    const PESGeometry ref=pes.Reference();
    auto P0=pes.ConvergeAt(ref.cell, ref.R); ASSERT_TRUE(P0) << P0.Error().details;
    const auto [A1, R1]=Strained(ref, Hydrostatic(0.01));
    auto P1=pes.ConvergeAt(A1, R1); ASSERT_TRUE(P1) << P1.Error().details;
    EXPECT_GT(std::fabs(P1->energy-P0->energy), 1e-5) << "a 1% strain moves the energy";
    auto P2=pes.ConvergeAt(ref.cell, ref.R); ASSERT_TRUE(P2) << P2.Error().details;
    EXPECT_NEAR(P2->energy, P0->energy, 2e-8) << "back at the reference, on the same sets: the same answer (each stop good to ~1e-8)";

    const auto [A4, R4]=Strained(ref, Hydrostatic(0.04));
    auto P4=pes.ConvergeAt(A4, R4);
    ASSERT_FALSE(P4); EXPECT_EQ(P4.Error().why, PESFailure::Why::CellChange);
    EXPECT_NE(P4.Error().details.find("allowance"), std::string::npos) << P4.Error().details;
}

// THE WHOLE-CELL ORACLE's precondition: the re-converged E(e) is differentiable, so central differences at h and 2h agree to the Taylor
// term E'''h^2/2 plus the SCF noise.  A G set, raster or mesh count re-decided at each cell jumps E inside the scan, and a jump d moves
// the difference by ~d/(4h): at h = 5e-4 the gate sees d >~ 6e-8 Ha.  Measured 2026-11-10 at h = 2e-3: dE/de -1.41076 / -1.41093,
// d2E/de2 13.495 / 13.482 (E''' ~ -86, so its Taylor term is ~1e-5 at h = 5e-4).
TEST(SolidSession, TheReconvergedEnergyIsDifferentiableInTheStrain)
{
    SolidSession pes(SiStrainSetup(0.03), {SCFStage{Tight(1e-9), SCFAccelerators::Type::Auto}});
    const PESGeometry ref=pes.Reference();
    const double h=5e-4;
    auto E=[&](double e)
    {
        const auto [A, R]=Strained(ref, Hydrostatic(e));
        auto p=pes.ConvergeAt(A, R);
        EXPECT_TRUE(p) << p.Error().details;
        return p ? p->energy : 0.0;
    };
    const double Em2=E(-2*h), Em1=E(-h), E0=E(0.0), Ep1=E(h), Ep2=E(2*h);
    const double d1=(Ep1-Em1)/(2*h), d2=(Ep2-Em2)/(4*h);
    const double c1=(Ep1-2*E0+Em1)/(h*h), c2=(Ep2-2*E0+Em2)/(4*h*h);
    std::cout<<std::setprecision(12)<<"[strain scan] E(-2h..2h) = "<<Em2<<" "<<Em1<<" "<<E0<<" "<<Ep1<<" "<<Ep2<<"  dE/de "<<d1<<" / "<<d2
             <<"  d2E/de2 "<<c1<<" / "<<c2<<std::endl;
    EXPECT_NEAR(d1, d2, 3e-5) << "dE/de at h and 2h";
    EXPECT_NEAR(c1, c2, 0.01*std::fabs(c1)) << "d2E/de2 at h and 2h";
    EXPECT_GT(c1, 0.0) << "Si's energy is convex in a hydrostatic strain near its reference";
}
