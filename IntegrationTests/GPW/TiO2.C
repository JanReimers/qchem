// File: IntegrationTests/GPW/TiO2.C  rutile TiO2 -- the CONVERGENCE PANEL's little-group REFERENCE case (doc/SCFRobustnessPlan.md §4).
//
// littleGroupK measures the ground state's group on the reference orbitals to 1e-6, so its SCF must reach drho/N 1e-9
// (`a2e664f1`, deck::Vet refuses less).  The case is the SCF of `decks/campaigns/TiO2/tio2_k222_q222_littleGroupGate.json` with its
// response post-actions and its state file stripped: can the recipe deliver a 1e-9 reference within budget?  Full tier.
#include "gtest/gtest.h"
#include <fstream>
#include <string>
#include <nlohmann/json.hpp>
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>

import qchem.Tests.Panel;   // panel::Run, panel::BasisDataDir
import qchem.SolidCalculation;
import qchem.GeometryOptimizer;
import qchem.Matrix3D;                           // Matrix3D<double> (the rutile start cell)
import qchem.Lattice_3D;
import qchem.UnitCell;
import qchem.BasisSet.Gaussian.Point.Factory;   // Angular
import qchem.SCFParams;
import qchem.SCFAccelerator.Factory;            // SCFAccelerators::Type
import qchem.Types;
import qchem.Reporting;                         // report::FormatGuard

TEST(GPW_TiO2, k222_Pol_Smear_Kerker_Pulay_Converges_Panel_Long)
{
    std::ifstream in(qchem::panel::DeckDir().parent_path()/"campaigns"/"TiO2"/"tio2_k222_q222_littleGroupGate.json");
    ASSERT_TRUE(in) << "the TiO2 little-group gate deck";
    nlohmann::json deck=nlohmann::json::parse(in);
    deck.erase("postSCF"); deck.erase("state");
    deck["basis"]["file"]=(qchem::panel::BasisDataDir()/"valence_lowq_tm.bsd").string();   // the deck names /home/janr/Code/qchem6/...
    const auto o=qchem::panel::Run(deck, 49, "TiO2 k222 little-group reference to drho/N 1e-9 (regime 1, kT 5 mHa; 39 measured, 2026-10-29 panel baseline)");
    ASSERT_TRUE(o.run.converged) << o.run.summary;
}

//------------------------------------------------------------------------------------------------------------------------------------
// RUTILE RELAXES END TO END ON GPW (doc/OpenWork.md §2 "Real relaxations end to end", step 3; lane B 2026-11-20).  The cell (a, c) and
// the O 4f coordinate u of P4₂/mnm, on the WHOLE GPW force and stress (every term's pullback and virial, XC on the free Becke mesh).
// The oracle is CP2K CELL_OPT on the identical model (`~/Code/cp2k-runs/TiO2/tio2_rutile_k222_cellopt.{inp,out}`): the Γ-centred
// 2x2x2 mesh (eight TRIM points, listed), PADE (the GTH-PADE PPs' own functional, stated: the solid default is VWN5 LDA), the same
// basis file (valence_lowq_tm.bsd, spherical), GTH q4 Ti / q6 O, no smearing.  CP2K's minimum: a 8.209552, c 5.422002 bohr, u 0.301844,
// E -71.337813131 Ha (pressure 0.94 bar).
//   WHY k222, NOT Γ: at Γ alone CP2K's CELL_OPT left the rutile sites (a +7.6 %, c -6.7 %, 0.31 Ha down; `tio2_rutile_G_cellopt.out`).
//   WHY THE GROUP IS HELD BY THE MAP: rutile is UNSTABLE at k222 in this model too.  A free CP2K CELL_OPT from rutile + 0.002 fractional
// noise (`tio2_rutile_k222_free_cellopt.out`) falls to -71.388 Ha, and our free relaxation (2026-11-20, 55 steps, 7.7 h) reached CP2K's
// rutile minimum (E -71.3378282 at step ~14, 15 μHa from CP2K) and then slid to the same -71.3885 basin, until the cell stretched past
// the strain allowance.  CP2K's KEEP_SYMMETRY run stays rutile only because a symmetric start has symmetric forces.  So, as CP2K does,
// the DENSITY is free (all eight k-points, the plain Becke mesh) and the GEOMETRY keeps the space group: a CellWyckoffMap over the
// group detected on the reference (3 coordinates: a, c, u).  An IMPOSED session is not used: it refuses a Becke force (the folded
// mesh, row "Analytic forces"), and on the uniform raster its raw forces break the rutile mirror (OpenWork row "imposed IBZ forces").
namespace
{
//! Rutile at a 2x2x2 Γ-centred mesh at (a, c, u): the free geometry run and the map over the detected group (see above), handed to \c body.
void WithRutile(double a0, double c0, double u0, const char* label,
                const std::function<void(qchem::PotentialEnergySurface&, const qchem::CellWyckoffMap&)>& body)
{
    using namespace qchem;
    UnitCell cell(Matrix3D<double>(a0,0,0, 0,a0,0, 0,0,c0));
    for (const rvec3_t f : {rvec3_t(0,0,0), rvec3_t(0.5,0.5,0.5)}) cell.AddAtom(22, f);
    for (const rvec3_t f : {rvec3_t(u0,u0,0), rvec3_t(0.5+u0,0.5-u0,0.5), rvec3_t(0.5-u0,0.5+u0,0.5), rvec3_t(1-u0,1-u0,0)}) cell.AddAtom(8, f);
    SCFParams p; p.NMaxIter=100; p.MinΔρPerElectron=1e-7; p.MinΔE=std::nullopt; p.MinΔFD=std::nullopt; p.MinVirial=std::nullopt; p.MinFD=std::nullopt;
    RunOptions o; o.Nelec=32; o.species={{"Ti",4},{"O",6}}; o.label=label; o.scf=p;
    o.model=Hamiltonian::Model::PADE;   // the GTH-PADE PPs' own functional, and CP2K's &PADE
    o.solid.cellStrainAllowance=0.05;   // CP2K's minimum is 5.4 % (a) and 3.0 % (c) COMPRESSED from the start: the frozen sizes over-converge
    // STATED GRIDS (accuracy 0: the rasters at the basis floor, 72 Ha = 2 alpha_max; Becke nR 40 x L 29).  The pilot audit on this cell
    // needs 8-16 GB (measured 2026-11-20: 'pilot reference (x + 8) + Becke (nR=100, L=65)', estimate 10.4 GB) and refuses on resources on
    // this box (OpenWork G2 "Mesh sizing for geometry runs").
    o.accuracy.eps=0.0;
    o.solid.xcMeshSymmetric=false;   // a FREE geometry run: one plain Becke mesh at every point (SolidSession refuses a group-adapted one)
    GaussianBasisRecipe recipe{.file="valence_lowq_tm.bsd", .angular=BasisSet::Gaussian::Angular::Spherical};
    const auto setup=SolidSetup::Resolve(Lattice_3D(cell, ivec3_t(2,2,2)), recipe, o);
    auto run=MakeSolidGeometryRun(setup, {SCFStage{p, SCFAccelerators::Type::Auto}}, nullptr, SolidSession::Carry::AO, /*relaxCell*/true);
    // The free session's own map keeps only the identity; hold rutile's group instead (the detection tolerance, as the factory states it).
    const PESGeometry ref=run.surface->Reference();
    const CellWyckoffMap map(ref, setup->lattice.GetSpaceGroup(), ElasticModel{}, 1e-4*(*ref.cell)(1,1));
    ASSERT_EQ(map.Dimension(), 3u) << "rutile's free coordinates are a, c and u";
    body(*run.surface, map);
}

//! Rutile relaxed from the start (a, c, u) to its minimum, checked against CP2K's (see above).
void RelaxRutileToCP2KsCell(double a0, double c0, double u0, const char* label)
{
    using namespace qchem;
    WithRutile(a0, c0, u0, label, [](PotentialEnergySurface& pes, const CellWyckoffMap& map)
    {
        auto r=GeometryOptimizer(pes, map, OptimizerOptions::FromEps(1e-6)).Run();
        ASSERT_TRUE(bool(r)) << r.Error().details;
        report::FormatGuard keep(std::cout);
        std::cout<<"[relax] rutile: "<<r->Steps()<<" steps, "<<r->history.size()<<" evaluations"<<std::endl;
        for (const auto& h : r->history)
            std::cout<<"[relax]   "<<std::setw(3)<<h.index<<"  E "<<std::setprecision(10)<<std::fixed<<h.value<<std::defaultfloat
                     <<"  max|g| "<<std::setprecision(3)<<std::scientific<<h.maxGradient<<"  |s| "<<h.stepLength
                     <<std::defaultfloat<<(h.accepted ? "" : "  rejected")<<std::endl;
        const rmat3d_t& A=*r->geometry.cell;
        const double a=A(1,1), c=A(3,3), u=r->geometry.R[2].x/a;
        std::cout<<std::setprecision(10)<<"[relax] rutile a "<<a<<"  c "<<c<<"  u "<<u<<"  E "<<r->value
                 <<"  (CP2K a 8.209552, c 5.422002, u 0.301844, E -71.337813131)"<<std::endl;
        EXPECT_NEAR(A(1,2), 0.0, 1e-8); EXPECT_NEAR(A(2,2), a, 1e-8);   // still tetragonal (the map holds it)
        // MEASURED 2026-11-20, cold from the Materials cell (6 steps, 55 min on 4 threads): a 8.208446, c 5.431751, u 0.301768, E -71.3378274.  a and u sit on CP2K's;
        // c is 0.18 % long, and that is the STATED Becke mesh: at CP2K's geometry the (40, 29) mesh gives σzz = -3.9e-5 Ha/bohr³ (~1.1 GPa,
        // ~0.2 % of c), the (80, 41) mesh -8.4e-7 (σxx 1.7e-6), i.e. CP2K's cell is our minimum on the finer mesh.  So c's tolerance is the
        // stated mesh's shift; an eps-sized mesh closes it (OpenWork G2 "Mesh sizing for geometry runs").  The energy's 14 μHa is partly
        // the stated mesh's luck: (80, 41) puts the same point 152 μHa above CP2K's.
        EXPECT_NEAR(a, 8.209552, 3e-3);
        EXPECT_NEAR(c, 5.422002, 1.5e-2);
        EXPECT_NEAR(u, 0.301844, 2e-4);
        EXPECT_NEAR(r->value, -71.337813131, 1e-4);
    });
}
} // namespace

// COLD, from the Materials cell (a 8.681591, c 5.591509, u 0.3057): the PERIODIC tier (55 min on 4 threads; root CMakeLists
// QCHEM_PERIODIC), its cheap twin below carries the claim in the long tier.
TEST(GPW_TiO2, k222_RelaxesToCP2KsCell_Long) { RelaxRutileToCP2KsCell(8.681591, 5.591509, 0.3057, "rutile relax"); }

// THE TWIN: ONE point, at CP2K's minimum, is stationary on our surface (the gradient over a, c, u that the optimizer steps on).
// A relaxation STARTED there is no twin: it pays the stated mesh's own c shift in 6 evaluations (measured 2026-11-20: 3300 s, the
// cold run's 3320 s).  MEASURED 2026-11-20 at this point: max|g| 2.0e-3 (the stated mesh's σzz, see above), E -71.3378876, 75 μHa
// below CP2K's; the cold start's max|g| is 1.1e-1, and the optimizer stops at 4.5e-4.  So the gate holds the stationary point to the
// stated mesh's shift: a force or stress error beyond that mesh's goes red.
TEST(GPW_TiO2, k222_CP2KsCellIsStationary_Long)
{
    using namespace qchem;
    WithRutile(8.209552, 5.422002, 0.301844, "rutile twin", [](PotentialEnergySurface& pes, const CellWyckoffMap& map)
    {
        const rvec_t q=map.Start();
        const PESGeometry g=map.At(q);
        auto point=pes.ConvergeAt(g.cell, g.R, map.PointRequest());
        ASSERT_TRUE(bool(point)) << point.Error().details;
        const rvec_t grad=map.Project(q, map.Gradient(q, *point));
        double gmax=0.0;
        for (size_t i=0;i<grad.size();++i) gmax=std::max(gmax, std::fabs(grad[i]));
        report::FormatGuard keep(std::cout);
        std::cout<<std::setprecision(10)<<"[stationary] rutile at CP2K's cell: E "<<point->energy<<" (CP2K -71.337813131)  max|g| "
                 <<std::setprecision(3)<<std::scientific<<gmax<<"  g";
        for (size_t i=0;i<grad.size();++i) std::cout<<" "<<grad[i];
        std::cout<<std::defaultfloat<<std::endl;
        EXPECT_LT(gmax, 3e-3);
        // KNOWN RED (root CMakeLists QCHEM_KNOWN_RED, controller ruling 2026-11-24): the absolute pair routing (2026-11-24) removed 68 uHa of
        // routing error, E -71.3378876 -> -71.33795527 (routing converged to 1e-8 from ln(1/eps) 16 to 32), 142 uHa below CP2K.  This line
        // is MESH-limited: the stated (40,29) Becke mesh was luck, (80,41) puts the point 152 uHa above CP2K.  The gate stays; the claim
        // moves to CP2K's own quadrature (a Uni twin), then this comes off the list (OpenWork §3 "TiO2 k222 energy vs CP2K is MESH-limited").
        EXPECT_NEAR(point->energy, -71.337813131, 1e-4);
    });
}
