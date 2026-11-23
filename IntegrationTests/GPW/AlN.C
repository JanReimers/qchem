// File: IntegrationTests/GPW/AlN.C  wurtzite AlN -- the geometry easy case with a free internal coordinate (u) and two cell lengths.
//
// doc/OpenWork.md §2 "Real relaxations end to end", the easy cases (controller plan 2026-11-20 (c)).  The question these gates answer is
// the user's (1) of 2026-11-22: "Does CP2K and SOLID^2 converge to same (wrong) geometry for a given Exc,{k},orbital basis etc."  Whether
// that geometry is close to experiment is the separate MODEL question (2), OpenWork §3 "AlN on the low-q basis relaxes far from experiment".
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>

import qchem.SolidCalculation;
import qchem.GeometryOptimizer;
import qchem.Lattice_3D;
import qchem.UnitCell;
import qchem.SCFParams;
import qchem.Mesh.XCPolicy;                     // BeckeXCParams (the stated XC mesh)
import qchem.SCFAccelerator.Factory;            // SCFAccelerators::Type
import qchem.Types;
import qchem.Reporting;                         // report::FormatGuard

//------------------------------------------------------------------------------------------------------------------------------------
// WURTZITE AlN ON GPW: the cell (a, c) and the N 2b coordinate u of P6₃mc, on the whole GPW force and stress (XC on the free Becke mesh).
// The oracle is CP2K CELL_OPT on the identical model (`~/Code/cp2k-runs/AlN/cellopt/aln_k222_free_cellopt.{inp,out}`): the Γ-centred
// 2x2x2 mesh (the 8 TRIM points, listed), PADE, OUR basis valence_lowq_aln.bsd (ALN-BASIS by scripts/bsd2cp2k; s and p only, so
// spherical = Cartesian), GTH-PADE q3 Al / q5 N, CUTOFF 800 Ry (1200 Ry moves E 1.0 μHa at the minimum).  CP2K's FREE run (KEEP_SYMMETRY
// F, from a perturbed 120° wurtzite: u 0.38, the cell skewed 1e-3) stays wurtzite: a 5.9153946, c 8.9044783 bohr, u 0.401717,
// E -24.2074417 Ha (37 steps; its symmetric run reached the same E to 2e-7).  Al-N 3 x 1.866 + 1.893 Å, 4-coordinated, direct Γ gap 3.71 eV.
//   ⚠ CP2K's `SYMMETRY HEXAGONAL` rewrites a 120° cell to γ = 60° and re-reads SCALED (1/3, 2/3) sites in it: a 3-coordinated start that
// KEEP_SYMMETRY then holds.  The oracle above is the free run, which never met that keyword.
//   As for rutile and NaF the DENSITY is free (all 8 k-points, the plain Becke mesh) and the GEOMETRY keeps the group: a CellWyckoffMap
// over the group detected on the reference.  P6₃mc is POLAR: the map keeps the z of BOTH 2b orbits (4 coordinates: a, c, z_Al, z_N), so
// the uniform z translation is one of them, with zero gradient by construction; u = z_N - z_Al.
//   ⚠ k222 is a POOR MODEL for AlN (the run WARNS "k mesh COARSE", 0.62 > 0.5 bohr⁻¹).  QE pw.x with these same GTH PPs (gth2upf,
// 250 Ry, `~/Code/qe-runs/AlN/aln_vcr_k{2,4}`): at k222 wurtzite collapses to the layered u = 0.5 phase (c/a 1.186, as CP2K DZVP),
// at k444 it holds at c/a 1.598, u 0.3825 (the LDA literature).  This basis (no d) holds wurtzite at k222 only as a trap.  So the TWIN
// below is a same-model claim (our stationary point = CP2K's, controller 2026-11-22), not a structure claim, and there is NO cold
// relaxation at k222: that gate goes to k444 on DZVP after the geometry hold (OpenWork §3 AlN row).
namespace
{
//! Wurtzite AlN at (a, c, u) on the k222 mesh (Al q3 / N q5, PADE, FREE density) with the P6₃mc map over its group; runs \a body on the
//! run's surface and the map.
void WithAlN(double a0, double c0, double u0, const char* label,
             const std::function<void(qchem::PotentialEnergySurface&, const qchem::CellWyckoffMap&)>& body)
{
    using namespace qchem;
    UnitCell cell(a0, a0, c0, 90.0, 90.0, 120.0);   // a along x, b at 120° in the xy plane: (1/3, 2/3) is a 3-fold site
    cell.AddAtom(13, rvec3_t(1.0/3.0, 2.0/3.0, 0.0));
    cell.AddAtom(13, rvec3_t(2.0/3.0, 1.0/3.0, 0.5));
    cell.AddAtom( 7, rvec3_t(1.0/3.0, 2.0/3.0, u0));
    cell.AddAtom( 7, rvec3_t(2.0/3.0, 1.0/3.0, 0.5+u0));
    SCFParams p; p.NMaxIter=100; p.MinΔρPerElectron=1e-7; p.MinΔE=std::nullopt; p.MinΔFD=std::nullopt; p.MinVirial=std::nullopt; p.MinFD=std::nullopt;
    RunOptions o; o.Nelec=16; o.species={{"Al",3},{"N",5}}; o.label=label; o.scf=p;
    o.model=Hamiltonian::Model::PADE;   // the GTH-PADE PPs' own functional, and CP2K's &PADE
    o.solid.cellStrainAllowance=0.07;   // CP2K's minimum is 5.4 % compressed in c from the experimental start
    o.solid.xcMeshSymmetric=false;      // a FREE geometry run: one plain Becke mesh at every point
    // STATED GRIDS: the default accuracy's raster pilot at CP2K's cell (2026-11-22, its "a deck restates these grids as" line: the ball
    // doubled to cutoffFactor 16, Becke (100, 53); its audit reference needed 3.05 GB, "it may not fit"), so the twin and the cold run skip it.
    o.accuracy.eps=0.0;
    o.solid.cutoffFactor=16.0;
    o.solid.tolerances.rasterX=o.solid.tolerances.xcRasterX=19.35139721329466;
    o.solid.tolerances.vlocRasterX=21.583068172270572;
    o.solid.tolerances.vlocEps=1e-6;
    o.xcMesh=qcMesh::BeckeXCParams(100, 2.0, 53);
    const auto setup=SolidSetup::Resolve(Lattice_3D(cell, ivec3_t(2,2,2)), {.file="valence_lowq_aln.bsd"}, o);
    auto run=MakeSolidGeometryRun(setup, {SCFStage{p, SCFAccelerators::Type::Auto}}, nullptr, SolidSession::Carry::AO, /*relaxCell*/true);
    const PESGeometry ref=run.surface->Reference();
    const CellWyckoffMap map(ref, setup->lattice.GetSpaceGroup(), ElasticModel{}, 1e-4*a0);
    ASSERT_EQ(map.Dimension(), 4u) << "wurtzite's free coordinates: a, c, and the z of each 2b orbit (u is their difference)";
    body(*run.surface, map);
}

} // namespace

// THE TWIN: CP2K's minimum is a stationary point of OUR surface -- one point, the map gradient (a, c, u) at CP2K's cell.
// Measured (lane B 2026-11-22, 4 threads, stated grids): E -24.20743513 (6.6 μHa above CP2K), |g| ≤ 1.8e-5 Ha/coordinate, 3205 s -- so it
// runs in the PERIODIC tier (root CMakeLists QCHEM_PERIODIC), not the long one.
TEST(GPW_AlN, k222_CP2KsCellIsStationary_Long)
{
    using namespace qchem;
    WithAlN(5.9153946, 8.9044783, 0.401717, "AlN twin", [](PotentialEnergySurface& pes, const CellWyckoffMap& map)
    {
        const rvec_t q=map.Start();
        const PESGeometry g=map.At(q);
        auto point=pes.ConvergeAt(g.cell, g.R, map.PointRequest());
        ASSERT_TRUE(bool(point)) << point.Error().details;
        const rvec_t grad=map.Project(q, map.Gradient(q, *point));
        report::FormatGuard keep(std::cout);
        std::cout<<std::setprecision(10)<<"[stationary] AlN at CP2K's cell: E "<<point->energy<<" (CP2K -24.2074417)  g";
        double gmax=0.0;
        for (size_t i=0; i<grad.size(); ++i) { std::cout<<" "<<std::setprecision(3)<<std::scientific<<grad[i]; gmax=std::max(gmax, std::fabs(grad[i])); }
        std::cout<<std::defaultfloat<<std::endl;
        EXPECT_LT(gmax, 1e-3);
        EXPECT_NEAR(point->energy, -24.2074417, 1e-4);
    });
}
