// File: Calculation/tests/SolidRelaxation.C  REAL SOLID RELAXATIONS (doc/OpenWork.md §2 "Real relaxations end to end", step 2): the
// geometry optimizer over a SolidSession on PLANE WAVES, whose whole force is gated (every PW term's force at fixed D, and the session's
// force against its own re-converged FD: ForceFacade.C).  The internal coordinates only, at a fixed cell: a FREE setup, so the map is
// P1 (q = all 3N Cartesian coordinates, the rigid translations held) and nothing but the forces brings the atoms home.  The oracle is
// SYMMETRY: Si₂ distorted off its diamond site and F off its rocksalt site each relax back to the site the space group fixes, where
// every force vanishes by symmetry.  Si₂ only, and no cell: plane waves are Γ-only so far, and at Γ NaF (Na q1) and Si₈ are UNSTABLE on
// their symmetric sites (F walks to the tetrahedral site, 77 mHa lower; Si₈ leaves diamond, QE's own relax too, from the same energy to
// 2e-7 Ha), so neither can carry a "returns to its site" claim (doc/GeometryOptimizationPlan.md §15).
#include <gtest/gtest.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
import qchem.SolidCalculation;
import qchem.GeometryOptimizer;
import qchem.Lattice_3D;
import qchem.UnitCell;
import qchem.BasisSet.Gaussian.Point.Factory;   // BasisSetData (the Gaussian recipe a PW setup must still name)
import qchem.SCFParams;
import qchem.SCFAccelerator.Factory;            // SCFAccelerators::Type
import qchem.Types;
import qchem.Reporting;                         // report::FormatGuard

using namespace qchem;

namespace
{
SCFParams Params()
{
    // Δρ/N 1e-7 (the PW force gates' gate): the energy is settled to 1e-15 by then, and Δρ/N floats at 2-4e-8 below it on these cells
    // (measured, Si₈: 1e-8 was never met).
    SCFParams p; p.NMaxIter=100; p.MinΔρPerElectron=1e-7; p.MinΔE=std::nullopt; p.MinΔFD=std::nullopt; p.MinVirial=std::nullopt; p.MinFD=std::nullopt;
    return p;
}

//! Relax the free setup's atoms at its fixed cell; print the history; hand back the result (asserting it converged).
OptimizedGeometry Relax(std::shared_ptr<const SolidSetup> setup, const std::string& what)
{
    auto run=MakeSolidGeometryRun(setup, {SCFStage{Params(), SCFAccelerators::Type::Auto}});
    OptimizerOptions o=OptimizerOptions::FromEps(1e-6);   // max|g| 4.5e-5 Ha/bohr
    auto r=GeometryOptimizer(*run.surface, *run.map, o).Run();
    EXPECT_TRUE(bool(r)) << what << ": " << (r ? std::string() : r.Error().details);
    if (!r) return {};
    report::FormatGuard keep(std::cout);
    std::cout<<"[relax] "<<what<<": "<<r->Steps()<<" steps, "<<r->history.size()<<" evaluations"<<std::endl;
    for (const auto& h : r->history)
        std::cout<<"[relax]   "<<std::setw(3)<<h.index<<"  E "<<std::setprecision(10)<<std::fixed<<h.value<<std::defaultfloat
                 <<"  max|g| "<<std::setprecision(3)<<std::scientific<<h.maxGradient<<"  |s| "<<h.stepLength
                 <<std::defaultfloat<<(h.accepted ? "" : "  rejected")<<std::endl;
    return r.TakeValue();
}

//! The fcc cell matrix of cubic constant \a a (its columns the primitive vectors; symmetric).
Matrix3D<double> Fcc(double a) {const double h=0.5*a; return Matrix3D<double>(0.0,h,h,  h,0.0,h,  h,h,0.0);}
} // namespace

// Si₂ on plane waves (Ecut 4 Ha, density ball 20 Ha, Γ), the second atom started at fractional (0.27, 0.24, 0.26), relaxed back to the
// diamond site (1/4, 1/4, 1/4) relative to the first.
TEST(SolidRelaxation, OnPlaneWavesSiReturnsToTheDiamondSite)
{
    const double a=10.26;
    UnitCell cell(Fcc(a));
    cell.AddAtom(14, rvec3_t(0.0,0.0,0.0));
    cell.AddAtom(14, rvec3_t(0.27,0.24,0.26));
    SolidCalcOptions o; o.Nelec=8; o.species={{"Si",4}}; o.densityEcut=20.0; o.label="pw Si relax";
    o.orbitalBasis=SolidCalcOptions::OrbitalBasis::PlaneWave; o.planeWaveEcut=4.0;
    const auto setup=SolidSetup::Resolve(Lattice_3D(cell, ivec3_t(1,1,1)), {.data=BasisSet::Gaussian::BasisSetData::SIPP_SR}, o);
    const auto r=Relax(setup, "Si2 PW");
    ASSERT_EQ(r.geometry.R.size(), 2u);
    const rvec3_t d=r.geometry.R[1]-r.geometry.R[0], ideal=Fcc(a)*rvec3_t(0.25,0.25,0.25);
    std::cout<<"[relax]   R1-R0 "<<d<<"  diamond "<<ideal<<"  off by "<<norm(d-ideal)<<" bohr"<<std::endl;
    EXPECT_LT(norm(d-ideal), 1e-3) << "back on the diamond site";
}
