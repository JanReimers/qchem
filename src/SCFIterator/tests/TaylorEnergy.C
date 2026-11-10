// File: SCFIterator/tests/TaylorEnergy.C  The first-order energy E_T of a D-mixing step is SECOND ORDER (controller ruling 2026-11-09).
//
// A D-mixer's next Fock is built from its mixed D_in, not from D_out, so E_KS[D_out] would cost a second density refresh (a second
// collocation per iteration on GPW: +33 % per iteration on MnO PulD, measured 2026-11-09).  The trace reports instead
//     E_T = E[D_in] + Re Tr(F[D_in](D_out - D_in)),     E_KS[D_out] - E_T = O(|D_out - D_in|^2),
// and refreshes ONCE at the end for the exact E_KS[D_out].  This file is the claim's gate: with the diagnostic energy_taylor on (it
// prints the exact E_KS[D_out] beside E_T every iteration, changing no number), the gap must fall as the SQUARE of the step, on a
// GPW solid (Si k222, Pulay on D: at Gamma the step direction rotates between iterations, so gap/step^2 wanders 0.0006-0.056
// and a log-log fit is not a clean line; at k222 it holds at 0.05-0.058 over eight iterations) and on a POLARIZED molecule through exact exchange (O2 HF triplet: spin-resolved F and D, and
// the whole-system HF terms), and the run's final energy must be the exact E_KS[D_out].
#include "gtest/gtest.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
import qchem.Diagnostics;
import qchem.SolidCalculation;
import qchem.MoleculeCalculation;
import qchem.Structure;
import qchem.Lattice_3D;
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.Types;
import qchem.SCFIterator;
using namespace qchem;

namespace {

std::string Captured(const std::function<void()>& body)
{
    std::ostringstream sink;
    std::streambuf* saved = std::cout.rdbuf(sink.rdbuf());
    try { body(); } catch (...) { std::cout.rdbuf(saved); throw; }
    std::cout.rdbuf(saved);
    return sink.str();
}

struct Point { double gap, step; };   // step = |D_out - D_in|_F, the Taylor expansion's own step
struct Parsed { std::vector<Point> pts; double finalE=NAN, finalGap=NAN; };

Parsed Parse(const std::string& out)
{
    Parsed p;
    std::istringstream is(out);
    for (std::string line; std::getline(is, line); )
    {
        if (line.rfind("[energy_taylor]",0)==0)
        {
            std::istringstream ls(line);
            std::string tag, k; int it; double ET, EKS, gap, step;
            ls >> tag >> k >> it >> k >> ET >> k >> EKS >> k >> gap >> k >> step;
            p.pts.push_back({gap, step});
        }
        if (auto i=line.find("[energy] final E_KS[D_out] = "); i!=std::string::npos)
        {
            p.finalE=std::stod(line.substr(i+29));
            p.finalGap=std::stod(line.substr(line.find("|E_T - E_KS| = ")+15));
        }
    }
    return p;
}

//! The least-squares slope of log(gap) on log(step), over the points whose gap sits between roundoff and the pre-asymptotic start.
double Slope(const std::vector<Point>& pts, double gapMin, double gapMax, size_t& n)
{
    double sx=0, sy=0, sxx=0, sxy=0; n=0;
    for (const Point& q : pts)
        if (q.gap>gapMin && q.gap<gapMax && q.step>0)
        {
            const double x=std::log(q.step), y=std::log(q.gap);
            sx+=x; sy+=y; sxx+=x*x; sxy+=x*y; ++n;
        }
    return (n*sxy-sx*sy)/(n*sxx-sx*sx);
}

} // anon

TEST(TaylorEnergy, GPW_Si_k222_PulayOnD_GapIsSecondOrder)
{
    Diagnostics::Scoped on("energy_taylor");
    const Parsed p=Parse(Captured([]{
        FCCUnitCell cell(10.26);
        cell.AddAtom(14, {0,0,0});
        cell.AddAtom(14, {0.25,0.25,0.25});
        Lattice_3D lat(cell, ivec3_t(2,2,2));
        auto mol = std::shared_ptr<const BasisSet::Real_BS>(
            BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::SIPP_SR, &cell,
                                        BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
        SCFParams par;
        par.NMaxIter=40; par.MinΔρPerElectron=1e-8; par.MinΔE=std::nullopt; par.MinΔFD=std::nullopt;
        par.MinVirial=std::nullopt; par.MinFD=std::nullopt; par.StartingRelaxRo=0.3; par.MergeTol=1e-4; par.Verbose=true;
        SolidCalculation calc(SolidSetup::Resolve(lat, mol, {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0}), par);
        (void)calc.Result();
    }));
    ASSERT_GE(p.pts.size(), 6u) << "too few D-mixing iterations printed their energy_taylor line";
    size_t n=0;
    // The window: above the 1e-10 Ha floor where E's own roundoff and the screening tolerances take over (measured: gap/step^2 holds at
    // 0.05-0.058 from step 0.26 down to 5e-4, then the gap sits at 1e-12..1e-15 Ha).
    const double s=Slope(p.pts, 1e-10, 1e-2, n);
    std::cerr << "[TaylorEnergy] Si k222 PulD: slope " << s << " over " << n << " iterations" << std::endl;
    ASSERT_GE(n, 4u) << "too few iterations between roundoff and the pre-asymptotic start to fit a slope";
    EXPECT_NEAR(s, 2.0, 0.3) << "|E_T - E_KS[D_out]| should fall as the square of the step (" << n << " points)";
    for (const Point& q : p.pts)   // DIRECTION-FREE: a first-order error (gap ~ c*step) breaks this once step < c; a remainder cannot.
        // Above the roundoff floor of E itself (1e-13 |E|: O2's 150 Ha sits at 2e-13 Ha once converged, ~8 ulp).
        if (q.step<1e-2 && q.gap>1e-13*std::fabs(p.finalE)) EXPECT_LE(q.gap, q.step*q.step) << "gap " << q.gap << " at step " << q.step << " is not second order";
    ASSERT_FALSE(std::isnan(p.finalE)) << "the run did not end on the exact energy refresh";
    EXPECT_LT(p.finalGap, 1e-10) << "the converged E_T should equal E_KS[D_out] to second order in a converged step";
}

TEST(TaylorEnergy, O2_HF_Triplet_Polarized_GapIsSecondOrder)
{
    Diagnostics::Scoped on("energy_taylor");
    const Parsed p=Parse(Captured([]{
        Molecule o2;
        o2.Insert(Site{8, Vector3D<double>(0, 0,  1.141)});
        o2.Insert(Site{8, Vector3D<double>(0, 0, -1.141)});
        MoleculeCalculation calc(MoleculeSetup::Resolve(o2, {.basis = "dzvp", .multiplicity = 3}),
                         {.NMaxIter = 80, .MinΔρPerElectron = 1e-8, .MinΔFD = std::nullopt, .MinVirial = std::nullopt,
                          .MinFD = 1e-7, .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = true});
        (void)calc.Result();
    }));
    ASSERT_GE(p.pts.size(), 6u);
    size_t n=0;
    const double s=Slope(p.pts, 1e-10, 1e-2, n);
    std::cerr << "[TaylorEnergy] O2 HF triplet: slope " << s << " over " << n << " iterations" << std::endl;
    ASSERT_GE(n, 4u);
    EXPECT_NEAR(s, 2.0, 0.3) << "the spin-resolved |E_T - E_KS[D_out]| should fall as the square of the step (" << n << " points)";
    for (const Point& q : p.pts)   // DIRECTION-FREE: a first-order error (gap ~ c*step) breaks this once step < c; a remainder cannot.
        // Above the roundoff floor of E itself (1e-13 |E|: O2's 150 Ha sits at 2e-13 Ha once converged, ~8 ulp).
        if (q.step<1e-2 && q.gap>1e-13*std::fabs(p.finalE)) EXPECT_LE(q.gap, q.step*q.step) << "gap " << q.gap << " at step " << q.step << " is not second order";
    ASSERT_FALSE(std::isnan(p.finalE));
    EXPECT_LT(p.finalGap, 1e-10);
}
