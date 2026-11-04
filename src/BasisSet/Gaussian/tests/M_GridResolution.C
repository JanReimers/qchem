// File: M_GridResolution.C  The grid builders' question asked of a Gaussian basis (V-BS-RASTER, user 2026-11-02): GridResolution on
// every molecular Gaussian block and its whole-basis reduction tBasisSet::MakeResolution.
//
// Each tail is checked against its DEFINITION, evaluated independently of the formula that produced it: the spectrum of the sharpest
// pair product (a Gaussian of exponent 2*alpha_max) at the returned cutoff, and the most diffuse primitive at the returned radius,
// must both sit at eps.
#include "gtest/gtest.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>
import qchem.BasisSet;
import qchem.BasisSet.GridResolution;
import qchem.BasisSet.AoShellSource;                  // the shells (their l and exponents) the margin is checked over
import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;   // GaussianSharpness: the exponent range the answer must agree with
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.Structure;
import qchem.Types;
import qchem.Blaze;   // range-for over rvec_t (the exported iterator operator!=)

using namespace qchem;
using namespace qchem::BasisSet::Gaussian;

namespace
{
//! Water over DZVP: H and O blocks with different exponent ranges, so the whole-basis reduction has something to reduce.
std::unique_ptr<BasisSet::Real_BS> Water(const Molecule& m) {return std::unique_ptr<BasisSet::Real_BS>(Factory(BasisSetData::DZVP, &m));}
Molecule MakeWater()
{
    Molecule m;
    m.Insert(Site{8, Vector3D(0.0, 0.0, 0.0)});
    m.Insert(Site{1, Vector3D(0.0, 1.431, 1.107)});
    m.Insert(Site{1, Vector3D(0.0,-1.431, 1.107)});
    return m;
}
}

//! Every block answers, and each answer IS its definition: the sharpest pair's spectrum exp(-G^2/4p), p = 2 alpha_max, is eps at
//! 1/2 G^2 = cutoff; the most diffuse pair's at coarseCutoff; the most diffuse primitive exp(-alpha_min r^2) is eps at the radius.
TEST(GridResolution, AGaussianBlocksTailsSitAtEpsByDefinition)
{
    const Molecule m=MakeWater();
    const auto bs=Water(m);
    size_t n=0;
    for (const auto* ibs : bs->Iterate())
    {
        const auto* g=dynamic_cast<const BasisSet::GridResolution*>(ibs);
        const auto* s=dynamic_cast<const GaussianSharpness*>(ibs);
        ASSERT_NE(g, nullptr) << "a molecular Gaussian block must answer the grid builders' question";
        ASSERT_NE(s, nullptr);
        for (double eps : {1e-4, 1e-8, 1e-12})
        {
            const BasisSet::Resolution r=g->MakeResolution(eps);
            auto spectrum=[](double p, double halfG2){ return std::exp(-2.0*halfG2/(4.0*p)); };   // G^2 = 2*(1/2 G^2)
            EXPECT_NEAR(spectrum(2.0*s->MaxExponent(), r.cutoff)/eps,       1.0, 1e-12) << "cutoff, eps=" << eps;
            EXPECT_NEAR(spectrum(2.0*s->MinExponent(), r.coarseCutoff)/eps, 1.0, 1e-12) << "coarseCutoff, eps=" << eps;
            EXPECT_NEAR(std::exp(-s->MinExponent()*r.radius*r.radius)/eps,   1.0, 1e-12) << "radius, eps=" << eps;
            EXPECT_LE(r.coarseCutoff, r.cutoff);
        }
        n++;
    }
    EXPECT_GT(n, 0u);
}

//! The whole basis reduces over its blocks: the finest cutoff, the coarsest coarse cutoff, the widest radius -- each attained by some
//! block, none exceeded.  And a tighter eps only ever asks the grid to reach further.
TEST(GridResolution, AWholeBasisReducesOverItsBlocksAndTightensWithEps)
{
    const Molecule m=MakeWater();
    const auto bs=Water(m);
    const double eps=1e-10;
    const BasisSet::Resolution R=bs->MakeResolution(eps);
    bool hitCut=false, hitCoarse=false, hitRad=false;
    for (const auto* ibs : bs->Iterate())
    {
        const BasisSet::Resolution b=dynamic_cast<const BasisSet::GridResolution&>(*ibs).MakeResolution(eps);
        EXPECT_LE(b.cutoff, R.cutoff);
        EXPECT_GE(b.coarseCutoff, R.coarseCutoff);
        EXPECT_LE(b.radius, R.radius);
        hitCut   |= b.cutoff==R.cutoff;
        hitCoarse|= b.coarseCutoff==R.coarseCutoff;
        hitRad   |= b.radius==R.radius;
    }
    EXPECT_TRUE(hitCut && hitCoarse && hitRad) << "a reduced value no block attains";
    const BasisSet::Resolution T=bs->MakeResolution(eps*1e-4);
    EXPECT_GT(T.cutoff, R.cutoff);
    EXPECT_GT(T.radius, R.radius);
    EXPECT_GT(T.coarseCutoff, R.coarseCutoff) << "the diffuse pair's spectrum, too, reaches further at a tighter eps";
}

//! eps outside (0,1) has no tail to answer: a caller-reachable choice, refused (not an assert).
TEST(GridResolution, AnEpsOutsideTheUnitIntervalIsRefused)
{
    const Molecule m=MakeWater();
    const auto bs=Water(m);
    EXPECT_THROW((void)bs->MakeResolution(0.0), std::invalid_argument);
    EXPECT_THROW((void)bs->MakeResolution(1.0), std::invalid_argument);
    EXPECT_THROW((void)bs->MakeResolution(-1e-6), std::invalid_argument);
}

//! The s-formula answer and an l>0 shell (user note, 2026-11-02): r^l e^{-a r^2}, measured against ITS peak (r0^2 = l/2a), sits at
//! eps (2eL/l)^{l/2} at the returned radius -- above eps, but within that stated margin, for every shell of every block.
TEST(GridResolution, AnLShellsTailStaysWithinItsPrefactorMargin)
{
    const Molecule m=MakeWater();
    const auto bs=Water(m);
    const double eps=1e-10, L=std::log(1.0/eps);
    size_t nl=0;
    for (const auto* ibs : bs->Iterate())
    {
        const BasisSet::Resolution r=dynamic_cast<const BasisSet::GridResolution&>(*ibs).MakeResolution(eps);
        const auto* src=dynamic_cast<const BasisSet::AoShellSource*>(ibs);
        ASSERT_NE(src, nullptr);
        for (const auto& sh : src->GetAoShells())
        {
            const int l=sh.rep->L();
            for (double a : sh.exponents)
            {
                const double rr=std::max(r.radius, std::sqrt(l/(2.0*a)));   // past the peak (the radius always is, for the diffuse end)
                const double r0=std::sqrt(l/(2.0*a));
                const double rel = l==0 ? std::exp(-a*rr*rr) : std::pow(rr/r0,l)*std::exp(-a*(rr*rr-r0*r0));
                const double margin = l==0 ? 1.0 : std::pow(2.0*std::exp(1.0)*L/l, 0.5*l);
                EXPECT_LE(rel, eps*margin*(1.0+1e-12)) << "l=" << l << " a=" << a;
                if (l>0) nl++;
            }
        }
    }
    EXPECT_GT(nl, 0u) << "DZVP water has p and d shells: the margin must have been exercised";
}
