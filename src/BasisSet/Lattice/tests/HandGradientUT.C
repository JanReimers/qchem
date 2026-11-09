// File: HandGradientUT.C  THE HAND-GRADIENT AUDIT, lattice families (lane B 2026-11-09, ahead of the user's Enzyme decision): every
// lattice basis's hand Gradient(r) against the central difference of its own operator()(r), every function and axis.  Plane waves, LAPW
// and APW (their interstitial plane waves), the GPW Bloch sum of Cartesian Gaussians, and the spherical lattice view's blocks (the
// harmonic congruence under GPW).  Each at Γ and off it where the family has a k.  Step h = 1e-5: a plane wave's h²/6 |K|³ is ≲ 1e-9
// for the cutoffs here, a Gaussian's as M_GeometryDerivatives GATE 2, so 1e-7 (relative above 1).  The largest derivative probed must be
// 0.01 or more (10⁵ × the tolerance; LAPW and APW reach only ~0.07, their plane waves carry 1/√V with V = 512), so a sign or a missing i
// would show.
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <string>
#include <vector>
import qchem.BasisSet;
import qchem.BasisSet.Orbital_1E_IBS;
import qchem.BasisSet.IrrepBasisSet;          // APW_IBS is an IrrepBasisSet (its VectorFunction face)
import qchem.BasisSet.PlaneWave.PlaneWave_IBS;
import qchem.BasisSet.Lattice.LAPW_IBS;
import qchem.BasisSet.Lattice.APW_IBS;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.BasisSet.Gaussian.Lattice.GPW_IBS;
import qchem.BasisSet.Gaussian.Lattice.SphericalLatticeView;
import qchem.Structure;
import qchem.UnitCell;
import qchem.Lattice_3D;
import qchem.VectorFunction;
import qchem.Types;
import qchem.Blaze;
using namespace qchem;

namespace
{
double Mag(double x) {return std::fabs(x);}
double Mag(const dcmplx& x) {return std::abs(x);}
template <class T> T Comp(const vec3_t<T>& v, int k) {return k==0 ? v.x : k==1 ? v.y : v.z;}
//! The gate: returns the largest |FD| seen, so the caller can check it means something.
template <class T, class F> double CheckHandGradient(const F& f, const std::vector<rvec3_t>& pts, const std::string& who)
{
    const double h=1e-5;
    double largest=0.0;
    for (const rvec3_t& r : pts)
    {
        const vec3vec_t<T> g=f.Gradient(r);
        const vec_t<T> v=f(r);
        EXPECT_EQ(g.size(), v.size()) << who << ": one gradient per function";
        if (g.size()!=v.size()) return largest;
        for (int k=0;k<3;k++)
        {
            rvec3_t p=r, m=r; (k==0 ? p.x : k==1 ? p.y : p.z)+=h; (k==0 ? m.x : k==1 ? m.y : m.z)-=h;
            const vec_t<T> vp=f(p), vm=f(m);
            for (size_t i=0;i<v.size();i++)
            {
                const T fd=(vp[i]-vm[i])/(2*h);
                EXPECT_LT(Mag(Comp(g[i],k)-fd), 1e-7*std::max(1.0,Mag(fd))) << who << ": function " << i << " axis " << k << " at " << r
                    << ", hand " << Comp(g[i],k) << " vs FD " << fd;
                largest=std::max(largest, Mag(fd));
            }
        }
    }
    return largest;
}
} //anon

TEST(HandGradient, EveryLatticeFamilysGradientIsTheFiniteDifferenceOfItsValues)
{
    const std::vector<rvec3_t> pts{rvec3_t(0.31,-0.42,0.27), rvec3_t(1.7,2.3,-0.9), rvec3_t(-2.6,0.8,3.1)};
    {
        const double a=8.0; UnitCell cell(a); const ivec3_t N(4,4,4); Lattice_3D lat(cell,N);
        for (const ivec3_t k : {ivec3_t(0,0,0), ivec3_t(1,0,0)})
        {
            const std::string at=" k=("+std::to_string(k.x)+","+std::to_string(k.y)+","+std::to_string(k.z)+")";
            BasisSet::PlaneWave::PlaneWave_IBS pw(lat.Reciprocal(), N, k, 4.0);
            EXPECT_GT(CheckHandGradient<dcmplx>(pw, pts, "PlaneWave_IBS"+at), 0.01);
            BasisSet::Lattice::LAPW_IBS lapw(lat.Reciprocal(), N, k, 2.0, 2.0, 8, 0.2);
            EXPECT_GT(CheckHandGradient<dcmplx>(lapw, pts, "LAPW_IBS"+at), 0.01);
            BasisSet::Lattice::APW_IBS apw(lat.Reciprocal(), N, k, 2.0, 2.0, 8);
            EXPECT_GT(CheckHandGradient<dcmplx>(apw, pts, "APW_IBS"+at), 0.01);
        }
    }
    {
        const double a=8.40;
        FCCUnitCell cell(a);
        cell.AddAtom(25, {0,0,0});
        std::shared_ptr<const BasisSet::Real_BS> cart(
            BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::VALENCE_LOWQ_SR, &cell, BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
        auto view=BasisSet::Gaussian::PG_Spherical::MakeSphericalLatticeView(cart);
        for (auto b : view->Iterate<BasisSet::Real_OIBS>())
            EXPECT_GT(CheckHandGradient<double>(*b, pts, "SphericalLatticeView "+b->Name()), 0.01);
        for (const ivec3_t k : {ivec3_t(0,0,0), ivec3_t(1,0,0)})
            for (const auto& [name, src] : {std::pair<std::string,std::shared_ptr<const BasisSet::Real_BS>>{"Cartesian", cart}, {"spherical view", view}})
            {
                BasisSet::Gaussian::GPW_IBS gpw(cell, ivec3_t(4,4,4), k, src, /*densityEcut=*/0.0);
                const BasisSet::Complex_OIBS& g=gpw;
                EXPECT_GT(CheckHandGradient<dcmplx>(g, pts, "GPW_IBS "+name+" k=("+std::to_string(k.x)+",0,0)"), 0.01);
            }
    }
}
