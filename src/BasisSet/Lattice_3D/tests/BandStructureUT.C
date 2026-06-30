// file: BandStructureUT.C  k-point sweep / band structure across the Brillouin zone.
//
// Validates the shared k-layer (KPath + SolveBands) against the exact empty-lattice free-electron band
// structure along Gamma-X-M-Gamma, for both lineages (PlaneWave and LAPW share one SolveBands).

#include <vector>
#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include "gtest/gtest.h"

import qchem.BasisSet.Lattice_3D.BandStructure;   // SolveBands, KPath
import qchem.BasisSet.Lattice_3D.PlaneWave_IBS;
import qchem.BasisSet.Lattice_3D.LAPW_IBS;
import qchem.Lattice_3D;                          // UnitCell, Lattice_3D, ReciprocalLattice
import qchem.Types;
import qchem.Math;                                // Pi
using namespace qchem;

using BasisSet::Lattice_3D::SolveBands;
using BasisSet::Lattice_3D::KPath;
using BasisSet::Lattice_3D::PlaneWave_IBS;
using BasisSet::Lattice_3D::LAPW_IBS;

namespace
{
// Lowest free-electron energies 1/2|k+G|^2 (cubic, B=(2pi/a)I) at k = kIndex/N, ascending.
std::vector<double> FreeEnergies(double a, ivec3_t N, ivec3_t kIndex, int mr)
{
    double b=2*Pi/a;
    rvec3_t kf(kIndex.x/double(N.x), kIndex.y/double(N.y), kIndex.z/double(N.z));
    std::vector<double> e;
    for (int mx=-mr;mx<=mr;mx++) for (int my=-mr;my<=mr;my++) for (int mz=-mr;mz<=mr;mz++)
    {
        rvec3_t q((kf.x+mx)*b,(kf.y+my)*b,(kf.z+mz)*b);
        e.push_back(0.5*(q.x*q.x+q.y*q.y+q.z*q.z));
    }
    std::sort(e.begin(),e.end());
    return e;
}

std::vector<double> ascending(const rvec_t& e)
{
    std::vector<double> v; for (size_t i=0;i<e.size();i++) v.push_back(e[i]);
    std::sort(v.begin(),v.end());
    return v;
}

// Separable cubic cosine V = V0 (cos bx + cos by + cos bz): Vtilde(dm) = V0/2 for a unit +/- e_i step.
std::function<dcmplx(const ivec3_t&)> CosineVtilde(double V0)
{
    return [V0](const ivec3_t& dm)->dcmplx
    { return (dm.x*dm.x+dm.y*dm.y+dm.z*dm.z==1) ? dcmplx(0.5*V0) : dcmplx(0.0); };
}
} // namespace

class BandStructureTests : public ::testing::Test {};

// Plane waves: the empty-lattice band structure along Gamma-X-M-Gamma is exactly the folded free-electron
// ladder 1/2|k+G|^2 at every k on the path.
TEST_F(BandStructureTests, FreeElectronBandsPlaneWave)
{
    double a=8.0, Ecut=2.0; ivec3_t N(8,8,8);          // N also sets the k-path resolution (k=kIndex/N)
    UnitCell cell(a); Lattice_3D lat(cell,N);
    auto recip=lat.Reciprocal();
    std::vector<ivec3_t> path=KPath({ivec3_t(0,0,0), ivec3_t(4,0,0), ivec3_t(4,4,0), ivec3_t(0,0,0)}, 4);
    ASSERT_GT(path.size(),8u);

    for (const ivec3_t& k : path)
    {
        PlaneWave_IBS pw(recip,N,k,Ecut);
        std::vector<double> bands=ascending(SolveBands(pw,&cell));   // empty cell -> Nuclear = 0
        std::vector<double> fe=FreeEnergies(a,N,k,2);
        ASSERT_GE(bands.size(),4u);
        for (int n=0;n<4;n++) EXPECT_NEAR(bands[n],fe[n],1e-9);
    }
}

// The very same SolveBands drives LAPW (overlap != I): empty-lattice (Znuc=0) bands reproduce the
// free-electron ladder across the path, up to the LAPW linearization error.  That error is set by the
// distance of each band from the single linearization energy E_l, so a band structure spanning E=0..E_l
// shows a ~1e-2 spread with one E_l (real LAPW uses per-l E_l / local orbitals to tighten this).
TEST_F(BandStructureTests, FreeElectronBandsLAPW)
{
    double a=8.0, Ecut=2.0, Rmt=2.0, Elin=0.5; size_t lmax=6; ivec3_t N(8,8,8);
    UnitCell cell(a); Lattice_3D lat(cell,N);
    auto recip=lat.Reciprocal();
    std::vector<ivec3_t> path=KPath({ivec3_t(0,0,0), ivec3_t(4,0,0), ivec3_t(4,4,0)}, 3);

    for (const ivec3_t& k : path)
    {
        LAPW_IBS b(recip,N,k,Ecut,Rmt,lmax,Elin,0.0);
        std::vector<double> bands=ascending(SolveBands(b,&cell));
        std::vector<double> fe=FreeEnergies(a,N,k,2);
        ASSERT_GE(bands.size(),3u);
        for (int n=0;n<3;n++) EXPECT_NEAR(bands[n],fe[n],1e-2);      // single-E_l linearization spread
    }
}

TEST_F(BandStructureTests, DISABLED_CosineGapCalibration)
{
    double a=8.0, Ecut=2.0; ivec3_t N(8,8,8), X(4,0,0);
    UnitCell cell(a); Lattice_3D lat(cell,N);
    auto recip=lat.Reciprocal();
    for (double V0 : {0.0,0.05,0.1,0.2})
    {
        PlaneWave_IBS pw(recip,N,X,Ecut);
        std::vector<double> b=ascending(SolveBands(pw, pw.MakePotential(CosineVtilde(V0))));
        printf("V0=%.2f  gap(X)=b1-b0=%.5f  (NFE: V0)\n",V0,b[1]-b[0]);
    }
}

// Cosine potential -> a band gap.  Along Gamma->X the lowest two bands are degenerate at X for V=0;
// turning on V0 cos(b.r) opens a gap there.  Nearly-free-electron theory: the zone-boundary gap is
// 2|Vtilde(b)| = 2*(V0/2) = V0 (to first order; higher order ~ V0^2).
TEST_F(BandStructureTests, CosinePotentialOpensZoneBoundaryGap)
{
    double a=8.0, Ecut=2.0, V0=0.1; ivec3_t N(8,8,8), X(4,0,0);
    UnitCell cell(a); Lattice_3D lat(cell,N);
    auto recip=lat.Reciprocal();

    // Empty lattice: the two lowest bands are degenerate at the zone boundary X.
    PlaneWave_IBS pw0(recip,N,X,Ecut);
    std::vector<double> empty=ascending(SolveBands(pw0, pw0.MakePotential(CosineVtilde(0.0))));
    EXPECT_NEAR(empty[1]-empty[0], 0.0, 1e-9);

    // Cosine on: a gap of ~V0 opens between them.
    PlaneWave_IBS pw(recip,N,X,Ecut);
    std::vector<double> bands=ascending(SolveBands(pw, pw.MakePotential(CosineVtilde(V0))));
    double gap=bands[1]-bands[0];
    EXPECT_GT(gap, 0.5*V0);                 // a real gap opened
    EXPECT_NEAR(gap, V0, 0.1*V0);           // nearly-free-electron value 2|Vtilde(b)| = V0 (10% for O(V0^2))
}
