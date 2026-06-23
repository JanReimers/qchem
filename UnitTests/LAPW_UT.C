// file: LAPW_UT.C  Empty-lattice validation of the Linearized APW (LAPW) basis (lineage B).
//
// LAPW is a true band solver: a single generalized eigenproblem H c = e O c gives the band.  For the
// empty lattice the eigenvalues reproduce the free-electron ladder 1/2|k+G|^2 up to the LAPW
// linearization error (smallest for states near the linearization energy E_l).

#include <vector>
#include <algorithm>
#include <cmath>
#include <complex>
#include "gtest/gtest.h"

import qchem.BasisSet.Lattice_3D.LAPW_IBS;
import qchem.Lattice_3D;
import qchem.LASolver;
import qchem.Types;
import qchem.Blaze;
import qchem.Math;           // Pi

using BasisSet::Lattice_3D::LAPW_IBS;

namespace
{
std::vector<double> Bands(const LAPW_IBS& b, const Structure* cl)
{
    // Honour Orbital_1E_IBS: assemble the Hamiltonian as H = 1/2 Kinetic + Nuclear, exactly as the
    // Hamiltonian layer does for atoms/molecules/plane waves.
    chmat_t O=b.MakeOverlap();
    chmat_t H=0.5*b.MakeKinetic() + b.MakeNuclear(cl);
    LASolver<dcmplx>* las=LASolver<dcmplx>::Factory(qchem::Eigen);
    las->SetBasisOverlap(O);
    auto [U,e]=las->Solve(H);
    delete las;
    std::vector<double> ev;
    for (size_t i=0;i<e.size();i++) ev.push_back(e[i]);
    std::sort(ev.begin(),ev.end());
    return ev;
}

std::vector<double> FreeEnergies(double a, rvec3_t kf, int mr)
{
    double b=2*Pi/a;
    std::vector<double> e;
    for (int mx=-mr;mx<=mr;mx++) for (int my=-mr;my<=mr;my++) for (int mz=-mr;mz<=mr;mz++)
    {
        rvec3_t q((kf.x+mx)*b,(kf.y+my)*b,(kf.z+mz)*b);
        e.push_back(0.5*(q.x*q.x+q.y*q.y+q.z*q.z));
    }
    std::sort(e.begin(),e.end());
    return e;
}
} // namespace

class LAPWTests : public ::testing::Test {};

TEST_F(LAPWTests, DISABLED_Calibration)
{
    double a=8.0, R=2.0, Ecut=2.0; size_t lmax=8;
    ivec3_t N(4,4,4), k(1,0,0);  rvec3_t kf(0.25,0,0);
    UnitCell cell(a); Lattice_3D lat(cell,N);
    std::vector<double> fe=FreeEnergies(a,kf,3);
    for (double El : {0.1,0.2,0.4})
    {
        LAPW_IBS b(lat.Reciprocal(),N,k,Ecut,R,lmax,El);
        std::vector<double> bands=Bands(b,&cell);
        printf("--- Elin=%.2f  nPW=%zu ---\n",El,b.GetNumFunctions());
        for (int i=0;i<6 && i<(int)bands.size();i++)
            printf("  band[%d]=% .6f  free=% .6f  d=% .2e\n",i,bands[i],fe[i],bands[i]-fe[i]);
    }
}

// LAPW is a band solver: a single generalized eigenproblem reproduces the free-electron ladder.
TEST_F(LAPWTests, EmptyLatticeBandsMatchFreeElectron)
{
    double a=8.0, R=2.0, Ecut=2.0; size_t lmax=8; double El=0.2;
    ivec3_t N(4,4,4), k(1,0,0); rvec3_t kf(0.25,0,0);
    UnitCell cell(a); Lattice_3D lat(cell,N);
    LAPW_IBS b(lat.Reciprocal(),N,k,Ecut,R,lmax,El);

    std::vector<double> bands=Bands(b,&cell), fe=FreeEnergies(a,kf,3);
    ASSERT_GE(bands.size(),5u);
    for (int i=0;i<5;i++) EXPECT_NEAR(bands[i],fe[i],1e-4);  // linearised, so not exact -- but close
}

// The LAPW linearization is most accurate for states near E_l: putting E_l on the 2nd level makes that
// band's error beat the (further) lowest band's error.
TEST_F(LAPWTests, LinearizationMostAccurateNearElin)
{
    double a=8.0, R=2.0, Ecut=2.0; size_t lmax=8;
    ivec3_t N(4,4,4), k(1,0,0); rvec3_t kf(0.25,0,0);
    UnitCell cell(a); Lattice_3D lat(cell,N);
    std::vector<double> fe=FreeEnergies(a,kf,3);
    LAPW_IBS b(lat.Reciprocal(),N,k,Ecut,R,lmax,fe[1]); // E_l on the 2nd level
    std::vector<double> bands=Bands(b,&cell);

    double errNear=std::abs(bands[1]-fe[1]);            // at E_l
    double errFar =std::abs(bands[0]-fe[0]);            // below E_l
    EXPECT_LT(errNear,errFar);
}
