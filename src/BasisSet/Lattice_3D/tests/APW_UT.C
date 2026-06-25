// file: APW_UT.C  Empty-lattice validation of the Augmented Plane Wave (APW) basis (lineage B).
//
// The APW secular matrix Gamma(E) = H(E) - E O(E) must be SINGULAR exactly at the free-electron
// energies E = 1/2|k+G|^2 (the empty-lattice limit), and non-singular between them.  This exercises
// the augmentation + value matching + interstitial/sphere assembly against the exact answer.

#include <vector>
#include <algorithm>
#include <cmath>
#include <complex>
#include "gtest/gtest.h"

import qchem.BasisSet.Lattice_3D.APW_IBS;
import qchem.Lattice_3D;     // UnitCell, Lattice_3D, ReciprocalLattice
import qchem.Types;
import qchem.Blaze;
import qchem.Math;           // Pi

using BasisSet::Lattice_3D::APW_IBS;

namespace
{
double MinAbsEig(const chmat_t& G)
{
    rvec_t d; mat_t<dcmplx> U;
    blazem::eigen(G,d,U);
    double m=1e300;
    for (size_t i=0;i<d.size();i++) m=std::min(m,std::abs(d[i]));
    return m;
}

// Free-electron energies 1/2|k+G|^2 for a cubic cell (B = (2pi/a) I), sorted ascending.
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

class APWTests : public ::testing::Test {};

TEST_F(APWTests, DISABLED_Calibration)
{
    double a=8.0, R=2.0, Ecut=2.0;
    ivec3_t N(4,4,4), k(1,0,0);
    rvec3_t kf(0.25,0,0);
    UnitCell cell(a);
    Lattice_3D lat(cell,N);
    std::vector<double> fe=FreeEnergies(a,kf,2);

    for (size_t lmax : {4u,6u,8u,12u})
    {
        APW_IBS apw(lat.Reciprocal(),N,k,Ecut,R,lmax);
        printf("--- lmax=%2zu  nPW=%zu ---\n",lmax,apw.GetNumFunctions());
        for (int n=0;n<4;n++)
        {
            double Es=fe[n];
            double Eoff=0.5*(fe[n]+fe[n+1]);
            printf("  E*=%.5f  minEig(E*)=%.3e   Eoff=%.5f minEig(Eoff)=%.3e\n",
                   Es,MinAbsEig(apw.MakeSecular(Es)),Eoff,MinAbsEig(apw.MakeSecular(Eoff)));
        }
    }
}

// The APW secular matrix is singular exactly at the free-electron energies and non-singular between
// them -- the empty-lattice limit, recovered to machine precision at lmax=8.
TEST_F(APWTests, EmptyLatticeSecularSingularAtFreeElectronEnergies)
{
    double a=8.0, R=2.0, Ecut=2.0; size_t lmax=8;
    ivec3_t N(4,4,4), k(1,0,0); rvec3_t kf(0.25,0,0);
    UnitCell cell(a); Lattice_3D lat(cell,N);
    APW_IBS apw(lat.Reciprocal(),N,k,Ecut,R,lmax);
    std::vector<double> fe=FreeEnergies(a,kf,2);

    EXPECT_LT(MinAbsEig(apw.MakeSecular(fe[0])),1e-9);  // singular at 1/2|k+G|^2 ...
    EXPECT_LT(MinAbsEig(apw.MakeSecular(fe[1])),1e-9);
    EXPECT_GT(MinAbsEig(apw.MakeSecular(0.5*(fe[0]+fe[1]))),1e-2); // ... non-singular between levels
}

// Raising l_max drives the residual singularity at a free-electron energy toward zero (convergence of
// the augmentation).
TEST_F(APWTests, SingularityConvergesWithLmax)
{
    double a=8.0, R=2.0, Ecut=2.0;
    ivec3_t N(4,4,4), k(1,0,0); rvec3_t kf(0.25,0,0);
    UnitCell cell(a); Lattice_3D lat(cell,N);
    double Es=FreeEnergies(a,kf,2)[1];                 // 2nd level (1st is ~machine-zero already)
    APW_IBS lo(lat.Reciprocal(),N,k,Ecut,R,4);
    APW_IBS hi(lat.Reciprocal(),N,k,Ecut,R,8);
    EXPECT_LT(MinAbsEig(hi.MakeSecular(Es)),MinAbsEig(lo.MakeSecular(Es)));
}
