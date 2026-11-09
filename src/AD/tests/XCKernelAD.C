// File: AD/tests/XCKernelAD.C  The XC kernel's AD gradients against its own finite differences, on a small synthetic system (two atoms,
// two s shells and a p shell, a few mesh points): dE/dR and dE/dD.  The facade-level checks against the library's term are in
// src/Calculation/tests/ForceFacade.C.
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
import qchem.AD.XCKernel;
import qchem.Types;

using namespace qchem;

namespace
{
AD::XCKernelData Small()
{
    AD::XCKernelData d;
    d.natom=2; d.beckeOrder=3; d.nBlocks=3;
    const double u[][3]={{0.3,0.1,-0.2},{-0.7,0.9,0.4},{1.1,-0.5,0.8},{0.2,0.2,0.2},{-0.4,-0.6,0.1}};
    for (auto& p : u) { d.u.push_back(p[0]); d.u.push_back(p[1]); d.u.push_back(p[2]); d.w.push_back(0.37); }
    d.shells.push_back({0, 0, {{0,0,0}}, {1.3}, {0.8, 2.5}, {0.6, 0.4}});
    d.shells.push_back({1, 1, {{0,0,0}}, {0.9}, {1.1}, {1.0}});
    d.shells.push_back({0, 2, {{1,0,0},{0,1,0},{0,0,1}}, {1.7,1.7,1.7}, {0.9}, {1.0}});
    d.nbf=5;
    const double Dm[5][5]={{0.9,0.2,0.1,0.0,0.05},{0.2,0.7,0.0,0.1,0.0},{0.1,0.0,0.4,0.0,0.0},{0.0,0.1,0.0,0.3,0.02},{0.05,0.0,0.0,0.02,0.2}};
    for (auto& r : Dm) for (double v : r) d.D.push_back(v);
    return d;
}
const std::vector<double> R0={0.0,0.0,0.0, 0.4,0.3,1.6};
}

TEST(XCKernelAD, ItsGradientsAreItsFiniteDifferences)
{
    const AD::XCKernelData d=Small();
    const auto gR=AD::XCEnergyGradientR(d, R0);
    const auto gD=AD::XCEnergyGradientD(d, R0);
    const double h=1e-5;
    for (size_t i=0;i<R0.size();i++)
    {
        auto p=R0, m=R0; p[i]+=h; m[i]-=h;
        EXPECT_NEAR(gR[i], (AD::XCEnergy(d,p)-AD::XCEnergy(d,m))/(2*h), 1e-8) << "dE/dR " << i;
    }
    for (size_t i=0;i<d.nbf;i++) for (size_t j=0;j<d.nbf;j++)
    {
        auto p=d, m=d; p.D[i*d.nbf+j]+=h; m.D[i*d.nbf+j]-=h;
        EXPECT_NEAR(gD(i,j), (AD::XCEnergy(p,R0)-AD::XCEnergy(m,R0))/(2*h), 1e-8) << "dE/dD " << i << j;
    }
}
