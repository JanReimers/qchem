// File: AD/tests/BeckeAD.C  The AD kernel's Becke partition against qchem.Mesh.BeckePartition: its weight equals the library's, and its
// AD gradient equals the library's (which is AD of the library's own weight, in place: src/Mesh/tests checks THAT against the hand
// oracle).  Without Enzyme both gradients refuse with NoGradient.
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
import qchem.AD.XCKernel;
import qchem.Mesh.BeckePartition;
import qchem.Types;

using namespace qchem;

namespace
{
// A bent triatomic (water-like, bohr) and a spread of points riding with each owner.
std::vector<rvec3_t> Atoms() {return {rvec3_t(0,0,0.2217), rvec3_t(0,1.4309,-0.8867), rvec3_t(0,-1.4309,-0.8867)};}
std::vector<double> Flat(const std::vector<rvec3_t>& R) {std::vector<double> x; for (auto& r : R) {x.push_back(r.x); x.push_back(r.y); x.push_back(r.z);} return x;}
std::vector<rvec3_t> Offsets() {return {rvec3_t(0.3,0.1,-0.2), rvec3_t(-0.7,0.9,0.4), rvec3_t(1.1,-0.5,0.8), rvec3_t(0.05,0.02,0.01), rvec3_t(-1.6,0.3,-0.9)};}
}

TEST(BeckeAD, TheKernelsWeightIsTheLibrarysWeight)
{
    const auto R=Atoms();
    for (int k : {2,3})
        for (size_t A=0;A<R.size();A++)
            for (const rvec3_t& u : Offsets())
                EXPECT_NEAR(AD::XCKernelBeckeWeight(Flat(R), A, u, k), qcMesh::BeckePartitionWeight(R, A, R[A]+u, k), 1e-15) << "owner " << A;
}

TEST(BeckeAD, ItsADGradientIsTheHandDerivedOne)
{
    const auto R=Atoms();
    for (int k : {2,3})
        for (size_t A=0;A<R.size();A++)
            for (const rvec3_t& u : Offsets())
            {
                const auto g=AD::XCKernelBeckeWeightGradient(Flat(R), A, u, k);
                const auto hand=qcMesh::BeckePartitionGradient(R, A, R[A]+u, k);
                rvec3_t sum(0,0,0);
                for (size_t D=0;D<R.size();D++)
                {
                    EXPECT_NEAR(g[3*D],   hand[D].x, 1e-13) << "owner " << A << " atom " << D;
                    EXPECT_NEAR(g[3*D+1], hand[D].y, 1e-13) << "owner " << A << " atom " << D;
                    EXPECT_NEAR(g[3*D+2], hand[D].z, 1e-13) << "owner " << A << " atom " << D;
                    sum+=rvec3_t(g[3*D], g[3*D+1], g[3*D+2]);
                }
                EXPECT_LT(norm(sum), 1e-14) << "translation invariance: the weight's gradient sums to zero";
            }
}
