// File: BeckePartitionTests.C  The molecular Becke partition weight and its geometry gradient (qchem.Mesh.BeckePartition): the gradient
// against a central finite difference of the SAME weight, with the point riding with its owner, and translation invariance.
#include "gtest/gtest.h"
#include <cmath>
#include <vector>
import qchem.Mesh.BeckePartition;

using namespace qchem;

namespace
{
double& C(rvec3_t& v, size_t k) {return k==0 ? v.x : k==1 ? v.y : v.z;}
double  C(const rvec3_t& v, size_t k) {return k==0 ? v.x : k==1 ? v.y : v.z;}
// A bent triatomic plus a fourth atom off its plane: no symmetry to hide a sign error behind
const std::vector<rvec3_t> kAtoms={rvec3_t(0.0,0.0,0.0), rvec3_t(1.43,0.0,1.11), rvec3_t(-1.43,0.0,1.11), rvec3_t(0.4,2.1,-0.7)};
const std::vector<rvec3_t> kOffsets={rvec3_t(0.3,-0.2,0.1), rvec3_t(-0.9,0.6,0.8), rvec3_t(1.4,1.1,-0.5), rvec3_t(0.05,0.02,-0.04), rvec3_t(-2.0,0.3,1.7)};
}

TEST(BeckePartition, TheCutoffDerivativeIsTheIteratedChainRule)
{
    for (int k=0;k<=3;++k)
        for (double mu : {-0.95,-0.6,-0.2,0.0,0.15,0.5,0.9})
        {
            const double h=1e-6, fd=(qcMesh::BeckeCutoff(mu+h,k)-qcMesh::BeckeCutoff(mu-h,k))/(2*h);
            EXPECT_NEAR(qcMesh::BeckeCutoffDerivative(mu,k), fd, 1e-8) << "k " << k << " mu " << mu;
        }
}

TEST(BeckePartition, TheWeightsAreAPartitionOfUnity)
{
    for (const auto& off : kOffsets)
        for (int k : {2,3})
        {
            const rvec3_t r=kAtoms[1]+off;
            double sum=0.0;
            for (size_t A=0;A<kAtoms.size();++A) sum+=qcMesh::BeckePartitionWeight(kAtoms, A, r, k);
            EXPECT_NEAR(sum, 1.0, 1e-14);
        }
    EXPECT_DOUBLE_EQ(qcMesh::BeckePartitionWeight({rvec3_t(0,0,0)}, 0, rvec3_t(1,2,3), 2), 1.0);
}

TEST(BeckePartition, TheGradientIsTheFiniteDifferenceOfTheWeightWithThePointRidingWithItsOwner)
{
    const double h=1e-5;
    double largest=0.0;   // the points must sit where the partition actually moves, or the test proves nothing
    for (int k : {2,3})
        for (size_t A=0;A<kAtoms.size();++A)
            for (const auto& off : kOffsets)
            {
                const auto g=qcMesh::BeckePartitionGradient(kAtoms, A, kAtoms[A]+off, k);
                ASSERT_EQ(g.size(), kAtoms.size());
                rvec3_t sum(0,0,0);
                for (size_t D=0;D<kAtoms.size();++D)
                {
                    sum+=g[D]; largest=std::max(largest, norm(g[D]));
                    for (size_t x=0;x<3;++x)
                    {
                        auto Rp=kAtoms, Rm=kAtoms; C(Rp[D],x)+=h; C(Rm[D],x)-=h;
                        const double fd=(qcMesh::BeckePartitionWeight(Rp, A, Rp[A]+off, k)-qcMesh::BeckePartitionWeight(Rm, A, Rm[A]+off, k))/(2*h);
                        EXPECT_NEAR(C(g[D],x), fd, 1e-8) << "k " << k << " owner " << A << " atom " << D << " axis " << x;
                    }
                }
                EXPECT_LT(norm(sum), 1e-13) << "translation invariance: the gradient sums to zero over the atoms";
            }
    EXPECT_GT(largest, 0.1) << "the sampled points are in the partition's transition regions";
    const auto one=qcMesh::BeckePartitionGradient({rvec3_t(0,0,0)}, 0, rvec3_t(1,0,0), 2);
    ASSERT_EQ(one.size(), 1u);
    EXPECT_EQ(norm(one[0]), 0.0) << "one atom: the weight is 1 everywhere";
}
