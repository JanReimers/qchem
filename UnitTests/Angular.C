// File: Angular.C  Test some identities for angular ERI integrals.

#include "gtest/gtest.h"
#include <iostream>
#include <blaze/Math.h>
import qchem.BasisSet.Atom.Internal.AngularIntegrals;

using std::cout;
using std::endl;

class AngularTests : public ::testing::Test
{
    public:
    using rvec11_t=AngularIntegrals::rvec11_t;
    AngularTests() {}

    rvec11_t msum_direct(size_t la, size_t lc, const std::set<int> m1,const std::set<int> m2)
    {
        rvec11_t Ak(0.0);
        for (auto ma:m1)
        for (auto mc:m2)
        {
            assert(abs(ma)<=la);
            assert(abs(mc)<=lc);
            Ak+=AngularIntegrals::Coulomb(la,lc,ma,mc);
        }
        return Ak;
    }
    rvec11_t msum_exchange(size_t la, size_t lc, const std::set<int> m1,const std::set<int> m2)
    {
        rvec11_t Ak(0.0);
        for (auto ma:m1)
        for (auto mc:m2)
        {
            assert(abs(ma)<=la);
            assert(abs(mc)<=lc);
            Ak+=AngularIntegrals::Exchange(la,lc,ma,mc);
        }
        return Ak;
    }

};

TEST_F(AngularTests,FullSums)
{
    size_t LMax=4;
    for (size_t la=0;la<=LMax;la++)
    for (size_t lc=0;lc<=LMax;lc++)
    {
        rvec11_t dac=AngularIntegrals::Coulomb(la,lc);
        rvec11_t eac=AngularIntegrals::Exchange(la,lc);

        rvec11_t dac_sum(0.0),eac_sum(0.0);
        for (int ma=-(int)la;ma<=(int)la;ma++)
            for (int mc=-(int)lc;mc<=(int)lc;mc++)
            {
                dac_sum+=AngularIntegrals::Coulomb(la,lc,ma,mc);
                eac_sum+=AngularIntegrals::Exchange(la,lc,ma,mc);

            }

        dac_sum/=(2*la+1)*(2*lc+1);
        eac_sum/=(2*la+1)*(2*lc+1);
        EXPECT_NEAR(max(abs(dac-dac_sum)),0.0,1e-12);
        EXPECT_NEAR(max(abs(eac-eac_sum)),0.0,3e-14);
        // cout << la << " " << lc << "  |  " << Max(fabs(dac-dac_sum))  << " " << Max(fabs(eac-eac_sum)) << endl;
    
    }
}

TEST_F(AngularTests,p_orbitals)
{
    size_t la=1,lc=1;
    std::set<int> m11({-1,0}),m12({1});
    std::set<int> m21({-1,0,1});

    rvec11_t d11_21=msum_direct  (la,lc,m11,m21);
    rvec11_t d12_21=msum_direct  (la,lc,m12,m21);
    rvec11_t e11_21=msum_exchange(la,lc,m11,m21);
    rvec11_t e12_21=msum_exchange(la,lc,m12,m21);

    double de=947.48202250457894;
    EXPECT_NEAR(d11_21[0],de,0.0);
    EXPECT_EQ  (d11_21[1],0.0);
    EXPECT_NEAR(d12_21[0],de/2,6e-14);
    EXPECT_EQ  (d12_21[1],0.0);
    EXPECT_NEAR(e11_21[0],de/3,6e-14);
    EXPECT_EQ  (e11_21[1],0.0);
    EXPECT_NEAR(e11_21[2],de/7.5,1e-13);
    EXPECT_EQ  (e11_21[3],0.0);
    EXPECT_NEAR(e12_21[0],de/6,3e-14);
    EXPECT_EQ  (e12_21[1],0.0);
    EXPECT_NEAR(e12_21[2],de/15,4e-14);
    EXPECT_EQ  (e12_21[3],0.0);

}
