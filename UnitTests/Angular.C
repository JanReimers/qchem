// File: Angular.C  Test some identities for angular ERI integrals.

#include "gtest/gtest.h"
#include <iostream>
import qchem.BasisSet.Atom.Internal.AngularIntegrals;

using std::cout;
using std::endl;

class AngularTests : public ::testing::Test
{
    public:
    AngularTests() {StreamableObject::SetToPretty();}

    RVec msum_direct(size_t la, size_t lc, const std::set<int> m1,const std::set<int> m2)
    {
        RVec Ak(1);
        for (auto ma:m1)
        for (auto mc:m2)
        {
            assert(abs(ma)<=la);
            assert(abs(mc)<=lc);
            Ak+=AngularIntegrals::Coulomb(la,lc,ma,mc);
        }
        return Ak;
    }
    RVec msum_exchange(size_t la, size_t lc, const std::set<int> m1,const std::set<int> m2)
    {
        RVec Ak(AngularIntegrals::Exchange(la,lc).GetLimits());
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
        RVec dac=AngularIntegrals::Coulomb(la,lc);
        RVec eac=AngularIntegrals::Exchange(la,lc);
        RVec dac_sum(dac.GetLimits(),0.0),eac_sum(eac.GetLimits(),0.0);
        for (int ma=-(int)la;ma<=(int)la;ma++)
            for (int mc=-(int)lc;mc<=(int)lc;mc++)
            {
                dac_sum+=AngularIntegrals::Coulomb(la,lc,ma,mc);
                eac_sum+=AngularIntegrals::Exchange(la,lc,ma,mc);

            }

        dac_sum/=(2*la+1)*(2*lc+1);
        eac_sum/=(2*la+1)*(2*lc+1);
        EXPECT_NEAR(Max(fabs(dac-dac_sum)),0.0,1e-12);
        EXPECT_NEAR(Max(fabs(eac-eac_sum)),0.0,3e-14);
        // cout << la << " " << lc << "  |  " << Max(fabs(dac-dac_sum))  << " " << Max(fabs(eac-eac_sum)) << endl;
    
    }
}

TEST_F(AngularTests,p_orbitals)
{
    size_t la=1,lc=1;
    std::set<int> m11({-1,0}),m12({1});
    std::set<int> m21({-1,0,1});

    RVec d11_21=msum_direct  (la,lc,m11,m21);
    RVec d12_21=msum_direct  (la,lc,m12,m21);
    RVec e11_21=msum_exchange(la,lc,m11,m21);
    RVec e12_21=msum_exchange(la,lc,m12,m21);

    cout << d11_21 << " " << d12_21 << endl;
    cout << e11_21 << " " << e12_21 << endl;
}
