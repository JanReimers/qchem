// File: DiracIntegral.C  Test the Dirac basis sets and integral engine.


#include "gtest/gtest.h"
#include "nlohmann/json.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <blaze/Math.h>
import qchem.LAParams;


import qchem.Factory;
import qchem.BasisSet.Internal.ERI4;
import Common.Constants;
import qchem.Cluster;
import qchem.Mesh.Integrator;
import qchem.Streamable;
import qchem.BasisSet.Internal.Orbital_DHF_IBS;
using BasisSet::Real_OIBS;
using Real_IBS=Real_OIBS;
using RKB_OIBS=BasisSet::Orbital_RKB_IBS_Imp<double>;
using RKBL_OIBS=BasisSet::Orbital_RKBL_IBS<double>;
using RKBS_OIBS=BasisSet::Orbital_RKBS_IBS<double>;

using std::cout;
using std::endl;


//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class DiracIntegralTests : public ::testing::Test
{
public:
    DiracIntegralTests()
    : Lmax(0   )
    , Z(1)
    , sbs(0)
    , gbs(0)
    , cl(new Atom(Z,0.0,Vector3D(0,0,0)))
    {
        nlohmann::json js = {{"type",abs_t::Slater_RKB}, {"N", 3}, {"emin", 0.1}, {"emax", 10.0} };
        sbs=BasisSet::Atom::Factory(js,2);
        js = {{"type",abs_t::Gaussian_RKB}, {"N", 3}, {"emin", 0.1}, {"emax", 10.0} };
        gbs=BasisSet::Atom::Factory(js,2);
        MeshParams mp({qchem::MHL,200,3,2.0,qchem::Gauss,1,0,0,3});
        mintegrator=new MeshIntegrator<double>(cl->CreateMesh(mp));
    }
    
    static const RKBL_OIBS* GetLarge(const Real_IBS* ibs)
    {
        assert(ibs);
        const RKB_OIBS* dirbs=dynamic_cast<const RKB_OIBS*>(ibs);
        assert(dirbs);
        return dirbs->itsRKBL;
    }
    static const RKBS_OIBS* GetSmall(const Real_IBS* ibs)
    {
        assert(ibs);
        const RKB_OIBS* dirbs=dynamic_cast<const RKB_OIBS*>(ibs);
        assert(dirbs);
        return dirbs->itsRKBS;
    }

    
    
    static rsmat_t merge_diag(const rsmat_t& l,const rsmat_t& s)
    {
        return RKB_OIBS::merge_diag(l,s);
    }

    static rsmat_t merge_off_diag(const rmat_t& l)
    {
        return RKB_OIBS::merge_off_diag(l);
    }

    int Lmax, Z;
    Real_BS* sbs;
    Real_BS* gbs;
    Cluster* cl;
    MeshIntegrator<double>* mintegrator;
};

TEST_F(DiracIntegralTests, SlaterBasisSet)
{
    
    cout << *sbs << endl;
}

TEST_F(DiracIntegralTests, GaussianBasisSet)
{
    
    cout << *gbs << endl;
}


TEST_F(DiracIntegralTests, SlaterOverlap)
{
    
    for (auto oi:sbs->Iterate<Real_OIBS >())
    {
        rsmat_t S=oi->Overlap();
        {
            rvec_t d=blaze::diagonal(S);
            for (size_t i=0;i<d.size()/2;i++) EXPECT_NEAR(d[i],1.0,1e-15);
        }
        // cout << std::fixed << std::setprecision(3) << std::setw(6) << S << S1 << endl;
        rsmat_t SLnum = mintegrator->Overlap(*GetLarge(oi));
        rsmat_t SSnum = mintegrator->Overlap(*GetSmall(oi));
        rsmat_t Snum=merge_diag(SLnum,SSnum);
        // cout << Snum << S << endl;

        EXPECT_NEAR(max(abs(S-Snum)),0.0,1e-14);
    }
}

TEST_F(DiracIntegralTests, GaussianOverlap)
{
    
    for (auto oi:gbs->Iterate<Real_OIBS >())
    {
        rsmat_t S=oi->Overlap();
        {
            rvec_t d=blaze::diagonal(S);
            for (size_t i=0;i<d.size()/2;i++) EXPECT_NEAR(d[i],1.0,1e-15);
        }
        // cout << std::fixed << std::setprecision(3) << std::setw(6) << S << S1 << endl;
        rsmat_t SLnum = mintegrator->Overlap(*GetLarge(oi));
        rsmat_t SSnum = mintegrator->Overlap(*GetSmall(oi));
//        cout << SLnum << SSnum << endl;
        rsmat_t Snum=merge_diag(SLnum,SSnum);
        // cout << Max(fabs(S-Snum)) << endl;
        EXPECT_NEAR(max(abs(S-Snum)),0.0,1e-14);
    }
}


TEST_F(DiracIntegralTests, SlaterNuclear)
{
    
    int Z=cl->GetNuclearCharge();
    for (auto oi:sbs->Iterate<Real_OIBS >())
    {
        rsmat_t Ven=oi->Nuclear(cl);
        rsmat_t VenLnum = -Z*mintegrator->Inv_r1(*GetLarge(oi));
        rsmat_t VenSnum = -Z*mintegrator->Inv_r1(*GetSmall(oi));
        rsmat_t Vennum=merge_diag(VenLnum,VenSnum);
        //cout << "Ven=" << Ven << endl << "Ven num=" << Vennum << endl;
        //Because of the singularity at the origin, the error is larger than the other integrals.
        EXPECT_NEAR(max(abs(Ven-Vennum)),0.0,1e-11);        
    }
}


TEST_F(DiracIntegralTests, GaussianNuclear)
{
    
    int Z=cl->GetNuclearCharge();
    for (auto oi:gbs->Iterate<Real_OIBS >())
    {
        rsmat_t Ven=oi->Nuclear(cl);
        rsmat_t VenLnum = -Z*mintegrator->Inv_r1(*GetLarge(oi));
        rsmat_t VenSnum = -Z*mintegrator->Inv_r1(*GetSmall(oi));
        rsmat_t Vennum=merge_diag(VenLnum,VenSnum);
        //cout << "Ven=" << Ven << endl << "Ven num=" << Vennum << endl;
        // cout << "Ven=" << Ven << endl << "Ven1=" << Ven1 << endl;
        //Because of the singularity at the origin, the error is larger than the other integrals.
        EXPECT_NEAR(max(abs(Ven-Vennum)),0.0,1e-11);        
    }
}


TEST_F(DiracIntegralTests, SlaterKinetic)
{
    
    for (auto oi:sbs->Iterate<Real_OIBS >())
    {
        rsmat_t K=oi->Kinetic();
        rvec_t d=blaze::diagonal(K);    
        for (size_t i=0;i<d.size();i++) EXPECT_NEAR(d[i],0.0,1e-15);       
        //cout << std::fixed << std::setprecision(3) << std::setw(6) << K << endl;
         rmat_t KnumL = mintegrator->Grada_b(*GetLarge(oi),*GetSmall(oi));
        rsmat_t Knum=merge_off_diag(KnumL);
        // cout << "K=" << K << endl << "Knum=" << Knum << endl;
        EXPECT_NEAR(max(abs(K-Knum)),0.0,1e-11);      
    }
}
TEST_F(DiracIntegralTests, GaussianKinetic)
{
    
    for (auto oi:gbs->Iterate<Real_OIBS >())
    {
        rsmat_t K=oi->Kinetic();
        rvec_t d=blaze::diagonal(K);    
        for (size_t i=0;i<d.size();i++) EXPECT_NEAR(d[i],0.0,1e-15);       
        //cout << std::fixed << std::setprecision(3) << std::setw(6) << K << endl;
         rmat_t KnumL = mintegrator->Grada_b(*GetLarge(oi),*GetSmall(oi));
        rsmat_t Knum=merge_off_diag(KnumL);
        // cout << "K=" << K << endl << "Knum=" << Knum << endl;
        EXPECT_NEAR(max(abs(K-Knum)),0.0,1e-11);         
    }
}

 //  For symmetric matricies, Sum() should know to double all the off diagonal elements.
TEST_F(DiracIntegralTests, Contraction)
{
    size_t N=5;
    rsmat_t H(N),D(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
        {
            H(i,j)=std::rand();
            D(i,j)=std::rand();
        }
    rmat_t Hf(H),Df(D);
    double c =sum(H %D );
    double cf=sum(Hf%Df);
    EXPECT_NEAR(c,cf,5e-15);
}

