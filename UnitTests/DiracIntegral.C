// File: DiracIntegral.C  Test the Dirac basis sets and integral engine.


#include "gtest/gtest.h"
#include "nlohmann/json.hpp"
#include <iostream>
#include <valarray>
#include <cmath>
#include <iomanip>
#include <blaze/Math.h>
import qchem.LAParams;

import qchem.BasisSet.Internal.IrrepBasisSet;

import qchem.Factory;
import qchem.BasisSet;
import qchem.IrrepBasisSet;
import qchem.BasisSet.Internal.ERI4;
import Common.Constants;
import qchem.Cluster;
import qchem.Mesh.Integrator;
import qchem.Atom;
import qchem.Molecule;
import qchem.Streamable;

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
    , lap({qchem::SVD,1e-6})
    , sbs(0)
    , gbs(0)
    , cl(new Molecule())
    {
        nlohmann::json js = {{"type",BasisSetAtom::Type::Slater_RKB}, {"N", 3}, {"emin", 0.1}, {"emax", 10.0} };
        sbs=BasisSetAtom::Factory(js,2);
        sbs->Set(lap);
        js = {{"type",BasisSetAtom::Type::Gaussian_RKB}, {"N", 3}, {"emin", 0.1}, {"emax", 10.0} };
        gbs=BasisSetAtom::Factory(js,2);
        gbs->Set(lap);
        cl->Insert(new Atom(Z,0.0,Vector3D(0,0,0)));
        MeshParams mp({qchem::MHL,200,3,2.0,qchem::Gauss,1,0,0,3});
        mintegrator=new MeshIntegrator<double>(cl->CreateMesh(mp));
    }
    
    static const Real_IBS* GetLarge(const Real_IBS* ibs)
    {
        assert(ibs);
        const Orbital_RKB_IBS_Common<double>* dirbs=dynamic_cast<const Orbital_RKB_IBS_Common<double>*>(ibs);
        assert(dirbs);
        return dirbs->itsRKBL;
    }
    static const Real_IBS* GetSmall(const Real_IBS* ibs)
    {
        assert(ibs);
        const Orbital_RKB_IBS_Common<double>* dirbs=dynamic_cast<const Orbital_RKB_IBS_Common<double>*>(ibs);
        assert(dirbs);
        return dirbs->itsRKBS;
    }

    
    
    static rsmat_t merge_diag(const rsmat_t& l,const rsmat_t& s)
    {
        return Orbital_RKB_IBS_Common<double>::merge_diag(l,s);
    }

    static rsmat_t merge_off_diag(const rmat_t& l)
    {
        return Orbital_RKB_IBS_Common<double>::merge_off_diag(l);
    }

    int Lmax, Z;
    LAParams lap;
    BasisSet* sbs;
    BasisSet* gbs;
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

