// File: DiracIntegral.C  Test the Dirac basis sets and integral engine.


#include "gtest/gtest.h"
#include "nlohmann/json.hpp"
#include <iostream>
#include <valarray>
#include <cmath>
#include <iomanip>

import qchem.LAParams;

import qchem.BasisSet.IBS_Common;

import qchem.Factory;
import qchem.BasisSet;
import qchem.Irrep_BS;
import qchem.BasisSet.ERI4;
import Common.Constants;
import qchem.Cluster;
import Mesh.Integrator;
import qchem.Atom;
import qchem.Molecule;

using std::cout;
using std::endl;


//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class DiracIntegralTests : public ::testing::Test
{
public:
    typedef SMatrix<double> SMat;
    typedef  Matrix<double>  Mat;
    
    DiracIntegralTests()
    : Lmax(0   )
    , Z(1)
    , lap({qchem::Lapack,qchem::SVD,1e-6,1e-12})
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
    
    static const TIrrepBasisSet<double>* GetLarge(const IrrepBasisSet* ibs)
    {
        assert(ibs);
        const Orbital_RKB_IBS_Common<double>* dirbs=dynamic_cast<const Orbital_RKB_IBS_Common<double>*>(ibs);
        assert(dirbs);
        return dirbs->itsRKBL;
    }
    static const TIrrepBasisSet<double>* GetSmall(const IrrepBasisSet* ibs)
    {
        assert(ibs);
        const Orbital_RKB_IBS_Common<double>* dirbs=dynamic_cast<const Orbital_RKB_IBS_Common<double>*>(ibs);
        assert(dirbs);
        return dirbs->itsRKBS;
    }

    
    
    static SMat merge_diag(const SMat& l,const SMat& s)
    {
        return Orbital_RKB_IBS_Common<double>::merge_diag(l,s);
    }

    static SMat merge_off_diag(const Mat& l)
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
    StreamableObject::SetToPretty();
    cout << *sbs << endl;
}

TEST_F(DiracIntegralTests, GaussianBasisSet)
{
    StreamableObject::SetToPretty();
    cout << *gbs << endl;
}


TEST_F(DiracIntegralTests, SlaterOverlap)
{
    StreamableObject::SetToPretty();
    for (auto oi:sbs->Iterate<TOrbital_IBS<double> >())
    {
        SMatrix<double> S=oi->Overlap();
        for (auto d:Vector<double>(S.GetDiagonal())) EXPECT_NEAR(d,1.0,1e-15);
        // cout << std::fixed << std::setprecision(3) << std::setw(6) << S << S1 << endl;
        const TIrrepBasisSet<double>* l=GetLarge(oi);
        const TIrrepBasisSet<double>* s=GetSmall(oi);
        SMatrix<double> SLnum = mintegrator->Overlap(*l);
        SMatrix<double> SSnum = mintegrator->Overlap(*s);
        SMat Snum=merge_diag(SLnum,SSnum);
        // cout << Snum << S << endl;

        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,1e-14);
    }
}

TEST_F(DiracIntegralTests, GaussianOverlap)
{
    StreamableObject::SetToPretty();
    for (auto oi:gbs->Iterate<TOrbital_IBS<double> >())
    {
        SMatrix<double> S=oi->Overlap();
        for (auto d:Vector<double>(S.GetDiagonal())) EXPECT_NEAR(d,1.0,1e-15);
        // cout << std::fixed << std::setprecision(3) << std::setw(6) << S << S1 << endl;
        const TIrrepBasisSet<double>* l=GetLarge(oi);
        const TIrrepBasisSet<double>* s=GetSmall(oi);
        SMatrix<double> SLnum = mintegrator->Overlap(*l);
        SMatrix<double> SSnum = mintegrator->Overlap(*s);
//        cout << SLnum << SSnum << endl;
        SMat Snum=merge_diag(SLnum,SSnum);
        // cout << Max(fabs(S-Snum)) << endl;
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,1e-14);
    }
}


TEST_F(DiracIntegralTests, SlaterNuclear)
{
    StreamableObject::SetToPretty();
    int Z=cl->GetNuclearCharge();
    for (auto oi:sbs->Iterate<TOrbital_IBS<double> >())
    {
        SMatrix<double> Ven=oi->Nuclear(cl);

        const TIrrepBasisSet<double>* l=GetLarge(oi);
        const TIrrepBasisSet<double>* s=GetSmall(oi);
        SMatrix<double> VenLnum = -Z*mintegrator->Inv_r1(*l);
        SMatrix<double> VenSnum = -Z*mintegrator->Inv_r1(*s);
        SMat Vennum=merge_diag(VenLnum,VenSnum);
        //cout << "Ven=" << Ven << endl << "Ven num=" << Vennum << endl;
        //Because of the singularity at the origin, the error is larger than the other integrals.
        EXPECT_NEAR(Max(fabs(Ven-Vennum)),0.0,1e-11);        
    }
}


TEST_F(DiracIntegralTests, GaussianNuclear)
{
    StreamableObject::SetToPretty();
    int Z=cl->GetNuclearCharge();
    for (auto oi:gbs->Iterate<TOrbital_IBS<double> >())
    {
        SMatrix<double> Ven=oi->Nuclear(cl);

        const TIrrepBasisSet<double>* l=GetLarge(oi);
        const TIrrepBasisSet<double>* s=GetSmall(oi);
        SMatrix<double> VenLnum = -Z*mintegrator->Inv_r1(*l);
        SMatrix<double> VenSnum = -Z*mintegrator->Inv_r1(*s);
        SMat Vennum=merge_diag(VenLnum,VenSnum);
        //cout << "Ven=" << Ven << endl << "Ven num=" << Vennum << endl;
        // cout << "Ven=" << Ven << endl << "Ven1=" << Ven1 << endl;
        //Because of the singularity at the origin, the error is larger than the other integrals.
        EXPECT_NEAR(Max(fabs(Ven-Vennum)),0.0,1e-11);        
    }
}


TEST_F(DiracIntegralTests, SlaterKinetic)
{
    StreamableObject::SetToPretty();
    for (auto oi:sbs->Iterate<TOrbital_IBS<double> >())
    {
        SMatrix<double> K=oi->Kinetic();
        for (auto d:Vector<double>(K.GetDiagonal())) EXPECT_NEAR(d,0.0,1e-15);
        //cout << std::fixed << std::setprecision(3) << std::setw(6) << K << endl;
        const TIrrepBasisSet<double>* l=GetLarge(oi);
        const TIrrepBasisSet<double>* s=GetSmall(oi);
        Matrix<double> KnumL = mintegrator->Grada_b(*l,*s);
        SMat Knum=merge_off_diag(KnumL);
        // cout << "K=" << K << endl << "Knum=" << Knum << endl;
        EXPECT_NEAR(Max(fabs(K-Knum)),0.0,1e-11);      
    }
}
TEST_F(DiracIntegralTests, GaussianKinetic)
{
    StreamableObject::SetToPretty();
    for (auto oi:gbs->Iterate<TOrbital_IBS<double> >())
    {
        SMatrix<double> K=oi->Kinetic();
        // for (auto d:Vector<double>(K.GetDiagonal())) EXPECT_NEAR(d,0.0,1e-15);
        //cout << std::fixed << std::setprecision(3) << std::setw(6) << K << endl;
        const TIrrepBasisSet<double>* l=GetLarge(oi);
        const TIrrepBasisSet<double>* s=GetSmall(oi);
        Matrix<double> KnumL = mintegrator->Grada_b(*l,*s);
        SMat Knum=merge_off_diag(KnumL);
        //cout << "K=" << K << endl << "Knum=" << Knum << endl;
        EXPECT_NEAR(Max(fabs(K-Knum)),0.0,1e-11);      
        // EXPECT_NEAR(Max(fabs(K-K1)),0.0,1e-14);      
    }
}

 //  For symmetric matricies, Sum() should know to double all the off diagonal elements.
TEST_F(DiracIntegralTests, Contraction)
{
    size_t N=5;
    SMatrix<double> H(N),D(N);
    FillRandomPositive(H);
    FillRandomPositive(D);
    Matrix<double> Hf(H),Df(D);
    double c=Sum(DirectMultiply(H,D));
    double cf=Sum(DirectMultiply(Hf,Df));
    EXPECT_NEAR(c,cf,5e-15);
}

/*
TEST_F(DiracIntegralTests, Overlap3C)
{
    for (auto i=bs->beginT();i!=bs->end();i++)
    {
        ERI3 Sabc=ie->MakeOverlap3C(*i,*i);
        
        auto c=i->beginT();
        for (auto sab:Sabc)
        {
            SMatrix<double> Sabcnum = mintegrator->Overlap3C(**i,**c);
            EXPECT_NEAR(Max(fabs(sab-Sabcnum)),0.0,1e-8);
            c++;
        }
    }
}

TEST_F(DiracIntegralTests, Repulsion)
{
    for (auto i=bs->beginT();i!=bs->end();i++)
    {
        SMatrix<double> S=ie->MakeRepulsion(*i);
        for (auto j=i;j!=bs->end();j++)
        {
            Matrix<double> Sx=ie->MakeRepulsion(*i,*j);
            
        }
    }
}

TEST_F(DiracIntegralTests, Repulsion3C)
{
    for (auto i=bs->beginT();i!=bs->end();i++)
    {
        ERI3 Sabc=ie->MakeRepulsion3C(*i,*i);
    }
}


TEST_F(DiracIntegralTests, CoulombExchange)
{
    for (auto iabt=bs->beginT();iabt!=bs->end();iabt++)
    for (auto icdt=bs->beginT();icdt!=bs->end();icdt++)
    {
        const Slater::IrrepBasisSet* iab=dynamic_cast<const Slater::IrrepBasisSet*>(*iabt);
        const Slater::IrrepBasisSet* icd=dynamic_cast<const Slater::IrrepBasisSet*>(*icdt);
        int Nab=iab->GetNumFunctions(), Ncd=icd->GetNumFunctions();
        ERI4 J=ie->MakeDirect(*iabt,*icdt);
       
        for (int ia=1 ;ia<=Nab;ia++)
        for (int ib=ia;ib<=Nab;ib++)
        {
            SMatrix<double> Jab=J(ia,ib);
            for (int ic=1 ;ic<=Ncd;ic++)
            for (int id=ic;id<=Ncd;id++)
            {
                double norm=iab->ns(ia)*iab->ns(ib)*icd->ns(ic)*icd->ns(id);
                if (supported(*iab,*icd,ia,ib,ic,id))
                {
                        
                    double jv=Jab(ic,id)/norm, r0=R0(*iab,*icd,ia,ib,ic,id);
                    if (fabs(jv-r0)/jv>1e-12)
                    {
                        cout << "(a,b,c,d)=(" << ia << "," << ib << "," << ic << "," << id << ")" << endl;
                        cout << iab->GetSymmetry() << " " << icd->GetSymmetry() << endl; 
                        cout << "j,r=" << jv << " " << r0 << endl;
                        assert(false);                 
                    }
                    // cout << Jview(ia,ib,ic,id)/norm << " " << R0(*iab,*icd,ia,ib,ic,id) << endl;
                    double rerr=fabs(jv-r0)/jv;
                    EXPECT_NEAR(rerr,0.0,1e-13);
                }
            }
        }
    }
}
*/
