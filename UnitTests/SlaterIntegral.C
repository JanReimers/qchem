// File: ERIList.C  Test the DFT persistance classes


#include "gtest/gtest.h"
#include "Imp/Integrals/SlaterIntegrals.H"

#include "Imp/BasisSet/Slater/IntegralEngine.H"
#include "Imp/BasisSet/Slater/BasisSet.H"
#include "Imp/BasisSet/Slater/BasisFunction.H"
#include "Imp/BasisSet/SphericalGaussian/QuantumNumber.H"

//#include "Imp/BasisSet/Slater_m/IrrepBasisSet.H"
//#include "Imp/BasisSet/Slater_m/IntegralEngine.H"
//#include "Imp/BasisSet/Slater_m/BasisSet.H"
//#include "Imp/BasisSet/Slater_m/BasisFunction.H"

#include "Imp/BasisSet/SphericalGaussian/QuantumNumber.H"
#include "Imp/BasisSet/Slater/IrrepBasisSet.H"
#include "Imp/Integrals/MeshIntegrator.H"
#include "Imp/Misc/DFTDefines.H"
#include "Imp/Cluster/Atom.H"
#include "Imp/Cluster/Molecule.H"
#include "Imp/Containers/ERI4.H"
#include "Imp/Containers/ptr_vector.h"

#include <MeshParams.H>
#include <Cluster.H>
#include <BasisSet.H>
#include "oml/imp/ran250.h"
#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;


//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class SlaterRadialIntegralTests : public ::testing::Test
{
public:
    SlaterRadialIntegralTests()
    : Lmax(4    )
    , Z(1)
    , lap({qchem::Lapack,qchem::SVD,1e-6,1e-12})
    , ie(new Slater::IntegralEngine())
    , bs(new Slater::BasisSet(lap,6,0.1,10,Lmax))
    , cl(new Molecule())
    {
        StreamableObject::SetToPretty();
        cl->Insert(new Atom(Z,0.0,Vector3D(0,0,0)));
        MeshParams mp({qchem::MHL,200,3,2.0,qchem::Gauss,1,0,0,3});
        mintegrator=new MeshIntegrator<double>(cl->CreateMesh(mp));
        MeshParams rmp({qchem::MHL,200,3,2.0,qchem::Gauss,32,0,0,3});
        rmintegrator=new MeshIntegrator<double>(cl->CreateMesh(rmp));
        //cout << *bs << endl;
    }
    
    bool   supported(const Slater::IrrepIEClient&,const Slater::IrrepIEClient&,int ia, int ib, int ic, int id) const;
    double R0(const Slater::IrrepIEClient&,const Slater::IrrepIEClient&,int ia, int ib, int ic, int id) const;
    
    typedef AnalyticIE<double>::ERI3 ERI3;
    
    int Lmax, Z;
    LAParams lap;
    AnalyticIE<double>* ie;
    Slater::BasisSet* bs;
//    Slater_m::BasisSet* bsm;
    Cluster* cl;
    MeshIntegrator<double>* mintegrator;
    MeshIntegrator<double>* rmintegrator;
};

bool SlaterRadialIntegralTests::supported(const Slater::IrrepIEClient& ab, const Slater::IrrepIEClient& cd,int ia, int ib, int ic, int id) const
{
    int nab=ab.Ns(ia)+ab.Ns(ib);
    int ncd=cd.Ns(ic)+cd.Ns(id);
    return nab<=6 && ncd<=6;
}


double SlaterRadialIntegralTests::R0(const Slater::IrrepIEClient& ab, const Slater::IrrepIEClient& cd,int ia, int ib, int ic, int id) const
{
    double a=ab.es(ia)+ab.es(ib);
    double b=cd.es(ic)+cd.es(id);
    int nab=ab.Ns(ia)+ab.Ns(ib);
    int ncd=cd.Ns(ic)+cd.Ns(id);
    double f=4*4*Pi*Pi/(a*b*(a+b));
    if (nab==2 && ncd==2)
        return 2*f*( 
                    1/(pow(a,1)*pow(b,1)*pow(a+b,0))
                    +1/pow(a+b,2)
                    );
    if (nab==4 && ncd==2)
        return 12*f*( 
                        2/(pow(a,3)*pow(b,1)*pow(a+b,0))
                      + 2/(pow(a,0)*pow(b,0)*pow(a+b,4))
                      + 1/(pow(a,1)*pow(b,0)*pow(a+b,3))
                      - 1/(pow(a,3)*pow(b,0)*pow(a+b,1))
                        );
    if (nab==2 && ncd==4)
        return 12*f*( 
                       2/(pow(a,1)*pow(b,3)*pow(a+b,0))
                     + 2/(pow(a,0)*pow(b,0)*pow(a+b,4))
                     + 1/(pow(a,0)*pow(b,1)*pow(a+b,3))
                     - 1/(pow(a,0)*pow(b,3)*pow(a+b,1))
                    );
    if (nab==4 && ncd==4)
         return 144*f*( 
           1/(pow(a,2)*pow(b,2)*pow(a+b,2))
         + 2/(pow(a,1)*pow(b,1)*pow(a+b,4))
         + 5/(pow(a,0)*pow(b,0)*pow(a+b,6))
         + 1/(pow(a,3)*pow(b,3)*pow(a+b,0))
          );
          
    if (nab==2 && ncd==6)
        return 240*f*( 
                       1/(pow(a,1)*pow(b,5)*pow(a+b,0))
                     + 1/(pow(a,2)*pow(b,4)*pow(a+b,0))
                     + 3/(pow(a,0)*pow(b,0)*pow(a+b,6))
                     - 2/(pow(a,1)*pow(b,0)*pow(a+b,5))
                     - 1/(pow(a,2)*pow(b,0)*pow(a+b,4))
                    );
     if (nab==6 && ncd==2)
        return 240*f*( 
                       1/(pow(a,5)*pow(b,1)*pow(a+b,0))
                     + 1/(pow(a,4)*pow(b,2)*pow(a+b,0))
                     + 3/(pow(a,0)*pow(b,0)*pow(a+b,6))
                     - 2/(pow(a,0)*pow(b,1)*pow(a+b,5))
                     - 1/(pow(a,0)*pow(b,2)*pow(a+b,4))
                    );
    if (nab==4 && ncd==6)
        return 1440*f*( 
                        2/(pow(a,3)*pow(b,5)*pow(a+b,0))
                     +  2/(pow(a,4)*pow(b,4)*pow(a+b,0))
                     + 28/(pow(a,0)*pow(b,0)*pow(a+b,8))
                     -  7/(pow(a,1)*pow(b,0)*pow(a+b,7))
                     - 12/(pow(a,2)*pow(b,0)*pow(a+b,6))
                     -  7/(pow(a,3)*pow(b,0)*pow(a+b,5))
                     -  2/(pow(a,4)*pow(b,0)*pow(a+b,4))
                    );
    if (nab==6 && ncd==4)
        return 1440*f*( 
                        2/(pow(a,5)*pow(b,3)*pow(a+b,0))
                     +  2/(pow(a,4)*pow(b,4)*pow(a+b,0))
                     + 28/(pow(a,0)*pow(b,0)*pow(a+b,8))
                     -  7/(pow(a,0)*pow(b,1)*pow(a+b,7))
                     - 12/(pow(a,0)*pow(b,2)*pow(a+b,6))
                     -  7/(pow(a,0)*pow(b,3)*pow(a+b,5))
                     -  2/(pow(a,0)*pow(b,4)*pow(a+b,4))
                    );
    
    if (nab==6 && ncd==6)
        return 86400*f*( 
                        1/(pow(a,5)*pow(b,5)*pow(a+b,0))
                     +  1/(pow(a,6)*pow(b,4)*pow(a+b,0))
                     + 42/(pow(a,0)*pow(b,0)*pow(a+b,10))
                     - 14/(pow(a,2)*pow(b,0)*pow(a+b,8))
                     - 14/(pow(a,3)*pow(b,0)*pow(a+b,7))
                     -  9/(pow(a,4)*pow(b,0)*pow(a+b,6))
                     -  4/(pow(a,5)*pow(b,0)*pow(a+b,5))
                     -  1/(pow(a,6)*pow(b,0)*pow(a+b,4))
                    );
    
    
                    
    assert(false);
    return 0;

}

TEST_F(SlaterRadialIntegralTests, Overlap)
{
    for (auto i=bs->beginT();i!=bs->end();i++)
    {
        SMatrix<double> S=ie->MakeOverlap(*i);
        for (auto d:Vector<double>(S.GetDiagonal())) EXPECT_NEAR(d,1.0,1e-15);
        //cout << S << endl;
        SMatrix<double> Snum = mintegrator->Overlap(**i);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,1e-8);

    }
}

TEST_F(SlaterRadialIntegralTests, Nuclear)
{
    for (auto i=bs->beginT();i!=bs->end();i++)
    {
        SMatrix<double> Hn=ie->MakeNuclear(*i,*cl);
        SMatrix<double> Hnnum = -1*mintegrator->Nuclear(**i);
        EXPECT_NEAR(Max(fabs(Hn-Hnnum)),0.0,1e-7);

    }
}

TEST_F(SlaterRadialIntegralTests, Kinetic)
{
    for (auto i=bs->beginT();i!=bs->end();i++)
    {
        SMatrix<double> K=ie->MakeKinetic(*i);
        //cout << S << endl;
        SMatrix<double> Knum = 0.5*mintegrator->Grad(**i);
            // We need to add the l*(l+1) term that comes from the angular integrals.
        // Lost of dynamic cast just to get at L!
        const QuantumNumber& qn=i->GetQuantumNumber();
        const SphericalSymmetryQN& sqn=dynamic_cast<const SphericalSymmetryQN& >(qn);
        int l=sqn.GetL();
        const Slater::IrrepBasisSet* sg=dynamic_cast<const Slater::IrrepBasisSet*>(*i);
        assert(sg);
        int n=2*l+2;
        for (auto i:Knum.rows())
            for (auto j:Knum.cols(i))
                Knum(i,j)+=0.5*(l*(l+1))*SlaterIntegral(sg->es(i)+sg->es(j),n-2)*sg->ns(i)*sg->ns(j);
            
        EXPECT_NEAR(Max(fabs(K-Knum)),0.0,1e-10);

    }
}

TEST_F(SlaterRadialIntegralTests, Overlap3C)
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

TEST_F(SlaterRadialIntegralTests, Repulsion)
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

TEST_F(SlaterRadialIntegralTests, Repulsion3C)
{
    for (auto i=bs->beginT();i!=bs->end();i++)
    {
        ERI3 Sabc=ie->MakeRepulsion3C(*i,*i);
    }
}

struct Vf : public VectorFunction<double>
{
    typedef VectorFunction<double> Base;
     virtual size_t  GetVectorSize() const
     {
         return bfs.size();
     }
    virtual Vec operator()(const RVec3& r) const
    {
        Vec ret(bfs.size());
        int i=1;
        for (auto bf:bfs) ret(i++)=(*bf)(r);
        return ret;
    }
    virtual Vec3Vec  Gradient  (const RVec3&) const
    {
        return Vec3Vec();
    }
    typedef TBasisFunction<double> bf_t;
    optr_vector1<bf_t*> bfs;
};
//            vfab.bfs.push_back(new Slater::BasisFunction(i->es(ia)+i->es(ib),i->Ns(ia)+i->Ns(ib),i->Ls(ia)+i->Ls(ib),i->ns(ia)*i->ns(ib)));
//                vfcd.bfs.push_back(new Slater::BasisFunction(i->es(ic)+i->es(id),i->Ns(ic)+i->Ns(id),i->Ls(ic)+i->Ls(id),i->ns(ic)*i->ns(id)));
//                Matrix<double> R=mintegrator->Repulsion(vfab,vfcd);


TEST_F(SlaterRadialIntegralTests, CoulombExchange)
{
    ERI4 J,K;
    J.SetSize(bs->size(),0.0);
    K.SetSize(bs->size(),0.0);
    ie->Make4C(J,K,bs);
    //int L=0;
    for (auto iabt=bs->beginT();iabt!=bs->end();iabt++)
    for (auto icdt=bs->beginT();icdt!=bs->end();icdt++)
    {
        const Slater::IrrepBasisSet* iab=dynamic_cast<const Slater::IrrepBasisSet*>(*iabt);
        const Slater::IrrepBasisSet* icd=dynamic_cast<const Slater::IrrepBasisSet*>(*icdt);
        int Nab=iab->GetNumFunctions(), Ncd=icd->GetNumFunctions();
        ERI4view Jview(J,iab->GetStartIndex(),icd->GetStartIndex());
       
        for (int ia=1 ;ia<=Nab;ia++)
        for (int ib=ia;ib<=Nab;ib++)
        {
            for (int ic=1 ;ic<=Ncd;ic++)
            for (int id=ic;id<=Ncd;id++)
            {
                double norm=iab->ns(ia)*iab->ns(ib)*icd->ns(ic)*icd->ns(id);
                if (supported(*iab,*icd,ia,ib,ic,id))
                {
                        
                    double jv=Jview(ia,ib,ic,id)/norm, r0=R0(*iab,*icd,ia,ib,ic,id);
                    if (fabs(jv-r0)/jv>1e-12)
                    {
                        cout << "(a,b,c,d)=(" << ia << "," << ib << "," << ic << "," << id << ")" << endl;
                        cout << iab->GetQuantumNumber() << " " << icd->GetQuantumNumber() << endl; 
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

TEST_F(SlaterRadialIntegralTests, YlmExchange)
{
     double eps=1e-14;
     SlaterRadialIntegrals S(0.5,0.5);
     int la=0,lb=0,lc=0,ld=0;

     double ssss=S.DoExchangeSum(la,lb,lc,ld);
     double ssss_m=S.DoExchangeSum(la,lb,lc,ld,0,0,0,0);
     EXPECT_NEAR(ssss,ssss_m,ssss*eps);
     
     lb=ld=1;
     double spsp=S.DoExchangeSum(la,lb,lc,ld);
     double spsp_m=0.0;
     for (int mb=-1;mb<=1;mb++)
        spsp_m+=S.DoExchangeSum(la,lb,lc,ld,0,mb,0,mb);
    
     spsp_m*=1.0/(2*la+1)/(2*lb+1);
     EXPECT_NEAR(spsp,spsp_m,spsp*eps);
     
     la=lc=1;
     lb=ld=0;
     double psps=S.DoExchangeSum(la,lb,lc,ld);
     double psps_m=0.0;
     for (int ma=-1;ma<=1;ma++)
        psps_m+=S.DoExchangeSum(la,lb,lc,ld,ma,0,ma,0);
     
     psps_m*=1.0/(2*la+1)/(2*lb+1);
     EXPECT_NEAR(psps,psps_m,psps*eps);

     la=ld=0;
     lb=lc=1;
     double spps=S.DoExchangeSum(la,lb,lc,ld);
     double spps_m=0.0;
     for (int mb=-1;mb<=1;mb++)
        spps_m+=S.DoExchangeSum(la,lb,lc,ld,0,mb,mb,0);
     spps_m*=1.0/(2*la+1)/(2*lb+1);
    
     EXPECT_NEAR(spps,spps_m,spps*eps);
     
     la=ld=1;
     lb=lc=0;
     double pssp=S.DoExchangeSum(la,lb,lc,ld);
     double pssp_m=0.0;
     for (int ma=-1;ma<=1;ma++)
        pssp_m+=S.DoExchangeSum(la,lb,lc,ld,ma,0,0,ma);
     pssp_m*=1.0/(2*la+1)/(2*lb+1);
    
     EXPECT_NEAR(pssp,pssp_m,pssp*eps);
     
     la=lb=lc=ld=1;
     double pppp=S.DoExchangeSum(la,lb,lc,ld);
     double pppp_mac=0.0,pppp_mad=0.0;
     for (int ma=-1;ma<=1;ma++)
     for (int mb=-1;mb<=1;mb++)
     {
        pppp_mac+=S.DoExchangeSum(la,lb,lc,ld,ma,mb,ma,mb);
        pppp_mad+=S.DoExchangeSum(la,lb,lc,ld,ma,mb,mb,ma);         
     }
     pppp_mac*=1.0/(2*la+1)/(2*lb+1);
     pppp_mad*=1.0/(2*la+1)/(2*lb+1);
       
     EXPECT_NEAR(pppp,pppp_mac,pppp*eps);
     EXPECT_NEAR(pppp,pppp_mad,pppp*eps);
//
//  d-s
// 
{
    la=lc=0;
    lb=ld=2;
    double sdsd=S.DoExchangeSum(la,lb,lc,ld);
    double sdsd_m=0.0;
    for (int mb=-lb;mb<=lb;mb++)
        sdsd_m+=S.DoExchangeSum(la,lb,lc,ld,0,mb,0,mb);

    sdsd_m*=1.0/(2*la+1)/(2*lb+1);
    EXPECT_NEAR(sdsd,sdsd_m,sdsd*eps);

    la=lc=2;
    lb=ld=0;
    double dsds=S.DoExchangeSum(la,lb,lc,ld);
    double dsds_m=0.0;
    for (int ma=-la;ma<=la;ma++)
        dsds_m+=S.DoExchangeSum(la,lb,lc,ld,ma,0,ma,0);

    dsds_m*=1.0/(2*la+1)/(2*lb+1);
    EXPECT_NEAR(dsds,dsds_m,dsds*eps);

    la=ld=0;
    lb=lc=2;
    double sdds=S.DoExchangeSum(la,lb,lc,ld);
    double sdds_m=0.0;
    for (int mb=-lb;mb<=lb;mb++)
        sdds_m+=S.DoExchangeSum(la,lb,lc,ld,0,mb,mb,0);
    sdds_m*=1.0/(2*la+1)/(2*lb+1);

    EXPECT_NEAR(sdds,sdds_m,sdds*eps);

    la=ld=2;
    lb=lc=0;
    double dssd=S.DoExchangeSum(la,lb,lc,ld);
    double dssd_m=0.0;
    for (int ma=-la;ma<=la;ma++)
        dssd_m+=S.DoExchangeSum(la,lb,lc,ld,ma,0,0,ma);
    dssd_m*=1.0/(2*la+1)/(2*lb+1);

    EXPECT_NEAR(dssd,dssd_m,dssd*eps);
}
//
//  d-p
//
    {
    
    la=lc=1;
    lb=ld=2;
    double pdpd=S.DoExchangeSum(la,lb,lc,ld);
    double pdpd_m=0.0;
    for (int ma=-la;ma<=la;ma++)
        for (int mb=-lb;mb<=lb;mb++)
            pdpd_m+=S.DoExchangeSum(la,lb,lc,ld,ma,mb,ma,mb);

    pdpd_m*=1.0/(2*la+1)/(2*lb+1);
    EXPECT_NEAR(pdpd,pdpd_m,pdpd*eps);

    la=lc=2;
    lb=ld=1;
    double dpdp=S.DoExchangeSum(la,lb,lc,ld);
    double dpdp_m=0.0;
    for (int ma=-la;ma<=la;ma++)
        for (int mb=-lb;mb<=lb;mb++)
            dpdp_m+=S.DoExchangeSum(la,lb,lc,ld,ma,mb,ma,mb);

    dpdp_m*=1.0/(2*la+1)/(2*lb+1);
    EXPECT_NEAR(dpdp,dpdp_m,dpdp*eps);

    la=ld=1;
    lb=lc=2;
    double pddp=S.DoExchangeSum(la,lb,lc,ld);
    double pddp_m=0.0;
    for (int ma=-la;ma<=la;ma++)
        for (int mb=-lb;mb<=lb;mb++)
            pddp_m+=S.DoExchangeSum(la,lb,lc,ld,ma,mb,mb,ma);
    pddp_m*=1.0/(2*la+1)/(2*lb+1);

    EXPECT_NEAR(pddp,pddp_m,pddp*eps);

    la=ld=2;
    lb=lc=1;
    double dppd=S.DoExchangeSum(la,lb,lc,ld);
    double dppd_m=0.0;
    for (int ma=-la;ma<=la;ma++)
        for (int mb=-lb;mb<=lb;mb++)
            dppd_m+=S.DoExchangeSum(la,lb,lc,ld,ma,mb,mb,ma);
    dppd_m*=1.0/(2*la+1)/(2*lb+1);

    EXPECT_NEAR(dppd,dppd_m,dppd*eps);
    }

//
//  d-d
//
    la=ld=2;
    lb=lc=2;
    double dddd=S.DoExchangeSum(la,lb,lc,ld);
    double dddd_mac=0.0,dddd_mad=0.0;
    for (int ma=-la;ma<=la;ma++)
        for (int mb=-lb;mb<=lb;mb++)
        {
            dddd_mac+=S.DoExchangeSum(la,lb,lc,ld,ma,mb,ma,mb);
            dddd_mad+=S.DoExchangeSum(la,lb,lc,ld,ma,mb,mb,ma);
            
        }
    dddd_mac*=1.0/(2*la+1)/(2*lb+1);
    dddd_mad*=1.0/(2*la+1)/(2*lb+1);

    EXPECT_NEAR(dddd,dddd_mac,dddd*eps);
    EXPECT_NEAR(dddd,dddd_mad,dddd*eps);
    
}

TEST_F(SlaterRadialIntegralTests, YlmCoulomb)
{
    double eps=1e-14;
    SlaterRadialIntegrals S(0.5,0.5);
    int la=0,lb=0,lc=0,ld=0;
    double ssss=S.Coulomb(la,lb,lc,ld);
    double ssss_m=S.Coulomb(la,lb,lc,ld,0,0,0,0);
    EXPECT_NEAR(ssss,ssss_m,ssss*eps);

    lc=ld=1;
    double sspp=S.Coulomb(la,lb,lc,ld);
    double sspp_m=0.0;
    for (int mb=-1;mb<=1;mb++)
        sspp_m+=S.Coulomb(la,lb,lc,ld,0,0,mb,mb);
    sspp_m*=1.0/(2*la+1)/(2*lc+1);

    EXPECT_NEAR(sspp,sspp_m,sspp*eps);

    la=lb=1;
    lc=ld=0;
    double ppss=S.Coulomb(la,lb,lc,ld);
    double ppss_m=0.0;
    for (int ma=-1;ma<=1;ma++)
        ppss_m+=S.Coulomb(la,lb,lc,ld,ma,ma,0,0);
    ppss_m*=1.0/(2*la+1)/(2*lc+1);

    EXPECT_NEAR(ppss,ppss_m,ppss*eps);


    lc=ld=1;
    double pppp=S.Coulomb(la,lb,lc,ld);
    double pppp_m=0.0;
    for (int ma=-la;ma<=la;ma++)
    for (int mc=-lc;mc<=lc;mc++)
        pppp_m+=S.Coulomb(la,lb,lc,ld,ma,ma,mc,mc);

    pppp_m*=1.0/(2*la+1)/(2*lc+1);       

    EXPECT_NEAR(pppp,pppp_m,pppp*eps);

    la=lb=0;
    lc=ld=2;
    double ssdd=S.Coulomb(la,lb,lc,ld);
    double ssdd_m=0.0;
    for (int mc=-lc;mc<=lc;mc++)
        ssdd_m+=S.Coulomb(la,lb,lc,ld,0,0,mc,mc);
    ssdd_m*=1.0/(2*la+1)/(2*lc+1);

    EXPECT_NEAR(ssdd,ssdd_m,ssdd*eps);
   
    la=lb=2;
    lc=ld=0;
    double ddss=S.Coulomb(la,lb,lc,ld);
    double ddss_m=0.0;
    for (int ma=-la;ma<=la;ma++)
        ddss_m+=S.Coulomb(la,lb,lc,ld,ma,ma,0,0);
    ddss_m*=1.0/(2*la+1)/(2*lc+1);

    EXPECT_NEAR(ddss,ddss_m,ddss*eps);
    
    la=lb=1;
    lc=ld=2;
    double ppdd=S.Coulomb(la,lb,lc,ld);
    double ppdd_m=0.0;
    for (int ma=-la;ma<=la;ma++)
        for (int mc=-lc;mc<=lc;mc++)
            ppdd_m+=S.Coulomb(la,lb,lc,ld,ma,ma,mc,mc);
    ppdd_m*=1.0/(2*la+1)/(2*lc+1);

    EXPECT_NEAR(ppdd,ppdd_m,ppdd*eps);
   
    la=lb=2;
    lc=ld=1;
    double ddpp=S.Coulomb(la,lb,lc,ld);
    double ddpp_m=0.0;
    for (int ma=-la;ma<=la;ma++)
        for (int mc=-lc;mc<=lc;mc++)
            ddpp_m+=S.Coulomb(la,lb,lc,ld,ma,ma,mc,mc);
    ddpp_m*=1.0/(2*la+1)/(2*lc+1);

    EXPECT_NEAR(ddpp,ddpp_m,ddpp*eps);
    
    la=lb=2;
    lc=ld=2;
    double dddd=S.Coulomb(la,lb,lc,ld);
    double dddd_m=0.0;
    for (int ma=-la;ma<=la;ma++)
        for (int mc=-lc;mc<=lc;mc++)
            dddd_m+=S.Coulomb(la,lb,lc,ld,ma,ma,mc,mc);
    dddd_m*=1.0/(2*la+1)/(2*lc+1);

    EXPECT_NEAR(dddd,dddd_m,dddd*eps);
}
//
//
//#ifndef DEBUG 
//#include "oml/io3d.h"
//TEST_F(SlaterRadialIntegralTests, Numerical)
//{
//    TIrrepBasisSet<double>* vf=*bs->beginT();
//    TBasisFunction<double>* sf=*vf->beginT();
//
//    Vector<double> cnum=mintegrator->Integrate(*vf);
//    Vector<double> c=ie->MakeCharge(vf);
//    EXPECT_NEAR(Max(fabs(cnum-c)),0.0,1e-13);
//
//    Vector<double> nnum=mintegrator->Normalize(*vf);
//    EXPECT_NEAR(Max(fabs(nnum-1.0)),0.0,1e-12);
//    {
//        Vector <double> onum=mintegrator->Overlap(*sf,*vf);
//        SMatrix<double> o=ie->MakeOverlap(vf);
//        EXPECT_NEAR(Max(fabs(onum-o.GetRow(1))),0.0,1e-12);        
//    }
//    {
//        Matrix<double> onum=mintegrator->Overlap(*vf,*vf);
//        SMatrix<double> o=ie->MakeOverlap(vf);
//        EXPECT_NEAR(Max(fabs(onum-o)),0.0,1e-12); 
//    }
//    {
//        SMatrix<double> onum=mintegrator->Overlap3C(*vf,*sf);
//        ERI3 o=ie->MakeOverlap3C(vf,vf);
//        EXPECT_NEAR(Max(fabs(onum-o[0])),0.0,1e-12); 
//    }
//    {
//        SMatrix<double> rnum=rmintegrator->Repulsion(*vf);
//        SMatrix<double> r=ie->MakeRepulsion(vf);
//        EXPECT_NEAR(Max(fabs(DirectDivide(rnum-r,r))),0.0,0.1); 
//    }
//    {
//        Vector<double> rnum=rmintegrator->Repulsion(*sf,*vf);
//        Vector<double> r=ie->MakeRepulsion(vf).GetRow(1);
//        EXPECT_NEAR(Max(fabs(DirectDivide(rnum-r,r))),0.0,0.1); 
//    }
//    {
//        SMatrix<double> rnum=rmintegrator->Repulsion(*vf,*vf);
//        SMatrix<double> r=ie->MakeRepulsion(vf,vf);
//        EXPECT_NEAR(Max(fabs(DirectDivide(rnum-r,r))),0.0,0.1); 
//    }
//    {
//        SMatrix<double> rnum=rmintegrator->Repulsion3C(*vf,*sf);
//        ERI3 r=ie->MakeRepulsion3C(vf,vf);
//        EXPECT_NEAR(Max(fabs(DirectDivide(rnum-r[0],r[0]))),0.0,0.1); 
//    }
//}
// #endif
