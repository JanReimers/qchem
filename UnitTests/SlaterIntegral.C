// File: ERIList.C  Test the DFT persistance classes


#include "gtest/gtest.h"
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/BasisSet/Slater/IntegralEngine.H"
#include "Imp/BasisSet/Slater/BasisSet.H"
#include "Imp/BasisSet/Slater/BasisFunction.H"
#include "Imp/BasisSet/SphericalGaussian/QuantumNumber.H"

//#include "Imp/BasisSet/Slater/IEClient.H"
#include "Imp/BasisSet/Slater/IrrepBasisSet.H"
#include "DFTDataBase/HeapDB/HeapDB.H"
//#include "Imp/Integrals/SlaterRadialIntegrals.H"
//#include "Imp/Integrals/Wigner3j.H"
#include "Mesh/RadialMesh/MHLRadialMesh.H"
#include "Mesh/AngularMesh/GaussAngularMesh.H"
#include "Mesh/AtomMesh.H"
#include "Mesh/MeshIntegrator.H"
#include "Misc/DFTDefines.H"
//#include "Cluster/Atom.H"
#include "Cluster/Molecule.H"
#include "Cluster.H"
#include "BasisSet.H"
#include "Imp/Containers/ERI4.H"
#include "oml/imp/ran250.h"
#include "Imp/Containers/ptr_vector.h"
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
    , db(new HeapDB<double>())
    , ibs(new Slater::IrrepBasisSet(lap,db,6,0.1,10,Lmax))
    , bs(new Slater::BasisSet(lap,6,0.1,10,Lmax))
    , mesh(0)
    , cl(new Molecule())
    , mintegrator()
    {
        StreamableObject::SetToPretty();
        {
            RadialMesh*  rm=new MHLRadialMesh(200,2U,1.0); //mem leak
            AngularMesh* am=new GaussAngularMesh(1);      //mem leak
            mesh=new AtomMesh(*rm,*am); 
            mintegrator=new MeshIntegrator<double>(mesh);
        }
        {
            RadialMesh*  rm=new MHLRadialMesh(50,1U,1.0); //mem leak
            AngularMesh* am=new GaussAngularMesh(1);      //mem leak
            rmesh=new AtomMesh(*rm,*am); 
            rmintegrator=new MeshIntegrator<double>(rmesh);
       }
        cl->Insert(new Atom(Z,0.0,Vector3D(0,0,0)));
    }
    
    bool   supported(const Slater::IrrepIEClient&,const Slater::IrrepIEClient&,int ia, int ib, int ic, int id) const;
    double R0(const Slater::IrrepIEClient&,const Slater::IrrepIEClient&,int ia, int ib, int ic, int id) const;
    
    typedef AnalyticIE<double>::ERI3 ERI3;
    
    int Lmax, Z;
    LAParams lap;
    AnalyticIE<double>* ie;
    IntegralDataBase<double>* db;
    Slater::IrrepBasisSet* ibs;
    Slater::BasisSet* bs;
    Mesh* mesh;
    Mesh* rmesh;
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
        //cout << S << endl;
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

