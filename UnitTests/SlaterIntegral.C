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
    : Lmax(1    )
    , Z(1)
    , lap({qchem::Lapack,qchem::SVD,1e-6,1e-12})
    , ie(new Slater::IntegralEngine())
    , db(new HeapDB<double>())
    , ibs(new Slater::IrrepBasisSet(lap,db,2,0.5,1,Lmax))
    , bs(new Slater::BasisSet(lap,2,0.5,1,Lmax))
    , mesh(0)
    , cl(new Molecule())
    , mintegrator()
    {
        StreamableObject::SetToPretty();
        RadialMesh*  rm=new MHLRadialMesh(50,2U,1.0); //mem leak
        AngularMesh* am=new GaussAngularMesh(8);      //mem leak
        mesh=new AtomMesh(*rm,*am); 
        mintegrator=new MeshIntegrator<double>(mesh);
        cl->Insert(new Atom(Z,0.0,Vector3D(0,0,0)));
    }
    
    double R0(const Slater::IrrepIEClient&,const Slater::IrrepIEClient&,int ia, int ib, int ic, int id) const;
    
    typedef AnalyticIE<double>::ERI3 ERI3;
    
    int Lmax, Z;
    LAParams lap;
    AnalyticIE<double>* ie;
    IntegralDataBase<double>* db;
    Slater::IrrepBasisSet* ibs;
    Slater::BasisSet* bs;
    Mesh* mesh;
    Cluster* cl;
    MeshIntegrator<double>* mintegrator;
};


//inline double r(double a, double b, double apb,int )
//{
//    
//}
double SlaterRadialIntegralTests::R0(const Slater::IrrepIEClient& ab, const Slater::IrrepIEClient& cd,int ia, int ib, int ic, int id) const
{
    double a=ab.es(ia)+ab.es(ib);
    double b=cd.es(ic)+cd.es(id);
    int nab=ab.Ns(ia)+ab.Ns(ib);
    int ncd=cd.Ns(ic)+cd.Ns(id);
    double f=1.0/(a*b*(a+b));
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
    assert(false);
    return 0;

}
//  double R02222=144/(a*b*(a+b))
//                    *(
//                           1/(pow(a,2)*pow(b,2)*pow(a+b,2))
//                         + 2/(pow(a,1)*pow(b,1)*pow(a+b,4))
//                         + 5/(pow(a,0)*pow(b,0)*pow(a+b,6))
//                         + 1/(pow(a,3)*pow(b,3)*pow(a+b,0))
//                      );
//                    double R01122=12/(a*b*(a+b))
//                    *(
//                           2/(pow(a,1)*pow(b,3)*pow(a+b,0))
//                         + 2/(pow(a,0)*pow(b,0)*pow(a+b,4))
//                         + 1/(pow(a,0)*pow(b,1)*pow(a+b,3))
//                         - 1/(pow(a,0)*pow(b,3)*pow(a+b,1))
//                      );
//                    double R02211=12/(a*b*(a+b))
//                    *(
//                           2/(pow(a,3)*pow(b,1)*pow(a+b,0))
//                         + 2/(pow(a,0)*pow(b,0)*pow(a+b,4))
//                         + 1/(pow(a,1)*pow(b,0)*pow(a+b,3))
//                         - 1/(pow(a,3)*pow(b,0)*pow(a+b,1))
//                      );
//
//TEST_F(SlaterRadialIntegralTests, Overlap)
//{
//    for (auto i=bs->beginT();i!=bs->end();i++)
//    {
//        SMatrix<double> S=ie->MakeOverlap(*i);
//        for (auto d:Vector<double>(S.GetDiagonal())) EXPECT_NEAR(d,1.0,1e-15);
//        //cout << S << endl;
//        SMatrix<double> Snum = mintegrator->Overlap(**i);
//        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,1e-8);
//
//    }
//}
//
//TEST_F(SlaterRadialIntegralTests, Nuclear)
//{
//    for (auto i=bs->beginT();i!=bs->end();i++)
//    {
//        SMatrix<double> Hn=ie->MakeNuclear(*i,*cl);
//        //cout << S << endl;
//        SMatrix<double> Hnnum = -1*mintegrator->Nuclear(**i);
//        EXPECT_NEAR(Max(fabs(Hn-Hnnum)),0.0,1e-8);
//
//    }
//}
//
//TEST_F(SlaterRadialIntegralTests, Kinetic)
//{
//    for (auto i=bs->beginT();i!=bs->end();i++)
//    {
//        SMatrix<double> K=ie->MakeKinetic(*i);
//        //cout << S << endl;
//        SMatrix<double> Knum = 0.5*mintegrator->Grad(**i);
//            // We need to add the l*(l+1) term that comes from the angular integrals.
//        // Lost of dynamic cast just to get at L!
//        const QuantumNumber& qn=i->GetQuantumNumber();
//        const SphericalSymmetryQN& sqn=dynamic_cast<const SphericalSymmetryQN& >(qn);
//        int l=sqn.GetL();
//        const Slater::IrrepBasisSet* sg=dynamic_cast<const Slater::IrrepBasisSet*>(*i);
//        assert(sg);
//        int n=2*l+2;
//        for (auto i:Knum.rows())
//            for (auto j:Knum.cols(i))
//                Knum(i,j)+=0.5*(l*(l+1))*SlaterIntegral(sg->es(i)+sg->es(j),n-2)*sg->ns(i)*sg->ns(j);
//            
//        EXPECT_NEAR(Max(fabs(K-Knum)),0.0,1e-11);
//
//    }
//}
//
//TEST_F(SlaterRadialIntegralTests, Overlap3C)
//{
//    for (auto i=bs->beginT();i!=bs->end();i++)
//    {
//        ERI3 Sabc=ie->MakeOverlap3C(*i,*i);
//        
//        auto c=i->beginT();
//        for (auto sab:Sabc)
//        {
//            SMatrix<double> Sabcnum = mintegrator->Overlap3C(**i,**c);
//            EXPECT_NEAR(Max(fabs(sab-Sabcnum)),0.0,1e-8);
//            c++;
//        }
//    }
//}

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
                //double a=iab->es(ia)+iab->es(ib) ,b=icd->es(ic)+icd->es(id);
                //int nab=i->Ns(ia)+i->Ns(ib) ,ncd=i->Ns(ic)+i->Ns(id);
                //double R01111=2.0/(a*b*(a+b))*(1/(a*b)+1/pow(a+b,2));
              
                //cout << "(a,b,c,d)=(" << ia << "," << ib << "," << ic << "," << id << ")" << endl;
               // cout << Jview(ia,ib,ic,id)/norm << " " << R0(*iab,*icd,ia,ib,ic,id) << endl;
                EXPECT_NEAR(Jview(ia,ib,ic,id)/norm,R0(*iab,*icd,ia,ib,ic,id),1e-13);
            }
        }
    }
}

