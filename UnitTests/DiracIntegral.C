// File: DiracIntegral.C  Test the Dirac basis sets and integral engine.


#include "gtest/gtest.h"

#include "Imp/BasisSet/Slater_mj/IntegralEngine.H"
#include "Imp/BasisSet/Slater_mj/BasisSet.H"
#include "Imp/BasisSet/Slater_mj/BasisFunction.H"
#include "Imp/BasisSet/Slater_mj/IrrepBasisSet.H"
#include "Imp/Symmetry/OkmjQN.H"

#include "Imp/Misc/DFTDefines.H"
#include "Imp/Cluster/Atom.H"
#include "Imp/Cluster/Molecule.H"
#include "Imp/Containers/ERI4.H"
#include "Imp/Containers/ptr_vector.h"

#include <Cluster.H>
#include <BasisSet.H>
#include <iostream>

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
    : Lmax(0    )
    , Z(1)
    , lap({qchem::Lapack,qchem::SVD,1e-6,1e-12})
    , bs(new Slater_mj::Dirac_BasisSet(lap,3,0.1,10,Lmax))
    , ie(bs->itsIE)
    , cl(new Molecule())
    {
        cl->Insert(new Atom(Z,0.0,Vector3D(0,0,0)));
    }
    
    int Lmax, Z;
    LAParams lap;
    Slater_mj::Dirac_BasisSet* bs;
    AnalyticIE<double>* ie;
    Cluster* cl;
};

TEST_F(DiracIntegralTests, BasisSet)
{
    StreamableObject::SetToPretty();
    cout << *bs << endl;
}

/*
TEST_F(DiracIntegralTests, Overlap)
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

TEST_F(DiracIntegralTests, Nuclear)
{
    for (auto i=bs->beginT();i!=bs->end();i++)
    {
        SMatrix<double> Hn=ie->MakeNuclear(*i,*cl);
        SMatrix<double> Hnnum = -1*mintegrator->Nuclear(**i);
        EXPECT_NEAR(Max(fabs(Hn-Hnnum)),0.0,1e-7);

    }
}

TEST_F(DiracIntegralTests, Kinetic)
{
    for (auto i=bs->beginT();i!=bs->end();i++)
    {
        SMatrix<double> K=ie->MakeKinetic(*i);
        //cout << S << endl;
        SMatrix<double> Knum = 0.5*mintegrator->Grad(**i);
            // We need to add the l*(l+1) term that comes from the angular integrals.
        // Lost of dynamic cast just to get at L!
        const QuantumNumber& qn=i->GetQuantumNumber();
        const YlQN& sqn=dynamic_cast<const YlQN& >(qn);
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
*/
