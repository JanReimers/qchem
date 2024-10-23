// File MeshIntegrals.C  Run through varios mesh types and parameters

#include "gtest/gtest.h"

#include "Imp/BasisSet/PolarizedGaussian/IntegralEngine.H"
#include "Imp/BasisSet/PolarizedGaussian/BasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/IrrepBasisSet.H"
#include "Imp/Integrals/MeshIntegrator.H"
#include "Imp/Cluster/Molecule.H"
#include "Imp/Cluster/Atom.H"
#include <LAParams.H>
#include <MeshParams.H>
#include <QuantumNumber.H>
#include "oml/smatrix.h"
#include "oml/matrix.h"
#include <iomanip>


using std::cout;
using std::endl;
using std::setw;

double norm1(const SMatrix<double>& m)
{
    return sqrt(Sum(DirectMultiply(m,m)));
}
//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class MeshIntegralsTests : public ::testing::Test
{
public:
    MeshIntegralsTests()
    : Lmax(4   )
    , Z(1)
    , lap({qchem::Lapack,qchem::SVD,1e-6,1e-12})
    , ie(new PolarizedGaussian::IntegralEngine())
    , bs(0)
    , cl(new Molecule())
    {
        
        
    }
    
    void InitAtom()
    {
        cl->Insert(new Atom(Z,0.0,Vector3D(0,0,0)));
        bs=new PolarizedGaussian::BasisSet(lap,3,.1,10.0,Lmax,cl);
    }
    
    void InitMolecule()
    {
        cl->Insert(new Atom(Z,0.0,Vector3D( 1.,0.,0.)));
        cl->Insert(new Atom(Z,0.0,Vector3D(-1.,0.,0.)));
        bs=new PolarizedGaussian::BasisSet(lap,3,.1,10.0,Lmax,cl);        
    }
    
    int Lmax, Z;
    LAParams lap;
    AnalyticIE<double>* ie;
    PolarizedGaussian::BasisSet* bs;
    Cluster* cl;
};

TEST_F(MeshIntegralsTests, GaussAngles)
{
    StreamableObject::SetToPretty();
    InitAtom();
    cout << *bs << endl;
    //cout.precision(2);
    size_t Nradial=100;
    cout << "alpha m      s        p         d"  << endl;
    for (size_t mMHL=2;mMHL<=4;mMHL++)
    for (int aMHL=1;aMHL<=4;aMHL++)
    {
        double alpha=0.5*aMHL;
        MeshParams mp({qchem::MHL,Nradial,mMHL,alpha,qchem::Gauss,12,0,0});
        MeshIntegrator<double> mi(cl->CreateMesh(mp));
        for (auto ibs=bs->beginT();ibs!=bs->end();ibs++)
        {
            SMatrix<double> o= ie->MakeOverlap(*ibs);
            SMatrix<double> on= mi.Overlap(**ibs);
            double err=norm1(o-on);
            cout << std::fixed << setw(3) << alpha << " " << setw(2) << mMHL << " " << setw(6) << log10(err) << ",   " ;
//            EXPECT_NEAR(err,0.0,1e-7/Nradial);
        }
    cout << endl;
        
    }
}

TEST_F(MeshIntegralsTests, GObritals)
{
    StreamableObject::SetToPretty();
    InitMolecule();
    cout << *bs << endl;
    auto ibs=bs->beginT();
    SMatrix<double> o= ie->MakeOverlap(*ibs);

    MeshParams mp({qchem::MHL,100,3,2,qchem::GaussLegendre,12,4,3,3});
    MeshIntegrator<double> mi(cl->CreateMesh(mp));
    SMatrix<double> on= mi.Overlap(**ibs);
   // cout << o << on << o-on << endl;
    double err=norm1(o-on);
    cout <<  std::setprecision(12) << err << endl ;
    
}
