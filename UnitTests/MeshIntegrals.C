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
    : Lmax(2    )
    , Z(1)
    , lap({qchem::Lapack,qchem::SVD,1e-6,1e-12})
    , ie(new PolarizedGaussian::IntegralEngine())
    , bs(new PolarizedGaussian::BasisSet(lap,3,.1,10.0,Lmax))
    , cl(new Molecule())
    {
        cl->Insert(new Atom(Z,0.0,Vector3D(0,0,0)));
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

