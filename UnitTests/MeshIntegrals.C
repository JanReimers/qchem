// File MeshIntegrals.C  Run through varios mesh types and parameters

#include "gtest/gtest.h"
#include <iomanip>
#include <blaze/Math.h>
import qchem.LAParams;
import qchem.BasisSet;
import qchem.BasisSet.Molecule.PolarizedGaussian;
import qchem.Mesh.Integrator;
import qchem.Cluster;
import qchem.Symmetry;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Readers.Gaussian94;

using namespace BasisSet::Molecule;
using BasisSet::Real_OIBS;

using std::cout;
using std::endl;
using std::setw;

// double norm1(const SMatrix<double>& m)
// {
//     return sqrt(Sum(DirectMultiply(m,m)));
// }
//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class MeshIntegralsTests : public ::testing::Test
{
public:
    MeshIntegralsTests()
    : Z(5)
    , reader("../../../BasisSetData/dzvp.bsd")
    , bs(0)
    , cl(new Atom(Z,0.0,Vector3D(0,0,0)))
    {
        
        
    }
    
    void InitAtom()
    {
        bs=new PolarizedGaussian::BasisSet(&reader,cl);
    }
    
    void InitMolecule()
    {
        Molecule* m=new Molecule();
        m->Insert(new Atom(Z,0.0,Vector3D( 1.,0.,0.)));
        m->Insert(new Atom(Z,0.0,Vector3D(-1.,0.,0.)));
        cl=m;
        bs=new PolarizedGaussian::BasisSet(&reader,cl);        
    }
    
    int Z;
    PolarizedGaussian::Gaussian94Reader reader;
    PolarizedGaussian::BasisSet* bs;
    Cluster* cl;
};

TEST_F(MeshIntegralsTests, PolGaussianOverlap)
{
    
    InitAtom();
    
    //cout.precision(2);
    size_t Nradial=100;
    double af=0.1;
    cout << "alpha1     alpha2       m      norm     max(abs)"  << endl;
    for (size_t mMHL=1;mMHL<=5;mMHL++)
    {
        int a1=0,a2=0;
        double min1=100.0,min2=100.0;
        for (int aMHL=1;aMHL<=50;aMHL++)
        {
            double alpha=af*aMHL;
            MeshParams mp({qchem::MHL,Nradial,mMHL,alpha,qchem::Gauss,12,0,0});
            MeshIntegrator<double> mi(cl->CreateMesh(mp));
            for (auto ibs:bs->Iterate<Real_OIBS>())
            {
                // cout << *ibs << endl;
                rsmat_t delta= ibs->Overlap()-mi.Overlap(*ibs);
                double err=norm(delta);
                double merr=max(abs(delta));
                if (err<min1)
                {
                    min1=err;
                    a1=aMHL;
                }
                if (merr<min2)
                {
                    min2=merr;
                    a2=aMHL;
                }
    //            EXPECT_NEAR(err,0.0,1e-7/Nradial);
            }
            
        }
        cout << std::fixed << setw(3) << af*a1 << " " << af*a2 << " " << setw(2) << mMHL << " " << setw(6) << log10(min1) << " "<< setw(6) << log10(min2);
        cout << endl;
    }
}

