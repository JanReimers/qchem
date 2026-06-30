// File MeshIntegrals.C  Run through varios mesh types and parameters

#include "gtest/gtest.h"
#include <iomanip>
#include <filesystem>
#ifndef BASISSET_DATA_PATH
#error "BASISSET_DATA_PATH must be defined by CMake"
#endif

import qchem.BasisSet.Molecule.Readers.Gaussian94;
import qchem.BasisSet.Molecule.PG_Cart;
import qchem.BasisSet;
import qchem.Structure;
import qchem.Structure.MolecularMesh;   // MakeMolecularMesh (qcMesh mesh)
import qchem.Mesh.Quadrature;           // qcMesh::Overlap + BasisField
import qchem.VectorFunction;
import qchem.Symmetry;
import qchem.Blaze;
using namespace qchem;

using namespace qchem::BasisSet::Molecule;
using BasisSet::Real_OIBS;

namespace
{
class BFView : public qcMesh::BasisField<double>
{
    const VectorFunction<double>& its;
public:
    explicit BFView(const VectorFunction<double>& v) : its(v) {}
    size_t     size()                       const override {return its.GetVectorSize();}
    rvec_t     operator()(const rvec3_t& r) const override {return its(r);}
    rvec3vec_t Gradient  (const rvec3_t& r) const override {return its.Gradient(r);}
};
} //anon

using std::cout;
using std::endl;
using std::setw;
static const std::filesystem::path basisset_data_dir = BASISSET_DATA_PATH;

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
    , reader(basisset_data_dir / "dzvp.bsd")
    , bs(0)
    , cl(new Atom(Z,0.0,Vector3D(0,0,0)))
    {
        
        
    }
    
    void InitAtom()
    {
        bs=new PG_Cart::BasisSet(&reader,cl);
    }
    
    void InitMolecule()
    {
        Molecule* m=new Molecule();
        m->Insert(new Atom(Z,0.0,Vector3D( 1.,0.,0.)));
        m->Insert(new Atom(Z,0.0,Vector3D(-1.,0.,0.)));
        cl=m;
        bs=new PG_Cart::BasisSet(&reader,cl);        
    }
    
    int Z;
    ::qchem::BasisSet::Molecule::Gaussian94Reader reader;
    PG_Cart::BasisSet* bs;
    Structure* cl;
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
            qcMesh::Mesh mesh = MakeMolecularMesh(*cl,
                {.radial=qcMesh::RadialKind::MHL, .nRadial=int(Nradial), .mhl_m=int(mMHL),
                 .mhl_alpha=alpha, .angular=qcMesh::AngularKind::Gauss, .nAngular=12});
            for (auto ibs:bs->Iterate<Real_OIBS>())
            {
                // cout << *ibs << endl;
                rsmat_t delta= ibs->Overlap()-qcMesh::Overlap(mesh,BFView(*ibs));
                double err=blazem::norm(delta);
                double merr=blazem::max(blazem::abs(delta));
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

