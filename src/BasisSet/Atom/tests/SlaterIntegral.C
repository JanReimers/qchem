// File: UnitTests/SlaterIntegral.C 


#include "gtest/gtest.h"
#include "nlohmann/json.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

import qchem.BasisSet.Atom.Factory;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet;
import qchem.Constants;
import qchem.Structure;
import qchem.Structure.MolecularMesh;   // MakeMolecularMesh (qcMesh mesh)
import qchem.Mesh.Quadrature;           // qcMesh quadrature + ScalarField/BasisField
import qchem.VectorFunction;
import qchem.Symmetry.Atom.Spherical;
import qchem.Blaze;
using namespace qchem;

using std::cout;
using std::endl;
using Real_BS  =BasisSet::Real_BS;
using Real_OIBS=BasisSet::Real_OIBS;

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
struct OneOverR  : qcMesh::ScalarField<double>
{
    double  operator()(const rvec3_t& r) const override {double m=norm(r); return m==0.0?0.0:1.0/m;}
    rvec3_t Gradient  (const rvec3_t&)   const override {return rvec3_t(0,0,0);}
};
struct OneOverR2 : qcMesh::ScalarField<double>
{
    double  operator()(const rvec3_t& r) const override {double m=norm(r); return m==0.0?0.0:1.0/(m*m);}
    rvec3_t Gradient  (const rvec3_t&)   const override {return rvec3_t(0,0,0);}
};
qcMesh::Mesh AtomMesh(const Structure& st, int nRadial, int mhl_m, double alpha, int nAngular)
{
    return MakeMolecularMesh(st, {.radial=qcMesh::RadialKind::MHL, .nRadial=nRadial, .mhl_m=mhl_m,
                                  .mhl_alpha=alpha, .angular=qcMesh::AngularKind::Gauss, .nAngular=nAngular});
}
} //anon
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
    , bs(0)
    , cl(new Atom(Z,0.0,Vector3D(0,0,0)))
    {
        nlohmann::json js = {
        {"type",BasisSet::Atom::Type::Slater},
        {"N", 6}, {"emin", 0.1}, {"emax", 10.0},
        };
        bs=BasisSet::Atom::Factory(js,75);
        itsMesh=AtomMesh(*cl,200,3,2.0,1);
    }

    int Lmax, Z;
    Real_BS* bs;
    Structure* cl;
    qcMesh::Mesh itsMesh;
};


TEST_F(SlaterRadialIntegralTests, Overlap)
{
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t S=oi->Overlap();
        for (auto d:blazem::diagonal(S)) EXPECT_NEAR(d,1.0,1e-15);
        rsmat_t Snum = qcMesh::Overlap(itsMesh,BFView(*oi));
        EXPECT_NEAR(blazem::max(blazem::abs(S-Snum)),0.0,1e-8);
    }
}

TEST_F(SlaterRadialIntegralTests, Nuclear)
{
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t Hn=oi->Nuclear(cl);
        rsmat_t Hnnum = -1*qcMesh::WeightedOverlap(itsMesh,BFView(*oi),OneOverR());
        EXPECT_NEAR(blazem::max(blazem::abs(Hn-Hnnum)),0.0,1e-7);

    }
}

TEST_F(SlaterRadialIntegralTests, Kinetic)
{
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t K=oi->Kinetic();   // the <p^2>=<-nabla^2> block (no 1/2); see BasisSet/Orbital_1E_IBS.C
        //cout << S << endl;
        int l=qchem::Symmetry::Atom::Getl(oi->GetSymmetry());;
        // ...which equals Grad2 (radial) + centrifugal l(l+1)<r^-2>, confirming the no-1/2 convention.
        rsmat_t Knum = qcMesh::KineticGrad2(itsMesh,BFView(*oi))
                     + l*(l+1)*qcMesh::WeightedOverlap(itsMesh,BFView(*oi),OneOverR2());
        EXPECT_NEAR(blazem::max(blazem::abs(K-Knum)),0.0,1e-10);
        
        // cout << "K=" << K << endl;
        // cout << "Knum=" << Knum << endl;
    }
}


