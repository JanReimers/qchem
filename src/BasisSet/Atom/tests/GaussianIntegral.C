// File: ERIList.C  Test the DFT persistance classes
#include <blaze/math/expressions/DMatDMatEqualExpr.h>
#include "gtest/gtest.h"
#include "nlohmann/json.hpp"
#include <iostream>
#include <fstream>

import qchem.BasisSet.Atom.Factory;
import qchem.BasisSet.Atom.Evaluators.Gaussian.IBS;
import qchem.BasisSet.Atom.Evaluators.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators;
import qchem.BasisSet;

import qchem.Constants;
import qchem.Structure;
import qchem.Structure.MolecularMesh;   // MakeMolecularMesh (qcMesh mesh)
import qchem.Mesh.Quadrature;           // qcMesh::Overlap/WeightedOverlap/KineticGrad2 + views
import qchem.VectorFunction;             // VectorFunction<double> (the basis evaluator interface)
import qchem.Symmetry.Spherical;
import qchem.Symmetry.Factory;
import qchem.Blaze;

using std::cout;
using std::endl;
using Real_BS=BasisSet::Real_BS;
using Real_OIBS=BasisSet::Real_OIBS;

namespace
{
// Adapters to drive qcMesh's free-function quadrature from the old VectorFunction basis evaluators.
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
// The atom integration mesh: MHL radial x Gauss angular, single-center (natom=1, no Becke).
qcMesh::Mesh AtomMesh(const Structure& cl, int nRadial, int mhl_m, double alpha, int nAngular)
{
    return MakeMolecularMesh(cl, {.radial=qcMesh::RadialKind::MHL, .nRadial=nRadial, .mhl_m=mhl_m,
                                  .mhl_alpha=alpha, .angular=qcMesh::AngularKind::Gauss, .nAngular=nAngular});
}
} //anon

//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class GaussianRadialIntegralTests : public ::testing::Test
{
public:
    GaussianRadialIntegralTests()
    : Lmax(4    )
    , Z(1)
    , bs(0)
    , cl(new Atom(Z,0.0,Vector3D(0,0,0)))
    {
        nlohmann::json js = {
        {"type",BasisSet::Atom::Type::Gaussian},
        {"N", 5}, {"emin", 0.01}, {"emax", 100.0},
        };
        bs=BasisSet::Atom::Factory(js,75);
        itsMesh=AtomMesh(*cl,200,3,2.0,1);
    }

    int Lmax, Z;
    Real_BS* bs;
    Structure* cl;
    qcMesh::Mesh itsMesh;
};




TEST_F(GaussianRadialIntegralTests, Overlap)
{
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t S=oi->Overlap();

        for (auto d:blazem::diagonal(S)) EXPECT_NEAR(d,1.0,1e-15);
        //cout << S << endl;
        rsmat_t Snum = qcMesh::Overlap(itsMesh,BFView(*oi));
        EXPECT_NEAR(blazem::max(blazem::abs(S-Snum)),0.0,1e-8);
       
    }
}

TEST_F(GaussianRadialIntegralTests, Nuclear)
{
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t Hn=oi->Nuclear(cl);
        //cout << S << endl;
        rsmat_t Hnnum = -1*qcMesh::WeightedOverlap(itsMesh,BFView(*oi),OneOverR());
        EXPECT_NEAR(blazem::max(blazem::abs(Hn-Hnnum)),0.0,1e-8);

    }
}

TEST_F(GaussianRadialIntegralTests, Kinetic)
{
    
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t K=oi->Kinetic();   // the <p^2>=<-nabla^2> block (no 1/2); see BasisSet/Orbital_1E_IBS.C
        //cout << S << endl;
        int l=Getl(oi->GetSymmetry());;
        // ...which equals Grad2 (radial) + centrifugal l(l+1)<r^-2>, confirming the no-1/2 convention.
        rsmat_t Knum = qcMesh::KineticGrad2(itsMesh,BFView(*oi))
                     + l*(l+1)*qcMesh::WeightedOverlap(itsMesh,BFView(*oi),OneOverR2());
        EXPECT_NEAR(blazem::max(blazem::abs(K-Knum)),0.0,1e-12);
        
    }
}

TEST_F(GaussianRadialIntegralTests,RkSymmetry_l0)
{
    using namespace BasisSet::Atom::Evaluators::Gaussian;
    typedef rvec11_t rvec11_t; 
    auto eval=new NR_Evaluator(15,.03,20.0,Symmetry::YFactory(0));
    auto cache4=eval->MakeCache4();
    cache4->Register(eval);
    auto ns=eval->Norm();

    rvec11_t Ak({1,1,1,1,1,1,1,1,1,1,1});

    {
        {
            auto rk0001=dynamic_cast<const Rk*>(cache4->Create(0,0,0,1));
            auto rk0100=dynamic_cast<const Rk*>(cache4->Create(0,1,0,0));
            auto rk0010=dynamic_cast<const Rk*>(cache4->Create(0,0,1,0));
            cout << " 0001=" << std::setprecision(12) << rk0001->DirectR0  (0,0   )*ns[0]*ns[0]*ns[0]*ns[1]*FourPi2 << endl;
            cout << " 0100=" << std::setprecision(12) << rk0100->DirectR0  (0,0   )*ns[0]*ns[0]*ns[0]*ns[1]*FourPi2 << endl;
            cout << " 0010=" << std::setprecision(12) << rk0010->DirectR0  (0,0   )*ns[0]*ns[0]*ns[1]*ns[0]*FourPi2 << endl;
            cout << "J0001=" << std::setprecision(12) << eval->direct(rk0001,0,0,Ak)*ns[0]*ns[0]*ns[0]*ns[1]*FourPi2 << endl;
            cout << "J0100=" << std::setprecision(12) << eval->direct(rk0100,0,0,Ak)*ns[0]*ns[1]*ns[0]*ns[0]*FourPi2 << endl;
            cout << "J0010=" << std::setprecision(12) << eval->direct(rk0010,0,0,Ak)*ns[0]*ns[0]*ns[1]*ns[0]*FourPi2 << endl;

            EXPECT_NEAR(rk0001->DirectR0  (0,0   ),rk0100->DirectR0  (0,0   ),1e-15);
            EXPECT_NEAR(rk0001->DirectR0  (0,0   ),rk0100->DirectR0  (0,0   ),1e-15);
            EXPECT_NEAR(rk0001->DirectRk  (0,0,Ak),rk0100->DirectRk  (0,0,Ak),1e-15);
            EXPECT_NEAR(rk0001->DirectRk  (0,0,Ak),rk0100->DirectRk  (0,0,Ak),1e-15);
            EXPECT_NEAR(rk0001->ExchangeRk(0,0,Ak),rk0100->ExchangeRk(0,0,Ak),1e-15);
            EXPECT_NEAR(rk0001->ExchangeRk(0,0,Ak),rk0100->ExchangeRk(0,0,Ak),1e-15);
        }
    }
}



