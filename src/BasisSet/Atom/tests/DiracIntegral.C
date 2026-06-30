// File: DiracIntegral.C  Test the Dirac basis sets and integral engine.


#include "gtest/gtest.h"
#include "nlohmann/json.hpp"
#include <iostream>
#include <iomanip>

import qchem.BasisSet.Atom.Factory;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet;
import qchem.Math;
import qchem.Structure;
import qchem.Structure.MolecularMesh;   // MakeMolecularMesh (qcMesh mesh)
import qchem.Mesh.Quadrature;           // qcMesh quadrature + ScalarField/BasisField
import qchem.VectorFunction;
import qchem.Streamable;
import qchem.BasisSet.Internal.Orbital_DHF_IBS;
import qchem.Blaze;
using namespace qchem;

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
qcMesh::Mesh AtomMesh(const Structure& st, int nRadial, int mhl_m, double alpha, int nAngular)
{
    return MakeMolecularMesh(st, {.radial=qcMesh::RadialKind::MHL, .nRadial=nRadial, .mhl_m=mhl_m,
                                  .mhl_alpha=alpha, .angular=qcMesh::AngularKind::Gauss, .nAngular=nAngular});
}
} //anon

using BasisSet::Real_BS;
using BasisSet::Real_OIBS;
using Real_IBS=Real_OIBS;
using RKB_OIBS=BasisSet::Orbital_RKB_IBS_Imp<double>;
using RKBL_OIBS=BasisSet::Orbital_RKBL_IBS<double>;
using RKBS_OIBS=BasisSet::Orbital_RKBS_IBS<double>;

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
    : Lmax(0   )
    , Z(1)
    , sbs(0)
    , gbs(0)
    , cl(new Atom(Z,0.0,Vector3D(0,0,0)))
    {
        nlohmann::json js = {{"type",BasisSet::Atom::Type::Slater_RKB}, {"N", 3}, {"emin", 0.1}, {"emax", 10.0} };
        sbs=BasisSet::Atom::Factory(js,2);
        js = {{"type",BasisSet::Atom::Type::Gaussian_RKB}, {"N", 3}, {"emin", 0.1}, {"emax", 10.0} };
        gbs=BasisSet::Atom::Factory(js,2);
        itsMesh=AtomMesh(*cl,200,3,2.0,1);
    }
    
    static const RKBL_OIBS* GetLarge(const Real_IBS* ibs)
    {
        assert(ibs);
        const RKB_OIBS* dirbs=dynamic_cast<const RKB_OIBS*>(ibs);
        assert(dirbs);
        return dirbs->itsRKBL;
    }
    static const RKBS_OIBS* GetSmall(const Real_IBS* ibs)
    {
        assert(ibs);
        const RKB_OIBS* dirbs=dynamic_cast<const RKB_OIBS*>(ibs);
        assert(dirbs);
        return dirbs->itsRKBS;
    }

    
    
    static rsmat_t merge_diag(const rsmat_t& l,const rsmat_t& s)
    {
        return RKB_OIBS::merge_diag(l,s);
    }

    static rsmat_t merge_off_diag(const rmat_t& l)
    {
        return RKB_OIBS::merge_off_diag(l);
    }

    int Lmax, Z;
    Real_BS* sbs;
    Real_BS* gbs;
    Structure* cl;
    qcMesh::Mesh itsMesh;
};

TEST_F(DiracIntegralTests, SlaterBasisSet)
{
    
    cout << *sbs << endl;
}

TEST_F(DiracIntegralTests, GaussianBasisSet)
{
    
    cout << *gbs << endl;
}


TEST_F(DiracIntegralTests, SlaterOverlap)
{
    
    for (auto oi:sbs->Iterate<Real_OIBS >())
    {
        rsmat_t S=oi->Overlap();
        {
            rvec_t d=blazem::diagonal(S);
            for (size_t i=0;i<d.size()/2;i++) EXPECT_NEAR(d[i],1.0,1e-15);
        }
        // cout << std::fixed << std::setprecision(3) << std::setw(6) << S << S1 << endl;
        rsmat_t SLnum = qcMesh::Overlap(itsMesh,BFView(*GetLarge(oi)));
        rsmat_t SSnum = qcMesh::Overlap(itsMesh,BFView(*GetSmall(oi)));
        rsmat_t Snum=merge_diag(SLnum,SSnum);
        // cout << Snum << S << endl;

        EXPECT_NEAR(blazem::max(blazem::abs(S-Snum)),0.0,1e-14);
    }
}

TEST_F(DiracIntegralTests, GaussianOverlap)
{
    
    for (auto oi:gbs->Iterate<Real_OIBS >())
    {
        rsmat_t S=oi->Overlap();
        {
            rvec_t d=blazem::diagonal(S);
            for (size_t i=0;i<d.size()/2;i++) EXPECT_NEAR(d[i],1.0,1e-15);
        }
        // cout << std::fixed << std::setprecision(3) << std::setw(6) << S << S1 << endl;
        rsmat_t SLnum = qcMesh::Overlap(itsMesh,BFView(*GetLarge(oi)));
        rsmat_t SSnum = qcMesh::Overlap(itsMesh,BFView(*GetSmall(oi)));
//        cout << SLnum << SSnum << endl;
        rsmat_t Snum=merge_diag(SLnum,SSnum);
        // cout << Max(fabs(S-Snum)) << endl;
        EXPECT_NEAR(blazem::max(blazem::abs(S-Snum)),0.0,1e-14);
    }
}


TEST_F(DiracIntegralTests, SlaterNuclear)
{
    
    int Z=cl->GetNuclearCharge();
    for (auto oi:sbs->Iterate<Real_OIBS >())
    {
        rsmat_t Ven=oi->Nuclear(cl);
        rsmat_t VenLnum = -Z*qcMesh::WeightedOverlap(itsMesh,BFView(*GetLarge(oi)),OneOverR());
        rsmat_t VenSnum = -Z*qcMesh::WeightedOverlap(itsMesh,BFView(*GetSmall(oi)),OneOverR());
        rsmat_t Vennum=merge_diag(VenLnum,VenSnum);
        //cout << "Ven=" << Ven << endl << "Ven num=" << Vennum << endl;
        //Because of the singularity at the origin, the error is larger than the other integrals.
        EXPECT_NEAR(blazem::max(blazem::abs(Ven-Vennum)),0.0,1e-11);        
    }
}


TEST_F(DiracIntegralTests, GaussianNuclear)
{
    
    int Z=cl->GetNuclearCharge();
    for (auto oi:gbs->Iterate<Real_OIBS >())
    {
        rsmat_t Ven=oi->Nuclear(cl);
        rsmat_t VenLnum = -Z*qcMesh::WeightedOverlap(itsMesh,BFView(*GetLarge(oi)),OneOverR());
        rsmat_t VenSnum = -Z*qcMesh::WeightedOverlap(itsMesh,BFView(*GetSmall(oi)),OneOverR());
        rsmat_t Vennum=merge_diag(VenLnum,VenSnum);
        //cout << "Ven=" << Ven << endl << "Ven num=" << Vennum << endl;
        // cout << "Ven=" << Ven << endl << "Ven1=" << Ven1 << endl;
        //Because of the singularity at the origin, the error is larger than the other integrals.
        EXPECT_NEAR(blazem::max(blazem::abs(Ven-Vennum)),0.0,1e-11);        
    }
}


// NOTE: the Dirac relativistic-kinetic numerical cross-check (SlaterKinetic/GaussianKinetic) used
// MeshIntegrator::Grada_b -- the l=0-only <grad a|b> gradient cross-term that qcMesh deliberately
// dropped.  The analytic Kinetic() it cross-checked is validated end-to-end by A_DHF (UTMain).

 //  For symmetric matricies, Sum() should know to double all the off diagonal elements.
TEST_F(DiracIntegralTests, Contraction)
{
    size_t N=5;
    rsmat_t H(N),D(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
        {
            H(i,j)=std::rand();
            D(i,j)=std::rand();
        }
    rmat_t Hf(H),Df(D);
    double c =blazem::sum(H %D );
    double cf=blazem::sum(Hf%Df);
    EXPECT_NEAR(c,cf,5e-15);
}

