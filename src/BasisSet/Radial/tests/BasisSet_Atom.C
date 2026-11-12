// File: UnitTests/BasisSet_Atom.C  Unit test the Atom IBS Evaluators
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
using std::cout;
using std::endl;

import qchem.BasisSet.Radial.Factory;
import qchem.BasisSet.Radial.Evaluators.Slater.IBS;
import qchem.BasisSet.Radial.Evaluators.Gaussian.IBS; 
import qchem.BasisSet.Radial.Evaluators.BSpline.IBS;
import qchem.BasisSet.Radial.Evaluators;
import qchem.BasisSet.Radial.IBS;          // Radial::RadialAngularID (the atom radial|angular identity face)
import qchem.BasisSet.Internal.Orbital_ERI4_IBS;   // tests may cheat: the ERI4 substrate face
import qchem.BasisSet;
import qchem.BasisSet.Orbital_DHF_IBS;     // Real_ORKB: the Dirac blocks (the hand-gradient audit)

import qchem.Structure;
import qchem.Mesh.Quadrature;           // qcMesh quadrature + ScalarField/BasisField
import qchem.VectorFunction;
import qchem.Symmetry.Factory;
import qchem.ElectronConfiguration.AtomNR;
import qchem.Blaze;
import qchem.Reporting;   // the announced basis.exponents rows (the heat-map gate)
using namespace qchem;

namespace
{
struct OneOverR  : ScalarFunction<double>
{
    double  operator()(const rvec3_t& r) const override {double m=norm(r); return m==0.0?0.0:1.0/m;}
    rvec3_t Gradient  (const rvec3_t&)   const override {return rvec3_t(0,0,0);}
};
struct OneOverR2 : ScalarFunction<double>
{
    double  operator()(const rvec3_t& r) const override {double m=norm(r); return m==0.0?0.0:1.0/(m*m);}
    rvec3_t Gradient  (const rvec3_t&)   const override {return rvec3_t(0,0,0);}
};
qcMesh::Mesh AtomMesh(const Structure& st, int nRadial, int mhl_m, double alpha, int angularDegree)
{
    return st.CreateIntegrationMesh({.radial=qcMesh::RadialKind::MHL, .nRadial=nRadial, .mhl_m=mhl_m,
                                  .mhl_alpha=alpha, .angular=qcMesh::AngularKind::Lebedev, .angularDegree=angularDegree});
}
} //anon

using rohfbs_t =BasisSet::Orbital_ERI4_IBS<double>;   // ERI4 substrate: these tests build blocks
using BasisSet::Real_BS;
using BasisSet::Real_OIBS;
using namespace qchem::BasisSet::Radial;
using namespace Evaluators;
using enum BasisSetAccuracy;

// The tiny fixed test bases that BasisSetAccuracy::N3/N5 used to build (N3 = 3 even-tempered exponents
// {0.5,1,2} / BSpline N=3; N5 = 5 exponents {0.25..4} / BSpline N=5).  N3/N5 were REMOVED from the production
// BasisSetAccuracy pool (ill-conditioned, invalid for SCF), but these low-level construction/ID/cache tests
// legitimately want a small deterministic basis, so build it here from the equivalent JSON (test-local).
static Real_BS* N3Basis(Type t, size_t Z)
{
    nlohmann::json js = (t==Type::BSpline6||t==Type::BSpliner6)
        ? nlohmann::json{{"type",t},{"rmin",0.5},{"rmax",2},{"N",3}}
        : nlohmann::json{{"type",t},{"ltrim",0},{"exponents",{0.5,1,2.0}}};
    return Factory(js, Z);
}
static Real_BS* N5Basis(Type t, size_t Z)
{
    nlohmann::json js = (t==Type::BSpline6||t==Type::BSpliner6)
        ? nlohmann::json{{"type",t},{"rmin",0.25},{"rmax",4},{"N",5}}
        : nlohmann::json{{"type",t},{"ltrim",0},{"exponents",{0.25,0.5,1,2.0,4.0}}};
    return Factory(js, Z);
}

//----------------------------------------------------------------------------------------
//
//  Testing common to all atom basis set evaluators
//
template <is1E_HF_Evaluator E> class BasisSet_Common : public ::testing::Test
{
public:
    BasisSet_Common() 
        : es{0.5,1.0,2.0}
        , N(es.size())
        , LMax(3)
        , cl(new Atom(1,0.0,Vector3D(0,0,0)))
        , bs(0)
    {
        
        itsMesh=AtomMesh(*cl,500,3,2.0, 0);
    }
    ~BasisSet_Common()
    {
        for (auto ev:evals) delete ev;
        delete cl;
        delete bs;
    }
    void Insert(E* eval)
    {
        evals.push_back(eval);
    }

    void TestOverlap(double eps) const;
    void TestGrad2  (double eps) const;
    void TestInv_r1 (double eps) const;
    void TestInv_r2 (double eps) const;
    // void TestCharge (double eps) const;


    rvec_t es;
    size_t N,LMax;
    std::vector<E*> evals;
    Structure* cl;
    qcMesh::Mesh itsMesh;
    Real_BS* bs;
};
template <is1E_HF_Evaluator E> void BasisSet_Common<E>::TestOverlap(double eps) const
{
    EXPECT_GT(evals.size(),0);
    size_t index=0;
    for (auto ev:evals)
    {
        size_t N=ev->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= ev->Overlap(i,j);
        for (auto d:blazem::diagonal(S)) EXPECT_NEAR(d,1.0,1e-15);
        rsmat_t Snum = qcMesh::Overlap(itsMesh,*ev);
        cout.precision(2);
        // cout << S-Snum << endl;
        cout << "l=" << index++ << std::endl;
        EXPECT_NEAR(blazem::max(blazem::abs(S-Snum)),0.0,eps);
    }
}
template <is1E_HF_Evaluator E> void BasisSet_Common<E>::TestGrad2  (double eps) const
{
    for (auto ev:evals)
    {
        size_t N=ev->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= ev->Grad2(i,j);
        rsmat_t Snum = qcMesh::KineticGrad2(itsMesh,*ev);
        EXPECT_NEAR(blazem::max(blazem::abs(S-Snum)),0.0,eps);
    }
        
}
template <is1E_HF_Evaluator E> void BasisSet_Common<E>::TestInv_r1 (double eps) const
{
    for (auto ev:evals)
    {
        size_t N=ev->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= ev->Inv_r1(i,j);
        rsmat_t Snum = qcMesh::MatrixOverlap(itsMesh,*ev,OneOverR());
        EXPECT_NEAR(blazem::max(blazem::abs(S-Snum)),0.0,eps);
    }
        
}
template <is1E_HF_Evaluator E> void BasisSet_Common<E>::TestInv_r2 (double eps) const
{
    for (auto ev:evals)
    {
        size_t N=ev->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= ev->Inv_r2(i,j);
        rsmat_t Snum = qcMesh::MatrixOverlap(itsMesh,*ev,OneOverR2());
        EXPECT_NEAR(blazem::max(blazem::abs(S-Snum)),0.0,eps);
    }
        
}
// template <isFit_Evaluator E> void BasisSet_Common<E>::TestCharge (double eps) const
// {
//     for (auto ev:evals)
//     {
//         size_t N=ev->size();
//         rvec_t S(N);
//         for (auto i:iv_t(0,N))
//                 S[i]= ev->Charge(i);
//         rvec_t Snum = mintegrator->Integrate(*ev);
//         EXPECT_NEAR(max(abs(S-Snum)),0.0,eps);
//     }
        
// }

//----------------------------------------------------------------------------------------
//
//  Testing atom Slater basis set evaluators
//
class BasisSet_SL: public BasisSet_Common<Slater::NR_Evaluator>
{
public:

    BasisSet_SL() : BasisSet_Common<Slater::NR_Evaluator>()
    {
        Atom_EC ec(86); //Radon has f orbtials with no magnetic splitting.
        for (auto ir:ec.GetIrreps())
            Insert(new Slater::NR_Evaluator(es,ir)); 
        bs=N3Basis(Type::Slater,86);
        cout << es << endl << *bs << endl;
    }
    
    static double R0(double a, double b, int la, int lb);

   
};

double BasisSet_SL::R0(double a, double b, int la, int lb) 
{
    double f=FourPi2/(a*b*(a+b));
    if (la==0 && lb==0)
        return 2*f*( 
                    1/(pow(a,1)*pow(b,1)*pow(a+b,0))
                    +1/pow(a+b,2)
                    );
    if (la==1 && lb==1)
         return 144*f*( 
           1/(pow(a,2)*pow(b,2)*pow(a+b,2))
         + 2/(pow(a,1)*pow(b,1)*pow(a+b,4))
         + 5/(pow(a,0)*pow(b,0)*pow(a+b,6))
         + 1/(pow(a,3)*pow(b,3)*pow(a+b,0))
          );
          
    if (la==2 && lb==2)
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

TEST_F(BasisSet_SL,Overlap) {TestOverlap(1e-15);}
TEST_F(BasisSet_SL,Grad2  ) {TestGrad2  (3e-15);}
TEST_F(BasisSet_SL,Inv_r1 ) {TestInv_r1 (3e-14);}
TEST_F(BasisSet_SL,Inv_r2 ) {TestInv_r2 (4e-9);}
// TEST_F(BasisSet_SL,Charge ) {TestCharge (9e-14);}

TEST_F(BasisSet_SL,AnalyticOverlap)
{
    for (auto ev:evals)
    {
        int l=ev->Getl();
        size_t N=ev->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= ev->Overlap(i,j);
        // cout << S << endl;
        for (auto i:iv_t(0,S.rows()))
            for (auto j:iv_t(i,S.rows()))
            {
                double a=es[i],b=es[j];
                EXPECT_NEAR(S(i,j),pow(2.0/(sqrt(a/b)+sqrt(b/a)),2*l+3),1e-15);
            }
    }
        
}
TEST_F(BasisSet_SL,AnalyticRepulsion)
{
    for (auto ev:evals)
    {
        int l=ev->Getl();
        if (l>=3) continue;
        size_t N=ev->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= ev->Repulsion(i,j);
        // cout << "l=" << l << " S=" << S << endl;
        rvec_t   ns=ev->Norm();
        for (auto i:iv_t(0,S.rows()))
            for (auto j:iv_t(i,S.rows()))
            {
                double a=es[i],b=es[j];
                double rerr=(S(i,j) - R0(a,b,l,l)*ns[i]*ns[j])/S(i,j);
                EXPECT_NEAR(rerr,0.0,3e-14);
            }
    }
        
}

TEST_F(BasisSet_SL,HF_ERIs)
{
    auto a=evals.begin();
    for (auto aibs:bs->Iterate<rohfbs_t>())
    {
        auto c=evals.begin();
        for (auto cibs:bs->Iterate<rohfbs_t>())
        {
            if (aibs->BasisSetID() <= cibs->BasisSetID())   // ERI4 cache is canonical-only (doc/ERI4Rework.md §5.2)
            {
                ERI4 J2=aibs->Direct(*cibs);
                ERI4 K2=aibs->Exchange(*cibs);
            }
            ++c;
        }
        ++a;
    }
}

std::string angularIDs[]={"0 {}","1 {}","2 {}","3 {}"};
TEST_F(BasisSet_SL,IDs)
{
    size_t index=0;
    for (auto ibs:bs->Iterate<RadialAngularID>())
    {
        EXPECT_EQ(ibs->Name(),"SL");
        EXPECT_EQ(ibs->RadialID(),"SL N=3 {0.5 1 2 }");   // EXACT since 2026-11-12: every exponent, %.17g (a cache key)
        EXPECT_EQ(ibs->AngularID(),angularIDs[index++]);
    }
}
//----------------------------------------------------------------------------------------
//
//  Testing atom Gaussian basis set evaluators
//

class BasisSet_SG: public BasisSet_Common<Gaussian::NR_Evaluator>
{
public:

    BasisSet_SG() : BasisSet_Common()
    {
        Atom_EC ec(86); //Radon has f orbtials with no magnetic splitting.
        for (auto ir:ec.GetIrreps())
            Insert(new Gaussian::NR_Evaluator(es,ir)); 
        bs=N3Basis(Type::Gaussian,86);
    }
    static double R0(double a, double b, int la, int lb);
};

double BasisSet_SG::R0(double a, double b, int la, int lb) 
{
    double f=FourPi2*Pi12/8*(1/(a*b*sqrt(a+b)));
    if (la==0 && lb==0)
        return f;
    if (la==1 && lb==1)
         return f*3/4*(2/(a*b)+1/((a+b)*(a+b)));
          
    if (la==2 && lb==2)
        return f*15/(4*pow(a+b,4))*(2*(a+b)*(b/(a*a)+a/(b*b)) + 7*(b/a+a/b)+63./4.);
    
    
                    
    assert(false);
    return 0;

}


TEST_F(BasisSet_SG,Overlap) {TestOverlap(1e-15);}
TEST_F(BasisSet_SG,Grad2  ) {TestGrad2  (6e-15);}
TEST_F(BasisSet_SG,Inv_r1 ) {TestInv_r1 (3e-14);}
TEST_F(BasisSet_SG,Inv_r2 ) {TestInv_r2 (4e-9);}
// TEST_F(BasisSet_SG,Charge ) {TestCharge (5e-14);}

TEST_F(BasisSet_SG,AnalyticOverlap)
{
    for (auto ev:evals)
    {
        int l=ev->Getl();
        size_t N=ev->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= ev->Overlap(i,j);

         for (auto i:iv_t(0,S.rows()))
            for (auto j:iv_t(i,S.rows()))
            {
                double a=es[i],b=es[j];
                EXPECT_NEAR(S(i,j),pow(2.0/(sqrt(a/b)+sqrt(b/a)),(2.0*l+3.0)/2.0),1e-15);
            }
    }
        
}
TEST_F(BasisSet_SG,AnalyticRepulsion)
{
    for (auto ev:evals)
    {
        int l=ev->Getl();
        if (l>=3) continue;
        size_t N=ev->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= ev->Repulsion(i,j);
        // cout << "l=" << l << " S=" << S << endl;
        rvec_t   ns=ev->Norm();
         for (auto i:iv_t(0,S.rows()))
            for (auto j:iv_t(i,S.rows()))
            {
                double a=es[i],b=es[j];
                double rerr=(S(i,j) - R0(a,b,l,l)*ns[i]*ns[j])/S(i,j);
                EXPECT_NEAR(rerr,0.0,3e-14);
            }
    }
        
}
TEST_F(BasisSet_SG,HF_ERIs)
{
    auto a=evals.begin();
    for (auto aibs:bs->Iterate<rohfbs_t>())
    {
        auto c=evals.begin();
        for (auto cibs:bs->Iterate<rohfbs_t>())
        {
            if (aibs->BasisSetID() <= cibs->BasisSetID())   // ERI4 cache is canonical-only (doc/ERI4Rework.md §5.2)
            {
                ERI4 J2=aibs->Direct(*cibs);
                ERI4 K2=aibs->Exchange(*cibs);
            }
            ++c;
        }
        ++a;
    }
}

TEST_F(BasisSet_SG,IDs)
{
    size_t index=0;
    for (auto ibs:bs->Iterate<RadialAngularID>())
    {
        EXPECT_EQ(ibs->Name(),"SG");
        EXPECT_EQ(ibs->RadialID(),"SG N=3 {0.5 1 2 }");   // EXACT since 2026-11-12: every exponent, %.17g (a cache key)
        EXPECT_EQ(ibs->AngularID(),angularIDs[index++]);
    }
}

//----------------------------------------------------------------------------------------
//
//  Testing atom BSpline basis set evaluators
//

class BasisSet_BS: public BasisSet_Common<BSpline::Evaluator<6>>
{
public:

    BasisSet_BS() : BasisSet_Common()
    {
        for (size_t l=0;l<=3;l++)
            Insert(new BSpline::Evaluator<6>(5,0.01,20.0,Symmetry::YFactory(l)));

        bs=N5Basis(Type::BSpline6,86);
    }
   
};

TEST_F(BasisSet_BS,Overlap) {TestOverlap(3e-11);}
TEST_F(BasisSet_BS,Grad2  ) {TestGrad2  (1e-4);}
TEST_F(BasisSet_BS,Inv_r1 ) {TestInv_r1 (5e-9);}
TEST_F(BasisSet_BS,Inv_r2 ) {TestInv_r2 (5e-3);}
// TEST_F(BasisSet_BS,Charge ) {TestCharge (8e-7);}

TEST_F(BasisSet_BS,IDs)
{
    size_t index=0;
    for (auto ibs:bs->Iterate<RadialAngularID>())
    {
        // cout << "index=" << index << " ibs=" << *ibs << endl;
        EXPECT_EQ(ibs->Name(),"BSpline<6>");
        EXPECT_EQ(ibs->RadialID(),"BSpline<6> grid: N=6 {0,0.25,0.5,1,2,4,}");   // EXACT since 2026-11-12: the whole knot vector (a cache key)
        EXPECT_EQ(ibs->AngularID(),angularIDs[index++]);
    }
}

//----------------------------------------------------------------------------------------
//
//  THE HAND-GRADIENT AUDIT (lane B 2026-11-09, ahead of the user's Enzyme decision): every radial family's hand Gradient(r), through
//  the blocks the Factory builds, against the central difference of the same block's operator()(r), every function and axis.  Slater,
//  Gaussian, BSpline<6> (r^l and plain), and the two Dirac (RKB) families, whose blocks stack the large and small components.  Points
//  0.3 to 2.1 bohr off the nucleus, off every BSpline knot.  Step h = 1e-5: the central difference's h²/6 |φ‴| is ≲ 1e-10 |φ‴| for the
//  exponents here (≤ 10), so 1e-7 (relative above 1).  The largest derivative probed must be O(1e-1) or more, so a sign would show.
//  (The kinetic Grad2 gates above use Gradient quadratically, so they cannot see a sign, and never reached the RKB blocks.)
//
TEST(RadialHandGradient, EveryRadialFamilysGradientIsTheFiniteDifferenceOfItsValues)
{
    using enum Type;
    const std::vector<std::pair<Type,nlohmann::json>> cases{
        {Slater,    {{"type",Slater},   {"ltrim",0},{"exponents",{0.5,1.0,2.0,4.0}}}},
        {Gaussian,  {{"type",Gaussian}, {"ltrim",0},{"exponents",{0.5,1.0,2.0,4.0}}}},
        {BSpline6,  {{"type",BSpline6}, {"rmin",0.25},{"rmax",4},{"N",5}}},
        {BSpliner6, {{"type",BSpliner6},{"rmin",0.25},{"rmax",4},{"N",5}}},
        {Gaussian_RKB, {{"type",Gaussian_RKB},{"N",4},{"emin",0.5},{"emax",8}}},
        {Slater_RKB,   {{"type",Slater_RKB},  {"N",4},{"emin",0.5},{"emax",8}}}};
    const std::vector<rvec3_t> pts{rvec3_t(0.21,-0.13,0.17), rvec3_t(0.9,0.55,-0.7), rvec3_t(-1.3,1.1,1.2)};
    const double h=1e-5;
    for (const auto& [t,js] : cases)
    {
        std::unique_ptr<Real_BS> bs(Factory(js, 30));   // Zn: s, p and d blocks
        const bool rkb = t==Gaussian_RKB || t==Slater_RKB;
        std::vector<const VectorFunction<double>*> blocks;
        if (rkb) for (auto b : bs->Iterate<BasisSet::Real_ORKB>()) blocks.push_back(b);
        else     for (auto b : bs->Iterate<Real_OIBS>())           blocks.push_back(b);
        ASSERT_GE(blocks.size(), 3u) << js.dump();
        double largest=0.0;
        for (size_t ib=0;ib<blocks.size();ib++)
        {
            const VectorFunction<double>& f=*blocks[ib];
            for (const rvec3_t& r : pts)
            {
                const rvec3vec_t g=f.Gradient(r);
                const rvec_t v=f(r);
                ASSERT_EQ(g.size(), v.size()) << js.dump() << " block " << ib << ": one gradient per function";
                for (int k=0;k<3;k++)
                {
                    rvec3_t p=r, m=r; (k==0 ? p.x : k==1 ? p.y : p.z)+=h; (k==0 ? m.x : k==1 ? m.y : m.z)-=h;
                    const rvec_t vp=f(p), vm=f(m);
                    for (size_t i=0;i<v.size();i++)
                    {
                        const double fd=(vp[i]-vm[i])/(2*h), an=(k==0 ? g[i].x : k==1 ? g[i].y : g[i].z);
                        EXPECT_NEAR(an, fd, 1e-7*std::max(1.0,std::fabs(fd))) << js.dump() << " block " << ib << " function " << i << " axis " << k << " at " << r;
                        largest=std::max(largest, std::fabs(fd));
                    }
                }
            }
        }
        EXPECT_GT(largest, 0.1) << js.dump() << ": the derivatives probed are not trivially small";
    }
}

// OpenWork §3 "valgen's usage heat map shows Gaussian exponents DOUBLED" (2026-11-02; closed 2026-11-11): the heat map labels its rows
// from basis.exponents, which each exponential shell ANNOUNCES as it is built.  The announced values must be the exponents the caller
// stated (alpha in exp(-alpha r^2) / exp(-zeta r)), for every l of a Gaussian and a Slater basis -- so the heat map's labels are the .bsd's.
TEST(BasisSet_Atom, EachShellAnnouncesTheExponentsItWasGiven)
{
    const std::vector<double> given={0.25,0.5,1,2.0,4.0};
    for (Type t : {Type::Gaussian, Type::Slater})
    {
        qchem::report::Begin("announce test");
        {
            qchem::report::Scope basis("basis");
            std::unique_ptr<Real_BS> bs(Factory(nlohmann::json{{"type",t},{"ltrim",0},{"exponents",given}}, size_t{18}));   // Ar: s and p shells
            const auto& rows=qchem::report::CurrentRunReport().at("basis").at("exponents");
            size_t n=0;
            for (const auto& r : rows)
            {
                EXPECT_EQ(r.at("values").get<std::vector<double>>(), given) << r.at("irrep").get<std::string>() << " type " << int(t);
                ++n;
            }
            EXPECT_GE(n, 2u) << "an s and a p shell announced";
        }
        qchem::report::End();
    }
}
