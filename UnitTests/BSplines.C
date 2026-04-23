// File: BSplines.C  Test the BSplinebasis package included as a submodule


#include "gtest/gtest.h"
#include "nlohmann/json.hpp"
#include <bspline/Core.h>
#include <iostream>
#include <iomanip>
#include <blaze/Math.h>

import qchem.Basisset.Atom.BSpline.GLQuadrature;
import qchem.Symmetry.Angular;
import BasisSet.Atom.BSpline.NR.IBS_Evaluator;


import qchem.Factory;
import qchem.BasisSet;
import qchem.IrrepBasisSet;
import qchem.Mesh.Integrator;
import qchem.Cluster;
import qchem.Symmetry;
import qchem.stl_io;
import qchem.Streamable;
using std::cout;
using std::endl;

class BSplineTests : public ::testing::Test
{
public:
    static constexpr size_t K=9; //Spline order.
    typedef bspline::Spline<double, K> spline_t;
    typedef Real_OIBS ibs_t;
    BSplineTests() 
        : LMax(4)
        , cl(new Atom(1,0.0,Vector3D(0,0,0)))
    {
        
        MeshParams mp({qchem::MHL,500,3,2.0,qchem::Gauss,1,0,0,3});
        mintegrator=new MeshIntegrator<double>(cl->CreateMesh(mp));
    }
    void Init(int N, double rmin, double rmax)
    {
        nlohmann::json js = {
        {"type",BasisSetAtom::Type::BSpline9},
        {"N", N}, {"rmin", rmin}, {"rmax", rmax},
        };
        bs=BasisSetAtom::Factory(js,75);
        for (auto io:bs->Iterate<ibs_t>()) itsIBSs.push_back(io);

        std::vector<double> knots=MakeLogKnots(rmin,rmax,K,N);
        splines=bspline::generateBSplines<K>(knots);
        // splines.erase(splines.begin());
    }

    static const spline_t& GetSpline(const BSpline_IBS<K>* eval,size_t index)
    {
        return (*eval)[index];
    }

    std::vector<double> MakeLogKnots(double rmin, double rmax, size_t SPLINE_ORDER, size_t numberOfGridPoints);
    template <class S,class B> rsmat_t MakeSMat(const std::vector<S>& splines, const B&);
   
    size_t LMax;
    BasisSet* bs;
    std::vector<const ibs_t*> itsIBSs;
    Cluster* cl;
    MeshIntegrator<double>* mintegrator;
    std::vector<double> knots;
    std::vector<spline_t> splines;
};

std::vector<double> BSplineTests::MakeLogKnots(double rmin, double rmax, size_t k, size_t Ngrid)
{
    assert(Ngrid>1);
    std::vector<double> knots;
    size_t numberOfZeros = 1,L=0;

    if (k + 1 > L)  numberOfZeros = k + 1 - L;

    for (size_t i = 0; i < numberOfZeros; i++) knots.push_back(0.0);

    // logarithmic step
    const double step =  pow(rmax / rmin, 1 / static_cast<double>(Ngrid-1));
    for (size_t i = 0; i < Ngrid-1; i++) 
        knots.push_back(rmin * pow(step, i));
    knots.push_back(rmax); //Make the last one exact
    // for (size_t i = 0; i < numberOfZeros; i++) knots.push_back(rmax);
    return knots;
}
template <class S,class B> rsmat_t BSplineTests::MakeSMat(const std::vector<S>& splines, const B& integrator)
{
    size_t N=splines.size();
    rsmat_t A(N);
    for (auto i:iv_t(0,A.rows()))
        for (auto j:iv_t(i,A.rows()))
            A(i,j)=integrator(splines[i],splines[j]);
    return A;
    
}

TEST_F(BSplineTests, Example1)
{
    static constexpr size_t SPLINE_ORDER = 3;
    using Spline = bspline::Spline<double, SPLINE_ORDER>;

    // Define knots vector.
    const std::vector<double> knots{0.0, 1.0, 2.0, 3.0, 4.0, 5.0};

    // Generate Splines.
    const std::vector<Spline> splines = bspline::generateBSplines<SPLINE_ORDER>(knots);
}

using namespace bspline::operators;
using namespace bspline::integration;

TEST_F(BSplineTests, Knots)
{
    
    static constexpr size_t SPLINE_ORDER = 3;
    std::vector<double> knots=MakeLogKnots(0.01,2000.0,SPLINE_ORDER,10);
    // cout << knots << endl;
    using Spline = bspline::Spline<double, SPLINE_ORDER>;
    std::vector<Spline> splines=bspline::generateBSplines<SPLINE_ORDER>(knots);
    // for (double r=0.01;r<2000;r*=2.0)
    //     cout << r << " " << splines[5](r) << endl;

    const BilinearForm bilinearForm{IdentityOperator{}};
    rsmat_t S=MakeSMat(splines,bilinearForm);
    // cout << "overlap=" << S << endl;
}

#include <map>
template <class T, size_t K> struct cmpSplines {
    bool operator()(const bspline::Spline<T,K>& a, const bspline::Spline<T,K>& b) const
    {
        return a.getSupport().getStartIndex()<b.getSupport().getStartIndex();
    }
};

TEST_F(BSplineTests, SplineMap)
{
    static constexpr size_t K = 3;
    std::vector<double> knots=MakeLogKnots(0.01,2000.0,K,10);
    using Spline = bspline::Spline<double, K>;
    std::vector<Spline> splines=bspline::generateBSplines<K>(knots);
    std::map<Spline,size_t,cmpSplines<double,K>> smap;
    size_t index=0;
    for (auto& s:splines)
    {
        auto i=smap.find(s);
        if (i==smap.end())
        {
            smap[s]=index;
        }
        index++;
    }
}


TEST_F(BSplineTests,GLQIntegration)
{
    Init(100,.0001,40);
    GLCache cache0(splines[0].getSupport().getGrid(),K+1); //K+1 is the minimum that works for w=x^0
    GLCache cache2(splines[0].getSupport().getGrid(),K+2); //K+2 is the minimum that works for w=x^2
    std::function< double (double)> w0 = [](double x){return 1.0;};
    std::function< double (double)> w2 = [](double x){return x*x;};
    
    cout.precision(3);
    cout << "Grid = ";
    auto& grid=splines[0].getSupport().getGrid();
    for (auto r:grid) cout << r << ",";
    cout << endl;
    cout << "spline 0 =" << splines[0].front() << "," << splines[0].back() << "B(0)" << splines[0](0.0) << endl;
    cout << "spline 1 =" << splines[1].front() << "," << splines[1].back() << "B(0)" << splines[1](0.0) << endl;
//
//  Test all combos of indefinite integrals.
//    
    size_t nerr0=0,nerr2=0;
    double maxerr0=0.0, maxerr2=0.0;
    for (auto spa:splines)
        for (auto spb:splines)
        {
            if (!spa.getSupport().calcIntersection(spb.getSupport()).containsIntervals()) continue;
            double Sab2=BilinearForm{IdentityOperator{}}(spa,spb);
            double Sab3=cache0.Integrate(w0,spa,spb);
            assert(Sab2!=0.0);
            assert(Sab3!=0.0);
            // EXPECT_NEAR(Sab2/Sab3,1.0,7e-14);
            double err0=fabs(Sab2/Sab3-1.0);
            if (err0>1e-14)
            {
                nerr0++;
                if (err0>maxerr0) maxerr0=err0;
            }

            Sab2=BilinearForm{X<2>{}}(spa,spb);
            Sab3=cache2.Integrate(w2,spa,spb);
            assert(Sab2!=0.0);
            assert(Sab3!=0.0);
            // EXPECT_NEAR(Sab2/Sab3,1.0,2e-13);
            double err2=fabs(Sab2/Sab3-1.0);
            if (err2>1e-14)
            {
                nerr2++;
                if (err2>maxerr2) maxerr2=err2;
            }
        }
    cout << "weight x^0: " << nerr0 << " error>1e-14, max=" << maxerr0 << endl;
    cout << "weight x^2: " << nerr2 << " error>1e-14, max=" << maxerr2 << endl;

    double rmin=0.5,rmax=5.0;
    for (auto spa:splines)
        for (auto spb:splines)
        {
            double Sdef=cache0.Integrate(w0,spa,spb,rmin,rmax);
            double Sind=cache0.Integrate(w0,spa,spb);
            EXPECT_LE(Sdef,Sind);
        }
}
TEST_F(BSplineTests, Overlap)
{
    cout << "Overlap ";
    cout.precision(3);
    Init(4+K,.1,40.);
    for (auto ibs:bs->Iterate<Real_OIBS >())
    {
        cout << ibs->GetSymmetry();
        rsmat_t S=ibs->Overlap();
        for (auto d:blaze::diagonal(S)) EXPECT_NEAR(d,1.0,1e-15);
        for (auto i:iv_t(0,S.rows()-K-1)) //Check banded
            for (auto j:iv_t(i+K+1,S.rows())) 
                EXPECT_EQ(S(i,j),0.0);
            
        rsmat_t Snum = mintegrator->Overlap(*ibs);
        EXPECT_NEAR(max(abs(S-Snum)),0.0,3e-6);

        // cout << "S=" << S << endl;
        // cout << "Snum=" << Snum << endl;
        // Now try GLQ integration.
        const BSpline_IBS<K>* eval=dynamic_cast<const BSpline_IBS<K>*>(ibs);
        auto ns=eval->Norm();
        auto grid=GetSpline(eval,0).getSupport().getGrid();
        GLCache cache2(grid,K+3);
        std::function< double (double)> w2 = [](double x){return x*x;};
        for (auto ia:iv_t(0,S.rows()))
            for (auto ib:iv_t(ia,S.rows())) 
            {
                double nab=ns[ia]*ns[ib]*4*M_PI;
                auto a=GetSpline(eval,ia),b=GetSpline(eval,ib);
                double Sab=cache2.Integrate(w2,a,b)*nab;
                EXPECT_NEAR(Sab,S(ia,ib),1e-14);
                Sab=cache2.Integrate(w2,a,b,grid.front(),grid.back())*nab;
                EXPECT_NEAR(Sab,S(ia,ib),1e-14);
                Sab=0.0;
                for (size_t ig=1;ig<grid.size();ig++)
                    Sab+=cache2.Integrate(w2,a,b,grid[ig-1],grid[ig]);
                Sab*=nab;
                EXPECT_NEAR(Sab,S(ia,ib),1e-14);
            }
        
        
    }
     cout << endl;
}
TEST_F(BSplineTests, Nuclear)
{
    cout << "Nuclear ";
    cout.precision(3);
    Init(4+K,.1,40.);
    for (auto ibs:bs->Iterate<const Real_OIBS>())
    {
        cout << ibs->GetSymmetry();
        // const Real_OIBS* ibs1=ibs;
        rsmat_t Ven=ibs->Nuclear(cl);
        for (auto i:iv_t(0,Ven.rows()-K-1)) //Check banded
            for (auto j:iv_t(i+K+1,Ven.rows())) EXPECT_EQ(Ven(i,j),0.0);
        
        rsmat_t Vennum = -cl->GetNuclearCharge()*mintegrator->Inv_r1(*ibs);
        EXPECT_NEAR(max(abs(Ven-Vennum)),0.0,1e-7);

        // cout << "Ven=" << Ven << endl;
        // cout << "Vennum=" << Vennum << endl;
    }
    cout << endl;
}
TEST_F(BSplineTests, Kinetic)
{
    cout << "Kinetic ";
    cout.precision(3);
    Init(4+K,.1,40.);
    for (auto ibs:bs->Iterate<const Real_OIBS>())
    {
        cout << ibs->GetSymmetry();
        rsmat_t T=ibs->Kinetic();
        for (auto i:iv_t(0,T.rows()-K-1)) //Check banded
            for (auto j:iv_t(i+K+1,T.rows())) EXPECT_EQ(T(i,j),0.0);
        
        int l=ibs->CastSymmetry<Angular_Sym>().GetL();
        rsmat_t Tnum = mintegrator->Grad2(*ibs);
        rsmat_t Cen  = mintegrator->Inv_r2(*ibs);
        Tnum+=l*(l+1)*Cen;
        EXPECT_NEAR(max(abs(T-Tnum)),0.0,3e-5);
        
        // cout << "T=" << T << endl;
        // cout << "Tnum=" << Tnum << endl;
    }
    cout << endl;
}

