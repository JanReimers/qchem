// File: BSplines.C  Test the BSplinebasis package included as a submodule


#include "gtest/gtest.h"
#include "Imp/BasisSet/Atom/l/BSpline_BS.H"
#include "Imp/BasisSet/Atom/l/BSpline_IBS.H"
#include "Imp/Containers/stl_io.h"
#include "Imp/Integrals/MeshIntegrator.H"
#include "Imp/Cluster/Atom.H"
#include "Imp/Cluster/Molecule.H"
#include <MeshParams.H>
#include <Cluster.H>

#include "oml/smatrix.h"
#include <bspline/Core.h>
#include <iostream>

using std::cout;
using std::endl;

class BSplineTests : public ::testing::Test
{
public:
    BSplineTests() 
        : LMax(0)
        , cl(new Molecule())
    {
        StreamableObject::SetToPretty();
        cl->Insert(new Atom(1,0.0,Vector3D(0,0,0)));
        MeshParams mp({qchem::MHL,200,3,2.0,qchem::Gauss,1,0,0,3});
        mintegrator=new MeshIntegrator<double>(cl->CreateMesh(mp));
    }
    void Init(int N, double rmin, double rmax)
    {
        bs=new Atoml::BSpline::BasisSet<K>(N,rmin,rmax,LMax);
    }

    std::vector<double> MakeLogKnots(double rmin, double rmax, size_t SPLINE_ORDER, size_t numberOfGridPoints);
    template <class S,class B> SMatrix<double> MakeSMat(const std::vector<S>& splines, const B&);
   
    static constexpr size_t K=6;
    size_t LMax;
    BasisSet* bs;
    Cluster* cl;
    MeshIntegrator<double>* mintegrator;
};

std::vector<double> BSplineTests::MakeLogKnots(double rmin, double rmax, size_t k, size_t Ngrid)
{
    std::vector<double> knots;
    size_t numberOfZeros = 1,L=1;

    if (k + 1 > L)  numberOfZeros = k + 1 - L;

    for (size_t i = 0; i < numberOfZeros; i++) knots.push_back(0.0);

    // logarithmic step
    const double step =  pow(rmax / rmin, 1 / static_cast<double>(Ngrid));
    for (int i = 0; i <= Ngrid; i++) 
        knots.push_back(rmin * pow(step, i));
    return knots;
}
template <class S,class B> SMatrix<double> BSplineTests::MakeSMat(const std::vector<S>& splines, const B& integrator)
{
    size_t N=splines.size();
    SMatrix<double> A(N);
    for (size_t i:A.rows())
        for (size_t j:A.cols(i))
            A(i,j)=integrator(splines[i-1],splines[j-1]);
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
    SMatrix<double> S=MakeSMat(splines,bilinearForm);
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
TEST_F(BSplineTests, Overlap)
{
    cout.precision(3);
    Init(10,.1,40.);
    for (auto ibs:bs->Iterate<Atoml::BSpline::Orbital_IBS<K> >())
    {
        const TOrbital_IBS<double>* ibs1=ibs;
        SMatrix<double> S=ibs1->Overlap();
        for (auto d:Vector<double>(S.GetDiagonal())) EXPECT_NEAR(d,1.0,1e-15);
        for (auto i:S.rows()) //Check banded
            for (auto j:S.cols(i+K+1)) EXPECT_EQ(S(i,j),0.0);
        
        SMatrix<double> Snum = mintegrator->Overlap(*ibs);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,1e-8);

        // cout << "S=" << S << endl;
        // cout << "Snum=" << Snum << endl;
    }
}
TEST_F(BSplineTests, Nuclear)
{
    cout.precision(3);
    Init(10,.1,40.);
    for (auto ibs:bs->Iterate<Atoml::BSpline::Orbital_IBS<K> >())
    {
        const TOrbital_IBS<double>* ibs1=ibs;
        SMatrix<double> Ven=ibs1->Nuclear(cl);
        for (auto i:Ven.rows()) //Check banded
            for (auto j:Ven.cols(i+K+1)) EXPECT_EQ(Ven(i,j),0.0);
        
        SMatrix<double> Vennum = -cl->GetNuclearCharge()*mintegrator->Nuclear(*ibs);
        EXPECT_NEAR(Max(fabs(Ven-Vennum)),0.0,1e-7);

        // cout << "Ven=" << Ven << endl;
        // cout << "Vennum=" << Vennum << endl;
    }
}
TEST_F(BSplineTests, Kinetic)
{
    cout.precision(3);
    Init(10,.1,40.);
    for (auto ibs:bs->Iterate<Atoml::BSpline::Orbital_IBS<K> >())
    {
        const TOrbital_IBS<double>* ibs1=ibs;
        SMatrix<double> T=ibs1->Grad2();
        for (auto i:T.rows()) //Check banded
            for (auto j:T.cols(i+K+1)) EXPECT_EQ(T(i,j),0.0);
        
        SMatrix<double> Tnum = mintegrator->Grad(*ibs);
        EXPECT_NEAR(Max(fabs(T-Tnum)),0.0,3e-5);

        // cout << "T=" << T << endl;
        // cout << "Tnum=" << Tnum << endl;
    }
}