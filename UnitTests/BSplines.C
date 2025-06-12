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
#include <Symmetry.H>

#include "oml/smatrix.h"
#include <bspline/Core.h>
#include <iostream>

using std::cout;
using std::endl;

class BSplineTests : public ::testing::Test
{
public:
    static constexpr size_t K=6; //Spline order.
    typedef bspline::Spline<double, K> spline_t;
    typedef TOrbital_IBS<double> ibs_t;
    typedef SMatrix<double> smat_t;
    BSplineTests() 
        : LMax(4)
        , cl(new Molecule())
    {
        StreamableObject::SetToPretty();
        cl->Insert(new Atom(1,0.0,Vector3D(0,0,0)));
        MeshParams mp({qchem::MHL,500,3,2.0,qchem::Gauss,1,0,0,3});
        mintegrator=new MeshIntegrator<double>(cl->CreateMesh(mp));
    }
    void Init(int N, double rmin, double rmax)
    {
        bs=new Atoml::BSpline::BasisSet<K>(N,rmin,rmax,LMax);
        for (auto io:bs->Iterate<ibs_t>()) itsIBSs.push_back(io);

        std::vector<double> knots=MakeLogKnots(0.01,2000.0,K,10);
        splines=bspline::generateBSplines<K>(knots);
    }

    std::vector<double> MakeLogKnots(double rmin, double rmax, size_t SPLINE_ORDER, size_t numberOfGridPoints);
    template <class S,class B> SMatrix<double> MakeSMat(const std::vector<S>& splines, const B&);
   
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
    for (size_t i = 0; i < Ngrid; i++) 
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

#include "Imp/BasisSet/Atom/radial/BSpline/GLQuadrature.H"
#include <Imp/BasisSet/Atom/Angular.H>
TEST_F(BSplineTests,GLQIntegration)
{
    Init(10,.1,10);
    GLCache cache0(splines[0].getSupport().getGrid(),K+1);
    GLCache cache2(splines[0].getSupport().getGrid(),K+3);
    std::function< double (double)> w0 = [](double x){return 1.0;};
    std::function< double (double)> w2 = [](double x){return x*x;};
//
//  Test all combos of indefinite integrals.
//    
    double max_error0=0.0,max_error2=0.0;
    for (auto spa:splines)
        for (auto spb:splines)
        {
            double Sab2=BilinearForm{IdentityOperator{}}(spa,spb);
            double Sab3=cache0.Integrate(w0,spa,spb);
            if (Sab3!=0)
                EXPECT_NEAR(Sab2/Sab3,1.0,4e-14);
            else
                EXPECT_NEAR(Sab2,0.0,1e-16);
            max_error0=std::max(max_error0,fabs(Sab2-Sab3)/Sab2);
            Sab2=BilinearForm{X<2>{}}(spa,spb);
            Sab3=cache2.Integrate(w2,spa,spb);
            if (Sab3!=0)
                EXPECT_NEAR(Sab2/Sab3,1.0,2e-13);
            else
                EXPECT_NEAR(Sab2,0.0,1e-16);
            max_error2=std::max(max_error2,fabs(Sab2-Sab3)/Sab2);
        }
    // cout << "max_error 0,2=" << max_error0 << "," << max_error2 << endl;
    // Now try some Definite integrals.
    cout.precision(3);
    cout << "Grid = ";
    auto& grid=splines[0].getSupport().getGrid();
    for (auto r:grid) cout << r << ",";
    cout << endl;
    // size_t Nsp=splines.size(), Ng=grid.size();
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
    Init(10,.1,40.);
    for (auto ibs:bs->Iterate<TOrbital_IBS<double> >())
    {
        cout << *ibs->GetSymmetry();
        SMatrix<double> S=ibs->Overlap();
        for (auto d:Vector<double>(S.GetDiagonal())) EXPECT_NEAR(d,1.0,1e-15);
        for (auto i:S.rows()) //Check banded
            for (auto j:S.cols(i+K+1)) EXPECT_EQ(S(i,j),0.0);
        
        SMatrix<double> Snum = mintegrator->Overlap(*ibs);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,3e-6);

        // cout << "S=" << S << endl;
        // cout << "Snum=" << Snum << endl;
        // Now try GLQ integration.
        const BSpline::IrrepIEClient<K>* iec=dynamic_cast<const BSpline::IrrepIEClient<K>*>(ibs);
        auto grid=iec->splines[0].getSupport().getGrid();
        GLCache cache2(grid,K+3);
        std::function< double (double)> w2 = [](double x){return x*x;};
        for (auto ia:S.rows())
            for (auto ib:S.cols(ia)) 
            {
                double nab=iec->ns(ia)*iec->ns(ib)*4*M_PI;
                auto a=iec->splines[ia-1],b=iec->splines[ib-1];
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
    Init(10,.1,40.);
    for (auto ibs:bs->Iterate<Atoml::BSpline::Orbital_IBS<K> >())
    {
        cout << *ibs->GetSymmetry();
        const TOrbital_IBS<double>* ibs1=ibs;
        SMatrix<double> Ven=ibs1->Nuclear(cl);
        for (auto i:Ven.rows()) //Check banded
            for (auto j:Ven.cols(i+K+1)) EXPECT_EQ(Ven(i,j),0.0);
        
        SMatrix<double> Vennum = -cl->GetNuclearCharge()*mintegrator->Inv_r1(*ibs);
        EXPECT_NEAR(Max(fabs(Ven-Vennum)),0.0,1e-7);

        // cout << "Ven=" << Ven << endl;
        // cout << "Vennum=" << Vennum << endl;
    }
    cout << endl;
}
TEST_F(BSplineTests, Kinetic)
{
    cout << "Kinetic ";
    cout.precision(3);
    Init(10,.1,40.);
    for (auto ibs:bs->Iterate<Atoml::BSpline::Orbital_IBS<K> >())
    {
        cout << *ibs->GetSymmetry();
        const TOrbital_IBS<double>* ibs1=ibs;
        SMatrix<double> T=ibs1->Kinetic();
        for (auto i:T.rows()) //Check banded
            for (auto j:T.cols(i+K+1)) EXPECT_EQ(T(i,j),0.0);
        
        int l=dynamic_cast<const Angular_Sym* >(ibs->GetSymmetry().get())->GetL();
        SMatrix<double> Tnum = mintegrator->Grad(*ibs);
        SMatrix<double> Cen = mintegrator->Inv_r2(*ibs);
        Tnum+=l*(l+1)*Cen;
        EXPECT_NEAR(Max(fabs(T-Tnum)),0.0,3e-5);
        
        // cout << "T=" << T << endl;
        // cout << "Tnum=" << Tnum << endl;
    }
    cout << endl;
}

