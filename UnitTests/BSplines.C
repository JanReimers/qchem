// File: BSplines.C  Test the BSplinebasis package included as a submodule


#include "gtest/gtest.h"
#include <bspline/Core.h>
#include "Imp/Containers/stl_io.h"
#include "oml/smatrix.h"
#include <iostream>

using std::cout;
using std::endl;

class BSplineTests : public ::testing::Test
{
public:
    BSplineTests() {};

    std::vector<double> MakeLogKnots(double rmin, double rmax, size_t SPLINE_ORDER, size_t numberOfGridPoints);
    template <class S,class B> SMatrix<double> MakeSMat(const std::vector<S>& splines, const B&);
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
    StreamableObject::SetToPretty();
    static constexpr size_t SPLINE_ORDER = 3;
    std::vector<double> knots=MakeLogKnots(0.01,2000.0,SPLINE_ORDER,10);
    cout << knots << endl;
    using Spline = bspline::Spline<double, SPLINE_ORDER>;
    std::vector<Spline> splines=bspline::generateBSplines<SPLINE_ORDER>(knots);
    for (double r=0.01;r<2000;r*=2.0)
        cout << r << " " << splines[5](r) << endl;

    const BilinearForm bilinearForm{IdentityOperator{}};
    double S11=bilinearForm(splines[0],splines[0]);
    cout << "overlap=" << S11 << endl;

    SMatrix<double> S=MakeSMat(splines,bilinearForm);
    cout << "overlap=" << S << endl;
}