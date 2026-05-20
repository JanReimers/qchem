// File: src/BasisSet1/Atom/Evaluators/BSpline/Imp/IBS_Evaluator.C
module;
#include <bspline/Core.h>
#include <cmath>
#include <cassert>
#include <iostream>
#include <functional>
#include <sstream>

module qchem.BasisSet1.Atom.Evaluators.BSpline.IBS;
import qchem.BasisSet1.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet1.Atom.Evaluators.BSpline.Internal.SplineGrouper;
import Common.Constants;
// import Common.IntPow;
using namespace bspline::operators; 
using namespace bspline::integration; 

template<size_t K> using spline_t = bspline::Spline<double, K>;

    // Alternate Overlap without operators.
    // std::function< double (double)> x2 = [](double r) {return r*r;};
    // return gl.Integrate(x2,a,b)*FourPi;

    // Alternate Grad2
    // std::function< double (double)> x1 = [](double r) {return r;};    
    // static const auto T = -X<2>{} * Dx<2>{};
    // assert(la==lb);
    // auto dbdx=transformSpline(bspline::operators::Dx<1>{},b);
    // double Iadb=gl.Integrate(x1,a,dbdx);
    // return (BilinearForm{T}(a,b) - 2*Iadb)*FourPi;

    // Alternate Inv_r1
    // std::function< double (double)> x1 = [](double r) {return r;};
    // return gl.Integrate(x1,a,b)*FourPi;

    // Alternate Inv_r2
    // std::function< double (double)> x0 = [](double r){return 1.0;};
    // return gl.Integrate(x0,a,b)*FourPi;    

//---------------------------------------------------------------------------
//
//  Start member functions.
//
template <size_t K> void BSpline_IBS_Evaluator<K>::Register(Grouper* _grouper)
{
    assert(_grouper);
    auto grouper=static_cast<SplineGrouper<K>*>(_grouper);
    assert(grouper);
    for (auto s:splines) es_indices.push_back(grouper->Insert(s,l));
}

template <size_t K> BSpline_IBS_Evaluator<K>::BSpline_IBS_Evaluator(size_t Ngrid, double _rmin, double _rmax,const Irrep_QNs::sym_t& ylm) 
: IBS_Evaluator(ylm), rmin(_rmin), rmax(_rmax) , itsGrid({0,1})
{
    knots=MakeLogKnots(Ngrid,rmin,rmax);
    // std::cout << "Knots=" << knots << std::endl;
    splines=bspline::generateBSplines<K>(knots);
    // splines.erase(splines.begin()); //First spline has B(0)=1.0 with violates B(0)=0 boundary condition for 1/r prefactor.
    for (size_t n=0;n<=3-l;n++) splines.pop_back(); //For s orbital the last spline has B(R)=1.0 with violates B(R)=0 boundary condition for 1/r prefactor.
    itsGrid=splines[0].getSupport().getGrid();
    // std::cout << "Grid = " << grid.size() << "    ";
    // for (auto r:grid) std::cout << r << ",";
    // std::cout << std::endl;
    ns=norms();
    // std::cout << "BSpline_IBS_Evaluator<K>::BSpline_IBS_Evaluator size=" << size() << std::endl;
    assert(size()==splines.size());
    assert(size()==ns.size());
};


template <size_t K> std::vector<double> BSpline_IBS_Evaluator<K>::MakeLogKnots(size_t Ngrid, double rmin, double rmax)
{
    assert(Ngrid>1);
    std::vector<double> knots;
    size_t numberOfZeros = 1;

    if (K + 1 > l)  numberOfZeros = K + 1 - l;

    for (size_t i = 0; i < numberOfZeros; i++) knots.push_back(0.0);

    // logarithmic step
    const double step =  pow(rmax / rmin, 1 / static_cast<double>(Ngrid-1));
    for (size_t i = 0; i < Ngrid-1; i++) //Skip 0.0 and rmax
        knots.push_back(rmin * pow(step, i));
    
    // std::cout << Ngrid << " " << l << " " << numberOfZeros << " ";
    // if (numberOfZeros>Ngrid-numberOfZeros) numberOfZeros=Ngrid-numberOfZeros;
    // if (numberOfZeros<1) numberOfZeros=1;
    //     std::cout << numberOfZeros << std::endl;

     for (size_t i = 0; i < numberOfZeros; i++) knots.push_back(rmax);
    // std::cout << knots << std::endl;
    return knots;
}

template <size_t K> std::string BSpline_IBS_Evaluator<K>::Name () const
{
    std::ostringstream os;
    os << "BSpline<" << K << "> ";
    return os.str();
}

template <size_t K> std::string BSpline_IBS_Evaluator<K>::RadialID () const
{
    std::ostringstream os;
    os << Name() << " {";
    for (auto k:knots) os << k << " ";
    os << "}";
    return os.str();
}

template <size_t K> std::string BSpline_IBS_Evaluator<K>::RadialType() const
{
    std::ostringstream os;
    os << "BS<" << K << "> grid=" << splines[0].getSupport().getGrid();;
    return os.str();
}

template <size_t K> Cache41*    BSpline_IBS_Evaluator<K>::MakeCache4() const
{
    return new BSpline_Cache4<K>(itsGrid);
}

template <size_t K> rvec_t BSpline_IBS_Evaluator<K>::norms() const
{
    size_t N=splines.size();
    // std::cout << "BSpline_IBS_Evaluator<K>::norms() N=" << N << std::endl;
    rvec_t ret(N);
    for (size_t i=0;i<N;i++) ret[i]=1.0/sqrt(BilinearForm{X<2>{}}(splines[i],splines[i])*FourPi); 
    return ret;
}


template <size_t K> rvec_t BSpline_IBS_Evaluator<K>::operator() (const rvec3_t& r) const
{
    rvec_t ret(size());
    double mr=norm(r);
    size_t i=0;
    for (auto s:splines) 
    {
        ret[i]=ns[i]*s(mr);
        ++i;
    }
    return ret;
}

template <size_t K> rvec3vec_t BSpline_IBS_Evaluator<K>::Gradient(const rvec3_t& r) const
{
    rvec3vec_t ret(size());
    double mr=norm(r);
    if (mr==0.0) 
    {
        
        ret=rvec3_t(0,0,0);
        return ret; //Cusp at the origin so grad is undefined.
    }
    assert(mr>0);
    ret=r/mr;
    size_t i=0;
    for (auto s:splines) 
    {
        auto dsdx=transformSpline(bspline::operators::Dx<1>{},s);
        ret[i]*=ns[i]*dsdx(mr);
        ++i;
    }
    return ret;
}

template <size_t K> std::ostream&  BSpline_IBS_Evaluator<K>::Write(std::ostream& os) const
{
    return os << " N= " << size() << " basis functions, {" << rmin << " ... " << rmax << "}" << std::endl;
}


#define INSTANCEk(k) template class BSpline_IBS_Evaluator<k>;
#include "../Internal/Instance.hpp"