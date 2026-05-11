// File: BasisSet1/Atom/Evaluators/BSpline/Imp/IBS_r_Evaluator.C
module;
#include <bspline/Core.h>
#include <InvPosition.H> // 1/x^n operator are not provided in bspline package, so a roll our own.
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

template <size_t K> double Repulsion(const spline_t<K>& ab , const spline_t<K>& c,size_t la,size_t lc)
{    
    assert(false);
    return 0.0;
}

template <size_t K1,size_t K2> double Overlap(const spline_t<K1>& a , const spline_t<K2>& b,size_t l_total,const GLCache& gl)
{
    std::function< double (double)> x0 = [](double r)
    {
        return 1.0;
    };
    return gl.Integrate(x0,a,b)*FourPi;
}

template <size_t K> double Grad2(const spline_t<K>& a , const spline_t<K>& b,size_t la, size_t lb,const GLCache& gl)
{
    // static const auto T = -X<2>{} * Dx<2>{} - 2 * X<1>{} * Dx<1>{};
    static const auto T = -Dx<2>{};
    assert(la==lb);
    return (BilinearForm{T}(a,b))*FourPi;
}

template <size_t K> double Inv_r1(const spline_t<K>& a , const spline_t<K>& b,size_t l_total,const GLCache& gl)
{
    std::function< double (double)> xm1 = [](double r)
    {
        assert(r!=0.0);
        return 1.0/r;
    };
    return gl.Integrate(xm1,a,b)*FourPi;
}
template <size_t K> double Inv_r2(const spline_t<K>& a , const spline_t<K>& b,size_t l_total,const GLCache& gl)
{
    std::function< double (double)> xm2 = [](double r)
    {
        assert(r!=0.0);
        return 1.0/(r*r);
    };
    return gl.Integrate(xm2,a,b)*FourPi; 
}

template <size_t K> double Charge(const spline_t<K>& a , size_t l)
{
    return LinearForm{IdentityOperator{}}(a)*FourPi;
}

//---------------------------------------------------------------------------
//
//  Start member functions.
//
template <size_t K> void BSpline_r_IBS_Evaluator<K>::Register(Grouper* _grouper)
{
    assert(_grouper);
    auto grouper=static_cast<SplineGrouper<K>*>(_grouper);
    assert(grouper);
    for (auto s:splines) es_indices.push_back(grouper->Insert(s,l));
    grouper->itsGLs[l]=itsGL.get();
}

template <size_t K> BSpline_r_IBS_Evaluator<K>::BSpline_r_IBS_Evaluator(size_t Ngrid, double _rmin, double _rmax,const Irrep_QNs::sym_t& ylm) 
: IBS_Evaluator(ylm), rmin(_rmin), rmax(_rmax) 
{
    knots=MakeLogKnots(Ngrid,rmin,rmax);
    // std::cout << "Knots=" << knots << std::endl;
    splines=bspline::generateBSplines<K>(knots);
    splines.erase(splines.begin()); //First spline has B(0)=1.0 with violates B(0)=0 boundary condition for 1/r prefactor.
    splines.pop_back(); //Last spline has B(R)=1.0 with violates B(R)=0 boundary condition for 1/r prefactor.
    auto grid=splines[0].getSupport().getGrid();
    // std::cout << "Grid = " << grid.size() << "    ";
    // for (auto r:grid) std::cout << r << ",";
    // std::cout << std::endl;
    itsGL.reset(new GLCache(grid,K+1));
    ns=norms();
    assert(size()==splines.size());
};

 template <size_t K> std::vector<double> BSpline_r_IBS_Evaluator<K>::MakeLogKnots(size_t Ngrid, double rmin, double rmax)
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
    
    // cout << Ngrid << " " << l << " " << numberOfZeros << " ";
    if (numberOfZeros>Ngrid-numberOfZeros) numberOfZeros=Ngrid-numberOfZeros;
    if (numberOfZeros<1) numberOfZeros=1;
    // cout << numberOfZeros << endl;

    for (size_t i = 0; i < numberOfZeros; i++) knots.push_back(rmax);
    // std::cout << knots << std::endl;
    return knots;
}


template <size_t K> std::string BSpline_r_IBS_Evaluator<K>::Name () const
{
    std::ostringstream os;
    os << "BSpline<" << K << "> 1/r ";
    return os.str();
}

template <size_t K> std::string BSpline_r_IBS_Evaluator<K>::RadialID () const
{
    std::ostringstream os;
    os << Name() << " {";
    for (auto k:knots) os << k << " ";
    os << "}";
    return os.str();
}

template <size_t K> rvec_t BSpline_r_IBS_Evaluator<K>::norms() const
{
    size_t N=splines.size();
    rvec_t ret(N);
    for (size_t i=0;i<N;i++) ret[i]=1.0/sqrt(::Overlap(splines[i],splines[i],2*l,*itsGL)); 
    return ret;
}

template <size_t K> rsmat_t BSpline_r_IBS_Evaluator<K>::Overlap() const
{
    size_t N=size();
    rsmat_t S(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
            S(i,j)= ::Overlap(splines[i],splines[j],2*l,*itsGL)*ns[i]*ns[j];

    return S;
}

template <size_t K> rsmat_t BSpline_r_IBS_Evaluator<K>::Grad2() const
{
    size_t N=size();
    rsmat_t S(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
            S(i,j)= ::Grad2(splines[i],splines[j],l,l,*itsGL)*ns[i]*ns[j];

    return S;
}

template <size_t K> rsmat_t BSpline_r_IBS_Evaluator<K>::Inv_r1() const
{
    size_t N=size();
    rsmat_t S(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
            S(i,j)= ::Inv_r1(splines[i],splines[j],2*l,*itsGL)*ns[i]*ns[j];

    return S;
}

template <size_t K> rsmat_t BSpline_r_IBS_Evaluator<K>::Inv_r2() const
{
    size_t N=size();
    rsmat_t S(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
            S(i,j)= ::Inv_r2(splines[i],splines[j],2*l,*itsGL)*ns[i]*ns[j];

    return S;
}

template <size_t K> rsmat_t BSpline_r_IBS_Evaluator<K>::Repulsion() const
{
    size_t N=size();
    rsmat_t S(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
            S(i,j)= ::Repulsion(splines[i],splines[j],l,l)*ns[i]*ns[j];

    return S;
}

template <size_t K> rvec_t BSpline_r_IBS_Evaluator<K>::Charge() const
{
    rvec_t V(size());
    for (auto i:iv_t(0,size()))
            V[i]= ::Charge(splines[i],l)*ns[i];

    return V;
}

template <size_t K> rmat_t BSpline_r_IBS_Evaluator<K>::XRepulsion(const IBS_Evaluator& _b) const
{
    const BSpline_r_IBS_Evaluator<K>& b=dynamic_cast<const BSpline_r_IBS_Evaluator<K>&>(_b);
    size_t Nr=size(), Nc=b.size();
    rmat_t M(Nr,Nc);
    for (auto i:iv_t(0,Nr))
            for (auto j:iv_t(0,Nc))
                M(i,j)=::Repulsion(splines[i],b.splines[j],l,b.l)*ns[i]*b.ns[j];
    return M;
}

template <size_t K> rmat_t BSpline_r_IBS_Evaluator<K>::XKinetic(const IBS_Evaluator* _b) const
{
    const BSpline_r_IBS_Evaluator<K>* b=dynamic_cast<const BSpline_r_IBS_Evaluator<K>*>(_b);
    assert(b);
    assert(l==b->l);
    size_t Nr=size(), Nc=b->size();
    rmat_t M(Nr,Nc);
    for (auto i:iv_t(0,Nr))
            for (auto j:iv_t(0,Nc))
                M(i,j)=(::Grad2(splines[i],b->splines[j],l,l,*itsGL) + l*(l+1)*::Inv_r2(splines[i],b->splines[j],2*l,*itsGL))*ns[i]*b->ns[j];
    return M;
}

template <size_t K> dERI3 BSpline_r_IBS_Evaluator<K>::Overlap(const IBS_Evaluator& _c) const
{
    const BSpline_r_IBS_Evaluator<K>& c=dynamic_cast<const BSpline_r_IBS_Evaluator<K>&>(_c);
    dERI3 S3;
    size_t N=size();
    for (size_t ic=0;ic<c.size();ic++) 
    {
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
            {
                auto ab=splines[i-1]+splines[j-1];
                S(i,j)=::Overlap(ab,c.splines[ic],l+l+c.l,*itsGL)*ns[i-1]*ns[j-1]*c.ns[ic];  
            }
        
        S3.push_back(S);
    }
    return S3;
}
template <size_t K> dERI3 BSpline_r_IBS_Evaluator<K>::Repulsion(const IBS_Evaluator& _c) const
{
    const BSpline_r_IBS_Evaluator<K>& c=dynamic_cast<const BSpline_r_IBS_Evaluator<K>&>(_c);
    dERI3 S3;
    size_t N=size();
    for (size_t ic=0;ic<c.size();ic++) 
    {
        rsmat_t S(size());
        // for (auto i:iv_t(0,N))
        //     for (auto j:iv_t(i,N))
        //         S(i,j)=::Repulsion(splines[i]*splines[j],c.splines[ic],l,c.l)*ns[i]*ns[j]*c.ns[ic];  
        
        S3.push_back(S);
    }
    return S3;
}

template <size_t K> rvec_t BSpline_r_IBS_Evaluator<K>::operator() (const rvec3_t& r) const
{
    rvec_t ret(size());
    double mr=norm(r);
    size_t i=0;
    for (auto s:splines) 
    {
        ret[i]=ns[i]*s(mr)/mr;
        ++i;
    }
    return ret;
}

template <size_t K> rvec3vec_t BSpline_r_IBS_Evaluator<K>::Gradient(const rvec3_t& r) const
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

template <size_t K> std::ostream&  BSpline_r_IBS_Evaluator<K>::Write(std::ostream& os) const
{
    return os << " with " << size() << " basis functions, {" << rmin << " ... " << rmax << "}" << std::endl;
}


#define INSTANCEk(k) template class BSpline_r_IBS_Evaluator<k>;
#include "../Internal/Instance.hpp"