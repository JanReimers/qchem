// File: BasisSet/Atom/radial/Imp/BSpline_IBS.C
module;
#include <bspline/Core.h>
#include <valarray>
#include <cmath>
#include <cassert>
#include <iostream>
module BasisSet.Atom.BSpline.NR.IBS_Evaluator;
import qchem.BasisSet.Atom.BSpline.Rk;
import qchem.BasisSet.Atom.BSpline.SplineGrouper;
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

template <size_t K1,size_t K2> double Overlap(const spline_t<K1>& a , const spline_t<K2>& b,size_t l_total)
{
    return BilinearForm{X<2>{}}(a,b)*FourPi;
    // return BilinearForm{IdentityOperator{}}(a,b);
}

template <size_t K> double Grad2(const spline_t<K>& a , const spline_t<K>& b,size_t la, size_t lb)
{
    static const auto T = -X<2>{} * Dx<2>{} - 2 * X<1>{} * Dx<1>{};
    assert(la==lb);
    return BilinearForm{T}(a,b)*FourPi;
}

template <size_t K> double Inv_r1(const spline_t<K>& a , const spline_t<K>& b,size_t l_total)
{
    return BilinearForm{X<1>{}}(a,b)*FourPi; 
}
template <size_t K> double Inv_r2(const spline_t<K>& a , const spline_t<K>& b,size_t l_total)
{
    return BilinearForm{IdentityOperator{}}(a,b)*FourPi; 
}

template <size_t K> double Charge(const spline_t<K>& a , size_t l)
{
    return LinearForm{X<2>{}}(a)*FourPi;
}

//---------------------------------------------------------------------------
//
//  Start member functions.
//
template <size_t K> void BSpline_IBS<K>::Register(Grouper* _grouper)
{
    assert(_grouper);
    auto grouper=static_cast<SplineGrouper<K>*>(_grouper);
    assert(grouper);
    for (auto s:splines) es_indices.push_back(grouper->Insert(s,l));
    grouper->itsGLs[l]=itsGL.get();
}

template <size_t K> BSpline_IBS<K>::BSpline_IBS(size_t Ngrid, double _rmin, double _rmax, int l, const is_t& mls) 
: IBS_Evaluator(l,mls), rmin(_rmin), rmax(_rmax) 
{
    std::vector<double> knots=MakeLogKnots(Ngrid,rmin,rmax);
    // std::cout << "Knots=" << knots << std::endl;
    splines=bspline::generateBSplines<K>(knots);
    auto grid=splines[0].getSupport().getGrid();
    // std::cout << "Grid = " << grid.size() << "    ";
    // for (auto r:grid) std::cout << r << ",";
    // std::cout << std::endl;
    itsGL.reset(new GLCache(grid,K+2*l));
    ns=norms();
    assert(size()==splines.size());
};

 template <size_t K> std::vector<double> BSpline_IBS<K>::MakeLogKnots(size_t Ngrid, double rmin, double rmax)
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

template <size_t K> BSpline_IBS<K>::ds_t BSpline_IBS<K>::norms() const
{
    size_t N=splines.size();
    ds_t ret(N);
    for (size_t i=0;i<N;i++) ret[i]=1.0/sqrt(::Overlap(splines[i],splines[i],2*l)); 
    return ret;
}

template <size_t K> rsmat_t BSpline_IBS<K>::Overlap() const
{
    size_t N=size();
    rsmat_t S(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
            S(i,j)= ::Overlap(splines[i],splines[j],2*l)*ns[i]*ns[j];

    return S;
}

template <size_t K> rsmat_t BSpline_IBS<K>::Grad2() const
{
    size_t N=size();
    rsmat_t S(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
            S(i,j)= ::Grad2(splines[i],splines[j],l,l)*ns[i]*ns[j];

    return S;
}

template <size_t K> rsmat_t BSpline_IBS<K>::Inv_r1() const
{
    size_t N=size();
    rsmat_t S(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
            S(i,j)= ::Inv_r1(splines[i],splines[j],2*l)*ns[i]*ns[j];

    return S;
}

template <size_t K> rsmat_t BSpline_IBS<K>::Inv_r2() const
{
    size_t N=size();
    rsmat_t S(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
            S(i,j)= ::Inv_r2(splines[i],splines[j],2*l)*ns[i]*ns[j];

    return S;
}

template <size_t K> rsmat_t BSpline_IBS<K>::Repulsion() const
{
    size_t N=size();
    rsmat_t S(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
            S(i,j)= ::Repulsion(splines[i],splines[j],l,l)*ns[i]*ns[j];

    return S;
}

template <size_t K> rvec_t BSpline_IBS<K>::Charge() const
{
    rvec_t V(size());
    for (auto i:iv_t(0,size()))
            V[i]= ::Charge(splines[i],l)*ns[i];

    return V;
}

template <size_t K> rmat_t BSpline_IBS<K>::XRepulsion(const Fit_IBS& _b) const
{
    const BSpline_IBS<K>& b=dynamic_cast<const BSpline_IBS<K>&>(_b);
    size_t Nr=size(), Nc=b.size();
    rmat_t M(Nr,Nc);
    for (auto i:iv_t(0,Nr))
            for (auto j:iv_t(0,Nc))
                M(i,j)=::Repulsion(splines[i],b.splines[j],l,b.l)*ns[i]*b.ns[j];
    return M;
}

template <size_t K> rmat_t BSpline_IBS<K>::XKinetic(const Orbital_RKBS_IBS<double>* _b) const
{
    const BSpline_IBS<K>* b=dynamic_cast<const BSpline_IBS<K>*>(_b);
    assert(b);
    assert(l==b->l);
    size_t Nr=size(), Nc=b->size();
    rmat_t M(Nr,Nc);
    for (auto i:iv_t(0,Nr))
            for (auto j:iv_t(0,Nc))
                M(i,j)=(::Grad2(splines[i],b->splines[j],l,l) + l*(l+1)*::Inv_r2(splines[i],b->splines[j],2*l))*ns[i]*b->ns[j];
    return M;
}

template <size_t K> dERI3 BSpline_IBS<K>::Overlap(const Fit_IBS& _c) const
{
    const BSpline_IBS<K>& c=dynamic_cast<const BSpline_IBS<K>&>(_c);
    dERI3 S3;
    for (size_t ic=0;ic<c.size();ic++) 
    {
        omls_t S(size());
        for (auto i:S.rows())
            for (auto j:S.cols(i))
            {
                auto ab=splines[i-1]+splines[j-1];
                S(i,j)=::Overlap(ab,c.splines[ic],l+l+c.l)*ns[i-1]*ns[j-1]*c.ns[ic];  
            }
        
        S3.push_back(S);
    }
    return S3;
}
template <size_t K> dERI3 BSpline_IBS<K>::Repulsion(const Fit_IBS& _c) const
{
    const BSpline_IBS<K>& c=dynamic_cast<const BSpline_IBS<K>&>(_c);
    dERI3 S3;
    for (size_t ic=0;ic<c.size();ic++) 
    {
        omls_t S(size());
        // for (auto i:S.rows())
        //     for (auto j:S.cols(i))
        //         S(i,j)=::Repulsion(splines[i-1]*splines[j-1],c.splines[ic],l,c.l)*ns[i-1]*ns[j-1]*c.ns[ic];  
        
        S3.push_back(S);
    }
    return S3;
}

template <size_t K> BSpline_IBS<K>::Vec    BSpline_IBS<K>::operator() (const RVec3& r) const
{
    Vec ret(size());
    double mr=norm(r);
    size_t i=0;
    for (auto s:splines) 
    {
        ++i;
        ret(i)=ns[i-1]*s(mr);
    }
    return ret;
}

template <size_t K> BSpline_IBS<K>::Vec3Vec BSpline_IBS<K>::Gradient(const RVec3& r) const
{
    Vec3Vec ret(size());
    double mr=norm(r);
    if (mr==0.0) 
    {
        
        Fill(ret,RVec3(0,0,0));
        return ret; //Cusp at the origin so grad is undefined.
    }
    assert(mr>0);
    Fill(ret,r/mr);
    size_t i=0;
    for (auto s:splines) 
    {
        auto dsdx=transformSpline(bspline::operators::Dx<1>{},s);
        ++i;
        ret(i)*=ns[i-1]*dsdx(mr);
    }
    return ret;
}

template <size_t K> std::ostream&  BSpline_IBS<K>::Write(std::ostream& os) const
{
    return os << " with " << size() << " basis functions, {" << rmin << " ... " << rmax << "}" << std::endl;
}


#define INSTANCEk(k) template class BSpline_IBS<k>;
#include "../../Instance.hpp"