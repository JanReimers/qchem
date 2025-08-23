// File: BasisSet/Atom/radial/Imp/BSpline_IBS.C
module;
#include <bspline/Core.h>
#include <valarray>
#include <cmath>
#include <cassert>
using namespace bspline::operators; 
using namespace bspline::integration; 
module BasisSet.Atom.BSpline_IBS;
import qchem.Basisset.Atom.radial.BSpline.Rk;
import Common.Constants;
// import Common.IntPow;

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
template <size_t K> void BSpline_IBS<K>::Register(ExponentGrouper* _grouper)
{
    assert(_grouper);
    grouper=_grouper;
    for (auto s:splines) es_indices.push_back(_grouper->Insert(s.getSupport().front(),l));
}

template <size_t K> BSpline_IBS<K>::BSpline_IBS(size_t Ngrid, double _rmin, double _rmax, int _l, const is_t& _mls) 
: rmin(_rmin), rmax(_rmax), l(_l), mls(_mls),ns(norms()) 
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
    ds_t ret(size());
    for (size_t i=0;i<size();i++) ret[i]=1.0/sqrt(::Overlap(splines[i],splines[i],2*l)); 
    return ret;
}

template <size_t K> BSpline_IBS<K>::omls_t BSpline_IBS<K>::Overlap() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Overlap(splines[i-1],splines[j-1],2*l)*ns[i-1]*ns[j-1];

    return S;
}

template <size_t K> BSpline_IBS<K>::omls_t BSpline_IBS<K>::Grad2() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Grad2(splines[i-1],splines[j-1],l,l)*ns[i-1]*ns[j-1];

    return S;
}

template <size_t K> BSpline_IBS<K>::omls_t BSpline_IBS<K>::Inv_r1() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Inv_r1(splines[i-1],splines[j-1],2*l)*ns[i-1]*ns[j-1];

    return S;
}

template <size_t K> BSpline_IBS<K>::omls_t BSpline_IBS<K>::Inv_r2() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Inv_r2(splines[i-1],splines[j-1],2*l)*ns[i-1]*ns[j-1];

    return S;
}

template <size_t K> BSpline_IBS<K>::omls_t BSpline_IBS<K>::Repulsion() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Repulsion(splines[i-1],splines[j-1],l,l)*ns[i-1]*ns[j-1];

    return S;
}

template <size_t K> BSpline_IBS<K>::omlv_t BSpline_IBS<K>::Charge() const
{
    omlv_t V(size());
    for (auto i:V.indices())
            V(i)= ::Charge(splines[i-1],l)*ns[i-1];

    return V;
}

template <size_t K> IBS_Evaluator::omlm_t BSpline_IBS<K>::XRepulsion(const Fit_IBS& _b) const
{
    const BSpline_IBS<K>& b=dynamic_cast<const BSpline_IBS<K>&>(_b);
    omlm_t M(size(),b.size());
    for (auto i:M.rows())
            for (auto j:M.cols())
                M(i,j)=::Repulsion(splines[i-1],b.splines[j-1],l,b.l)*ns[i-1]*b.ns[j-1];
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

template <size_t K> Rk* BSpline_IBS<K>::CreateRk(size_t ia,size_t ic,size_t ib,size_t id) const
{
    assert(false);
    return 0;
    // assert(grouper);
    // assert(itsRkCache);
    // size_t lmax=grouper->LMax(ia,ib,ic,id);
    // const GLCache* gl=this->GetGL(lmax);
    // return new BSpline::RkEngine(grouper->unique_spv,ia,ib,ic,id,lmax,*gl,*itsRkCache);
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



#define INSTANCEk(k) template class BSpline_IBS<k>;
#include "../BSpline/Instance.hpp"