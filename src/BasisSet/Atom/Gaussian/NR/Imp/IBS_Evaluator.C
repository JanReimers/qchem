// File: BasisSet/Atom/radial/Imp/Gaussian_IBS.C
module;
#include <valarray>
#include <cmath>
#include <cassert>
#include <iostream>
module BasisSet.Atom.Gaussian.NR.IBS_EValuator;
import qchem.BasisSet.Atom.Gaussian.Rk;
import qchem.BasisSet.Atom.GaussianIntegrals;
import Common.Constants;
import Common.IntPow;


inline double Overlap(double ea , double eb,size_t l_total)
{
    return Gaussian::Integral(ea+eb,l_total); //Already has 4*Pi and r^2 from dr.
}

inline double Grad2(double ea , double eb,size_t la, size_t lb)
{
    assert(la==lb);
    double t=ea+eb;
    size_t l1=la+1;
    return  l1*l1     * Gaussian::Integral(t,2*la-2)
            -2*l1 * t * Gaussian::Integral(t,2*la  )
            +4*ea*eb  * Gaussian::Integral(t,2*la+2);
}

inline double Inv_r2(double ea , double eb,size_t l_total)
{
    return Gaussian::Integral(ea+eb,l_total-2); //Already has 4*Pi
}

inline double Repulsion(double eab, double ec,size_t la,size_t lc)
{    
    Gaussian::RkEngine cd(eab,ec,std::max(la,lc));
    return FourPi2*cd.Coulomb_R0(la,lc);
}

inline double Charge(double ea, size_t l)
{
    return ::Gaussian::Integral(ea,l);
}
//---------------------------------------------------------------------------
//
//  Start member functions.
//

// This need overridability.
double Gaussian_IBS::Inv_r1(double ea , double eb,size_t l_total) const
{
    return Gaussian::Integral(ea+eb,l_total-1); //Already has 4*Pi
}


Gaussian_IBS::ds_t Gaussian_IBS::norms() const
{
    size_t N=es.size();    
    ds_t ret(N);
    for (size_t i=0;i<N;i++) ret[i]=1.0/sqrt(::Overlap(es[i],es[i],2*l)); 
    return ret;
}

Gaussian_IBS::omls_t Gaussian_IBS::Overlap() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Overlap(es[i-1],es[j-1],2*l)*ns[i-1]*ns[j-1];

    return S;
}

Gaussian_IBS::omls_t Gaussian_IBS::Grad2() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Grad2(es[i-1],es[j-1],l,l)*ns[i-1]*ns[j-1];

    return S;
}

Gaussian_IBS::omls_t Gaussian_IBS::Inv_r1() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= Inv_r1(es[i-1],es[j-1],2*l)*ns[i-1]*ns[j-1];

    return S;
}

Gaussian_IBS::omls_t Gaussian_IBS::Inv_r2() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Inv_r2(es[i-1],es[j-1],2*l)*ns[i-1]*ns[j-1];

    return S;
}

Gaussian_IBS::omls_t Gaussian_IBS::Repulsion() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Repulsion(es[i-1],es[j-1],l,l)*ns[i-1]*ns[j-1];

    return S;
}

Gaussian_IBS::omlv_t Gaussian_IBS::Charge() const
{
    omlv_t V(size());
    for (auto i:V.indices())
            V(i)= ::Charge(es[i-1],l)*ns[i-1];

    return V;
}

IBS_Evaluator::omlm_t Gaussian_IBS::XRepulsion(const Fit_IBS& _b) const
{
    const Gaussian_IBS& b=dynamic_cast<const Gaussian_IBS&>(_b);
    omlm_t M(size(),b.size());
    for (auto i:M.rows())
            for (auto j:M.cols())
                M(i,j)=::Repulsion(es[i-1],b.es[j-1],l,b.l)*ns[i-1]*b.ns[j-1];
    return M;
}

IBS_Evaluator::omlm_t Gaussian_IBS::XKinetic(const Orbital_RKBS_IBS<double>* _b) const
{
    const Gaussian_IBS* b=dynamic_cast<const Gaussian_IBS*>(_b);
    assert(b);
    assert(l==b->l);
    omlm_t M(size(),b->size());
    for (auto i:M.rows())
            for (auto j:M.cols())
                M(i,j)=(::Grad2(es[i-1],b->es[j-1],l,l) + l*(l+1)*::Inv_r2(es[i-1],b->es[j-1],2*l))*ns[i-1]*b->ns[j-1];
    return M;
}


dERI3 Gaussian_IBS::Overlap  (const Fit_IBS& _c) const
{
    const Gaussian_IBS& c=dynamic_cast<const Gaussian_IBS&>(_c);
    dERI3 S3;
    for (size_t ic=0;ic<c.size();ic++) 
    {
        omls_t S(size());
        for (auto i:S.rows())
            for (auto j:S.cols(i))
                S(i,j)=::Overlap(es[i-1]+es[j-1],c.es[ic],l+l+c.l)*ns[i-1]*ns[j-1]*c.ns[ic];  
        
        S3.push_back(S);
    }
    return S3;
}
dERI3 Gaussian_IBS::Repulsion(const Fit_IBS& _c) const
{
    const Gaussian_IBS& c=dynamic_cast<const Gaussian_IBS&>(_c);
    dERI3 S3;
    for (size_t ic=0;ic<c.size();ic++) 
    {
        omls_t S(size());
        for (auto i:S.rows())
            for (auto j:S.cols(i))
                S(i,j)=::Repulsion(es[i-1]+es[j-1],c.es[ic],l,c.l)*ns[i-1]*ns[j-1]*c.ns[ic];  
        
        S3.push_back(S);
    }
    return S3;
}


template <class v> v gaussian(double r,size_t l,const v& e, const v& n)
{
    return n*uintpow(r,l)*exp(-e*r*r);
}
template <class v> v grad_gaussian(double r,size_t l,const v& e, const v& n)
{
    double lr= r==0 ? 0 : l/r;
    return (lr-2*r*e)*gaussian(r,l,e,n);
}

template <class T> Vector<T> convert(const std::valarray<T>& v) 
{
    Vector<T> ret(v.size());
    size_t i=0;
    for (auto vi:v) ret(++i)=vi;
    return ret;
}


Gaussian_IBS::Vec    Gaussian_IBS::operator() (const RVec3& r) const
{
    return convert(gaussian(norm(r),l,es,ns));
}

Gaussian_IBS::Vec3Vec Gaussian_IBS::Gradient(const RVec3& r) const
{
    Vec3Vec ret(size());
    double mr=norm(r);
    if (mr==0.0)
    {
        Fill(ret,RVec3{0,0,0});
        return ret;
    }
    ds_t grad=grad_gaussian(norm(r),l,es,ns);
    RVec3 rhat=r/norm(r);
    size_t i=0;
    for (auto& g:grad) ret(++i)=g*rhat;
    return ret;
}

std::ostream&  Gaussian_IBS::Write(std::ostream& os) const
{
    return os << " with " << size() << " basis functions, alpha={" << es[0] << " ... " << es[size()-1] << "}" << std::endl;
}



Gaussian_IBS::ds_t Gaussian_RKBS_IBS::norms() const
{
    size_t N=es.size();
    ds_t ret(N);
    for (size_t i=0;i<N;i++) 
    {
        double k=::Grad2(es[i],es[i],l,l) + l*(l+1)*::Inv_r2(es[i],es[i],2*l);
        ret[i]=1.0/sqrt(k); 
    }
    return ret;
}

double Gaussian_RKBS_IBS::Inv_r1(double ea , double eb,size_t l_total) const
{
    return 4*ea*eb*::Gaussian::Integral(ea+eb,l_total+1); //Don't count the r^2 in dr^3
}


Gaussian_IBS::ds_t Gaussian_RKBS_IBS::eval(const RVec3& r) const
{
    double mr=norm(r);
    ds_t f=-2*es*mr;
    if (kappa >0) 
        f+=(2*kappa+1)/mr;
        
    return f*gaussian(mr,l,es,ns);
}

Gaussian_IBS::Vec    Gaussian_RKBS_IBS::operator() (const RVec3& r) const
{
   return convert(eval(r)); 
}

Gaussian_IBS::Vec3Vec Gaussian_RKBS_IBS::Gradient(const RVec3& r) const
{
    Vec3Vec ret(size());
    double mr=norm(r);
    if (mr==0.0)
    {
        Fill(ret,RVec3{0,0,0});
        return ret;
    }
    ds_t grad=eval(r)*(l/mr-2*mr*es);
   
    // std::cout << "grad=";
    // for (auto g:grad) std::cout << g << " ";
    // std::cout << std::endl;
    RVec3 rhat=r/norm(r);
    size_t i=0;
    for (auto& g:grad) ret(++i)=g*rhat;
    // std::cout << "ret=";
    // for (auto g:ret) std::cout << g << " ";
    // std::cout << std::endl;
    return ret;
}
