// File: BasisSet/Atom/radial/Imp/Slater_IBS.C
module;
#include <valarray>
#include <cmath>
#include <cassert>
#include <iostream>
module BasisSet.Atom.Slater_IBS;
import qchem.BasisSet.Atom.Internal.radial.Slater.Rk;
import qchem.BasisSet.Atom.Internal.radial.Slater.Integrals;
import Common.Constants;
import Common.IntPow;
import qchem.stl_io;


inline double Overlap(double ea , double eb,size_t l_total)
{
    return Slater::Integral(ea+eb,l_total); //Already has 4*Pi and r^2 from dr.
}

inline double Grad2(double ea , double eb,size_t la, size_t lb)
{
    assert(la==lb);
    double ab=ea+eb;
    int l=la; //Safer to do formulas with int.
    // int ll=l*(l+1);
    double Term1=(l+1)*(l+1)*Slater::Integral(ab,2*l-2); //SlaterIntegral already has 4*Pi
    double Term2=-(l+1)*ab* Slater::Integral(ab,2*l-1);
    double Term3=ea*eb*Slater::Integral(ab,2*l);
    return Term1+Term2+Term3;
}

inline double Inv_r2(double ea , double eb,size_t l_total)
{
    return Slater::Integral(ea+eb,l_total-2); //Already has 4*Pi
}

inline double Repulsion(double eab, double ec,size_t la,size_t lc)
{    
    Slater::RkEngine cd(eab,ec,std::max(la,lc));
    return FourPi2*cd.Coulomb_R0(la,lc);
}

inline double Charge(double ea, size_t l)
{
    return ::Slater::Integral(ea,l);
}

//---------------------------------------------------------------------------
//
//  Start member functions.
//

// This need overridability.
double Slater_IBS::Inv_r1(double ea , double eb,size_t l_total) const
{
    return Slater::Integral(ea+eb,l_total-1); //Already has 4*Pi
}

Slater_IBS::ds_t Slater_IBS::norms() const
{
    size_t N=es.size();    
    ds_t ret(N);
    for (size_t i=0;i<N;i++) ret[i]=1.0/sqrt(::Overlap(es[i],es[i],2*l)); 
    return ret;
}

Slater_IBS::omls_t Slater_IBS::Overlap() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Overlap(es[i-1],es[j-1],2*l)*ns[i-1]*ns[j-1];

    return S;
}

Slater_IBS::omls_t Slater_IBS::Grad2() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Grad2(es[i-1],es[j-1],l,l)*ns[i-1]*ns[j-1];

    return S;
}

Slater_IBS::omls_t Slater_IBS::Inv_r1() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= Inv_r1(es[i-1],es[j-1],2*l)*ns[i-1]*ns[j-1];

    return S;
}

Slater_IBS::omls_t Slater_IBS::Inv_r2() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Inv_r2(es[i-1],es[j-1],2*l)*ns[i-1]*ns[j-1];

    return S;
}

Slater_IBS::omls_t Slater_IBS::Repulsion() const
{
    omls_t S(size());
    for (auto i:S.rows())
        for (auto j:S.cols(i))
            S(i,j)= ::Repulsion(es[i-1],es[j-1],l,l)*ns[i-1]*ns[j-1];

    return S;
}

IBS_Evaluator::omlm_t Slater_IBS::XRepulsion(const Fit_IBS& _b) const
{
    const Slater_IBS& b=dynamic_cast<const Slater_IBS&>(_b);
    omlm_t M(size(),b.size());
    for (auto i:M.rows())
            for (auto j:M.cols())
                M(i,j)=::Repulsion(es[i-1],b.es[j-1],l,b.l)*ns[i-1]*b.ns[j-1];
    return M;
}

IBS_Evaluator::omlm_t Slater_IBS::XKinetic(const Orbital_RKBS_IBS<double>* _b) const
{
    const Slater_IBS* b=dynamic_cast<const Slater_IBS*>(_b);
    assert(b);
    assert(l==b->l);
    omlm_t M(size(),b->size());
    for (auto i:M.rows())
            for (auto j:M.cols())
                M(i,j)=(::Grad2(es[i-1],b->es[j-1],l,l) + l*(l+1)*::Inv_r2(es[i-1],b->es[j-1],2*l))*ns[i-1]*b->ns[j-1];
    return M;
}


Slater_IBS::omlv_t Slater_IBS::Charge() const
{
    omlv_t V(size());
    for (auto i:V.indices())
            V(i)= ::Charge(es[i-1],l)*ns[i-1];

    return V;
}

dERI3 Slater_IBS::Overlap  (const Fit_IBS& _c) const
{
    const Slater_IBS& c=dynamic_cast<const Slater_IBS&>(_c);
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
dERI3 Slater_IBS::Repulsion(const Fit_IBS& _c) const
{
    const Slater_IBS& c=dynamic_cast<const Slater_IBS&>(_c);
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

template <class v> v slater(double r,size_t l,const v& e, const v& n)
{
    return n*uintpow(r,l)*exp(-e*r);
}
template <class v> v grad_slater(double r,size_t l,const v& e, const v& n)
{
    double lr= r==0 ? 0 : l/r;
    return (lr-e)*slater(r,l,e,n);
}

Slater_IBS::Vec    Slater_IBS::operator() (const RVec3& r) const
{
    return convert(slater(norm(r),l,es,ns));
}

Slater_IBS::Vec3Vec Slater_IBS::Gradient(const RVec3& r) const
{
    Vec3Vec ret(size());
    double mr=norm(r);
    if (mr==0.0)
    {
        Fill(ret,RVec3{0,0,0});
        return ret;
    }
    ds_t grad=grad_slater(norm(r),l,es,ns);
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


std::ostream&  Slater_IBS::Write(std::ostream& os) const
{
    return os << " with " << size() << " basis functions, alpha={" << es[0] << " ... " << es[size()-1] << "}" << std::endl;
}



Slater_IBS::ds_t Slater_RKBS_IBS::norms() const
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

double Slater_RKBS_IBS::Inv_r1(double ea , double eb,size_t l_total) const
{
    return ea*eb*::Slater::Integral(ea+eb,l_total-1); //Already has 4*Pi
}

Slater_IBS::ds_t Slater_RKBS_IBS::eval(const RVec3& r) const
{
    double mr=norm(r);
    ds_t f=-es;
    if (kappa >0) 
        f+=(2*kappa+1)/mr;
        
    return f*slater(mr,l,es,ns);
}

Slater_IBS::Vec    Slater_RKBS_IBS::operator() (const RVec3& r) const
{
   return convert(eval(r)); 
}

Slater_IBS::Vec3Vec Slater_RKBS_IBS::Gradient(const RVec3& r) const
{
    Vec3Vec ret(size());
    double mr=norm(r);
    if (mr==0.0)
    {
        Fill(ret,RVec3{0,0,0});
        return ret;
    }
    ds_t grad=eval(r)*(l/mr-es);
   
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
