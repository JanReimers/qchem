// File: BasisSet/Atom/radial/Imp/Slater_IBS.C
module;
#include <cmath>
#include <cassert>
#include <iostream>
#include <blaze/math/DynamicVector.h>

module BasisSet.Atom.Slater.NR.IBS_Evaluator;
import qchem.BasisSet.Atom.Slater.Rk;
import qchem.BasisSet.Atom.Slater.Integrals;
import Common.Constants;



inline double Overlap(double ea , double eb,size_t l_total)
{
    return Slater::Integral(ea+eb,l_total); //Already has 4*Pi and r^2 from dr.
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
double Slater_IBS::Grad2(double ea , double eb,size_t la, size_t lb)
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

double Slater_IBS::Inv_r2(double ea , double eb,size_t l_total)
{
    return Slater::Integral(ea+eb,l_total-2); //Already has 4*Pi
}

double Slater_IBS::Inv_r1(double ea , double eb,size_t l_total) const
{
    return Slater::Integral(ea+eb,l_total-1); //Already has 4*Pi
}

rvec_t Slater_IBS::norms() const
{
    size_t N=es.size();    
    rvec_t ret(N);
    for (size_t i=0;i<N;i++) ret[i]=1.0/sqrt(::Overlap(es[i],es[i],2*l)); 
    return ret;
}

rsmat_t Slater_IBS::Overlap() const
{
    size_t N=size();
    rsmat_t S(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
            S(i,j)= ::Overlap(es[i],es[j],2*l)*ns[i]*ns[j];

    return S;
}

rsmat_t Slater_IBS::Grad2() const
{
    size_t N=size();
    rsmat_t S(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
            S(i,j)= Grad2(es[i],es[j],l,l)*ns[i]*ns[j];

    return S;
}

rsmat_t Slater_IBS::Inv_r1() const
{
    size_t N=size();
    rsmat_t S(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
            S(i,j)= Inv_r1(es[i],es[j],2*l)*ns[i]*ns[j];

    return S;
}

rsmat_t Slater_IBS::Inv_r2() const
{
    size_t N=size();
    rsmat_t S(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
            S(i,j)= Inv_r2(es[i],es[j],2*l)*ns[i]*ns[j];

    return S;
}

rsmat_t Slater_IBS::Repulsion() const
{
    size_t N=size();
    rsmat_t S(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
            S(i,j)= ::Repulsion(es[i],es[j],l,l)*ns[i]*ns[j];

    return S;
}

rmat_t Slater_IBS::XRepulsion(const Fit_IBS& _b) const
{
    const Slater_IBS& b=dynamic_cast<const Slater_IBS&>(_b);
    size_t Nr=size(), Nc=b.size();
    rmat_t M(Nr,Nc);
    for (auto i:iv_t(0,Nr))
            for (auto j:iv_t(0,Nc))
                M(i,j)=::Repulsion(es[i],b.es[j],l,b.l)*ns[i]*b.ns[j];
    return M;
}

rmat_t Slater_IBS::XKinetic(const Orbital_RKBS_IBS<double>* _b) const
{
    const Slater_IBS* b=dynamic_cast<const Slater_IBS*>(_b);
    assert(b);
    assert(l==b->l);
    size_t Nr=size(), Nc=b->size();
    rmat_t M(Nr,Nc);
    for (auto i:iv_t(0,Nr))
            for (auto j:iv_t(0,Nc))
                M(i,j)=(Grad2(es[i],b->es[j],l,l) + l*(l+1)*Inv_r2(es[i],b->es[j],2*l))*ns[i]*b->ns[j];
    return M;
}


rvec_t Slater_IBS::Charge() const
{
    rvec_t V(size());
    for (auto i:iv_t(0,size()))
            V[i]= ::Charge(es[i],l)*ns[i];

    return V;
}

dERI3 Slater_IBS::Overlap  (const Fit_IBS& _c) const
{
    const Slater_IBS& c=dynamic_cast<const Slater_IBS&>(_c);
    dERI3 S3;
    size_t N=size();
    for (size_t ic=0;ic<c.size();ic++) 
    {
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)=::Overlap(es[i]+es[j],c.es[ic],l+l+c.l)*ns[i]*ns[j]*c.ns[ic];  
        
        S3.push_back(S);
    }
    return S3;
}
dERI3 Slater_IBS::Repulsion(const Fit_IBS& _c) const
{
    const Slater_IBS& c=dynamic_cast<const Slater_IBS&>(_c);
    dERI3 S3;
    size_t N=size();
    for (size_t ic=0;ic<c.size();ic++) 
    {
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)=::Repulsion(es[i]+es[j],c.es[ic],l,c.l)*ns[i]*ns[j]*c.ns[ic];  
        
        S3.push_back(S);
    }
    return S3;
}



rvec_t Slater_IBS::operator() (const rvec3_t& r) const
{
    return slater(norm(r),l,es,ns);
}

rvec3vec_t Slater_IBS::Gradient(const rvec3_t& r) const
{
    rvec3vec_t ret(size());
    double mr=norm(r);
    if (mr==0.0)
    {
        ret=rvec3_t{0,0,0};
        return ret;
    }
    rvec_t grad=grad_slater(norm(r),l,es,ns);
    // std::cout << "grad=";
    // for (auto g:grad) std::cout << g << " ";
    // std::cout << std::endl;
    rvec3_t rhat=r/norm(r);
    size_t i=0;
    for (auto& g:grad) ret[i++]=g*rhat;
    // std::cout << "ret=";
    // for (auto g:ret) std::cout << g << " ";
    // std::cout << std::endl;
    return ret;
}


std::ostream&  Slater_IBS::Write(std::ostream& os) const
{
    return os << " with " << size() << " basis functions, alpha={" << es[0] << " ... " << es[size()-1] << "}" << std::endl;
}
