// File: BasisSet/Atom/Gaussian/NR/Imp/IBS_Evaluator.C
module;
#include <cmath>
#include <cassert>
#include <iostream>
#include <blaze/math/DynamicVector.h>

module qchem.BasisSet1.Atom.Evaluators.Gaussian.IBS; 
import qchem.BasisSet1.Atom.Evaluators.Gaussian.Internal.Rk; 
import qchem.BasisSet1.Atom.Evaluators.Gaussian.Internal.GaussianIntegrals; 
import qchem.BasisSet1.Atom.Evaluators.Gaussian.Internal.ExponentScaler; 
import Common.Constants;



inline double Overlap(double ea , double eb,size_t l_total)
{
    return Gaussian::Integral(ea+eb,l_total); //Already has 4*Pi and r^2 from dr.
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
std::string Gaussian_IBS_Evaluator::Name() const
{
    return "Spherical Gaussian ";
}


double Gaussian_IBS_Evaluator::Grad2(double ea , double eb,size_t la, size_t lb)
{
    assert(la==lb);
    double t=ea+eb;
    size_t l1=la+1;
    return  l1*l1     * Gaussian::Integral(t,2*la-2)
            -2*l1 * t * Gaussian::Integral(t,2*la  )
            +4*ea*eb  * Gaussian::Integral(t,2*la+2);
}

double Gaussian_IBS_Evaluator::Inv_r2(double ea , double eb,size_t l_total)
{
    return Gaussian::Integral(ea+eb,l_total-2); //Already has 4*Pi
}

// This need overridability.
double Gaussian_IBS_Evaluator::Inv_r1(double ea , double eb,size_t l_total) const
{
    return Gaussian::Integral(ea+eb,l_total-1); //Already has 4*Pi
}

rvec_t Gaussian_IBS_Evaluator::exponents(size_t N, double emin, double emax, const Irrep_QNs::sym_t& ir)
{
    size_t LMax=3; //TODO how do we get the real LMax(Z) into this?
    ::Gaussian::ExponentScaler ss(N,emin,emax,LMax);
    return ss.Get_es(ir);
}

rvec_t Gaussian_IBS_Evaluator::norms() const
{
    size_t N=es.size();    
    rvec_t ret(N);
    for (size_t i=0;i<N;i++) ret[i]=1.0/sqrt(::Overlap(es[i],es[i],2*l)); 
    return ret;
}


rsmat_t Gaussian_IBS_Evaluator::Repulsion() const
{
    size_t N=size();
    rsmat_t S(N);
    for (auto i:iv_t(0,N))
        for (auto j:iv_t(i,N))
            S(i,j)= ::Repulsion(es[i],es[j],l,l)*ns[i]*ns[j];

    return S;
}

rvec_t Gaussian_IBS_Evaluator::Charge() const
{
    rvec_t V(size());
    for (auto i:iv_t(0,size()))
            V[i]= ::Charge(es[i],l)*ns[i];

    return V;
}

rmat_t Gaussian_IBS_Evaluator::XRepulsion(const IBS_Evaluator& _b) const
{
    const Gaussian_IBS_Evaluator& b=dynamic_cast<const Gaussian_IBS_Evaluator&>(_b);
    size_t Nr=size(), Nc=b.size();
    rmat_t M(Nr,Nc);
    for (auto i:iv_t(0,Nr))
            for (auto j:iv_t(0,Nc))
                M(i,j)=::Repulsion(es[i],b.es[j],l,b.l)*ns[i]*b.ns[j];
    return M;
}

rmat_t Gaussian_IBS_Evaluator::XKinetic(const IBS_Evaluator& _b) const
{
    const Gaussian_IBS_Evaluator& b=dynamic_cast<const Gaussian_IBS_Evaluator&>(_b);
    assert(l==b.l);
    size_t Nr=size(), Nc=b.size();
    rmat_t M(Nr,Nc);
    for (auto i:iv_t(0,Nr))
            for (auto j:iv_t(0,Nc))
                M(i,j)=(Grad2(es[i],b.es[j],l,l) + l*(l+1)*Inv_r2(es[i],b.es[j],2*l))*ns[i]*b.ns[j];
    return M;
}



rvec_t Gaussian_IBS_Evaluator::operator() (const rvec3_t& r) const
{
    return gaussian(norm(r),l,es,ns); 
}

rvec3vec_t Gaussian_IBS_Evaluator::Gradient(const rvec3_t& r) const
{
    rvec3vec_t ret(size());
    double mr=norm(r);
    if (mr==0.0)
    {
        ret=rvec3_t{0,0,0};
        return ret;
    }
    rvec_t grad=grad_gaussian(norm(r),l,es,ns);
    rvec3_t rhat=r/norm(r);
    size_t i=0;
    for (auto& g:grad) ret[i++]=g*rhat;
    return ret;
}

std::ostream&  Gaussian_IBS_Evaluator::Write(std::ostream& os) const
{
    return os << " with N=" << es.size() << " basis functions, alpha={" << es[0] << " ... " << es[size()-1] << "}" << std::endl;
}


