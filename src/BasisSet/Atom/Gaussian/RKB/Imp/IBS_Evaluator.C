// File: BasisSet/Atom/Gaussian/RKB/Imp/IBS_Evaluator.C
module;
#include <valarray>
#include <cmath>
#include <cassert>
module BasisSet.Atom.Gaussian.RKB.IBS_EValuator;
import qchem.BasisSet.Atom.GaussianIntegrals;

Gaussian_IBS::ds_t Gaussian_RKBS_IBS::norms() const
{
    size_t N=es.size();
    ds_t ret(N);
    for (size_t i=0;i<N;i++) 
    {
        double k=Grad2(es[i],es[i],l,l) + l*(l+1)*Inv_r2(es[i],es[i],2*l);
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
