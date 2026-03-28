// File: BasisSet/Atom/Gaussian/RKB/Imp/IBS_Evaluator.C
module;
#include <valarray>
#include <cmath>
#include <cassert>
module BasisSet.Atom.Gaussian.RKB.IBS_EValuator;
import qchem.BasisSet.Atom.GaussianIntegrals;
import Common.Constants;
import qchem.Conversions;

Gaussian_IBS::ds_t Gaussian_RKBS_IBS::norms() const
{
    return Gaussian_IBS::norms()/(2*c_light);
}

double Gaussian_RKBS_IBS::Inv_r1(double ea , double eb,size_t l_total) const
{
    return 4*ea*eb*::Gaussian::Integral(ea+eb,l_total+1); //Don't count the r^2 in dr^3
}


Gaussian_IBS::ds_t Gaussian_RKBS_IBS::eval(const rvec3_t& r) const
{
    double mr=norm(r);
    ds_t f=-2*es*mr;
    if (kappa >0) 
        f+=(2*kappa+1)/mr;
        
    return f*gaussian(mr,l,es,ns);
}

rvec_t Gaussian_RKBS_IBS::operator() (const rvec3_t& r) const
{
   return convert1(eval(r)); //valarray -> rvec_t 
}

rvec3vec_t Gaussian_RKBS_IBS::Gradient(const rvec3_t& r) const
{
    rvec3vec_t ret(size());
    double mr=norm(r);
    if (mr==0.0)
    {
        ret=rvec3_t{0,0,0};
        return ret;
    }
    ds_t grad=eval(r)*(l/mr-2*mr*es);
   
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
