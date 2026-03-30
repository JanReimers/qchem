// File: BasisSet/Atom/radial/Imp/Slater_IBS.C
module;
#include <cmath>
#include <cassert>
#include <blaze/math/DynamicVector.h>

module BasisSet.Atom.Slater.RKB.IBS_Evaluator;
import qchem.BasisSet.Atom.Slater.Integrals;
import Common.Constants;


rvec_t Slater_RKBS_IBS::norms() const
{
    return Slater_IBS::norms()/(2*c_light);
}

double Slater_RKBS_IBS::Inv_r1(double ea , double eb,size_t l_total) const
{
    return ea*eb*::Slater::Integral(ea+eb,l_total-1); //Already has 4*Pi
}

rvec_t Slater_RKBS_IBS::eval(const rvec3_t& r) const
{
    double mr=norm(r);
    rvec_t f=-es;
    if (kappa >0) 
        f+=(2*kappa+1)/mr;
        
    return f*slater(mr,l,es,ns);
}

rvec_t Slater_RKBS_IBS::operator() (const rvec3_t& r) const
{
   return eval(r); 
}

rvec3vec_t Slater_RKBS_IBS::Gradient(const rvec3_t& r) const
{
    rvec3vec_t ret(size());
    double mr=norm(r);
    if (mr==0.0)
    {
        ret=rvec3_t{0,0,0};
        return ret;
    }
    rvec_t grad=eval(r)*(l/mr-es);
   
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
